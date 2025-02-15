import pandas as pd
import numpy as np
import time
import pandas as pd
import os
from tqdm import tqdm
import re
import sys


# Caclulating CIGAR 
# From query to target

def parse_cigar(cigar):
    """Parses a CIGAR string into a list of (operation, length) tuples."""
    return [(m[-1], int(m[:-1])) for m in re.findall(r'\d+[MIDNSHP=X]', cigar)]

def query_to_target_position(query_pos, query_start, target_start, cigar):
    """Converts a query position to a target position based on the PAF CIGAR string."""
    if cigar is None:
        raise ValueError("CIGAR string not found in PAF entry")
    
    cigar_ops = parse_cigar(cigar)
    q_pos = query_start
    t_pos = target_start
    
    for op, length in cigar_ops:
        if op in 'MD=X':  # Query consumes sequence
            if q_pos + length > query_pos:
                if op in 'M=X':  # Match or mismatch
                    return t_pos + (query_pos - q_pos)
                elif op == 'D':  # Deletion in query (gap in target)
                    return t_pos  # Stays at the same target position
            q_pos += length
        if op in 'M=D=X':  # Target consumes sequence
            t_pos += length
    
    raise ValueError("Query position out of range")

# READ inputs 
# READ GAF 
def readGAF(gaf_file, region):
    # print("read gaf file...")
    gaf = pd.read_csv(gaf_file, sep="\t", header=None)
    gaf.columns = ["qname", "qlen", "qstart", "qend", "strand", "path", "path_len", "path_start", "path_end", 
                "num_match", "alignment_block", "mapq", "NM", "AS", "dv", "id", "cg"]

    # region
    # print(f"Finding node alignment for {region}")
    parts = re.split(r"[:\-]", region)
    chr = parts[0]
    start, end = map(int, parts[1:3])

    # Ensure correct filtering for overlapping alignments
    gaf_chr_db = gaf[(gaf['qname'] == chr) & ((gaf['qstart'] <= end) & (gaf['qend'] >= start))]
    gaf_chr_db.loc[:,'path_split'] = gaf_chr_db['path'].str.split(r'[><]').apply(lambda x: [i for i in x if i])
    # gaf_chr_db = pd.DataFrame(gaf_chr_db)
    return gaf_chr_db, chr, start, end

def readGraph(graph_file):
    # Read and parse graph data
    # print("read graph file...")
    graph = pd.read_csv(graph_file, sep="\t", header=None)

    # Extract segment and link information
    segment = graph.loc[graph[0] == "S", [0, 1, 3]]
    link = graph[graph[0] == "L"]

    # Rename columns
    segment.columns = ["s", "node", "len"]
    link.columns = ["l", "from", "fromOrient", "to", "toOrient", "overlap"]

    # Extract segment lengths safely
    segment['len'] = segment['len'].str.extract(r'LN:i:(\d+)')[0]
    segment['len'] = pd.to_numeric(segment['len'], errors='coerce').fillna(0).astype(int)

    # Extract overlap values safely
    link['overlap'] = link['overlap'].str.extract(r'(\d+)M$')[0]
    link['overlap'] = pd.to_numeric(link['overlap'], errors='coerce').fillna(0).astype(int)
    return segment, link


def finding_nodes(gaf_chr_db, start, end,segment, link):
    final_nodes = []
    total_bed_start = None
    total_bed_end = None
    node_space = None
      # Initialize to avoid UnboundLocalError
    # check the dtype of get_chr_db, if it is tuple, then print "it is tuple"
    num_path = gaf_chr_db.shape[0]
    # print(f"{num_path} paths found for the region")
    
    for dimension in range(gaf_chr_db.shape[0]):
        gaf_chr = gaf_chr_db.iloc[dimension, :].copy()  # Ensure it's a copy, not a view
        
        # Compute normalized start and end positions
        if gaf_chr['qstart'] < start:
            total_bed_start = query_to_target_position(start, gaf_chr['qstart'], gaf_chr['path_start'], gaf_chr['cg'])
        else:
            total_bed_start = 1  # Start from the beginning if within range
        
        if gaf_chr['qend'] > end:
            total_bed_end = query_to_target_position(end, gaf_chr['qstart'], gaf_chr['path_start'], gaf_chr['cg'])
        else: 
            total_bed_end = gaf_chr['path_len']  # End at the last position if within range

        if not isinstance(gaf_chr['path_split'], list) or len(gaf_chr['path_split']) == 0:
            print("No node alignment found for the region")
            continue  
        
        elif len(gaf_chr['path_split']) == 1:
            final_nodes.append(gaf_chr['path_split'][0])
        
        else:
            # Create a DataFrame for node space
            idx_start = 0
            node_space = pd.DataFrame()

            for i in range(1, len(gaf_chr['path_split'])):
                pre = gaf_chr['path_split'][i-1]
                suf = gaf_chr['path_split'][i]

                # Extract values as native Python numbers
                overlapM = link.loc[(link['from'] == pre) & (link['to'] == suf), 'overlap']
                pre_len = segment.loc[segment['node'] == pre, 'len']
                suf_len = segment.loc[segment['node'] == suf, 'len']

                # Convert Pandas Series to a single value (float -> int)
                overlapM = int(overlapM.iloc[0]) if not overlapM.empty else 0
                pre_len = int(pre_len.iloc[0]) if not pre_len.empty else 0
                suf_len = int(suf_len.iloc[0]) if not suf_len.empty else 0

                # print(f"overlapM: {overlapM}, pre_len: {pre_len}, suf_len: {suf_len}")

                # Append row to node_space DataFrame
                node_space = pd.concat([
                    node_space,
                    pd.DataFrame({'node': [pre], 'start_coor': [idx_start], 'end_coor': [idx_start + pre_len]})
                ], ignore_index=True)

                # Update idx_start
                idx_start = idx_start + pre_len - overlapM

            # Append last node
            node_space = pd.concat([
                node_space,
                pd.DataFrame({'node': [suf], 'start_coor': [idx_start], 'end_coor': [idx_start + suf_len]})
            ], ignore_index=True)
            # Find nodes that overlap with the region
            nodes = node_space[(node_space['end_coor'] >=total_bed_start) & (node_space['start_coor'] <= total_bed_end)]
            # nodes = node_space[(node_space['end_coor'].between(total_bed_start, total_bed_end) & (node_space['start_coor'].between(total_bed_start, total_bed_end)))]
            final_nodes.extend(nodes['node'].tolist())
    
    # Flatten and remove duplicates
    final_nodes = list(set(final_nodes))
    return final_nodes, total_bed_start, total_bed_end, node_space

def getNodes_from_unHPCregion(gaf_file, graph_file, regions_list):
    """
    getNodes_from_unHPCregion reads a GAF file, a graph file, and a list of regions, and returns a DataFrame with the nodes that overlap with the regions.
    
    Parameters
    ----------
    gaf_file : str
        The path to the GAF file.
    graph_file : str
        The path to the graph file.
    regions_list : list
        A list of regions in the format "chr:start-end".
    
    Returns
    -------
    regions_node_db : DataFrame
        A DataFrame with the regions and the nodes that overlap with them.
    """
    regions_node_db = pd.DataFrame()
    for i in tqdm(range(len(regions_list)), desc="Finding nodes for regions"): 
        region = regions_list[i]
        # print(f"Finding nodes for {region}")
        gaf_chr_db, chr, start, end = readGAF(gaf_file, region)
        segment, link = readGraph(graph_file)
        final_nodes, total_bed_start, total_bed_end, node_space = finding_nodes(gaf_chr_db, start, end,segment, link)
        regions_node_db = pd.concat([regions_node_db, pd.DataFrame({"region": [region], "nodes": [final_nodes]})], ignore_index=True)
    return regions_node_db

def bed_to_regionsList(bed_file):
    """
    bed_to_regionsList reads a BED file and returns a list of regions in the format "chr:start-end".

    Parameters
    ----------
    bed_file : str
        The path to the BED file.

    Returns
    -------
    regions_list : list
        A list of regions in the format "chr:start-end".
    """
    bed = pd.read_csv(bed_file, sep="\t", header=None)
    bed.columns = ["chrom", "start", "end"]
    regions_list = list(bed['chrom'] + ":" + bed['start'].astype(str) + "-" + bed['end'].astype(str))
    return regions_list
