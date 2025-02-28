import pandas as pd
import numpy as np
import time
from tqdm import tqdm
import pandas as pd
import os

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
def readGAF_extractRegion(gaf_file, region):
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
            print(f"final_nodes: {final_nodes}")
            print(gaf_chr)
            # Create a DataFrame for node space
            node_space = pd.concat([
                node_space,
                pd.DataFrame({'node': gaf_chr_db['path_split'].values[0], 
                              'start_coor': 1, 'end_coor': gaf_chr['path_len']})
            ], ignore_index=True)
        
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
        A list of regions in the format "chr:start-end". compressed coordinates
    
    Returns
    -------
    regions_node_db : DataFrame
        A DataFrame with the regions and the nodes that overlap with them.
    """
    regions_node_db = pd.DataFrame()
    regions_node_coor = pd.DataFrame()
   
    for i in tqdm(range(len(regions_list)), desc="Finding nodes for regions"): 
        region = regions_list[i]
        # print(f"Finding nodes for {region}")
        gaf_chr_db, chr, start, end = readGAF_extractRegion(gaf_file, region)
        
        segment, link = readGraph(graph_file)
        
        final_nodes, total_bed_start, total_bed_end, node_space = finding_nodes(gaf_chr_db, start, end,segment, link)
        # print(f"final_nodes: {final_nodes}")
        # print(f"total_bed_start: {total_bed_start}")
        # print(f"total_bed_end: {total_bed_end}")    
        # print(f"node_space: {node_space}")
        # get the node bed
        # sub_db = node_space.loc[(node_space['start_coor'].between(total_bed_start, total_bed_end) | (node_space['end_coor'].between(total_bed_start,total_bed_end))),]
        sub_db = node_space[(node_space['end_coor'] >=total_bed_start) & (node_space['start_coor'] <= total_bed_end)]
        sub_db['len_node'] = sub_db['end_coor'] - sub_db['start_coor']
        sub_db['start_coor_on_node'] = total_bed_start - sub_db['start_coor']
        sub_db.loc[sub_db['start_coor_on_node'] < 0, 'start_coor_on_node'] = 0

        sub_db['end_coor_on_node'] = total_bed_end - sub_db['end_coor']
        idx = sub_db['end_coor_on_node'] < 0
        sub_db.loc[idx, 'end_coor_on_node'] = sub_db.loc[idx, 'len_node']


        sub_db['region'] = region
        sub_db.columns= ['node', 'start_coor_path_comp', 'end_coor_path_comp', 'len_node', 'start_coor_on_node','end_coor_on_node','region']
        regions_node_coor = pd.concat([regions_node_coor, sub_db], ignore_index=True)
        regions_node_db = pd.concat([regions_node_db, pd.DataFrame({"region": [region], "nodes": [final_nodes]})], ignore_index=True)
    
    return regions_node_db, regions_node_coor

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


def read_untig_Scfmap(file_path = "6-layoutContigs/unitig-popped.layout.scfmap"):
    """
    read_untig_Scfmap reads a unitig scfmap file and returns a DataFrame with the contig, unitig, and piece information.

    Parameters
    ----------
    file_path : str
        The path to the unitig scfmap file.

    Returns
    -------
    scfmap : DataFrame
        A DataFrame with the contig, unitig, and piece information.
    """

    data = []
    with open(file_path, "r") as f:
        lines = f.readlines()
        
        for i in range(0, len(lines), 3):  # Process every 3 lines as a block
            path_line = lines[i].strip().split()  # Extract path details
            piece_line = lines[i+1].strip()  # Extract piece
            # Extract contig and utig values
            contig = path_line[1]
            utig = path_line[2]
            piece = piece_line
            
            data.append([contig, utig, piece])

    # Create DataFrame
    scfmap = pd.DataFrame(data, columns=["contig", "node", "piece"])

    # Display result
    return scfmap


def read_hapAssignRead(file):

    hapmer = pd.read_csv(file, header=None, sep ='\t')
    hapmer.columns = ['read','len_read','hap1_total_kmer','hap1_found_kmer','hap2_total_kmer','hap2_found_kmer','hap_read']
    hapmer.drop_duplicates(inplace=True)
    return hapmer


def readNodeInfo(csv = 'assembly.colors.csv',
                 graph = 'assembly.homopolymer-compressed.noseq.gfa',
                 scfamp_node = "6-layoutContigs/unitig-popped.layout.scfmap",
                 layout_node = '6-layoutContigs/unitig-popped.layout'):
    """
    readNodeInfo reads a CSV file, a graph file, and a scfmap file, and returns a DataFrame with the node information.

    Parameters
    ----------
    csv : str
        The path to the CSV file.
    graph : str
        The path to the graph file.
    scfamp_node : str
        The path to the scfmap file.

    Returns
    -------
    nodeinfo : DataFrame
        A DataFrame with the node information.
    """
    
    scfmap = read_untig_Scfmap(file_path = scfamp_node)
    del scfmap['contig']

    nodeinfo = pd.read_csv(csv, header = 0, sep = '\t')
    color_dict = dict({'#8888FF' : 'pat',
                    '#AAAAAA' : 'lowcov',
                    '#FFFF00' : 'ambiguous',
                    '#FF8888' : 'mat'
                    })
    nodeinfo['hap_node'] = nodeinfo['color'].map(color_dict)

    segment, link = readGraph(graph)

    del color_dict
    del segment['s']

    nodeinfo = nodeinfo.merge(segment, on = 'node', how = 'outer')
    nodeinfo = nodeinfo.merge(scfmap, on = 'node', how = 'outer')

    # read node layout 
    node_layout = pd.read_csv(layout_node, header=None)
    node_layout.columns = ['layout_info']
    node_layout_idx = node_layout.loc[node_layout['layout_info'].fillna('').str.startswith('tig'), 'layout_info'].str.split(r'\s+', expand=True)

    node_layout_idx.columns = ['tig','piece']
    del node_layout_idx['tig']
    node_layout_idx = node_layout_idx.reset_index()
    # del node_layout['index']

    return nodeinfo, node_layout, node_layout_idx


def get_hap_ratio(obj, node, hifi_read, node_layout, node_layout_idx, nodeinfo, regions_node_coor):
    color = nodeinfo.copy()
    contigName = obj.paths[obj.paths['path'].str.contains(node+'+') | obj.paths['path'].str.contains(node+'-')]['name'].values[0]
    # len_node = regions_node_coor.loc[regions_node_coor['node'] == node, 'len_node'].values[0]
    start_coor = regions_node_coor.loc[regions_node_coor['node'] == node, 'start_coor_on_node'].values[0]
    end_coor = regions_node_coor.loc[regions_node_coor['node'] == node, 'end_coor_on_node'].values[0]

    hap_ratio = color.loc[color['node'] == node,]['mat:pat'].values[0]

    hap = color.loc[color['node'] == node,]['hap_node'].values[0]

    piece = nodeinfo.loc[nodeinfo['node'] == node,]['piece'].values[0]
    idx = node_layout_idx.loc[node_layout_idx['piece'] == piece,'index'].values[0] + 4
    idx_next = node_layout_idx.loc[node_layout_idx['piece'] == piece].index[0] + 1
    idx_next = node_layout_idx.loc[idx_next,'index'] - 2

    layout_sub = node_layout.loc[idx:idx_next,:]['layout_info'].str.split(r'\t', expand=True)
    layout_sub.columns = ['read','start_on_node','end_on_node']

    mergedb = hifi_read.merge(layout_sub, on='read', how='right')
    # mergedb.loc[mergedb['hap'].isna(), 'hap'] = 'unassigned'

    mergedb['platform']='ont'
    mergedb.loc[mergedb['read'].str.endswith('ccs'), 'platform'] = 'ccs'

    mergedb['start_on_node'] = mergedb['start_on_node'].astype(int)
    mergedb = mergedb.sort_values(by=['start_on_node']).reset_index(drop=True)

    # print(f"Node: {node}")
    # print(f"hap: {hap}")
    # print(f"Piece: {piece}")
    # print(f"Start on utig: {start_coor}")
    # print(f"ENd on utig: {end_coor}")

    return mergedb, hap, contigName, start_coor, end_coor, hap_ratio


def getNodeCoor(obj, regions_node_db, hifi_read, node_layout, node_layout_idx, nodeinfo, regions_node_coor):
    loc_on_node = pd.DataFrame()
    mergedb_all = pd.DataFrame()

    for i in tqdm(range(len(regions_node_db)), desc="regions of interest"):
        region = regions_node_db['region'].values[i]
        for node_num in range(len(regions_node_db.loc[i, 'nodes'])):
            node = regions_node_db.loc[i, 'nodes'][node_num]
            mergedb, hap, contigName, start_coor, end_coor, hap_ratio = get_hap_ratio(obj, node, hifi_read, node_layout, node_layout_idx, nodeinfo, regions_node_coor)
            mergedb['node'] = node
            mergedb['contig'] = contigName
            mergedb['hap_node'] = hap
            loc_on_node = pd.concat([loc_on_node, pd.DataFrame(dict(node=[node], region = [region], contig=[contigName], start=[start_coor], end=[end_coor], hap_node=[hap]))])
            mergedb_all = pd.concat([mergedb_all, mergedb])

            
    mergedb_all['start_on_node'] = mergedb_all['start_on_node'].astype(int)
    mergedb_all['end_on_node'] = mergedb_all['end_on_node'].astype(int)
    mergedb_all['mid_on_node'] = (mergedb_all['end_on_node'].astype(int) - mergedb_all['start_on_node'].astype(int))/2 + mergedb_all['start_on_node'].astype(int)
    mergedb_all['mid_on_node'] = mergedb_all['mid_on_node'].astype(int)

    return loc_on_node, mergedb_all


def nodeExtract(node):
    if node.endswith("-"):
        strand = '-'
        node = node.rstrip('-')
    elif node.endswith("+"):
        strand = '+'
        node = node.rstrip('+')
    else:
        strand = '.'
    return node, strand


def update_seg_link_withGap(gap_list, segment, link):
    segment_sub = segment.copy()
    link_sub = link.copy()

    for z in range(len(gap_list)):
        gap = gap_list[z]
        gap_pre , gap_pre_strand= nodeExtract(gap[0])
        gap_sur , gap_sur_strand = nodeExtract(gap[2])
        gap_gapName = gap[1]
        gap_len = gap_gapName.split('N')[1]

        segment_sub = pd.concat([segment_sub, pd.DataFrame({'s' : ['s'],'node': [gap_gapName], 'len': [gap_len]})])
        link_sub = pd.concat([link_sub, pd.DataFrame({'l' : ['L'],'from': [gap_pre], 'fromOrient' : [gap_pre_strand], 'to': [gap_gapName], 'toOrient':[gap_sur_strand], 'overlap': [0]})])
        link_sub = pd.concat([link_sub, pd.DataFrame({'l' : ['L'],'from': [gap_gapName], 'fromOrient' : [gap_pre_strand], 'to': [gap_sur], 'toOrient':[gap_sur_strand], 'overlap': [0]})])

    return segment_sub, link_sub


def getNodeSpace_from_onePath(node_list, segment, link):
    idx_start = 0
    node_space = pd.DataFrame()

    if len(node_list) == 1:
        node, strand = nodeExtract(node_list[0])
        node_len = segment.loc[segment['node'] == node, 'len']
        node_len = int(node_len.iloc[0]) if not node_len.empty else 0

        node_space = pd.concat([
            node_space,
            pd.DataFrame({'node': [node], 'start_coor': [idx_start], 'end_coor': [idx_start + node_len], 'strand': [strand]})
        ], ignore_index=True)

        return node_space

    for i in range(1,len(node_list)):
        pre , pre_strand= nodeExtract(node_list[i-1])
        suf, suf_strand = nodeExtract(node_list[i])
        

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
            pd.DataFrame({'node': [pre], 'start_coor': [idx_start], 'end_coor': [idx_start + pre_len], 'strand': [pre_strand]})
        ], ignore_index=True)

        # Update idx_start
        idx_start = idx_start + pre_len - overlapM

    # Append last node
    node_space = pd.concat([
        node_space,
        pd.DataFrame({'node': [suf], 'start_coor': [idx_start], 'end_coor': [idx_start + suf_len], 'strand': [suf_strand]})
    ], ignore_index=True)

    return node_space


def getNodeSpace_from_allPath(obj, segment, link, file):
    node_space = pd.DataFrame()

    for j in tqdm(range(len(obj.paths))):
        node_list = obj.paths['path'].values[j].split(',')
        contig_name = obj.paths['name'].values[j]
        if '[' in obj.paths['path'].values[j]:
            gap_list = list(obj.gaps[obj.gaps['name'] == contig_name]['gaps'])
            segment, link = update_seg_link_withGap(gap_list, segment, link)
            
        node_space_sub = getNodeSpace_from_onePath(node_list, segment, link)
        node_space_sub['chr'] = contig_name
        node_space_sub['score'] = 100
        node_space_sub = node_space_sub[['chr', 'start_coor','end_coor','node','score','strand']]
        node_space = pd.concat([node_space, node_space_sub])
        
    node_space.columns = ['chrom','chromStart','chromEnd','name','score','strand']
    
    node_space.to_csv(file, index=False, sep = '\t', header = True)
    
    return node_space