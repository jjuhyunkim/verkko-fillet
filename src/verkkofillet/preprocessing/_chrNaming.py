import re
import networkx as nx
from collections import Counter
import pandas as pd
import sys
import shlex
import subprocess
import os
import shutil

def make_unique(column):
    counts = {}
    result = []
    for name in column:
        if name in counts:
            counts[name] += 1
            result.append(f"{name}_{counts[name]}")
        else:
            counts[name] = 0
            result.append(name)
    return result
    
def create_pairs(input_str):
    # Split the string into a list
    split_values = input_str.split(',')
    # Generate pairs of consecutive elements
    return [(split_values[i], split_values[i+1]) for i in range(len(split_values)-1)]

def remove_elements_starting_with_bracket(lst):
    return [item for item in lst if not item.startswith('[')]

def remove_ignore_nodes(lst, ignore_lst):
    return [item for item in lst if item not in ignore_lst]

def find_multi_used_node(obj):
    path = obj.paths.copy()
    scfmap = obj.scfmap.copy()
    stats = obj.stats.copy()
    
    path['nodeSet'] = path['path'].apply(
        lambda x: set(word.rstrip('-+') for word in x.replace(',', ' ').split())
    )
    path = pd.merge(scfmap, path, how = 'outer',left_on = "pathName",right_on="name")
    path['nodeSet'] = path['nodeSet'].apply(remove_elements_starting_with_bracket)
    
    assignedContig = list(obj.stats['contig']) 
    unassignedPath = path.loc[~path['contig'].isin(assignedContig)].reset_index()
    assignedPath = path.loc[path['contig'].isin(assignedContig)].reset_index()
    
    assignedContig = list(stats['contig'])
    assignedPath = pd.merge(stats, assignedPath, on = 'contig', how = 'outer')
    
    path_grouped = assignedPath.groupby('ref_chr')['nodeSet'].agg(
        lambda x: set([item for sublist in x if isinstance(sublist, (list, set)) for item in sublist])
    ).reset_index()
    
    list_of_lists = path_grouped['nodeSet'].tolist()
    
    # Example list of lists
    
    # Flatten the list of lists and count the occurrences of each element
    flat_list = [item for sublist in list_of_lists for item in sublist]
    element_counts = Counter(flat_list)
    
    # Find the elements that appear in more than one list
    duplicates = [item for item, count in element_counts.items() if count > 1]
    
    return duplicates, path_grouped

def naming_contigs(obj, node_database, duplicate_nodes , 
                   gfa = "assembly.homopolymer-compressed.noseq.gfa", dam = "mat", sire = "pat", fai = "assembly.fasta.fai"):
    stats = obj.stats.copy()
    path = obj.paths.copy()
    scfmap = obj.scfmap.copy()
    fai_chr = pd.read_csv(fai, sep='\t', header=None, usecols=[0])[0].tolist()

    scfmap = scfmap.loc[scfmap['contig'].isin(fai_chr)]
    
    gfa_link = pd.read_csv(gfa, sep = '\t', comment = "S", header = None, usecols= [1,3])
    gfa_link.columns = ['start','end']
    gfa_link = gfa_link[~gfa_link['start'].isin(duplicate_nodes)]
    gfa_link = gfa_link[~gfa_link['end'].isin(duplicate_nodes)]
    
    path['nodeSet'] = path['path'].apply(
        lambda x: set(word.rstrip('-+') for word in x.replace(',', ' ').split())
    )
    path['path'] = path['nodeSet'].apply(lambda x: ','.join(x))
    path = pd.merge(scfmap, path, how = 'left',left_on = "pathName",right_on="name")
    path['nodeSet'] = path['nodeSet'].apply(remove_elements_starting_with_bracket)
    path['nodeSet'] = path['nodeSet'].apply(lambda lst: remove_ignore_nodes(lst, duplicate_nodes))
    
    assignedContig = list(obj.stats['contig'])
    assignedContig = [re.sub(r':.*', '', string) for string in assignedContig]
    unassignedPath = path.loc[~path['contig'].isin(assignedContig)].reset_index()
    assignedPath = path.loc[path['contig'].isin(assignedContig)].reset_index()
    
    # Exploding the list into separate rows
    df_exploded = node_database.explode('nodeSet').reset_index(drop=True)
    
    # Renaming columns to match the desired output
    df_exploded.columns = ['chr', 'node']

    result_df = gfa_link.copy()

    # Create a graph
    G = nx.Graph()
    G.add_edges_from(result_df.values)
    
    # Find connected components
    connected_components = [set(component) for component in nx.connected_components(G)]
    len(connected_components)
    
    # assign chromosome
    dict1 = {}
    
    for i in range(0,len(connected_components)):
        dict2 = {}
        chr_assign = df_exploded.loc[df_exploded['node'].isin(connected_components[i]), "chr"].unique()
        if len(chr_assign) == 1:
            dict2 = {chr_assign[0] : connected_components[i]}
            dict1.update(dict2)
            #print("component_" + str(i) + " : " + chr_assign.astype(str))
        if len(chr_assign) >1: 
            chrName = "_".join(chr_assign)
            dict2 = {chrName : connected_components[i]}
            dict1.update(dict2)
            #print("component_" + str(i) + " : " + chr_assign)
        if len(chr_assign) < 1:
            print("component_" + str(i) + " : empty")
    
    # Assuming dict1 has sets or lists as values and we want to check if 'nodeSet' is a subset of any of those sets
    for i in range(0, unassignedPath.shape[0]):
        # Access the nodeSet for the current row and convert it to a set
        node_set = set(unassignedPath.loc[i, "nodeSet"])
        
        # Find the key that corresponds to a matching value in dict1
        some_key = None  # Initialize the key
        
        for key, value in dict1.items():
            if node_set.issubset(set(value)):  # Check if node_set is a subset of value
                some_key = key  # If a match is found, assign the key
                break  # No need to continue checking after the first match
        
        # If a match was found, assign the key to the 'assignChr' column
        if some_key is not None:
            unassignedPath.loc[i, "assignChr"] = some_key
        
        # Print if it's a subset and which key was assigned
        # print(f"Row {i}: Is subset? {some_key is not None}, Assigned Key: {some_key}")
    
    # update unassigned
    unassignedPath.loc[unassignedPath['assignChr'].isna(), 'assignChr'] = "chrUn"
    unassignedPath = unassignedPath[~unassignedPath['contig'].isna()].reset_index()
    del unassignedPath['index']
    
    # Split the 'contig' column and extract the first part
    unassignedPath['hap'] = unassignedPath['contig'].str.split(pat="_", expand=True)[0]
    
    # Find rows where 'hap' starts with 'unassigned'
    unassigned_rows = unassignedPath['hap'].str.startswith('unassigned')
    
    # Update the 'hap' column for these rows
    unassignedPath.loc[unassigned_rows, 'hap'] = "hapUn"
    unassignedPath['hap'] = (
        unassignedPath['hap']
        .str.replace('sire', sire, regex=False)
        .str.replace('dam', dam, regex=False)
    )
    unassignedPath['path_id'] = unassignedPath['name'].apply(lambda x: re.sub(r".*_utig", "utig", x))
    unassignedPath['new_contig_name'] = (
        unassignedPath['assignChr'] + "_" + unassignedPath['hap'] + "_random_" + unassignedPath['path_id']
    )

    assignedPath = pd.merge(assignedPath,stats, on='contig', how = 'left')
    assignedPath['new_contig_name'] = assignedPath['ref_chr'].astype(str) + "_" + assignedPath['hap'].astype(str)
    assignedPath= assignedPath[['contig','new_contig_name']]

    unassignedPath= unassignedPath[['contig','new_contig_name']]

    final_contigNaming = pd.concat([assignedPath, unassignedPath])
    final_contigNaming['new_contig_name'] = make_unique(final_contigNaming['new_contig_name'])
    # final_contigNaming.to_csv(out_mapFile, sep ='\t', header = None, index=False)
    # Display the updated DataFrame
    return final_contigNaming