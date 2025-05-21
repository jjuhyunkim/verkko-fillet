import pandas as pd
import re
import copy
import os
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import subprocess


def readGaf(obj, gaf="graphAlignment/verkko.graphAlign_allONT.gaf", force = False):
    """
    Reads a GAF file and stores it as a pandas DataFrame in the provided object.

    Parameters
    ----------
    obj
        An object where the parsed GAF data will be stored (as `obj.gaf`).
    gaf
        Path to the GAF file to be loaded.

    Returns
    -------
    obj: object
        The updated object with the `gaf` attribute containing the DataFrame.
    """
    # Check if obj.gaf already exists and stop if it does
    obj = copy.deepcopy(obj)
    
    if not hasattr(obj, 'gaf'):
        obj.gaf = None
    
    if obj.gaf is not None and not force:
        print("GAF data already loaded, skipping loading process.")
        return obj
    
    if force : 
        print("Force is set to True, regenerating `obj.gaf`.")

    
    gaf_path = os.path.abspath(gaf)  # Ensure absolute path for compatibility
    print(f"Looking for GAF file at: {gaf_path}")
    

    if os.path.exists(gaf_path):
        print("Loading ONT alignment GAF file...")
        try:
            # Load the GAF file into a pandas DataFrame
            gaf = pd.read_csv(gaf,
                    header=None, sep='\t', low_memory=False, index_col=None)
            names=['qname', 'qlen', 'qstart', 'qend', 'strand',
                    'pname', 'plen', 'pstart', 'pend', 'nmatch', 'blocksize',
                    'mapq', 'nm','as','dv','id']
            if gaf.shape[1] == 17:
                names.append('cg')
            
            gaf.columns = names

            # Step 1: Modify 'path_modi' column to replace special characters
            gaf['path_modi'] = gaf['pname'].str.replace(r'[><\[\]\$]', '@', regex=True).str.replace(r'$', '@', regex=True)
            
            # Step 2: Clean 'identity' column by removing "id:f:" prefix and convert to float
            gaf['nm'] = gaf['nm'].str.replace(r'^NM:i:', '', regex=True)
            gaf['nm'] = pd.to_numeric(gaf['nm'], errors='coerce')  # Handle non-numeric values gracefully
            
            # Step 2: Clean 'identity' column by removing "id:f:" prefix and convert to float
            gaf['as'] = gaf['as'].str.replace(r'^AS:f:', '', regex=True)
            gaf['as'] = pd.to_numeric(gaf['as'], errors='coerce')  # Handle non-numeric values gracefully
            
            # Step 2: Clean 'identity' column by removing "id:f:" prefix and convert to float
            gaf['id'] = gaf['id'].str.replace(r'^id:f:', '', regex=True)
            gaf['id'] = pd.to_numeric(gaf['id'], errors='coerce')  # Handle non-numeric values gracefully

            # Step 2: Clean 'identity' column by removing "id:f:" prefix and convert to float
            gaf['dv'] = gaf['dv'].str.replace(r'^dv:f:', '', regex=True)
            gaf['dv'] = pd.to_numeric(gaf['dv'], errors='coerce')  # Handle non-numeric values gracefully
            
            # Attach the DataFrame to the object
            obj.gaf = gaf
            print("GAF file successfully loaded.")
            return obj
        except Exception as e:
            print(f"Error loading GAF file: {e}")
            return obj
    else:
        print(f"GAF file not found at: {gaf_path}")
        return obj

def searchNodes(obj, node_list_input, multimap_filter = 'mapq', force = False):
    """
    Extracts and filters paths containing specific nodes from the graph alignment file (GAF).
    
    Parameters
    ----------
    obj
        An object containing graph alignment data (obj.gaf) and path frequency data (obj.paths_freq).
    node_list_input
        A list of node identifiers to search for.
    multimap_filter
        The column name used to filter the GAF data (default is 'mapq').
    force
        If True, forces the regeneration of obj.paths_freq even if it already exists.

    Returns
    -------
        A styled pandas DataFrame with paths containing the specified nodes and associated frequencies.
    """
    # Prepare node markers
    # obj = copy.deepcopy(obj)
    node_list = [f"@{node}@" for node in node_list_input]
    
    if force:
        print(f"Force is set to True, regenerating `obj.paths_freq`.")

    if obj.paths_freq is not None and not force:
        print(f"`obj.paths_freq` already exists.")
        print(f"skip generating `obj.paths_freq`.")
    # Check if path frequency data exists, otherwise generate it
    elif obj.paths_freq is None or force:
        print(f"Path frequency is empty, generating `obj.paths_freq`.")
        print(f"Filter by best {multimap_filter}")
        
        gaf = obj.gaf.copy()  # Assume obj.gaf is a DataFrame
        # gaf.head()
        # group by 'qname' and pick max of 'mapq' for each group
        if multimap_filter != None:
            gaf = gaf.loc[gaf.groupby('qname')[multimap_filter].idxmax()]
            gaf = gaf.loc[gaf.groupby('qname')[multimap_filter].idxmax()]
        gaf = gaf.reset_index(drop=True)

       
        gaf_size = pd.DataFrame(gaf.groupby('pname').size().reset_index())
        gaf_size.columns = ['pname', 'nsupport']
        
        # Modify path column by adding '@' around key graph elements
        gaf_size['path_modi'] = (
            gaf_size['pname']
            .str.replace(r'[><]', '@', regex=True)  # Replace '>' and '<' with '@'
            .str.replace(r'(?<=\[)', '@', regex=True)  # Add '@' after '['
            .str.replace(r'(?=\])', '@', regex=True)   # Add '@' before ']'
            .str.replace(r'($)', '@', regex=True)      # Add '@' at the start and end
        )
        # gaf_size['path_modi'] = ["@".join(list(x) + ['']) for x in gaf_size['path_modi'].str.split("@")]    
        
        obj.paths_freq = gaf_size
    else:
        print("`obj.paths_freq` already exists.")
    
    # Debug information
    print(f"Extracting paths containing nodes: {node_list_input}")
    
    # Input data
    # node_list_input = ['utig4-2658', 'utig4-2656']
    
    pattern = '|'.join(node_list)  # Create regex pattern for filtering

    # Copy dataset
    paths_freq = obj.paths_freq.copy()

    # Generate reversed path
    paths_freq['reversed_path'] = paths_freq['path_modi'].apply(lambda x: "@".join(reversed(x.split("@"))))

    # Combine paths and remove duplicates
    paths_freq['combined_paths'] = paths_freq.apply(lambda row: "_".join(sorted(set([row['path_modi'], row['reversed_path']]))), axis=1)

    # Filter rows containing specified nodes
    filtered_df = paths_freq[paths_freq['path_modi'].str.contains(pattern, regex=True)].copy()  # Ensure we copy
    
    if filtered_df.empty:
        print("No paths found")
        return
    # Add presence columns for each node
    for node in node_list_input:
        filtered_df[node] = filtered_df['path_modi'].str.contains(node, regex=False).map({True: 'Y', False: ''})

    # Handle NaN values in nsupport before grouping
    filtered_df['nsupport'] = filtered_df['nsupport'].fillna(0)

    # Group and reshape nsupport
    grouped = filtered_df.groupby(['combined_paths'])['nsupport'].apply(list)
    df_result = grouped.apply(lambda x: x[:2] + [None] * (2 - len(x))).apply(pd.Series)
    df_result.columns = ["fw", "rv"]
    df_result = df_result.reset_index()
    df_result['fw'] = df_result['fw'].fillna(0).astype(int)
    df_result['rv'] = df_result['rv'].fillna(0).astype(int)

    # Keep the row with the highest nsupport per group
    filtered_fullinfo = filtered_df.loc[filtered_df.groupby('combined_paths')['nsupport'].idxmax()].drop(columns=['nsupport'])

    # Merge full info with grouped data
    cleaned = filtered_fullinfo.merge(df_result, on='combined_paths')

    # Compute total support
    cleaned['total_support'] = cleaned['fw'] + cleaned['rv']

    # Drop unnecessary columns
    cleaned = cleaned.drop(columns=['combined_paths', 'reversed_path'])

    # Rename columns for clarity
    cleaned.columns = ['path', 'path_modi'] + node_list_input + ['fw', 'rv', 'total_support']

    # Sorting logic
    cleaned['sort_index'] = cleaned[node_list_input].sum(axis=1)  # Sorting key
    cleaned = cleaned.sort_values(by=['sort_index', 'total_support'], ascending=False).drop(columns=['sort_index'])

    # Replace HTML special characters
    cleaned['path'] = cleaned['path'].replace({'<': '&lt;', '>': '&gt;'}, regex=True)
    del cleaned['path_modi']
    cleaned.reset_index(drop=True, inplace=True)
    # Final cleaned DataFrame


    # Styling for display
    headers = {
        'selector': 'th.col_heading',
        'props': 'background-color: #5E17EB; color: white;'
    }
    styled_df = (
        cleaned.style
        .set_table_styles([headers])
        .bar(color='#FFCFC9', subset=['total_support'])
        .set_properties(subset=['path'], **{'width': '500px'})
        .set_properties(subset=node_list_input, **{'width': '30'})
        .set_properties(subset=['total_support'], **{'width': '50px'})
    )
    
    return styled_df

def searchSplit(obj, node_list_input, min_mapq = 0 , min_qlen=5000, min_mapStart = 5000):
    """\
    Searches for paths containing all specified nodes with a minimum mapping quality and length.

    Parameters
    ----------
    obj
        The VerkkoFillet object containing the GAF data.
    node_list_input
        A list of node identifiers to search for.
    min_mapq
        The minimum mapping quality required for a path to be considered. Default is 0.
    min_qlen
        The minimum query length required for a path to be considered. Default is 5000.
    min_mapStart
        The minimum distance from the start of the query or path for a path to be considered. Default is 5000.

    Returns
    -------
        A DataFrame containing the Qname and path_modi columns of paths that meet the criteria.
    """
    # Create the regex pattern from the node list
    obj = copy.deepcopy(obj)
    gaf = obj.gaf.copy()

    node_pattern = '|'.join(node_list_input)  # Creates 'utig4-2329|utig4-2651'
    contains_nodes = (
    gaf['path_modi'].str.contains(node_pattern, na=False) &
    (gaf['mapq'] > min_mapq ) &
    (gaf['qlen'] > min_qlen) &
    ((gaf['qstart'] < min_mapStart) | (gaf['qlen'] - gaf['qend'] < min_mapStart)) & 
    ((gaf['pstart'] < min_mapStart) | (gaf['plen'] - gaf['pend'] < min_mapStart))
    )
    filtered_gaf = obj.gaf.loc[contains_nodes, :]
    result = filtered_gaf.groupby("qname")['path_modi'].agg(set).reset_index()
    target_elements = set([f"@{node}@" for node in node_list_input])
    rows_with_both = result[result['path_modi'].apply(lambda x: target_elements.issubset(x))].reset_index(drop=True)

    num_rows = rows_with_both.shape[0]
    print(f"{num_rows} reads were found that contain both nodes {node_list_input}")
    
    return rows_with_both

# Use subprocess to run the grep command

def read_Scfmap(scfmap_file = "assembly.scfmap"):
    """\
    Read the scfmap file and return a DataFrame with the 'fasta_name' and 'path_name' columns.

    Parameters
    ----------
    scfmap_file
        The path to the scfmap file. Default is "assembly.scfmap".

    Returns
    -------
        A DataFrame containing the 'fasta_name' and 'path_name' columns
    """
    command = f'grep "^path" {scfmap_file} | cut -d" " -f 2,3'
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    
    # Process the output into a list of lines
    lines = result.stdout.strip().split('\n')
    
    # Convert lines into a DataFrame
    # Assuming space-separated values (adjust delimiter if needed)
    scf = pd.DataFrame([line.split() for line in lines], columns=["fasta_name","path_name"])  # Replace with actual column names
    return scf

def get_NodeChr(obj): 
    """\
    Get the node and chromosome mapping from the VerkkoFillet object.
    """
    df = obj.paths[['name','path']]
    df['path'] = df['path'].str.split(',')
    df = df.explode('path')
    df['path'] = df['path'].str.rstrip('+-')
    df = df.reset_index(drop=True)
    return df

def find_hic_support(obj, node, 
                     hic_support_file = "8-hicPipeline/hic.byread.compressed", 
                     max_print = 20, 
                     scfmap_file = "assembly.scfmap", 
                     exclude_chr = ['chrX_mat', 'chrY_pat']):
    """\
    Find HiC support for a specific node.

    Parameters
    ----------
    obj
        The VerkkoFillet object to be used.
    node
        The node for which to find HiC support.
    hic_support_file
        The path to the HiC support file. Default is "8-hicPipeline/hic.byread.compressed".
    max_print
        The maximum number of results to display. Default is 20.
    scfmap_file
        The path to the scfmap file. Default is "assembly.scfmap".
    exclude_chr
        A list of chromosomes to exclude from the results. Default is None.

    Returns
    -------
        dot plot of HiC support for the specified node.
    """
    # read data
    
    if node[-1] in ['+','-']:
        node = node[:-1]
    
    print("Finding HiC support for node:", node)

    # chech the hic file exist
    if not os.path.exists(hic_support_file):
        print("HiC support file not found at:", hic_support_file)
        return
    # if node have + or - at the end, remove it
    hic = pd.read_csv(hic_support_file, sep =' ', header = None)
    hic.columns = ['X','node1','node2','num']
    # filter the desiring node
    filtered_hic = hic[(hic['node1'] == node) | (hic['node2'] == node)]
    if filtered_hic.shape[0] == 0:
        print("No HiC support found for node:", node)
        return
    
    obj = copy.deepcopy(obj)
    stat = obj.stats[['contig','ref_chr','hap']]
    

    scf = read_Scfmap(scfmap_file)
    nodeChr = get_NodeChr(obj)
    
    # read HiC data and parsing
    
    filtered_hic['searchNode'] = node
    filtered_hic['counterpart'] = filtered_hic['node2'].copy()
    filtered_hic.loc[filtered_hic['counterpart'] == node, 'counterpart'] = filtered_hic['node1']
    # merge datasets to map between chromosome naming
    merge = pd.merge(stat,scf, how = 'inner', left_on="contig",right_on= "fasta_name")
    merge = pd.merge(nodeChr,merge, how = 'inner', left_on="name",right_on= "path_name")
    merge = merge[['ref_chr','path','hap']]
    merge.columns = ['ref_chr','node','hap']
    merge = merge.groupby('node').agg(
        hap=('hap', lambda x: set(x)),  # Aggregate 'hap' into a list for each 'node'
        chr=('ref_chr', 'first')       # Keep the first 'chr' value for each 'node'
    ).reset_index()
    merge = merge[merge['node'].str.startswith('utig')]
    # print(merge.head())
    
    merge['hap'] = merge['hap'].apply(lambda x: '-'.join(map(str, x)) if isinstance(x, (set, list)) else x)
    # print(merge.head())
    
    # merge with hic data
    data = pd.merge(merge,filtered_hic,how = 'right', left_on="node", right_on = "counterpart")
    
    # excluding chr if user gave list.
    if exclude_chr is not None:
        data = data[~data['chr'].isin(exclude_chr)]
    
    # sort the data and make index and cut 
    data = data.drop_duplicates()
    data = data.sort_values(by = "num", ascending=False)
    data['index'] = range(1,data.shape[0]+1)
    data = data.head(max_print)
    
    # Sort data to find the top 5 by 'Value'
    data['Label'] = ''  # Initialize empty label column
    data = data.drop_duplicates()
    label_num=10
    
    # Update labels for the top 5
    # Assign values from 'counterpart' to 'Label' for the specified range
    data['Label'] = ""  # Initialize the column
    data = data.reset_index(drop=True)  # Reset index after filtering or merging
    data.loc[:label_num - 1, 'Label'] = data.loc[:label_num - 1, 'counterpart']
    
    
    layout = go.Layout(hovermode=False)
    fig = px.scatter(
        data,
        x='index',
        y='num',
        title='HiC support for ' + node,
        text='Label',  # Add labels from the 'Label' column
        color='chr',
        hover_data={'chr': True, 
                    'node': False, 'X': False, 'node1': False, 
                    'node2': False, 'num': True, 'index': False, 
                    'searchNode': False, 'counterpart': True, 'Label': False, 'hap': True})
    
    fig.update_layout(
        plot_bgcolor='white',  # Set the plot background to white
        xaxis=dict(
            showgrid=True,  # Show gridlines on the x-axis
            gridcolor='lightgrey'  # Set gridline color to light grey
        ),
        yaxis=dict(
            title="num. of HiC link",  # Custom title for the y-axis
            tickformat="d",  # Format ticks as integers (optional)
            showgrid=True  # Optionally show gridlines
        ),
        width=600,  # Figure width
        height=500  # Figure height
    )
    
    
    # Show the plot
    fig.show()