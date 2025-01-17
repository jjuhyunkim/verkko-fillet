import pandas as pd
import re
import os
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import subprocess

def readGaf(obj, gaf="graphAlignment/verkko.graphAlign_allONT.gaf"):
    """
    Reads a GAF file and stores it as a pandas DataFrame in the provided object.

    Args:
        obj: An object where the parsed GAF data will be stored (as `obj.gaf`).
        gaf (str): Path to the GAF file to be loaded.

    Returns:
        obj: The updated object with the `gaf` attribute containing the DataFrame.
    """
    # Check if obj.gaf already exists and stop if it does
    if obj.gaf is not None:
        print("GAF data already loaded, skipping loading process.")
        return obj
    
    gaf_path = os.path.abspath(gaf)  # Ensure absolute path for compatibility
    print(f"Looking for GAF file at: {gaf_path}")
    
    if os.path.exists(gaf_path):
        print("Loading ONT alignment GAF file...")
        try:
            # Load the GAF file into a pandas DataFrame
            gaf = pd.read_csv(gaf_path,
                              header=None, usecols=[0, 1, 5, 11, 15], sep='\t', low_memory=False, index_col=None,
                              names=['Qname', 'len','path', 'mapQ', 'identity'])
            # Step 1: Modify 'path_modi' column to replace special characters
            gaf['path_modi'] = gaf['path'].str.replace(r'[><\[\]\$]', '@', regex=True).str.replace(r'$', '@', regex=True)
            
            # Step 2: Clean 'identity' column by removing "id:f:" prefix and convert to float
            gaf['identity'] = gaf['identity'].str.replace(r'^id:f:', '', regex=True)
            gaf['identity'] = pd.to_numeric(gaf['identity'], errors='coerce')  # Handle non-numeric values gracefully
            
            # Attach the DataFrame to the object
            obj.gaf = gaf
            print("GAF file successfully loaded.")
            return obj
        except Exception as e:
            print(f"Error loading GAF file: {e}")
            return None
    else:
        print(f"GAF file not found at: {gaf_path}")
        return None

def searchNodes(obj, node_list_input):
    """
    Extracts and filters paths containing specific nodes from the graph alignment file (GAF).
    
    Args:
        obj: An object containing graph alignment data (obj.gaf) and path frequency data (obj.paths_freq).
        node_list_input: A list of node identifiers to search for.

    Returns:
        A styled pandas DataFrame with paths containing the specified nodes and associated frequencies.
    """
    # Prepare node markers
    node_list = [f"@{node}@" for node in node_list_input]
    
    # Check if path frequency data exists, otherwise generate it
    if obj.paths_freq is None:
        print("Path frequency is empty, generating `obj.paths_freq`.")
        
        gaf = obj.gaf  # Assume obj.gaf is a DataFrame
        gaf_size = pd.DataFrame(gaf.groupby('path').size().reset_index())
        gaf_size.columns = ['path', 'size']
        
        # Modify path column by adding '@' around key graph elements
        gaf_size['path_modi'] = (
            gaf_size['path']
            .str.replace(r'[><]', '@', regex=True)  # Replace '>' and '<' with '@'
            .str.replace(r'(?<=\[)', '@', regex=True)  # Add '@' after '['
            .str.replace(r'(?=\])', '@', regex=True)   # Add '@' before ']'
            .str.replace(r'($)', '@', regex=True)      # Add '@' at the start and end
        )
        
        obj.paths_freq = gaf_size
    else:
        print("`obj.paths_freq` already exists.")
    
    # Debug information
    print(f"Extracting paths containing nodes: {node_list_input}")
    
    # Build regex pattern for filtering
    pattern = '|'.join(map(re.escape, node_list))  # Escape special characters
    
    # Filter rows based on presence of nodes in path
    filtered_df = obj.paths_freq[obj.paths_freq['path_modi'].str.contains(pattern, regex=True)]
    
    # Add presence columns for each node
    for node in node_list:
        filtered_df.loc[:, node] = filtered_df['path_modi'].str.contains(node).map({True: 'Y', False: ''})
    
    # Sorting logic
    filtered_df['sort_index'] = filtered_df[node_list].sum(axis=1)
    filtered_df = filtered_df.sort_values(['sort_index', 'size'], ascending=False)
    
    # Drop intermediate columns
    filtered_df.drop(columns=['path_modi', 'sort_index'], inplace=True)
    
    # Escape special HTML characters in the path for better visualization
    filtered_df['path'] = filtered_df['path'].str.replace('<', '&lt;').str.replace('>', '&gt;')
    filtered_df=filtered_df.reset_index()
    del filtered_df['index']
    
    # Styling for display
    headers = {
        'selector': 'th.col_heading',
        'props': 'background-color: #5E17EB; color: white;'
    }
    styled_df = (
        filtered_df.style
        .set_table_styles([headers])
        .bar(color='#FFCFC9', subset=['size'])
        .set_properties(subset=['path'], **{'width': '500px'})
        .set_properties(subset=['size'], **{'width': '50px'})
    )
    
    return styled_df

def searchSplit(obj, node_list_input, min_mapq=0, min_len=50000):
    node_list = node_list_input
    # Create the regex pattern from the node list
    node_pattern = '|'.join(node_list)  # Creates 'utig4-2329|utig4-2651'
    contains_nodes = (
    obj.gaf['path_modi'].str.contains(node_pattern, na=False) &
    (obj.gaf['mapQ'] > 0) &
    (obj.gaf['len'] > 50000)
    )
    filtered_gaf = obj.gaf.loc[contains_nodes, :]
    result = filtered_gaf.groupby("Qname")['path_modi'].agg(set).reset_index()
    target_elements = set([f"@{node}@" for node in node_list_input])
    rows_with_both = result[result['path_modi'].apply(lambda x: target_elements.issubset(x))]
    return rows_with_both

# Use subprocess to run the grep command

def read_Scfmap(scfmap_file = "assembly.scfmap"):
    command = f'grep "^path" {scfmap_file} | cut -d" " -f 2,3'
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    
    # Process the output into a list of lines
    lines = result.stdout.strip().split('\n')
    
    # Convert lines into a DataFrame
    # Assuming space-separated values (adjust delimiter if needed)
    scf = pd.DataFrame([line.split() for line in lines], columns=["fasta_name","path_name"])  # Replace with actual column names
    return scf

def get_NodeChr(obj): 
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
                     exclude_chr = None):
    # read data
    stat = obj.stats[['contig','ref_chr','hap']]
    scf = read_Scfmap(scfmap_file)
    nodeChr = get_NodeChr(obj)
    
    # read HiC data and parsing
    hic = pd.read_csv(hic_support_file, sep =' ', header = None)
    hic.columns = ['X','node1','node2','num']
    # filter the desiring node
    filtered_hic = hic[(hic['node1'] == node) | (hic['node2'] == node)]
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
    if exclude_chr != None:
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