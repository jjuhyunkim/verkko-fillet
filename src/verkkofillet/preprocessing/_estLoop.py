import subprocess
import pandas as pd
import plotly.graph_objects as go
import os

import itertools
import copy
import pandas as pd
from collections import Counter
import seaborn as sns
import matplotlib.pyplot as plt
import os

def calNodeDepth(obj, hap_info = "assembly.colors.csv", gfa = "assembly.homopolymer-compressed.noseq.gfa", width = 8, height = 5, save = True, figName = None, dpi = 300):
    """
    Plot the depth of nodes in the graph.

    Parameters
    -----------
    obj : VerkkoFillet object
        The VerkkoFillet object.
    width : int
        The width of the plot. Default is 8.
    height : int
        The height of the plot. Default is 5.
    save : bool
        Whether to save the plot. Default is True.
    figName : str
        The name of the figure file. Default is None.
    dpi : int
        The resolution of the figure. Default is 300.
    """


    obj = copy.deepcopy(obj)

    hap_info = os.path.abspath(hap_info)
    gfa = os.path.abspath(gfa)

    if not os.path.exists(hap_info):
        raise FileNotFoundError(f"The hap_info file does not exist: {hap_info}")
    if not os.path.exists(gfa):
        raise FileNotFoundError(f"The GFA file does not exist: {gfa}")

    if obj.paths_freq is None:
        raise ValueError("No paths_freq found. Please run `vf.tl.graphIdx` and `vf.tl.graphAlign` first.")

    # Count the frequency of each node in the paths
    print("Counting the frequency of each node in the paths...")
    paths_freq = obj.paths_freq.copy()
    nested_list = paths_freq['path_modi'].str.split("@").tolist()
    flatlist = list(filter(None, list(itertools.chain(*nested_list))))
    del nested_list

    count_dict = Counter(flatlist)
    count_df = pd.DataFrame(count_dict.items(), columns=['node', 'count'])
    count_df = count_df.sort_values(by='count', ascending=False)
    median = count_df['count'].median()
    print(f"Median: {median}")


    nodeLen = pd.read_csv(gfa, sep='\t', header=None)
    nodeLen = nodeLen[nodeLen[0] == "S"]
    nodeLen = nodeLen.iloc[:,[1,3]]
    nodeLen.columns = ['node', 'len']
    nodeLen['len'] = nodeLen['len'].str.replace(r'^LN:i:', '', regex=True)
    nodeLen['len'] = pd.to_numeric(nodeLen['len'], errors='coerce')  # Handle non-numeric values gracefully
    # nodeLen

    hapdb = pd.read_csv(hap_info, sep='\t', header=0)

    color_map = {
        "#AAAAAA" : "ambiguous",
        "#FF8888" : "hap1",
        "#8888FF" : "hap2",
        "#FFFF00" : "ambiguous"
    }
    hapdb['hap'] = hapdb['color'].map(color_map)
    # hapdb

    merged = pd.merge(nodeLen, hapdb, on='node')
    merged = merged[['node', 'len', 'hap']]
    merged = pd.merge(merged, count_df, on='node')
    

    hapSpecificMed = merged.loc[merged['hap'].isin(["hap1", "hap2"]),"count"].median()
    merged['cov_hap'] = merged['count'] / hapSpecificMed
    print(f"Median of haplotype specific nodes : {hapSpecificMed}")
    
    count_df = merged.copy()
    count_df = count_df.sort_values(by='count', ascending=False)

    # plot line plots
    fig = plt.figure(figsize=(8, 5))
    gs = fig.add_gridspec(1, 2, width_ratios=[4, 1])  # Line plot wider than box plot

    # Line plot
    ax0 = fig.add_subplot(gs[0])
    sns.lineplot(data=count_df, x='node', y='count', marker="o", ax=ax0, markersize=3, alpha=0.5)
    ax0.set_title("Line Plot")
    ax0.set_yscale("log")
    ax0.set(xticklabels=[], xlabel=None, xticks=[])
    ax0.axhline(median, color='r', linestyle='--', label='Median', linewidth=.5, alpha=0.5)

    # Box plot (on the right)
    ax1 = fig.add_subplot(gs[1])
    sns.boxplot(y=count_df["count"], ax=ax1)
    ax1.set_title("Box Plot")
    ax1.set_yscale("log")

    plt.title("Depth of nodes in the graph")
    plt.tight_layout()
    
    if figName is None:
        figName = f"figs/nodeDepth.png"

    if save:
        if not os.path.exists("figs"):
            os.makedirs("figs")

        if os.path.exists(figName):
            print(f"File {figName} already exists")
            print("Please remove the file or change the name")

        elif not os.path.exists(figName):
            plt.savefig(figName, dpi=dpi)
            print(f"File {figName} saved")

    plt.show()

    count_df.reset_index(drop=True, inplace=True)
    count_df['cov'] = count_df['count']/median
    return count_df

def estLoops(obj, nodeList, gaf="graphAlignment/verkko.graphAlign_allONT.gaf"):
    """
    Estimate the number of loops between two nodes in the graph.

    Parameters
    ----------
    obj
        The VerkkoFillet object to be used.
    nodeList
        A list of two nodes to be compared.
    gaf
        Path to the GAF file containing the graph alignment information. Default is "graphAlignment/verkko.graphAlign_allONT.gaf".   
    """
    if obj is None:
        raise ValueError("The 'obj' parameter is required.")
    if nodeList is None or len(nodeList) < 2:
        raise ValueError("The 'nodeList' parameter must contain two nodes.")

    
    node1, node2 = nodeList
    gaf = os.path.abspath(gaf)

    if not os.path.exists(gaf):
        raise FileNotFoundError(f"The GAF file does not exist: {gaf}")
    
    # Step 1: Run the first grep command to search for node2 in the file
    grep1 = subprocess.run(['grep', '-w', node2, gaf], capture_output=True, text=True)
    
    if grep1.returncode != 0:  # Check if grep1 failed
        print(f"No reads detected for node2: {node2}")
        return
    
    # Step 2: Run the second grep command to search for node1 in the output of the previous command
    grep2 = subprocess.run(['grep', '-w', node1], input=grep1.stdout, capture_output=True, text=True)
    
    if grep2.returncode != 0:  # Check if grep2 failed
        print(f"No reads detected for node1: {node1}")
        return
    
    # Proceed with the rest of the pipeline if both grep commands succeed
    try:
        sed1 = subprocess.run(['sed', 's/id:f://g'], input=grep2.stdout, capture_output=True, text=True)
        awk1 = subprocess.run(['awk', '{if ($NF > 0.99 && $3 < 100 && $4 + 100 > $2) print $6}'], input=sed1.stdout, capture_output=True, text=True)
        tr1 = subprocess.run(['tr', '>', ','], input=awk1.stdout, capture_output=True, text=True)
        tr2 = subprocess.run(['tr', '<', ','], input=tr1.stdout, capture_output=True, text=True)
        sed2 = subprocess.run(['sed', f's/{node1}//g'], input=tr2.stdout, capture_output=True, text=True)
        sed3 = subprocess.run(['sed', f's/{node2}//g'], input=sed2.stdout, capture_output=True, text=True)
        sed4 = subprocess.run(['sed', 's/,,//g'], input=sed3.stdout, capture_output=True, text=True)
        awk2 = subprocess.run(['awk', '-F', ',', '{print NF-1}'], input=sed4.stdout, capture_output=True, text=True)
        sort1 = subprocess.run(['sort'], input=awk2.stdout, capture_output=True, text=True)
        uniq = subprocess.run(['uniq', '-c'], input=sort1.stdout, capture_output=True, text=True)
        awk3 = subprocess.run(['awk', '{print $2, $1}'], input=uniq.stdout, capture_output=True, text=True)
        final_output = subprocess.run(['sed', '-e', f"1i {node1}_to_{node2}"], input=awk3.stdout, capture_output=True, text=True)
    except Exception as e:
        print(f"Error during processing: {e}")
        return
    
    # Process and plot results
    lines = final_output.stdout.strip().split('\n')
    header = [lines[0]]  # First line contains the header
    rows = [line.split() for line in lines[1:]]  # Remaining lines are data
    
    df = pd.DataFrame(rows, columns=['nLoop', 'freq'])
    df['nLoop'] = pd.to_numeric(df['nLoop'])
    df['freq'] = pd.to_numeric(df['freq'])
    df = df.sort_values(by="nLoop")
    
    all_nLoops = pd.DataFrame({'nLoop': range(df['nLoop'].min(), df['nLoop'].max() + 1)})
    df = pd.merge(all_nLoops, df, on='nLoop', how='left').fillna({'freq': 0})
    
    # Display DataFrame
    # print(df)
    
    # Plotting
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=df['nLoop'], 
        y=df['freq'],
        mode='lines+markers',
        marker=dict(symbol='circle', size=10, color='blue'),
        name='Freq'
    ))
    fig.update_layout(
        title=str(header),
        xaxis_title='Number of loops',
        yaxis_title='Frequency',
        showlegend=False,
        width=500,
        height=400,
        plot_bgcolor='white',
        xaxis=dict(showgrid=True, gridcolor='grey'),
        yaxis=dict(showgrid=True, gridcolor='grey')
    )
    fig.show()
