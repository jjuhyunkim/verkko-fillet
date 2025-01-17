import subprocess
import pandas as pd
import plotly.graph_objects as go
import os

def estLoops(obj, nodeList, gaf="graphAlignment/verkko.graphAlign_allONT.gaf"):
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
