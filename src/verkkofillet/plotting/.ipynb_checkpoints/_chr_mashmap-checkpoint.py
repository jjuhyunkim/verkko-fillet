import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from .._run_shell import run_shell
import os
    
def showMashmapOri(obj, mashmap_out = "assembly.mashmap.out.filtered.out", by = "chr_hap"):
    """
    Generates a bar plot showing the covered regions of the assembly for each reference.

    Parameters:
    -----------
    obj (verko-fillet object):
        An object that contains a .stats attribute, which should be a pandas DataFrame.
    mashmap_out (str, optional):
        The Mashmap output file, aligning the assembly to the reference. Default is assembly.mashmap.out.filtered.out.
    by (str, optional):
        Specifies the grouping method for the plot. Default is "chr_hap". Available options are ['contig', 'all', 'chr_hap'].
    """
    working_dir = os.path.abspath(obj.verkko_fillet_dir)  # Ensure absolute path for the working directory
    
    mashmap = pd.read_csv(working_dir + "/" + mashmap_out , header = None, sep ='\t')
    
    mashmap.columns = ['qname','qlen','qstart','qend','strand','tname','tlen','tstart','tend','nmatch','blocklen','mapQ','id','kc']
    
    mashmap['block_q'] =  mashmap['qend'] - mashmap['qstart']
    # mashmap.head(2)
    # Group the data by 'qname', 'tname', and 'strand'
    grouped = mashmap.groupby(['qname', 'tname', 'strand'])
    data = grouped.agg(
        qlen=('qlen', 'first'),  # Take the first value of qlen as representative
        qcover=('block_q','sum')  # Calculate coverage
    ).reset_index()
    
    # Calculate percentage coverage
    data['qcover_perc'] = data['qcover'] / data['qlen'] * 100
    # Copy the stats DataFrame from obj
    stats = obj.stats.copy()
    # Filter rows based on 'qname' matching contig names in obj.stats['contig']
    contig_list = list(stats['contig'])  # Assuming this is a list of contig names
    data = data.loc[data['qname'].str.contains('|'.join(contig_list)), :]
    
    
    # Create a new column 'name' by concatenating 'ref_chr' and 'hap'
    if by == 'chr_hap':
        stats['by'] = stats['ref_chr'].astype(str) + "_" + stats['hap']
    if by == 'contig':
        stats['by'] = stats['contig']
    if by == 'all':
        stats['by'] = stats['contig'] + '_' + stats['ref_chr'].astype(str) + "_" + stats['hap']
        
    # Display the first two rows of the DataFrame
    data = pd.merge(data,stats,how= 'left',left_on = 'qname', right_on = 'contig')
    data.loc[data['strand'] == "-", 'qcover_perc'] *= -1
    
    # Separate positive and negative values for clarity
    data['positive_qcover'] = data['qcover_perc'].where(data['strand'] == '+', 0)
    data['negative_qcover'] = data['qcover_perc'].where(data['strand'] == '-', 0)
    
    # Sort the data by qname for better visualization
    data = data.sort_values(by = "qcover_perc", ascending = False)
    
    # Prepare the plot
    fig, ax = plt.subplots(figsize=(5, 5))
    
    # Add horizontal bars for positive and negative values
    ax.barh(data['by'], data['positive_qcover'], color='purple', label='Positive Strand')
    ax.barh(data['by'], data['negative_qcover'], color='skyblue', label='Negative Strand')
    
    # Add labels, legend, and gridlines
    ax.set_xlabel('qcover_perc (%)')
    ax.set_ylabel('Contig')
    ax.set_title('Horizontal Stacked Bar Plot with Positive and Negative Strands')
    ax.axvline(0, color='black', linewidth=0.8, linestyle='--')  # Zero line for reference
    ax.legend()
    
    # Adjust layout and show the plot
    plt.tight_layout()
    plt.show()
