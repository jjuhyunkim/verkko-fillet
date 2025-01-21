import pandas as pd
import os
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from natsort import natsorted


def qvPlot(obj):
    """
    Generates a bar plot showing QV stats by haplotype.

    Parameters:
    -----------
    obj : verko-fillet object
        An object that contains a `.stats` attribute, which should be a pandas DataFrame.
    """
    # Create figure and axes
    qvTab=obj.qv
    
    fig, ax1 = plt.subplots(figsize=(3, 3))
    
    # Create barplot for QV
    barplot = sns.barplot(x='asmName', y='QV', data=qvTab, ax=ax1, color='grey')
    
    # Add labels to the bars
    for index, bar in enumerate(barplot.patches):
        height = bar.get_height()
        ax1.text(
            bar.get_x() + bar.get_width() / 2,  # X-coordinate (center of bar)
            height-5,                        # Y-coordinate (middle of bar)
            f'{qvTab["QV"][index]:.0f}',       # Label (formatted as an integer)
            ha='center',                       # Horizontal alignment
            va='center',                       # Vertical alignment
            color='white',                     # Text color
            fontsize=10                        # Font size
        )
    
    # Create the second y-axis
    ax2 = ax1.twinx()
    
    # Create a line plot for ErrorRate
    sns.lineplot(x='asmName', y='ErrorRate', data=qvTab, ax=ax2, color='black', label='ErrorRate', marker='o')
    
    # Set axis labels
    ax1.set_ylabel('QV', color='grey')
    ax2.set_ylabel('Error Rate', color='black')
    
    # Set x-axis label
    ax1.set_xlabel('Name of assembly')
    
    # Adjust colors for visibility
    ax1.tick_params(axis='y', colors='grey')
    ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45, ha='right')
    ax2.tick_params(axis='y', colors='black')
    
    # Remove legend from the second axis (if not needed)
    ax2.legend_.remove() if ax2.legend_ else None
    
    # Show the plot
    plt.show()

def completePlot(obj):
    """
    Generates a bar plot showing contig completeness grouped by reference chromosome and haplotype. The completeness of each chromosome is calculated by comparing it to the reference length. A completeness value greater than 100 indicates that the contig length exceeds the original reference length.

    Parameters:
    -----------
    obj : verko-fillet object
        An object that contains a `.stats` attribute, which should be a pandas DataFrame 
        with the following columns:
        - `ref_chr` (str): Reference chromosome identifier.
        - `hap` (str): Haplotype information.
        - `contig_len` (int): Length of the contigs.
    """
    stat_db = obj.stats
    plt.figure(figsize=(5, 4)) 
    sns.barplot(stat_db.groupby(['ref_chr','hap'])['completeness'].sum().reset_index(),
                x="ref_chr", y="completeness", hue="hap")
    plt.title("completeness", fontsize=14)
    plt.xticks(rotation=45)


def contigLenPlot(obj):
    """
    Generates a bar plot showing length of contig by haplotype.

    Parameters:
    -----------
    obj : verko-fillet object
        An object that contains a `.stats` attribute, which should be a pandas DataFrame.
    """
    stat_db = obj.stats
    plt.figure(figsize=(5, 4)) 
    sns.barplot(stat_db.groupby(['ref_chr','hap'])['contig_len'].sum().reset_index(),
                x="ref_chr", y="contig_len", hue="hap")
    plt.title("len(contig)", fontsize=14)
    plt.xticks(rotation=45)

def contigPlot(obj):
    """
    Generates a heatmap of statistics for each haplotype and contig. Brick color represents T2T contigs without gaps, salmon color indicates T2T contigs with gaps, and beige color denotes non-T2T contigs.

    Parameters:
    -----------
    obj : verko-fillet object
        An object that contains a `.stats` attribute, which should be a pandas DataFrame 
    """
    stat_db = obj.stats.copy()
    stat_db.loc[stat_db['t2tStat'] == "not_t2t", "scf_ctg"] = 0
    stat_db.loc[stat_db['t2tStat'] == "scf", "scf_ctg"] = 1
    stat_db.loc[stat_db['t2tStat'] == "ctg", "scf_ctg"] = 2
    
    # Create the pivot table
    ctg = pd.pivot_table(stat_db,values='scf_ctg',index='ref_chr',columns='hap',aggfunc='max')
    
    plt.figure(figsize=(3, 5))  # Adjust the figure size as needed
    sns.heatmap(ctg, cmap="Reds", linecolor="white", linewidths=0.005, cbar=False, vmin=0, vmax=2)
    
    # Display the plot
    plt.show()

    