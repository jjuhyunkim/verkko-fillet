import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import copy


def barCovKmer(obj, node_list):
    """"
    Plot the coverage and trio kmers by nodes"

    Parameters
    ----------
    obj: object
        A verkkofillet object
    node_list: list
        A list of nodes to plot
    """
    obj= copy.deepcopy(obj)

    subset = obj.node.loc[obj.node['node'].isin(node_list),:]
    subset = subset.reset_index(drop=True)
    subset['norm_len'] = subset['norm_len'].apply(lambda x: f'{x:.4f}').astype(float)
    subset['mat_len'] = subset['mat']/subset['len'].apply(lambda x: f'{x:.4f}').astype(float)
    subset['pat_len'] = subset['pat']/subset['len'].apply(lambda x: f'{x:.4f}').astype(float)
    subset = subset.sort_values(by=['norm_len'], ascending=False)
    subset = subset.reset_index(drop=True)
    subset

        # Sample Data
    data = subset
    fig, ax1 = plt.subplots(figsize=(6, 4))

    # Bar width and x-ticks adjustment
    bar_width = 0.3
    x = np.arange(len(data['node']))  # Numeric positions for bars

    # Primary bar plot (First Y-Axis) with transparency
    norm_bars = ax1.bar(x, data['norm_len'], width=bar_width, color=data['color'], 
                        label='norm_len', hatch='....', edgecolor='black', alpha=0.8)

    # Secondary y-axis
    ax2 = ax1.twinx()

    # Stack "mat" and "pat" on the secondary axis with different alpha levels
    mat_bars = ax2.bar(x + bar_width + bar_width/10, data['mat_len'], width=bar_width, 
                    color='red', alpha=0.5, label='mat')
    pat_bars = ax2.bar(x + bar_width + bar_width/10, data['pat_len'], width=bar_width, 
                    color='blue', alpha=0.3, label='pat', bottom=data['mat_len'])

    # Add value labels to primary y-axis bars
    for p in ax1.patches:
        ax1.annotate(f'{p.get_height()}', 
                    (p.get_x() + p.get_width() / 2., p.get_height()), 
                    ha='center', va='bottom', fontsize=10, color='black', 
                    xytext=(0, 5), textcoords='offset points', 
                    rotation=70)

    # Set labels and titles
    ax1.set_xlabel('nodes')
    ax1.set_ylabel('Coverage (norm by length)')
    ax2.set_ylabel('Trio kmers (Stacked mat & pat)')
    ax1.set_title('Coverage and trio kmers by nodes')

    # Adjust x-axis labels to match the categories
    ax1.set_xticks(x)
    ax1.set_xticklabels(data['node'], rotation=45)

    # Add merged legend to ax1 (removing individual legends)
    # ax1.legend(handles, labels, loc='upper left', bbox_to_anchor=(1.2, 1), frameon=False)
    legend1 = [mpatches.Patch(facecolor='white', hatch='....', edgecolor='black', label='Coverage'),
               mpatches.Patch(facecolor='white', edgecolor='grey', label='Kmers')]
    legend2 = [
        mpatches.Patch(color='red', label='mat', alpha = .7),
        mpatches.Patch(color='blue', label='pat', alpha = .7),
        mpatches.Patch(color='grey', label='ambiguous', alpha = .7)
    ]
    ax2.legend(handles=legend1 + legend2, loc='upper left', bbox_to_anchor=(1.15, 1), frameon=False)
    
    # Show the plot
    plt.show()
