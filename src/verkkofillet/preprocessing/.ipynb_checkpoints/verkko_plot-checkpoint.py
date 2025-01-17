import pandas as pd
import os
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from natsort import natsorted


def complete_plot(obj):
    stat_db = obj.stats
    plt.figure(figsize=(5, 4)) 
    sns.barplot(stat_db.groupby(['chr','hap'])['completeness'].sum().reset_index(),
                x="chr", y="completeness", hue="hap")
    plt.title("completeness", fontsize=14)
    plt.xticks(rotation=45)


def contigLen_plot(obj):
    stat_db = obj.stats
    plt.figure(figsize=(5, 4)) 
    sns.barplot(stat_db.groupby(['chr','hap'])['contig_len'].sum().reset_index(),
                x="chr", y="contig_len", hue="hap")
    plt.title("len(contig)", fontsize=14)
    plt.xticks(rotation=45)


def contig_plot(obj):
    stat_db = obj.stats
    ctg = pd.pivot_table(stat_db, values='scf_ctg', index=['chr'], columns=['hap'], aggfunc="max")
    ctg['chr_num'] = [s.replace('chr', '') for s in ctg.index.tolist()]
    
    chr_num = []
    for s in ctg.index:
        if s.replace('chr', '').isdigit():  # Check if the remaining part is numeric
            chr_num.append(int(s.replace('chr', '')))
        else:
            chr_num.append(float('inf'))  # Temporarily assign a large value for 'X' or other characters
    
    ctg['chr_num'] = chr_num
    
    # Step 3: Replace 'inf' or non-numeric characters with max + 1
    max_num = max([num for num in chr_num if num != float('inf')])
    ctg['chr_num'] = [max_num + 1 if num == float('inf') else num for num in chr_num]
    
    # Step 4: Sort by 'chr_num' and clean up
    ctg = ctg.sort_values(by='chr_num').drop(columns=['chr_num'])
    ctg.fillna(0, inplace=True)
    # Step 5: Display the result
    plt.figure(figsize=(2, 4)) 
    sns.heatmap(ctg,cmap="Reds", linecolor="white",linewidths = 0.005,cbar=False, vmin=0, vmax=2)
    # sns.heatmap(ctg,cmap="Reds", linecolor="white",linewidths = 0.005,cbar=False)
    # plt.title(, fontsize=14)
    
    # Display the plot
    plt.show()

    