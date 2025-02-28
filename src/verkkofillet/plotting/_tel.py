import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from .._run_shell import run_shell
import os
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


def percTel(intra_telo , showContig=None,
            width = 5, height = 7, save = True, figName = None):
    
    if showContig is None:
        showContig = ['contig']

    check_columns = ['internal-left', 'internal-right', 'non-internal-left', 'non-internal-right']
    intra_telo = intra_telo.copy()
    intra_telo['by'] = intra_telo[showContig].astype(str).agg('_'.join, axis=1)
    heatmapDb = intra_telo.loc[:,['by']+ check_columns]
    heatmapDb = heatmapDb.set_index('by')
    heatmapDb = heatmapDb.dropna()
    heatmapDb = heatmapDb.sort_values(by=['by'])

    plt.figure(figsize=(width, height))

    ax = sns.heatmap(heatmapDb, annot=True, fmt=".2f", cmap = "Reds", 
                    cbar_kws={'label': 'Telomere Percentage'}, vmin=0, vmax=1,  linewidth=.3)
    ax.set(xlabel="", ylabel="")
    ax.set(title="Telomere Percentage by Contig")

    if figName is None:
        figName = f"figs/intra_telo.heatmap.png"

    if save:
        if not os.path.exists("figs"):
            print("Creating figs directory")
            os.makedirs("figs")

        if os.path.exists(figName):
            print(f"File {figName} already exists")
            print("Please remove the file or change the name")

        elif not os.path.exists(figName):
            plt.savefig(figName)
            print(f"File {figName} saved")

    plt.show()