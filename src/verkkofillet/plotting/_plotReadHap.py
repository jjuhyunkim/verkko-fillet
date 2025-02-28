import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

def plotHist_readOnNode(nodeinfo, mergedb_all, loc_on_node, node, **kwargs):
    mergedb = mergedb_all.loc[mergedb_all['node'] == node]
    hap = nodeinfo.loc[nodeinfo['node'] == node, 'hap_node'].values[0]
    hap_ratio = nodeinfo.loc[nodeinfo['node'] == node, 'mat:pat'].values[0]

    sns.histplot(data = mergedb, bins = 100, kde =True, **kwargs).set_title(f"""{node} ({hap} {hap_ratio})""")

    start_coor = loc_on_node.loc[loc_on_node['node'] == node, 'start'].values[0] 
    end_coor = loc_on_node.loc[loc_on_node['node'] == node, 'end'].values[0]

    plt.axvline(x=start_coor , color='red', linestyle='--')
    plt.axvline(x=end_coor , color='blue', linestyle='--')

    plt.show()