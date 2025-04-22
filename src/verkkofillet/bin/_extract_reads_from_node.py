import pandas as pd
import sys
import os


node = sys.argv[1] # "utig4-0"
out = os.path.realpath(sys.argv[2]) # "test.out"
layout_node = os.path.realpath(sys.argv[3]) # "unitig-popped.layout"
# scfmap= os.path.realpath(sys.argv[3]) # "unitig-popped.layout.scfmap"
scfmap = layout_node + ".scfmap"

print("input node : " + node)
print("output file : " + out)
print("input scfmap : " + scfmap)
print("input layout : " + layout_node)



# read scfmap
data = []
with open(scfmap, "r") as f:
    lines = f.readlines()
    
    for i in range(0, len(lines), 3):  # Process every 3 lines as a block
        path_line = lines[i].strip().split()  # Extract path details
        piece_line = lines[i+1].strip()  # Extract piece
        # Extract contig and utig values
        contig = path_line[1]
        utig = path_line[2]
        piece = piece_line

        if utig == node:
            data.append([contig, utig, piece])

scfmap = pd.DataFrame(data, columns=["contig", "node", "piece"])

# read layout
node_layout = pd.read_csv(layout_node, header=None)
node_layout.columns = ['layout_info']
node_layout_idx = node_layout.loc[node_layout['layout_info'].fillna('').str.startswith('tig'), 'layout_info'].str.split(r'\s+', expand=True)

node_layout_idx.columns = ['tig','piece']
del node_layout_idx['tig']
node_layout_idx = node_layout_idx.reset_index()

# get piece
piece = scfmap.loc[scfmap['node'] == node, 'piece'].values[0]

# parse layout
idx = node_layout_idx.loc[node_layout_idx['piece'] == piece,'index'].values[0] + 4
idx_next = node_layout_idx.loc[node_layout_idx['piece'] == piece].index[0] + 1
idx_next = node_layout_idx.loc[idx_next,'index'] - 2

layout_sub = node_layout.loc[idx:idx_next,:]['layout_info'].str.split(r'\t', expand=True)
layout_sub.columns = ['read','start_on_node','end_on_node']
layout_sub['start_on_node'] = layout_sub['start_on_node'].astype(int)
layout_sub['end_on_node'] = layout_sub['end_on_node'].astype(int)

# min should be start on node and add direction in + or - at the end
layout_sub['strand'] = "+"
layout_sub.loc[layout_sub['start_on_node'] > layout_sub['end_on_node'], 'strand'] = "-"

layout_sub['start_on_node_fix'] = layout_sub[['start_on_node', 'end_on_node']].min(axis=1)
layout_sub['end_on_node_fix'] = layout_sub[['start_on_node', 'end_on_node']].max(axis=1)

layout_sub = layout_sub[['read', 'start_on_node_fix', 'end_on_node_fix', 'strand']]
layout_sub.columns = ['read', 'start_on_node', 'end_on_node', 'strand']

layout_sub = layout_sub.reset_index(drop=True)

# writing output
layout_sub.to_csv(out, header = True, index = False, sep = '\t')