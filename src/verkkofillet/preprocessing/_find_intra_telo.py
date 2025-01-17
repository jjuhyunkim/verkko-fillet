import pandas as pd
import os
import re

def find_intra_telo(obj, file="internal_telomere/assembly_1/assembly.windows.0.5.bed", 
                    fai_file = "assembly.fasta.fai" , loc_from_end=15000):
    
    working_dir = os.getcwd()
    file = os.path.abspath(file)
    fai_file = os.path.abspath(fai_file)
    
    tel = pd.read_csv(file, sep='\t', header=None)
    tel.columns = ['contig', 'start', 'end']
    
    # Check if the fai file exists, if not, run the command to generate it
    if not os.path.exists(fai_file):  # Corrected the file existence check
        cmd = f"samtools faidx {fai_file}"
        os.system(cmd)  # This will run the command to create the fai file
    
    fai = pd.read_csv(fai_file, sep='\t', header=None, usecols=[0, 1])
    fai.columns = ['contig', 'len_fai']
    
    tel = tel.loc[tel['contig'].str.startswith(('dam', 'sire'))]
    tel = tel.groupby('contig', as_index=False).agg({
        'start': 'min',  # Take the minimum start value
        'end': 'max',    # Take the maximum end value
    })
    
    tel_merged = pd.merge(tel, fai, on='contig', how='inner')
    
    # Filter for regions where the start and end conditions are met
    int_telo = tel_merged[~((tel_merged['len_fai'] - tel_merged['end'] > loc_from_end) |
                            (tel_merged['start'] < loc_from_end))]
    int_telo = pd.merge(obj.stats, int_telo, how = 'inner',on='contig').loc[:,['contig','old_chr','ref_chr','hap','start','end','len_fai']]
    return int_telo


def find_reads_intra_telo(intra_telo, pos,scfmap = "assembly.scfmap",layout = "6-layoutContigs/unitig-popped.layout"):
    print("Finding the reads support for the additional artifical sequences outside of the telomere...")
    contig = str(intra_telo['contig'][0])
    if pos == 'start':
        bp = int(intra_telo['start'][0])
    elif pos == 'end' :
        bp = int(intra_telo['end'][0])
    else :
        print ("the pos argument should be either start or end")
        return
    len_fai = int(intra_telo['len_fai'][0])

    with open(scfmap, 'rb') as f:
        data = f.read().decode('utf-8')  # Decode bytes to string
    
    # Regular expression to match 'path' to 'end'
    pattern = r'(path.*?end)'
    
    # Find all matches
    matches = re.findall(pattern, data, re.DOTALL)
    filtered_matches = [match for match in matches if contig in match]
    
    # Regular expression to match all pieces
    pattern = r'piece\d{6}'
    
    # Extract pieces from the data string
    pieces = re.findall(pattern, filtered_matches[0])
    
    # Get the first and last piece
    first_piece = pieces[0] if pieces else None
    last_piece = pieces[-1] if pieces else None
    
    # Output the results
    # print(f"First piece: {first_piece}")
    # print(f"Last piece: {last_piece}")
    
    if pos == "start":
        piece = first_piece
    elif pos == "end":
        piece = last_piece
    else:
        print("pos should be either start or end")

    print("Looking for the reads from " + piece)
    
    with open(layout, 'rb') as f:
        data = f.read().decode('utf-8')  # Decode bytes to string
    
    # Regular expression to match 'path' to 'end'
    pattern = r'(tig.*?end)'
    
    # Find all matches
    matches = re.findall(pattern, data, re.DOTALL)
    filtered_matches = [match for match in matches if piece in match]
    filtered_matches = filtered_matches[0].split("\n")

    filtered_matches_body = filtered_matches[4:-1]
    filtered_matches_body = [entry.split("\t") for entry in filtered_matches_body]
    df = pd.DataFrame(filtered_matches_body, columns=["readName", "5prime", "3prime"])
    
    df['5prime'] = df['5prime'].astype(int)
    df['3prime'] = df['3prime'].astype(int)
    
    df['start'] = df[['3prime', '5prime']].min(axis=1) * 1.5  # multiply 1.5 cuz this is baesd on HPC coordinates
    df['end'] = df[['3prime', '5prime']].max(axis=1) * 1.5    # multiply 1.5 cuz this is baesd on HPC coordinates
    
    pieceinfo = filtered_matches[0:4]
    pieceinfo = [entry.split("\t") for entry in pieceinfo]

    if pos == "start":
        df_sub = df.loc[(df['start'] < bp)|(df['end'] < bp)]
    elif pos == "end":
        bp_new = int(pieceinfo[1][1]) - (len_fai - bp)
        df_sub = df.loc[(df['start'] > bp_new)|(df['end'] > bp_new)]
    else:
        print("pos should be either start or end")
    df_sub['type'] = df_sub['readName'].apply(lambda x: 'ont' if ';' in x else 'hifi')
    df_sub['type'] = pd.Categorical(df_sub['type'], categories=['ont','hifi'], ordered=True)
    df_sub_count = df_sub.groupby('type')['5prime'].count().reset_index()
    
    print("Summary : ")
    print("   Num of ONT reads : " + str(df_sub_count.iloc[0,1]))
    print("   Num of HiFi reads : " + str(df_sub_count.iloc[1,1]))
    
    return df_sub