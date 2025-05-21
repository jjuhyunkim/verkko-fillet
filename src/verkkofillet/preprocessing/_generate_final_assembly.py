import pandas as pd
import copy
from collections import defaultdict
import subprocess
import os
import sys
import re
import numpy as np
from natsort import natsorted
from itertools import chain



def anotateContig(path, nodeList, screenName):
    path = copy.deepcopy(path)  # Create a deep copy of the DataFrame
    excludelst = []
    for i in range(len(path)):
        path_list = path.loc[i, 'path'].split(",")
        path_list = [re.sub(r"[+-]$", "", s) for s in path_list]
        path_name = path.loc[i, 'name']

        if set(path_list).issubset(set(nodeList)):
            excludelst.append(path_name)

    path.loc[path['name'].isin(excludelst), 'cat'] = screenName
    return path  # ✅ return the updated DataFrame

def annoteContigDict(path, nodeListDict):
    path = copy.deepcopy(path)  # Create a deep copy of the DataFrame
    for key, value in nodeListDict.items():
        value = value.replace(" ", "").split(",")
        path = anotateContig(path, value, key)
    return path  # ✅ return the final annotated DataFrame


def make_cat_column_unique(df, col='cat'):
    if col not in df.columns:
        raise ValueError(f"Column '{col}' not found in DataFrame.")
    df[col+'_unique'] = ""  # Create a new column for unique categories
    counter = defaultdict(int)
    new_cats = []

    for val in df[col]:
        if pd.isna(val):  # Keep NaNs as-is
            new_cats.append(val)
        else:
            counter[val] += 1
            new_cats.append(f"{val}_{counter[val]}")  # Always append _1, _2, etc.

    df[col+'_unique']  = new_cats
    return df

# catList = ['chr14_pat_rDNA', 'chr14_mat_rDNA', 'chr8_pat_rDNA']


def generateJointPathFile(obj, nodeListDict):
    """
    Generate a joint path file from the given object and node list dictionary.
    """
    obj = copy.deepcopy(obj)  # Create a deep copy of the object
    
    test = annoteContigDict(obj.paths, nodeListDict)
    catList = nodeListDict.keys()
    
    paths = pd.merge(test, obj.scfmap, left_on='name', how='left', right_on="pathName")
    paths = paths[paths['cat'].isin(catList)]
    paths = make_cat_column_unique(paths, col='cat')
    paths['contig'] = paths['contig'].fillna("NA")   
    return paths


contigList = ['chr14_pat_rDNA', 'chr14_mat_rDNA', 'chr8_pat_rDNA']
outFasta  = '/data/Phillippy/projects/giraffeT2T/assembly/verkko2.2_hifi-duplex_trio-hic/verkko-thic/giraffe_rDNA.fasta'
oriFasta = '/vf/users/Phillippy/projects/giraffeT2T/assembly/verkko2.2_hifi-duplex_trio-hic/verkko-thic/assembly.fasta'

def writeSeparateFastaFileWithNewName(contigList, oriFasta, outFasta, paths):
    """
    Write separate FASTA files for each contig in contigList.
    """
    if os.path.exists(outFasta):
        print(f"Output file {outFasta} already exists. Please remove it before running the script.")
        sys.exit(1)
    for contigs in contigList:
        print(contigs)
        contigPaths = paths[(paths['cat'] == contigs) & (paths['contig'] != "NA")]
        for contig in contigPaths['cat_unique']:
            contigName =  paths.loc[paths['cat_unique'] == contig,"contig"].values[0]
            # print(contig,contigName)
            try : 
                subprocess.call(f"samtools faidx {oriFasta} {contigName} | sed -e '1d' | sed -e \"1i >{contig}\" >> {outFasta}", shell=True)
            except subprocess.CalledProcessError as e:
                print(f"Error occurred while processing {contigName} {contig}: {e}")
    print("Finished writing separate fasta files.")
    print("Indexing the output fasta file...")
    subprocess.call(f"samtools faidx {outFasta}", shell=True)


def pickPrimaryContigs(qv, tel, gap):
    pd.set_option('future.no_silent_downcasting', True)
    # tel = f"pattern/{name}.telo.bed"
    
    tel = pd.read_csv(tel, sep="\t", header=None)
    tel.columns = ["chrom", "start", "end", "totalLength"]
    tel.loc[tel['start'] == 0, "loc"] = "p"
    tel.loc[tel['end'] == tel['totalLength'], "loc"] = "q"
    tel = tel.groupby("chrom")['loc'].apply(lambda x: ''.join(map(str, x))).reset_index(name = "telomere")
    # tel.head()
    
    # gap = f"pattern/{name}.exclude.bed"
    gap = pd.read_csv(gap, sep="\t", header=None)
    gap.columns = ["chrom", "start", "end"]
    gap = gap['chrom'].value_counts().reset_index(name='gapCount').rename(columns={'index': 'chrom'})
    gap.head()

    

    qv = pd.read_csv(qv, sep="\t", header=None, usecols=[0, 3])
    qv.columns = ["chrom", "qv"]
    qv = qv.drop_duplicates()
    # qv

    df = pd.merge(tel, gap, on="chrom", how="left")
    df = pd.merge(df, qv, on="chrom", how="left")
    df = df.fillna(0)
    df['chrom'] = df['chrom'].str.replace("_withNewGap", "")
    df['gapCount'] = df['gapCount'].astype(int)
    df['qv'] = df['qv'].astype(float)
    df['qv'] = df['qv'].round(2)

    df = df.loc[~df['chrom'].str.contains("random"),]
    df[['chromosome', 'hap']] = df['chrom'].str.split("_", expand=True)

    df.drop_duplicates(subset=['chrom'], inplace=True)
    sorted_chroms = natsorted(df['chromosome'].unique())
    df = df.set_index('chromosome').loc[sorted_chroms].reset_index()
    df.reset_index(drop=True, inplace=True)
    # df.groupby('chromosome').agg({'gapCount': 'sum', 'qv': 'mean'}).reset_index()
    # Group and get first row per chromosome + hap
    df = df.groupby(['chromosome', 'hap'])[['telomere', 'gapCount', 'qv']].first().reset_index()

    # Pivot so each hap becomes a column group
    df = df.pivot(index='chromosome', columns='hap', values=['telomere', 'gapCount', 'qv'])

    # Swap levels so 'hap' is on top
    df.columns = df.columns.swaplevel(0, 1)

    # Sort columns by haplotype (if needed)
    df = df.sort_index(axis=1, level=0)

    # Optional: sort index naturally
    df = df.loc[natsorted(df.index)]

    from itertools import chain

    df_mat = df['mat']
    df_pat = df['pat']
    df_mat.loc[:,'qv'] = df_mat['qv'].fillna(0).astype(float)

    df_pat.loc[:,'qv'] = df_pat['qv'].fillna(0).astype(float)


    df['mainContig_byNumGap'] = ""
    df.loc[df_mat.gapCount > df_pat.gapCount, "mainContig_byNumGap"] = "pat"
    df.loc[df_mat.gapCount < df_pat.gapCount, "mainContig_byNumGap"] = "mat"

    df['mainContig_byQV'] = ""
    df.loc[df_mat.qv < df_pat.qv, "mainContig_byQV"] = "pat"
    df.loc[df_mat.qv > df_pat.qv, "mainContig_byQV"] = "mat"

    df['mainContig'] = ""
    df['mainContig'] = df.apply(lambda row: [row['mainContig_byNumGap'], row['mainContig_byQV']], axis=1)
    df['mainContig'] = df['mainContig'].apply(
        lambda x: ','.join(
            sorted(set(i for i in chain.from_iterable(x) if i))
        )
    )

    primary = ' '.join(list(df.index + "_" + df.mainContig))
    alternate = ' '.join(list(df.index + "_"  + df['mainContig'].replace({'mat': 'pat', 'pat': 'mat'})))

    print(f"primary contigs: {primary}")
    print(f"alternate contigs: {alternate}")

    return df