import sys
import shlex
import subprocess
import os
import re
import shutil
import pandas as pd
from IPython.display import Image, display
from datetime import datetime
from .._run_shell import run_shell

script_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../bin/'))

def addPadding_to_bed(bed_file, fai, pad = 10000, out_bed_file = None, force = False):
    """
    addPadding_to_bed reads a BED file and adds padding to the start and end positions.
    
    Parameters
    ----------
    bed_file : str
        The path to the BED file.
    pad : int
        The padding size to add to the start and end positions. Default is 10000.
    out_bed_file : str
        The path to the output BED file. If None, the output file will be saved in the same directory as the input file with the suffix ".pad_{pad}.bed".

    Returns
    -------
    None
    """
    if force:
        print("Force mode is enabled. Existing output files will be overwritten.")
        if out_bed_file is None:
            out_bed_file = os.path.splitext(bed_file)[0] + f".pad_{pad}.bed"
    else:
        if out_bed_file is None:
            out_bed_file = os.path.splitext(bed_file)[0] + f".pad_{pad}.bed"
        if os.path.exists(out_bed_file):
            print(f"Output file already exists: {out_bed_file}")
            return

    if not os.path.exists(bed_file):
        print(f"BED file not found: {bed_file}")
        return
    
    if not os.path.exists(fai):
        print(f"FAI file not found: {fai}")
        return
    
    fai = pd.read_csv(fai, sep="\t", header=None)
    fai.columns = ["chrom", "length", "start", "bases", "line_length"]
    
    bed = pd.read_csv(bed_file, sep="\t", header=None)
    bed.columns = ["chrom", "start", "end", "platform"]
    bed["start"] = bed["start"].apply(lambda x: max(x - pad, 0))
    bed["end"] = bed.apply(lambda row: min(row["end"] + pad, fai[fai["chrom"] == row["chrom"]]["length"].values[0]), axis=1)

    print(f"Saving padded BED file to {out_bed_file}, with adding {pad} padding to the start and end positions.")
    bed.to_csv(out_bed_file, sep="\t", header=None, index=None)

def build_sparse_compression_map(uncomp_fasta, comp_fasta = None, mapJson_file = None, dist = 100000, showOnly=False, force = False):
    """
    build_sparse_compression_map builds a sparse compression map for a compressed assembly.

    Parameters
    ----------
    uncomp_fasta : str
        The path to the uncompressed assembly in FASTA format.
    comp_fasta : str
        The path to the compressed assembly in FASTA format. If None, the compressed assembly will be saved in the same directory as the uncompressed assembly with the suffix ".comp.fasta". Default is None.
    mapJson_file : str
        The path to the output JSON file that contains the compression map. If None, the output file will be saved in the same directory as the uncompressed assembly with the suffix ".map.json". Default is None.
    dist : int
    
    Returns
    -------
    None
    """
    if comp_fasta is None:
        comp_fasta = os.path.splitext(uncomp_fasta)[0] + ".comp.fasta"
    if mapJson_file is None:
        mapJson_file = os.path.splitext(uncomp_fasta)[0] + ".map.json"
    if not os.path.exists(uncomp_fasta):
        print(f"Uncompressed assembly not found: {uncomp_fasta}")
        return
    
    if force:
        print("Force mode is enabled. Existing output files will be overwritten.")
        if mapJson_file is None:
            mapJson_file = os.path.splitext(uncomp_fasta)[0] + ".map.json"
    else:
        if mapJson_file is None:
            mapJson_file = os.path.splitext(uncomp_fasta)[0] + ".map.json"
        if os.path.exists(mapJson_file):
            print(f"Output file already exists: {mapJson_file}")
            return
    
    if not os.path.exists(f"{script_path}/build_sparse_compression_map.py"):
        print(f"Script not found: {script_path}/build_sparse_compression_map.py")
        return
    
    cmd=f"python {script_path}/build_sparse_compression_map.py -c {uncomp_fasta} {comp_fasta} {mapJson_file} -d {dist}"

    run_shell(cmd, wkDir=os.getcwd(), functionName="build_sparse_compression_map", longLog=False, showOnly=showOnly)


    
def lift_seqs(uncomp_fasta, mapJson_file, uncomp_bed, comp_bed = None, mode = "uncompressed_to_compressed", showOnly=False, force = False):
    """
    lift_seqs lifts the coordinates of a BED file from uncompressed to compressed space.

    Parameters
    ----------
    uncomp_fasta : str
        The path to the uncompressed assembly in FASTA format.
    mapJson_file : str
        The path to the JSON file that contains the compression map.
    uncomp_bed : str
        The path to the uncompressed BED file.
    uncompressed_to_compressed : str
        The mode of lifting the coordinates. Default is "uncompressed_to_compressed". Other option is "compressed_to_uncompressed".
    comp_bed : str
        The path to the output compressed BED file.
    
    Returns
    -------
    None
    """
    if mode not in ["uncompressed_to_compressed", "compressed_to_uncompressed"]:
        print(f"Invalid mode: {mode}")
        return
    
    if force:
        print("Force mode is enabled. Existing output files will be overwritten.")
        if comp_bed is None:
            comp_bed = os.path.splitext(uncomp_bed)[0] + ".comp.bed"
    else:
        if comp_bed is None:
            comp_bed = os.path.splitext(uncomp_bed)[0] + ".comp.bed"
        if os.path.exists(comp_bed):
            print(f"Output file already exists: {comp_bed}")
            return
    
    if not os.path.exists(uncomp_fasta):
        print(f"Uncompressed assembly not found: {uncomp_fasta}")
        return
    
    if not os.path.exists(mapJson_file):
        print(f"Compression map not found: {mapJson_file}")
        return
    
    if not os.path.exists(f"{script_path}/lift_seqs.py"):
        print(f"Script not found: {script_path}/lift_seqs.py")
        return
    
    if mode == "uncompressed_to_compressed":
        print("Lifting coordinates from uncompressed to compressed space.")
    else:
        print("Lifting coordinates from compressed to uncompressed space.")
    
    print(f"output bed file will be saved in {comp_bed}")

    cmd=f"python {script_path}/lift_seqs.py --{mode} {uncomp_fasta} {mapJson_file} {uncomp_bed} {comp_bed}"
    run_shell(cmd, wkDir=os.getcwd(), functionName="lift_seqs", longLog=False, showOnly=showOnly)


def make_bandage_csv(finalnodes, uncomp_bed, comp_bed, force = False, df_flattened_file = None):
    """
    make_bandage_csv creates a CSV file for Bandage visualization.

    Parameters
    ----------
    finalnodes : DataFrame
        The DataFrame with the final nodes. Output from `getNodes_from_unHPCregion`.
    df_flattened_file : str
        The path to the flattened DataFrame file to be used for Bandage visualization.
    uncomp_bed : str
        The path to the uncompressed BED file.
    comp_bed : str
        The path to the compressed BED file.

    Returns
    -------
    None
    """
    if force:
        print("Force mode is enabled. Existing output files will be overwritten.")
        if df_flattened_file is None:
            df_flattened_file = os.path.splitext(comp_bed)[0] + ".csv"
    else:
        if df_flattened_file is None:
            df_flattened_file = os.path.splitext(comp_bed)[0] + ".csv"
        if os.path.exists(df_flattened_file):
            print(f"Output file already exists: {df_flattened_file}")
            return
    
    if not os.path.exists(uncomp_bed):
        print(f"Uncompressed BED file not found: {uncomp_bed}")
        return
    
    if not os.path.exists(comp_bed):
        print(f"Compressed BED file not found: {comp_bed}")
        return
    
    uncomp = pd.read_csv(uncomp_bed, sep = '\t', header = None, usecols=[0,1,2])
    comp = pd.read_csv(comp_bed, sep = '\t', header = None)

    uncomp.columns = ["chrom","start","end"]
    comp.columns = ["chrom","start","end"]

    uncomplist = list(uncomp['chrom'] + ":"+ uncomp['start'].astype(str) + "-" + uncomp['end'].astype(str))
    complist = list(comp['chrom'] + ":"+ comp['start'].astype(str) + "-" + comp['end'].astype(str))

    if  len(uncomplist) != len(complist):
        print("lengths are not equal")
        sys.exit(1)

    nameMap = pd.DataFrame({"uncomp":uncomplist, "comp":complist})

    df_flattened = finalnodes.explode('nodes', ignore_index=True)
    df_flattened = df_flattened.loc[:,['nodes',"region"]]
    df_flattened.columns = ["node","comp"]

    df_flattened = df_flattened.merge(nameMap, left_on = "comp", right_on = "comp", how = "left")
    df_flattened = df_flattened.loc[:,['node',"uncomp", "comp"]]

    print(f"writing to {df_flattened_file}")
    df_flattened.to_csv(df_flattened_file, sep = '\t', header = True, index = None)