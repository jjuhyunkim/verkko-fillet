import sys
import shlex
import subprocess
import os
import pandas as pd
from .._run_shell import run_shell


script_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../bin/'))

def mkMeryl(obj, fofn, working_directory="kmer", prefix="assembly", fasta = "assembly.fasta", k = 31):
    """\
    Create meryl database for k-mer counting.
    
    Parameters
    ----------
    obj
        Object containing the `verkkoDir` attribute.
    fofn
        Path to the file of filenames (FOFN).
    prefix
        Prefix for the meryl output.
    """
    # Define paths
    working_dir = os.path.abspath(working_directory)
    fasta = os.path.abspath(fasta)  # Ensure correct path format
    script= os.path.join(script_path,"_sumit_meryl.sh")
    
    meryl_db_path = os.path.abspath(os.path.join(working_dir,f"{prefix}_meryl.k{k}.meryl"))
    
    if os.path.exists(meryl_db_path):
        print(f"Meryl database for {prefix} already exists.")
        return
    
    # Check and create the working directory if it does not exist
    os.makedirs(working_dir, exist_ok=True)
    print(f"Working directory: {working_dir}")
        
    # Log file path
    log_file_path = os.path.join(working_dir, "meryl.logs")
    
    # Kmer calculation
    merqury_cmd = f"sh {script} {k} {fofn} {prefix}_meryl"
    try:
        with open(log_file_path, "w") as log_file:
            subprocess.run(
                merqury_cmd,
                shell=True,
                check=True,
                stdout=log_file,  # Redirect standard output to the log file
                stderr=subprocess.STDOUT,  # Redirect standard error to the log file
                cwd=working_dir  # Set the working directory
            )
        print(f"Kmer calculation completed. Logs are saved to {log_file_path}.")
    except subprocess.CalledProcessError as e:
        print(f"Error during kmer calculation: {e}. Check the log file at {log_file_path}.")

def calQV(obj, 
          working_directory = "kmer" , prefix="assembly", fasta = "assembly.fasta", k = 31, showOnly = False, trio = False):
    """
    Perform quality evaluation using Merqury's qv.sh script.

    Parameters
    ----------
    obj
        Object containing the `verkkoDir` attribute.
    working_directory
        The working directory for the analysis. Default is "kmer".
    prefix
        Prefix for the meryl output. Default is "assembly".

    Returns
    -------
        {prefix}.qv_cal.qv
    """
    # Define paths
    print(f"This step might take a while and need lots of memory.")
    print(f"Please make sure you have enough resources.")
    print(f"Or please run this step on a cluster")
    print(f"You can see the command line using `showOnly=True`")

    working_dir = os.path.abspath(working_directory)
    asm = os.path.abspath(fasta)
    
    # Ensure the working directory exists
    os.makedirs(working_dir, exist_ok=True)
    # print(f"Working directory: {working_dir}")
    
    if not os.path.exists(f"{working_dir}/{prefix}_meryl.k{k}.meryl"):
        print(f"There's no meryl db for {prefix} prepared, Please run Meryl first.")
        return
        
    if not os.path.exists(f"{asm}.fai"):
        print(f"There's no FAI file for {asm}, generating one.")
        subprocess.run(f"samtools faidx {asm}", shell=True, check=True, cwd=working_dir)
    
    if os.path.exists(f"{working_dir}/{prefix}.qv_cal.qv"):
        print(f"The QV output file already exists.")
        return
    
    # Extract assembly for the maternal and paternal haplotypes
    if trio:
        print(f"trio mode is on")
        haplotypes = {
            "dam": "assembly_dam.fasta",
            "sire": "assembly_sire.fasta"
        }
        dam=os.path.join(working_dir,"assembly_dam.fasta")
        sire=os.path.join(working_dir,"assembly_sire.fasta")
        
        if not (os.path.exists(dam) and os.path.exists(sire)):
            print(f"Extracting haplotypes from the assembly.")
            for hap, output in haplotypes.items():
                samtools_extract = f"samtools faidx {asm} $(grep '^{hap}' {asm}.fai | cut -f 1 | tr '\\n' ' ') > {output}"
                run_shell(cmd_qv, wkDir=working_dir, functionName=f"samtools_extract_{output}", longLog=False, showOnly=showOnly)
                samtools_idx = f"samtools faidx {output}"
                run_shell(cmd_qv, wkDir=working_dir, functionName=f"samtools_index_{output}", longLog=False, showOnly=showOnly)
        
        # Run Merqury qv.sh for quality evaluation
        hap1 = os.path.abspath(haplotypes["dam"])
        hap2 = os.path.abspath(haplotypes["sire"])
        script = os.path.abspath(os.path.join(script_path, "qv.sh"))

        cmd_qv=f"sh {script} {prefix}_meryl.k{k}.meryl {hap1} {hap2} {prefix}.qv_cal"

        run_shell(cmd_qv, wkDir=working_dir, functionName="qvCal", longLog=False, showOnly=showOnly)
    else:
        print(f"trio mode is off")
        cmd_qv=f"sh {script} {prefix}_meryl.k{k}.meryl {asm} {prefix}.qv_cal"
        run_shell(cmd_qv, wkDir=working_dir, functionName="qvCal", longLog=False, showOnly=showOnly)
