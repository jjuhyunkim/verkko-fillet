import sys
import shlex
import subprocess
import os
import pandas as pd


script_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../bin/'))

def mkMeryl(obj, fofn, working_directory="kmer", prefix="assembly", fasta = "assembly.fasta", k = 31):
    """
    Create meryl database for k-mer counting.
    
    Args:
        obj: Object containing the `verkkoDir` attribute.
        fofn: Path to the file of filenames (FOFN).
        prefix: Prefix for the meryl output.
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
          working_directory = "kmer" , prefix="assembly", fasta = "assembly.fasta", k = 31,
         ):
    """
    Perform quality evaluation using Merqury's qv.sh script.

    Args:
        obj: Object containing the `verkkoDir` attribute.
        prefix: Prefix for the meryl output.
    """
    # Define paths
    script = os.path.abspath(os.path.join(script_path, "qv.sh"))
    
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
    haplotypes = {
        "dam": "assembly_dam.fasta",
        "sire": "assembly_sire.fasta"
    }
    dam=os.path.join(working_dir,"assembly_dam.fasta")
    sire=os.path.join(working_dir,"assembly_sire.fasta")
    
    if not (os.path.exists(dam) and os.path.exists(sire)):
        for hap, output in haplotypes.items():
            samtools_cmd = f"samtools faidx {asm} $(grep '^{hap}' {asm}.fai | cut -f 1 | tr '\\n' ' ') > {output}"
            subprocess.run(samtools_cmd, shell=True, check=True, cwd=working_dir)
            subprocess.run(f"samtools faidx {output}", shell=True, check=True, cwd=working_dir)
    
    # Run Merqury qv.sh for quality evaluation
    log_file_path = os.path.join(working_dir, "qv_cal.logs")
    qv_cmd = [
        script,
        f"{prefix}_meryl.k{k}.meryl",
        haplotypes["dam"],
        haplotypes["sire"],
        f"{prefix}.qv_cal"
    ]
    with open(log_file_path, "w") as log_file:
        subprocess.run(
            qv_cmd,
            stdout=log_file,  # Redirect standard output
            stderr=subprocess.STDOUT,  # Redirect standard error to the same log file
            check=True,
            cwd=working_dir
        )
    
    print(f"Quality evaluation completed. Logs are saved to {log_file_path}.")