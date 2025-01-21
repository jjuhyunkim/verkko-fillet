import sys
import shlex
import subprocess
import os
import re
import shutil
from IPython.display import Image, display
from datetime import datetime
from .._run_shell import run_shell

script_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../bin/'))

def chrAssign(obj, ref, working_directory = "chromosome_assignment", fasta="assembly.fasta", chr_name="chr", idx=99, showOnly = False):
    """
    Run the script to align the assembly to the given reference using mashmap and obtain the chromosome assignment results.

    Parameters:
    -----------
    obj (verko-fillet object):
        An object that contains a .stats attribute, which should be a pandas DataFrame.
    ref (str) :
        Existing reference
    fasta (str):
        verkko assembly. [default: `assembly.fasta`]
    working_directory (str):
        output directory [default : `./stats/`]
    chr_name (str):
        prefix of the chromosome name in the previous reference. [default : "chr"]
    idx (int):
        Identity threshold to filter mashmap result [defualt : 99]
    showOnly (bool): 
        If set to True, the script will not be executed; it will only display the intended operations. [default : FALSE]

    Return:
    -----------
    {working_directory}/assembly.mashmap.out
    {working_directory}/assembly.mashmap.out.filtered.out
    {working_directory}/chr_completeness_max_hap1
    {working_directory}/chr_completeness_max_hap2
    {working_directory}/translation_hap1
    {working_directory}/translation_hap2
    """
    # Ensure absolute paths
    working_dir = os.path.abspath(working_directory)
    script = os.path.abspath(os.path.join(script_path, "getChrNames.sh"))  # Assuming script_path is defined elsewhere
    fasta = os.path.abspath(fasta)
    ref = os.path.abspath(ref)
    
    # Check if the script exists
    if not os.path.exists(script):
        print(f"Script not found: {script}")
        return
    
    # Check if the working directory exists
    if not os.path.exists(working_dir):
        os.mkdir(working_dir)
    
    output_files = [
        "translation_hap1",
        "translation_hap2",
        "chr_completeness_max_hap1",
        "chr_completeness_max_hap2",
        "assembly.mashmap.out"
    ]

    # Check if all output files already exist
    if all(os.path.exists(os.path.join(working_dir, file)) for file in output_files):
        print("All output files already exist. Skipping chromosome assignment.")
        return
        
    if os.path.exists(f"{working_dir}/assembly.mashmap.out"):
        print(f"The [assembly.mashmap.out] file is already exists")
        print(f"If you want to re-run this job, please detete {working_directory}/[assembly.mashmap.out]")
        return
        
    # Construct the shell command
    cmd = f"bash {shlex.quote(script)} {shlex.quote(ref)} {shlex.quote(str(idx))} {shlex.quote(fasta)} {shlex.quote(chr_name)}"
    
    run_shell(cmd, wkDir=os.getcwd(), functionName = "chrAssign" ,longLog = False, showOnly = showOnly)

    for output in output_files :
        shutil.move(output, f"{working_dir}/{output}")

def convertRefName(fasta, map_file, out_fasta=None, showOnly = False):
    """
    Replace the name in the given FASTA file.
    
    Parameters:
    -----------    
    fasta (str):
        FASTA file in which the contig name is to be replaced
    map_file (str):
        A two-column file, delimited by tabs, containing the old and new contig names.
    showOnly (bool): 
        If set to True, the script will not be executed; it will only display the intended operations. [default : FALSE]
        
    Return:
    -----------
    out_fasta (str):
        output fasta file [default : {prefix}.rename.fa]
    
    """
    # Default out_fasta if not provided
    ref_fasta = os.path.abspath(fasta)
    map_file = os.path.abspath(map_file)

    
    if out_fasta is None:
        # Extract the base name of the file and directory
        basename_fasta = os.path.basename(ref_fasta)
        dir_fasta = os.path.dirname(ref_fasta)
        
        # Always remove the last extension (e.g., .gz, .fa, .fasta)
        if basename_fasta.endswith(".gz"):
            basename_fasta = os.path.splitext(basename_fasta)[0]  # Remove .gz
        
        basename_fasta = os.path.splitext(basename_fasta)[0]  # Remove the actual file extension
        out_fasta = os.path.join(dir_fasta, f"{basename_fasta}.rename.fa")

    # Check if the output file already exists
    if os.path.exists(out_fasta):
        print(f"The renamed reference fasta already exists: {out_fasta}")
        return
    
    # Construct the awk command to replace headers
    cmd = f"awk 'NR==FNR {{map[$1]=$2; next}} /^>/ {{header=substr($1,2); if (header in map) $1=\">\" map[header];}} {{print}}' {shlex.quote(map_file)} {shlex.quote(ref_fasta)} > {shlex.quote(out_fasta)}"
    run_shell(cmd, wkDir=working_dir, functionName = "convertRefName" ,longLog = False, showOnly = showOnly)


def showPairwiseAlign(obj, 
                      size="large", 
                      working_directory = "chromosome_assignment",
                      mashmap_out="chromosome_assignment/assembly.mashmap.out", 
                      prefix="refAlign", 
                      idx=0.99, 
                      minLen=50000, 
                      showOnly = False):
    # Ensure paths are absolute
    working_dir = os.path.abspath(working_directory)
    script = os.path.abspath(os.path.join(script_path, "generateDotPlot"))
    log_file = os.path.join(working_dir, "logs", "showPairwiseAlign.log")
    mashmap_out = os.path.abspath(mashmap_out)
    
    # Check if gnuplot is available
    gnuplot_path = shutil.which("gnuplot")
    if not gnuplot_path:
        print(f"Command 'gnuplot' is not available.")
        return

    # Filtering command
    cmd1 = (
        f"awk -F'\t' '{{ split($13, arr, \":\"); "
        f"if ((arr[3] > {idx}) && ($4 - $3 > {minLen})) print }}' "
        f"{shlex.quote(mashmap_out)} > {shlex.quote(mashmap_out)}.filtered.out"
    )
    run_shell(cmd1, wkDir=working_dir, functionName = "showPairwiseAlign_1" ,longLog = False, showOnly = showOnly)

    # Generate plot command
    cmd2 = f"perl {shlex.quote(script)} png {shlex.quote(size)} {shlex.quote(mashmap_out)}.filtered.out"
    run_shell(cmd2, wkDir=working_dir, functionName = "showPairwiseAlign_2" ,longLog = False, showOnly = showOnly)

    # Rename output files
    output_files = ['out.fplot', 'out.rplot', 'out.gp', 'out.png']
    for file in output_files:
        old_path = os.path.join(working_dir, file)
        new_path = os.path.join(working_dir, f"{prefix}.{file.split('.')[1]}")
        if os.path.exists(old_path):
            os.rename(old_path, new_path)

    # Display the PNG image
    image_path = os.path.join(working_dir, f"{prefix}.png")
    if os.path.exists(image_path):
        display(Image(filename=image_path, width=500))
    else:
        print(f"Image {image_path} not found.")