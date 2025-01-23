import sys
import shlex
import pandas as pd
import subprocess
import os

script_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../bin/'))

def save_list_to_file(data_list, file_path="insertONTsupport.list"):
    """
    Saves a list to a file with one column. If the folder does not exist, it creates it.

    Parameters:
    - data_list: List of items to be saved.
    - file_path: Path to the file where the data will be saved.
    """
    # Extract directory from the file path and handle empty cases
    directory = os.path.dirname(file_path) if os.path.dirname(file_path) else '.'

    # Create the directory if it doesn't exist
    os.makedirs(directory, exist_ok=True)
    
    # Open the file in append mode and write the data
    with open(file_path, 'a') as f:
        for item in data_list:
            f.write(f"{item}\n")


def insertGap(obj, gapid, split_reads,
              outputDir="missing_edge",
              alignGAF="graphAlignment/verkko.graphAlign_allONT.gaf",
              graph="assembly.homopolymer-compressed.gfa"):
    """
    Inserts a gap into the graph using split reads.

    Parameters:
    - obj: verkko fillet obj.
    - gapid: Identifier for the gap.
    - split_reads: Pandas DataFrame containing reads information.
    - outputDir: Output directory (default: 'missing_edge').
    - alignGAF: Path to alignment GAF file (default: 'graphAlignment/verkko.graphAlign_allONT.gaf').
    - graph: Path to graph file (default: 'assembly.homopolymer-compressed.gfa').
    """

    # Ensure absolute paths
    outputDir = os.path.abspath(outputDir)
    alignGAF = os.path.abspath(alignGAF)
    graph = os.path.abspath(graph)
    
    # Check if the working directory exists
    os.makedirs(outputDir, exist_ok=True)

    try:
        # Extract Verkko script path
        script_path_proc = subprocess.run(
            ["verkko", "-h"], 
            text=True, 
            capture_output=True, 
            check=True
        )
        script_path_output = script_path_proc.stdout
        script_path_line = [line for line in script_path_output.splitlines() if "Verkko module path" in line]
        
        if not script_path_line:
            raise ValueError("Verkko module path not found in output.")
        
        verkko_path = script_path_line[0].split()[-1]
        script = os.path.abspath(os.path.join(verkko_path, "scripts", "insert_aln_gaps.py"))
        

        
    except (subprocess.CalledProcessError, ValueError) as e:
        script = os.path.join(script_path,"insert_aln_gaps.py")
        
    
    # Check if the script exists
    if not os.path.exists(script):
        print(f"Script not found: {script}")
        return

    print("Extracting reads...")

    # Ensure the column exists in split_reads
    if 'Qname' not in split_reads.columns:
        print("Error: 'Qname' column not found in input data.")
        return

    reads = list(set(split_reads['Qname']))
    file_path = os.path.abspath(os.path.join(outputDir, f"{gapid}.missing_edge.ont_list.txt"))
    
    save_list_to_file(reads, file_path)
    print(f"The split reads for {gapid} were saved to {file_path}")

    subset_gaf = os.path.abspath(os.path.join(outputDir, f"{gapid}.missing_edge.gaf"))

    # Grep reads from GAF file
    cmd_grep = f"grep -w -f {shlex.quote(file_path)} {shlex.quote(alignGAF)} > {shlex.quote(subset_gaf)}"
    try:
        result = subprocess.run(
            cmd_grep, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE, 
            shell=True, 
            check=True, 
            cwd=outputDir
        )
    except subprocess.CalledProcessError as e:
        print(f"Command failed: {cmd_grep}")
        print(f"Error code: {e.returncode}")
        print(f"Error output: {e.stderr.decode().strip()}")
        return

    # Run Verkko gap insertion script
    patch_nogap = os.path.join(outputDir, f"patch.nogap.{gapid}.gaf")
    patch_gaf = os.path.join(outputDir, f"patch.{gapid}.gaf")
    patch_gfa = os.path.join(outputDir, f"patch.{gapid}.gfa")

    cmd_insert = f"python { shlex.quote(script)} {shlex.quote(graph)} {shlex.quote(subset_gaf)} 1 50000 {shlex.quote(patch_nogap)} {shlex.quote(patch_gaf)} gapmanual y > {shlex.quote(patch_gfa)}"
    
    try:
        result = subprocess.run(
            cmd_insert, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE, 
            shell=True, 
            check=True, 
            cwd=outputDir
        )

        print(f"The gap filling was completed for {gapid}!")
        
        # Display the final path contents
        final_paths = set(pd.read_csv(patch_gaf, header=None, usecols=[5], sep='\t')[5])
        print("The final path looks like:")
        print(final_paths)
        
    except subprocess.CalledProcessError as e:
        print(f"Command failed: {cmd_insert}")
        print(f"Error code: {e.returncode}")
        print(f"Error output: {e.stderr.decode().strip()}")
