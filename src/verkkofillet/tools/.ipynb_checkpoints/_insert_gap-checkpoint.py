import sys
import shlex
import pandas as pd
import subprocess
import os

script_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../bin/'))

def save_list_to_file(data_list, file_path="insertONTsupport.list"):
    """
    Saves a list to a file with one column. If the folder does not exist, it creates it.

    Parameters:
    - data_list: List of items to be saved.
    - file_path: Path to the file where the data will be saved.
    """
    # Extract directory from the file path
    directory = os.path.dirname(file_path)
    
    # Create the directory if it doesn't exist
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    # Open the file in append mode and write the data
    with open(file_path, 'a') as f:
        for item in data_list:
            f.write(f"{item}\n")


def insertGap(obj,gapid, split_reads,
              outputDir="missing_edge",
              alignGAF="graphAlignment/verkko.graphAlign_allONT.gaf",
              graph="assembly.homopolymer-compressed.gfa"):
    # Ensure absolute paths
    outputDir=os.path.abspath(alignGAF)
    alignGAF=os.path.abspath(alignGAF)
    graph=os.path.abspath(graph)

    
    script_path = subprocess.run(
        "verkko -h | grep 'Verkko module path' | cut -d' ' -f 6",
        shell=True,  # Enables shell commands
        text=True,   # Ensures the output is in text format
        capture_output=True,  # Captures stdout and stderr
        check=True   # Raises an exception for non-zero exit codes
    )
    script_path = script_path.stdout.strip()
    script = os.path.abspath(os.path.join(script_path, "scripts", "insert_aln_gaps.py"))
    
    # Check if the script exists
    if not os.path.exists(script):
        print(f"Script not found: {script}")
        return
    
    # Check if the working directory exists
    if not os.path.exists(outputDir):
        print(f"Working directory not found: {outputDir}")
        return
    
    print("Extracting reads...")
    split_reads = split_reads
    reads = list(set(list(split_reads['Qname'])))
    file_path =  os.path.abspath(os.path.join(outputDir,gapid+".missing_edge.ont_list.txt"))
    
    save_list_to_file(reads, file_path)
    print(f"The split reads for {gapid} was saved to {file_path}")

    subset_gaf = os.path.abspath(os.path.join(outputDir,gapid+".missing_edge.gaf"))

    cmd=f"grep -w -f {file_path} {alignGAF} > {subset_gaf}"
     # print(cmd)
    try:
        # Run the command
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,  # Capture stdout
            stderr=subprocess.PIPE,  # Capture stderr
            shell=True,  # Allow shell-specific syntax
            check=True,  # Raise an exception if the command fails
            cwd=outputDir  # Set working directory
        )
        # Debugging output
        
        # print(f"The fill the edge was done for {gapid}!")
        # print("Standard Output:", result.stdout.decode().strip())
    except subprocess.CalledProcessError as e:
        # Handle errors
        print(f"Command failed: {cmd}")
        print(f"Error code: {e.returncode}")
        print(f"Error output: {e.stderr.decode().strip()}")
        print(f"Standard output: {e.stdout.decode().strip()}")
    
    cmd=f"{script} {graph} {subset_gaf} 1 50000 {outputDir}/patch.nogap.{gapid}.gaf {outputDir}/patch.{gapid}.gaf gapmanual y > {outputDir}/patch.{gapid}.gfa"
    # print(cmd)
    try:
        # Run the command
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,  # Capture stdout
            stderr=subprocess.PIPE,  # Capture stderr
            shell=True,  # Allow shell-specific syntax
            check=True,  # Raise an exception if the command fails
            cwd=outputDir  # Set working directory
        )
        # Debugging output
        
        print(f"The fill the edge was done for {gapid}!")
        # print("Standard Output:", result.stdout.decode().strip())
        print("The final path looks like:")
        print(set(list(pd.read_csv(f"{outputDir}/patch.{gapid}.gaf", header =None, usecols=[5], sep ='\t')[5])))
    except subprocess.CalledProcessError as e:
        # Handle errors
        print(f"Command failed: {cmd}")
        print(f"Error code: {e.returncode}")
        print(f"Error output: {e.stderr.decode().strip()}")
        print(f"Standard output: {e.stdout.decode().strip()}")