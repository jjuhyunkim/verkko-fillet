import subprocess
import os
import pandas as pd
import pickle
import glob
import sys
import shutil
from .._default_func import check_user_input, print_directory_tree,addHistory
from .._run_shell import run_shell
from datetime import datetime

# --------------------------------------------------------------------------------
# Reading and Writing data files and AnnData objects
# --------------------------------------------------------------------------------

# Function to check user input
script_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../bin/'))

class FilletObj:
    def __init__(self):
        self.verkkoDir = None
        self.verkko_fillet_dir = None
        self.paths = None
        self.version = None
        self.species = None
        self.stats = None
        self.gaps = None
        self.gaf = None
        self.paths_freq = None
        self.qv = None
        self.history = None
        self.scfmap = None

    def __repr__(self):
        attributes = vars(self)
        
        # Filter out None values (optional, or customize as needed)
        existing_attributes = {key: value for key, value in attributes.items() if value is not None}
        
        # Start forming the string for representation
        repr_str = f"{self.__class__.__name__}\n"
        
        # Add each attribute and its value to the string
        for attribute, value in existing_attributes.items():
            # Check if the value is a pandas DataFrame
            if isinstance(value, pd.DataFrame):
                value = ', '.join(value.columns.tolist())  # Convert the column names to a list
            
            repr_str += f"  {attribute}: {value}\n"
        
        return repr_str

def read_Verkko(verkkoDir, 
                verkko_fillet_dir=None, 
                paths_path=None,
                force = False,
                scfmap_path = None, 
                version=None, 
                species=None, 
                lock_original_folder = True, showOnly = False, longLog = False):
    """
    Prepares the Verkko environment by creating necessary directories, locking the original directory, 
    and loading the paths file for further processing.

    Parameters:
    - verkkoDir (str): Base directory of Verkko data.
    - verkko_fillet_dir (str, optional): Target directory for fillet data. Defaults to None.
    - paths_path (str, optional): Path to 'assembly.paths.tsv' file. Defaults to 'assembly.paths.tsv'.
    - version (str, optional): Version of the data. Defaults to None.
    - species (str, optional): Species name. Defaults to None.
    - lock_original_folder (bool, optional): Whether to lock the original directory. Defaults to True.

    Returns:
    - FilletObj: A FilletObj instance with the configured directories and loaded paths data.
    """
    # make filletObj
    obj = FilletObj()
    
    verkkoDir = os.path.realpath(verkkoDir)
    
    # set verkko_fillet output dir
    if verkko_fillet_dir == None:
        verkko_fillet_dir = os.path.join(verkkoDir+"_verkko_fillet")

    # check the verkko fillet output dir
    if os.path.exists(verkko_fillet_dir) and force == False:
        print(f"The Verkko fillet target directory already exists: {verkko_fillet_dir}")
        print(f"If you didn't mean this, please set another directory or for overwirting, please use force= True")
    else:
        # Create the directory if it does not exist
        print(f"The Verkko fillet target directory has been created and set to: {verkko_fillet_dir}")
        print("All temporary and output files will be written to this directory.")
        script = os.path.abspath(os.path.join(script_path, "make_verkko_fillet_dir.sh"))
        cmd=f"sh {script} {verkkoDir} {verkko_fillet_dir}"
        run_shell(cmd, wkDir=verkko_fillet_dir, functionName = "make_verkko_fillet_dir" ,longLog = longLog, showOnly = showOnly)
        
    working_dir=verkko_fillet_dir
    
    # lock original verkko folder to prevent mess up
    if lock_original_folder :
        print(f"Lock the original Verkko folder to prevent it from being modified.")
        script = os.path.abspath(os.path.join(script_path, "lock_folder.sh"))
        cmd=f"sh {script} {verkkoDir}"
        run_shell(cmd, wkDir=working_dir, functionName = "lock_original_folder" ,longLog = longLog, showOnly = showOnly)

    # Set the additional attributes on the object
    obj.species = species
    obj.verkkoDir = verkkoDir
    obj.verkko_fillet_dir = verkko_fillet_dir
    obj.version = version
    obj.history = pd.DataFrame({
    "timestamp": [datetime.now()],
    "activity": [f"verkko-fillet obj is generated. from : {verkkoDir}, outdir : {verkko_fillet_dir}"],
        "function" : "read_Verkko"
})
    # Check and set the paths_path
    if paths_path == None:
        paths_path = os.path.abspath(os.path.join(verkko_fillet_dir, "assembly.paths.tsv"))
    
    # Load paths file if it exists
    if paths_path is not None and os.path.exists(paths_path):
        print(f"Path file loading...from {paths_path}")
        try:
            # Read CSV file with pandas
            obj.paths = pd.read_csv(paths_path, header=0, sep='\t', index_col=None)
            print("Path file loaded successfully.")
            # obj.history = addHistory(f"path file is loaded from {paths_path}")
        except Exception as e:
            print(f"Error loading paths file: {e}")
    else:
        print("Paths file not found or path is None.")
    
    if scfmap_path == None:
        scfmap_path = os.path.abspath(os.path.join(verkko_fillet_dir, "assembly.scfmap"))
    
    if scfmap_path is not None and os.path.exists(scfmap_path):
        print(f"scfmap file loading...from {scfmap_path}")
        try:
            scfmap = pd.read_csv(scfmap_path, sep = ' ', header = None)
            scfmap.columns = ['info','contig','pathName']
            scfmap= scfmap.loc[scfmap['info']=='path']
            del scfmap['info']
            obj.scfmap = scfmap
            print("scfmap file loaded successfully.")
            # obj.history = addHistory(f"path file is loaded from {paths_path}")
        except Exception as e:
            print(f"Error loading paths file: {e}")
    else:
        print("scfmap file not found or path is None.")
        
    return obj


def save_Verkko(obj, fileName):
    print("save verkko fllet obj to -> " + fileName)
    obj = addHistory(obj,f"Writing verkko-fillet obj to {fileName}", 'save_Verkko')
    with open(fileName, "wb") as f:
        pickle.dump(obj, f)

def load_Verkko(fileName):
    print("load verkko fllet obj from <- " + fileName)
    # Open the file in read-binary mode
    with open(fileName, "rb") as f:
        # Load the object from the file using pickle
        obj = pickle.load(f)
        
    obj = addHistory(obj,f"Reading verkko-fillet obj from {fileName}", 'load_Verkko')
    return obj
    
def hard_copy_symlink(symlink_path, destination_path):
    """
    Creates a hard copy of the file pointed to by the symbolic link.
    
    Parameters:
    - symlink_path: Path to the symbolic link.
    - destination_path: Path where the file should be copied.
    """
    if os.path.islink(symlink_path):
        # Get the target of the symbolic link
        target_path = os.readlink(symlink_path)
        #print(f"Symbolic link points to: {target_path}")

        # Copy the actual file to the destination
        shutil.copy(target_path, destination_path)
        #print(f"Hard copy of the symlink created at: {destination_path}")
    else :
        shutil.copy(symlink_path, destination_path)

def mkCNSdir(obj, new_folder_path, final_gaf = "final_rukki_fixed.paths.gaf"):
    newFolder = os.path.abspath(new_folder_path)
    verkkoDir = os.path.abspath(obj.verkkoDir)  # Define oriDir only once
    
    # Check if the new folder exists
    if os.path.exists(newFolder):
        print("New verkko folder for CNS already exists!")
    else:
        # Create the new folder
        os.makedirs(newFolder)  # Using os.makedirs to handle any intermediate directories
        
        try:
            # Create symbolic links for each directory/file
            subprocess.run(["ln", "-s", os.path.join(verkkoDir, "1-buildGraph"), os.path.join(newFolder, "1-buildGraph")], check=True)
            subprocess.run(["ln", "-s", os.path.join(verkkoDir, "2-processGraph"), os.path.join(newFolder, "2-processGraph")], check=True)
            subprocess.run(["ln", "-s", os.path.join(verkkoDir, "3-align"), os.path.join(newFolder, "3-align")], check=True)
            subprocess.run(["ln", "-s", os.path.join(verkkoDir, "3-alignTips"), os.path.join(newFolder, "3-alignTips")], check=True)
            subprocess.run(["ln", "-s", os.path.join(verkkoDir, "4-processONT"), os.path.join(newFolder, "4-processONT")], check=True)
            subprocess.run(["ln", "-s", os.path.join(verkkoDir, "hifi-corrected.fasta.gz"), os.path.join(newFolder, "hifi-corrected.fasta.gz")], check=True)
            subprocess.run(["ln", "-s", os.path.join(verkkoDir, "6-rukki"), os.path.join(newFolder, "6-rukki")], check=True)
            subprocess.run(["ln", "-s", os.path.join(verkkoDir, "5-untip"), os.path.join(newFolder, "5-untip")], check=True)

            # Copy the "6-layoutContigs" directory and ensure path joins correctly
            layoutContigs_folder = os.path.join(newFolder, "6-layoutContigs")
            os.makedirs(layoutContigs_folder, exist_ok=True)  # Ensure directory exists
            source_files = glob.glob(os.path.join(verkkoDir, "6-layoutContigs",'**', '*'), recursive=True)
            for source in source_files:
                try:
                    hard_copy_symlink(source, layoutContigs_folder)
                except subprocess.CalledProcessError as e:
                    print(f"Error copying {source} to {consensus_folder}: {e}")
            subprocess.run(["cp", final_gaf, os.path.join(newFolder, "6-layoutContigs", "consensus_paths.txt")],check=True)
            
            # Create the 7-consensus directory
            consensus_folder = os.path.join(newFolder, "7-consensus")
            os.makedirs(consensus_folder, exist_ok=True)  # Ensure directory exists
            source_files = glob.glob(os.path.join(verkkoDir, "7-consensus", "ont_subset.*"))
            
            for source in source_files:
                try:
                    hard_copy_symlink(source, consensus_folder)
                except subprocess.CalledProcessError as e:
                    print(f"Error copying {source} to {consensus_folder}: {e}")
            
            print(f"Symbolic links and files created from {verkkoDir} to {newFolder}")
            print(" ")
            print_directory_tree(new_folder_path, max_depth=1, prefix="", is_root=True)
        except subprocess.CalledProcessError as e:
            print(f"Error creating symbolic link or copying files: {e}")


def updateCNSdir_missingEdges(obj, new_folder_path):
    """
    Updates the CNS directory by handling missing edges and creating necessary symbolic links or files.
    
    Parameters:
    - obj: Object containing the original verkko directory path.
    - new_folder_path: Path to the new folder to be created.
    """
    newFolder = os.path.abspath(new_folder_path)
    filletDir = os.path.abspath(obj.verkko_fillet_dir)  # Define oriDir only once
    verkkoDir = os.path.abspath(obj.verkkoDir)
    
    # Check if the new folder exists
    if not os.path.exists(newFolder):
        print("New verkko folder for CNS is not exists!")
        return
    
    try:
        # Create symbolic links or handle files in `7-consensus`
        consensus_folder = os.path.join(newFolder, "7-consensus")
        os.makedirs(consensus_folder, exist_ok=True)
        subprocess.run(
            ["bash", "-c", f"cat {filletDir}/missing_edge/patch.*.gaf | awk '{{print $1}}' >> ont_subset.id"], 
            check=True, 
            cwd=consensus_folder
        )
        subprocess.run(["gunzip", "ont_subset.fasta.gz"], check=True, cwd=consensus_folder)
        subprocess.run(
            ["bash", "-c", f"zcat {verkkoDir}/3-align/split/ont*.fasta.gz | seqtk subseq - ont_subset.id >> ont_subset.fasta"], 
            check=True, 
            cwd=consensus_folder
        )
        subprocess.run(["bgzip", "ont_subset.fasta"], check=True, cwd=consensus_folder)

        # Handle files in `6-layoutContigs`
        layout_folder = os.path.join(newFolder, "6-layoutContigs")
        os.makedirs(layout_folder, exist_ok=True)
        subprocess.run(["rm", "consensus_paths.txt"], check=True, cwd=layout_folder)
        subprocess.run(
            ["bash", "-c", f"cat {filletDir}/missing_edge/patch.*.gaf >> combined-alignments.gaf"], 
            check=True, 
            cwd=layout_folder
        )
        subprocess.run(
            ["bash", "-c", f"cat {filletDir}/missing_edge/patch.*.gfa | grep '^L' | grep gap >> combined-edges.gfa"], 
            check=True, 
            cwd=layout_folder
        )
        subprocess.run(
            ["bash", "-c", f"cat {filletDir}/missing_edge/patch.*.gfa | awk 'BEGIN {{ FS=\"[ \\t]+\"; OFS=\"\\t\"; }} ($1 == \"S\") && ($3 != \"*\") {{ print $2, length($3); }}' >> nodelens.txt"], 
            check=True, 
            cwd=layout_folder
        )
        subprocess.run(
            ["bash", "-c", f"tail -n 2 ../7-consensus/ont_subset.id >> ont-gapfill.txt"], 
            check=True, 
            cwd=layout_folder
        )

        # Fetch the Verkko module path and construct script path
        script_path_proc = subprocess.run(
            "verkko -h | grep 'Verkko module path' | cut -d' ' -f 6",
            shell=True,  # Enables shell commands
            text=True,   # Ensures the output is in text format
            capture_output=True,  # Captures stdout and stderr
            check=True   # Raises an exception for non-zero exit codes
        )
        script_path = script_path_proc.stdout.strip()
        script = os.path.abspath(os.path.join(script_path, "scripts", "replace_path_nodes.py"))
        
        # Run the script
        subprocess.run(
            ["bash", "-c", f"{script} ../4-processONT/alns-ont-mapqfilter.gaf ../6-layoutContigs/combined-nodemap.txt | grep -F -v -w -f ../6-layoutContigs/ont-gapfill.txt > ../6-layoutContigs/ont.alignments.gaf || true"], 
            check=True, 
            cwd=layout_folder
        )
        print("All files are updated! the new folder is ready for verkko-cns")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while executing a subprocess command: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

testDir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../data/'))
def loadGiraffe():
    fileName = f"{testDir}/test_giraffe/giraffe_before_gap_filling.pkl"
    obj = load_Verkko(fileName)
    return obj
