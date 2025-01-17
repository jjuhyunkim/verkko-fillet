import re
import networkx as nx
from collections import Counter
import pandas as pd
import sys
import shlex
import subprocess
import os
import shutil
from .._run_shell import run_shell

script_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../bin/'))

def renameContig(obj, mapFile = "assembly.final.mapNaming.txt", original_fasta= "assembly.fasta", output_fasta = None, showOnly = False):
    working_dir = os.path.abspath(obj.verkko_fillet_dir)
    script = os.path.abspath(os.path.join(script_path, "changeChrName.sh"))  # Assuming script_path is defined elsewhere
    
    if output_fasta is None:  # Use 'is None' for comparison
        prefix = re.sub(r"\.gz|\.fasta", "", original_fasta)
        outFasta = prefix + "_rename.fasta"
    
    # Check if the script exists
    if not os.path.exists(script):
        print(f"Script not found: {script}")
        return
    
    # Check if the working directory exists
    if not os.path.exists(working_dir):
        print(f"Working directory not found: {working_dir}")
        return
        
    if not os.path.exists(mapFile):
        print(f"chromosome map file not found : {mapFile}")
        return
        
    # Construct the shell command
    cmd = f"sh {shlex.quote(script)} {shlex.quote(mapFile)} {shlex.quote(str(original_fasta))} {shlex.quote(outFasta)}"
    
    run_shell(cmd, wkDir=working_dir, functionName = "chrRename" ,longLog = False, showOnly = showOnly)

    print(f"Final renamed fasta file : {outFasta}")
