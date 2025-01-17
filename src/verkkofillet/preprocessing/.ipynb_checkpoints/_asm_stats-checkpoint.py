import sys
import shlex
import subprocess
import os

bindir = "/vf/users/Phillippy/projects/giraffeT2T/assembly/script/post_verkko_pkg/src/bin/"

def getT2T(obj, fasta="assembly.fasta"):
    script = os.path.abspath(os.path.join(bindir, "getT2T.sh"))  # Ensure absolute path
    working_dir = os.path.abspath(obj.verkkoDir)  # Ensure absolute path for the working directory
    
    if not os.path.exists(script):
        print(f"Script not found: {script}")
        return
    
    if not os.path.exists(working_dir):
        print(f"Working directory not found: {working_dir}")
        return
    
    cmd = f"sh {shlex.quote(script)} {shlex.quote(fasta)}"
    
    try:
        subprocess.run(
            cmd,
            stdout=subprocess.PIPE,  # Capture stdout for debugging
            stderr=subprocess.PIPE,  # Capture stderr for debugging
            shell=True,  # Allowing shell-specific syntax
            check=True,  # Raises CalledProcessError for non-zero return code
            cwd=working_dir  # Run command in the specified working directory
        )
        print("getT2T was done!")
    except subprocess.CalledProcessError as e:
        print(f"Command failed: {cmd}")
        print(f"Error code: {e.returncode}")
        print(f"Error output: {e.stderr.decode().strip()}")
        print(f"Standard output: {e.stdout.decode().strip()}")
