import sys
import shlex
import subprocess
import os

script_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../bin/'))

def findRevChr(obj, size = "large", mashmap_out = "assembly.mashmap.out"):
    working_dir = os.path.abspath(obj.verkkoDir)
    script = os.path.abspath(os.path.join(script_path, "generateDotPlot"))

# Construct the shell command
    cmd = f"perl {shlex.quote(script)} png {shlex.quote(size)} {shlex.quote(mashmap_out)}"
    
    try:
        # Run the command
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,  # Capture stdout
            stderr=subprocess.PIPE,  # Capture stderr
            shell=True,  # Allow shell-specific syntax
            check=True,  # Raise an exception if the command fails
            cwd=working_dir  # Set working directory
        )
        # Debugging output
        print("Plotting chromosome alignment is done!")
        print("Standard Output:", result.stdout.decode().strip())
    except subprocess.CalledProcessError as e:
        # Handle errors
        print(f"Command failed: {cmd}")
        print(f"Error code: {e.returncode}")
        print(f"Error output: {e.stderr.decode().strip()}")
        print(f"Standard output: {e.stdout.decode().strip()}")