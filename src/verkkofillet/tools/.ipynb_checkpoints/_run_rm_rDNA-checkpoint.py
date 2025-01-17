import sys
import shlex
import subprocess
import os

script_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../bin/'))

def rmrDNA(obj):
    print("Starting removing rDNA nodes in the graph")

    script = os.path.abspath(os.path.join(script_path, "removeRDNA.sh"))  # Ensure absolute path
    working_dir = os.path.abspath(obj.verkko_fillet_dir)  # Ensure absolute path for the working directory

    # Define expected output files
    output_files = [
        "target.screennodes.out",
        "assembly.homopolymer-compressed.noseq.telo_rdna.gfa",
        "assembly.colors.telo_rdna.csv"
    ]

    # Check if all output files already exist
    if all(os.path.exists(os.path.join(working_dir, file)) for file in output_files):
        print("All output files already exist. Skipping rDNA removal.")
        return

    # Ensure script and working directory exist
    if not os.path.exists(script):
        print(f"Script not found: {script}")
        return

    if not os.path.exists(working_dir):
        print(f"Working directory not found: {working_dir}")
        return

    cmd = f"sh {shlex.quote(script)}"

    try:
        subprocess.run(
            cmd,
            stdout=subprocess.PIPE,  # Capture stdout for debugging
            stderr=subprocess.PIPE,  # Capture stderr for debugging
            shell=True,  # Allowing shell-specific syntax
            check=True,  # Raises CalledProcessError for non-zero return code
            cwd=working_dir  # Run command in the specified working directory
        )
        print("remove rDNA was done!")
        print("Output files: ")
        for file in output_files:
            print(file)
    except subprocess.CalledProcessError as e:
        print(f"Command failed: {cmd}")
        print(f"Error code: {e.returncode}")
        print(f"Error output: {e.stderr.decode().strip()}")
        print(f"Standard output: {e.stdout.decode().strip()}")
