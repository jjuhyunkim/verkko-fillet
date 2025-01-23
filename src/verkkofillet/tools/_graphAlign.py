import sys
import os
import subprocess
import tqdm
import shlex
import pandas as pd
from tqdm import tqdm

def graphIdx(obj, threads=1,
             GraphAligner_path="GraphAligner", 
             working_directory = "graphAlignment", 
             prefix="verkko.graphAlign", 
             graph = "assembly.homopolymer-compressed.gfa"):
    
    folder_path = os.path.abspath(obj.verkko_fillet_dir)
    working_directory = os.path.abspath(working_directory)
    graph = os.path.abspath(graph)
    folder_path_9 = os.path.abspath(working_directory)
    done_file_path = os.path.join(folder_path_9, "graphIndex.done")
    

    # Ensure the output directory exists
    if not os.path.exists(folder_path_9):
        os.makedirs(folder_path_9)
        print(f"Folder {folder_path_9} created.")
    
    if not os.path.exists(done_file_path):
        print(f"Done file already exists")
        return 
        
    # Proceed only if the index file doesn't exist
    if not os.path.exists(done_file_path):
        print("The diploid index will be stored in:", folder_path_9)
        
        cmd = (
            f"touch {os.path.join(folder_path_9, 'empty.fasta')} && "
            f"{GraphAligner_path} -t {threads} -g {graph} "
            f"-f {os.path.join(folder_path_9, 'empty.fasta')} "
            f"-a {os.path.join(folder_path_9, 'empty.gaf')} "
            f"--diploid-heuristic 21 31 "
            f"--diploid-heuristic-cache {os.path.join(folder_path_9, 'diploid.index')} "
            f"--seeds-mxm-cache-prefix {os.path.join(folder_path_9, prefix)} "
            f"--bandwidth 15 --seeds-mxm-length 30 --mem-index-no-wavelet-tree "
            f"--seeds-mem-count 10000 > {os.path.join(folder_path_9, 'graph_index.log')} && "
            f"touch {done_file_path}"
        )

        try:
            subprocess.run(
                cmd,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,  # Captures errors for debugging
                shell=True,
                check=True,
                cwd= folder_path_9
            )
            print("Indexing completed successfully!")
        except subprocess.CalledProcessError as e:
            print(f"Command failed with error: {e.stderr.decode().strip()}")


def graphAlign(obj, threads=1, GraphAligner_path="GraphAligner",
               prefix="verkko.graphAlign",
               graph = "assembly.homopolymer-compressed.gfa", 
               ontReadList = None,
               working_directory = "graphAlignment"):
    
    # Construct the graph path within the function
    working_directory = os.path.abspath(working_directory)
    graph = os.path.abspath(graph)
    done_file_path = os.path.join(working_directory, "graphAlignment.done")
    
    # Step 1: Create alignment folder if it doesn't exist
    if not os.path.exists(working_directory):
        os.makedirs(working_directory)
        print(f"Folder {working_directory} created.")
    
    # Step 2: Generate ONT read list
    if ontReadList==None:
        ontReadList = os.path.join(working_directory, f"{prefix}.ontReadList.txt")
        ontReads = os.path.join(obj.verkkoDir, "3-align/split/")
        cmd = f"ls {ontReads}*.fasta.gz > {ontReadList}"
        print(f"Command: {cmd}")
        try:
            subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error generating ONT read list: {e}")
            return

    ontReadList = os.path.abspath(ontReadList)
    print(ontReadList)
    # Step 3: Align reads
    gaf_path = os.path.join(working_directory, f"{prefix}_allONT.gaf")
    
    if not os.path.exists(gaf_path):
        print(f"Aligning reads to graph: {graph}")
        read_list = pd.read_csv(ontReadList, header=None)[0].tolist()
        
        for i in tqdm(range(len(read_list)), desc="Processing reads"):
            read_file = read_list[i]
            gaf_file = f"{prefix}_ont{i}.gaf"
            log_file = f"{prefix}_ont{i}.log"
        
            # Safely construct the command
            cmd = (
                f"{shlex.quote(GraphAligner_path)} -t {threads} -g {shlex.quote(graph)} "
                f"-f {shlex.quote(read_file)} -a {shlex.quote(gaf_file)} "
                f"--diploid-heuristic 21 31 "
                f"--diploid-heuristic-cache diploid.index "
                f"--seeds-mxm-cache-prefix {shlex.quote(prefix)} "
                f"--seeds-mxm-windowsize 5000 "
                f"--seeds-mxm-length 30 --seeds-mem-count 10000 "
                f"--bandwidth 15 --multimap-score-fraction 0.99 "
                f"--precise-clipping 0.85 --min-alignment-score 5000 "
                f"--hpc-collapse-reads --discard-cigar "
                f"--clip-ambiguous-ends 100 --overlap-incompatible-cutoff 0.15 "
                f"--max-trace-count 5 --mem-index-no-wavelet-tree > {shlex.quote(log_file)}"
            )
            print(f"Align command for read {i}: {cmd}")
            
            try:
                subprocess.run(cmd, shell=True, check=True, cwd= working_directory)
            except subprocess.CalledProcessError as e:
                print(f"Error during alignment of read {i}: {e}")
                return
        
        # Concatenate GAF files
        concat_cmd = (
            f"cat {prefix}_ont*.gaf > {gaf_path} && "
            f"rm {prefix}_ont* "
        )
        print(f"Concatenation command: {concat_cmd}")
        
        try:
            subprocess.run(concat_cmd, shell=True, check=True, cwd= working_directory)
        except subprocess.CalledProcessError as e:
            print(f"Error during concatenation: {e}")
            return
        
        print(f"Alignment completed. Final GAF file: {gaf_path}")