## How to Run Verkko Consensus Using a Fixed Path

After fixing the paths and filling the gaps with new paths, you can run Verkko consensus using the updated path file based on the original Verkko assembly directory. The `vf.pp.mkCNSdir(obj, new_folder_path)` function uses `obj.verkkoDir` as the base directory and copies the necessary files and folders to `new_folder_path`. The final fixed path (GAF) will be saved in `new_folder_path/6-layoutContigs/consensus_paths.txt`. You can then run the Verkko consensus command as shown below:

```bash
baseDir=base_folder_path
outDir=new_folder_path
path="/path/to/verkko-thic_verkko_fillet/assembly.fixed.paths.gaf"
hifi="/path/to/hifi1.fastq /path/to/hifi2.fastq " # Provide a space-separated list of HiFi read file paths
nano="/path/to/nanopore_reads1.fastq /path/to/nanopore_reads2.fastq"

verkko \
--snakeopts "--dry-run" \
--slurm -d $outDir \
--hifi $hifi --nano $nano \
--paths $path --assembly $baseDir
```

The `--dry-run` option in Verkko will simulate the steps to be executed with the newly generated directory. Ensure that the `layoutContigs` step is included in the list of steps to be executed. Below is an example of the expected job statistics for a Verkko consensus run:

```bash
Job stats:
job                 count    min threads    max threads
----------------  -------  -------------  -------------
buildPackages           1              8              8
cnspath                 1              1              1
combineConsensus        1              8              8
layoutContigs           1              1              1
total                   4              1              8
```

Once verified, you can proceed with the actual run by removing the `--dry-run` option, as shown below:

```bash
verkko \
--slurm -d $outDir \
--hifi $hifi --nano $nano \
--paths $path --assembly $baseDir
```
