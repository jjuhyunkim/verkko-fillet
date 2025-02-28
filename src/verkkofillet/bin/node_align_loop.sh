#! /bin/bash
# set -x

# check the input args, print usage if not correct

if [ "$#" -lt 4 ]; then
    echo "Usage: $0 --hpc_fa <fa> --fai <fai> --mode <mode> --verkko_dir <verkko_dir> --path <path>"

    echo -e "\t--hpc_fa <fa>: HPC fasta file (required)"
    echo -e "\t--fai <fai>: FAI file (optional)"
    echo -e "\t--mode <mode>: OPTIONAL. Must be either 'nohup' or 'slurm'. Default is 'nohup'."
    echo -e "\t--verkko_dir <verkko_dir>: Verkko directory (required)"
    echo -e "\t--path <path>: Path to the verkko directory (required)"

    exit 1
fi

mkdir -p logs

# get arguments
while [ "$#" -gt 0 ]; do
    case $1 in
        --hpc_fa) hpc_fa=$2; shift 2;;
        --fai) fai=$2; shift 2;;
        --mode) mode=$2; shift 2;;
        --verkko_dir) verkko_dir=$2; shift 2;;
        --path) path=$2; shift 2;;
        *) echo "Unknown parameter passed: $1"; exit 1;;
    esac
done

# check mode is either "nohup" or "slurm", if not set, default is "nohup"
if [ -z $mode ]; then
    mode="nohup"
fi

# check mode is either "nohup" or "slurm"
if [ $mode != "nohup" ] && [ $mode != "slurm" ]; then
    echo "Error: mode must be either nohup or slurm"
    exit 1
fi

# check verkko directory exists
if [ ! -d $verkko_dir ]; then
    echo "Error: $verkko_dir does not exist"
    exit 1
fi

# check if fai file exists
if [ ! -f $fai ]; then
    echo "Error: $fai does not exist"
    exit 1
fi

for contig in `cut -f 1 $fai`
do echo $contig
name=$contig.node_to_noHPC
log=logs/$name.%A_%a.log
script=$(dirname $(readlink -f $0))/node_align.sh
partition=norm
walltime=1-00:00:00
args="--hpc_fa $hpc_fa --contig $contig --verkko_dir $verkko_dir --path $path"
mem=100g
cpus=50
# if mode is "slum", then use slurm
if [ $mode == "slurm" ]; then
    echo -e "slurm"
    sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus --time=$walltime --error=$log --output=$log $script $args
elif [ $mode == "nohup" ]; then
    echo -e "nohup"
    nohup $script --hpc_fa $fa --contig $contig --verkko_dir $verkko_dir --path $path 1> $log 2> $log &
fi
done
