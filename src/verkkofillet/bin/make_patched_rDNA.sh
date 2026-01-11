#!/bin/bash
# usage: ./make_patched_rDNA.sh <input_asm> <output_fasta> <contig> <rDNA_consensus> <gapSize> <force>

set -e
set -o pipefail
set -u

if [ "$#" -lt 12 ]; then
    echo "Usage: $0 --input_asm <input_asm> --output_fasta <output_fasta> --contig <contig> --rDNA_consensus <rDNA_consensus> --gapSize <gapSize> --orderGap <orderGap> [ --force <force> --tmp_file <tmp_file> ]"
    echo "  --input_asm: input assembly fasta file"
    echo "  --output_fasta: output patched fasta file"
    echo "  --contig: contig name to be patched"
    echo "  --rDNA_consensus: rDNA consensus fasta file"
    echo "  --gapSize: size of the gap to be inserted"
    echo "  --orderGap: order of the gap to be inserted"
    echo "  --force: whether to overwrite existing files (true/false), default false"
    echo "  --tmp_file: temporary file prefix, default is <input_asm>.tmp"
    exit 1
fi
# Arguments

while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --input_asm)
            input_asm="$2"
            shift; shift
            ;;
        --output_fasta)
            output_fasta="$2"
            shift; shift
            ;;
        --contig)
            contig="$2"
            shift; shift
            ;;
        --rDNA_consensus)
            rDNA_consensus="$2"
            shift; shift
            ;;
        --gapSize)
            gapSize="$2"
            shift; shift
            ;;
        --orderGap)
            orderGap="$2"
            shift; shift
            ;;
        --force)
            force="$2"
            shift; shift
            ;;
        --tmp_file)
            tmp_file="$2"
            shift; shift
            ;;
        *)
            echo "Unknown option $1"
            exit 1
            ;;
    esac
done

# check input
if [ -z "${input_asm:-}" ] || [ -z "${output_fasta:-}" ] || [ -z "${contig:-}" ] || [ -z "${rDNA_consensus:-}" ] || [ -z "${gapSize:-}" ] || [ -z "${orderGap:-}" ]; then
    echo "Error: Missing required arguments."
    exit 1
fi

# set default values
force=${force:-false}
tmp_file=${tmp_file:-""}
if [ -z "$tmp_file" ]; then
    tmp_file="${input_asm}.tmp"
fi

# check inputs
echo "Input assembly: $input_asm"
echo "Output fasta: $output_fasta"
echo "Contig to patch: $contig"
echo "rDNA consensus: $rDNA_consensus"
echo "Gap size: $gapSize"
echo "Order of gap: $orderGap"
echo "Force overwrite: $force"
echo "Temporary file prefix: $tmp_file"
echo " "


# get the prefix without suffix whether .fa or .fasta and fa.gz and fasta.gz
prefix=$(basename "$input_asm")
prefix=${prefix%.fa}
prefix=${prefix%.fasta}
prefix=${prefix%.fa.gz}
prefix=${prefix%.fasta.gz}

# get the directory of the script
script_dir=$(dirname "$(realpath "$0")")

# Extract the contig from assembly
echo -e "Extracting contig $contig from assembly... \n"
samtools faidx $input_asm
samtools faidx $input_asm $contig > $prefix.$contig.fa

# Cleaning the rDNA consensus sequence
echo -e "Cleaning rDNA consensus sequence... \n"

# vm=$(verkko 2>&1 | grep "module path" | awk '{print $4}')
verkkoDir=$(which verkko | sed -e "s/bin\/verkko//")
script=$verkkoDir/lib/verkko/scripts/circularize_ctgs.py
cleanrDNA=${rDNA_consensus%.fa*}.circ.fasta
python $script -o $cleanrDNA $rDNA_consensus

# Run inserting the patch
echo -e "Inserting rDNA patch... \n"
script=$script_dir/insert_rDNA.py
python $script $prefix.$contig.fa $output_fasta $cleanrDNA $gapSize $orderGap $force