#! /bin/bash

bed=$1
original_fa=$2
output_fa=$3


# Check the bed file.
if [ ! -f "$bed" ]; then
    echo "Error: Bed file $bed does not exist."
    exit 1
fi
if [ ! -f "$original_fa" ]; then
    echo "Error: Original fasta file $original_fa does not exist."
    exit 1
fi
if [ -f "$output_fa" ]; then
    echo "Error: Output fasta file $output_fa does already exist."
    exit 1
fi

