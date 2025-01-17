#!/bin/bash

# Get the command-line arguments
mapFile=$1
inputFasta=$2
outputFasta=$3

# Ensure correct number of arguments
if [[ $# -ne 3 ]]; then
  echo "Usage: $0 <mapFile> <inputFasta> <outputFasta>"
  exit 1
fi

# Debugging: Echo input parameters
echo "Map File: $mapFile"
echo "Input FASTA: $inputFasta"
echo "Output FASTA: $outputFasta"

# Create the awk command
cmd="awk 'NR==FNR{a[\$1]=\$2; next} /^>/{header=\$0; for (i in a) if (index(header, i) > 0) {gsub(i, a[i], header)} print header; next} !/^>/{print}' $mapFile $inputFasta > $outputFasta"

# Debugging: Echo the final command
echo "Executing command: $cmd"

# Execute the command
eval $cmd

samtools faidx $outputFasta

