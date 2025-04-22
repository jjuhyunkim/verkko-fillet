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
# grep -w "True" $mapFile > ${mapFile}.tmp
cp $mapFile $mapFile.tmp
# Debugging: Echo input parameters
echo "Map File: $mapFile"
echo "Input FASTA: $inputFasta"
echo "Output FASTA: $outputFasta"

# Create the awk command
cmd1="awk 'NR==FNR{a[\$1]=\$2; next} /^>/{header=\$0; for (i in a) if (index(header, i) > 0) {gsub(i, a[i], header)} print header; next} !/^>/{print}' ${mapFile}.tmp $inputFasta > ${outputFasta}.tmp"
cmd2="samtools faidx  ${outputFasta}.tmp"
cmd3="samtools faidx ${outputFasta}.tmp $(cut -f 2 ${mapFile}.tmp  | tr '\n' ' ') > $outputFasta"
cmd4="samtools faidx $outputFasta"
cmd5="rm ${outputFasta}.tmp*"

echo "Executing command 1: $cmd1" && eval $cmd1 &&
echo "Executing command 2: $cmd2" && eval $cmd2 &&
echo "Executing command 3: $cmd3" && eval $cmd3 &&
echo "Executing command 4: $cmd4" && eval $cmd4 &&
echo "Executing command 5: $cmd5" && eval $cmd5 &&
echo "All commands executed successfully."
