#!/bin/bash

# Inputs
original_fasta=$1
trimmed_fasta=$2
trim_bed=$3

# Check if the input files are provided
if [ -z "$original_fasta" ] || [ -z "$trimmed_fasta" ] || [ -z "$trim_bed" ]; then
    echo "Usage: $0 <original_fasta> <trimmed_fasta> <trim_bed>"
    exit 1
fi

# Step 1: If .fai index file doesn't exist for the original fasta, create it
if [ ! -f "$original_fasta.fai" ]; then
    echo "Indexing the original FASTA file..."
    samtools faidx "$original_fasta"
fi

# Step 2: Identify sequences in original_fasta that are not in the trimmed bed file
# Create a temporary file to store the sequence names from the BED file
temp_bed=$(mktemp)
cut -f 1 "$trim_bed" > "$temp_bed"

# Use grep to find sequences in the original FASTA that are not in the BED file
notTrim=$(cut -f 1 "$original_fasta.fai" | grep -w -v -f "$temp_bed" | tr '\n' ' ')

# Clean up the temporary BED file
rm "$temp_bed"

# If no sequences are identified, print a message and exit
if [ -z "$notTrim" ]; then
    echo "No sequences to extract. Exiting."
    exit 1
fi

# Step 3: Create a new FASTA file for the sequences that aren't trimmed
echo "Creating FASTA file with untrimmed sequences..."
cmd="samtools faidx $original_fasta $notTrim > assembly_trimmed.fasta"
echo "Running command: $cmd"
$cmd

# Step 4: Trim the original fasta based on the trim_bed file
if [ ! -f "$trimmed_fasta" ]; then
    echo "Trimming original FASTA with the provided BED file..."
    cmd="samtools faidx $original_fasta -r $trim_bed > $trimmed_fasta"
    echo "Running command: $cmd"
    $cmd
else
    echo "Trimmed FASTA file $trimmed_fasta already exists. Skipping trimming step."
fi

# Step 5: Index the newly trimmed FASTA file
echo "Indexing the trimmed FASTA file..."
cmd="samtools faidx $trimmed_fasta"
echo "Running command: $cmd"
$cmd

echo "Process completed successfully!"

