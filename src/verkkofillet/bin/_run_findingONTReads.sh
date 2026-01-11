#! /bin/bash
verkkofilletdir=$1
verkkodir=$2
threads=$3

if [ ! -d $verkkodir/3-align/split ]; then
    echo "Error: $verkkodir/3-align/split does not exist"
    exit 1
fi 

if [ ! -d $verkkofilletdir/missing_edge ]; then
    echo "Error: $verkkofilletdir/missing_edge does not exist"
    exit 1
fi
# check if seqtk is enable
if ! command -v seqkit &> /dev/null
then
    echo "seqkit could not be found"
    exit 1
fi

if [ ! -d $verkkofilletdir ]; then
    echo "Error: $verkkofilletdir does not exist"
    exit 1
fi

cat $verkkofilletdir/missing_edge/patch.*.gaf |awk '{print $1}' | sort | uniq >> $verkkofilletdir/missing_edge/ont_subset.tmp.id &&
echo -e "$(wc -l ont_subset.tmp.id | awk '{print $1}') ont reads are used for gapfilling"


echo "Extract ont reads from 3-align/split/ont*.fasta.gz"
echo "This step may take a while"
echo " " 
zcat $verkkodir/3-align/split/ont*.fasta.gz | seqkit grep -j $threads -I -n -f $verkkofilletdir/missing_edge/ont_subset.tmp.id > $verkkofilletdir/missing_edge/ont_subset.tmp.fasta &&
shopt -s nullglob
files=("$verkkofilletdir/missing_edge/"*patching.fasta)
if [ ${#files[@]} -eq 1 ] && [ "${files[0]}" = "$verkkofilletdir/missing_edge/*patching.fasta" ]; then
    echo "No patching.fasta files found, skipping."
else
    for f in "${files[@]}"; do
        if [ -s "$f" ]; then
            cat "$f" >> "$verkkofilletdir/missing_edge/ont_subset.tmp.fasta"
        fi
    done
fi
shopt -u nullglob
echo "Extracted ont reads are saved in $verkkofilletdir/missing_edge/ont_subset.tmp.fasta"
echo " "
echo "Done"
echo " "