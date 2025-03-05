#! /bin/bash

verkkofilletdir=$1
verkkodir=$2
newfolder=$3
finalGaf=$4

finalGaf=$(realpath $finalGaf)
newfolder=$(realpath $newfolder)
verkkofilletdir=$(realpath $verkkofilletdir)
verkkodir=$(realpath $verkkodir)


# check if seqtk is enable
if ! command -v seqkit &> /dev/null
then
    echo "seqkit could not be found"
    exit 1
fi

echo -e "Update CNS dir with missing edges"
echo " "


if [ ! -d $verkkodir ]; then
    echo "Error: $verkkodir does not exist"
    exit 1
fi

if [ ! -d $verkkofilletdir ]; then
    echo "Error: $verkkofilletdir does not exist"
    exit 1
fi

if [ ! -d $newfolder ]; then
    echo "Error: $newfolder does not exist"
    exit 1
fi

if [ ! -f $finalGaf ]; then
    echo "Error: $finalGaf does not exist"
    exit 1
fi


echo "processing 7-consensus directory..."
mkdir -p $newfolder/7-consensus &&
cd $newfolder/7-consensus/ &&
cp $verkkodir/7-consensus/ont_subset.* ./ &&
chmod a+w * &&
touch ont_subset.tmp.id
cat $verkkofilletdir/missing_edge/patch.*.gaf |awk '{print $1}' | sort | uniq >> ont_subset.tmp.id &&
echo -e "$(wc -l ont_subset.tmp.id | awk '{print $1}') ont reads are used for gapfilling"
echo " "

# check if the fasta file is compressed
gunzip ont_subset.fasta.gz &&
cat ont_subset.tmp.id >> ont_subset.id &&
echo "Extract ont reads from 3-align/split/ont*.fasta.gz"
echo "This step may take a while"
echo " " 

if [ ! -d $verkkodir/3-align/split ]; then
    echo "Error: $verkkodir/3-align/split does not exist"
    exit 1
fi  
# check if the ontReadsList.txt exists
# for i in `ls $verkkodir/3-align/split/ont*gz`; do  echo $i; zcat $i | grep "^>" | awk -v fileName=$i '{print fileName,$0}' | grep -f ont_subset.tmp.id >> ontReadsList.txt;done
# fastaList=$(grep -f ont_subset.tmp.id ontReadsList.txt | awk '{print $1}'| tr '\n' ' ' )
# zcat $fastaList |seqtk subseq - ont_subset.tmp.id >> ont_subset.fasta  || true

# Calculate total size of all gzipped files for accurate progress
touch ont_subset.tmp.fasta
zcat $verkkodir/3-align/split/ont*.fasta.gz | seqkit grep -j "$(nproc)" -I -n -f ont_subset.tmp.id > ont_subset.tmp.fasta &&
cat ont_subset.tmp.fasta >> ont_subset.fasta &&
bgzip ont_subset.fasta &&

echo "7-consensus directory is updated"
cd ..

echo "processing 6-layoutContigs directory..."
cp -r $verkkodir/6-layoutContigs/ .
chmod -R a+w 6-layoutContigs/ &&
cd 6-layoutContigs/
chmod a+w * &&
rm consensus_paths.txt
cat $verkkofilletdir/missing_edge/patch.*.gaf >> combined-alignments.gaf &&
cat $verkkofilletdir/missing_edge/patch.*.gaf |grep "^L" |grep gap >> combined-edges.gfa &&
cat $verkkofilletdir/missing_edge/patch.*.gaf| grep gap | awk 'BEGIN { FS="[ \t]+"; OFS="\t"; } ($1 == "S") && ($3 != "*") { print $2, length($3); }' >> nodelens.txt &&
cp $finalGaf ./consensus_paths.txt &&
cat ../7-consensus/ont_subset.tmp.id >> ont-gapfill.txt &&


echo "running replace_path_nodes.py"
verkkoLib=$(verkko | grep "Verkko module"| awk '{print $4'})
$verkkoLib/scripts/replace_path_nodes.py ../4-processONT/alns-ont-mapqfilter.gaf combined-nodemap.txt |grep -F -v -w -f ont-gapfill.txt > ont.alignments.gaf &&
cd ..

echo "6-layoutContigs directory is updated"