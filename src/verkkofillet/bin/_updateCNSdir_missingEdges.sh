#! /bin/bash

ml verkko
ml seqtk


verkkofilletdir=$1
verkkodir=$2
newfolder=$3
finalGaf=$4
missing_edge_dir=$5

finalGaf=$(realpath $finalGaf)
newfolder=$(realpath $newfolder)
verkkofilletdir=$(realpath $verkkofilletdir)
verkkodir=$(realpath $verkkodir)
subsetid=$(realpath $missing_edge_dir)
subsetfasta="$missing_edge_dir/ont_subset.tmp.fasta"

echo -e "finalGaf : $finalGaf"
echo -e "newfolder : $newfolder"
echo -e "verkkofilletdir : $verkkofilletdir"
echo -e "verkkodir : $verkkodir"
echo -e "missing_edge_dir : $missing_edge_dir"


for dir in $newfolder $verkkodir $verkkofilletdir $missing_edge_dir; do
    if [ ! -d $dir ]; then
        echo "Error: $dir does not exist or is not a directory"
        exit 1
    fi
done

for file in $finalGaf $subsetfasta; do
    if [ ! -f $file ]; then
        echo "Error: $file does not exist"
        exit 1
    fi
done


if [ ! -f $missing_edge_dir/ont_subset.tmp.fasta ]; then
    echo "Extracting fasta sequences for missing edges..."
    cat $missing_edge_dir/patch.gapid*id > $missing_edge_dir/ont_subset.tmp.id
    cat $missing_edge_dir/patch.gapid*gaf > $missing_edge_dir/ont_subset.tmp.gaf
    zcat $verkkodir/3-align/split/ont*.fasta.gz | seqtk subseq - $missing_edge_dir/ont_subset.tmp.id > $missing_edge_dir/ont_subset.tmp.fasta
else 
    echo "ont_subset.tmp.fasta already exists, skipping fasta extraction"
fi

echo "processing 7-consensus directory..."
mkdir -p $newfolder/7-consensus &&
cd $newfolder/7-consensus/ &&
cp $verkkodir/7-consensus/ont_subset.* ./ &&
chmod a+w * &&
cat $missing_edge_dir/ont_subset.tmp.id >> ont_subset.id &&
gunzip ont_subset.fasta.gz &&
cat $missing_edge_dir/ont_subset.tmp.fasta >> ont_subset.fasta &&
bgzip ont_subset.fasta &&
echo "7-consensus directory is updated"
cd ..

echo "processing 6-layoutContigs directory..."
cp -r $verkkodir/6-layoutContigs/ .
chmod -R a+w 6-layoutContigs/ &&
cd 6-layoutContigs/
chmod a+w * &&
rm consensus_paths.txt &&
cat $verkkofilletdir/missing_edge/patch.gapid_*.gaf >> hifi.alignments.gaf &&
cat $verkkofilletdir/missing_edge/patch.gapid_*.gfa | grep "^L" | grep gap >> combined-edges.gfa &&
cat $verkkofilletdir/missing_edge/patch.gapid_*.gfa | grep gap | awk 'BEGIN { FS="[ \t]+"; OFS="\t"; } ($1 == "S") && ($3 != "*") { print $2, length($3); }' >> nodelens.txt &&
cp $finalGaf ./consensus_paths.txt &&
cat $missing_edge_dir/ont_subset.tmp.id >> ont-gapfill.txt &&
rm ./unitig*

echo " "
echo "running replace_path_nodes.py"
echo " "
verkkoLib=$(verkko | grep "Verkko module"| awk '{print $4'})
$verkkoLib/scripts/replace_path_nodes.py ../4-processONT/alns-ont-mapqfilter.gaf combined-nodemap.txt |grep -F -v -w -f ont-gapfill.txt > ont.alignments.gaf &&
cd ..

echo "6-layoutContigs directory is updated!"
