#! /bin/bash

ml verkko

verkkofilletdir=$1
verkkodir=$2
newfolder=$3
finalGaf=$4
subsetid=$5
subsetfasta=$6

finalGaf=$(realpath $finalGaf)
newfolder=$(realpath $newfolder)
verkkofilletdir=$(realpath $verkkofilletdir)
verkkodir=$(realpath $verkkodir)
subsetid=$(realpath $subsetid)
subsetfasta=$(realpath $subsetfasta)


echo -e "finalGaf : $finalGaf"
echo -e "newfolder : $newfolder"
echo -e "verkkofilletdir : $verkkofilletdir"
echo -e "verkkodir : $verkkodir"
echo -e "subsetid : $subsetid"
echo -e "subsetfasta : $subsetfasta"


if ! command -v verkko &> /dev/null
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
cat $subsetid >> ont_subset.id &&
gunzip ont_subset.fasta.gz &&
cat $subsetfasta >> ont_subset.fasta &&
bgzip ont_subset.fasta &&
echo "7-consensus directory is updated"
cd ..

echo "processing 6-layoutContigs directory..."
cp -r $verkkodir/6-layoutContigs/ .
chmod -R a+w 6-layoutContigs/ &&
cd 6-layoutContigs/
chmod a+w * &&
rm consensus_paths.txt &&
cat $verkkofilletdir/missing_edge/patch.gapid_*.gaf >> combined-alignments.gaf &&
cat $verkkofilletdir/missing_edge/patch.gapid_*.gfa |grep "^L" |grep gap >> combined-edges.gfa &&
cat $verkkofilletdir/missing_edge/patch.gapid_*.gfa| grep gap | awk 'BEGIN { FS="[ \t]+"; OFS="\t"; } ($1 == "S") && ($3 != "*") { print $2, length($3); }' >> nodelens.txt &&
cp $finalGaf ./consensus_paths.txt &&
cat $subsetid >> ont-gapfill.txt &&

echo " "
echo "running replace_path_nodes.py"
echo " "
verkkoLib=$(verkko | grep "Verkko module"| awk '{print $4'})
$verkkoLib/scripts/replace_path_nodes.py ../4-processONT/alns-ont-mapqfilter.gaf combined-nodemap.txt |grep -F -v -w -f ont-gapfill.txt > ont.alignments.gaf &&
cd ..

echo "6-layoutContigs directory is updated!"