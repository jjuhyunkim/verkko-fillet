#! /bin/bash

echo "Usage: _node_to_assembly.sh <assembly> <mapping> <path_file> <gfa> <name of contig>"
echo "assembly: assembly file in fasta format(uncompressed, final fasta)"
echo "mapping: mapping file with contig names, new name for each contig(eg. chr1) and name of path for each contig"
echo "path_file: file with path information"
echo "gfa: graph with sequences of unitigs, with compressed"
echo "name of contig: name of the contig to be processed"
echo " "

# example of map file 
# column 1 : name of the contig
# column 2 : name of the output chromosome
# column 3 : name of the path for each contig
#dam_compressed.k31.hapmer-0000001       CM029953.1      dam_compressed.k31.hapmer_from_utig4-1003
#dam_compressed.k31.hapmer-0000002       CM029959.1      dam_compressed.k31.hapmer_from_utig4-1104
#dam_compressed.k31.hapmer-0000003       CM029950.1      dam_compressed.k31.hapmer_from_utig4-1243
#dam_compressed.k31.hapmer-0000004       CM029954.1      dam_compressed.k31.hapmer_from_utig4-1282
#dam_compressed.k31.hapmer-0000005       CM029958.1      dam_compressed.k31.hapmer_from_utig4-1315
#dam_compressed.k31.hapmer-0000006       CM029952.1      dam_compressed.k31.hapmer_from_utig4-1425


echo "Loading modules"
ml minimap2
ml bedtools
ml samtools

asm=$(realpath $1)
map=$(realpath $2)
path_file=$(realpath $3)
gfa=$(realpath $4)
name=$5

exc="random"
threads=10

printf "Assembly: $asm\nMapping: $map\nPath file: $path_file\ngfa: $gfa\nthreads: $threads\n"

# check the inputs
if [ ! -f $asm ]; then
echo "Assembly file not found"
exit 1
fi

if [ ! -f $map ]; then
echo "Mapping file not found"
exit 1
fi

if [ ! -f $path_file ]; then
echo "Path file not found"
exit 1
fi

if [ ! -f $gfa ]; then
echo "GFA file not found"
exit 1
fi

if [ ! -f $name.alignment.sorted.bed ]; then
echo "Processing $name"

samtools faidx $asm $name > $name.fa &&
samtools faidx $name.fa &&

path=$(fgrep -w $name $map | cut -f 3)
echo -e "$name\t$path\n"
# grep $path $path_file | cut -f 2 | tr ',' '\n' | sed -e 's/-$//'| sed -e 's/+$//' | grep -v "\["  > $name.utig4.list &&
grep $path $path_file | cut -f 2 | sed 's/>/\n/g' | sed 's/,/\n/g' | sed 's/</\n/g' | sed -e 's/\[.*//g' | sed -e '/^$/d' | sed 's/[-+]$//' > $name.utig4.list &&
awk 'NR==FNR {list[$1]; next} $1 == "S" && $2 in list {print ">"$2; print $3}' "$name.utig4.list" "$gfa" > $name.unitigs.hpc.subset.fasta &&
minimap2 -t $threads -H -d $name.mmi $name.fa &&
minimap2 -t $threads -a $name.mmi $name.unitigs.hpc.subset.fasta > $name.alignment.sam &&
rm $name.unitigs.hpc.subset.fasta &&
rm $name.mmi &&
rm $name.fa &&
samtools view -bS $name.alignment.sam > $name.alignment.bam && 
rm $name.alignment.sam &&
samtools sort -@ $threads $name.alignment.bam -o $name.alignment.sorted.bam && 
rm $name.alignment.bam &&
samtools index $name.alignment.sorted.bam && 
bedtools bamtobed -i $name.alignment.sorted.bam > $name.alignment.sorted.bed &&
rm $name.alignment.sorted.bam
fi

echo "Done"