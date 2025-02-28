#! /bin/bash
set -x 
# contig_com=sire_compressed.k31.hapmer-0000228
# verkkodir=/vf/users/Phillippy/projects/giraffeT2T/assembly/verkko2.2_hifi-duplex_trio-hic/verkko-thic_v0.1.0

if [ "$#" -lt 3 ]; then
    echo "Usage: $0 --hpc_fa <fa> --contig <contig> --verkko_dir <verkko_dir> --path <path>"

    echo -e "\t--hpc_fa <fa>: HPC fasta file (required)"
    echo -e "\t--contig <contig>: Contig name (required)"
    echo -e "\t--verkko_dir <verkko_dir>: Verkko directory (required)"
    echo -e "\t--path <path>: Path to the verkko directory (required)"

    exit 1
fi

mkdir -p logs

# get arguments
while [ "$#" -gt 0 ]; do
    case $1 in
        --hpc_fa) hpc_fa=$2; shift 2;;
        --contig) contig=$2; shift 2;;
        --verkko_dir) verkko_dir=$2; shift 2;;
        --path) path=$2; shift 2;;
        *) echo "Unknown parameter passed: $1"; exit 1;;
    esac
done


if [ ! -d $verkko_dir ]; then
    echo "Error: $verkko_dir does not exist"
    exit 1
fi



# print the input args
echo "hpc_fa: $hpc_fa"
echo "contig_com: $contig_com"
echo "verkko_dir: $verkko_dir"
echo "path: $path"

samtools faidx $hpc_fa $contig > $contig.fasta &&
seqtk hpc $contig.fasta > $contig.hpc.fasta &&
rm $contig.fasta


grep "$contig" $path| cut -f 2 | tr ',' '\n' | sed -e 's/-$//'| sed -e 's/+$//' | grep -v "\[" | grep -v "gapmanual"| sort | uniq | tr '\n' ' ' > $contig.utig4.list
samtools faidx $verkko_dir/8-hicPipeline/unitigs.hpc.fasta $(cat $contig.utig4.list) > $contig.unitigs.subset.fasta
seqtk hpc $contig.unitigs.subset.fasta > $contig.unitigs.hpc.subset.fasta &&
rm $contig.unitigs.subset.fasta
samtools faidx $contig.unitigs.hpc.subset.fasta

minimap2 -t 50 -d $contig.hpc.fasta.mmi $contig.hpc.fasta
minimap2 -t 10 -a $contig.hpc.fasta.mmi $contig.unitigs.hpc.subset.fasta > $contig.alignment.sam
samtools view -@ 50 -bS $contig.alignment.sam > $contig.alignment.bam
samtools sort -@ 50 $contig.alignment.bam -o $contig.alignment.sorted.bam &&
rm $contig.alignment.bam
samtools index $contig.alignment.sorted.bam
bedtools bamtobed -i $contig.alignment.sorted.bam > $contig.alignment.sorted.bed