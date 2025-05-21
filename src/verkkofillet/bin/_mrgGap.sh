#! /bin/bash

# Inputs --------------------------------------------------------------------
bed=$1
original_fa=$2
output_fa=$3
gaplen=$4

# Check if the input files are provided --------------------------------------------------------------------
if [ $# -ne 4 ]; then
    echo " "
    echo -e "verkko-fillet/src/verkkofillet/bin/_mrgGap.sh: Merge gaps in a fasta file using a bed file.\n"
    echo -e "Usage: $0 <bed_file> <input_fa> <output_fa> <gap_length>\n"
    echo -e "\tExample: $0 input.bed original.fasta output.fasta gaplength.txt\n"

    echo -e "\tbed_file\tA bed file with the gaps to be merged. The first three columns should be the chromosome, start, and end positions of the gaps."
    echo -e "\tinput_fa\tThe original fasta file."
    echo -e "\toutput_fa\tThe output fasta file with merged gaps."
    echo -e "\tgap_length\tA file with the gap lengths. The first column should be the gap ID, the second column should be the gap length.\n"

    #echo -e "\tExample gap_length_file content:"
    #echo -e "\tmrg1\t1000000\nmrg2\t\nmrg3\t\nmrg4\t4000000\n"
    # echo "The second column of the gap_length_file can be empty, in which case it will be set to 1000000."

    #echo -e "\tExample of the bed file content:"
    #echo -e "\tchr14_mat\t0\t1000000\tmrg1\nchr14_mat\t1000000\t2000000\tmrg2\nchr14_mat\t2000000\t3000000\tmrg3"
    # echo "The bed file should contain the chromosome, start, end positions, and gap IDs."
    exit 1
fi


# Load the modules --------------------------------------------------------------------
ml samtools
ml bedtools
ml seqtk

# Check files --------------------------------------------------------------------
fai=$original_fa.fai
if [ ! -f "$fai" ]; then
    samtools faidx $original_fa
fi
# Check the inputs
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

# Check the bed file and Join the files -------------------------------------------------------------------- 
gaps=$(awk '$4!=""{print $4}' $bed | sort | uniq | tr '\n' ', ' | sed 's/,$/\n/')
echo -e "Contigs: $gaps\n"
contigs=$(awk '$4!=""{print $1}' $bed | sort | uniq | tr '\n' ', ' | sed 's/,$/\n/')
echo -e "Contigs: $contigs\n"

awk '{
    if ($4 != "") {
        group = $4
        chr[group] = $1  # store chromosome
        if (!(group in min2) || $2 < min2[group]) min2[group] = $2
        if (!(group in max3) || $3 > max3[group]) max3[group] = $3
    }
}
END {
    for (g in min2)
        print chr[g], min2[g], max3[g], g
}' OFS='\t' "$bed" | bedtools sort > gap.mrg.bed

join -a 1 -a 2 -1 4 -2 1 -e "1000000" -o 1.1 1.2 1.3 1.4 2.2 gap.mrg.bed gapLength.txt  | tr ' ' '\t' > gapName.gapLength.outer_joined.txt &&
rm gap.mrg.bed

# Extract unplaced regions --------------------------------------------------------------------
awk '{print $1,$2,$3 }' OFS='\t' gapName.gapLength.outer_joined.txt > gapName.gapLength.outer_joined.tmp.bed
awk '$4!=""{print $1,$2,$3 }' OFS='\t' $bed > $bed.tmp.bed

rm -rf unplaced.bed
bedtools subtract -a gapName.gapLength.outer_joined.tmp.bed -b $bed.tmp.bed | awk '{print $0,$1"_"$2"_"$3}' OFS='\t' > unplaced.bed
rm $bed.tmp.bed gapName.gapLength.outer_joined.tmp.bed

# Make not gap regions --------------------------------------------------------------------
for contig in `cut -f 1 gapName.gapLength.outer_joined.txt | sort -u`; do
echo "Processing contig: $contig"


grep -w $contig $fai | awk '{print $1, 0, $2}' OFS='\t'> fai.tmp.bed
grep -w $contig gapName.gapLength.outer_joined.txt > gapName.gapLength.outer_joined.tmp.txt
bedtools subtract -a fai.tmp.bed -b gapName.gapLength.outer_joined.tmp.txt > $contig.noGap.bed
bedtools getfasta -fi $original_fa -bed $contig.noGap.bed -fo $contig.noGap.fa
samtools faidx $contig.noGap.fa


numGap=$(wc -l < "$contig.noGap.fa.fai" | awk '{print $1 - 1}')
echo "Number of gaps: $numGap"
for ((gapnum = 1; gapnum <= numGap; gapnum++)); do
echo "Processing gap: $gapnum"
if [ $gapnum -eq 1 ]; then
    samtools faidx $contig.noGap.fa $(sed -n "${gapnum}p" $contig.noGap.fa.fai | awk '{print $1}') > $contig.noGap_wNewGap_1.fa
else
    mv $contig.noGap_wNewGap.fa $contig.noGap_wNewGap_1.fa
fi
line=$((gapnum + 1))
samtools faidx "$contig.noGap.fa" $(sed -n "${line}p" "$contig.noGap.fa.fai" | awk '{print $1}') > "$contig.noGap_wNewGap_2.fa"
lenGap=$(sed -n "${gapnum}p" gapName.gapLength.outer_joined.txt | awk '{print $5}')
{
  echo ">${contig}_withNewGap"
  seqtk seq $contig.noGap_wNewGap_1.fa | sed -e '1d'
  printf "%0.sN" $(seq 1 $lenGap)
  seqtk seq $contig.noGap_wNewGap_2.fa | sed -e '1d'
} | seqtk seq | fold -w 60 > $contig.noGap_wNewGap.fa
done
samtools faidx $contig.noGap_wNewGap.fa 
rm $contig.noGap_wNewGap_1.fa* $contig.noGap_wNewGap_2.fa* $contig.noGap.fa*

awk -v contig=${contig}_withNewGap 'BEGIN {FS=""; OFS="\t"} 
    /^>/ {if (seq) {for (i=1; i<=length(seq); i++) if (substr(seq, i, 1) == "N") print contig, start+i-1, start+i} 
          chrom=$1; start=0; seq=""} 
    {if ($0 !~ /^>/) seq=seq$0} 
    END {if (seq) {for (i=1; i<=length(seq); i++) if (substr(seq, i, 1) == "N") print contig, start+i-1, start+i}}' $contig.noGap_wNewGap.fa > $contig.noGap_wNewGap.tmp.bed
bedtools merge -i  $contig.noGap_wNewGap.tmp.bed > $contig.noGap_wNewGap.bed && rm  $contig.noGap_wNewGap.tmp.bed
done

# Finalize --------------------------------------------------------------------
echo -e "Merging the final fasta"

rm -rf $output_fa
for contig in `cut -f 1 $original_fa.fai`; do
echo "Processing contig: $contig"
if [ ! -f $contig.noGap_wNewGap.fa ]; then
    samtools faidx $original_fa $contig | seqtk seq | fold -w 60  >> $output_fa
else
    cat $contig.noGap_wNewGap.fa | seqtk seq | fold -w 60 >> $output_fa
    rm $contig.noGap_wNewGap.fa*
fi
done
bedtools getfasta -fi $original_fa -bed unplaced.bed -name | seqtk seq | fold -w 60 >> $output_fa && samtools faidx $output_fa
samtools faidx $output_fa
rm *.noGap.bed 