

$script/vgp-assembly/find assembly.fasta > tmp.telomere

java FindTelomereWindows tmp.telomere 99.9 |awk '{if ($NF > 0.5) print $2"\t"$4"\t"$5"\t"$3"\t"$NF}'|sed s/\>//g|bedtools merge -d -500 -i - -c 4 -o distinct

$script/vgp-assembly/telomere/telomere_analysis.sh giraffe 0.5 50000 assembly.fasta giraffe
