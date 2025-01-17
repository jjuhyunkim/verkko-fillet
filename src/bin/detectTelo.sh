#! /bin/bash

fasta=$1
prefix=$2

# /data/korens/devel/utils/telomere/find $fasta > tmp.telomere

# java FindTelomereWindows tmp.telomere 99.9 |awk '{if ($NF > 0.5) print $2"\t"$4"\t"$5"\t"$3"\t"$NF}'|sed s/\>//g|bedtools merge -d -500 -i - -c 4 -o distinct

$VGP_PIPELINE/telomere/telomere_analysis.sh $prefix 0.5 50000 $fasta

