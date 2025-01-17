#! /bin/bash

input_fofn=$1

LEN=`wc -l $input_fofn | awk '{print $1}'`

### Count
$MERQURY/build/count.sh


$MERQURY/build/union_sum.sh


### Union

for line_num in $(seq 1 $LEN)
do
  input=`sed -n ${line_num}p $input_fofn`
  name=`echo $input | sed 's/.fastq.gz$//g' | sed 's/.fq.gz$//g' | sed 's/.fasta$//g' | sed 's/.fa$//g' | sed 's/.fasta.gz$//g' | sed 's/.fa.gz$//g'`
  name=`basename $name`
  echo "$name.k$k.$line_num.meryl" >> $out_prefix.meryl_count.meryl.list
done







