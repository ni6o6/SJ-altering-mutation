#!/bin/bash

#Naoko Iida

ref="s3://niida-tokyo/lung_wg/hg19_0chr.fa"
fin="26merge_fil.annot.ass.adj.freq.AI-mut.cano.txt"
fo='26merge_fil.annot.ass.adj.freq.AI-mut.cano.mpileup.txt'

sed 's/ /_/g' ${fin} > tmp.txt

while IFS= read -r line; do

F=($(echo ${line} | tr '\t' ' '))
echo ${F[3]}
input="s3://niida-tokyo/lung_wg/"${F[3]}".markdup.bam"
pos=${F[10]}":"${F[11]}"-"${F[11]}
echo ${pos}
samtools mpileup -r ${pos} -f ${ref} ${input} > mpileup.txt
    while IFS= read -r row ; do
    echo ${F[3]}\t${row} >> ${fo}
    done < mpileup.txt
done < tmp.txt

python SJ_annot_mpileup.py
