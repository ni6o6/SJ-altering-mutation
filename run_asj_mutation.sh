#!/bin/bash

#python3

mkdir alterativeSJ_filtered_annot
mkdir alterativeSJ_adjustment
mkdir alterativeSJ_freq
mkdir alterativeSJ_tabixAI
mkdir alterativeSJ_add_mut

cntl1=./reference/$1
#cntl1="./reference/SJ_control_2_4.bed.gz"
#cntl2="./reference/SJ_26common2.bed.gz" #savnetSJ
cntl2=./reference/$2
#cntl2="./reference/26common2_hg19_20190606v2.bed.gz"
#cntl2="./reference/26common2_recount2hg19.bed.gz" #recount2SJ

for file in `\find ./junction -name '*.SJ.out.tab'`; do

arr=(`echo $file | tr -s '.|/' ' '`)
pr=${arr[1]}
echo ${arr[1]}
out1="./alterativeSJ_filtered_annot/"${pr}".SJ.filx1.txt"
out2="./alterativeSJ_filtered_annot/"${pr}".SJ.filx2.txt"
out3="./alterativeSJ_filtered_annot/"${pr}".SJ.filx2.annot.txt"

junc_utils filter --pooled_control_file ${cntl1} ${file} ${out1}

#junc_utils filter --pooled_control_file ${cntl2} ${out1} ${out2}

junc_utils annotate ${out2} ${out3} --genome_id hg19

done

python SJ_AssAdj2bed.py

python SJ_Freq.py

python SJ_tabix_hg19.py

python SJ_extract_canonicalMut.py
		
#after mpileup
#python SJ_annot_mpileup.py

