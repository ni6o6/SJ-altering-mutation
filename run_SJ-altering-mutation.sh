#!/bin/bash

#mkdir alterativeSJ_filtered_annot
#mkdir alterativeSJ_adjustment
#mkdir alterativeSJ_freq
#mkdir alterativeSJ_tabixAI
#mkdir alterativeSJ_add_mut

#source activate py27
#cntl1="./reference/SJ_control_2_4.bed.gz"
#cntl2="./reference/SJ_26common2.bed.gz"

#for file in `\find ./junction -name '*.SJ.out.tab'`; do

#arr=(`echo $file | tr -s '.|/' ' '`)
#pr=${arr[1]}
#echo ${arr[1]}
#out1="./alterativeSJ_filtered_annot/"${pr}".SJ.filx1.txt"
#out2="./alterativeSJ_filtered_annot/"${pr}".SJ.filx2.txt"
#out3="./alterativeSJ_filtered_annot/"${pr}".SJ.filx2.annot.txt"

#junc_utils filter
#junc_utils filter --pooled_control_file ${cntl1} ${file} ${out1}
#junc_utils filter
#junc_utils filter --pooled_control_file ${cntl2} ${out1} ${out2}

#junc_utils annotate *.SJ.filtered.out.tab
#junc_utils annotate ${out2} ${out3} --genome_id hg19

#done

#python SJ_AssAdj2bed.py

#sh SJ_IntersectDepth.sh #need bedtools

#python SJ_Freq.py

#python SJ_tabix.py

python SJ_add_mut.py

sh_merge26celllines.sh
