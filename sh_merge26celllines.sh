#!/bin/bash

echo "chr\tstart\tend\tsample\tclass\tstrand\treads\treads.1\ttotal\tfreq\tCHR\tPOS\tID\tREF_x\tMUT_x\tSCORE1\tSCORE2\tSYMBOL\tSTRAND\tTYPE\tDIST\tDS_AG\tDS_AL\tDS_DG\tDS_DL\tDP_AG\tDP_AL\tDP_DG\tDP_DL\tSAVnet_SpliceAI\tREF_y\tMUT_y"> merge26.SJ.fil.annot.ass.adj.freq.AI-mut_merge.txt

for file in `\find . -name 'A427.SJ.fil.annot.ass.adj.freq.AI-mut.txt'`; do
	echo ${file}
	pr=($(echo ${file} | tr '.' ' '))
	echo ${pr[0]}
	awk '{if($9>60 && $10>=0.05) print $0}' ${file} >> merge26.SJ.fil.annot.ass.adj.freq.AI-mut_merge.txt
	
done 
