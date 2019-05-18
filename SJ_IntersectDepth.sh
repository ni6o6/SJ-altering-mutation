#!/bin/bash

#cd ./alterativeSJ_adjustment
for file in `\find ./alterativeSJ_adjustment -name '*SJ.filtered.annot.ass.adj.bed'`; do
	echo ${file}
	pr=($(echo ${file} | tr '.|/' ' '))
	echo ${pr[1]}
	b="./junction/${pr[1]}.SJ.out.tab"
	c="${pr[1]}.SJ.filtered.annot.ass.adj.intersect.bed"
	intersectBed -wao -a ${file} -b ${b} > ./alterativeSJ_adjustment/${c} #result.bed
	
	
done 
