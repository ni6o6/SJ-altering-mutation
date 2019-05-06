#!/bin/bash

#mkdir alterativeSJ_filtered_annot
#mkdir alterativeSJ_adjustment
#mkdir alterativeSJ_freq
#mkdir alterativeSJ_tabixAI
mkdir alterativeSJ_add_mut

source activate py27
#python SJ_FilterAnno.py

#python SJ_AssAdj2bed.py
#done

#sh SJ_IntersectDepth.sh #need bedtools

#python SJ_Freq.py

#python SJ_tabix.py

python SJ_add_mut.py
