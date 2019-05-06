#! /usr/bin/env python
"""
Naoko Iida
#source activate py27
"""

import glob
import os
import subprocess
import pandas as pd
import numpy as np

cntl="./reference/SJ_control_2_4.bed.gz"
in1=glob.glob("./junction/*.SJ.out.tab")

for filepath in in1:
    basename = os.path.basename(filepath)
    pr=basename.split('.', 4)[0]
    print(pr)
    out1='./alterativeSJ_filtered_annot2/%s.SJ.filtered.out.tab' %(pr)
    out2='./alterativeSJ_filtered_annot2/%s.SJ.filtered.out.tab' %(pr)
    #junc_utils filter
    #filter_commands = ["junc_utils", "filter", "--pooled_control_file", cntl, filepath, pr+".SJ.filtered.out.tab"]
    filter_commands = ["junc_utils", "filter", "--pooled_control_file", cntl, filepath, out1]
    subprocess.check_call(filter_commands)
    #junc_utils annotate *.SJ.filtered.out.tab
    #annotate_commands = ["junc_utils", "annotate", pr+".SJ.filtered.out.tab", pr+".SJ_filtered.annot.txt", "--genome_id", "hg19"]
    annotate_commands = ["junc_utils", "annotate", out1, out2, "--genome_id", "hg19"]
    subprocess.call(annotate_commands)
    
    #Add new colmun
    #df = pd.read_table(pr+'.SJ_filtered.annot.txt', sep='\t')
    #n=len(df.index)
    #df['Samplename'] = np.array([pr]*n)
    #df.to_csv('26lungcell.SJ_filtered.annot.txt', index=None, sep='\t', mode='a', header=False)
    #df.to_csv('catall.SJ_filtered.annot.txt', index=None, sep='\t', mode='a', header=df.columns)




