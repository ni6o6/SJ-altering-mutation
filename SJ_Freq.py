""""""
#SJ_20190501_Freq.py
""""""
import numpy as np
import pandas as pd
import glob
import os

fo1='./alterativeSJ_adjustment/'
fo2='./alterativeSJ_freq/'

in1=glob.glob(fo1+"*.SJ.filtered.annot.ass.adj.intersect.bed")

for infile1 in in1:
    basename = os.path.basename(infile1)
    pr=basename.split('.', 8)[0]
    print(pr)

    df1 = pd.read_csv(infile1,sep='\t',header=None, index_col=None,names=('A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q') )

    ta1=df1.groupby(["A","B","C"]).agg({'G': np.max, 'N': np.sum}).reset_index()
    ratio=ta1["G"]/ta1["N"]
    ta2 = pd.concat([ta1,ratio], axis=1, join='inner')

    #out1='A427.SJ.filtered.annot.ass.adj.freq.txt'
    infile2 ='%s.SJ.filtered.annot.ass.adj.bed' %(pr)
    df2 = pd.read_csv(fo1+infile2,sep='\t',header=None, index_col=None,names=('A','B','C','D','E','F','G'))

    rec1=pd.merge(df2, ta2, on=['A', 'B', 'C'])
    col = ["chr","start","end","sample","class","strand","reads","reads","total","freq"] #keep two "reads" because the multiple rows of same SSs are present.
    rec1.columns = col
    out1='%s.SJ.filtered.annot.ass.adj.freq.txt' %(pr)
    rec1.to_csv(fo2+out1, index=False, sep='\t')
