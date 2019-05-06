#! /usr/bin/env python
"""
Naoko Iida

"""
import numpy as np
import pandas as pd
import glob
import os

fo1='./reference/'
fo2='./alterativeSJ_tabixAI/'
fo3='./alterativeSJ_add_mut/'
#mutations file from SAVnet.result 
f1='lung_cellline.savnet.result.txt'

df1 = pd.read_csv(fo1+f1,sep='\t',header=0, index_col=None)

in1=glob.glob(fo2+"*.SJ.fil.annot.ass.adj.freq.AI.txt")
for f in in1:
    basename = os.path.basename(f)
    pr=basename.split('.', 9)[0]
    print(pr)

#Extract row, sample=pr
#   #open file2 
    df2 = pd.read_csv(f,sep='\t' ,header=0, index_col=None,dtype = 'object')
    df1s=df1[df1['Sample_Name'].isin([pr])]
    mt_key=df1s.iloc[:,2]
    mt=mt_key.str.split(',', expand=True)
    mt.rename(columns={0: 'CHR', 1: 'POS', 2: 'REF', 3: 'MUT'}, inplace=True)
    #df2['CHR']=df2['CHR'].astype(str)
    rec1=pd.merge(df2, mt, on=['CHR', 'POS'],how='left') 
#output
    out1=pr+'.SJ.fil.annot.ass.adj.freq.AI-mut.txt'    
    rec1.to_csv(fo3+out1, index=False, sep='\t')
