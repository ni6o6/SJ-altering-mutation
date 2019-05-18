#! /usr/bin/env python
"""
Naoko Iida

"""
import numpy as np
import pandas as pd
import glob
import os

ann='./alterativeSJ_add_mut/26merge_fil.annot.ass.adj.freq.AI-mut.cano.mpileup.txt'
fi1='./alterativeSJ_add_mut/26merge_fil.annot.ass.adj.freq.AI-mut.cano.txt'
fo1='./alterativeSJ_add_mut/26merge_fil.annot.ass.adj.freq.AI-mut.cano.mpileup2_20190515.txt'

df1 = pd.read_csv(ann,sep='\t',header=None, index_col=None, dtype = 'object')
df1.columns = ['sample', 'CHR', 'POS', 'ref', 'reads', 'info1', 'info2']
print(df1.head())
#merge
# 
df2 = pd.read_csv(fi1,sep='\t' ,header=0, index_col=None,dtype = 'object')
#print(df2.head())
rec1=pd.merge(df2, df1, on=['sample', 'CHR', 'POS'],how='left') 
#output
rec1.to_csv(fo1, index=False, sep='\t')
