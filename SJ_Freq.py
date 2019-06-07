#! /usr/bin/env python

import numpy as np
import pandas as pd
import glob
import os
import time

fo1='./alterativeSJ_adjustment/'
fo2='./alterativeSJ_freq/'
fo3='./junction/'

in1=glob.glob(fo1+"*.SJ.filtered.annot.ass.adj.bed")
#in1=glob.glob(fo1+"DRR016718_RERF-LC-OK.SJ.filtered.annot.ass.adj.bed")

for infile1 in in1:
    start = time.time()
    basename = os.path.basename(infile1)
    pr=basename.split('.', 7)[0]
    print(pr)
    outfile='./alterativeSJ_freq/%s.SJ.filtered.annot.ass.adj.freq.txt' %(pr)
    jfile= fo3 + '%s.SJ.out.tab' %(pr)
    jf = pd.read_csv(jfile, sep='\t', header=None, index_col=None,names=('CHR','START','END','A','B','C','reads','D','E'),
                 dtype = {'CHR':'object', 'START':'object', 'END':'object','A':'object','B':'object','C':'object','reads':'int','D':'object','E':'object'})
    with open(infile1, 'r') as df1:
        with open(outfile, 'w') as out1:
            for line in df1:
                F = line.rstrip('\n').split('\t')
                ln = line.rstrip('\n')
                if '5' in F[6] and '+' in F[7]:
                    c = (F[0],)
                    cl = list(c)
                    p = (F[4], F[2])
                    pl = list(p)
                    rec = jf[(jf['CHR'].isin(cl) ) & (jf['END'].isin(pl))]  
                    #total = rec.groupby('END').sum()['reads'] <--XX
                    total = int(rec.sum()['reads'])
                    freq = round(int(F[8])/total, 3)
                    out1.write(ln+'\t'+str(total)+'\t'+str(freq)+'\n')
                elif '5' in F[6] and '-' in F[7]:
                    c = (F[0],)
                    cl = list(c)
                    p = (F[3], F[1])
                    pl = list(p)
                    rec = jf[(jf['CHR'].isin(cl) ) & (jf['START'].isin(pl))]  
                    total = int(rec.sum()['reads'])
                    freq = round(int(F[8])/total, 3)
                    out1.write(ln+'\t'+str(total)+'\t'+str(freq)+'\n')
                elif '3' in F[6] and '+' in F[7]:
                    c = (F[0],)
                    cl = list(c)
                    p = (F[3], F[1])
                    pl = list(p)
                    rec = jf[(jf['CHR'].isin(cl) ) & (jf['START'].isin(pl))]  
                    total = int(rec.sum()['reads'])
                    freq = round(int(F[8])/total, 3)
                    out1.write(ln+'\t'+str(total)+'\t'+str(freq)+'\n')
                elif '3' in F[6] and '-' in F[7]:
                    c = (F[0],)
                    cl = list(c)
                    p = (F[4], F[2])
                    pl = list(p)
                    rec = jf[(jf['CHR'].isin(cl) ) & (jf['END'].isin(pl))]  
                    total = int(rec.sum()['reads'])
                    freq = round(int(F[8])/total, 3)
                    out1.write(ln+'\t'+str(total)+'\t'+str(freq)+'\n')
    end = time.time()
    print(round(end-start,2))
