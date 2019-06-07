#! /usr/bin/env python
"""
Naoko Iida
#source activate py27
"""
import glob
import os
import pysam
import sys, subprocess
import pandas as pd
import numpy as np
import re

tbx = pysam.TabixFile("./reference/whole_genome_filtered_spliceai_scores.vcf.gz")

in1=glob.glob("./alterativeSJ_freq/*.SJ.filtered.annot.ass.adj.freq.txt")

for infile in in1:
    basename = os.path.basename(infile)
    pr=basename.split('.', 8)[0]
    print(pr)
    outfile='./alterativeSJ_tabixAI/%s.SJ.fil.annot.ass.adj.freq.AI.txt' %(pr)

    with open(infile, 'r') as in2:
        with open(outfile, 'w') as out1:
            #he = 'chr\tstart\tend\tstart_ori\tend_ori\tsample\tclass\tstrand\treads\ttotal\tfreq\tCHR\tPOS\tID\tREF\tMUT\tSCORE1\tSCORE2\tSYMBOL\tSTRAND\tTYPE\tDIST\tDS_AG\tDS_AL\tDS_DG\tDS_DL\tDP_AG\tDP_AL\tDP_DG\tDP_DL\tSAVnet_SpliceAI\n' #keep two "reads" because the multiple rows of same SSs are present.
            #out1.write(he)

            for line in in2:
                F = line.rstrip('\n').split('\t')
                ln = line.rstrip('\n')
                
                if F[0] in ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]:
                    if "3" in F[6] and "+" in F[7]:
                        c=F[0]
                        s=int(F[2])-5
                        e=int(F[2])+5
                        rows = tbx.fetch(c, s, e)
                        num=0
                        for row in rows:
                            #out1.write(ln+'\t'+str(row)+'\n')
                            line=ln+'\t'+str(row)
                            G = line.split('\t')
                            cl = G[4] #class
                            ints = int(G[1]) #SAVnet intron start
                            inte = int(G[2])
                            aip = int(G[12]) #Splice AI postion
                            inf = G[18]
                            inf2 = re.split('[=;]', inf)
                            strand = G[7]
                            ag = int(inf2[17])
                            dg = int(inf2[21])
                            del G[-1]
                            line= '\t'.join(G)
                            G_list=(inf2[1],inf2[3],inf2[5],inf2[7],inf2[9],inf2[11],inf2[13],inf2[15],inf2[17],inf2[19],inf2[21],inf2[23])
                            line2= '\t'.join(G_list)
                            if inte == aip+ag-1:
                                num+=1
                                out1.write(line+'\t'+line2+'\tc'+str(num)+'\n')
                            else:
                                out1.write(line+'\t'+line2+'\t-\n')
                    elif "3" in F[6] and "-" in F[7]:
                        c=F[0]
                        s=int(F[1])-5
                        e=int(F[1])+5
                        rows = tbx.fetch(c, s, e)
                        num=0
                        for row in rows:
                            #out1.write(ln+'\t'+str(row)+'\n')
                            line=ln+'\t'+str(row)
                            G = line.split('\t')
                            cl = G[4] #class
                            ints = int(G[1]) #SAVnet intron start
                            inte = int(G[2])
                            aip = int(G[12]) #Splice AI postion
                            inf = G[18]
                            inf2 = re.split('[=;]', inf)
                            strand = G[7]
                            ag = int(inf2[17])
                            dg = int(inf2[21])
                            del G[-1]
                            line= '\t'.join(G)
                            G_list=(inf2[1],inf2[3],inf2[5],inf2[7],inf2[9],inf2[11],inf2[13],inf2[15],inf2[17],inf2[19],inf2[21],inf2[23])
                            line2= '\t'.join(G_list)
                            if ints == aip+ag+1:
                                num+=1
                                out1.write(line+'\t'+line2+'\tc'+str(num)+'\n')
                            else:
                                out1.write(line+'\t'+line2+'\t-\n')
                    elif "5" in F[6] and "-" in F[7]:
                        c=F[0]
                        s=int(F[2])-5
                        e=int(F[2])+5
                        rows = tbx.fetch(c, s, e)
                        num=0
                        for row in rows:
                            #out1.write(ln+'\t'+str(row)+'\n')
                            line=ln+'\t'+str(row)
                            G = line.split('\t')
                            cl = G[4] #class
                            ints = int(G[1]) #SAVnet intron start
                            inte = int(G[2])
                            aip = int(G[12]) #Splice AI postion
                            inf = G[18]
                            inf2 = re.split('[=;]', inf)
                            strand = G[7]
                            ag = int(inf2[17])
                            dg = int(inf2[21])
                            del G[-1]
                            line= '\t'.join(G)
                            G_list=(inf2[1],inf2[3],inf2[5],inf2[7],inf2[9],inf2[11],inf2[13],inf2[15],inf2[17],inf2[19],inf2[21],inf2[23])
                            line2= '\t'.join(G_list)
                            if inte == aip+dg-1:
                                num+=1
                                out1.write(line+'\t'+line2+'\tc'+str(num)+'\n')
                            else:
                                out1.write(line+'\t'+line2+'\t-\n')
                    elif "5" in F[6] and "+" in F[7]:
                        c=F[0]
                        s=int(F[1])-5
                        e=int(F[1])+5
                        rows = tbx.fetch(c, s, e)
                        num=0
                        for row in rows:
                            #out1.write(ln+'\t'+str(row)+'\n')
                            line=ln+'\t'+str(row)
                            G = line.split('\t')
                            cl = G[4] #class
                            ints = int(G[1]) #SAVnet intron start
                            inte = int(G[2])
                            aip = int(G[12]) #Splice AI postion
                            inf = G[18]
                            inf2 = re.split('[=;]', inf)
                            strand = G[5]
                            ag = int(inf2[17])
                            dg = int(inf2[21])
                            del G[-1]
                            line= '\t'.join(G)
                            G_list=(inf2[1],inf2[3],inf2[5],inf2[7],inf2[9],inf2[11],inf2[13],inf2[15],inf2[17],inf2[19],inf2[21],inf2[23])
                            line2= '\t'.join(G_list)
                            if ints == aip+dg+1:
                                num+=1
                                out1.write(line+'\t'+line2+'\tc'+str(num)+'\n')
                            else:
                                out1.write(line+'\t'+line2+'\t-\n')
#https://pysam.readthedocs.io/en/latest/usage.html#working-with-tabix-indexed-files


