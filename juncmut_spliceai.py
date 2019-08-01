#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 2019

@author: genome

#def asj_tabix_hg19(input_pr, folder):

"""

import pysam
import argparse
import re


parser = argparse.ArgumentParser() #make a parser

parser.add_argument("input", metavar = "./junction/sample.SJ.out.tab", default = None, type = str,
                        help = "Path to input file") 
    
parser.add_argument("folder", metavar = "group", default = "my_samples", type = str,
                        help = "Path to input file") 
    
args = parser.parse_args()

pr = args.input
print(pr) 
file44 ='./alterativeSJ_assadjfreq/%s/%s.SJ.fil.annot.assadjunifreqT.txt' %(args.folder,pr)
file5 = './alterativeSJ_tabixAI/%s/%s.SJ.fil.annot.assadjunifreqT.AI.txt' %(args.folder,pr)
e_file = './alterativeSJ_tabixAI/%s/%s.SJ.fil.annot.assadjunifreqT.noAI.txt' %(args.folder,pr) 
tbx = pysam.TabixFile("./reference/whole_genome_filtered_spliceai_scores.vcf.gz")


with open(file44, 'r') as in1:
    with open(file5, 'w') as out1:
        for line in in1:
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
                        line=ln+'\t'+str(row)
                        G = line.split('\t')
                        ints = int(G[1]) #SAVnet intron start
                        inte = int(G[2])
                        aip = int(G[12]) #Splice AI postion
                        inf = G[18]
                        inf2 = re.split('[=;]', inf)
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
                        line=ln+'\t'+str(row)
                        G = line.split('\t')
                        ints = int(G[1]) #SAVnet intron start
                        inte = int(G[2])
                        aip = int(G[12]) #Splice AI postion
                        inf = G[18]
                        inf2 = re.split('[=;]', inf)
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
                        line=ln+'\t'+str(row)
                        G = line.split('\t')
                        ints = int(G[1]) #SAVnet intron start
                        inte = int(G[2])
                        aip = int(G[12]) #Splice AI postion
                        inf = G[18]
                        inf2 = re.split('[=;]', inf)
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
                        line=ln+'\t'+str(row)
                        G = line.split('\t')
                        #cl = G[4] #class
                        ints = int(G[1]) #SAVnet intron start
                        inte = int(G[2])
                        aip = int(G[12]) #Splice AI postion
                        inf = G[18]
                        inf2 = re.split('[=;]', inf)
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
 
with open(file44, 'r') as in1:
    with open(e_file, 'w') as eout:
        for line in in1:
            F = line.rstrip('\n').split('\t')
            ln = line.rstrip('\n')
            if F[0] in ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]:
                if "3" in F[6] and "+" in F[7]:
                    c=F[0]
                    s=int(F[2])-5
                    e=int(F[2])+5
                    rows = tbx.fetch(c, s, e)                   
                    if not list(rows):
                        eout.write(line)
                elif "3" in F[6] and "-" in F[7]:
                    c=F[0]
                    s=int(F[1])-5
                    e=int(F[1])+5
                    rows = tbx.fetch(c, s, e)
                    if not list(rows):
                        eout.write(line)
                elif "5" in F[6] and "-" in F[7]:
                    c=F[0]
                    s=int(F[2])-5
                    e=int(F[2])+5
                    rows = tbx.fetch(c, s, e)
                    if not list(rows):
                        eout.write(line)                    
                elif "5" in F[6] and "+" in F[7]:
                    c=F[0]
                    s=int(F[1])-5
                    e=int(F[1])+5
                    rows = tbx.fetch(c, s, e)
                    if not list(rows):
                        eout.write(line)                    
                