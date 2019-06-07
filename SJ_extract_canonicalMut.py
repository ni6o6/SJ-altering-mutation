#! /usr/bin/env python
"""
Naoko Iida

"""
import numpy as np
import pandas as pd
import glob
import os
import subprocess

outfile="./alterativeSJ_tabixAI/26merge_fil.annot.ass.adj.freq.AI.cano.txt"

subprocess.call(["rm", "-f", outfile])

#he = 'chr\tstart\tend\tstart_ori\tend_ori\tsample\tclass\tstrand\treads\ttotal\tfreq\tCHR\tPOS\tID\tREF_x\tMUT_x\tSCORE1\tSCORE2\tSYMBOL\tSTRAND\tTYPE\tDIST\tDS_AG\tDS_AL\tDS_DG\tDS_DL\tDP_AG\tDP_AL\tDP_DG\tDP_DL\tSAVnet_SpliceAI\tREF_y\tMUT_y\tcanonicalMut\tDP\tDS\n'

infile=glob.glob("./alterativeSJ_tabixAI/*.SJ.fil.annot.ass.adj.freq.AI.txt") #cd
for f in infile:
    basename = os.path.basename(f)
    pr=basename.split('.', 8)[0]
    print(pr)

    with open(f, 'r') as in1:
            
        with open(outfile, 'a') as out1:
            for line in in1:
                F = line.rstrip('\n').split('\t')
                ln = line.rstrip('\n')
                if int(F[8])>=3 and float(F[10])>=0.05 and "c" in F[30]:
                #if int(F[9])>=60 and float(F[10])>=0.05 and "c" in F[30]:
        #3'SS strand=+
                    if "3" in F[6] and "+" in F[7]:
                        a =str(F[22])
                        DP =str(F[26])
                        if F[26] is "1" and F[15] is "G": #mt
                            m= "*"+F[14]+"|_"+F[15]
                            out1.write(ln+'\t'+m+'\t'+str(DP)+'\t'+str(a)+'\n')
                        elif F[26] is "2" and F[15] is "A":
                            m= F[14]+"*|_"+F[15]
                            out1.write(ln+'\t'+m+'\t'+str(DP)+'\t'+str(a)+'\n')
        #3'SS strand=- 
                    elif "3" in F[6] and "-" in F[7]:
                        if F[14] is "G":
                            rb="C"
                        elif F[14] is "C":
                            rb="G"
                        elif F[14] is "T":
                            rb="A"
                        elif F[14] is "A":
                            rb="T"
                        if F[15] is "G":
                            mb="C"
                        elif F[15] is "C":
                            mb="G"
                        elif F[15] is "T":
                            mb="A"
                        elif F[15] is "A":
                            mb="T"
                        a =str(F[22])
                        DP =str(F[26])
                        if F[26] is "-1" and F[15] is "C":
                            m= "*"+rb+"|_"+mb
                            out1.write(ln+'\t'+m+'\t'+str(DP)+'\t'+str(a)+'\n')
                        elif F[26] is "-2" and F[15] is "T":
                            m= rb+"*|_"+mb
                            out1.write(ln+'\t'+m+'\t'+str(DP)+'\t'+str(a)+'\n')
        #5'SS strand=+
                    elif "5" in F[5] and "+" in F[7]:
                        a =str(F[24])
                        DP =str(F[28])
                        if F[28] is "-1" and F[15] is "G":
                            m= "|"+F[14]+"*_"+F[15]
                            out1.write(ln+'\t'+m+'\t'+str(DP)+'\t'+str(a)+'\n')
                        elif F[28] is "-2" and F[15] is "T":
                            m= "|*"+F[14]+"_"+F[15]
                            out1.write(ln+'\t'+m+'\t'+str(DP)+'\t'+str(a)+'\n')
        #5'SS strand=-
                    elif "5" in F[6] and "-" in F[7]:
                        a =str(F[24])
                        DP =str(F[28])
                        if F[14] is "G":
                            rb="C"
                        elif F[14] is "C":
                            rb="G"
                        elif F[14] is "T":
                            rb="A"
                        elif F[14] is "A":
                            rb="T"
                        if F[15] is "G":
                            mb="C"
                        elif F[15] is "C":
                            mb="G"
                        elif F[15] is "T":
                            mb="A"
                        elif F[15] is "A":
                            mb="T"
                        if F[28] is "1" and F[15] is "C":
                            m= "|"+rb+"*_"+mb
                            out1.write(ln+'\t'+m+'\t'+str(DP)+'\t'+str(a)+'\n')
                        elif F[28] is "2" and F[15] is "A":
                            m= "|*"+rb+"_"+mb
                            out1.write(ln+'\t'+m+'\t'+str(DP)+'\t'+str(a)+'\n')


