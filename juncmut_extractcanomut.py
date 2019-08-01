#!/usr/bin/env python3

"""
Created on Wed Jul 31 13:56:22 2019

@author: genome
#def asj_cmut #extraction of_canonical mutation
"""

import argparse

parser = argparse.ArgumentParser() #make a parser

parser.add_argument("input", metavar = "./junction/sample.SJ.out.tab", default = None, type = str,
                        help = "Path to input file") 
    
parser.add_argument("folder", metavar = "group", default = "my_samples", type = str,
                        help = "Path to input file") 
    
args = parser.parse_args()

pr = args.input
print(pr) 

    
def complement_dna(string):
    comp=''
    for char in string:
        if   char == 'A': comp += 'T'
        elif char == 'T': comp += 'A'
        elif char == 'G': comp += 'C'
        elif char == 'C': comp += 'G'
        else:             comp += char
    return comp[::-1]

file5 = './alterativeSJ_tabixAI/%s/%s.SJ.fil.annot.assadjunifreqT.AI.txt' %(args.folder,pr)    
file6 = './alterativeSJ_cmut/%s/%s.SJ.fil.annot.assadjunifreqT.AI.cmut.txt' %(args.folder,pr)
he = 'chr\tstart\tend\tstart_ori\tend_ori\tsample\tclass\tstrand\treads\ttotal\tfreq\tCHR\tPOS\tID\tREF\tMUT\tSCORE1\tSCORE2\tSYMBOL\tSTRAND\tTYPE\tDIST\tDS_AG\tDS_AL\tDS_DG\tDS_DL\tDP_AG\tDP_AL\tDP_DG\tDP_DL\tSAVnet_SpliceAI\tMOTIF\tDP\tDS\n'
with open(file6, 'w') as out1:
    out1.write(he)
    
    with open(file5, 'r') as in1:   
        
        for line in in1:
            F = line.rstrip('\n').split('\t')
            ln = line.rstrip('\n')
            if "c" in F[30]:
        #3'SS strand=+
                if "3" in F[6] and "+" in F[7]:
                    a =str(F[22])
                    DP =str(F[26])
                    if F[26] == "1" and F[15] == "G":
                        m= "*"+F[14]+"|_"+F[15]
                        out1.write(ln+'\t'+m+'\t'+str(DP)+'\t'+str(a)+'\n')
                    elif F[26] == "2" and F[15] == "A":
                        m= F[14]+"*|_"+F[15]
                        out1.write(ln+'\t'+m+'\t'+str(DP)+'\t'+str(a)+'\n')
        #3'SS strand=- 
                elif "3" in F[6] and "-" in F[7]:
                    rb=complement_dna(F[14])
                    mb=complement_dna(F[15])
                    a =str(F[22])
                    DP =str(F[26])
                    if F[26] == "-1" and F[15] == "C":
                        m= "*"+rb+"|_"+mb
                        out1.write(ln+'\t'+m+'\t'+str(DP)+'\t'+str(a)+'\n')
                    elif F[26] == "-2" and F[15] == "T":
                        m= rb+"*|_"+mb
                        out1.write(ln+'\t'+m+'\t'+str(DP)+'\t'+str(a)+'\n')
        #5'SS strand=+
                elif "5" in F[6] and "+" in F[7]:
                    a =str(F[24])
                    DP =str(F[28])
                    if F[28] == "-1" and F[15] == "G":
                        m= "|"+F[14]+"*_"+F[15]
                        out1.write(ln+'\t'+m+'\t'+str(DP)+'\t'+str(a)+'\n')
                    elif F[28] == "-2" and F[15] == "T":
                        m= "|*"+F[14]+"_"+F[15]
                        out1.write(ln+'\t'+m+'\t'+str(DP)+'\t'+str(a)+'\n')
        #5'SS strand=-
                elif "5" in F[6] and "-" in F[7]:
                    rb=complement_dna(F[14])
                    mb=complement_dna(F[15])
                    a =str(F[24])
                    DP =str(F[28])
                    if F[28] == "1" and F[15] == "C":
                        m= "|"+rb+"*_"+mb
                        out1.write(ln+'\t'+m+'\t'+str(DP)+'\t'+str(a)+'\n')
                    elif F[28] == "2" and F[15] == "A":
                        m= "|*"+rb+"_"+mb
                        out1.write(ln+'\t'+m+'\t'+str(DP)+'\t'+str(a)+'\n')
                            