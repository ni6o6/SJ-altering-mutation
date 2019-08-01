#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 13:01:49 2019

@author: genome
#def asj_ass_adj_bed(input_file, output_file, pr):
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

file2 = './alterativeSJ_fil_annot/%s/%s.SJ.fil.annot.txt' %(args.folder,pr)
file3='./alterativeSJ_assadjfreq/%s/%s.SJ.fil.annot.assadj.txt' %(args.folder,pr)
    
with open(file2, 'r') as in1:
        with open(file3, 'w') as out1:

            for line in in1:
                F = line.rstrip('\n').split('\t')
                r=F[6]
                if F[0] in ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]:
                    if "Alternative" in F[9] or "alternative" in F[9]:
                        if "e" in F[13]: #strand=+ 3'SS
                            a =str(F[14])
                            b = a.split(';',2)
                            if b[0] == '*':
                                b[0] = '0'
                                c=F[0]
                                si=int(F[1])-1*int(b[0])
                                ei=int(F[2])-1*int(b[0])
                                out1.write(str(c)+'\t'+str(si)+'\t'+str(ei)+str(F[1])+'\t'+str(F[2])+'\t'+str(pr)+'\t'+str(F[9])+'\t+\t'+str(r)+'\n')
                            else:
                                c=F[0]
                                si=int(F[1])-1*int(b[0])
                                ei=int(F[2])-1*int(b[0])
                                out1.write(str(c)+'\t'+str(si)+'\t'+str(ei)+'\t'+str(F[1])+'\t'+str(F[2])+'\t'+str(pr)+'\t'+str(F[9])+'\t+\t'+str(r)+'\n')
                        elif "s" in F[13]: #strand=- 5'SS
                            a =str(F[14])
                            b = a.split(';',2)
                            if b[0] == '*':
                                b[0] = '0'
                                c=F[0]
                                si=int(F[1])-1*int(b[0])
                                ei=int(F[2])-1*int(b[0])
                                out1.write(str(c)+'\t'+str(si)+'\t'+str(ei)+'\t'+str(F[1])+'\t'+str(F[2])+'\t'+str(pr)+'\t'+str(F[9])+'\t-\t'+str(r)+'\n')
                            else:
                                c=F[0]
                                si=int(F[1])-1*int(b[0])
                                ei=int(F[2])-1*int(b[0])
                                out1.write(str(c)+'\t'+str(si)+'\t'+str(ei)+'\t'+str(F[1])+'\t'+str(F[2])+'\t'+str(pr)+'\t'+str(F[9])+'\t-\t'+str(r)+'\n')
                        elif "s" in F[17]: #strand=+ 5'SS
                            a =str(F[18])
                            b = a.split(';',2)
                            if b[0] == '*':
                                b[0] = '0'
                                c=F[0]
                                si=int(F[1])-1*int(b[0])
                                ei=int(F[2])-1*int(b[0])
                                out1.write(str(c)+'\t'+str(si)+'\t'+str(ei)+'\t'+str(F[1])+'\t'+str(F[2])+'\t'+str(pr)+'\t'+str(F[9])+'\t-\t'+str(r)+'\n')
                            else:
                                c=F[0]
                                si=int(F[1])-1*int(b[0])
                                ei=int(F[2])-1*int(b[0])
                                out1.write(str(c)+'\t'+str(si)+'\t'+str(ei)+'\t'+str(F[1])+'\t'+str(F[2])+'\t'+str(pr)+'\t'+str(F[9])+'\t+\t'+str(r)+'\n')
                        elif "e" in F[17]: #strand=- 3'SS
                            a =str(F[18])
                            b = a.split(';',2)
                            if b[0] == '*':
                                b[0] = '0'
                                c=F[0]
                                si=int(F[1])-1*int(b[0])
                                ei=int(F[2])-1*int(b[0])
                                out1.write(str(c)+'\t'+str(si)+'\t'+str(ei)+'\t'+str(F[1])+'\t'+str(F[2])+'\t'+str(pr)+'\t'+str(F[9])+'\t-\t'+str(r)+'\n')
                            else:
                                c=F[0]
                                si=int(F[1])-1*int(b[0])
                                ei=int(F[2])-1*int(b[0])
                                out1.write(str(c)+'\t'+str(si)+'\t'+str(ei)+'\t'+str(F[1])+'\t'+str(F[2])+'\t'+str(pr)+'\t'+str(F[9])+'\t-\t'+str(r)+'\n')
