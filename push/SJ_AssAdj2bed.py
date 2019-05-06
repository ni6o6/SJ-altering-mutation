#! /usr/bin/env python

"""
Naoko Iida
#source activate py27
#python SJ_20190501_AssAdj2bed.py
#input *.SJ.filtered.annot.txt
#Filter with Alterative SS (ass.
#Ajustment of intron start end (adj.
#next, intersectBed 
"""

import glob
import os
import pandas as pd
import numpy as np

in1=glob.glob("./alterativeSJ_filtered_annot/*.SJ.filtered.annot.txt")

for infile in in1:
    basename = os.path.basename(infile)
    pr=basename.split('.', 5)[0]
    print(pr)
    outfile='./alterativeSJ_adjustment/%s.SJ.filtered.annot.ass.adj.bed' %(pr)

    with open(infile, 'r') as in1:
        with open(outfile, 'w') as out1:
            #he = 'chr\tstart\tSJ_3\tSJ_4\tSJ_5\tSJ_6\tSJ_7\tSJ_8\tSJ_9\tSplicing_Class\tIs_Inframe\tGene_1\tExon_Num_1\tIs_Boundary_1\tOffset_1\tGene_2\tExon_Num_2\tIs_Boundary_2\tOffset_2\tchr\tintron_start_adj\tintron_end_adj\n'
            #out1.write(he)

            for line in in1:
                F = line.rstrip('\n').split('\t')
                ln = line.rstrip('\n')
                r=F[6]
                if F[0] in ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]:
                    if "Alternative" in F[9] or "alternative" in F[9]:
                        if "e" in F[13]: #strand=+ 3'SS
                            a =str(F[14])
                            b = a.split(';',2)
                            if b[0] is '*':
                                b[0] = '0'
                                c=F[0]
                                si=int(F[1])-1*int(b[0])
                                ei=int(F[2])-1*int(b[0])
                                out1.write(str(c)+'\t'+str(si)+'\t'+str(ei)+'\t'+str(pr)+'\t'+str(F[9])+'\t+\t'+str(r)+'\n')
                            else:
                                c=F[0]
                                si=int(F[1])-1*int(b[0])
                                ei=int(F[2])-1*int(b[0])
                                out1.write(str(c)+'\t'+str(si)+'\t'+str(ei)+'\t'+str(pr)+'\t'+str(F[9])+'\t+\t'+str(r)+'\n')
                        elif "s" in F[13]: #strand=- 5'SS
                            a =str(F[14])
                            b = a.split(';',2)
                            if b[0] is '*':
                                b[0] = '0'
                                c=F[0]
                                si=int(F[1])-1*int(b[0])
                                ei=int(F[2])-1*int(b[0])
                                out1.write(str(c)+'\t'+str(si)+'\t'+str(ei)+'\t'+str(pr)+'\t'+str(F[9])+'\t-\t'+str(r)+'\n')
                            else:
                                c=F[0]
                                si=int(F[1])-1*int(b[0])
                                ei=int(F[2])-1*int(b[0])
                                out1.write(str(c)+'\t'+str(si)+'\t'+str(ei)+'\t'+str(pr)+'\t'+str(F[9])+'\t-\t'+str(r)+'\n')
                        elif "s" in F[17]: #strand=+ 5'SS
                            a =str(F[18])
                            b = a.split(';',2)
                            if b[0] is '*':
                                b[0] = '0'
                                c=F[0]
                                si=int(F[1])-1*int(b[0])
                                ei=int(F[2])-1*int(b[0])
                                out1.write(str(c)+'\t'+str(si)+'\t'+str(ei)+'\t'+str(pr)+'\t'+str(F[9])+'\t-\t'+str(r)+'\n')
                            else:
                                c=F[0]
                                si=int(F[1])-1*int(b[0])
                                ei=int(F[2])-1*int(b[0])
                                out1.write(str(c)+'\t'+str(si)+'\t'+str(ei)+'\t'+str(pr)+'\t'+str(F[9])+'\t+\t'+str(r)+'\n')
                        elif "e" in F[17]: #strand=- 3'SS
                            a =str(F[18])
                            b = a.split(';',2)
                            if b[0] is '*':
                                b[0] = '0'
                                c=F[0]
                                si=int(F[1])-1*int(b[0])
                                ei=int(F[2])-1*int(b[0])
                                out1.write(str(c)+'\t'+str(si)+'\t'+str(ei)+'\t'+str(pr)+'\t'+str(F[9])+'\t-\t'+str(r)+'\n')
                            else:
                                c=F[0]
                                si=int(F[1])-1*int(b[0])
                                ei=int(F[2])-1*int(b[0])
                                out1.write(str(c)+'\t'+str(si)+'\t'+str(ei)+'\t'+str(pr)+'\t'+str(F[9])+'\t-\t'+str(r)+'\n')


