#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 15:30:10 2019

@author: genome
"""

import pysam
import argparse

parser = argparse.ArgumentParser() #make a parser

parser.add_argument("input", metavar = "./junction/sample.SJ.out.tab", default = None, type = str,
                        help = "Path to input file") 
    
parser.add_argument("folder", metavar = "group", default = "my_samples", type = str,
                        help = "Path to input file") 
    
args = parser.parse_args()

pr = args.input
print(pr) 

db = "./database/gnomad.genomes.r2.1.1.sites.vcf.bgz"
#db = "/Volumes/NIIDA_SSD1R/gnomad/gnomad.exomes.r2.1.1.sites.vcf.bgz"
tb = pysam.TabixFile(db)

file = './alterativeSJ_cmut/%s/%s.SJ.fil.annot.assadjunifreqT.AI.cmut.txt' %(args.folder,pr)   

he = 'chr\tstart\tend\tstart_ori\tend_ori\tsample\tclass\tstrand\treads\ttotal\tfreq\tCHR\tPOS\tID\tREF\tMUT\tSCORE1\tSCORE2\tSYMBOL\tSTRAND\tTYPE\tDIST\tDS_AG\tDS_AL\tDS_DG\tDS_DL\tDP_AG\tDP_AL\tDP_DG\tDP_DL\tSAVnet_SpliceAI\tMOTIF\tDP\tDS\n'
out_file = './alterativeSJ_cmut/%s/%s.SJ.fil.annot.assadjunifreqT.AI.cmut.snp.txt' %(args.folder,pr)
    
with open(file, 'r') as hin:
    with open(out_file,'w') as hout:
        for line in hin:
            line = line.rstrip('\n')
            junc_record = line
            F = line.split('\t')
            allele = "na"
            snp = "na" 
            gnomad_snp = []
            if F[11] in ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"]:
    #c = tb.fetch(F[11], int(F[12]) - 1, int(F[12]))
                for record_line in tb.fetch(F[11], int(F[12]) - 1, int(F[12])):
                    record = record_line.split('\t')
                    #print('ref',F[14], F[15])
                    #print('gnomad', record[3],record[4])
                        
                    if F[14] == record[3] and F[15] == record[4]:
                        allele = record[3]+">"+record[4]
                        infos = record[7].split(';')
                        for info in infos:
                            if info.startswith("AF="):
                                freq = float(info.replace("AF=", ''))
                        gnomad_snp.append(str(freq))
                            
                snp = ','.join(gnomad_snp)
                out_record = junc_record + "\t" + allele + "\t" + snp 
                print(out_record, file = hout)
                                       
            elif F[11] in ["Y"]:
                out_record = junc_record + "\tna\t" 
                print(out_record, file = hout)
            else:
                out_record = junc_record + "\tsnp\tfreq" 
                print(out_record, file = hout)

