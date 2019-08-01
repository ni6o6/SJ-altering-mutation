#!/usr/bin/env python3

"""
Created on Wed Jul 31 2019

@author: naokoIida
"""
import argparse
import subprocess
import shutil
from junc_utils import utils


parser = argparse.ArgumentParser() #make a parser

parser.add_argument("input", metavar = "./junction/sample.SJ.out.tab", default = None, type = str,
                        help = "Path to input file") 
    
parser.add_argument("folder", metavar = "group", default = "my_samples", type = str,
                        help = "Path to input file") 
    
parser.add_argument('--control_file', nargs='*', type = str,
                        help = "Path to control data created by merge_control (default: %(default)s), reads filter>=1")
    
parser.add_argument("--genome_id", choices = ["hg19", "hg38", "mm10"], default = "hg19",
                          help = "Genome id used for selecting UCSC-GRC chromosome name corresponding files (default: %(default)s)")
args = parser.parse_args()

pr = args.input   
print(pr)    

file1 = './alterativeSJ_fil_annot/%s/%s.SJ.fil.txt' %(args.folder,pr)
file2 = './alterativeSJ_fil_annot/%s/%s.SJ.fil.annot.txt' %(args.folder,pr)
    
if not args.control_file:
    shutil.copy(args.input_file, file1)

else:
    f = args.input_file
    cont_list=args.control_file
    n=1
    for cont in cont_list:
        out = 'tmp_out'+str(n)+'.txt'
        #junc_utils filter --pooled_control_file *.bed.gz input(./junction/A427.SJ.out.tab) output
        utils.proc_star_junction(f, out, cont, 1, 10, True, False) #<--reads>=1
        f = 'tmp_in'+str(n)+'.txt'
        shutil.copy(out, f)
        n=n+1 
    shutil.copy(out, file1)
    
    annotate_commands = ["junc_utils", "annotate", file1, file2, "--genome_id", args.genome_id]
    subprocess.call(annotate_commands)
    