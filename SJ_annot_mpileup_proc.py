#! /usr/bin/env python

import re
import sys

#input_file = sys.argv[1]
input_file = "./alterativeSJ_add_mut/26merge_fil.annot.ass.adj.freq.AI-mut.cano.mpileup2_20190515.txt"
fo = "./alterativeSJ_add_mut/26merge_fil.annot.ass.adj.freq.AI-mut.cano.mpileup3_20190515.txt"

def remove_indel_bases(bases):
    
    # remaining = bases.replace('$', '')
    remaining = bases
    proc = ""
    while len(remaining) > 0:
        match = re.search(r'([\+\-])(\d+)', remaining)
        if match is None:
            proc = proc + remaining
            remaining = ""
        else:
            del_num = int(match.group(2))
            del_pos = match.start()

            proc = proc + remaining[0:del_pos]
            remaining = remaining[(del_pos + len(str(del_num)) + del_num + 1):len(remaining)]

    remaining = proc
    proc = ""
    while len(remaining) > 0:
        match = re.search(r'\^\S', remaining)
        if match is None:
            proc = proc + remaining
            remaining = ""
        else:
            pos = match.start()
            proc = proc + remaining[0:pos]
            remaining = remaining[(pos + 2):len(remaining)]

    proc = proc.replace('$', '')
    return proc



with open(input_file, 'r') as hin:
    with open(fo, 'w') as out1:
        for line in hin:
        
            F = line.rstrip('\n').split('\t')

            F37 = remove_indel_bases(F[37])
        
            depth = 0
            base2num = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0, 'a': 0, 'c': 0, 'g': 0, 't': 0, 'n': 0}
            base2pos = {'A': [], 'C': [], 'G': [], 'T': [], 'N': [], 'a': [], 'c': [], 'g': [], 't': [], 'n': []}
        
            for i in range(len(F37)):
                if F37[i] in ['>', '<', '*']: continue
                depth = depth + 1
                if F37[i] == '.':
                    base2num[F[35].upper()] = base2num[F[35].upper()] + 1
                    #base2pos[F[3].upper()].append(pos_vector[i])
                elif F37[i] == ',':
                    base2num[F[35].lower()] = base2num[F[35].lower()] + 1 
                    #base2pos[F[3].lower()].append(pos_vector[i])
                else:
                    base2num[F37[i]] = base2num[F37[i]] + 1
                    #base2pos[F5[i]].append(pos_vector[i])

            if depth == 0: continue
            
            alt = F[14]
            nuc = alt.upper()
            alt_rate = float(base2num[nuc] + base2num[nuc.lower()]) / depth
        
            depth_p = base2num['A'] + base2num['C'] + base2num['G'] + base2num['T']
            depth_n = base2num['a'] + base2num['c'] + base2num['g'] + base2num['t']
        
            var_p = base2num[alt.upper()]
            var_n = base2num[alt.lower()]
            #var_info = str(depth_p) + ',' + str(var_p) + ';' + str(depth_n) + ',' + str(var_n)
            strand_ratio =  str(base2num[nuc]) + ":" + str(base2num[nuc.lower()])

            rec = '\t'.join(F) + '\t' +  alt + '\t' + str(depth_p + depth_n) + '\t' + \
                str(var_p + var_n) + '\t' + str(round(alt_rate, 4)) + '\t' + str(strand_ratio)+'\n'
        
            out1.write(rec)
			
