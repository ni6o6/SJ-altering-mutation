#!/usr/bin/env python3

"""
Created on Wed Jul 31 2019

@author: naokoIida
"""
import sys
import os

fo = sys.argv[1]

if not os.path.exists('alterativeSJ_fil_annot'):
    os.mkdir('alterativeSJ_fil_annot')
if not os.path.exists('alterativeSJ_fil_annot/'+fo):
    os.mkdir('alterativeSJ_fil_annot/'+fo)
if not os.path.exists('alterativeSJ_assadjfreq'):
    os.mkdir('alterativeSJ_assadjfreq')
if not os.path.exists('alterativeSJ_assadjfreq/'+fo):
    os.mkdir('alterativeSJ_assadjfreq/'+fo)
if not os.path.exists('alterativeSJ_tabixAI'):
    os.mkdir('alterativeSJ_tabixAI')
if not os.path.exists('alterativeSJ_tabixAI/'+fo):
    os.mkdir('alterativeSJ_tabixAI/'+fo)
if not os.path.exists('alterativeSJ_cmut'):
    os.mkdir('alterativeSJ_cmut')
if not os.path.exists('alterativeSJ_cmut/'+fo):
    os.mkdir('alterativeSJ_cmut/'+fo)
    
