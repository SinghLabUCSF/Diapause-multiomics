#!/usr/env python3

'''
Example Command Line Execution
python3 TE_Motif.py
'''

##########################################################
### Set environmental variables
##########################################################
import subprocess

CURRENT_DIR = './'
PACK = ''

file01 = CURRENT_DIR + 'Data/Interest_Motif_sort.bed'
file02 = CURRENT_DIR + 'Data//Nfur_TE.bed'
file1 = CURRENT_DIR + 'Data/TE_Motif_Overlap.bed'

##########################################################
### Overlap TEs and Motifs
##########################################################
with open(file1, 'wt') as temp:
	subprocess.call([PACK + 'bedtools', 'intersect', '-wa', '-wb', '-a', file01,  '-b', file02], stdout=temp)
   