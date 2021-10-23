#!/usr/env python3

'''
Example Command Line Entry
python3 Motif_Correction.py
'''

#######################################################################################
##### Define Constants
#######################################################################################
import sys
import os

CURRENT_DIR = './'

home = CURRENT_DIR
dictor = CURRENT_DIR + 'Data/Master_DE_Dia_Up_Peak.bed'
DES = []

#######################################################################################
##### Load all DE peak numbers into list
#######################################################################################
with open(dictor, 'r') as in1:
	for line in in1:
		line = line.strip('\n').split('\t')
		peak = line[3]
		DES.append(peak)

#######################################################################################
##### Correct DE status of found in each Motif Count File
#######################################################################################
for file in os.listdir(home + 'Data/Raw/'):
	namer = file[:-4]
	with open(home + 'Data/Raw/' + file, 'r') as temp_in, open(home + 'Data/DEC/' + file, 'w') as temp_out:
		for line in temp_in:
			line = line.rstrip('\n').split('\t')
			start = line[0]
			check = line[1]
			fix1 = line[3]
			fix2 = line[2]
			end = line[4:]
			
			if check in DES:
				fix1 = 'Diapause'
				fix2 = 'DE'
			else:
				fix1 = 'Other'
				fix2 = 'OT'
				
			rejon = '\t'.join(end)
			writ = start + '\t' + check + '\t' + fix2 + '\t' + fix1 + '\t' + rejon + '\n'
			temp_out.write(writ)
			