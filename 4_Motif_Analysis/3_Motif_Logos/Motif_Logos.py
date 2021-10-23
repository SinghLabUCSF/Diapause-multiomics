#!/usr/env python3

'''
Example Command Line Input
python3 Motif_Logos.py
'''

##########################################################
### Set environmental variables
##########################################################
import sys
import os
CURRENT_DIR = './'

##########################################################
### Generate Ideal Motif LOGO sets
##########################################################
#search data files for all inputs
for file in os.listdir(CURRENT_DIR + 'Data/Matrix/'):
	#verify input status
	if '_matrix.txt' in file:
		#name input and output
		header = file[:-11]
		Input = CURRENT_DIR + 'Data/Matix/' + file
		Output = CURRENT_DIR + 'Data/Manuals/' + header + '_Manual.txt'
		lineup = {}
		Converse = {1:'A', 2:'C', 3:'G', 4:'T'}
		base = 0
		
		#Read in input file
		with open(Input, 'r') as in1:
			for line in in1:
				base =  base + 1
				line = line.rstrip('\n').split('\t')
				#Store data in list-embedded dictionary
				lineup[base] = []
				freq = 0
				
				#Convert base percentages to counts out of 1000
				for entry in line:
					freq = freq + 1
					putter = float(entry)*1000
					step = 0
			
					while step < putter:
						step =  step + 1
						lineup[base].append(Converse[freq])
		
		#Open output file
		with open(Output, 'w') as out1:
			shift = 0			
			while shift < 999:
				namer = ''
				#Traverse data count structure to generate output
				for item in lineup.keys():
					namer = namer + lineup[item][shift]
				out1.write(namer + '\n')
				shift = shift + 1
	