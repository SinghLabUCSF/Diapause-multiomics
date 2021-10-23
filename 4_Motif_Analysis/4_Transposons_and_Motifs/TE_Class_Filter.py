#!/usr/env python3

'''
Example Command Line Execution:
python3 TE_Class_Filter.py
'''

##########################################################
### Set environmental variables
##########################################################
CURRENT_DIR = './'

file1 = CURRENT_DIR + 'Data/Nfur_TE.bed'
file2 = CURRENT_DIR + 'Data/TE_Family_Data.txt'
file3 = CURRENT_DIR + 'Data/TE_Family_Data_Filtered.txt'

Checker = {}

##########################################################
### Read in and parse the names of TE families
##########################################################
with open(file1, 'r') as in1:
	for line in in1:
		line = line.rstrip('\n').split('\t')
		TE = line[4]
		if '/' in TE:
			TE = TE.split('/')
			Checker[TE[1]] = TE[0]
		else:
			Checker[TE] = TE

##########################################################
### Generate output using corrected TE family names
##########################################################		
with open(file2, 'r') as in2, open(file3, 'w') as out1:
	for line in in2:
		line = line.rstrip('\n').split('\t')
		check = line[0]
		outer = '\t'.join(line)
		if check == 'TE' or check == 'Total':
			out1.write(outer + '\tClass\n')
		else:
			out1.write(outer + '\t' + Checker[check] + '\n')
	