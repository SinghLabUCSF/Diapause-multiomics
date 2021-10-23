#!/usr/env python3

'''
Example Command Line Execution
python3 NeoF_DE.py
'''

###############################################################
#### Set Environmental Variables
###############################################################
CURRENT_DIR = './'

file1 = CURRENT_DIR + 'Data/Paralog_List.txt'
file2 = CURRENT_DIR + 'Data/Nfur_Up_Master.bed'
file3 = CURRENT_DIR + 'Data/Nfur_Up_NeoN1.bed'
file4 = CURRENT_DIR + 'Data/Nfur_Up_NeoN2.bed'

N1 = []
N2 = []

###############################################################
#### Read all Neo-functionalized paralogs into python list
###############################################################
with open(file1, 'r') as in1:
	for line in in1:
		line = line.rstrip('\n').split('\t')
		if line[2] == 'NeoF':
			N1.append(line[0])
			N2.append(line[1])

###############################################################
#### Sort peaks based on closest genes paralog status
###############################################################			
with open(file2, 'r') as in2, open(file3, 'w') as out1, open(file4, 'w') as out2:
	for line in in2:
		temp = line.rstrip('\n').split('\t')
		if temp[6] in N1:
			out1.write(line)
		if temp[6] in N2:
			out2.write(line)
			
