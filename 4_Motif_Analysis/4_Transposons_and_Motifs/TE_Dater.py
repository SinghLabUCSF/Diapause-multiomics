#!/usr/env python3

'''
Example Command Line Execution:
python3 TE_Dater.py
'''

##########################################################
### Set environmental variables
##########################################################
CURRENT_DIR = './'

file1 = CURRENT_DIR + 'Data/Nfur_overlap.bed'
file2 = CURRENT_DIR + 'Data/Aaus_overlap.bed'
file3 = CURRENT_DIR + 'Data/Alim_overlap.bed'
file4 = CURRENT_DIR + 'Data/Olat_overlap.bed'
file5 = CURRENT_DIR + 'Data/Drer_overlap.bed'
file6 = CURRENT_DIR + 'Data/Peak_Master_0.0.txt'
file7 = CURRENT_DIR + 'Data/TE_Conservation.txt'

Nfur = {}
Aaus = {}
Alim = {}
Olat = {}
Drer = {}


##########################################################
### Define Functions
##########################################################

'''
This function parses a TE/Peak overlapped bed file for a given
species and saves the data within a python dictionary
'''

def TE_load(a_file, a_dict):
	with open(a_file, 'r') as in1:
		for line in in1:
			line = line.rstrip('\n').split('\t')
			TE = line[4]
			Peak = line[8]
			if Peak in a_dict.keys():
				if TE not in a_dict[Peak]:
					a_dict[Peak] = a_dict[Peak] + ',' + TE
			else:
				a_dict[Peak] = TE


##########################################################
### Read in TEs from each species
##########################################################				
TE_load(file1,Nfur)
TE_load(file2,Aaus)
TE_load(file3,Alim)
TE_load(file4,Olat)
TE_load(file5,Drer)

##########################################################
### Date TE sites by species conservation inference
##########################################################
with open(file6, 'r') as in2, open(file7, 'w') as out1:
	for line in in2:
		line = line.rstrip('\n').split('\t')
		if line[0] != 'Nfur':
			hold = line[0]
			
			if line[0] in Nfur.keys():
				hold = hold + '\t' + Nfur[line[0]]
			else:
				hold = hold + '\tNA'
			
			aak = line[1].split('X')
			aak = aak[0]	
			if aak in Aaus.keys():
				hold = hold + '\t' + Aaus[aak]
			else:
				hold = hold + '\tNA'
			
			alk = line[3].split('X')
			alk = alk[0]   
			if alk in Alim.keys():
				hold = hold + '\t' + Alim[alk]
			else:
				hold = hold + '\tNA'
			
			olk = line[4].split('X')
			olk = olk[0]	
			if olk in Olat.keys():
				hold = hold + '\t' + Olat[olk]
			else:
				hold = hold + '\tNA'
			
			drk = line[5].split('X')
			drk = drk[0]	
			if drk in Drer.keys():
				hold = hold + '\t' + Drer[drk]
			else:
				hold = hold + '\tNA'
				
			out1.write(hold + '\n')
				
  
			
			
			  
			