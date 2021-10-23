#!/usr/env python3

'''
Example Command Line Input
python3 PC_Loading_Bed.py
'''

##########################################################
### Set environmental variables
##########################################################
CURRENT_DIR = './'

file0 = CURRENT_DIR + 'Data/Peak_Master_0.0.txt'
file1 = CURRENT_DIR + 'Data/Alim_annotations.txt'
file2 = CURRENT_DIR + 'Data/Nfur_annotations.txt'
file3 = CURRENT_DIR + 'Data/Nfur_Alim_PC2.txt'
file4 = CURRENT_DIR + 'Data/Nfur_PC2.bed'
file5 = CURRENT_DIR + 'Data/Alim_PC2.bed'

##########################################################
### Create empty data structures for input
##########################################################
nf_2_al = {}
pk_2_if = {}
alpk_2_alif = {}

##########################################################
### Read in Peak overlap data
##########################################################
with open(file0, 'r') as in0:
	#Parse peak conservation data
    for line in in0:
        line = line.rstrip('\n').split('\t')
        if line[0] != 'Nfur':
            Nfur = line[0]
            Alim = line[3]
            #Remove name correction from peaks
            Alim = Alim.split('X')
            Alim = Alim[0]
            #Store data in dictionary
            nf_2_al[Nfur] = Alim

##########################################################
### Read in South American Killifish peak metadata
##########################################################       
with open(file1, 'r') as in1:
	#Parse bed file data
    for line in in1:
        line = line.rstrip('\n').split('\t')
        #Store data in dictionary
        if line[0] != 'names':
            alpk_2_alif[line[0]] = line[3] + '\t' + line[4] + '\t' + line[5] + '\t' + line[0] + '\t' + line[11] + '\t' + line[17]            

##########################################################
### Read in African Turquoise Killifish peak metadata
##########################################################            
with open(file2, 'r') as in2:
	#Parse bed file data
    for line in in2:
        line = line.rstrip('\n').split('\t')
        #Store data in dictionary
        if line[0] != 'names':
            pk_2_if[line[0]] = line[3] + '\t' + line[4] + '\t' + line[5] + '\t' + line[0] + '\t' + line[19] + '\t' + line[25]
            

##########################################################
### Combine metadata from peaks found in PC2 loadings
##########################################################       
with open(file3, 'r') as in3, open(file4, 'w') as out1, open(file5, 'w') as out2:
	#Parse PC2 loading list
    for line in in3:
        line = line.rstrip('\n').split('\t')
        if line[0] != 'peak':
        	
        	#Check for peak metadata in dictionaries
            peaker = 'peak_' + str(line[0])
            info = pk_2_if[peaker]
            peaker2 = nf_2_al[peaker]
            if peaker2 in alpk_2_alif.keys():
                info2 = alpk_2_alif[peaker2]
                load = line[1]
        		
        		#Output metadata as bed file for peaks
                out1.write(info + '\t' + load + '\n') 
                out2.write(info2 + '\t' + load + '\n') 
                            