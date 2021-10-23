#!/usr/env python3

'''
Example Command Line Execution
python3 Peak_Cons_Bed_Generator_0.py
'''

###############################################################
#### Set Environmental Variables
###############################################################

CURRENT_DIR = './'
ALLER = ['singletons_Comb.txt', 'lowCV_GENES.txt', 'Singletons_Dev.txt', 'Singletons_Dia.txt']

for THINGER in ALLER:
	SPECIAL = THINGER.split('.')
	SPECIAL = SPECIAL[0] + '_DE.' + SPECIAL[1]
	file1 = CURRENT_DIR + 'Data/Nfur_annotations.txt'
	file2 = CURRENT_DIR + 'Data/lowCV_PEAKS.bed'
	file3 = CURRENT_DIR + 'Data/lowCV_PEAKS_final.bed'
	file4 = CURRENT_DIR + 'Data/' + THINGER
	file5 = CURRENT_DIR + 'Data/' + THINGER[:-4] + '.bed'
	file6 = CURRENT_DIR + 'Data/' + SPECIAL[:-4] + '.bed'
	file7 = CURRENT_DIR + 'Data/Nfur_Peak_All.bed'

	hold1 = {}
	hold2 = {}
	hold3 = []

	###############################################################
	#### Read in Nfur Annotation Data
	###############################################################
	with open(file1, 'r') as in1:
		for line in in1:
			line = line.rstrip('\n').split('\t')
			peak = line[0]
			gene = line[25]
			exer = line[3] + '\t' + line[4] + '\t' + line[5] + '\t' + line[0] + '\t' + line[25] + '\n'
			hold1[peak] = gene
		
			if gene in hold2.keys():
				hold2[gene].append(exer)
		
			else:
				hold2[gene] = []
				hold2[gene].append(exer)						

	###############################################################
	#### Read and Clean Test set bed file
	###############################################################			
	with open(file2, 'r') as in2, open(file3, 'w') as out1:
		for line in in2:
			line = line.rstrip('\n')
			temp = line.split('\t')
			tempU = temp[3]
			if tempU in hold1.keys():
				fix = hold1[tempU]
				outer = line + '\t' + fix + '\n'
			else:
				fix = 'NA'
				outer = line + '\t' + fix + '\n'
			out1.write(outer)

	###############################################################
	#### Read in Gene coordinate file
	###############################################################		
	with open(file7, 'r') as in21:
		for line in in21:
			line = line.rstrip('\n').split('\t')
			pick = line[3]
			hold3.append(pick)

	###############################################################
	#### Read in gene list of interest to combine peaks and genes
	###############################################################		
	with open(file4, 'r') as in3, open(file5, 'w') as out2, open(file6, 'w') as out3:
		for line in in3:
			TempW = line.rstrip('\n')
			if TempW in hold2.keys():
				for item in hold2[TempW]:
					Tempx = item.split('\t')
					Tempy = Tempx[3]
					if Tempy in hold3:
						out2.write(item)
						out3.write(item)
					else:
						out2.write(item)

		