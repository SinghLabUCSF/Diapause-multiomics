#!/usr/env python3

'''
Example Command Line Input
python3 PS_Merger.py
'''

##########################################################
### Set environmental variables
##########################################################
CURRENT_DIR = './'

file1 = CURRENT_DIR + 'Data/Chunk_All_PS_Final.txt'
file2 = CURRENT_DIR + 'Data/Chunk_All_KS_Final.txt'
file3 = CURRENT_DIR + 'Data/Chunk_All_PSKS_Merge.txt'

inner1 = {}
inner2 = {}
fate = {}

##########################################################
### Read in Positive selection files and save data
##########################################################
with open(file1, 'r') as in1, open(file2, 'r') as in2: 
	for line in in1:
		temp = line.rstrip('\n')
		temp2 = temp.split('\t')
		peak = temp2[3]
		if peak != 'ChunkName':
			fdr = temp2[18]
			Group1 = '\t'.join(temp2[0:5])
			Group2 = '\t'.join(temp2[6:12])
			GroupT = Group1 +'\t'+ Group2 +'\t'+ temp2[17]
		
		
			if fdr != 'fdr':
				if float(fdr) <= 0.1:
					fate[peak] = 1
				else:
					fate[peak] = 0
			else:
				fate[peak] = fdr
			inner1[peak] = [GroupT, temp2[5], temp2[12], temp2[13], temp2[14], temp2[15], temp2[16], temp2[18]]
		
	for line in in2:
		temp = line.rstrip('\n')
		temp2 = temp.split('\t')
		peak = temp2[3]
		if peak != 'ChunkName':
			fdr = temp2[18]
			Group1 = '\t'.join(temp2[0:5])
			Group2 = '\t'.join(temp2[6:12])
			GroupT = Group1 +'\t'+ Group2 +'\t'+ temp2[17]
		
			if fdr != 'fdr':
				if float(fdr) <= 0.1 and peak in fate.keys() and fate[peak] == 1:
					fate[peak] = 'Both_Sig'
				elif float(fdr) <= 0.1:
					fate[peak] = 'KS_Sig'
				elif peak in fate.keys() and fate[peak] == 1:
					fate[peak] = 'PS_Sig'
				else:
					fate[peak] = 'No_Sig' 
			else:
				fate[peak] = fdr	 
			inner2[peak] = [GroupT, temp2[5], temp2[12], temp2[13], temp2[14], temp2[15], temp2[16], temp2[18]]

##########################################################
### Check for positive selection state of a give peak
##########################################################		
	for thing in fate.keys():
		if fate[thing] == 0:
			fate[thing] = 'No_Sig'
		elif fate[thing] == 1:
			fate[thing] = 'PS_Sig'

##########################################################
### Output merged set of peaks
##########################################################
with open(file3, 'w') as out1:
	out1.write('Chromosome\tStart\tEnd\tChunkName\tlength\tConservation\ttype\tDEStat\tgene\tparalog\tMotifs\tPeakName\tPS_PeakPercent\tKS_PeakPercent\tPS_Score\tKS_Score\tPS_Mismatches\tKS_Mismatches\tPS_PercentMM\tKS_PercentMM\tPS_PercentG\tKS_PercentG\tPS_Pvalue\tKS_Pvalue\tPS_fdr\tKS_fdr\tSUM_fdr\tUnion\n')
	for entry in inner1.keys():
		if fate[entry] != 'No_Sig':
			Final = 'Union'
		else:
			Final = 'Exclude'
		if entry in inner2.keys():
			out1.write(inner1[entry][0] + '\t' + inner1[entry][1] + '\t' + inner2[entry][1] + '\t' + inner1[entry][2] + '\t' + inner2[entry][2] + '\t' + inner1[entry][3] + '\t' + inner2[entry][3] + '\t' + inner1[entry][4] + '\t' + inner2[entry][4] + '\t' + inner1[entry][5] + '\t' + inner2[entry][5] + '\t' + inner1[entry][6] + '\t' + inner2[entry][6] + '\t' + inner1[entry][7] + '\t' + inner2[entry][7] + '\t' +fate[entry] + '\t' + Final + '\n')
		else:
			out1.write(inner1[entry][0] + '\t' + inner1[entry][1] + '\tNA\t' + inner1[entry][2] + '\tNA\t' + inner1[entry][3] + '\tNA\t' + inner1[entry][4] + '\tNA\t' + inner1[entry][5] + '\tNA\t' + inner1[entry][6] + '\tNA\t' + inner1[entry][7] + '\tNA\t' +fate[entry] + '\t' + Final + '\n')
			
	for entry in inner2.keys():
		if entry not in inner1.keys():
			out1.write(inner2[entry][0] + '\tNA\t' + inner2[entry][1] + '\tNA\t' + inner2[entry][2] + '\tNA\t' + inner2[entry][3] + '\tNA\t' + inner2[entry][4] + '\tNA\t' + inner2[entry][5] + '\tNA\t' + inner2[entry][6] + '\tNA\t' + inner2[entry][7] + '\t' +fate[entry] + '\t' + Final + '\n')

