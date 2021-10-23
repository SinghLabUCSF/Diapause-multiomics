#!/usr/env python3

'''
Example Command Line Execution:
python3 TE_Bubble_Prep.py
'''

##########################################################
### Set environmental variables
##########################################################
import subprocess

CURRENT_DIR = './'
PACK = ''

file1 = CURRENT_DIR + 'Data/Nfur_True_DE.bed' #DE Peak file
file2 = CURRENT_DIR + 'Data/Nfur_Fake_DE.bed' #Fake Peak file
file3 = CURRENT_DIR + 'Data/Nfur_True_TE.bed' #TE File
file4 = CURRENT_DIR + 'Data/Nfur_True.bed' #genome-wide peak file

file5 = CURRENT_DIR + 'Data/TE_DE_Peak.bed' #DE-TE
file6 = CURRENT_DIR + 'Data/TE_All_Peak.bed' #GW-TE
file7 = CURRENT_DIR + 'Data/TE_DE_Fake.bed' #FK-TE

file8 = CURRENT_DIR + 'Data/TE_Class_Data.txt' #Final Class
file9 = CURRENT_DIR + 'Data/TE_Family_Data.txt' #Final Family

DEF = {}
GWF = {}
FKF = {}
TTF = {}
DEC = {}
GWC = {}
FKC = {}
TTC = {}

##########################################################
### Define Functions
##########################################################
def Count_TEs(a_file,a_fam,a_cal):

	'''
	This function counts the total number of TEs in a given family
	or class that appear in an overlap file between TEs and Peaks
	'''
	
	with open(a_file, 'r') as temper:
		counter = 0
		for line in temper:
			counter = counter + 1
			line = line.rstrip('\n').split('\t')
			if '/' in line[4]:
				TET = line[4].split('/')
				TEF = TET[1]
				TEC = TET[0]
			else:
				TEF = line[4]
				TEC = line[4]
			if TEF in a_fam:
				a_fam[TEF] = a_fam[TEF] + 1
			else:
				a_fam[TEF] = 1
				
			if TEC in a_cal:
				a_cal[TEC] = a_cal[TEC] + 1
			else:
				a_cal[TEC] = 1 
			
		a_fam['Total'] = counter
		a_cal['Total'] = counter	


##########################################################
### Generate a fake TE control region bed file
##########################################################
with open(file1, 'r') as in1, open(file2, 'w') as out1:
	for line in in1:
		line = line.rstrip('\n').split('\t')
		one = line[0]
		two = str(int(line[1]) + 10000)
		three = str(int(line[2]) + 10000)
		four = line[3]
		out1.write(one + '\t' + two + '\t' + three + '\t' + four + '\n')


##########################################################
### Intersect TEs with All,DE,and Control coordinates
##########################################################		
with open(file5, 'wt') as temp:
	subprocess.call([PACK + 'bedtools', 'intersect', '-a', file3, '-b', file1], stdout=temp)
	
with open(file6, 'wt') as temp:
	subprocess.call([PACK + 'bedtools', 'intersect', '-a', file3, '-b', file4], stdout=temp)
	
with open(file7, 'wt') as temp:
	subprocess.call([PACK + 'bedtools', 'intersect', '-a', file3, '-b', file2], stdout=temp)


##########################################################
### Count Family and Class overlap between TEs and Peaks
##########################################################	
Count_TEs(file5,DEF,DEC)
Count_TEs(file6,GWF,GWC)
Count_TEs(file7,FKF,FKC)
Count_TEs(file3,TTF,TTC)

for entry in TTF.keys():
	if entry not in GWF.keys():
		GWF[entry] = 0
	if entry not in DEF.keys():
		DEF[entry] = 0
	if entry not in FKF.keys():
		FKF[entry] = 0

for entry in TTC.keys():
	if entry not in GWC.keys():
		GWC[entry] = 0
	if entry not in DEC.keys():
		DEC[entry] = 0
	if entry not in FKC.keys():
		FKC[entry] = 0

##########################################################
### Summarizes TE counts as percentages for output file
##########################################################		
with open(file8, 'w') as classer, open(file9, 'w') as fammer:
	header = 'TE\tTot.Count\tTot.Per\tGW.Count\tGW.Per\tDE.Count\tDE.Per\tFK.Count\tFK.Per\n'
	classer.write(header)
	fammer.write(header)
	
	for item in TTC.keys():
		step0 = item + '\t' + str(TTC[item]) + '\t' + str(round((TTC[item]/TTC['Total']),3)) + '\t'
		step1 = step0 + str(GWC[item]) + '\t' + str(round((GWC[item]/GWC['Total']),3)) + '\t'
		step2 = step1 + str(DEC[item]) + '\t' + str(round((DEC[item]/DEC['Total']),3)) + '\t'
		step3 = step2 + str(FKC[item]) + '\t' + str(round((FKC[item]/FKC['Total']),3)) + '\n'
		classer.write(step3)
	
	for item in TTF.keys():
		step0 = item + '\t' + str(TTF[item]) + '\t' + str(round((TTF[item]/TTF['Total']),3)) + '\t'
		step1 = step0 + str(GWF[item]) + '\t' + str(round((GWF[item]/GWF['Total']),3)) + '\t'
		step2 = step1 + str(DEF[item]) + '\t' + str(round((DEF[item]/DEF['Total']),3)) + '\t'
		step3 = step2 + str(FKF[item]) + '\t' + str(round((FKF[item]/FKF['Total']),3)) + '\n'
		fammer.write(step3)	
			   