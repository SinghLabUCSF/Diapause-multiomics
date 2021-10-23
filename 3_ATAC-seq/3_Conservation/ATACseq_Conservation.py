#!/usr/env python3

'''
Example Command Line Input
python3 ATACseq_Conservation.py ScriptX

Script1: Convert Other species in multiple whole-genome alignments to the coordinates of the reference species (Lift-Over).
Script2: Generate simple reference file for the coordinates of consensus peaks of reference species (Bed-Format)
Script3: Remove reference peaks not in multi-alignment and determine the conservation status of peaks/peak type/and peak intensity across species.
Script4: Determine the amount of each genome covered in pairwise/multiple genome alignment and the block size break down of each.
Script5: Generate Single file with information on peak conservation, paralogs status, etc.
'''

##########################################################
### Set environmental variables
##########################################################

import sys
import os
import subprocess
RUNNER = sys.argv[1]
Current_Dir = './'
PACK = ''

##########################################################
### Set output file headers and variables
##########################################################

#Used in Script3
Nfur = Current_Dir + 'Data/Nfur_Nfur_merged_sorted.bed'
#Amount of peak overlap required to call conservation
require_list = [0.0, 0.25, 0.5]
#Number of Samples per Species
Checkerdex = {'Aaus':4, 'Astr':5, 'Alim':3, 'Olat':4, 'Drer':4} 
#Header for General Output
Header = 'Nfur\tAaus\tAstr\tAlim\tOlat\tDrer\n' 
#Header for Species One Samples
NfurO = 'Nfur_Peak\tNfur_Dia6D1\tNfur_Dia6D2\tNfur_Escape1\tNfur_Escape2\tNfur_kvOld2\tNfur_kvYoung1\tNfur_kvYoung2\tNfur_Dia1m1\tNfur_Dia1m2\tNfur_Dia1m3\tNfur_kvOld3\t'
#Header for Species Two Samples
AausO = 'Aaus_KV1\tAaus_KV2\tAaus_HB1\tAaus_KV3\t'
#Header for Species Three Samples
AstrO = 'Astr_HB1\tAstr_KVO1\tAstr_KVO2\tAstr_KVY1\tAstr_KVY2\t'
#Header for Species Four Samples
AlimO = 'Alim_1M_1\tAlim_Dev_1\tAlim_KV_1\t'
#Header for Species Five Samples
OlatO = 'Olat_stage_19.1\tOlat_stage_19.2\tOlat_stage_25.1\tOlat_stage_25.2\t'
#Header for Species Six Samples
DrerO = 'Drer_X8som1\tDrer_X8som2\tDrer_X48h1\tDrer_X48h2'
#Order for combining headers
binder = NfurO + AausO + AstrO + AlimO + OlatO + DrerO + '\n'

#Data Names to be Ignored (Set by Script3)
banned = []

#Effective Genome Sizes (Used by Script 4)
NF_true = 856808511
AA_true = 838719370
AL_true = 483220538
OL_true = 1368780147
DR_true = 733566086



##########################################################
### Choose script to run
##########################################################

if RUNNER == 'Script1':

	##########################################################
	### Script 1
	##########################################################
	##########################################################
	### Define Constants and Indexing Data Structure Categories
	##########################################################
	NF_SF = {'Aaus':{}, 'Astr':{}, 'Alim':{}, 'Olat':{}, 'Drer':{}}
	OT_SF = {'Aaus':{}, 'Astr':{}, 'Alim':{}, 'Olat':{}, 'Drer':{}}
	
	for file in  os.listdir(Current_Dir + 'Data/'):
		if 'rpkm' in file:
			Pfile = file
			spec = file[:-9]
			counter = 1
			
			if spec != 'Nfur':
				
				if spec == 'Astr':
					tru_spec = 'Aaus'
				else:
					tru_spec = spec

				#####################################################################################
				### Create NFUR to other fish Dictionaries
				#####################################################################################

				with open(Current_Dir + 'Data/Nfur.' + tru_spec + '.sing.maf', 'r') as refer, open(Current_Dir + 'Data/' + spec + '_maf_peaks.bed', 'w') as comper:
					#parase maf file entries for each pairwise alignment
					for line in refer:
						line = line.rstrip('\n')
						line = '\t'.join(line.split())
						line = line.split('\t')
						if line[0] == 's':
							name = line[1].split('.')
							#Save all Nfur entries from the alignment
							if name[0] == 'Nfur':
								NF_SF[spec]['maf' + str(counter)] = []
								NF_SF[spec]['maf' + str(counter)].append(name[1] + '.' + name[2])
								NF_SF[spec]['maf' + str(counter)].append(int(line[2]) + 1)
								NF_SF[spec]['maf' + str(counter)].append(int(line[2]) + int(line[3]))
								NF_SF[spec]['maf' + str(counter)].append(line[6])
					
							#Save other species entries from the alignment
							else:
								OT_SF[spec]['maf' + str(counter)] = []
								#Special Case for Alim or Alim TE file
								if spec == 'Alim' or spec == 'AlimTE':
									OT_SF[spec]['maf' + str(counter)].append(name[1] + '.' + name[2])
								else:
									OT_SF[spec]['maf' + str(counter)].append(name[1])	
								
								OT_SF[spec]['maf' + str(counter)].append(int(line[2]) + 1)
								OT_SF[spec]['maf' + str(counter)].append(int(line[2]) + int(line[3]))
								OT_SF[spec]['maf' + str(counter)].append(line[6])
							
								#Special Case for Alim or Alim TE file
								if spec == 'Alim' or spec == 'AlimTE':
									temp = (name[1] + '.' + name[2] + '\t' + str((int(line[2])+1)) + '\t' + str((int(line[2]) + int(line[3]))) + '\t' + 'maf' + str(counter) + '\n')
								else:
									temp = (name[1] + '\t' + str((int(line[2])+1)) + '\t' + str((int(line[2]) + int(line[3]))) + '\t' + 'maf' + str(counter) + '\n')
								comper.write(temp)
								counter = counter + 1
				print('Finished cataloging blocks for Nfur vs ' + spec)

				#####################################################################################
				### Convert Narrow Peak to BED files
				#####################################################################################

				with open(Current_Dir + 'Data/' + Pfile, 'r') as trans1, open(Current_Dir + 'Data/' + spec + '.bed', 'w') as trans2:
					#Parse peak files for use
					for line in trans1:
						line = line.rstrip('\n')
						line = '\t'.join(line.split())
						line = line.split('\t')
						if line[0] == 'chromosome': #Ignore header line
							pass
						else:
							line = str(line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + line[3])
							trans2.write(line + '\n')
				print('Finished aligning peaks for Nfur vs. ' + spec)
				
	print('Step 1 & 2 Complete')
			
	#####################################################################################
	### Sort Peaks Files
	#####################################################################################

	for file in os.listdir(Current_Dir + 'Data/'):
		if '.bed' in file:
			with open(Current_Dir + 'Data/' + file[:-4] + '_sorted.bed', 'wt') as hold:			   
				subprocess.call([PACK + 'bedtools', 'sort', '-i', Current_Dir + 'Data/' + file], stdout=hold)

	print('Step 3 Complete')

	#####################################################################################
	### Overlap Peaks and Multi-Aligned Genome Fragments
	#####################################################################################

	for file in os.listdir(Current_Dir + 'Data/'):
		if 'maf_peaks_sorted.bed' in file:
			with open(Current_Dir + 'Data/' + file[:-11] + '_intersects.bed', 'wt') as sitter:
				subprocess.call([PACK + 'bedtools', 'intersect', '-wa', '-wb', '-a', Current_Dir + 'Data/' + file, '-b', Current_Dir + 'Data/' + file[:-21] + '_sorted.bed'], stdout=sitter)

	print('Step 4 Complete')

	#####################################################################################
	### Calculate NFUR Coordinates of Peaks
	#####################################################################################

	for file in os.listdir(Current_Dir + 'Data/'):
		if '_intersects.bed' in file:
			print(file)
			with open(Current_Dir + 'Data/' + file, 'r') as close, open(Current_Dir + 'Data/' + file[:-21] + '_Nfur.bed', 'w') as output:
				#Parse data from intersected bed files
				for line in close:
					line = line.rstrip('\n').split('\t')
					srch = line[3]
					pstr = line[5]
					pend = line[6]
					pidx = line[7]
					sstr = line[1]
					send = line[2]
					checker = OT_SF[file[:-25]][srch][3]
					verify = NF_SF[file[:-25]][srch][3]
		
					########################
					### Find Start Index ###
					########################
					if int(pstr) <= int(sstr):
						true_str = 0
						mover1 = int(pstr) - int(sstr) + 1
					else:
						mover1 = int(pstr) - int(sstr) + 1
						step = 0
						limit = 0
						for base in checker:
							step = step + 1
							if base == '-':
								pass
							else:
								limit = limit + 1
				
							if limit == mover1:
								break
						true_str = step - 1
			
			
					########################	
					### Find End Index   ###
					########################
					if int(pend) >= int(send):
						true_end = 0
					else:
						mover2 = (int(pend) - int(pstr)) + mover1 + 1
						step = 0
						limit = 0
						for base in checker:
							step = step + 1
							if base == '-':
								pass
							else: 
								limit = limit + 1
				
							if limit == mover2:
								break
						true_end = step - 1
			   
					########################	
					### Verify Start	 ###
					########################
					marcher = 0
					NT_str = 0
					while marcher < true_str:
						if NF_SF[file[:-25]][srch][3][marcher] == '-':
							pass
						else:
							NT_str = NT_str + 1
						marcher = marcher + 1
				
					########################	
					### Verify End	   ###
					########################
					if true_end == 0:
						pass
					else:
						marcher2 = 0
						NT_end = 0
						while marcher2 < true_end:
							if NF_SF[file[:-25]][srch][3][marcher2] == '-':
								pass
							else:
								NT_end = NT_end + 1
							marcher2 = marcher2 + 1
			
					########################	
					### Build Bed Entry  ###
					########################		
					texter1 = NF_SF[file[:-25]][srch][0]
					texter2 = str(int(NF_SF[file[:-25]][srch][1]) + NT_str)
					if true_end == 0:
						texter3 = str(NF_SF[file[:-25]][srch][2])
					else:
						texter3 = str(int(NF_SF[file[:-25]][srch][1]) + NT_end)
					texter4 = pidx
					output.write(texter1 + '\t' + texter2 + '\t' + texter3 + '\t' + texter4 + '\n')
		
	print('Step 5 Complete')

	#####################################################################################
	### Merge Split Peaks
	#####################################################################################
	
	#Sort bed files before merging
	for file in os.listdir(Current_Dir + 'Data/'):
		if '_Nfur.bed' in file:
			with open(Current_Dir + 'Data/' + file[:-4] + '_sorted.bed', 'wt') as presort:
				subprocess.call([PACK + 'bedtools', 'sort', '-i', Current_Dir + 'Data/' + file], stdout=presort)

	for file in os.listdir(Current_Dir + 'Data/'):
		if 'Nfur_sorted.bed' in file: 
			merger = {}
			lister = {}
			with open(Current_Dir + 'Data/' + file, 'r') as trans3, open(Current_Dir + 'Data/' + file[:-11] + '_merged.bed', 'w') as trans4:
				#Parse data from the bed file
				for line in trans3:
					line = line.rstrip('\n').split('\t')
					peakN = line[3] 
					CHR = line[0]
					Cstart = int(line[1])
					Cend = int(line[2])
		
					if peakN in lister.keys():
						holder = lister[peakN]
						PNU = peakN + 'X' + str(holder)
						ICHR = merger[PNU][0]
						ISTR = merger[PNU][1]
						IEND = merger[PNU][2]
						#Add peak prefix allowing for multiple entries for split peaks across blocks
						if CHR == ICHR:
							if Cstart - int(IEND) < 501:
								merger[PNU] = [ICHR, ISTR, line[2], PNU]
							else:
								lister[peakN] = lister[peakN] + 1
								PNU = peakN + 'X' + str(holder + 1)
								merger[PNU] = [line[0], line[1], line[2], PNU]
						else:
							lister[peakN] = lister[peakN] + 1
							PNU = peakN + 'X' + str(holder + 1)
							merger[PNU] = [line[0], line[1], line[2], PNU]
					else:
						lister[peakN] = 1
						PNU = peakN + 'X1'
						merger[PNU] = [line[0], line[1], line[2], PNU]
				
				#Write all peaks to file
				for entry in merger.keys():
					texterF = '\t'.join(merger[entry])
					trans4.write(texterF + '\n')

	print('Step 6 Complete')
   
	#####################################################################################
	### Sort Other Fish Peaks by NFUR Chromosome
	#####################################################################################

	for file in os.listdir(Current_Dir + 'Data/'):
		if '_Nfur_merged.bed' in file:
			with open(Current_Dir + 'Data/' + file[:-4] + '_sorted.bed', 'wt') as final:
				subprocess.call([PACK + 'bedtools', 'sort', '-i', Current_Dir + 'Data/' + file], stdout=final)			  

	print('Step 7 Complete')
	
	#####################################################################################
	### Remove intermediate files
	#####################################################################################
	
	for file in os.listdir(Current_Dir + 'Data/'):
		if '.bed' in file and 'Nfur_merged_sorted.bed' not in file and 'Master' not in file:
			subprocess.call(['rm', Current_Dir + 'Data/' + file])
	
	print('Step 8 Complete')
	


if RUNNER == 'Script2':
	
	##########################################################
	### Script 2
	##########################################################
	##########################################################
	### NFUR TRUE BED
	##########################################################
	
	with open(Current_Dir + '/Data/Nfur_rpkm.txt', 'r') as in1, open(Current_Dir + 'Data/Nfur_Peak_All.bed', 'w') as out1:
		for line in in1:
			line = line.rstrip('\n').split('\t')
			Peak = line[3]
			Chr = line[0]
			Start = line[1]
			End = line[2]
			if line[0] != 'chromosome':
				out1.write(Chr + '\t' + Start + '\t' + End + '\t' + Peak + '\n')
	
	print('Step 1 Complete')
	

if RUNNER == 'Script3':	

	##########################################################
	### Script 3
	##########################################################
	##########################################################
	### Convert MAF blocks to bed coordinate file
	##########################################################

	with open(Current_Dir + 'Data/fish4_WGA.maf', 'r') as in1, open(Current_Dir + 'Data/Peak_Exclusion.bed', 'w') as out1:
		counter = 0
		for line in in1:
			if line[0] != '#' and 'Nfur' in line:
				line = line.rstrip('\n')
				line = '\t'.join(line.split())
				line = line.split('\t')
				chr = line[1]
				chr = chr[5:]
				start = str(line[2])
				finish = str(int(line[2]) + int(line[3]))
				counter = counter + 1
				out1.write(chr + '\t' + start + '\t' + finish + '\ttemp' + str(counter) + '\n')
	
	print('Step 1 Complete')
			
	##########################################################
	### Overlap Coordinate Files
	##########################################################		 
			
	with open(Current_Dir + 'Data/Peak_Exclusion_Merged.bed', 'wt') as out2:
		subprocess.call([PACK + 'bedtools', 'intersect', '-wo', '-a', Current_Dir + 'Data/Nfur_Peak_All.bed', '-b', Current_Dir + 'Data/Peak_Exclusion.bed'], stdout=out2)
		
	print('Step 2 Complete')

	##########################################################
	### Finalize Coordinate Converstion
	########################################################## 
 
	with open(Current_Dir + 'Data/Nfur_Peak_All.bed', 'r') as in2, open(Current_Dir + 'Data/Peak_Exclusion_Merged.bed', 'r') as in3:
		over = []
		for line in in3:
			line = line.rstrip('\n').split('\t')
			peak = line[3]
			over.append(peak)
		
		for line in in2:
			line = line.rstrip('\n').split('\t')
			peak = line[3]
			if peak not in over:
				banned.append(peak)
				
	
	print('Step 3 Complete')
	
	##########################################################
	### Generate Libraires Holding Peak Data
	##########################################################
	
	#Species Included data structures
	for require in require_list:
		Codex = {'Nfur':{}, 'Aaus':{}, 'Astr':{}, 'Alim':{}, 'Olat':{}, 'Drer':{}}
		Minidex = {'Nfur':{},  'Aaus':{}, 'Astr':{}, 'Alim':{}, 'Olat':{}, 'Drer':{}} 
		Ultradex = {'Aaus':{}, 'Astr':{}, 'Alim':{}, 'Olat':{}, 'Drer':{}}	
	
		with open(Nfur, 'w') as input:
			for file in os.listdir(Current_Dir + 'Data/'):
				if 'rpkm' in file: 
					with open(Current_Dir + 'Data/' + file, 'r') as temp:
						locator = file[:4]
						for line in temp:
							line = line.rstrip('\n')
							line = '\t'.join(line.split())
							line = line.split('\t')

							if line[0] == 'chromosome':
								pass
							else:
								ider = line[3]
								looper = len(line)
								hopper = 11
								Codex[locator][ider] = []
	
								while hopper < looper:
									Codex[locator][ider].append(line[hopper])
									hopper = hopper + 1
			 
							if locator == 'Nfur' and line[0] != 'chromosome':
								input.write(line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + line[3] + '\n')	   
					
				if 'annotations' in file:
					with open(Current_Dir + 'Data/' + file, 'r') as temp:
						locator = file[:4]
						for line in temp:
							line = line.rstrip('\n').split('\t')
				
							if line[0] == 'names':
								truth = line.index('annotation')
							else:
								locator = file[:4]
								ider = line[0]
								Minidex[locator][ider] = line[truth]

	##########################################################
	### Relate Nfur Peaks to all Other species
	##########################################################
	
		for file in os.listdir(Current_Dir + 'Data/'):
			if Current_Dir + 'Data/' + file != Nfur and 'Nfur_merged_sorted.bed' in file:
				#Find overlap 
				with open(Current_Dir + 'Data/' + file[:4] + '_inter.bed', 'wt') as sitter:
					subprocess.call([PACK + 'bedtools', 'intersect', '-wo', '-a', Current_Dir + 'Data/Nfur_Peak_All.bed', '-b', Current_Dir + 'Data/' + file], stdout=sitter)
				#Parse overlap file data
				with open(Current_Dir + 'Data/' + file[:4] + '_inter.bed', 'r') as newly:
					local = file[:4]
					for line in newly:
						line = line.rstrip('\n').split('\t')
						killi = line[3]
						fish = line[7]
						over = float(line[-1])
						start = float(line[1])
						finish = float(line[2])
						size = finish - start
						percent = over/size
						
						#Include peak info if overlap between peaks meets criteria
						if percent >= require:
							if killi in Ultradex[local].keys():
								if over > Ultradex[local][killi][1]:
									Ultradex[local][killi] = [fish, over]
							else:
								Ultradex[local][killi] = [fish, over]

		print('Peaks compared across Nfur and ' + file[:4] + ' at ' + str(require) + ' overlap')

		##########################################################
		### Write all peak data to master files
		##########################################################
		
		#Generate master conservation files and write their headers			
		with open(Nfur, 'r') as last, open(Current_Dir + 'Data/Peak_Master_' + str(require) + '.txt', 'w') as out1, open(Current_Dir + 'Data/RPKM_Master_' + str(require) + '.txt', 'w') as out2, open(Current_Dir + 'Data/Annotate_Master_' + str(require) + '.txt', 'w') as out3:
			out1.write(Header)
			out2.write(binder)
			out3.write(Header)
	
			for line in last:
				line = line.rstrip('\n').split('\t')
				peak = line[3]
		
				peaker = []
				rpkmer = []
				annotator = []
				
				#Begin entry for each file with Nfur peak name
				peaker.append(peak)
				rpkmer.append(peak)
				annotator.append(peak)
				annotator.append('NA')
				for item in Codex['Nfur'][peak]:
					rpkmer.append(item)
				
				#Bridge across adjacent split peaks
				for entry in Ultradex:
					if peak in Ultradex[entry].keys():
						bridge = Ultradex[entry][peak][0]
						peaker.append(bridge)
						bridge = bridge.split('X')
						bridge = bridge[0]
						
						#Verify bridge and remove banned peaks
						if bridge in Minidex[entry].keys():
							ant = Minidex[entry][bridge]
							annotator.append(ant)
						elif peak in banned:
							annotator.append('NA')
						else:
							annotator.append('None')
				
						for item in Codex[entry][bridge]:
							rpkmer.append(item)
			
					#Remove non bridged banned peaks
					elif peak in banned:
						peaker.append('NA')
						annotator.append('NA')
						repper = Checkerdex[entry]
						lenner = 1			   
						while lenner <= repper:
							rpkmer.append('NA')
							lenner = lenner + 1
					
					#Make peak as not conserved if all other options exhausted	
					else:	  
						peaker.append('None')
						annotator.append('None')
						repper = Checkerdex[entry]
						lenner = 1			   
						while lenner <= repper:
							rpkmer.append('None')
							lenner = lenner + 1
				
				#Join together data and write to master files
				one = '\t'.join(peaker)
				two = '\t'.join(rpkmer)
				three = '\t'.join(annotator)
				out1.write(one + '\n')
				out2.write(two + '\n')
				out3.write(three + '\n')
			
		print('All species intergrated at ' + str(require) + ' overlap')
		
	print('Step 4, 5, & 6 Complete')
	
	
	
	
if RUNNER == 'Script4':	

	##########################################################
	### Script 4
	##########################################################
	##########################################################
	### Set up data structures for script
	##########################################################
	
	NF_tot = 0
	AA_tot = 0
	AL_tot = 0
	OL_tot = 0
	DR_tot = 0
	BINS = {'BP':0, 'TBP':0, 'HBP':0, 'KB':0, 'TKB':0, 'HKB':0, 'MB':0, 'TMB':0}
	
	##########################################################
	### Catalog the size of each alignment block and total genome coverage
	##########################################################

	with open(Current_Dir + './Data/fish4_WGA.maf', 'r') as in1:
		for line in in1:
			#Parse Nfur entries and add to genome coverage
			if line[0] != '#' and 'Nfur' in line:
				line = line.rstrip('\n')
				line = '\t'.join(line.split())
				line = line.split('\t')
				length = int(line[3])
				NF_tot = NF_tot + length
				
				#Check which size bin each block fits on
				if 0 < length < 10:
					BINS['BP'] = BINS['BP'] + 1
				elif 11 < length < 100:
					BINS['TBP'] = BINS['TBP'] + 1
				elif 101 < length < 1000:
					BINS['HBP'] = BINS['HBP'] + 1
				elif 1001 < length < 10000:
					BINS['KB'] = BINS['KB'] + 1
				elif 10001 < length < 100000:
					BINS['TKB'] = BINS['TKB'] + 1
				elif 100001 < length < 1000000:
					BINS['HKB'] = BINS['HKB'] + 1
				elif 1000001 < length < 10000000:
					BINS['MB'] = BINS['MB'] + 1
				elif 10000001 < length:
					BINS['TMB'] = BINS['TMB'] + 1		
			
			#Parse Aaus entries and add to genome coverage
			elif line[0] != '#' and 'Aaus' in line:
				line = line.rstrip('\n')
				line = '\t'.join(line.split())
				line = line.split('\t')
				length = int(line[3])
				AA_tot = AA_tot + length
			
			#Parse Alim entries and add to genome coverage
			elif line[0] != '#' and 'Alim' in line:
				line = line.rstrip('\n')
				line = '\t'.join(line.split())
				line = line.split('\t')
				length = int(line[3])
				AL_tot = AL_tot + length
			
			#Parse Olat entries and add to genome coverage
			elif line[0] != '#' and 'Olat' in line:
				line = line.rstrip('\n')
				line = '\t'.join(line.split())
				line = line.split('\t')
				length = int(line[3])
				OL_tot = OL_tot + length
			
			#Parse Drer entries and add to genome coverage
			elif line[0] != '#' and 'Drer' in line:
				line = line.rstrip('\n')
				line = '\t'.join(line.split())
				line = line.split('\t')
				length = int(line[3])
				DR_tot = DR_tot + length
		
	print('Step 1 Complete')
				
	##########################################################
	### Output all cataloged data
	##########################################################
		
	with open(Current_Dir + 'Data/Coverage_Stats.txt', 'w') as out1, open(Current_Dir + 'Data/Size_Stats.txt', 'w') as out2:
		out1.write('Nfur' + '\t' + str(NF_tot/NF_true) + '\n')
		out1.write('Aaus' + '\t' + str(AA_tot/AA_true) + '\n')
		out1.write('Alim' + '\t' + str(AL_tot/AL_true) + '\n')
		out1.write('Olat' + '\t' + str(OL_tot/OL_true) + '\n')
		out1.write('Drer' + '\t' + str(DR_tot/DR_true) + '\n')
	
		for item in BINS.keys():
			out2.write(item + '\t' + str(BINS[item]) + '\n')

	print('Step 2 Complete')



if RUNNER == 'Script5':	
	
	##########################################################
	### Script 5
	##########################################################
	##########################################################
	### Set up data structures
	##########################################################

	peak_2_type = {}
	peak_2_gene = {}
	classer = {}
	ager = {}
	UPPER = []
	DOWNER = []

	##########################################################
	### Catalog feature types for each peak
	##########################################################

	with open(Current_Dir + '/Data/Nfur_annotations.txt', 'r') as in1:
		#Parse data from annotation files to get feature data
		for line in in1:
			line = line.rstrip('\n').split('\t')
			peak = line[0]
			gene = line[25]
			tip = line[19]
			tip = tip.split(' ')
			tip = tip[0]
			if tip == '5\'' or tip == '3\'':
				tip = 'UTR'
			elif tip == 'Distal' or tip == 'Downstream':
				tip = 'Intergenic'
			#Save feature data
			peak_2_type[peak] = tip
			peak_2_gene[peak] = gene
	
	print('Step 1 Complete')
	
	##########################################################
	### Catalog paralog status for each peak
	##########################################################
		
	with open(Current_Dir + 'Data/Paralog_List.txt', 'r') as in2:
		for line in in2:
			#Parse paralog data for specialization information
			line = line.rstrip('\n').split('\t')
			gene1 = line[0]
			gene2 = line[1]
			stat = line[2]
			age = line[3]
			
			#Check for specialization
			if stat == 'NeoF':
				classer[gene1] = 'N1'
				classer[gene2] = 'N2'
			
			#Check for paralog duplication age
			if age == 'Ancient vertebrate' and gene1 not in ager.keys() and gene2 not in ager.keys():
				ager[gene1] = 'V'
				ager[gene2] = 'V'
		
			elif age == 'Recent all fish':
				if gene1 in ager.keys() and ager[gene1] == 'Ancient vertebrate':
					ager[gene1] = 'F'
					ager[gene2] = 'F'
			
				elif gene2 in ager.keys() and ager[gene2] == 'Ancient vertebrate':
					ager[gene1] = 'F'
					ager[gene2] = 'F'
			
				elif gene1 not in ager.keys() and gene2 not in ager.keys():
					ager[gene1] = 'F'
					ager[gene2] = 'F'
			
			elif age == 'Recent all killifish':
				ager[gene1] = 'K'
				ager[gene2] = 'K' 
				
	print('Step 2 Complete')
			
	##########################################################
	### Catalog DE status for each peak
	##########################################################
		
	with open(Current_Dir + 'Data/Master_DE_Dia_Up_Peak.bed', 'r') as in3:
		for line in in3:
			line = line.rstrip('\n').split('\t')
			peak = line[3]
			UPPER.append(peak)
		
	with open(Current_Dir + 'Data/Master_DE_Dia_Down_Peak.bed', 'r') as in4:
		for line in in4:
			line = line.rstrip('\n').split('\t')
			peak = line[3]
			DOWNER.append(peak)
			
	print('Step 3 Complete')
			
	##########################################################
	### Generate files for alignment and peak conservation
	##########################################################
	
	#Generate output for both aligned region and peak conservation		
	method = ['Peak', 'Align']
	for TYPER in method:
	
		with open(Current_Dir + 'Data/Peak_Master_0.0.txt', 'r') as in5, open(Current_Dir + 'Data/Final_Con_List_' + TYPER + '_0.0.txt', 'w') as out1:
			header = ('Peak\tCon\tType\tDE\tGene\tStat\tage\n') 
			out1.write(header)
		
			for line in in5:
				line = line.rstrip('\n').split('\t')
				if line[0] != 'Nfur':
					Peak = line[0]
					checker = Peak.split('_')
					checker = int(checker[1])
					#Do not keep non-chromosomal peaks in file as they are not in the mutli-alignment
					if checker <= 46951:
					
						##########################################################
						### Process all metadata for peak conservation
						##########################################################
						
						if TYPER == 'Peak':
							if 'peak' in line[1] or 'peak' in line[2]:
								if 'peak' in line[3]:
									if 'peak' in line[4] or 'peak' in line[5]:
										Stat = '3.Broadly-Conserved'
									else:
										Stat = '2.Killifish-Specific'
								elif 'peak' in line[4] or 'peak' in line[5]:
									Stat = '3.Broadly-Conserved'
								else:
									Stat = '2.Killifish-Specific'		
							elif 'peak' in line[3] and 'peak' not in line[4] and 'peak' not in line[5]:
								Stat = '1.Nfur-Specific'
							elif 'peak' not in line[3] and ('peak' in line[4] or 'peak' in line[5]):
								Stat = '3.Broadly-Conserved'
							elif line[1] == 'NA' and line[2] == 'NA' and line[3] == 'NA' and line[4] == 'NA' and line[5] == 'NA':
								Stat = '1.Nfur-Specific'
							else:
								Stat = '1.Nfur-Specific'
			
						##########################################################
						### Process all metadata for alignment conservation
						##########################################################
						
						elif TYPER == 'Align':
							if line[1] != 'NA' or line[2] != 'NA':
								if line[3] != 'NA':
									if line[4] != 'NA' or line[5] != 'NA':
										Stat = '3.Broadly-Conserved'
									else:
										Stat = '2.Killifish-Specific'
								elif line[4] != 'NA' or line[5] != 'NA':
									Stat = '3.Broadly-Conserved'
								else:
									Stat = '2.Killifish-Specific'
							elif line[3] != 'NA' and line[4] == 'NA' and line[5] == 'NA':
								Stat = '1.Nfur-Specific'
							elif line[3] == 'NA' and (line[4] != 'NA' or line[5] != 'NA'):
								Stat = '3.Broadly-Conserved'
							elif line[1] == 'NA' and line[2] == 'NA' and line[3] == 'NA' and line[4] == 'NA' and line[5] == 'NA':
								Stat = '1.Nfur-Specific'
							else:
								Stat = '1.Nfur-Specific'
								
						##########################################################
						### If in missing data and output all data
						##########################################################
				
						if Peak in peak_2_type.keys():	
							tipper = peak_2_type[Peak]
						else:
							tipper = 'Unknown'
			
						if Peak in peak_2_gene.keys():
							genie = peak_2_gene[Peak]
						else:
							genie = 'Unknown'
				
						if Peak in UPPER:
							DER = 'DE_Dia'
						elif Peak in DOWNER:
							DER = 'DE_Dev'
						else:
							DER = 'Non_DE'
			
						if genie in classer.keys():
							grouper = classer[genie]
						else:
							grouper = 'Other'
			
						if genie in ager.keys():
							timer = ager[genie]
						else:
							timer = 'U'
			
						out1.write(Peak + '\t' + Stat + '\t' + tipper + '\t' + DER + '\t' + genie + '\t' + grouper + '\t' + timer + '\n')
						
	print('Step 4 Complete')

	