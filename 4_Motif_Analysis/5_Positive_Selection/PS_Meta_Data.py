#!/usr/env python3

'''
Example Command Line Input
python3 Infer_Ancestor_Chunks.py
'''

##########################################################
### Set environmental variables
##########################################################
import subprocess
import os

CURRENT_DIR = './'
PACK = ''
INPUTS = ['Chunk_All_KS.txt', 'Chunk_All_PS.txt']


##########################################################
### Set environmental variables for each ancestral node
##########################################################
for INSTANT in INPUTS:

	outer = CURRENT_DIR + 'Data/'
	file1 = outer + INSTANT
	file2 = outer + 'Final_Con_List_peak_Edit_R.txt'
	file3 = outer + 'Nfur_Peak_All.bed'
	file4 = outer + 'fish4_WGA.maf'
	file5 = outer + 'Interest_Motif_sort.bed'
	file6 = outer + 'temp_in.bed'
	file7 = outer + 'temp_out.bed' 

	seqer_home = outer + 'Peak_Seqs/'
	OUTPUT_F = outer + INSTANT[:-4] + '_Final.txt'

	NFUR_C = {}
	NFUR_P = {}
	Align = {}
	POS = {}
	Motifs = {}
	
	
	##########################################################
	### Parse positive selection data
	##########################################################
	with open(file1, 'r') as in1:
		for line in in1:
			line = line.rstrip('\n').split('\t')
			POS[line[0]] = (str(line[1]) + '\t' + str(line[2]) + '\t' + str(line[3]))

	print(INSTANT + ' Step 1 Complete')
	
	##########################################################
	### Parse peak conservation data
	##########################################################
	with open(file2, 'r') as in2:
			for line in in2:
				line = line.rstrip('\n').split('\t')
				peak = line[0]
				cons = line[1]
				typer = line[2]
				DES = line[3]
				gene = line[4]
				stat = line[5]
				Align[peak] = cons + '\t' + typer + '\t' + DES + '\t' + gene + '\t' + stat

	print(INSTANT + ' Step 2 Complete')
	
	##########################################################
	### Parse all peak bed-formatted coordinates
	##########################################################
	with open(file3, 'r') as in3:
		for line in in3:
			line = line.rstrip('\n').split('\t')
			chr = line[0]
			start = line[1]
			end = line[2]
			peak = line[3]
			if chr in NFUR_P.keys():
				NFUR_P[chr][peak] = [int(start), int(end)]
			else:
				NFUR_P[chr] = {}
				NFUR_P[chr][peak] = [int(start), int(end)]
			
	print(INSTANT + ' Step 3 Complete')
	
	##########################################################
	### Parse maf blocks and write all peak metadata
	##########################################################
	counter = 1
	with open(file4, 'r') as in4, open(OUTPUT_F, 'w') as final:
		final.write('Chromosome\tStart\tEnd\tName\tlength\tPeakPercent\tConservation\ttype\tDEStat\tgene\tparalog\tMotifs\tScore\tMismatches\tPercentMM\tPercentG\tPvalue\n')
		for line in in4:
			line = line.rstrip('\n')
			line = '\t'.join(line.split())
			line = line.split('\t')
			if line[0] == 's':
				name = line[1].split('.')
				if name[0] == 'Nfur':
					user = 0
					chr = name[1] + '.' + name[2]
					start = int(line[2])
					end = int(line[2]) + int(line[3])
					if chr in NFUR_P.keys():
						for peak in NFUR_P[chr].keys():
							if (NFUR_P[chr][peak][0] <= start <= NFUR_P[chr][peak][1]) or (NFUR_P[chr][peak][0] <= end <= NFUR_P[chr][peak][1]) or (start <=  NFUR_P[chr][peak][0] <= end) or (start <=  NFUR_P[chr][peak][1] <= end):
								if peak in NFUR_C:
									NFUR_C[peak] = NFUR_C[peak] + 1
								else:
									NFUR_C[peak] = 1
							
								########################################
								#Switch comment for chunks vs peaks
								########################################
								namer = peak
								#namer = peak + '.' + str(NFUR_C[peak])
								if NFUR_P[chr][peak][0] <= start:
									O_start = start
								else:
									O_start = NFUR_P[chr][peak][0]
								if NFUR_P[chr][peak][1] <= end:
									O_end = NFUR_P[chr][peak][1]
								else:
									O_end = end
								percent = str((O_end-O_start)/(NFUR_P[chr][peak][1]-NFUR_P[chr][peak][0]))
								user = 1
							
								tempy = chr + '\t' +  str(start) + '\t' + str(end) + '\n'
						
								with open(file6, 'w') as in6:
									in6.write(tempy)
						
								with open(file7, 'wt') as in71:
									subprocess.call([PACK + 'bedtools', 'intersect', '-a', file6, '-b', file5, '-wa', '-wb'], stdout = in71)
							
								with open(file7, 'r') as in72:
									for line in in72:
										line = line.rstrip('\n').split('\t')
										motif = line[6]
										if namer not in Motifs.keys():
											Motifs[namer] = []
											Motifs[namer].append(motif)
										else:
											Motifs[namer].append(motif)
										
									if namer not in Motifs.keys():
										Motifs[namer] = []
										Motifs[namer].append('NA')
								
								##########################################################
								### Clean workspace and condense motif data
								##########################################################	
								subprocess.call(['rm', file6])
								subprocess.call(['rm', file7])
								MOT = ','.join(Motifs[namer])
								Trunc = namer.split('.')
								Trunc = Trunc[0]
								ALLY = Align[Trunc]
								if namer in POS.keys():
									DATA = POS[namer]
									temp_data = DATA.split('\t')
									DATA1 = temp_data[0]
									DATA2 = temp_data[1]
									DATA3 = temp_data[2]
									DATA2H = int(DATA2)/(end-start)
								
								else:
									DATA1 = 'NA'
									DATA2 = 'NA'
									DATA3 = 'NA'
									DATA2H = 'NA'
								
								gappy = 0
								if os.path.isfile(seqer_home + namer + '.fa'):
									with open(seqer_home + namer + '.fa', 'r') as refer:
										for thing in refer:
											if '>nfur' in thing:
												gappy = 1
											elif gappy == 1:
												thing =  thing.rstrip('\n')
												thing = list(thing)
												denom = len(thing)
												numor = 0
												for entry in thing:
													if entry == '-':
														numor = numor + 1
												if denom != 0:
													gapper = str(numor/denom)
												else:
													gapper = 0
												break
								else:
									gapper = 0
									
								##########################################################
								### Combine data for writing to output file
								##########################################################
								nfur = chr + '\t' + str(start) + '\t' + str(end) + '\t' + namer + '\t' + str(end-start) + '\t' + percent + '\t' + ALLY + '\t' + MOT + '\t' + str(DATA1) + '\t' + str(DATA2) + '\t' + str(DATA2H) + '\t' + str(gapper) + '\t' + str(DATA3) + '\n'
								final.write(nfur)
								print('page ' + str(counter) + ' done')
								counter =  counter + 1

	print(INSTANT + ' Step 5 Complete')
