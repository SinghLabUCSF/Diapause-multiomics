#!/usr/env python3

'''
Example Command Line Entry
python3 Motif_Master_Cleanup.py
'''

#######################################################################################
##### Define Constants
#######################################################################################

import sys
import os

CURRENT_DIR = ''

file01 = CURRENT_DIR + 'Data/Peak_Master_0.0.txt'
file02 = CURRENT_DIR + 'Data/TE_Conservation.txt'
file03 = CURRENT_DIR + 'Data/TE_Motif_Overlap.bed'
Proccesses = ['Diapause', 'Other']

#######################################################################################
##### Define Functions
#######################################################################################
'''
This function converts all infinite fold change data to 100 as a maximum
'''
def Inf_correct(score):
	if score == 'Inf':
		return(100)
	elif score == 'NA':
		return('NA')
	else:
		return(score)

'''
This function generates a unique number id that denotes if a TE(s) overlaps with a given peak
and which species the TE(s) originate in.
'''		
def TE_Con_Count(a_peak, a_TE, a_dict):
	temp = 0
	if a_peak in a_dict.keys():
		if a_TE in a_dict[a_peak][1]:
			temp = temp + 10
		if a_TE in a_dict[a_peak][2]:
			temp = temp + 100
		if a_TE in a_dict[a_peak][3]:
			temp = temp + 1000
		if a_TE in a_dict[a_peak][4]:
			temp = temp + 10000	  
	return(temp)
	
'''
This function returns a list of aligned bases in Nfur as compared to other species, specifically
skipping over alignment gaps and unknown alignment states
'''	
def count_nfdex(a_seq):
	counter = 0
	temp = []
	tseq = ''
	for base in a_seq:
		if base != '-':
			temp.append(counter)
			tseq = tseq + base
		counter = counter + 1

	data = [tseq, temp]
	return(data)
		   
'''
This function extrat bases in non-Nfur species according to the exact base pair coordinates 
determine to not be gaps in NFUR in the function COUNT_NFDEX.
'''					
def apply_nfdex(a_seq, a_list):
	temp = ''
	for place in a_list:
		if a_seq[place] != '-':
			temp =  temp + a_seq[place]
		elif a_seq[place] == '-':
			temp = temp + 'N'
	
	return(temp)

#######################################################################################
##### Iterate through each
#######################################################################################
for file in os.listdir(CURRENT_DIR + 'Data/DEC/'):
	MOT = file[:-4]
	for Formula in Proccesses:
		file1 = CURRENT_DIR + 'Data/DEC/' + MOT + '.txt'
		file2 = CURRENT_DIR + 'Data/Final/' + MOT + '_' + Formula + '_Final.txt'
		file3 = CURRENT_DIR + 'Data/Logos/' + MOT + '_' + Formula + '_NF.txt' 
		file4 = CURRENT_DIR + 'Data/Logos/' + MOT + '_' + Formula + '_AA.txt'
		file5 = CURRENT_DIR + 'Data/Logos/' + MOT + '_' + Formula + '_AL.txt'
		file6 = CURRENT_DIR + 'Data/Logos/' + MOT + '_' + Formula + '_OL.txt'
		file7 = CURRENT_DIR + 'Data/Logos/' + MOT + '_' + Formula + '_DR.txt'

		CONS = {}
		TES = {}
		COMB = {}
		TIME = ['5.TE.Nfur', '6.TE.African', '7.TE.Killifish']
		Faster = 0
		begin = 0
	
		#######################################################################################
		##### Genrate Peak Conservation Dictionary
		#######################################################################################	

		with open(file01, 'r') as in0:
			for line in in0:
				line = line.rstrip('\n').split('\t')
				if line[0] != 'Nfur':
		
					if 'peak' in line[4]:
						CONS[line[0]] = 'Broad1'
					elif 'peak' in line[5]:
						CONS[line[0]] = 'Broad2'
					else:
						if 'peak' in line[3] and ('peak' in line[2] or 'peak' in line[1]):
							CONS[line[0]] = '4.Killifish'
						elif 'peak' in line[3]:
							CONS[line[0]] = '3.Annual'  
						elif 'peak' in line[2] or 'peak' in line[1]:
							CONS[line[0]] = '2.African'
						elif line[1] == 'NA' and line[2] == 'NA' and line[3] == 'NA' and line[4] == 'NA' and  line[5] == 'NA':
							CONS[line[0]] = '6.Unaligned'
						else:
							CONS[line[0]] = '1.Nfur-Specific'
					
		print(MOT + ' Volume 1 Complete')
 
 
		#######################################################################################
		##### Generate Motif-TE Overlap Dictionary
		#######################################################################################
				   
		with open(file03, 'r') as in0:
			for line in in0:
				line = line.rstrip('\n').split('\t')
				Peak = line[6]
				Motif = line[3]
				leng = int(line[2]) - int(line[1])
		
				if Motif == 'HRE':
					if leng == 15:
						Motif = 'HRE_ST'
					elif leng == 20:
						Motif = 'HRE_HP'

				elif Motif == 'COUP-TFII':
					if leng == 12:
						Motif = 'COUP-TFII_K5'
					elif leng == 8:
						Motif = 'COUP-TFII_AR'

				elif Motif == 'c-Myc':
					if leng == 8:
						Motif = 'c-Myc_LN'
					elif leng == 10:
						Motif = 'c-Myc_mE'

				elif Motif == 'PAX5':
					if leng == 14:
						Motif = 'PAX5_co'
					elif leng == 16:
						Motif = 'PAX5_no'
				TE =  line[12]
				if Motif == MOT:
					if int(line[9]) < int(line[1]) and int(line[2]) < int(line[10]):
						if Peak in COMB.keys():
							if TE not in COMB[Peak]:
								COMB[Peak].append(TE)
						else:
							COMB[Peak] = []
							COMB[Peak].append(TE)
					
		print(MOT + ' Volume 2 Complete')
		
		
		#######################################################################################
		##### Generate TE Site Conservation Dictionary
		#######################################################################################
					
		with open(file02, 'r') as in0:
			for line in in0:
				line = line.rstrip('\n').split('\t')
				TES[line[0]] = [line[1], line[2], line[3], line[4], line[5]]
		
		print(MOT + ' Volume 3 Complete')
		

		#######################################################################################
		##### Parse Each Motif Site
		#######################################################################################

		with open(file1, 'r') as in1, open(file2, 'w') as out1:
			out1.write('Motif\tPeak\tGene\tMethod\tPara\tConserve\tConfound\tTranslocation\tNfur\tAaus/Astr\tAlim\tOlat\tDrer\n')
			with open(file3, 'w') as seq1, open(file4, 'w') as seq2, open(file5, 'w') as seq3, open(file6, 'w') as seq4, open(file7, 'w') as seq5:
				for line in in1:
					line = line.rstrip('\n').split('\t')
					if line[0] != 'Motif' and line[3] == Formula:
						if 'NF' in line[5] or 'NF' in line[6]:
							Method = 'NA'
							mot = line[0]
							peak = line[1]
							gene = line[4]
							Trans = line[9]
							sco1 = Inf_correct(line[16])
							sco2 = Inf_correct(line[17])
							sco3 = Inf_correct(line[18])
							sco4 = Inf_correct(line[19])
							sco5 = Inf_correct(line[20])
					
							if 'NF' in line[5] and 'NF' in line[6]:
								ParStat = 'Both'
							elif 'NF' in line[5]:
								ParStat = 'N1'
							elif 'NF' in line[6]:
								ParStat = 'N2'
					
							Cont = 0
							if 'D' in line[11]:
								Cont = Cont + 10000
							if 'D' in line[12]:
								Cont = Cont + 1000
							if 'D' in line[13]:
								Cont = Cont + 100
							if 'D' in line[14]:
								Cont = Cont + 10
							if 'D' in line[15]:
								Cont = Cont + 1
							if Cont != 0:
								ConStat = 'Confound' + str(Cont)
							else:
								ConStat = 'Safe'
					 
						
							#################################################
							##### Identify Translocation Events
							################################################# 
							if peak in CONS.keys():
								SetStat = CONS[peak]
						
								if Trans == 'Detected' and line[24] == 'Bound' and SetStat == 'Broad1':
									TTran = 'Translocated'
									Method = '4.Translocated'  
								elif Trans == 'Detected' and line[24] == 'Bound':
									TTran = 'Translocated'
									Method = '4.Translocated'
								elif Trans == 'Detected':
									TTran = 'Translocated'
								else:
									TTran = 'None'
						
								if SetStat == 'Broad1' or SetStat == 'Broad2':
									SetStat = '5.Broad'
						
							else:
								SetStat = 'NA '
						
						
							#################################################
							##### Identify De Novo Mutation Events
							#################################################
							if 'Bound' not in line[22:26]:
								Method = '1.Nfur.Novel'
							elif 'Bound' in line [22] and 'Bound' not in line[23:26]:
								Method = '2.African.Novel'
							elif 'Bound' in line[22] and 'Bound' in line[23] and 'Bound' not in line[24:26]:
								Method = '3.Killifish.Novel'
						
					
							#################################################
							##### Identify Transposition Events
							#################################################   
							great = 0	
							if peak in COMB.keys():
								if len(COMB[peak]) > 1:
									for item in COMB[peak]:
										grabber = TE_Con_Count(peak, item, TES)
										great = great + grabber			 
								else:
									great = TE_Con_Count(peak, COMB[peak][0], TES)
						
								if len(str(great)) < 4:
									Method =  TIME[len(str(great))-1]
								else:
									Method = '8.TE.Conserved'
					
							if Method == 'NA':
								Method = '9.Conserved'
								
							tot1 = mot + '\t' + peak + '\t' + gene + '\t' + Method + '\t'
							tot2 = ParStat + '\t' + SetStat + '\t' + ConStat + '\t' + TTran + '\t'
							tot3 = str(sco1) + '\t' + str(sco2) + '\t' + str(sco3) + '\t' + str(sco4) + '\t' + str(sco5) + '\n'
							out1.write(tot1 + tot2 + tot3)
					
					
					
							#################################################
							##### Generate LOGO Flat Files
							#################################################
							if line[26] != 'NA' and line[26] != 'Unaligned' and ('NF' in line[5] or 'NF' in line[6]):
								master = count_nfdex(line[26])
								if begin == 0:
									begin = 1
									Storer = len(master[0])
								Faster = Faster + 1
								if len(master[0]) == Storer:
									seq1.write('>Motif' + str(Faster) + '\n')
									seq1.write(master[0] + '\n')
					
							if line[27] != 'NA' and line[27] != 'Unaligned' and ('NF' in line[5] or 'NF' in line[6]):
								aaus = apply_nfdex(line[27], master[1])
								if len(aaus) == Storer:
									seq2.write('>Motif' + str(Faster) + '\n')
									seq2.write(aaus + '\n') 
						  
							if line[28] != 'NA' and line[28] != 'Unaligned' and ('NF' in line[5] or 'NF' in line[6]):
								alim = apply_nfdex(line[28], master[1])
								if len(alim) == Storer:
									seq3.write('>Motif' + str(Faster) + '\n')
									seq3.write(alim + '\n')
						
							if line[29] != 'NA' and line[29] != 'Unaligned' and ('NF' in line[5] or 'NF' in line[6]):
								olat = apply_nfdex(line[29], master[1])
								if len(olat) == Storer:
									seq4.write('>Motif' + str(Faster) + '\n')
									seq4.write(olat + '\n')
						
							if line[30] != 'NA' and line[30] != 'Unaligned' and ('NF' in line[5] or 'NF' in line[6]):
								drer = apply_nfdex(line[30], master[1])
								if len(drer) == Storer:
									seq5.write('>Motif' + str(Faster) + '\n') 
									seq5.write(drer + '\n')	
						
		print(MOT + ' Volume 4 Complete')

								 