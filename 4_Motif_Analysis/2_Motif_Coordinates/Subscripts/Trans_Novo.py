#!/usr/env python3

'''
This Script is initiated by the script TN_Submission.sh which in turn
is initiated by the script MasterNT.sh
'''

###############################################################
#### Set Constant Variables
###############################################################
import sys
import math

CURRENT_DIR = './'

temp_mot = sys.argv[1]

file1 = CURRENT_DIR + 'Data/fish4_WGA.maf'
file2 = CURRENT_DIR + 'Data/Nfur_Para.txt'
file3 = CURRENT_DIR + 'Data/Nfur_Gene_Peak_Overlap.bed'
file4 = CURRENT_DIR + 'Data/master_motif_file.motif'
file5 = CURRENT_DIR + 'Data/Interest_Motif_sort.bed'
file6 = CURRENT_DIR + 'Data/Peak_Master_0_Final.txt'
file7 = CURRENT_DIR + 'Data/AllMotifs_ConsensusPeaks_Aaul.txt'
file8 = CURRENT_DIR + 'Data/AllMotifs_ConsensusPeaks_Astr.txt'
file9 = CURRENT_DIR + 'Data/AllMotifs_ConsensusPeaks_Alim.txt'
file10 = CURRENT_DIR + 'Data/AllMotifs_ConsensusPeaks_Olat.txt'
file11 = CURRENT_DIR + 'Data/AllMotifs_ConsensusPeaks_Drer.txt'


mastery = {'Aaus':{}, 'Astr':{}, 'Alim':{}, 'Olat':{}, 'Drer':{}}
forge = {'Aaus':{}, 'Astr':{}, 'Alim':{}, 'Olat':{}, 'Drer':{}}
IND_2_MOT = {}
maffer = {}
maccer = {}
P1_P2 = {}
P2_P1 = {}
P1_TY = {}
P2_TY = {}
PK_GN = {}
GN_CH = {}
Scorer = {}
final_hold = {}
Interests = {temp_mot:0}
COMP = {'A':'T', 'a':'t', 'T':'A', 't':'a', 'C':'G', 'c':'g', 'G':'C', 'g':'c', 'N':'N', 'n':'n', '-':'-'}
boop = 4
True_move = 0


###############################################################
#### Define Functions
###############################################################
'''
This function traverses the bases of a motif sequence or aligned region sequence and calculates
a likelihood of binding score. This score is then compared to the minimum binding score to determine
if the sequence can or cannot bind the given transcription factor
'''
def score_maker(a_case,a_memory,a_score,a_Scorer,a_motif,a_number):
	if a_case[a_number] == 'NA':
		a_score = 'NA'
		return(['NA','NA'])
	else:
		badder = 0
		trace = 0
		fix = 0
		for letter in a_case[a_number]:
			if trace in a_memory:
				fix = fix + 1
			elif letter == '-':
				badder =  badder + 1
			else:
				a_score = a_score * a_Scorer[a_motif][trace-fix][letter]
			trace = trace + 1
		
		if a_score == 1:
			return('NA', 'NA')
		elif a_score == 0:
			return('Unbound', 'Inf')
			
		else:
			if badder != 0:
				a_score = a_score/(10*badder)
			if len(memory) != 0:
				a_score = a_score/(10*len(a_memory))
		
			if a_score == 1:
				a_score = 'NA'
				return('NA', 'NA')
			elif a_score == 0:
				a_score = 'Inf'
				return('Unbound', 'Inf')
			else: 
				a_score = -1 * math.log(a_score/(1-a_score))

				if a_score < float(a_Scorer[a_motif]['Thresh']):
					return('Bound', a_score)
				elif a_score >= float(a_Scorer[a_motif]['Thresh']):
					return('Unbound', a_score)


'''
This function searches through the combined motif list and extracts just the current motif
of interest for a given ATACseq peak and also adjust the name if the motif contains a 
duplicate name within the file.
'''
def Assign_Motifs(A_spec, A_file, A_dic):
	with open(A_file, 'r') as A_inner:
		for line in A_inner:
			line = line.rstrip('\n').split('\t')
			if 'PeakID' in line[0]:
				counter = 0
				for item in line[21:]:
					temp = item.split('(')
					temp = temp[0]
					if temp == 'HRE':
						if 'Striatum' in item:
							temp = 'HRE_ST'
						elif 'HepG2' in item:
							temp = 'HRE_HP'
					elif temp == 'COUP-TFII':
						if 'K5' in item:
							temp = 'COUP-TFII_K5'
						elif 'Artia' in item:
							temp = 'COUP-TFII_AR'
					elif temp == 'c-Myc':
						if 'LNCAP' in item:
							temp = 'c-Myc_LN'
						elif 'mES' in item:
							temp = 'c-Myc_mE'
					elif temp == 'PAX5':
						if 'condensed' in item:
							temp = 'PAX5_co'
						elif 'condensed' not in item:
							temp = 'PAX5_no'
					elif temp == 'THRb':
						if 'Liver' in item:
							temp = 'THRb_LV'
						elif 'Hep' in item:
							temp = 'THRb_HP'
					elif temp == 'FOXA1':
						if 'MC' in item:
							temp = 'FOXA1_MC'
						elif 'LNCAP' in item:
							temp = 'FOXA1_LN'

					IND_2_MOT[counter + 21] = temp
					counter = counter + 1
			else:
				anchor = line[0]
				for item in line[21:]:
					if item != '':
						site = IND_2_MOT[line.index(item)]
						if anchor not in A_dic[A_spec].keys():
							A_dic[A_spec][anchor] = []
							A_dic[A_spec][anchor].append(site)
						else:
							A_dic[A_spec][anchor].append(site)


###############################################################
#### Generate Maf Sequence Dictionary
###############################################################
with open(file1, 'r') as in1:
	for line in in1:
		if line[0] != '#' and 'Nfur' in line:
			if boop != 4:
				if boop != 3:
					if boop != 2:
						if boop != 1:
							maffer[Tchr][bottle].append('NA')
							maffer[Tchr][bottle].append('NA')
							maffer[Tchr][bottle].append('NA')
							maffer[Tchr][bottle].append('NA')
							maccer[Tchr][bottle].append('NA')
							maccer[Tchr][bottle].append('NA')
							maccer[Tchr][bottle].append('NA')
							maccer[Tchr][bottle].append('NA')
						else:
							maffer[Tchr][bottle].append('NA')
							maffer[Tchr][bottle].append('NA')
							maffer[Tchr][bottle].append('NA')
							maccer[Tchr][bottle].append('NA')
							maccer[Tchr][bottle].append('NA')
							maccer[Tchr][bottle].append('NA')
					else:
						maffer[Tchr][bottle].append('NA')
						maffer[Tchr][bottle].append('NA')
						maccer[Tchr][bottle].append('NA')
						maccer[Tchr][bottle].append('NA')
				else:
					maffer[Tchr][bottle].append('NA')
					maccer[Tchr][bottle].append('NA') 
		
			line = line.rstrip('\n')
			line = '\t'.join(line.split())
			line = line.split('\t')
			Tchr = line[1]
			Tchr = Tchr[5:]
			if Tchr not in maffer.keys():
				maffer[Tchr] = {}
			if Tchr not in maccer.keys():
				maccer[Tchr] = {}
			start = str(line[2])
			seq = line[6]
			finish = str(int(line[2]) + int(line[3]))
			bottle = (Tchr, start, finish)
			maffer[Tchr][bottle] = []
			maffer[Tchr][bottle].append(seq)
			maccer[Tchr][bottle] = []
			boop = 0
			
		elif line[0] != '#' and 'Aaus' in line:
			line = line.rstrip('\n')
			line = '\t'.join(line.split())
			line = line.split('\t')
			chr = line[1]
			chr = chr[5:]
			start = str(line[2])
			seq = line[6]
			finish = str(int(line[2]) + int(line[3]))
			maffer[Tchr][bottle].append(seq)
			maccer[Tchr][bottle].append(chr)
			boop = 1
			
		elif line[0] != '#' and 'Alim' in line:
			if boop != 1:
				maffer[Tchr][bottle].append('NA')
				maccer[Tchr][bottle].append('NA')
				
			line = line.rstrip('\n')
			line = '\t'.join(line.split())
			line = line.split('\t')
			chr = line[1]
			chr = chr[5:]
			start = str(line[2])
			seq = line[6]
			finish = str(int(line[2]) + int(line[3]))
			maffer[Tchr][bottle].append(seq)
			maccer[Tchr][bottle].append(chr)
			boop = 2
			
		elif line[0] != '#' and 'Olat' in line:
			if boop != 2:
				if boop != 1:
					maffer[Tchr][bottle].append('NA')
					maffer[Tchr][bottle].append('NA')
					maccer[Tchr][bottle].append('NA')
					maccer[Tchr][bottle].append('NA')
				else:
					maffer[Tchr][bottle].append('NA')
					maccer[Tchr][bottle].append('NA')
					
			line = line.rstrip('\n')
			line = '\t'.join(line.split())
			line = line.split('\t')
			chr = line[1]
			chr = chr[5:]
			start = str(line[2])
			seq = line[6]
			finish = str(int(line[2]) + int(line[3]))
			maffer[Tchr][bottle].append(seq)
			maccer[Tchr][bottle].append(chr)
			boop = 3
			
		elif line[0] != '#' and 'Drer' in line:
			if boop != 3:
				if boop != 2:
					if boop != 1:
						maffer[Tchr][bottle].append('NA')
						maffer[Tchr][bottle].append('NA')
						maffer[Tchr][bottle].append('NA')
						maccer[Tchr][bottle].append('NA')
						maccer[Tchr][bottle].append('NA')
						maccer[Tchr][bottle].append('NA')
					else:
						maffer[Tchr][bottle].append('NA')
						maffer[Tchr][bottle].append('NA')
						maccer[Tchr][bottle].append('NA')
						maccer[Tchr][bottle].append('NA')
				else:
					maffer[Tchr][bottle].append('NA')
					maccer[Tchr][bottle].append('NA')
							   
			line = line.rstrip('\n')
			line = '\t'.join(line.split())
			line = line.split('\t')
			chr = line[1]
			chr = chr[5:]
			start = str(line[2])
			seq = line[6]
			finish = str(int(line[2]) + int(line[3]))
			maffer[Tchr][bottle].append(seq)
			maccer[Tchr][bottle].append(chr)
			boop = 4
		
print('Volume 1 Complete')


###############################################################
#### Generate Paralog Dictionaries
###############################################################
with open(file2, 'r') as in2:
	for line in in2:
		line = line.rstrip('\n').split('\t')
		P1 = line[0]
		P2 = line[1]
		TY = line[2]
		if P1 not in P1_P2.keys():
			P1_P2[P1] = []
			P1_P2[P1].append(P2)
			P1_TY[P1] = []
			P1_TY[P1].append(TY)
		else:
			P1_P2[P1].append(P2)
			P1_TY[P1].append(TY)
		
		if P2 not in P2_P1.keys():
			P2_P1[P2] = []
			P2_P1[P2].append(P1)
			P2_TY[P2] = []
			P2_TY[P2].append(TY)
		else:
			P2_P1[P2].append(P1)
			P2_TY[P2].append(TY)

print('Volume 2 Complete')


###############################################################
#### Generate Gene Dictionaries
###############################################################
with open(file3, 'r') as in3:
	for line in in3:
		line = line.rstrip('\n').split('\t')
		gene = line[4]
		gst = line[1]
		ged = line[2]
		peak = line[8]
		chr = line[0]
		PK_GN[peak] = gene
		GN_CH[gene] = (chr, gst, ged)

print('Volume 3 Complete')


###############################################################
#### Generate Score Dictionaries
###############################################################
with open(file4, 'r') as in4:
	for line in in4:
		if line[0] == '>':
			line = line.split('\t')
			THRESH = line[2]
			lline = line[1].split('(')
			KEEP = lline[0]

			if KEEP == 'HRE':
				if len(line[0]) == 16:
					KEEP = 'HRE_ST'
				elif len(line[0]) == 21:
					KEEP = 'HRE_HP'

			elif KEEP == 'COUP-TFII':
				if len(line[0]) == 13:
					KEEP = 'COUP-TFII_K5'
				elif len(line[0]) == 9:
					KEEP = 'COUP-TFII_AR'

			elif KEEP == 'c-Myc':
				if len(line[0]) == 9:
					KEEP = 'c-Myc_LN'
				elif len(line[0]) == 11:
					KEEP = 'c-Myc_mE'

			elif KEEP == 'PAX5':
				if len(line[0]) == 15:
					KEEP = 'PAX5_co'
				elif len(line[0]) == 17:
					KEEP = 'PAX5_no'
			elif KEEP == 'THRb':
				if len(line[0]) == 16:
					KEEP = 'THRb_HP'
				elif len(line[0]) == 9:
					KEEP = 'THRb_LV'
			elif KEEP == 'FOXA1':
				if len(line[0]) == 11:
					KEEP = 'FOXA1_MC'
				elif len(line[0]) == 21:
					KEEP = 'FOXA1_LN'


			jumper = 'x'
			numer = 0
			if KEEP in Interests.keys():
				jumper = Interests[KEEP]
				Scorer[KEEP] = {}
				Scorer[KEEP]['Thresh'] = THRESH
		else:
			for item in Interests.keys():
				if Interests[item] == jumper:
					line = line.rstrip('\n').split('\t')
					Scorer[item][numer] = {'A':float(line[0]), 'C':float(line[1]), 'G':float(line[2]), 'T':float(line[3]), 'N':0}
					numer = numer + 1

print('Volume 4 Complete')


###############################################################
#### Generate Alternate Species Site Dictionary
###############################################################

with open(file6, 'r') as mref:
	for line in mref:
		line = line.rstrip('\n').split('\t')
		if line[0] != 'Nfur':
			mastery['Aaus'][line[0]] = line[1]
			mastery['Astr'][line[0]] = line[2]
			mastery['Alim'][line[0]] = line[3]
			mastery['Olat'][line[0]] = line[4]
			mastery['Drer'][line[0]] = line[5]


Assign_Motifs('Aaus', file7, forge)
Assign_Motifs('Astr', file8, forge)
Assign_Motifs('Alim', file9, forge)
Assign_Motifs('Olat', file10, forge)
Assign_Motifs('Drer', file11, forge)

print('Volume 5 Complete')


###############################################################
#### Gather Motif Info
###############################################################
with open(file5, 'r') as in5:

	for line in in5:
		line = line.rstrip('\n').split('\t')
		#Info from Motif file
		chr = line[0]
		start = int(line[1]) - 1
		end = int(line[2]) - 1
		motif = line[3]
		sign = line[7]
		DEID = line[4]
		leng = end - start
		
		if motif == 'HRE':
			if leng == 15:
				motif = 'HRE_ST'
			elif leng == 20:
				motif = 'HRE_HP'
		elif motif == 'COUP-TFII':
			if leng == 12:
				motif = 'COUP-TFII_K5'
			elif leng == 8:
				motif = 'COUP-TFII_AR'
		elif motif == 'c-Myc':
			if leng == 8:
				motif = 'c-Myc_LN'
			elif leng == 10:
				motif = 'c-Myc_mE'
		elif motif == 'PAX5':
			if leng == 14:
				motif = 'PAX5_co'
			elif leng == 16:
				motif = 'PAX5_no'
		elif motif == 'THRb':
			if leng == 15:
				motif = 'THRb_HP'
			elif leng == 8:
				motif = 'THRb_LV'
		elif motif  == 'FOXA1':
			if len(line[0]) == 10:
				motif = 'FOXA1_MC'
			elif len(line[0]) == 20:
				motif = 'FOXA1_LN'


		if motif in Scorer.keys():
			if line[4] == 'NON_DE':
				move = 'ST'
			else:
				move = 'DE'
			ASP = line[6]
		
			#Info on Closest Gene and Chromosome
			if ASP in PK_GN.keys():
				ASG =  PK_GN[ASP]
			else:
				ASG = 'NONE'
			if ASG in GN_CH.keys():
				ASC = GN_CH[ASG]
			else: 
				ASC = 'NA'
			
			#Info on Paralog Grouping
		
			if ASG in P1_TY:
				status = 'P1'
				temp1 = P1_TY[ASG]
				if 'NeoF' in temp1:
					status = status + '/NF'  
				if 'Both up' in temp1:
					status = status + '/BU'  
				if 'Both down' in temp1:
					status = status + '/BD'
				if 'Unclassified' in temp1:
					status = status + '/UC'
			else:
				status = 'NA'
				
			if ASG in P2_TY:
				status2 = 'P2'
				temp2 = P2_TY[ASG]
				if 'NeoF' in temp2:
					status2 = status2 + '/NF'  
				if 'Both up' in temp2:
					status2 = status2 + '/BU'  
				if 'Both down' in temp2:
					status2 = status2 + '/BD'
				if 'Unclassified' in temp2:
					status2 = status2 + '/UC'
			else:
				status2 = 'NA'
		
			if ASG not in P1_TY and ASG not in P2_TY:
				status1 = 'Singleton'
				status2 = 'Singleton'
			
			t_comp = []
			for fish in mastery.keys():
				decomp = 'X'
				if ASP in mastery[fish].keys():
					OF_ASP = mastery[fish][ASP]
					OF_ASP = OF_ASP.split('X')
					OF_ASP = OF_ASP[0]
					if OF_ASP in forge[fish].keys():
						if motif in forge[fish][OF_ASP]:
							decomp = 'D'
				t_comp.append(decomp)
			t_comp = '\t'.join(t_comp)


			###############################################################
			#### Pull Sequences for De Novo Analysis
			###############################################################
		
			case = []
			score1 = 0
			score2 = 0
			score3 = 0
			score4 = 0
			scroe5 = 0
			breaker = []
			CHROM_ARRAY1 = ['NA', 'NA', 'NA', 'NA']
			CHROM_ARRAY2 = ['NA', 'NA', 'NA', 'NA']
			
			if chr in maffer.keys():
				for bottle in maffer[chr].keys():
					if chr == bottle[0]:
						if int(bottle[1]) <= start and int(bottle[2]) >= end:
							dif = start - int(bottle[1])
							CHROM_ARRAY1 = maccer[chr][bottle]
					
							home = maffer[chr][bottle][0]
							home = list(home)
							walker = 0
							bridge = 0
							counter = 0
							hold = []
							while counter < dif:
								if home[walker] != '-':
									counter = counter + 1
								walker = walker + 1
							while counter < dif + leng:
								if home[walker + bridge] != '-':
									counter = counter + 1
								bridge = bridge + 1
					
							for entry in maffer[chr][bottle]:
								if entry == 'NA':
									case.append('NA')
								else:
									if sign == '-':
										temp1 = entry[walker:walker + bridge]
										backer = len(temp1) - 1
										temp2 = []
										while backer > -1:
											temp2.append(COMP[temp1[backer]])
											backer = backer - 1
										case.append(temp2)
									else:
										case.append(entry[walker:walker + bridge])
					
							memory = []
							trace = 0
							fix = 0
							score1 = 1
							score2 = 1
							score3 = 1
							score4 = 1
							score5 = 1
							case[0] = [x.upper() for x in case[0]]
							case[1] = [x.upper() for x in case[1]]
							case[2] = [x.upper() for x in case[2]]
							case[3] = [x.upper() for x in case[3]]
							case[4] = [x.upper() for x in case[4]]

							for letter in case[0]:
								if letter == '-':
									memory.append(trace)
									fix = fix + 1
								else:
									score1 = score1 * Scorer[motif][trace - fix][letter]
								trace = trace + 1
							score1 = -1 * math.log(score1/(1-score1))
							if score1 < float(Scorer[motif]['Thresh']):
								Thresh1 = 'Bound'
							elif score1 >= float(Scorer[motif]['Thresh']):
								Thresh1 = 'Unbound'


							Thresh2 = score_maker(case,memory,score2,Scorer,motif,1)
							Thresh3 = score_maker(case,memory,score3,Scorer,motif,2)
							Thresh4 = score_maker(case,memory,score4,Scorer,motif,3)
							Thresh5 = score_maker(case,memory,score5,Scorer,motif,4)
							
							score2 = Thresh2[1]
							Thresh2 = Thresh2[0]
							score3 = Thresh3[1]
							Thresh3 = Thresh3[0]
							score4 = Thresh4[1]
							Thresh4 = Thresh4[0]
							score5 = Thresh5[1]
							Thresh5 = Thresh5[0]


							breaker.append('thing1') 

					if ASC[0] == bottle[0]:
						if (int(ASC[1]) <= int(bottle[1]) <= int(ASC[2])) or (int(ASC[1]) <= int(bottle[2]) <= int(ASC[2])):
							moving = 0
							for entry in CHROM_ARRAY2:
								if entry == 'NA':
									CHROM_ARRAY2[moving] = maccer[chr][bottle][moving]
								moving = moving + 1
							breaker.append('thing2')
			
					if 'thing1' in breaker and 'thing2' in breaker:
						break
					

			if 'thing1' not in breaker:
				case = ['Unaligned', 'Unaligned', 'Unaligned', 'Unaligned', 'Unaligned',]
				CHROM_ARRAY1 = ['NA', 'NA', 'NA', 'NA',]
				CHROM_ARRAY2 = ['NA', 'NA', 'NA', 'NA',]
				Thresh1 = 'NA'
				score1 = 'NA'
				Thresh2 = 'NA'
				score2 = 'NA'
				Thresh3 = 'NA'
				score3 = 'NA'
				Thresh4 = 'NA'
				score4 = 'NA'
				Thresh5 = 'NA'
				score5 = 'NA'
			
			elif 'thing2' not in breaker:
				CHROM_ARRAY2 = ['NA', 'NA', 'NA', 'NA',]

			step = 0
			transit = []
			while step < 4:
				if CHROM_ARRAY1[step] == 'NA' or CHROM_ARRAY2[step] == 'NA':
					transit.append('NA')
				elif CHROM_ARRAY1[step] == CHROM_ARRAY2[step]:
					transit.append('None')
				else:
					transit.append('Detected')
				step = step + 1
			transit = '\t'.join(transit)


			total1 =  motif + '\t' + ASP + '\t' + move + '\t' + DEID + '\t' + ASG + '\t' + status + '\t' + status2+ '\t' + transit + '\t' + t_comp + '\t'
			total2 = str(score1) + '\t' + str(score2) + '\t' + str(score3) + '\t' + str(score4) + '\t' + str(score5) + '\t'
			total2b = Thresh1 + '\t' + Thresh2 + '\t' + Thresh3 + '\t' + Thresh4 + '\t' + Thresh5 + '\t'
			total3 = ''.join(case[0]) + '\t' + ''.join(case[1]) + '\t' + ''.join(case[2]) + '\t' + ''.join(case[3]) + '\t' + ''.join(case[4]) + '\n'
		
			if motif in final_hold.keys():
				final_hold[motif].append(total1 + total2 + total2b + total3)
			else:
				final_hold[motif] = []
				final_hold[motif].append(total1 + total2 + total2b + total3)
			True_move = True_move + 1
			print(True_move)

	print('Volume 6 Complete')


###############################################################
#### Write all motif files
###############################################################

head = 'Motif\tPeak\tDE_Stat\tUP_DW\tGene\tP1_Stat\tP2_Stat\tAA_T\tAL_T\tOL_T\tDR_T\tAA_AS\tAS_AS\tAL_AS\tOL_AS\tDR_AS\tNF_SS\tAA_SS\tAL_SS\tOL_SS\tDR_SS\tNF_Call\tAA_Call\tAL_Call\tOL_Call\tDR_Call\tNFseq\tAAseq\tALseq\tOLseq\tDRseq\n'

for item in final_hold.keys():
	with open(CURRENT_DiR + 'Data/Raw/' + item + '.txt', 'w') as out1:
		out1.write(head)
		for tbline in final_hold[item]:
			out1.write(tbline)
			

print('Volume 7 Complete')
