#!/usr/env python3

'''
Example Command Line Execution
python3 Peak_Cons_Bed_Generator_1.py
'''

###############################################################
#### Set Environmental Variables
###############################################################

CURRENT_DIR = './'
ALLER = ['Singletons_Comb_DE.bed','Singletons_Comb.bed','Singletons_Dia_DE.bed', 'Singletons_Dia.bed', 'Singletons_Dev_DE.bed', 'Singletons_Dev.bed', 'lowCV_GENES.bed','lowCV_GENES_DE.bed','lowCV_PEAKS_Final.bed', 'Dev_Relaxed.bed', 'Dev_Strict.bed', 'Dia_Relaxed.bed']
USER = ['0.0', '0.25', '0.5']

###############################################################
#### Define Functions
###############################################################
def bed_dic(a_file,a_dic):
	'''
	This function reads a given bedfile into a dictionary containing each be entry
	'''
	with open(a_file, 'r') as intemp:
		for line in intemp:
			line = line.rstrip('\n').split('\t')
			if line[0] != 'chromosome':
				a_dic[line[3]] = [line[0], int(line[1]), int(line[2])]

###############################################################		  
def para_list(a_file, a_list1, a_list2):
	'''
	This function reads in a list of paralog pairs and saves them within a python list
	'''
	with open(a_file, 'r') as intemp:
		for line in intemp:
			line = line.rstrip('\n').split('\t')
			if line[2] == 'NeoF':
				a_list1.append(line[0])
				a_list2.append(line[1])
					
###############################################################
def cons_dic(a_file, a_dic1, a_dic2, a_dic3, a_dic4, a_dic5):
	'''
	This function reads in peak conservation data and simplifies each entry for storage in a python dictionary
	'''
	with open(a_file, 'r') as intemp:
		for line in intemp:
			line = line.rstrip('\n').split('\t')
			if 'peak' in line[1]:
				temp = line[1].split('X')
				temp = temp[0]
				a_dic1[line[0]] = temp
			elif 'None' in line[1]:
				a_dic1[line[0]] = 'None'
			else:
				a_dic1[line[0]] = 'NA'
				
			if 'peak' in line[2]:
				temp = line[2].split('X')
				temp = temp[0]
				a_dic2[line[0]] = temp
			elif 'None' in line[2]:
				a_dic2[line[0]] = 'None'
			else:
				a_dic2[line[0]] = 'NA'
				
			if 'peak' in line[3]:
				temp = line[3].split('X')
				temp = temp[0]
				a_dic3[line[0]] = temp
			elif 'None' in line[3]:
				a_dic3[line[0]] = 'None' 
			else:
				a_dic3[line[0]] = 'NA'
				
			if 'peak' in line[4]:
				temp = line[4].split('X')
				temp = temp[0]
				a_dic4[line[0]] = temp
			elif 'None' in line[4]:
				a_dic4[line[0]] = 'None'
			else:
				a_dic4[line[0]] = 'NA'
				
			if 'peak' in line[5]:
				temp = line[5].split('X')
				temp = temp[0]
				a_dic5[line[0]] = temp
			elif 'None' in line[5]:
				a_dic5[line[0]] = 'None' 
			else:
				a_dic5[line[0]] = 'NA'

###############################################################
def Define_Achors(a_chr, a_start, a_end):
	'''
	This function defines a central anchor base between conserved peaks between species
	'''
	flood = 0
	if a_chr in Nfur.keys():
		for entry in Nfur[a_chr]:
			if Nfur[a_chr][entry][0] <= a_start <= Nfur[a_chr][entry][1]:
				BA = entry
				flood = 1
				Boff = Nfur[a_chr][entry][0] - a_start
				stepper = 0
				truer = 0
				while truer < Boff:
						base = Nfur[a_chr][entry][2][stepper]
						if base != '-':
							truer = truer + 1
							stepper =  stepper + 1
						else:
							stepper = stepper + 1
				TBoff = stepper
							
			if flood == 1 and (Nfur[a_chr][entry][0] <= a_end <= Nfur[a_chr][entry][1]):
				AA = entry
				flood = 2
				Aoff = Nfur[a_chr][entry][0] - a_end
				stepper = 0
				truer = 0
				while truer < Aoff:
						base = Nfur[a_chr][entry][2][stepper]
						if base != '-':
							truer = truer + 1
							stepper =  stepper + 1
						else:
							stepper = stepper + 1
				TAoff = stepper
				
			if flood == 2:
				break
		
		if flood ==2:
			return([BA, TBoff, AA, TAoff])
		else:
			return(1,0,1,0)


###############################################################
def Grab_Coord(a_chr, a_list, a_spec):
	'''
	This function runs through a sequence entry align across species and produces a bed entry of the 
	aligned chromosomal coordinates in all other species of interest.
	''' 
	if chr in a_spec.keys():
		follow = 0
		if a_list[0] in a_spec[chr].keys():
			holder = 0
			jumper = 0
			while holder < a_list[1]:
				if a_spec[chr][a_list[0]][3][holder] != '-':
					holder = holder + 1
					jumper = jumper + 1
				else:
					holder = holder + 1
			spec_chr1 = a_spec[chr][a_list[0]][0]
			spec_start = a_spec[chr][a_list[0]][1]
			spec_start = spec_start + jumper
			follow = 1
			
		else:
			inner = a_list[0]
			while inner < a_list[2]:
				inner = inner + 1
				if inner in a_spec[chr].keys():
					spec_chr1 = a_spec[chr][inner][0]
					spec_start = a_spec[chr][inner][1]
					follow = 1
					break
		
		if follow == 0:
			spec_chr1 = 'wrong'
			spec_start = 'wrong2'
		
		follow = 0
		if a_list[2] in a_spec[chr].keys():
			holder = 0
			jumper = 0
			while holder < a_list[3]:
				if a_spec[chr][a_list[2]][3][holder] != '-':
					holder = holder + 1
					jumper = jumper + 1
				else:
					holder = holder + 1
					
			spec_chr2 = a_spec[chr][a_list[2]][0]
			spec_end = a_spec[chr][a_list[2]][1]
			spec_end = spec_end + jumper
			follow = 1
			
		else:
			outer = a_list[2]
			while outer > a_list[0]:
				outer = outer - 1
				if outer in a_spec[chr].keys():
					spec_chr2 = a_spec[chr][outer][0]
					spec_end = a_spec[chr][outer][2]
					follow = 1
		
		if follow == 0:
			spec_chr2 = 'wrong3'
			spec_end = 'wrong4'
		
		if spec_chr1 == spec_chr2:
			if spec_start <= spec_end:
				return([spec_chr1, spec_start, spec_end])
			else:
				return([spec_chr1, spec_end, spec_start])
		else:
			return(['chimeric', 'chimeric', 'chimeric'])

###############################################################
#### Iterate over all possible file combinations for naming
###############################################################

for THINGER in ALLER:
	for ITEM in USER:
		SPECIAL = THINGER[:-4] + '_' + str(ITEM) + '_DIAPAUSE.txt'
		geneid1 = 4 #normal = 18, trunc = 4
		geneid2 = 4 #normal = 19, trunc = 4

		file1 = CURRENT_DIR + 'Data/fish4_WGA.maf'
		file2 = CURRENT_DIR + 'Data/' + THINGER
		file3 = CURRENT_DIR + 'Data/Nfur_rpkm.txt'
		file4 = CURRENT_DIR + 'Data/Aaus_rpkm.txt'
		file5 = CURRENT_DIR + 'Data/Astr_rpkm.txt'
		file6 = CURRENT_DIR + 'Data/Alim_rpkm.txt'
		file7 = CURRENT_DIR + 'Data/Olat_rpkm.txt'
		file8 = CURRENT_DIR + 'Data/Drer_rpkm.txt'
		file9 = CURRENT_DIR + 'Data/Peak_Master_' + ITEM + '.txt'
		file10 = ''#should remain empty
		file11 = CURRENT_DIR + 'Data/Paralog_List.txt'
		file12 = CURRENT_DIR + 'Data/' + THINGER[:-4] + '_' + str(ITEM) + '/' + SPECIAL

		Nfur = {}
		Aaus = {}
		Alim = {}
		Olat = {}
		Drer = {}
		nfur_p = {}
		aaus_p = {}
		astr_p = {}
		alim_p = {}
		olat_p = {}
		drer_p = {}
		nfur_c = {}
		aaus_c = {}
		astr_c = {}
		alim_c = {}
		olat_c = {}
		drer_c = {}
		neo1 = []
		neo2 = []				 
		 
		###############################################################
		#### Generate Maf Dictionaries
		###############################################################
		with open(file1, 'r') as in1:
			current = 0
			for line in in1:
				if line[0] != '#' and 'Nfur' in line:
					line = line.rstrip('\n')
					line = '\t'.join(line.split())
					line = line.split('\t')
					chr = line[1]
					chr = chr[5:]
					start = int(line[2])
					finish = int(line[2]) + int(line[3])
					seq = str(line[6])
					if chr in Nfur.keys():
						current = current + 1
						Nfur[chr][current] = []
						Nfur[chr][current].append(start)
						Nfur[chr][current].append(finish)
						Nfur[chr][current].append(seq)
					else:
						current = 1
						Nfur[chr] = {}
						Nfur[chr][current] = []
						Nfur[chr][current].append(start)
						Nfur[chr][current].append(finish)
						Nfur[chr][current].append(seq)
					

				elif line[0] != '#' and 'Aaus' in line:
					line = line.rstrip('\n')
					line = '\t'.join(line.split())
					line = line.split('\t')
					temp = line[1]
					temp = temp[5:]
					start = int(line[2])
					finish = int(line[2]) + int(line[3])
					seq = str(line[6])
					if chr in Aaus.keys():
						Aaus[chr][current] = []
						Aaus[chr][current].append(temp)
						Aaus[chr][current].append(start)
						Aaus[chr][current].append(finish)
						Aaus[chr][current].append(seq)
					else: 
						Aaus[chr] = {}
						Aaus[chr][current] = []
						Aaus[chr][current].append(temp)
						Aaus[chr][current].append(start)
						Aaus[chr][current].append(finish)
						Aaus[chr][current].append(seq)
		
		
				elif line[0] != '#' and 'Alim' in line:
					line = line.rstrip('\n')
					line = '\t'.join(line.split())
					line = line.split('\t')
					temp = line[1]
					temp = temp[5:]
					start = int(line[2])
					finish = int(line[2]) + int(line[3])
					seq = str(line[6])
			
					if chr in Alim.keys():
						Alim[chr][current] = []
						Alim[chr][current].append(temp)
						Alim[chr][current].append(start)
						Alim[chr][current].append(finish)
						Alim[chr][current].append(seq)
			
					else: 
						Alim[chr] = {}
						Alim[chr][current] = []
						Alim[chr][current].append(temp)
						Alim[chr][current].append(start)
						Alim[chr][current].append(finish)
						Alim[chr][current].append(seq)
		
				
				elif line[0] != '#' and 'Olat' in line:
					line = line.rstrip('\n')
					line = '\t'.join(line.split())
					line = line.split('\t')
					temp = line[1]
					temp = temp[5:]
					start = int(line[2])
					finish = int(line[2]) + int(line[3])
					seq = str(line[6])
			
					if chr in Olat.keys():
						Olat[chr][current] = []
						Olat[chr][current].append(temp)
						Olat[chr][current].append(start)
						Olat[chr][current].append(finish)
						Olat[chr][current].append(seq)
			
					else: 
						Olat[chr] = {}
						Olat[chr][current] = []
						Olat[chr][current].append(temp)
						Olat[chr][current].append(start)
						Olat[chr][current].append(finish)
						Olat[chr][current].append(seq)
		
			
				elif line[0] != '#' and 'Drer' in line:
					line = line.rstrip('\n')
					line = '\t'.join(line.split())
					line = line.split('\t')
					temp = line[1]
					temp = temp[5:]
					start = int(line[2])
					finish = int(line[2]) + int(line[3])
					seq = str(line[6])
			
					if chr in Drer.keys():
						Drer[chr][current] = []
						Drer[chr][current].append(temp)
						Drer[chr][current].append(start)
						Drer[chr][current].append(finish)
						Drer[chr][current].append(seq)
			
					else: 
						Drer[chr] = {}
						Drer[chr][current] = []
						Drer[chr][current].append(temp)
						Drer[chr][current].append(start)
						Drer[chr][current].append(finish)
						Drer[chr][current].append(seq)
				
		print(THINGER +' Volume 1 Complete')
				
		###############################################################
		#### Generate Meta Data
		###############################################################	

		bed_dic(file3, nfur_p)
		bed_dic(file4, aaus_p)
		bed_dic(file5, astr_p)
		bed_dic(file6, alim_p)
		bed_dic(file7, olat_p)
		bed_dic(file8, drer_p)
		cons_dic(file9, aaus_c,astr_c,alim_c,olat_c,drer_c)			
		para_list(file11, neo1, neo2)

		print(THINGER + ' Volume 2 Complete')

		###############################################################
		#### Generate Output File
		###############################################################	

		with open(file2, 'r') as in3, open(file12, 'w') as out1:
			head1 = 'Nf_chr\tNF_st\tNf_ed\tNf_pk\tNf_ne'
			head2 = 'Aa_chr\tAa_st\tAa_ed\tAa_mt'
			head3 = 'As_chr\tAs_st\tAs_ed\tAs_mt'
			head4 = 'Al_chr\tAl_st\tAl_ed\tAl_mt'
			head5 = 'Ol_chr\tOl_st\tOl_ed\tOl_mt'
			head6 = 'Dr_chr\tDr_st\tDr_ed\tDr_mt'
			header = head1 + '\t' + head2 + '\t' + head3 + '\t' + head4 + '\t' + head5 + '\t' + head6 + '\n'
			out1.write(header)
			checker = 0
			track = 0
	
			for line in in3:
				liner = ''
				line = line.rstrip('\n').split('\t')
				if line[geneid1] in neo1 or line[geneid2] in neo1:
					status = 'NeoF1'
				elif line[geneid1] in neo2 or line[geneid2] in neo2:
					status = 'NeoF2'
				else:
					status = 'Other'
			
				liner = liner + line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + line[3] + '\t' + status + '\t'
		
				chrom = line[0]
				sanchor = int(line[1])
				eandcor = int(line[2])
				peaker = line[3]
		
				if chrom[1] != 'W':
					if aaus_c[peaker] != 'NA' and aaus_c[peaker] != 'None':
						temper = aaus_p[aaus_c[peaker]]
						liner = liner + temper[0] + '\t' + str(temper[1]) + '\t' + str(temper[2]) + '\t' + aaus_c[peaker] + '\t'
					elif aaus_c[peaker] != 'NA':
						Anchor = Define_Achors(line[0], int(line[1]), int(line[2]))
						Cord = Grab_Coord(line[0], Anchor, Aaus)
						liner = liner + Cord[0] + '\t' + str(Cord[1]) + '\t' + str(Cord[2]) + '\t' + 'Aligned' + '\t'
						checker = 1
					else:
						liner = liner + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'Unaligned' + '\t'
		
			
					if astr_c[peaker] != 'NA' and astr_c[peaker] != 'None':
						temper = astr_p[astr_c[peaker]]
						liner = liner + temper[0] + '\t' + str(temper[1]) + '\t' + str(temper[2]) + '\t' + astr_c[peaker] + '\t'
					elif checker == 1 and astr_c[peaker] != 'NA':
						Cord = Grab_Coord(line[0], Anchor, Aaus)
						liner = liner + Cord[0] + '\t' + str(Cord[1]) + '\t' + str(Cord[2]) + '\t' + 'Aligned' + '\t'   
					elif astr_c[peaker] != 'NA':
						Anchor = Define_Achors(line[0], int(line[1]), int(line[2]))
						Cord = Grab_Coord(line[0], Anchor, Aaus)
						liner = liner + Cord[0] + '\t' + str(Cord[1]) + '\t' + str(Cord[2]) + '\t' + 'Aligned' + '\t'
						checker = 1
					else:
						liner = liner + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'Unaligned' + '\t'
		
			
					if alim_c[peaker] != 'NA' and alim_c[peaker] != 'None':
						temper = alim_p[alim_c[peaker]]
						liner = liner + temper[0] + '\t' + str(temper[1]) + '\t' + str(temper[2]) + '\t' + alim_c[peaker] + '\t'
					elif checker == 1 and alim_c[peaker] != 'NA':
						Cord = Grab_Coord(line[0], Anchor, Alim) 
						liner = liner + Cord[0] + '\t' + str(Cord[1]) + '\t' + str(Cord[2]) + '\t' + 'Aligned' + '\t'  
					elif alim_c[peaker] != 'NA':
						Anchor = Define_Achors(line[0], int(line[1]), int(line[2]))
						Cord = Grab_Coord(line[0], Anchor, Alim)
						liner = liner + Cord[0] + '\t' + str(Cord[1]) + '\t' + str(Cord[2]) + '\t' + 'Aligned' + '\t'
						checker = 1
					else:
						liner = liner + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'Unaligned' + '\t'
		
		
					if olat_c[peaker] != 'NA' and olat_c[peaker] != 'None':
						temper = olat_p[olat_c[peaker]]
						liner = liner + temper[0] + '\t' + str(temper[1]) + '\t' + str(temper[2]) + '\t' + olat_c[peaker] + '\t'
					elif checker == 1 and olat_c[peaker] != 'NA':
						Cord = Grab_Coord(line[0], Anchor, Olat)
						liner = liner + Cord[0] + '\t' + str(Cord[1]) + '\t' + str(Cord[2]) + '\t' + 'Aligned' + '\t'	
					elif olat_c[peaker] != 'NA':
						Anchor = Define_Achors(line[0], int(line[1]), int(line[2]))
						Cord = Grab_Coord(line[0], Anchor, Olat)
						liner = liner + Cord[0] + '\t' + str(Cord[1]) + '\t' + str(Cord[2]) + '\t' + 'Aligned' + '\t'
						checker = 1
					else:
						liner = liner + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'Unaligned' + '\t'
		
		
					if drer_c[peaker] != 'NA' and drer_c[peaker] != 'None':
						temper = drer_p[drer_c[peaker]]
						liner = liner + temper[0] + '\t' + str(temper[1]) + '\t' + str(temper[2]) + '\t' + drer_c[peaker] + '\n'
					elif checker == 1 and drer_c[peaker] != 'NA':
						Cord = Grab_Coord(line[0], Anchor, Drer)
						liner = liner + Cord[0] + '\t' + str(Cord[1]) + '\t' + str(Cord[2]) + '\t' + 'Aligned' + '\n'	
					elif drer_c[peaker] != 'NA':
						Anchor = Define_Achors(line[0], int(line[1]), int(line[2]))
						Cord = Grab_Coord(line[0], Anchor, Drer)
						liner = liner + Cord[0] + '\t' + str(Cord[1]) + '\t' + str(Cord[2]) + '\t' + 'Aligned' + '\n'
						checker = 1
					else:
						liner = liner + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'Unaligned' + '\n'
		
					out1.write(liner)
					checker = 0
					track = track + 1
					print('finished page ' + str(track))
			
		print(THINGER + ' Volume 3 Complete')
		
			