#!/usr/env python3

'''
Example Command Line Input
python3 ATACseq_IVG_Coordinates.py ScriptX

Script1: Generate coordinates list for the center of differential ATACseq peaks in reference species
Script2: Generate coordinates list for the promoter of genes in the reference species
Script3: Generate coordinate link between paralogs in reference species
Script4: Use coordinates to generate anchor base pair for features across species using maf blocks
Script5: Expand anchors in species to width of IGV window and pair peak/gene/paralog entries
'''

##########################################################
### Set environmental variables
##########################################################

import sys
import os
import subprocess
RUNNER = sys.argv[1]
Current_Dir = './'

##########################################################
### Set output file headers and variables
##########################################################

file1 = Current_Dir + 'Data/Nfur_rpkm.txt'
file2 = Current_Dir + 'Data/Nfur_Peak_All.bed'
file3 = Current_Dir + 'Data/Nfur.gtf'
file4 = Current_Dir + 'Data/Nfur_gene.bed'
file5 = Current_Dir + 'Data/Promoter_List.txt'
file6 = Current_Dir + 'Data/Partner_List.txt'
file7 = Current_Dir + 'Data/fish4_WGA.maf'
file8 = Current_Dir + 'Data/Peak_Master_0.0.txt'
file9 = Current_Dir + 'Data/Manual_Cans.txt'

hold = {} 
SIGN = {}
FINAL = {}
SIGN = {}
maf_lib = {}
pek_con = {}

Nfur = {}
Aaus = {}
Alim = {}
Olat = {}
Drer = {}

Nfur2 = {}
Aaus2 = {}
Alim2 = {}
Olat2 = {}
Drer2 = {}


block = 10000000 #Maximum number of blocks in maf file
filer = ['Promoter_Pre_Anchors.txt', 'DE_Pre_Anchors.txt', 'Partner_Pre_Anchors.txt']



###############################################################
#### Define Functions
###############################################################	   
			
def Species_Anchor(a_peak, a_abr, a_spec):

	'''
	Used by Script 4 to relate anchor coordinates between species
	'''
	
	#Evaluate if peak exist in data and if conserved anchor directly on peak in other species
	if a_peak in pek_con.keys() and a_abr in pek_con[a_peak]:
		if a_spec in maf_lib[step].keys():
			true_spec = str(maf_lib[step][a_spec][0]) + '\t' + str(maf_lib[step][a_spec][1])
		
	#If the peak is not conserved find exact basepair match to reference species anchor
		else:
			full = 0
			temp = 0
			sign = 1
			#Traverse along the length of other species maf block
			while full == 0:
				temp = temp + 1
				sign = sign*(-1)
				guess = step + (temp*sign)
				
				if guess > 0 and guess < len(maf_lib.keys()):
					#If desired coordinate of bp falls on a gap in alignment alternate moving left and right to find closest base
					if a_spec in maf_lib[guess].keys():
						full = 1
						#If base found is to the left of gap
						if sign < 0:
							tac = maf_lib[step]['Nfur'][1] - maf_lib[guess]['Nfur'][1]
							true_spec = str(maf_lib[guess][a_spec][0]) + '\t' + str(maf_lib[guess][a_spec][1] + tac)
					
						#if base found is to the right of the gap
						elif sign > 0:
							tac = maf_lib[guess]['Nfur'][1] - maf_lib[step]['Nfur'][1]
							true_spec = str(maf_lib[guess][a_spec][0]) + '\t' + str(maf_lib[guess][a_spec][1] - tac)
				else:
					full = 1
					true_spec = 'NA\tNA'
	
	#If peak coordinate does not exist, write a null entry
	else:
		true_spec = 'NA\tNA'
	return(true_spec)
	
	
def build_entry(a_P1, a_P2, a_lib, a_Alt, a_alt_lib):

	'''
	Used by Script 5 to to expand single base anchors for ATACseq peaks to a coordinate window for IGV track files
	'''
	
	#generate entry for peak 1 if it exists
	if a_P1 in a_lib.keys() and a_lib[a_P1][0] != 'NA':
		#make sure that the anchor is at least 2000bp away from a chromosome/scaffold end
		if a_lib[a_P1][1] < 2000:
			step1 = a_lib[a_P1][0] + ':' + str(0) + '-' + str(a_lib[a_P1][1] + 2000)
		else:
			step1 = a_lib[a_P1][0] + ':' + str(a_lib[a_P1][1] - 2000) + '-' + str(a_lib[a_P1][1] + 2000)
			
	#If peak does not exist find the promoter entry anchor for the gene
	elif a_Alt in a_alt_lib.keys() and a_alt_lib[a_Alt][0] != 'NA':
		#make sure that the anchor is at least 2000bp away from a chromosome/scaffold end
		if a_alt_lib[a_Alt][1] < 2000:
			step1 = a_alt_lib[a_Alt][0] + ':' + str(0) + '-' + str(a_alt_lib[a_Alt][1] + 2000)
		else:
			step1 = a_alt_lib[a_Alt][0] + ':' + str(a_alt_lib[a_Alt][1] - 2000) + '-' + str(a_alt_lib[a_Alt][1] + 2000)
	
	#If peak and promoter are missing, write a null entry
	else:
		step1 = 'NA:NA-NA'
	
	#generate entry for peak 2 if it exists
	if a_P2 in a_lib.keys() and a_lib[a_P2][0] != 'NA':
		#make sure that the anchor is at least 2000bp away from a chromosome/scaffold end
		if a_lib[a_P2][1] < 2000:
			step2 = a_lib[a_P2][0] + ':' + str(0) + '-' + str(a_lib[a_P2][1] + 2000)
		else:
			step2 = a_lib[a_P2][0] + ':' + str(a_lib[a_P2][1] - 2000) + '-' + str(a_lib[a_P2][1] + 2000)
			
	#If peak does not exist find the promoter entry anchor for the gene
	elif a_Alt in a_alt_lib.keys() and a_alt_lib[a_Alt][0] != 'NA':
		#make sure that the anchor is at least 2000bp away from a chromosome/scaffold end
		if a_alt_lib[a_Alt][1] < 2000:
			step2 = a_alt_lib[a_Alt][0] + ':' + str(0) + '-' + str(a_alt_lib[a_Alt][1] + 2000)
		else:
			step2 = a_alt_lib[a_Alt][0] + ':' + str(a_alt_lib[a_Alt][1] - 2000) + '-' + str(a_alt_lib[a_Alt][1] + 2000)
	
	#If peak and promoter are missing, write a null entry
	else:
		step2 = 'NA:NA-NA'
	
	#Put together peak 1 and peak 2 entry
	step3 =  step1 + ' ' + step2 + '\n'
	return(step3)



if RUNNER == 'Script1':

	###############################################################
	#### Script 1
	###############################################################	
	###############################################################
	#### Generate Differential Pre-Coordinate Anchors
	###############################################################	

	with open(file1, 'r') as in1:
		#Read in Nfur ATACseq peaks
		for line in in1:
			temp = line.rstrip('\n').split('\t')
			use = temp[3]
			hold[use] = line

	with open(file1, 'r') as in2, open(Current_Dir + 'Data/DE_Pre_Anchors.txt', 'w') as out1:
		#Keep peak data for all differential peaks
		for line in in2:
			line = line.rstrip('\n').split('\t')
			peak = line[3]
			if peak in hold.keys() and peak != 'Peak':
				out1.write(hold[peak])



if RUNNER == 'Script2':

	###############################################################
	#### Script 2
	###############################################################	          
	###############################################################
	#### Generate Promoter Pre-Coordinate Anchors
	###############################################################	           

	with open(file3, 'r') as in1:
		#Read in GTF file for NFUR
		for line in in1:
			line = line.rstrip('\n').split('\t')
			#Identify the promoter of each gene
			if (line[2] == 'exon' or 'CDS') and 'gene_name' in line[8]:
				temp = line[8].split(';')
				for entry in temp:
					if 'gene_name' in entry:
						temp2 = entry.split('"')
						GENE = temp2[1]
				WAY = line[6]
				SIGN[GENE] = WAY
			
	with open(file4, 'r') as in2:
		for line in in2:
			line = line.rstrip('\n').split('\t')
			gene = line[4]
			#Determine gene orientation and start/stop coordinates
			if SIGN[gene] == '-':
				FINAL[gene] = (str(line[0]), str(line[2]))
			elif SIGN[gene] == '+':
				FINAL[gene] = (str(line[0]), str(line[1]))
	
	#Read in list of need promoters for paralogs	
	with open(file5, 'r') as in3, open(Current_Dir + 'Data/Promoter_Pre_Anchors.txt', 'w') as out1:
		for line in in3:
			line = line.rstrip('\n').split('\t')
			line = line[0]
			#Output all data on promoters
			out1.write(FINAL[line][0] + '\t' + FINAL[line][1] + '\t' + FINAL[line][1] + '\t' + line + '\n')
					  
					  
					  
if RUNNER == 'Script3':
			  
	###############################################################
	#### Script 3
	###############################################################						  
	###############################################################
	#### Generate pre-coordinate anchors paired ATACseq peaks
	###############################################################					  				  
	
	with open(file1, 'r') as in1:
		#Read in Nfur ATACseq peaks
		for line in in1:
			temp = line.rstrip('\n').split('\t')
			use = temp[3]
			hold[use] = temp
	
	with open(file6, 'r') as in1:
		for line in in1:
			#Read in ATACseq peak data
			line = line.rstrip('\n').split('\t')
			peak = line[0]
			if peak in hold.keys():
				SIGN[peak] = hold[peak][0] + '\t' + hold[peak][1] + '\t' + hold[peak][2] + '\t' + hold[peak][3] + '\n' 	
	
	with open(Current_Dir + 'Data/Partner_Pre_Anchors.txt', 'w') as out1:
		for entry in SIGN.keys():
			line = SIGN[entry]
			out1.write(line)
		
		
		
if RUNNER == 'Script4':
			  
	###############################################################
	#### Script 4
	###############################################################						  
	###############################################################
	#### Graph species coordinates using maf alignment blocks
	###############################################################

	with open(file7, 'r') as in1:
		current = 0
		for line in in1:
			#Sort Nfur block entry data for storage
			if line[0] != '#' and 'Nfur' in line:
				current = current + 1
				line = line.rstrip('\n')
				line = '\t'.join(line.split())
				line = line.split('\t')
				chr = line[1]
				chr = chr[5:]
				start = str(line[2])
				finish = str(int(line[2]) + int(line[3]))
				maf_lib[current] = {}
				maf_lib[current]['Nfur'] = (chr, int(start), int(finish))
			
			#Sort Aaus block data for storage
			elif line[0] != '#' and 'Aaus' in line:
				line = line.rstrip('\n')
				line = '\t'.join(line.split())
				line = line.split('\t')
				chr = line[1]
				chr = chr[5:]
				start = str(line[2])
				finish = str(int(line[2]) + int(line[3]))
				maf_lib[current]['Aaus'] = (chr, int(start), int(finish))
			
			#Sort Alim block data for storage
			elif line[0] != '#' and 'Alim' in line:				
				line = line.rstrip('\n')
				line = '\t'.join(line.split())
				line = line.split('\t')
				chr = line[1]
				chr = chr[5:]
				start = str(line[2])
				finish = str(int(line[2]) + int(line[3]))
				maf_lib[current]['Alim'] = (chr, int(start), int(finish))
			
			#Sort Olat block data for storage
			elif line[0] != '#' and 'Olat' in line:
				line = line.rstrip('\n')
				line = '\t'.join(line.split())
				line = line.split('\t')
				chr = line[1]
				chr = chr[5:]
				start = str(line[2])
				finish = str(int(line[2]) + int(line[3]))
				maf_lib[current]['Olat'] = (chr, int(start), int(finish))

			#Sort Drer block data for storage
			elif line[0] != '#' and 'Drer' in line:			   
				line = line.rstrip('\n')
				line = '\t'.join(line.split())
				line = line.split('\t')
				chr = line[1]
				chr = chr[5:]
				start = str(line[2])
				finish = str(int(line[2]) + int(line[3]))
				maf_lib[current]['Drer'] = (chr, int(start), int(finish))
	
	print('Step 1 Complete')

	###############################################################
	#### Generate Conservation Dictionaries
	###############################################################

	with open(file8, 'r') as in3:
		for line in in3:
			line = line.rstrip('\n').split('\t')
			pek_con[line[0]] = []
			if line[1] != 'NA' or line[2] != 'NA':
				pek_con[line[0]].append('Aa')
			if line[2] != 'NA':
				pek_con[line[0]].append('Al')
			if line[3] != 'NA':
				pek_con[line[0]].append('Ol')
			if line[4] != 'NA':
				pek_con[line[0]].append('Dr')

	print('Step 2 Complete')

	###############################################################
	#### Output Peak Anchors
	###############################################################	   

	for sets in filer:
		#Run each of the generated pre-coordinate anchor lists
		with open(Current_Dir + 'Data/' + sets, 'r') as in2, open(Current_Dir + 'Data/' + sets[:-15] + 'Final_Anchors.txt', 'w') as out1:
			#Generate header for output file
			header = 'Peak\tNfurC\tNfurA\tAausC\tAausA\tAlimC\tAlimA\tOlatC\tOlatA\tDrerC\tDrerA\n'
			out1.write(header)
			place = 0
	
			for line in in2:
				#Parse data for pre-coordinate anchor lists
				line = line.rstrip('\n').split('\t')
				if line[0] != 'chromosome':
					peak = line[3]
					chrom = line[0]
					starter = int(line[1])
					diff = int((int(line[2]) - int(line[1]))/2)
					cent = int(starter) + int(diff)
					fin = 0
					step = 0
					#Traverse the maf block containing peak until the center is found
					while fin == 0:
						step = step + 1
						#If your reach the end of the block without finding the anchor generate null entry for peak
						#if step >= block:
						if step > len(maf_lib.keys()):
							fin = 1
							out1.write(peak + '\t' + chrom + '\t' + str(cent) + '\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n')
							place = place + 1
							print('Page ' + str(place) + ' Complete')
						#Check if you have found the right chromosome and are centered on the peak
						elif maf_lib[step]['Nfur'][0] == chrom:
							if maf_lib[step]['Nfur'][1] < cent < maf_lib[step]['Nfur'][2]:
								fin = 1
								true_Nf = peak + '\t' + str(maf_lib[step]['Nfur'][0]) + '\t' + str(maf_lib[step]['Nfur'][1])
								true_Aa = Species_Anchor(peak, 'Aa', 'Aaus')
								true_Al = Species_Anchor(peak, 'Al', 'Alim')
								true_Ol = Species_Anchor(peak, 'Ol', 'Olat')
								true_Dr = Species_Anchor(peak, 'Dr', 'Drer')
								#Write the correct anchor coordinates for each species for each pre-coordinate anchors
								out1.write(true_Nf + '\t' + true_Aa + '\t' + true_Al + '\t' + true_Ol + '\t' + true_Dr + '\n')
								place = place + 1
								print('Page ' + str(place) + ' Complete')
				
					
		print('Step 3 Complete')
	
	

if RUNNER == 'Script5':
			  
	###############################################################
	#### Script 5
	###############################################################						  
	###############################################################
	#### Expand anchor coordinate to IGV window width
	###############################################################	

	#Save anchors of differential ATACseq peaks
	with open(Current_Dir + 'Data/DE_Final_Anchors.txt', 'r') as in1:
		for line in in1:
			line = line.rstrip('\n').split('\t')
			if line[0] != 'Peak':
				Peak = line[0]
				if line[2] != 'NA':
					Nfur[Peak] = (line[1], int(line[2]))
				else:
					Nfur[Peak] = (line[1], 'NA')
			
				if line[4] != 'NA':
					Aaus[Peak] = (line[3], int(line[4]))
				else:
					Aaus[Peak] = (line[3], 'NA')
			
				if line[6] != 'NA':
					Alim[Peak] = (line[5], int(line[6]))
				else:
					Alim[Peak] = (line[5], 'NA')
			
				if line[8] != 'NA':
					Olat[Peak] = (line[7], int(line[8]))
				else:
					Olat[Peak] = (line[7], 'NA')
			
				if line[10] != 'NA':
					Drer[Peak] = (line[9], int(line[10]))
				else:
					Drer[Peak] = (line[9], 'NA')
	
	#Save anchors for promoters of Nfur genes			
	with open(Current_Dir + 'Data/Promoter_Final_Anchors.txt', 'r') as in13:
		for line in in13:
			line = line.rstrip('\n').split('\t')
			if line[0] != 'Peak':
				Peak = line[0]
				if line[2] != 'NA':
					Nfur[Peak] = (line[1], int(line[2]))
				else:
					Nfur[Peak] = (line[1], 'NA')
			
				if line[4] != 'NA':
					Aaus[Peak] = (line[3], int(line[4]))
				else:
					Aaus[Peak] = (line[3], 'NA')
			
				if line[6] != 'NA':
					Alim[Peak] = (line[5], int(line[6]))
				else:
					Alim[Peak] = (line[5], 'NA')
			
				if line[8] != 'NA':
					Olat[Peak] = (line[7], int(line[8]))
				else:
					Olat[Peak] = (line[7], 'NA')
			
				if line[10] != 'NA':
					Drer[Peak] = (line[9], int(line[10]))
				else:
					Drer[Peak] = (line[9], 'NA')
	
	#Save the anchors for paralogs of Nfur genes		
	with open(Current_Dir + 'Data/Partner_Final_Anchors.txt', 'r') as in12:
		for line in in12:
			line = line.rstrip('\n').split('\t')
			if line[0] != 'Peak':
				Peak = line[0]
				if line[2] != 'NA':
					Nfur2[Peak] = (line[1], int(line[2]))
				else:
					Nfur2[Peak] = (line[1], 'NA')
			
				if line[4] != 'NA':
					Aaus2[Peak] = (line[3], int(line[4]))
				else:
					Aaus2[Peak] = (line[3], 'NA')
			
				if line[6] != 'NA':
					Alim2[Peak] = (line[5], int(line[6]))
				else:
					Alim2[Peak] = (line[5], 'NA')
			
				if line[8] != 'NA':
					Olat2[Peak] = (line[7], int(line[8]))
				else:
					Olat2[Peak] = (line[7], 'NA')
			
				if line[10] != 'NA':
					Drer2[Peak] = (line[9], int(line[10]))
				else:
					Drer2[Peak] = (line[9], 'NA')
				
	#Generate output file to hold IGV window entry data
	with open(file9, 'r') as in2, open(Current_Dir + 'Data/Manual_IVGC.txt', 'w') as out1:
		for line in in2:
			line = line.rstrip('\n').split('\t')
			P1 = line[0]
			P2 = line[1]
			Alt = line[2]
			
			#Genrate the output set for each species (See Functions above)
			true1 = build_entry(P1, P2, Nfur, Alt, Nfur2)
			true2 = build_entry(P1, P2, Aaus, Alt, Aaus2)
			true3 = build_entry(P1, P2, Alim, Alt, Alim2)
			true4 = build_entry(P1, P2, Olat, Alt, Olat2)
			true5 = build_entry(P1, P2, Drer, Alt, Drer2)
			
			#Write output if proper pair exists for peak
			if line[1] != 'NA':
				P2 = line[1]
			else:
				P2 = line[2]
			#Name each entry based on the saved IGV track that will be manually generated
			true0 = P1 + '_' + P2 + '.svg' + '\n'
			out1.write(true0 + true1 + true2 + true3 + true4 + true5 + '\n\n\n')
	
