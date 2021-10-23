#!/usr/env python3

'''
Example Command-line Initiation:
python3 Multi_Alignment.py ScriptX Species_List.txt

Choose script to run
Script1: Split genome into smaller fasta files and covert them to 2bit format. Generate scaffold size files. Generate whole-genome 2bit file.
Script2: Submit jobs to perform pairwise alignment between species 2bit scaffolds (All-by-All)
Script3: Submit jobs to convert LAV format pairwise alignments to PSL format
Script4: Submit jobs to generate pairwise chain files for all PSL pairwise alignments
Script5: Merge all pairwise chain. Generate pairwise pre-nets and nets. Use merged-chains and nets to generate AXT pairwise alignments. Convert AXT files to MAF format.
Script6: Generate multiple whole-genome alignment using all pairwise whole-genome alignments
'''

import sys
import subprocess
import os
import time
import datetime

SCRIPT = sys.argv[1]
Species = sys.argv[2]
Current_Dir = './'

#####################################
### Define Functions
#####################################

def No_Chromosome_Sorter(input, output):
	
	'''
	Utilized by Script1 to split non-chromosomal-level genome assemblies into equally sized fasta
	files to prevent the need for millions of parallel jobs during SCRIPT2, 3, &4
	'''
	#Define line counting and holding data structures
	counter = 0
	track = 1
	switch = 0
	fil_lin = {}
	
	#Open fasta input file and a renamed fasta file to replace it for whole-genome 2bit conversion
	with open(GENOME + 'whole/' + input, 'r') as inner, open(GENOME + 'whole/' + output, 'w') as edited:
		for line in inner:
			counter = counter + 1
			line = line.rstrip('\n')
			
			#Parse fasta headers
			if line[0] == '>':
				edit = line.rstrip('\n').split('|')
				edit = '>' + edit[3] + '\n'
				edited.write(edit)
				
				#Generate name for new divided fasta files
				if switch == 0:
					switch = 1
					tracker = 'Scaf' + str(track)
					fil_lin[tracker] = []
				
				#Check if the end of the current scaffold extends to or over desired file size
				if counter >= 600000: #Number of total fasta file lines to keep. Will complete current scaffold if over limit
					counter = 0
					track = track + 1
					tracker = 'Scaf' + str(track)
					fil_lin[tracker] = []
					line = line.strip('\n').split('|')
					contig = '>' + line[3]
					fil_lin[tracker].append(contig)
				
				#If file size not yet reached, start on next scaffold
				else:
					tracker =  'Scaf' + str(track)
					line = line.strip('\n').split('|')
					contig = '>' + line[3]
					fil_lin[tracker].append(contig)
			
			#Add line of sequence to dictionary for scaffold
			else:
				edited.write(line + '\n')
				tracker =  'Scaf' + str(track)
				fil_lin[tracker].append(line)
	
	#Once all scaffolds are divided for smaller file use, generate files		
	for unit in fil_lin.keys():
		with open(GENOME + 'split/' + input[:4] + '/' + unit + '.fa', 'w') as temp:
			for entry in fil_lin[unit]:
				temp.write(entry + '\n')
	
	#Rename input fasta with 'old' suffix and rename output fasta as input fasta (for 2bit conversion)
	subprocess.call(['mv', GENOME + 'whole/' + input, GENOME + 'whole/' + input[:4] + '_old.fa'])
	subprocess.call(['mv', GENOME + 'whole/' + output, GENOME + 'whole/' + input])			
		
				
def Submit_Jobs(group1,group2):

	'''
	Utilized by Script2 to submit jobs to cluster to conduct pairwise alignment of scaffolds from each genome comparison
	'''
	
	#Iterate over all 2bit scaffold files for species A and B to generate a pairwise LAV file
	for scaf1 in os.listdir(TWOB + group1[:4] + '/'):
		for scaf2 in os.listdir(TWOB + group2[:4] + '/'):
			if (scaf1[:-5] + '_' + scaf2[:-5] + '.lav') not in os.listdir(LAV + group1[:4] + '_' + group2[:4] + '/'): #Only produce LAV files for comparisons not yet made
				
				#Generate variables for cluster job
				var1 = (TWOB + group1[:4] + '/' + scaf1)
				var2 = (TWOB + group2[:4] + '/' + scaf2)
				var3 = (LAV + group1 + '_' + group2 + '/' + scaf1[:-5] + '_' + scaf2[:-5] + '.lav')
				var4 = (scaf1[:-5] + '_' + scaf2[:-5])
				#Submit job
				subprocess.call(['sbatch', '--job-name=' + var4, '-o', LAV + group1 + '_' + group2 + '/out/out.%j', '-e', LAV + group1 + '_' + group2 + '/err/err.%j', '--export=ONE=' + var1 + ',TWO=' + var2 + ',THREE=' + var3 + ',ZERO=' + Current_Dir, Current_Dir + 'Sub_Scripts/Align_Submission.sh'])


def Submit_Jobs2(groupx):

	'''
	Utilized by Script4 to submit jobs to generate chains for pairwise PSL files
	'''
	
	#Iterate over all PSL files for job submission
	for align in os.listdir(PSL + groupx + '/'):
		if len(align) > 4 and align[:-4] + '.chain' not in os.listdir(CHAINS + groupx + '/'): #Do not submit job for already generated files
			#Generate variables for cluster job
			tmp_name = align[:-4].split('_')
			if len(tmp_name) == 2:
				var1 = TWOB + groupx[:4] + '/' + tmp_name[0] + '.2bit'
				var2 = TWOB + groupx[5:] + '/' + tmp_name[1] + '.2bit'
				var3 = PSL + groupx + '/' + align
				var4 = CHAINS + groupx + '/' + align[:-4] + '.chain'
				var5 = align[:-4]
			elif len(tmp_name) == 3:
				var1 = TWOB + groupx[:4] + '/' + tmp_name[0] + '_' + tmp_name[1] + '.2bit'
				var2 = TWOB + groupx[5:] + '/' + tmp_name[2] + '.2bit'
				var3 = PSL + groupx + '/' + align
				var4 = CHAINS + groupx + '/' + align[:-4] + '.chain'
				var5 = align[:-4]
			#Submit job
			subprocess.call(['sbatch', '--job-name=' + var5, '-o', CHAINS + groupx + '/out/out.%j', '-e', CHAINS + groupx + '/err/err.%j', '--export=THREE=' + var3 + ',ONE=' + var1 + ',TWO=' + var2 + ',FOUR=' + var4 + ',ZERO=' + Current_Dir, Current_Dir + 'Sub_Scripts/Chain_Submission.sh'])
			


#####################################
### Step 0.1 Define Directories
#####################################

GENOME = Current_Dir + 'Data/1_genomes/'
TWOB = Current_Dir + 'Data/2_2bit/'
LAV = Current_Dir + 'Data/3_lav/'
PSL = Current_Dir + 'Data/4_psl/'
CHAINS = Current_Dir + 'Data/5_chain/'
MERCHAINS = Current_Dir + 'Data/6_mchain/'
PNETS = Current_Dir + 'Data/7_pnet/'
NETS = Current_Dir + 'Data/8_net/'
AXTS = Current_Dir + 'Data/9_axt/'
MAFS = Current_Dir + 'Data/10_maf/'
ROAST = Current_Dir + 'Data/11_Final/'
PACK = Current_Dir + '/packages/'
REFER = 'Nfur'


#####################################
### Step 0.2 Make Directories for Data
##################################### 
	
if not os.path.isfile(GENOME + 'sizes/'):
	subprocess.call(['mkdir', GENOME + 'sizes/'])
if not os.path.isfile(GENOME + 'split/'):	
	subprocess.call(['mkdir', GENOME + 'split/'])
if not os.path.isfile(GENOME + 'w2b/'):
	subprocess.call(['mkdir', GENOME + 'w2b/'])
if not os.path.isfile(TWOB):
	subprocess.call(['mkdir', TWOB])
if not os.path.isfile(LAV):
	subprocess.call(['mkdir', LAV])
if not os.path.isfile(PSL):
	subprocess.call(['mkdir', PSL])
if not os.path.isfile(CHAINS):
	subprocess.call(['mkdir', CHAINS])
if not os.path.isfile(MERCHAINS):
	subprocess.call(['mkdir', MERCHAINS])
if not os.path.isfile(PNETS):
	subprocess.call(['mkdir', PNETS])
if not os.path.isfile(NETS):
	subprocess.call(['mkdir', NETS])
if not os.path.isfile(AXTS):
	subprocess.call(['mkdir', AXTS])
if not os.path.isfile(MAFS):
	subprocess.call(['mkdir', MAFS])
if not os.path.isfile(ROAST):
	subprocess.call(['mkdir', ROAST])


#####################################
### Step 0.3 Define Species
#####################################

tracker = 1
#Open Species_List.txt for evaluation
with open(Species, 'r') as inner:
	for line in inner:
		if tracker == 1:
			users = line.rstrip('\n').split(',')
		if tracker == 2:
			line = line.rstrip('\n').split('_')
			TYPE = line[0]
			TREE = line[1]
		tracker = tracker + 1
		
		
#####################################
### Choose SCRIPT to run
#####################################

if SCRIPT == 'Script1':
	#####################################
	### SCRIPT1
	#####################################
	#####################################
	### Step 1 Split Fastas and Reference
	#####################################
	
	#Check if data files match names provide in species_list.txt
	for sample in os.listdir(GENOME + 'whole/'):
		if sample[:4] in users:
			subprocess.call(['mkdir', GENOME + 'split/' + sample[:4]])
			
			#Generate split scaffold/chromosome files based on species parameters
			if sample[:4] == 'Aaus':
				#See function 1
				No_Chromosome_Sorter('Aaus.fa', 'Aaus_new.fa')
			elif sample[:4] == 'Alim':
				#See function 1
				No_Chromosome_Sorter('Alim.fa', 'Alim_new.fa')
			else:
				subprocess.call([PACK + 'faSplit', 'byName', GENOME + 'whole/' + sample, GENOME + 'split/' + sample[:4] + '/'])
			
			#Generate scaffold/chromosome size files	
			with open(GENOME + 'sizes/' + sample[:4] + '.sizes', 'wt') as outer1:
				subprocess.call([PACK + 'faSize', GENOME + 'whole/' + sample, '-detailed',], stdout=outer1)
			
			#Convert split scaffold/chromosome size files from fasta to 2bit format
			subprocess.call([PACK + 'faToTwoBit', GENOME + 'whole/' + sample, GENOME + 'w2b/' + sample[:-3] + '.2bit'])
			
			#Remove all non-chromosomal scaffolds from species with chromosome-level assemblies based on species naming conventions
			Sp_Group = {'Nfur':'NW', 'Aaus':'ZZ', 'Alim':'ZZ', 'Olat':'len', 'Gacu':'scaffold', 'Drer':'len', 'Blan':'xfSc', 'Test':'NW', 'Exam':'NW', 'Samp':'NW'} #Add additional species naming convention here
			subprocess.call(['mkdir', GENOME + 'split/' + sample[:4] + '/temp'])
			for scaf in os.listdir(GENOME + 'split/' + sample[:4] + '/'):
				#eliminate non-chromosomal scaffolds based on length of file name
				if Sp_Group[sample[:4]] == 'len':
					if len(scaf[:-3]) >= 3:
						subprocess.call(['mv', GENOME + 'split/' + sample[:4] + '/' + scaf, GENOME + 'split/' + sample[:4] + '/temp/'])
					else:
						pass
				#eliminate non-chromosomal scaffolds by regular expression criteria
				else:
					if Sp_Group[sample[:4]] in scaf:
						subprocess.call(['mv', GENOME + 'split/' + sample[:4] + '/' + scaf, GENOME + 'split/' + sample[:4] + '/temp/'])
					else:
						pass
			
			#Once separated, combined non-chromosomal scaffolds
			if len(os.listdir(GENOME + 'split/' + sample[:4] + '/temp/')) != 0:
				with open(GENOME + 'split/' + sample[:4] + '/NonChr.fa', 'wt') as joiner:
					for cont in os.listdir(GENOME + 'split/' + sample[:4] + '/temp/'):
						with open(GENOME + 'split/' + sample[:4] + '/temp/' + cont, 'r') as temper:
							for line in temper:
								joiner.write(line)
				#Delete non-chromosomal scaffolds
				subprocess.call(['rm', '-rf', GENOME + 'split/' + sample[:4] + '/NonChr.fa'])
			#Delete temp file
			subprocess.call(['rm', '-rf', GENOME + 'split/' + sample[:4] + '/temp'])
			print(sample[:4] + ' split by chromosome')
		
		else:
			pass

	print('\n\nStep 1.1 Complete\n\n')

	#####################################
	### Step 2 Generate 2bit files
	#####################################		
	
	for folder in os.listdir(GENOME + 'split/'):
		if folder[:4] in users:
			subprocess.call(['mkdir', TWOB + folder[:4]])
			for scaffold in os.listdir(GENOME + 'split/' + folder[:4] + '/'):
				subprocess.call([PACK + 'faToTwoBit', GENOME + 'split/' + folder[:4] + '/' + scaffold, TWOB + folder[:4] + '/' + scaffold[:-3] + '.2bit'])
	
			print('Converted ' + folder[:4] + 'to 2bit format' )
		
	print('\n\nStep 1.2 Complete\n\n')
	


elif SCRIPT == 'Script2':
	#####################################
	### SCRIPT2
	#####################################
	#####################################
	### Step 3 Pairwise Alignment
	#####################################
	
	#Check if data files match names provide in species_list.txt
	for folder1 in os.listdir(TWOB):
		for folder2 in os.listdir(TWOB):
			if (folder1[:4] == REFER) and (folder2[:4] in users) and folder1[:4] != folder2[:4]:
				if (folder1[:4] + '_' + folder2[:4] + '/') not in os.listdir(LAV):
					print('Alinging ' + folder1[:4] + ' and ' + folder2[:4])
					#Generate directories for job submission
					subprocess.call(['mkdir', LAV + folder1[:4] + '_' + folder2[:4]])
					subprocess.call(['mkdir', LAV + folder1[:4] + '_' + folder2[:4] + '/out'])
					subprocess.call(['mkdir', LAV + folder1[:4] + '_' + folder2[:4] + '/err'])
					#See function 2
					Submit_Jobs(folder1,folder2)
			
	print('\n\nStep 2 Complete\n\n')



elif SCRIPT == 'Script3':
	#####################################
	### SCRIPT3
	#####################################
	#####################################
	### Step 4 PSL Conversion
	#####################################	
	
	#Check if data files match names provide in species_list.txt
	for folder in os.listdir(LAV):
		if len(folder) > 3 and folder[:4] in users and folder[5:] in users:
			#Generate directories for job submission
			subprocess.call(['mkdir', PSL + folder])
			subprocess.call(['mkdir', PSL + folder + '/out/'])
			subprocess.call(['mkdir', PSL + folder + '/err/'])
			for align in os.listdir(LAV + folder + '/'):
				if '.lav' in align:
					#Generate variables for job submission
					var1 = LAV + folder + '/' + align
					var2 = PSL + folder + '/' + align[:-4] + '.psl'
					#Submit job
					subprocess.call(['sbatch', '--job-name=' + var1, '-o', PSL + folder + '/out/out.%j', '-e', PSL + folder + '/err/err.%j', '--export=ONE=' + var1 + ',TWO=' + var2 + ',ZERO=' + Current_Dir, Current_Dir + 'Sub_Scripts/PSL_Submission.sh'])
			print('Converted ' + folder + ' to psl format')
			
	print('\n\nStep 3 Complete\n\n')




elif SCRIPT == 'Script4':
	#####################################
	### SCRIPT4
	#####################################
	#####################################
	### Step 5 Chain Generation
	#####################################
	
	#Check if data files match names provide in species_list.txt
	for folder in os.listdir(PSL):
		if folder[:4] in users and folder[5:] in users:
			print('Chaining ' + folder[:4] + ' and ' + folder[5:])
			#Generate directories for job submission
			subprocess.call(['mkdir', CHAINS + folder])
			subprocess.call(['mkdir', CHAINS + folder + '/out'])
			subprocess.call(['mkdir', CHAINS + folder + '/err'])
			#See function 3
			counter = Submit_Jobs2(folder)

	print('\n\nStep 4 Complete\n\n')



elif SCRIPT == 'Script5':
	#####################################
	### SCRIPT5
	#####################################
	#####################################
	### Step 6 Chain Merging and Sorting
	#####################################
	
	#Check if data files match names provide in species_list.txt
	for folder in os.listdir(CHAINS):
		if len(folder) > 3 and folder[:4] in users and folder[5:] in users:
			with open(MERCHAINS + folder + '.chain', 'wt') as outer2:
				crystal = []
				#Avoid counting out/err files from chain generation jobs
				for shard in os.listdir(CHAINS + folder + '/'):
					if 'out' in shard or 'err' in shard:
						pass
					else:
						crystal.append(CHAINS + folder + '/' + shard)
				#Generate list of all chains to be merged
				with open('./chain_list.txt', 'wt') as temper1:
					temper1.write(str('\n'.join(crystal)))
				#Use list of chains to merge them
				subprocess.call([PACK + 'chainMergeSort', '-inputList=./chain_list.txt'], stdout=outer2)
				#Delete list of all chains
				subprocess.call(['rm', './chain_list.txt'])
				print('Merged and Sorted Chains of ' + folder[:4] + ' and ' + folder[5:])
			
	print('\n\nStep 6 Complete\n\n')

	#####################################
	### Step 7 Pre-Netting Chains
	#####################################

	#Check if data files match names provide in species_list.txt
	for chain in os.listdir(MERCHAINS):
		if chain[:4] in users and chain[5:-6] in users:
			#Find chromosomal/scaffold size files for each species
			back1 = GENOME + 'sizes/' + chain[:4] + '.sizes'
			back2 = GENOME + 'sizes/' + chain[5:-6] + '.sizes'
			#Genrate pre-net files
			subprocess.call([PACK + 'chainPreNet', MERCHAINS + chain, back1, back2, PNETS + 'pn_' + chain])
			print('Pre-net Generated for ' + chain[:4] + ' and ' + chain[5:-6])
		

	print('\n\nStep 7 Complete\n\n')

	#####################################
	### Step 8 Net Generation
	#####################################
	
	#Check if data files match names provide in species_list.txt
	for prenet in os.listdir(PNETS):
		if prenet[3:7] in users and prenet[8:-6] in users:
			#Find chromosomal/scaffold size files for each species
			back1 = GENOME + 'sizes/' + prenet[3:7] + '.sizes'
			back2 = GENOME + 'sizes/' + prenet[8:-6] + '.sizes'
			#Generate net files
			subprocess.call([PACK + 'chainNet', PNETS + prenet, '-minSpace=1', back1, back2, NETS + 'temp.net', './dump.net'])
			subprocess.call([PACK + 'netSyntenic', NETS + 'temp.net', NETS + prenet[3:7] + '_' + prenet[8:-6] + '.net'])
			#Remove temporary files
			subprocess.call(['rm', NETS + 'temp.net'])
			subprocess.call(['rm', './dump.net'])
			print('Net Generated for ' + prenet[3:7] + ' and ' + prenet[8:-6])
		
	print('\n\nStep 8 Complete\n\n')

	#####################################
	### Step 9 AXT Conversion and Sorting
	#####################################		
	
	#Check if data files match names provide in species_list.txt
	for net in os.listdir(NETS):
		if net[:4] in users and net[5:-4] in users:
			#Find whole-genome 2bit files to search sequences information from
			ref1 = GENOME + 'w2b/' + net[:4] + '.2bit'
			ref2 = GENOME + 'w2b/' + net[5:-4] + '.2bit'
			#Generate AXT files
			subprocess.call([PACK + 'netToAxt', NETS + net, PNETS + 'pn_' + net[:-4] + '.chain', ref1, ref2, AXTS + 'temp.axt'])
			subprocess.call([PACK + 'axtSort', AXTS + 'temp.axt', AXTS + net[:-4] + '.axt'])
			subprocess.call(['rm', AXTS + 'temp.axt'])
			print('Converted and Sorted ' + net[:4] + ' and ' + net[5:-4] + ' to axt format')
		
	print('\n\nStep 9 Complete\n\n')

	#####################################
	### Step 10 Maf Generation
	#####################################		
	
	#Check if data files match names provide in species_list.txt	
	for axt in os.listdir(AXTS):
		if axt[:4] in users and axt[5:-4] in users:
			#Find chromosomal/scaffold size files for each species
			back1 = GENOME + 'sizes/' + axt[:4] + '.sizes'
			back2 = GENOME + 'sizes/' + axt[5:-4] + '.sizes'
			#Convert AXT files to Maf files
			subprocess.call([PACK + 'axtToMaf', AXTS + axt, back1, back2, MAFS + axt[:4] + '.' + axt[5:-4] + '.sing.maf', '-tPrefix=' + axt[:4] + '.', '-qPrefix=' + axt[5:-4] + '.'])
			print('Maf Generated for' + axt[:4] + ' and ' + axt[5:-4])	
		
	print('\n\nStep 10 Complete\n\n')


	
elif SCRIPT == 'Script6':
	#####################################
	### SCRIPT6
	#####################################
	
	#Genrate multi-alignment from pairwise MAFs
	subprocess.call([PACK + 'roast', '+', 'E=Nfur', TREE, '*.*.maf', 'fish4_WGA.maf' ])

	print('\n\nStep 11 Complete\n\n')

	
