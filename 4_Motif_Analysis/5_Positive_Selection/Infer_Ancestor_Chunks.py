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

home = CURRENT_DIR + 'Data/'
#PAML/baseml-formatted instruction inputs
control = '        noisy = 0   * 0,1,2,3: how much rubbish on the screen\n      verbose = 0   * 1: detailed output, 0: concise output\n      runmode = 0   * 0: user tree;  1: semi-automatic;  2: automatic\n                    * 3: StepwiseAddition; (4,5):PerturbationNNI \n\n        model = 0   * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85\n                    * 5:T92, 6:TN93, 7:REV, 8:UNREST, 9:REVu; 10:UNRESTu\n\n        Mgene = 0   * 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff\n\n*        ndata = 100\n        clock = 0   * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis\n    fix_kappa = 0   * 0: estimate kappa; 1: fix kappa at value below\n        kappa = 5  * initial or fixed kappa\n\n    fix_alpha = 0   * 0: estimate alpha; 1: fix alpha at value below\n        alpha = 0.5   * initial or fixed alpha, 0:infinity (constant rate)\n       Malpha = 0   * 1: different alpha\'s for genes, 0: one alpha\n        ncatG = 5   * # of categories in the dG, AdG, or nparK models of rates\n        nparK = 0   * rate-class models. 1:rK, 2:rK&fK, 3:rK&MK(1/K), 4:rK&MK \n\n        nhomo = 0   * 0 & 1: homogeneous, 2: kappa for branches, 3: N1, 4: N2\n        getSE = 0   * 0: don\'t want them, 1: want S.E.s of estimates\n RateAncestor = 1   * (0,1,2): rates (alpha>0) or ancestral states\n\n   Small_Diff = 7e-6\n    cleandata = 0  * remove sites with ambiguity data (1:yes, 0:no)?\n*        icode = 0  * (with RateAncestor=1. try "GC" in data,model=4,Mgene=4)\n*  fix_blength = -1  * 0: ignore, -1: random, 1: initial, 2: fixed\n       method = 0  * Optimization method 0: simultaneous; 1: one branch a time\n'

out1 = home + 'Nfur_Chunk.fa'
out2 = home + 'Teleost_Ancestor_Chunk.fa'
out3 = home + 'Pre_Kill_Ancestor_Chunk.fa'
out4 = home + 'Killi_Ancestor_Chunk.fa'
out5 = home + 'African_Ancestor_Chunk.fa'

Nfur_s = []
Telo_s = []
Pkil_s = []
Kill_s = []
Afri_s = []


##########################################################
### Generate baseml control file
##########################################################
for file in os.listdir(home + 'Chunk_Seqs/'):
	seqs = '      seqfile = ' + home + 'Chunk_Seqs/' + file + ' \n'
	tree = '     treefile = ' + home + 'fish.trees' + '\n'
	out = '\n      outfile = ' + home + file + '       * main result file\n'
	temper = seqs + tree + out + control
	with open(home + 'temp_baseml.ctl', 'w') as temp:
		temp.write(temper)
		
	
	##########################################################
	### Run baseml
	##########################################################
	subprocess.call([PACK + 'baseml', home + 'temp_baseml.ctl'])
	subprocess.call(['rm', home + 'temp_baseml.ctl'])
	subprocess.call(['rm', home + file])
	
	##########################################################
	### Extract Ancestral Sequence Information
	##########################################################
	with open('./rst', 'r') as final:
		counter = 0
		for line in final:
			line = line.rstrip('\n')
			
			if counter == 4:
				counter = counter + 1
				line = line.rstrip(' ')
				line = line[18:]
				line = line.split(' ')
				line = ''.join(line)
				line = line.split('-')
				line = ''.join(line)
				
				iter = '>' + file[:-3] + '\n' + line + '\n'
				Nfur_s.append(iter)
				
			
			elif 12 >= counter >= 9:
				line = line.rstrip(' ')
				line = line[18:]
				line = line.split(' ')
				line = ''.join(line)
				line = line.split('-')
				line = ''.join(line)
				
				iter = '>' + file[:-3] + '\n' + line + '\n'
				if counter == 9:
					Telo_s.append(iter)
				elif counter == 10:
					Pkil_s.append(iter)
				elif counter == 11:
					Kill_s.append(iter)
				elif counter == 12:
					Afri_s.append(iter)
				counter = counter + 1
					
			elif counter >= 13:
				break
				
			elif counter > 0:
				counter = counter + 1
				
			if line[0:4] == 'List':
				counter = 1
				

##########################################################
### Merge all ancestral sequence information
##########################################################
with open(out1, 'w') as o1, open(out2, 'w') as o2, open(out3, 'w') as o3, open(out4, 'w') as o4, open(out5, 'w') as o5:
	for entry in Nfur_s:
		o1.write(entry)
	for entry in Telo_s:
		o2.write(entry)
	for entry in Pkil_s:
		o3.write(entry)
	for entry in Kill_s:
		o4.write(entry)
	for entry in Afri_s:
		o5.write(entry)
	
	