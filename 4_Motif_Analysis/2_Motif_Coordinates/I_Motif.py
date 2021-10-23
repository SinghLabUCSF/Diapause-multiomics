#!/usr/env python3

'''
Example command line execution
python3 I_Motif.py
'''

#################################################################################		   
#### Define Constants
################################################################################# 
import subprocess
CURRENT_DIR = './'
PACK = ''

file01 = CURRENT_DIR + 'Data/Para_UP.txt'
file02 = CURRENT_DIR + 'Data/Para_DW.txt'
file03 = CURRENT_DIR + 'Data/TE_DE_UP.bed' 
file04 = CURRENT_DIR + 'Data/TE_DE_DW.bed' 
file1 = CURRENT_DIR + 'Data/AllMotifs_ConsensusPeaks_Nfur.txt'
file2 = CURRENT_DIR + 'Data/Interest_Motif.bed'
file3 = CURRENT_DIR + 'Data/Nfur_TE.bed'
file4 = CURRENT_DIR + 'Data/Interest_Motif_sort.bed'
file5 = CURRENT_DIR + 'Data/Motif_TE.bed'
file6 = CURRENT_DIR + 'Data/Motif_TE.txt'

mark1 = 'Neo_P1'
mark2 = 'Neo_P2'
mark3 = 'Both_U'
mark4 = 'Both_D'
mark5 = 'Unchar'
multi = '),'
counter = 0

IPEAK1 = []
IPEAK2 = []
GPEAK1 = []
GPEAK2 = []
GPEAK3 = []
GPEAK4 = []
GPEAK5 = []

IND_2_MOT = {}
Interests = ['AMYB', 'AP-1', 'Ap4', 'Arnt:Ahr', 'Ascl1', 'Atf1', 'Atf2', 'Atf3', 'Atf4', 'Atf7', 'Atoh1', 'Bach1', 'Bach2', 'Bapx1', 'Barx1', 'BATF', 'BHLHA15', 'bHLHE40', 'bHLHE41', 'BMAL1', 'BMYB', 'BORIS', 'Brn1', 'Brn2', 'bZIP:IRF', 'c-Jun-CRE', 'c-Myc', 'c-Myc', 'CArG', 'Cdx2', 'CDX4', 'CEBP', 'Chop', 'CLOCK', 'COUP-TFII', 'COUP-TFII', 'CRE', 'CRX', 'CTCF-SatelliteElement', 'CTCF', 'CUX1', 'Cux2', 'Dlx3', 'Duxbl', 'E-box', 'E2F3', 'E2F6', 'EAR2', 'Egr1', 'Egr2', 'Erra', 'ERRg', 'Esrrb', 'Fosl2', 'Fox:Ebox', 'Foxa2', 'Foxa3', 'FoxD3', 'Foxf1', 'FOXK1', 'FOXK2', 'FoxL2', 'FOXM1', 'Foxo1', 'Foxo3', 'FOXP1', 'Fra1', 'Fra2', 'FXR', 'GATA:SCL', 'Gfi1b', 'GFX', 'GFY-Staf', 'GSC', 'Hand2', 'HIF-1b', 'HIF2a', 'HLF', 'HNF4a', 'HNF6', 'Hnf6b', 'HOXA1', 'Hoxa10', 'Hoxa11', 'Hoxa13', 'HOXA2', 'Hoxa9', 'HOXB13', 'Hoxb4', 'Hoxc9', 'Hoxd10', 'Hoxd11', 'Hoxd12', 'Hoxd13', 'HRE', 'HRE', 'IRF1', 'IRF2', 'IRF8', 'Isl1', 'ISRE', 'Jun-AP1', 'JunB', 'JunD', 'KLF10', 'KLF14', 'KLF3', 'Klf4', 'KLF5', 'KLF6', 'Klf9', 'LEF1', 'Lhx1', 'Lhx2', 'Lhx3', 'LXH9', 'LXRE', 'MafA', 'MafB', 'Max', 'Maz', 'Meis1', 'MITF', 'MNT', 'MYB', 'MyoD', 'n-Myc', 'Nanog', 'NeuroD1', 'NeuroG2', 'NF-E2', 'NF1', 'NFE2L2', 'NFIL3', 'NFY', 'Nkx3.1', 'Nkx6.1', 'NPAS', 'NPAS2', 'Nrf2', 'Nur77', 'OCT:OCT', 'Oct11', 'Oct2', 'OCT4-SOX2-TCF-NANOG', 'Oct4:Sox17', 'Oct4', 'Oct6', 'Olig2', 'Otx2', 'PAX3:FKHR-fusion', 'PAX5', 'PAX5', 'PAX6', 'Pax8', 'PBX2', 'Pbx3', 'Pdx1', 'Phox2a', 'Pit1', 'Pitx1:Ebox', 'Pitx1', 'Pknox1', 'PPARa', 'PPARE', 'PR', 'PRDM1', 'PRDM14', 'Prop1', 'PU.1-IRF', 'PU.1', 'RAR:RXR', 'RARa', 'RBPJ:Ebox', 'Rbpj1', 'REST-NRSF', 'Reverb', 'RFX', 'Rfx1', 'Rfx2', 'Rfx5', 'RORa', 'RORg', 'RORgt', 'RUNX2', 'RXR', 'SCL', 'Six1', 'Six2', 'Six4', 'Smad2', 'Smad3', 'Smad4', 'Sox10', 'Sox15', 'Sox17', 'Sox2', 'Sox3', 'Sox4', 'Sox6', 'Sox9', 'Sp1', 'Sp2', 'Sp5', 'SpiB', 'Srebp1a', 'Srebp2', 'STAT1', 'TATA-Box', 'Tbet', 'Tbx20', 'Tcf12', 'Tcf21', 'Tcf3', 'TCF4', 'Tcf7', 'TEAD', 'TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 'TFE3', 'Tgif1', 'Tgif2', 'THRa', 'THRb', 'TR4', 'Twist2', 'Unknown-ESC-element', 'Unknown', 'USF1', 'Usf2', 'VDR', 'WT1', 'X-box', 'YY1', 'Zac1', 'ZBTB18', 'Zfp281', 'Zic', 'Zic3', 'ZKSCAN1', 'ZNF143|STAF', 'ZNF189', 'Znf263', 'ZNF382', 'ZNF41', 'ZNF467', 'ZNF519', 'ZNF669', 'ZNF675', 'ZNF692', 'ZNF7', 'ZNF768', 'ZSCAN22']


#################################################################################		   
#### Define Functions
################################################################################# 
'''
This function take a single peak and check whether it is differentially expressing in 
diapause, development, or neither. It then evaluates whether the closest gene to the given
peak is a paralog, it neo-functionalized status, or lack there of.
'''

def determine_setting(peak):
	setting = ['x','x']
	if peak in IPEAK1:
		setting[0] = 'Diapause'
	elif peak in IPEAK2:
		setting[0] = 'Development'
	else:
		setting[0] = 'NON_DE'
				
	if peak in GPEAK1:
		setting[1] = mark1
	elif peak in GPEAK2:
		setting[1] = mark2
	elif peak in GPEAK3:
		setting[1] = mark3	
	elif peak in GPEAK4:
		setting[1] = mark4
	elif peak in GPEAK5:
		setting[1] = mark5
	else:
		setting[1] = 'Other'
	
	return(setting)

#################################################################################		   
#### Generate Peak TYPE and STATE lists
################################################################################# 

with open(file01, 'r') as in01, open(file02, 'r') as in02:
	for line in in01:
		line = line.rstrip('\n').split('\t')
		stepper = 0 
		type = line[0]
		for item in line[2:1679]:
			stepper = stepper + 1
			if stepper%3 == 1:
				temp = item
				if item != 'NA':
					if type == 'NeoF':
						GPEAK1.append(temp)
					elif type == 'Both up':
						GPEAK3.append(temp)
					elif type == 'Both down':
						GPEAK4.append(temp)
					elif type == 'Uncategorized':
						GPEAK5.append(temp)
						
		stepper = 0
		for item in line[1680:]:
			stepper = stepper + 1
			if stepper%3 == 1:
				temp = item
				if item != 'NA':
					if type == 'NeoF':
						GPEAK2.append(temp)
					elif type == 'Both up':
						GPEAK3.append(temp)
					elif type == 'Both down':
						GPEAK4.append(temp)
					elif type == 'Uncategorized':
						GPEAK5.append(temp)				  
						
	for line in in02:
		line = line.rstrip('\n').split('\t')
		stepper = 0 
		type = line[0]
		for item in line[2:1679]:
			stepper = stepper + 1
			if stepper%3 == 1:
				temp = item
				if item != 'NA':
					if type == 'NeoF':
						GPEAK1.append(temp)
					elif type == 'Both up':
						GPEAK3.append(temp)
					elif type == 'Both down':
						GPEAK4.append(temp)
					elif type == 'Uncategorized':
						GPEAK5.append(temp)
						
		stepper = 0
		for item in line[1680:]:
			stepper = stepper + 1
			if stepper%3 == 1:
				temp = item
				if item != 'NA':
					if type == 'NeoF':
						GPEAK2.append(temp)
					elif type == 'Both up':
						GPEAK3.append(temp)
					elif type == 'Both down':
						GPEAK4.append(temp)
					elif type == 'Uncategorized':
						GPEAK5.append(temp)	  


with open(file03, 'r') as in03, open(file04, 'r') as in04:
	for line in in03:
		line = line.rstrip('\n').split('\t')
		upper = line[3]
		IPEAK1.append(upper)
			   
	for line in in04:
		line = line.rstrip('\n').split('\t')
		downer = line[3]
		IPEAK2.append(downer)  
		
print('Volume 1 Finished')

			  
#################################################################################		   
#### Generate Motif Bed File
################################################################################# 

with open(file1, 'r') as in1, open(file2, 'w') as out1:
	for line in in1:
		line = line.rstrip('\n').split('\t')
		IDER = line[0]
		CHR = line[1]
		pst = line[2]
		pen = line[3]
		
		if 'PeakID' in IDER:
			for item in line[21:]:
				temp = item.split('(')
				temp = temp[0]
				IND_2_MOT[counter + 21] = temp
				counter = counter + 1
				
		else:
			held = determine_setting(IDER)
			state = held[0]
			set = held[1]
			for item in line[21:]:
				if item != '' and (IND_2_MOT[line.index(item)] in Interests):
					site = IND_2_MOT[line.index(item)]
					if multi in item:
						item = item.split('),')
						for mot in item:
							mot = mot.split('(')
							tstart = int(pst) + int(mot[0])
							mini = mot[1].split(',')
							tend = int(tstart) + int(len(mini[0]))
							sign = mini[1]
							out1.write(CHR + '\t' + str(tstart) + '\t' + str(tend) + '\t' + site + '\t' + state + '\t' + set + '\t' + IDER + '\t' + sign + '\n')
							
					else:
						item = item.split('(')
						tstart = int(pst) + int(item[0])
						mini = item[1].split(',')
						tend = int(tstart) + int(len(mini[0]))
						sign = mini[1]
						out1.write(CHR + '\t' + str(tstart) + '\t' + str(tend) + '\t' + site + '\t' + state + '\t' + set + '\t' + IDER + '\t' + sign + '\n')

print('Volume 2 Finished')


#################################################################################		   
#### Sort Bed File
#################################################################################
		   
with open(file4, 'wt') as sort:
	subprocess.call([PACK + 'bedtools', 'sort', '-i', file2], stdout = sort)
print('Volume 3 Finished')


#################################################################################		   
#### Merge Bed Files
#################################################################################
	
with open(file5, 'wt') as temp:
	subprocess.call([PACK + 'bedtools', 'intersect', '-a', file4, '-b', file3, '-wao'], stdout = temp)
print('Volume 4 Finished')	

#################################################################################		   
#### Generate Stats Document
#################################################################################
	
with open(file5, 'r') as temp, open(file6, 'w') as final:
	for line in temp:
		line = line.rstrip('\n').split('\t')
		peak = line[6]
		motif = line[3]
		if line[7] != '.':
			te = line[11]
		else:
			te = 'NA'
		state = line[4]
		set = line[5]
		final.write(peak + '\t' + motif + '\t' + te + '\t' + state + '\t' + set + '\n')
		
print('Volume 5 Finished')		

				