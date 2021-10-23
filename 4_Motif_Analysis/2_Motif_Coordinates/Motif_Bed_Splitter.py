#/usr/env python3

'''
Example Command Line Entry

python3 Motif_Bed_Splitter.py
'''

#######################################################################################
##### Define Constants
#######################################################################################
CURRENT_DIR = './'

home = CURRENT_DIR
input = 'Data/Interest_Motif_sort.bed'

Motifs = {}

#######################################################################################
##### Read in combined motif data from master file
#######################################################################################
with open(home + input, 'r') as in1:
	for line in in1:
		line = line.rstrip('\n').split('\t')
		chrom = line[0]
		start = line[1]
		end = line[2]
		mot = line[3]
		holder = chrom + '\t' + start + '\t' + end + '\t' + mot + '\n'
		
		if mot in Motifs.keys():
			Motifs[mot].append(holder)
		else:
			Motifs[mot] = []
			Motifs[mot].append(holder)

#######################################################################################
##### Generate individual bed files for each motif
#######################################################################################
for motter in Motifs.keys():
	with open(home + 'Data/Interest_Beds/' + motter + '.bed', 'w') as temp:
		for entry in Motifs[motter]:
			temp.write(entry)
			