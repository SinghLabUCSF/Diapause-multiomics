#!/usr/env python3

'''
Example Command Line Execution
python3 Peak_Cons_Bed_Generator_2.py
'''

###############################################################
#### Set Environmental Variables
###############################################################
import sys
import subprocess

CURRENT_DIR = './'
TYPER = ['N', 'A']
ALLER = ['Dev_Relaxed_0.0/', 'Dev_Relaxed_0.25/', 'Dev_Relaxed_0.5/', 'Dev_Strict_0.0/', 'Dev_Strict_0.25/', 'Dev_Strict_0.5/', 'Dia_Relaxed_0.0/', 'Dia_Relaxed_0.25/', 'Dia_Relaxed_0.5/', 'LowCV_GENES_0.0/', 'LowCV_GENES_0.25/', 'LowCV_GENES_0.5/', 'LowCV_GENES_DE_0.0/', 'LowCV_GENES_DE_0.25/', 'LowCV_GENES_DE_0.5/', 'LowCV_PEAKS_Final_0.0/', 'LowCV_PEAKS_Final_0.25/', 'LowCV_PEAKS_Final_0.5/', 'Singletons_Comb_0.0/', 'Singletons_Comb_0.25/', 'Singletons_Comb_0.5/',  'Singletons_Comb_DE_0.0/', 'Singletons_Comb_DE_0.25/', 'Singletons_Comb_DE_0.5/', 'Singletons_Dev_0.0/', 'Singletons_Dev_0.25/', 'Singletons_Dev_0.5/',  'Singletons_Dev_DE_0.0/', 'Singletons_Dev_DE_0.25/', 'Singletons_Dev_DE_0.5/', 'Singletons_Dia_0.0/', 'Singletons_Dia_0.25/', 'Singletons_Dia_0.5/',  'Singletons_Dia_DE_0.0/', 'Singletons_Dia_DE_0.25/', 'Singletons_Dia_DE_0.5/']

for THINGER in ALLER:
	home = CURRENT_DIR + 'Data/'+ THINGER
	file1 = home + THINGER[:-1] + '_DIAPAUSE.txt'


	###############################################################
	#### Generate All Output files for given input
	###############################################################
	for decid in TYPER:
		if decid == 'A':
			plate = '_All.bed'
		elif decid == 'N':
			plate = '_Neo_Combine.bed'
		
		file2 = home + 'Nfur' + plate
		file3 = home + 'Nfur_Neo1.bed'
		file4 = home + 'Nfur_Neo2.bed'
		file5 = home + 'Z_trash1.txt'
		file6 = home + 'Z_trash2.txt'
		file7 = home + 'Z_trash3.txt'
		file8 = home + 'Z_trash4.txt'
		file9 = home + 'Z_trash5.txt'
		file10 = home + 'Z_trash6.txt'
		
		subprocess.call(['mkdir', home + 'Align'])
		subprocess.call(['mkdir', home + 'Either'])
		subprocess.call(['mkdir', home + 'Peak'])
		file11 = home + 'Align/Aaus' + plate
		file12 = home + 'Align/Aaus_Neo1.bed'
		file13 = home + 'Align/Aaus_Neo2.bed'
		file14 = home + 'Either/Aaus' + plate
		file15 = home + 'Either/Aaus_Neo1.bed'
		file16 = home + 'Either/Aaus_Neo2.bed'
		file17 = home + 'Peak/Aaus' + plate
		file18 = home + 'Peak/Aaus_Neo1.bed'
		file19 = home + 'Peak/Aaus_Neo2.bed'

		file20 = home + 'Align/Astr' + plate
		file21 = home + 'Align/Astr_Neo1.bed'
		file22 = home + 'Align/Astr_Neo2.bed'
		file23 = home + 'Either/Astr' + plate
		file24 = home + 'Either/Astr_Neo1.bed'
		file25 = home + 'Either/Astr_Neo2.bed'
		file26 = home + 'Peak/Astr' + plate
		file27 = home + 'Peak/Astr_Neo1.bed'
		file28 = home + 'Peak/Astr_Neo2.bed'

		file29 = home + 'Align/Alim' + plate
		file30 = home + 'Align/Alim_Neo1.bed'
		file31 = home + 'Align/Alim_Neo2.bed'
		file32 = home + 'Either/Alim' + plate
		file33 = home + 'Either/Alim_Neo1.bed'
		file34 = home + 'Either/Alim_Neo2.bed'
		file35 = home + 'Peak/Alim' + plate
		file36 = home + 'Peak/Alim_Neo1.bed'
		file37 = home + 'Peak/Alim_Neo2.bed'

		file38 = home + 'Align/Olat' + plate
		file39 = home + 'Align/Olat_Neo1.bed'
		file40 = home + 'Align/Olat_Neo2.bed'
		file41 = home + 'Either/Olat' + plate
		file42 = home + 'Either/Olat_Neo1.bed'
		file43 = home + 'Either/Olat_Neo2.bed'
		file44 = home + 'Peak/Olat' + plate
		file45 = home + 'Peak/Olat_Neo1.bed'
		file46 = home + 'Peak/Olat_Neo2.bed'

		file47 = home + 'Align/Drer' + plate
		file48 = home + 'Align/Drer_Neo1.bed'
		file49 = home + 'Align/Drer_Neo2.bed'
		file50 = home + 'Either/Drer' + plate
		file51 = home + 'Either/Drer_Neo1.bed'
		file52 = home + 'Either/Drer_Neo2.bed'
		file53 = home + 'Peak/Drer' + plate
		file54 = home + 'Peak/Drer_Neo1.bed'
		file55 = home + 'Peak/Drer_Neo2.bed'

		file56 = home + 'Align/AaB_Nfur' + plate
		file57 = home + 'Align/AaB_Nfur_Neo1.bed'
		file58 = home + 'Align/AaB_Nfur_Neo2.bed'
		file59 = home + 'Either/AaB_Nfur' + plate
		file60 = home + 'Either/AaB_Nfur_Neo1.bed'
		file61 = home + 'Either/AaB_Nfur_Neo2.bed'
		file62 = home + 'Peak/AaB_Nfur' + plate
		file63 = home + 'Peak/AaB_Nfur_Neo1.bed'
		file64 = home + 'Peak/AaB_Nfur_Neo2.bed'

		file65 = home + 'Align/AsB_Nfur' + plate
		file66 = home + 'Align/AsB_Nfur_Neo1.bed'
		file67 = home + 'Align/AsB_Nfur_Neo2.bed'
		file68 = home + 'Either/AsB_Nfur' + plate
		file69 = home + 'Either/AsB_Nfur_Neo1.bed'
		file70 = home + 'Either/AsB_Nfur_Neo2.bed'
		file71 = home + 'Peak/AsB_Nfur' + plate
		file72 = home + 'Peak/AsB_Nfur_Neo1.bed'
		file73 = home + 'Peak/AsB_Nfur_Neo2.bed'

		file74 = home + 'Align/AlB_Nfur' + plate
		file75 = home + 'Align/AlB_Nfur_Neo1.bed'
		file76 = home + 'Align/AlB_Nfur_Neo2.bed'
		file77 = home + 'Either/AlB_Nfur' + plate
		file78 = home + 'Either/AlB_Nfur_Neo1.bed'
		file79 = home + 'Either/AlB_Nfur_Neo2.bed'
		file80 = home + 'Peak/AlB_Nfur' + plate
		file81 = home + 'Peak/AlB_Nfur_Neo1.bed'
		file82 = home + 'Peak/AlB_Nfur_Neo2.bed'

		file83 = home + 'Align/OlB_Nfur' + plate
		file84 = home + 'Align/OlB_Nfur_Neo1.bed'
		file85 = home + 'Align/OlB_Nfur_Neo2.bed'
		file86 = home + 'Either/OlB_Nfur' + plate
		file87 = home + 'Either/OlB_Nfur_Neo1.bed'
		file88 = home + 'Either/OlB_Nfur_Neo2.bed'
		file89 = home + 'Peak/OlB_Nfur' + plate
		file90 = home + 'Peak/OlB_Nfur_Neo1.bed'
		file91 = home + 'Peak/OlB_Nfur_Neo2.bed'

		file92 = home + 'Align/DrB_Nfur' + plate
		file93 = home + 'Align/DrB_Nfur_Neo1.bed'
		file94 = home + 'Align/DrB_Nfur_Neo2.bed'
		file95 = home + 'Either/DrB_Nfur' + plate
		file96 = home + 'Either/DrB_Nfur_Neo1.bed'
		file97 = home + 'Either/DrB_Nfur_Neo2.bed'
		file98 = home + 'Peak/DrB_Nfur' + plate
		file99 = home + 'Peak/DrB_Nfur_Neo1.bed'
		file100 = home + 'Peak/DrB_Nfur_Neo2.bed'


		###############################################################
		#### Create Master Dictionary and list to hold all data
		###############################################################
		Filer = [file2, file3, file4, file5, file6, file7, file8, file9, file10, file11, file12, file13, file14, file15, file16, file17, file18, file19, file20, file21, file22, file23, file24, file25, file26, file27, file28, file29, file30, file31, file32, file33, file34, file35, file36, file37, file38, file39, file40, file41, file42, file43, file44, file45, file46, file47, file48, file49, file50, file51, file52, file53, file54, file55, file56, file57, file58, file59, file60, file61, file62, file63, file64, file65, file66, file67, file68, file69, file70, file71, file72, file73, file74, file75, file76, file77, file78, file79, file80, file81, file82, file83, file84, file85, file86, file87, file88, file89, file90, file91, file92, file93, file94, file95, file96, file97, file98, file99, file100]
		Datar = {2:[], 3:[], 4:[], 11:[], 12:[], 13:[], 17:[], 18:[], 19:[], 20:[], 21:[], 22:[], 26:[], 27:[], 28:[], 29:[], 30:[], 31:[], 35:[], 36:[], 37:[], 38:[], 39:[], 40:[], 44:[], 45:[], 46:[], 47:[], 48:[], 49:[], 53:[], 54:[], 55:[], 56:[], 57:[], 58:[], 62:[], 63:[], 64:[], 65:[], 66:[], 67:[], 71:[], 72:[], 73:[], 74:[], 75:[], 76:[], 80:[], 81:[], 82:[], 83:[], 84:[], 85:[], 89:[], 90:[], 91:[], 92:[], 93:[], 94:[], 98:[], 99:[], 100:[]} 

		###############################################################
		#### Parse input file
		###############################################################
		with open(file1, 'r') as inner:
			for line in inner:
				line = line.rstrip('\n').split('\t')
				Nfur = line[0] + '\t' + line[1] + '\t' +line[2] + '\t' + line[3] + '\n'
				Aaus = line[5] + '\t' + line[6] + '\t' +line[7] + '\t' + line[8] + '\n'
				Astr = line[9] + '\t' + line[10] + '\t' +line[11] + '\t' + line[12] + '\n'
				Alim = line[13] + '\t' + line[14] + '\t' +line[15] + '\t' + line[16] + '\n'
				Olat = line[17] + '\t' + line[18] + '\t' +line[19] + '\t' + line[20] + '\n'
				Drer = line[21] + '\t' + line[22] + '\t' +line[23] + '\t' + line[24] + '\n'
				typer = {8:Aaus,12:Astr,16:Alim,20:Olat,24:Drer}
		
		
				trip = 0
				back = 49
				spec = 4
				agar = 43
				agor = -2
			
				###############################################################
				#### Sort NEOF1 sites
				###############################################################
				if line[4] == 'NeoF1':
					for item in typer.keys():
						back = back + 5
						spec = spec + 5
						agar = agar + 5
						agor = agor + 5
						if 'peak' in line[item]:
							if trip == 0: 
								Datar[2].append(Nfur) #All
								Datar[3].append(Nfur) #Neo1
								trip = 1
							Datar[back + item].append(Nfur) #B_All
							Datar[back + item + 1].append(Nfur) #B_Neo1
							Datar[spec + item].append(typer[item]) #S_All
							Datar[spec + item + 1].append(typer[item]) #S_Neo1
					
						elif line[item] == 'Aligned' and line[item-1] != 'chimeric':
							if trip == 0:
								Datar[2].append(Nfur)
								Datar[3].append(Nfur)
								trip = 1
							Datar[agar + item].append(Nfur)
							Datar[agar + item + 1].append(Nfur)
							Datar[agor + item].append(typer[item])
							Datar[agor + item + 1].append(typer[item])
					
						elif line[item] == 'Aligned':
							if trip == 0:
								Datar[2].append(Nfur)
								Datar[3].append(Nfur)
								trip = 1  
				
						elif line[item] == 'Unaligned':
							if trip == 0:
								Datar[2].append(Nfur)
								Datar[3].append(Nfur)
								trip = 1
					
						
				###############################################################
				#### Sort NEOF2 sites
				###############################################################
				elif line[4] == 'NeoF2':
					for item in typer.keys():
						back = back + 5
						spec = spec + 5
						agar = agar + 5
						agor = agor + 5
						if 'peak' in line[item]:
							if trip == 0: 
								Datar[2].append(Nfur) #All
								Datar[4].append(Nfur) #Neo2
								trip = 1
							Datar[back + item].append(Nfur) #B_All
							Datar[back + item + 2].append(Nfur) #B_Neo2
							Datar[spec + item].append(typer[item]) #S_All
							Datar[spec + item + 2].append(typer[item]) #S_Neo2
				
						elif line[item] == 'Aligned' and line[item-1] != 'chimeric':
							if trip == 0:
								Datar[2].append(Nfur)
								Datar[4].append(Nfur)
								trip = 1
							Datar[agar + item].append(Nfur)
							Datar[agar + item + 2].append(Nfur)
							Datar[agor + item].append(typer[item])
							Datar[agor + item + 2].append(typer[item])
					
						elif line[item] == 'Aligned':
							if trip == 0:
								Datar[2].append(Nfur)
								Datar[3].append(Nfur)
								trip = 1  
				
						elif line[item] == 'Unaligned':
							if trip == 0:
								Datar[2].append(Nfur)
								Datar[4].append(Nfur)
								trip = 1
			
				###############################################################
				#### Sort all other sites
				###############################################################
				else:
					for item in typer.keys():
						back = back + 5
						spec = spec + 5
						agar = agar + 5
						agor = agor + 5
						if 'peak' in line[item]:
							if decid == 'A':
								if trip == 0: 
									Datar[2].append(Nfur) #All
									trip = 1
								Datar[back + item].append(Nfur)
								Datar[spec + item].append(typer[item]) #S_All
				
						elif line[item] == 'Aligned' and line[item-1] != 'chimeric':
							if decid == 'A':
								if trip == 0:
									Datar[2].append(Nfur) #All
									trip = 1
								Datar[agar + item].append(Nfur)
								Datar[agor + item].append(typer[item])
					
						elif line[item] == 'Aligned':
							if decid == 'A':
								if trip == 0:
									Datar[2].append(Nfur)
									trip = 1
				
						elif line[item] == 'Unaligned':
							if decid == 'A':
								if trip == 0:
									Datar[2].append(Nfur)
									trip = 1
	
		###############################################################
		#### Write master dictionary to output files
		###############################################################					
		for outer in Filer:
			with open(outer, 'w') as temp:
				hold = (Filer.index(outer) + 2)
				if hold in Datar.keys():
					for entry in Datar[hold]:
						temp.write(entry)
				elif hold not in Datar and hold > 10:
					hold1 = hold - 3
					hold2 = hold + 3
					for entry in Datar[hold1]:
						temp.write(entry)
					for entry in Datar[hold2]:
						temp.write(entry)
				