#!/usr/env python3

'''
Example Command Line Execution
python3 Motif_Meta_Summary.py
'''

#######################################################################################
##### Define Constants
#######################################################################################

import sys
import os

CURRENT_DIR = './'

home = CURRENT_DIR + '/Data/Final/'
file1 = CURRENT_DIR + 'Data/DE_Novo_Summary.txt'
file2 = CURRENT_DIR + 'Data/TE_Mediate_Summary.txt'

#######################################################################################
##### Define Functions
#######################################################################################
'''
This function summarizes the number of Motif/peak overlap sites fall into each of provided 
conservation categories.
'''
def Count_Summary(a_file):
    DNhold = [0,0,0,0,0,0,0,0]
    TEhold = [0,0,0,0,0,0,0,0] 
    with open(a_file, 'r') as temp:
        for line in temp:
            line = line.rstrip('\n').split('\t')
            if line[4] == 'N1' and line[6] == 'Safe':
                if line[3] == '1.Nfur.Novel':
                    DNhold[0] = DNhold[0] + 1
                elif line[3] == '2.African.Novel':
                    DNhold[1] = DNhold[1] + 1
                elif line[3] == '3.Killifish.Novel':
                    DNhold[2] = DNhold[2] + 1
                elif line[3] == '9.Conserved' or line[4] == '4.Translocated':
                    DNhold[3] = DNhold[3] + 1
                elif line[3] == '5.TE.Nfur':
                    TEhold[0] = TEhold[0] + 1
                elif line[3] == '6.TE.African':
                    TEhold[1] = TEhold[1] + 1
                elif line[3] == '7.TE.Killifish':
                    TEhold[2] = TEhold[2] + 1
                elif line[3] == '8.TE.Conserved':
                    TEhold[3] = TEhold[3] + 1
                
            elif line[4] == 'N2' and line[6] == 'Safe':
                if line[3] == '1.Nfur.Novel':
                    DNhold[4] = DNhold[4] + 1
                elif line[3] == '2.African.Novel':
                    DNhold[5] = DNhold[5] + 1
                elif line[3] == '3.Killifish.Novel':
                    DNhold[6] = DNhold[6] + 1
                elif line[3] == '9.Conserved' or line[4] == '4.Translocated':
                    DNhold[7] = DNhold[7] + 1
                elif line[3] == '5.TE.Nfur':
                    TEhold[4] = TEhold[4] + 1
                elif line[3] == '6.TE.African':
                    TEhold[5] = TEhold[5] + 1
                elif line[3] == '7.TE.Killifish':
                    TEhold[6] = TEhold[6] + 1
                elif line[4] == '8.TE.Conserved':
                    TEhold[7] = TEhold[7] + 1
        fini = [DNhold,TEhold]
        return(fini)

'''
This function takes the counts of motif/peak overlaps generated in COUNT_SUMMARY and divides
then by the total number of sites to generate a percentage for each category.
'''
def Percent_Summary(a_list, a_sum1, a_sum2):
    hold = []
    track = 0
    for item in a_list:
        if track < 4:
            if a_sum1 != 0:
                per = item/a_sum1
                hold.append(str(per))
            else:
                per = 0
                hold.append(str(per))
        else:
            if a_sum2 != 0:
                per = item/a_sum2
                hold.append(str(per))
            else:
                per = 0
                hold.append(str(per))
        track = track + 1
        
    return(hold)  


#######################################################################################
##### Generate Summary Table
#######################################################################################

with open(file1, 'w') as out1, open(file2, 'w') as out2:
    head = 'Motif\tGroup\tN1NF\tN1AF\tN1KF\tN1CS\tN2NF\tN2AF\tN2KF\tN2CS\tP_N1NF\tP_N1AF\tP_N1KF\tP_N1CS\tP_N2NF\tP_N2AF\tP_N2KF\tP_N2CS\n'
    out1.write(head)
    out2.write(head)
    for file in os.listdir(home):
        if '_Final.txt' in file:
            name = file.split('_')
            if name[1] == 'Diapause' or name[1] == 'Development':
                name = name[0] + '\t' + name[1]
            else:
                name = name[0] + '_' + name[1] + '\t' + name[2]
            hold = Count_Summary(home + file)
            DNL = hold[0]
            TEL = hold[1]
        
            DNtotn1 = 0
            DNtotn2 = 0
            TEtotn1 = 0
            TEtotn2 = 0
            
            dnc = 0
            tec = 0
            
            for item in DNL:
                if dnc < 4:
                    DNtotn1 = DNtotn1 + item
                else:
                    DNtotn2 = DNtotn2 + item
                dnc = dnc + 1
                    
            for item in TEL:
                if tec < 4:
                    TEtotn1 = TEtotn1 + item
                else:
                    TEtotn2 = TEtotn2 + item
                tec = tec + 1
            
            DNL2 = Percent_Summary(DNL, DNtotn1, DNtotn2)
            TEL2 = Percent_Summary(TEL, TEtotn1, TEtotn2)
            
            writ1 = '\t' + str(DNL[0]) + '\t' + str(DNL[1]) + '\t' + str(DNL[2]) + '\t' + str(DNL[3]) + '\t' + str(DNL[4]) + '\t' + str(DNL[5]) + '\t' + str(DNL[6]) + '\t' + str(DNL[7]) + '\t'
            writ2 = '\t' + str(TEL[0]) + '\t' + str(TEL[1]) + '\t' + str(TEL[2]) + '\t' + str(TEL[3]) + '\t' + str(TEL[4]) + '\t' + str(TEL[5]) + '\t' + str(TEL[6]) + '\t' + str(TEL[7]) + '\t'
        
            out1.write(name + writ1 + '\t'.join(DNL2) + '\n')
            out2.write(name + writ2 + '\t'.join(TEL2) + '\n')
            
