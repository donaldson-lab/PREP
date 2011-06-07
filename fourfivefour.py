'''
Created on Apr 11, 2011

@author: virushunter2
'''

from Bio import SeqIO

def main(infile, infile2, threshold, length):
    length = int(length)
    threshold = int(threshold)
    
    good_quality = []
    
    for record in SeqIO.parse(open(infile), "qual"):
        qscore = record.letter_annotations["phred_quality"]
        flag = False
        i = 0
        for num in qscore:
            if num >= threshold and flag == False:
                i = i + 1
            else:
                flag = True
        good_quality.append(i)
    
    i = 0
    string = ""
    for record in SeqIO.parse(open(infile2), 'fasta'):
        if good_quality[i] >= length:
            string = string + ">" + str(record.id) + '\n' + str(record.seq[:good_quality[i]])+'\n'
            print i
            i = i + 1
        else:
            i = i + 1
    
    with open('454quality_trimmed.txt', 'w') as file:
        file.write(str(string))
        
        


                
        