'''
Created on Jun 2, 2010

@author: virushunter2
'''

import os
import sys
from Bio import SeqIO

def main(in_file, out_file, threshold, length):
    length = int(length)
    threshold = int(threshold)
    num_seq = 0
    trim_seq = 0
    
    for record in SeqIO.parse(open(in_file), "fastq-illumina"):
        qscore = record.letter_annotations["phred_quality"]
        i = 0
        j = 0
        for number in qscore:
            
            if number > threshold and j == 0:
                i +=1
            else:
                j = 1
        new_seq = record.seq[:i]
        
        if len(new_seq) >=length:
            
            if os.path.exists(out_file):
                file = open(out_file, 'a')
                file.write('>'+ record.id + '\n')
                file.write(str(new_seq) + '\n')
                trim_seq +=1
            else:
                file = open(out_file, 'w')
                file.write('>'+ record.id + '\n')
                file.write(str(new_seq) + '\n')
                trim_seq +=1
        num_seq +=1
        
    print "Total Records: %d" %num_seq
    print "%d sequences " %trim_seq + "of at least %d nt " %length + 'were added to the new file'
        
if __name__ == "__main__":
        main(*sys.argv[1:])
        print "Completed"
        

