'''
Created on Apr 11, 2011

@author: virushunter2
'''
from Bio.SeqIO.QualityIO import FastqGeneralIterator

def main(infile, min_length):
    min_length = min_length
    i = 1
    for title, seq, qual in FastqGeneralIterator(open(infile)):
        print i
        i = i + 1
        qual = qual.rstrip("B") #Remove any trailing B characters
        length = len(qual)
        if length >= min_length:
            seq = seq[:length] #trim to match
            handle = open("trimmed.fastq", "a")
            handle.write(">%s\n%s\n" % (title, seq))
            handle.close