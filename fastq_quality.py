'''
Created on Mar 30, 2011

@author: virushunter2
'''
from Bio import SeqIO

def main(infile, outfile, threshold, length):
    threshold = int(threshold)
    length = len(threshold)
    good_reads = (rec for rec in SeqIO.parse(open(infile), "fastq") if min (rec.letter_annotations["phred_quality"]) >= threshold and len(rec) >= length)
    count = SeqIO.write(good_reads, outfile, "fastq")
    print "Saved %i reads" %count