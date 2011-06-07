'''
Created on Apr 11, 2011

@author: virushunter2
'''

import os
import wx
from Bio import SeqIO

def main(infile, threshold, length):
    length = int(length)
    threshold = int(threshold)
    filters = 'FASTA files (*.fasta; *.fastq; *.fa; *.fna; *.seq; *.txt)|*.fasta;*.fastq;*.fa;*.fna;*.seq;*.txt'
    dialog = wx.FileDialog(None, style = wx.OPEN, wildcard = filters)
    if dialog.ShowModal() == wx.ID_OK:
        infile = dialog.GetPath()
    for record in SeqIO.parse(open(infile), "fastq-illumina"):
        qscore = record.letter_annotations["phred_quality"]
        print qscore