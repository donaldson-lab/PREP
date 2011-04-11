'''
Created on Jun 8, 2010

@author: virushunter2
'''
import sys
from Bio.Seq import Seq
from Bio import SeqIO

class contig():
    def __init__self(self, filename):
        self.filename = filename
        

    def filter_contigs_by_length(self, contig_file, output_file, min_length):
        out_file = open('%s' %output_file, 'a')
        sort_id = {}
        sort_seq = {}
        alist= []
        with open(contig_file) as file:
            for rec in SeqIO.parse(file, 'fasta'):
                length = len(rec.seq)
                if length >= min_length:
                    sort_id[str(rec.id)] = length
                    sort_seq[str(rec.id)] = str(rec.seq)
                    alist = sorted(sort_id.iteritems(), key=lambda (k,v):(v,k), reverse=True)
            for item in alist:
                out_file.write(">" + item[0] + '\n' + sort_seq[item[0]] + '\n')
                    
        print "Filtering Completed"
        print "Sorting Completed"
def main(input_file, output_file, min_length):
    min_length = int(min_length)
    c = contig()
    c.filter_contigs_by_length(input_file, output_file, min_length)


if __name__ == "__main__":
    main(*sys.argv[1:])
   



        