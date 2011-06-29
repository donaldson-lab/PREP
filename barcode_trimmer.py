import sys, os, re
import pp #@UnresolvedImport
from Bio.Seq import Seq
from Bio import SeqIO

class Barcode:
    start_list = []

    def __init__self(self,filename=[]):
        self.filename = filename

    def read_Barcode_List(self, filename):
        new_contents = []
        with open(filename) as file:
            contents = file.readlines()
            for line in contents:
                if not line.strip():
                    continue
                else:
                    new_contents.append(line)
            new_contents = [seq.partition('\n') for seq in new_contents]
            blist = []
            for tuple in new_contents:
                blist.append(tuple[0])
            return blist

    def read_Barcode_Entries(self, filename):
        blist = [seq.partition('\t') for seq in self.read_Barcode_List(filename)]
        return blist

    def match_barcodes(self, seq, target, mismatches, trunc, trim, directory, remove):
        current_dir = directory

        sequence_list = []
        os.chdir(current_dir)
        seq_number = 0

        pattern = target[0]
        bcode = pattern
        barcode = list(bcode)
        r_barcode = list(str(Seq(pattern).reverse_complement()))
        b_dict = {}
        f = target[1]
        f = f.split('\t')[0]
        core_pattern = ""
        r_core_pattern = ""
        if trunc == 1:
            start_pattern = '(%(first)s?)'
            end_pattern = '(%last)s)?'
            for char in barcode:
                core_pattern = core_pattern + char
            core_pattern = core_pattern[1:-1]
            pattern = re.compile(start_pattern % {'first':barcode[0]} + core_pattern + end_pattern % {'last':barcode[len(barcode)-1]})
            for char in r_barcode:
                r_core_pattern = r_core_pattern + char
            r_core_pattern = r_core_pattern[1:-1]
            r_pattern = re.compile(start_pattern % {'first':r_barcode[0]} + r_core_pattern + end_pattern % {'last':r_barcode[len(r_barcode)-1]})

        elif trunc == 2:
            start_pattern = '(%(first)s)?(%(second)s)?'
            end_pattern = '(%(penultimate)s)?(%(last)s)?'
            for char in barcode:
                core_pattern = core_pattern + char
            core_pattern = core_pattern[2:-2]
            pattern = re.compile(start_pattern % {'first' :barcode[0], 'second': barcode[1]} + core_pattern + end_pattern % {'penultimate':barcode[len(barcode)-2], 'last':barcode[len(barcode)-1]})
            for char in r_barcode:
                r_core_pattern = r_core_pattern + char
            r_core_pattern = r_core_pattern[2:-2]
            r_pattern = re.compile(start_pattern % {'first' :r_barcode[0], 'second': r_barcode[1]} + r_core_pattern + end_pattern % {'penultimate':r_barcode[len(r_barcode)-2], 'last':r_barcode[len(r_barcode)-1]})
                
        elif trunc == 3:
            start_pattern = '(%(first)s)?(%(second)s)?(%(third)s)?'
            end_pattern = '(%(third_last)s)?(%(penultimate)s)?(%(last)s)?'
            for char in barcode:
                core_pattern = core_pattern + char
            core_pattern = core_pattern[3:-3]
            pattern = re.compile(start_pattern % {'first' :barcode[0], 'second': barcode[1], 'third': barcode[2]} + core_pattern + end_pattern % {'third_last': barcode[len(barcode)-3], 'penultimate':barcode[len(barcode)-2], 'last':barcode[len(barcode)-1]})
            for char in r_barcode:
                r_core_pattern = r_core_pattern + char
            r_core_pattern = r_core_pattern[3:-3]
            r_pattern = re.compile(start_pattern % {'first' :r_barcode[0], 'second': r_barcode[1], 'third': r_barcode[2]} + r_core_pattern + end_pattern % {'third_last': r_barcode[len(r_barcode)-3], 'penultimate':r_barcode[len(r_barcode)-2], 'last':r_barcode[len(r_barcode)-1]})
        else:
            pattern = re.compile(bcode)
            r_pattern = re.compile(r_barcode)

        target = str(seq)
        g = re.split(pattern, target)
        for item in g:
            if item == None:
                item = ''
            h = re.split(r_pattern,item)

        if len(g) == 1 and len(h) == 1:
            sequence_list.append(seq_number)
        else:
            for fragment in h:
                if fragment == None:
                    fragment = ''
                if len(fragment) >=20:
                    if not f in b_dict:
                        b_dict[f] = ""
                        b_dict[f] = b_dict[f] + ">" + f + '\n' + fragment + '\n'

        seq_number += 1
        return sequence_list, target, b_dict

    def remove_adaptor(self, seq, region, trim):
        try:
            pos = seq.find(region)
        except AttributeError:
            pos = seq.seq.find(region)
        if (pos-trim < 0):
            place = 0
        else:
            place = pos - trim
            
        return seq[:place], seq[pos+len(region)+trim:], place
    
    def hamming_distance(self, s1, s2):
        return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))  
  
    
    def trim_adaptor(self, seq, adaptor, mismatches, trim):
        barcode = Barcode()
        distance = barcode.hamming_distance(seq, adaptor)
        i = 0
        while i < len(seq) - len(adaptor):
            if distance <= mismatches:
                return barcode.remove_adaptor(seq, adaptor, trim)
            else:
                i = i+1
                seq = seq[i:]
                
    def num_seqs(self, directory):
        orig_dir = os.getcwd()
        for file in os.walk(directory):
            for item in file[2]:
                if '.txt' in item:
                    if os.getcwd() is not directory:
                        os.chdir(directory)
                    number = 1
                    handle = open(item)
                    outfile = open('tmp.txt', 'a')
                    for line in handle:
                        if line.startswith('>'):
                            line = line.rstrip('\n')
                            outfile.write(line + ".%d" %number + '\n')
                            number += 1
                        else:
                            outfile.write(line)
                    os.rename('tmp.txt', item)
                    os.chdir(orig_dir)

def main(fasta_file, barcode_file, mismatches, trunc, trim, output, remove):
    bcode = Barcode()
    trunc = int(trunc)
    mismatches = int(mismatches)
    trim = int(trim)
    current_dir = os.getcwd()
    num = 0
    files = []
    targets = bcode.read_Barcode_Entries(barcode_file)
    inputs = []
    for target in targets:
        inputs.append([target[0],target[2].split('\t')[0]])
    if output != "":
        os.chdir(output)
   
    while num < len(targets):
            
        for pattern in targets[num]:
            if '\t' in pattern and len(pattern) >1:
                tmp = pattern.split('\t')
                pattern = tmp[0]
                file = open("%s.txt" %pattern, 'a')
                files.append(file)
            num += 1
    
    os.chdir(current_dir)
    inputs = tuple(inputs)
    ppservers = ()
    job_server = pp.Server(ppservers=ppservers)
    print "Starting pp with", job_server.get_ncpus(), "workers"
    with open(fasta_file) as in_handle:
        seq_number = 0
        seq_list = []
        b_dict = {}
        no_string = ""
        for rec in SeqIO.parse(in_handle, "fasta"):
            jobs = [(input, (job_server.submit(bcode.match_barcodes,(rec.seq,input,mismatches,trunc,trim,output,remove,),(),('sys','os','re','from Bio.Seq import *','from Bio.Seq import _dna_complement_table','from Bio import *')))) for input in inputs]
            for input, job in jobs:
                if job()[0] == []:
                    seq_list.append(seq_number)
                    key = ""
                    value = ""
                    for k, v in zip(job()[2].iterkeys(),job()[2].itervalues()):
                        key = k
                        value = v
                    if key not in b_dict:
                        b_dict[key] = ""
                    b_dict[key] = b_dict[key] + value
            if seq_number not in seq_list:
                no_string = no_string + ">No Barcode" + '\n' + str(rec.seq) + '\n'
                os.chdir(current_dir)
            print seq_number
            seq_number += 1

    os.chdir(output)
    for barcode in b_dict.iterkeys():
        file = open('%s.txt' %barcode, 'w')
        file.write(b_dict[barcode])
        file.close()
            
    file = open("no_barcode.txt", 'w')
    file.write(no_string)
    file.close()           
    bcode.num_seqs(output)
    print "Complete"
           

if __name__ == "__main__":
    main(*sys.argv[1:])

