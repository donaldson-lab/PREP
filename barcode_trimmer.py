'''
Created on May 28, 2010

@author: Doug Crandell
'''
import sys
import os
import re

from Bio.Seq import Seq #@UnresolvedImport @UnusedImport
from Bio import SeqIO #@UnresolvedImport

class Barcode:
    start_list = []
     
    def __init__self(self, filename = []):
        self.filename = filename;
    
    # Read in the list of barcodes   
    def read_Barcode_List(self, filename):
        with open(filename) as file:
            return file.read().split('\n')[:]
     
    # Read data from barcode list    
    def read_Barcode_Entries(self, filename):
        barcode = Barcode()
        blist = [seq.partition('\t') for seq in barcode.read_Barcode_List(filename)]
        return blist
    
    def match_barcodes(self, seq, filename2, mismatches, trunc, trim, directory, remove):
        current_dir = directory
        match = Barcode()
        targets = match.read_Barcode_Entries(filename2)
        
        sequence_list = []
        barcode_number = 0
        barcode = []
        seq_number = 0
        
        # Create a regular expression to from the barcode to search the sequence
        for pattern in targets:
            pattern = targets[barcode_number][0]
            bcode = targets[barcode_number][0] 
            barcode = list(pattern)
            string = ""
            problem_string = ""
            filename = targets[barcode_number][1]
            filename1 = filename.split('\t')
            filename = filename1[0]
            
            file = open("%s.txt" %filename, 'a')
            core_pattern = ""
            
            if trunc == 1:
                start_pattern = '(%(first)s?)'
                end_pattern = '(%(last)s)?'
                for char in barcode:
                    core_pattern = core_pattern + char
                core_pattern = core_pattern[1:-1]
                pattern = re.compile(start_pattern % {'first' :barcode[0]} + core_pattern + end_pattern % {'last':barcode[21]})
            
            elif trunc == 2:
                start_pattern = '(%(first)s)?(%(second)s)?'
                end_pattern = '(%(penultimate)s)?(%(last)s)?'
                for char in barcode:
                    core_pattern = core_pattern + char
                core_pattern = core_pattern[2:-2]
                pattern = re.compile(start_pattern % {'first' :barcode[0], 'second': barcode[1]} + core_pattern + end_pattern % {'penultimate':barcode[20], 'last':barcode[21]})
            
            elif trunc == 3:
                start_pattern = '(%(first)s)?(%(second)s)?(%(third)s)?'
                end_pattern = '(%(third_last)s)?(%(penultimate)s)?(%(last)s)?'
                for char in barcode:
                    core_pattern = core_pattern + char
                core_pattern = core_pattern[3:-3]
                pattern = re.compile(start_pattern % {'first' :barcode[0], 'second': barcode[1], 'third': barcode[2]} + core_pattern + end_pattern % {'third_last': barcode[19], 'penultimate':barcode[20], 'last':barcode[21]})
            
            else:
                pattern = re.compile(bcode)     
            
            target = (str(seq))           
            g = re.search(pattern, target)
                
            if (g):
                
                match.start_list.append(g.start())
                if remove:
                    trimmed = match.remove_adaptor(str(seq), bcode, trim)
                    
                    if len(trimmed[1]) >=20 and len(trimmed[0]) >=20:
                        string = string +">" + filename + '\n' + trimmed[0] + '\n'
                        string = string +">" + filename + '\n' + trimmed[1] + '\n'
                        sequence_list.append(seq_number)
                        
                    elif len(trimmed[1]) >=20:
                        string = string +">" + filename + '\n' + trimmed[1] + '\n'
                        sequence_list.append(seq_number)
                     
                    elif len(trimmed[0]) >=20:
                        string = string +">" + filename + '\n' + trimmed[0] + '\n'
                        sequence_list.append(seq_number)
                    
                    elif 0 < len(trimmed[0]) < 20 and 0 < len(trimmed[1]) < 20:
                        problem_string = problem_string +">" + filename + '\n' + trimmed[0] + '\n'
                        problem_string = problem_string +">" + filename + '\n' + trimmed[1] + '\n'
                        sequence_list.append(seq_number)
                        
                    elif 0 < len(trimmed[0]) <20:
                        problem_string = problem_string +">" + filename + '\n' + trimmed[0] + '\n'
                        sequence_list.append(seq_number)
                    
                    elif 0 < len(trimmed[1]) <20:
                        problem_string = problem_string +">" + filename + '\n' + trimmed[1] + '\n'
                        sequence_list.append(seq_number)
                else:
                    string = string + ">" + filename + '\n'
                    string = string + str(seq) + '\n'
                    

            barcode_number = barcode_number + 1
            file.write(string)
            file = open('problem_barcodes.txt', 'a')
            file.write(problem_string)
            file.close()
            os.chdir(current_dir)                         
        seq_number = seq_number +1
        
        return sequence_list, targets, match.start_list
    
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
    
    with open(fasta_file) as in_handle:
        seq_number = 0
        seq_list = []
        no_list = []
        start_list = []
        no_string = ""
        for rec in SeqIO.parse(in_handle, "fasta"):
            matches = bcode.match_barcodes(rec.seq, barcode_file, mismatches, trunc, trim, output, remove)
                        
            if matches[2] not in start_list:
                start_list.append(matches[2])
                            
            if matches[0] != []:
                seq_list.append(seq_number)
            
            if seq_number not in seq_list and seq_number not in no_list:
                no_string = no_string + ">No Barcode" + '\n' + str(rec.seq) + '\n'
                os.chdir(current_dir)
            #print seq_number
            seq_number += 1
       
    list = []
    under_30 = 0
    under_50 = 0
    other = 0
           
    for item in start_list[0]:
        list.append(item)
    for item in start_list[1:]:
        list.append(item)
          
            
    for number in list:
        if number <= 30:
            under_30 +=1
        elif number >30 and number <=50:
            under_50 +=1
        else:
            other += 1
    file = open("no_barcode.txt", 'w')
    file.write(no_string)
    file.close()           
    bcode.num_seqs(output)
        
    file = open("Summary.txt", 'w')
    file.write("Number of barcodes beginning before nt #30: %d" %under_30 + '\n')
    file.write("Number of barcodes beginning after nt #30 and before nt #50: %d" % under_50 + '\n')
    file.write("Number of barcodes beginning after nt #50: %d" %other + '\n')
    file.close()
    print "Complete"
           

if __name__ == "__main__":
    main(*sys.argv[1:])
    