'''
Created on Jun 11, 2010

@author: virushunter2
'''
import re, sys, os
from Bio import SeqIO
import wx
import sqlite3

class Barcode:
    def __init__self(self, filename=[]):
        self.filename = filename

    def read_Barcode_List(self, filename):
        new_contents, blist = [], []
        with open(filename) as file:
            for line in file:
                if not line.strip():
                    continue
                else:
                    new_contents.append(line)
            new_contents = [seq.partition('\n') for seq in new_contents]
            for tuple in new_contents:
                blist.append(tuple[0])
            return blist

    def read_Barcode_Entries(self, filename, listbox):
        bdict = {}
        name_list = []
        blist = [seq.partition('\t') for seq in self.read_Barcode_List(filename)]
        for item in blist:
            bdict[item[2]] = item[0]
        connection = sqlite3.connect('barcodes.db')
        cursor = connection.cursor()
        cursor.execute('''CREATE TABLE IF NOT EXISTS barcodes (id integer primary key, sequence, barcode_name)''')
        cursor.execute('''SELECT barcode_name FROM barcodes''')
        for row in cursor:
            name_list.append(str(row[0]))
        for name, bcode in bdict.items():
            if name not in name_list:
                t = bcode, name
                cursor.execute('''INSERT INTO barcodes (id, sequence, barcode_name) VALUES (NULL, ?, ?)''', t)
                insert = len(filter(lambda x: x<name, listbox.GetItems()))
                listbox.Insert(name, insert)
        connection.commit()
        return bdict

class SuffixArray():
    def lce(self, string1, int1, string2, int2, mismatches):
            res, mm = 0, 0
            for i,j in zip(string1[int1:], string2[int2:]):
                if i == j:
                    res += 1
                else:
                    if mm < mismatches:
                        mm += 1
                        res += 1
                    else:
                        return res
            return res

    def create_array(self, sequence):
        ta = []
        i = 0
        sa = range(len(sequence))
        sa.sort(key=lambda a: sequence[a:])
        for item in sa: #@UnusedVariable
            ta.append(sequence[sa[i]:])
            i += 1
        return ta

    def findPatterns(self, fasta_file, patterns, truncation, mismatches, trim, index, remove, min_length):
        num_seqs, progressMax = 0, 0
        output = {}
        output['problem_barcodes'] = []
        output['no_barcode'] = []
        inv_map = dict(zip(patterns.values(), patterns.keys()))
        for pattern in patterns.values():
            filename = inv_map[pattern]
            output[filename.rstrip('r')] = []
        with open(fasta_file) as f:
            for rec in SeqIO.parse(f, "fasta"):
                progressMax += 1
        dialog = wx.ProgressDialog("Barcode Trimming","Time remaining :", progressMax, style=wx.PD_ELAPSED_TIME | wx.PD_REMAINING_TIME)
        with open(fasta_file) as in_handle:
            for rec in SeqIO.parse(in_handle, "fasta"):
                results, starts, ends = [], [], []
                num_seqs += 1
                if num_seqs < progressMax:
                    dialog.Update(num_seqs)
                else:
                    dialog.Destroy()
                sequence = str(rec.seq)
                ta = self.create_array(sequence)
                n = len(ta) #@UnusedVariable
                for pattern in patterns.values():
                    filename = inv_map[pattern]
                    a = self.findPattern(sequence, ta, pattern, filename, truncation, mismatches, index)
                    if a[0] != None:
                        if filename.rstrip('r') not in results:
                            results.append(filename.rstrip('r'))
                        for s, e in zip(a[1], a[2]):
                            starts.append(s)
                            ends.append(e)
                cuts = zip(starts,ends)
                cuts = sorted(cuts)
                length = len(sequence)
                for start in starts:
                    if start > 50 and start < length - 50:
                        results.append('middle')
                if len(results) == 1:
                    diff, num, prev_end = 0, 0, 0
                    if remove:
                        for start, end in cuts:
                            if len(sequence[:start-diff]) >= min_length and len(sequence[:end-diff:]) >= min_length:
                                sequence = sequence[:start-diff-trim] + sequence[end-diff+trim:]
                            elif len(sequence[:start-diff]) >=min_length:
                                sequence = sequence[:start-diff-trim]
                            elif len(sequence[end-diff:]) >= min_length:
                                sequence = sequence[end-diff+trim:]
                            else:
                                sequence = ""
                            diff += (end-start)
                            if num >= 1:
                                if start - prev_end < min_length:
                                    sequence = ""
                            num += 1
                            prev_end = end
                        if len(sequence) >= 20:
                            tmp = output[results[0]]
                            tmp.append(sequence + '\n')
                            output[results[0]] = tmp
                    else:
                        tmp = output[results[0]]
                        tmp.append(sequence + '\n')
                        output[results[0]] = tmp
                elif len(results) == 0:
                    tmp = output['no_barcode']
                    tmp.append(sequence + '\n')
                    output['no_barcode'] = tmp
                elif len(results) >= 1:
                    tmp = output['problem_barcodes']
                    tmp.append(sequence + '\n')
                    output['problem_barcodes'] = tmp
            return output
            
    def findPattern(self, sequence, suffix_array, pattern, filename, truncation, mismatches, index):
        n = len(suffix_array)
        bcodes, starts, ends = [], [], []
        result = False
        if pattern <= suffix_array[1]:
            ans = 1
        elif pattern > suffix_array[n-1]:
            ans = n-1
        else:
            L = 1
            R = n
            while R - L > 1:
                M = (L + R)/2
                if pattern <= suffix_array[M]:
                    R = M
                else:
                    L = M
            ans = R
        if len(pattern) > len(suffix_array[ans]):
            b = suffix_array[ans]
        else:
            b = suffix_array[ans][:len(pattern)]
        if pattern == b:
            result = True
            target = b
        else:
            s2 = suffix_array[ans-1]
            a = self.lce(pattern, 0, s2, 0, mismatches)
            if a == len(pattern):
                result = True
                target = s2
            else:
                if index < truncation:
                    index += 1
                    self.findPattern(sequence, suffix_array, pattern[1:], filename, truncation, mismatches, index)
                else:
                    result = False
        if result == True:
            for m in re.finditer(target, sequence):
                bcodes.append(filename)
                starts.append(m.start())
                ends.append(m.start()+len(pattern))
            return bcodes, starts, ends
        else:
            return None, None, None

    def writeOutput(self, dict):
        for key in dict.keys():
            filename = key + '.txt'
            with open(filename, 'w') as f:
                num = 0
                for seq in dict[key]:
                    num += 1
                    f.write('>' + str(key) + '.' + str(num) + '\n' + seq)
                    
def main(fasta_file, barcode_file, mismatches, truncation, trim, output, remove, min_length, listbox):
    SA = SuffixArray()
    B = Barcode()
    sequence = fasta_file
    patterns = B.read_Barcode_Entries(barcode_file, listbox)
    index = 0
    a = SA.findPatterns(sequence, patterns, truncation, mismatches, trim, index, remove, min_length)
    if output != "":
        os.chdir(output)
    SA.writeOutput(a)
    
if __name__ == "__main__":
    main(*sys.argv[1:])