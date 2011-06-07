import wx
import os
import fourfivefour
import barcode_trimmer
import quality_trim
import fastq_quality
import illumina_quality
import sqlite3
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter

from threading import Thread
from wx.lib.pubsub import Publisher


class TestThread(Thread):

    def __init__(self):
        Thread.__init__(self)
        self.start()

    def run(self):
        MainFrame.trim_barcodes(MainFrame(None, -1, 'PREP: Pyrosequence Read Extractor and Processor'))


class MainFrame(wx.Frame):
    def __init__(self, parent, id, title):
        wx.Frame.__init__(self, parent, id, title, size=(725,370)) 
        self.mainPanel = mainPanel = wx.Panel(self) 
        sizer = wx.BoxSizer(wx.VERTICAL)

        toolbar = self.CreateToolBar()
        toolbar.AddLabelTool(1, 'Barcode Trimming', wx.Bitmap('barcode.png'))
        toolbar.AddLabelTool(2, 'File Processing', wx.Bitmap('folder.png'))
        toolbar.Realize()
                
        #create five very simple panels and add them to the sizer 
        barcode_panel = wx.Panel(mainPanel) 
        sizer.Add(barcode_panel, 1, wx.EXPAND)

        
        vbox = wx.BoxSizer(wx.VERTICAL)
        self.file_button = wx.Button(barcode_panel, -1, 'Select Sequence File')
        self.barcode_button = wx.Button(barcode_panel, -1, 'Select Barcode File')
        self.add_bc_button = wx.Button(barcode_panel, -1, 'Add Barcode to Database')
        self.submit_button = wx.Button(barcode_panel, -1, 'Submit')
        self.fasta_path_label = wx.StaticText(barcode_panel, -1, 'FASTA filepath:')
        self.barcode_path_label = wx.StaticText(barcode_panel, -1, 'Barcode filepath:')
        self.output_label = wx.StaticText(barcode_panel, -1, 'Output location:')
        self.mismatch_label = wx.StaticText(barcode_panel, -1, "Mismatches:")
        self.truncation_label = wx.StaticText(barcode_panel, -1, 'Truncation:')
        self.trim_label = wx.StaticText(barcode_panel, -1, 'Trim:')
        self.remove = wx.CheckBox(barcode_panel, -1, "Remove Barcodes")

        self.fasta_path = wx.TextCtrl(barcode_panel, -1, size = (300,20))
        self.barcode_path = wx.TextCtrl(barcode_panel, -1, size = (300,20))
        self.output_path = wx.TextCtrl(barcode_panel, -1, size = (300,20))
        self.mismatches = wx.TextCtrl(barcode_panel, -1, size = (50,20))
        self.truncation = wx.TextCtrl(barcode_panel, -1, size = (50,20))
        self.trim = wx.TextCtrl(barcode_panel, -1, size = (50, 20))

        self.barcode_list = []
        connection = sqlite3.connect('barcodes.db')
        cursor = connection.cursor()
        cursor.execute('''SELECT barcode_name from barcodes''')
        for row in cursor:
            self.barcode_list.append(row[0])
        self.barcode_list = sorted(self.barcode_list)
        self.barcode_listbox = wx.ListBox(barcode_panel, 26, wx.DefaultPosition, (125, 150), self.barcode_list, wx.LB_EXTENDED)

        mismatch_box = wx.BoxSizer(wx.HORIZONTAL)
        truncation_box = wx.BoxSizer(wx.HORIZONTAL)
        trim_box = wx.BoxSizer(wx.HORIZONTAL)
        fasta_path_box = wx.BoxSizer(wx.HORIZONTAL)
        barcode_path_box = wx.BoxSizer(wx.HORIZONTAL)
        output_box = wx.BoxSizer(wx.HORIZONTAL)
        button_box = wx.BoxSizer(wx.HORIZONTAL)

        fasta_path_box.AddMany([self.fasta_path_label, (self.fasta_path, 0, wx.LEFT, 19)])
        barcode_path_box.AddMany([self.barcode_path_label, (self.barcode_path, 0, wx.LEFT, 10), (self.remove, 0, wx.LEFT, 10)])
        output_box.AddMany([self.output_label, (self.output_path, 0, wx.LEFT, 12)])
        mismatch_box.AddMany([self.mismatch_label, (self.mismatches, 0, wx.LEFT, 10)])
        truncation_box.AddMany([self.truncation_label, (self.truncation, 0, wx.LEFT, 10)])
        trim_box.AddMany([self.trim_label, (self.trim, 0, wx.LEFT, 10)])
        button_box.AddMany([self.file_button, (self.barcode_button, 0, wx.LEFT, 25), (self.add_bc_button, 0, wx.LEFT, 25), (self.barcode_listbox, 0, wx.LEFT, 50)])
        
        vbox.Add(button_box, 0, wx.ALL, 10)
        vbox.Add(fasta_path_box, 0 , wx.TOP|wx.EXPAND, -120)
        vbox.Add(barcode_path_box, 0, wx.TOP|wx.EXPAND, 10)
        vbox.Add(output_box, 0, wx.TOP|wx.EXPAND, 10)
        vbox.Add(mismatch_box, 0, wx.TOP|wx.EXPAND, 10)
        vbox.Add(truncation_box, 0, wx.TOP|wx.EXPAND, 10)
        vbox.Add(trim_box, 0, wx.TOP|wx.EXPAND, 10)
        vbox.Add(self.submit_button, 0, wx.TOP|wx.LEFT, 10)

        self.file_button.Bind(wx.EVT_LEFT_DOWN, self.select_fasta)
        self.barcode_button.Bind(wx.EVT_LEFT_DOWN, self.select_barcodes)
        self.add_bc_button.Bind(wx.EVT_LEFT_DOWN, self.on_add_barcode)

        self.file_button.Bind(wx.EVT_ENTER_WINDOW, self.enter_fasta_button)
        self.barcode_button.Bind(wx.EVT_ENTER_WINDOW, self.enter_barcode_button)
        self.add_bc_button.Bind(wx.EVT_ENTER_WINDOW, self.enter_add)
        self.fasta_path_label.Bind(wx.EVT_ENTER_WINDOW, self.enter_fasta)
        self.fasta_path.Bind(wx.EVT_ENTER_WINDOW, self.enter_fasta)
        self.barcode_path_label.Bind(wx.EVT_ENTER_WINDOW, self.enter_barcode)
        self.barcode_path.Bind(wx.EVT_ENTER_WINDOW, self.enter_barcode)
        self.output_label.Bind(wx.EVT_ENTER_WINDOW, self.enter_output)
        self.output_path.Bind(wx.EVT_ENTER_WINDOW, self.enter_output)
        self.mismatch_label.Bind(wx.EVT_ENTER_WINDOW, self.enter_mismatch)
        self.mismatches.Bind(wx.EVT_ENTER_WINDOW, self.enter_mismatch)
        self.truncation.Bind(wx.EVT_ENTER_WINDOW, self.enter_truncation)
        self.trim.Bind(wx.EVT_ENTER_WINDOW, self.enter_trim)
        self.truncation_label.Bind(wx.EVT_ENTER_WINDOW, self.enter_truncation)
        self.trim_label.Bind(wx.EVT_ENTER_WINDOW, self.enter_trim)
        self.submit_button.Bind(wx.EVT_LEFT_DOWN, self.onButton)
        self.submit_button.Bind(wx.EVT_ENTER_WINDOW, self.enter_submit)
        self.barcode_listbox.Bind(wx.EVT_ENTER_WINDOW, self.enter_listbox)
        self.remove.Bind(wx.EVT_ENTER_WINDOW, self.enter_remove)
        self.remove.SetValue(True)

        sizer.Hide(0)
        
        processing_panel = wx.Panel(mainPanel)
        vbox2 = wx.BoxSizer(wx.VERTICAL)
        threshold_box = wx.BoxSizer(wx.HORIZONTAL)
        length_box = wx.BoxSizer(wx.HORIZONTAL)
        quality_box = wx.BoxSizer(wx.HORIZONTAL)
        fff_box = wx.BoxSizer(wx.HORIZONTAL)
        self.threshold_label = wx.StaticText(processing_panel, -1, 'Trimming Threshold:')
        self.threshold = wx.TextCtrl(processing_panel, -1, size = (50,20))
        self.convert_button = wx.Button(processing_panel, -1, 'Convert Illumina')
        self.quality_trim_button = wx.Button(processing_panel, -1, 'FASTQ Sanger Trim')
        #self.quality_filter_button = wx.Button(processing_panel, -1, 'Quality Filter')
        self.convert_button.Bind(wx.EVT_LEFT_DOWN, self.convert_ill_fasta)
        self.convert_button.Bind(wx.EVT_ENTER_WINDOW, self.enter_convert)
        self.quality_trim_button.Bind(wx.EVT_LEFT_DOWN, self.fastq_trim)
        self.quality_trim_button.Bind(wx.EVT_ENTER_WINDOW, self.enter_quality_trim)
        '''self.quality_filter_button.Bind(wx.EVT_LEFT_DOWN, self.quality_filter)
        self.quality_filter_button.Bind(wx.EVT_ENTER_WINDOW, self.enter_quality_filter)'''
        self.fff_button = wx.Button(processing_panel, -1, '454 Trim')
        self.fff_button.Bind(wx.EVT_LEFT_DOWN, self.fourfivefour_quality)
        self.illumina_quality_button = wx.Button(processing_panel, -1, 'Illumina 1.5+ Trim')
        self.illumina_quality_button.Bind(wx.EVT_LEFT_DOWN, self.illumina_quality)
        self.illumina_quality_button.Bind(wx.EVT_ENTER_WINDOW, self.enter_illumina)
        self.threshold.Bind(wx.EVT_ENTER_WINDOW, self.enter_threshold)
        self.threshold_label.Bind(wx.EVT_ENTER_WINDOW, self.enter_threshold)
        self.length_label = wx.StaticText(processing_panel, -1, 'Minimum Sequence Length:')
        self.length = wx.TextCtrl(processing_panel, -1, size = (50,20))
        self.length_label.Bind(wx.EVT_ENTER_WINDOW, self.enter_length)
        self.length.Bind(wx.EVT_ENTER_WINDOW, self.enter_length)
        self.fff_button.Bind(wx.EVT_ENTER_WINDOW, self.enter_fff)
        self.length.SetValue(str(20))
        self.threshold.SetValue(str(25))
        
        
        threshold_box.AddMany([self.threshold_label, (self.threshold, 0, wx.LEFT, 10)])
        length_box.AddMany([(self.length_label), (self.length, 0, wx.LEFT, 10)])
        quality_box.AddMany([self.quality_trim_button,  (threshold_box, 0, wx.LEFT, 10)])
        fff_box.AddMany([self.illumina_quality_button, (length_box, 0, wx.LEFT, 10)])
        
        vbox2.Add(self.convert_button, 0, wx.TOP|wx.LEFT, 10)
        vbox2.Add(quality_box, 0, wx.TOP|wx.LEFT, 10)
        vbox2.Add(fff_box, 0, wx.TOP|wx.LEFT, 10)
        vbox2.Add(self.fff_button, 0, wx.TOP|wx.LEFT, 10)
        sizer.Add(processing_panel, 1, wx.EXPAND)
        
        
        sizer.Hide(1)

        self.Bind(wx.EVT_TOOL, self.onBarcode, id=1)               
        self.Bind(wx.EVT_TOOL, self.onProcessing, id=2)
                
        #show the first panel 
        sizer.Show(0) 
        mainPanel.Sizer = sizer
        barcode_panel.Sizer = vbox
        processing_panel.Sizer = vbox2
        

        self.sb = self.CreateStatusBar()
        self.Centre()
        self.Show(True)

    def onButton(self, event):
        TestThread()

    def onBarcode(self, event):
        self.mainPanel.Sizer.Hide(1)
        self.mainPanel.Sizer.Show(0)
        self.mainPanel.Layout()
        event.Skip()
        
                
    def onProcessing(self, event):
        self.mainPanel.Sizer.Hide(0)
        self.mainPanel.Sizer.Show(1)
        self.mainPanel.Layout() 
        event.Skip()

    def on_add_barcode(self, event):
        add = AddBarcode(None, -1, 'Add Barcode')
        add.ShowModal()
        if isinstance(add.bc_name, str):
            name = add.bc_name
        else:
            name = str(add.bc_name.GetValue())
        add.Destroy()
        if name != 'Barcode Name':
            insert = len(filter(lambda x: x<name, self.barcode_listbox.GetItems()))
            self.barcode_listbox.Insert(name, insert)

    def select_fasta(self, event):
        filters = 'FASTA files (*.fasta; *.fastq; *.fa; *.fna; *.seq; *.txt)|*.fasta;*.fastq;*.fa;*.fna;*.seq;*.txt'
        dialog = wx.FileDialog(None, style = wx.OPEN, wildcard = filters)
        if dialog.ShowModal() == wx.ID_OK:
            self.fasta = dialog.GetPath()
            self.fasta_path.SetValue(str(self.fasta))

    def select_barcodes(self, event):
        filters = 'Text files (*.txt|*.txt'
        dialog = wx.FileDialog(None, style = wx.OPEN, wildcard = filters)
        if dialog.ShowModal() == wx.ID_OK:
            self.barcodes = dialog.GetPath()
            self.barcode_path.SetValue(str(self.barcodes))

    def trim_barcodes(self):
        if str(self.barcode_path.GetValue()) != '':
            barcode_trimmer.main(str(self.fasta_path.GetValue()), str(self.barcode_path.GetValue()), int(self.mismatches.GetValue()), int(self.truncation.GetValue()), int(self.trim.GetValue()), str(self.output_path.GetValue()), int(self.remove.GetValue()))
        else:
            names = {}
            name_list = []
            connection = sqlite3.connect('barcodes.db')
            cursor = connection.cursor()
            cursor.execute('''SELECT sequence, barcode_name, barcode_id from barcodes''')
            for row in cursor:
                name_list.append(row[1])
                names[row[1]] = [row[0], row[2]]
            sorted_names = self.sorted_dict_values(names)
            name_list = sorted(name_list)
            file = open('tmp_barcodes.txt', 'w')
            for num in self.barcode_listbox.GetSelections():
                file.write(str(sorted_names[num][0]) + '\t' + str(name_list[num]) + '\t' + str(sorted_names[num][1]) + '\n')
            file.close()
            barcode_file = os.getcwd + '/tmp_barcodes.txt'
            self.barcode_path.SetValue(barcode_file)
            barcode_trimmer.main(str(self.fasta_path.GetValue()), str(self.barcode_path.GetValue()), int(self.mismatches.GetValue()), int(self.truncation.GetValue()), int(self.trim.GetValue()), str(self.output_path.GetValue()), int(self.remove.GetValue()))
            os.remove('/tmp_barcodes.txt')

    def sorted_dict_values(self, adict):
        items = adict.items()
        items.sort()
        return [value for key, value in items]

    def enter_mismatch(self, event):
        self.sb.SetStatusText('Enter the maximum number of mismatches between the barcode and the sequence')
        event.Skip()

    def enter_truncation(self, event):
        self.sb.SetStatusText('Enter the maximum number of nucleotides which can be truncated from either end of the barcode')
        event.Skip()

    def enter_trim(self, event):
        self.sb.SetStatusText('Enter the number of additional bases to trim off the sequence after the barcode has been removed')
        event.Skip()
        
    def enter_fasta(self, event):
        self.sb.SetStatusText('Type in the filepath of the FASTA file or select "Select Sequence File" to select the filepath via dropdown menu')
        event.Skip()
        
    def enter_fasta_button(self, event):
        self.sb.SetStatusText('Select FASTA file via dropdown menu')
        event.Skip()

    def enter_barcode_button(self, event):
        self.sb.SetStatusText('Select barcode file via dropdown menu')
        event.Skip()

    def enter_add(self, event):
        self.sb.SetStatusText('Add barcode information to the database for permanent storage')
        event.Skip()

    def enter_output(self, event):
        self.sb.SetStatusText('Type in the filepath of the output directory')
        event.Skip()

    def enter_barcode(self, event):
        self.sb.SetStatusText('Type in the filepath of the barcode file or select "Select Barcode File" to select the filepath via dropdown menu')
        event.Skip()

    def enter_listbox(self, event):
        self.sb.SetStatusText('Select barcodes to trim from FASTA file')
        event.Skip()

    def enter_remove(self, event):
        self.sb.SetStatusText('Option to remove barcodes when trimming or bin sequences by barcode')
        event.Skip()

    def enter_submit(self, event):
        self.sb.SetStatusText('Submit parameters for barcode trimming')
        event.Skip()
    
    def enter_quality_trim(self, event):
        self.sb.SetStatusText('Trim Sanger style FASTQ files using PHRED scores with ASCII offset 33, allowing PHRED scores from 0 to 93. Trim nucleotides with PHRED scores below trimming threshold')
        event.Skip()
    
    '''def enter_quality_filter(self, event):
        self.sb.SetStatusText('Filter Sanger style FASTQ files using PHRED scores with ASCII offset 33, allowing PHRED scores from 0 to 93. Filters out sequences with bases below trimming threshold.')
        event.Skip()'''
        
    def enter_threshold(self, event):
        self.sb.SetStatusText('Set threshold value for quality trimming')
        event.Skip()

    def enter_convert(self, event):
        self.sb.SetStatusText('Convert Illumina to FASTA format')
        event.Skip()
        
    def enter_length(self, event):
        self.sb.SetStatusText('Set minimum length for trimmed sequences')
        event.Skip()
        
    def enter_fff(self, event):
        self.sb.SetStatusText('Trim 454 data based off of quality file: select quality file first, then sequence file')
        event.Skip()
        
    def enter_illumina(self, event):
        self.sb.SetStatusText('Trim Illumina 1.5+ files based on quality scores')
    
    def convert_ill_fasta(self, event):
        filters = 'Text files (*.txt)|*.txt'
        dialog = wx.FileDialog(None, style = wx.OPEN, wildcard = filters)
        if dialog.ShowModal() == wx.ID_OK:
            self.illumina = dialog.GetPath()
        self.filename = self.illumina.split('/')[-1]
        self.filename = self.filename + '.fasta'
        records = SeqIO.parse(open(self.illumina), "fastq-illumina")
        handle = open(self.filename, "w")
        count = FastaWriter(handle, wrap=80).write_file(records)
        handle.close()
        print "Converted %i records" %count
        
    '''def quality_filter(self, event):
        filters = 'FASTA files (*.fasta; *.fastq; *.fa; *.fna; *.seq; *.txt)|*.fasta;*.fastq;*.fa;*.fna;*.seq;*.txt'
        dialog = wx.FileDialog(None, style = wx.OPEN, wildcard = filters)
        if dialog.ShowModal() == wx.ID_OK:
            in_file = dialog.GetPath()
        threshold = int(self.threshold.GetValue())
        length = int(self.length.GetValue()) 
        out_file = in_file + '.quality_trimmed.fasta'
        if in_file != "":
            fastq_quality.main(in_file, out_file, threshold, length)
        else:
            dlg = wx.MessageDialog(self, "", 'Please enter an input file')
            dlg.ShowModal()
            dlg.Destroy()'''
            
    def fastq_trim(self,event):
        filters = 'FASTA files (*.fasta; *.fastq; *.fa; *.fna; *.seq; *.txt)|*.fasta;*.fastq;*.fa;*.fna;*.seq;*.txt'
        dialog = wx.FileDialog(None, style = wx.OPEN, wildcard = filters)
        if dialog.ShowModal() == wx.ID_OK:
            in_file = dialog.GetPath()
        out_file = in_file + '.good_quality.fasta'
        threshold = int(self.threshold.GetValue())
        if in_file != "":
            quality_trim.main(in_file, out_file, threshold)
        else:
            dlg = wx.MessageDialog(self, "", "Please enter an input file")
            dlg.ShowModal()
            dlg.Destroy()
            
    def fourfivefour_quality(self, event):
        filters = 'FASTA files (*.fasta; *.fastq; *.fa; *.fna; *.seq; *.txt)|*.fasta;*.fastq;*.fa;*.fna;*.seq;*.txt'
        dialog = wx.FileDialog(None, style = wx.OPEN, wildcard = filters)
        if dialog.ShowModal() == wx.ID_OK:
            in_file = dialog.GetPath()
        dialog = wx.FileDialog(None, style = wx.OPEN, wildcard = filters)
        if dialog.ShowModal() == wx.ID_OK:
            in_file2 = dialog.GetPath()
        threshold = int(self.threshold.GetValue())
        length = int(self.length.GetValue())
        if in_file != "":
            fourfivefour.main(in_file, in_file2, threshold, length)
            
    def illumina_quality(self, event):
        filters = 'FASTA files (*.fasta; *.fastq; *.fa; *.fna; *.seq; *.txt)|*.fasta;*.fastq;*.fa;*.fna;*.seq;*.txt'
        dialog = wx.FileDialog(None, style = wx.OPEN, wildcard = filters)
        if dialog.ShowModal() == wx.ID_OK:
            in_file = dialog.GetPath()
        length = int(self.length.GetValue())
        if in_file != "":
            illumina_quality.main(in_file, length)
            
class db():
    
    def __init__(self):
        self = self
        
    def create_db(self):
        connection = sqlite3.connect('barcodes.db')
        cursor = connection.cursor()
        cursor.execute('''CREATE TABLE IF NOT EXISTS barcodes (id integer primary key, sequence, barcode_name, barcode_id)''')
        connection.commit()
        
class AddBarcode(wx.Dialog):
    def __init__(self, parent, id, title):
        wx.Dialog.__init__(self, parent, id, title, size=(300, 200))

        vbox = wx.BoxSizer(wx.VERTICAL)
        self.sequence = wx.TextCtrl(self, -1, 'Barcode Sequence', size = (300,25))
        self.bc_name = wx.TextCtrl(self, -1, 'Barcode Name', size = (300,25))
        self.bc_id = wx.TextCtrl(self, -1, 'Barcode ID', size = (300,25))
        vbox.AddMany([self.sequence, (self.bc_name, 0, wx.TOP, 10), (self.bc_id, 0, wx.TOP, 10)])
        
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        okButton = wx.Button(self, -1, 'Ok', size=(70, 30))
        okButton.Bind(wx.EVT_LEFT_DOWN, self.add)
        closeButton = wx.Button(self, -1, 'Close', size=(70, 30))
        closeButton.Bind(wx.EVT_LEFT_DOWN, self.close)
        hbox.Add(okButton, 1)
        hbox.Add(closeButton, 1, wx.LEFT, 5)

        vbox.Add(hbox, 1, wx.ALIGN_CENTER | wx.TOP | wx.BOTTOM, 10)

        self.SetSizer(vbox)
        self.Centre()
        
    def close(self, event):
        self.bc_name = "Barcode Name"
        self.Close()
        
    def add(self, event):
        connection = sqlite3.connect('barcodes.db')
        cursor = connection.cursor()    
        if str(self.sequence.GetValue()) == 'Barcode Sequence':
            dlg = wx.MessageDialog(self, '', 'Please enter a nucleotide sequence', wx.OK)
            dlg.ShowModal()
            dlg.Destroy()
        elif str(self.bc_id.GetValue()) == "Barcode Name" or str(self.bc_id.GetValue()) == 'Barcode ID':
            dlg = wx.MessageDialog(self, '', 'Please enter a unique Barcode ID and Barcode Name', wx.OK)
            dlg.ShowModal()
        else:
            t = str(self.sequence.GetValue().upper()), str(self.bc_name.GetValue()), str(self.bc_id.GetValue())
            cursor.execute('''INSERT INTO barcodes (id, sequence, barcode_name, barcode_id) VALUES (NULL, ?, ?, ?)''', t)
            connection.commit()
            self.Close()
            
if __name__ == "__main__":
    if not os.path.exists('barcodes.db'):
        db.create_db(db())          
    app = wx.App(False)
    MainFrame(None, -1, 'PREP: Pyrosequence Read Extractor and Processor')
    app.MainLoop()
