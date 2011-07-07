import wx
import os
import fourfivefour
import suffix_array #@UnresolvedImport
import quality_trim
#import fastq_quality
import illumina_quality
import sqlite3
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter

class MainFrame(wx.Frame):
    def __init__(self, parent, id, title):
        wx.Frame.__init__(self, parent, id, title, size=(725,370)) 
        self.mainPanel = mainPanel = wx.Panel(self) 
        sizer = wx.BoxSizer(wx.VERTICAL)

        toolbar = self.CreateToolBar()
        toolbar.AddLabelTool(1, 'Barcode Trimming', wx.Bitmap('barcode.png'))
        toolbar.AddLabelTool(2, 'File Processing', wx.Bitmap('folder.png'))
        toolbar.Realize()
                
        barcode_panel = wx.Panel(mainPanel) 
        sizer.Add(barcode_panel, 1, wx.EXPAND)

        vbox = wx.BoxSizer(wx.VERTICAL)
        self.file_button = wx.Button(barcode_panel, -1, 'Select Sequence File')
        self.barcode_button = wx.Button(barcode_panel, -1, 'Select Barcode File')
        self.add_bc_button = wx.Button(barcode_panel, -1, 'Add Barcode to Database')
        self.del_bc_button = wx.Button(barcode_panel, -1, 'Delete Barcodes')
        self.submit_button = wx.Button(barcode_panel, -1, 'Submit')
        self.fasta_path_label = wx.StaticText(barcode_panel, -1, 'FASTA filepath:')
        self.barcode_path_label = wx.StaticText(barcode_panel, -1, 'Barcode filepath:')
        self.output_label = wx.StaticText(barcode_panel, -1, 'Output location:')
        self.mismatch_label = wx.StaticText(barcode_panel, -1, "Mismatches:")
        self.truncation_label = wx.StaticText(barcode_panel, -1, 'Truncation:')
        self.trim_label = wx.StaticText(barcode_panel, -1, 'Trim:')
        self.remove = wx.CheckBox(barcode_panel, -1, "Remove Barcodes")
        self.length_label = wx.StaticText(barcode_panel, -1, "Minimum Sequence Length:")

        self.fasta_path = wx.TextCtrl(barcode_panel, -1, size = (300,20))
        self.barcode_path = wx.TextCtrl(barcode_panel, -1, size = (300,20))
        self.output_path = wx.TextCtrl(barcode_panel, -1, size = (300,20))
        self.mismatches = wx.TextCtrl(barcode_panel, -1, size = (20,20))
        self.truncation = wx.TextCtrl(barcode_panel, -1, size = (20,20))
        self.trim = wx.TextCtrl(barcode_panel, -1, size = (20, 20))
        self.min_length = wx.TextCtrl(barcode_panel, -1, size = (40,20))
        self.min_length.SetValue(str(20))

        self.barcode_list = []
        connection = sqlite3.connect('barcodes.db')
        cursor = connection.cursor()
        cursor.execute('''SELECT barcode_name from barcodes''')
        for row in cursor:
            self.barcode_list.append(row[0])
        self.barcode_list = sorted(self.barcode_list)
        self.barcode_label = wx.StaticText(barcode_panel, -1, "Barcodes")
        self.barcode_listbox = wx.ListBox(barcode_panel, 26, wx.DefaultPosition, (125, 300), self.barcode_list, wx.LB_EXTENDED)

        mismatch_box = wx.BoxSizer(wx.HORIZONTAL)
        truncation_box = wx.BoxSizer(wx.HORIZONTAL)
        trim_box = wx.BoxSizer(wx.HORIZONTAL)
        fasta_path_box = wx.BoxSizer(wx.HORIZONTAL)
        barcode_path_box = wx.BoxSizer(wx.HORIZONTAL)
        length_box = wx.BoxSizer(wx.HORIZONTAL)
        output_box = wx.BoxSizer(wx.HORIZONTAL)
        button_box = wx.BoxSizer(wx.HORIZONTAL)
        barcode_box = wx.BoxSizer(wx.VERTICAL)
        barcode_box.AddMany([self.barcode_label, self.barcode_listbox])
        self.submit_button.Bind(wx.EVT_BUTTON, self.trim_barcodes)
        self.submit_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.submit_sizer.AddMany([self.submit_button, (self.del_bc_button, 0, wx.LEFT, 375)])

        fasta_path_box.AddMany([self.fasta_path_label, (self.fasta_path, 0, wx.LEFT, 19)])
        barcode_path_box.AddMany([self.barcode_path_label, (self.barcode_path, 0, wx.LEFT, 10), (self.remove, 0, wx.LEFT, 10)])
        output_box.AddMany([self.output_label, (self.output_path, 0, wx.LEFT, 12)])
        mismatch_box.AddMany([self.mismatch_label, (self.mismatches, 0, wx.LEFT, 10)])
        truncation_box.AddMany([self.truncation_label, (self.truncation, 0, wx.LEFT, 10)])
        trim_box.AddMany([self.trim_label, (self.trim, 0, wx.LEFT, 10)])
        length_box.AddMany([self.length_label, (self.min_length, 0, wx.LEFT, 10)])
        button_box.AddMany([self.file_button, (self.barcode_button, 0, wx.LEFT, 25), (self.add_bc_button, 0, wx.LEFT, 25), (barcode_box, 0, wx.LEFT, 50)])
        
        vbox.Add(button_box, 0, wx.ALL|wx.EXPAND, 10)
        vbox.Add(fasta_path_box,0, wx.TOP|wx.EXPAND, -280)
        vbox.Add(barcode_path_box, 1, wx.TOP|wx.EXPAND, 10)
        vbox.Add(output_box, 1, wx.TOP|wx.EXPAND, 10)
        vbox.Add(mismatch_box, 1, wx.TOP|wx.EXPAND, 10)
        vbox.Add(truncation_box, 1, wx.TOP|wx.EXPAND, 10)
        vbox.Add(trim_box, 1, wx.TOP|wx.EXPAND, 10)
        vbox.Add(length_box, 1, wx.TOP|wx.EXPAND, 10)
        vbox.Add(self.submit_sizer, 1, wx.TOP|wx.EXPAND, 10)

        self.file_button.Bind(wx.EVT_BUTTON, self.select_fasta)
        self.barcode_button.Bind(wx.EVT_BUTTON, self.select_barcodes)
        self.add_bc_button.Bind(wx.EVT_BUTTON, self.on_add_barcode)
        self.del_bc_button.Bind(wx.EVT_BUTTON, self.delete_barcodes)

        self.remove.SetValue(True)
        self.truncation.SetValue(str(2))
        self.mismatches.SetValue(str(2))
        self.trim.SetValue(str(0))
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
        self.convert_button.Bind(wx.EVT_BUTTON, self.convert_ill_fasta)
        self.quality_trim_button.Bind(wx.EVT_BUTTON, self.fastq_trim)
        #self.quality_filter_button.Bind(wx.EVT_BUTTON, self.quality_filter)
        self.fff_button = wx.Button(processing_panel, -1, '454 Trim')
        self.fff_button.Bind(wx.EVT_BUTTON, self.fourfivefour_quality)
        self.illumina_quality_button = wx.Button(processing_panel, -1, 'Illumina 1.5+ Trim')
        self.illumina_quality_button.Bind(wx.EVT_BUTTON, self.illumina_quality)
        self.length_label = wx.StaticText(processing_panel, -1, 'Minimum Sequence Length:')
        self.length = wx.TextCtrl(processing_panel, -1, size = (50,20))
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
        self.SetMinSize(self.GetSize())
        self.SetMaxSize(self.GetSize())
        self.Centre()
        self.Show(True)

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

    def trim_barcodes(self, event):
        if str(self.barcode_path.GetValue()) != '':
            try:
                suffix_array.main(str(self.fasta_path.GetValue()), str(self.barcode_path.GetValue()), int(self.mismatches.GetValue()), int(self.truncation.GetValue()), int(self.trim.GetValue()), str(self.output_path.GetValue()), int(self.remove.GetValue()), int(self.min_length.GetValue()))
            except IOError:
                dlg = wx.MessageDialog(None, "Please enter a valid filepath", "Invalid Filepath")
                dlg.ShowModal()
                dlg.Destroy()
        else:
            names = {}
            name_list = []
            connection = sqlite3.connect('barcodes.db')
            cursor = connection.cursor()
            cursor.execute('''SELECT sequence, barcode_name FROM barcodes''')
            for row in cursor:
                name_list.append(row[1])
                names[row[1]] = row[0]
            sorted_names = self.sorted_dict_values(names)
            name_list = sorted(name_list)
            file = open('tmp_barcodes.txt', 'w')
            for num in self.barcode_listbox.GetSelections():
                file.write(str(sorted_names[num]) + '\t' + str(name_list[num]) + '\n')
            file.close()
            barcode_file = os.getcwd() + '/tmp_barcodes.txt'
            orig_dir = os.getcwd()
            self.barcode_path.SetValue(barcode_file)
            try:
                suffix_array.main(str(self.fasta_path.GetValue()), str(self.barcode_path.GetValue()), int(self.mismatches.GetValue()), int(self.truncation.GetValue()), int(self.trim.GetValue()), str(self.output_path.GetValue()), int(self.remove.GetValue()), int(self.min_length.GetValue()))
            except IOError:
                dlg = wx.MessageDialog(None, "Invalid Filepath")
                dlg.ShoModal()
                dlg.Destroy()
            os.remove(str(orig_dir) + '/tmp_barcodes.txt')
            
    def delete_barcodes(self, event):
        name_list = []
        names = {}
        connection = sqlite3.connect('barcodes.db')
        cursor = connection.cursor()
        cursor.execute('''SELECT id, barcode_name FROM barcodes''')
        for row in cursor:
            name_list.append(row[1])
            names[row[1]] = row[0]
        sorted_names = sorted(name_list)
        for selection in self.barcode_listbox.GetSelections():
            self.barcode_listbox.Delete(selection)
            cursor.execute('''DELETE FROM barcodes WHERE id = %d'''%names[sorted_names[selection]])
        connection.commit()
        self.Update()

    def sorted_dict_values(self, adict):
        items = adict.items()
        items.sort()
        return [value for key, value in items] #@UnusedVariable
    
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
        cursor.execute('''CREATE TABLE IF NOT EXISTS barcodes (id integer primary key, sequence, barcode_name)''')
        connection.commit()
        
class AddBarcode(wx.Dialog):
    def __init__(self, parent, id, title):
        wx.Dialog.__init__(self, parent, id, title, size=(300, 200))

        vbox = wx.BoxSizer(wx.VERTICAL)
        self.sequence = wx.TextCtrl(self, -1, 'Barcode Sequence', size = (300,25))
        self.bc_name = wx.TextCtrl(self, -1, 'Barcode Name', size = (300,25))
        vbox.AddMany([self.sequence, (self.bc_name, 0, wx.TOP, 10)])
        
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        okButton = wx.Button(self, -1, 'Ok', size=(70, 30))
        okButton.Bind(wx.EVT_BUTTON, self.add)
        closeButton = wx.Button(self, -1, 'Close', size=(70, 30))
        closeButton.Bind(wx.EVT_BUTTON, self.close)
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
        cursor.execute('''SELECT barcode_name FROM barcodes''')
        name_list = []
        for row in cursor:
            name_list.append(row[0])
        if str(self.sequence.GetValue()) == 'Barcode Sequence':
            dlg = wx.MessageDialog(self, '', 'Please enter a nucleotide sequence', wx.OK)
            dlg.ShowModal()
            dlg.Destroy()
        elif str(self.bc_name.GetValue()) == "Barcode Name":
            dlg = wx.MessageDialog(self, '', 'Please enter a unique Barcode Name', wx.OK)
            dlg.ShowModal()
        else:
            if str(self.bc_name.GetValue()) not in name_list:
                t = str(self.sequence.GetValue().upper()), str(self.bc_name.GetValue())
                cursor.execute('''INSERT INTO barcodes (id, sequence, barcode_name) VALUES (NULL, ?, ?)''', t)
                connection.commit()
                self.Close()
            else:
                dlg = wx.MessageDialog(self, '', 'Please enter a UNIQUE Barcode Name', wx.OK)
                dlg.ShowModal()
            
if __name__ == "__main__":
    if not os.path.exists('barcodes.db'):
        db.create_db(db())          
    app = wx.App(False)
    MainFrame(None, -1, 'PREP: Pyrosequence Read Extractor and Processor')
    app.MainLoop()
        


        



    
