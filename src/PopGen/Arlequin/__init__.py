#2007 by Tiago Antao <tiagoantao@gmail.com>.  All rights reserved.

"""
This module provides code to work with Arlequin3.


See http://cmpg.unibe.ch/software/arlequin3/ .

Classes:
Record           Holds Arlequin data (arp).
Profile          Holds Arlequin profile data (part of Record)
RecordParser     Parses an Arlequin record (arp file) into a Record object.
_Scanner         Scans a Arlequin record.
_RecordConsumer  Consumes Arlequin data to a Record object.


"""

import re
from PopGen.Common import Marker

class Record:
    """Holds information from a Arlequin (arp) record.

    Members:
    profile   Meta information about the file, this includes:
              (This info is in the Profile class)
      title         Title
      num_pops      Number of sampled populations
      is_genotypic  Genotypic data?
      gam_phase     Gametic Phase
      is_recessive  Recessive data?
      data_type     Data Type
      locus_sep     Locus Separator
      missing_tok   Missing data token
      markers_type  Types of markers (calculated)
    pop_data        A list with pop data, each element (tuple) has a pop:
      pop_name      Population name
      samp_size     Sample size
      data          A list with individuals, composed of (tuple)
        name        Name
        number      A number (unknown function)
        alleles     A list or pair of lists (diploidy case) of alleles

        STRUCTURE INFORMATION IS MISSING!
    """
    def __init__(self):
        self.profile = Profile()
        self.pop_data = []

    #incomplete
    def __str__(self):
        rep = ''
        rep = self.profile.title + '\n'
        rep = str(self.profile.markers_type) + '\n'
        rep += str(self.pop_data)
        return rep

class Profile:
    """Profile of an Arlequin file.
    
    For details please see the Record class.
    """
    def __init__(self):
        self.title        = ''
        self.num_pops     = 0
        self.is_genotypic = True
        self.gam_phase    = 1
        self.is_recessive = False
        self.data_type    = 'STANDARD'
        self.locus_sep    = 'WHITESPACE'
        self.missing_tok  = '?'
        self.markers_type = []
        
class RecordParser:
    """Parses Arlequin data into a Record object.

    """
    def __init__(self):
        self._scanner = _Scanner()
        self._consumer = _RecordConsumer()

    def parse(self, handle):
        self._scanner.feed(handle, self._consumer)
        return self._consumer.data

class _Scanner:
    """Scans a Arlequin record.
    
    There is only one record per file.
    
    """
    
    def report_param(self, name, line, rep_fun, prep_fun, is_str = False):
        if is_str:
            rexp = re.compile('^[^#]*' + name + '="(.*)"')
        else:
            rexp = re.compile('^[^#]*' + name + '= *([^# ]+)')
        m = rexp.search(line)
        if m <> None:
            rep_fun(prep_fun(m.group(1)))

    def report_allele(self, token, consumer):
        token = token.replace(' ', '')
        #this is incomplete
        try:
            msat = int(token)
            consumer.MSAT(msat)
        except ValueError:
            if len(token) == 1:
                consumer.SNP(token)
            else:
                consumer.seq(token)

    def feed(self, handle, consumer):
        """feed(self, handle, consumer)

        Feed in a FDist unit record for scanning.  handle is a file-like
        object that contains a FDist record.  consumer is a
        Consumer object that will receive events as the report is scanned.

        """

        consumer.start_record()
        l = handle.readline()
        on_profile = False
        on_data = False
        on_structure = False
        while l<>'':
            l = l.rstrip()
            if l.find('[Profile]') <> -1: #Should control for comment
                on_profile   = True
                on_data      = False
                on_structure = False
            if l.find('[Data]') <> -1: #Should control for comment
                on_profile   = False
                on_data      = True
                on_structure = False
            if l.find('[Structure]') <> -1: #Should control for comment
                on_profile   = False
                on_data      = False
                on_structure = True
            if on_profile:
                self.report_param('Title', l, consumer.title, str, True)
                #if l.find('NbSamples'): consumer.()
                #if l.find('GameticPhase'): consumer.()
                #if l.find('RecessiveData'): consumer.()
                #if l.find('DataType'): consumer.()
                #if l.find('LocusSeparator'): consumer.()
                self.report_param('GenotypicData',l,consumer.ploidy,int,False)
            if on_data:
                self.report_param('SampleName',l, consumer.pop_name,str,True)
                self.report_param('SampleSize',l, consumer.pop_size,int,False)
                tokens = l.split('\t')
                if tokens[0].find('_') <> -1:
                    pop_i, indiv_name = tokens[0].split('_')
                    consumer.new_indiv(indiv_name)
                    consumer.new_chromatid()
                    #skipping tokens[1] - the told unk number
                    for tok in tokens[2:]:
                        self.report_allele(tok, consumer) 
                    consumer.end_chromatid()
                    if consumer.data.is_genotypic:
                        l = handle.readline().rstrip()
                        consumer.new_chromatid()
                        tokens = l.split('\t')
                        for tok in tokens[2:]:
                            self.report_allele(tok, consumer) 
                        consumer.end_chromatid()
                    consumer.end_indiv()
                elif l.find('}') <> -1:
                    consumer.end_pop()
            l = handle.readline()
        consumer.end_record()
        
class _RecordConsumer:
    """Consumer that converts a Arlequin record to a Record object.

    Members:
    data    Record with Arlequin data.

    """
    def __init__(self):
        self.data = None

    def start_record(self):
        self.data = Record()
        self._detecting_markers = True
        self._markers = []

    def start_profile(self):
        pass
        
    def start_pop_data(self):
        pass

    def end_record(self):
        pass

    def ploidy(self, ploidy):
        if ploidy == 1:
            self.data.is_genotypic = True
        else:
            self.data.is_genotypic = False

    def pop_name(self, pop_name):
        self._current_pop = [pop_name, None, []]

    def end_pop(self):
        self.data.pop_data.append(tuple(self._current_pop))

    def pop_size(self, pop_size):
        self._current_pop[1] = pop_size

    def new_indiv(self, name):
        self._current_indiv = [name, None, []]
        self._chrs = []

    def end_indiv(self):
        if self.data.profile.is_genotypic:
            self._current_indiv[2] = tuple(self._chrs)
        else:
            self._current_indiv[2] = self._chrs[0]
        self._current_pop[2].append(tuple(self._current_indiv))

    def new_chromatid(self):
        self._chrs.append([])

    def end_chromatid(self):
        if self._detecting_markers:
            self.data.profile.markers_type = self._markers
        self._detecting_markers = False

    def SNP(self, SNP):
        if self._detecting_markers:
            self._markers.append(Marker.SNP)
        self._chrs[-1].append(SNP)

    def MSAT(self, MSAT):
        if self._detecting_markers:
            self._markers.append(Marker.MSAT)
        self._chrs[-1].append(Marker.MSAT)

    def seq(self, seq):
        if self._detecting_markers:
            self._markers.append('seq')
        self._chrs[-1].append(Marker.SEQ)

    def title(self, title):
        self.data.profile.title = title
