# Copyright 2008 by Tiago Antao.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
This module provides code to work with FStat.

Classes:
Record           Holds FStat data.
RecordParser     Parses a FStat record (file) into a Record object.

_Scanner         Scans a FStat record.
_RecordConsumer  Consumes FStat data to a Record object.


"""
from copy import deepcopy
from types import *

from Bio import File
from Bio.ParserSupport import *
from Bio.PopGen import GenePop


class Record:
    """Holds information from a FStat record.

    Members:
    sub_pops           Number of sub populations.

    num_loci           Number of loci.

    marker_len         Maximum number of alleles (normally 99).    

    ploidy             Ploidy (Currently we only support diploid data).
    
    loci_list          List of loci names.
    
    populations        List of population data.
    
    populations has one element per population. Each element is itself
    a list of individuals, each individual is a
    list of alleles (2 per marker): Example
    [
        [
            [(1,2),    (3,3), (200,201)],
            [(2,None), (3,3), (None,None)],
        ],
        [
            [(1,1),  (4,3), (200,200)],
        ]
    ]

    
    """
    def __init__(self):
        self.sub_pops        = 0
        self.num_loci        = 0
        self.marker_len      = 0
        self.ploidy          = 2
        self.loci_list       = []
        self.populations     = []

    def __str__(self):
        rep  = [str(self.sub_pops) + ' ' + str(self.num_loci) + ' ' +
                str(self.marker_len) + ' ' + str(self.ploidy) + '\n']
        rep.append('\n'.join(self.loci_list) + '\n')
        for pop_i in range(len(self.populations)):
            for indiv in self.populations[pop_i]:
                rep.append(str(pop_i+1) + ' ')
                for marker in indiv:
                    for al in marker:
                        if al == None:
                            al = '0'
                        aStr = str(al)
                        while len(aStr)<len(str(self.marker_len)):
                            aStr = "".join(['0', aStr])
                        rep.append(aStr)
                    rep.append(' ')
                rep.append('\n')
        return "".join(rep)

    def convert_to_genepop(self):
        '''Converts to GenePop format.
        '''
        gp = GenePop.Record()
        gp.marker_len = len(str(self.marker_len))
        gp.comment_line = 'Converted from FStat format'
        gp.loci_list = deepcopy(self.loci_list)
        gp.populations = []
        for mypop in self.populations:
            indIdx = 1
            pop = []
            for indiv in mypop:
                pop.append((str(indIdx), deepcopy(indiv)))
                indIdx += 1
            gp.populations.append(pop)
        return gp
    

class RecordParser(AbstractParser):
    """Parses FStat data into a Record object.

    """
    def __init__(self):
        self._scanner = _Scanner()
        self._consumer = _RecordConsumer()

    def parse(self, handle):
        self._scanner.feed(handle, self._consumer)
        return self._consumer.data

def parse(handle):
   """Parses a handle containing a FStat file.
   """
   parser = RecordParser()
   return parser.parse(handle)

class _Scanner:
    """Scans a FStat record.
    
    There is only one record per file.
    """

    def feed(self, handle, consumer):
        """feed(self, handle, consumer)

        Feed in a FStat unit record for scanning.  handle is a file-like
        object that contains a FStat record.  consumer is a
        Consumer object that will receive events as the report is scanned.

        """
        if isinstance(handle, File.UndoHandle):
            uhandle = handle
        else:
            uhandle = File.UndoHandle(handle)


        consumer.start_record()
        
        boot_line = uhandle.readline().rstrip()
        pops, num_loci, marker_len, ploidy = map(
                lambda x: int(x), boot_line.split(' '))
        consumer.sub_pops(pops)
        consumer.num_loci(num_loci)
        consumer.marker_len(marker_len)
        consumer.ploidy(ploidy)

        for i in range(num_loci): 
            loci_name = uhandle.readline().rstrip()
            consumer.loci_name(loci_name)

        line = uhandle.readline()
        while line<>'':
            line = line.rstrip()
            toks = line.split(' ')
            consumer.individual(toks[0], toks[1:])
            line = uhandle.readline()
        consumer.end_record()

class _RecordConsumer(AbstractConsumer):
    """Consumer that converts a FStat record to a Record object.

    Members:
    data    Record with FStat data.

    """
    def __init__(self):
        self.data = None

    def start_record(self):
        self.data = Record()

    def end_record(self):
        pass

    def loci_name(self, locus):
        self.data.loci_list.append(locus)

    def sub_pops(self, sub_pops):
        self.data.populations = []
        for i in range(sub_pops):
            self.data.populations.append([])
        self.data.sub_pops = sub_pops

    def num_loci(self, num_loci):
        self.data.num_loci = num_loci

    def marker_len(self, marker_len):
        self.data.marker_len = marker_len

    def ploidy(self, ploidy):
        self.data.ploidy = ploidy

    def individual(self, pop, allele_list):
        alleles = []
        for a in allele_list:
            alleles.append((int(a[:len(a)/2]),int(a[len(a)/2:])))
        self.data.populations[int(pop)-1].append(alleles)
    

