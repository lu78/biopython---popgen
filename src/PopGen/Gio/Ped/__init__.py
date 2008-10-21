#!/usr/bin/env python
# Copyright 2008 by Giovanni Dall'Olio.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
PED file parser.

The format is described here:
- http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped

Every line in the PED file corresponds to an individual. 

Example of PED file:
 # population id father mother char1 phenotype locus1_allele1 locus1_allele2 locus2_allele1 .....
 Mandenka HGDP00912 0 0 1 2 C C G A G G T T T T A A G G
 Mandenka HGDP01283 0 0 1 2 T C G A G G C T C C G G G A
 Yoruba HGDP00928 0 0 2 2 C C G G G G C T T T G A G G
 Yoruba HGDP00937 0 0 1 2 C C A A A G C C T C A A G A


Classes:

Record            Holds PED file data
RecordParser      Parses a PED record (file) into a Record object.

_Scanner         Scans a PED record.
_RecordConsumer  Consumes PED data to a Record object.


"""


from copy import deepcopy
from types import *

from Bio import File
from Bio.ParserSupport import *     # overwriting previous import?
from PopGen.AbstractPopRecord import AbstractPopRecord 
        


class Record(AbstractPopRecord):
    """Holds information from a PEP record.

    Members:
    marker_len         The marker length (2 or 3 digit code per allele).    
    
    comment_line       Comment line.

    loci_list          List of loci names.
    
    populations        List of population data.
    
    populations has one element per population. Each element is itself
    a list of individuals, each individual is a pair composed by individual
    name and a list of alleles (2 per marker): Example
    [
        [
            ('Ind1', [(1,2),    (3,3), (200,201)],
            ('Ind2', [(2,None), (3,3), (None,None)],
        ],
        [
            ('Other1', [(1,1),  (4,3), (200,200)],
        ]
    ]

    
    """
    
    def __init__(self):
        self.marker_len      = 0
        self.comment_line    = ""
        self.loci_list       = []
        self.populations     = []
        
#    def __str__(self):
#        return 



class RecordParser(AbstractParser):
    """Parses GenePop data into a Record object.

    """
    def __init__(self):
        self._scanner = _Scanner()
        self._consumer = _RecordConsumer()

    def parse(self, handle):
        self._scanner.feed(handle, self._consumer)
        return self._consumer.data
    
def parse(handle):
    """Parses a handle containing a PED file.
    """
    parser = RecordParser()
    return parser.parse(handle)

class _Scanner:
    """Scans a GenePop record.
    
    There is only one record per file.    # no there are multiple records per file?????
    
    """

    def feed(self, handle, consumer):
        """feed(self, handle, consumer)

        Feed in a PED unit record for scanning.  handle is a file-like
        object that contains a Genepop record.  consumer is a
        Consumer object that will receive events as the report is scanned.

        """
        # Check whetherfile exists????
        if isinstance(handle, File.UndoHandle):
            uhandle = handle
        else:
            uhandle = File.UndoHandle(handle)
            
        consumer.start_record()
        
        comment_line = uhandle.readline().rstrip()
        consumer.comment(comment_line)
        
        # Here it goes the code to parse the single PED line
        pass


class _RecordConsumer(AbstractConsumer):
    """Consumer that converts a GenePop record to a Record object.

    Members:
    data    Record with GenePop data.

    """
    def __init__(self):
        self.data = None

    def start_record(self):
        self.data = Record()

    def end_record(self):
        # Here it goes the code to create a Population object and return it to the Record handler????
        pass
    
    def comment(self, comment_line):
        self.data.comment_line = comment_line       # what if there are 2 comment lines or more?


