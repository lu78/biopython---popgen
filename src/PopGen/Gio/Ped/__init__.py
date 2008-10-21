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
from logging import debug

from Bio import File
from Bio.ParserSupport import *     # overwriting previous import?

#from PopGen.AbstractPopRecord import AbstractPopRecord 
from PopGen.GenePop import Record       # Use the same Record object from GenePop.


class RecordParser(AbstractParser):
    """Parses GenePop data into a Record object.

    >>> import tempfile
    >>> sample_PED_file = tempfile.TemporaryFile()
    >>> sample_PED_file.write('''
    ... # population id father mother char1 phenotype locus1_allele1 locus1_allele2 locus2_allele1 .....
    ... Mandenka HGDP00912 0 0 1 2 C C G A G G T T T T A A G G
    ... Mandenka HGDP01283 0 0 1 2 T C G A G G C T C C G G G A
    ... 
    ... Yoruba HGDP00928 0 0 2 2 C C G G G G C T T T G A G G
    ... Yoruba HGDP00937 0 0 1 2 C C A A A G C C T C A A G A
    ... ''')
    >>> sample_PED_file.seek(0)
    >>> rp = RecordParser()
    >>> pops = rp.parse(sample_PED_file)
    
    # >>> pops
    
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
    
    There is only one record per file.
    
    """

    def feed(self, handle, consumer):
        """feed(self, handle, consumer)

        Feed in a PED unit record for scanning.  handle is a file-like
        object that contains a Genepop record.  consumer is a
        Consumer object that will receive events as the report is scanned.

        """
        # uhandle is the handler for the current file. I have not clear what UndoHandle does.
        if isinstance(handle, File.UndoHandle):
            uhandle = handle
        else:
            uhandle = File.UndoHandle(handle)
            
        consumer.start_record()
        
        
        for line in uhandle.readlines():
            # Start parsing the PED file, line per line!
            line = line.strip()    
            if line.startswith('#'):
                # line starting with # are comments
                consumer.comment(line)
#               debug("comment ", line)
            elif line != '':
                # parse a valid PED line and put its content in Record.
                ped_fields = line.strip().split()       # not sure strip is needed
                print consumer.data.populations # ops! I need to change Record object and transfomr populations in a dictionary
        


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

def _test():
    import doctest
    doctest.testmod()
    
if __name__ == '__main__':
    _test()