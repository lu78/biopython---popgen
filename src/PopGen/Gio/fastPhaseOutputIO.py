#!/usr/bin/env python
# Copyright 2008 by Giovanni Dall'Olio.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
fastPHASE output file parser.
This module contains only a 'fastPhaseOutputIterator', because you are not supposed to write PHASE output files  

The format is described here:

Example of PHASE output file:
>>> from StringIO import StringIO
>>> phasefile = StringIO('''
... ********************************************
... *                                          *
... *      Output from fastPHASE 1.2.3         *
... *      Code by P Scheet                    *
... *                                          *
... ********************************************
... BEGIN COMMAND_LINE
... /usr/bin/fastPHASE -T10 -K2 -p -usamplefile-fastphase.popinfo -osamplefile samplefile.inp 
... END COMMAND_LINE
... 
... BEGIN COMMAND_EXPLAIN
...  K no. clusters (chosen or supplied): 2
...  S seed for random numbers (chosen or supplied): 1224252023
... END COMMAND_EXPLAIN
... 
... BEGIN DESCRIBE_TASKS
... minimize switch error
... END DESCRIBE_TASKS
... 
... BEGIN GENOTYPES
... Ind1  # subpop. label: 6  (internally 1)
... T T T T T G A A A C C A A A G A C G C T G C G T C A G C C T G C A A T C T G T G T T A A G A C T C G
... T T T T T G C C C C C A A A A G C G C G T C G T C A G T C T A A G A C C T A T G C T A A G G C T T G
... Ind2  # subpop. label: 6  (internally 1)
... C T T T T G C C C T C A A A A G T G C T G T G C C A G T C T A C G G C C T G C A T T A A G A T T C G
... T T T T T G A A A C C A A A G A C G C T T C G T C A G T A T A C G A T C T A T G C T A A T G C T T G
... END GENOTYPES
... ''')  
>>> for seq1 in fastPhaseOutputIterator(phasefile):
...     print seq1.seq
TTTTTGAAACCAAAGACGCTGCGTCAGCCTGCAATCTGTGTTAAGACTCG
TTTTTGCCCCCAAAAGCGCGTCGTCAGTCTAAGACCTATGCTAAGGCTTG
CTTTTGCCCTCAAAAGTGCTGTGCCAGTCTACGGCCTGCATTAAGATTCG
TTTTTGAAACCAAAGACGCTTCGTCAGTATACGATCTATGCTAATGCTTG
"""

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def fastPhaseOutputIterator(handle):
    """
    Iterates over a Phase file output handler
    Returns a SeqRecord object.
    """
    while True:
        line = handle.readline()
        if line == "": return   # premature end of file or just empty
        if line.startswith('BEGIN GENOTYPES'):
            break
    
    while True:
        line = handle.readline().strip()
        if not line: break
        if line == "END GENOTYPES": return
        descr = line.strip().split('#')
        id1 = descr[0].strip() + '_all1'
        id2 = descr[0].strip() + '_all2'
        name1 = id1
        name2 = id2
        seq1 = handle.readline().replace(" ", "").replace("\r", "").strip()
        seq2 = handle.readline().replace(" ", "").replace("\r", "").strip()
        yield SeqRecord(Seq(seq1), id = id1, name = name1, description = descr) 
        yield SeqRecord(Seq(seq2), id = id2, name = name2, description = descr)
#        print "line", line
        
    assert False, "should not reach this line"
                                            

def _test():
    """tests current module"""
    import doctest
    doctest.testmod()
    
if __name__ == '__main__':
    _test()
