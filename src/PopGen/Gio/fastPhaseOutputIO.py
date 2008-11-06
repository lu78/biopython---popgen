#!/usr/bin/env python
# Copyright 2008 by Giovanni Dall'Olio.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
fastPHASE output file parser.
This module contains only a 'fastPhaseOutputIterator', because you are not 
supposed to write PHASE output files  

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
... /usr/bin/fastPHASE -T10 -K20 -p -u sample.popinfo -o sample sample.inp 
... END COMMAND_LINE
... 
... BEGIN COMMAND_EXPLAIN
...  K no. clusters (chosen or supplied): 20
...  S seed for random numbers (chosen or supplied): 1224252023
... END COMMAND_EXPLAIN
... 
... BEGIN DESCRIBE_TASKS
... minimize switch error
... END DESCRIBE_TASKS
... 
... BEGIN GENOTYPES
... Ind1  # subpop. label: 6  (internally 1)
... T T T T T G A A A C C A A A G A C G C T G C G T C A G C C T G C A A T C T G
... T T T T T G C C C C C A A A A G C G C G T C G T C A G T C T A A G A C C T A
... Ind2  # subpop. label: 6  (internally 1)
... C T T T T G C C C T C A A A A G T G C T G T G C C A G T C T A C G G C C T G
... T T T T T G A A A C C A A A G A C G C T T C G T C A G T A T A C G A T C T A
... END GENOTYPES
... ''')  
>>> for record in fastPhaseOutputIterator(phasefile):
...    print record.seq
TTTTTGAAACCAAAGACGCTGCGTCAGCCTGCAATCTG
TTTTTGCCCCCAAAAGCGCGTCGTCAGTCTAAGACCTA
CTTTTGCCCTCAAAAGTGCTGTGCCAGTCTACGGCCTG
TTTTTGAAACCAAAGACGCTTCGTCAGTATACGATCTA
"""

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from PopGenExceptions import InvalidInputFile
import re 
import logging

def fastPhaseOutputIterator(handle):
    """
    Iterates over a Phase file output handler
    Returns a SeqRecord object.
    """
    while True:
        line = handle.readline()
        if line == "": 
            return   # premature end of file or just empty
        if line.startswith("BEGIN GENOTYPES"):
            break
    
    while True:
        line = handle.readline().strip()
        if line == "END GENOTYPES": 
            return      # exit cycle, file has been parsed
        if line.strip() == "": 
            break
        
        descr = line.split('#')
        id1 = descr[0].strip() + '_all1'
        id2 = descr[0].strip() + '_all2'
        name1 = id1
        name2 = id2
        
        # TO FIX: if there are blank lines in the file (there shouldn't), 
        # an error exception should be thrown.
        seqs = []
        while True:
            line = handle.readline()    # re-defining line var, but it doesn't matter               
            
            if re.match('^(\w\s)+', line):  
                seq = line.replace(" ", "").replace("\r", "").strip()
                
                seqs.append(seq)
                if len(seqs) == 2:
                    break
                
            elif not line or re.match('^\w+', line):
                # the file has been parsed, but no other sequence has been found.
                raise InvalidInputFile("Missing sequence in input file")
            
#        logging.debug(seqs)
        
        # check that len(seq1) == len(seq2)
        if len(seqs[0]) != len(seqs[1]):
            raise InvalidInputFile("Two chromosomes with different length")
                    
        yield SeqRecord(Seq(seqs[0]), id = id1, name = name1, description = descr[1]) 
        yield SeqRecord(Seq(seqs[1]), id = id2, name = name2, description = descr[1])
#        print "line", line
        
    assert False, "should not reach this line"
                                            

def _test():
    """tests current module"""
    import doctest
    logging.basicConfig(level=logging.DEBUG)
    doctest.testmod()
    testfilesdir = '../../Tests/SeqIO/'
    testfilesdir = '.'
    doctest.testfile(testfilesdir + 'test_fastPhaseOutputIO.py')
    
if __name__ == '__main__':
    _test()
