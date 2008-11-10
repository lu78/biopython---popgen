#!/usr/bin/env python
# Copyright 2008 by Giovanni Dall'Olio.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
fastPHASE is software for haplotype reconstruction and missing genotype
estimation from population genetic SNP data.

Scheet, P and Stephens, M (2006) "A fast and flexible statistical model for
large-scale population genotype data: applications to inferring missing
genotypes and haplotypic phase." Am J Hum Genet 78(4):629-44.


This module contains only a 'fastphaseoutputIterator', because you are not 
supposed to write PHASE output files  

The format is described here:
- http://stephenslab.uchicago.edu/software.html


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
>>> for record in fastphaseoutputIterator(phasefile):
...    print record
IUPACUnambiguousDNA() alignment with 4 rows and 38 columns
TTTTTGAAACCAAAGACGCTGCGTCAGCCTGCAATCTG Ind1_all1
TTTTTGCCCCCAAAAGCGCGTCGTCAGTCTAAGACCTA Ind1_all2
CTTTTGCCCTCAAAAGTGCTGTGCCAGTCTACGGCCTG Ind2_all1
TTTTTGAAACCAAAGACGCTTCGTCAGTATACGATCTA Ind2_all2
"""

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from PopGenExceptions import InvalidInputFile
import re 
import logging
from Bio.Align.Generic import Alignment
from Bio.Alphabet import IUPAC

def fastphaseoutputIterator(handle, alphabet = None):
    """
    Iterates over a fastPhase file output handler
    Returns an Alignment object.
    
    inputs:
    o handle:     a file handler for a fastPhase output file
    o alphabet:   a biopython alphabet object(if None, UnambiguousDNa is used)
      
    """
    if alphabet is None:
        alphabet = IUPAC.IUPACUnambiguousDNA()
        
    align = Alignment(alphabet)
    records = align._records    # hack to append SeqRecord objects to Alignment
                                # see bugs 2553/2554
    
    while True:
        line = handle.readline()
        if line == "": 
            return   # premature end of file or just empty
        if line.startswith("BEGIN GENOTYPES"):
            break
    
    while True:
        line = handle.readline().strip()
        if line == "END GENOTYPES": 
            yield align
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
                   
        seq1 = SeqRecord(Seq(seqs[0]), id = id1, name = name1, description = descr[1])
        seq2 = SeqRecord(Seq(seqs[1]), id = id2, name = name2, description = descr[1])
        records.append(seq1)
        records.append(seq2)
#        print "line", line
        
    assert False, "should not reach this line"
                                            

def _test():
    """tests current module"""
    import doctest
    logging.basicConfig(level=logging.DEBUG)
    doctest.testmod()
    testfilesdir = '../../Tests/SeqIO/'
    testfilesdir = './'
    doctest.testfile(testfilesdir + 'test_fastPhaseOutputIO.py')
    
if __name__ == '__main__':
    _test()
