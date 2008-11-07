#!/usr/bin/env python
# Copyright 2008 by Giovanni Dall'Olio.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
Parse, writes, and generates ldhat files.

LdHatGenerator:
    generates a random ldhat imput file, given:
        - number of sequences
        - lenght of sequences
        - random seed (optional)
        - list of frequencies per site
        - list of alleles per site
    the output sets are supposed to be used for testing purposes.
    e.g.: you want an input ldhat file in which  

"""
import random
import logging
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def LdHatGenerator(nseq, seqlen, freqs_per_site, alleles_per_site, seed = None):
    """
    generates a random ldhat input file (sites.txt)
    
#    >>> import LdHatGenerator
    >>> ldhat = LdHatGenerator(nseq = 3, seqlen = 4, 
    ...        freqs_per_site = [1.0, 1.0, 1.0, 1.0], 
    ...        alleles_per_site = ['AT', 'CT', 'GA', 'GT']    
    ...        )
    >>> print ldhat
    3 4 1
    >seq1
    ACGG
    >seq2
    ACGG
    >seq3
    ACGG
    
#    >>> ldhat = LdHatGenerator(nseq = 3, seqlen = 4,
#    ...        freqs_per_site = [0.2, 0.5, 1.0, 1.0],
#    ...        alleles_per_site = ['AT', 'CT', 'GA', 'GT'])
#    >>> print ldhat
#    
    """
    if seed is None:
        random.seed()
    # Add a check for arguments here
    
    seqs = [SeqRecord('')] * seqlen
    
    for seq in seqs:
        for pos in xrange(seqlen):
#        logging.debug((pos, alleles_per_site[pos]))          
            
            if random.random() <= freqs_per_site:
                nt = alleles_per_site[pos][0]
            else:
                nt = alleles_per_site[pos][1]
            logging.debug(nt)
            seq.seq += nt
            
    logging.debug([seq.seq.tostring() for seq in seqs])
        
            
    
    
    
def _test():
    import doctest
    doctest.testmod()
    
if __name__ == '__main__':
    logging.basicConfig(level = logging.DEBUG)
    _test()
    
    