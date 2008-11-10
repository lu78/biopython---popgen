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
from Bio.Align.Generic import Alignment
from Bio.Alphabet import IUPAC

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
    <BLANKLINE>
    
    >>> ldhat = LdHatGenerator(nseq = 100, seqlen = 10,
    ...        freqs_per_site = [0.2, 0.5, 0.1, 0.4, 0.3, 0.6, 0.8, 0.1, 1.0, 0.42],
    ...        alleles_per_site = ['AT', 'CT', 'GA', 'GT', 'TG', 'GT', 'TC', 'CA', 'TA', 'GT'])
    >>> print ldhat # this test will always fail, until I'll fix it
    >>> 
    
    """
    
    align = Alignment(IUPAC.IUPACUnambiguousDNA())
                      
    if seed is None:
        random.seed()
    # Add a check for arguments here
    
    if not (len(freqs_per_site) == len(alleles_per_site) == seqlen):
        raise ValueError("note that (len(freqs_per_site) != len(allele_per_site) != seqlen) ")

    
    for n in xrange(nseq):
        seqrecord = SeqRecord(Seq(''), id = 'seq%d' % (n+1), description = '')
        align._records.append(seqrecord)
        
        for pos in xrange(seqlen):
#        logging.debug((pos, alleles_per_site[pos]))          
            
            if random.random() <= freqs_per_site[pos]:
                nt = alleles_per_site[pos][0]
            else:
                nt = alleles_per_site[pos][1]
                
            seqrecord.seq += nt
        logging.debug(seqrecord.seq.tostring())
            
    logging.debug([seqrecord.align._records.tostring() for seqrecord in align._records])
    
    
    # write output. This could be refactored to avoid to repeat the cycle.
    output = '%s %s 1\n' %(nseq, seqlen)
    for seq in align._records:
        output += seq.format('fasta')
    
    return output
    
    
        
            
    
    
    
def _test():
    import doctest
    doctest.testmod()
    
if __name__ == '__main__':
    logging.basicConfig(level = logging.DEBUG)
    _test()
    
    