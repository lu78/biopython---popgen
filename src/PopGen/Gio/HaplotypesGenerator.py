#!/usr/bin/env python
# Copyright 2008 by Giovanni Dall'Olio.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
Parse, writes, and generates ldhat files.

HaplotypesGenerator:
    generates a random haplotypes seq, given:
        - number of sequences
        - lenght of sequences
        - random seed (optional)
        - list of frequencies per site
        - list of alleles per site
        
    this script is meant to produce sets to be used for testing purposes.

"""
import random
import logging
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Generic import Alignment

class NotImplementedException(Exception): pass

class HaplotypesGenerator(object):
    """
    generates a random biopython Alignment object, 
     containing a set of haplotypes generated randomly.
    
    parameters:
    o nseq              -> number of sequences
    o seqlen            -> lenght of sequences
    o freqs_per_site    -> a list indicating the frequency in every position.
    o alleles_per_site  -> a list with the two possible alleles for every site. 
                           only biallelic alleles are supported.
    
#    >>> import HaplotypesGenerator 
    >>> ag = HaplotypesGenerator(nseq = 3, seqlen = 10, 
    ...        freqs_per_site = [1.0] * 10, 
    ...        alleles_per_site = ['AT'] * 10,    
    ...        )
    >>> (ldhat, alignment) = ag.generate()
    >>> print ldhat
    3 10 1
    >seq1 
    AAAAAAAAAA
    >seq2 
    AAAAAAAAAA
    >seq3 
    AAAAAAAAAA
    <BLANKLINE>
    
    >>> ag = HaplotypesGenerator(nseq = 2, seqlen = 10,
    ...        freqs_per_site = [0.2, 0.5, 0.1, 0.4, 0.3, 0.6, 0.8, 0.1, 1.0, 0.42],
    ...        alleles_per_site = ['AT', 'CT', 'GA', 'GT', 'TG', 'GT', 'TC', 'CA', 'TA', 'GT'])
    
    >>> (ldhat, alignment) = ag.generate(seed=10)      # be careful with this test
    >>> print alignment.format('fasta')
    >seq1 
    TCAGGTTATG
    >seq2 
    TTAGGTTATT
    <BLANKLINE>
    
    >>> (ldhat, alignment) = ag.generate(seed=20)
    >>> print alignment.format('fasta')
    >seq1 
    TTATTTCATG
    >seq2 
    TTAGGGTATG
    <BLANKLINE>
    
    """
    
    def __init__(self, nseq, seqlen, freqs_per_site, alleles_per_site, alphabet = None):
        """
        Initialize an AlignmentGenerator object.
    
        defaults:
        alphabet -> IUPACUnambiguousDNA
        """
                          
        # Check for arguments
        if not (len(freqs_per_site) == len(alleles_per_site) == seqlen):
            raise ValueError("note that (len(freqs_per_site) != len(allele_per_site) != seqlen) ")
    
        self.nseq = nseq
        self.seqlen = seqlen
        self.freqs_per_site = freqs_per_site
        self.alleles_per_site = alleles_per_site
        
        if alphabet is None:
            from Bio.Alphabet import IUPAC
            self.alphabet = IUPAC.IUPACUnambiguousDNA()
        
        
    def generate(self, seed = None):
        """
        Generate a random set with the given parameters.
        Returns an Alignment object.
        
        """
        alignment = Alignment(self.alphabet)
        
        # initialize random generator
        if seed is not None:
            random.seed(seed)

        records = alignment._records    # hack to add SeqRecord objects to Alignment until # fixed
        
        for n in xrange(self.nseq):
            seqrecord = SeqRecord(Seq(''), id = 'seq%d' % (n+1), description = '')
            records.append(seqrecord)
            
            for pos in xrange(self.seqlen):
    #        logging.debug((pos, alleles_per_site[pos]))          
                
                if random.random() <= self.freqs_per_site[pos]:
                    nt = self.alleles_per_site[pos][0]
                else:
                    nt = self.alleles_per_site[pos][1]
                    
                seqrecord.seq += nt
            logging.debug(seqrecord.seq.tostring())
                
        logging.debug([seqrecord.seq.tostring() for seqrecord in records])
    
    
        # write output. It could use the alignment._format function when it will
        # support ldhat files.
        output = '%s %s 1\n' %(self.nseq, self.seqlen)
        for seq in records:
            output += seq.format('fasta')
        
        return output, alignment
        

def paramsGenerator(mode = None, seqlen = 20, nseq = 10):
    """Create values to be used with LDHAT Generator
    For example, if you want a set in which all the nucleotides have a 
    frequency of '1.0', you can use this function with the parameter
    mode = 'equals'.
    
    Possible values for mode:
    o equals
    o onehotspot
    o example1
    
    Note: 
        sometimes the parameters seqlen, nseq are ignored.
    
    examples:
    >>> params = paramsGenerator('equals', 3, 4)
    >>> print params
    [4, 3, [1.0, 1.0, 1.0], ['AG', 'AG', 'AG']]
    
    # these parameters could be used as inputs to LdHatGenerator: 
    >>> ag = HaplotypesGenerator(params[0], params[1], params[2], params[3])
    >>> print ag.generate()[1]
    >seq1 
    AAA
    >seq2 
    AAA
    >seq3 
    AAA
    >seq4 
    AAA
    <BLANKLINE>
    
    >>> params = paramsGenerator('onehotspot', 3, 4)
    >>> print params[2] == [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
    ... 0.8, 0.8, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
    True
    
    """
    if mode == 'equals':
        freqs = [1.0] * seqlen
        alleles = ['AG'] * seqlen
    elif mode == 'onehotspot':  # to implement
        freqs = [0.5] * 23
        freqs[10] = freqs[11] = 0.8
        alleles = ['AG'] * len(freqs)
    elif mode == 'example1':
        freqs = [0.2, 0.5, 0.1, 0.4, 0.3, 0.6, 0.8, 0.1, 1.0, 0.42]
        alleles = ['AT', 'CT', 'GA', 'GT', 'TG', 'GT', 'TC', 'CA', 'TA', 'GT']
    elif mode == 'example2':
        freqs = [0.4, 0.1, 0.1, 0.1, 0.4, 0.2, 0.5, 0.7, 0.4, 0.3]
        alleles = ['AT', 'CT', 'GA', 'GT', 'TG', 'GT', 'TC', 'CA', 'TA', 'GT']
    elif mode == 'all0.5':
        freqs = [0.5] * seqlen
        alleles = ['AT'] * seqlen
    elif mode == 'all0.3':
        freqs = [0.3] * seqlen
        alleles = ['AT'] * seqlen
    else:
        raise NotImplementedException
        
    output = [nseq, seqlen, freqs, alleles] 
    return output
        
            
    
    
    
def _test():
    import doctest
    doctest.testmod()
    
if __name__ == '__main__':
#    logging.basicConfig(level = logging.DEBUG)
    logging.basicConfig(level = None)
    _test()
    
    