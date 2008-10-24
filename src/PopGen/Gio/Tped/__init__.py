#!/usr/bin/env python

"""
Parses a Tped file.

>>> import tempfile
>>> sample_Tped_file = TemporaryFile()
>>> sample_Tped_file.write(\
... '''# Chromosome SNP_Id ??? position locus1_allele1 locus1_allele2 locus2_allele1
... 4 rs10000543 0 30979886 C C T C C C C C C C C C T C C C C C C C C C T C C C 
... 4 rs10000929 0 131516474 A A A A G A G A A A G A G A G A A A G A G A G A G G 
... 4 rs10002472 0 159087423 A G G G G G A A G G G G A G G G A G A G G G G G G G 
... 4 rs10001548 0 166098831 T T C T C C C T C T C C T T C T C C C C T T C T C T 
... 4 rs10001378 0 182579995 C C T T C C T T T C T T T C T T T T T T T T C C T T
... 4 rs10004399 0 183794360 A A A A A A A A G A G A G A G G G A A A A A G G G A
... 4 rs10000918 0 186505570 G A G G G G G G A A A A G G G A 0 0 G G G G G A G G
... '''
>>> sample_Tped_file.seek(0)

"""


def TpedIterator(handle):
    """
    Iterates on an TPed file handler.
    Returns Marker objects.
    
    >>> ti = TPedIterator(sample_Tped_file)
    >>> for marker in ti:
    ...     print marker
    """
    while True:
        line = handle.readline()
        if line.strip() == "":
            return # premature end of file?
        if line.startswith('#'):
            comment = line