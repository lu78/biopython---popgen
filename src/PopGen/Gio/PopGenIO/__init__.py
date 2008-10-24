#!/usr/bin/env python
# Copyright 2008 by Giovanni Dall'Olio.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
PopGenIO parser

Generic parser for population genetics, especially SNP data, like GenePop, Ped and Tped formats.

>>> from PopGen.Gio import PopGenIO
>>> import tempfile
>>> sample_tped_file = tempfile.TemporaryFile()
>>> sample_tped_file.write('''
... # Chromosome SNP_Id ??? position locus1_allele1 locus1_allele2 locus2_allele1
... 4 rs10000543 0 30979886 C C T C C C C C C C C C T C C C C C C C C C T C
... 4 rs10000929 0 131516474 A A A A G A G A A A G A G A G A A A G A G A G A
... 4 rs10002472 0 159087423 A G G G G G A A G G G G A G G G A G A G G G G G
... 4 rs10001548 0 166098831 T T C T C C C T C T C C T T C T C C C C T T C T 
... 4 rs10001378 0 182579995 C C T T C C T T T C T T T C T T T T T T T T C C
... 4 rs10004399 0 183794360 A A A A A A A A G A G A G A G G G A A A A A G G
... 4 rs10000918 0 186505570 G A G G G G G G A A A A G G G A 0 0 G G G G G A
... ''')
>>> sample_tped_file.seek(0)
>>> pi = PopGenIO.parse(sample_tped_file)
>>> for marker in pi:
...    print marker


"""


from PopGen.Gio import TpedIO
from PopGen.Gio import PedIO
import logging

_FormatToIterator = {
                     # for every format, define which iterator should be used, and which kind of object is returned
                     'tped' : (TpedIO.TpedIterator, 'Marker'),
                     'ped'  : (PedIO.PedIterator, 'Individual')
                     }        
_FormatToWriter = {}

def parse(handle, format = 'tped') :
    """Turns a TPED file into an iterator returning Marker objects.

    handle   - handle to the file.
#    format   - lower case string describing the file format.
    
    If you have the file name in a string 'filename', use:

#    >>> from PopGen.Gio import Tped
#    >>> my_iterator = Tped.parse(open(filename,"rU"))

    Note that file will be parsed with default settings. For more control, you
    must use the format specific iterator directly... (no)
    """
    
    #Try and give helpful error messages:
    if isinstance(handle, basestring) :
        raise TypeError("Need a file handle, not a string (i.e. not a filename)")
    if not isinstance(format, basestring) :
        raise TypeError("Need a string for the file format (lower case)")
    if not format :
        raise ValueError("Format required (lower case string)")
    if format != format.lower() :
        raise ValueError("Format string '%s' should be lower case" % format)
    
    #Map the file format to a sequence iterator:    
    if format in _FormatToIterator :
        iterator_generator = _FormatToIterator[format][0]
        return iterator_generator(handle)        
    else :
        raise ValueError("Unknown format '%s'" % format)
    

def _test():
    """Test using doctest"""
    import doctest
    doctest.testmod(verbose=False)
    
if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    _test()
