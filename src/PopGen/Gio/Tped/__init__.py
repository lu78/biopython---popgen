#!/usr/bin/env python

"""
Parses a Tped file.

>>> import tempfile
>>> sample_Tped_file = tempfile.TemporaryFile()
>>> sample_Tped_file.write('''
... # Chromosome SNP_Id ??? position locus1_allele1 locus1_allele2 locus2_allele1
... 4 rs10000543 0 30979886 C C T C C C C C C C C C T C C C C C C C C C T C C
... 4 rs10000929 0 131516474 A A A A G A G A A A G A G A G A A A G A G A G A G 
... 4 rs10002472 0 159087423 A G G G G G A A G G G G A G G G A G A G G G G G G 
... 4 rs10001548 0 166098831 T T C T C C C T C T C C T T C T C C C C T T C T C 
... 4 rs10001378 0 182579995 C C T T C C T T T C T T T C T T T T T T T T C C T
... 4 rs10004399 0 183794360 A A A A A A A A G A G A G A G G G A A A A A G G G
... 4 rs10000918 0 186505570 G A G G G G G G A A A A G G G A 0 0 G G G G G A G
... ''')
>>> sample_Tped_file.seek(0)

"""

from PopGen.Gio.Marker import Marker
import logging
from PopGen.Gio.PopGenExceptions import *


def TpedIterator(handle):
    """
    Iterates on an TPed file handler.
    Returns Marker objects.
    
    >>> import tempfile
    >>> sample_Tped_file = tempfile.TemporaryFile()
    >>> sample_Tped_file.write('''# Chromosome Marker_Id ??? position locus1_allele1 locus1_allele2 locus2_allele1
    ... 4 rs10000543 0 30979886 C C T C 
    ... 4 rs10000929 0 131516474 A A A A
    ... 4 rs10002472 0 159087423 A G G G
    ... ''')
    >>> sample_Tped_file.seek(0)
    >>> ti = TpedIterator(sample_Tped_file)
    >>> for marker in ti:
    ...     print marker
    """
    for line in handle:
#        if line.strip() == "":
#            debug('premature end of file?')
#            return
        if line.startswith('#'):
            comment = line
        else:
            # should add a check for line length and syntax here.
            tped_fields = line.split()
            if len(tped_fields) < 4:
                raise InvalidInputFile('line too short - check input file\n(%s)' % line)
            # should check that current line has the same length than the previous
            # should check that characters after [4] are rigth
            
            chromosome = tped_fields[0]
            marker_name = tped_fields[1]
            unknown = tped_fields[2]
            position = tped_fields[3]
            current_marker = Marker(marker_name)
            current_marker.position = "%s - chrom %s" % (position, chromosome)
            logging.debug(current_marker)

            genotypes = [(tped_fields[i], tped_fields[i+1]) for i in xrange(4, len(tped_fields), 2)]
            logging.debug(genotypes)
            current_marker.genotypes = genotypes
            
            yield current_marker

_FormatToIterator = {
                     'tped' : TpedIterator,
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
        iterator_generator = _FormatToIterator[format]
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