

import TpedIO
import PedIO

_FormatToIterator = {
                     'tped' : TpedIO.TpedIterator,
                     'ped'  : PedIO.PedIterator
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
