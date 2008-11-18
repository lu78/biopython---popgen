#!/usr/bin/env python
# biopython license
"""

Parses files in HGDP (I don't know how to call it) format.

Example:
>>> import StringIO
>>> hgdp_file = StringIO('''
... [Header]
... BSGT Version    2.3.41.16318
... Processing Date    8/20/2006 1:40 PM
... Content        BDCHP-1X12-HUMANHAP650Y_11226216_A.csv
... Num SNPs    19
... Total SNPs    19
... Num Samples    30
... Total Samples    159
... [Data]
...     Ind1    Ind2    Ind3    Ind4    Ind5    Ind6    Ind7    Ind8    Ind9    Ind10
... MitoA10045G    AA    AA    AA    AA    AA    AA    AA    AA    AA    AA   
... MitoA10551G    AA    AA    AA    AA    AA    AA    AA    AA    AA    AA
... MitoA11252G    AA    AA    AA    AA    AA    AA    AA    AA    AA    AA
... MitoA11468G    AA    AA    AA    AA    AA    AA    AA    AA    AA    AA
... MitoA11813G    AA    AA    AA    AA    AG    AA    AA    AA    AA    AA
... MitoA12309G    AA    AA    AA    AA    AA    AA    AA    AA    AA    AA
... MitoA13106G    AA    AA    AA    AA    AA    AA    AA    AA    AA    AA
... MitoA13264G    GG    AA    AA    AA    GG    AG    AA    AA    AA    AA
... MitoA13781G    AA    AA    AA    AA    AA    AA    --    AA    AA    AA
... MitoA14234G    AA    AA    AA    AA    AA    AA    AA    AA    AA    AA
... MitoA14583G    AA    AA    AA    AA    AA    AA    AA    AA    AA    AA
... MitoA14906G    GG    GG    GG    GG    GG    GG    GG    GG    GG    GG
... MitoA15219G    AA    AA    AA    GG    AA    AA    AA    AA    AA    AA
... MitoA15245G    AA    AA    AA    AA    AA    AA    AA    AA    AA    AA
... MitoA15302G    AA    AA    AA    AA    AA    GG    GG    AA    AA    GG
... MitoA15759G    AA    AA    AA    AA    AA    AA    AA    AA    AA    AA
... MitoA15908G    AA    AA    AA    AA    AA    AA    AA    AA    AA    AA
... rs1112391    TC    CC    CC    CC    CC    CC    CC    CC    CC    CC
... rs11124185    TC    TT    TT    TT    TT    TT    TT    TT    TT    TT
... ''')
>>> h = HGDPIterator(hgdp_file)     # what HGDP iterator should return?
>>> 
"""

def HGDPIterator(handle, markers_filter = None, samples_filter = None):
    pass




def _test():
    import doctest
    doctest.testmod()
    
if __name__ == '__main__':
    _test()