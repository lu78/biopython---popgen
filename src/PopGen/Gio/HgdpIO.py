#!/usr/bin/env python
# biopython license
"""

Parses files in HGDP (I don't know how to call it) format.

HGDP results consists in three different file: 
- one with the genotypes
- one with gene informations
- one with individual informations

Example of Genotype file
>>> import StringIO
>>> genotypes_file = StringIO('''
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
>>> h = HGDPIterator(genotypes_file)     # what HGDP iterator should return?
>>> 

Example of file with SNP annotations
>>> gene_annotation = StringIO('''
... SNP    Chromosome    Coordinate    GenomeBuild    GeneSymbol    Gene    Location    LocationRelativeToGene    Coding_Status    AminoAcidChange    IDwithMouse    PhastConservation
... rs3094315    1    792430    35    FLJ22639    NM_024796    flanking_3UTR    -9020    -1    -1    -1    -1
... rs12562034    1    808312    35    FLJ22639    NM_024796    flanking_5UTR    -5563    -1    -1    -1    -1
... rs3934834    1    1045730    35    C1orf159    NM_017891    flanking_3UTR    -11398    -1    -1    -1    -1
... rs9442372    1    1058628    35    C1orf159    NM_017891    intron    -338    -1    -1    -1    -1
... rs3737728    1    1061339    35    C1orf159    NM_017891    intron    -24    -1    -1    -1    -1
... rs11260588    1    1061582    35    C1orf159    NM_017891    intron    -267    -1    -1    -1    -1
... rs9442398    1    1061618    35    C1orf159    NM_017891.2    intron    -303    -1    -1    -1    -1
... rs6687776    1    1070489    35    C1orf159    NM_017891    intron    -3083    -1    -1    -1    -1
... rs9651273    1    1071464    35    C1orf159    NM_017891    intron    -4058    -1    -1    -1    -1
... rs4970405    1    1088879    35    C1orf159    NM_017891    intron    -2484    -1    -1    -1    -1
... rs12726255    1    1089874    35    C1orf159    NM_017891    intron    -1489    -1    -1    -1    -1
... rs7540009    1    1100158    35    C1orf159    NM_017891.2    flanking_5UTR    -8766    -1    -1    -1    -1
... rs11807848    1    1101090    35    C1orf159    NM_017891    flanking_5UTR    -9698    -1    -1    -1    -1
... ''')

Example of file with Individuals annotations
>>> 
"""

def HGDPgenotypesIterator(handle, markers_filter = None, samples_filter = None):
    pass




def _test():
    import doctest
    doctest.testmod()
    
if __name__ == '__main__':
    _test()