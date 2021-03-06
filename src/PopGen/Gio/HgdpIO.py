#!/usr/bin/env python
# biopython license
"""
Parse files in HGDP (I don't know how its proper name) format.

HGDP results consists in three different file: 
- one with the genotypes
- one with gene informations
- one with individual informations

Example of Genotypes file
>>> from StringIO import StringIO
>>> genotypes_file = StringIO(
... '''  Ind1    Ind2    Ind3    Ind4    Ind5    Ind6    Ind7    Ind8    Ind9    Ind10
... MitoA10045G    AA    GG    AG    AA    AA    AA    AA    AA    AA    AA   
... rs1112391    TT    TC    CC    CC    CC    CC    CC    CC    CC    CC
... MitoA11252G    AA    AA    AA    AA    AA    AA    AA    AA    AA    AA
... rs11124185    TC    TT    TT    TT    TT    TT    TT    TT    TT    TT
... MitoA11468G    AA    AA    AA    AA    AA    AA    AA    AA    AA    AA
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
... ''')
>>> individuals = hgdpgenotypesParser(genotypes_file)     # what HGDP iterator should return?
>>> 

Example of file with SNP annotations
>>> gene_annotation = StringIO(
... '''SNP    Chromosome    Coordinate    GenomeBuild    GeneSymbol    Gene    Location    LocationRelativeToGene    Coding_Status    AminoAcidChange    IDwithMouse    PhastConservation
... rs3094315    1    792430    35    FLJ22639    NM_024796    flanking_3UTR    -9020    -1    -1    -1    -1
... rs12562034    1    808312    35    FLJ22639    NM_024796    flanking_5UTR    -5563    -1    -1    -1    -1
... rs3934834    1    1045730    35    C1orf159    NM_017891    flanking_3UTR    -11398    -1    -1    -1    -1
... rs9442372    1    1058628    35    C1orf159    NM_017891    intron    -338    -1    -1    -1    -1
... rs3737728    1    1061339    35    C1orf159    NM_017891    intron    -24    -1    -1    -1    -1
... rs11260588    1    1061582    35    C1orf159    NM_017891    intron    -267    -1    -1    -1    -1
... rs9442398    1    1061618    35    C1orf159    NM_017891.2    intron    -303    -1    -1    -1    -1
... rs6687776    1    1070489    35    C1orf159    NM_017891    intron    -3083    -1    -1    -1    -1
... ''')

Example of file with Individuals annotations
>>> samples_file = StringIO(
... '''"sample"
... "code"  "sex"    "population"    "region"        "continent"     "Unidad"
... "HGDP00001"    "M"    "Brahui"    "Pakistan"      "Asia"  "Brahui"
... "HGDP00003"    "M"    "Brahui"    "Pakistan"      "Asia"  "Brahui"
... "HGDP01362"    "M"    "French Basque"    "France"    "Europe test"    "Basque"
... "HGDP00151"    "F"    "Makrani"    "Pakistan"    "Asia"    "Makrani"
... "HGDP00013"    "M"    "Brahui test"    "Pakistan"      "Asia"  "Brahui"
... "HGDP00015"    "M"    "Brahui"    "Pakistan test"      "Asia"  "Brahui"
... ''') 
>>> hgdpsamplesfileParser(samples_file)
[Mr. HGDP00001 (Brahui), Mr. HGDP00003 (Brahui), Mr. HGDP01362 (French Basque), Mr. HGDP00151 (Makrani), Mr. HGDP00013 (Brahui test), Mr. HGDP00015 (Brahui)]
"""

import logging
from PopGen.Gio.Individual import Individual
#from PopGen.Gio.Genotype import Genotype
from PopGen.Gio.Marker import Marker
import csv
import re

def hgdpsamplesfileParser(handle, ):
    """
    parse a file with descriptions of Individuals (samples) in hgdp
    
    >>> from StringIO import StringIO
    >>> samples_file = StringIO(
    ... '''"sample"
    ... "code"  "sex"    "population"    "region"        "continent"     "Unidad"
    ... "HGDP00001"    "M"    "Brahui Test"    "Pakistan"      "Asia"  "Brahui"
    ... "HGDP00003"    "M"    "Brahui"    "Pakistan"      "Asia"  "Brahui"
    ... "HGDP01362"    "M"    "French Basque"    "France"    "Europe"    "Basque"
    ... "HGDP00151"    "F"    "Makrani"    "Pakistan"    "Asia"    "Makrani"''')
    >>> samples = hgdpsamplesfileParser(samples_file)
    >>> print [sample for sample in samples if sample.region == "Pakistan"]
    [Mr. HGDP00001 (Brahui Test), Mr. HGDP00003 (Brahui), Mr. HGDP00151 (Makrani)]
    """
    splitter = re.compile('"\s+"')
    handle.readline()   # skip headers
    header = handle.readline()
    if header is None:
        raise ValueError('Empty file!!')
    
    individuals = []
    
#    for line in handle.readlines():    
    line = handle.readline()

    while line: 
        row = splitter.split(line)
        
        if row is None: break   # FIXME: optimize
        if len(row) != 6:
            raise ValueError("wrong number of columns in current line")
        
        ind_id = row[0].replace('"', '')
#        logging.debug(row)
        sex = row[1]    # TODO: translate this to 1/2 
        pop = row[2]    # TODO: use Population object
        region = row[3]
        continent = row[4]
        unit = row[5].replace('"', '')
        
        # create an Individual object
        Ind = Individual(ind_id, pop, region=region, continent=continent, 
                        working_unit=unit, sex=sex)
        individuals.append(Ind)
        
        line = handle.readline()

    return individuals
    

def hgdpgenotypesParser(handle, individuals_filter = None, markers_filter = None):
    """
    Parse a genotypes file handler.
    
    It returns a Marker object for every line of the file
    
    >>> from StringIO import StringIO
    >>> genotypes_file = StringIO(
    ... '''  HGDP00001    HGDP00002    HGDP00003    HGDP00004    HGDP00005    HGDP00006    HGDP00007    HGDP00008    HGDP00009    HGDP000010
    ... rs1112390    AA    GG    AG    AA    AA    AA    AA    AA    AA    AA   
    ... rs1112391    TT    TC    CC    CC    CC    CC    CC    CC    CC    CC
    ... MitoA11252G    AA    AA    AA    AA    AA    AA    AA    AA    AA    AA
    ... rs11124185    TC    TT    TT    TT    TT    TT    TT    TT    TT    TT
    ... MitoA13265G    AA    AA    AA    AA    AA    AA    AA    AA    AA    AA
    ... MitoA13264G    GG    AA    AA    AA    GG    AG    AA    AA    AA    AA
    ... MitoA13781G    AA    AA    AA    AA    AA    AA    --    AA    AA    AA
    ... MitoA14234G    AA    AA    AA    AA    AA    AA    AA    AA    AA    AA
    ... MitoA14583G    AA    AA    AA    AA    AA    AA    AA    AA    AA    AA
    ... MitoA14906G    GG    GG    GG    GG    GG    GG    GG    GG    GG    GG
    ... MitoA15219G    AA    AA    AA    GG    AA    AA    AA    AA    AA    AA''')
    
    >>> individuals_filter = ['HGDP00001', 'HGDP00004', ]  
    >>> markers = hgdpgenotypesParser(genotypes_file, individuals_filter)
    >>> print '\t' + '\t'.join(markers[0].individuals)
     HGDP00001   HGDP00004
    >>> for marker in markers:
    ...    print marker.to_geno_format()    #doctest: +NORMALIZE_WHITESPACE
    rs1112390    AA    AA    
    rs1112391    TT    CC    
    MitoA11252G    AA    AA    
    rs11124185    TC    TT    
    MitoA13265G    AA    AA    
    MitoA13264G    GG    AA    
    MitoA13781G    AA    AA    
    MitoA14234G    AA    AA    
    MitoA14583G    AA    AA    
    MitoA14906G    GG    GG    
    MitoA15219G    AA    GG    

    """
    # initialize output var
    markers = []
    
    # read the header, containing the Individuals names
#    handle.readline()       # first line is empty??
    header = handle.readline()
    if header is None:
        raise ValueError('Empty file!!')
    individuals = [Individual(ind_id) for ind_id in header.split()]
    if individuals_filter is None:      # TODO: ugly 
        individuals_filter = [ind.individual_id for ind in individuals]
    
    columns_to_filter = []
    for ind in header.split():
        if ind in individuals_filter:
            columns_to_filter.append(1)
        else:
            columns_to_filter.append(0)
    
    # Read the remaining lines of genotypes file, containin genotypes info.
    for line in handle.readlines():
        fields = line.split()   # TODO: add more rigorous conditions
        if fields is None:
            break
        # Initialize a Genotype object 
        marker = Marker(name = fields[0], individuals = individuals_filter)
        markers.append(marker)
        
        for n in range(1, len(fields)):
            current_individual = individuals[n-1]
            if current_individual in individuals_filter:    #TODO: this consumes CPU time
#            if columns_to_filter[0] == 1:
                current_genotype = fields[n]
                marker.add_genotype(current_genotype)
#                print current_individual
            else:
                pass
            
    return markers


def test_doc():
    """test the module"""
    import doctest
    doctest.testmod()
    
if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    test_doc()
