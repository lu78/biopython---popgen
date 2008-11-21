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
... "HGDP01362"    "M"    "French_Basque"    "France"    "Europe"    "Basque"
... "HGDP00151"    "F"    "Makrani"    "Pakistan"    "Asia"    "Makrani"
... "HGDP00013"    "M"    "Brahui"    "Pakistan"      "Asia"  "Brahui"
... "HGDP00015"    "M"    "Brahui"    "Pakistan"      "Asia"  "Brahui"
... ''') 
>>> hgdpsamplesfileIterator(samples_file)
{'Pakistan': [Mr. HGDP00001 (Brahui), Mr. HGDP00003 (Brahui), Mr. HGDP00151 (Makrani), Mr. HGDP00013 (Brahui), Mr. HGDP00015 (Brahui)], 'France': [Mr. HGDP01362 (French_Basque)]}
"""

import logging
from PopGen.Gio.Individual import Individual
import csv
import re

def hgdpsamplesfileIterator(handle, ):
    """
    parse a file with descriptions onn Individuals (samples) in hgdp
    """
    handle.readline()   # skip headers
    header = handle.readline()
    if header is None:
        raise ValueError('Empty file!!')
    
#    individuals_by_population = {}
    individuals_by_region = {}
#    individuals_by_continent = {}
    
    for line in handle.readlines(): 
        row = line.split()
        if row is None: break   # FIXME: optimize
        id = row[0].replace('"', '')
#        logging.debug(row)
        sex = row[1].replace('"', '')   # TODO: translate this to 1/2 
        pop = row[2].replace('"', '')
        region = row[3].replace('"', '')
        continent = row[4].replace('"', '')
        unit = row[5].replace('"', '')
        
        # create an Individual object
        Ind = Individual(id, pop, sex = sex)
        
#        individuals_by_population.setdefault(pop, [])
#        individuals_by_population[pop].append(id)
        individuals_by_region.setdefault(region, [])
        individuals_by_region[region].append(Ind)
#        individuals_by_continent.setdefault(continent, [])
#        individuals_by_continent[continent].append(id)
    logging.debug(individuals_by_region)
#    logging.debug(individuals_by_population)
    
    return individuals_by_region
    

def hgdpgenotypesParser(handle, samples_filter = None, markers_filter = None):
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
    
    >>> samples_filter = ['HGDP00001', 'HGDP00004']  
    >>> individuals = hgdpgenotypesParser(genotypes_file, samples_filter)
    >>> for ind in individuals:
    ...     print ind, ind.markers
    HGDP00001 ['AA', 'TT', 'AA', 'TC', 'AA', 'GG', 'AA', 'AA', 'AA', 'GG', 'AA']
    HGDP00004 ['AA', 'CC', 'AA', 'TT', 'AA', 'AA', 'AA', 'AA', 'AA', 'GG', 'GG']
    """
    # read the header, containing the Individuals names
#    handle.readline()       # first line is empty??
    header = handle.readline()
    if header is None:
        raise ValueError('Empty file!!')
    individuals = [Individual(ind_id) for ind_id in header.split()]
    
    for line in handle.readlines():
        fields = line.split()
        if fields is None:
            break
        for n in range(1, len(fields)):
            current_individual = individuals[n-1]
            if current_individual in samples_filter:
                current_marker = fields[n]
                current_individual.markers.append(current_marker) # or append?
#                logging.debug(current_individual)
#                print current_individual
            else:
                pass
            
        # for every line in the file, return a list of 'Individual' object
    filtered_individuals = [ind for ind in individuals if ind in samples_filter]
    return filtered_individuals




def _test():
    """test the module"""
    import doctest
    doctest.testmod()
    
if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    _test()