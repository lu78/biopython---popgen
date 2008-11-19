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
>>> genotypes_file = StringIO('''
...   Ind1    Ind2    Ind3    Ind4    Ind5    Ind6    Ind7    Ind8    Ind9    Ind10
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
>>> h = hgdp_genotypesIterator(genotypes_file)     # what HGDP iterator should return?
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
... ''')

Example of file with Individuals annotations
>>> samples_file = StringIO('''
... HGDPIndividualNumber    PopulationCode    PopulationName    SamplingLocation    GeographicRegionOfPopulation    Sex    InHGDP.CEPHpanel.CannEtAl2002..    AnalyzedInRosenbergEtAl2002.    PopulationLabelBelievedToBeCorrect.    HasNoDuplicatesInPanel.    AnalyzedInRosenbergEtAl2005.datasetH1048..    HasNoKnown1stDegreeRelativesInPanel.    HasNoKnown1stOr2ndDegreeRelativesInPanel.    HasAParentOrOffspringInPanel.    IncludedInDataset971.No1stDegreeRelatives..    IncludedInDataset952.No1stOr2ndDegreeRelatives..    DuplicationConjugate.OnlyDiffersFromSampleNumberIfSampleIsDuplicated.    AlternatePopulationCode.OnlyDiffersFromPrimaryPopulationCodeForBantuSouthAfrica.    AlternatePopulationName.OnlyDiffersFromPrimaryPopulationNameForBantuSouthAfrica.    OrangeQInRosenbergEtAl2002Fig1K.5    BlueQInRosenbergEtAl2002Fig1K.5    PinkQInRosenbergEtAl2002Fig1K.5    GreenQInRosenbergEtAl2002Fig1K.5    PurpleQInRosenbergEtAl2002Fig1K.5    hgdp.id    population    studySetLiNorel    continent    Li
... 1    51    Brahui    Pakistan    CENTRAL_SOUTH_ASIA    m    1    1    1    1    1    1    1    0    1    1    1    51    Brahui    0.001    0.995    0.002    0.001    0.001    HGDP00001    Brahui    TRUE    CSASIA    TRUE
... 3    51    Brahui    Pakistan    CENTRAL_SOUTH_ASIA    m    1    1    1    1    1    1    1    0    1    1    3    51    Brahui    0.002    0.99    0.003    0.001    0.004    HGDP00003    Brahui    TRUE    CSASIA    TRUE
... 5    51    Brahui    Pakistan    CENTRAL_SOUTH_ASIA    m    1    1    1    1    1    1    1    0    1    1    5    51    Brahui    0.003    0.865    0.096    0.031    0.006    HGDP00005    Brahui    TRUE    CSASIA    TRUE
... 7    51    Brahui    Pakistan    CENTRAL_SOUTH_ASIA    m    1    1    1    1    1    1    1    0    1    1    7    51    Brahui    0.001    0.889    0.107    0.002    0.002    HGDP00007    Brahui    TRUE    CSASIA    TRUE
... 9    51    Brahui    Pakistan    CENTRAL_SOUTH_ASIA    m    1    1    1    1    1    1    1    0    1    1    9    51    Brahui    0.006    0.902    0.053    0.006    0.033    HGDP00009    Brahui    TRUE    CSASIA    TRUE
... 11    51    Brahui    Pakistan    CENTRAL_SOUTH_ASIA    m    1    1    1    1    1    1    1    0    1    1    11    51    Brahui    0.004    0.972    0.013    0.003    0.008    HGDP00011    Brahui    TRUE    CSASIA    TRUE
... 13    51    Brahui    Pakistan    CENTRAL_SOUTH_ASIA    m    1    1    1    1    1    1    1    0    1    1    13    51    Brahui    0.216    0.761    0.014    0.006    0.003    HGDP00013    Brahui    TRUE    CSASIA    TRUE
... 15    51    Brahui    Pakistan    CENTRAL_SOUTH_ASIA    m    1    1    1    1    1    1    1    0    1    1    15    51    Brahui    0.001    0.878    0.109    0.004    0.008    HGDP00015    Brahui    TRUE    CSASIA    TRUE
... 17    51    Brahui    Pakistan    CENTRAL_SOUTH_ASIA    m    1    1    1    1    1    1    1    0    1    1    17    51    Brahui    0.002    0.984    0.008    0.005    0.001    HGDP00017    Brahui    TRUE    CSASIA    TRUE
... 19    51    Brahui    Pakistan    CENTRAL_SOUTH_ASIA    m    1    1    1    1    1    1    1    0    1    1    19    51    Brahui    0.062    0.912    0.008    0.015    0.002    HGDP00019    Brahui    TRUE    CSASIA    TRUE
... 21    51    Brahui    Pakistan    CENTRAL_SOUTH_ASIA    m    1    1    1    1    1    1    1    0    1    1    21    51    Brahui    0.002    0.994    0.001    0.001    0.002    HGDP00021    Brahui    TRUE    CSASIA    TRUE
... 23    51    Brahui    Pakistan    CENTRAL_SOUTH_ASIA    m    1    1    1    1    1    1    1    0    1    1    23    51    Brahui    0.004    0.934    0.025    0.008    0.029    HGDP00023    Brahui    TRUE    CSASIA    TRUE
... 25    51    Brahui    Pakistan    CENTRAL_SOUTH_ASIA    m    1    1    1    1    1    1    1    0    1    1    25    51    Brahui    0.002    0.867    0.098    0.028    0.004    HGDP00025    Brahui    TRUE    CSASIA    TRUE
... ''') 
"""

import logging
from PopGen.Gio.Individual import Individual

def hgdp_genotypesIterator(handle, markers_filter = None, samples_filter = None):
    """
    Parse a genotypes file handler.
    
    It returns a Marker object for every line of the file  
    #>>> for snp in hgdp_genotypes_iterator(handle):
    #...     print snp
    <snp object at ...>
    """
    # read the header, containing the Individuals names
    handle.readline()       # first line is empty??
    header = handle.readline()
    if header is None:
        raise ValueError('Empty file!!')
    individuals = [Individual(id) for id in header.split()]
    
    for line in handle.readlines():
        fields = line.split()
        if fields is None:
            break
        for n in range(1, len(fields)):
            individuals[n-1].markers.append(fields[n])
            logging.debug(individuals[n-1])
            logging.debug(individuals[n-1].markers)




def _test():
    """test the module"""
    import doctest
    doctest.testmod()
    
if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    _test()