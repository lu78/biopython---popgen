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
... 27    51    Brahui    Pakistan    CENTRAL_SOUTH_ASIA    m    1    1    1    1    1    1    1    0    1    1    27    51    Brahui    0.01    0.957    0.024    0.004    0.005    HGDP00027    Brahui    TRUE    CSASIA    TRUE
... 29    51    Brahui    Pakistan    CENTRAL_SOUTH_ASIA    m    1    1    1    1    1    1    1    0    1    1    29    51    Brahui    0.175    0.688    0.049    0.086    0.002    HGDP00029    Brahui    TRUE    CSASIA    TRUE
... 31    51    Brahui    Pakistan    CENTRAL_SOUTH_ASIA    m    1    1    1    1    1    1    1    0    1    1    31    51    Brahui    0.001    0.975    0.009    0.012    0.003    HGDP00031    Brahui    TRUE    CSASIA    TRUE
... 33    51    Brahui    Pakistan    CENTRAL_SOUTH_ASIA    m    1    1    1    1    1    1    1    0    1    1    33    51    Brahui    0.001    0.991    0.003    0.002    0.003    HGDP00033    Brahui    TRUE    CSASIA    TRUE
... 35    51    Brahui    Pakistan    CENTRAL_SOUTH_ASIA    m    1    1    1    1    1    1    1    0    1    1    35    51    Brahui    0.001    0.903    0.085    0.008    0.003    HGDP00035    Brahui    TRUE    CSASIA    TRUE
... 37    51    Brahui    Pakistan    CENTRAL_SOUTH_ASIA    m    1    1    1    1    1    1    1    0    1    1    37    51    Brahui    0.001    0.95    0.045    0.001    0.003    HGDP00037    Brahui    TRUE    CSASIA    TRUE
... ''') 
"""

def HGDPgenotypesIterator(handle, markers_filter = None, samples_filter = None):
    pass




def _test():
    import doctest
    doctest.testmod()
    
if __name__ == '__main__':
    _test()