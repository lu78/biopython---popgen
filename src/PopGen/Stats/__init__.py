# Copyright 2007 by Tiago Antao.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


"""
  Init for PopGen stats.
  
  Some shared code is here.
"""


class RequiresGenotypeException(Exception):
    """A certain statistic requires a Genotype.
    
    The exception is raised if a certain statistic requires genotype
    information (i.e., both alleles for an individual), as opposed to
    only have to supply allele counts per population/whole sample.
    
    Most statistics don't require genotypical information, one exception
    is observed heterozygosity (Ho).
    """
    pass

