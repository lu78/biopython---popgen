# Copyright 2007 by Tiago Antao.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
'Simple' statistics. That is statistics that operate on a single
population and have no notion of substructure.

If you need statistics to work with population structure please
have a look at Structural.py .

Both simple and structural statistics are single locus (a sequence
can be seen as a single locus, or split on several SNPs, it all
depends on your objective). Multi-locus statistics (like the
ones related to linkage desiquilibrium) will be dealt elsewhere.
"""

from PopGen.Exceptions import *

class Simple:
    """'Abstract' Simple stats class.

    Any concrete subclass must implement calc_stat(self).

    Concrete statistics can either used aggregated data (e.g.,
    50 counts of allele A, 10 of C) or genotypical data (e.g.,
    Individual 1 is C/C, ... individual 30 is A/C). Some
    statistics, like Observed Heterozigosity, REQUIRE genotypical
    data, other, like Expected Heterozigosity, can live with allele
    counts only.

    Concrete classes are expected to override __str__ in order to
    provide a descriptive statistical name.

    For future meta-programming purposes all statistics are expected
    to have a small acronym on acronym.

    Statistics that require genotypical data must set
    counts_acceptable to False

    """

    def __init__(self):
        self.counts_acceptable = True   # ?
        self.removable = ['000', '00', '0', 'NA']
        self.counts = {}
        self.indiv_data = {}

    def set_data(self, indiv_data):
        if self.counts_acceptable is True: #We convert to counts
            self.counts = self._convert_to_counts(indiv_data)
        else:
            self.indiv_data = indiv_data

    def set_counts(self, counts):
        if self.counts_acceptable is False:
            raise RequiresGenotypeException
        else:
            self.counts = counts

    def _convert_to_counts(self, indiv_data):
        count_data = {}
        for indiv in indiv_data:
            for allele in indiv:
                if allele is not None:
                    allele_count = count_data.get(str(allele), 0)
                    count_data[str(allele)] = allele_count + 1
        return count_data


class NumberOfAlleles(Simple):
    """ Number of alleles.

    A simple statistic that returns the number of existing alleles
    of a certain locus in a population.

    As an example, an SNP, can't have more than 4 alleles, and the
    most typical value is 2.
    """
    def __init__(self):
        Simple.__init__(self)
        self.acronym = 'NAll'

    def __str__(self):
        return "Number of alleles"

    def calc_stat(self):
        return len(self.counts.keys())

class ExpectedHeterozygosity(Simple):
    def __init__(self):
        Simple.__init__(self)
        self.acronym = 'ExHe'

    def __str__(self):
        return "Expected heterozygosity"

    def calc_stat(self):
        sumLocus = 0
        for count in self.counts.values():
            sumLocus += count
        exHe = 1.0
        for count in self.counts.values():
            exHe -= (float(count)/float(sumLocus))**2
        return exHe
