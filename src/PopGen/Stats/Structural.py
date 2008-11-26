# Copyright 2007 by Tiago Antao.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
  Implementation of Statistics dependent on population structure.
  
  Here can be found statistics that work with population structure
  (e.g., Fst). These statistics require several populations, and
  cannot operate only on (unstructured) sets of individuals.
  
  An 'abstract' class is available, with all methods that need
  to be implemented
  
  Obviously other Structural statistics can be implemented
  outside of this module.
"""

from copy import deepcopy
from PopGen.Exceptions import PopulationExistsException, RequiresGenotypeException
import logging

class Structural:
    """'Abstract' Structural class.

    A constructor is available. When overriden on subclasses
    one parameter might have to be changed:
    counts_acceptable - Which is True (the default)
    if the statistic only needs allele counts (and not genotypical) information,
    i.e., which (two) alleles exist per individual
    
    The following members are made available:
    
    add_pop(pop_name, indiv_data) - Adds a population,
      indiv_data is a list of pairs
      [('a','a'), ('a','b'), ('b','a'), ...], with the two
      alleles per individual. If you only have allele counts
      then use add_pop_counts
      
    add_pop_counts(pop_name, counts) - Adds a population,
      counts is a dictionary of allele_name : count.
      For some statistics counts might not be enough, in that
      case only add_pop can be used
      
    Each subclass needs to implement:
        
    - calc_stat() - Calculates the statistic
    - __str__
    - acronym        # short name for the test type (Fst, FstB...)
    
    >>> s = Structural()      # Usually this will remain as an abstract class
    >>> s.counts_acceptable = False     
    >>> Ind1 = (1, 2)         # only one locus per individual
    >>> Ind2 = (2, None)
    >>> Other1 = (1, 2)
    >>> s.add_pop('Vulcanians', [Ind1, Ind2])
    >>> s.add_pop('Martians', [Other1])
    >>> s.pop_names
    ['Vulcanians', 'Martians']
    >>> s.pop_indivs
    {'Vulcanians': [(1, 2), (2, None)], 'Martians': [(1, 2)]}
    
    TODO: test
"""

    def __init__(self):
        self.counts_acceptable = True   # ???
        self.pop_counts = {}
        self.pop_names = []
        self.pop_indivs = {}
        self.removable = ['000', '00', '0', 'NA']   # What does it means?
    
    def add_pop(self, pop_name, indiv_data):
        """
        Use this if the statistics needs allele data
        
        for every individual, a list of tuples is provided
        every tuple represent the two values in the locus
        >>> s = Structural()
        >>> s.add_pop('Martians', [('A', 'G'), ('G', 'C')])
        """
        if pop_name in self.pop_names:
            raise PopulationExistsException(pop_name)   # or maybe append?

        self.pop_counts[pop_name] = self._convert_to_counts(indiv_data) # We convert to counts
        if self.counts_acceptable is False: 
            self.pop_names.append(pop_name)
            self.pop_indivs[pop_name] = indiv_data
    
    def add_pop_counts(self, pop_name, counts):
        """
        Use this if the statistic only needs allele counts.
        
        counts is a dictionary of the occurrencies for every allele
        >>> s = Structural()
        >>> s.add_pop_counts('Martians', {'A': 10, 'C': 5, 'G': 6, 'T': 3})
        """
        if pop_name in self.pop_names:
            raise PopulationExistsException(pop_name)
        if self.counts_acceptable is False:
            raise RequiresGenotypeException
        else:
            self.pop_names.append(pop_name)
            self.pop_counts[pop_name] = counts
            
    def calc_stat(self):
        raise NotImplementedError

    def _convert_to_counts(self, indiv_data):
        """ 
        >>> s = Structural()
        >>> print s._convert_to_counts([('A', 'C'), ('C', 'G')])
        {'A': 1, 'C': 2, 'G': 1}
        """
        count_data = {}
        for indiv in indiv_data:
            for allele in indiv:
                if allele is not None:
                    allele_count = count_data.get(str(allele), 0)
                    count_data[str(allele)] = allele_count + 1
        return count_data

class FstBeaumont(Structural):
    """Implements Fst a la Beaumont.
    
    This is a variation of Weir and Cockerham, should only be
    used with FDist. A Weir and Cockerham Theta is be implemented.
    
    TODO: test
    """
    
    def __init__(self):
        Structural.__init__(self)
        self.acronym = 'FstB'
        
    def __str__(self):
        return "Fst, Beaumont Style"

    def calc_stat(self):
        x0 = 0.0
        skip = 0
        for pop_count in self.pop_counts:
            alleles = self.pop_counts[pop_count]
            sample_size = 0
            for allele in alleles:
                sample_size += alleles[allele]
            if sample_size == 0:
                skip += 1
            else:
                x2 = 0.0
                for allele in alleles.keys():
                  x2 += alleles[allele]**2
                x0 += (x2-sample_size)/(sample_size*(sample_size-1))
        yy = 0.0
        pops = self.pop_counts.keys()
        for j in range(0, len(pops)):
            alleles = self.pop_counts[pops[j]]
            sample_size = 0
            for allele in alleles:
                sample_size += alleles[allele]
            if sample_size != 0:
                for k in range(j+1, len(pops)):
                    alleles2 = self.pop_counts[pops[k]]
                    sample_size2 = 0
                    for allele2 in alleles2:
                        sample_size2 += alleles2[allele2]
                    if sample_size2 != 0:
                        y1 = 0.0
                        for all in alleles.keys():
                            if alleles2.has_key(all):
                                y1 += alleles[all]*alleles2[all]
                        yy += y1/(sample_size*sample_size2)
        real_samples = len(self.pop_counts)-skip
        q2 = x0/real_samples
        q3 = 2*yy/(real_samples*(real_samples-1))

        het0 = 1.0 - q2;
        het1 = 1.0 - q3;
        if het1 < 1.0e-10:
            self.fst = -100.0;
        else:
            self.fst = 1.0 - het0/het1;
        self.he = het1
        self.h0 = het0
        return self.fst

def _calc_refs(counts):
    """Returns the sum of all referenced alleles.
    
    If the population is diploid then
       num_indivs = calc_refs / 2
    if is haploid then
       num_indivs = calc_refs

    This is only needed if only count data is available
    """
    refs = 0
    for allele in counts.keys():
        refs += counts[allele]
    return refs

def _get_all_alleles(pop_counts):
    alleles = []
    for pop in pop_counts.keys():
        alleles += _get_all_alleles_pop(pop_counts[pop])
    return list(set(alleles))

def _get_all_alleles_pop(pop_counts):
    alleles = []
    for allele in pop_counts.keys():
        alleles += allele
    return list(set(alleles))

def _get_allele_freq(pop_count, allele):
    all_alleles = _get_all_alleles_pop(pop_count)
    all_alleles_count = 0.0 # integer, but double for float division
    logging.debug(all_alleles)
    for this_allele in all_alleles:        
        all_alleles_count += pop_count[this_allele]
    return 1.0 * pop_count[allele] / all_alleles_count

def _get_het_allele_freq(indivs, allele):
    het_count = 0
    for indiv in indivs:
        if indiv[0] == allele or indiv[1] == allele:
            if indiv[0] != indiv[1]:
                het_count += 1
    return 1.0 * het_count / len(indivs)

def _calcFs(pop_indivs, pop_counts):
    """Calculates F statistics.
    
    Requires diploidy for now. "Stolen" and converted from simuPop C++ code
    """
    pop_names = pop_counts.keys()
    r = len(pop_names)
    n_i = []
    for pop_name in pop_names:
        n_i.append(len(pop_indivs[pop_name]))
    n = reduce(lambda x, y: x+y, n_i)
    n_bar = 1.0 * n / r
    n_c = n
    for ni in n_i:
        n_c -= 1.0*(ni**2)/n
    n_c = n_c / (r-1)

    alleles = _get_all_alleles(pop_counts)
    a = 0.0
    b = 0.0
    c = 0.0
    for allele in alleles:
        p_i = []
        for pop_name in pop_names:
            p_i.append(_get_allele_freq(pop_counts[pop_name], allele))
        p_bar = 0.0
        for i in range(len(p_i)):
            p_bar += n_i[i] * p_i[i]
        p_bar = 1.0 * p_bar / n
        s_2 = 0.0
        for i in range(len(p_i)):
            s_2 += n_i[i] * (p_i[i] - p_bar) * (p_i[i] - p_bar)
        h_bar = 0.0
        for i in range(len(p_i)):
            h_bar += _get_het_allele_freq(pop_indivs[pop_name], allele) *n_i[i]
        h_bar = 1.0 * h_bar / n
        a += n_bar / n_c * (s_2 - (p_bar * (1-p_bar) - (r - 1.0) / r * s_2 - h_bar / 4.0) / (n_bar - 1.0) )
        b += n_bar / (n_bar - 1) * (p_bar * (1-p_bar) - (r - 1.0) / r * s_2 - (2 * n_bar - 1) / (4.0 * n_bar) * h_bar )
        c += h_bar / 2.0
    if a + b + c == 0:
        fst = 0.0
    else:
        fst = a / (a + b + c)
    if a + b + c == 0:
        fit = 1.0
    else:
        fit = (1.0 - c) / (a + b + c)
    if b + c == 0:
        fis = 1.0
    else:
        fis = (1.0 - c) / (b + c)

    return fst, fit, fis


class Fst(Structural):
    """Implements Fst (plus Fis and Fit).
    
    Weir and Cockerham.
    
    TODO: test
    
    >>> s = Fst()
    >>> Ind1 = (1, 2)
    >>> Ind2 = (2, None)
    >>> Other1 = (1, 2)        # BUG: if fails if there are not heterozygotes 
    ...                        #     in a population (e.g. Other = (1, 1)).
    
    >>> s.add_pop('Vulcanians', [Ind1, Ind2])
    >>> s.add_pop('Martians', [Other1, ])
    >>> s.pop_names
    ['Vulcanians', 'Martians']
    >>> s.pop_indivs
    {'Vulcanians': [(1, 2), (2, None)], 'Martians': [(1, 2)]}
    >>> s.pop_counts
    {'Vulcanians': {'1': 1, '2': 2}, 'Martians': {'1': 1, '2': 1}}
    >>> print s.calc_stat()
    
    """
    
    def __init__(self):
        Structural.__init__(self)
        self.acronym = 'Fst'
        self.counts_acceptable = False
        
    def __str__(self):
        return "Fst"

    def calc_stat(self):
        fst, fit, fis = _calcFs(self.pop_indivs, self.pop_counts)
        return fst
    
    def add_pop_counts(self, indiv_data):
        raise RequiresGenotypeException()

class Fk(Structural):
    """Implements Fk.

    See Krimbas & Tsakas 1971 and Pollak 1983.

    Should only have 2 populations (if more are needed, do pairs)
    This only implements the sum part.
    A full Fk requires all loci Fk averaged by number of alleles.
    
    TODO: test
    """

    def __init__(self):
        Structural.__init__(self)
        self.acronym = 'Fk'
        self.counts_acceptable = True
        
    def __str__(self):
        return "Fk"

    def calc_stat(self):
        p1, p2 = self.pop_counts.keys()
        alleles = list(set(p1.keys() + p2.keys()))
        sum = 0.0
        for a in alleles:
            sum += 1.0*(p1.get(a, 0)-p2.get(b, 0))**2 / (       # undefined var b
                    (p1.get(a,0)+p2.get(a,0))/2.0)


def _test():
    """Test the module using doctest"""
    import doctest
    doctest.testmod()
    
if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    _test()