# Copyright 2008 by Tiago Antao.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
Hardy-Weinberg
"""

from numpy import asarray
from scipy.stats import chisquare
from Bio.PopGen.Stats import Simple

def return_expected(indiv_data):
    '''Data should be filtered for nulls'''
    naStat = Simple.NumberOfAlleles()
    naStat.set_data(indiv_data)
    numberAlleles = naStat.calc_stat()
    alleles = naStat.counts.keys()
    counts = {}
    expected_counts = {}
    freqs = {}
    for k in naStat.counts.keys():
        freqs[k] = naStat.counts[k]/(2.0*len(indiv_data))
    for a1 in range(len(alleles)):
       for a2 in range(a1,len(alleles)):
           la = [alleles[a1], alleles[a2]]
           la.sort()
           counts[tuple(la)] = 0
           if la[0] == la[1]:
               expected_counts[tuple(la)] = freqs[alleles[a1]]**2
           else:
               expected_counts[tuple(la)] = 2*freqs[alleles[a1]]*freqs[alleles[a2]]
    for indiv in indiv_data:
        la = list(indiv)
        la.sort()
        counts[tuple(la)] += 1
    f_obs = []
    f_exp = []
    for k in counts.keys():
      f_obs.append(counts[k]/(1.0*len(indiv_data)))
      f_exp.append(expected_counts[k])
    #print f_obs, f_exp,
    cs = chisquare(f_obs, asarray(f_exp))
    return cs
