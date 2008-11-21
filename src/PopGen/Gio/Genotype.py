#!/usr/bin/env python
# Copyright 2008 by Giovanni Dall'Olio.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
Genotype object

>>> import Individual, Marker
>>> g = Genotype('snr10000')
>>> g.samples = ['AA', 'AC', 'CC', 'AC']
"""


class Genotype(object):
    """
    # TODO: rewrite everything :(
    >>> g = Genotype('snr10000')
    >>> g.samples = ['AA', 'AC', 'CC', 'AC']
    """
    def __init__(self, name, individuals = None, values = None):   # careful you are messing with synonimous
        self.genotype_id = name
        
        if values is None:
            self.values = []
        else:
            self.samples = values
        
        if individuals is None:
            individuals = []
        else:
            self.individuals = individuals
    