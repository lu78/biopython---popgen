#!/usr/bin/env python
# Copyright 2008 by Giovanni Dall'Olio.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
Genotype object

>>> import Individual, Marker
>>> g = Genotype('snr10000')
>>> g.genotypes = ['AA', 'AC', 'CC', 'AC']
"""


class Genotype(object):
    def __init__(self, name):
        self.genotype_id = name
        self.genotypes = []
    