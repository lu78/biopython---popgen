#!/usr/bin/env python
# Copyright 2008 by Giovanni Dall'Olio.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


class AbstractPopRecord(object):
    """Abstract class for GenePop/PEP/etc parsers

    Holds information from a GenePop/PEP record.

    Members:
    marker_len         The marker length (2 or 3 digit code per allele).    
    
    comment_line       Comment line.

    loci_list          List of loci names.
    
    populations        List of population data.
    
    populations has one element per population. Each element is itself
    a list of individuals, each individual is a pair composed by individual
    name and a list of alleles (2 per marker): Example
    [
        [
            ('Ind1', [(1,2),    (3,3), (200,201)],
            ('Ind2', [(2,None), (3,3), (None,None)],
        ],
        [
            ('Other1', [(1,1),  (4,3), (200,200)],
        ]
    ]
    
    """

    def __init__(self):
        self.marker_len      = 0
        self.comment_line    = ""
        self.loci_list       = []
        self.populations     = []