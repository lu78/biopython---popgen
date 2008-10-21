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
        
    def __str__(self):
        """
        Represent the record
        
        example:
        # put an example here
        """
        rep  = [self.comment_line + '\n']
        rep.append('\n'.join(self.loci_list) + '\n')
        for pop in self.populations:
            rep.append('Pop\n')
            for indiv in pop:
                name, markers = indiv
                rep.append(name)
                rep.append(',')
                for marker in markers:
                    rep.append(' ')
                    for al in marker:
                        if al == None:
                            al = '0'
                        aStr = str(al)
                        while len(aStr) < self.marker_len:
                            aStr = "".join(['0', aStr])
                        rep.append(aStr)
                rep.append('\n')
        return "".join(rep)

    def split_in_pops(self, pop_names):
        """Splits a GP record in a dictionary with 1 pop per entry.

            Given a record with n pops and m loci returns a dictionary
            of records (key pop_name) where each item is a record
            with a single pop and m loci.

            Parameters:
            pop_names - Population names
        """
        gp_pops = {}
        for i in range(len(self.populations)):
            gp_pop = GenePop.Record()
            gp_pop.marker_len = self.marker_len
            gp_pop.comment_line = self.comment_line
            gp_pop.loci_list = deepcopy(self.loci_list)
            gp_pop.populations = [deepcopy(self.populations[i])]
            gp_pops[pop_names[i]] = gp_pop
        return gp_pops

    def split_in_loci(self, gp):
        """Splits a GP record in a dictionary with 1 locus per entry.

            Given a record with n pops and m loci returns a dictionary
            of records (key locus name) where each item is a record
            with a single locus and n pops.
        """
        gp_loci = {}
        for i in range(len(self.loci_list)):
            gp_pop = GenePop.Record()
            gp_pop.marker_len = self.marker_len
            gp_pop.comment_line = self.comment_line
            gp_pop.loci_list = [self.loci_list[i]]
            gp_pop.populations = []
            for pop in self.populations:
                my_pop = []
                for indiv in pop:
                    my_pop.append((indiv[0], [indiv[1][i]]))
                gp_pop.populations.append(my_pop)
            gp_loci[gp_pop.loci_list[0]] = gp_pop
        return gp_loci


    def remove_population(self, pos):
        """Removes a population (by position).
        """
        del self.populations[pos]
    
    def remove_locus_by_position(self, pos):
        """Removes a locus by position.
        """
        del self.loci_list[pos]
        for pop in self.populations:
            for indiv in pop:
                name, loci = indiv
                del loci[pos]

    def remove_locus_by_name(self, name):
        """Removes a locus by name.
        """
        for i in range(len(self.loci_list)):
            if self.loci_list[i] == name:
                self.remove_locus_by_position(i)
                return
        #If here than locus not existent... Maybe raise exception?
        #   Although it should be Ok... Just a boolean return, maybe?