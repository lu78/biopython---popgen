    # Copyright 2008 by Tiago Antao.  All rights reserved.
    # This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
This module provides code to work with LDNe.

Classes:
Record           Holds LDNe OUTPUT data.
RecordParser     Parses a LDNe record (file) into a Record object.

_Scanner         Scans a LDNe record.
_RecordConsumer  Consumes LDNe data to a Record object.

"""

from Bio import File
from Bio.ParserSupport import *
from os import sep, access, F_OK
from sys import path
import re
import decimal
#from copy import deepcopy
#from types import *

#This is a workaround to work with the test system
#In any case the problem is with the test system
for instance in path:
    test_path = instance + sep + sep.join(['Bio', 'PopGen', 'LDNe', 'data'])
    if access(test_path, F_OK):
        builtin_tpl_dir = test_path
        break



class Record:
    """Holds information from a LDNe OUTPIUT file.

    Members:
    mating             Mating model (random, bound).    

    freqs_used         Frequencies used as a cutoff.
    
    populations        List of population data.
                       Composed of (popid, [cutoff_cases])
    
                       For each cutoff the following applies (presented
                       as a tuple for each cutoff)
                       XXX NOTE: ONE JACKNIFE IS MISSING ON LINUX
                       Harmonic mean
                       Independent comparisons
                       Overall r^2
                       Expected r^2 sample
                       Estimated Ne (plus 2 vals for 95% CIs)
                       IGNORING PARAMETRIC AND JACKNIFE

    """
    def __init__(self):
        self.mating      = None
        self.freqs_used  = []
        self.populations = []

    def __str__(self):
        rep  = ['Mating: ' + self.mating + '\nCutoffs: ']
        rep.append(', '.join(self.freqs_used) + '\n')
        pop_count = 1
        #for pop in self.populations:
        for id, fcases in self.populations:
            rep.append('Population ' + str(id) + '\n')
            for hm, ic, or2, er2, (ne, (ne95, ne105), (j95, j105)) in fcases:
                rep.append('  Harmonic Mean             ' + str(hm)  + '\n')
                rep.append('  Independent comparisons   ' + str(ic)  + '\n')
                rep.append('  Overall r^2               ' + str(or2) + '\n')
                rep.append('  Expected r^2 sample       ' + str(er2) + '\n')
                rep.append('  NE (Para 95% CI)  (JN CI) ' + str(ne)  +
                            ' (' + str(ne95) + ', ' + str(ne105) + ') ' +
                            ' (' + str(j95) + ', ' + str(j105) + ')\n')
            pop_count += 1
        return "".join(rep)


class RecordParser(AbstractParser):
    """Parses LDNe data into a Record object.

    """
    def __init__(self):
        self._scanner = _Scanner()
        self._consumer = _RecordConsumer()

    def parse(self, handle):
        self._scanner.feed(handle, self._consumer)
        return self._consumer.data

def parse(handle):
   """Parses a handle containing a LDNe file.
   """
   parser = RecordParser()
   return parser.parse(handle)

class _Scanner:
    """Scans a LDNe record.
    
    There is only one record per file.
    
    """

    def feed(self, handle, consumer):
        """feed(self, handle, consumer)

        Feed in a LDNe unit record for scanning.  handle is a file-like
        object that contains a LDNe record.  consumer is a
        Consumer object that will receive events as the report is scanned.

        """
        if isinstance(handle, File.UndoHandle):
            uhandle = handle
        else:
            uhandle = File.UndoHandle(handle)


        consumer.start_record()
        l = uhandle.readline()
        match = None
        while (match == None):
            if l=='':
                raise ValueError('Mating mode not found')
            match = re.search('Mating Model is ([^ \n]+)', l)
            if (match == None):
                l = uhandle.readline()
        consumer.mating(match.group(1))
        match = None
        while (match == None):
            if l=='':
                raise ValueError('Cutoffs not found')
            match = re.search('Lowest Allele Fr', l)
            if (match == None):
                l = uhandle.readline()
        consumer.cutoff(filter(lambda x : x<>'',
                l.rstrip().split('=')[1].split(' ')))
        while l<>'':
            match =re.search('Population ([0-9])+',l)
            if (match <> None):
                consumer.start_pop(int(match.group(1)))
            self.fetch_values('Harmonic Mean', l, consumer.harmonic_mean)
            self.fetch_values('Independent Comparisons', l, consumer.independent_comparison)
            self.fetch_values('OverAll', l, consumer.overall_r2)
            self.fetch_values('Expected r', l, consumer.expected_r2)
            self.fetch_values('Estimated Ne', l, consumer.estimated_ne)
            if l.startswith('95% CI'):
                l = uhandle.readline().rstrip()
                consumer.parametric_bottom(
                        filter(lambda x: x<>'' and x<>'*' and not x.startswith('P'),
                            [l[30:40], l[40:50], l[50:60]]
                        )
                )
                l = uhandle.readline().rstrip()
                consumer.parametric_top([l[30:40], l[40:50], l[50:60]])
                l = uhandle.readline().rstrip()
                l = uhandle.readline().rstrip()
                consumer.jacknife_bottom(
                        filter(lambda x: x<>'' and x<>'*' and not x.startswith('J') and not x.startswith('o') and not x.startswith('L'),
                            [l[30:40], l[40:50], l[50:60]]
                        )
                )
                l = uhandle.readline().rstrip()
                consumer.jacknife_top(
                        filter(lambda x: x<>'',
                            [l[30:40], l[40:50], l[50:60]]
                        )
                )
            l = uhandle.readline()
        consumer.end_record()

    def fetch_values(self, hook, l, fun):
        l = l.rstrip()
        if (re.search(hook, l) <> None):
            fun(filter(lambda x: x<>'',
                [l[30:40], l[40:50], l[50:60]]))
        

class _RecordConsumer(AbstractConsumer):
    """Consumer that converts a LDNe record to a Record object.

    Members:
    data    Record with LDNe data.

    """
    def __init__(self):
        self.data = None
        self.pop = None

    def start_record(self):
        self.data = Record()

    def attach_pop(self):
        for i in range(len(self.data.freqs_used)):
            self.pop.append((
                    self.hm[i],
                    self.ic[i],
                    self.or2[i],
                    self.er2[i],
                    (self.ne[i],
                        (self.pb[i], self.pt[i]),
                        (self.jb[i], self.jt[i]))
            ))
        self.data.populations.append((self.pop_id, self.pop))

    def end_record(self):
        if self.pop <> None: self.attach_pop()

    def mating(self, mating):
        self.data.mating = mating

    def cutoff(self, cutoff):
        self.data.freqs_used = max(lambda x: float(x) ,cutoff)

    def start_pop(self, id):
        if self.pop <> None: self.attach_pop()
        self.pop    = []
        self.pop_id = id

    def harmonic_mean(self, means):
        self.hm = map(lambda x: float(x), means)

    def independent_comparison(self, ic):
        self.ic = map(lambda x: int(x), ic)

    def overall_r2(self, or2):
        self.or2 = map(lambda x: float(x), or2)

    def expected_r2(self, er2):
        self.er2 = map(lambda x: float(x), er2)

    def estimated_ne(self, ne):
        self.ne = map(lambda x: float(x), ne)

    def parametric_bottom(self, pb):
        self.pb = map(lambda x: float(x), pb)

    def _number_or_inf(self, x):
        try:
            return float(x)
        except ValueError:
            return float(decimal.Inf)
    def parametric_top(self, pt):
        self.pt = map(lambda x: self._number_or_inf(x) , pt)

    def jacknife_bottom(self, jb):
        #print "AAAA"
        self.jb = map(lambda x: float(x), jb)

    def jacknife_top(self, jt):
        #print "AbAA"
        self.jt = map(lambda x: self._number_or_inf(x), jt)
