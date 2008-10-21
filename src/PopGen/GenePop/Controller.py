# Copyright 2008 by Tiago Antao <tiagoantao@gmail.com>.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.



"""
This module allows to control GenePop.

"""

import os
import tempfile
import re
#from sys import platform, maxint
#from shutil import copyfile
#from random import randint, random
#from time import strftime, clock
#from logging import debug

class GenePopController:
    def __init__(self, genepop_dir = '', ext = None):
        """Initializes the controller.
        
        genepop_dir is the directory where GenePop is.
        ext is the extension of binaries (.exe on windows, 
          none on Unix)

        The binary should be called Genepop (capital G)
        
        """
        self.tmp_idx = 0
        self.genepop_dir = genepop_dir
        self.os_name = os.name
        if self.os_name == 'nt':
            py_ext = '.exe'
        else:
            py_ext = ''
        if ext is None:
            self.ext = py_ext
        else:
            self.ext = ext

    def _get_path(self):
        """Returns the path to the GenePop application.

           Includes path where GenePop can be found plus executable extension.
        """
        if self.genepop_dir == '':      # if genepop_dir is not defined, return what?
            return app + self.ext       # undefined variable app
        else:
            return os.sep.join([self.genepop_dir, 'Genepop']) + self.ext

    def calc_fst_all(self, fname):
        """Executes GenePop and gets Fst/Fis/Fit (all populations)
        
        Parameters:
        fname - file name

        Returns:
        List of tuples
          (Locus name, Fis, Fst, Fit)
        
        Creates a file named fname.FST
        """
        os.system(
            self._get_path() + ' Mode=Batch MenuOptions=6.1 InputFile='+fname +
            '> /dev/null 2> /dev/null')
        f = open(fname+'.FST')
        fsts = []
        l = f.readline()
        while l != '':
            l = l.rstrip()
            if l.startswith('  Locus:'):
                locus = l.split(':')[1].lstrip()
            elif l.startswith('Fis^='):
                fis = l.split(' ')[1]
            elif l.startswith('Fst^='):
                fst = l.split(' ')[1]
            elif l.startswith('Fit^='):
                fit = l.split(' ')[1]
                fsts.append((locus, fis, fst, fit))
            l = f.readline()
        f.close()
        return fsts

    def get_loci_genotype_counts(self, fname):
        """Executes GenePop and gets genotype counts (all populations)
        
        Parameters:
        fname - file name

        Returns:
        List (1 item per pop) of List (1 item per locus) of tuples:
          (Locus name, [(a1,a2,count), ...], expHomo, obsHomo, expHet, obsHet)
          locus might be repeated (once per subpop)
        
        Creates a file named fname.FST
        """
        os.system(
            self._get_path() + ' Mode=Batch MenuOptions=5.1 InputFile='+fname +
            '> /dev/null 2> /dev/null') # this would only work on Unix?
        locus = None
        onInfo = False
        f = open(fname+'.INF')
        l = f.readline()
        doneLocus = []
        counts = []
        allPopCounts = []
        
        while l != '':
            l = l.rstrip()
            match = re.match(".*Pop: .* Locus: (.+)", l)
            
            if match is not None:
                locus = match.group(1)
                if locus in doneLocus:
                    doneLocus=[locus]
                    allPopCounts.append(counts)
                    counts = []
                else:
                    doneLocus.append(locus)
                genoCounts = []
                
            if locus is not None:
                if l.find("Genotypes  Obs.")>-1:
                    onInfo = True
                elif onInfo:
                    m2 = re.match(" +([0-9]+) , ([0-9]+) * ([0-9]+)",l)
                    if m2 is not None:
                        genoCounts.append((m2.group(1), m2.group(2), m2.group(3)))
                    else:
                      onInfo = False
                elif l.find("Expected number of ho")>-1:
                    expHo =  float(l[38:])
                elif l.find("Observed number of ho")>-1:
                    obsHo =  int(l[38:])
                elif l.find("Expected number of he")>-1:
                    expHe =  float(l[38:])
                elif l.find("Observed number of he")>-1:
                    obsHe =  int(l[38:])
                    counts.append((locus, genoCounts, expHo, obsHo, expHe, obsHe))
                    locus = None

            l = f.readline()
        f.close()
        allPopCounts.append(counts)
        return allPopCounts

