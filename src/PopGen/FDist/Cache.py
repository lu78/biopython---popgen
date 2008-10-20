# Copyright 2007 by Tiago Antao <tiagoantao@gmail.com>.  All rights reserved.

"""
This module allows to cache fdist results, and return on the fly
in case the calculation was done.

http://www.rubic.rdg.ac.uk/~mab/software.html
"""

from bz2 import BZ2File
from logging import debug
from math import log10
import os
import tempfile
from PopGen.FDist.Controller import FDistController

class Cache:
    def __init__(self, fdist_dir, cache_dir, limit=0.001):
        """Initializes the cache.
        
        cache_dir - Where the cache will be stored.
        temp_dir - Where temporary calculations will be done.
        """
        self.fdist_dir = fdist_dir
        self.cache_dir = cache_dir
        self.limit = limit

    def run_fdist_force_fst(self, npops, nsamples, fst, sample_size,
        mut = 0, num_sims = 50000, data_dir='.', try_runs = 5000):
        """Executes fdist (if not in cache).
        
        Parameters
        npops - Number of populations
        nsamples - Number of populations sampled
        fst - expected Fst
        sample_size - Sample size per population
        mut - 0=Stepwise, 1=Infinite allele
        num_sims - number of simulations
        dir - directory where fdist will be executed (must be rw)
        try_runs - number of simulations on the part trying to get
                   Fst correct.        
        Important Note: This can take quite a while to run!
        """
        cache_fst = str(round(fst, -int(log10(self.limit))))
        cache_name = '_'.join([str(npops), str(nsamples), str(sample_size),
            str(mut), cache_fst]) + '.bz2'
        full_cache_name = self.cache_dir + os.sep + cache_name
        full_final_name = data_dir + os.sep + 'out.dat'
        
        if os.access(full_cache_name, os.R_OK):
            bzf = BZ2File(full_cache_name, 'r')
            of = open(full_final_name, 'w')
            of.writelines(bzf.readlines())
            of.close()
            bzf.close()
            return float(cache_fst)
        else:
            fdc = FDistController(self.fdist_dir)
            real_fst = fdc.run_fdist_force_fst(npops, nsamples, fst,
                sample_size, mut, num_sims, data_dir, try_runs, self.limit)
            cache_name = '_'.join([str(npops), str(nsamples), str(sample_size),
                str(mut), str(round(real_fst, -int(log10(self.limit))))]) \
                + '.bz2'
            full_cache_name = self.cache_dir + os.sep + cache_name
            if not os.access(full_cache_name, os.R_OK):
                bzf = BZ2File(full_cache_name, 'w')
                inf = open(data_dir + os.sep + 'out.dat')
                bzf.writelines(inf.readlines())
                inf.close()
                bzf.close()
            if real_fst < fst:
                begin = real_fst
                end   = fst
            else:
                begin = fst
                end   = real_fst
            while begin <= end + (self.limit/2):
                interval_full_cache_name = self.cache_dir + os.sep + \
                    '_'.join([str(npops), str(nsamples), str(sample_size),
                    str(mut), str(round(begin, -int(log10(self.limit))))]) \
                    + '.bz2'
                try:
                    os.symlink(full_cache_name, interval_full_cache_name)
                except OSError:
                    pass #safe to ignore
                begin += self.limit
            return real_fst

    def run_datacal(self, data_dir='.'):
        fdc = FDistController(self.fdist_dir)
        return fdc.run_datacal(data_dir)

    def run_cplot(self, ci = 0.95, data_dir='.'):
        fdc = FDistController(self.fdist_dir)
        return fdc.run_cplot(ci, data_dir)

    def run_pv(self, out_file='probs.dat', data_dir='.'):
        fdc = FDistController(self.fdist_dir)
        return fdc.run_pv(out_file, data_dir)

