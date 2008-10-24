#!/usr/bin/env python
# Marker Object (a SNP, a gene, etc..)

from PopGen.Gio.PopGenExceptions import InvalidGenotype

class Marker(object):
    '''
    A Marker object(like a SNP, a gene, etc...)
    # Note: would it be better to use 'Marker' or 'Locus'?
    
    >>> C10G = Marker('A130G')
    >>> C10G.genotypes = (('A', 'A'), ('G', 'A'))    # should be done via add_genotype
    >>> C10G.genotypes.append(('A', 'A'))
    Traceback (most recent call last):
        ... 
    AttributeError: 'tuple' object has no attribute 'append'
    '''      
          
    def __init__(self, name = None):
        # should add a check for parameters type (e.g. missing -> should be integer)
        if name is None:
            name = 'Un-named Marker'
        self.name = name
        self.position = ""  # should be a 'position' object. For now, just a description (e.g. chromosome 11 pos 23131)
        self.genotypes = () # list of genotypes object (e.g.: [('A', 'A'), ('G', 'A')]) Should be an object
        self.original_strand = ''    # should be '+' or '-'
        self.references = ''         # gene name, associated diseases, etc..
        
    def _check_genotype_input(self, genotype):
        """
        check if a genotype input is correct (a tuple with two elements)
        
        >>> m = Marker()
        >>> m._check_genotype_input(['A', 'C'])
        Traceback (most recent call last):
            ... 
        InvalidGenotype: not a tuple of two elements
        >>> m._check_genotype_input(('A'))
        Traceback (most recent call last):
            ... 
        InvalidGenotype: not a tuple of two elements
        >>> m._check_genotype_input(('A', 'C'))
        """
        if not isinstance(genotype, tuple):
            raise InvalidGenotype('not a tuple of two elements')
        elif len(genotype) != 2:
            raise InvalidGenotype('not a tuple of two elements') 
        
    def add_multiple_genotypes(self, genotypes):
        """
        Please use this method when adding multiple genotypes (e.g. when parsing a file), 
        instead of overwriting Marker.genotypes
        """
        for genotype in genotypes:
            self._check_genotype_input(genotype)
        self.genotypes = genotypes  # avoid calling self.genotype multiple times        
        self._recalculate_freqs()
        
    def add_genotype(self, genotype):
        """
        use this method when adding a genotype.
        """
        # should check genotype format here
        
        self.genotypes = (self.genotypes, genotype)
        
        self._recalculate_freqs() 
        
    def _recalculate_freqs(self):
        """
        recalculate allele frequency, total population, missing_data_count, etc...
        """ 
        
def _test():
    import doctest
    doctest.testmod(verbose = False)
    
if __name__ == '__main__':
    _test()
     
        
    