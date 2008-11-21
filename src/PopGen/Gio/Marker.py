#!/usr/bin/env python
# Marker Object (a SNP, a gene, etc..)

from PopGen.Gio.PopGenExceptions import InvalidGenotype

        
def _check_genotype_input(genotype):
    """
    check if a genotype input is correct (a tuple with two elements)
    
    >>> m = Marker()
    >>> _check_genotype_input(['A', 'C'])
    Traceback (most recent call last):
        ... 
    InvalidGenotype: not a tuple of two elements
    >>> _check_genotype_input(('A'))
    Traceback (most recent call last):
        ... 
    InvalidGenotype: not a tuple of two elements
    >>> _check_genotype_input(('A', 'C'))
    """
    if not isinstance(genotype, tuple):
        raise InvalidGenotype('not a tuple of two elements')
    elif len(genotype) != 2:
        raise InvalidGenotype('not a tuple of two elements') 

class Marker(object):
    '''
    A Marker object(like a SNP, a gene, etc...)
    
    Attributes:
    - name                    -> name/id of the marker (e.g. sn1334)
    - position                -> a description of the position, e.g. 'chromosome 11 pos. 12331'
    - genotypes               -> list of genotypes (e.g.: (('A', 'A'), ('None', 'None')))
    - reference_allele_freq
    - derived_allele_freq  
    - minor_allele_freq
    - total_individual_count
    
    >>> C10G = Marker('A130G')
    >>> C10G.genotypes = (('A', 'A'), ('G', 'A'))    # should be done via add_genotype
    >>> C10G.genotypes.append(('A', 'A'))
    Traceback (most recent call last):
        ... 
    AttributeError: 'tuple' object has no attribute 'append'
    >>> C10G.populations = {'Vesuvians': (0, 1), 'Martians': (2)}
    '''      
          
    def __init__(self, name = None, individuals = None):
        # should add a check for parameters type (e.g. missing -> should be integer)
        if name is None:
            name = 'Un-named Marker'
        self.name = name
        self.position = ""  # should be a 'position' object. For now, just a description (e.g. chromosome 11 pos 23131)
        self.genotypes = () # list of genotypes object (e.g.: (('A', 'A'), ('G', 'A'))) Should be an object
        self.populations = {}       # should define populations. which are the positions in self.genotypes which correspond to populations.
        self.individuals = individuals 
        self.individual_count = 0        # Individuals for which the Marker is genotyped
        self.missing_data_count = 0      # Individuals for which data is not available
        self.reference_allele_freq = 0.0 # frequency of the reference allele (the first!)
        self.derived_allele_freq = 0.0   # frequency of the derived allele (the second!)
        self.original_strand = ''        # should be '+' or '-'
        self.references = ''             # gene name, associated diseases, etc..
    
    def __repr__(self):
        return "Marker %s, %s individuals" % (self.name, self.individual_count)

    def to_geno_format(self):
        """
        Temporary function to print a marker in .geno (e.g. HGDP) format
        """
        output = self.name + '\t'
        for genotype in self.genotypes:
            output += ''.join(genotype) + '\t'
#        output += '\n'
        return output
        
    def add_multiple_genotypes(self, genotypes):
        """
        Please use this method when adding multiple genotypes (e.g. when parsing a file), 
        instead of overwriting Marker.genotypes
        """
        for genotype in genotypes:
            _check_genotype_input(genotype)
        self.genotypes = genotypes  # avoid calling self.genotype multiple times        
        self._recalculate_freqs()
        
    def add_genotype(self, genotype):
        """
        use this method when adding a genotype.
        """
        # should check genotype format here
        if isinstance(genotype, basestring) and len(genotype) == 2:
            # genotype is probably a string like 'AA'
            if self.genotypes == ():  # first genotype added
                self.genotypes = (genotype[0], genotype[1])
            else:
                self.genotypes = (self.genotypes, (genotype[0], genotype[1]))
        else:        
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
     
        
    