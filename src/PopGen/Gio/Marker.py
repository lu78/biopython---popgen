#!/usr/bin/env python
# Marker Object (a SNP, a gene, etc..)

class wrong_add_individual_input(Exception): 
    def __string__(self): "Wrong individual input"
    

class Marker():
    '''
    A Marker object(like a SNP, a gene, etc...)
    # Note: would it be better to use 'Marker' or 'Locus'?
    
    >>> C10G = Marker('A130G')
    >>> C10G.genotypes = (('A', 'A'), ('G', 'A'))    # should be done via add_genotype
    >>> C10G.genotypes.append(('A', 'A'))
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
     
        
    