#!/usr/bin/env python
# Marker Object (a SNP, a gene, etc..)

class wrong_add_individual_input(Exception): 
    def __string__(self): "Wrong individual input"
    

class Marker(object):
    '''
    A Marker object(like a SNP, a gene, etc...)
    # Note: would it be better to use 'Marker' or 'Locus'?
    
    >>> C10G = Marker('A130G')
    >>> C10G.add_individual(1, 1, 0) # heterozygote individual
    >>> C10G.add_individual(2, 0, 0) # individual with two G or A
    >>> C10G.add_individual(0, 0, 2) # individual with missing data
    >>> C10G.heterozygotes
    1
    >>> C10G.Pur_count
    3
    >>> C10G.missing_data_count
    2
    >>> C10G.total_population
    2
    '''      
          
    def __init__(self, name = None):
        # should add a check for parameters type (e.g. missing -> should be integer)
        if name is None:
            name = 'Un-named Marker'
        self.name = name
        self.position = ""  # should be a 'position' object. For now, just a description (e.g. chromosome 11 pos 23131)
        self.reference_allele_freq = 0.0
        self.derived_allele_freq = 0.0
        self.minor_allele_freq = 0.0
        self.original_strand = ''    # should be '+' or '-'
        self.references = ''         # gene name, associated diseases, etc..
        
        
    def add_individual(self, individual):
        '''
        '''
        pass
        
        
def _test():
    import doctest
    doctest.testmod(verbose = False)
    
if __name__ == '__main__':
    _test()
     
        
    