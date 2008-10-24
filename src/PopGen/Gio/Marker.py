#!/usr/bin/env python
# Marker Object (a SNP, a gene, etc..)

class wrong_add_individual_input(Exception): 
    def __string__(self): "Wrong individual input"
    

class Marker(object):
    '''
    A Marker object(like a SNP, a gene, etc...)
    # Note: would it be better to use 'Marker' or 'Locus'?
    
    >>> C10G = Marker('A130G')
    >>> C10G.genotypes = [('A', 'A'), ('G', 'A')])
    >>> C10G.genotypes.append(('A', 'A'))

    '''      
          
    def __init__(self, name = None):
        # should add a check for parameters type (e.g. missing -> should be integer)
        if name is None:
            name = 'Un-named Marker'
        self.name = name
        self.position = ""  # should be a 'position' object. For now, just a description (e.g. chromosome 11 pos 23131)
        self.genotypes = [] # list of genotypes object (e.g.: [('A', 'A'), ('G', 'A')])
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
     
        
    