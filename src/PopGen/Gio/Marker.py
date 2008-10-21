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
    
    heterozygotes_count = 0
    Pur_count = 0
    Pyr_count = 0
    missing_data_count = 0
    total_population = 0 # people with missing data are excluded from this count
      
    def __init__(self, name = None):
        # should add a check for parameters type (e.g. missing -> should be integer)
        if name is None:
            name = 'Un-named Marker'
        self.name = name
        
    def add_individual(self, Pur_count, Pyr_count, missing_data_count):
        '''
        '''
        
        if Pur_count + Pyr_count + missing_data_count != 2:
            raise wrong_add_individual_input()
        
        if missing_data_count != 1:
            self.missing_data_count += missing_data_count
            
        else:
            self.total_population += 1
            
            self.Pur_count += Pur_count
            self.Pyr_count += Pyr_count
            
            if Pyr_count == Pur_count == 1:
                # Current individual is heterozygote.
                self.heterozygotes_count += 1  
        
        
def _test():
    import doctest
    doctest.testmod(verbose = False)
    
if __name__ == '__main__':
    _test()
     
        
    