#!/usr/bin/env python
# Population object

class Population:
    """
    Describe a Population objects.    
    
    >>> Vulcanians = Population("Vulcanians")
    >>> from Marker import Marker
    >>> cox2 = Marker('cox2')
    >>> cox2.heterozygotes = 6
    >>> cox2.total_population = 30)
    >>> Vulcanians.markers.append(cox2)
    
    """
    name = ''
    markers = []
    
    def __init__(self, population_name):
        self.name = population_name
        
    def __repr__(self):
        return "%s" % (self.name)
    
    
def _test():
    import doctest
    doctest.testmod(verbose = True)
    
if __name__ == '__main__':
    _test()
    
