#!/usr/bin/env python
# Population object

class Population(object):
    """
    Describe a Population objects.    
    
    >>> Vulcanians = Population("Vulcanians")
    >>> from Marker import Marker        # I wrote this when I was handling PED files
    >>> cox2 = Marker('cox2')
    >>> cox2.heterozygotes = 6
    >>> cox2.total_population = 30
    >>> Vulcanians.markers.append(cox2)
    
    >>> Vulcanians.region = ''
    >>> Vulcanians.continent = ''
    
    """
    
    def __init__(self, name, region = None, continent = None, unit = None):
        self.name = name
        self.region = region
        self.continent = continent
        self.unit = unit
        self.markers = []       # not sure I will put markers here
        self.individuals = []   # not sure I will use this
        
    def __repr__(self):
        return "%s" % (self.name)
    
    
def _test():
    import doctest
    doctest.testmod(verbose = False)
    
if __name__ == '__main__':
    _test()
    
