#!/usr/bin/env python
# Population object
class Continent(object):
    """
    >>> Europe = Continent("Europe")
    """
    def __init__(self, name):
        self.name = name

class Population(object):
    """
    Describe a Population objects.    
    
    >>> Vulcanians = Population("Vulcanians")
    >>> from Marker import Marker        # I wrote this when I was handling PED files
    >>> cox2 = Marker('cox2')
    >>> cox2.heterozygotes = 6
    >>> cox2.total_population = 30
    >>> Vulcanians.markers.append(cox2)
    
    >>> pops = (Population("Brahui", "Pakistan", "Asia", "Brahui"), 
    ...         Population("French_Basque", "France", "Europe", "Basque"),
    ...          )
    >>> Europe = Continent("Europe")
    >>> print [pop for pop in pops if pop == Europe]
    [French_Basque]
    """
    
    def __init__(self, name, region = None, continent = None, working_unit = None):
        self.name = name
        self.region = region
        self.continent = continent
        self.working_unit = working_unit
        
        self.markers = []       # not sure I will put markers here
        self.individuals = []   # not sure I will use this
        
    def __repr__(self):
        return "%s" % (self.name)
    
    def __eq__(self, other):
        if isinstance(other, Continent):
            return self.continent == other.name
        # TODO: implement other cases
        else:
            return self.name == other
    
def _test():
    import doctest
    doctest.testmod(verbose = False)
    
if __name__ == '__main__':
    _test()
    
