class Individual(object):
    """
    Represents an Individual
    It is suggested to instantiate Individuals with the .from_ped_line method.
    
    >>> Einstein = Individual().from_ped_line(["Vulcanians", "Einstein", "0", "0", "1", "2", "C", "T", "C", "C"])
    >>> Einstein                              # Test __repr__ method
    Mr. Einstein (Vulcanians)
    >>> print Einstein + ' Albert'            # Test __add__ method
    Einstein Albert
    >>> print Einstein in ('Einstein', )      # Test __eq__ method
    True
    """    
    
    def __init__(self, id = None, population = 'unkn. population',
                 region = None, continent = None, area = None, working_unit = None, 
                 father = None, mother = None, sex = '0', phenotype = None, markers = []):
        if id is not None:
            self.individual_id = id
        else:
            self.individual_id = None       # is this ok?
        self.population = population    # better a 'population' object?
        self.region = region
        self.continent = continent
        self.area = area
        self.working_unit = working_unit
        self.father = ''
        self.mother = ''
        self.sex = '0'
        self.phenotype = '0'
        self.markers = []
    
    def __repr__(self):
        if self.sex in ('0', '1'):
            r = "Mr. %s (%s)" %(self.individual_id, self.population)
        else:
            r = "Mrs. %s (%s)" %(self.individual_id, self.population)
        return r
    
    def __str__(self):
        """
        """
        return self.individual_id
    
    def __add__(self, other):
        return str(self.individual_id) + other
    
    def __eq__(self, other):
        return self.individual_id == other
    
    def __ne__(self, other):
        return self.individual_id != other
    
    def from_ped_line(self, ped_line):
        self.population = ped_line[0]
        self.individual_id = ped_line[1]
        self.father = ped_line[2]
        self.mother = ped_line[3]
        self.sex = ped_line[4]
        return self
    
def _test():
    import doctest
    doctest.testmod(verbose = False)
    
if __name__ == '__main__':
    _test()