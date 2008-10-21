class Individual(object):
    """
    Represents an Individual
    It is suggested to instantiate Individuals with the .from_ped_line method.
    
    >>> Einstein = Individual().from_ped_line(["Vulcanians", "Einstein", "0", "0", "1", "2", "C", "T", "C", "C"])
    >>> Einstein
    Mr. Einstein (Vulcanians)
    """
    
    population = ''
    name = ''
    father = ''
    mother = ''
    sex = '0'
    phenotype = '0'
    markers = []
    
    def __init__(self):
        pass
    
    def __repr__(self):
        if self.sex == '0' or self.sex == '1':
            r = "Mr. %s (%s)" %(self.name, self.population)
        else:
            r = "Mrs. %s (%s)" %(self.name, self.population)
        return r
    
    def from_ped_line(self, ped_line):
        self.population = ped_line[0]
        self.name = ped_line[1]
        self.father = ped_line[2]
        self.mother = ped_line[3]
        self.sex = ped_line[4]
        return self
    
def _test():
    import doctest
    doctest.testmod(verbose = True)
    
if __name__ == '__main__':
    _test()