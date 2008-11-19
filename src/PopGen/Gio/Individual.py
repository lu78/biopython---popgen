class Individual(object):
    """
    Represents an Individual
    It is suggested to instantiate Individuals with the .from_ped_line method.
    
    >>> Einstein = Individual().from_ped_line(["Vulcanians", "Einstein", "0", "0", "1", "2", "C", "T", "C", "C"])
    >>> Einstein
    Mr. Einstein (Vulcanians)
    """    
    
    def __init__(self, id):
        self.individual_id = id
        self.population = 'unkn. population'
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
    
    def from_ped_line(self, ped_line):
        self.population = ped_line[0]
        self.individual_id = ped_line[1]
        self.father = ped_line[2]
        self.mother = ped_line[3]
        self.sex = ped_line[4]
        return self
    
def _test():
    import doctest
    doctest.testmod(verbose = True)
    
if __name__ == '__main__':
    _test()