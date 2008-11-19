class GenericPopGenException(Exception):
    def __init__(self, error_message):
        self.error_message = error_message
        return
    
    def __str__(self):
        return self.error_message

class InvalidInputFile(GenericPopGenException):
    pass

class InvalidGenotype(GenericPopGenException):
    pass

class NotImplementedError(GenericPopGenException): pass

class PopulationExistsException(GenericPopGenException):
    """A certain population already exists.
    """

    def __init__(self, pop_name):
        self.pop_name = pop_name

    def __str__(self):
        return self.pop_name + ' already exists'
    
class RequiresGenotypeException(GenericPopGenException):
    pass
