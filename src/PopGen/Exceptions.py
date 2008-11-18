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
        

class RequiresGenotypeException(GenericPopGenException): 
    pass