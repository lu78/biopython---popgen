import os

directory = None

def prepare_dir():
    try:
        os.mkdir(directory)
    except OSError:
        pass #Probably OK, dir already exists

def init(hapmap_dir = '.'):
    '''Initializes the HapMap module.

    It doesn't have to be called directly, as long as any init method
    in the module is called, this root one will be indirectly called also.
    '''
    global directory
    directory = hapmap_dir
    prepare_dir()
