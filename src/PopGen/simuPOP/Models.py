#2007 by Tiago Antao <tiagoantao@gmail.com>. All rights reserved.

from simuPOP import *

"""
Provides data structures for several demographic models.

"""

def _prepareArgs(args, kwargs):
    """Prepares argument lists for use (reporting size).
    """
    kwargs.pop('subPop', 0)
    for i in range(len(args), 0, -1):
        pos = i-1
        name = args[pos]
        if name == 'subPop':
            del args[pos]

def island(demes, mig, size, *args, **kwargs):
    """Creates an island model representation.

    Parameters are the same for the population object plus number
        of demes, migration rate and less subPop.
    Returns a population AND an operator.
    """
    _prepareArgs(args, kwargs)
    pop = population(size = 0, subPop = demes * [size], *args, **kwargs)
    migMat = []
    for i in range(demes):
        migLine = []
        for j in range(demes):
            if i == j: migLine.append(0.0)
            else: migLine.append(mig)
        migMat.append(migLine)
    op = migrator(rate = migMat)
    return pop, op

def ssm1D(demes, mig, size, *args, **kwargs):
    """Creates a 1D stepping stone model representation.

    Parameters are the same for the population object plus number
        of demes, migration rate and less subPop.
    Returns a population AND an operator.
    """
    _prepareArgs(args, kwargs)
    pop = population(size = 0, subPop = demes * [size], *args, **kwargs)
    migMat = []
    for i in range(demes):
        migLine = []
        for j in range(demes):
            if abs(i, j) == 1: migLine.append(mig)
            else: migLine.append(0.0)
        migMat.append(migLine)
    op = migrator(rate = migMat, begin = 1)
    return pop, op

def ssm2D(xSize, ySize, mig, size, *args, **kwargs):
    """Creates a 2D stepping stone model representation.

    Parameters are the same for the population object plus x and y
        dimension sizes, migration rate and less subPop.
    Returns a population AND an operator.
    """
    size = _prepareArgs(args, kwargs)
    pop = population(size = 0, subPop = (xSize*ySize) * [size], *args, **kwargs)
    migMat = []
    for y in range(ySize):
        for x in range(xSize):
            migLine = []
            for y2 in range(ySize):
                for x2 in range(xSize):
                    if abs(x, x2) == 1 and y == y2: migLine.append(mig)
                    if abs(y, y2) == 1 and x == x2: migLine.append(mig)
                    else: migLine.append(0.0)
            migMat.append(migLine)
    op = migrator(rate = migMat)
    return pop, op
