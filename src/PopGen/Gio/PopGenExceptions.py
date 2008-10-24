#!/usr/bin/env python
# Copyright 2008 by Giovanni Dall'Olio.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

class GenericPopGenException(Exception):
    def __init__(self, error_message):
        self.error_message = error_message
        
    def __repr__(self):
        return self.error_message

class InvalidInputFile(GenericPopGenException):
    pass
        