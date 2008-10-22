#!/usr/bin/env python
# test unit for Ped Parser

import unittest
from pprint import pprint
import os
import logging

from PopGen.Gio import Ped

class PedParserTest(unittest.TestCase):
    """A test class for the Ped parser module"""
    
    def setUp(self):
        self.parser = Ped.RecordParser() 
        self.testfilespath = '../../data/ped/'
        self.testfiles = os.listdir(self.testfilespath)
        
    def testParser(self):
        for filename in self.testfiles:
            pop = self.parser.parse(file(self.testfilespath + filename, 'r'))
            # how to test now?
    
    def tearDown(self):
        """
        tear down any data used in tests
        tearDown is called after each test function execution.
        """

        pass
        
        
if __name__ == '__main__':
    logging.basicConfig(level = logging.DEBUG)
    unittest.main()
    