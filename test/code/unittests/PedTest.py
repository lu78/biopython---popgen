#!/usr/bin/env python
# test unit for Ped Parser

import unittest
from pprint import pprint
import os

import Ped

class PedTest(unittest.TestCase):
    """A test class for the Ped parser module"""
    
    def setUp(self):
        self.testfiles = os.listdir('test/data/ped')