#!/usr/bin/env python
# Copyright 2008 by Giovanni Dall'Olio.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
PED file parser.

The format is described here:
- http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped

Every line in the PED file corresponds to an individual. 

Example of PED file:
# population id father mother char1 phenotype locus1_allele1 locus1_allele2 locus2_allele1 .....
Mandenka HGDP00912 0 0 1 2 C C G A G G T T T T A A G G
Mandenka HGDP01283 0 0 1 2 T C G A G G C T C C G G G A
Yoruba HGDP00928 0 0 2 2 C C G G G G C T T T G A G G
Yoruba HGDP00937 0 0 1 2 C C A A A G C C T C A A G A


Classes:

Record            Holds PED file data
RecordParser      Parses a PED record (file) into a Record object.

_Scanner         Scans a PED record.
_RecordConsumer  Consumes PED data to a Record object.


"""


from copy import deepcopy
from types import *

from Bio import File
from Bio.ParserSupport import *

class Record:
    pass

