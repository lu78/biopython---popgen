#!/usr/bin/env python
# biopython license
"""

Parses files in HGDP (I don't know how to call it) format.

Example:
>>> import StringIO
>>> hgdp_file = StringIO('''
... [Header]
... BSGT Version    2.3.41.16318
... Processing Date    8/20/2006 1:40 PM
... Content        BDCHP-1X12-HUMANHAP650Y_11226216_A.csv
... Num SNPs    655352
... Total SNPs    655352
... Num Samples    30
... Total Samples    159
... [Data]
...     WG0001019-DNAA03_NA18547    WG0001019-DNAD03_NA18593    WG0001019-DNAE03_NA18550    WG0001019-DNAF02_NA18577    WG0001019-DNAF03_NA18542    WG0001019-DNAG02_NA18537    WG0001019-DNAH02_NA18529    WG0001019-DNAH03_NA18576    WG0001019-DNAB03_NA18558    WG0001019-DNAC03_NA18612
... MitoA10045G    AA    AA    AA    AA    AA    AA    AA    AA    AA    AA   
... MitoA10551G    AA    AA    AA    AA    AA    AA    AA    AA    AA    AA
... MitoA11252G    AA    AA    AA    AA    AA    AA    AA    AA    AA    AA
... MitoA11468G    AA    AA    AA    AA    AA    AA    AA    AA    AA    AA
... MitoA11813G    AA    AA    AA    AA    AA    AA    AA    AA    AA    AA
... MitoA12309G    AA    AA    AA    AA    AA    AA    AA    AA    AA    AA
... MitoA13106G    AA    AA    AA    AA    AA    AA    AA    AA    AA    AA
... MitoA13264G    GG    AA    AA    AA    GG    AA    AA    AA    AA    AA
... MitoA13781G    AA    AA    AA    AA    AA    AA    --    AA    AA    AA
... MitoA14234G    AA    AA    AA    AA    AA    AA    AA    AA    AA    AA
... MitoA14583G    AA    AA    AA    AA    AA    AA    AA    AA    AA    AA
... MitoA14906G    GG    GG    GG    GG    GG    GG    GG    GG    GG    GG
... MitoA15219G    AA    AA    AA    GG    AA    AA    AA    AA    AA    AA
... MitoA15245G    AA    AA    AA    AA    AA    AA    AA    AA    AA    AA
... MitoA15302G    AA    AA    AA    AA    AA    GG    GG    AA    AA    GG
... MitoA15759G    AA    AA    AA    AA    AA    AA    AA    AA    AA    AA
... MitoA15908G    AA    AA    AA    AA    AA    AA    AA    AA    AA    AA
... MitoA15925G    AA    AA    AA    AA    AA    AA    AA    AA    AA    AA
... MitoA16163G    AA    AA    AA    AA    AA    AA    AA    AA    AA    AA

... ''')
"""