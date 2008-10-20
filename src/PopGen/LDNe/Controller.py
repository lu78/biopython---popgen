# Copyright 2008 by Tiago Antao <tiagoantao@gmail.com>.  All rights reserved.

"""
This module allows to control LNDe.

"""

import os
import tempfile
from shutil import copyfile
from logging import debug

class LDNeController:
    def __init__(self, ldne_dir):
        """Initializes the controller.
        
        ldne_dir is the directory where LDNe is.
        
        The initializer checks for existance and executability of binaries.
        """
        self.ldne_dir = ldne_dir
        self.os_name = os.name
        if self.os_name=='nt':
            self.bin_name = 'LDNe.exe'
            #this is wrong (the exe name), most probably
        else:
            self.bin_name = 'ldne'
        dir_contents = os.listdir(self.ldne_dir)
        if self.bin_name in dir_contents:
            if not os.access(self.ldne_dir + os.sep +
                self.bin_name, os.X_OK):
                raise IOError, "LDNe not executable"
        else:
            raise IOError, "LDNe not available"

    def run_ldne(self, gen_file, out_file):
        """Executes LDNe.
        """
        in_name = 'in.le'
        inf = open(in_name,'w')
        inf.write(gen_file + "\n")
        inf.write(out_file + "\n")
        inf.close()
        os.system(self.ldne_dir + os.sep + self.bin_name + 
          ' < ' + in_name + ' >/dev/null 2>&1')
    

