import csv
from urllib2 import urlopen
import gzip
import os
from Bio.PopGen import HapMap


def get_parent_info():
    '''Get parent info.
       The HapMap module has to be initialized
    '''
    file_location = HapMap.directory + os.sep + "parent.csv"
    if (not os.path.exists(file_location)):
        httpDir = 'http://www.hapmap.org/downloads/samples_individuals'
        f = urlopen(httpDir + '/' + 'pedinfo2sample_YRI.txt.gz')
        file(HapMap.directory + os.sep + 'pedYRI.gz', 'wb').write(f.read())
        f.close()
        f = urlopen(httpDir + '/' + 'pedinfo2sample_CEU.txt.gz')
        file(HapMap.directory + os.sep + 'pedCEU.gz', 'wb').write(f.read())
        f.close()
        parentf = open(file_location, "wb")
        gf = gzip.open(HapMap.directory + os.sep + 'pedYRI.gz')
        parentf.write("YRI\n")
        parentf.write(gf.read())
        gf.close()
        gf = gzip.open(HapMap.directory + os.sep + 'pedCEU.gz')
        parentf.write("CEU\n")
        parentf.write(gf.read())
        gf.close()
        parentf.close()
        

    reader = csv.reader(open(file_location, "rb"), delimiter='\t')
    hasFamily = 0
    parents = {}
    for row in reader:
        if len(row) == 1:
            pop=row[0]
            continue
        if int(row[3]) > 0: # is offspring
            ofs = row[-1].split(':')[-2]
        elif row[4] == '2': #father
            father = row[-1].split(':')[-2]
        else:
            mother = row[-1].split(':')[-2]
        hasFamily += 1
        if hasFamily==3:
            hasFamily = 0
            parents[pop, ofs] = (father, mother)
    return parents
