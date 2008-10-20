from ftplib import FTP
import gzip
import os
from Bio.PopGen import HapMap
from Bio.PopGen.Stats import Simple

def prepare_geno_dir():
    global geno_dir
    geno_dir = HapMap.directory + os.sep +  'genotypes'
    try:
        os.mkdir(geno_dir)
    except OSError:
        pass #Probably OK, dir already exists

def init(hapmap_dir = '.'):
    HapMap.init(hapmap_dir)
    prepare_geno_dir()

def getPopsForChr(self, chr):
    '''Check if database has a certain chromosome. Returns list of pops.
    '''
    pops = []
    c = self.dbConn.cursor()
    c.execute('''
	SELECT DISTINCT population
	  FROM freq_snp
	  WHERE chromosome = ?''', chr)
    for pop in c:
	pops.append(pop)
    return pops

def getChrsForPop(self, pop):
    '''Check if database has a certain population. Returns list of chrs.
    '''
    chrs = []
    c = self.dbConn.cursor()
    c.execute('''
	SELECT chromosome 
	  FROM freq_snp
	  WHERE population = ?''', pop)
    for chr in c:
	chrs.append(chr)
    return chrs

def has_chr_pop(chr, pop):
    '''Check if database has a certain chromosome and population.
    '''
    return os.path.exists(geno_dir + os.sep + pop + "-" + str(chr))

def infer_from_list(lst):
    for f in lst:
        f = f.split('/')[-1]
        if (f.find('17') > -1 and f.find('YRI') > -1):
            return f.replace('17', 'CHRNUM').replace('YRI','PNAME')

def infer_file_names():
    '''Infer the file names.
    Sometimes HapMap changes the file names, so we try to be robust
    '''
    ftp = FTP('www.hapmap.org')
    ftpDir = '/genotypes/latest/fwd_strand/non-redundant'
    ftp.login()
    flist = ftp.nlst(ftpDir)
    ftp.close()
    return infer_from_list(flist)



def require_chr_pop(chrom, pop):
    '''Requires a chromosome and population.

    This will retrieve data from HapMap and load it on the database.
    Hogs network and takes time and disk space!
    Load might fail, repeating might be a possibility.
    '''
    if not has_chr_pop(chrom, pop):
        tempHMFile = geno_dir + os.sep + 'tmp_hp.gz'
        try:
            os.remove(tempHMFile)
        except OSError:
            pass #Doesn't exist -> OK
        print 'Will load', chrom, pop, 'NOTE: This will take a lot of time'
        #tpl_name = infer_file_names()
        #name = tpl_name.replace('PNAME', pop).replace('CRHNUM', str(chrom))
        #ftp = FTP('www.hapmap.org')
        #ftpDir = 'frequencies/latest/fwd_strand/non-redundant'
        #ftp.login()
        #ftp.retrbinary('RETR ' + ftpDir + '/' + name,
        #    open(tempHMFile, 'wb').write)
        #ftp.close()
        #os.rename(geno_dir + os.sep + 'tmp_hp.gz', geno_dir + os.sep + pop + "-" + str(chrom))

def require_pop(pop):
    '''Requires a population, all chromosomes.

    This will retrieve data from HapMap and load it on the database.
    Hogs network and takes time and disk space!
    Load might fail, repeating might be a possibility.
    '''
    for chrom in range(1, 23):
        require_chr_pop(str(chrom), pop)
    require_chr_pop('X', pop)
    require_chr_pop('Y', pop)

def require_chr(chrom, pop_list):
    '''Requires a chromosome for specified populations.

    This will retrieve data from HapMap and load it on the database.
    We let the user decide which populations to load, mainly because
	we don't know between (JPT and CHB) or JPT+CHB
    Hogs network and takes time and disk space!
    Load might fail, repeating might be a possibility.
    '''
    for pop in pop_list:
	require_chr_pop(chrom, pop)

def cleanDB(self):
    '''Cleans the database. Caution...
    '''
    self.dbConn.execute('DELETE FROM freq_snp')

def getRSsForInterval(self, chrom, begin, end):
    '''Gets a list of rs_ids for an interval (limits included)
    '''
    c = self.dbConn.cursor()
    c.execute('''
	SELECT rs_id
	  FROM freq_snp
	WHERE chromosome = ?
	  AND position >= ?
	  AND position <= ?''', (chrom, begin, end))
    rs_list = []
    for rs in c:
	rs_list.append(rs[0])
    return rs_list

def report_dups(chrom, pop):
    '''Reports indivuduals with duplicates.'''
    dups = []
    require_chr_pop(chrom, pop)
    gf = gzip.open(geno_dir + os.sep + pop + "-" + str(chrom))
    for name in gf.readline().split(' ')[11:]:
        if name.endswith('.dup'): dups.append(name[:-4])
    gf.close()
    return dups

def start_chr_pop(chrom, pop):
    '''Starts the processing of a chrom/pop case.

       Returns an handle.
    '''
    require_chr_pop(chrom, pop)
    dups = report_dups(chrom, pop)
    gf = gzip.open(geno_dir + os.sep + pop + "-" + str(chrom))
    indivs = gf.readline().rstrip().split(' ')[11:]
    indivPos = {}
    for i in range(len(indivs)):
        indivPos[indivs[i]] = i + 11
    return gf, dups, indivPos

def get_dups(handle):
    return handle[1]

def get_indivs(handle):
    return filter(lambda x : not x.endswith('.dup'), handle[2].keys())


def next_record(handle):
    '''Returns the next record.

       The record should be processed by record processing functions.
    '''
    line = handle[0].readline()
    if line=='': return None
    return line.rstrip().split(' ')

def end_chr_pop(handle):
    handle[0].close()

def get_snp(record):
    return record[0][2:]

def get_alleles(record):
    return record[1][0]+record[1][2]

def get_indiv(handle, record, indiv, dup = False, dup_constraint = True):
    if dup_constraint and handle[2].has_key(indiv+'.dup'):
        dup = record[handle[2][indiv+'.dup']]
        orig = record[handle[2][indiv]]
        if dup=='NN': return orig
        elif orig=='NN': return dup
        elif dup<>orig:
            return None
        else:
            return orig # == dup
    if dup:
        return record[handle[2][indiv+'.dup']]
    else:
        return record[handle[2][indiv]]

def has_same_dup(handle, record, indiv):
    if record[handle[2][indiv+'.dup']] == record[handle[2][indiv]]:
      return True
    if record[handle[2][indiv+'.dup']] == (
        record[handle[2][indiv]][1] + record[handle[2][indiv]][0]):
      return True

    return False
