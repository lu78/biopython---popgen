'''
Accessing UCSC Table Browser.
'''
__docformat__ = 'javadoc en'
from ftplib import FTP
import sqlite3
from sys import exit
import gzip
import os

from PopGen import Config

ftp_site = 'hgdownload.cse.ucsc.edu'
ftp_root = 'goldenPath'

def prepareDataDir(db):
    tableDir = Config.dataDir + os.sep +  'TableBrowser'
    try:
        os.mkdir(tableDir)
    except OSError:
        pass #Probably OK, dir already exists
    dbConn = sqlite3.connect(tableDir + os.sep + db + '.db')
    try:
        dbConn.execute('''CREATE TABLE known_gene (
                            asc_id     VARCHAR(20),
                            chromosome VARCHAR(2),
                            strand     CHAR(1),
                            tx_start   INTEGER,
                            tx_end     INTEGER,
                            cds_start  INTEGER,
                            cds_end    INTEGER,
                            prot_id    VARCHAR(100)
                            )''')
        dbConn.execute('CREATE INDEX kg_asc_id ON known_gene(asc_id)')
        dbConn.execute('CREATE INDEX kg_prot_id ON known_gene(prot_id)')
        dbConn.execute('''CREATE TABLE gene_exons (
                            asc_id     VARCHAR(20),
                            start      INTEGER,
                            finish     INTEGER
                            )''')
        dbConn.execute('CREATE INDEX ge_asc_id ON gene_exons(asc_id)')
    except sqlite3.OperationalError:
        pass
    return dbConn


def loadFile(db, dir, file, unzip = False):
    tableDir = Config.dataDir + os.sep +  'TableBrowser'
    tempTBFile = tableDir + os.sep + 'tmp_tb'
    tempTBFileGz = tempTBFile + '.gz'

    ftp = FTP(ftp_site)
    ftpDir = ftp_root + '/' + db + '/' + dir
    fname = ftpDir + '/' + file
    ftp.login()
    if unzip:
        ftp.retrbinary('RETR ' + fname, open(tempTBFileGz, 'wb').write)
        gz = gzip.open(tempTBFileGz, 'rb')
        uncomp = open(tempTBFile, 'w')
        l = gz.readline()
        while l<>'':
            l = gz.readline()
            uncomp.write(l)
        gz.close()
        uncomp.close()
    else:
        ftp.retrbinary('RETR ' + fname, open(tempTBFile, 'wb').write)
    ftp.close()
    return tempTBFile

def countLines(dbConn, table):
    c = dbConn.cursor()
    c.execute('SELECT count(*) FROM ' + table)
    cnt = iter(c).next()
    return cnt[0]


class KnownGene:
    '''Wrapper class for KnownGene table
    '''
    def __init__(self, db):
        '''Inits the object. NOTE: will load the database, if not loaded!

        @param db UCSC table, like hg18 or bosTau1
        '''
        self.dbConn = prepareDataDir(db)
        self.db     = db

    def loadDB(self):
        '''Loads the database (also does cleanup).
        '''
        self.cleanDB()
        fName = loadFile(self.db, 'database', 'knownGene.txt.gz' , True)
        f = open(fName, 'r')
        l = f.readline()
        c = self.dbConn.cursor()
        c.execute('DELETE FROM known_gene')
        c.execute('DELETE FROM gene_exons')
        while l<>'':
            toks       = l.rstrip().split('\t')
            ascId      = toks[0]
            chr        = toks[1][3:]
            strand     = toks[2]
            txStart    = int(toks[3])
            txEnd      = int(toks[4])
            cdsStart   = int(toks[5])
            cdsEnd     = int(toks[6])
            numExons   = int(toks[7])
            exonsStart = toks[8].split(',')
            exonsEnd   = toks[9].split(',')
            protId     = toks[10]
            c.execute('''
                INSERT INTO known_gene (
                    asc_id, chromosome, strand, tx_start, tx_end,
                    cds_start, cds_end, prot_id)
                VALUES (?,?,?,?,?,?,?,?)''',
                (ascId, chr, strand, txStart, txEnd,
                  cdsStart, cdsEnd, protId))
            for i in range(len(exonsStart)):
                exonStart = exonsStart[i]
                if exonStart == '': continue
                exonEnd = exonsEnd[i]
                print 'es', exonStart
                c.execute('''
                    INSERT INTO gene_exons (asc_id, start, finish)
                         VALUES (?, ?, ?)
                ''', (ascId, exonStart, exonEnd))
            l = f.readline()
        f.close()
        self.dbConn.commit()

    def cleanDB(self):
        """Cleans the database.
 
        Deletes all data from gene_exons and known_gene.
        """
        self.dbConn.execute('DELETE FROM gene_exons')
        self.dbConn.execute('DELETE FROM known_gene')
        
    def close(self):
        '''Closes the database.
        '''
        print 'close'
        self.dbConn.close()

    def getAscIdsFromProtId(self, protId):
        '''Returns a list of ascIds for a certain protId.

        @param protId Protein Id, like P53_HUMAN.

        @return List of ascIds.
        '''
        c = self.dbConn.cursor()
        c.execute('''
            SELECT asc_id
              FROM known_gene
             WHERE prot_id = ?''', (protId,))
        asc_list = []
        for asc in c:
            asc_list.append(asc[0])
        return asc_list

    def getAscId(self, ascId):
        '''Gets all info (except exons) for a certain ascId.

        @param ascId ascId.
        @return a tuple with all data.
        '''
        c = self.dbConn.cursor()
        c.execute('''
            SELECT *
              FROM known_gene
             WHERE asc_id = ?''', (ascId,))
        for asc in c:
            return asc


