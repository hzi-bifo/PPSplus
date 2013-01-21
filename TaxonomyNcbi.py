#!/usr/bin/env python

#/AM/home-0/shared/python/Python-2.7.1/python

import os
import sys
import re
import sqlite3
from sets import Set


#Represents an interface to the sqlite3 database in which the NCBI database is stored.
#(NOTE that methods and variables starting with "_" are local and shouldn`t be used from the outside)
class TaxonomyNcbi():

    #Constructor
    #@param databaseFile: usually file named "ncbitax_sqlite.db"
    #@param allowedRanks: taxonomic ranks that will be considered (where 'root' is the root of the taxonomy)
    #@param considerNoRank: consider ranks 'no rank' if true
    def __init__(self, databaseFile, allowedRanks=['root','superkingdom','phylum','class','order','family','genus','species'],
                 considerNoRank=False):

        self._allowedRanks = Set(allowedRanks)
        if considerNoRank:
            self._allowedRanks.add('no rank')
        try:
            self.conn = sqlite3.connect(os.path.normpath(databaseFile))
            self.cursor = self.conn.cursor()
        except Exception:
            sys.stderr.write(str('TaxonomyNcbi: Failed to create connection to database: ' + databaseFile))
            raise


    #@return scientific name or None
    def getScientificName(self, ncbid, checkRank = False):

        if ncbid == -1:
            ncbid = 1
            sys.stderr.write('ncbid -1 converted to 1')

        if checkRank and (not self.isRankNcbidAllowed(ncbid)):
            return None

        self.cursor.execute(str('SELECT TN.name FROM taxon_name TN, taxon T WHERE T.ncbi_taxon_id=?' +
                                ' AND T.taxon_id = TN.taxon_id AND TN.name_class="scientific name"'),(ncbid,))
        result = self.cursor.fetchall()
        if len(result) == 1:
            return result[0][0]
        else:
            sys.stderr.write(str('TaxonomyNcbi: Cannot find name for ncbi: ' + str(ncbid)))
            return None


    #@return ncbid or None
    def getNcbid(self, scientificName, checkRank = False):
        self.cursor.execute(str('SELECT T.ncbi_taxon_id FROM taxon_name TN, taxon T ' +
                                'WHERE TN.name_class="scientific name" AND TN.name=? AND TN.taxon_id=T.taxon_id'),
                                (scientificName,))
        result = self.cursor.fetchall()
        if len(result) == 1:
            ncbid = int(result[0][0])
            if checkRank and (not self.isRankNcbidAllowed(ncbid)):
                return None
            return ncbid
        else:
            if len(result) == 0:
                sys.stderr.write(str('TaxonomyNcbi: Cannot find scientific name "' + scientificName + '" in the database.\n'))
            else:
                sys.stderr.write(str('TaxonomyNcbi: scientific name "' + scientificName + '" is ambiguous!\n'))
            return None

    #@return ncbid or None; name doesn`t have to be a scientific name
    def getNcbid2(self, name, checkRank = False):
        self.cursor.execute(str('SELECT T.ncbi_taxon_id FROM taxon_name TN, taxon T ' +
                                'WHERE TN.name=? AND TN.taxon_id=T.taxon_id'),
                                (name,))
        result = self.cursor.fetchall()
        if len(result) == 1:
            ncbid = int(result[0][0])
            if checkRank and (not self.isRankNcbidAllowed(ncbid)):
                return None
            return ncbid
        else:
            if len(result) == 0:
                sys.stderr.write(str('TaxonomyNcbi: Cannot find scientific name "' + name + '" in the database.\n'))
            else:
                sys.stderr.write(str('TaxonomyNcbi: scientific name "' + name + '" is ambiguous!\n'))
            return None


    def childrenNcbids(self, ncbid): #SELECT T1.ncbi_taxon_id from taxon T1 where T1.parent_taxon_id=818;
        self.cursor.execute(str('SELECT T1.ncbi_taxon_id from taxon T1 where T1.parent_taxon_id=?'),(ncbid,))
        result = self.cursor.fetchall()
        if len(result) == 0:
            return None
        else:
            resultList = []
            for i in result:
                resultList.append(i[0])
            return resultList


    #@return ncbid or None
    def getParentNcbid(self, ncbid):
        if ncbid == 1:
            return None
        taxonId = self._getTaxonId(ncbid)

        while True:
            if ncbid == 1 or ncbid == None or taxonId == None: #the root of the taxonomy reached
                return None
            ncbid = self._getParentNcbid(taxonId)

            taxonId = self._getTaxonId(ncbid)
            rank = self._getRank(taxonId)
            if (rank in self._allowedRanks) or (ncbid == 1 and 'root' in self._allowedRanks):
                return ncbid


    #@return rank or None
    def getRank(self, ncbid, checkRank = False):
        if checkRank and (not self.isRankNcbidAllowed(ncbid)):
            return None
        return self._getRank(self._getTaxonId(ncbid))


    #@return: True or False
    def isRankNcbidAllowed(self, ncbid):
        taxonId = self._getTaxonId(ncbid)
        rank = self._getRank(taxonId)
        if rank in self._allowedRanks:
            return True
        else:
            return False


    #@return True or False
    def isRankAllowed(self, rank):
        if rank in self._allowedRanks:
            return True
        else:
            return False


    #close the database after you stop using it
    def close(self):
        self.cursor.close()
        self.conn.close()


    def _getTaxonId(self, ncbid):
        if ncbid == None:
            return None
        if ncbid == -1:
            ncbid = 1
            sys.stderr.write('ncbid=(-1) converted to ncbid=(1)\n')
        self.cursor.execute('SELECT taxon_id FROM taxon T WHERE T.ncbi_taxon_id=?',(ncbid,))
        result = self.cursor.fetchall()
        if len(result) != 1:
            sys.stderr.write('TaxonomyNcbi: Cannot find taxon_id for ncbi:' + str(ncbid) + ' result:' + str(result) + ' \n')
            return None
        return int(result[0][0])


    def _getParentNcbid(self, taxonId):
        if taxonId == None:
            return None
        self.cursor.execute('SELECT parent_taxon_id FROM taxon T WHERE T.taxon_id=?', (taxonId,))
        result = self.cursor.fetchall()
        if len(result) != 1:
            sys.stderr.write(str('TaxonomyNcbi: Cannot find parent for taxon_id' + str(taxonId)))
            return None
        return int(result[0][0])


    def _getRank(self, taxonId):
        if taxonId == None:
            return None
        self.cursor.execute('SELECT node_rank FROM taxon T WHERE T.taxon_id=?', (taxonId,))
        result = self.cursor.fetchall()
        if len(result) != 1:
            sys.stderr.write(str('TaxonomyNcbi: Cannot find rank for taxon_id: ' + str(taxonId)))
            return None
        return str(result[0][0])


def test():
    databaseFile = "/Users/ivan/Documents/work/binning/taxonomy/ncbi_taxonomy_20110629/ncbitax_sqlite.db"
    taxonomy = TaxonomyNcbi(databaseFile)

    print 'Scientific name for ncbid 286730 (Alkaliflexus imshenetskii) is:', taxonomy.getScientificName(286730)

    print 'Ncbid of Lachnospiraceae (186303) is:', str(taxonomy.getNcbid('Lachnospiraceae'))

    print 'Get Parent of 167965 which is 74152:', str(taxonomy.getParentNcbid(167965))

    print 'Get rank of Lachnosipraceae (family):', str(taxonomy.getRank(186803))


    taxonomy.close()

def test2():
    databaseFile = "/Users/ivan/Documents/work/binning/taxonomy/ncbi_taxonomy_20110629/ncbitax_sqlite.db"
    taxonomy = TaxonomyNcbi(databaseFile)


    ncbid = 1
    parent = taxonomy.getParentNcbid(ncbid)
    print ncbid, taxonomy.getScientificName(ncbid)
    print parent, taxonomy.getScientificName(parent)
    print parent, taxonomy.getRank(parent)




    taxonomy.close()


if __name__ == "__main__":
  test2()