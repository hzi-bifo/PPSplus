#!/usr/bin/env python

#/AM/home-0/shared/python/Python-2.7.1/python

import sqlite3
import re
from Config import Config
from sets import Set
import os
import glob
import re


class DBData():

    def __init__(self, ncbiProcessDir, databaseFile):
        self._databaseFile = databaseFile
        #buffer genome/WGS info
        count = 0
        self._ncbiBuffer = Set([]) #set of ncbids for which there is a genome/WGS in the database
        for filePath in glob.glob(os.path.join(os.path.normpath(ncbiProcessDir),'*.[0-9]*.f[an][sa]')):
            count += 1
            ncbid = int(re.sub(os.path.join(r'.*[^0-9]([0-9]+).[0-9]+.f[an][sa]$'),r'\1' ,filePath))
            self._ncbiBuffer.add(ncbid)
            #print ncbid
        print count, 'genomes/WGS found in directory ', os.path.normpath(ncbiProcessDir)

    #get the number of genomes/wgs available in the directory
    #
    #param threshold: it returns max. threshold genomes/wgs (after this threshold is reached, it returns)
    #param dir: directory that contains genomes/wgs in the form: "ncbid.[0-9]*.f[an][sa]"
    #
    #return: the number of genomes/wgs from different species that are subclades of the input ncbid
    def getGenomeWgsCount(self, ncbid, threshold):

        try:
            conn = sqlite3.connect(os.path.normpath(self._databaseFile))
            cursor = conn.cursor()

            speciesIdsList = []
            self._collectSpecies(speciesIdsList, cursor, ncbid, threshold)
            return len(speciesIdsList)
        except Exception:
            print "Failed to create connection to a database:", self._databaseFile
            raise
        finally:
            cursor.close()
            conn.close()


    #get the first ncbid of species/subspecies from the input list that is contained in the directory
    #
    #param dir: directory that contain genomes/wgs files
    #
    #return None or the first ncbid for which there is a genome or a draft genome in the directory
    def _genomeExists(self, listOfNcbids):

        #for ncbid in listOfNcbids:
        #    if glob.glob(os.path.join(os.path.normpath(dir), str(str(ncbid) + '.[0-9]*.f[an][sa]'))):
        #        return ncbid
        for ncbid in listOfNcbids:
            if ncbid in self._ncbiBuffer:
                return ncbid
        return None


    #param speciesIds: output list of ncbids (species or subspecies) for which there are genomes/wgs in the directory
    def _collectSpecies(self, speciesIds, cursor, root, threshold):
        if len(speciesIds) >= threshold:
            return
        cursor.execute('SELECT node_rank FROM taxon T WHERE T.ncbi_taxon_id=?',(root,))
        result = cursor.fetchall()
        if len(result) == 0:
            return
        assert len(result) == 1
        if result[0][0] == 'species':
            list = []
            list.append(root)
            self._collectSubSpecies(list, cursor, root)
            ncbid = self._genomeExists(list) #now check if at least one genome exists from the list
            if ncbid != None:
                speciesIds.append(ncbid)
        else:
            cursor.execute('SELECT ncbi_taxon_id FROM taxon T WHERE T.parent_taxon_id=?',(root,))
            result = cursor.fetchall()
            for item in result:
                self._collectSpecies(speciesIds, cursor, int(item[0]), threshold)


    #get ncbids of all subspecies of the root
    #
    #list: output list with all ncbids of the subspecies
    def _collectSubSpecies(self, list, cursor, root):
        cursor.execute('SELECT ncbi_taxon_id FROM taxon T WHERE T.parent_taxon_id=?',(root,))
        result = cursor.fetchall()
        if len(result) > 0:
            for item in result:
                id = int(item[0])
                list.append(id)
                self._collectSubSpecies(list, cursor, id)


def test(ncbid):
    config = Config(open(os.path.normpath('D://A_Phylo//A_Metagenomic//pPPS//workspace//pPPS//config01.cfg')), 'pPPS')
    databaseFile = os.path.normpath(config.get('databaseFile'))
    ncbiProcessDir = os.path.normpath('D://A_Phylo//A_Metagenomic//pPPS//workspace//pPPS//wdir02//ncbiProcDir')
    dbData = DBData(ncbiProcessDir, databaseFile)
    threshold = 3
    print dbData.getGenomeWgsCount(ncbid, threshold)

    #config = Config(open(os.path.normpath('//AM//metagenomic//work//projects//pPPS//tests//TW//TW01//config.cfg')), 'pPPS')

    #threshold = 3
    #dir = 'D://A_Phylo//A_Metagenomic//pPPS//workspace//pPPS//genomes'
    #dir = '//AM//metagenomic//work//projects//pPPS//tests//TW//TW01//ncbiProcDir'

    #databaseFile = os.path.normpath(config.get('databaseFile'))
    #taxonomicRanks = config.get('taxonomicRanks').split(',')

    #count = getGenomeWgsCount(ncbid, threshold, dir, databaseFile, taxonomicRanks)

    #print count, 'genomes/wgs for ncbid:', ncbid


if __name__ == "__main__":
  #test(122)
  #haveData(126)
  #haveData(84999) #Coriobacteriales
  #haveData(171549) #Bacteroidales
  #haveData(815) #Bacteriodaceae
  #haveData(171551) #Porphyromonadaceae
  #haveData(171552) #Prevotellaceae
  #haveData(171550) #Rikenellaceae
  test(976) #Bacteroidetes
  #haveData(200666) #Sphingobacteriales
  #haveData(768503) #Cytophagia
  #haveData(117743) #Flavobacteria
  #haveData(475963) #Caldilineales
  #haveData(292625) #Anaerolineae
  #haveData(200795) #Chloroflexi
  #haveData(204431) #Fibrobacteraceae (59374, 834)
  #haveData(186803) #Lachnospiraceae
  #haveData(541000) #Ruminococcaceae
  #haveData(186802) #Clostridiales
  #haveData(31979) #Clostridiaceae
  #haveData(186806) #Eubacteriaceae
  #haveData(186807) #Peptococcaceae
  #haveData(31977) #Veillonellaceae
  #haveData(186801) #Clostridia
  #haveData(128827) #Erysipelotrichaceae
  #haveData(1239) #Firmicutes
  #haveData(91061) #Bacilli
  #haveData(255528) #Victivallaceae (340101)
  #haveData(126) #Planctomycetaceae
  #haveData(481) #Neisseriaceae
  #haveData(213121) #Desulfobulbaceae (577650, 177439, 589865)
  #haveData(213421) #Desulfuromonaceae
  #haveData(69541) #Desulfuromonadales
  #haveData(72294) #Campylobacteraceae
  #haveData(1224) #Proteobacteria
  #haveData(28211) #Alphaproteobacteria
  #haveData(1236) #Gammaproteobacteria
  #haveData(137) #Spirochaetaceae
  #haveData(186333) #Anaeroplasmataceae
  #haveData(186332) #Anaeroplasmatales
  #haveData(31969) #Mollicutes

  #haveData(278082) #Victivallales
  #haveData(256845) #Lentisphaerae

