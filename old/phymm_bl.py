#!/usr/bin/env python

import os

from com.csv import getMapping
from com.csv import OutFileBuffer
from com.fasta import fastaFileToDict
from com.taxonomy_ncbid import TaxonomyNcbi


class DummyOutBuffer():

    def writeText(self, text):
        pass
    def close(self):
        pass

class PhymmBLConfig():

    def __init__(self, fastaFilePath, seqIdToNcbid, databaseFile, ranksToConsider, outDir, outConfigFile, pathPrefix):
        """
            @param: referenceFastaFile: seq_id -> sequence
            @param: map: seq_id -> ncbid
            @param: database
            @param: taxonomic ranks: phylum, class, order, family, genus, species, (not defined at a rank: NO_VALUE)
            @param: output directory
            @param: outConfigFile: customGenomicData configuration file
            @param: relative path prefix
        """
        self._seqNameToSeq = fastaFileToDict(fastaFilePath)
        self._seqIdToNcbid = getMapping(seqIdToNcbid, 0, 1, sep=None, comment = '#')
        self._taxonomy = TaxonomyNcbi(databaseFile)
        self._ranksToConsider = ranksToConsider
        self._outDir = outDir
        self._pathPrefix = pathPrefix
        self._outConfigFile = outConfigFile


    def generateCustomData(self, genConfigFile=True, genCustomData=True):

        if genConfigFile:
            configOutFile = OutFileBuffer(self._outConfigFile)
        else:
            configOutFile = DummyOutBuffer()

        ncbidToSeqNameList = dict([])
        for seqName in self._seqNameToSeq:

            assert seqName in self._seqIdToNcbid
            assert len(self._seqIdToNcbid[seqName]) == 1
            ncbid = int(self._seqIdToNcbid[seqName][0])
            if ncbid not in ncbidToSeqNameList:
                ncbidToSeqNameList[ncbid] = [seqName]
            else:
                ncbidToSeqNameList[ncbid].append(seqName)

        for ncbid in ncbidToSeqNameList:
            seqNameList = ncbidToSeqNameList[ncbid]

            fileName = str(str(ncbid) + '.fna')
            if genCustomData:
                outFastaFile = OutFileBuffer(os.path.join(self._outDir, fileName))
            else:
                outFastaFile = DummyOutBuffer()

            for seqName in seqNameList:
                seq = self._seqNameToSeq[seqName]
                outFastaFile.writeText(str('>' + str(ncbid) + '_' + str(seqName) + '\n'))
                outFastaFile.writeText(str(str(seq) + '\n'))
            outFastaFile.close()

            configOutFile.writeText(str(os.path.normpath(os.path.join(self._pathPrefix, fileName)) + '\t1\t'))

            rankToNcbid = dict([])
            rankToNcbid['species'] = ncbid
            current = ncbid
            while current != 1 and current != None:
                current = self._taxonomy.getParentNcbid(current)
                rank = self._taxonomy.getRank(current, checkRank=True)
                rankToNcbid[rank] = current

            for rank in self._ranksToConsider:
                if rank in rankToNcbid:
                    configOutFile.writeText(str(self._taxonomy.getScientificName(rankToNcbid[rank]) + '\t'))
                else:
                    configOutFile.writeText('NO_VALUE\t')

            configOutFile.writeText('NO_VALUE\n') # strain label

        configOutFile.close()


    def close(self):
        self._taxonomy.close()




def main():
    fastaFilePath = '/local/data/work/johdro/an_paper_new/results_simmc-oldrefdata/genus-mask/reference.fna'
    seqIdToNcbid = '/local/data/work/johdro/an_paper_new/results_simmc-oldrefdata/genus-mask/gi_taxid_dna_species-ncbi_taxonomy_20110629.dmp'
    db = '/AM/metagenomic/work/users/johdro/swap/ncbi_taxonomy_20110629/ncbitax_sqlite.db'
    ranksToConsider = ['phylum','class','order','family','genus','species']
    outDir = '/AM/metagenomic/work/projects/pPPS/data/PhymmBL/simMC/genusExcluded/customData'
    outConfigFile = '/AM/metagenomic/work/projects/pPPS/data/PhymmBL/simMC/genusExcluded/configCustomData.txt'
    pathPrefix = '../simMCGenusExcluded'
    genConfigFile=True
    genCustomData=True
    phymm = PhymmBLConfig(fastaFilePath, seqIdToNcbid, db, ranksToConsider, outDir, outConfigFile, pathPrefix)
    phymm.generateCustomData(genConfigFile, genCustomData)
    phymm.close()

#120604516


if __name__ == "__main__":
    main()