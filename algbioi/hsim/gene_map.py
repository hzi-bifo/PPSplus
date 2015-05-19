#!/usr/bin/env python

"""
    Copyright (C) 2015  Ivan Gregor

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


    ***********************************************************************

    Manage mapping of (joined pair-end) reads onto reference sequences based on reference genes.
"""

import os
# import sys
import re
import gzip

import numpy as np
# import multiprocessing as mp

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

import comh
from algbioi.com import fq
from algbioi.com import fasta as fas
from algbioi.com import parallel


class _RefAccVersRec():
    def __init__(self, accVersion):
        """
            Record for one accession version, i.e. one reference sequence
        """
        self._geneArray = []  # tuple array (gPos, gLen, gStrain, gName)
        self._genePosArray = None  # numpy position array
        self._accVersion = accVersion
        self._sorted = False

    def update(self, gPos, gLen, gStrain, gName):
        """
            Add a gene to the list.
        """
        assert not self._sorted
        self._geneArray.append((gPos, gLen, gStrain, gName))

    def sort(self):
        """
            Sort the genes according to the start positions.
        """
        self._geneArray.sort(key=lambda x: x[0])
        # numpy array of start positions
        self._genePosArray = np.zeros(len(self._geneArray), np.int32)
        for i in range(len(self._geneArray)):
            self._genePosArray[i] = self._geneArray[i][0]
        self._sorted = True

    def getIntersect(self, rPos, rLen):
        """
            Given coordinates of a read, get all intersections with the genes.

            Return a list of: (rPosOnRead, rPosOnGene, overlapLen, gStrain, gName)
        """
        retList = None
        assert self._sorted
        if len(self._genePosArray) == 0:
            return None

        # do binary search to find this
        i = np.searchsorted(self._genePosArray, rPos)
        i = int(i)

        if i - 1 >= 0:  # there is a gene before
            gPos, gLen, gStrain, gName = self._geneArray[i - 1]

            if rPos < gPos + gLen:  # got an overlap
                rPosOnRead = 0
                rPosOnGene = rPos - gPos
                overlapLen = min(gPos + gLen - rPos, rLen)
                retList = [(rPosOnRead, rPosOnGene, overlapLen, gStrain, gName)]

                assert 0 <= rPosOnRead < rLen and 0 <= rPosOnGene < gLen and 0 < overlapLen <= rLen

        if i < len(self._genePosArray):
            gPos, gLen, gStrain, gName = self._geneArray[i]
            if rPos + rLen > gPos:
                rPosOnRead = gPos - rPos
                rPosOnGene = 0
                if gPos + gLen < rPos + rLen:
                    overlapLen = gLen
                else:
                    overlapLen = rPos + rLen - gPos
                entry = (rPosOnRead, rPosOnGene, overlapLen, gStrain, gName)
                if retList is None:
                    retList = [entry]
                else:
                    retList.append(entry)

                assert 0 <= rPosOnRead < rLen and 0 <= rPosOnGene < gLen and 0 < overlapLen <= rLen

        if i + 1 < len(self._geneArray):
            gPos, gLen, gStrain, gName = self._geneArray[i + 1]
            if rPos < gPos + gLen and gPos < rPos + rLen:
                if gPos <= rPos:
                    rPosOnRead = 0
                    rPosOnGene = rPos - gPos
                    overlapLen = min(rPos + rLen, gPos + gLen) - rPos
                else:
                    rPosOnRead = gPos - rPos
                    rPosOnGene = 0
                    overlapLen = min(rPos + rLen, gPos + gLen) - gPos
                entry = (rPosOnRead, rPosOnGene, overlapLen, gStrain, gName)

                if retList is None:
                    retList = [entry]
                else:
                    retList.append(entry)

                assert 0 <= rPosOnRead < rLen and 0 <= rPosOnGene < gLen and 0 < overlapLen <= rLen

        return retList


class _RefAccRec():
    def __init__(self, accession, refFasFile):
        """
            Record for one accession, i.e. a reference file containing a list of accession version ids.
        """
        self._accession = accession
        self._refFasFile = refFasFile
        self._seqD = fas.fastaFileToDictWholeNames(refFasFile)
        self._accVerToAccVerRec = {}

    def update(self, headerDict, geneSeq):
        """
            Add one gene entry, update the records
        """
        accVersion = headerDict['accVersion']
        if accVersion in self._seqD:

            # get the position of the gene within the reference sequence
            refSeq = self._seqD[accVersion]
            pos = refSeq.find(geneSeq)
            strain = 1
            if pos == -1:
                revGeneSeq = str(Seq(geneSeq, generic_dna).reverse_complement())
                pos = refSeq.find(revGeneSeq)
                strain = -1
            assert pos != -1

            # add the accVersion record if it doesn't exist
            if accVersion not in self._accVerToAccVerRec:
                self._accVerToAccVerRec[accVersion] = _RefAccVersRec(accVersion)

            self._accVerToAccVerRec[accVersion].update(pos, len(geneSeq), strain, headerDict['geneName'])

        else:
            print 'skipping record %s not found in %s' % (accVersion, self._refFasFile)

    def sort(self):
        for accVerRec in self._accVerToAccVerRec.values():
            accVerRec.sort()

    def getAccVtoR(self):
        return self._accVerToAccVerRec


def mapReadsToGenes(accVtoR, inSamFile, outSamFile):
    """


    """
    # get SAM read generator
    if inSamFile.endswith('.gz'):
        samGen = gzip.open(inSamFile, mode='r')
    else:
        samGen = open(inSamFile)

    # get a SAM writer
    samOut = fq.WriteFq(outSamFile)
    mappedReadCount = 0

    # for each SAM read entry
    for line in samGen:
        if line.startswith('@'):
            continue
        line = line.strip()
        tokens = line.split('\t')

        assert len(tokens) >= 9, str(len(tokens))

        # read SAM entries
        # flag = int(tokens[1])
        accVersion = tokens[2]
        readPos = int(tokens[3]) - 1  # 0 based now!
        readLen = abs(int(tokens[8]))
        # strain?
        # if flag & 0x10 != 0:
        #     revCompl = True
        # else:
        #     revCompl = False

        # sequence with this accession version contain genes
        record = ''
        if accVersion in accVtoR:
            entryList = accVtoR[accVersion].getIntersect(readPos, readLen)

            if entryList is not None:
                mappedReadCount += 1
                e = 0
                for rPosOnRead, rPosOnGene, overlapLen, gStrain, gName in entryList:
                    e += 1
                    record += '%s,%s,%s,%s,%s' % (rPosOnRead, rPosOnGene, overlapLen, gStrain, gName)
                    if e < len(entryList):
                        record += '|'

        samOut.write('%s\t%s\n' % (line, record))

    samOut.close()
    return mappedReadCount


def getReadsToGenesMap(refList, samList, genesDir):
    """
        Map reads onto genes, update the SAM file annotations for the reads.
    """
    # get map: accession -> reference file
    accToRefPath = {}
    for refFna in refList:
        acc = re.sub('(.*).fna.gz', r'\1', os.path.basename(refFna).split('_', 1)[1])
        # for one accession, get one reference sequence (get the longest if there are more)
        if acc not in accToRefPath or os.path.getsize(accToRefPath[acc]) < os.path.getsize(refFna):
            accToRefPath[acc] = refFna

    # get map: accession -> reference record
    accToRefRecord = {}
    for acc, refPath in accToRefPath.iteritems():
        accToRefRecord[acc] = _RefAccRec(acc, refPath)

    # map individual genes onto the reference sequences
    for f in os.listdir(genesDir):
        # for all gene sequences
        for header, seq in fas.fastaFileToDictWholeNames(os.path.join(genesDir, f)).iteritems():
            # get gene sequence header
            headerD = dict(zip(map(lambda x: x.split(':')[0], header.split(';')),
                               map(lambda x: x.split(':')[1], header.split(';'))))
            # get accession
            acc = headerD['accession']
            # gene sequence belong to one of the reference sequences being considered
            if acc in accToRefRecord:
                accToRefRecord[acc].update(headerD, seq)

    # sort all acc records
    for refRecord in accToRefRecord.values():
        refRecord.sort()

    # get some statistics
    # print 'ref sequences:', len(accToRefRecord)
    # c = 0
    # g = 0
    # for refRecord in accToRefRecord.values():
    #     c += len(refRecord._accVerToAccVerRec)
    #     for e in refRecord._accVerToAccVerRec.values():
    #         g += len(e._geneArray)
    # print 'acc versions: ', c
    # print 'genes', g

    # get mapping: accVersion -> _RefAccVersRec
    # man = mp.Manager()
    # accVtoR = man.dict()  # too slow with the manager
    accVtoR = {}
    for refRecord in accToRefRecord.values():  # for all ref records
        refRecordD = refRecord.getAccVtoR()
        sizeI = len(refRecordD)
        sizeA = len(accVtoR)
        accVtoR.update(refRecordD)
        assert len(accVtoR) == sizeI + sizeA, ' %s %s != %s' % (len(accVtoR), sizeI, sizeA)

    # define tasks: map reads onto genes
    taskList = []
    mappedReadCount = 0
    for inSamFile, outSamFile in samList:
        # mappedReadCount += mapReadsToGenes(accVtoR, inSamFile, outSamFile)  # to run in serial
        taskList.append(parallel.TaskThread(mapReadsToGenes, (accVtoR, inSamFile, outSamFile)))
    retList = parallel.runThreadParallel(taskList, comh.MAX_PROC)
    for e in retList:
        mappedReadCount += e

    print 'Mapped read count: %s' % mappedReadCount


if __name__ == "__main__":
    pass