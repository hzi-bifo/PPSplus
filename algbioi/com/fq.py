#!/usr/bin/env python

"""
    Copyright (C) 2014  Ivan Gregor

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

    Contains basic functionality to work with FASTAQ files.
"""

import gzip
import parallel
import multiprocessing as mp
import numpy as np
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


class ReadFqGen():
    def __init__(self, fqFilePath):
        """
            A generator for reading from a FASTAQ file. A file is considered to be compressed if it ends with ".gz".
            The generator generates tuples: (@name, DNA string, + ,QS string)

            @param fqFilePath: FASTAQ file (.fq or .fq.gz)
        """
        if fqFilePath.strip().endswith('.gz'):
            # reading from a compressed file
            self._fqOpen = gzip.open(fqFilePath)
        else:
            self._fqOpen = open(fqFilePath)

    def _readLine(self):
        """
            @return: a next stripped line
        """
        e = self._fqOpen.readline()
        if not e:
            raise StopIteration()
        return e.strip()

    def next(self):
        """
            @raise: StopIteration
            @return: tuple (@name, DNA string, +, QS string)
        """
        try:
            return self._readLine(), self._readLine(), self._readLine(), self._readLine()
        except StopIteration as e:
            self._fqOpen.close()
            raise e

    def __iter__(self):
        return self


class WriteFq():
    def __init__(self, fqFileOut, compressLevel=1, blockSize=10000):
        """
            To write to a FASTQ file. Call the close method at the end !!!

            @param fqFileOut: a file path, a file is gzipped if it ends with .gz
            @param compressLevel: int between 1 ~ fastest and 9 ~ the best compression
            @param blockSize: writing large blocks to a gzipped file is much faster than writing entry by entry
        """
        if fqFileOut.endswith('.gz'):
            self._outFq = gzip.open(fqFileOut, 'wb', compressLevel)
        else:
            self._outFq = open(fqFileOut, 'w')
        self._blocks = []
        self._blockSize = blockSize

    def _writeBlocks(self):
        self._outFq.write(''.join(self._blocks))
        self._blocks = []

    def writeFqEntry(self, name, dna, qs):
        self.write(name + '\n' + dna + '\n+\n' + qs + '\n')

    def write(self, text):
        self._blocks.append(text)
        if len(self._blocks) == self._blockSize:
            self._writeBlocks()

    def close(self):
        self._writeBlocks()
        self._outFq.close()


def joinPairEnd(fileTupleList, minOverlap=0.05, minOverlapIdentity=0.9, maxCpu=mp.cpu_count()):
    """
        Join pair end reads

        Insert must be max: 2 * readLen - int(round(readLen * minOverlap))

        @param fileTupleList: list of (fq1, fq2, fqJoin, readLen, insert, sd, qsMax); but sd not used!?
        @param minOverlap: percentage of the readLen, minimum that must overlap
        @param minOverlapIdentity: the overlapping region must have at least this identity
    """
    # define a task list
    taskList = []
    for fq1, fq2, fqJoin, readLen, insert, sd, qsMax in fileTupleList:
        taskList.append(parallel.TaskThread(_joinPairEnd, (fq1, fq2, fqJoin, readLen, insert,
                                                           qsMax, minOverlap, minOverlapIdentity)))
    # run tasks in parallel
    retList = parallel.runThreadParallel(taskList, maxCpu)

    # collect statistics
    allPairs = 0
    allSkipped = 0
    for stored, skipped in retList:
        allPairs += stored + skipped
        allSkipped += skipped
    return round(float(allSkipped) / float(allPairs) * 100, 3)


def _joinPairEnd(fq1, fq2, fqJoin, readLen, insert, qsMax, minOverlap, minOverlapIdentity, trace=False):
    """
        Join pair end reads.

        @param fq1: in FASTAQ file 1
        @param fq2: in FASTAQ file 2
        @param fqJoin: out FASTAQ file with joined reads
        @param readLen:
        @param insert:
        @param qsMax:
        @param minOverlap:
        @param minOverlapIdentity:
        @param trace:
    """
    # defines Qs for overlapping positions that match
    dnaConsensus = QsMultMatrix(qsMax)

    # min and max insert sizes to be explored
    minInsert = readLen
    maxInsert = 2 * readLen - int(round(readLen * minOverlap))
    assert insert <= maxInsert

    # possible insert sizes to be tested for the reads overlap
    possibleInsertList = _getPossibleInsertSizes(minInsert, maxInsert, insert)

    # define output file
    out = WriteFq(fqJoin)
    stored = 0
    skipped = 0

    # for each read
    for r1, r2 in zip(ReadFqGen(fq1), ReadFqGen(fq2)):
        n1, dna1, p, qs1 = r1
        n2, dna2, p, qs2 = r2

        # get reverse complement of dna2
        dna2Rev = str(Seq(dna2, generic_dna).reverse_complement())
        # reverse qs2
        qs2Rev = qs2[::-1]

        # get overlap insert size
        overlapInsert = _getOverlap(dna1, dna2Rev, readLen, possibleInsertList, minOverlapIdentity=minOverlapIdentity)

        # no reliable overlap found
        if overlapInsert is None:
            skipped += 1
            continue

        # get consensus dna and Qs of th overlapping region
        dnaCons, qsCons = dnaConsensus.getConsensus(dna1[(overlapInsert - readLen):readLen],
                                                    dna2Rev[:(2*readLen - overlapInsert)],
                                                    qs1[(overlapInsert - readLen):readLen],
                                                    qs2Rev[:(2*readLen - overlapInsert)], qsMax)

        # resulting dna and qs
        dna = str(dna1[:(overlapInsert - readLen)] + dnaCons + dna2Rev[(2*readLen - overlapInsert):])  # result dna
        qs = str(qs1[:(overlapInsert - readLen)] + qsCons + qs2Rev[(2*readLen - overlapInsert):])  # result Qs

        # store the joined read to a file
        out.writeFqEntry(n1[:-2], dna, qs)
        stored += 1

        if trace:
            print overlapInsert  # overlap
            print dna1  # dna1
            print dna1[:(overlapInsert - readLen)]  # dna1 prefix
            print str(' ' * (overlapInsert - readLen) + dna1[(overlapInsert - readLen):readLen])  # dna1 overlap
            print str(' ' * (overlapInsert - readLen) + dna2Rev[:(2*readLen - overlapInsert)])  # dna2 overlap
            print str(' ' * readLen + dna2Rev[(2*readLen - overlapInsert):])  # dna2 suffix
            print str(' ' * (overlapInsert - readLen) + dna2Rev)  # dna2
            print dna  # consensus dna
            print qs
            print qs1  # qs1
            print str(' ' * (overlapInsert - readLen) + qs2Rev)  # qs2 (rev)
            print map(lambda x: f(ord(x) - 33), list(qs))  # qs as list
            print map(lambda x: f(ord(x) - 33), list(qs1))  # qs1 as list
            print str('[' + str("'__', " * (overlapInsert - readLen)).rstrip() + str(map(lambda x: f(ord(x) - 33),
                                                                                         list(qs2Rev))))  # qs2 as list

    # close the file !!!
    out.close()
    return stored, skipped


def f(i):
    if i <= 9:
        return '0%s' % i
    else:
        return str(i)


class QsMultMatrix():
    def __init__(self, qsMax=62):
        # probability multiplication matrix Qs x Qs -> char(Qs + 33)
        self._qsMulMatrix = np.zeros((qsMax, qsMax), dtype=np.uint8)

        # Qs to probability
        qsToProb = np.zeros(qsMax, dtype=np.float64)
        for i in range(qsMax):
            qsToProb[i] = np.power(10, ((-1.) * i) / 10.)

        for i in range(qsMax):
            for j in range(qsMax):
                self._qsMulMatrix[i][j] = min(qsMax - 1, self._probToQs(qsToProb[i] * qsToProb[j])) + 33

    def _probToQs(self, prob):
        return int((-10.) * np.log10(prob))

    def _qsMul(self, qs1, qs2, qsMax):
        return min(qsMax + 32, self._qsMulMatrix[ord(qs1) - 33][ord(qs2) - 33])

    def getConsensus(self, dna1, dna2, qs1, qs2, qsMax):
        assert len(dna1) == len(dna2) == len(qs1) == len(qs2)

        consDna = np.zeros(len(dna1), dtype=np.uint8)

        consQs = np.zeros(len(dna1), dtype=np.uint8)
        for i in range(len(dna1)):
            if dna1[i] == dna2[i]:
                consDna[i] = ord(dna1[i])
                consQs[i] = self._qsMul(qs1[i], qs2[i], qsMax)
            elif ord(qs1[i]) >= ord(qs2[i]):
                consDna[i] = ord(dna1[i])
                consQs[i] = ord(qs1[i])
                # print 'mismatch ', i
            else:
                # print 'mismatch_', i
                consDna[i] = ord(dna2[i])
                consQs[i] = ord(qs2[i])

        return ''.join(map(lambda x: chr(x), consDna)), ''.join(map(lambda x: chr(x), consQs))


def _getOverlap(dna1, dna2, readLen, possibleInsertList, minOverlapIdentity):
    """"
        Given a list of possible insert sizes and required minimum overlap, return the first insert size satisfying
        the condition for an overlap.

        @type minOverlapIdentity: float
        @return: the insert size of an allowed overlap or None
    """
    for insert in possibleInsertList:
        tryDifferentInsert = False
        mismatch = 0
        overlapSize = 2 * readLen - insert
        maxError = int(round(overlapSize * (1. - minOverlapIdentity)))
        for i in range(overlapSize):
            if dna1[insert - readLen + i] != dna2[i]:
                mismatch += 1
            if mismatch > maxError:
                # print insert, mismatch
                tryDifferentInsert = True
                break
        if not tryDifferentInsert:
            if mismatch <= maxError:
                # print 'mismatch', mismatch
                return insert
    return None


def _getPossibleInsertSizes(minInsert, maxInsert, insert):
    """
        @return: a list of possible insert sizes, starting from the most probable
    """
    assert minInsert <= insert <= maxInsert
    exploreChain = [insert]
    i = 1
    while True:
        add = False
        if minInsert <= insert + i <= maxInsert:
            exploreChain.append(insert + i)
            add = True
        if minInsert <= insert - i <= maxInsert:
            exploreChain.append(insert - i)
            add = True
        if add:
            i += 1
        else:
            break
    return exploreChain


def _testJoin():
    fileTupleList = (('/Users/ivan/Documents/nobackup/hsim01/562/samples/0/NZ_AKLB00000000/0_pair1.fq.gz',
                      '/Users/ivan/Documents/nobackup/hsim01/562/samples/0/NZ_AKLB00000000/0_pair2.fq.gz',
                      '/Users/ivan/Documents/nobackup/hsim01/562/samples/0/NZ_AKLB00000000/0_join.fq.gz',
                      100, 150, 60),)
    print joinPairEnd(fileTupleList)
    # a = DnaOverlapConcensus()


def qsFilter(fileTupleList):
    """
        Filter reads according to the QS cutoffs.
        TODO: not implemented yet, not in production stage

        @param fileTupleList: list of tuples (pair1.fq, pair2.fq, filtered.fq, qsArray, readLen, trimRemain)
    """
    for fq1, fq2, fqFilter, qsArray, readLen, trimRemain in fileTupleList:

        _gsFilterTask(fq1, fq2, fqFilter, qsArray, readLen, trimRemain)
        break  # !!!
    # print fileTupleList[0]


def _getQualityChunk(qs, readLen, qsArray):
    """
        Get the longest part of the read satisfying min quality scores for respective positions.

        TODO: not in production stage

        @param qs: array of quality scores
        @param qsArray: min QS for each position
        @return: (start position, chunk length)
    """
    start = None
    count = 0
    startMax = None
    countMax = 0
    # check QS for each read position
    for i in range(readLen):
        if ord(qs[i]) - 33 >= qsArray[i]:
            if count == 0:
                start = i
            count += 1
        else:
            if count > countMax:
                startMax = start
                countMax = count
            count = 0
    if count > countMax:
        countMax = count
        startMax = start
    if startMax is None:
        startMax = 0
        countMax = 0

    return startMax, countMax


def _gsFilterTask(fq1, fq2, fqFilter, qsArray, readLen, trimRemain):
    """
        TODO: not in production stage
        @param fq1:
        @param fq2:
        @param fqFilter:
        @param qsArray:
        @param readLen:
        @param trimRemain:
        @return:
    """
    # print fq1, fq2, fqFilter, qsArray, readLen, trimRemain

    # out = WriteFq(fqFilter)
    okPair = 0
    allPair = 0
    countListMin = []
    countListMax = []
    countListAvg = []
    avgStartList = []

    threashold = int(round(trimRemain * readLen))  # !!!
    # print threashold
    for r1, r2 in zip(ReadFqGen(fq1), ReadFqGen(fq2)):
        n1, dna1, p, qs1 = r1
        n2, dna2, p, qs2 = r2
        start1, count1 = _getQualityChunk(qs1, readLen, qsArray)
        start2, count2 = _getQualityChunk(qs2, readLen, qsArray)

        if count1 >= threashold and count2 >= threashold:
            # keep else throw away both pairs!
            okPair += 1
        allPair += 1

        countListMin.append(min(count1, count2))
        countListMax.append(max(count1, count2))
        countListAvg.append(count1)
        countListAvg.append(count2)
        avgStartList.append(start1)
        avgStartList.append(start2)
        # countList.append(count2)
        # print zip(map(lambda x: ord(x) - 33, qs1), range(100)) ,'\n'

        # print qsArray
        # break
    print round((float(okPair) / float(allPair)) * 100., 3), '%', okPair, allPair
    print np.mean(countListMin), np.std(countListMin)
    print np.mean(countListMax), np.std(countListMax)
    print np.mean(countListAvg), np.std(countListAvg)
    print np.mean(avgStartList), np.std(avgStartList)

    # out.close()
    # return  # TODO: report


def _testFilter():
    fileTupleList = [('/Users/ivan/Documents/nobackup/hsim01/562/samples/0/NZ_AKLB00000000/0_pair1.fq.gz',
                      '/Users/ivan/Documents/nobackup/hsim01/562/samples/0/NZ_AKLB00000000/0_pair2.fq.gz',
                      '/Users/ivan/Documents/nobackup/hsim01/562/samples/0/NZ_AKLB00000000/0_filtered_qs.fq.gz',
                      np.array([11, 11, 11,  8, 11,  9, 10, 10, 10,  9,  9,  9,  9,  9,  9,  9,  9,
                                9,  9,  9,  9,  9,  9,  9,  9,  9,  9, 10, 10, 10, 10, 10, 10, 10,
                                10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,  9, 10,
                                9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,
                                9,  9,  9,  9,  9,  9, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11,
                                11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11], dtype=np.uint16),
                      100, 0.4),
                     ('/Users/ivan/Documents/nobackup/hsim01/562/samples/0/NZ_AKLB00000000/1_pair1.fq.gz',
                      '/Users/ivan/Documents/nobackup/hsim01/562/samples/0/NZ_AKLB00000000/1_pair2.fq.gz',
                      '/Users/ivan/Documents/nobackup/hsim01/562/samples/0/NZ_AKLB00000000/1_filtered_qs.fq.gz',
                      np.array([11, 11, 11,  8, 11,  9, 10, 10, 10,  9,  9,  9,  9,  9,  9,  9,  9,
                                9,  9,  9,  9,  9,  9,  9,  9,  9,  9, 10, 10, 10, 10, 10, 10, 10,
                                10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,  9, 10,
                                9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,
                                9,  9,  9,  9,  9,  9, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11,
                                11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11], dtype=np.uint16),
                      100, 0.4)]
    qsFilter(fileTupleList)
    # pass
    # fq1 = '/Users/ivan/Documents/nobackup/hsim01/562/samples/small.fq.gz'
    # fq1 = '/Users/ivan/Documents/nobackup/hsim01/562/samples/0/NZ_ANLR00000000/0_pair1.fq.gz'
    # fq2 = '/Users/ivan/Documents/nobackup/hsim01/562/samples/0/NZ_ANLR00000000/0_pair2.fq.gz'
    # fq1 = '/Users/ivan/Documents/nobackup/hsim01/562/samples/0_pair1.fq.gz'
    #

    # out = WriteFq('/Users/ivan/Documents/nobackup/hsim01/562/samples/test.fq.gz')

    # out.close()


def _fqReadWriteTest():
    fq1 = '/Users/ivan/Documents/nobackup/hsim01/562/samples/0_pair1.fq.gz'
    # fq1 = '/Users/ivan/Documents/nobackup/hsim01/562/samples/0/NZ_ANLR00000000/0_pair1.fq.gz'

    out = WriteFq('/Users/ivan/Documents/nobackup/hsim01/562/samples/test.fq.gz')
    c = 0
    for e in ReadFqGen(fq1):
        # print i
        out.writeFqEntry(e[0], e[1], e[3])
        c += 1
    print("%s, %s" % (c, c*4))
    out.close()


if __name__ == "__main__":
    # _fqReadWriteTest()
    # _testFilter()
    _testJoin()