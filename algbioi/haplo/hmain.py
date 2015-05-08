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

    Implements the main function of the package.
"""

# import os
# import gzip
import numpy as np

from algbioi.haplo import hio
from algbioi.haplo.snowball import runSnowball


def getHmmCovArray(recList):
    """
        Returns an array representing the coverage, i.e. how the HMM "dna" alignment is covered by the HMM annot. reads.

        @type recList list[algbioi.haplo.read_rec.ReadRec]

        @rtype: list[int]
    """
    # find out the alignment length considering all read-records
    aliLen = 0
    for rec in recList:
        aliLen = max(aliLen, 3 * rec.hmmCoordStart + rec.annotLen)

    # alignment coverage array
    aliCovArray = np.zeros(aliLen, np.uint32)

    # go over all records
    for rec in recList:
        hmmCoordStart = rec.hmmCoordStart * 3

        # go over all alignment positions, count the coverage
        j = 0
        for i in range(rec.annotStart, rec.annotStart + rec.annotLen):
            aliCovArray[hmmCoordStart + j] += rec.posCovArray[i]
            j += 1

    return aliCovArray


def findSeed(recList, aliCovArray):
    """
        Finds the seed, i.e. records that covers the most alignment positions with the highest coverage.

        For all annotated positions of the read, sums up the alignment coverage at the corresponding positions.
        Returns a list sorted according this sum (the biggest first).

        @type recList: list[algbioi.haplo.read_rec.ReadRec]
        @type aliCovArray: list[int]

        @return: seed list, seeds with the highest support first
        @rtype: list[algbioi.haplo.read_rec.ReadRec]
    """
    recSumList = []
    # compute the overlap sum of each read-record
    for rec in recList:
        hmmCoordStart = rec.hmmCoordStart * 3
        covSum = 0
        j = 0
        for i in range(rec.annotStart, rec.annotStart + rec.annotLen):
            covSum += aliCovArray[hmmCoordStart + j]
            j += 1

        recSumList.append((rec, covSum))

    # sort according to the longest annotation
    # recSumList.sort(key=lambda x: x[0].annotLen, reverse=True)

    # sort according to the overlap sum
    recSumList.sort(key=lambda x: x[1], reverse=True)

    # return only the records
    retList = []
    for rec in recSumList:
        retList.append(rec[0])
    return retList


def buildSuperReads(inFq, inDomtblout, pairEndReadLen, inProtFna=None, outFile=None, translTable=11,
                    maxMismatchQSAllowed=9, minScoreReqiured=75, minAnnotOverlapScore=40, scoreStopSearch=100,
                    maxQS=94):

    # TODO: make prot sequence voluntary

    # 1. read in read records
    recList = hio.parse(inFq, inDomtblout, inProtFna, pairEndReadLen)

    # 2. sort read record list according to the HMM start positions
    recList.sort(key=lambda x: x.hmmCoordStart)

    # 3. build the alignment coverage array
    aliCovArray = getHmmCovArray(recList)

    # 4. find the hotspot (starting seed)
    seedList = findSeed(recList, aliCovArray)

    # print("Record list input length: %s" % len(recList))

    # 5. run the snowball alg.
    recSet = runSnowball(recList, seedList, translTable, maxMismatchQSAllowed, minScoreReqiured,
                         minAnnotOverlapScore, scoreStopSearch, maxQS)

    if outFile is not None:
        hio.storeReadRec(list(recSet), outFile)

    return recSet


# TODO: 6. clean up discart low support contigs, output results

# TODO: 7. evaluate results !!!

if __name__ == "__main__":
    pass