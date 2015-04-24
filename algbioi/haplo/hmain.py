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

import os

import numpy as np

from algbioi.haplo import in_parse


def getHmmCovArray(recList):
    """
        Returns an array representing the coverage, i.e. how the HMM "dna" alignment is covered by the HMM annot. reads.

        @type recList list[in_parse.ReadRec]

        @rtype: list[int]
    """
    # find out the alignment length considering all read-records
    aliLen = 0
    for rec in recList:
        aliLen = max(aliLen, (rec.hmmCoordStart + rec.hmmCoordLen) * 3)

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
        Finds the seed, i.e. a read that covers the most alignment positions with the highest coverage.

        For all annotated positions of the read, sums up the alignment coverage at the corresponding positions.
        Returns a list sorted according this sum (the biggest first.)

        @type recList: list[in_parse.ReadRec]
        @type aliCovArray: list[int]

        @return: seed list, seeds with the highest support first
        @rtype: list[in_parse.ReadRec]
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

    # sort according to the overlap sum
    recSumList.sort(key=lambda x: x[1], reverse=True)

    # return only the records
    retList = []
    for rec in recSumList:
        retList.append(rec[0])
    return retList


def runSnowball(recList, recSeedList):
    """
        Run the main algorithm.

        @param recList: read-record list sorted according to the start alignment positions
        @param recSeedList: read-record list sorted according to the biggest overlap with all the annotations

        @type recList: list[in_parse.ReadRec]
        @type recSeedList: list[in_parse.ReadRec]

        @return: TODO:
        @rtype: list[in_parse.ReadRec]

        # TODO: move it snowball.py
    """
    pass
    # TODO: continue here!


def buildSuperReads(inFq, inDomtblout, inProtFna, pairEndReadLen):

    # 1. read in read records
    recList = in_parse.parse(inFq, inDomtblout, inProtFna, pairEndReadLen)

    # 2. sort read record list according to the HMM start positions
    recList.sort(key=lambda x: x.annotStart)

    # 3. build the alignment coverage array
    aliCovArray = getHmmCovArray(recList)

    # 4. find the hotspot (starting seed)
    seedList = findSeed(recList, aliCovArray)

    # 5. run the snowball alg. TODO: implement
    runSnowball(recList, seedList)

    # TODO: 6. output results

    # TODO: 7. evaluate results !!!


if __name__ == "__main__":
    baseDir = '/Users/ivan/Documents/nobackup/hsim01/562/samples/0/sample_partitioned'
    buildSuperReads(inFq=os.path.join(baseDir, 'o_aminotran_3_1_join.fq.gz'),
                    inDomtblout=os.path.join(baseDir, 'o_aminotran_3_1_join_prot.domtblout.gz'),
                    inProtFna=os.path.join(baseDir, 'o_aminotran_3_1_join_prot.fna.gz'),
                    pairEndReadLen=150)

    # TODO: explore tha "translation table 11", print out all combinations ! what are the stop codons ???

    pass