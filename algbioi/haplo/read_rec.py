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

    Represents a (super-)read record.
"""
import copy
import numpy as np

from Bio.Alphabet import generic_dna
from Bio.Seq import Seq


class ReadRec(object):
    def __init__(self, recordId, dnaSeq, qsArray, protSeq, readingFrameTag, annotStart, annotLen, protStart, protLen,
                 hmmCoordStart, hmmCoordLen, posCovArray, translTable=11):
        """
            A record representing one read (also a super-read or a joined pair-end read).
            The constructor is used to create an initial record, the function "merge records" of this module
            is then used to create super-read-records by merging of read-records.

            @param recordId: record identifier
            @param dnaSeq: dna sequence of the (super-)read (oriented as +strain)
            @param qsArray: QS array
            @param protSeq: protein sequence (may drop this?)
            @param readingFrameTag: reading frame indicator (1..6), (for super-reads: 1..3)

            @param annotStart: start of the HMM annotation within the read (0 based)
            @param annotLen: length of the HMM annotation
            @param protStart: start of the HMM annotation on the prot sequence (0 based) (or None)
            @param protLen: length of the HMM annotation on the prot sequence (or None)

            @param hmmCoordStart: start of the HMM annotation within the Gene family (0 based)
            @param hmmCoordLen: length of the annotation withing the Gene family (or None)

            @param posCovArray: position coverage array (i.e. for each position, how many pair-end reads cover it,
                                                        e.g. 111222111 for a joined pair-end read)
            @param translTable: a translation table used to translate dna to prot sequences

            @type recordId: str
            @type dnaSeq: str
            @type qsArray: str
            @type protSeq: str
            @type readingFrameTag: int
            @type annotStart: int
            @type annotLen: int
            @type protStart: int
            @type protLen: int
            @type hmmCoordStart: int
            @type hmmCoordLen: int
            @type posCovArray: ndarray
            @type translTable: int
        """
        assert len(dnaSeq) == len(qsArray)
        assert annotLen == 3 * protLen

        # the gene family annotation is on the reverse strain, orient it as if it was on the positive strain
        if 4 <= readingFrameTag <= 6:
            dnaSeq = str(Seq(dnaSeq, generic_dna).reverse_complement())  # dna reverse complement
            qsArray = qsArray[::-1]  # reverse the QS array
            annotStart = len(dnaSeq) - annotStart - annotLen  # reverse the start annotation position

        self.qsArray = qsArray
        # QS array as a numpy array if needed.. uncomment
        # self.qsArray = np.zeros(len(qsArray), dtype=np.uint8)
        # for i in range(len(qsArray)):
        #     self.qsArray[i] = ord(qsArray[i]) - 33
        self.recordId = recordId
        self.dnaSeq = dnaSeq
        self.protSeq = protSeq
        self.readingFrameTag = readingFrameTag
        self.annotStart = annotStart
        self.annotLen = annotLen
        self.protStart = protStart
        self.protLen = protLen
        self.hmmCoordStart = hmmCoordStart
        self.hmmCoordLen = hmmCoordLen
        self.posCovArray = posCovArray

        # a list of ReadRec merging of which this read-record originated
        # None - if this read-record did not originate by the composition of other read-records
        # (may remove in the future)
        self.readRecList = None

        # a list of (overlapIdx, score) used to merge the read-records (list of (idx, score), i.e. list[(int, int)]
        self.overlapIdxList = None

        # the position of this (super-)read within a bigger contig
        self.posWithinContig = 0

        # TODO: remove this checking
        # p = self.protSeq[self.protStart:self.protStart + self.protLen]
        # d = str(Seq(self.dnaSeq[self.annotStart:self.annotStart + self.annotLen], generic_dna).translate(translTable))
        # assert p == d

        # TODO: remove: checking of the stop codons withing the annotation
        # for i in range(self.annotStart, self.annotStart + self.annotLen, 3):
        #     if (self.dnaSeq[i:i+3] == 'TAG' or self.dnaSeq[i:i+3] == 'TAA' or self.dnaSeq[i:i+3] == 'TGA'):
        #         print self.dnaSeq[i:i+3], i, self.annotStart, self.annotStart + self.annotLen, \
        #             self.protSeq[self.protStart: self.protStart + self.protLen], \
        #             str(self.qsArray[i:i+3])

    def getMostCoveredRec(self):
        """
            Get read-records (of which this read-record consist of) sorted according to how they are covered
            by the other reads, the most covered first.

            @return: list of the read-records, most covered first (or None if this read-record consists only of itself)
            @rtype: list[ReadRec]
        """
        if self.readRecList is None:
            return None
        else:
            recCovList = []
            for rec in self.readRecList:
                covSum = 0
                offset = rec.posWithinContig
                for i in range(len(rec.dnaSeq)):
                    covSum += self.posCovArray[offset + i]
                recCovList.append((rec, covSum))
            recCovList.sort(key=lambda x: x[1], reverse=True)
        return map(lambda x: x[0], recCovList)


def mergeRecords(rec1, rec2, overlapIdx, consSeqGen, score=None):
    """
        Merges record-1 into record-2, returns a merged object.

        @param rec1: (super-)read-record
        @param rec2: read-record (that hasn't been joined with any other so far)
        @param overlapIdx: position within rec2, start pos of rec1 mapped onto rec2 (it can be a negative value)
        @param consSeqGen: generator of consensus sequences and consensus QS arrays
        @param score: score used for this overlap

        @type rec1: ReadRec
        @type rec2: ReadRec
        @type overlapIdx: int
        @type consSeqGen: algbioi.com.fq.QsMultMatrix
        @type score: int
        @rtype: ReadRec

        @return: a merged read record
    """
    # make a shallow copy of the object
    mRec = copy.copy(rec1)

    # concatenate record ids
    mRec.recordId = rec1.recordId + '|' + rec2.recordId

    # lengths of the sequences
    dnaLen1 = len(rec1.dnaSeq)
    dnaLen2 = len(rec2.dnaSeq)

    # get start positions withing the dna sequences and corresponding lengths
    if overlapIdx <= 0:

        # part before overlap rec1
        boPos1 = 0
        boLen1 = - overlapIdx

        # part before overlap rec2
        boPos2 = 0
        boLen2 = 0

        # overlap start pos within rec1, rec2
        oPos1 = - overlapIdx
        oPos2 = 0

    else:

        # part before overlap rec1
        boPos1 = 0
        boLen1 = 0

        # part before overlap rec2
        boPos2 = 0
        boLen2 = overlapIdx

        # overlap start pos within rec1, rec2
        oPos1 = 0
        oPos2 = overlapIdx

    # length of the dna overlap
    oLen = min(dnaLen1 - oPos1, dnaLen2 - oPos2)

    # after overlap rec1
    aPos1 = oPos1 + oLen
    aLen1 = dnaLen1 - aPos1

    # after overlap rec2
    aPos2 = oPos2 + oLen
    aLen2 = dnaLen2 - aPos2

    dnaCons, qsCons = consSeqGen.getConsensus(rec1.dnaSeq[oPos1:oPos1+oLen],
                                              rec2.dnaSeq[oPos2:oPos2+oLen],
                                              rec1.qsArray[oPos1:oPos1+oLen],
                                              rec2.qsArray[oPos2:oPos2+oLen])

    mRec.dnaSeq = rec1.dnaSeq[boPos1:boPos1+boLen1] + rec2.dnaSeq[boPos2:boPos2+boLen2] + dnaCons \
                  + rec1.dnaSeq[aPos1:aPos1+aLen1] + rec2.dnaSeq[aPos2:aPos2+aLen2]

    mRec.qsArray = rec1.qsArray[boPos1:boPos1+boLen1] + rec2.qsArray[boPos2:boPos2+boLen2] + qsCons \
                   + rec1.qsArray[aPos1:aPos1+aLen1] + rec2.qsArray[aPos2:aPos2+aLen2]

    assert len(mRec.dnaSeq) == len(mRec.qsArray)

    # position coverage array
    mRec.posCovArray = np.zeros(len(mRec.dnaSeq), dtype=np.uint32)

    # copy before overlap coverage
    if boLen1 > 0:
        assert boLen2 == 0
        for i in range(boLen1):
            mRec.posCovArray[i] = rec1.posCovArray[i]
    else:
        assert boLen1 == 0
        for i in range(boLen2):
            mRec.posCovArray[i] = rec2.posCovArray[i]

    # sum up overlap coverage
    offset = max(boLen1, boLen2)
    assert min(boLen1, boLen2) == 0
    for i in range(oLen):
        mRec.posCovArray[offset + i] = rec1.posCovArray[oPos1 + i] + rec2.posCovArray[oPos2 + i]

    # copy after overlap coverage
    offset += oLen
    if aLen1 > 0:
        assert aLen2 == 0
        for i in range(aLen1):
            mRec.posCovArray[offset + i] = rec1.posCovArray[aPos1 + i]
    else:
        assert aLen1 == 0
        for i in range(aLen2):
            mRec.posCovArray[offset + i] = rec2.posCovArray[aPos2 + i]

    if overlapIdx > 0:
        mRec.readingFrameTag = rec2.readingFrameTag
    else:
        mRec.readingFrameTag = rec1.readingFrameTag
    if mRec.readingFrameTag > 3:
        mRec.readingFrameTag -= 3

    # Hmm annotation start positions
    anPos1 = rec1.annotStart + boLen2
    anPos2 = rec2.annotStart + boLen1

    # Hmm annotation end + 1 positions
    anPosEndA1 = anPos1 + rec1.annotLen
    anPosEndA2 = anPos2 + rec2.annotLen

    # the annotations must have an overlap
    assert max(anPos1, anPos2) < min(anPosEndA1, anPosEndA2)

    # get the union of the annotations
    mRec.annotStart = min(anPos1, anPos2)
    mRec.annotLen = max(anPosEndA1, anPosEndA2) - mRec.annotStart

    # Hmm coord within gene family, take coordinates of the one at the interval boundary or the one with highest support
    if anPos1 < anPos2:
        mRec.hmmCoordStart = rec1.hmmCoordStart
    elif anPos1 > anPos2:
        mRec.hmmCoordStart = rec2.hmmCoordStart
    elif rec1.annotLen >= rec2.annotLen:
        mRec.hmmCoordStart = rec1.hmmCoordStart
    else:
        mRec.hmmCoordStart = rec2.hmmCoordStart

    # a list of all read-records this merged read-record consist of, also the overlap index and score
    assert rec2.readRecList is None and rec2.overlapIdxList is None
    if rec1.readRecList is None:
        assert rec1.overlapIdxList is None
        mRec.readRecList = [rec1, rec2]
        mRec.overlapIdxList = [(overlapIdx, score)]
        # set the positions of the read-records within the merged record
        rec1.posWithinContig = boLen2
    else:
        mRec.readRecList = rec1.readRecList + [rec2]
        mRec.overlapIdxList = rec1.overlapIdxList + [(overlapIdx, score)]
        # set the positions of the read-records within the merged record
        if boLen2 > 0:
            for r in rec1.readRecList:
                r.posWithinContig += boLen2

    rec2.posWithinContig = boLen1

    # TODO: may drop these !!!
    mRec.protSeq = None
    mRec.protStart = None
    mRec.protLen = None
    mRec.hmmCoordLen = None

    return mRec
