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

    Parsing of the input files.
"""

import os
import gzip

import numpy as np

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

from algbioi.com import fq
from algbioi.com import fasta as fas

from algbioi.hsim import pfam


class ReadRec(object):
    def __init__(self, name, dnaSeq, qsArray, protSeq, readingFrameTag, annotStart, annotLen, hmmCoordStart,
                 hmmCoordLen, posCovArray, readRecList):
        """
            A record representing one input read (also a super read or a joined pair-end read).

            @param name: record description
            @param dnaSeq: dna sequence of the read (oriented as +strain)
            @param qsArray: QS array
            @param protSeq: protein sequence
            @param readingFrameTag: reading frame indicator (1..6)

            @param annotStart: start of the HMM annotation (0 based)
            @param annotLen: length of the HMM annotation
            @param hmmCoordStart: start of the HMM annotation within the Gene family (0 based)
            @param hmmCoordLen: length of the annotation withing the Gene family

            @param posCovArray: position coverage array
                (i.e. for each position, how many pair-end reads it cover, e.g. 111222111 for a joined pair-end read)

            @param readRecList: list of ReadRec of which this read consists of (may remove it in the future)
            @type readRecList: list[ReadRec]
        """
        assert annotLen == hmmCoordLen * 3
        assert len(dnaSeq) == len(qsArray)
        self.name = name
        if 4 <= readingFrameTag <= 6:
            self.dnaSeq = str(Seq(dnaSeq, generic_dna).reverse_complement())
            qsArray = qsArray[::-1]  # reverse the quality scores
        else:
            self.dnaSeq = dnaSeq
        self.qsArray = np.zeros(len(qsArray), dtype=np.uint16)
        for i in list(qsArray):
            self.qsArray = ord(i)
        self.protSeq = protSeq
        self.readingFrameTag = readingFrameTag
        self.annotStart = annotStart
        self.annotLen = annotLen
        self.hmmCoordStart = hmmCoordStart
        self.hmmCoordLen = hmmCoordLen
        self.posCovArray = posCovArray
        self.readRecList = readRecList


def readDomblout(inDomblout):
    """
        Read in a dom file, allow only one hit per read.

        @return: map: readName -> (list of tokens, frame-tag)
        @rtype: dict[str,(list[str], int)]
    """
    nameToHit = {}
    for line in gzip.open(inDomblout):
        line = line.strip()
        if line.startswith('#'):
            continue
        assert line.startswith('@')
        tokens = line.split()
        name = tokens[0][1:-2]
        tag = int(tokens[0][-1])

        assert name not in nameToHit
        nameToHit[name] = (tokens, tag)

    return nameToHit


def parse(inFq, inDomtblout, inProtFna, pairEndReadLen):
    """
        Read in joined pair-end reads and its HMM annotation from the FASTQ, prot FASTA, and DOMTBLOUT files.

        @param inFq: FASTQ file containing joined pair-end reads
        @param inDomtblout: HMM annotation file
        @param inProtFna: corresponding prot sequence
        @param pairEndReadLen: length of one pair-end read (i.e. its length, not the insert size)

        @return: a list of read records
        @rtype: list[ReadRec]
    """
    recList = []

    # read in prot sequences
    nameToProtSeq = fas.fastaFileToDictWholeNames(inProtFna)

    # read in dom file
    nameToDom = readDomblout(inDomtblout)

    assert len(nameToProtSeq) == len(nameToDom)

    # read in pair end reads, create ReadRec
    for readName, dna, p, qs in fq.ReadFqGen(inFq):

        protSeq = nameToProtSeq[readName]

        hit, frameTag = nameToDom[readName]

        annotStart, annotLen, strain, score, acc = pfam.dnaHitInfo(hit, dna, protSeq)

        # TODO: consider also strain, score, acc ?

        hmmCoordStart = int(hit[19]) - 1
        hmmCoordLen = int(hit[20]) - hmmCoordStart

        assert annotLen == 3 * hmmCoordLen

        # create the coverage array
        posCovArray = np.zeros(len(dna), dtype=np.uint16)
        start = len(dna) - pairEndReadLen  # incl.
        end = pairEndReadLen - 1  # incl.
        for i in range(len(posCovArray)):
            if start <= i <= end:
                posCovArray[i] = 2
            else:
                posCovArray[i] = 1

        readRecList = None  # this read record did not originate by the composition of several read records
        recList.append(ReadRec(readName, dna, qs, protSeq, frameTag, annotStart, annotLen, hmmCoordStart, hmmCoordLen,
                               posCovArray, readRecList))

    return recList


if __name__ == "__main__":
    # test !
    baseDir = '/Users/ivan/Documents/nobackup/hsim01/562/samples/0/sample_partitioned'
    baseName = 'duf2158_3'
    # baseName = 'aminotran_3_1'
    parse(inFq=os.path.join(baseDir, 'o_' + baseName + '_join.fq.gz'),
                    inDomtblout=os.path.join(baseDir, 'o_' + baseName + '_join_prot.domtblout.gz'),
                    inProtFna=os.path.join(baseDir, 'o_' + baseName + '_join_prot.fna.gz'),
                    pairEndReadLen=150)