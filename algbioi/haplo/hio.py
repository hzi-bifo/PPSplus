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
import cPickle

import numpy as np

from algbioi.com import fq
from algbioi.com import fasta as fas
from algbioi.haplo import read_rec
from algbioi.hsim import pfam


def readDomblout(inDomblout):
    """
        Read in a dom file, allow only one hit per read.

        @return: map: readName -> (list of line tokens, frame-tag)
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


def parse(inFq, inDomtblout, inProtFna, pairEndReadLen):  # TODO: remove pair-end read length, inProtFna=None
    """
        Read in joined pair-end reads and its HMM annotation from the FASTQ, prot FASTA, and DOMTBLOUT files.

        @param inFq: FASTQ file containing joined pair-end reads
        @param inDomtblout: HMM annotation file
        @param inProtFna: corresponding prot sequence (can be None)
        @param pairEndReadLen: length of one pair-end read (i.e. its length, not the insert size)

        @return: a list of read records
        @rtype: list[read_rec.ReadRec]
    """
    recList = []

    # read in prot sequences
    nameToProtSeq = fas.fastaFileToDictWholeNames(inProtFna)  # TODO: inProtFna can be None

    # read in dom file
    nameToDom = readDomblout(inDomtblout)

    assert len(nameToProtSeq) == len(nameToDom)

    # read in pair end reads, create ReadRec
    for readName, dna, p, qs in fq.ReadFqGen(inFq):

        readName = readName[1:]

        protSeq = nameToProtSeq[readName]  # TODO: con be None

        hit, frameTag = nameToDom[readName]

        annotStart, annotLen, strain, score, acc = pfam.dnaHitInfo(hit, dna, protSeq)

        if strain == 1:
            assert 1 <= frameTag <= 3
        else:
            assert 4 <= frameTag <= 6

        # TODO: consider also strain, score, acc ?

        protStart = int(hit[19]) - 1
        protLen = int(hit[20]) - protStart

        assert annotLen == 3 * protLen


        # TODO: correct the coordinate start !!!
        hmmCoordStart = int(hit[15]) - 1
        hmmCoordLen = int(hit[16]) - hmmCoordStart

        # create the coverage array
        posCovArray = np.zeros(len(dna), dtype=np.uint8)
        start = len(dna) - pairEndReadLen  # incl.
        end = pairEndReadLen - 1  # incl.
        for i in range(len(posCovArray)):
            if start <= i <= end:
                posCovArray[i] = 2
            else:
                posCovArray[i] = 1

        # TODO: replace qs and posCovArray !!!

        recList.append(read_rec.ReadRec(readName, dna, qs, protSeq, frameTag, annotStart, annotLen, protStart, protLen,
                                        hmmCoordStart, hmmCoordLen, posCovArray))

    return recList


def storeReadRec(recList, outFilePath, compressLevel=1):
    """
        Store a list of read-records to a file.

        @param recList: list of read-records to be stored
        @param outFilePath: output dmp file
        @type recList: list[read_rec.ReadRec]
        @type outFilePath: str
    """
    # open file for writing
    out = gzip.open(outFilePath, 'wb', compressLevel)
    # write the list to the file
    cPickle.dump(recList, out, cPickle.HIGHEST_PROTOCOL)

    out.close()


def loadReadRec(srcFilePath):
    """
        Load a list of records from a file.

        @param srcFilePath: a file containing stored read-records
        @type srcFilePath: str
        @return: list of read-records or None
        @rtype: list[read_rec.ReadRec]
    """
    out = gzip.open(srcFilePath, 'rb')
    try:
        return cPickle.load(out)
    except EOFError:
        return None


def _test():
    # test !
    baseDir = '/home/igregor/Documents/work/hsim/562/samples/0/sample_partitioned'
    # baseName = 'duf2158_3'
    baseName = 'aminotran_3_1'
    recList = parse(inFq=os.path.join(baseDir, 'r_' + baseName + '_join.fq.gz'),
                    inDomtblout=os.path.join(baseDir, 'r_' + baseName + '_join_prot.domtblout.gz'),
                    inProtFna=os.path.join(baseDir, 'r_' + baseName + '_join_prot.fna.gz'),
                    pairEndReadLen=150)
    f = '/home/igregor/Documents/work/hsim/562/samples/0/tmp_test.gz'
    storeReadRec(recList, f)

    recList2 = loadReadRec(f)

    for r1, r2 in zip(recList, recList2):
        assert r1.dnaSeq == r2.dnaSeq
        for i in range(len(r1.posCovArray)):
            assert r1.posCovArray[i] == r2.posCovArray[i]

    print len(recList)
    print len(recList2)

    if __name__ == "__main__":
        _test()