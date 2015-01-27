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

    Contains basic functionality to work with SAM files.
"""
import os
import sys
import gzip
import numpy as np
import parallel
import fasta as fas
import csv


def getErrorStatAndQSCutoff(inFileTupleList, outFilePathProfile, outFilePathQSCutoff=None, maxCpu=2):
    """
        Get error statistics and QS cutoffs from a SAM and reference FASTA files.

        http://samtools.github.io/hts-specs/SAMv1.pdf

        @param inFileTupleList: list of tuples (SAM file path, reference FASTA file path, readLen, qsArrayLen)
        @type inFileTupleList: list
        @param outFilePathProfile: the error profile will be stored in this file
        @param outFilePathQSCutoff: for each cutoff, contains a QS for each position
        @param maxCpu: max cpu to be used
    """
    # count errors and total counts, read vs. reference - define tasks
    countErrTaskList = []
    qsArrayLenD = {}
    readLenSet = set()
    for samFile, fastaFile, readLen, qsArrayLen in inFileTupleList:
        # define task
        countErrTaskList.append(parallel.TaskThread(_getErrorMAllMReadCount, (samFile, fastaFile, readLen, qsArrayLen)))
        # store all read lengths
        readLenSet.add(readLen)
        # store QS array lengths for individual read lengths
        if readLen not in qsArrayLenD:
            qsArrayLenD[readLen] = qsArrayLen
        elif qsArrayLenD[readLen] < qsArrayLen:
            qsArrayLenD[readLen] = qsArrayLen

    # define the error count and all count matrices for different read lengths (row ~ QS, col ~ read position)
    pqeD = {}  # error count (for a particular read length)
    pqaD = {}  # all count (..)
    readCountD = {}  # read count (..)
    insertStatListD = {}
    cigarReportListD = {}
    # get the matrices
    for readLen in readLenSet:
        pqeD[readLen], pqaD[readLen] = getMatrixList(qsArrayLenD[readLen], readLen, matrixCount=2)
        readCountD[readLen] = 0
        insertStatListD[readLen] = []
        cigarReportListD[readLen] = []

    # derive statistics from the SAM files, run in parallel
    rl = parallel.runThreadParallel(countErrTaskList, maxCpu)

    # for each read length, add up all values: error count, all count, read count
    # collect insert size and std for each SAM file
    for pqe, pqa, readCount, readLen, insertSizeMean, insertSizeStd, cigarReport in rl:
        insertStatListD[readLen].append(' %s:%s' % (round(insertSizeMean, 1), round(insertSizeStd, 1)))
        cigarReportListD[readLen].append(cigarReport)
        for i in range(qsArrayLenD[readLen]):
            for j in range(readLen):
                pqeD[readLen][i][j] += pqe[i][j]
                pqaD[readLen][i][j] += pqa[i][j]
        readCountD[readLen] += readCount

    # output profile file
    out = csv.OutFileBuffer(outFilePathProfile)

    # write insert sizes and std for each SAM file (and readLen)
    out.writeText('# For each SAM file: insert : std\n')
    for readLen in readLenSet:
        out.writeText('rl: %s, %s\n' % (readLen, ','.join(insertStatListD[readLen])))
    out.writeText('-\n')

    # write cigar report for each SAM file (and readLen)
    out.writeText('# For each SAM file: cigar report\n')
    for readLen in readLenSet:
        out.writeText('rl: %s, %s\n' % (readLen, ', '.join(cigarReportListD[readLen])))
    out.writeText('-\n')

    # output for each read length
    for readLen in readLenSet:
        pqe = pqeD[readLen]
        pqa = pqaD[readLen]
        # compute overall error
        error = 0
        allEntry = 0
        for i in range(qsArrayLenD[readLen]):
            for j in range(readLen):
                error += pqe[i][j]
                allEntry += pqa[i][j]

        # write overall error
        out.writeText('# Overall Substitution error (readLen: %s)\n' % readLen)
        out.writeText('%s\n' % round((float(error) / float(allEntry)) * 100., 3))
        out.writeText('-\n')

        # write per QS error
        out.writeText('# Per QS error (readLen: %s)\n' % readLen)
        out.writeText(_getProfilePrint(pqe, pqa, qsArrayLenD[readLen], readLen, readCountD[readLen]))
        out.writeText('\n-\n')

        # compute cumulative values
        for i in range(qsArrayLenD[readLen] - 2, -1, -1):
            for j in range(readLen):
                pqe[i][j] += pqe[i+1][j]
                pqa[i][j] += pqa[i+1][j]

        # write cumulative error
        out.writeText('# Cumulative error (readLen: %s)\n' % readLen)
        out.writeText(_getProfilePrint(pqe, pqa, qsArrayLenD[readLen], readLen, readCountD[readLen]))
        out.writeText('\n-\n')

    out.close()

    # get QS cutoffs for different errors cutoffs
    if outFilePathQSCutoff is not None:
        out = csv.OutFileBuffer(outFilePathQSCutoff)
        out.writeText('# readLen, cutoff, min support, lowest QS pos 0,\n')
        # different read lengths
        for readLen in readLenSet:
            pqe = pqeD[readLen]
            pqa = pqaD[readLen]
            # compute qs cutoff for different error cutoffs
            for cutoff in np.arange(0.01, 1.002, 0.002):
                line = [str(readLen), str(cutoff), str(None)]
                support = sys.maxint
                # different positions
                for j in range(readLen):
                    # entry format: readLen, cutoff (error), min support
                    # pos 0 lowest QS, pos 1 lowest QS, ..., pos readLen - 1 lowest QS
                    # where: qs = 0 .. max QS - 1 or None
                    qs = None
                    for i in range(len(pqe)):
                        if pqa[i][j] > 0:
                            if (float(pqe[i][j]) / float(pqa[i][j])) * 100. <= cutoff:
                                qs = i
                                support = min(support, pqa[i][j])
                                break
                    if qs is None:
                        support = None
                    line.append(str(qs))
                if support is not None and support != sys.maxint:
                    line[2] = str(round((float(support) / float(readCountD[readLen])) * 100., 2))
                out.writeText(', '.join(line) + '\n')
        out.close()


def readCutoffArray(cutoffFile, cutoff, readLen):
    """
        For each read position, get min QS, s.t. the probability of error at this position is at most cutoff (in %).

        @param cutoffFile: file containing QS cutoffs for different errors and read lengths
        @param cutoff: max error (%) allowed per read position
        @return: array of min QS for each read position corresponding to the given error (or None)
    """
    lastLine = None
    for line in open(cutoffFile):
        line = line.strip()
        # line not a comment or empty string
        if not line.startswith('#') and len(line) != 0:
            readLenF, cutoffF, supportF, valF = line.split(',', 3)
            cutoffF = float(cutoffF)
            # correct readLen, support for all read positions defined
            if int(readLenF) == readLen and supportF.strip() != 'None':
                # the entry corresponding to the cutoff found
                if np.isclose(float(cutoffF), cutoff):
                    return np.array(map(lambda x: int(x), valF.split(',')), dtype=np.uint16)
                # the entry corresponding to the exact cutoff not found, take the last (more strict) one
                elif cutoffF > cutoff:
                    # take the previous one
                    if lastLine is None:
                        return None
                    else:
                        readLenF, cutoffF, supportF, valF = lastLine.split(',', 3)
                        sys.stderr.write('sam.readCutoffArray: cutoff %s instead of %s taken' % (cutoffF, cutoff))
                        return np.array(map(lambda x: int(x), valF.split(',')), dtype=np.uint16)
                else:
                    lastLine = line
    # no suitable entry found
    return None


def _getErrorMAllMReadCount(samFilePath, fastaFilePath, readLen, qsArrayLen):
    """
        Get the matrix for error counts and for all counts of particular QS and position within a read,
        and a number of reads, (row ~ QS, col ~ read positions).

        @param qsArrayLen: length of the QS array (QSs are mapped onto 0.. qsArrayLen -1)
        @return: (matrix error count, matrix all count, readCount, readLen, insertSizeMean, insertSizeStd, cigarReport)
        @rtype: tuple
    """
    try:
        assert os.path.isfile(samFilePath)
        assert os.path.isfile(fastaFilePath)

        # get two zeroed matrices
        pqe, pqa = getMatrixList(qsArrayLen, readLen, matrixCount=2)

        # read in reference sequences
        seqNameToSeq = fas.fastaFileToDictWholeNames(fastaFilePath)

        # read the SAM file read by read and count error positions and all positions for a given QS
        readCount = 0

        # list of insert sizes
        insertSizeList = []

        # initialization for the cigar string parsing
        x = i = d = 0
        digitsCharList = map(lambda y: str(y), range(10))
        cigarSymbolCharList = ['=', 'X', 'I', 'D']

        if samFilePath.endswith('.gz'):
            samOpen = gzip.open(samFilePath)
        else:
            samOpen = open(samFilePath)
        for line in samOpen:
            if not line.startswith('@'):
                tokens = line.split('\t')
                if len(tokens) >= 11:
                    # qname = tokens[0]  # query read name
                    # check the flag
                    flag = int(tokens[1])
                    assert flag & 0x900 == 0  # it's the primary line of the read
                    assert flag & 0x1 != 0  # template having multiple segments in sequencing
                    assert flag & 0x2 != 0  # each segment properly aligned according to the aligner
                    assert flag & 0x4 == 0 and flag & 0x8 == 0  # mapping exists (also for the next segment)
                    # in a read pair, only one read is a reverse complement
                    assert (flag & 0x10 != 0 and flag & 0x20 == 0) or (flag & 0x10 == 0 and flag & 0x20 != 0)
                    # start position within the reference (zero based)
                    startPos = int(tokens[3]) - 1
                    # reference sequence name
                    rname = tokens[2]
                    # corresponding reference sequence
                    refSeq = seqNameToSeq[rname]
                    refSeq = refSeq[startPos:(startPos + readLen)]
                    # the cigar string, i.e. substitutions, deletions, insertions
                    cigar = tokens[5]
                    buff = ''
                    for e in list(cigar):
                        if e in cigarSymbolCharList:
                            n = int(buff)
                            buff = ''
                            if e == '=':
                                pass  # eq += n
                            elif e == 'X':
                                x += n
                            elif e == 'I':
                                i += n
                            elif e == 'D':
                                d += n
                        else:
                            assert e in digitsCharList
                            buff += e

                    # read - always as if it was on the positive strain
                    read = tokens[9]
                    # quality score array
                    qual = tokens[10]  # qual = np.fromstring(tokens[10], dtype=np.uint8)
                    readCount += 1

                    # read is a reverse complement
                    if flag & 0x10 != 0:
                        revCompl = True
                        # add the insert size of the pair (where only one in the pair is the reverse complement)
                        insertSizeList.append(abs(int(tokens[8])))
                    else:
                        revCompl = False

                    # count errors and all positions
                    for j in range(len(read)):
                        i = ord(qual[j]) - 33
                        if i > qsArrayLen - 1:
                            sys.stderr('Read from SAM: QS "%i" too high, changed to max="%s" \n' % (i, qsArrayLen - 1))
                            i = qsArrayLen - 1
                        # read is reverse complement (but in the SAM file it is annotated as if it was on the + strain)
                        if revCompl:
                            posErrProfile = readLen - 1 - j
                        else:
                            posErrProfile = j
                        if read[j] != refSeq[j]:
                            pqe[i][posErrProfile] += 1
                        pqa[i][posErrProfile] += 1

        allPositions = float(readLen * readCount)

        cigarReport = str("S:%s%% I%s D%s" % (round((float(x) / allPositions) * 100., 2),
                                              round((float(i) / allPositions), 4),
                                              round((float(d) / allPositions), 4)))
    except Exception as e:
        print type(e)
        print e.message
        print e.args
        raise e

    return pqe, pqa, readCount, readLen, np.mean(insertSizeList), np.std(insertSizeList), cigarReport


def _getProfilePrint(pqe, pqa, qsArrayLen, readLen, readCount):
    """
        Get a string representing an error profile (row ~ QS, col ~ position), (error %, overall %)

        @param pqe: error matrix (a number of reads having an error at this position and QS)
        @param pqa: all matrix (a number of all reads having this QS at this position)
        @param qsArrayLen: length of quality scores (i.e. QS ~ 0 .. max quality score - 1)
        @rtype: str
    """
    lineList = []
    for i in range(qsArrayLen):
        line = ""
        for j in range(readLen):
            if pqa[i][j] > 0:
                line += '%s : %s, ' % (round((float(pqe[i][j]) / float(pqa[i][j])) * 100., 3),
                                       round((float(pqa[i][j]) / float(readCount)) * 100., 1))
            else:
                line += ', '
                assert pqe[i][j] == pqa[i][j] == 0
        lineList.append(line)
    return '\n'.join(lineList)


def getMatrixList(rowNum, colNum, matrixCount=1):
    """
        Get a list of zeroed matrices of the same given dimensions.
    """
    matrixList = []
    for e in range(matrixCount):
        matrixList.append(np.zeros((rowNum, colNum), dtype=np.int64))
    return matrixList
