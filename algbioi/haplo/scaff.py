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

    Handle scaffolding.
"""

import os
import numpy as np

from algbioi.hsim import comh
from algbioi.haplo import hio
from algbioi.haplo import heval
from algbioi.haplo import snowball
from algbioi.haplo import join_rec
from algbioi.haplo import read_rec
# from algbioi.com import common as com
# from algbioi.com import fasta as fas


def calculatePairOverlaps(readFq, readDomtblout, contigRecFile, outFile=None,
                          minPOverlap=0.8, minScoreFraction=0.8, minAnnotFraction=0.5,
                          readProtFna=None, readSam=None, tTable=11):
    """
        Calculate the overlap probabilities between assembled contigs (super-reads) and reads from the second library.
        If output-file defined, store, else return.

        @param readFq: FASTAQ file containing reads
        @param readDomtblout: domtblout file containing HMM read alignments
        @param contigRecFile: file containing contigs as read-records
        @param outFile: path to an output file or None.
        @param readSam: read mapping SAM file (or None)

        @type readFq: str
        @type readDomtblout: str
        @type contigRecFile: str
        @type outFile: str | None
        @type minPOverlap: float
        @type minScoreFraction: float
        @type minAnnotFraction: float
        @type readProtFna: str | None
        @type readSam: str | None
        @type tTable: int
        @return: if outFile not defined, return (contig-list, read-list, contig-read-mapping), else None
        @rtype (list[read_rec.ReadRec], list[read_rec.ReadRec], dict[(int,int), (float, float, float, int)]) | None
    """
    # read in read and contig lists
    readList = hio.parse(readFq, readDomtblout, readProtFna)
    contigList = hio.loadReadRec(contigRecFile)

    # get true read mapping (if available)
    if readSam is None:
        readMap = None
    else:
        readMap = heval.getReadTrueMapSimple(readSam)

    # define the triplet mapping if PROT sequences are used for comparison
    # tripletMap = join_rec.getTripletMap(tTable)
    tripletMap = None

    # contig-read overlap mapping
    overlapMap = {}

    i = -1
    for contig in contigList:
        i += 1
        if not contig.isSuperRead():
            continue

        j = -1
        for read in readList:
            j += 1

            # minimum length score and annotation len required for the overlap
            minScore = int(float(read.getLen()) * minScoreFraction)
            minAnnotLen = int(float(read.getLen()) * minAnnotFraction)

            # overlaps to be inspected
            overlapList = snowball.getPossibleOverlaps(contig, read)
            if overlapList is not None:
                offset = 0
                oList = []
                for overlap in overlapList:

                    # count the overlap probability
                    overlapProb, overlapScore, annoLen = snowball.inspectOverlap2(contig, read, overlap, tripletMap,
                                                                                  False, minScore, minAnnotLen,
                                                                                  stopOverlapMaxMismatch=0.05,
                                                                                  continuousOverlap=True)
                    # overlap is sufficient
                    if overlapProb is not None and overlapProb > minPOverlap:
                        oList.append((overlapProb, overlapScore, annoLen, offset))
                    offset += 1

                # get the best overlap
                if len(oList) > 0:
                    oList.sort(key=lambda x: x[0], reverse=True)
                    assert (i, j) not in overlapMap
                    overlapMap[(i, j)] = oList[0]  # overlapProb, overlapScore, annoLen, offset

    # store the record labels
    if readMap is not None:
        for read in readList:
            strainIdR, seqIdR = readMap[read.recordId[:-2]]
            read.labelEval = [heval.LabelR(None, None, strainIdR, seqIdR, None)]

    # store contig-list, read-list, overlap-mapping to an output file
    rEntry = (contigList, readList, overlapMap)
    if outFile is not None:
        hio.storeObj(rEntry, outFile)
        return None
    else:
        return rEntry


def test1():
    gf = 'mfs_2_3'

    # gf = 'acr_tran_3'
    # 5 !!!
    tTable = 11
    baseDir = os.path.join(comh.REFERENCE_DIR_ROOT, '562/samples/5/sample_partitioned')
    readFq = os.path.join(baseDir, 'r_%s_pair.fq.gz' % gf)
    readProtFna = os.path.join(baseDir, 'r_%s_pair_prot.fna.gz' % gf)
    readDomtblout = os.path.join(baseDir, 'r_%s_pair_prot.domtblout.gz' % gf)
    readSam = os.path.join(baseDir, 'r_%s_pair_gmap.sam.gz' % gf)
    contigRecFile = os.path.join(baseDir, 'r_%s_join_read_rec.pkl.gz' % gf)

    outFile = os.path.join(baseDir, 'r_%s_scaff_overlap.pkl.gz' % gf)

    print outFile

    minPOverlap = 0.8
    minScoreFraction = 0.8
    minAnnotFraction = 0.5

    calculatePairOverlaps(readFq, readDomtblout, contigRecFile, outFile, minPOverlap, minScoreFraction,
                          minAnnotFraction, readProtFna, readSam, tTable)


READ_UNIQUE_ARRAY_LEN = 9


def readUniqueMap(scaffOverlapFile, pFrom=0.8, pTo=0.99, step=0.01):
    """
        @param scaffOverlapFile:
        @type scaffOverlapFile: str
        @rtype: ndarray
    """
    # an array of probability thresholds
    pArray = np.arange(pFrom, pTo + step, step, dtype=np.float64)
    arrayLen = len(pArray)

    rArray = np.zeros((READ_UNIQUE_ARRAY_LEN, arrayLen), dtype=np.int64)

    tmpArray = np.zeros(arrayLen, dtype=np.int64)

    # read entries from the file
    contigList, readList, overlapMap = hio.loadObj(scaffOverlapFile)

    # reads
    for j in range(len(readList)):

        tmpArray.fill(0)

        readLenMinScore = 0.9 * float(readList[j].getLen())

        # contigs
        for i in range(len(contigList)):

            if contigList[i].getAvgCov() < 9:  # from 7?
                continue

            if (i, j) in overlapMap:
                overlapProb, overlapScore, annoLen, offset = overlapMap[(i, j)]

                for k in range(arrayLen):

                    if overlapProb > pArray[k] and readLenMinScore > overlapScore:

                        tmpArray[k] += 1

        for k in range(arrayLen):
            if tmpArray[k] > 0:
                found = False
                for l in range(1, READ_UNIQUE_ARRAY_LEN - 1):
                    if tmpArray[k] == l:
                        rArray[l-1][k] += 1
                        found = True
                        break

                if not found:
                    assert tmpArray[k] >= READ_UNIQUE_ARRAY_LEN - 1
                    rArray[READ_UNIQUE_ARRAY_LEN - 2][k] += 1

                rArray[READ_UNIQUE_ARRAY_LEN - 1][k] += 1

                # tmpArray[k]

                # rArray[0][k] += tmpArray[k]  # number of times it matched to a contig
                # rArray[1][k] += 1.  # number of entries (reads matching to a contig)
                # if tmpArray[k] > 1.:
                #     rArray[2][k] += 1.  # ambiguously mapped reads

    return rArray


def readUniqueMapReport(rArrayList, pFrom=0.8, pTo=0.99, step=0.01):
    if len(rArrayList) > 0:
        arrayLen = len(rArrayList[0][0])
        pArray = np.arange(pFrom, pTo + step, step, dtype=np.float64)
        rArray = np.zeros((READ_UNIQUE_ARRAY_LEN, arrayLen), dtype=np.int64)
        for r in rArrayList:
            rArray += r

        buff = ''
        for k in range(arrayLen):
            if rArray[READ_UNIQUE_ARRAY_LEN - 1][k] > 0:

                buff += '%s  ' % (pArray[k] * 100.)
                for l in range(1, READ_UNIQUE_ARRAY_LEN - 1):
                    buff += '%s  ' % round(((float(rArray[l-1][k]) / float(rArray[READ_UNIQUE_ARRAY_LEN - 1][k])) * 100.), 3)
                buff += '%s  ' % (rArray[READ_UNIQUE_ARRAY_LEN - 1][k])
                buff += '\n'
                # buff += '%s%%\t%s\t%s\t%s\n' % (int(pArray[k] * 100.), round(rArray[0][k] / rArray[1][k], 5),
                #                                 int(round(rArray[1][k])), int(round(rArray[2][k])))
        return buff


def checkReadLabels(scaffOverlapFile):

    contigList, readList, overlapMap = hio.loadObj(scaffOverlapFile)

    # contigs
    for i in range(len(contigList)):

        contig = contigList[i]

        if contig.getAvgCov() < 10:
            continue

        # reads
        for j in range(len(readList)):

            read = readList[j]

            if (i, j) in overlapMap:

                readLenMinScore = 0.9 * float(read.getLen())

                overlapProb, overlapScore, annoLen, offset = overlapMap[(i, j)]

                if overlapProb > 0.9 and overlapScore > readLenMinScore:
                    readLabel = read.labelEval[0]
                    k = 0
                    found = False
                    for contigLabel in contig.labelEval:
                        if readLabel.seqId == contigLabel.seqId and readLabel.strainId == contigLabel.strainId:
                            print k, contig.labelEval[0].error, contig.labelEval[0].codonError, contigLabel.error, contigLabel.codonError
                            found = True
                            break
                        k += 1
                    if not found:
                        print 'NF', k
        print '-----'

            # print read.labelEval[0].strainId, read.labelEval[0].seqId


def checkContigLabels(scaffOverlapFile):
    tripletMap = join_rec.getTripletMap(11)
    contigList, readList, overlapMap = hio.loadObj(scaffOverlapFile)
    # contigs
    for contig in contigList:
        try:
            if contig.labelEval is None:
                continue
        except AttributeError:
            continue

        # labels
        print contig.getLen()
        cSet = set()

        for contigLabel in contig.labelEval:

            print contigLabel.error, contigLabel.codonError, contigLabel.strainId, contigLabel.seqId
            cSet.add(contigLabel.getRefDna())

            # print contig.dnaSeq
            # print contigLabel.getRefDna()
            # errors = 0
            # for i, j in zip(list(contig.dnaSeq), list(contigLabel.getRefDna())):
            #     if i != j:
            #         errors += 1
            # print 'err', errors
            # print heval.getError(contig.dnaSeq, contigLabel.getRefDna(), contig.getPosCovArray(), 10, contig.annotStart, contig.annotLen, tripletMap)

        print 'setLen:', len(cSet)
        print '----------'


def test2():

    gf = 'acr_tran_3'
    # gf = 'usher_1'
    # gf = 'aldedh_1'
    # gf = 'aa_permease_2'
    # gf = 'dde_tnp_is66_6'
    # gf = 'peptidase_s6_2'  # differences..
    # gf = 'molybdopterin_1'
    # gf = 'bpd_transp_1_3'
    # gf = 'hth_1_3'
    # gf = 'okr_dc_1_5'
    # gf = 'pfl_3'
    # gf = 'response_reg_1'
    # gf = 'xan_ur_permease_1'
    # gf = 'terminase_gpa_3'
    # gf = 'eal_3'
    # gf = 'amp_binding_4'
    # gf = 'lysr_substrate_2'
    # gf = 'aminotran_3_1'
    # gf = 'mfs_1_3'
    # gf = 'gntp_permease_2'
    # gf = 'sbp_bac_5_3'
    # gf = 'sigma54_activat_1'
    # gf = 'hlyd_2'
    # gf = 'hatpase_c_4'
    # gf = 'deorc_3'
    # gf = 'fe_adh_5'

    scaffOverlapFile = os.path.join(comh.REFERENCE_DIR_ROOT, '562/samples/5/sample_partitioned/r_%s_scaff_overlap.pkl.gz' % gf)

    # checkReadLabels(scaffOverlapFile)

    checkContigLabels(scaffOverlapFile)

    # recListFile = os.path.join(comh.REFERENCE_DIR_ROOT, '562/samples/5/sample_partitioned/r_%s_join_read_rec.pkl.gz' % gf)
    # refGmap = os.path.join(comh.REFERENCE_DIR_ROOT, '562/samples/5/sample_partitioned/r_%s_join_gmap.sam.gz' % gf)
    #
    # refSeqBuff = fas.getSequenceBuffer([os.path.join(comh.REFERENCE_DIR_ROOT, '562', comh.FASTA_GENOMES_DIR_NAME),
    #                    os.path.join(comh.REFERENCE_DIR_ROOT, '562', comh.FASTA_GENOMES_DRAFT_DIR_NAME)])
    # recSet = set(hio.loadReadRec(recListFile))
    # heval.getPerBaseError(recSet, refGmap, refSeqBuff, maxCov=20, translTable=11, allLabels=True)

    # print readUniqueMapReport([readUniqueMap(scaffOverlapFile)])


def test():

    gf = 'mfs_2_3'

    # gf = 'acr_tran_3'
    # 5 !!!
    translTable = 11
    baseDir = os.path.join(comh.REFERENCE_DIR_ROOT, '562/samples/5/sample_partitioned')
    inFq = os.path.join(baseDir, 'r_%s_pair.fq.gz' % gf)
    inProtFna = os.path.join(baseDir, 'r_%s_pair_prot.fna.gz' % gf)
    inDomtblout = os.path.join(baseDir, 'r_%s_pair_prot.domtblout.gz' % gf)
    refGmapR = os.path.join(baseDir, 'r_%s_pair_gmap.sam.gz' % gf)
    refGmapC = os.path.join(baseDir, 'r_%s_join_gmap.sam.gz' % gf)
    readRecFile = os.path.join(baseDir, 'r_%s_join_read_rec.pkl.gz' % gf)

    contigList = hio.loadReadRec(readRecFile)
    readList = hio.parse(inFq, inDomtblout, inProtFna)
    readMap = heval.getReadTrueMapSimple(refGmapR)


    # for c in contigList:
    #     try:
    #         label = c.labelEval
    #         print label.strainId
        # except AttributeError:
        #     label = None
        #     print c.readRecList
    #
    # return

    # print len(contigList), len(readList), len(readMap)
    tripletMap = join_rec.getTripletMap(translTable)

    okCount = 0
    noCount = 0

    for contig in contigList:
        if contig.readRecList is None:
            continue

        strainIdC = contig.labelEval.strainId
        seqIdC = contig.labelEval.seqId
        errorC = contig.labelEval.error
        codonErrorC = contig.labelEval.codonError
        startPosC = contig.labelEval.startPos
        print('LEN: %s' % contig.getLen())

        for read in readList:
            minScore = int(float(read.getLen()) * 0.95)
            minAnnotScore = int(float(read.getLen()) * 0.5)

            strainIdR, seqIdR = readMap[read.recordId[:-2]]

            overlapList = snowball.getPossibleOverlaps(contig, read)
            oList = []
            if overlapList is not None:
                for overlap in overlapList:
                    overlapProb, overlapScore, annotScore = snowball.inspectOverlap2(contig, read, overlap, tripletMap,
                                                                                     False, minScore, minAnnotScore,
                                                                                     stopOverlapMaxMismatch=0.05,
                                                                                     continuousOverlap=True)
                    if overlapProb is not None and overlapProb > 0.9:
                        oList.append((overlapProb, overlapScore, annotScore))
            if len(oList) > 0:

                oList.sort(key=lambda x: x[0], reverse=True)

                overlapProb, overlapScore, annotScore = oList[0]

                if strainIdR != strainIdC or seqIdC != seqIdR:
                    report = 'Different '
                    if strainIdR != strainIdC:
                        report += 'strain %s %s ' % (strainIdR, strainIdC)
                    if seqIdC != seqIdR:
                        report += 'seq %s %s ' % (seqIdC, seqIdR)
                    report += 'error: %s %s %s' % (errorC, codonErrorC, startPosC)
                    noCount += 1
                else:
                    report = 'OK', strainIdC, seqIdC
                    okCount += 1

                print overlapProb, report
        print 'ok, no', okCount, noCount
        print '----------'
        okCount = noCount = 0


    # print okCount, noCount
    heval.assemblyStat(set(contigList))
    heval.assemblyPurity(set(contigList), refGmapC)


if __name__ == "__main__":
    # test()
    # test1()
    test2()