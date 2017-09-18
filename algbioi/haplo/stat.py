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

    Compute assembly statistics.
"""

import os
# import sys
import numpy as np

from algbioi.hsim import comh
from algbioi.haplo import hio, qs

SNP_ERROR_ARRAY_LEN = 10


class StatInfo(object):
    def __init__(self, uniqueA, totalCorrect, snpErrorA, totalNoErrA, strainToBp, strainToBpU):

        self.uniqueA = uniqueA
        self.totalCorrect = totalCorrect
        self.snpErrorA = snpErrorA
        self.totalNoErrA = totalNoErrA
        self.strainToBp = strainToBp
        self.strainToBpU = strainToBpU

    def __str__(self):
        buff = '%s Mbp (%s Mbp) # assembled (0 error)\n' \
               % (round((float(self.uniqueA[0][1]) / 1000000.), 3), round(float(self.totalCorrect) / 1000000., 3))

        #
        l = []
        for i in self.totalNoErrA:
            l.append(str(round((i[1] / float(self.uniqueA[0][1])) * 100., 2)) + '%')
            if i[1] == self.uniqueA[0][1]:
                break

        buff += ' '.join(l) + ' # % bp of total with error at coverage <= position\n'

        #
        buff += '%sbp (%s) # avgLen (contigs)\n' % (round(float(self.uniqueA[0][1]) / float(self.uniqueA[0][0]), 1), self.uniqueA[0][0])


        a = self.uniqueA

        comment = ('count', 'bp', 'cov', 'snp array')

        for j in range(2):
            for i in range(len(a))[1:]:

                buff += '%s%%\t' % round((a[i][j] / float(a[0][j])) * 100., 2)

            buff += '# %s\n' % comment[j]

        for i in range(len(a))[1:]:

            buff += '%s\t' % round((a[i][2] / float(a[i][1])), 2)

        buff += '# %s\n' % comment[2]
        buff += '%s # %s\n' % (self.snpErrorA, comment[3])

        for s, bp in self.strainToBp.iteritems():
            if s in self.strainToBpU:
                u = round(float(self.strainToBpU[s]) / 1000000., 3)
            else:
                u = 'NA'
            buff += '%s: %s Mbp (%s Mbp)\n' % (s, round(float(bp) / 1000000., 3), u)

        return buff


def sumUpStatInfo(statInfoList, maxStrains=3, maxCovT=10):

    uniqueA = np.zeros((maxStrains + 1, 3), dtype=np.int64)
    totalCorrect = 0

    snpErrorA = np.zeros(SNP_ERROR_ARRAY_LEN, dtype=np.int64)

    totalNoErrA = np.zeros((maxCovT, 2), dtype=np.int64)  # cov X (count, bp)

    strainToBp = {}
    strainToBpU = {}

    for si in statInfoList:

        if len(uniqueA) == len(si.uniqueA):
            uniqueA += si.uniqueA
        else:
            for i in range(len(si.uniqueA)):
                for j in range(3):
                    uniqueA[i][j] += si.uniqueA[i][j]

        totalCorrect += si.totalCorrect
        snpErrorA += si.snpErrorA

        totalNoErrA += si.totalNoErrA

        for s, bp in si.strainToBp.iteritems():
            if s not in strainToBp:
                strainToBp[s] = bp
            else:
                strainToBp[s] += bp

        for s, bp in si.strainToBpU.iteritems():
            if s not in strainToBpU:
                strainToBpU[s] = bp
            else:
                strainToBpU[s] += bp

    return StatInfo(uniqueA, totalCorrect, snpErrorA, totalNoErrA, strainToBp, strainToBpU)


def getStat(readRecFile, maxStrains=3, maxCovT=10):
    """
        @type readRecFile: str
        @type maxStrains: int
    """
    try:
        recList = hio.loadReadRec(readRecFile)

        # (total, strain 1, strain2 ...) x (count, bp, coverage)
        uniqueA = np.zeros((maxStrains + 1, 3), dtype=np.int64)

        # total (count, bp) of contigs with no error at the coverage >= index
        totalNoErrA = np.zeros((maxCovT, 2), dtype=np.int64)  # cov X (count, bp)

        snpErrorA = np.zeros(SNP_ERROR_ARRAY_LEN, dtype=np.int64)

        # cumulative length of contigs with no error
        totalCorrect = 0

        strainToBp = {}
        strainToBpU = {}  # strainId -> bp (unique for this strain)

        # for all super-reads
        for rec in recList:
            if rec.isSuperRead():
                recLen = rec.getLen()
                covCumul = rec.getAvgCov() * recLen

                labels = rec.getLabels()
                if labels is None:
                    print('No label for record: %s' % rec.recordId)
                if labels is not None:
                    bestLabel = labels[0]
                    bestLabels = [bestLabel]

                    # contig has no error
                    if bestLabel.error == 0:
                        totalCorrect += recLen
                        maxCov = 0
                    else:
                        # print '---'
                        # print bestLabel.error, bestLabel.codonError
                        maxCov = -1
                        err = 0
                        for i in range(1, len(bestLabel.totalErrorA[1])):
                            err += bestLabel.totalErrorA[1][i]
                            if err == bestLabel.error:
                                maxCov = i
                                break
                        assert maxCov >= 1
                        # print maxCov, bestLabel.totalErrorA[1]

                    for i in range(maxCov, maxCovT):
                        totalNoErrA[i][0] += 1
                        totalNoErrA[i][1] += recLen

                    # collect labels with the same minimal error
                    for label in labels[1:]:
                        if label.error == bestLabel.error and label.codonError == bestLabel.codonError:
                            bestLabels.append(label)
                        else:
                            break

                    # get the set of reference sequences and strain-ids
                    refSeqSet = set()
                    strainSet = set()
                    tmpSet = set()
                    for label in bestLabels:
                        refSeqSet.add(label.getRefDna())
                        strainSet.add(label.strainId)

                        if label.strainId not in tmpSet:
                            tmpSet.add(label.strainId)
                            if label.strainId not in strainToBp:
                                strainToBp[label.strainId] = recLen
                            else:
                                strainToBp[label.strainId] += recLen

                    if len(tmpSet) == 1:  # the contig maps only to one strain
                        s = bestLabels[0].strainId
                        if s not in strainToBpU:
                            strainToBpU[s] = recLen
                        else:
                            strainToBpU[s] += recLen

                    # there is more than one reference sequence
                    if len(refSeqSet) != 1:

                        if bestLabel.error == 0:
                            # there can be different sequences with 0 error, the difference is due to the N-characters
                            nFound = False
                            for r in refSeqSet:
                                if 'N' in r:
                                    nFound = True
                            if not nFound:
                                print('stat: 0 error but different references! %s' % rec.recordId)
                        else:
                            # print len(refSeqSet), bestLabel.error, bestLabel.codonError, recLen, len(strainSet), round(rec.getAvgCov(), 1)
                            # print rec.dnaSeq
                            # for r in refSeqSet:
                            #     print r

                            # the error is not 0, there are several reference sequences

                            # for each reference sequence, get the set of positions at which it differs from teh rec.dnaSeq
                            posSetList = []
                            for r in refSeqSet:
                                # diffPos = []
                                diffPosSet = set()
                                for i in range(len(rec.dnaSeq)):
                                    if r[i] != rec.dnaSeq[i] and qs.isDNAChar(r[i]):
                                        # diffPos.append(i)
                                        diffPosSet.add(i)
                                # print ','.join(map(lambda x: str(x), diffPos))
                                posSetList.append(diffPosSet)

                            # for each position, get how many times there was an error at the respective position
                            unionMap = {}
                            for s in posSetList:
                                for e in s:
                                    if e in unionMap:
                                        unionMap[e] += 1
                                    else:
                                        unionMap[e] = 1

                            # at how many positions there was an error that wasn't in all reference sequences
                            c = 0
                            for v in unionMap.values():
                                if v < len(refSeqSet):
                                    c += 1

                            snpError = c / len(refSeqSet)

                            snpErrorA[min(snpError - 1, len(snpErrorA) - 1)] += 1
                            # print snpErrorA

                            # print rec.getPosCovArray()
                            # print ''

                    # update total values
                    uniqueA[0][0] += 1
                    uniqueA[0][1] += recLen
                    uniqueA[0][2] += covCumul

                    # update values for a particular number of strains having the same error
                    strainSetLen = min(len(strainSet), maxStrains)
                    uniqueA[strainSetLen][0] += 1
                    uniqueA[strainSetLen][1] += recLen
                    uniqueA[strainSetLen][2] += covCumul

        si = StatInfo(uniqueA, totalCorrect, snpErrorA, totalNoErrA, strainToBp, strainToBpU)
        return si
    except Exception as e:
        print 'Exception: getStat(readRecFile=%s, maxStrains=%s)' % (readRecFile, maxStrains)
        print type(e)
        print e.message
        print e.args
        raise e


def test1():
    gf1 = 'acr_tran_3'
    gf2 = 'usher_1'
    # gf = 'aldedh_1'
    # gf = 'aa_permease_2'
    # gf = 'dde_tnp_is66_6'
    # gf = 'peptidase_s6_2'  # differences..
    # gf = 'molybdopterin_1'
    # gf = 'bpd_transp_1_3'
    # gf = 'hth_1_3'

    sample = 5

    readRecFile1 = os.path.join(comh.REFERENCE_DIR_ROOT,
                                '562/samples/%s/sample_partitioned/r_%s_join_read_rec.pkl.gz' % (sample, gf1))
    readRecFile2 = os.path.join(comh.REFERENCE_DIR_ROOT,
                                '562/samples/%s/sample_partitioned/r_%s_join_read_rec.pkl.gz' % (sample, gf2))

    print sumUpStatInfo([getStat(readRecFile1), getStat(readRecFile2)])


if __name__ == "__main__":
    # getStat(readRecFile='/net/metagenomics/projects/PPSmg/hsim/hsim01/562/samples/45/sample_partitioned/r_phage_min_tail_1_join_read_rec.pkl.gz', maxStrains=4)
    test1()
