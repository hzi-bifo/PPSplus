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

    Gene assembly evaluation and statistics.
"""
import os
import gzip

import numpy as np

from Bio.Alphabet import generic_dna
from Bio.Seq import Seq

from algbioi.com import fasta as fas
from algbioi.com import common as com
from algbioi.hsim import comh
from algbioi.haplo import hmain
from algbioi.haplo import hio
from algbioi.haplo import join_rec
from algbioi.com import qs


def assemblyStat(recSet):
    """
        Get the basic assembly statistics.

        @type recSet: set[read_rec.ReadRec]
    """
    print('Contigs total: %s' % len(recSet))

    print('Contig distribution: (#reads, count, contigLengths)')

    readsToCount = {}

    for e in recSet:
        if e.readRecList is not None:
            reads = len(e.readRecList)
        else:
            reads = 1
        if reads not in readsToCount:
            readsToCount[reads] = []
        readsToCount[reads].append(len(e.dnaSeq))
    t = []
    for k, v in readsToCount.iteritems():
        t.append((k, v))
    t.sort(key=lambda x: x[0], reverse=True)
    for k, v in t:
        v.sort(reverse=True)
        print('%s\t%s\t%s' % (k, len(v), v))


def getReadTrueMap(refGmap):
    """
        Read the true read-mapping.
        @param refGmap: sam mapping file, where the last column is the strain identifier, last but one gene annotation
        @type refGmap: str
        @return: map: readId -> (strainId, seqId, seqPos, geneAnnot)
        @rtype: dict[str,(str,str,int,str)]
    """
    mapping = {}
    for line in gzip.open(refGmap):
        line = line.strip()
        if line.startswith('#'):
            continue
        tokens = line.split('\t')
        readId = tokens[0]
        strainId = tokens[12]
        seqId = tokens[2]
        seqPos = int(tokens[3])
        geneAnnot = tokens[11]
        assert readId not in mapping, readId
        mapping[readId] = (strainId, seqId, seqPos, geneAnnot)
    return mapping


def getReadTrueMapSimple(refGmap):
    """
        Read the true read-mapping.
        @param refGmap: sam mapping file
        @type refGmap: str
        @return: map: readId -> (strainId, seqId)
        @rtype: dict[str,(str,str)]
    """
    mapping = {}
    for line in gzip.open(refGmap):
        line = line.strip()
        if line.startswith('#'):
            continue
        tokens = line.split('\t')
        readId = tokens[0]
        seqId = tokens[2]
        strainId = tokens[12]

        if readId not in mapping:
            mapping[readId] = (strainId, seqId)
        else:
            assert mapping[readId] == (strainId, seqId)
    return mapping


def assemblyPurity(recSet, refGmap):
    """
        For each assembled super-read, get from which strains its constituent reads come from and corresponding counts.
        @param recSet: set of read-records
        @param refGmap: sam mapping file
        @type recSet:  set[read_rec.ReadRec]
        @type refGmap: str
    """
    print('#reads, #strains, list[(sequenceId, #reads)]')
    # get true read mapping
    gmap = getReadTrueMap(refGmap)

    resultList = []
    for rec in recSet:
        if rec.readRecList is not None:
            strainSet = set()
            seqSet = set()
            seqCountDict = {}
            for r in rec.readRecList:
                strainSet.add(gmap[r.recordId][0])
                seqName = gmap[r.recordId][1]
                seqSet.add(seqName)

                if seqName not in seqCountDict:
                    seqCountDict[seqName] = 0
                seqCountDict[seqName] += 1

            t = []
            for k, v in seqCountDict.iteritems():
                t.append((k, v))
            t.sort(key=lambda x: x[1], reverse=True)

            resultList.append((len(rec.readRecList), len(strainSet), str(t),
                               ' | '.join((map(lambda x: ', '.join(x.getContent().split(',')[:4]), rec.labelEval)))))

    resultList.sort(key=lambda x: x[0], reverse=True)

    for e1, e2, e3, e4 in resultList:
        print('%s\t%s\t%s' % (e1, e2, e3))

    print('#reads, #strains, seq, strain, error, codonError, startPos, strain')
    resultList.sort(key=lambda x: x[0], reverse=True)
    for e1, e2, e3, e4 in resultList:
        print('%s\t%s\t%s' % (e1, e2, e4))


# def getContigLabel(rec, gmap):
#     """
#         Get the label of the contig, based on the labels of the constituent reads,
#         s.t. most of the reads have this label.
#
#         @param rec: read-record
#         @param gmap: map: readId -> (strainId, seqId, seqPos, geneAnnot)
#         @type rec: read_rec.ReadRec
#         @type gmap: dict[str,(str,str,int,str)]
#         @return: (strainId, seqId)
#         @rtype: (str, str)
#     """
#     the record consists only of itself
    # if rec.readRecList is None:
    #     return gmap[rec.recordId][:2]
    # else:
    #     countDict = {}  # map: (strainId, seqId) -> count
    #     for r in rec.readRecList:
    #         strainId, seqId = gmap[r.recordId][:2]
    #         if (strainId, seqId) not in countDict:
    #             countDict[(strainId, seqId)] = 0
    #         countDict[(strainId, seqId)] += 1
    #     tList = []
    #     for k, v in countDict.iteritems():
    #         tList.append((k, v))
    #     tList.sort(key=lambda x: x[1], reverse=True)
    #     return tList[0][0]
#


class LabelR(object):
    def __init__(self, error, codonError, strainId, seqId, startPos, refContig=None, strand=None, totalErrorA=None):
        self.error = error
        self.codonError = codonError
        self.strainId = strainId
        self.seqId = seqId
        self.startPos = startPos
        self.refContig = refContig
        self.strand = strand
        self.totalErrorA = totalErrorA

    def getRefDna(self):
        if self.strand is None:
            return None
        if self.strand == 1:
            return self.refContig
        else:
            assert self.strand == -1
            return str(Seq(self.refContig, generic_dna).reverse_complement())

    def getContent(self):
        return '(%s, %s, %s, %s, %s, %s, %s)' \
               % (self.seqId, self.strainId, self.error, self.codonError, self.startPos, self.strand, self.refContig)


def getRefContig(refSeq, startPos, seqLen):
    """
        Gets a substring of the reference sequence corresponding to the start position and length.
        If the start position is negative or the sequence is shorter than seqLen, it's padded by non-DNA characters (?).

        @param refSeq: reference dna sequence
        @param startPos: a start position within the sequence (can be negative)
        @param seqLen: the length of the substring sequence
        @type refSeq: str
        @type startPos: int
        @type seqLen: int
        @rtype: str
        @return: a substring sequence corresponding to the startPos and seqLen
    """
    if 0 <= startPos and (startPos + seqLen) <= len(refSeq):
        refContig = refSeq[startPos:startPos+seqLen]
    else:
        # the reference sequence doesn't match, needs to be padded from the left or right by non-DNA characters
        padLeft = ''
        padRight = ''
        if startPos + seqLen > len(refSeq):
            padRight = '?' * (startPos + seqLen - len(refSeq))
            endPos = len(refSeq)
        else:
            endPos = startPos + seqLen
        if startPos < 0:
            padLeft = '?' * (-startPos)
            startPos = 0
        refContig = padLeft + refSeq[startPos:endPos] + padRight

    assert seqLen == len(refContig)
    return refContig


def getError(dnaSeq, refSeq, covArray, maxCov, annotStart, annotLen, tripletMap):
    """
        Compute substitution errors between two DNA sequences.
        A comparison to a non-DNA character in the reference sequence is considered to be "correct" (except for '?').

        @param dnaSeq: a sequence
        @param refSeq: a reference sequence
        @param covArray: represents a coverage value for each position
        @param maxCov: maximum coverage considered (higher values will be mapped to this one
        @param annotStart: position within the dnaSeq at which the Hmm annotation starts
        @param annotLen: length of the annotation
        @param tripletMap: map: triplet -> codon
        @type dnaSeq: str
        @type refSeq: str
        @type covArray: ndarray
        @type maxCov: int
        @type annotStart: int
        @type tripletMap: dict[str,str]
        @rtype: (ndarray, int, int, int)
        @return: the result error array (1 row ~ total count, 2 row error count, 3 row codon-error, length ~ coverage),
            error, correct, codon-error
    """
    assert len(dnaSeq) == len(refSeq) == len(covArray)
    totalErrorA = np.zeros((3, maxCov + 1), dtype=np.uint64)
    # count the per base error for each position
    error = 0
    correct = 0
    codonError = 0
    for i in range(len(dnaSeq)):
        cov = min(int(round(covArray[i])), maxCov)
        totalErrorA[0][cov] += 1
        # count non-DNA in the reference as true
        if dnaSeq[i] == refSeq[i] or (refSeq[i] != '?' and (not qs.isDNAChar(refSeq[i]))):  # TODO: consider all the ambiguous characters ?
            correct += 1
        else:
            totalErrorA[1][cov] += 1
            error += 1
            codonErr = True
            if annotStart <= i < annotStart + annotLen:  # does it encode the same codon?
                offset = (i - annotStart) % 3
                start = i - offset
                codon1 = tripletMap.get(dnaSeq[start:start+3])
                codon2 = tripletMap.get(refSeq[start:start+3])
                if codon1 is not None and codon1 == codon2:
                    codonErr = False
            if codonErr:
                totalErrorA[2][cov] += 1
                codonError += 1

    return totalErrorA, error, correct, codonError


def getPerBaseErrorPkl(recListFile, refGmap, refSeqBuff, maxCov, translTable, storeLabels=True):
    try:
        recSet = set(hio.loadReadRec(recListFile))
        error = getPerBaseError(recSet, refGmap, refSeqBuff, maxCov, translTable, allLabels=storeLabels, tag=recListFile)
        if storeLabels:
            hio.storeReadRec(list(recSet), recListFile)
        return error

    except Exception as e:
        print "Exception in heval.getPerBaseErrorPkl"
        print recListFile, refGmap, 'buffer here', maxCov, translTable
        print e.message
        print type(e)
        return None


def getPerBaseError(recSet, refGmap, refSeqBuff, maxCov=20, translTable=11, allLabels=False, tag=''):
    """
        Get the per position, per coverage error.

        Skip singletons - not assembled read-records that consist only of itself.

        @param recSet: a set of records for which the per-base error is computed
        @param refGmap: reference file containing the true mapping for each read
        @param refSeqBuff: mapping for reference sequenceing map: (strainId, seqId) -> sequence
        @param maxCov: higher coverages will be mapped to this one
        @param allLabels: compute the error for all labels and set the labels (incl. errors) to the read-record

        @type recSet: set[read_rec.ReadRec]
        @type refGmap: str
        @type refSeqBuff: dict[(str,str),str]
        @type maxCov: int
        @type allLabels: bool

        @return: a two dimensional array, first row ~ total counts, second row ~ error counts (array length ~ coverage)
        @rtype: ndarray
    """
    # get the true mapping for each read (joined-pair-read)
    gmap = getReadTrueMap(refGmap)

    # define map: triplet -> codon
    tripletMap = join_rec.getTripletMap(translTable)

    # the result error array (1st row ~ total count, 2nd row error count, 3rd row codon error); array length ~ coverage
    errorA = np.zeros((3, maxCov + 1), dtype=np.uint64)
    errorSum = 0

    # second buffer for frequently retrieved entries
    refSeqBuff2 = {}

    # sum up the error for all read-records
    for rec in recSet:
        # skip not assembled reads
        if rec.readRecList is None:
            continue

        # examined settings (strainId, seqId, startPos, strain{1,-1})
        examinedSet = set()

        # count the error with respect to each read (also consider reverse complements)
        errList = []
        for recP in rec.readRecList:

            # start position within the reference sequence (position within the contig)
            startPos = gmap[recP.recordId][2] - 1 - recP.posWithinContig

            strainId, seqId = gmap[recP.recordId][:2]

            if (strainId, seqId, startPos, 1) not in examinedSet:

                # get the reference sequence
                refSeq = refSeqBuff2.get((strainId, seqId))
                if refSeq is None:
                    refSeq = refSeqBuff[(strainId, seqId)]
                    refSeqBuff2[(strainId, seqId)] = refSeq

                # get the reference contig
                refContig = getRefContig(refSeq, startPos, len(rec.dnaSeq))

                # compute error
                totalErrorA, error, correct, codonError = getError(rec.dnaSeq, refContig, rec.getPosCovArray(), maxCov,
                                                                   rec.annotStart, rec.annotLen, tripletMap)
                # store error
                errList.append((totalErrorA, error, codonError, correct, refContig, refSeq, seqId, strainId, startPos, 1))
                if error == 0 and not allLabels:
                    assert codonError == 0
                    break
                examinedSet.add((strainId, seqId, startPos, 1))
            else:
                refSeq = None

            # try reverse complement
            startPos = gmap[recP.recordId][2] - 1 - (len(rec.dnaSeq) - (len(recP.dnaSeq) + recP.posWithinContig))

            if (strainId, seqId, startPos, -1) not in examinedSet:

                if refSeq is None:
                    refSeq = refSeqBuff2.get((strainId, seqId))
                    if refSeq is None:
                        refSeq = refSeqBuff[(strainId, seqId)]
                        refSeqBuff2[(strainId, seqId)] = refSeq

                refContig = getRefContig(refSeq, startPos, len(rec.dnaSeq))

                totalErrorA, error, correct, codonError = getError(rec.dnaSeq,
                                                                   str(Seq(refContig, generic_dna).reverse_complement()),
                                                                   rec.getPosCovArray(), maxCov, rec.annotStart,
                                                                   rec.annotLen, tripletMap)

                errList.append((totalErrorA, error, codonError, correct, refContig, refSeq, seqId, strainId, startPos, -1))
                if error == 0 and not allLabels:
                    assert codonError == 0
                    break
                examinedSet.add((strainId, seqId, startPos, -1))

        if len(errList) == 0:
            print('Cannot determine error for: %s %s' % (rec.recordId, len(rec.dnaSeq)))
            continue

        # get the results from the smallest error (first sort according to the error, then codon-error)
        errList.sort(key=lambda x: (x[1], x[2], x[6], x[7], x[8], x[9]))

        totalErrorA, error, codonError, correct, refContig, refSeq, seqId, strainId, startPos, s = errList[0]

        # if the start position is outside of the reference sequence, try a different reference seq. with the same error
        # if startPos < 0:  # TODO: is it necessary?
        #     for err in errList[1:]:
        #         totalErrorA2, error2, codonError2, correct2, refContig2, refSeq2, seqId2, strainId2, startPos2, s2 = err
        #         if error2 > error or codonError2 > codonError:
        #             break
        #         elif startPos2 >= 0:
        #             totalErrorA = totalErrorA2
        #             break

        errorA += totalErrorA
        errorSum += error

        # store all the possible labels with corresponding errors
        if allLabels:
            examinedSet2 = set()  # (strainId, seqId)  # TODO: consider also the straind?
            labelList = []
            for err in errList:
                totalErrorA2, error2, codonError2, correct2, refContig2, refSeq2, seqId2, strainId2, startPos2, s2 = err
                if (strainId2, seqId2) not in examinedSet2:
                    labelList.append(LabelR(error2, codonError2, strainId2, seqId2, startPos2, refContig2, s2, totalErrorA2))
                    examinedSet2.add((strainId2, seqId2))

            # all possible labels of the read-record
            # if rec.labelEval is not None:
            #
            #     if len(rec.labelEval) != len(labelList):
            #         print("Diff label list len: %s (%s, %s) %s" % (rec.recordId, len(rec.labelEval), len(labelList)), tag)
            #     else:
            #         for lb1, lb2 in zip(rec.labelEval, labelList):
            #             if lb1.error != lb2.error or lb1.codonError != lb2.codonError:
            #                 print('Diff err %s (%s, %s) (%s, %s) %s' % (rec.recordId, lb1.error, lb2.error, lb1.codonError,
            #                                                    lb2.codonError, tag))
                        # else:
                        #     print '.',
                # pass

            rec.labelEval = labelList

        # if errList[0][1] > 0 and len(rec.readRecList) > 30:
        #     print 'Error "%s" at contig len: "%s" reads: "%s" matrix: %s' % (errList[0][1], len(rec.dnaSeq), len(rec.readRecList), errList[0][0][1])

    return errorA, errorSum  # TODO: remove errorSum !!!


def getAssemblyReport(errorAList, maxCov=20):
    """

        @type errorAList: list[ndarray]
        @return:
    """
    errorASum = np.zeros((3, maxCov + 1), dtype=np.uint64)
    eSum = 0
    for e, es in errorAList:
        if e is not None:
            errorASum += e
            eSum += es
        else:
            print('Warning: None in heval.getAssemblyReport')

    buff = 'err: %s\n' % eSum
    buff += str('Error per coverage total:\n%s\n\n' % errorASum)
    buff += str('@01:(' + ', '.join(map(lambda x: str(x), errorASum[0])) + ')\n')
    buff += str('@02:(' + ', '.join(map(lambda x: str(x), errorASum[1])) + ')\n')
    buff += str('@03:(' + ', '.join(map(lambda x: str(x), errorASum[2])) + ')\n\n')
    buff += str('Error per coverage %%\n')
    eList = []
    for i in range(1, maxCov + 1):
        if not com.isclose(0, float(errorASum[0][i])):
            v = (float(errorASum[1][i]) / float(errorASum[0][i])) * 100.
            buff += str('Cov: %s\tErr: %s%%\n' % (i, v))
            eList.append(str(v))
    buff += str('@04:(' + ', '.join(eList) + ')\n\n')
    buff += str('\nError-codon per coverage\n')
    eList = []
    for i in range(1, maxCov + 1):
        if not com.isclose(0, float(errorASum[0][i])):
            v = (float(errorASum[2][i]) / float(errorASum[0][i])) * 100.
            buff += str('Cov: %s\tErr: %s%%\n' % (i, v))
            eList.append(str(v))
    buff += str('@05:(' + ', '.join(eList) + ')\n\n')

    # compute cumulative error values:
    for i in range(maxCov, 1, -1):
        errorASum[0][i - 1] += errorASum[0][i]
        errorASum[1][i - 1] += errorASum[1][i]
        errorASum[2][i - 1] += errorASum[2][i]

    buff += str('Error per coverage cumulative total:\n%s\n\n' % errorASum)
    buff += str('@06:(' + ', '.join(map(lambda x: str(x), errorASum[0])) + ')\n')
    buff += str('@07:(' + ', '.join(map(lambda x: str(x), errorASum[1])) + ')\n')
    buff += str('@08:(' + ', '.join(map(lambda x: str(x), errorASum[2])) + ')\n\n')
    buff += str('Error per coverage cumulative %%\n')
    eList = []
    for i in range(1, maxCov + 1):
        if not com.isclose(0, float(errorASum[0][i])):
            v = (float(errorASum[1][i]) / float(errorASum[0][i])) * 100.
            buff += str('Cov: %s\tErr: %s%%\n' % (i, v))
            eList.append(str(v))
    buff += str('@09:(' + ', '.join(eList) + ')\n\n')
    buff += str('\nError-codon per coverage cumulative %%\n')
    eList = []
    for i in range(1, maxCov + 1):
        if not com.isclose(0, float(errorASum[0][i])):
            v = (float(errorASum[2][i]) / float(errorASum[0][i])) * 100.
            buff += str('Cov: %s\tErr: %s%%\n' % (i, v))
            eList.append(str(v))
    buff += str('@10:(' + ', '.join(eList) + ')\n\n')

    buff += str('@11:(' + ', '.join(map(lambda x: str(x), range(1, maxCov + 1))) + ')')

    return buff


def _test():

    # gf = 'pfl_3'
    # gf = 'aminotran_3_1'

    # gf = 'acr_tran_3'
    # gf = 'abc_tran_3'  # too many pieces..
    # gf = 'usher_1'
    gf = 'aldedh_1'  # TODO: explore two long contigs belonging to one reference !!!
    # gf = 'molybdopterin_1'
    # gf = 'aa_permease_2'
    # gf = 'hth_1_3'
    # gf = 'response_reg_1'
    # gf = 'lysr_substrate_2'
    # gf = 'xan_ur_permease_1'
    # gf = 'bpd_transp_1_3'
    # gf = 'feccd_3'
    # gf = 'terminase_gpa_3'

    baseDir = os.path.join(comh.REFERENCE_DIR_ROOT, '562/samples/0/sample_partitioned')
    recSet = hmain.buildSuperReads(inFq=os.path.join(baseDir, 'r_%s_join.fq.gz' % gf),
                                   inDomtblout=os.path.join(baseDir, 'r_%s_join_prot.domtblout.gz' % gf),
                                   inProtFna=os.path.join(baseDir, 'r_%s_join_prot.fna.gz' % gf))

    refGmap = os.path.join(baseDir, 'r_%s_join_gmap.sam.gz' % gf)

    refSeqBuff = fas.getSequenceBuffer([os.path.join(comh.REFERENCE_DIR_ROOT, '562', comh.FASTA_GENOMES_DIR_NAME),
                                        os.path.join(comh.REFERENCE_DIR_ROOT, '562', comh.FASTA_GENOMES_DRAFT_DIR_NAME)])

    print 'assembly finished'
    errA = getPerBaseError(recSet, refGmap, refSeqBuff, maxCov=10, allLabels=True)

    print getAssemblyReport([errA], maxCov=10)

    assemblyStat(recSet)
    assemblyPurity(recSet, refGmap)


if __name__ == "__main__":
    _test()

# gf = 'okr_dc_1_5'
    # gf = 'mfs_1_3'
    # gf = 'eal_3'
    # gf = 'sigma54_activat_1'
    #### gf = 'pfkb_2'
    # gf = 'sbp_bac_5_3'
    # gf = 'aminotran_3_1'
    # gf = 'phage_integrase_1'
    # gf = 'fhipep_6'
    # gf = 'hatpase_c_4'
    # gf = 'aldo_ket_red_1'
    # gf = 'hlyd_2'
    # gf = 'big_1_1'
    # gf = 'sdf_3'
    # gf = 'fe_adh_5'
    # gf = 'mfs_2_3'
    # gf = 'na_sulph_symp_1'
    # gf = 'ggdef_5'
    # gf = 'papd_n_3'
    # gf = 'phosphorylase_1'  # all contigs to the same label..
    # gf = 'papc_n_3'
    # gf = 'f_bp_aldolase_2'  # high cov. error
    # gf = 'ald_xan_dh_c2_3'
    # gf = 'gntp_permease_2'  # two contigs with the same label
    # gf = 'dde_tnp_is66_6'
    # gf = 'bpd_transp_2_3'
    # gf = 'sugar_tr_1'
    # gf = 'deorc_3'
    # gf = 'pep_utilizers_c_5'  # three contigs same label
    #gf = 'sbp_bac_3_3'
    # gf = 'asma_2'

    # gf = 'dnag_1'
    # gf = 'infc_1'
    # gf = 'pgk_0'
    # gf = 'tsf_0'
    # gf = 'frr_0'
    # gf = 'nusa_1'
    # gf = 'pyrg_1'
    # gf = 'rpma_1'
    # gf = 'smpb_1'
    # gf = 'rpsc_1'
    # gf = 'rpsi_1'
    # gf = 'rpsk_1'
    # gf = 'rpss_1'
    # gf = 'rpsb_1'
    # gf = 'rpse_1'
    # gf = 'rpsj_1'
    # gf = 'rpsm_1'
    # gf = 'rpla_1'
    # gf = 'rplb_1'
    # gf = 'rplc_1'
    # gf = 'rpld_1'
    # gf = 'rple_1'
    # gf = 'rplf_1'  # two variants !
    # gf = 'rplk_1'
    # gf = 'rpll_1'
    # gf = 'rplm_1'
    # gf = 'rpln_1'
    # gf = 'rplp_1'
    # gf = 'rpls_1'
    # gf = 'rplt_1'

    # now skipping several records

    # gf = 'big_3_4_1'  # all just from one?

    # 7 mlst
    # gf = 'adk_3'  # too small
    # gf = 'lyase_1_1'  # ok
    # gf = 'ter_1'  # small, ok
    # gf = 'fumarasec_c_3'  # small
    # gf = 'iso_dh_1'  # too small
    # gf = 'adenylsucc_synt_1'  # only once..
    # gf = 'hatpase_c_4'  # mixture of too many
    # gf = 'dna_gyraseb_4'
    # gf = 'toprim_1'  # small
    # gf = 'dna_gyraseb_c_5' # small but ok
    # gf = 'reca_2'  # too small..
    # gf = 'ldh_1_n_2' # only one strain..
    # gf = 'ldh_1_c_2'  # only one strain..

    # amphora
    # gf = 'dnag_dnab_bind_4'
    # gf = 'toprim_n_2'
    # gf = 'toprim_4_1'
    # gf = 'zf_chc2_4'
    # gf = 'dnab_bind_2'
    # gf = 'pgk_3'
    # gf = 'rna_pol_rpb2_6_4'
    # gf = 'ef_ts_4'
    # gf = 'rrf_3'
    # gf = 's1_1'
    # gf = 'nusa_n_3'
    # gf = 'gatase_3'
    # gf = 'ctp_synth_n_4'
    # gf = 'ribosomal_l27_2'
    # gf = 'ribosomal_l27_2'
    # gf = 'smpb_2'
    # gf = 'kh_2_2'
    # gf = 'ribosomal_s9_2'
    # gf = 'ribosomal_s2_2'
    # gf = 'ribosomal_s2_2'
    # gf = 'ribosomal_l1_2'
    # gf = 'ribosomal_l2_c_3'
    # gf = 'ribosomal_l4_2'
    # gf = 'ribosomal_l3_2'
    # gf = 'ribosomal_l4_2'
    # gf = 'pdt_3'
    # gf = 'ribosomal_s8_2'
    # gf = 'ribosomal_s9_2'
    # gf = 'ribosomal_l13_2'
    # gf = 'ribosomal_l14_2'
    # gf = 'ribosomal_l16_2'
    # gf = 'ribosomal_l19_2'
    # gf = 'ribosomal_l20_2'
