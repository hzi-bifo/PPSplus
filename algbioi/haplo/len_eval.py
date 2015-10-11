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

    Length evaluation.
"""

import os
# import sys
import numpy as np
import gzip

from algbioi.com import csv
from algbioi.com import fasta as fas
from algbioi.com import parallel
from algbioi.haplo import hio
from algbioi.haplo import heval
from algbioi.hsim import comh

from Bio.Alphabet import generic_dna
from Bio.Seq import Seq

UINT8_MAX = np.iinfo(np.uint8).max


def hmmAnnotContigsSnowball(samplePartDir, tmpDir):
    """

    """
    hmmProfileDir = os.path.join(comh.HMM_PROFILE_DIR, comh.HMM_PROFILE_FILES_PARTITIONED_DIR)
    if not os.path.isdir(tmpDir):
        os.mkdir(tmpDir)

    taskList = []
    for f in os.listdir(samplePartDir):
        fPath = os.path.join(samplePartDir, f)
        if f.endswith('join_read_rec.pkl.gz') and os.path.isfile(fPath):

            hmmProfile = os.path.join(hmmProfileDir, f[2:-21] + '.hmm')
            assert os.path.isfile(hmmProfile)

            # annotLen(fPath, tmpDir, hmmProfile)
            taskList.append(parallel.TaskThread(annotLen, (fPath, tmpDir, hmmProfile)))

    parallel.runThreadParallel(taskList, comh.MAX_PROC, keepRetValues=False)


def hmmAnnotContigsSat(contPklDir, tmpDir):
    """

    """
    hmmProfileDir = os.path.join(comh.HMM_PROFILE_DIR, comh.HMM_PROFILE_FILES_PARTITIONED_DIR)
    if not os.path.isdir(tmpDir):
        os.mkdir(tmpDir)

    # buffer pfNames, map: pfName -> hmmProfileFile
    pfNameToProfile = {}
    for f in os.listdir(hmmProfileDir):
        fPath = os.path.join(hmmProfileDir, f)
        for line in open(fPath):
            line = line.strip()
            if line.startswith('ACC'):
                pfName = line.split()[1].split('.')[0]
                pfNameToProfile[pfName] = fPath
                break

    # for each read record
    taskList = []
    for f in os.listdir(contPklDir):
        fPath = os.path.join(contPklDir, f)
        if f.endswith('.pkl.gz') and os.path.isfile(fPath):
            pfName = f[8:-7]
            assert pfName in pfNameToProfile, pfName

            taskList.append(parallel.TaskThread(annotLen, (fPath, tmpDir, pfNameToProfile[pfName])))

    parallel.runThreadParallel(taskList, comh.MAX_PROC, keepRetValues=False)


def annotLen(recFile, tmpDir, hmmProfileFile):
    """
        Get the length of the covered model for each sequence.
        Record attribute (evalHmmCoord) contains tuple (hmmCoordLen, modelLen, hmmCoordStart, tag)

        @param recFile: record file
        @param tmpDir: temporary directory
        @param hmmProfileFile: corresponding hmm model
        @type recFile: str
        @type tmpDir: str
        @type hmmProfileFile: str
    """
    try:
        # get the DNA model length
        modelLen = None
        for line in open(hmmProfileFile):
            line = line.strip()
            if line.startswith('LENG'):
                modelLen = int(line.split()[1])
                break
        assert modelLen is not None

        # write record sequences into a fasta file
        recSet = hio.loadObj(recFile)
        counter = 0
        idToRec = {}
        contigFile = os.path.join(tmpDir, os.path.basename(hmmProfileFile)[:-3] + 'fna')
        contigProtFile = contigFile[:-4] + '_prot.fna'
        assert not os.path.isfile(contigFile) and not os.path.isfile(contigProtFile)
        out = csv.OutFileBuffer(contigFile)
        for rec in recSet:
            if not hasattr(rec, 'evalHmmCoord'):  # add the attribute if it does not exist !!!
                rec.evalHmmCoord = None
            assert counter not in idToRec
            idToRec[counter] = rec
            out.writeText('>%s\n%s\n' % (counter, rec.dnaSeq))
            counter += 1
        out.close()

        # translate DNA into PROT sequences
        fas.fastaToProt(contigFile, contigProtFile, comh.TRANSLATION_TABLE)

        # run HMM search on the PROT fasta file
        domOut = contigProtFile[:-3] + 'dom'
        cmd = '%s/hmmsearch -o /dev/null --noali --domtblout %s -E 0.01 --cpu 5 %s %s' % (comh.HMMER_BINARY, domOut,
                                                                                          hmmProfileFile,
                                                                                          contigProtFile)

        parallel.reportFailedCmd(parallel.runCmdSerial([parallel.TaskCmd(cmd, os.path.dirname(contigProtFile))]))

        # parse the DOM file, store covered model length for each sequence
        for line in open(domOut):
            line = line.strip()
            if not line.startswith('#') and len(line) > 0:
                hit = line.split()

                counter, tag = map(lambda x: int(x), hit[0].split('_'))
                rec = idToRec[counter]

                # alignment env coord
                protStart = int(hit[19]) - 1
                protLen = int(hit[20]) - protStart

                # alignment coord
                protStartAli = int(hit[17]) - 1
                protLenAli = int(hit[18]) - protStartAli

                # hmm coordinates
                hmmCoordStart = int(hit[15]) - 1
                hmmCoordLen = int(hit[16]) - hmmCoordStart

                # the env coordinates start (and end) often before the alignment coordinates and end after, get the offsets
                offsetEnv = protStartAli - protStart
                offsetEnvE = protLen - offsetEnv - protLenAli

                assert offsetEnv >= 0
                assert offsetEnvE >= 0

                hmmCoordStart -= offsetEnv
                hmmCoordLen += offsetEnv + offsetEnvE

                # update the record
                if rec.evalHmmCoord is None or rec.evalHmmCoord[0] < hmmCoordLen:
                    rec.evalHmmCoord = (hmmCoordLen, modelLen, hmmCoordStart, tag)

        # clean up
        os.remove(contigFile)
        os.remove(contigProtFile)
        parallel.reportFailedCmd(parallel.runCmdSerial([parallel.TaskCmd('gzip -f %s' % domOut, tmpDir)]))

        # for c, rec in idToRec.iteritems():
        #     print c, rec.evalHmmCoord

        # store to a file
        hio.storeObj(recSet, recFile)
    except Exception as e:
        print recFile, tmpDir, hmmProfileFile
        print e.message
        print e.args
        raise e


def contLenEval(srcDir):
    """

    """
    errA = np.zeros((2, 100), dtype=np.int64)

    for f in os.listdir(srcDir):
        fPath = os.path.join(srcDir, f)

        if os.path.isfile(fPath) \
                and (f.endswith('join_read_rec.pkl.gz') or (f.startswith('contigs_PF') and f.endswith('.pkl.gz'))):
            recSet = hio.loadObj(fPath)
            for rec in recSet:
                if hasattr(rec, 'evalHmmCoord') and rec.evalHmmCoord is not None:
                    hmmCoordLen, modelLen = rec.evalHmmCoord[:2]
                    idx = max(1, min(int((hmmCoordLen / float(modelLen)) * 100.), 100)) - 1
                    errA[0][idx] += len(rec.dnaSeq)
                    errA[1][idx] += 1

    for i in range(98, -1, -1):
        errA[0][i] += errA[0][i+1]
        errA[1][i] += errA[1][i+1]

    return errA


def contLenEvalReport(errAList, tag):
    """

    """
    errA = np.zeros((2, 100), dtype=np.int64)
    for e in errAList:
        errA += e

    buff = ''
    buff += '@%s\tb\tc(%s)\n' % (tag, ', '.join(map(lambda x: str(x), errA[0])))
    buff += '@%s\tc\tc(%s)\n' % (tag, ', '.join(map(lambda x: str(x), errA[1])))
    buff += '@%s\tb50\tc(%s)\n' % (tag, ', '.join(map(lambda x: str(x), errA[0][50:])))
    buff += '@%s\tc50\tc(%s)\n' % (tag, ', '.join(map(lambda x: str(x), errA[1][50:])))
    return buff


class RefRec(object):
    def __init__(self, strainId, seqId, seq):
        self.strainId = strainId
        self.seqId = seqId
        self.seq = seq
        self.covA = np.zeros((100, len(seq)), dtype=np.uint8)
        # self._cumul = False

    def addSeq(self, idx, recSeq, lab, checkErrFree=True):
        """
            @param checkErrFree: can be set to False for faster processing
            @type idx: int
            @type recSeq: str
            @type lab: heval.LabelR
            @type checkErrFree: bool
        """
        # extra checking for error free contig-to-reference mapping
        if checkErrFree:
            if lab.error == 0:
                if lab.strand == -1:
                    recSeq = str(Seq(recSeq, generic_dna).reverse_complement())

                refSeq = self.seq[lab.startPos: lab.startPos + len(recSeq)]
                if not (recSeq == refSeq or 'N' in refSeq or 'Y' in refSeq):
                    print('Different sequences:\n%s\n%s\n' % (recSeq, refSeq))

        if lab.startPos < 0:
            print('Label start position is not positive (%s), sequence considered from idx=0' % lab.startPos)

        for i in range(max(0, lab.startPos), min(lab.startPos + len(recSeq), len(self.seq))):
            for j in range(idx + 1):
                try:
                    self.covA[j][i] = min(self.covA[j][i] + 1, UINT8_MAX)
                except Exception as e:
                    print j, i
                    raise e

    def updateCumul(self, cumulA):
        """
            @type cumulA: ndarray
            @return sequence length
            @rtype int
        """
        for i in range(len(self.seq)):
            for j in range(100):
                if self.covA[j][i] > 0:
                    cumulA[j] += 1

        return len(self.seq)


def getRefCover(srcContigDir, refDirList, tag, samDir):

    try:
        # src = '/Users/ivan/Documents/nobackup/hsim01/562/samples10/0/tmp_len/tmp_smaz_hned_my.pkl.gz'  # test

        # get strains contained in the dataset
        strainSet = set()
        for f in os.listdir(samDir):
            fPath = os.path.join(samDir, f)
            if f.endswith('sam.gz') and os.path.isfile(fPath):
                for line in gzip.open(fPath):
                    line = line.strip()
                    strainSet.add(line.rsplit('\t', 1)[-1])

        # strainSet = set(['NZ_ANLR00000000', 'NZ_AKLX00000000', 'NZ_AKLB00000000'])  # test

        # buffer sequences
        strainToFile = {}
        for strain in strainSet:
            for refDir in refDirList:
                fPath = os.path.join(refDir, '%s.fna' % strain)
                if os.path.isfile(fPath):
                    strainToFile[strain] = fPath
                    break
        assert len(strainSet) == len(strainToFile)

        # map: (strainId, seqId) -> seqRec
        idToRefRec = {}
        for strainId, fPath in strainToFile.iteritems():

            for seqId, seq in fas.fastaFileToDictWholeNames(fPath).iteritems():
                rid = (strainId, seqId)
                assert rid not in idToRefRec
                idToRefRec[rid] = RefRec(strainId, seqId, seq)

        # for each gene domain
        for f in os.listdir(srcContigDir):
            fPath = os.path.join(srcContigDir, f)
            if os.path.isfile(fPath) \
                    and (f.endswith('join_read_rec.pkl.gz') or (f.startswith('contigs_PF') and f.endswith('.pkl.gz'))):

                # for each contig update the reference coverage
                for rec in hio.loadObj(fPath):
                    # if rec.labelEval is not None:
                    #     print hasattr(rec, 'evalHmmCoord'), hasattr(rec, 'labelEval'), rec.evalHmmCoord, rec.labelEval

                    # if the relevant attributes are defined, process the contig
                    if hasattr(rec, 'evalHmmCoord') and hasattr(rec, 'labelEval') and rec.evalHmmCoord is not None \
                            and rec.labelEval is not None:

                        hmmCoordLen, modelLen = rec.evalHmmCoord[:2]
                        idx = max(1, min(int((hmmCoordLen / float(modelLen)) * 100.), 100)) - 1

                        error = None
                        codonError = None
                        for lab in rec.labelEval:
                            if error is None:
                                error = lab.error
                                codonError = lab.codonError
                            elif error != lab.error or codonError != lab.codonError:
                                # process label !!!
                                # print lab.getContent()
                                # idToRefRec[(lab.strainId, lab.seqId)].addSeq(idx, rec.dnaSeq, lab)
                                break
                            try:
                                idToRefRec[(lab.strainId, lab.seqId)].addSeq(idx, rec.dnaSeq, lab)
                            except Exception as e:
                                print e.message
                                print e.args
                                print lab.strainId, lab.seqId, strainSet, '\n--------'

        # collect overall reference coverage
        cumulA = np.zeros(100, dtype=np.uint64)
        strainBp = 0
        # for each position, was it covered by a contig with length at least XXX
        for r in idToRefRec.values():
            strainBp += r.updateCumul(cumulA)

        return (cumulA, strainBp, tag)

    except Exception as e:
        print srcContigDir, refDirList, tag
        print e.message
        print e.args
        raise e


def getRefCoverReport(listOfCumulABp, tag):
    """
        @type listOfCumulABp: list[(ndarray,int)]
    """
    cumulA = np.zeros(100, dtype=np.uint64)
    strainBp = 0
    for cumulAPart, strainBpPart in listOfCumulABp:
        cumulA += cumulAPart
        strainBp += strainBpPart

    buff = '@0%s\n' % tag
    buff += '@1\tRA\t%s\n' % (', '.join(map(lambda x: str(x), cumulA)))
    buff += '@2\tRBP\t%s\n' % strainBp
    buff += '@3\tPA\t%s\n\n' % (', '.join(map(lambda x: str(round(((float(x) / float(strainBp)) * 100.), 2)), cumulA)))

    return buff


def _test1():
    samplePart = '/home/igregor/Documents/work/hsim/562/samples/0/sample_partitioned'
    tmpDir = '/home/igregor/Documents/work/hsim/562/samples/0/tmp_len'
    hmmAnnotContigsSnowball(samplePart, tmpDir)


def _test2():
    contPklDir = '/home/igregor/Documents/work/hsim/562/samples/0/sat/cont_pkl'
    tmpDir = '/home/igregor/Documents/work/hsim/562/samples/0/sat/tmp_len'
    hmmAnnotContigsSat(contPklDir, tmpDir)


def _test3():
    srcDir = '/home/igregor/Documents/work/hsim/562/samples/0/sample_partitioned'
    # srcDir = '/home/igregor/Documents/work/hsim/562/samples/0/sat/cont_pkl'
    errA = contLenEval(srcDir)
    print contLenEvalReport([errA], 'tag')


def _test4():
    # srcDir = '/Users/ivan/Documents/nobackup/hsim01/562/samples10/0/sample_partitioned'
    srcDir = '/Users/ivan/Documents/nobackup/hsim01/562/samples10/0/sat/cont_pkl'
    refDirList = [os.path.join(comh.REFERENCE_DIR_ROOT, '562', comh.FASTA_GENOMES_DIR_NAME),
                  os.path.join(comh.REFERENCE_DIR_ROOT, '562', comh.FASTA_GENOMES_DRAFT_DIR_NAME)]

    getRefCover(srcDir, refDirList)

if __name__ == "__main__":
    # _test4()

    # print getRefCover('/net/metagenomics/projects/PPSmg/hsim/hsim01/562/samples/5/sample_partitioned',
    #             ['/net/metagenomics/projects/PPSmg/hsim/hsim01/562/fasta_genomes',
    #              '/net/metagenomics/projects/PPSmg/hsim/hsim01/562/fasta_draft_genomes'],
    #             '5_snow', '/net/metagenomics/projects/PPSmg/hsim/hsim01/562/samples/5/sample_partitioned')
    # annotLen('/net/metagenomics/projects/PPSmg/hsim/hsim01/562/samples/8/sample_partitioned/r_eii_sor_5_join_read_rec.pkl.gz',
    #          '/net/metagenomics/projects/PPSmg/hsim/hsim01/562/samples/8/tmp_len',
    #          '/net/metagenomics/projects/PPSmg/database/pfamV27/nobackup/pub/databases/Pfam/releases/Pfam27.0/pfam_a_and_amphora2/eii_sor_5.hmm')
    pass




