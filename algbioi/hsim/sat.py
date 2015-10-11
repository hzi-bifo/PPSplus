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

    Wrapping scripts to run the SAT assembler and post-process its output.
"""

# ./SAT-Assembler.sh -m <HMM file> -f <fasta file> [options] -o:  output file name

import os
import gzip
import numpy as np

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

from algbioi.com import csv
from algbioi.com import fq
from algbioi.com import fasta as fas
from algbioi.com import parallel
from algbioi.hsim import comh
from algbioi.haplo import hio
from algbioi.haplo import heval


class CRec(object):
    def __init__(self, recordId, dnaSeq, posWithinContig=0):
        """
            SAT Contig record.

            @type recordId: str
            @type dnaSeq: str
            @type posWithinContig: int
        """
        self.recordId = recordId
        self.readRecList = None
        self.posWithinContig = posWithinContig
        self.dnaSeq = dnaSeq

        self.annotStart = 0  # TODO init currently doesn't consider
        self.annotLen = 0  # TODO init

        self.labelEval = None

        # hmm coordinates for evaluation purposes
        self.evalHmmCoord = None

        self._posCovArray = None

    def addRead(self, rRec):
        """
            Adds a read that well aligns to the contig.
            @type rRec: CRec
        """
        if self.readRecList is None:
            self.readRecList = []
        self.readRecList.append(rRec)

    def getPosCovArray(self):
        """
            Returns a position covarege array.
            @rtype: ndarray
        """
        if self._posCovArray is None:
            self._posCovArray = np.zeros(len(self.dnaSeq), dtype=np.float64)
            if self.readRecList is not None:
                for rec in self.readRecList:
                    for i in range(len(rec.dnaSeq)):
                        idx = i + rec.posWithinContig
                        if 0 <= idx < len(self.dnaSeq):
                            self._posCovArray[idx] += 1
                        elif idx >= 0:
                            break

        return self._posCovArray


def fqToFasta(fqPairList, outFile, checkUniqueNames=True):
    """
        Creates an input FASTA file for the SAT assembler.

        @param fqPairList: list of pairs of files containing pair end reads (fq1, fq2)
        @type fqPairList: list[(str,str)]
        @param outFile: output FASTA file path
    """
    out = csv.OutFileBuffer(outFile)
    nameSet = set()
    count = 1
    for fq1File, fq2File in fqPairList:
        for fq1, fq2 in zip(list(fq.ReadFqGen(fq1File)), list(fq.ReadFqGen(fq2File))):

            for fqEntry, pairEnd in zip([fq1, fq2], ['1', '2']):

                name, dna = fqEntry[:2]
                name = name[1:]
                assert name.endswith(pairEnd)
                if checkUniqueNames:
                    assert name not in nameSet
                    nameSet.add(name)
                name = '%s.%s.%s' % (name, count, pairEnd)
                out.writeText('>%s\n%s\n' % (name, dna))

            count += 1
    out.close()


def getSatRunTask(inFastaFile, outDir, hmmProfiles, satInstallDir, hmmerBinDir, options=''):
    """
        Get the task to run SAT.
        @param outDir: this directory will be created, must not exist, will contain all SAT output
    """
    if os.path.isdir(outDir):
        print('Directory "%s" already exists, remove the directory before running SAT!' % outDir)

    cmd = 'export PATH=%s:$PATH;%s -m %s -f %s -o %s %s' % (hmmerBinDir,
                                                            os.path.join(satInstallDir, 'SAT-Assembler.sh'),
                                                            hmmProfiles, inFastaFile, outDir, options)
    cwd = satInstallDir
    task = parallel.TaskCmd(cmd, cwd)  # parallel.reportFailedCmd(parallel.runCmdSerial([task]))
    return task


def getBwaIdxTask(fastaFile, aliIndexDir, idxName):
    """
        Get the task to build the BWA index.
    """
    cmd = '%s -a -q --seed 0 %s %s' \
          % (os.path.join(comh.BWA_INSTALL_DIR, 'bowtie2-build'), fastaFile, idxName)

    return parallel.TaskCmd(cmd, cwd=aliIndexDir)


def getBwaRunTask(aliIndexDir, idxName, fq1, fq2, samFile):  # --very-sensitive
    """
        Get the task to run BWA
    """
    cmd = '%s --all --no-unal --quiet --reorder --seed 0 --no-head --no-sq -x %s ' \
          '-1 %s -2 %s -S %s; gzip -f %s' \
          % (os.path.join(comh.BWA_INSTALL_DIR, 'bowtie2'), idxName, fq1, fq2, samFile, samFile)

    return parallel.TaskCmd(cmd, cwd=aliIndexDir)


def satMapReads(satOutDir, fqPairList, fqPartitionDir, aliIndexDir, samDir, maxProc=comh.MAX_PROC):
    """
        Map the reads onto the assembled contigs.

        @param satOutDir: output sat directory
        @param fqPairList: list of (fq1, fq2)
        @param fqPartitionDir: directory containing reads partitioned in pfam families
        @param aliIndexDir: directory containing alignment indices
        @param samDir: output directory containing SAM files

        @type satOutDir: str
        @type fqPairList list[(str,str)]
        @type fqPartitionDir: str
        @type aliIndexDir: str
        @type samDir: str
    """
    assert os.path.isdir(satOutDir)

    # create directories
    for d in [fqPartitionDir, aliIndexDir, samDir]:
        if not os.path.isdir(d):
            os.mkdir(d)

    # collect Pfam families for which contigs were assembled
    pfList = []
    for f in os.listdir(satOutDir):
        if f.endswith('contigs.fa'):
            fPath = os.path.join(satOutDir, f)
            if os.path.getsize(fPath) > 0:
                pfName = f.split('_', 1)[0]
                pfList.append(pfName)

    if True:
        # read-in, buffer all reads for a sample
        readNameToPairDict = {}
        for fq1File, fq2File in fqPairList:
            for fq1, fq2 in zip(list(fq.ReadFqGen(fq1File)), list(fq.ReadFqGen(fq2File))):
                name1 = fq1[0][1:-2]
                name2 = fq2[0][1:-2]
                assert name1 == name2
                assert name1 not in readNameToPairDict
                readNameToPairDict[name1] = (fq1, fq2)

        # for each Pfam for which contigs were assembled, get reads from which they were assembled
        for pfName in pfList:
            hmmFile = os.path.join(satOutDir, '_%s.hmmer' % pfName)
            if not os.path.isfile(hmmFile):
                hmmFile = os.path.join(satOutDir, '0_sat.fna_%s.hmmer' % pfName)  # TODO derive the name in a better way

            # collect read names
            readNameSet = set()
            for line in open(hmmFile):
                line = line.strip()
                readNameSet.add(line.split('/', 1)[0])

            # write reads belonging to this contig to fq files
            if len(readNameSet) > 0:
                fw1 = fq.WriteFq(os.path.join(fqPartitionDir, '%s_1.fq.gz' % pfName))
                fw2 = fq.WriteFq(os.path.join(fqPartitionDir, '%s_2.fq.gz' % pfName))

                for name in readNameSet:
                    fq1, fq2 = readNameToPairDict[name]
                    name1, dna1, p, qs1 = fq1
                    name2, dna2, p, qs2 = fq2
                    fw1.writeFqEntry(name1, dna1, qs1)
                    fw2.writeFqEntry(name2, dna2, qs2)

                fw1.close()
                fw2.close()

    if True:
        # building the BWA index
        taskList = []
        for pfName in pfList:
            taskList.append(getBwaIdxTask(os.path.join(satOutDir, '%s_contigs.fa' % pfName), aliIndexDir, pfName))

        parallel.reportFailedCmd(parallel.runCmdParallel(taskList, maxProc))

    if True:
        # running BWA
        taskList = []
        for pfName in pfList:
            fq1 = os.path.join(fqPartitionDir, '%s_1.fq.gz' % pfName)
            fq2 = os.path.join(fqPartitionDir, '%s_2.fq.gz' % pfName)
            taskList.append(getBwaRunTask(aliIndexDir, pfName, fq1, fq2, os.path.join(samDir, '%s.sam' % pfName)))

        parallel.reportFailedCmd(parallel.runCmdParallel(taskList, maxProc))

    # removing the BWA index
    if True:
        for f in os.listdir(aliIndexDir):
            if 'bt' in f:
                os.remove(os.path.join(aliIndexDir, f))

    # removing FQ files
    # if False:  # it is used later !!!
    #     for f in os.listdir(fqPartitionDir):
    #         if f.endswith('fq.gz'):
    #             os.remove(os.path.join(fqPartitionDir, f))


def diff(str1, str2):
    """
        @return: Number of differences in two strings.
        @rtype: int
    """
    str1Len = len(str1)
    str2Len = len(str2)

    if str1Len < str2Len:
        str1 += '?' * (str2Len - str1Len)
    elif str1Len > str2Len:
        str2 += '?' * (str1Len - str2Len)
    assert len(str1) == len(str2)

    c = 0
    for i, j in zip(list(str1), list(str2)):
        if i != j and i != '?' and j != '?':
            c += 1
    return c


def satToContigRecords(samDir, fqPartitionDir, satOutDir, dstDir=None, maxDiff=10):
    """
        @type samDir: str
        @type fqPartitionDir: str
        @type satOutDir: str
        @type dstDir: str
        @type maxDiff: int
    """
    try:
        over = 0
        overDiffSize = 0
        total = 0
        pfToRecSet = {}

        # for each sam file (i.e. pfam)
        for f in os.listdir(samDir):
            samPath = os.path.join(samDir, f)
            if f.endswith('.sam.gz') and os.path.isfile(samPath) and os.path.getsize(samPath) > 0:
                pfName = f.split('.')[0]

                # buffer reads of the sam
                readNameToPair = {}  # map: readName -> ((dna1, qs1), (dna2, qs2))
                fq1File = os.path.join(fqPartitionDir, '%s_1.fq.gz' % pfName)
                fq2File = os.path.join(fqPartitionDir, '%s_2.fq.gz' % pfName)
                for fq1, fq2 in zip(list(fq.ReadFqGen(fq1File)), list(fq.ReadFqGen(fq2File))):
                    name1, dna1, p, qs1 = fq1
                    name2, dna2, p, qs2 = fq2
                    name = name1[1:-2]
                    assert name1.endswith('/1')
                    assert name2.endswith('/2')
                    assert name1[:-2] == name2[:-2]
                    assert name not in readNameToPair
                    readNameToPair[name] = (dna1, qs1, dna2, qs2)

                # buffer contigs
                contigNameToSeq = fas.fastaFileToDictWholeNames(os.path.join(satOutDir, '%s_contigs.fa' % pfName))

                cNameToRec = {}  # map: cName -> CRec

                # read in a SAM file
                for line in gzip.open(samPath):
                    total += 1
                    line = line.strip()
                    tokens = line.strip().split('\t')
                    readName, flag, contigName, pos, m, cigar, rnext, pnext, tlen, seq, qs = tokens[:11]
                    cName = '%s_%s' % (pfName, contigName)
                    flag = int(flag)
                    pos = int(pos)

                    # parse flags
                    if flag & 0x4 == 0:  # segment mapped

                        contigSeq = contigNameToSeq[contigName]
                        # contigSeqRev = str(Seq(contigSeq, generic_dna).reverse_complement())

                        if flag & 0x10 == 0:  # reverse complement
                            rc = False
                        else:
                            rc = True
                            # revSeq = str(Seq(seq, generic_dna).reverse_complement())
                            # qsRev = qs[::-1]

                        if flag & 0x40 != 0:  # first segment
                            firstSeg = True
                        else:
                            assert flag & 0x80 != 0  # last segment
                            firstSeg = False

                        # get a read
                        dna1, qs1, dna2, qs2 = readNameToPair[readName]
                        assert len(dna1) == len(dna2) == len(seq)

                        # corresponding contig sequence
                        cSeq = contigSeq[pos - 1:len(seq) + pos - 1]

                        if firstSeg:
                            nameSuffix = '\\1'
                            if rc:
                                dna = str(Seq(dna1, generic_dna).reverse_complement())
                            else:
                                dna = dna1
                            assert dna == seq

                            if diff(dna, cSeq) > maxDiff and ('I' not in cigar and 'D' not in cigar):
                                over += 1
                                if len(dna) != len(cSeq):
                                    overDiffSize += 1
                                # print readName
                                # print dna
                                # print cSeq
                                # print cigar
                                # print diff(dna, cSeq)
                        else:
                            nameSuffix = '\\2'
                            if rc:
                                dna = str(Seq(dna2, generic_dna).reverse_complement())
                            else:
                                dna = dna2
                            assert dna == seq

                            if diff(dna, cSeq) > maxDiff and ('I' not in cigar and 'D' not in cigar):
                                over += 1
                                if len(dna) != len(cSeq):
                                    overDiffSize += 1

                        # creates a read record
                        rRec = CRec(readName + nameSuffix, dna, pos - 1)

                        # creates a contig record
                        if cName not in cNameToRec:
                            cRec = CRec(cName, contigSeq)
                            cNameToRec[cName] = cRec
                        else:
                            cRec = cNameToRec[cName]

                        # add the read record to the contig record
                        cRec.addRead(rRec)
                    else:
                        print 'unmapped'

                # initialize the position coverage array
                for cRec in cNameToRec.values():
                    cRec.getPosCovArray()

                assert pfName not in pfToRecSet
                pfToRecSet[pfName] = set(cNameToRec.values())

                # for debugging
                # print pfToRecSet[pfName]
                # for c in pfToRecSet[pfName]:
                #     print c.recordId
                #     print c.dnaSeq
                #     for i in c.readRecList:
                #         print str(' ' * i.posWithinContig) + i.dnaSeq + '   %s  %s' % (i.recordId, i.posWithinContig)
                #     print c.getPosCovArray()
                # return

        # store records to files
        if dstDir is not None:
            if not os.path.isdir(dstDir):
                os.mkdir(dstDir)

            for pfName, recSet in pfToRecSet.iteritems():
                dstFile = os.path.join(dstDir, 'contigs_%s.pkl.gz' % pfName)
                hio.storeObj(recSet, dstFile)

        print 'STAT: over: %s overDiffSize: %s total: %s' % (over, overDiffSize, total)
    except Exception as e:
        print samDir, fqPartitionDir, satOutDir, dstDir, maxDiff
        print e.message
        print e.args
        raise e


def getReadGmap(sampleDir, dstGmap):
    """
        Get the mapping of reads onto reference sequences of individual strains.
        Map: read -> (refSeq, refStrain)
        The strainId is the last entry in the sam file.

        @type sampleDir: str
        @type dstGmap: str
    """
    # NZ_KB001324.1-231       99      NZ_KB001324.1   6627    99      45=1X60=1X21=1X10=1X10= =
    # 2353    -4424   CAGAATGTGTAATCAGGTAACTGGCAAGCTTTTGCGTAAAGTAGCGTTCGCCGCGCAGGACGCCTTGCAATCCATTGACAACGCGTTCTTGATCCTCCATGGCATATAAAACGCCGTTGATATGAGGCTAGTTTTCAATCTCGCGGTACG
    # CCCGCCG=GGGGGJJGJJJJGJCGJCJJGJJGGGGJJGCGGGJJG$G=JCJCGJCCJJJJJCGG=C$GGGGGGCGCGC1GGGCCGCCGGCGCGGGGCGGGJCCGG1$CGCGGGGGCGCGGGGGGGCCG$GG=G=GCGGG$GCGG$GCGCC
    # *       NZ_ANLR00000000
    out = fq.WriteFq(dstGmap)

    for strainId in os.listdir(sampleDir):
        samPath = os.path.join(sampleDir, strainId, '0_pair.sam.gz')
        fq1Path = os.path.join(sampleDir, strainId, '0_pair1.fq.gz')
        fq2Path = os.path.join(sampleDir, strainId, '0_pair2.fq.gz')
        if os.path.isfile(samPath) and os.path.isfile(fq1Path) and os.path.isfile(fq2Path):
            # print samPath
            # for each sam record
            for line in gzip.open(samPath):
                line = line.strip()
                if not line.startswith('@'):  # skip comments

                    tokens = line.split('\t')
                    readName, flag, contigName, pos, m, cigar, rnext, pnext, tlen, seq, qs = tokens[:11]
                    flag = int(flag)
                    if flag & 0x40 != 0:  # first segment
                        readName += '\\1'
                    else:
                        assert flag & 0x80 != 0  # last segment
                        readName += '\\2'
                    seq = qs = geneAnnot = '*'

                    entry = '\t'.join(map(lambda x: str(x), [readName, flag, contigName, pos, m, cigar, rnext, pnext,
                                                             tlen, seq, qs, geneAnnot, strainId])) + '\n'
                    out.write(entry)

    out.close()


def countError(contRecDir, refGmap, refSeqBuff, storeLab=True, maxCov=30, maxProc=4):
    """
        Compute the per-base error.

        @type contRecDir: str
        @type refGmap: str
        @type refSeqBuff: dict[(str,str),str]
        @type storeLab: bool
        @type maxCov: int
        @type maxProc: int
        @rtype: list[ndarray]
    """
    assert os.path.isdir(contRecDir)
    assert os.path.isfile(refGmap)
    taskList = []
    # c = 0
    for f in os.listdir(contRecDir):
        readRecPkl = os.path.join(contRecDir, f)
        if os.path.isfile(readRecPkl) and readRecPkl.endswith('pkl.gz'):
            # c += 1
            # if c < 200:
            #     continue

            taskList.append(parallel.TaskThread(heval.getPerBaseErrorPkl, (readRecPkl, refGmap, refSeqBuff, maxCov,
                                                                           comh.TRANSLATION_TABLE, storeLab)))
            # error = heval.getPerBaseErrorPkl(readRecPkl, refGmap, refSeqBuff, maxCov, comh.TRANSLATION_TABLE, storeLab)
            # if c > 1000:
            #     break

    return parallel.runThreadParallel(taskList, maxProc)


def tmpF(contRecDir, refGmap, refSeqBuff, storeLab=True, maxCov=30, maxProc=4, pfToRecSet=None, labContigsPklFile=None):
    # TODO: remove
    # contigs = 0
    # total = 0
    # errAList = []
    taskList = []
    # c = 0
    # count = 0
    # maxi = 1000000
    for pfName, recSet in pfToRecSet.iteritems():
        # print pfName, len(recSet)
        # contigs = len(recSet)
        # total += 1
        recSet2 = set()

        for rec in recSet:
            # if len(rec.dnaSeq) > 0:
                # print len(rec.dnaSeq)
                # count += 1
                # maxi = min(maxi, len(rec.dnaSeq))

            if rec.readRecList is not None:  # dnaSeq
                recSet2.add(rec)

        if len(recSet2) > 0:
            # print contigs
            # errA = heval.getPerBaseError(recSet2, refGmap, refSeqBuff, maxCov, allLabels=True)
            taskList.append(parallel.TaskThread(heval.getPerBaseError, (recSet2, refGmap, refSeqBuff, maxCov, True)))
            # c += 1
            # if c > 100:
            #     break
            # print errA
            # errAList.append(errA)
            break  # TODO: !!!
    # print '1000 contigs: %s %s' % (count, maxi)

    print 'computing'
    errAList = parallel.runThreadParallel(taskList, maxProc)
    # print report
    print heval.getAssemblyReport(errAList, maxCov)

    hio.storeObj(pfToRecSet, labContigsPklFile)

    print 'loading'
    d2 = hio.loadObj(labContigsPklFile)
    for k, v in d2.iteritems():
        print k
        for i in v:
            if i.labelEval is not None:
                for e in i.labelEval:
                    print e.getContent()

    # print len(pfToRecSet)

    # print 'all', contigs, total, float(contigs) / total


def _test1():

    fnaFile = '/home/igregor/Documents/work/hsim/562/samples/0/sat/0_sat.fna'
    gmapFile = '/home/igregor/Documents/work/hsim/562/samples/0/sat/0_gmap.sam.gz'
    contigsPklFile = '/home/igregor/Documents/work/hsim/562/samples/0/sat/0_sat_contigs.pkl.gz'
    labContigsPklFile = '/home/igregor/Documents/work/hsim/562/samples/0/sat/0_sat_contigs_lab.pkl.gz'
    sampleDir = '/home/igregor/Documents/work/hsim/562/samples/0'
    fqPairList = [
        ('/home/igregor/Documents/work/hsim/562/samples/0/NZ_ANLR00000000/0_pair1.fq.gz',
    '/home/igregor/Documents/work/hsim/562/samples/0/NZ_ANLR00000000/0_pair2.fq.gz'),
    ('/home/igregor/Documents/work/hsim/562/samples/0/NZ_AKLX00000000/0_pair1.fq.gz',
    '/home/igregor/Documents/work/hsim/562/samples/0/NZ_AKLX00000000/0_pair2.fq.gz'),
    ('/home/igregor/Documents/work/hsim/562/samples/0/NZ_AKLB00000000/0_pair1.fq.gz',
    '/home/igregor/Documents/work/hsim/562/samples/0/NZ_AKLB00000000/0_pair2.fq.gz')]
    # fqToFasta(fqPairList, fnaFile)  # FASTA for SAT

    satOutDir = '/home/igregor/Documents/work/hsim/562/samples/0/sat/0_sat'
    hmmProfiles = os.path.join(comh.HMM_PROFILE_DIR, comh.HMM_PROFILE_FILE)
    satInstallDir = comh.SAT_INSTALL_DIR
    hmmerBinDir = comh.HMMER_BINARY
    # runSat(fnaFile, satOutDir, hmmProfiles, satInstallDir, hmmerBinDir, options='')  # SAT

    fqPartitionDir = '/home/igregor/Documents/work/hsim/562/samples/0/sat/fq_part'
    aliIndexDir = '/home/igregor/Documents/work/hsim/562/samples/0/sat/ali_idx'
    samDir = '/home/igregor/Documents/work/hsim/562/samples/0/sat/sam'
    # satMapReads(satOutDir, fqPairList, fqPartitionDir, aliIndexDir, samDir)  # BWA

    # satToContigRecords(samDir, fqPartitionDir, satOutDir, contigsPklFile)

    # getReadGmap(sampleDir, gmapFile)

    refSeqBuff = fas.getSequenceBuffer([os.path.join(comh.REFERENCE_DIR_ROOT, '562', comh.FASTA_GENOMES_DIR_NAME),
                                        os.path.join(comh.REFERENCE_DIR_ROOT, '562', comh.FASTA_GENOMES_DRAFT_DIR_NAME)])

    countError(contigsPklFile, gmapFile, refSeqBuff, labContigsPklFile)


def test_sat_out():
    satOutDir = '/Users/ivan/Documents/nobackup/hsim01/562/samples10/_0/sat/0_sat'

    pfSet = set()
    for f in os.listdir(satOutDir):
        if f.endswith('contigs.fa'):
            pfSet.add(f.split('_')[0])

    c = 0
    t = 0
    lim = 1024 * 3
    for f in os.listdir(satOutDir):
        if f.endswith('hmmer'):
            pf = f[1:].split('.')[0]
            fPath = os.path.join(satOutDir, f)
            size = os.path.getsize(fPath)
            if size > lim and (pf not in pfSet):
                # print 'attention'
                c += 1
            if size > lim:
                t += 1
    print c, t


def _test():
    f = '/home/igregor/Documents/work/hsim/562/samples/0/sat/cont_pkl/contigs_PF13558.pkl.gz'
    s = hio.loadObj(f)
    for r in s:
        print r.labelEval

if __name__ == "__main__":
    satToContigRecords('/home/igregor/Documents/work/hsim/562/samples/3/sat/sam',
                       '/home/igregor/Documents/work/hsim/562/samples/3/sat/fq_part',
                       '/home/igregor/Documents/work/hsim/562/samples/9/sat/0_sat',
                       '/home/igregor/Documents/work/hsim/562/samples/3/sat/cont_pkl',
                       10)
    # _test()
    # test_sat_out()
