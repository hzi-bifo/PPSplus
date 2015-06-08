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

    Handle read Hmm alignments.
"""
import os
import gzip
import numpy as np
import multiprocessing as mp
from algbioi.com import csv
from algbioi.com import fq
from algbioi.com import fasta as fas
from algbioi.com import parallel
from algbioi.hsim import comh
from algbioi.haplo import hio


def alignAllPartitioned(partDir, hmmProfileDir, hmmerBinDir, maxProc=mp.cpu_count()):
    """
        For each tag-domain file in the read partitioned directory, hmm-align all the corresponding prot-reads.
        The ending of the space-separated (seq-name align-seq) alignment files is 'read_ali_hmm.gz'.
        Moreover alignments are sorted according to the hmm coordinates.

        @param partDir: directory containing reads partitioned into domains
        @param hmmProfileDir: directory containing domain-hmm-profiles
        @param hmmerBinDir: directory containing HMMER binaries
        @type partDir: str
        @type hmmProfileDir: str
        @type hmmerBinDir: str
    """
    assert os.path.isdir(partDir)
    assert os.path.isdir(hmmProfileDir)
    assert os.path.isdir(hmmerBinDir)

    # read in all dom-tag names
    domTagSet = set()
    for f in os.listdir(hmmProfileDir):
        domTagSet.add(f[:-4])

    # read in all src files
    fileSet = set()
    for f in os.listdir(partDir):
        if f.endswith('_prot.domtblout.gz') or f.endswith('_prot.fna.gz'):
            fileSet.add(f)

    # for each dom-tag check whether there are corresponding files
    taskList = []
    for domTag in domTagSet:
        fList = None
        domJ = 'r_%s_join_prot.domtblout.gz' % domTag
        if domJ in fileSet:
            protJ = 'r_%s_join_prot.fna.gz' % domTag
            assert protJ in fileSet
            fList = [(os.path.join(partDir, protJ), os.path.join(partDir, domJ))]
        domP = 'r_%s_pair_prot.domtblout.gz' % domTag
        if domP in fileSet:
            protP = 'r_%s_pair_prot.fna.gz' % domTag
            assert protP in fileSet
            entry = (os.path.join(partDir, protP), os.path.join(partDir, domP))
            if fList is None:
                fList = [entry]
            else:
                fList.append(entry)

        # define task for this domain-tag
        if fList is not None:
            profileHmm = os.path.join(hmmProfileDir, '%s.hmm' % domTag)
            outAliFilePath = os.path.join(partDir, 'r_%s_read_ali_hmm.gz' % domTag)
            taskList.append(parallel.TaskThread(alignHmm, (fList, profileHmm, hmmerBinDir, outAliFilePath)))

    # run all tasks (alignments) in parallel
    parallel.runThreadParallel(taskList, maxProc)


def alignHmm(fList, profileHmm, hmmerBinDir, outAliFilePath, verbose=False):
    """
        Align reads given a list of pairs (prot_fasta_file, domtblout_file).

        @param fList: list of (prot_fasta_file, domtblout_file)
        @param profileHmm: file containing a HMM profile
        @param hmmerBinDir: directory containing HMMER binaries
        @param outAliFilePath: output alignment file (name SPACE alignment-sequence)

        @type fList: list[(str,str)]
        @type profileHmm: str
        @type hmmerBinDir: str
        @type outAliFilePath: str
    """
    hmmAlignBin = os.path.join(hmmerBinDir, 'hmmalign')
    assert os.path.isfile(profileHmm)
    assert os.path.isfile(hmmAlignBin)
    assert os.path.isdir(os.path.dirname(outAliFilePath))

    # read in prot sequences
    protSeqDict = {}
    protSeqDictLen = 0
    for prot in map(lambda x: x[0], fList):
        assert os.path.isfile(prot), '%s' % (prot)
        s = fas.fastaFileToDictWholeNames(prot)
        protSeqDictLen += len(s)
        protSeqDict.update(s)
    assert len(protSeqDict) == protSeqDictLen

    # read in dom files
    domDict = {}
    domDictLen = 0
    for dom in map(lambda x: x[1], fList):
        assert os.path.isfile(dom), '%s %s' % (dom)
        d = hio.readDomblout(dom)
        domDictLen += len(d)
        domDict.update(d)
    assert len(domDict) == domDictLen

    # create a FASTA file for the aligner
    buff = []
    for name, seq in protSeqDict.iteritems():
        hit = domDict[name][0]
        protStart = int(hit[19]) - 1
        protLen = int(hit[20]) - protStart
        # take only the substring that has been annotated
        buff.append('>%s\n%s\n' % (name, seq[protStart:protStart + protLen]))

    # randomly shuffle the fasta file records (to remove the order imposed by sequence names as stored in a dictionary)
    rand = np.random.RandomState(domDictLen)
    rand.shuffle(buff)
    protFileTmp = outAliFilePath + '_prot_tmp'
    # store records to a file
    out = csv.OutFileBuffer(protFileTmp)
    out.writeText(''.join(buff))
    out.close()

    # align the sequences
    aliFileTmp = outAliFilePath + '.ali_tmp'
    cmd = '%s -o %s --amino --informat FASTA --outformat Pfam %s %s; rm %s' \
        % (hmmAlignBin, os.path.basename(aliFileTmp), profileHmm, os.path.basename(protFileTmp),
           os.path.basename(protFileTmp))  # --allcol
    cwd = os.path.dirname(protFileTmp)
    parallel.reportFailedCmd(parallel.runCmdSerial([parallel.TaskCmd(cmd, cwd)]))
    # print cmd

    # read in alignments
    nameToAli = aliToDict(aliFileTmp)

    # remove the tmp alignment file
    os.remove(aliFileTmp)

    # store alignments to a file sorted according to the
    buff = []
    maxNameLen = 0
    for name, alignment in nameToAli.iteritems():
        hit = domDict[name][0]
        hmmCoordStart = int(hit[15]) - 1
        buff.append((name, alignment, hmmCoordStart))
        if len(name) > maxNameLen:
            maxNameLen = len(name)
    buff.sort(key=lambda x: x[2])

    maxNameLen += 3
    out = fq.WriteFq(outAliFilePath)
    for name, alignment, h in buff:
        padding = (maxNameLen - len(name)) * ' '
        out.write('%s%s%s\n' % (name, padding, alignment))
    out.close()

    # print the alignment with additional info
    if verbose:
        buff = []
        for name, ali in nameToAli.iteritems():
            hit = domDict[name][0]
            protStart = int(hit[19]) - 1
            protLen = int(hit[20]) - protStart

            protStartAli = int(hit[17]) - 1
            protLenAli = int(hit[18]) - protStartAli

            hmmCoordStart = int(hit[15]) - 1
            hmmCoordLen = int(hit[16]) - hmmCoordStart

            cut = protSeqDict[name][protStart: protStart + protLen]
            assert len(cut) == protLen

            header = name + ' hmm: ' + str(hmmCoordStart) + ' ' + str(hmmCoordLen) + ' prot: ' + str(protStart) + ' ' \
                     + str(protLen) + ' ali: ' + str(protStartAli) + ' ' + str(protLenAli)
            buff.append((header + (55 - len(header)) * ' ' + '\t' + ali, hmmCoordStart, name))
            # buff.sort(key=lambda x: x[2])
            buff.sort(key=lambda x: x[1])

        print '\n'.join(map(lambda x: x[0], buff))


def aliToDict(aliFile):
    """
        Reads alignments in a Pfam (Stockholm) format into a dictionary.

        @param aliFile: file in Pfam (Stockholm) format
        @type: str
        @rtype: dict[str,str]
        @return: map: seqName -> alignment sequence
    """
    assert os.path.isfile(aliFile)
    d = {}
    if aliFile.endswith('.gz'):
        f = gzip.open(aliFile)
    else:
        f = open(aliFile)
    for line in f:
        line = line.strip()
        if line.startswith('#') or line == '' or line.startswith('//'):
            continue
        name, ali = line.split()
        assert name not in d
        d[name] = ali
    return d


def test1():
    # geneTag = 'rna_pol_rpb6_4'
    geneTag = 'phage_antiter_q_2'
    d = '/home/igregor/Documents/work/hsim/562/samples/0/sample_partitioned'
    profileHmm = os.path.join(comh.HMM_PROFILE_DIR, comh.HMM_PROFILE_FILES_PARTITIONED_DIR, '%s.hmm' % geneTag)
    fList = [(os.path.join(d, 'r_%s_join_prot.fna.gz' % geneTag), os.path.join(d, 'r_%s_join_prot.domtblout.gz' % geneTag)),
         (os.path.join(d, 'r_%s_pair_prot.fna.gz' % geneTag), os.path.join(d, 'r_%s_pair_prot.domtblout.gz' % geneTag))]

    outAliFilePath = os.path.join(d, 'r_%s_read_ali_hmm.gz' % geneTag)

    alignHmm(fList, profileHmm, comh.HMMER_BINARY, outAliFilePath, verbose=True)

    print outAliFilePath

if __name__ == "__main__":
    pass
    # test1()

