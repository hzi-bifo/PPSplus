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

    Note that we could have written some parts of this code in a nicer way,
    but didn't have time. Be careful when reusing the source code.


    Miscellaneous functions to test the Illumina ART simulator.
"""
import os
import re
import sys
import subprocess
import time

ART_BINARY = '/net/metagenomics/projects/PPSmg/tools/art/art_bin_VanillaIceCreamLinux/art_illumina'


def runCommand(cmdList, cwd):
    counter = 0
    for cmd in cmdList:
        counter += 1
        process = subprocess.Popen(cmd, shell=True, bufsize=-1, cwd=cwd)
        print('run "%s" cmd: %s' % (counter, cmd))
        process.wait()
        if process.returncode != 0:
            print('Process: "%s" ended with return code: "%s' % (cmd, process.returncode))


def noNewLine(line):
    """
        Delete all '\n' and '\r' characters in a string.
    """
    return line.replace('\n', '').replace('\r', '')


def fastaFileToDictWholeNames(filePath):
    """
        Reads a fasta file and returns mapping: seqName -> sequence the whole sequence name is used
        as seqName!!! (even if it contains space)
    """
    seqIdToSeq = {}
    f = None
    try:
        f = open(os.path.normpath(filePath), 'r')
    except Exception:
        print "Cannot open file:", filePath
        raise
    else:
        name = ''
        seq = ''
        for line in f:
            line = noNewLine(line)
            if re.match('>', line):
                if seq != '':
                    assert name != ''
                    seqIdToSeq[name] = seq
                    seq = ''
                name = line.replace('>', '')
            else:
                seq += line
        if seq != '':
            assert name != ''
            seqIdToSeq[name] = seq
    finally:
        if f is not None:
            f.close()
    return seqIdToSeq


class OutFileBuffer():
    """
        To append text to a file.
    """
    def __init__(self, outFilePath, bufferText=False, fileOpenMode='w'):
        self.outFilePath = outFilePath
        self.empty = True
        self.bufferText = bufferText
        self.textBuffer = ''
        self.opened = False
        try:
            self.outFile = open(os.path.normpath(self.outFilePath), fileOpenMode)
            self.opened = True
        except Exception:
            sys.stderr.write('Cannot open a file for writing: ' + outFilePath)
            raise

    def writeText(self, text):
        try:
            if not self.opened:  # reopen to append
                self.outFile = open(os.path.normpath(self.outFilePath), 'a')
                self.opened = True
            self.outFile.write(text)
            if self.bufferText:
                self.textBuffer += text
            if self.empty:
                self.empty = False
        except Exception:
            sys.stderr.write('Cannot write to file: ' + self.outFilePath)
            self.close()
            raise

    def getTextBuffer(self):
        return self.textBuffer

    def isEmpty(self):
        return self.empty

    def close(self):
        self.outFile.close()
        self.opened = False


class ProcHandler():
    """
        Stores context of one process.
    """
    def __init__(self, cmd, cwd, timeout=None):
        self.process = subprocess.Popen('exec ' + cmd, shell=True, bufsize=-1, cwd=cwd)
        self.cmd = cmd
        self.runtime = 0.
        self.timeout = timeout

    def incRuntime(self, timeStep):
        self.runtime += timeStep

    def isTimeOut(self):
        if self.timeout is not None and self.runtime > self.timeout:
            return True
        return False

    def getPid(self):
        return self.process.pid


def runCommandParallel(cmdList, cwdList, maxProc=2, timeout=None, timeStep=0.5):
    """
        Run several commands in parallel.

        @param cmdList: list of commands
        @param cwdList: corresponding list of current working directories in which commands will be run
        @param maxProc: maximum number of processes that will be run in parallel
        @param timeout: after this number of seconds, the process will be killed, (None if no timeout set)
        @param timeStep: time interval in which processes will be actively checked whether they are running
        @return: list of failed commands, tuple (command, return code, runtime)
    """
    counter = 0
    failList = []
    cmdInGen = True
    cmdGen = iter(zip(cmdList, cwdList))  # generator of commands
    procArray = {}
    for i in range(maxProc):
        procArray[i] = None

    # loop until all processes finish (or are killed)
    while True:

        # run commands
        if cmdInGen:
            for i in range(maxProc):
                if procArray[i] is None:
                    try:
                        cmd, cwd = cmdGen.next()
                        procArray[i] = ProcHandler(cmd, cwd, timeout)  # run process
                        counter += 1
                        print('Running "%s" cmd: %s' % (counter, cmd))
                    except StopIteration:
                        cmdInGen = False  # there are no processes to be run

        # sleep for a while
        time.sleep(timeStep)

        # check for finished processes and processes passed timeout
        for i in range(maxProc):
            ph = procArray[i]
            if ph is not None:  # there is a process in the slot
                if ph.process.poll() is None:
                    # process running
                    ph.incRuntime(timeStep)
                    if ph.isTimeOut():
                        # process over time, kill it!
                        ph.process.kill()
                        print("Process (%s): %s killed! (after %ss)" % (ph.getPid(), ph.cmd, ph.runtime))
                        failList.append((ph.cmd, 9, ph.runtime, ph.getPid()))
                        procArray[i] = None  # free slot
                else:
                    # process finished
                    ph.process.wait()
                    if ph.process.returncode != 0:
                        print('Process(%s): "%s" ended with return code: "%s'
                              % (ph.getPid(), ph.cmd, ph.process.returncode))
                        failList.append((ph.cmd, ph.process.returncode, ph.runtime, ph.getPid()))
                    procArray[i] = None  # free slot

        # finish if no process is running and there is not process to be run
        if len(set(procArray.values())) == 1 and not cmdInGen:
            break

    return failList


def getSimsets(inDir, outDir):
    """


    """
    assert os.path.isfile(ART_BINARY), 'Binnary file does not exist: %s' % ART_BINARY
    # count = 0
    seed = 123
    coverage = 10
    minLen = 350
    cmdListA = []
    cwdList = []
    for fileName in os.listdir(inDir):
        filePath = os.path.join(inDir, fileName)
        cwd = os.path.join(outDir, fileName.rsplit('.', 1)[0])
        if not os.path.isdir(cwd):
            os.mkdir(cwd)
        fna = os.path.join(cwd, fileName)
        eraseShortSeq(filePath, fna, minLen)

        cmd = "%s -i %s -f %s -rs %s  -l 100 -m 450 -s 25 -p -o pair_ -sam " % (ART_BINARY, fna, coverage, seed)
        cmdListA.append(cmd)
        cwdList.append(cwd)

    failList = runCommandParallel(cmdListA, cwdList, maxProc=50, timeout=600, timeStep=0.2)
    print('Failed: ')
    for f in failList:
        print('fail %s' % str(f))


def eraseShortSeq(inFile, outFile, minLen):
    out = OutFileBuffer(outFile)
    for k, v in fastaFileToDictWholeNames(inFile).iteritems():
        if len(v) >= minLen:
            out.writeText('>%s\n%s\n' % (k, v))
    out.close()


def test1():
    # inDir = '/net/metagenomics/projects/PPSmg/hsim/nobackup/562/fasta_genomes'
    inDir = '/net/metagenomics/projects/PPSmg/hsim/nobackup/562/fasta_draft_genomes'
    # inDir = '/net/metagenomics/projects/PPSmg/hsim/arttest/chunks'
    outDir = '/net/metagenomics/projects/PPSmg/hsim/nobackup/562/test_simset'
    # outDir = '/net/metagenomics/projects/PPSmg/hsim/arttest/chunks_sim'
    getSimsets(inDir, outDir)

# def test2():
#     f = '/net/metagenomics/projects/PPSmg/hsim/arttest/NZ_AMUV00000000.fna'
    # f = '/Users/ivan/Documents/work/tools/art/test/NZ_AMUV00000000.fna'
    # dst = '/Users/ivan/Documents/work/tools/art/test/chunks'
    # from algbioi.com import fasta as fas
    # from algbioi.com import csv as csv
    #
    # for k, v in fas.fastaFileToDictWholeNames(f).iteritems():
    #     out = csv.OutFileBuffer(os.path.join(dst, k))
    #     out.writeText('>%s\n%s' % (k,v))
    #     out.close()


def test3():
    from algbioi.com import fasta as fas

    s = set()
    for line in open('/Users/ivan/Documents/work/tools/art/test/fail_list.txt'):
        s.add(line.strip())

    chunks = '/Users/ivan/Documents/work/tools/art/test/chunks'
    c = set(os.listdir(chunks))

    intersect = set()
    for i in c:
        if i not in s:
            intersect.add(i)

    lengths = []
    for n in intersect:
        for v in fas.fastaFileToDictWholeNames(os.path.join(chunks, n)).values():
            lengths.append(len(v))

    print min(lengths)
    print max(lengths)
    print float(sum(lengths)) / float(len(lengths))


def test4():
    from algbioi.com import fasta as fas
    allG = '/Users/ivan/Documents/nobackup/hsim01/562/phylo_acc_all_genes'
    l = []
    longS = 0
    shortS = 0
    for f in os.listdir(allG):
        filePath = os.path.join(allG, f)
        if os.path.isfile(filePath):
            for v in fas.fastaFileToDictWholeNames(filePath).values():
                l.append(len(v))
                if len(v) > 300:
                    longS += 1
                else:
                    shortS += 1

    print min(l)
    print max(l)
    print float(sum(l)) / float(len(l))
    print longS
    print shortS


def test5():
    from algbioi.com import fasta as fas
    genes = '/Users/ivan/Documents/nobackup/hsim01/562/phylo_acc_all_genes'
    draftGenomes = '/Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes'

    # accessions in genes
    genesAccSet = set()
    for f in os.listdir(genes):
        filePath = os.path.join(genes, f)
        for k in fas.fastaFileToDictWholeNames(filePath).keys():
            acc = k.split(':', 1)[1].split(';')[0]
            genesAccSet.add(acc)

    # accessions in draft genomes shorter than XXX
    minLen = 400
    dgAccSet = set()
    d = {}
    for f in os.listdir(draftGenomes):
        filePath = os.path.join(draftGenomes, f)
        for k, v in fas.fastaFileToDictWholeNames(filePath).iteritems():
            if len(v) < minLen:
                acc = k.lstrip('>')
                d[acc] = len(v)
                dgAccSet.add(acc)

    intersect = set()
    for i in dgAccSet:
        if i in genesAccSet:
            intersect.add(i)

    print len(intersect)
    for i in intersect:
        print d[i]


def testART(specDir, testDir):
    assert False  # mess !!! anything useful ???
    # src = const.FASTA_GENOMES_DIR_NAME
    src = None  # !!!! ???

    cmdList = []
    avoidSet = set()
    avoidSet.add('NZ_ABHO00000000.fna')
    avoidSet.add('NZ_ABHT00000000.fna')
    avoidSet.add('NZ_ABHW00000000.fna')
    avoidSet.add('NZ_ACGN00000000.fna')

    for f in os.listdir(os.path.join(specDir, src)):
        if f in avoidSet:
            continue
        srcFna = os.path.join(os.path.join(specDir, src, f))
        cmd = "%s -i %s -f %s -rs %s  -l 100 -m 450 -s 25 -p -o pair_ -sam " % (ART_BINARY, srcFna, 1, 123)
        cmdList.append(cmd)
    # cmdList = cmdList[16:]  # !!!!!!!!!!!!!!!!!!!!!!!!!
    # com.runCommand(cmdList, testDir)
    cmdList = ['/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_ADBT00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_ADUK00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_ADUL00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_ADUM00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_ADUN00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_ADUO00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_ADUP00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_ADUQ00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_ADUR00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_ADUX00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_ADUY00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_ADUZ00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_ADVB00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_ADVC00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AETX00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AETY00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AETZ00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AEUA00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AEUB00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AEUC00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AEZL00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AEZM00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AEZN00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AEZP00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AEZR00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AEZV00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AEZY00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AFDT00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AFDX00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AFDY00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AFOB00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AFPS00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AFVR00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AFWC00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AFWO00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AFWP00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AGSG00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AGTG00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AGTH00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AGTI00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AGTJ00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AGTK00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AGTL00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AHAV00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AIEX00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AIFA00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AIFE00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AIFH00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AIFL00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AIFN00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AIFQ00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AIFU00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AIFV00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AIFW00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AIFX00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AIGA00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AIGB00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AIGE00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AIGG00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AIGH00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AIGR00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AIGW00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AIHA00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AIHC00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AIHD00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AIHE00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AIHL00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AIHM00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AJMB00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AJQW00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AJVQ00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AJVR00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AJVS00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AJVT00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AJVU00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AJVV00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AJVW00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AJVX00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AJWO00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AJWP00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AJWQ00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AJWR00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AJWS00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AJWT00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AJWU00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AJWV00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AKAV00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AKAW00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AKAX00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AKAY00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AKAZ00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AKBA00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AKBB00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AKBC00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ',
               '/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_AKKW00000000.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ']

    # cmdList = cmdList[1:2]

    strains = '/Users/ivan/Documents/nobackup/hsim01/562/phylo_strains_cluster/strains.fna'
    strainSet = set()
    for line in fas.fastaFileToDictWholeNames('/Users/ivan/Documents/nobackup/hsim01/562/phylo_strains_cluster/strains.fna').keys():
        strainSet.add(line)
    print len(strainSet)

    print len(cmdList)
    l = []
    for i in cmdList:
        e = os.path.basename(i.split()[2]).rsplit('.')[0]
        if e in strainSet:
            l.append(e)
            print e
    print len(l)
    return

    testFileOrig = '/Users/ivan/Documents/nobackup/hsim01/562/fasta_draft_genomes/NZ_ADBT00000000.fna'  # ???
    testFileCp = os.path.join(testDir, 'test.fna')
    out = csv.OutFileBuffer(testFileCp)
    first = True
    for name, seq in fas.fastaFileToDictWholeNames(testFileOrig).iteritems():

        if first:
            out.writeText('>%s' % name)
            first = False
        else:
            out.writeText('\n>%s' % name)

        s = seq
        l = len(s)
        i = 0
        range = 80
        while i < l:
            j = i + range
            out.writeText('\n' + s[i:j])
            i += range

    out.close()

    # cmdList = ['/Users/ivan/Documents/work/tools/art/art_bin_VanillaIceCreamOSX/art_illumina -i /Users/ivan/Documents/nobackup/hsim01/562/samples/test.fna -f 1 -rs 123  -l 100 -m 450 -s 25 -p -o pair_ -sam ']
    print('size %s' % len(cmdList))
    counter = 0
    maxWait = 10.0
    killedList = []
    for cmd in cmdList:
        counter += 1
        wait = 0
        process = subprocess.Popen(cmd, shell=True, bufsize=-1, cwd=testDir)
        print('run "%s" cmd: %s' % (counter, cmd))
        while True:
            r = process.poll()
            if r is not None:
                break
            if wait > maxWait:
                process.kill()
                print('Process KILLED: %s' % cmd)
                killedList.append(cmd)
                break
            wait += 0.5
            time.sleep(0.5)
        process.wait()
        if process.returncode != 0:
            print('Process: "%s" ended with return code: "%s' % (cmd, process.returncode))
    print('Killed list:')
    print killedList


# are there some genes on the short sequences ?????


if __name__ == "__main__":
    test1()
    # test2()
    # test3()
    # test4()
    # test5()

