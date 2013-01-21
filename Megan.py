#!/usr/bin/env python

import re
import os
import glob
import Common
from sets import Set

def writeSequences(outFile, contigsList, namesToSeq, maxLineLen):
    try:
        f = open(os.path.normpath(outFile), 'w')
        k = 0
        range = maxLineLen
        for name in contigsList:
            if name not in namesToSeq:
                print str('Contig name "' + name + '" is not in the fasta file')
                continue
            s = namesToSeq[name]
            if k == 0:
                f.write('>' + name)
                k += 1
            else:
                f.write('\n>' + name)
            l = len(s)
            i = 0
            while i < l:
                j = i + range
                f.write('\n' + s[i:j])
                i += range
    except Exception:
        print "Cannot create a file or write to it:", outFile
        raise
    finally:
        f.close()


def fillInSeq(inDir, outDir, mapFastaFile, ignoreContigsFile, scaffoldContigsFile=None):
    """
        For each file in the inDir that contains contig/sequence names produces one file in the outDir that
        contains the appropriate sequences (i.e. >sequence_name \n sequence).

        @param inDir: directory that contains *.fasta files that contain only contig names
        @param outDir: directory to store output *.fasta files with the corresponding sequences
        @param mapFastaFile: file that contains all sequences (sequence/contig names + sequences) that can be contained in the inDir
        @param ignoreContigsFile: file that contains contigsNames that won`t be contained in any file in the outDir
        @param scaffoldContigsFile: file that contains: scaffold -> contigs mapping
    """
    nameToSeq = dict([])
    try:
        f = open(os.path.normpath(mapFastaFile),'r')
    except Exception:
        print "Cannot open file:", mapFastaFile
        raise
    else:
        name = ''
        seq = ''
        maxLineLen = 0
        for line in f:
            line = Common.noNewLine(line)
            if maxLineLen < len(line):
                maxLineLen = len(line)
            if re.match('>', line):
                if seq != '':
                    assert name != ''
                    nameToSeq[name] = seq #store seq
                    seq = ''
                name = line.replace('>','')
            else:
                seq += line
        if seq != '':
            assert name != ''
            nameToSeq[name] = seq #store seq
    finally:
        f.close()

    ignoreContigsSet = None
    if ignoreContigsFile != None:
        try:
            f = open(os.path.normpath(ignoreContigsFile),'r')
        except Exception:
            print "Cannot open file:", ignoreContigsFile
            raise
        else:
            ignoreContigsSet = set([])
            for line in f:
                ignoreContigsSet.add(Common.noNewLine(line))

    for filePath in glob.glob(os.path.join(os.path.normpath(inDir),r'*.fna')):
        try:
            f = open(os.path.normpath(filePath),'r')
        except Exception:
            print "Cannot open file:", filePath
            raise
        else:
            contigsList = []
            for line in f:
                line = Common.noNewLine(line)
                if (len(line) != 0) and not re.match(' ', line):
                    contigsList.append(line.replace('>',''))
                    if not re.match('>', line):
                        print 'The line doesn`t start with ">"', line

            if ignoreContigsSet != None:
                tempList = []
                for c in contigsList:
                    if c not in ignoreContigsSet:
                        tempList.append(c)
                    else:
                        print c, 'ignored'
                contigsList = tempList

            if scaffoldContigsFile != None:
                contigsList = addContigsFromTheSameScaffold(contigsList, scaffoldContigsFile)

            outFile = os.path.normpath(os.path.join(outDir, os.path.basename(filePath)))
            writeSequences(outFile, contigsList, nameToSeq, maxLineLen)


def addContigsFromTheSameScaffold(contigsList, scaffoldContigsFile):
    conToScaf = dict([])
    scafToContigs = dict([])
    try:
        f = open(os.path.normpath(scaffoldContigsFile),'r')
    except Exception:
        print "Cannot open file:", scaffoldContigsFile
        raise
    else:
        for line in f:
            line = Common.noNewLine(line)
            scaffold = re.sub(r'^[ ]*([^ ,]+),[^,]*',r'\1', line)
            contig = re.sub(r'^[^,]*,([^ ,]+)[ ]*',r'\1', line)
            conToScaf[contig] = scaffold
            if scaffold in scafToContigs:
                scafToContigs[scaffold].append(contig)
            else:
                temp = []
                temp.append(contig)
                scafToContigs[scaffold] = temp
            #print str('"' + scaffold + '" "' + contig + '"')

        cSet = Set([])
        for contig in contigsList:
            cSet.add(contig)
            if contig in conToScaf:
                scaffold = conToScaf[contig]
                conList = scafToContigs[scaffold]
                for c in conList:
                    cSet.add(c)
        temp = []
        for c in cSet:
            temp.append(c)
    return temp


if __name__ == "__main__":
    inDir = 'D:/A_Phylo/A_Metagenomic/reindeer/training_data02/names3'
    outDir = 'D:/A_Phylo/A_Metagenomic/reindeer/training_data02/names3_f_tr'
    mapFastaFile = 'D:/A_Phylo/A_Metagenomic/reindeer/data/SRM_Large_Contigs_namesOnly.fna'
    #scaffoldContigFile = 'D:/A_Phylo/A_Metagenomic/reindeer/contigs-scaffolds.csv'
    scaffoldContigFile = None
    ignoreContigsList = 'D:/A_Phylo/A_Metagenomic/reindeer/data/SRM-1_contigs_names'
    fillInSeq(inDir, outDir, mapFastaFile, ignoreContigsList, scaffoldContigFile)