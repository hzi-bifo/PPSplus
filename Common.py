#!/usr/bin/env python

#/AM/home-0/shared/python/Python-2.7.1/python

import os
import re
from sets import Set

#replace stretches of nonDNA characters with one N
def removeNonDna(seq):
    return re.sub(r'[^ATGCatgc]+', 'N', seq).upper()

#delete all '\n' and '\r' characters in a string
def noNewLine(str):
    return str.replace('\n', '').replace('\r','')

#return a path that results from the concatenation of dstDir and a file from filePath where the name of the file
#changes from e.g.: input.n.fas" to "input.n._ids.fas"
def createTagFilePath(dstDir, fileNameFromPath, tag):
    lastDotIdx = fileNameFromPath.rfind('.')
    return os.path.join(os.path.normpath(dstDir),
                        os.path.basename(os.path.normpath(str(fileNameFromPath[0:lastDotIdx] +  fileNameFromPath[lastDotIdx:]
                                                            + '.' + tag))))


#returns true if two files contain the same sequences regardless of their names (and empty spaces)
def seqFileCmp(file1, file2):
    seqList1 = seqFileToSeqList(file1)
    seqList2 = seqFileToSeqList(file2)
    if (len(seqList1) != len(seqList2)):
        print "The files contain different number of sequences", len(seqList1), len(seqList2)
        return False
    else:
        print "Number of sequences: ", len(seqList1)

    seqFoundInS2 = Set([])
    seqAF = Set([])
    idx1 = 0
    for s1 in seqList1:
        idx1 += 1
        idx2 = 0
        for s2 in seqList2:
            idx2 += 1
            if (s1 == s2):
                #print "same:", idx1, idx2
                if s2 not in seqAF:
                    #print idx1, idx2
                #else:
                    seqAF.add(s2)
                if idx2 not in seqFoundInS2:
                    #print "One sequence is in one file more than once", idx1, idx2
                    seqFoundInS2.add(idx2)
                    continue

    s1 = Set([])
    s2 = Set([])
    for s in seqList1:
        s1.add(s)
    for s in seqList2:
        s2.add(s)
    if len(s1) != len(s2):
        print "The length of unique sequences differ S1:", len(s1), "S2:", len(s2)
    else:
        print "The number of unique sequences is: ", len(s1)


    if len(seqFoundInS2) == len(seqList1):
        print "Both files contain the same sequences"
        return True
    else:
        print "Sequences matches: ", len(seqList1), " Sequences found: ", len(seqFoundInS2)
        return False

def seqFileToSeqList(file):
    seqList = []
    try:
        f = open(os.path.normpath(file),'r')
    except Exception:
        print "Cannot open file:", file1
        raise
    else:
        name = ''
        seq = ''
        for line in f:
            line = noNewLine(line)
            if re.match('>', line):
                if seq != '':
                    assert name != ''
                    seqList.append(seq) #store seq
                    seq = ''
                name = line.replace('>','')
            else:
                seq += line
        if seq != '':
            assert name != ''
            seqList.append(seq) #store seq
        return seqList
    finally:
        f.close()



if __name__ == "__main__":
    #file = "200643.1.fna"
    file1="/Users/ivan/Documents/work/binning/database/silva/align/LSURef_106_tax_silva_trunc.fasta"
    file2="/Users/ivan/Documents/work/binning/database/silva/lsuparc_silva106_ncbitax.bacteria+archaea.fna"
    print seqFileCmp(file1, file2)
    #seqFileCmp(str("D:\\A_Phylo\\A_Metagenomic\\pPPS\\workspace\\pPPS\\toKaustubh\\variant01\\trainingData\\" + file),
    #           str("D:\\A_Phylo\\A_Metagenomic\\pPPS\\workspace\\pPPS\\toKaustubh\\variant02\\trainingData\\" + file))
    #filePath = "dir/input.n.fas"
    #dstDir = "D:/A_Phylo/A_Metagenomic/pPPS/workspace/pPPS/toKaustubh"
    #tag = "ids"
    #print createTagFilePath(dstDir, filePath, tag)