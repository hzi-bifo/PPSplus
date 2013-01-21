#!/usr/bin/env python

"""
    Defines functions that implement common operations with Fasta files.
"""

import sys
import os
import re
import types

from Bio import SeqIO

from TabSepFileFunctions import OutFileBuffer
from Common import removeNonDna
from Common import noNewLine


def filterOutSequences(inFileName, outFileName, allowedNamesSet, formatName="fasta", seqNameModifyFunction = None):
    """
        From the input fasta file filter out sequences their names are not contained in the allowedNamesSet.

        @param allowedNamesSet: the set of entries that are allowed as a sequence names
        @param seqNameModifyFunction: a sequence`s name is modified by this function and then compared to the allowedNamesSet
    """
    outFileBuffer = OutFileBuffer(outFileName)
    recordCondition = RecordConditionFilterOutSequences(allowedNamesSet, seqNameModifyFunction)
    parser = RecordFilter(outFileBuffer, formatName, recordCondition)
    _forEachRecord(inFileName, parser)

def filterOutNonDna(inFileName, outFileName):
    outFileBuffer = OutFileBuffer(outFileName)
    parser = RemoveNonDnaParser(outFileBuffer)
    _forEachRecord(inFileName, parser)


def getSequenceToBpDict(fastaFilePath):
    """
        Reads a fasta file and returns mapping: sequenceName -> sequenceLength.
    """
    return _forEachRecord(fastaFilePath, SeqToBpParser()).getSeqToBpDict()


class SeqToBpParser():
    def __init__(self):
        self._seqToBp = dict([])
    def parse(self, record):
        self._seqToBp[record.id] = len(str(record.seq))
    def getSeqToBpDict(self):
        return self._seqToBp
    def getFormatName(self):
        return "fasta"


def getSequencesToList(fastaFilePath):
    """
        Reads a fasta file and returns a list of: (sequenceName, sequence).
    """
    return _forEachRecord(fastaFilePath, SeqToListParser()).getSeqToList()


class SeqToListParser():
    def __init__(self):
        self._seqToList = []
    def parse(self, record):
        self._seqToList.append((str(record.id), noNewLine(str(record.seq))))
    def getSeqToList(self):
        return self._seqToList
    def getFormatName(self):
        return "fasta"


class RemoveNonDnaParser():
    def __init__(self, outFileBuffer):
        self._outFileBuffer = outFileBuffer

    def finalize(self):
        self._outFileBuffer.close()

    def getFormatName(self):
        return "fasta"

    def parse(self, record):
        self._outFileBuffer.writeText(str('>' + str(record.id) + '\n'))
        self._outFileBuffer.writeText(str(removeNonDna(str(record.seq)) + '\n'))


class RecordConditionFilterOutSequences():
    def __init__(self, allowedNamesSet, seqNameModifyFunction = None):
        self.allowedNamesSet = allowedNamesSet
        self.seqNameModifyFunction = seqNameModifyFunction

    def takeRecord(self, record):
        """
            If the record.id (modified by the function) is in the allowedNamesSet then the entry will be accepted.
        """
        id = record.id
        if self.seqNameModifyFunction != None:
            id = self.seqNameModifyFunction(id)
        if id in self.allowedNamesSet:
            return True
        else:
            return False


class RecordFilter():
    """
        Appends a record that is currently parsed to the outFileBuffer if it satisfies the condition.
    """
    def __init__(self, outFileBuffer, formatName, recordCondition):
        self.outFileBuffer = outFileBuffer
        self.formatName = formatName
        self.recordCondition = recordCondition

    def parse(self, record):
        if self.recordCondition.takeRecord(record):
            self.outFileBuffer.writeText(record.format(self.formatName))

    def getFormatName(self):
        return self.formatName

    def finalize(self):
        self.outFileBuffer.close()


def fastaFileToDict(fastaFilePath, formatName='fasta'):
    """
        Reads a fasta file and returns mapping: seqName -> sequence.
    """
    return _forEachRecord(fastaFilePath, _RecordStorage(formatName),formatName=formatName).getSeqNameToSeq()


def fastaFileToDictWholeNames(filePath):
    """
        Reads a fasta file and returns mapping: seqName -> sequence the whole sequence name is used
        as seqName!!! (even if it contains space)
    """
    seqIdToSeq = dict([])
    try:
        f = open(os.path.normpath(filePath),'r')
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
                name = line.replace('>','')
            else:
                seq += line
        if seq != '':
            assert name != ''
            seqIdToSeq[name] = seq
    finally:
        f.close()
    return seqIdToSeq


class _RecordStorage():
    def __init__(self, formatName='fasta'):
        self._seqNameToSeq = dict([])
        self._formatName = formatName

    def parse(self, record):
        self._seqNameToSeq[str(record.id)] = str(record.seq)

    def getFormatName(self):
        return self._formatName

    def getSeqNameToSeq(self):
        return self._seqNameToSeq


def _forEachRecord(filePath, parser, formatName = "fasta"):
    """
        Call the parser for each record in the file.
    """
    try:
        if isinstance(parser.getFormatName, types.MethodType):
            formatName = parser.getFormatName()
    except Exception:
        pass
    try:
        f = open(os.path.normpath(filePath), 'r')
    except Exception:
        sys.stderr.write('Cannot open a ' + formatName + ' file for reading: ' + filePath + '\n')
        raise
    else:
        try:
            readBuffer = SeqIO.parse(f, parser.getFormatName())
            for record in readBuffer:
                parser.parse(record)
        except Exception:
            sys.stderr.write('Cannot read from a ' + formatName + ' file: ' + filePath + '\n')
            raise
        finally:
            f.close()
    try:
        if isinstance(parser.finalize, types.MethodType):
            parser.finalize()
    except Exception:
        pass

    return parser


def main():
    pass

def test():
    #filterOutNonDna('/Users/ivan/Documents/work/binning/data/TW/TWexpertSSD_nonDna/538960.1.fas',
    #                '/Users/ivan/Documents/work/binning/data/TW/TWexpertSSD_nonDna/out/538960.1.fas');
    filterOutNonDna('/Users/ivan/Documents/work/binning/data/Reindeer/927658.1.fas',
        '/Users/ivan/Documents/work/binning/data/Reindeer/nonDnaRemoved/927658.1.fas')

if __name__ == "__main__":
    #main()
    test()