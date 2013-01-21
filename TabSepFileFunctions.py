#!/usr/bin/env python

#/AM/home-0/shared/python/Python-2.7.1/python


import os
import sys
import re
import types
from sets import Set
from Common import noNewLine
#from GeneDB import protNameEntryToGid
#from GeneDB import dnaNameEntryToGid

#from FastaFileFunctions import filterOutSequences



#CORE FUNCTIONS

#Returns a column of a file as a list.
#
#@colNum: starts with 0, all entries will be taken from this column (default: 0)
#@sep: the column separator (default: None ~ white space)
#@comment: lines that starts with this substring are not considered (default: #)
#
def getColumnAsList(fileName, entryModifyFunction=None, colNum=0, sep=None, comment='#'):
    lineParser = _ColumnEntryListBuffer(colNum, sep, comment, entryModifyFunction)
    return forEachLine(fileName, lineParser).retVal()


#From the input file filter out lines that does not contain entries from the allowedSet
#in a specific column and store the rest lines to the output file.
#
#@allowedEntriesSet: the set of entries that are allowed in the respective column
#@colNum: starts with 0, entries in this column will be considered (default: 0)
#@sep: the column separator (default: None ~ white space)
#@comment: lines that starts with this substring are just copied to the output file
#
def filterOutLines(inFileName, outFileName, allowedEntriesSet, entryModifyFunction=None, colNum=0, sep=None, comment='#'):
    outBuffer = OutFileBuffer(outFileName)
    lineCondition = _LineConditionFilterOutLines(allowedEntriesSet, colNum, sep, comment, entryModifyFunction)
    parser = _LineFilter(lineCondition, outBuffer)
    forEachLine(inFileName, parser)


#Transforms a tab separated file to a dictionary.
#
#@inFileName: tab separated file
#@keyColNum: the number of a column that represents keys (starting from 0)
#@valColNum: the number of a column that represent values (starting from 0)
#
def getMapping(inFileName, keyColNum, valColNum, sep=None, comment = '#'):
    return forEachLine(inFileName, _MappingParser(keyColNum, valColNum, sep, comment)).getDict()

class _MappingParser():
    def __init__(self, keyColNum, valColNum, sep, comment):
        self._dict = dict([])
        self._keyColNum = keyColNum
        self._valColNum = valColNum
        self._sep = sep
        self._comment = comment
    def getDict(self):
        return self._dict
    def parse(self, line):
        if not isComment(line, self._comment):
            lineList = line.split(self._sep)
            if len(lineList) > self._keyColNum and len(lineList) > self._valColNum:
                key = lineList[self._keyColNum]
                val = lineList[self._valColNum]
                if key not in self._dict:
                    tmp = []
                    tmp.append(val)
                    self._dict[key] = tmp
                else:
                    self._dict[key].append(val)
            else:
                if len(lineList) > 0:
                    print str('TabSepFileFunctions:_MappingParser: line skipped: ' + line + ' doesn`t have enough entries\n')




#SPECIFIC HELPER FUNCTIONS AND CLASSES DEFINITION

class _LineConditionFilterOutLines():
    def __init__(self, allowedEntriesSet, colNum, sep, comment, entryModifyFunction=None):
        self.allowedEntriesSet = allowedEntriesSet
        self.colNum = colNum
        self.sep = sep
        self.comment = comment
        self.entryModifyFunction = entryModifyFunction

    def takeLine(self, line):
        if isComment(line, self.comment):
            return True
        else:
            lineList = line.split(self.sep)
            if len(lineList) > self.colNum:
                entry = lineList[self.colNum]
                if self.entryModifyFunction != None:
                    entry = self.entryModifyFunction(entry)
                if entry in self.allowedEntriesSet:
                    return True
        return False

#GENERAL HELPER FUNCTIONS AND CLASSES

class _LineFilter():
    def __init__(self, lineCondition, outFileBuffer):
        self.lineCondition = lineCondition
        self.outFileBuffer = outFileBuffer

    def parse(self, line):
        if self.lineCondition.takeLine(line):
            self.outFileBuffer.writeText(str(line + '\n'))

    def finalize(self):
        self.outFileBuffer.close()

# reads predictions
def predToDict(predFilePath):
    return forEachLine(predFilePath, _PredParser()).getContigToPredDict()

#To parse the prediction file
class _PredParser():
    def __init__(self):
        self._dict = dict([])

    def parse(self, line):
        if not isComment(line, '#'):
            lineList = line.split()
            if len(lineList) >= 2:
                key = lineList[0]
                val = lineList[len(lineList) - 1]
                if key in self._dict:
                    sys.stderr.write('Consistency:PredParser: the contig "' + key + '" has already been assigned' )
                self._dict[key] = int(val)


    def getContigToPredDict(self):
        return self._dict



#for each line of the file call the parser, at the end call the finalize method of the parser if it`s defined
def forEachLine(filePath, parser):
    try:
        f = open(os.path.normpath(filePath), 'r')
    except Exception:
        sys.stderr.write('Cannot open a file for reading: ' + filePath)
        raise
    else:
        try:
            for line in f:
                parser.parse(noNewLine(line))
        except Exception:
            sys.stderr.write('Cannot read from file: ' + filePath)
            raise
        finally:
            f.close()
    try:
        if isinstance(parser.finalize, types.MethodType):
            parser.finalize()
    except Exception:
        pass

    return parser


#is the line a comment
def isComment(line, comment):
    if re.match(str('[ \t]*' + comment), line):
        return True
    else:
        return False


#Parse a line and collects an entry in the specific column
class _ColumnEntryListBuffer():
    def __init__(self, colNum, sep, comment, entryModifyFunction=None):
        self.list = []
        self.colNum = colNum
        self.sep = sep
        self.comment = comment
        self.entryModifyFunction = entryModifyFunction

    def parse(self, line):
        if not isComment(line, self.comment): #the line is not a comment
            lineList = line.split(self.sep)
            if len(lineList) > self.colNum:
                entry = lineList[self.colNum]
                if self.entryModifyFunction != None:
                    entry = self.entryModifyFunction(entry)
                self.list.append(entry)

    def retVal(self):
        return self.list


#to append text to a file
class OutFileBuffer():
    def __init__(self, outFilePath, bufferText = False, fileOpenMode='w'):
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
            if not self.opened: #reopen to append
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






def main():
    pass

if __name__ == "__main__":
    #main()
    test()