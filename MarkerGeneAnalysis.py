#!/usr/bin/env python

#/AM/home-0/shared/python/Python-2.7.1/python

import os
import sys
import re
import subprocess
import Common
from sets import Set
from Bio import SeqIO

from TabSepFileFunctions import forEachLine
from dnaToProt import dnaToProt
from TabSepFileFunctions import isComment
from TabSepFileFunctions import OutFileBuffer



#helper class to store all paths to the marker gene files
class _MgFiles():
    def __init__(self):
        self.mgDict = dict([])

    def parse(self, line):
        val = line.split('\t')
        if len(val) == 3:
            self.addFilePath(val[0], val[1], val[2])

    def addFilePath(self, geneName, fileType, filePath):
        if geneName not in self.mgDict:
            self.mgDict[geneName] = dict([])
        self.mgDict[geneName][fileType] = os.path.normpath(filePath)

    def getFilePath(self, geneName, fileType):
        if (geneName in self.mgDict) and (fileType in self.mgDict[geneName]):
            return self.mgDict[geneName][fileType]
        else:
            return None

    def getGeneNameList(self):
        resultList = list([])
        for name in self.mgDict:
            resultList.append(name)
        return resultList


#helper class to read the hmmer dom? file to get the regions that correspond to the Amphora marker genes
class _MgRegions():
    def __init__(self):
        self.entryDict = dict([])

    def parse(self, line):
        if not isComment(line, '#'):
            entryArray = line.split()
            if len(entryArray) == 23:
                aliFrom = int(entryArray[17])
                aliTo = int(entryArray[18])
                targetName = entryArray[0]
                #if targetName in self.entryDict:
                #    print str('MarkerGeneAnalysis: regions for the targetName were already set old: ' + str(self.entryDict[targetName])
                #    + ' new: ' + str([aliFrom, aliTo]))
                if targetName not in self.entryDict:
                    self.entryDict[targetName] = []
                self.entryDict[targetName].append([aliFrom, aliTo])
                #entry['target_name'] = entryArray[0]
                #entry['query_name'] = entryArray[3]
                #entry['dom_i_evalue'] = entryArray[12]
                #entry['ali_from'] = int(entryArray[17])
                #entry['ali_to'] = int(entryArray[18])
            else:
                print 'Line not considered: ', line

    def getEntryDict(self):
        return self.entryDict


#helper class to parse the mothur output file and creates a standardized output file for each marker gene
class _MothurOutFileParser():
    def __init__(self, outBuffer, source):
        self.outBuffer = outBuffer
        self.source = source

    def parse(self, line):
        lineArray = line.split()
        if len(lineArray) != 2:
            print '_MothurOutFileParser: wrong line', line
            return
        name = re.sub(r'^([0-9]+_[0-9]+)_[0-9]+_[0-9]+_[pr]+[0-2]$',r'\1', lineArray[0])
        tag = re.sub(r'^[0-9]+_[0-9]+_([0-9]+_[0-9]+_[pr]+[0-2])$',r'\1', lineArray[0])
        placementList = lineArray[1].replace('unclassified;', '').rsplit(';')
        if len(placementList) < 2:
            #print '_MothurOutFileParser: skip line', line
            return

        placement = placementList[-2]
        clade = int(re.sub('([0-9]+)\(.*', r'\1' , placement))
        weight = float(re.sub('[0-9]+\(([0-9\.]+)\)', r'\1' , placement))

        entry = str(str(name) + '\t' + str(clade) + '\t' + str(weight) + '\t' + str(self.source) + '\t' + str(tag))
        if self.outBuffer.isEmpty():
            self.outBuffer.writeText(entry)
        else:
            self.outBuffer.writeText(str('\n' + entry))

    def finalize(self):
        self.outBuffer.close()


#main class to perform the marker gene analysis based on the Amphora marker genes
class MarkerGeneAnalysis():

    def __init__(self, config, configMG, configRRNA16S, mgWorkingDir):
        self.markerGeneListFile = os.path.normpath(configMG.get('markerGeneListFile'))
        self.markerGeneWorkingDir = mgWorkingDir #os.path.normpath(configMG.get('markerGeneWorkingDir'))
        self.hmmInstallDir = os.path.normpath(configRRNA16S.get('rnaHmmInstallDir'))
        self.hmmerBinDir = os.path.normpath(configRRNA16S.get('hmmerBinDir'))
        self.mothurParam = configRRNA16S.get('mothurClassifyParamOther')
        self.workingDir = os.path.normpath(config.get('workingDir'))
        self.hmmerBinDir = os.path.normpath(configRRNA16S.get('hmmerBinDir'))
        self.mothur = os.path.join(os.path.normpath(configRRNA16S.get('mothurInstallDir')), 'mothur')


    #run hmmer HMM and mothur classify (bayesian), same param as for the 16S analysis
    def runMarkerGeneAnalysis(self, fastaFileDNA, outLog=None):
        #read list of marker genes
        mgFiles = forEachLine(self.markerGeneListFile, _MgFiles())

        #translate DNA to protein sequences
        fastaFileProt = os.path.join(self.markerGeneWorkingDir, str(os.path.basename(fastaFileDNA) + '.PROT'))
        dnaToProt(fastaFileDNA, fastaFileProt)

        #read DNA fasta file
        try:
            handle = open(fastaFileDNA, "rU")
            dnaSeqDict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
            handle.close()
        except Exception:
            sys.stderr.write(str('Cannot read file: ' + str(fastaFileDNA)))
            raise

        #to output all predictions in one file
        outPredAllFileName = os.path.join(self.markerGeneWorkingDir,
                                           str(os.path.basename(fastaFileDNA) + '_all.mP'))
        outAllBuffer = OutFileBuffer(outPredAllFileName)

        #run HMM search
        mgList = mgFiles.getGeneNameList()

        if outLog != None:
            stdoutLog = open(outLog,'w')
        else:
            stdoutLog = subprocess.STDOUT

        #for each gene perform the analysis separately
        for geneName in mgList:

            domFileArray = [os.path.join(self.markerGeneWorkingDir, str(geneName + '_1.dom')),
                            os.path.join(self.markerGeneWorkingDir, str(geneName + '_2.dom'))]
            outFileArray = [os.path.join(self.markerGeneWorkingDir, str(geneName + '_1.out')),
                            os.path.join(self.markerGeneWorkingDir, str(geneName + '_2.out'))]
            hmmFileArray = [mgFiles.getFilePath(geneName, 'hmmPROTPrim'),
                            mgFiles.getFilePath(geneName, 'hmmPROTSec')]
            cmdArray = list([])

            #define cmd
            for i in range(2):
                if hmmFileArray[i] != None:
                    cmdArray.append(str(os.path.join(self.hmmerBinDir, 'hmmsearch') + ' --domtblout ' + domFileArray[i] + ' -E 0.01'
                               + ' -o ' + outFileArray[i] + ' ' + hmmFileArray[i] + ' ' + fastaFileProt))
                else:
                    cmdArray.append(None)

            #run cmd
            for cmd in cmdArray:
                if cmd != None and os.name == 'posix':
                    hmmProc = subprocess.Popen(cmd, shell=True, bufsize=-1, cwd=self.hmmInstallDir, stdout=stdoutLog)
                    print 'run cmd:', cmd
                    hmmProc.wait()
                    print 'HMM  return code:', hmmProc.returncode
                else:
                    print 'Marker genes analysis, doesn`t run (no posix): ', cmd

            #get regions that match to the HMM profile ()
            entryDictList = []
            for i in range(2):
                if cmdArray[i] != None:
                    entryDictList.append(forEachLine(domFileArray[i], _MgRegions()).getEntryDict())
                else:
                    entryDictList.append(None)

            entryDict1 = entryDictList[0]
            entryDict2 = entryDictList[1]

            #extract regions found in the protein sequences that were found by the HMM and generate corresponding DNA sequences
            regionDnaFasta = os.path.join(self.markerGeneWorkingDir, str(geneName + '_dna.gff'))
            outFileBuffer = OutFileBuffer(regionDnaFasta)

            for seqName in entryDict1:
                i = -1
                for e in entryDict1[seqName]:
                    i += 1
                    from1 = entryDict1[seqName][i][0]
                    to1 = entryDict1[seqName][i][1]
                    assert ((from1 != None) and (to1 != None))
                    #compare the results found by the primary and secondary HMM profiles
                    if (entryDict2 != None) and (seqName in entryDict2):
                        if len(entryDict2[seqName]) >= (i+1):
                            from2 = entryDict2[seqName][i][0]
                            to2 = entryDict2[seqName][i][1]
                            #if from1 != from2 or to1 != to2:
                            #    print str('Different positions in' + seqName + ' from1:' + str(from1) + ' from2:' + str(from2)
                            #                + ' to1:' + str(to1) + ' to2:' + str(to2))

                    #extract regions from the DNA sequences (consider 3 ORF and reverse complements)

                    #name of the whole sequence
                    dnaSeqName = re.sub(r'([0-9]+_[0-9]+)_[pr]+[012]', r'\1', seqName)
                    #whole DNA sequence
                    dnaSeq = dnaSeqDict[dnaSeqName].seq

                    #reverse complement (contains "pr")
                    tagRev = 'p'
                    if re.match(r'[0-9]+_[0-9]+_pr[012]', seqName):
                        dnaSeq = dnaSeq.reverse_complement()
                        tagRev = 'pr'

                    #shift "0"
                    if re.match(r'[0-9]+_[0-9]+_[pr]+0', seqName):
                        tagFrom = ((from1 - 1)*3)
                        tagTo = (to1*3)
                        tagRev += '0'
                        dnaSeq = dnaSeq[tagFrom:tagTo]

                    #shift "1"
                    elif re.match(r'[0-9]+_[0-9]+_[pr]+1', seqName):
                        tagFrom = (((from1 - 1)*3) + 1)
                        tagTo = ((to1*3) + 1)
                        tagRev += '1'
                        dnaSeq = dnaSeq[tagFrom:tagTo]

                    #shift "2"
                    elif re.match(r'[0-9]+_[0-9]+_[pr]+2', seqName):
                        tagFrom = (((from1 - 1)*3) + 2)
                        tagTo = ((to1*3) + 2)
                        tagRev += '2'
                        dnaSeq = dnaSeq[tagFrom:tagTo]

                    #error
                    else:
                        sys.stderr.write('Wrong seq name: ' + seqName + ' \n')
                        dnaSeq = None

                    tag = str(str(tagFrom) + '_' + str(tagTo) + '_' + tagRev)
                    outFileBuffer.writeText(str('>' + dnaSeqName + '_' + tag + '\n' + dnaSeq + '\n'))

            outFileBuffer.close()

            #if no marker gene found
            if outFileBuffer.isEmpty():
                continue

            #run mothur classify (bayesian? the same as for the 16S analysis)
            templateFile = mgFiles.getFilePath(geneName, 'templateDNA')
            taxonomyFile = mgFiles.getFilePath(geneName, 'taxonomyDNA')
            assert ((templateFile != None) and (taxonomyFile != None))
            cmd = str('time ' + self.mothur + ' "#classify.seqs(fasta=' + regionDnaFasta + ', template=' + templateFile
                + ', taxonomy=' +  taxonomyFile + ', ' + self.mothurParam + ')"')
            if os.name == 'posix':
                mothurProc = subprocess.Popen(cmd, shell=True, bufsize=-1, cwd=self.markerGeneWorkingDir, stdout=stdoutLog)
                print 'run cmd:', cmd
                mothurProc.wait()
                print 'mothur return code:', mothurProc.returncode
            else:
                print 'Cannot run mothur since your system is not "posix" but', str('"' + os.name + '"'), '\n', cmd

            #transform the mothur output to a simple output (name, ncbid, weight)

            mothurPredFileName = os.path.join(self.markerGeneWorkingDir,
                                              str(geneName + '_dna.' + os.path.basename(taxonomyFile) + 'onomy'))  # taxonomy
            outPredFileName = os.path.join(self.markerGeneWorkingDir,
                                           str(os.path.basename(fastaFileDNA) + '_' + geneName + '.mP'))
            outBuffer = OutFileBuffer(outPredFileName, bufferText = True)
            forEachLine(mothurPredFileName, _MothurOutFileParser(outBuffer, geneName))

            if not outAllBuffer.isEmpty():
                outAllBuffer.writeText('\n')
            outAllBuffer.writeText(outBuffer.getTextBuffer())

        if outLog != None:
            stdoutLog.close()
        outAllBuffer.close()

    #set candidate placement according to the marker gene analysis !!!!!!!!!!!!!!!!!
    def setCandidatePlacement(self, sequences, taxonomy, fastaFileDNA):
        outPredAllFileName = os.path.join(self.markerGeneWorkingDir,
                                          str(os.path.basename(fastaFileDNA) + '_all.mP'))
        return forEachLine(outPredAllFileName, _SetCandidatePlacement(sequences, taxonomy)).getAssignedSeqCount()


#helper class to set the candidate placements
class _SetCandidatePlacement():
    def __init__(self, sequences, taxonomy):
        self.sequences = sequences
        self.taxonomy = taxonomy
        self.assignedIdList = []

    def parse(self, line):
        if re.match(r'^[0-9]+_[0-9]+\t[0-9]+\t[0-9\.]+\t[^\t]+\t[^\t]+$', line):
            scaffoldId = int(re.sub(r'^([0-9]+)_[0-9]+\t[0-9]+\t[0-9\.]+\t[^\t]+\t[^\t]+$',r'\1' ,line))
            contigId = int(re.sub(r'^[0-9]+_([0-9]+)\t[0-9]+\t[0-9\.]+\t[^\t]+\t[^\t]+$',r'\1' ,line))
            ncbid = int(re.sub(r'^[0-9]+_[0-9]+\t([0-9]+)\t[0-9\.]+\t[^\t]+\t[^\t]+$',r'\1' ,line))
            weight = float(re.sub(r'^[0-9]+_[0-9]+\t[0-9]+\t([0-9\.]+)\t[^\t]+\t[^\t]+$',r'\1' ,line))
            source = str(re.sub(r'^[0-9]+_[0-9]+\t[0-9]+\t[0-9\.]+\t([^\t]+)\t[^\t]+$',r'\1' ,line))
            tag = str(re.sub(r'^[0-9]+_[0-9]+\t[0-9]+\t[0-9\.]+\t[^\t]+\t([^\t]+)$',r'\1' ,line))

        if ncbid != 1:
            taxPathDict = self.taxonomy.getPathToRoot(ncbid)
            if taxPathDict != None and taxPathDict.keys() >= 1:
                self.sequences.setCandidateTaxonomyPath(contigId, scaffoldId, taxPathDict, weight, source, tag)
                self.assignedIdList.append(contigId)
            else:
                sys.stderr.write(str('No taxonomic path found for ncbid: ' + str(ncbid)))

    def getAssignedSeqCount(self):
        return len(Set(self.assignedIdList))



def test():
    print 'Marker gene analysis'

if __name__ == "__main__":
    test()