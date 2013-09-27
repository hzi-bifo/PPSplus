#!/usr/bin/env python

import os
import sys
import re
import subprocess

from com.csv import getColumnAsList
from com.csv import filterOutLines
from com.fasta import filterOutSequences
from com.config import Config
from com import common
from core.taxonomy import Taxonomy
from core.dna_to_prot import dnaToProt
from ref_seq.download_mg import protNameEntryToGid
from ref_seq.download_mg import dnaNameEntryToGid


class Amphora():
    def __init__(self, config, configAmphora):
        self._config = config
        self._configAmphora = configAmphora
        self._amphoraInstallDir = os.path.normpath(self._configAmphora.get('amphoraInstallDir'))
        self._workingDir = os.path.normpath(self._config.get('workingDir'))

    def predict(self, fastaFile):

        #transform the fasta file into a protein file
        fastaFileProt = str(fastaFile + '.prot')

        try:
            from Bio import SeqIO
        except:
            cmd = str('python ' + os.path.normpath('/AM/metagenomic/work/projects/pPPS/tools/dna_to_prot.py')
                      + ' ' + fastaFile + ' ' + fastaFileProt)
            reProc = subprocess.Popen(cmd, shell=True, bufsize=-1, cwd=self._workingDir)
            print 'run cmd:', cmd
            reProc.wait()
            print 'protein file creation return code:', reProc.returncode
        else:
            dnaToProt(fastaFile, fastaFileProt)

        #run Amphora
        amphoraOut = str(fastaFile + '.ah')

        phylotyping = os.path.normpath(os.path.join(self._amphoraInstallDir, 'Scripts', 'Phylotyping.pl'))

        cmd=str('time perl ' + phylotyping + ' ' + self._configAmphora.get('amphoraOptions') + ' ' + fastaFileProt + ' '
                + amphoraOut)

        if os.name == 'posix':
            amphoraProc = subprocess.Popen(cmd, shell=True, bufsize=-1, cwd=self._amphoraInstallDir)
            print 'run cmd:', cmd
            amphoraProc.wait()
            print 'Amphora return code:', amphoraProc.returncode
        else:
            print 'Cannot run Amphora since your system is not "posix" but', str('"' + os.name + '"'), '\n', cmd



    def setCandidatePlacementFromAmphora(self, sequences, taxonomy, fastaFile):
        assignedContigIdList = []
        amphoraOut = str(fastaFile + '.ah')
        amphoraOptions = self._configAmphora.get('amphoraOptions')
        try:
            f = open(os.path.normpath(amphoraOut),'r')
        except Exception:
            sys.stderr.write(str("Cannot open file:" + amphoraOut))
            raise
        else:
            for line in f:
                line = common.noNewLine(line)
                if re.match(r'^[0-9]+_[0-9]+_[pr012]+\t[^\t]+\t[^\t]+\tBootstrap:[0-9]*$', line):
                    scaffoldId = int(re.sub(r'^([0-9]+)_[0-9]+_[pr012]+\t[^\t]+\t[^\t]+\tBootstrap:[0-9]*$',r'\1' ,line))
                    contigId = int(re.sub(r'^[0-9]+_([0-9]+)_[pr012]+\t[^\t]+\t[^\t]+\tBootstrap:[0-9]*$',r'\1' ,line))
                    #assign->scientific name (sometimes bacteria)
                    #pred1 = re.sub(r'^[0-9]+_[0-9]+_[pr012]+\t([^\t]+)\t[^\t]+\tBootstrap:[0-9]*$',r'\1' ,line)
                    #outgroup->scientific name
                    pred2 = re.sub(r'^[0-9]+_[0-9]+_[pr012]+\t[^\t]+\t([^\t]+)\tBootstrap:[0-9]*$',r'\1' ,line)
                    bootstrap = re.sub(r'^[0-9]+_[0-9]+_[pr012]+\t[^\t]+\t[^\t]+\tBootstrap:([0-9]*)$',r'\1' ,line)

                    #ncbid = taxonomy.getNcbidFromScientificName(pred1)
                    ncbid = taxonomy.getNcbidFromScientificName(pred2)
                    if ncbid == None:
                        ncbid = 1

                    if bootstrap == '':
                        if re.match(r'.*-BootstrapCutoff[ ]+[0-9]+', amphoraOptions):
                            bootstrap = int(re.sub(r'.*-BootstrapCutoff[ ]+([0-9]+)',r'\1' , amphoraOptions))
                        else:
                            bootstrap = 70 #default
                    else:
                        bootstrap = int(bootstrap)
                    #the weights are considered to be at most 80
                    weight = float(bootstrap) * 80.0

                    if ncbid != 1:
                        taxPathDict = taxonomy.getPathToRoot(ncbid)
                        if taxPathDict.keys() >= 1:
                            sequences.setCandidateTaxonomyPath(contigId, scaffoldId, taxPathDict, weight)
                            assignedContigIdList.append(contigId)
                        else:
                            raise Exception('No taxonomic path found for ncbid: ', ncbid)


        finally:
            f.close()

        return len(set(assignedContigIdList))


#---------------------------

def filterOut5S():
    """
        Marker gene analysis - filter out reference sequences that hasn`t been found via corresponding HMM profiles for 5S
    """
    mgDir = os.path.normpath('/net/metagenomics/projects/PPSmg/data/nobackup/markerGenes2')
    print 'for 5S'
    domFileDna1 = os.path.join(mgDir,'mGenesValidated', '5S_arc_hmmsearch.dom')
    domFileDna2 = os.path.join(mgDir,'mGenesValidated', '5S_bac_hmmsearch.dom')
    set1 = set(getColumnAsList(domFileDna1, entryModifyFunction = dnaNameEntryToGid))
    set2 = set(getColumnAsList(domFileDna2, entryModifyFunction = dnaNameEntryToGid))
    allowedEntriesSet = set1 | set2

    inFileNameDna = os.path.join(mgDir,'mGenesExtracted','5S_bact+arch_dna.tax')
    inFileNameProt = os.path.join(mgDir,'mGenesExtracted','5S_bact+arch_prot.tax')
    outFileNameDna = os.path.join(mgDir,'db','5S_bact+arch_dnaV.tax')
    outFileNameProt = os.path.join(mgDir,'db','5S_bact+arch_protV.tax')
    filterOutLines(inFileNameDna, outFileNameDna, allowedEntriesSet, dnaNameEntryToGid)
    filterOutLines(inFileNameProt, outFileNameProt, allowedEntriesSet, protNameEntryToGid)

    inFileNameFastaDna = os.path.join(mgDir,'mGenesExtracted','5S_bact+arch_dna.noalign.fna')
    inFileNameFastaProt = os.path.join(mgDir,'mGenesExtracted','5S_bact+arch_prot.noalign.fna')
    outFileNameFastaDna = os.path.join(mgDir,'db','5S_bact+arch_dnaV.noalign.fna')
    outFileNameFastaProt = os.path.join(mgDir,'db','5S_bact+arch_protV.noalign.fna')
    filterOutSequences(inFileNameFastaDna, outFileNameFastaDna, allowedEntriesSet, "fasta", dnaNameEntryToGid)
    filterOutSequences(inFileNameFastaProt, outFileNameFastaProt, allowedEntriesSet, "fasta", protNameEntryToGid)


def filterOutAmphoraMG():
    """
        Marker gene analysis - as filterOut5S, but for all 31 marker genes
    """
    #for Amphora genes
    mgDir = os.path.normpath('/net/metagenomics/projects/PPSmg/data/nobackup/markerGenes2')

    geneNameArray = ['dnaG', 'infC', 'pgk', 'rpoB', 'tsf', 'frr', 'nusA', 'pyrG', 'rpmA', 'smpB',
                 'rpsC', 'rpsI', 'rpsK', 'rpsS', 'rpsB', 'rpsE', 'rpsJ', 'rpsM', 'rplA', 'rplB',
                 'rplC', 'rplD', 'rplE', 'rplF', 'rplK', 'rplL', 'rplM', 'rplN', 'rplP', 'rplS', 'rplT']

    for geneName in geneNameArray:
        print geneName
        domFileProt1 = os.path.join(mgDir,'mGenesValidated',str(geneName + '_hmmsearch.dom'))
        #domFileProt2 = str('D:/A_Phylo/A_Metagenomic/data/markerGenes/mGenesValidated/' + geneName + '_lshmmsearch.dom')
        set1 = set(getColumnAsList(domFileProt1, entryModifyFunction = protNameEntryToGid))
        #set2 = set(getColumnAsList(domFileProt2, entryModifyFunction = protNameEntryToGid))
        allowedEntriesSet = set1 #& set2
        #print str('gene:' + geneName + ' sw:' + str(len(set1)) + ' ls:' + str(len(set2)))

        inFileNameDna = os.path.join(mgDir,'mGenesExtracted',str(geneName + '_bact+arch_dna.tax'))
        inFileNameProt = os.path.join(mgDir,'mGenesExtracted',str(geneName + '_bact+arch_prot.tax'))
        outFileNameDna = os.path.join(mgDir,'db',str(geneName + '_bact+arch_dnaV.tax'))
        outFileNameProt = os.path.join(mgDir,'db',str(geneName + '_bact+arch_protV.tax'))
        filterOutLines(inFileNameDna, outFileNameDna, allowedEntriesSet, dnaNameEntryToGid)
        filterOutLines(inFileNameProt, outFileNameProt, allowedEntriesSet, protNameEntryToGid)

        inFileNameFastaDna = os.path.join(mgDir,'mGenesExtracted',str(geneName + '_bact+arch_dna.noalign.fna'))
        inFileNameFastaProt = os.path.join(mgDir,'mGenesExtracted',str(geneName + '_bact+arch_prot.noalign.fna'))
        outFileNameFastaDna = os.path.join(mgDir,'db',str(geneName + '_bact+arch_dnaV.noalign.fna'))
        outFileNameFastaProt = os.path.join(mgDir,'db',str(geneName + '_bact+arch_protV.noalign.fna'))
        filterOutSequences(inFileNameFastaDna, outFileNameFastaDna, allowedEntriesSet, "fasta", dnaNameEntryToGid)
        filterOutSequences(inFileNameFastaProt, outFileNameFastaProt, allowedEntriesSet, "fasta", protNameEntryToGid)

def test():
    #filterOut5S()
    filterOutAmphoraMG()

#-----------------

if __name__ == "__main__":
    test()
    #main()

def main():
    amphoraOut = os.path.normpath("D:\\A_Phylo\\A_Metagenomic\\pPPS\\workspace\\pPPS\\wdir02\\inputTW.fas.ids.ah")
    config = Config(open(os.path.normpath('D://A_Phylo//A_Metagenomic//pPPS//workspace//pPPS//config01.cfg')), 'pPPS')
    databaseFile = os.path.normpath(config.get('databaseFile'))
    taxonomicRanks = config.get('taxonomicRanks').split(',')
    taxonomy = Taxonomy(databaseFile, taxonomicRanks)
    #readAmphoraOutFile(amphoraOut, taxonomy) # commented out !!!