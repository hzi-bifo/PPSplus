#!/usr/bin/env python

#/AM/home-0/shared/python/Python-2.7.1/python

import os
import re
import subprocess
import Common
from Config import Config
from Config import Config2
from Sequences import Sequences
from Taxonomy import Taxonomy


class Johdro():
    def __init__(self, config, configJoh):
        self._config = config
        self._configJoh = configJoh

    #test whether the env. variables are set - test only
    def setEnvVarCmd(self):
        return str('export TAXATORTK_NCBI_ROOT=' + os.path.normpath(self._configJoh.get('johdroTaxdump')) + '; ')


    #run blast and johdro postfilter
    def runBlast(self, fastaFile):
        workingDir = self._config.get('workingDir')

        cmd = str('time ' + os.path.join(os.path.normpath(self._configJoh.get('blastInstallDir')), 'blastn')
                  + ' ' + self._configJoh.get('blastParam') + ' -query ' + fastaFile
                  + ' | python ' + os.path.join(os.path.normpath(self._configJoh.get('johdroScriptsDir')), 'blastxml2alignments.py') + ' --parse-refid ''no'' '
                  + ' | gzip > ' + Common.createTagFilePath(workingDir, fastaFile, 'bA') + '.gz'
                  )

        if os.name == 'posix':
            blastProc = subprocess.Popen(cmd, shell=True, bufsize=-1, cwd=workingDir) #stdout=subprocess.STDOUT, stderr=subprocess.STDOUT)
            print 'run cmd:', cmd
            blastProc.wait()
            print 'Blast return code:', blastProc.returncode
        else:
            print 'Cannot run Blast since your system is not "posix" but', str('"' + os.name + '"'), '\n', cmd


    def runLast(self, fastaFile):
        workingDir = self._config.get('workingDir')
        cmd = str('time ' + os.path.join(os.path.normpath(self._configJoh.get('lastInstallDir')), 'lastal')
                  + ' ' + self._configJoh.get('lastParam') + ' ' + fastaFile
                  + ' | python ' + os.path.join(os.path.normpath(self._configJoh.get('johdroScriptsDir')), 'lastmaf2alignments.py')
                  #+ ' | ' + os.path.join(os.path.normpath(self._configJoh.get('johdroBinDir')), 'extract-fastacomment-ncbifield') + ' -k gi '
                  + ' | sort '
                  + ' | gzip > ' + Common.createTagFilePath(workingDir, fastaFile, 'lA') + '.gz'
                  )
        if os.name == 'posix':
            lastProc = subprocess.Popen(cmd, shell=True, bufsize=-1, cwd=workingDir) #stdout=subprocess.STDOUT, stderr=subprocess.STDOUT)
            print 'run cmd:', cmd
            lastProc.wait()
            print 'Last return code:', lastProc.returncode
        else:
            print 'Cannot run Last since your system is not "posix" but', str('"' + os.name + '"'), '\n', cmd


    def runBlastLCA(self, fastaFile):
        workingDir = self._config.get('workingDir')
        cmd = str(self.setEnvVarCmd() + 'time zcat ' + Common.createTagFilePath(workingDir, fastaFile, 'bA') + '.gz'
                  + ' | ' + os.path.join(os.path.normpath(self._configJoh.get('johdroBinDir')), 'alignments-prefilter')
                  + ' ' + self._configJoh.get('prefilterParam')
                  + ' | ' + os.path.join(os.path.normpath(self._configJoh.get('johdroBinDir')), 'predictor')
                  + ' ' + self._configJoh.get('predictorParam')
                  + ' |  cat > ' + Common.createTagFilePath(workingDir, fastaFile, 'bP')
                  )
        if os.name == 'posix':
            blastProc = subprocess.Popen(cmd, shell=True, bufsize=-1, cwd=workingDir) #stdout=subprocess.STDOUT, stderr=subprocess.STDOUT)
            print 'run cmd:', cmd
            blastProc.wait()
            print 'Blast LCA return code:', blastProc.returncode
        else:
            print 'Cannot run Blast LCA since your system is not "posix" but', str('"' + os.name + '"'), '\n', cmd


    def runLastLCA(self, fastaFile):
        workingDir = self._config.get('workingDir')
        cmd = str(self.setEnvVarCmd() + 'time zcat ' + Common.createTagFilePath(workingDir, fastaFile, 'lA') + '.gz'
                  + ' | ' + os.path.join(os.path.normpath(self._configJoh.get('johdroBinDir')), 'alignments-prefilter')
                  + ' ' + self._configJoh.get('prefilterParam')
                  + ' | ' + os.path.join(os.path.normpath(self._configJoh.get('johdroBinDir')), 'predictor')
                  + ' ' + self._configJoh.get('predictorParam')
                  + ' | cat > ' + Common.createTagFilePath(workingDir, fastaFile, 'lP')
                  )
        if os.name == 'posix':
            lastProc = subprocess.Popen(cmd, shell=True, bufsize=-1, cwd=workingDir) #stdout=subprocess.STDOUT, stderr=subprocess.STDOUT)
            print 'run cmd:', cmd
            lastProc.wait()
            print 'Last LCA return code:', lastProc.returncode
        else:
            print 'Cannot run Last LCA since your system is not "posix" but', str('"' + os.name + '"'), '\n', cmd


    def setScaffCandidatePlacementFromBlastLCA(self, sequences, taxonomy, fastaFileScaff):
        workingDir = self._config.get('workingDir')
        johdroPlacementFile = Common.createTagFilePath(workingDir, fastaFileScaff, 'bP')
        return self._setScaffCandidatePlacementsToSequences(sequences, taxonomy, johdroPlacementFile)


    def setScaffCandidatePlacementFromLastLCA(self, sequences, taxonomy, fastaFileScaff):
        workingDir = self._config.get('workingDir')
        johdroPlacementFile = Common.createTagFilePath(workingDir, fastaFileScaff, 'lP')
        return self._setScaffCandidatePlacementsToSequences(sequences, taxonomy, johdroPlacementFile)


    def _setScaffCandidatePlacementsToSequences(self, sequences, taxonomy, johdroPlacementFile):
        #THE WEIGHT IS SET TO "None" SINCE IT IS NOT CONTAINED IN THE JOHDRO OUTPUT FILE !!!!!!!!!!!!
        count = 0
        try:
            f = open(os.path.normpath(johdroPlacementFile),'r')
        except Exception:
            print "Cannot open file:", johdroPlacementFile
            raise
        else:
            for line in f:
                line = Common.noNewLine(line)
                if re.match(r'^[0-9]+\t[0-9]+$', line):
                    scaffoldId = int(re.sub(r'^([0-9]+)\t[0-9]+$',r'\1' ,line))
                    ncbid = int(re.sub(r'^[0-9]+\t([0-9]+)$',r'\1' ,line))
                    if ncbid != 1:
                        taxPathDict = taxonomy.getPathToRoot(ncbid)
                        if taxPathDict.keys() >= 1:
                            count += sequences.setScaffCandidateTaxonomyPath(scaffoldId, taxPathDict, taxonomy, None) #weight ~ None
                        else:
                            raise Exception('No taxonomic path found for ncbid: ', ncbid)
        finally:
            f.close()
            return count


    def setCandidatePlacementFromBlastLCA(self, sequences, taxonomy, fastaFile):
        workingDir = self._config.get('workingDir')
        johdroPlacementFile = Common.createTagFilePath(workingDir, fastaFile, 'bP')
        return self._setCandidatePlacementsToSequences(sequences, taxonomy, johdroPlacementFile)


    def setCandidatePlacementFromLastLCA(self, sequences, taxonomy, fastaFile):
        workingDir = self._config.get('workingDir')
        johdroPlacementFile = Common.createTagFilePath(workingDir, fastaFile, 'lP')
        return self._setCandidatePlacementsToSequences(sequences, taxonomy, johdroPlacementFile)


    def _setCandidatePlacementsToSequences(self, sequences, taxonomy, johdroPlacementFile):
        #THE WEIGHT IS SET TO "None" SINCE IT IS NOT CONTAINED IN THE JOHDRO OUTPUT FILE !!!!!!!!!!!!
        count = 0
        try:
            f = open(os.path.normpath(johdroPlacementFile),'r')
        except Exception:
                print "Cannot open file:", johdroPlacementFile
                raise
        else:
            for line in f:
                line = Common.noNewLine(line)
                if re.match(r'^[0-9]+_[0-9]+\t[0-9]+$', line):
                    scaffoldId = int(re.sub(r'^([0-9]+)_([0-9]+)\t[0-9]+$',r'\1' ,line))
                    contigId = int(re.sub(r'^[0-9]+_([0-9]+)\t[0-9]+$',r'\1' ,line))
                    ncbid = int(re.sub(r'^[0-9]+_[0-9]+\t([0-9]+)$',r'\1' ,line))
                    if ncbid != 1:
                        #print line, ":", scaffoldId, contigId, ncbid
                        taxPathDict = taxonomy.getPathToRoot(ncbid)
                        if taxPathDict.keys() >= 1:
                            sequences.setCandidateTaxonomyPath(contigId, scaffoldId, taxPathDict, None) # weight ~ None
                            count += 1
                        else:
                            raise Exception('No taxonomic path found for ncbid: ', ncbid)
        finally:
            f.close()
        return count


def test():

    config = Config(open(os.path.normpath('D:\\A_Phylo\\A_Metagenomic\\pPPS\\workspace\\pPPS\\config01.cfg')), 'pPPS')
    configJoh = Config2(config, 'Johdro')
    johdro = Johdro(configJoh)
    #sequences = Sequences(config)
    #sequences.writeSequences(config.get('inputIdsFastaFile'))
    #taxonomy = Taxonomy(config.get('databaseFile'), config.get('taxonomicRanks').split(','))
    #placeSequences(sequences, taxonomy, configJoh.get('placementFile'))






if __name__ == "__main__":
  test()
