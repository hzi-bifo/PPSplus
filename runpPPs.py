#!/usr/bin/env python

#/AM/home-0/shared/python/Python-2.7.1/python

import sys
import os
import shutil
import argparse
import re
import subprocess
import datetime
from Config import Config
from MGCluster import MGCluster
from Config import Config2
from Sequences import Sequences
from Taxonomy import Taxonomy
from TrainingData import PPSInput
from Johdro import Johdro
import Common
import MlTreeMap
import PPS
import SSDCrossVal
import TrainingData
import TabSepFileFunctions
import OutProc
from MarkerGeneAnalysis import MarkerGeneAnalysis
from Analysis16S import RRNA16S
from Amphora import Amphora
from Consistency import Consistency
from TabSepFileFunctions import OutFileBuffer

def main():
    parser = argparse.ArgumentParser(
                                     description='''Pipeline that pre-processes the input for PhyloPythiaS and run it.''',
                                     epilog='''''')#Read "readme.txt" for more details.

    parser.add_argument('-c', '--config', nargs=1, type=file, required=True,
                        help='configuration file of the pipeline', metavar='config.cfg',
                        dest='config')

    #parser.add_argument('-f', '--generate-working-fasta-files', action='store_true',
    #                    help='generates working fasta files (contigs/scaffolds are numbered)',
    #                    dest='f')

    #parser.add_argument('-m', '--run-mltreemap', action='store_true',
    #                    help='runs MlTreeMap',
    #                    dest='m')

    parser.add_argument('-n', '--run-rrna16S', action='store_true',
                        help='run hidden markov model and classify according to the 16S, 23S, and 5S genes',
                        dest='n')

    parser.add_argument('-g', '--run-marker-gene-analysis', action='store_true',
                        help='run hidden markov model and classify according to the marker genes (from Amphora)',
                        dest='g')

    #parser.add_argument('-a', '--run-amphora', action='store_true',
    #                    help='runs Amphora',
    #                    dest='a')

    parser.add_argument('-j', '--run-johdro', action='store', nargs='+',
                        help='run blast (b), run last (l), run LCA for blast (lcab), run LCA for last (lcal)',
                        dest='j')

    parser.add_argument('-o', '--process-preprocessed-output', action='store', nargs='+',
                        help=str('process output of Taxator Blast LCA (jlcab), Taxator Last LCA (jlcal), '
                                 + ' 16S rRNA analysis (s16), marker gene Amphora analysis without (mg)' +
                                 ' push down predictions to more specific clades (sc) ' +
                                 'build new OTUs (otu) '),
                        dest='o') #MlTreeMap (m); Amphora (ah); exclude some data from SSD (ex) place sequences;
                        #(exSSD) exclude defined contigs from SSD; create training data'

    parser.add_argument('-t', '--pps-train', action='store_true',
                        help='run PhyloPythiaS "train script"',
                        dest='t')

    parser.add_argument('-p', '--pps-predict', action='store', nargs='+',
                        help='run PhyloPythiaS "predict script" (c) for contigs, (s) for scaffolds, (v) run cross validation',
                        dest='p')

    parser.add_argument('-r', '--read-pps-out', action='store_true',
                        help='Reads the output/placements of PPS',
                        dest='r')

    parser.add_argument('-s', '--summary', action='store_true',
                        help='Summary',
                        dest='s')


    args = parser.parse_args()

    #read pPPS configuration
    #print args
    config = Config(args.config[0], 'pPPS')
    #configMl = Config2(config, 'MLTreeMap')
    configJoh = Config2(config, 'Johdro')
    configPPS = Config2(config, 'PPS')
    configRRNA16S = Config2(config, 'RRNA16S')
    #configAmph = Config2(config, 'Amphora')
    configMG = Config2(config, 'MarkerGenes')
    sequences = None

    workingDir = os.path.normpath(config.get('workingDir'))
    mgWorkingDir = os.path.join(workingDir,'mgWorking')
    outputDir = os.path.normpath(config.get('outputDir'))
    inputFastaFile = os.path.normpath(config.get('inputFastaFile'))

    #create the following directories in the working directory if they don`t exist
    dirArray = ['projectDir','sampleSpecificDir', 'mgWorking']
    for dirName in dirArray:
        dirPath = os.path.join(workingDir, dirName)
        if not os.path.isdir(dirPath):
            os.mkdir(dirPath)

    inputFastaScaffoldsFile = os.path.normpath(config.get('inputFastaScaffoldsFile'))
    if not os.path.isfile(inputFastaScaffoldsFile):
        inputFastaScaffoldsFile = None

    scaffoldsToContigsMapFile = os.path.normpath(config.get('scaffoldsToContigsMapFile'))
    if not os.path.isfile(scaffoldsToContigsMapFile):
        scaffoldsToContigsMapFile = None

    fastaFileIds = Common.createTagFilePath(workingDir , inputFastaFile, 'ids')

    if inputFastaScaffoldsFile != None:
        fastaFileScaffoldsIds = Common.createTagFilePath(workingDir , inputFastaScaffoldsFile, 'ids')
    else:
        fastaFileScaffoldsIds = None

    seqNameSeqIdFile = Common.createTagFilePath(workingDir, inputFastaFile, 'cToIds')
    if inputFastaScaffoldsFile != None:
        scaffNameScaffIdFile = Common.createTagFilePath(workingDir, inputFastaScaffoldsFile , 'sToIds')
    else:
        scaffNameScaffIdFile = None
    scaffoldContigMapIdsFile = Common.createTagFilePath(workingDir, inputFastaFile, 'mapSCIds')

    summaryAllFile = os.path.join(outputDir, 'summaryAll.txt')
    summaryTrainFile = os.path.join(outputDir, 'summaryTrain.txt')
    taxonomicRanks = config.get('taxonomicRanks').split(',')
    databaseFile = os.path.normpath(config.get('databaseFile'))

    referencePlacementFilePPOut = os.path.normpath(config.get('referencePlacementFilePPOut'))

    referencePlacementFileOut = os.path.normpath(config.get('referencePlacementFileOut'))

    useScaffInSSDPred = eval(config.get('useScaffInSSDPred'))

    #the command will be executed in Shell in Unix/Linux
    if os.name == 'posix':
        inShell=True
    else:
        inShell=False

    #generates working fasta files always
    #if args.f:
    sequences = Sequences(config, inputFastaFile, inputFastaScaffoldsFile, scaffoldsToContigsMapFile)
    sequences.writeSequences(fastaFileIds)
    print 'Working contigs input fasta file created:', fastaFileIds
    sequences.writeSeqNameSeqId(seqNameSeqIdFile)
    print 'Ids mapping for the working contigs fasta file created:', seqNameSeqIdFile
    sequences.writeScaffoldContigMap(scaffoldContigMapIdsFile)
    print 'Scaffolds -> contigs map ids file created:', scaffoldContigMapIdsFile
    #assert Common.seqFileCmp(inputFastaFile, fastaFileIds), 'The fasta IDs file contains different sequences!'
    if fastaFileScaffoldsIds != None:
        sequences.writeScaffolds(fastaFileScaffoldsIds)
        print 'Working scaffolds input fasta file created:', fastaFileScaffoldsIds #(scaffold id -> sequence)
        sequences.writeScaffNameScaffId(scaffNameScaffIdFile)
        print 'Ids mapping for the working scaffolds fasta file created:', scaffNameScaffIdFile

    #logging
    logDir = os.path.join(outputDir,'log')
    if not os.path.exists(logDir):
        os.mkdir(logDir)

    #getLogFileName(logDir, description)


    #is it specified what to do?
    if  not args.n and not args.g and not args.j and not args.o and not args.t and not args.p and not args.r and not args.s:
        print 'Choose what do you want to do!'
        print parser.print_help()


    #run MlTreeMap
    #if args.m:
    #    mlCmd = str('perl ' + os.path.normpath(os.path.join(configMl.get('mlTreeMapInstallDir'), 'mltreemap.pl'))
    #                    + ' -i ' + fastaFileIds
    #                    + ' -o ' + os.path.join(workingDir,'mlTreeMapOut')
    #                    + ' ' + configMl.get('arguments'))
    #    if os.name == 'posix':
    #        mlProc = subprocess.Popen(mlCmd, shell=True, bufsize=-1, cwd=configMl.get('mlTreeMapInstallDir')) #stdout=subprocess.STDOUT, stderr=subprocess.STDOUT)
    #        mlProc.wait()
    #        print 'MLTreeMap return code:', mlProc.returncode
    #    else:
    #        print 'Cannot run MLTreeMap since your system is not "posix" but', str('"' + os.name + '"'), '\n', mlCmd

    # run 16S analysis
    if args.n:
        rrna = RRNA16S(config, configRRNA16S)
        #if 'hmm' in args.n: #run hidden markov model
        print 'run Hidden Markov Model (16S)'
        rrna.runHMM(fastaFileIds, outLog=getLogFileName(logDir, 'hmm16S'))
        #if 's16' in args.n:
        print 'run 16S classification'
        rrna.classify16S(fastaFileIds, outLog=getLogFileName(logDir, 'mothurClassify16S'))
        print 'run 23S classification'
        rrna.classify23S(fastaFileIds, outLog=getLogFileName(logDir, 'mothurClassify23S'))
        print 'run 5S classification'
        rrna.classify5S(fastaFileIds, outLog=getLogFileName(logDir, 'mothurClassify5S'))

    #marker gene analysis (Amphora)
    if args.g:
        mg = MarkerGeneAnalysis(config, configMG, configRRNA16S, mgWorkingDir)
        print 'run (Amphora) marker gene analysis'
        mg.runMarkerGeneAnalysis(fastaFileIds, outLog=getLogFileName(logDir, 'amphoraMG'))


    # run Amphora
    #if args.a:
    #    amphora = Amphora(config, configAmph)
    #    print "run Amphora predict"
    #    amphora.predict(fastaFileIds)


    #run johdro (blast|last + predict)
    if args.j:
        johdro = Johdro(config, configJoh)
        if 'b' in args.j: #run blast
            print 'run Blast for contigs'
            johdro.runBlast(fastaFileIds)
            if fastaFileScaffoldsIds != None:
                print 'run Blast for scaffolds'
                johdro.runBlast(fastaFileScaffoldsIds)

        if 'l' in args.j: #run last
            print 'run Last for contigs'
            johdro.runLast(fastaFileIds)
            if fastaFileScaffoldsIds != None:
                print 'run Last for scaffolds'
                johdro.runLast(fastaFileScaffoldsIds)

        if 'lcab' in args.j: #run LCA blast
            print 'run LCA Blast for contigs'
            johdro.runBlastLCA(fastaFileIds)
            if fastaFileScaffoldsIds != None and useScaffInSSDPred:
                print 'run LCA Blast for scaffolds'
                johdro.runBlastLCA(fastaFileScaffoldsIds)

        if 'lcal' in args.j: #run LCA last
            print 'run LCA Last for contigs'
            johdro.runLastLCA(fastaFileIds)
            if fastaFileScaffoldsIds != None and useScaffInSSDPred:
                print 'run LCA Last for scaffolds'
                johdro.runLastLCA(fastaFileScaffoldsIds)


    #remove non-DNA characters in the sequences of contigs and scaffolds
    #this is important before the sample specific data are created and before PPS predicts the data
    sequences.setRemoveNonDna(True)
    #the working fasta files of scaffolds and contigs are regenerated, so that they wouldn`t contain non-DNA characters
    sequences.writeSequences(fastaFileIds)
    print 'Working contigs input fasta file without non-DNA chars created:', fastaFileIds
    if fastaFileScaffoldsIds != None:
        sequences.writeScaffolds(fastaFileScaffoldsIds)
        print 'Working scaffolds input fasta file without non-DNA chars created:', fastaFileScaffoldsIds #(scaffold id -> sequence)


    #make predictions files (run blast LCA, last LCA, MLTreeMap) result: contig/scaffold name \t ncbid

    #process MLTreeMap or johdro output and create file with ncbids and training data
    if args.o:
        if sequences == None:
            sequences = Sequences(config, inputFastaFile, inputFastaScaffoldsFile, scaffoldsToContigsMapFile)

        taxonomy = Taxonomy(databaseFile, taxonomicRanks)

        #place sequences from MLTreeMap #need to be rewritten!!!! place sequences only as candidates!!!
        #if 'm' in args.o:
        #    print 'process MLTreeMap output'
        #    MlTreeMap.placeSequences(sequences, taxonomy, configMl.get('final_RAxML_outputs')) #rewrite! only set candidates!
        #    print 'Placed sequences:', len(sequences.placedSeqSet) #revrite!

        #place sequences from johdro Blast LCA (as candidate placements)
        if 'jlcab' in args.o:
            johdro = Johdro(config, configJoh)
            print 'process johdro Blast LCA output'
            count = johdro.setCandidatePlacementFromBlastLCA(sequences, taxonomy, fastaFileIds)
            print 'Placed sequences by johdro Blast LCA:', count
            if fastaFileScaffoldsIds != None and useScaffInSSDPred:
                count = johdro.setScaffCandidatePlacementFromBlastLCA(sequences, taxonomy, fastaFileScaffoldsIds)
                print 'Placed sequences by johdro Blast LCA based on scaffolds:', count

        #place sequences from johdro Last LCA (as candidate placements)
        if 'jlcal' in args.o:
            johdro = Johdro(config, configJoh)
            print 'process johdro Last LCA output'
            count = johdro.setCandidatePlacementFromLastLCA(sequences, taxonomy, fastaFileIds)
            print 'Placed sequences by johdro Last LCA:', count
            if fastaFileScaffoldsIds != None and useScaffInSSDPred:
                count = johdro.setScaffCandidatePlacementFromLastLCA(sequences, taxonomy, fastaFileScaffoldsIds)
                print 'Placed sequences by johdro Last LCA based on scaffolds:', count

        #place sequences according to the 16S analysis
        if 's16' in args.o:
            rrna = RRNA16S(config, configRRNA16S)
            print 'process the 16S, 23S, and 5S output'
            countList = rrna.setCandidatePlacementFrom16S23S5S(sequences, taxonomy, fastaFileIds)
            print str('Placed sequences by the\n16S:' + str(countList[0]) + '\n23S:' + str(countList[1]) + '\n5S:' + str(countList[2]) + '\nall 16S, 23S, and 5S:' + str(countList[3]))

        if 'mg' in args.o:
             mg = MarkerGeneAnalysis(config, configMG, configRRNA16S, mgWorkingDir)
             print 'process marker gene analysis (Amphora) output'
             num = mg.setCandidatePlacement(sequences, taxonomy, fastaFileIds)
             print 'Placed sequences by the marker gene analysis (Amphora):', num


        #if 'ah' in args.o:
        #    amphora = Amphora(config, configAmph)
        #    print 'process Amphora output'
        #    count = amphora.setCandidatePlacementFromAmphora(sequences, taxonomy, fastaFileIds)
        #    print 'Placed sequences by Amphora:', count


        #place sequences according to the candidate placements
        sequences.setTaxonomyPathsFromCandidatePaths(taxonomy, float(config.get('candidatePlTopPercentThreshold')))
        print 'Sequences placed based on candidate placements, placed:', len(sequences.placedSeqSet)

        #assign not placed contigs of one scaffold to the lowest common ancestor of assigned contigs
        if eval(config.get('placeContigsFromTheSameScaffold')):
            sequences.placeContigsFromTheSameScaffold(taxonomy, float(config.get('agThreshold')),
                                                      float(config.get('assignedPartThreshold')),
                                                      float(config.get('candidatePlTopPercentThreshold')))
            print 'Not assigned contigs of one scaffold are placed to the lowest common ancestor of assigned contigs, placed:', len(sequences.placedSeqSet)
        else:
            print 'Not assigned contigs from the same scaffold are NOT placed'

        #marker gene clustering
        clust = MGCluster(config, configRRNA16S, configMG, mgWorkingDir, fastaFileIds, sequences,
                          taxonomy, os.path.basename(inputFastaFile))
        if 'sc' in args.o: #construct specific predictions
            clust.refineSpecificPred()

        if 'otu' in args.o: #build new otus
            clust.reconstructOTU()

        taxonomy.close()

        fFilePath = Common.createTagFilePath(workingDir, fastaFileIds, 'ssd_forbidden')
        if ('ex' in args.o) and (os.path.isfile(fFilePath)):
            forbiddenDict = TrainingData.loadDictFromAFile(fFilePath)
            print 'read file of forbidden entries'
        else:
            forbiddenDict = None

        exSSDFilePath = os.path.join(workingDir, 'excludeContigsAsSSD.txt')
        if ('exSSD' in args.o) and ((os.path.isfile(exSSDFilePath))):
            exSSDContigNameList = TabSepFileFunctions.getColumnAsList(exSSDFilePath, entryModifyFunction=None, colNum=0,
                                                                      sep='\n', comment='#')
            print 'read file of contigs to be excluded from SSD'
        else:
            exSSDContigNameList = None

        #create training data
        pps = PPSInput(sequences, taxonomicRanks, forbiddenDict, exSSDContigNameList, summaryAllFile)

        pps.createPPSInputFiles(os.path.join(workingDir, 'ncbids.txt'),
                                os.path.join(workingDir, 'sampleSpecificDir'),
                        int(config.get('rankIdAll')), int(config.get('rankIdCut')), int(config.get('rankIdCutMinBp')),
                        float(config.get('minPercentInLeaf')),
                        int(config.get('maxLeafClades')),
                        int(config.get('minBpToModel')),
                        int(config.get('minGenomesWgs')),
                        os.path.normpath(config.get('genomesWGS')),
                        None, #forbiddenDict
                        databaseFile,
                        taxonomicRanks,
                        int(config.get('fastaLineMaxChar')),
                        int(config.get('minSSDfileSize')), int(config.get('maxSSDfileSize')), float(config.get('weightStayAll')),
                        summaryTrainFile)

        #

    #run PhyloPythiaS train
    if args.t:
        trainCmd = str(os.path.normpath(os.path.join(configPPS.get('ppsInstallDir'), 'scripts', 'train.rb'))
                        + ' ' + os.path.normpath(configPPS.get('configPPS')))
        if os.name == 'posix':
            logOut = open(getLogFileName(logDir, 'PPS_train'),'w')
            trainProc = subprocess.Popen(trainCmd, shell=True, bufsize=-1, cwd=configPPS.get('ppsInstallDir'), stdout=logOut) #subprocess.STDOUT, stderr=subprocess.STDOUT)
            trainProc.wait()
            logOut.close()

            print 'PPS "train" return code:', trainProc.returncode
        else:
            print 'Cannot run PhyloPythiaS train since your system is not "posix" but', str('"' + os.name + '"'), 'CMD:\n', trainCmd


    #run PhyloPythiaS predict
    if args.p:
        if 'c' in args.p:
            predictCmd = str(os.path.normpath(os.path.join(configPPS.get('ppsInstallDir'), 'scripts', 'predict.rb'))
                            + ' ' + fastaFileIds
                            + ' ' + os.path.normpath(configPPS.get('configPPS')))

            if os.name == 'posix':
                logOut = open(getLogFileName(logDir, 'PPS_predict_c'),'w')
                predictProc = subprocess.Popen(predictCmd, shell=True, bufsize=-1, cwd=configPPS.get('ppsInstallDir'), stdout=logOut) #stdout=subprocess.STDOUT, stderr=subprocess.STDOUT)
                predictProc.wait()
                logOut.close()
                shutil.copy2(str(fastaFileIds + '.nox.fna.out'), str(fastaFileIds + '.out'))
                shutil.copy2(str(fastaFileIds + '.nox.fna.PP.out'), str(fastaFileIds + '.PP.out'))

                print 'PPS "predict" contigs return code:', predictProc.returncode
            else:
                print 'Cannot run PhyloPythiaS predict contigs since your system is not "posix" but', str('"' + os.name + '"'), 'CMD:\n', predictCmd

        if inputFastaScaffoldsFile != None:
            scaffoldsIdsNewPath = os.path.normpath(os.path.join(workingDir,'crossVal', os.path.basename(fastaFileScaffoldsIds)))

        if 's' in args.p:
            if inputFastaScaffoldsFile == None:
                print 'Cannot predict a scaffold file, since it is not specified!'
            else:
                crossValDir = os.path.normpath(os.path.join(workingDir,'crossVal'))
                try:
                    os.mkdir(os.path.normpath(crossValDir))
                    print 'Directory created:', crossValDir
                except OSError:
                    print 'Directory already exists:', crossValDir

                src = fastaFileScaffoldsIds
                dst = os.path.join(crossValDir, os.path.basename(fastaFileScaffoldsIds))
                try:
                    os.symlink(src, dst)
                except:
                    sys.stderr.write(str('Can`t create symbolic link: ' + dst + '\n'))
                    shutil.copyfile(src, dst)
                    sys.stderr.write(str('Copy of the file created: ' + dst + '\n'))
                #str('ln -s ' + fastaFileScaffoldsIds + ' ' + crossValDir + ';')
                predictCmd = str(os.path.normpath(os.path.join(configPPS.get('ppsInstallDir'), 'scripts', 'predict.rb'))
                            + ' ' +  scaffoldsIdsNewPath
                            + ' ' + os.path.normpath(configPPS.get('configPPS')))
                if os.name == 'posix':
                    logOut = open(getLogFileName(logDir, 'PPS_predict_s'),'w')
                    predictProc = subprocess.Popen(predictCmd, shell=True, bufsize=-1, cwd=configPPS.get('ppsInstallDir'), stdout=logOut) #stdout=subprocess.STDOUT, stderr=subprocess.STDOUT)
                    predictProc.wait()
                    logOut.close()
                    shutil.copy2(str(scaffoldsIdsNewPath + '.nox.fna.out'), str(scaffoldsIdsNewPath + '.out'))
                    shutil.copy2(str(scaffoldsIdsNewPath + '.nox.fna.PP.out'), str(scaffoldsIdsNewPath + '.PP.out'))
                    print 'PPS "predict" scaffolds return code:', predictProc.returncode
                else:
                    print 'Cannot run PhyloPythiaS predict scaffolds since your system is not "posix" but', str('"' + os.name + '"'), 'CMD:\n', predictCmd
        if 'v' in args.p:
            contigPred = os.path.normpath(str(fastaFileIds + '.PP.out'))
            scaffPred =  os.path.normpath(str(scaffoldsIdsNewPath + '.PP.out'))
            if not os.path.isfile(scaffPred):
                print 'Cross-validation: scaffolds were not predicted!'
            else:
                if not os.path.isfile(contigPred):
                    print 'Cross-validation: contigs were not predicted!', contigPred
                else:
                    #create prediction file for contigs from a prediction file for scaffolds
                    scaffPredAsContigs = os.path.normpath(os.path.join(workingDir,'crossVal', 'scaffPredAsContigs.PP.out'))
                    print 'scaff to contigs:', scaffoldsToContigsMapFile, scaffPred, scaffPredAsContigs
                    OutProc.scafToContigOutput(scaffoldContigMapIdsFile, scaffPred, scaffPredAsContigs)
                    i = 0
                    for rank in taxonomicRanks:
                        if i > int(config.get('tablesGeneratorMaxRank')):
                            break
                        i += 1
                        #do cross validation: contigs` prediction vs scaffolds` prediction
                        mlCmd = str('java -jar tables.jar ' + fastaFileIds + ' ' +
                                    contigPred + ' ' +
                                    scaffPredAsContigs + ' ' +
                                    os.path.normpath(os.path.join(os.path.dirname(scaffPredAsContigs),str('crossVal_' + str(rank) + '.csv'))) + ' ' +
                                    str(rank) + ' ' + 'unassigned' + ' ' + str(config.get('taxonomicRanks'))
                                    )
                        mlProc = subprocess.Popen(mlCmd, shell=inShell, bufsize=-1, cwd=config.get('tablesGeneratorDir')) #stdout=subprocess.STDOUT, stderr=subprocess.STDOUT)
                        mlProc.wait()
                        print 'Cross-validation:', mlCmd, 'return code: ', mlProc.returncode

    #read output from PPS and place sequences
    if args.r:
        if sequences == None:
            sequences = Sequences(config, inputFastaFile, inputFastaScaffoldsFile, scaffoldsToContigsMapFile)

        taxonomy = Taxonomy(databaseFile, taxonomicRanks)

        PPS.readPPSOutput(sequences, taxonomy, fastaFileIds, eval(config.get('overrideWithPPSPlacements')))


        #discard inconsistent contigs or scaffolds
        #discardedContigs = sequences.discardInconsistentPlacements(float(config.get('discardContigThreshold')),
        #                                                           float(config.get('discardScaffThreshold')),
        #                                                           taxonomicRanks)
        #print "The number of discarded contigs while reading PPS output:", discardedContigs

        taxonomy.close()
        print 'Placed sequences after PPS:', len(sequences.placedSeqSet)




    if args.s:
        #create PP.pOUT file
        sequences.writePlacementsPPOut(str(Common.createTagFilePath(outputDir, inputFastaFile, 'PP.pOUT')),
                                  taxonomicRanks,
                                  config.get('outputFileContigSubPattern'))

        #create .pOUT file (name tab ncbid)
        sequences.writePlacementsOut(str(Common.createTagFilePath(outputDir, inputFastaFile, 'pOUT')),
                                  taxonomicRanks,
                                  config.get('outputFileContigSubPattern'))

        sequences.writeSequences(str(Common.createTagFilePath(workingDir, inputFastaFile, 'pOUT.fas')), writeIds=False,
                                 outputFileContigSubPattern=config.get('outputFileContigSubPattern'))

        i = 0

        if os.path.isfile(referencePlacementFilePPOut):
            for rank in taxonomicRanks:
                if i > int(config.get('tablesGeneratorMaxRank')):
                    break
                i += 1
                mlCmd = str('java -jar tables.jar ' + str(Common.createTagFilePath(workingDir, inputFastaFile, 'pOUT.fas')) + ' ' +
                             str(referencePlacementFilePPOut) + ' ' +
                             str(Common.createTagFilePath(outputDir, inputFastaFile, 'PP.pOUT')) + ' ' +
                             str(Common.createTagFilePath(outputDir, inputFastaFile, str('table.' + str(rank) + '.csv')))  + ' ' +
                             str(rank) + ' ' +
                             'unassigned' + ' ' +
                             str(config.get('taxonomicRanks'))

                            )
                mlProc = subprocess.Popen(mlCmd, shell=inShell, bufsize=-1, cwd=config.get('tablesGeneratorDir')) #stdout=subprocess.STDOUT, stderr=subprocess.STDOUT)
                mlProc.wait()

        #compute scaffold-contig consistency
        #(call rb script, read it)
        consOutBuff = OutFileBuffer(Common.createTagFilePath(outputDir, inputFastaFile, 'cons'))

        cons = None
        if scaffoldsToContigsMapFile != None:
            cons = Consistency(inputFastaFile,
                               str(Common.createTagFilePath(outputDir, inputFastaFile, 'pOUT')),
                               scaffoldsToContigsMapFile,
                               databaseFile,
                               minScaffContigCount=None, minScaffBpLen=None, cladesSet=None, considerContigWithNoScaff=True)
        else:
            if os.path.isfile(os.path.normpath(str(fastaFileIds + '.out'))):
                cons = Consistency(fastaFileIds,
                                   str(fastaFileIds + '.out'),
                                   scaffoldContigMapIdsFile,
                                   databaseFile,
                                   minScaffContigCount=None, minScaffBpLen=None, cladesSet=None, considerContigWithNoScaff=True)



        if cons != None:
            consOutBuff.writeText(cons.getScaffoldsPrint() + '\n')
            consOutBuff.writeText(cons.getGroupedScaffoldsPrint() + '\n')
            consOutBuff.close()
            cons.close()


        #scCmd = str(os.path.normpath(os.path.join(configPPS.get('ppsInstallDir'), 'scripts', 'scaffold_consistency.rb'))
        #                + ' -p ' + fastaFileIds + '.out'
        #                + ' -n ' + os.path.join(workingDir, 'projectDir', 'tree.newick')
        #                + ' -f ' + fastaFileIds
        #                + ' -s ' + databaseFile
        #                + ' -c ' + scaffoldContigMapIdsFile
        #                + ' --ncbi-taxon-id '
        #                + ' > '  + Common.createTagFilePath(outputDir, inputFastaFile, 'cons')
        #                )
        #if os.name == 'posix':
        #    scProc = subprocess.Popen(scCmd, shell=True, bufsize=-1, cwd=configPPS.get('ppsInstallDir')) #stdout=subprocess.STDOUT, stderr=subprocess.STDOUT)
        #    scProc.wait()
        #    print 'PPS "scaffold_consistency" return code:', scProc.returncode
        #else:
        #    print 'Cannot run scaffold_consistency since your system is not "posix" but', str('"' + os.name + '"'), 'CMD:\n', scCmd

        #cross validation: how were training data assigned
        taxonomy = Taxonomy(databaseFile, taxonomicRanks)
        if os.path.isfile(Common.createTagFilePath(workingDir, fastaFileIds, 'out')):
            placementPPS = SSDCrossVal.ppsOut2Placements(Common.createTagFilePath(workingDir, fastaFileIds, 'out'), None)
            SSDPlacement = SSDCrossVal.ssd2Placements(os.path.join(workingDir, 'sampleSpecificDir'), None)
            cmpList = SSDCrossVal.cmpPlacements(SSDPlacement, placementPPS, taxonomy, taxonomicRanks)
            SSDCrossVal.cmp2Summary(cmpList, Common.createTagFilePath(workingDir, fastaFileIds, 'ssd_cross'))

            #store the sequences that should not be used as sample specific data for recpective clades
            #param: rankCut=0~bacteria, maxDist=10 ~ not limited, True ~ output the forbidden list
            forbiddenList = SSDCrossVal.filterCmpList(cmpList, 0, 10, taxonomy, True)
            TrainingData.updateForbiddenList(forbiddenList, Common.createTagFilePath(workingDir, fastaFileIds, 'ssd_forbidden'))


        #for i in list:
        #    print i

        #comparison to ref assignment:
        placementPPSn = SSDCrossVal.ppsOut2Placements(Common.createTagFilePath(outputDir, inputFastaFile, 'pOUT'), None)
        if os.path.isfile(referencePlacementFileOut):
            placementPPSref = SSDCrossVal.ppsOut2Placements(referencePlacementFileOut, None)
            cmpListR = SSDCrossVal.cmpPlacements(placementPPSref, placementPPSn, taxonomy, taxonomicRanks)
            SSDCrossVal.cmp2Summary(cmpListR, Common.createTagFilePath(workingDir, fastaFileIds, 'cmp_ref'))



        taxonomy.close()




def getLogFileName(logDir, description):
    now = datetime.datetime.now()
    timeStamp = str(str(now.year) + str(now.month) + str(now.day) + '_' + str(now.hour) + str(now.minute) + str(now.second))
    fileName = timeStamp + '_' + str(description) + '.txt'
    return os.path.normpath(os.path.join(logDir,fileName))



    #update Sequences according to the PhyloPythiaS assignments
    #option override MlTree map assignments by PhyloPythiaS assignments?



    #discard assignments with low consistency
    #option consistency in % e.g. 80%

    #output the right assignments by PPS and scaffold-contig consistency
    #option: regexp that defines which part of the name should be output as a sequences name

    #compute one table for each rank up to rankIdCut in csv format
    #option input file reference placements


#os.path.normpath

    #print args.config
    #print args.renamed_fasta_file




    #print parser.print_help()

if __name__ == "__main__":
  main()