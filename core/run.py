#!/usr/bin/env python

"""
    Master script of PPSplus.
"""

import sys
import os
import shutil
import argparse
import subprocess
import datetime

from com import common
from com import csv
from com.config import Config
from com.config import Config2
from core import ssd_eval
from core.training_data import PPSInput
from core.analysis_mg import MarkerGeneAnalysis
from core.cluster import MGCluster
from core.taxonomy import Taxonomy
from core.sequences import Sequences
from core.analysis16s import RRNA16S

from eval.consistency import Consistency
from misc import out_proc

from old.johdro import Johdro

# hera:
# export PATH=/net/programs/Debian-6.0.3-x86_64/python-2.7joh/bin:/net/programs/Debian-6.0.3-x86_64/jdk-1.7.0.07/bin:$PATH
# export PYTHONPATH=/net/metagenomics/projects/PPSmg/scripts/scriptsR13
#
#
def main():
    # external commands will be executed in Shell in Unix/Linux
    assert os.name == 'posix', str('The pipeline runs only on "posix" systems (i.e. Unix/Linux compatible). ' +
                                   'Your system is "' + os.name + '"')

    parser = argparse.ArgumentParser(
        description='''Pipeline that pre-processes the input for PhyloPythiaS and run it.''',
        epilog='''Read the user documentation for more details.''')

    parser.add_argument('-c', '--config', nargs=1, type=file, required=True, help='configuration file of the pipeline',
                        metavar='config.cfg', dest='config')

    parser.add_argument('-n', '--run-rrna16S', action='store_true',
                        help='run hidden markov model and classify according to the 16S, 23S, and 5S marker genes',
                        dest='n')

    parser.add_argument('-g', '--run-marker-gene-analysis', action='store_true',
                        help='run hidden markov model and classify according to the marker genes (from Amphora)',
                        dest='g')

    #parser.add_argument('-j', '--run-johdro', action='store', nargs='+',
    #                    help='run blast (b), run last (l), run LCA for blast (lcab), run LCA for last (lcal)',
    #                    dest='j')

    parser.add_argument('-o', '--process-preprocessed-output', action='store', nargs='+',
                        help='process output from the 16S rRNA analysis (s16), marker gene Amphora analysis (mg)',
                        choices=["s16", "mg"],
                        dest='o')
                        # push down predictions to more specific clades (sc)
                        # build new OTUs (otu)
                        # Taxator Blast LCA (jlcab), Taxator Last LCA (jlcal)
                        # MlTreeMap (m); Amphora (ah); exclude some data from SSD (ex) place sequences;
                        # (exSSD) exclude defined contigs from SSD; create training data'

    parser.add_argument('-t', '--pps-train', action='store_true',
                        help='run PhyloPythiaS "training phase"',
                        dest='t')

    parser.add_argument('-p', '--pps-predict', action='store', nargs='+', choices=["c","s","v"],
                        help="""run PhyloPythiaS "predict phase" (c) for contigs, (s) for scaffolds,
                        (v) compare predictions of contigs and scaffolds""",
                        dest='p')

    parser.add_argument('-r', '--read-pps-out', action='store_true',
                        help="""Reads the output placements of PhyloPythiaS (otherwise just predictions based on the
                        marker gene analysis will appear in the results)""",
                        dest='r')

    parser.add_argument('-s', '--summary', action='store_true',
                        help='Summary, output the results.',
                        dest='s')

    args = parser.parse_args()

    #read configuration
    config = Config(args.config[0], 'pPPS')
    configJoh = Config2(config, 'Johdro')
    configPPS = Config2(config, 'PPS')
    configRRNA16S = Config2(config, 'RRNA16S')
    configMG = Config2(config, 'MarkerGenes')

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

    fastaFileIds = common.createTagFilePath(workingDir , inputFastaFile, 'ids')

    if inputFastaScaffoldsFile is not None:
        fastaFileScaffoldsIds = common.createTagFilePath(workingDir , inputFastaScaffoldsFile, 'ids')
    else:
        fastaFileScaffoldsIds = None

    seqNameSeqIdFile = common.createTagFilePath(workingDir, inputFastaFile, 'cToIds')
    if inputFastaScaffoldsFile is not None:
        scaffNameScaffIdFile = common.createTagFilePath(workingDir, inputFastaScaffoldsFile , 'sToIds')
    else:
        scaffNameScaffIdFile = None
    scaffoldContigMapIdsFile = common.createTagFilePath(workingDir, inputFastaFile, 'mapSCIds')

    summaryAllFile = os.path.join(outputDir, 'summaryAll.txt')
    summaryTrainFile = os.path.join(outputDir, 'summaryTrain.txt')
    taxonomicRanks = config.get('taxonomicRanks').split(',')
    databaseFile = os.path.join(os.path.normpath(config.get('databaseFile')), 'ncbitax_sqlite.db')

    referencePlacementFilePPOut = os.path.normpath(config.get('referencePlacementFilePPOut'))
    referencePlacementFileOut = os.path.normpath(config.get('referencePlacementFileOut'))
    useScaffInSSDPred = eval(config.get('useScaffInSSDPred'))
    ppsConfigFilePath = os.path.join(workingDir, 'PPS_config_generated.txt')

    # generates working fasta files always
    sequences = Sequences(config, inputFastaFile, inputFastaScaffoldsFile, scaffoldsToContigsMapFile)
    sequences.writeSequences(fastaFileIds)
    print('Working contigs input fasta file created: %s' % fastaFileIds)
    sequences.writeSeqNameSeqId(seqNameSeqIdFile)
    print('Ids mapping for the working contigs fasta file created: %s' % seqNameSeqIdFile)
    sequences.writeScaffoldContigMap(scaffoldContigMapIdsFile)
    print('Scaffolds -> contigs map ids file created: %s' % scaffoldContigMapIdsFile)
    #assert Common.seqFileCmp(inputFastaFile, fastaFileIds), 'The fasta IDs file contains different sequences!'
    if fastaFileScaffoldsIds is not None:
        sequences.writeScaffolds(fastaFileScaffoldsIds)
        print('Working scaffolds input fasta file created: %s' % fastaFileScaffoldsIds) #(scaffold id -> sequence)
        sequences.writeScaffNameScaffId(scaffNameScaffIdFile)
        print('Ids mapping for the working scaffolds fasta file created: %s' % scaffNameScaffIdFile)

    #set up logging
    logDir = os.path.join(outputDir,'log')
    if not os.path.exists(logDir):
        os.mkdir(logDir)

    #is it specified what to do?
    if not args.n and not args.g and not args.o and not args.t and not args.p and not args.r and not args.s:
        print('Choose what do you want to do!') # args.j
        print(parser.print_help())

    # run 16S analysis
    if args.n:
        rrna = RRNA16S(config, configRRNA16S)
        print('run Hidden Markov Model for (16S, 23S, 5S)')
        rrna.runHMM(fastaFileIds, outLog=getLogFileName(logDir, 'hmm16S'))
        print('run 16S classification')
        rrna.classify16S(fastaFileIds, outLog=getLogFileName(logDir, 'mothurClassify16S'))
        print('run 23S classification')
        rrna.classify23S(fastaFileIds, outLog=getLogFileName(logDir, 'mothurClassify23S'))
        print('run 5S classification')
        rrna.classify5S(fastaFileIds, outLog=getLogFileName(logDir, 'mothurClassify5S'))

    #marker gene analysis (Amphora)
    if args.g:
        mg = MarkerGeneAnalysis(config, configMG, configRRNA16S, mgWorkingDir)
        print 'run (Amphora) marker gene analysis'
        mg.runMarkerGeneAnalysis(fastaFileIds, outLog=getLogFileName(logDir, 'amphoraMG'))

# JOHDRO    #run johdro (blast|last + predict)
#    if args.j:
#        johdro = Johdro(config, configJoh)
#        if 'b' in args.j: #run blast
#            print('run Blast for contigs')
#            johdro.runBlast(fastaFileIds)
#            if fastaFileScaffoldsIds is not None:
#                print('run Blast for scaffolds')
#                johdro.runBlast(fastaFileScaffoldsIds)
#
#        if 'l' in args.j: #run last
#            print('run Last for contigs')
#            johdro.runLast(fastaFileIds)
#            if fastaFileScaffoldsIds is not None:
#                print('run Last for scaffolds')
#                johdro.runLast(fastaFileScaffoldsIds)
#
#        if 'lcab' in args.j: #run LCA blast
#            print('run LCA Blast for contigs')
#            johdro.runBlastLCA(fastaFileIds)
#            if fastaFileScaffoldsIds is not None and useScaffInSSDPred:
#                print('run LCA Blast for scaffolds')
#                johdro.runBlastLCA(fastaFileScaffoldsIds)
#
#        if 'lcal' in args.j: #run LCA last
#            print('run LCA Last for contigs')
#            johdro.runLastLCA(fastaFileIds)
#            if fastaFileScaffoldsIds is not None and useScaffInSSDPred:
#                print('run LCA Last for scaffolds')
#                johdro.runLastLCA(fastaFileScaffoldsIds)

    # remove non-DNA characters in the sequences of contigs and scaffolds
    # this is important before the sample specific data are created and before PPS predicts the data
    sequences.setRemoveNonDna(True)
    # the working fasta files of scaffolds and contigs are regenerated, so that they
    # wouldn`t contain non-DNA characters
    sequences.writeSequences(fastaFileIds)
    print('Working contigs input fasta file without non-DNA chars created: %s' % fastaFileIds)
    if fastaFileScaffoldsIds is not None:
        sequences.writeScaffolds(fastaFileScaffoldsIds) #(scaffold id -> sequence)
        print('Working scaffolds input fasta file without non-DNA chars created: %s' % fastaFileScaffoldsIds)

    # process output of the marker gene analysis and create output for PPS
    if args.o:
        if sequences is None:
            sequences = Sequences(config, inputFastaFile, inputFastaScaffoldsFile, scaffoldsToContigsMapFile)
        taxonomy = Taxonomy(databaseFile, taxonomicRanks)

        #place sequences according to the 16S analysis
        if 's16' in args.o:
            rrna = RRNA16S(config, configRRNA16S)
            print('process the 16S, 23S, and 5S output')
            countList = rrna.setCandidatePlacementFrom16S23S5S(sequences, taxonomy, fastaFileIds)
            print(str('Placed sequences by the\n16S:' + str(countList[0]) + '\n23S:' + str(countList[1]) +
                      '\n5S:' + str(countList[2]) + '\nall 16S, 23S, and 5S:' + str(countList[3])))

        if 'mg' in args.o:
            mg = MarkerGeneAnalysis(config, configMG, configRRNA16S, mgWorkingDir)
            print('process marker gene analysis (Amphora) output')
            num = mg.setCandidatePlacement(sequences, taxonomy, fastaFileIds)
            print('Placed sequences by the marker gene analysis (Amphora): %s' % num)

# JOHDRO
#        #place sequences from johdro Blast LCA (as candidate placements)
#        if 'jlcab' in args.o:
#            johdro = Johdro(config, configJoh)
#            print('process johdro Blast LCA output')
#            count = johdro.setCandidatePlacementFromBlastLCA(sequences, taxonomy, fastaFileIds)
#            print('Placed sequences by johdro Blast LCA: %s' % count)
#            if fastaFileScaffoldsIds is not None and useScaffInSSDPred:
#                count = johdro.setScaffCandidatePlacementFromBlastLCA(sequences, taxonomy, fastaFileScaffoldsIds)
#                print('Placed sequences by johdro Blast LCA based on scaffolds: %s' % count)
        #place sequences from johdro Last LCA (as candidate placements)
#        if 'jlcal' in args.o:
#            johdro = Johdro(config, configJoh)
#            print('process johdro Last LCA output')
#            count = johdro.setCandidatePlacementFromLastLCA(sequences, taxonomy, fastaFileIds)
#            print('Placed sequences by johdro Last LCA: %s' % count)
#            if fastaFileScaffoldsIds is not None and useScaffInSSDPred:
#                count = johdro.setScaffCandidatePlacementFromLastLCA(sequences, taxonomy, fastaFileScaffoldsIds)
#                print('Placed sequences by johdro Last LCA based on scaffolds: %s' % count)

        # place sequences according to the candidate placements
        sequences.setTaxonomyPathsFromCandidatePaths(taxonomy, float(config.get('candidatePlTopPercentThreshold')))
        print('Sequences placed based on candidate placements, placed: %s' % len(sequences.placedSeqSet))

        # assign not placed contigs of one scaffold to the lowest common ancestor of assigned contigs
        if eval(config.get('placeContigsFromTheSameScaffold')):
            sequences.placeContigsFromTheSameScaffold(taxonomy, float(config.get('agThreshold')),
                                                      float(config.get('assignedPartThreshold')),
                                                      float(config.get('candidatePlTopPercentThreshold')))
            print('Not assigned contigs of one scaffold are placed to the lowest common ancestor of assigned contigs, '
                  'placed: %s' % len(sequences.placedSeqSet))

        # marker gene clustering (NOT TESTED PROPERLY)
        clust = MGCluster(config, configRRNA16S, configMG, mgWorkingDir, fastaFileIds, sequences,
                          taxonomy, os.path.basename(inputFastaFile))
        if 'sc' in args.o: #construct specific predictions
            clust.refineSpecificPred()

        if 'otu' in args.o: #build new otus
            clust.reconstructOTU()

        taxonomy.close()

        fFilePath = common.createTagFilePath(workingDir, fastaFileIds, 'ssd_forbidden')
        if ('ex' in args.o) and (os.path.isfile(fFilePath)):
            forbiddenDict = PPSInput.loadDictFromAFile(fFilePath)
            print 'read file of forbidden entries'
        else:
            forbiddenDict = None

        exSSDFilePath = os.path.join(workingDir, 'excludeContigsAsSSD.txt')
        if ('exSSD' in args.o) and ((os.path.isfile(exSSDFilePath))):
            exSSDContigNameList = csv.getColumnAsList(exSSDFilePath, entryModifyFunction=None, colNum=0,
                                                                      sep='\n', comment='#')
            print 'read file of contigs to be excluded from SSD'
        else:
            exSSDContigNameList = None

        # create input for PPS
        pps = PPSInput(sequences, taxonomicRanks, forbiddenDict, exSSDContigNameList, summaryAllFile)

        pps.createPPSInputFiles(
            os.path.join(workingDir, 'ncbids.txt'),
            os.path.join(workingDir, 'sampleSpecificDir'),
            int(config.get('rankIdAll')),
            int(config.get('rankIdCut')),
            int(config.get('rankIdCutMinBp')),
            float(config.get('minPercentInLeaf')),
            int(config.get('maxLeafClades')),
            int(config.get('minBpToModel')),
            int(config.get('minGenomesWgs')),
            int(config.get('minBpPerSpecies')),
            os.path.normpath(config.get('genomesWGS')),
            None, #forbiddenDict
            databaseFile,
            taxonomicRanks,
            int(config.get('fastaLineMaxChar')),
            int(config.get('minSSDfileSize')),
            int(config.get('maxSSDfileSize')),
            float(config.get('weightStayAll')),
            summaryTrainFile)

    # generate PPS configuration file (before train or predict)
    if args.t or args.p:
        keyDict = {'NCBI_PROCESSED_DIR': str(os.path.normpath(config.get('genomesWGS'))),
                   'NCBI_TAX_DIR': os.path.dirname(databaseFile),
                   'PROJECT_DIR': os.path.join(workingDir, 'projectDir'),
                   'TREE_FILE': os.path.join(workingDir, 'ncbids.txt'),
                   'SAMPLE_SPECIFIC_DIR': os.path.join(workingDir, 'sampleSpecificDir'),
                   'PARALLEL_MODELS': 'FALSE', # optional
                   'GENOMES_EXCLUDE': ''} # optional
        createPPSConfig(ppsConfigFilePath, keyDict)

    # run PhyloPythiaS train
    if args.t:
        trainCmd = str(os.path.normpath(os.path.join(configPPS.get('ppsInstallDir'), 'scripts', 'train.rb')) +
                       ' ' + ppsConfigFilePath)
        logOut = open(getLogFileName(logDir, 'PPS_train'), 'w')
        trainProc = subprocess.Popen(trainCmd, shell=True, bufsize=-1, cwd=configPPS.get('ppsInstallDir'),
                                     stdout=logOut) #subprocess.STDOUT, stderr=subprocess.STDOUT)
        trainProc.wait()
        logOut.close()
        print('PPS "train" return code: %s' % trainProc.returncode)

    # run PhyloPythiaS predict
    if args.p:
        if 'c' in args.p: # predict contigs
            predictCmd = str(os.path.normpath(os.path.join(configPPS.get('ppsInstallDir'), 'scripts', 'predict.rb')) +
                             ' ' + fastaFileIds + ' ' + ppsConfigFilePath)

            logOut = open(getLogFileName(logDir, 'PPS_predict_c'),'w')
            predictProc = subprocess.Popen(predictCmd, shell=True, bufsize=-1, cwd=configPPS.get('ppsInstallDir'),
                                           stdout=logOut) #stdout=subprocess.STDOUT, stderr=subprocess.STDOUT)
            predictProc.wait()
            logOut.close()
            shutil.copy2(str(fastaFileIds + '.nox.fna.out'), str(fastaFileIds + '.out'))
            shutil.copy2(str(fastaFileIds + '.nox.fna.PP.out'), str(fastaFileIds + '.PP.out'))
            print('PPS "predict" contigs return code: %s' % predictProc.returncode)

        if inputFastaScaffoldsFile is not None:
            scaffoldsIdsNewPath = os.path.normpath(
                os.path.join(workingDir, 'crossVal', os.path.basename(fastaFileScaffoldsIds)))

        # predict scaffolds
        if 's' in args.p:
            if inputFastaScaffoldsFile is None:
                print('Cannot predict a scaffold file, since it is not specified!')
            else:
                crossValDir = os.path.normpath(os.path.join(workingDir, 'crossVal'))
                if not os.path.isdir(crossValDir):
                    os.mkdir(crossValDir)
                src = fastaFileScaffoldsIds
                dst = os.path.join(crossValDir, os.path.basename(fastaFileScaffoldsIds))
                try:
                    os.symlink(src, dst)
                except:
                    sys.stderr.write(str('Can`t create symbolic link: ' + dst + '\n'))
                    shutil.copyfile(src, dst)
                    sys.stderr.write(str('Copy of the file created: ' + dst + '\n'))
                #str('ln -s ' + fastaFileScaffoldsIds + ' ' + crossValDir + ';')
                predictCmd = str(os.path.normpath(
                    os.path.join(configPPS.get('ppsInstallDir'), 'scripts', 'predict.rb')) + ' ' +
                                 scaffoldsIdsNewPath + ' ' + os.path.normpath(configPPS.get('configPPS')))

                logOut = open(getLogFileName(logDir, 'PPS_predict_s'),'w')
                predictProc = subprocess.Popen(predictCmd, shell=True, bufsize=-1,
                                               cwd=configPPS.get('ppsInstallDir'),
                                               stdout=logOut) #stdout=subprocess.STDOUT, stderr=subprocess.STDOUT)
                predictProc.wait()
                logOut.close()
                shutil.copy2(str(scaffoldsIdsNewPath + '.nox.fna.out'), str(scaffoldsIdsNewPath + '.out'))
                shutil.copy2(str(scaffoldsIdsNewPath + '.nox.fna.PP.out'), str(scaffoldsIdsNewPath + '.PP.out'))
                print('PPS "predict" scaffolds return code: %s' % predictProc.returncode)

        # compare prediction of scaffolds and contigs
        if 'v' in args.p:
            contigPred = os.path.normpath(str(fastaFileIds + '.PP.out'))
            scaffPred =  os.path.normpath(str(scaffoldsIdsNewPath + '.PP.out'))
            if not os.path.isfile(scaffPred):
                print('Comparison of contig/scaffold predictions: scaffolds were not predicted!')
            else:
                if not os.path.isfile(contigPred):
                    print('Comparison of contig/scaffold predictions: contigs were not predicted!')
                else:
                    #create prediction file for contigs from a prediction file for scaffolds
                    scaffPredAsContigs = os.path.normpath(
                        os.path.join(workingDir, 'crossVal', 'scaffPredAsContigs.PP.out'))
                    print('scaff to contigs:', scaffoldsToContigsMapFile, scaffPred, scaffPredAsContigs)
                    out_proc.scafToContigOutput(scaffoldContigMapIdsFile, scaffPred, scaffPredAsContigs)
                    i = 0
                    for rank in taxonomicRanks:
                        if i > int(config.get('tablesGeneratorMaxRank')):
                            break
                        i += 1
                        # do comparison: contigs` prediction vs scaffolds` prediction
                        mlCmd = str('java -jar tables.jar ' + fastaFileIds + ' ' + contigPred + ' ' +
                                    scaffPredAsContigs + ' ' +
                                    os.path.normpath(os.path.join(
                                        os.path.dirname(scaffPredAsContigs), str('crossVal_' + str(rank) + '.csv'))) +
                                    ' ' + str(rank) + ' ' + 'unassigned' + ' ' + str(config.get('taxonomicRanks'))
                                    )
                        mlProc = subprocess.Popen(mlCmd, shell=True, bufsize=-1,
                                                  cwd=config.get('tablesGeneratorDir'))
                                                  #stdout=subprocess.STDOUT, stderr=subprocess.STDOUT)
                        mlProc.wait()
                        print('Comparison of contig/scaffold predictions: %s return code: %s' %
                              (mlCmd, mlProc.returncode))

    # Read output from PPS and place sequences
    if args.r:
        if sequences is None:
            sequences = Sequences(config, inputFastaFile, inputFastaScaffoldsFile, scaffoldsToContigsMapFile)
        taxonomy = Taxonomy(databaseFile, taxonomicRanks)

        pps.readPPSOutput(sequences, taxonomy, fastaFileIds, eval(config.get('overrideWithPPSPlacements')))

        #discard inconsistent contigs or scaffolds
        #discardedContigs = sequences.discardInconsistentPlacements(float(config.get('discardContigThreshold')),
        #                                                           float(config.get('discardScaffThreshold')),
        #                                                           taxonomicRanks)
        #print "The number of discarded contigs while reading PPS output:", discardedContigs

        taxonomy.close()
        print('Placed sequences after PPS predict: %s' % len(sequences.placedSeqSet))


    if args.s:
        #create PP.pOUT file
        sequences.writePlacementsPPOut(str(common.createTagFilePath(outputDir, inputFastaFile, 'PP.pOUT')),
                                       taxonomicRanks, config.get('outputFileContigSubPattern'))

        #create .pOUT file (name tab ncbid)
        sequences.writePlacementsOut(str(common.createTagFilePath(outputDir, inputFastaFile, 'pOUT')),
                                     taxonomicRanks, config.get('outputFileContigSubPattern'))

        sequences.writeSequences(str(common.createTagFilePath(workingDir, inputFastaFile, 'pOUT.fas')),
                                 writeIds=False, outputFileContigSubPattern=config.get('outputFileContigSubPattern'))

        # generate comparison tables
        i = 0
        if os.path.isfile(referencePlacementFilePPOut):
            for rank in taxonomicRanks:
                if i > int(config.get('tablesGeneratorMaxRank')):
                    break
                i += 1
                mlCmd = str('java -jar tables.jar ' +
                            str(common.createTagFilePath(workingDir, inputFastaFile, 'pOUT.fas')) + ' ' +
                            str(referencePlacementFilePPOut) + ' ' +
                            str(common.createTagFilePath(outputDir, inputFastaFile, 'PP.pOUT')) + ' ' +
                            str(common.createTagFilePath(outputDir, inputFastaFile,
                                                         str('table.' + str(rank) + '.csv')))  + ' ' +
                            str(rank) + ' ' + 'unassigned' + ' ' +
                            str(config.get('taxonomicRanks'))
                )
                mlProc = subprocess.Popen(mlCmd, shell=True, bufsize=-1, cwd=config.get('tablesGeneratorDir'))
                # stdout=subprocess.STDOUT, stderr=subprocess.STDOUT)
                mlProc.wait()

        #compute scaffold-contig consistency
        #(call rb script, read it ???)
        consOutBuff = csv.OutFileBuffer(common.createTagFilePath(outputDir, inputFastaFile, 'cons'))

        cons = None
        if scaffoldsToContigsMapFile is not None:
            cons = Consistency(inputFastaFile,
                               str(common.createTagFilePath(outputDir, inputFastaFile, 'pOUT')),
                               scaffoldsToContigsMapFile, databaseFile, minScaffContigCount=None,
                               minScaffBpLen=None, cladesSet=None, considerContigWithNoScaff=True)
        else:
            if os.path.isfile(os.path.normpath(str(fastaFileIds + '.out'))):
                cons = Consistency(fastaFileIds,
                                   str(fastaFileIds + '.out'),
                                   scaffoldContigMapIdsFile,
                                   databaseFile, minScaffContigCount=None,
                                   minScaffBpLen=None, cladesSet=None, considerContigWithNoScaff=True)

        if cons is not None:
            consOutBuff.writeText(cons.getScaffoldsPrint() + '\n')
            consOutBuff.writeText(cons.getGroupedScaffoldsPrint() + '\n')
            consOutBuff.close()
            cons.close()

#        # ODL VERSION (don`t use!) of scaffold contig consistency
#        scCmd = str(os.path.normpath(
#            os.path.join(configPPS.get('ppsInstallDir'), 'scripts', 'scaffold_consistency.rb')) +
#                    ' -p ' + fastaFileIds + '.out' +
#                    ' -n ' + os.path.join(workingDir, 'projectDir', 'tree.newick') +
#                    ' -f ' + fastaFileIds +
#                    ' -s ' + databaseFile +
#                    ' -c ' + scaffoldContigMapIdsFile +
#                    ' --ncbi-taxon-id ' + ' > ' + common.createTagFilePath(outputDir, inputFastaFile, 'cons'))
#        scProc = subprocess.Popen(scCmd, shell=True, bufsize=-1, cwd=configPPS.get('ppsInstallDir'))
#        #stdout=subprocess.STDOUT, stderr=subprocess.STDOUT)
#        scProc.wait()
#        print('PPS "scaffold_consistency" return code: %s' % scProc.returncode)

        # validation: how were training data assigned
        taxonomy = Taxonomy(databaseFile, taxonomicRanks)
        if os.path.isfile(common.createTagFilePath(workingDir, fastaFileIds, 'out')):
            placementPPS = ssd_eval.ppsOut2Placements(common.createTagFilePath(workingDir, fastaFileIds, 'out'), None)
            SSDPlacement = ssd_eval.ssd2Placements(os.path.join(workingDir, 'sampleSpecificDir'), None)
            cmpList = ssd_eval.cmpPlacements(SSDPlacement, placementPPS, taxonomy, taxonomicRanks)
            ssd_eval.cmp2Summary(cmpList, common.createTagFilePath(workingDir, fastaFileIds, 'ssd_cross'))

            #store the sequences that should not be used as sample specific data for respective clades
            #param: rankCut=0~bacteria, maxDist=10 ~ not limited, True ~ output the forbidden list
            forbiddenList = ssd_eval.filterCmpList(cmpList, 0, 10, taxonomy, True)
            PPSInput.updateForbiddenList(forbiddenList,
                                         common.createTagFilePath(workingDir, fastaFileIds, 'ssd_forbidden'))

        # comparison to ref assignment:
        placementPPSn = ssd_eval.ppsOut2Placements(common.createTagFilePath(outputDir, inputFastaFile, 'pOUT'), None)
        if os.path.isfile(referencePlacementFileOut):
            placementPPSref = ssd_eval.ppsOut2Placements(referencePlacementFileOut, None)
            cmpListR = ssd_eval.cmpPlacements(placementPPSref, placementPPSn, taxonomy, taxonomicRanks)
            ssd_eval.cmp2Summary(cmpListR, common.createTagFilePath(workingDir, fastaFileIds, 'cmp_ref'))

        taxonomy.close()


def getLogFileName(logDir, description):
    now = datetime.datetime.now()
    timeStamp = str(str(now.year) + str(now.month) + str(now.day) + '_' + str(now.hour) + str(now.minute) + str(now.second))
    fileName = timeStamp + '_' + str(description) + '.txt'
    return os.path.normpath(os.path.join(logDir,fileName))


def createPPSConfig(ppsConfigFilePath, keyDict):
    """
        Generates the up to date PPS configuration file

        @param keyDict: must specify: NCBI_PROCESSED_DIR, NCBI_TAX_DIR, PROJECT_DIR, TREE_FILE, SAMPLE_SPECIFIC_DIR,
            can be specified: PARALLEL_MODELS, GENOMES_EXCLUDE
        @type keyDict: dict
    """
    out = csv.OutFileBuffer(ppsConfigFilePath)
    out.writeText(
        """#######################################################################
        #configuration file (PPS+ GENERATED !!!)
        #please make sure that there is no space before or after ":"
        #lines starting with character "#" are treated as comments
        #please provide complete paths instead of only file or directory names
        #######################################################################
        #directory where processed NCBI data is stored, provide empty directory to create new
        #REUSABLE\n""")
    out.writeText('NCBI_PROCESSED_DIR:%s\n' % keyDict.get('NCBI_PROCESSED_DIR',''))
    out.writeText(
        """#Directory containing NCBI taxonomy in SQlite3 format with file name "ncbitax_sqlite.db"
        #provide empty directory to create new database
        #REUSABLE\n""")
    out.writeText('NCBI_TAX_DIR:%s\n' % keyDict.get('NCBI_TAX_DIR', ''))
    out.writeText('#project directory, the directory must be empty\n')
    out.writeText('PROJECT_DIR:%s\n' % keyDict.get('PROJECT_DIR',''))
    out.writeText(
        """#############################
        #!!!FOLLOWING ARE OPTIONAL!!!
        #############################
        ###### Output space options #####
        #a file containing a tree in newick format (see restrictions in INSTALL.txt)
        #OR a file with ncbi taxon ids (one id per line) to create a tree from\n""")
    out.writeText('TREE_FILE:%s\n' % keyDict.get('TREE_FILE',''))
    out.writeText(
        """#Taxonomic ranks (comma separated, no space) starting at the lowest rank. \
        Please make sure that "root" is there at the end.
        TAXONOMY_RANKS:species,genus,family,order,class,phylum,superkingdom,root
        #number of minimum genomes a clade must contain to be included in generic model
        #effective only if tree file is not provided
        N_MIN_GENOMES_GENERIC:3
        #action on loss 0:disabled, 1:invert
        LOSS_ACTION:0
        ###### Input space options #####
        #a directory with sample specific fasta files (file names must start with appropriate organism/species \
        ncbi taxonomic id)
        #leave empty if you don't have any\n""")
    out.writeText('SAMPLE_SPECIFIC_DIR:%s\n' % keyDict.get('SAMPLE_SPECIFIC_DIR',''))
    out.writeText(
        """#kmer feature space for multiple kmers use kmer_min-kmer_max
        KMER:4-6
        #Fragment lengths for different models (comma separated, no space)
        FRAGMENT_LEN:1000,3000,5000,10000,15000,50000
        #kmer feature
        #use reverse complement for computing kmer features?
        REV_COMPLEMENT:1
        #remove reverse complement features?
        RM_REV_COMPLEMENT:1
        #0:disabled, 1:sequence length, 2:sequence_length-k+1, 3:embedded monomer frequency
        KMER_NORMALIZATION:1
        #Number of examples per training file
        NUMBER_EXAMPLES:10000
        #step size for sample specific data; either a single number (for all fragment lengths) or an array separated with ","
        SAMPLE_SPECIFIC_STEP:1000,300,500,1000,1500,5000
        ###### Training options #####
        #C values for SVM, if single value is given then models will be build with that value.
        #If comma separated (no space) values are given then cross-validation will be performed.
        #If a single value is provided, all models will be built with it. Our experience shows that in general
        #values less than 1 (e.g. 0.01 and 0.1) do not provide good models.
        C_GRID:1000
        #clean-up the data (sampled_fasta and train_data directories) created after training? TRUE/FALSE
        CLEAN_UP_TRAIN:TRUE
        #kernel type 0:linear, 1:polynomial, 2:rbf (on-linear kernels are computationally expensive)
        KERNEL:0
        ##polynomial kernel degree
        KERNEL_POLYNOMIAL_DEGREE:2
        ##rbf kernel gamma
        KERNEL_RBF_GAMMA:1
        ##polynomial kernel s
        KERNEL_POLYNOMIAL_S:1
        ######  Predictions options #####
        #number of classifiers to use, keep this odd to avoid ties
        N_CLASSIFIERS:3
        #Create Pie charts for every taxonomic rank TRUE/FALSE (in prediction)
        #slice colors are determined automatically so no color consistency is guaranteed
        PIE_CHARTS:TRUE
        ###### Misc options #####
        #should the models be built in parallel (please make sure that you have enough number of \
        processors and main memory)\n""")
    out.writeText('PARALLEL_MODELS:%s\n' % keyDict.get('PARALLEL_MODELS','FALSE'))
    out.writeText(
        """#allowed file extensions
        EXTENSIONS:
        #genomes to exclude: file containing one ncbi tax_id per line\n""")
    out.writeText('GENOMES_EXCLUDE:%s\n' % keyDict.get('GENOMES_EXCLUDE',''))
    out.writeText(
        """#if the training data is already there then just build models (TRUE/FALSE)
        ONLY_MODELS:FALSE\n""")
    out.close()


if __name__ == "__main__":
    main()