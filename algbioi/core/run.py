#!/usr/bin/env python

"""
    Master script of the PPSplus.
"""

import os
import sys
import shutil
import subprocess
import datetime
import traceback
import argparse

from algbioi.com import csv
from algbioi.com import fasta as fas
from algbioi.com import common
from algbioi.com import taxonomy_ncbi
from algbioi.com.config import Config
from algbioi.core import ssd_eval
from algbioi.core.pps import ppsOut2ppOut
from algbioi.core.pps import readPPSOutput
from algbioi.core.pps import computeTrainingAccuracy
from algbioi.core.pps import generateCladesForGeneralModel
from algbioi.core.training_data import PPSInput
from algbioi.core.analysis_mg import MarkerGeneAnalysis
from algbioi.core.cluster import MGCluster
from algbioi.core.taxonomy import Taxonomy
from algbioi.core.sequences import Sequences
from algbioi.core.analysis16s import RRNA16S
from algbioi.eval import consistency
from algbioi.eval import accuracy
from algbioi.eval import confusion_matrix
from algbioi.misc import out_proc
from algbioi.ref import mask_db

# paths on hera/gaia:
# export PATH=/net/programs/Debian-6.0.3-x86_64/python-2.7joh/bin:$PATH
# export PYTHONPATH=/net/metagenomics/projects/PPSmg/scripts/scriptsR31
#

def main():
    # external commands will be executed in Shell in Unix/Linux
    assert os.name == 'posix', str('The pipeline runs only on "posix" systems (i.e. Unix/Linux compatible). ' +
                                   'Your system is "' + os.name + '"')

    parser = argparse.ArgumentParser(
        description='''PhyloPythiaS Plus is an extension of PhyloPythiaS.''',
        epilog='''Read the user documentation for more details.''')

    parser.add_argument('-c', '--config', nargs=1, type=file, required=True, help='configuration file of the pipeline',
                        metavar='config.cfg', dest='config')

    parser.add_argument('-n', '--run-rrna16S', action='store_true',
                        help='run hidden markov model and classify according to the 16S, 23S, and 5S marker genes',
                        dest='n')

    parser.add_argument('-g', '--run-marker-gene-analysis', action='store_true',
                        help='run hidden markov model and classify according to the "31" marker genes',
                        dest='g')

    parser.add_argument('-o', '--process-preprocessed-output', action='store', nargs='+',
                        help='process output from the 16S rRNA analysis (s16), marker gene "31" analysis (mg); '
                             'to build the general model (general)' ,
                        choices=["s16", "mg", "general"],
                        dest='o')
    # push down predictions to more specific clades (sc)
    # build new OTUs (otu)
    # Taxator Blast LCA (jlcab), Taxator Last LCA (jlcal)
    # MlTreeMap (m); Amphora (ah); exclude some data from SSD (ex) place sequences;
    # (exSSD) exclude defined contigs from SSD; create training data'

    parser.add_argument('-t', '--pps-train', action='store_true', help='run the PhyloPythiaS "training phase"',
                        dest='t')

    parser.add_argument('-a', '--train-accuracy', action='store_true', help='Compute the training accuracy.',
                        dest='a')

    parser.add_argument('-p', '--pps-predict', action='store', nargs='+', choices=["c", "s", "v"],
                        help="""run the PhyloPythiaS "predict phase" (c) for contigs, (s) for scaffolds,
                        (v) compare predictions of contigs and scaffolds""",
                        dest='p')

    parser.add_argument('-r', '--read-pps-out', action='store_true',
                        help="""Reads the output placements of PhyloPythiaS (otherwise just predictions based on the
                        marker gene analysis will appear in the results if run with option '-o' as well.)""",
                        dest='r')

    parser.add_argument('-s', '--summary', action='store_true',
                        help='Summary, outputs the results, compute precision and recall, compute comparison tables '
                             '(if reference placement is available), compute scaffold-contig consistency '
                             '(if the mapping between scaffold and contigs is available)',
                        dest='s')

    args = parser.parse_args()

    # read configuration
    config = Config(args.config[0], 'PhyloPythiaS_Plus')

    # pipeline directory
    pipelineDir = os.path.normpath(config.get('pipelineDir'))
    if not os.path.isdir(pipelineDir):
        print("Pipeline directory doesn't exist: ", pipelineDir)
        return

    # create the following directories in the working/output directory if they don't exist
    workingDir = os.path.join(pipelineDir, 'working')
    mgWorkingDir = os.path.join(workingDir, 'mgWorking')
    outputDir = os.path.join(pipelineDir, 'output')
    profilesDir = os.path.join(outputDir, 'profiles')
    dirArray = [workingDir, outputDir, os.path.join(workingDir, 'projectDir'),
                os.path.join(workingDir, 'sampleSpecificDir'), mgWorkingDir, profilesDir]
    for dirPath in dirArray:
        if not os.path.isdir(dirPath):
            try:
                os.mkdir(dirPath)
            except OSError:
                print("Can't create directory", dirPath)
                return

    # read input fasta: contigs, scaffolds, mapping
    inputFastaFile = os.path.normpath(config.get('inputFastaFile'))
    if (inputFastaFile is None) or ((inputFastaFile is not None) and (not os.path.isfile(inputFastaFile))):
        print("The input fasta file %s doesn't exist" % inputFastaFile)
        return
    if ('-' in inputFastaFile) or ('+' in inputFastaFile):
        print('The input fasta file path is not allowed to contain "+" or "-" characters ' +
              'due to the "Mothur" software, given path:\n' + inputFastaFile)
        return

    if ('-' in pipelineDir) or ('+' in pipelineDir):
        print('The project directory is not allowed to contain "+" or "-" characters ' +
              'due to the "Mothur" software, directory given:\n' + pipelineDir)
        return

    inputFastaScaffoldsFile = config.get('inputFastaScaffoldsFile')
    if (inputFastaScaffoldsFile is not None) and (not os.path.isfile(inputFastaScaffoldsFile)):
        print("The given scaffolds fasta file doesn't exist: ", inputFastaScaffoldsFile)
        return

    scaffoldsToContigsMapFile = config.get('scaffoldsToContigsMapFile')
    if (scaffoldsToContigsMapFile is not None) and (not os.path.isfile(scaffoldsToContigsMapFile)):
        print("The given scaffolds to contigs map file doesn't exist: ", scaffoldsToContigsMapFile)
        return

    # create input id files (contigs, scaffolds, mappings)
    fastaFileIds = common.createTagFilePath(workingDir, inputFastaFile, 'ids')
    seqNameSeqIdFile = common.createTagFilePath(workingDir, inputFastaFile, 'cToIds')

    fastaFileScaffoldsIds = None
    scaffNameScaffIdFile = None
    if inputFastaScaffoldsFile is not None:
        fastaFileScaffoldsIds = common.createTagFilePath(workingDir, inputFastaScaffoldsFile, 'ids')
        scaffNameScaffIdFile = common.createTagFilePath(workingDir, inputFastaScaffoldsFile, 'sToIds')

    scaffoldContigMapIdsFile = common.createTagFilePath(workingDir, inputFastaFile, 'mapSCIds')

    # output files
    summaryAllFile = os.path.join(outputDir, 'summary_all.txt')
    summaryTrainFile = os.path.join(outputDir, 'summary_train.txt')
    if config.get('configPPS') is None:
        ppsConfigFilePath = os.path.join(workingDir, 'PPS_config_generated.txt')  # config will be generated
    else:
        if not os.path.isfile(config.get('configPPS')):
            print("The PPS configuration file (configPPS) '%s' doesn't exist." % config.get('configPPS'))
            return
        ppsConfigFilePath = os.path.normpath(config.get('configPPS'))  # use custom config file

    taxonomicRanks = taxonomy_ncbi.TAXONOMIC_RANKS[1:]  # without root
    try:
        minSeqLen = int(config.get('minSeqLen'))
    except Exception as e:
        print("Can't parse configuration entry (minSeqLen), make sure that it's an integer number")
        raise e

    if config.get('databaseFile') is None:
        print("The taxonomy (databaseFile) is not specified.")
        return
    if os.path.isdir(config.get('databaseFile')):
        databaseFile = os.path.join(os.path.normpath(config.get('databaseFile')), 'ncbitax_sqlite.db')
        if not os.path.isfile(databaseFile):
            print("The directory '%s' doesn't contain taxonomy file 'ncbitax_sqlite.db'")
            return
    if not os.path.isfile(databaseFile):
        print("The database file doesn't exist:", databaseFile)
        return

    # reference predictions
    referencePlacementFileOut = config.get('referencePlacementFileOut')
    if referencePlacementFileOut is not None and not os.path.isfile(referencePlacementFileOut):
        print("The reference placement file doesn't exist: %s" % referencePlacementFileOut)
        return

    # leaf clades contained in the reference
    if referencePlacementFileOut is not None:
        try:
            refLeafCladesSet = set(map(int, csv.predToDict(referencePlacementFileOut).values()))
        except Exception as e:
            print("Can't get clade taxon ids from the reference placement file: %s" % referencePlacementFileOut)
            raise e
    else:
        refLeafCladesSet = None
    referencePlacementFilePPOut = None  # config.get('referencePlacementFilePPOut')

    #if (referencePlacementFilePPOut is not None) and (not os.path.isfile(referencePlacementFilePPOut)):
    #    print("The reference file doesn't exist: ", referencePlacementFilePPOut)
    #    return
    if (referencePlacementFilePPOut is None) and (referencePlacementFileOut is not None):
        # generate PP placement file
        try:
            referencePlacementFilePPOut = os.path.join(workingDir, os.path.basename(referencePlacementFileOut) + '_.PP.out')
            ppsOut2ppOut(referencePlacementFileOut, referencePlacementFilePPOut, taxonomicRanks, databaseFile)
        except Exception as e:
            print("An error '%s' occurred while file '%s' has been generated."
                  % (e.message, referencePlacementFilePPOut))

    # generates working fasta files always when the configuration file is newer than particular files to be generated
    sequences = Sequences(inputFastaFile, inputFastaScaffoldsFile, scaffoldsToContigsMapFile, taxonomicRanks, minSeqLen)
    mtimeConfig = os.path.getmtime(config.getConfigFile())
    if ((not os.path.isfile(fastaFileIds) or not os.path.isfile(seqNameSeqIdFile) or
             not os.path.isfile(scaffoldContigMapIdsFile)) or
            (mtimeConfig > min(os.path.getmtime(fastaFileIds), os.path.getmtime(seqNameSeqIdFile),
                               os.path.getmtime(scaffoldContigMapIdsFile)))):
        sequences.writeSequences(fastaFileIds)
        print('Working contigs input fasta file created: %s' % fastaFileIds)
        sequences.writeSeqNameSeqId(seqNameSeqIdFile)
        print('Ids mapping for the working contigs fasta file created: %s' % seqNameSeqIdFile)
        sequences.writeScaffoldContigMap(scaffoldContigMapIdsFile)
        print('Scaffolds -> contigs map ids file created: %s' % scaffoldContigMapIdsFile)
        # assert Common.seqFileCmp(inputFastaFile, fastaFileIds), 'The fasta IDs file contains different sequences!'
    if ((fastaFileScaffoldsIds is not None) and
            ((not os.path.isfile(fastaFileScaffoldsIds) or not os.path.isfile(scaffNameScaffIdFile)) or
                 (mtimeConfig > min(os.path.getmtime(fastaFileScaffoldsIds), os.path.getmtime(scaffNameScaffIdFile))))):
        sequences.writeScaffolds(fastaFileScaffoldsIds)
        print('Working scaffolds input fasta file created: %s' % fastaFileScaffoldsIds)  # (scaffold id -> sequence)
        sequences.writeScaffNameScaffId(scaffNameScaffIdFile)
        print('Ids mapping for the working scaffolds fasta file created: %s' % scaffNameScaffIdFile)

    # set up logging
    logDir = os.path.join(outputDir, 'log')
    if not os.path.exists(logDir):
        try:
            os.mkdir(logDir)
        except OSError:
            print("Can't create logging directory:", logDir)
            return

    try:
        parallelPPSmodels = eval(config.get('parallelPPSmodels'))
    except:
        parallelPPSmodels = False

    # is it specified what to do?
    if not (args.n or args.g or args.o or args.t or args.p or args.r or args.s or args.a):
        print('Choose what do you want to do!')
        print(parser.print_help())
        return

    # exclude mg from the reference ?
    excludeRefMgRank = config.get('excludeRefMgRank')
    s16Database = os.path.normpath(config.get('s16Database'))
    mgDatabase = os.path.normpath(config.get('mgDatabase'))
    if (excludeRefMgRank is not None) and (args.n or args.g):
        if not mask_db.isRankAllowed(excludeRefMgRank):
            print("It's not allowed to exclude rank '%s' from the reference sequence database" % excludeRefMgRank)
        else:
            maskRefMgS16Dir = os.path.join(workingDir, str('ref_mg_16S_mask_' + excludeRefMgRank))  # mask dir
            maskRefMgDir = os.path.join(workingDir, str('ref_mg_mask_' + excludeRefMgRank))  # mask dir
            maskRefOk = True
            if ((not os.path.isdir(maskRefMgS16Dir)) or (not os.path.isdir(maskRefMgDir)) or
                    (min(os.path.getmtime(maskRefMgS16Dir), os.path.getmtime(maskRefMgDir)) < mtimeConfig)):
                if refLeafCladesSet is None:
                    print("Can't exclude reference sequences since the reference file doesn't exist.")
                    maskRefOk = False
                else:
                    # remove the mask directories if they already exist
                    if os.path.isdir(maskRefMgS16Dir):
                        shutil.rmtree(maskRefMgS16Dir)
                    if os.path.isdir(maskRefMgDir):
                        shutil.rmtree(maskRefMgDir)
                    # create the directories
                    try:
                        os.mkdir(maskRefMgS16Dir)
                        os.mkdir(os.path.join(maskRefMgS16Dir, 'db'))
                        os.mkdir(maskRefMgDir)
                        os.mkdir(os.path.join(maskRefMgDir, 'db'))
                    except OSError:
                        print("Can't create one of the directories (%s, %s, 'db') for the masked reference marker gene "
                              "sequences." % (maskRefMgS16Dir, maskRefMgDir))
                        maskRefOk = False
                    else:
                        mask_db.maskDb('mg', os.path.join(s16Database, 'db'), os.path.join(maskRefMgS16Dir, 'db'),
                                       excludeRefMgRank, refLeafCladesSet, databaseFile)
                        mask_db.maskDb('mg', os.path.join(mgDatabase, 'db'), os.path.join(maskRefMgDir, 'db'),
                                       excludeRefMgRank, refLeafCladesSet, databaseFile)
                        shutil.copy2(os.path.join(s16Database, 'content.csv'), maskRefMgS16Dir)
                        shutil.copy2(os.path.join(mgDatabase, 'content.csv'), maskRefMgDir)
                        shutil.copytree(os.path.join(mgDatabase, 'hmmAmphora'),
                                        os.path.join(maskRefMgDir, 'hmmAmphora'))
            if maskRefOk:
                s16Database = maskRefMgS16Dir
                mgDatabase = maskRefMgDir
            else:
                print("Couldn't mask reference at rank '%s'." % excludeRefMgRank)

    # run 16S analysis
    if args.n:
        rrna = RRNA16S(config, s16Database, workingDir)
        print('run Hidden Markov Model for (16S, 23S, 5S)')
        rrna.runHMM(fastaFileIds, outLog=getLogFileName(logDir, 'hmm16S'))
        print('run (16S, 23S, 5S) classification')
        rrna.classify16S(fastaFileIds, outLog=getLogFileName(logDir, 'mothurClassify16S'))
        rrna.classify23S(fastaFileIds, outLog=getLogFileName(logDir, 'mothurClassify23S'))
        rrna.classify5S(fastaFileIds, outLog=getLogFileName(logDir, 'mothurClassify5S'))

    # marker gene analysis (Amphora)
    if args.g:
        mg = MarkerGeneAnalysis(config, mgDatabase, workingDir, mgWorkingDir)
        print 'run 31 marker gene analysis'
        mg.runMarkerGeneAnalysis(fastaFileIds, outLog=getLogFileName(logDir, 'amphoraMG'))

    # remove non-DNA characters in the sequences of contigs and scaffolds
    # this is important before the sample specific data are created and before PPS predicts the data
    sequences.setRemoveNonDna(True)
    # the working fasta files of scaffolds and contigs are regenerated, so that they
    # wouldn't contain non-DNA characters
    sequences.writeSequences(fastaFileIds)
    print('Working contigs input fasta file without non-DNA chars created: %s' % fastaFileIds)
    if fastaFileScaffoldsIds is not None:
        sequences.writeScaffolds(fastaFileScaffoldsIds)
        print('Working scaffolds input fasta file without non-DNA chars created: %s' % fastaFileScaffoldsIds)

    # exclude ref. sequences from the reference ?
    refSeq = os.path.normpath(config.get('refSeq'))  # reference sequences
    excludeRefSeqRank = config.get('excludeRefSeqRank')  # rank
    if (excludeRefSeqRank is not None) and (args.o or args.p or args.t):  # exclude rank ?
        if not mask_db.isRankAllowed(excludeRefSeqRank):  # check rank
            print("It's not allowed to exclude rank '%s' from the reference sequence database!" % excludeRefSeqRank)
        else:
            maskRefSeqDir = os.path.join(workingDir, str('ref_seq_mask_' + excludeRefSeqRank))  # mask dir
            if (not os.path.isdir(maskRefSeqDir)) or (os.path.getmtime(maskRefSeqDir) < mtimeConfig):  # create mask dir
                if refLeafCladesSet is None:
                    print("Can't exclude reference sequences since the reference file doesn't exist.")
                else:
                    if os.path.isdir(maskRefSeqDir):  # remove the mask dir if it already exists
                        shutil.rmtree(maskRefSeqDir)
                    try:
                        os.mkdir(maskRefSeqDir)
                    except OSError:
                        print("Can't create directory %s for the masked reference sequences." % maskRefSeqDir)
                    else:
                        mask_db.maskDb('mr', refSeq, maskRefSeqDir, excludeRefSeqRank, refLeafCladesSet, databaseFile)
                        refSeq = maskRefSeqDir

    # build the general model
    if args.o and ('general' in args.o):
        rankCut = taxonomy_ncbi.TAXONOMIC_RANKS[int(config.get('rankIdCut')) + 1]
        count = generateCladesForGeneralModel(refSeq,
                                              databaseFile,
                                              rankCut,
                                              int(config.get('minGenomesWgs')),
                                              int(config.get('minBpPerSpecies')),
                                              int(config.get('maxLeafClades')),
                                              os.path.join(workingDir, 'ncbids.txt'))
        # the sample specific dir must be empty and the sample specific dir must be set to '' in the PPS config file
        sampleSpecDir = os.path.join(workingDir, 'sampleSpecificDir')
        if os.path.isdir(sampleSpecDir):
            for f in os.listdir(sampleSpecDir):
                os.remove(os.path.join(sampleSpecDir, f))

        print('General Model: the list of "%s" clades at rank "%s" generated.' % (count, rankCut))

    # process output of the marker gene (31 and 16, 23, 5) analysis and create output for PPS
    elif args.o:
        if sequences is None:
            sequences = Sequences(inputFastaFile, inputFastaScaffoldsFile,
                                  scaffoldsToContigsMapFile, taxonomicRanks, minSeqLen)
        taxonomy = Taxonomy(databaseFile, taxonomicRanks)

        # place sequences according to the 16S analysis
        if 's16' in args.o:
            rrna = RRNA16S(config, s16Database, workingDir)
            countList = rrna.setCandidatePlacementFrom16S23S5S(sequences, taxonomy, fastaFileIds)
            print(str('Candidate placement of sequences by the\n16S:' +
                      str(countList[0]) + '\n23S:' + str(countList[1]) + '\n5S:' + str(countList[2]) +
                      '\nall 16S, 23S, and 5S:' + str(countList[3])))

        if 'mg' in args.o:
            mg = MarkerGeneAnalysis(config, mgDatabase, workingDir, mgWorkingDir)
            num = mg.setCandidatePlacement(sequences, taxonomy, fastaFileIds)
            print('Candidate placement of sequences by the marker gene analysis (31): %s' % num)

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
        if 'sc' in args.o or 'otu' in args.o:
            clust = MGCluster(config, mgWorkingDir, fastaFileIds, sequences,
                              taxonomy, os.path.basename(inputFastaFile))
        if 'sc' in args.o:  # construct specific predictions
            clust.refineSpecificPred()

        if 'otu' in args.o:  # build new otus
            clust.reconstructOTU()

        taxonomy.close()

        # fFilePath = common.createTagFilePath(workingDir, fastaFileIds, 'ssd_forbidden')
        # if ('ex' in args.o) and (os.path.isfile(fFilePath)):
        #    forbiddenDict = PPSInput.loadDictFromAFile(fFilePath)
        #    print 'read file of forbidden entries'
        #else:
        forbiddenDict = None

        #exSSDFilePath = os.path.join(workingDir, 'excludeContigsAsSSD.txt')
        #if ('exSSD' in args.o) and ((os.path.isfile(exSSDFilePath))):
        #    exSSDContigNameList = csv.getColumnAsList(exSSDFilePath, entryModifyFunction=None, colNum=0,
        #                                              sep='\n', comment='#')
        #    print 'read file of contigs to be excluded from SSD'
        #else:
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
            refSeq,
            None,  # forbiddenDict
            databaseFile,
            taxonomicRanks,
            100,  # fastaLineMaxChar
            int(config.get('minSSDfileSize')),
            int(config.get('maxSSDfileSize')),
            float(config.get('weightStayAll')),
            summaryTrainFile)

    # generate PPS configuration file (before train or predict) if custom configuration file is not available !
    if (args.t or args.p or args.a) and (config.get('configPPS') is None):
        if parallelPPSmodels:
            pm = 'TRUE'
        else:
            pm = 'FALSE'

        # if the sample specific dir is empty or doesn't exist, it must be set to '' in the PPS configuration file
        sampleSpecDir = os.path.join(workingDir, 'sampleSpecificDir')
        if not os.path.isdir(sampleSpecDir) or len(os.listdir(sampleSpecDir)) == 0:
            sampleSpecDir = ''

        keyDict = {'NCBI_PROCESSED_DIR': str(refSeq),
                   'NCBI_TAX_DIR': os.path.dirname(databaseFile),
                   'PROJECT_DIR': os.path.join(workingDir, 'projectDir'),
                   'TREE_FILE': os.path.join(workingDir, 'ncbids.txt'),
                   'SAMPLE_SPECIFIC_DIR': sampleSpecDir,
                   'PARALLEL_MODELS': pm,  # optional
                   'GENOMES_EXCLUDE': ''}  # optional
        createPPSConfig(ppsConfigFilePath, keyDict)


    # run PhyloPythiaS train (optionally compute training accuracy)
    ppsScripts = os.path.normpath(os.path.join(config.get('ppsInstallDir'), 'scripts'))
    if args.t:

        trainCmd = str(os.path.join(ppsScripts, 'train.rb') + ' ' + ppsConfigFilePath)
        logOut = open(getLogFileName(logDir, 'PPS_train'), 'w')
        trainProc = subprocess.Popen(trainCmd, shell=True, bufsize=-1, cwd=config.get('ppsInstallDir'),
                                     stdout=logOut, stderr=subprocess.STDOUT)
        trainProc.wait()
        logOut.close()
        print('PPS "train" return code: %s' % trainProc.returncode)
        if trainProc.returncode != 0:
            raise Exception("PPS train returned with non-zero %s status: %s" % (trainProc.returncode, trainCmd))

        # compute training accuracy (precision and recall) and confusion matrices
    if args.a:
        try:
            if not os.path.isdir(os.path.join(workingDir, 'projectDir', 'sampled_fasta')):
                print("Train accuracy: Can't compute the training accuracy since the training "
                      "phase hasn't been run yet.")
            else:
                taDir = os.path.join(workingDir, 'train_accuracy')
                if os.path.isdir(taDir):
                    shutil.rmtree(taDir)
                try:
                    os.mkdir(taDir)
                except IOError:
                    print("Train accuracy: Can't create directory: " + taDir)
                else:
                    taOutDir = os.path.join(outputDir, 'train_accuracy')
                    try:
                        if not os.path.isdir(taOutDir):
                            os.mkdir(taOutDir)
                    except IOError:
                        print("Train accuracy: Can't create directory: " + taDir)
                    else:
                        computeTrainingAccuracy(workingDir, taDir, os.path.join(workingDir, 'sampleSpecificDir'),
                                                os.path.join(workingDir, 'projectDir', 'sampled_fasta'), taOutDir,
                                                config.get('ppsInstallDir'), ppsScripts, ppsConfigFilePath,
                                                getLogFileName(logDir, 'PPS_predict_train_data'),
                                                os.path.join(workingDir, 'ncbids.txt'), databaseFile)
                        print("Train accuracy: done")
        except Exception as e:
            print("Training data accuracy won't be computed due to an error: %s" % e.message)

    # run PhyloPythiaS predict
    if args.p:
        if 'c' in args.p:  # predict contigs
            predictCmd = str(os.path.join(ppsScripts, 'predict.rb') + ' ' + fastaFileIds + ' ' + ppsConfigFilePath)
            logOut = open(getLogFileName(logDir, 'PPS_predict_c'), 'w')
            predictProc = subprocess.Popen(predictCmd, shell=True, bufsize=-1, cwd=config.get('ppsInstallDir'),
                                           stdout=logOut, stderr=subprocess.STDOUT)  # stdout=subprocess.STDOUT
            predictProc.wait()
            logOut.close()
            shutil.move(str(fastaFileIds + '.nox.fna.out'), str(fastaFileIds + '.out'))
            shutil.move(str(fastaFileIds + '.nox.fna.PP.out'), str(fastaFileIds + '.PP.out'))
            print('PPS "predict" contigs return code: %s' % predictProc.returncode)
            if predictProc.returncode != 0:
                raise Exception("PPS predict contigs returned with non-zero %s status: %s" %
                                (predictProc.returncode, predictCmd))

        if inputFastaScaffoldsFile is not None:
            scaffoldsIdsNewPath = os.path.join(workingDir, 'crossVal', os.path.basename(fastaFileScaffoldsIds))

        # predict scaffolds
        if 's' in args.p:
            if inputFastaScaffoldsFile is None:
                print('Cannot predict a scaffold file, since it is not specified!')
            else:
                crossValDir = os.path.join(workingDir, 'crossVal')
                if not os.path.isdir(crossValDir):
                    try:
                        os.mkdir(crossValDir)
                    except:
                        print("Can't create directory:", crossValDir)
                        return
                shutil.copy2(fastaFileScaffoldsIds, os.path.join(crossValDir, os.path.basename(fastaFileScaffoldsIds)))
                predictCmd = str(os.path.join(ppsScripts, 'predict.rb') + ' ' + scaffoldsIdsNewPath + ' ' +
                                 ppsConfigFilePath)

                logOut = open(getLogFileName(logDir, 'PPS_predict_s'), 'w')
                predictProc = subprocess.Popen(predictCmd, shell=True, bufsize=-1, cwd=config.get('ppsInstallDir'),
                                               stdout=logOut, stderr=subprocess.STDOUT)  # stdout=subprocess.STDOUT
                predictProc.wait()
                logOut.close()
                shutil.move(str(scaffoldsIdsNewPath + '.nox.fna.out'), str(scaffoldsIdsNewPath + '.out'))
                shutil.move(str(scaffoldsIdsNewPath + '.nox.fna.PP.out'), str(scaffoldsIdsNewPath + '.PP.out'))
                print('PPS "predict" scaffolds return code: %s' % predictProc.returncode)
                if predictProc.returncode != 0:
                    raise Exception("PPS predict scaffolds returned with non-zero %s status: %s" %
                                (predictProc.returncode, predictCmd))

        # compare prediction of scaffolds and contigs
        if 'v' in args.p:
            contigPred = str(fastaFileIds + '.out')
            scaffPred = str(scaffoldsIdsNewPath + '.out')
            if not os.path.isfile(scaffPred):
                print('Comparison of contig/scaffold predictions: scaffolds were not predicted!')
            elif not os.path.isfile(contigPred):
                print('Comparison of contig/scaffold predictions: contigs were not predicted!')
            else:
                cmpOutDir = os.path.join(outputDir, 'contigs_vs_scaff')
                try:
                    os.mkdir(cmpOutDir)
                except IOError:
                    print("Contigs vs. scaffolds: can't create directory: " + cmpOutDir)
                else:
                    # create prediction file for contigs from a prediction file for scaffolds
                    scaffPredAsContigs = os.path.join(workingDir, 'crossVal', 'scaffPredAsContigs.PP.out')
                    out_proc.scafToContigOutput(scaffoldContigMapIdsFile, scaffPred, scaffPredAsContigs)

                    ranks = taxonomy_ncbi.TAXONOMIC_RANKS[1:]
                    cm = confusion_matrix.ConfusionMatrix(fastaFileIds, contigPred,
                                                          scaffPredAsContigs, databaseFile, ranks)
                    for rank in ranks:
                        cm.generateConfusionMatrix(rank, os.path.join(cmpOutDir, 'contigs_vs_scaff'))
                    cm.close()

            createScaffoldAbundanceProfiles(scaffoldsIdsNewPath, scaffPred, profilesDir, databaseFile,
                                            taxonomicRanks, int(config.get('minSeqLen')))

    # Read output from PPS predict and place sequences
    if args.r:
        if sequences is None:
            sequences = Sequences(inputFastaFile, inputFastaScaffoldsFile,
                                  scaffoldsToContigsMapFile, taxonomicRanks, minSeqLen)
        taxonomy = Taxonomy(databaseFile, taxonomicRanks)

        readPPSOutput(sequences, taxonomy, fastaFileIds, eval(config.get('overrideWithPPSPlacements')))
        # one could discard inconsistent contigs here
        taxonomy.close()
        print('Placed sequences after PPS predict: %s' % len(sequences.placedSeqSet))

    predOutFilePath = common.createTagFilePath(outputDir, inputFastaFile, 'pOUT')
    predPPOutFilePath = common.createTagFilePath(outputDir, inputFastaFile, 'PP.pOUT')
    if args.s:
        # create PP.pOUT file
        sequences.writePlacementsPPOut(predPPOutFilePath, taxonomicRanks, config.get('outputFileContigSubPattern'))

        # create .pOUT file (name tab ncbid)
        sequences.writePlacementsOut(predOutFilePath, taxonomicRanks, config.get('outputFileContigSubPattern'))

        sequences.writeSequences(str(common.createTagFilePath(workingDir, inputFastaFile, 'pOUT.fas')),
                                 writeIds=False, outputFileContigSubPattern=config.get('outputFileContigSubPattern'))

        # generate comparison tables
        outCmpRefDir = os.path.join(outputDir, 'cmp_ref')
        try:
            if not os.path.isdir(outCmpRefDir):
                os.mkdir(outCmpRefDir)
        except IOError:
            print("Cmp to ref.: Can't create directory: " + outCmpRefDir)
        else:
            if (referencePlacementFileOut is not None) and os.path.isfile(referencePlacementFileOut):
                ranks = taxonomy_ncbi.TAXONOMIC_RANKS[1:]
                cm = confusion_matrix.ConfusionMatrix(inputFastaFile, predOutFilePath,
                                                      referencePlacementFileOut, databaseFile, ranks)
                for rank in ranks:
                    cm.generateConfusionMatrix(rank, os.path.join(outCmpRefDir, os.path.basename(inputFastaFile)))
                cm.close()

        # compute scaffold-contig consistency
        consOutBuff = csv.OutFileBuffer(common.createTagFilePath(outputDir, inputFastaFile, 'cons'))
        cons = None
        if scaffoldsToContigsMapFile is not None:
            cons = consistency.Consistency(inputFastaFile,
                               str(common.createTagFilePath(outputDir, inputFastaFile, 'pOUT')),
                               scaffoldsToContigsMapFile, databaseFile, minScaffContigCount=None,
                               minScaffBpLen=None, cladesSet=None, considerContigWithNoScaff=True)
        else:
            if os.path.isfile(os.path.normpath(str(fastaFileIds + '.out'))):
                cons = consistency.Consistency(fastaFileIds,
                                   str(fastaFileIds + '.out'),
                                   scaffoldContigMapIdsFile,
                                   databaseFile, minScaffContigCount=None,
                                   minScaffBpLen=None, cladesSet=None, considerContigWithNoScaff=True)

        if cons is not None:
            consOutBuff.writeText(cons.getScaffoldsPrint() + '\n')
            consOutBuff.writeText(cons.getGroupedScaffoldsPrint() + '\n')
            consOutBuff.close()
            cons.close()

        # validation: how the sample specific training data were assigned
        taxonomy = Taxonomy(databaseFile, taxonomicRanks)
        if os.path.isfile(common.createTagFilePath(workingDir, fastaFileIds, 'out')):
            placementPPS = ssd_eval.ppsOut2Placements(common.createTagFilePath(workingDir, fastaFileIds, 'out'), None)
            if os.path.isdir(os.path.join(workingDir, 'sampleSpecificDir')) and \
                    (len(os.listdir(os.path.join(workingDir, 'sampleSpecificDir'))) > 0):
                SSDPlacement = ssd_eval.ssd2Placements(os.path.join(workingDir, 'sampleSpecificDir'), None)
                cmpList = ssd_eval.cmpPlacements(SSDPlacement, placementPPS, taxonomy, taxonomicRanks)
                ssd_eval.cmp2Summary(cmpList, common.createTagFilePath(workingDir, fastaFileIds, 'ssd_cross'))

            # store the sequences that should not be used as sample specific data for respective clades
            # param: rankCut=0~bacteria, maxDist=10 ~ not limited, True ~ output the forbidden list
            #forbiddenList = ssd_eval.filterCmpList(cmpList, 0, 10, taxonomy, True)
            #PPSInput.updateForbiddenList(forbiddenList,
            #                             common.createTagFilePath(workingDir, fastaFileIds, 'ssd_forbidden'))

        # comparison to ref assignment:
        # This block works only if the taxon Ids in the reference are at the species level
        # try:
        #     placementPPSn = ssd_eval.ppsOut2Placements(common.createTagFilePath(outputDir, inputFastaFile, 'pOUT'), None)
        #     if os.path.isfile(referencePlacementFileOut):
        #         placementPPSref = ssd_eval.ppsOut2Placements(referencePlacementFileOut, None)
        #         cmpListR = ssd_eval.cmpPlacements(placementPPSref, placementPPSn, taxonomy, taxonomicRanks)
        #         ssd_eval.cmp2Summary(cmpListR, common.createTagFilePath(workingDir, fastaFileIds, 'cmp_ref'))
        # except AssertionError:
        #     pass

        taxonomy.close()

        # compute the precision and recall values
        if (referencePlacementFileOut is not None) and os.path.isfile(referencePlacementFileOut):
            buff = csv.OutFileBuffer(os.path.join(outputDir, 'precision_recall.csv'))
            buff2 = csv.OutFileBuffer(os.path.join(outputDir, 'precision_recall_corrections.csv'))
            acc = accuracy.Accuracy(inputFastaFile, common.createTagFilePath(outputDir, inputFastaFile, 'pOUT'),
                                    referencePlacementFileOut, databaseFile)
            acc2 = accuracy.Accuracy(inputFastaFile, common.createTagFilePath(outputDir, inputFastaFile, 'pOUT'),
                                    referencePlacementFileOut, databaseFile, float(config.get('correctLabelThreshold')))

            buff.writeText(acc.getAccuracyPrint(taxonomy_ncbi.TAXONOMIC_RANKS[1:],
                                                minFracClade=float(config.get('recallMinFracClade')),
                                                minFracPred=float(config.get('precisionMinFracPred')), overview=True))
            buff2.writeText(acc2.getAccuracyPrint(taxonomy_ncbi.TAXONOMIC_RANKS[1:],
                                                minFracClade=float(config.get('recallMinFracClade')),
                                                minFracPred=float(config.get('precisionMinFracPred')), overview=True))

            buff.close()
            buff2.close()
            acc.close()
            acc2.close()

            # compute precision and recall values for the data that were not used as the sample specific data
            # generates such a fasta file that can be then used for other binning methods!
            ssdDir = os.path.join(workingDir, 'sampleSpecificDir')
            if os.path.isdir(ssdDir):
                ssdFastaFileList = os.listdir(ssdDir)
                if len(ssdFastaFileList) > 0:
                    seqIdList = []
                    # get all sequences in the sample specific directory
                    for ssdFile in ssdFastaFileList:
                        for seqId in fas.fastaFileToDict(os.path.join(ssdDir, ssdFile)):
                            seqIdList.append(seqId)
                    contigIdList = map(lambda x: x.split('_', 1)[1], seqIdList)
                    contigIdToContigName = csv.getMapping(
                        os.path.join(workingDir, str(os.path.basename(inputFastaFile) + '.cToIds')),
                        keyColNum=1, valColNum=0, sep='\t')
                    # transform the seqIds to real contigNames
                    ssdContigNameSet = set(map(lambda x: contigIdToContigName[x][0], contigIdList))
                    out = csv.OutFileBuffer(os.path.join(outputDir, 'no_ssd.fas'))
                    seqIdToBp = {}
                    seqIdToSeq = fas.fastaFileToDict(inputFastaFile)
                    for seqId, seq in seqIdToSeq.iteritems():
                        if seqId not in ssdContigNameSet:
                            out.writeText(str(seqId) + '\n' + str(seq) + '\n')
                            seqIdToBp[seqId] = len(seq)
                    out.close()
                    buff = csv.OutFileBuffer(os.path.join(outputDir, 'precision_recall_no_ssd.csv'))
                    acc = accuracy.Accuracy(seqIdToBp, common.createTagFilePath(outputDir, inputFastaFile, 'pOUT'),
                                            referencePlacementFileOut, databaseFile)
                    buff.writeText(acc.getAccuracyPrint(taxonomy_ncbi.TAXONOMIC_RANKS[1:],
                                                        minFracClade=float(config.get('recallMinFracClade')),
                                                        minFracPred=float(config.get('precisionMinFracPred')),
                                                        overview=True))
                    buff.close()
                    acc.close()

        # generate abundance profiles
        createContigAbundanceProfiles(sequences, profilesDir, int(config.get('minSeqLen')))


def createScaffoldAbundanceProfiles(scaffFilePath, scaffPred, profilesDir, databaseFile, taxonomicRanks, minScaffLen):
    """
        Generates abundance profiles for scaffolds
        @param scaffFilePath:
        @param scaffPred:
        @param databaseFile:
        @param minScaffLen: shorter scaffolds won't be considered
        @return:
    """
    seqIdToBp = fas.getSequenceToBpDict(scaffFilePath)
    seqIdToTaxonId = csv.predToDict(scaffPred)
    taxonomy = Taxonomy(databaseFile, taxonomicRanks)
    errorIdToCount = {}
    entryList = []  # (taxPath, bp)
    for seqId, bp in seqIdToBp.iteritems():
        if bp < minScaffLen:
            continue
        taxonId = seqIdToTaxonId.get(seqId, None)
        if taxonId is None:
            continue

        try:
            taxPath = taxonomy.getPathToRoot(taxonId)
        except Exception:
            if taxonId not in errorIdToCount:
                errorIdToCount[taxonId] = 1
            else:
                errorIdToCount[taxonId] += 1
            continue
    if len(errorIdToCount) > 0:
        for errorId, count in errorIdToCount.iteritems():
            print('Error for taxonId: %s, records skipped: %s' % (errorId, count))

        entryList.append((taxPath, bp))

    taxonomy.close()
    createAbundanceProfiles(entryList, profilesDir, 'scaffold')


def createContigAbundanceProfiles(sequences, profilesDir, minSeqLen):
        """
            Generate abundance profiles for assigned sequences. One file for one rank.
            @param sequences:
            @param profilesDir:
            @param minSeqLen: shorter sequences won't be considered
            @return:
        """
        entryList = []
        for seq in sequences.sequences:
            bp = seq.seqBp
            if bp < minSeqLen:
                continue
            taxPath = seq.getTaxonomyPath()
            entryList.append((taxPath, bp))

        createAbundanceProfiles(entryList, profilesDir, 'contig')


def createAbundanceProfiles(entryList, profilesDir, label):

        for rank in taxonomy_ncbi.TAXONOMIC_RANKS[1:]:

            taxonIdToScientificName = {}
            taxonIdToBp = {}
            taxonIdToCount = {}
            for taxPath, bp in entryList:
                try:
                    if (taxPath is None) or (rank not in taxPath):
                        continue
                except Exception as e:
                    print taxPath
                    raise e

                taxonId = int(taxPath[rank].ncbid)
                if taxonId not in taxonIdToScientificName:
                    taxonIdToScientificName[taxonId] = taxPath[rank].name
                if taxonId in taxonIdToBp:
                    taxonIdToBp[taxonId] += bp
                else:
                    taxonIdToBp[taxonId] = bp
                if taxonId in taxonIdToCount:
                    taxonIdToCount[taxonId] += 1
                else:
                    taxonIdToCount[taxonId] = 1

            infoList = []
            bpSum = 0
            countSum = 0
            for taxonId, bp in taxonIdToBp.iteritems():
                count = taxonIdToCount[taxonId]
                name = taxonIdToScientificName[taxonId]
                infoList.append((taxonId, name, bp, count))
                bpSum += bp
                countSum += count

            infoList.sort(key=lambda x: x[2], reverse=True)

            outFile = csv.OutFileBuffer(os.path.join(profilesDir, str(rank + '_' + label + '_profile.csv')))
            outFile.writeText('bp%, count%, bp, count, taxonId, scientific name\n')

            bpPercentSum = 0.0
            countPercentSum = 0.0
            if countSum > 0:
                assert bpSum > 0
                for taxonId, name, bp, count in infoList:
                    if countSum > 0:
                        assert bpSum > 0
                        bpPercent = round((float(bp) / float(bpSum)) * 100.0, 1)
                        countPercent = round((float(count) / float(countSum)) * 100.0, 1)
                        outFile.writeText(str(bpPercent) + ' %, ' + str(countPercent) + ' %, ' + str(bp) + ', ' +
                                          str(count) + ', ' + str(taxonId) + ', ' + str(name) + '\n')
                        bpPercentSum += bpPercent
                        countPercentSum += countPercent


                outFile.writeText(str(round(bpPercentSum, 0)) + ', ' + str(round(countPercentSum)) + ', ' + str(bpSum) +
                                  ', ' + str(countSum) + ', ,\n')
            outFile.close()


def getLogFileName(logDir, description):
    now = datetime.datetime.now()
    timeStamp = str(
        str(now.year) + str(now.month) + str(now.day) + '_' + str(now.hour) + str(now.minute) + str(now.second))
    fileName = str(description) + '_' + timeStamp + '.txt'
    return os.path.normpath(os.path.join(logDir, fileName))


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
    out.writeText('NCBI_PROCESSED_DIR:%s\n' % keyDict.get('NCBI_PROCESSED_DIR', ''))
    out.writeText(
        """#Directory containing NCBI taxonomy in SQlite3 format with file name "ncbitax_sqlite.db"
#provide empty directory to create new database
#REUSABLE\n""")
    out.writeText('NCBI_TAX_DIR:%s\n' % keyDict.get('NCBI_TAX_DIR', ''))
    out.writeText('#project directory, the directory must be empty\n')
    out.writeText('PROJECT_DIR:%s\n' % keyDict.get('PROJECT_DIR', ''))
    out.writeText(
        """#############################
#!!!FOLLOWING ARE OPTIONAL!!!
#############################
###### Output space options #####
#a file containing a tree in newick format (see restrictions in INSTALL.txt)
#OR a file with ncbi taxon ids (one id per line) to create a tree from\n""")
    out.writeText('TREE_FILE:%s\n' % keyDict.get('TREE_FILE', ''))
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
    out.writeText('SAMPLE_SPECIFIC_DIR:%s\n' % keyDict.get('SAMPLE_SPECIFIC_DIR', ''))
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
CLEAN_UP_TRAIN:FALSE
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
PIE_CHARTS:FALSE
###### Misc options #####
#should the models be built in parallel (please make sure that you have enough number of
processors and main memory)\n""")
    out.writeText('PARALLEL_MODELS:%s\n' % keyDict.get('PARALLEL_MODELS', 'FALSE'))
    out.writeText(
        """#allowed file extensions
EXTENSIONS:
#genomes to exclude: file containing one ncbi tax_id per line\n""")
    out.writeText('GENOMES_EXCLUDE:%s\n' % keyDict.get('GENOMES_EXCLUDE', ''))
    out.writeText(
        """#if the training data is already there then just build models (TRUE/FALSE)
ONLY_MODELS:FALSE\n""")
    out.close()


if __name__ == "__main__":
    main()