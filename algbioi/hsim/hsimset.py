#!/usr/bin/env python

"""
    Copyright (C) 2015  Ivan Gregor

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Note that we could have written some parts of this code in a nicer way,
    but didn't have time. Be careful when reusing the source code.

    ***********************************************************************


    Manages the simulated data generation for the haplotype reconstruction.

"""
import os
import re
# import sys
# import time
import numpy as np
# from Bio.Seq import Seq
# from Bio.Alphabet import generic_dna
from algbioi.com import taxonomy_ncbi
from algbioi.com import csv
from algbioi.com import fasta as fas
from algbioi.com import parallel
from algbioi.com import rand
from algbioi.com import sam
from algbioi.com import fq
from algbioi.hsim import comh
from algbioi.hsim import gene_map
from algbioi.hsim import pfam
from algbioi.haplo import hmain
from algbioi.haplo import heval


def _main():
    """
        Main function, choose steps you want to execute
    """
    assert os.path.isfile(comh.NCBI_TAXONOMY_FILE)
    assert os.path.isdir(comh.REFERENCE_DIR_ROOT)
    taxonomy = taxonomy_ncbi.TaxonomyNcbi(comh.NCBI_TAXONOMY_FILE, considerNoRank=True)

    for spec in comh.SPECIES_LIST:
        specDir = os.path.join(comh.REFERENCE_DIR_ROOT, spec)
        assert comh.isSpeciesDirectory(specDir, taxonomy), "Not a species directory: %s" % specDir

        # get genes of strains containing all the phylo hause keeping genes
        if False:
            createPhyloGeneDir(specDir)

        # compute alignments for all genes or just the genes used for the phylotyping
        if False:
            getPhyloAlignments(specDir, onlyPhyloGenes=True)

        # concatenate phylo genes for clustering
        if False:
            concatenatePhyloGenesForClustering(specDir)

        # cluster concatenated phylo genes
        if False:
            phyloClustering(specDir)

        # define random samples from clusters.
        if False:
            defineRandomSamplesFromClusters(specDir, mode=1)  # mode 0 for testing

        # generate simulated reads for each sample
        if False:
            genSimulatedReadsForSamples(specDir)

        # get error profiles and QS cutoffs
        if False:
            getReadStatAndQSCutoff(specDir)

        # join overlapping pair-end reads
        if False:
            joinPairEndReads(specDir)

        # create a SAM file for the joined pair-end reads
        if False:
            createSamFileForJoinedPairEndReads(specDir)

        # compute error profile for the joined pair-end reads
        if False:
            getJoinedReadsStat(specDir)

        # Get gene mapping
        if False:
            getGeneMapping(specDir)

        if False:
            translateJoinedReadsToProt(specDir)

        # Get Pfam annotation (for joined reads)
        if False:
            getPfamAnnotation(specDir)

        # Get the precision/recall of the pfam-annotations.
        if False:
            getAnnotationAccuracy(specDir)

        # Partition reads into Pfam-domains
        if False:
            partitionReadsToPfamDom(specDir)

        # Get Stat for partitioned reads into Pfam-domains
        if False:
            getStatPartitionedReadsPfamDom(specDir)

        if False:
            assembleContigs(specDir)

        if True:
            computePerBaseAssemblyError(specDir)


        # filter reads given QS cutoffs
        # if False:
        #     filterReadsQS(specDir)


        # getGenesStat(specDir)  # check, put to a new module,


        ###
    taxonomy.close()


def createPhyloGeneDir(specDir):
    """
        Creates directory containing all genes for strains containing all phylo genes.

        The house keeping genes are defined in a file: phylo_genes_*

        Create directory: comh.FASTA_PULL_GENES_PHYLO_DIR_NAME
    """
    print("Getting phylo genes from genomes and draft genomes.")
    # get the file containing phylogenetic gene names (to build a phylogeny / clustering of the genes)
    phyloFile = comh.getPhyloGenesFilePath(specDir)

    # get accessions of genomes and draft genomes containing all the phylo genes
    pullGeneDir = os.path.join(specDir, comh.FASTA_PULL_GENES_DIR_NAME)
    phyloGeneList = csv.getColumnAsList(phyloFile)
    accList = _getAccIntersectAll(pullGeneDir, phyloGeneList)
    print('There are "%s" genomes or draft genomes accessions for phylo genes: "%s"' % (len(accList), phyloGeneList))
    accSet = set(accList)

    # create the destination directory
    dstDir = os.path.join(specDir, comh.FASTA_PULL_GENES_PHYLO_DIR_NAME)
    if not os.path.isdir(dstDir):
        os.mkdir(dstDir)

    # for each gene file, create a file containing only gene sequences of the respective accessions
    for f in os.listdir(pullGeneDir):
        filePath = os.path.join(pullGeneDir, f)
        if os.path.isfile(filePath):
            # sequences containing respective accessions
            seqDict = {}
            for name, seq in fas.fastaFileToDictWholeNames(filePath).iteritems():
                nameDict = comh.getSeqNameToDict(name)
                if nameDict['accession'] in accSet:
                    seqDict[name] = seq
            # store the sequences containing respective accessions to a file
            if len(seqDict) > 0:
                out = csv.OutFileBuffer(os.path.join(dstDir, f))
                first = True
                for k, v in seqDict.iteritems():
                    if first:
                        out.writeText('>%s\n%s' % (k, v))
                        first = False
                    else:
                        out.writeText('\n>%s\n%s' % (k, v))
                out.close()


def _getAccIntersectAll(genesDir, geneList):
    """
        Get a list of all accessions of genomes or draft genomes that contain all genes from the gene list.
        S.t. each genome or draft genome, the accession of which is returned by this function, contain all
        the genes from the gene list

        @param genesDir: directory containing all the genes (one file for one gene)
        @param geneList: list of gene names

        @return: a list of accessions
        @rtype: list
    """
    # get the file names of the phylo genes
    geneSet = set(map(lambda x: comh.getGeneNameToFileName(x), geneList))

    # get mapping: accession -> list of the phylo gene files
    accDict = {}
    for f in os.listdir(genesDir):
        if f in geneSet:
            for name in fas.fastaFileToDictWholeNames(os.path.join(genesDir, f)).keys():
                d = comh.getSeqNameToDict(name)
                assert comh.getGeneNameToFileName(d['geneName']) == f, "%s %s" % (d['geneName'], f)
                acc = d['accession']
                if acc not in accDict:
                    accDict[acc] = [f]
                else:
                    accDict[acc].append(f)

    # list of the accessions containing all the phylo genes
    allGeneAcc = []
    geneSetLen = len(geneSet)
    for acc, geneL in accDict.iteritems():
        l = len(set(geneL))
        if l == geneSetLen:
            allGeneAcc.append(acc)
        else:
            assert l < geneSetLen
    return allGeneAcc


def getPhyloAlignments(specDir, onlyPhyloGenes=True):
    """
        Computing alignments for all phylo genes.

        @param onlyPhyloGenes: if true, only files used for the phylogenetic reconstruction / clustering will be used

        Create directory: comh.FASTA_PULL_GENES_PHYLO_ALIGN_DIR_NAME
    """
    print("Computing alignments for the phylo genes")
    # create the destination directory for the aligned sequences
    pullAlignDir = os.path.join(specDir, comh.FASTA_PULL_GENES_PHYLO_ALIGN_DIR_NAME)
    if not os.path.isdir(pullAlignDir):
        os.mkdir(pullAlignDir)

    fileNameList = None
    if onlyPhyloGenes:
        # get gene names used of the phylo reconstruction / clustering
        fileNameList = map(lambda x: comh.getGeneNameToFileName(x),
                           csv.getColumnAsList(comh.getPhyloGenesFilePath(specDir)))

    # compute alignments
    comh.getAlignments(os.path.join(specDir, comh.FASTA_PULL_GENES_PHYLO_DIR_NAME),
                       os.path.join(specDir, comh.FASTA_PULL_GENES_PHYLO_ALIGN_DIR_NAME),
                       listOfInFileNames=fileNameList, reportFailedCmd=True)
    #


def concatenatePhyloGenesForClustering(specDir):
    """
        Concatenate phylo genes for the phylo clustering / reconstruction for all strains.

        Create file: comh.CLUSTER_PHYLO_STRAINS_FNA
        Create directory: comh.CLUSTER_PHYLO_DIR
    """
    print('Concatenating phylo genes.')
    # list of file paths of the phylo genes
    geneFilePathList = map(lambda x: os.path.join(specDir, comh.FASTA_PULL_GENES_PHYLO_ALIGN_DIR_NAME,
                                                  comh.getGeneNameToFileName(x)),
                           csv.getColumnAsList(comh.getPhyloGenesFilePath(specDir)))

    # map: accession -> concatenated sequence
    accToSeq = {}
    geneNameSet = set()
    # concatenate gene sequences
    for filePath in geneFilePathList:
        assert os.path.isfile(filePath)
        for name, seq in fas.fastaFileToDictWholeNames(filePath).iteritems():
            d = comh.getSeqNameToDict(name)
            acc = d['accession']
            geneName = d['geneName']
            geneNameSet.add(geneName)
            if acc in accToSeq:
                accToSeq[acc] = accToSeq[acc] + seq
            else:
                accToSeq[acc] = seq

    # checking the gene names
    s1 = set(map(lambda x: os.path.basename(x), geneFilePathList))
    s2 = set(map(lambda x: comh.getGeneNameToFileName(x), geneNameSet))
    assert s1 == s2, 'Different gene sets: %s %s' % (s1, s2)

    # create a directory for the gene clustering
    clusterDir = os.path.join(specDir, comh.CLUSTER_PHYLO_DIR)
    if not os.path.isdir(clusterDir):
        os.mkdir(clusterDir)

    # create the fasta file containing concatenated phylo genes
    out = csv.OutFileBuffer(os.path.join(clusterDir, comh.CLUSTER_PHYLO_STRAINS_FNA))
    lengthSet = set()
    first = True
    for name, seq in accToSeq.iteritems():
        lengthSet.add(len(seq))
        if first:
            out.writeText(">%s\n%s" % (name, seq))
            first = False
        else:
            out.writeText("\n>%s\n%s" % (name, seq))
    out.close()
    assert len(lengthSet) == 1, 'all sequences should be of the same length and not: %s' % lengthSet
    print("Alignment length: %s" % lengthSet)


def phyloClustering(specDir):
    """
        Cluster concatenated phylo genes using mothur, furthest neighbor.
    """
    print('Clustering concatenated phylo genes.')
    # create the mothur output directory
    cwd = os.path.join(specDir, comh.CLUSTER_PHYLO_DIR, comh.CLUSTER_MOTHUR_CWD)
    if not os.path.isdir(cwd):
        os.mkdir(cwd)

    # concatenated alignments
    alignFile = os.path.join(specDir, comh.CLUSTER_PHYLO_DIR, comh.CLUSTER_PHYLO_STRAINS_FNA)
    taskList = []
    # compute distance matrix
    taskList.append(parallel.TaskCmd('time %s "#dist.seqs(fasta=%s, processors=%s, countends=T, calc=onegap, cutoff=0.3, output=lt)"'
                % (comh.MOTHUR_BINARY, alignFile, comh.MAX_PROC), cwd))
    # dstance matrix file
    distFile = alignFile.rsplit('.', 1)[0] + '.phylip.dist'

    # clustering
    taskList.append(parallel.TaskCmd('time %s "#cluster(phylip=%s, method=%s, precision=1000)"'
                                     % (comh.MOTHUR_BINARY, distFile, comh.CLUSTER_MOTHUR_METHOD), cwd))
    # run commands
    parallel.reportFailedCmd(parallel.runCmdSerial(taskList))

    # mothurCmd = str('time /Users/ivan/Documents/work/tools/mothur/mothur/mothur "#dist.seqs(fasta=/Users/ivan/Documents/nobackup/hsim01/562/phylo_strains_cluster.fna, processors=2, countends=T, calc=onegap, cutoff=0.3, output=lt)" ')
    # mothurCmd = str('time /Users/ivan/Documents/work/tools/mothur/mothur/mothur "#cluster(phylip=/Users/ivan/Documents/nobackup/hsim01/562/phylo_strains_cluster.phylip.dist, method=nearest, precision=1000)"')


def defineRandomSamplesFromClusters(specDir, mode=1):
    """
        Define random samples from clusters.

        Samples definitions file: identifier(int) TAB comma concatenated accessions TAB comma concatenated abundances

        Mode 1: Generate samples from each cluster at a particular defined threshold

        Mode 0: for testing, take three strains from the first unique cluster that has at least three strains
        (take first genome found, first draft genome found, and first strain that is not equal to the first two taken).

        Creates a directory and a file with sample definitions: comh.SAMPLES_DIR and comh.SAMPLES_DEF_FILE

        @param mode: a strategy of choosing a sample (1 or 0 for testing)
    """
    print('Define random samples from clusters.')
    # define the random generator
    randGen = np.random.RandomState(comh.SAMPLES_DEF_RAND_SEED)
    # verifies that the clustering file is there and the clustering method was right
    assert comh.CLUSTER_MOTHUR_METHOD == 'furthest'
    clusterListFile = os.path.join(specDir, comh.CLUSTER_PHYLO_DIR,
                                   comh.CLUSTER_PHYLO_STRAINS_FNA.split('.', 1)[0] + '.phylip.fn.list')
    assert os.path.isfile(clusterListFile)
    assert mode in {0, 1}

    # create the samples directory
    samplesDir = os.path.join(specDir, comh.SAMPLES_DIR)
    if not os.path.isdir(samplesDir):
        os.mkdir(samplesDir)

    # list of defined samples
    sampleList = []

    if mode == 0:
        # creates a simple test sample with three strains
        for line in open(clusterListFile):
            # find line containing unique clusters
            if re.match(r'^unique', line):
                # for each unique cluster
                for c in line.split()[2:]:
                    genome = None
                    draftGenome = None
                    thirdStrain = None
                    cList = c.split(',')
                    # cluster has at least three elements
                    if len(cList) >= 3:
                        # get genome
                        for acc in cList:
                            if os.path.isfile(os.path.join(specDir, comh.FASTA_GENOMES_DIR_NAME, acc + '.fna')):
                                genome = acc
                                break
                        # get draft genome
                        for acc in cList:
                            if os.path.isfile(os.path.join(specDir, comh.FASTA_GENOMES_DRAFT_DIR_NAME, acc + '.fna')):
                                draftGenome = acc
                                break
                        # get third strain
                        for acc in cList:
                            if acc != genome and acc != draftGenome:
                                thirdStrain = acc
                                break
                    # define the sample
                    if genome is not None and draftGenome is not None and thirdStrain is not None:
                        sampleList.append([genome, draftGenome, thirdStrain])
                        break
                break
    elif mode == 1:
        # create many samples exploring each cluster at a defined threshold
        thresholdLine = None
        lastLine = None
        # find the line defining clusters at the required threshold
        for line in open(clusterListFile):
            line = line.strip()
            # the clustering on this line correspond to this threshold
            try:
                threshold = float(line.split('\t', 1)[0])
                # print line.split('\t', 1)[0], map(lambda x: len(x.split(',')), line.strip().split('\t')[2:])
            except ValueError:
                continue
            # this line correspond to the threshold
            if abs(threshold - comh.SAMPLES_DEF_IDENTITY_CUTOFF) < 0.000001:
                thresholdLine = line
                break
            # it may happen that the threshold is not defined, then the previous smaller is taken
            elif threshold > comh.SAMPLES_DEF_IDENTITY_CUTOFF:
                thresholdLine = lastLine
                break
            lastLine = line
        # the clustering with the required threshold has been found
        if thresholdLine is not None:
            # generate samples for each cluster
            for cluster in thresholdLine.split('\t')[2:]:
                samples = rand.getRandSubsets(cluster.split(','), comh.SAMPLES_DEF_STRAIN_COUNT_LIST,
                                              comh.SAMPLES_DEF_MAX_SAMPLES_PER_CLUSTER, rand=randGen)
                # add the samples to the resulting list
                if len(samples) > 0:
                    sampleList += samples
    # unknown mode
    else:
        assert False

    # no samples has been generated
    if len(sampleList) == 0:
        return

    # generate random strain abundances for the samples
    abundanceList = rand.getRandLognormNumbers(map(lambda x: len(x), sampleList),
                                               mean=comh.SAMPLES_DEF_LOGNORM_MEAN,
                                               sd=comh.SAMPLES_DEF_LOGNORM_SD,
                                               minVal=comh.SAMPLES_DEF_MIN_COVERAGE,
                                               maxVal=comh.SAMPLES_DEF_MAX_COVERAGE,
                                               rand=randGen)
    # store samples to a file
    out = csv.OutFileBuffer(os.path.join(specDir, comh.SAMPLES_DIR, comh.SAMPLES_DEF_FILE))
    counter = 0
    for sampleAccessions, sampleAbundances in zip(sampleList, abundanceList):
        out.writeText('%s\t%s\t%s\n' % (counter, ','.join(sampleAccessions),
                                        ','.join(map(lambda x: str(x), sampleAbundances))))
        counter += 1
    out.close()
    # return (sampleList, abundanceList) # @return: tuple (a list of lists of accessions, a list of lists of abundances)


def genSimulatedReadsForSamples(specDir):
    """
        Generate simulated reads for each sample.
    """
    print('Generating simulated reads.')
    taskList = []
    # define tasks for each sample
    for sampleDef in open(os.path.join(specDir, comh.SAMPLES_DIR, comh.SAMPLES_DEF_FILE)):
        # get a sample definition
        sampleId, accList, abundanceList = sampleDef.split('\t')
        accList = accList.split(',')
        abundanceList = map(lambda x: float(x.strip()), abundanceList.split(','))
        assert len(accList) == len(abundanceList)
        # create the sample directory
        sampleDir = os.path.join(specDir, comh.SAMPLES_DIR, str(sampleId))
        if not os.path.isdir(sampleDir):
            os.mkdir(sampleDir)

        # for each strain, define a run of the simulator
        for acc, coverage in zip(accList, abundanceList):
            # create a directory for this strain
            strainSimDir = os.path.join(sampleDir, acc)
            if not os.path.isdir(strainSimDir):
                os.mkdir(strainSimDir)
            # get the strain fasta file
            srcFna = None
            genomeFastaFile = os.path.join(specDir, comh.FASTA_GENOMES_DIR_NAME, acc + '.fna')
            draftGenomeFastaFile = os.path.join(specDir, comh.FASTA_GENOMES_DRAFT_DIR_NAME, acc + '.fna')
            if os.path.isfile(genomeFastaFile):
                srcFna = genomeFastaFile
            elif os.path.isfile(draftGenomeFastaFile):
                srcFna = draftGenomeFastaFile
            assert srcFna is not None

            # there should be data for all libraries
            assert len(comh.ART_READ_LEN) == len(comh.ART_INSERT_SIZE) == len(comh.ART_INSERT_SD) == \
                len(comh.ART_MIN_SEQ_LEN)

            count = 0
            for readLen, insertSize, insertSD, minSeqLen, qProfile in zip(comh.ART_READ_LEN, comh.ART_INSERT_SIZE,
                                                                comh.ART_INSERT_SD, comh.ART_MIN_SEQ_LEN,
                                                                comh.ART_Q_PROFILE):
                # filter out short sequences from the (draft) genome fasta files
                srcNoShortSeq = os.path.join(strainSimDir, str(count) + '_' + os.path.basename(srcFna))
                fas.cpSeqNoShortSeq(srcFna, srcNoShortSeq, minSeqLen)

                # use a particular error profile
                qP1, qP2 = qProfile
                if qP1 is not None:
                    qP1 = ' -1 %s ' % qP1
                else:
                    qP1 = ''
                if qP2 is not None:
                    qP2 = ' -2 %s ' % qP2
                else:
                    qP2 = ''

                # define the read generation command
                cmd = "%s -i %s -f %s -rs %s -l %s -m %s -s %s -p -o %s_pair -na -sam %s %s; " \
                      "gzip %s_pair1.fq; gzip %s_pair2.fq; gzip %s_pair.sam; gzip %s" \
                      % (comh.ART_ILLUMINA_BINARY, srcNoShortSeq, coverage,
                         rand.strToRandInt(str(comh.SAMPLES_DEF_RAND_SEED) + acc + str(round(coverage))
                                           + str(os.path.getsize(srcNoShortSeq))), readLen, insertSize, insertSD,
                         count, qP1, qP2, count, count, count, srcNoShortSeq)

                # define the logging file for stdout
                stdoutLog = os.path.join(strainSimDir, '%s_art_log.txt' % count)
                taskList.append(parallel.TaskCmd(cmd, strainSimDir, stdout=stdoutLog))
                # print cmd
                count += 1

    # generate all simulated reads, run the simulator in parallel
    parallel.reportFailedCmd(parallel.runCmdParallel(taskList, maxProc=comh.MAX_PROC))


def getReadStatAndQSCutoff(specDir):
    """
        Get read statistics and QS cutoff
    """
    print('Getting read statistics and QS cutoff')
    # list of SAM and reference FASTA files to compute the statistics
    fileTupleList = []
    # directory containing all samples
    samplesDir = os.path.join(specDir, comh.SAMPLES_DIR)
    # output profile file
    outFilePathProfile = os.path.join(samplesDir, comh.SAMPLES_ERROR_PROFILE)
    # output QS cutoff file
    outFilePathQSCutoff = os.path.join(samplesDir, comh.SAMPLES_ERROR_QS_CUTOFF)

    # collect all files to compute the statistics and QS cutoff
    for sample in os.listdir(samplesDir):
        if sample.isdigit():
            sampleDir = os.path.join(samplesDir, sample)
            if os.path.isdir(sampleDir):
                # for all strains of one sample
                for strain in os.listdir(sampleDir):
                    strainDir = os.path.join(sampleDir, strain)
                    if os.path.isdir(strainDir):
                        # get FASTA and SAM files for all libraries of one strain
                        fastaD = {}
                        samD = {}
                        for f in os.listdir(strainDir):
                            i = f.split('_', 1)[0]
                            if i.isdigit():
                                if f.endswith('.fna.gz'):
                                    fastaD[int(i)] = os.path.join(strainDir, f)
                                elif f.endswith('.sam.gz'):
                                    samD[int(i)] = os.path.join(strainDir, f)
                        assert len(samD) == len(fastaD)
                        for i, fastaFilePath in fastaD.iteritems():
                            samFilePath = samD[i]
                            # store the files and param (SAM file path, reference FASTA file path, readLen, qsArrayLen)
                            fileTupleList.append(
                                (samFilePath, fastaFilePath, comh.ART_READ_LEN[i], comh.ART_QS_MAX[i]))
    # compute the statistics (profile and QS cutoff)
    sam.getErrorStatAndQSCutoff(fileTupleList, outFilePathProfile, outFilePathQSCutoff, maxCpu=comh.MAX_PROC)


def joinPairEndReads(specDir):
    """
        Join overlapping pair end reads.
    """
    print('Joining overlapping pair end reads.')

    # collect all libraries with overlapping reads
    librarySet = set()
    for i in range(len(comh.ART_INSERT_SD)):
        if (comh.ART_INSERT_SIZE[i] + comh.ART_INSERT_SD[i]) < 2 * comh.ART_READ_LEN[i]:
            librarySet.add(i)

    # collect files for the joining
    fileTupleList = []

    # directory containing all samples
    samplesDir = os.path.join(specDir, comh.SAMPLES_DIR)

    # collect all files to join the reads
    for sample in os.listdir(samplesDir):
        if sample.isdigit():
            sampleDir = os.path.join(samplesDir, sample)
            if os.path.isdir(sampleDir):
                # for all strains of one sample
                for strain in os.listdir(sampleDir):
                    strainDir = os.path.join(sampleDir, strain)
                    if os.path.isdir(strainDir):
                        # get fq1 and fq2 files for all libraries of one strain
                        fq1D = {}
                        fq2D = {}
                        for f in os.listdir(strainDir):
                            i = f.split('_', 1)[0]
                            if i.isdigit() and (int(i) in librarySet):
                                if f.endswith('pair1.fq.gz'):
                                    fq1D[int(i)] = os.path.join(strainDir, f)
                                elif f.endswith('pair2.fq.gz'):
                                    fq2D[int(i)] = os.path.join(strainDir, f)
                        assert len(fq1D) == len(fq2D)
                        for i, fq1Path in fq1D.iteritems():
                            fq2Path = fq2D[i]
                            fqJoinPath = os.path.join(os.path.dirname(fq1Path), '%s_join.fq.gz' % i)
                            # store the files and param (fq1, fq2, fqJoin, readLen, insert, sd, qsMax)
                            fileTupleList.append(
                                (fq1Path, fq2Path, fqJoinPath, comh.ART_READ_LEN[i], comh.ART_INSERT_SIZE[i],
                                 comh.ART_INSERT_SD[i], comh.ART_QS_MAX[i]))
    # join pair end reads
    r = fq.joinPairEnd(fileTupleList,
                       minOverlap=comh.SAMPLES_PAIRED_END_JOIN_MIN_OVERLAP,
                       minOverlapIdentity=comh.SAMPLES_PAIRED_END_JOIN_MIN_OVERLAP_IDENTITY,
                       maxCpu=comh.MAX_PROC)
    print("Not joined: %s %%" % r)


def createSamFileForJoinedPairEndReads(specDir):
    """
        Create SAM files for the joined pair-end reads.
    """
    print('Creating SAM files for the joined pair-end reads')
    # all files for the SAM files creation
    fileTupleList = []

    # directory containing all samples
    samplesDir = os.path.join(specDir, comh.SAMPLES_DIR)

    # collect all files for the SAM files creation
    for sample in os.listdir(samplesDir):
        if sample.isdigit():
            sampleDir = os.path.join(samplesDir, sample)
            if os.path.isdir(sampleDir):
                # for all strains of one sample
                for strain in os.listdir(sampleDir):
                    strainDir = os.path.join(sampleDir, strain)
                    if os.path.isdir(strainDir):
                        for f in os.listdir(strainDir):
                            i, name = f.split('_', 1)
                            if i.isdigit() and name.startswith('join'):
                                fqJoinPath = os.path.join(strainDir, f)
                                pairEndSamPath = os.path.join(strainDir, '%s_pair.sam.gz' % i)
                                joinSamPath = os.path.join(strainDir, "%s_join.sam.gz" % i)
                                assert os.path.isfile(pairEndSamPath)
                                fileTupleList.append((fqJoinPath, pairEndSamPath, joinSamPath))
    # create SAM files for all joined pair end reads
    r = sam.createSamFileForJoinedPairEndReads(fileTupleList, maxCpu=comh.MAX_PROC)
    print('SAM files for "%s" joined pair end reads created' % r)


def getJoinedReadsStat(specDir):
    """
        Compute error profile for the joined pair-end reads.

    """
    print("Computing error profile for the joined pair-end reads.")
    # all files for the error profile computation
    fileTupleList = []

    # directory containing all samples
    samplesDir = os.path.join(specDir, comh.SAMPLES_DIR)

    # collect all files for the error profile computation (refPath, fqJoinedPath, samPath, readLen)
    for sample in os.listdir(samplesDir):
        if sample.isdigit():
            sampleDir = os.path.join(samplesDir, sample)
            if os.path.isdir(sampleDir):
                # for all strains of one sample
                for strain in os.listdir(sampleDir):
                    strainDir = os.path.join(sampleDir, strain)
                    if os.path.isdir(strainDir):
                        for f in os.listdir(strainDir):
                            i, name = f.split('_', 1)
                            if i.isdigit() and name.endswith('.fna.gz'):
                                refPath = os.path.join(strainDir, f)
                                fqJoinPath = os.path.join(strainDir, '%s_join.fq.gz' % i)
                                if os.path.isfile(fqJoinPath):
                                    joinSamPath = os.path.join(strainDir, '%s_join.sam.gz' % i)
                                    assert os.path.isfile(joinSamPath)
                                    fileTupleList.append((refPath, fqJoinPath, joinSamPath,
                                                          comh.ART_READ_LEN[int(i)], comh.ART_QS_MAX[int(i)]))
    # compute error profiles
    rc = sam.getJoinedReadStat(fileTupleList, os.path.join(samplesDir, comh.SAMPLES_JOIN_ERROR_PROFILE), comh.MAX_PROC)
    print('Statistics for "%s" joined pair-end reads counted' % rc)


def getGeneMapping(specDir):
    """
        Getting gene mapping for the reads.
    """
    print('Getting gene mapping for the reads.')

    # directory containing all genes
    genesDir = os.path.join(specDir, comh.FASTA_PULL_GENES_DIR_NAME)

    # directory containing all samples
    samplesDir = os.path.join(specDir, comh.SAMPLES_DIR)

    # collect all files: FNA references for the joined pair-end reads and the corresponding SAM files
    refList = []  # list of reference FNA files
    samList = []  # list of SAM file tuples (input SAM, output SAM)
    for sample in os.listdir(samplesDir):
        if sample.isdigit():
            sampleDir = os.path.join(samplesDir, sample)
            if os.path.isdir(sampleDir):
                # for all strains of one sample
                for strain in os.listdir(sampleDir):
                    strainDir = os.path.join(sampleDir, strain)
                    if os.path.isdir(strainDir):
                        for f in os.listdir(strainDir):
                            i, name = f.split('_', 1)
                            if i.isdigit() and name.endswith('.fna.gz'):
                                samPath = os.path.join(strainDir, str(i) + '_join.sam.gz')
                                if os.path.isfile(samPath):
                                    outSamPath = os.path.join(strainDir, str(i) + "_join_gmap.sam.gz")
                                    samList.append((samPath, outSamPath))
                                    refList.append(os.path.join(strainDir, f))
    gene_map.getReadsToGenesMap(refList, samList, genesDir)


def translateJoinedReadsToProt(specDir):
    """
        Translate reads to all 6 reading frames.

        TODO: not implemented yet !
    """
    print('Translating joined reads into 6 reading frames')
    # directory containing all samples
    samplesDir = os.path.join(specDir, comh.SAMPLES_DIR)

    # collect tasks
    taskList = []
    for sample in os.listdir(samplesDir):
        if sample.isdigit():
            sampleDir = os.path.join(samplesDir, sample)
            if os.path.isdir(sampleDir):
                # for all strains of one sample
                for strain in os.listdir(sampleDir):
                    strainDir = os.path.join(sampleDir, strain)
                    if os.path.isdir(strainDir):
                        for f in os.listdir(strainDir):
                            i, name = f.split('_', 1)
                            if i.isdigit() and name == 'join.fq.gz':
                                inFq = os.path.join(strainDir, f)
                                outFasta = os.path.join(strainDir, str(i) + '_join_prot.fna.gz')
                                taskList.append(parallel.TaskThread(fq.readsToProt, (inFq, outFasta)))
    # run tasks in parallel
    parallel.runThreadParallel(taskList, comh.MAX_PROC)


def getPfamAnnotation(specDir):
    """
        Get the Pfam annotation for all the joined pair-end reads.
    """
    print('Searching through PROT read sequences')
    # directory containing all samples
    samplesDir = os.path.join(specDir, comh.SAMPLES_DIR)

    # collect tasks
    taskList = []
    for sample in os.listdir(samplesDir):
        if sample.isdigit():
            sampleDir = os.path.join(samplesDir, sample)
            if os.path.isdir(sampleDir):
                # for all strains of one sample
                for strain in os.listdir(sampleDir):
                    strainDir = os.path.join(sampleDir, strain)
                    if os.path.isdir(strainDir):
                        for f in os.listdir(strainDir):
                            i, name = f.split('_', 1)
                            if i.isdigit() and name == 'join_prot.fna.gz':
                                protReadsInGzip = os.path.join(strainDir, f)
                                assert protReadsInGzip.endswith('fna.gz')
                                protReadsIn = protReadsInGzip[:-3]
                                domOut = os.path.join(strainDir, str(i) + '_join_prot.domtblout')
                                cmd = 'zcat %s > %s;%s/hmmsearch -o /dev/null --noali --domtblout %s -E 0.01 ' \
                                      '--cpu 1 %s %s;rm %s;gzip -f %s' \
                                      % (protReadsInGzip, protReadsIn, comh.HMMER_BINARY, domOut,
                                         os.path.join(comh.PFAM, 'Pfam-A.hmm'), protReadsIn, protReadsIn, domOut)
                                cwd = os.path.dirname(protReadsIn)
                                taskList.append(parallel.TaskCmd(cmd, cwd))

    # search reads again the Pfam database in parallel
    parallel.runCmdParallel(taskList, comh.MAX_PROC)


def getAnnotationAccuracy(specDir):
    """
        Getting the precision/recall of the Pfam annotations.
    """
    print('Getting the precision and recall of the Pfam read annotations.')
    # directory containing all samples
    samplesDir = os.path.join(specDir, comh.SAMPLES_DIR)

    # collect tasks
    taskList = []
    for sample in os.listdir(samplesDir):
        if sample.isdigit():
            sampleDir = os.path.join(samplesDir, sample)
            if os.path.isdir(sampleDir):
                # for all strains of one sample
                for strain in os.listdir(sampleDir):
                    strainDir = os.path.join(sampleDir, strain)
                    if os.path.isdir(strainDir):
                        for f in os.listdir(strainDir):
                            i, name = f.split('_', 1)
                            if i.isdigit() and name == 'join_prot.domtblout.gz':
                                taskList.append(parallel.TaskThread(pfam.getHmmAnnotationAccuracy,
                                    (os.path.join(strainDir, f),
                                     os.path.join(specDir, comh.FASTA_PULL_GENES_PHYLO_PFAM_TO_GENE_NAME),
                                     comh.SAMPLES_PFAM_EVAN_MIN_SCORE,
                                     comh.SAMPLES_PFAM_EVAN_MIN_ACCURACY)))

    rl = parallel.runThreadParallel(taskList, comh.MAX_PROC)
    pfam.mergeResultsHmmAnnotationAccuracy(rl, os.path.join(samplesDir, comh.SAMPLES_JOIN_PFAM_MAP_QUALITY),
                                           comh.SAMPLES_PFAM_EVAN_MIN_SCORE, comh.SAMPLES_PFAM_EVAN_MIN_ACCURACY)


def partitionReadsToPfamDom(specDir):
    """
        Partition reads into Pfam-domains.
    """
    print('Partitioning reads into the individual Pfam-domains.')
    samplesDir = os.path.join(specDir, comh.SAMPLES_DIR)

    # collect tasks
    taskList = []
    for sample in os.listdir(samplesDir):
        if sample.isdigit():
            sampleDir = os.path.join(samplesDir, sample)
            if os.path.isdir(sampleDir):
                taskList.append(parallel.TaskThread(pfam.partitionReads, (sampleDir, comh.SAMPLES_PFAM_EVAN_MIN_SCORE,
                                    comh.SAMPLES_PFAM_EVAN_MIN_ACCURACY, comh.SAMPLES_SHUFFLE_RAND_SEED,
                                    comh.SAMPLES_PFAM_PARTITIONED_DIR)))
    parallel.runThreadParallel(taskList, comh.MAX_PROC)


def getStatPartitionedReadsPfamDom(specDir):
    """
        Get statistics for the partitioned reads into the Pfam-domains
    """
    print('Getting statistics for the partitioned reads into Pfam-domains.')
    samplesDir = os.path.join(specDir, comh.SAMPLES_DIR)

    # collect tasks
    taskList = []
    for sample in os.listdir(samplesDir):
        if sample.isdigit():
            sampleDir = os.path.join(samplesDir, sample)
            if os.path.isdir(sampleDir):
                taskList.append(parallel.TaskThread(pfam.getStatPartitionedReads,
                                                    (sampleDir,
                                                     comh.SAMPLES_PFAM_PARTITIONED_DIR,
                                                     comh.SAMPLES_PFAM_PARTITIONED_STAT_FILE)))

    parallel.runThreadParallel(taskList, comh.MAX_PROC)


def assembleContigs(specDir):
    """
        Assemble contigs of each simulated dataset.
    """
    print('Running the contig assembly.')
    samplesDir = os.path.join(specDir, comh.SAMPLES_DIR)

    # collect tasks
    taskList = []
    for sample in os.listdir(samplesDir):
        if sample.isdigit():
            partitionedDir = os.path.join(samplesDir, sample, comh.SAMPLES_PFAM_PARTITIONED_DIR)
            if os.path.isdir(partitionedDir):
                for f in os.listdir(partitionedDir):
                    fPath = os.path.join(partitionedDir, f)
                    if f.endswith('fq.gz') and os.path.isfile(fPath):
                        base = fPath[:-6]
                        inFq = fPath
                        inDomtblout = '%s_prot.domtblout.gz' % base
                        inProtFna = '%s_prot.fna.gz' % base
                        outPath = '%s_read_rec.pkl.gz' % base
                        taskList.append(parallel.TaskThread(hmain.buildSuperReads, (inFq, inDomtblout,
                            comh.ART_READ_LEN[0], inProtFna, outPath, comh.TRANSLATION_TABLE,
                            comh.ASSEMBLY_MAX_MISMATCH_QS_ALLOWED, comh.ASSEMBLY_MIN_SCORE_REQIURED,
                            comh.ASSEMBLY_MIN_ANNOT_OVERLAP_SCORE, comh.ASSEMBLY_SCORE_STOP_SEARCH,
                            comh.ASSEMBLY_MAX_QS)))

    parallel.runThreadParallel(taskList, comh.MAX_PROC)


def computePerBaseAssemblyError(specDir):
    """
        Compute the per base assembly error.
    """
    print('Computing the per-base assembly error')
    samplesDir = os.path.join(specDir, comh.SAMPLES_DIR)

    # buffer the reference sequences
    refSeqBuff = fas.getSequenceBuffer([os.path.join(specDir, comh.FASTA_GENOMES_DIR_NAME),
                           os.path.join(specDir, comh.FASTA_GENOMES_DRAFT_DIR_NAME)])
    # collect tasks
    taskList = []
    # rList = []
    for sample in os.listdir(samplesDir):
        if sample.isdigit():
            partitionedDir = os.path.join(samplesDir, sample, comh.SAMPLES_PFAM_PARTITIONED_DIR)
            if os.path.isdir(partitionedDir):
                for f in os.listdir(partitionedDir):
                    fPath = os.path.join(partitionedDir, f)
                    if f.endswith('join_read_rec.pkl.gz') and os.path.isfile(fPath):
                        base = fPath[:-21]
                        readRecPkl = fPath
                        gmapSam = '%s_join_gmap.sam.gz' % base
                        assert os.path.isfile(gmapSam)

                        taskList.append(parallel.TaskThread(heval.getPerBaseErrorPkl,
                                                            (readRecPkl, gmapSam, refSeqBuff, 30,
                                                             comh.TRANSLATION_TABLE)))

                        # a = heval.getPerBaseErrorPkl(readRecPkl, gmapSam, refGenomesDraftGDirList, 30, comh.TRANSLATION_TABLE)
                        # rList.append(a)
                        # heval.getAssemblyReport([a], maxCov=30)
    rList = parallel.runThreadParallel(taskList, comh.MAX_PROC)
    report = heval.getAssemblyReport(rList, maxCov=30)
    print report
    out = csv.OutFileBuffer(os.path.join(samplesDir, comh.ASSEMBLY_SUPER_READ_EVAL_INIT))
    out.writeText(report)
    out.close()


def filterReadsQS(specDir):
    """
        Filter reads according to the QS cutoffs
        @attention: TODO: not ready yet !
    """
    assert False  # TODO: not ready yet !
    print("Filtering reads according to the QS cutoffs")

    # get the QS cutoff arrays for each library
    qsArrayD = {}
    for cutoff, readLen, i in zip(comh.SAMPLES_READ_TRIM_CUTOFFS, comh.ART_READ_LEN, range(len(comh.ART_READ_LEN))):
        qsArrayD[i] = sam.readCutoffArray(os.path.join(specDir, comh.SAMPLES_DIR, comh.SAMPLES_ERROR_QS_CUTOFF),
                                          cutoff, readLen)

    # list of (pair1.fq, pair2.fq, filtered.fq, qsArray, readLen, trimRemain)
    fileTupleList = []

    # directory containing all samples
    samplesDir = os.path.join(specDir, comh.SAMPLES_DIR)

    # collect all pairs of (fq1, fq2) files, define the output file
    for sample in os.listdir(samplesDir):
        if sample.isdigit():
            sampleDir = os.path.join(samplesDir, sample)
            if os.path.isdir(sampleDir):
                # for all strains of one sample
                for strain in os.listdir(sampleDir):
                    strainDir = os.path.join(sampleDir, strain)
                    if os.path.isdir(strainDir):
                        # get fq1 and fq2 files for all libraries of one strain
                        fq1D = {}
                        fq2D = {}
                        for f in os.listdir(strainDir):
                            i = f.split('_', 1)[0]
                            if i.isdigit():
                                if f.endswith('1.fq.gz'):
                                    fq1D[int(i)] = os.path.join(strainDir, f)
                                elif f.endswith('2.fq.gz'):
                                    fq2D[int(i)] = os.path.join(strainDir, f)
                        assert len(fq1D) == len(fq2D)
                        for i, fq1 in fq1D.iteritems():
                            fq2 = fq2D[i]
                            # store the files and param (pair1.fq, pair2.fq, filtered.fq, qsArray, readLen, trimRemain)
                            fileTupleList.append((fq1, fq2, os.path.join(os.path.dirname(fq1),
                                                                         '%s_filtered_qs.fq.gz' % i),
                                                                         qsArrayD[i], comh.ART_READ_LEN[i],
                                                                         comh.SAMPLES_READ_TRIM_REMAIN[i]))
    # filter reads according to the QS cutoffs
    fq.qsFilter(fileTupleList)


    # SAMPLES_READ_TRIM_CUTOFFS = [0.07, 0.07]
    # ART_READ_LEN = [100, 100]
    # SAMPLES_ERROR_QS_CUTOFF = 'samples_error_qs_cutoffs.csv'
    # implement in the SAM file

    # define it as a task.. create com. fq.py !!! implementing the reading from fq and filtering

    # output filtering report..


def _tmp():
    # test ..
    readLen = 100
    qsArrayLen = 64
    fq1 = '/Users/ivan/Documents/nobackup/hsim01/562/samples/4/NZ_AIEZ00000000/0_pair1.fq'
    fq2 = '/Users/ivan/Documents/nobackup/hsim01/562/samples/4/NZ_AIEZ00000000/0_pair2.fq'
    samFile = '/Users/ivan/Documents/nobackup/hsim01/562/samples/4/NZ_AIEZ00000000/0_pair.sam'
    # samFile = '/Users/ivan/Documents/nobackup/hsim01/562/samples/5/NZ_AIHP00000000/0_pair.sam'
    fasFile = '/Users/ivan/Documents/nobackup/hsim01/562/samples/4/NZ_AIEZ00000000/0_NZ_AIEZ00000000.fna'
    # _getStatSam(samFile)


def getGenesStat(specDir):
    srcDir = os.path.join(specDir, comh.FASTA_PULL_GENES_PHYLO_ALIGN_DIR_NAME)
    geneToAccSet = {}
    for f in os.listdir(srcDir):
        filePath = os.path.join(srcDir, f)
        for name in fas.fastaFileToDictWholeNames(filePath).keys():
            d = comh.getSeqNameToDict(name)
            acc = d['accession']
            geneName = d['geneName']
            if geneName not in geneToAccSet:
                geneToAccSet[geneName] = set()
            geneToAccSet[geneName].add(acc)

    counts = []
    for v in geneToAccSet.values():
        counts.append(len(v))

    total = max(counts)

    for i in range(100, 0, -1):
        count = 0
        threshold = total * float(i) / 100.
        for e in counts:
            if e >= threshold:
                count += 1

        print('%s%% -> %s' % (i, count))



    # count = 0
    # for f in os.listdir(srcDir):
    #     filePath = os.path.join(srcDir, f)
    #
    #     names = fas.fastaFileToDictWholeNames(filePath).keys()
    #     s = set(map(lambda x: getSeqNameToDict(x)['geneName'], names))
    #     if len(s) != 1:
    #         print f, len(names)
    #         count += 1

        # geneName = getSeqNameToDict(names[0])['geneName']
        # if len(names) > 290:
        #     geneCount[geneName] = len(names)
    # TODO !!!
    #
    #
    # print geneCount
    # print count

def fun():
    specList = ["562"]
    geneFileList = ['/Users/ivan/Documents/nobackup/hsim01/562/mg_31.txt',
                    '/Users/ivan/Documents/nobackup/hsim01/562/phylo_genes_mlst_7.txt',
                    '/Users/ivan/Documents/nobackup/hsim01/562/phylotyping_3.txt',
                    '/Users/ivan/Documents/nobackup/hsim01/562/virulent_3.txt',
                    '/Users/ivan/Documents/nobackup/hsim01/562/virulent_typing_20.txt']

    for f in geneFileList:
        geneList = csv.getColumnAsList(f)
        for species in specList:
            speciesDir = os.path.join(comh.REFERENCE_DIR_ROOT, species)
            pullGenesDir = os.path.join(speciesDir, comh.FASTA_PULL_GENES_DIR_NAME)
            accAllList = _getAccIntersectAll(pullGenesDir, geneList)  # list of accessions of genomes or draft genomes containing all
            print('spec: %s, accession count: %s, file: %s' % (species, len(accAllList), f))

    # Get genes based on which
    # Generates simulated dataset scenarios for one folder !!!
    # Escherichia coli, 562, species
    # assert os.path.isfile(NCBI_TAXONOMY_FILE), "The taxonomy file does not exists!"
    # assert os.path.isdir(REFERENCE_DIR_ROOT), "The root of the reference does not exist!"
    # taxonomy = taxonomy_ncbi.TaxonomyNcbi(NCBI_TAXONOMY_FILE, considerNoRank=True)
    # getAlignments('/Users/ivan/Documents/nobackup/hsim01/562/s', '/Users/ivan/Documents/nobackup/hsim01/562/d')
    # taxonomy.close()


if __name__ == "__main__":
    _main()
    # getAlignments('/Users/ivan/Documents/nobackup/hsim01/562/a', '/Users/ivan/Documents/nobackup/hsim01/562/b', maxCPU=1)