#!/usr/bin/env python

"""
    Copyright (C) 2014  Ivan Gregor

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
import sys
import time
import numpy as np
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from algbioi.com import taxonomy_ncbi
from algbioi.com import csv
from algbioi.com import fasta as fas
from algbioi.com import parallel
from algbioi.com import rand
from algbioi.com import sam
from algbioi.hsim import comh


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


        # getGenesStat(specDir)  # check, put to a new module, gzip everything !!!
        # TODO: Get mapping (check: gzip.open) !!!

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
            for readLen, insertSize, insertSD, minSeqLen in zip(comh.ART_READ_LEN, comh.ART_INSERT_SIZE,
                                                                comh.ART_INSERT_SD, comh.ART_MIN_SEQ_LEN):
                # filter out short sequences from the (draft) genome fasta files
                srcNoShortSeq = os.path.join(strainSimDir, str(count) + '_' + os.path.basename(srcFna))
                fas.cpSeqNoShortSeq(srcFna, srcNoShortSeq, minSeqLen)

                # define the read generation command
                cmd = "%s -i %s -f %s -rs %s -l %s -m %s -s %s -p -o %s_pair -na -sam; " \
                      "gzip %s_pair1.fq; gzip %s_pair2.fq; gzip %s_pair.sam; gzip %s" \
                      % (comh.ART_ILLUMINA_BINARY, srcNoShortSeq, coverage,
                         rand.strToRandInt(str(comh.SAMPLES_DEF_RAND_SEED) + acc + str(round(coverage))
                                           + str(os.path.getsize(srcNoShortSeq))), readLen, insertSize, insertSD,
                         count, count, count, count, srcNoShortSeq)

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


# TODO: get the mapping !!!

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