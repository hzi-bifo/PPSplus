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


    ***********************************************************************

    Manage Pfam-A, get mapping between bacterial genes and Pfam gene families.
"""
import os
import sys

from algbioi.com import parallel
from algbioi.com import fasta as fas
from algbioi.com import csv

from algbioi.hsim import comh


def _main():

    for spec in comh.SPECIES_LIST:
        specDir = os.path.join(comh.REFERENCE_DIR_ROOT, spec)

        # translate gene DNA to gene PROT sequences
        if False:
            dnaToProt(specDir)

        # run HMM search to match gene PROT sequences to PFAM HMMs
        if False:
            searchForGeneFam(specDir)

        # get mapping between geneName and domainName
        if True:
            getMapGeneNameAndPfamDomainName(specDir)


def dnaToProt(specDir):
    """
        Translate gene DNA sequences to PROT gene sequences, create one file per geneName.
        Uses translation table 11.
    """
    print('Translating gene DNA sequences to PROT sequences.')

    srcGeneDir = os.path.join(specDir, comh.FASTA_PULL_GENES_PHYLO_DIR_NAME)
    dstProtGeneDir = os.path.join(specDir, comh.FASTA_PULL_GENES_PHYLO_PROT_DIR_NAME)

    # create dst directory
    if not os.path.isdir(dstProtGeneDir):
        os.mkdir(dstProtGeneDir)

    # translate sequences
    for f in os.listdir(srcGeneDir):
        fas.dnaToProt(os.path.join(srcGeneDir, f), os.path.join(dstProtGeneDir, f), translTable=11)


def searchForGeneFam(specDir):
    """
        Running hmmsearch to match genes and Pfam domains.
        Very time consuming!
    """
    print('Searching for Pfam gene families')
    srcProtGeneDir = os.path.join(specDir, comh.FASTA_PULL_GENES_PHYLO_PROT_DIR_NAME)
    taskList = []
    # for each gene DNA fasta file
    for f in os.listdir(srcProtGeneDir):
        filePath = os.path.join(srcProtGeneDir, f)
        basePath = filePath
        if filePath.endswith('.fna'):
            basePath = filePath[:-4]
        # define hmmsearch command
        cmd = '%s/hmmsearch -o /dev/null --noali --domtblout %s -E %s --cpu %s %s %s' % \
              (comh.HMMER_BINARY,
               # str(pathBase + '.out'),
               # str(pathBase + '.aln'),  # -A %s
               # str(pathBase + '.tblout'),  # --tblout %s
               str(basePath + '.domtblout'),
               0.01,
               1,
               os.path.join(comh.PFAM, 'Pfam-A.hmm'),
               filePath)
        taskList.append(parallel.TaskCmd(cmd, cwd=srcProtGeneDir))

    # run commands in parallel
    parallel.reportFailedCmd(parallel.runCmdParallel(taskList, comh.MAX_PROC))


def getMapGeneNameAndPfamDomainName(specDir):
    """
        Map Pfam-domain names onto geneNames.
    """
    print("Mapping Pfam-domain names onto gene names.")
    srcDir = os.path.join(specDir, comh.FASTA_PULL_GENES_PHYLO_PROT_DIR_NAME)
    assert os.path.isdir(srcDir)

    # output mapping files
    geneNameToDomOut = csv.OutFileBuffer(os.path.join(specDir, comh.FASTA_PULL_GENES_PHYLO_GENE_NAME_TO_PFAM))
    domToGeneNameOut = csv.OutFileBuffer(os.path.join(specDir, comh.FASTA_PULL_GENES_PHYLO_PFAM_TO_GENE_NAME))

    geneNameToDomEntryList = {}

    # for each dom file
    for f in os.listdir(srcDir):
        if f.endswith('.domtblout'):
            seqNameToList = {}
            geneName = None

            # for each line in the dom file
            for line in open(os.path.join(srcDir, f)):
                # skip comments
                if line.startswith('#'):
                    continue

                # parse a line
                tokens = line.split()
                assert len(tokens) == 23
                targetName, accession1, tlen, queryName, accession2, qlen, eValue, score1, bias, r1, r2, cEvalue, \
                iEvalue, score2, bias, from1, to1, from2, to2, from3, to3, acc, descriptionOfTarget = tokens

                # get the gene name
                if geneName is None:
                    headerD = dict(zip(map(lambda x: x.split(':')[0], targetName.split(';')),
                                       map(lambda x: x.split(':')[1], targetName.split(';'))))
                    geneName = headerD['geneName']

                # store the entry: gene name ->  (score, dom name)
                entry = (float(score1), queryName)
                if targetName not in seqNameToList:  # first record for the gene sequence name
                    seqNameToList[targetName] = [entry]
                else:
                    seqNameToList[targetName].append(entry)  # there are more records for the gene sequence name

            # for each sequence name, find the match with the highest score
            for targetName in seqNameToList:
                seqNameToList[targetName].sort(key=lambda x: x[0], reverse=True)

            # get a list of dom names
            domNameList = []
            for entry in seqNameToList.values():
                domNameList.append(entry[0])

            # store map: gene-name -> list of (score, dom-name)
            if len(domNameList) > 0:
                geneNameToDomEntryList[geneName] = domNameList

    # map: dom-name -> list of (gene-name, score)
    domNameToGeneNameList = {}

    # distribution: geneName -> count of dom names
    geneNameDist = {}

    # store map: gene-name -> dom-name list
    for geneName, domEntryList in geneNameToDomEntryList.iteritems():
        entrySet = set(map(lambda x: x[1], domEntryList))

        # get the support for each dom-name
        domToSupport = {}
        scoreSum = 0.
        for e in domEntryList:
            domName = e[1]
            score = e[0]
            scoreSum += score
            if domName in domToSupport:
                domToSupport[domName] += score
            else:
                domToSupport[domName] = score
        for domName in domToSupport.keys():
            domToSupport[domName] = (domToSupport[domName] / scoreSum) * 100.
        domNameTupleList = []
        for dom, support in domToSupport.iteritems():
            domNameTupleList.append((dom, support))
        domNameTupleList.sort(key=lambda x: x[1], reverse=True)
        outList = []
        for k, v in domNameTupleList:
            outList.append(k)
            outList.append(str(round(v, 1)))
        # store map: gene-name -> list of (dom-name, support)
        geneNameToDomOut.writeText('%s\t%s\n' % (geneName, ','.join(outList)))
        # update distribution
        dLen = len(entrySet)
        if dLen not in geneNameDist:
            geneNameDist[dLen] = 1
        else:
            geneNameDist[dLen] += 1
        # update mapping: dom name -> gene name
        for scoreDom in domEntryList:
            if scoreDom[1] not in domNameToGeneNameList:
                domNameToGeneNameList[scoreDom[1]] = []
            domNameToGeneNameList[scoreDom[1]].append((geneName, scoreDom[0]))

    # map: dom-name -> list of (gene-name, score)
    domNameToGeneScoreList = {}
    for domName, entryList in domNameToGeneNameList.iteritems():
        geneToScore = {}
        scoreSum = 0.
        for entry in entryList:
            geneName = entry[0]
            score = entry[1]
            scoreSum += score
            if geneName in geneToScore:
                geneToScore[geneName] += score
            else:
                geneToScore[geneName] = score
        geneScoreTupleList = []
        for geneName in geneToScore.keys():
            geneScoreTupleList.append((geneName, (geneToScore[geneName] / scoreSum) * 100.))

        geneScoreTupleList.sort(key=lambda x: x[1], reverse=True)

        domNameToGeneScoreList[domName] = []
        for geneScore in geneScoreTupleList:
            domNameToGeneScoreList[domName].append(geneScore[0])
            domNameToGeneScoreList[domName].append(str(geneScore[1]))

    # distribution: dom name -> count of gene names
    domNameDist = {}

    # store map: dom name -> gene name
    for domName, geneNameScoreList in domNameToGeneScoreList.iteritems():
        domToGeneNameOut.writeText('%s\t%s\n' % (domName, ','.join(geneNameScoreList)))
        # update distribution
        gLen = len(geneNameScoreList) / 2
        if gLen not in domNameDist:
            domNameDist[gLen] = 1
        else:
            domNameDist[gLen] += 1

    geneNameToDomOut.close()
    domToGeneNameOut.close()
    print('geneName -> domName distr: %s' % str(geneNameDist))
    print('domName -> geneName distr: %s' % str(domNameDist))


if __name__ == "__main__":
    _main()


    # print 'targetName', targetName
    # print 'accession1', accession1
    # print 'tlen', tlen
    # print 'queryName', queryName
    # print 'accession2', accession2
    # print 'qlen', qlen
    # print 'eValue', eValue
    # print 'score1', score1
    # print 'bias', bias
    # print 'cEvalue', cEvalue
    # print 'iEvalue', iEvalue
    # print 'score2', score2
    # print 'bias', bias
    # print 'from1', from1
    # print 'to1', to1
    # print 'from2', from2
    # print 'to2', to2
    # print 'from3', from3
    # print 'to3', to3
    # print 'acc', acc
    # print 'descriptionOfTarget', descriptionOfTarget
    # print 'rest', r1, r2