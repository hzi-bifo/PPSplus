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
    Manage the strains - reference data for a list of given species.

    NCBI genomes and draft genomes downloaded using commands (for e.coli):
    wget -r --no-parent ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Escherichia_coli*
    wget -r --no-parent ftp://ftp.ncbi.nih.gov/genomes/Bacteria_DRAFT/Escherichia_coli*

    Other links:
    http://www.ncbi.nlm.nih.gov/books/NBK44863/
    http://www.ncbi.nlm.nih.gov/genome/?term=escherichia+coli

    Genetic code.
    http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi

    # Batch entrez
    http://www.ncbi.nlm.nih.gov/sites/batchentrez
"""
import os
# import sys
# import re
from Bio import SeqIO
# from Bio.Seq import Seq
# from Bio.Alphabet import generic_dna
# from Bio.SeqFeature import SeqFeature
# from Bio.SeqFeature import FeatureLocation
from algbioi.com import taxonomy_ncbi
from algbioi.com import csv
from algbioi.com import fasta
from algbioi.com import gbk
from algbioi.hsim import comh
from algbioi.com import parallel


# common directory and file names
GENOMES_DIR = 'bacteria_genomes'  # directory containing downloaded (decompressed) genomes gbk files
DRAFT_GENOMES_DIR = 'bacteria_draft'  # directory containing downloaded (decompressed) draft genomes gbk files

GENOMES_SOURCES_FILE_NAME = 'bacteria_genomes_sources.csv'  # contains a list of source files containing genomes


def _main():
    """
        Main function, choose steps you want to execute (change: False -> True)
    """
    assert os.path.isfile(comh.NCBI_TAXONOMY_FILE)
    assert os.path.isdir(comh.REFERENCE_DIR_ROOT)
    taxonomy = taxonomy_ncbi.TaxonomyNcbi(comh.NCBI_TAXONOMY_FILE, considerNoRank=True)

    for spec in comh.SPECIES_LIST:
        specDir = os.path.join(comh.REFERENCE_DIR_ROOT, spec)  # directory containing species data
        assert comh.isSpeciesDirectory(specDir, taxonomy), "Not a species directory: %s" % specDir

        # Get a list of relevant files for each GENOME strain
        if False:
            getBacteriaGenomesFileList(specDir, taxonomy)

        # Extract DNA sequences: for each strain one file
        if False:
            storeGenomesAsFastaDNA(specDir, taxonomy)

        # Extract DNA genes for genomes: for each gene one file
        if False:
            storeGenesAsFastaDNA(specDir, taxonomy)

        # Decompress draft genome (gbk.tgz) files
        if False:
            decompressDraftGenomeFiles(specDir, suffixList=['gbk.tgz'])

        # Extract draft genome DNA sequences: for each draft genome one file
        if False:
            storeDraftGenomesAsFastaDNA(specDir, taxonomy)

        # Extract DNA genes for draft genomes: for each gene one file
        if False:
            storeDraftGenomeGenesAsFastaDNA(specDir, taxonomy)

        # Pull genes from genomes and draft genomes together.
        if False:
            pullGenesFastaTogether(specDir)

    taxonomy.close()


def getBacteriaGenomesFileList(specDir, taxonomy):
    """
        Get a list of relevant files for each GENOME strain (importantly the gbk files).

        Create file: GENOMES_SOURCES_FILE_NAME
    """
    taxonId = int(os.path.basename(specDir))
    speciesName = taxonomy.getScientificName(taxonId)
    genomesDir = os.path.join(specDir, GENOMES_DIR)
    genomeCount = len(os.listdir(genomesDir))
    print("Getting relevant gbk files for genomes")
    print("Species id: %s name: %s" % (taxonId, speciesName))
    print("Genome count: %s" % genomeCount)

    # for each genome in a directory, get a gbk containing a "complete genome"
    gbkCompleteGenomeList = []
    count = 0.
    for genomeDir in os.listdir(genomesDir):
        genomeDirPath = os.path.join(genomesDir, genomeDir)
        # print('Processing: %s' % genomeDir)
        count += 1.
        print 'Progress: ' + str(round((float(count)/genomeCount) * 100, 1)) + '%\r',

        # for all files in the genome directory, search for a gbk file containing a "complete genome"
        gbkCompleteGenome = None
        for f in os.listdir(genomeDirPath):
            if f.split('.')[-1] == 'gbk':  # got a gbk file
                records = []  # get all records
                fPath = os.path.join(genomeDirPath, f)
                for record in SeqIO.parse(fPath, "genbank"):
                    records.append(record)
                if len(records) == 0:
                    print('No record contained in file "%s' % fPath)
                    continue
                # get a record with a complete genome
                record = None
                for r in records:
                    if 'complete genome' in str(r.description).lower():
                        if record is None:
                            record = r
                        else:
                            print('File "%s" contains more than one record containing a complete genome!' % fPath)
                # get a file containing a record with a complete genome
                if record is not None:
                    if gbkCompleteGenome is None:
                        gbkCompleteGenome = fPath
                    else:
                        print('There are at least two gbk files describing a "complete genome" for one strain: '
                              '"%s" and "%s" (the first one was taken)' % (gbkCompleteGenome, fPath))

        # there is a gbk file containing a complete genome, get all corresponding relevant files
        if gbkCompleteGenome is not None:
            d = {}
            for suffix in ['ptt', 'rpt', 'fna', 'gbk']:
                f = ".".join(gbkCompleteGenome.split('.')[0:-1]) + '.' + suffix
                if not os.path.isfile(f):
                    print('File "%s" does not exists!' % f)
                else:
                    d[suffix] = os.path.sep.join(f.split(os.path.sep)[-4:])
            gbkCompleteGenomeList.append(d)
        else:
            print('No file containing a complete genome found in "%s"' % genomeDirPath)

    # write paths to the relevant (importantly gbk) files to a source file (relative to the root)
    out = csv.OutFileBuffer(os.path.join(specDir, GENOMES_SOURCES_FILE_NAME))
    for d in gbkCompleteGenomeList:
        out.writeText(d.get('gbk') + ',' + d.get('ptt') + ',' + d.get('fna') + ',' + d.get('rpt') + '\n')
    out.close()


def storeGenomesAsFastaDNA(specDir, taxonomy):
    """
        Get a FASTA file for each genome.

        Create directory: comh.FASTA_GENOMES_DIR_NAME
    """
    print('Storing genomes as FASTA.')
    # make the destination directory (to store genomes as FASTA) if it doesn't exists
    dstDir = os.path.join(specDir, comh.FASTA_GENOMES_DIR_NAME)
    if not os.path.isdir(dstDir):
        os.mkdir(dstDir)

    # get all gbk files containing genomes
    gbkFilePathList = csv.getColumnAsList(os.path.join(specDir, GENOMES_SOURCES_FILE_NAME), colNum=0, sep=',')
    # define tasks
    taskList = []
    for gbkFilePath in gbkFilePathList:
        taskList.append(parallel.TaskThread(gbk.readFromGbkFile, (os.path.join(comh.REFERENCE_DIR_ROOT, gbkFilePath),
                                                                  taxonomy, ['accession', 'acc_version', 'seq'])))
    # extract genome sequences from gbk files
    itemsList = parallel.runThreadParallel(taskList, maxThreads=comh.MAX_PROC)

    # write genome sequences to files (one genome one FASTA file)
    for items in itemsList:
        assert len(items) == 1
        accession = items[0]['accession']
        accVersion = items[0]['acc_version']
        seq = items[0]['seq']
        out = csv.OutFileBuffer(os.path.join(dstDir, accession + '.fna'))
        out.writeText('>' + accVersion + '\n' + seq)
        out.close()


def storeGenesAsFastaDNA(specDir, taxonomy):
    """
        Get a FASTA file for each gene for genomes. (Do not store multiple copies !!!)

        Create directory: comh.FASTA_GENOMES_GENES_DIR_NAME
    """
    print("Storing genes as fasta for genomes")
    # make the destination directory (to store a FASTA file for each gene) if it doesn't exists
    dstDir = os.path.join(specDir, comh.FASTA_GENOMES_GENES_DIR_NAME)
    if not os.path.isdir(dstDir):
        os.mkdir(dstDir)

    # get all gbk files containing genomes
    gbkFilePathList = csv.getColumnAsList(os.path.join(specDir, GENOMES_SOURCES_FILE_NAME), colNum=0, sep=',')
    # define tasks
    taskList = []
    for gbkFilePath in gbkFilePathList:
        taskList.append(parallel.TaskThread(gbk.readFromGbkFile, (os.path.join(comh.REFERENCE_DIR_ROOT, gbkFilePath),
                                                                  taxonomy, ['genes_annotation', 'accession'])))
    # extract gene sequences from gbk files
    itemsList = parallel.runThreadParallel(taskList, maxThreads=comh.MAX_PROC)

    # for each genome, write gene sequences to files (one gene one FASTA file)
    for items in itemsList:
        assert len(items) == 1
        genesAnnotation = items[0]['genes_annotation']
        accession = items[0]['accession']

        # for each gene, write gene sequence to a file
        for annot in genesAnnotation.values():
            filePath = os.path.join(dstDir, comh.getGeneNameToFileName(annot.geneName))
            if os.path.isfile(filePath):
                out = csv.OutFileBuffer(filePath, fileOpenMode='a')
                out.writeText('\n')
            else:
                out = csv.OutFileBuffer(filePath)

            out.writeText('>accVersion:%s;geneName:%s;locusTag:%s;src:genome;accession:%s\n%s'
                          % (annot.sequenceAccessionVersion, annot.geneName, annot.locusTag, accession, annot.seqDNA))
            out.close()


def decompressDraftGenomeFiles(specDir, suffixList=['gbk.tgz']):
    """
        Decompresses all draft genome gbk files with the given suffix,
        does nothing in the case the decompressed directory already exists.
    """
    print("Decompressing draft genome gbk.tgz files")
    draftDirRoot = os.path.join(specDir, DRAFT_GENOMES_DIR)

    # for each draft genome
    for draft in os.listdir(draftDirRoot):
        draftGenomeDir = os.path.join(draftDirRoot, draft)
        assert os.path.isdir(draftGenomeDir)
        taskList = []
        # for each file in the draft genome directory
        for f in os.listdir(draftGenomeDir):
            # take files with the corresponding suffixes
            for suffix in suffixList:
                if f.endswith(suffix):
                    # add task to a list (file will be decompressed)
                    taskList.append(parallel.TaskThread(comh.extract, (os.path.join(draftGenomeDir, f),)))
                    break

        # decompress all files (a file won't be decompressed in the case it has already been decompressed)
        parallel.runThreadParallel(taskList, maxThreads=comh.MAX_PROC)


def storeDraftGenomesAsFastaDNA(specDir, taxonomy):
    """
        Get a FASTA file for each draft genome.

        Create directory: comh.FASTA_GENOMES_DRAFT_DIR_NAM
    """
    print('Storing draft genomes as FASTA.')
    # make the destination directory (to store draft genomes as FASTA) if it doesn't exists
    dstDir = os.path.join(specDir, comh.FASTA_GENOMES_DRAFT_DIR_NAME)
    if not os.path.isdir(dstDir):
        os.mkdir(dstDir)

    draftGenomesRoot = os.path.join(specDir, DRAFT_GENOMES_DIR)
    # for each draft genome
    for draftGenome in os.listdir(draftGenomesRoot):
        draftGenomeDir = os.path.join(draftGenomesRoot, draftGenome)
        # find a gbk file and gbk dir
        gbkFile = None
        gbkDir = None
        for f in os.listdir(draftGenomeDir):
            filePath = os.path.join(draftGenomeDir, f)
            if filePath.endswith('.gbk'):
                if os.path.isfile(filePath):
                    assert gbkFile is None, filePath
                    gbkFile = filePath
                elif os.path.isdir(filePath):
                    assert gbkDir is None, filePath
                    gbkDir = filePath
                else:
                    assert False, filePath

        if gbkFile is None or gbkDir is None:
            print("Entry skipped, draft genome reference not complete! gbkFile: %s gbkDir: %s" % (gbkFile, gbkDir))
            continue

        # get the accession of this draft genome
        item = gbk.readFromGbkFile(gbkFile, taxonomy, ['accession'])
        assert len(item) == 1, 'More records: %s' % item
        accession = item[0]['accession']

        # open a file to store DNA sequences (scaffolds) of this draft genome
        out = csv.OutFileBuffer(os.path.join(dstDir, accession + '.fna'))

        # for each scaffold of the draft genome, define task - read its sequence and accession version
        taskList = []
        for f in os.listdir(gbkDir):
            gbkFilePath = os.path.join(gbkDir, f)
            if gbkFilePath.endswith('.gbk'):
                taskList.append(parallel.TaskThread(gbk.readFromGbkFile,
                                                    (gbkFilePath, taxonomy, ['acc_version', 'seq'])))

        # read all scaffold sequences
        itemList = parallel.runThreadParallel(taskList, maxThreads=comh.MAX_PROC)

        # write all scaffolds to a file
        addNewLine = False
        for item, task in zip(itemList, taskList):
            assert len(item) == 1, 'More than one item: items: %s task args: %s' % (len(item), task.args)
            accVersion = item[0]['acc_version']
            seq = item[0]['seq']
            if addNewLine:
                out.writeText("\n")
            else:
                addNewLine = True
            out.writeText(">%s\n%s" % (accVersion, seq))
        out.close()


def storeDraftGenomeGenesAsFastaDNA(specDir, taxonomy):
    """
        Get a FASTA file for each gene for draft genomes. (Do not store multiple copies !!!)

        Create directory: comh.FASTA_GENOMES_DRAFT_GENES_DIR_NAME
    """
    print("Storing genes as fasta for draft genomes")
    # make the destination directory (to store a FASTA file for each gene) if it doesn't exists
    dstDir = os.path.join(specDir, comh.FASTA_GENOMES_DRAFT_GENES_DIR_NAME)
    if not os.path.isdir(dstDir):
        os.mkdir(dstDir)

    # for each draft genome
    draftGenomesRoot = os.path.join(specDir, DRAFT_GENOMES_DIR)
    for draftGenome in os.listdir(draftGenomesRoot):
        draftGenomeDir = os.path.join(draftGenomesRoot, draftGenome)
        # find a gbk file and gbk dir
        gbkFile = None
        gbkDir = None
        for f in os.listdir(draftGenomeDir):
            filePath = os.path.join(draftGenomeDir, f)
            if f.endswith('.gbk'):
                if os.path.isfile(filePath):
                    assert gbkFile is None, filePath
                    gbkFile = filePath
                elif os.path.isdir(filePath):
                    assert gbkDir is None, filePath
                    gbkDir = filePath
                else:
                    assert False, filePath

        if gbkFile is None or gbkDir is None:
            print("Entry skipped, draft genome reference not complete! gbkFile: %s gbkDir: %s" % (gbkFile, gbkDir))
            continue

        # get the accession of this draft genome
        item = gbk.readFromGbkFile(gbkFile, taxonomy, ['accession'])
        assert len(item) == 1, 'More records: %s' % item
        accession = item[0]['accession']

        # for each scaffold of the draft genome, define task - read all its gene annotations
        taskList = []
        for f in os.listdir(gbkDir):
            gbkFilePath = os.path.join(gbkDir, f)
            if gbkFilePath.endswith('.gbk'):
                taskList.append(parallel.TaskThread(gbk.readFromGbkFile, (gbkFilePath, taxonomy, ['genes_annotation'])))

        # read all gene annotations
        itemList = parallel.runThreadParallel(taskList, maxThreads=comh.MAX_PROC)

        # store mapping: gene name -> list of gene records
        geneToRecordList = {}
        for item, task in zip(itemList, taskList):
            assert len(item) == 1, 'More than one item: items: %s task args: %s' % (item, task.args)
            for geneName, geneRecord in item[0]['genes_annotation'].iteritems():
                if geneName in geneToRecordList:
                    geneToRecordList[geneName].append(geneRecord)
                else:
                    geneToRecordList[geneName] = [geneRecord]

        # check whether there is only one record for one gene, store single copy genes to files (one gene, one file)
        for v in geneToRecordList.values():
            if len(v) == 1:  # the gene is single copy
                annot = v[0]
                filePath = os.path.join(dstDir, comh.getGeneNameToFileName(annot.geneName))
                if os.path.isfile(filePath):
                    out = csv.OutFileBuffer(filePath, fileOpenMode='a')
                    out.writeText('\n')
                else:
                    out = csv.OutFileBuffer(filePath)
                out.writeText('>accVersion:%s;geneName:%s;locusTag:%s;src:draft;accession:%s\n%s'
                              % (annot.sequenceAccessionVersion, annot.geneName, annot.locusTag, accession,
                                 annot.seqDNA))
                out.close()
            else:
                assert len(v) > 0, draftGenomeDir


def pullGenesFastaTogether(specDir):
    """
        Pull genes from genomes and draft genomes together.

        Create directory: comh.FASTA_PULL_GENES_DIR_NAME
    """
    print("Pulling genes from genomes and draft genomes together.")
    # make the destination directory if it doesn't exists
    dstDir = os.path.join(specDir, comh.FASTA_PULL_GENES_DIR_NAME)
    if not os.path.isdir(dstDir):
        os.mkdir(dstDir)

    # get a set of all gene names contained in genomes or draft genomes (as the corresponding file names)
    allFileNameSet = set()
    srcDirList = [comh.FASTA_GENOMES_GENES_DIR_NAME, comh.FASTA_GENOMES_DRAFT_GENES_DIR_NAME]
    for srcDir in srcDirList:
        for f in os.listdir(os.path.join(specDir, srcDir)):
            if f.endswith('.fna'):
                allFileNameSet.add(f)
            else:
                print('file %s does not end with ".fna"' % f)

    # for each gene create a file and copy all sequences for this gene to this file
    for f in allFileNameSet:
        dstFile = csv.OutFileBuffer(os.path.join(dstDir, f))
        firstLine = True
        for srcDir in srcDirList:
            filePath = os.path.join(specDir, srcDir, f)
            if os.path.isfile(filePath):
                for name, seq in fasta.fastaFileToDictWholeNames(filePath).iteritems():
                    if firstLine:
                        dstFile.writeText(">%s\n%s" % (name, seq))
                        firstLine = False
                    else:
                        dstFile.writeText("\n>%s\n%s" % (name, seq))
        dstFile.close()


if __name__ == "__main__":
    _main()