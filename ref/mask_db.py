#!/usr/bin/env python

import re

from com.csv import getMapping
from com.csv import OutFileBuffer
from com.csv import getColumnAsList
from com.fasta import fastaFileToDict
from com.fasta import filterOutSequences
from com.taxonomy_ncbid import TaxonomyNcbi


def _getLabelsCreateFasta():
    """
        To process the original mercier dataset with 59 strains. Take only contigs that were mapped to the reference
        genomes. Output a fasta file and a mapping file.
    :rtype : None
    """
    # input fasta file
    fastaFilePath = '/net/metagenomics/projects/PPSmg/data/V35/contigs_1000.txt' #contigs_1000.txt
    seqIdToSeq = fastaFileToDict(fastaFilePath)

    # contigs mapped to genome names
    nameLabelsFilePath = '/net/metagenomics/projects/PPSmg/data/V35/contigs_1000_blast_labels.txt' #contigs_1000_blast_labels.txt
    seqIdToNameLabels = getMapping(nameLabelsFilePath, 0, 1, sep='\t', comment = '#')

    # mapping: genome name -> taxon id
    genomeListFilePath = '/net/metagenomics/projects/PPSmg/data/V35/genome_list2.txt' #genome_list.txt
    nameLabelToNcbid = getMapping(genomeListFilePath, 0, 2, sep=';', comment = '#')

    # to store mapped sequences
    outFastaFilePath = '/net/metagenomics/projects/PPSmg/data/V35/contigsMappedBlast1000.fna' #contigsMappedBlast1000.fna
    outFasta = OutFileBuffer(outFastaFilePath)
    # to stored taxonomic mapping of mapped sequences
    outLabelsFilePath = '/net/metagenomics/projects/PPSmg/data/V35/contigsMappedBlast1000Labels.txt' #contigsMappedBlast1000Labels.txt
    outLabels = OutFileBuffer(outLabelsFilePath)

    for seqId in seqIdToSeq:
        if seqId in seqIdToNameLabels:
            outFasta.writeText(str('>' + str(seqId) + '\n' + seqIdToSeq[seqId] + '\n'))

    outFasta.close()
    print 'fasta created'

    for seqId in seqIdToSeq:
        if seqId in seqIdToNameLabels:
            nameLabel = seqIdToNameLabels[seqId][0]
            ncbid = nameLabelToNcbid[nameLabel][0]
            outLabels.writeText(str(str(seqId) + '\t' + str(ncbid) + '\n'))

    outLabels.close()
    print 'labels created'


def getFirstLabelAtAllowedRank():
    rank='species' # !!!!!!!

    predFile1 = '/net/metagenomics/projects/PPSmg/data/V35/contigsMappedBlast1000Labels.txt'
    predFile2 = '/net/metagenomics/projects/PPSmg/data/V35/contigsMappedBlast1000LabelsSpecies.txt'
    seqIdToLabel = getMapping(predFile1, 0, 1, sep='\t', comment = '#')
    outPred = OutFileBuffer(predFile2)

    taxonomy = TaxonomyNcbi('/net/metagenomics/projects/PPSmg/data/nobackup/NCBI20120828/ncbiTax/ncbitax_sqlite.db')

    for seqId in seqIdToLabel:
        ncbid = int(seqIdToLabel[seqId][0])
        while not taxonomy.isRankNcbidAllowed(ncbid):
            ncbid = taxonomy.getParentNcbid(ncbid)
        outPred.writeText(str(seqId + '\t' + str(ncbid) + '\n'))

    taxonomy.close()
    outPred.close()


def removeLines(mg):
    removeListFilePath = '/net/metagenomics/projects/PPSmg/data/V35/genome_ncbids.txt'
    #removeListFilePath = '/net/metagenomics/projects/PPSmg/data/V35/genome_accession_silva.txt'
    srcFilePath = str('/net/metagenomics/projects/PPSmg/data/markerGenes/db/' + mg + '_bact+arch_dnaV.tax')
    dstFilePath = str('/net/metagenomics/projects/PPSmg/data/V35/genomesRemoved/markerGenes/db/' + mg + '_bact+arch_dnaV.tax')
    #srcFilePath = str('/net/metagenomics/projects/PPSmg/data/silva/' + mg + '_silva106_ncbitax.bacteria+archaea.tax' )
    #dstFilePath = str('/net/metagenomics/projects/PPSmg/data/V35/genomesRemoved/silva/' + mg + '_silva106_ncbitax.bacteria+archaea.tax' )
    pattern = r'.*ncbid:([0-9]+)$'
    #pattern = r'^([^\-]+)\-.*$'

    removeSet = set(getColumnAsList(removeListFilePath, colNum=0, comment='#'))
    col0 = getColumnAsList(srcFilePath, colNum=0, sep='\t', comment='#')
    col1 = getColumnAsList(srcFilePath, colNum=1, sep='\t', comment='#')
    out = OutFileBuffer(dstFilePath)
    removed = 0
    for col0,col1 in zip(col0,col1):
        if re.sub(pattern, r'\1', col0) not in removeSet:
            out.writeText(str(col0 + '\t' + col1 + '\n'))
        else:
            removed += 1

    out.close()
    print mg, 'removeLines', removed


def removeEntries(mg):
    """
        Removes sequences from the marker gene files at the level from species, genus, family etc.
    """
    removeListPath = '/net/metagenomics/projects/PPSmg/data/V35/genome_ncbids_species.txt'
    srcFilePath = str('/net/metagenomics/projects/PPSmg/data/markerGenes2/db/' + mg + '_bact+arch_dnaV.tax')
    dstFilePath = str('/net/metagenomics/projects/PPSmg/data/V35/mgScenarios/speciesRemoved/db/' + mg + '_bact+arch_dnaV.tax')
    out = OutFileBuffer(dstFilePath)
    removeSet = set(getColumnAsList(removeListPath, colNum=0, comment='#'))
    removeSetInt = set([])
    removeSetIds = set([])
    removed = 0
    for s in removeSet:
        if s != '':
            removeSetInt.add(int(s))
    col0 = getColumnAsList(srcFilePath, colNum=0, sep='\t', comment='#')
    col1 = getColumnAsList(srcFilePath, colNum=1, sep='\t', comment='#')
    for col0,col1 in zip(col0,col1):
        lineSetInt = set([])
        for s in col1.split(';'):
            if s != '':
                lineSetInt.add(int(s))
        if len(removeSetInt.intersection(lineSetInt)) > 0: #the intersection is not empty
            removed += 1
            removeSetIds.add(col0)
        else:
            out.writeText(str(col0 + '\t' + col1 + '\n'))
    out.close()

    print mg, 'removedEntries', removed

    srcFilePath = str('/net/metagenomics/projects/PPSmg/data/markerGenes2/db/' + mg + '_bact+arch_dnaV.noalign.fna')
    dstFilePath = str('/net/metagenomics/projects/PPSmg/data/V35/mgScenarios/speciesRemoved/db/' + mg + '_bact+arch_dnaV.noalign.fna')
    out = OutFileBuffer(dstFilePath)
    seqIdToSeq = fastaFileToDict(srcFilePath)
    removed=0
    for seqId in seqIdToSeq:
        if seqId in removeSetIds:
            removed += 1
        else:
            out.writeText(str('>' + str(seqId) + '\n' + str(seqIdToSeq[seqId]) + '\n'))

    out.close()

    print mg, 'removedSeq', removed


def removeSequences(mg):
    removeListFilePath = '/net/metagenomics/projects/PPSmg/data/V35/genome_ncbids.txt'
    #removeListFilePath = '/net/metagenomics/projects/PPSmg/data/V35/genome_accession_silva.txt'
    srcFilePath = str('/net/metagenomics/projects/PPSmg/data/markerGenes/db/' + mg + '_bact+arch_dnaV.noalign.fna')
    dstFilePath = str('/net/metagenomics/projects/PPSmg/data/V35/genomesRemoved/markerGenes/db/' + mg + '_bact+arch_dnaV.noalign.fna')
    #srcFilePath = str('/net/metagenomics/projects/PPSmg/data/silva/' + mg + '_silva106_ncbitax.bacteria+archaea.fna' )
    #dstFilePath = str('/net/metagenomics/projects/PPSmg/data/V35/genomesRemoved/silva/' + mg + '_silva106_ncbitax.bacteria+archaea.fna' )
    pattern = r'.*ncbid:([0-9]+)$'
    #pattern = r'^([^\-]+)\-.*$'

    removeSet = set(getColumnAsList(removeListFilePath, colNum=0, comment='#'))
    seqIdToSeq = fastaFileToDict(srcFilePath)
    out = OutFileBuffer(dstFilePath)
    removed = 0
    for seqId in seqIdToSeq:
        if re.sub(pattern, r'\1', str(seqId)) not in removeSet:
            out.writeText(str('>' + str(seqId) + '\n' + str(seqIdToSeq[seqId]) + '\n'))
        else:
            removed += 1

    out.close()
    print mg, 'removeSequences', removed


def filterSequences():
    """
        To filter sequences with a specific label.
    """
    inFileName = '/net/metagenomics/projects/PPSmg/data/V35/contigsMappedBlast1000.fna'
    outFileName = '/net/metagenomics/projects/PPSmg/data/V35/nostocRemoved/contigsMappedBlast1000NostocRm.fna'
    mapFileName = '/net/metagenomics/projects/PPSmg/data/V35/contigsMappedBlast1000Labels.txt'
    labelRemove = 103690
    #seq id -> label
    labelToIdsDict = getMapping(mapFileName, 1, 0, sep='\t', comment = '#')
    allowedNamesSet = set([])
    for i in labelToIdsDict:
        if int(i) != int(labelRemove):
            for j in labelToIdsDict[i]:
                allowedNamesSet.add(j)

    filterOutSequences(inFileName, outFileName, allowedNamesSet)


def _test1():
    pass


def _test2():

    #getFirstLabelAtAllowedRank()
    list = ['dnaG', 'infC', 'pgk', 'rplA', 'rplC', 'rplE', 'rplK', 'rplM', 'rplP', 'rplT', 'rpoB', 'rpsC', 'rpsI', 'rpsK', 'rpsS', 'tsf', 'frr', 'nusA', 'pyrG', 'rplB',
               'rplD', 'rplF', 'rplL', 'rplN', 'rplS', 'rpmA', 'rpsB', 'rpsE', 'rpsJ', 'rpsM', 'smpB', '5S']
    #list = ['ssuref','lsuparc']
    for mg in list:
        removeLines(mg)
        removeSequences(mg)

if __name__ == "__main__":
    _test1()

