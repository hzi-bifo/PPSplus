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


    Contains functionality to read from a gbk file.
"""

import re
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqFeature import FeatureLocation
# from algbioi.com import taxonomy_ncbi
# import os
# from Bio import Entrez
# from Bio.SeqFeature import SeqFeature
# from algbioi.com import csv


KEYWORDS = {'description', 'accession', 'gi', 'acc_version', 'taxonId', 'seq', 'seq_rev', 'gene_name_list',
            'genes_annotation'}


class GeneRecord():
    """
        Represents a gene within a sequence.
    """
    def __init__(self, geneName, locationList, seqDNA, seqPROT, translTable, locusTag, sequenceAccessionVersion):
        self.geneName = geneName
        self.locationList = locationList
        self.seqDNA = seqDNA
        self.seqPROT = seqPROT
        self.translTable = translTable
        self.locusTag = locusTag
        self.sequenceAccessionVersion = sequenceAccessionVersion


def getLocationList(locationAnnotation):
    """
        Get a list of feature locations from a join location annotation.

        Example input: 'join{[872761:872836](+), [872837:873860](+)}'

        @type locationAnnotation: str
        @rtype: list
    """
    assert locationAnnotation.startswith('join')
    retList = []
    l = locationAnnotation.replace('join{', '').replace('}', '').split(',')
    for entry in l:
        entry = entry.strip(' ')
        start = int(re.sub(r'^\[([0-9]+):[0-9]+\]\([\+\-]\)$', r'\1', entry))
        stop = int(re.sub(r'^\[[0-9]+:([0-9]+)\]\([\+\-]\)$', r'\1', entry))
        strand = re.sub(r'^\[[0-9]+:[0-9]+\]\(([\+\-])\)$', r'\1', entry)
        if strand == '+':
            strand = 1
        else:
            strand = -1
        retList.append(FeatureLocation(start, stop, strand=strand))
    return retList


def getRegions(seqStr, locationList):
    """
        Get a concatenated region from individual location annotations

        @param seqStr: DNA sequence
        @type seqStr: str
        @param locationList: list of locations to be extracted from the sequence
        @type locationList: list of FeatureLocation
        @return: Concatenated substrings corresponding to the location
        @rtype: str
    """
    if len(locationList) == 1:
        location = locationList[0]
        seq = seqStr[location.start:location.end]
        if location.strand < 0:
            seq = str(Seq(seq, generic_dna).reverse_complement())
        return seq
    else:
        return str(getRegions(seqStr, [locationList[0]]) + getRegions(seqStr, locationList[1:]))


def isLocationContinuous(locationList):
    """
        Is the list of locations within a genome continuous?

        @type locationList: list of FeatureLocation
        @return: True if the location described by the list is continuous
    """
    lastEnd = None
    for location in locationList:
        if (lastEnd is not None) and (lastEnd + 1 != location.start):
            return False
        lastEnd = location.end
    return True


def readFromGbkFile(gbkFile, taxonomy, listOfKeywords=list(KEYWORDS)):
    """
        For a list of input keywords, returns a list dictionaries (keyword -> value).
        One dictionary for one record (often, there is just one record in one gbk file)

        Don't consider genes of:
        -multiple copies
        -with no continuous location
        -if the translated region doesn't match the annotated protein

        Keywords:
        -description: record.description (DEFINITION)
        -accession: record.annotations["accessions"] (ACCESSION) take the first one if more available
        -gi: record.annotations["gi"] (GI)
        -acc_version: record.id (VERSION), accession version
        -taxonId: feature source, feature.qualifiers['db_xref']
        -seq: DNA sequence, str(record.seq)
        -seq_rev: reverse complement DNA sequence, str(record.reverse_complement().seq)
        -gene_name_list: list of gene names, feature gene, feature.qualifiers['gene'], (not including unnamed!)
        -genes_annotation: dict of genes, feature CDS (not including unnamed) (geneName -> GeneRecord);
        doesn't contain genes existing in multiple copies! (read this documentation!)

        Documentation biopython:
        http://biopython.org/DIST/docs/api/Bio.SeqFeature.SeqFeature-class.html

        @raise: Throws an exception when the taxon ID does not exist in the taxonomy.

        @param gbkFile: gbk annotation input file
        @type taxonomy: taxonomy_ncbi.TaxonomyNcbi
        @param listOfKeywords: list of keywords, values of which will be retrieved (described above)
        @return: list of dictionaries, one dictionary for one record, a dictionary contains values for input keywords
        @rtype: list of dict
    """
    retList = []
    # check all keywords
    keywordSet = set(listOfKeywords)
    for keyword in keywordSet:
        if keyword not in KEYWORDS:
            sys.stderr('Keyword %s not defined!\n' % keyword)

    # if gene annotation is required, gene accession will be retrieved as well
    if 'genes_annotation' in keywordSet:
        keywordSet.add('acc_version')

    # for each record
    for record in SeqIO.parse(gbkFile, "genbank"):
        recDict = {}

        if 'description' in keywordSet:
            recDict['description'] = record.description

        if 'accession' in keywordSet:
            accessionList = record.annotations['accessions']
            recDict['accession'] = accessionList[0]
            if len(accessionList) > 1:
                print('There are more than one accession: "%s", the first was taken!' % accessionList)

        if 'gi' in keywordSet:
            recDict['gi'] = record.annotations['gi']

        if 'acc_version' in keywordSet:
            recDict['acc_version'] = record.id

        if 'taxonId' in keywordSet:
            for feature in record.features:
                if feature.type == 'source':
                    dbxref = feature.qualifiers['db_xref']
                    assert len(dbxref) == 1, 'TaxonId: this list should be of length 1: "%s"' % dbxref
                    taxon, taxonId = dbxref[0].split(':')
                    assert taxon == 'taxon', 'This string should be "taxon", not "%s"' % taxon
                    if taxonomy.exists(int(taxonId)):
                        recDict['taxonId'] = int(taxonId)
                    else:
                        print("TaxonId does not exists in the taxonomy: %s" % taxonId)
                    break
            if 'taxonId' not in recDict:
                print('TaxonId not found for "%s" in "%s"' % (record.id, gbkFile))
                recDict['taxonId'] = None

        if 'seq' in keywordSet:
            recDict['seq'] = str(record.seq)

        if 'seq_rev' in keywordSet:
            recDict['seq_rev'] = str(record.reverse_complement().seq)

        if 'gene_name_list' in keywordSet:
            recDict['gene_name_list'] = []
            for feature in record.features:
                if feature.type == 'gene' and 'gene' in feature.qualifiers:  # !!!
                    recDict['gene_name_list'].append(feature.qualifiers['gene'].strip())

        if 'genes_annotation' in keywordSet:
            geneCount = 0
            genesDict = {}
            locSet = set()  # set of gene locations
            multipleCopyGeneSet = set()
            for feature in record.features:
                if feature.type == 'gene' and 'gene' in feature.qualifiers:
                    # print feature.qualifiers['locus_tag']
                    if 'locus_tag' not in feature.qualifiers:
                        print('No locus tag found for %s in %s' % (feature, gbkFile))
                        continue
                    assert len(feature.qualifiers['locus_tag']) == 1, '%s' % feature.qualifiers['locus_tag']
                    locSet.add(feature.qualifiers['locus_tag'][0])

            for feature in record.features:
                if feature.type == 'CDS' and 'locus_tag' in feature.qualifiers:
                    # print feature.qualifiers['locus_tag']
                    assert len(feature.qualifiers['locus_tag']) == 1, '%s' % feature.qualifiers['locus_tag']
                    locusTag = feature.qualifiers['locus_tag'][0]
                    if locusTag in locSet:
                        assert 'gene' in feature.qualifiers

                        # gene location
                        location = feature.location  # FeatureLocation: start, end, strand

                        # gene name
                        assert len(feature.qualifiers['gene']) == 1, 'file: %s' % (gbkFile)
                        geneName = feature.qualifiers['gene'][0].strip()
                        try:
                            assert len(feature.qualifiers['translation']) == 1, 'geneName: %s; file: %s' \
                                                                                % (geneName, gbkFile)
                        except:
                            print 'skipped: geneName: %s; in file: %s' % (geneName, gbkFile)
                            continue

                        # protein
                        protein = feature.qualifiers['translation'][0]
                        assert len(feature.qualifiers['transl_table']) == 1, 'file: %s' % (gbkFile)
                        translTable = int(feature.qualifiers['transl_table'][0])
                        seq = str(record.seq)

                        if str(location).startswith('join'):
                            locationList = getLocationList(str(location))
                            # print location
                        else:
                            locationList = [location]
                        seq = Seq(getRegions(seq, locationList), generic_dna)

                        try:
                            protein2 = seq.translate(table=translTable, stop_symbol='', cds=True)
                        except Exception as e:
                            print(e.message)
                            protein2 = 'X'

                        if not isLocationContinuous(locationList):
                            print('The location "%s" is not continuous for gene "%s" in file "%s"'
                                  % (location, geneName, gbkFile))
                        elif str(protein) != str(protein2):
                            print("Problem! translated proteins doesn't match: \n%s\n%s\n%s: %s\n"
                                  % (protein, protein2, geneName, gbkFile))
                        else:
                            geneCount += 1
                            # print str(round((float(geneCount)/len(locSet)) * 100, 0)) + '%\r',
                            # this gene exists in multiple copies
                            if geneName in genesDict:
                                multipleCopyGeneSet.add(geneName)
                            genesDict[geneName] = GeneRecord(geneName, locationList, seq, protein, translTable,
                                                             locusTag, recDict['acc_version'])

            # delete genes existing in multiple copies
            for geneName in multipleCopyGeneSet:
                del genesDict[geneName]

            recDict['genes_annotation'] = genesDict

        # print 'subannotations:', record.subannotations

        # append the record dictionary
        retList.append(recDict)

    return retList
