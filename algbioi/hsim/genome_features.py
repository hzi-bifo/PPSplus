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


    Manage the strains - reference data for a species.

    NCBI genomes and draft genomes downloaded using commands:
    wget -r --no-parent ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Escherichia_coli*
    wget -r --no-parent ftp://ftp.ncbi.nih.gov/genomes/Bacteria_DRAFT/Escherichia_coli*


    Other links:
    http://www.ncbi.nlm.nih.gov/books/NBK44863/
    http://www.ncbi.nlm.nih.gov/genome/?term=escherichia+coli

    # Batch entrez
    http://www.ncbi.nlm.nih.gov/sites/batchentrez
"""
import os
import sys
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from algbioi.com import taxonomy_ncbi
from algbioi.com import csv

NCBI_TAXONOMY_FILE = '/Users/ivan/Documents/work/binning/taxonomy/20140916/ncbitax_sqlite.db'
REFERENCE_DIR_ROOT = '/Volumes/VerbatimSSD/work/hsim01'

#


def getBacteriaGenomesFileList(taxonomy):
    """
        @type taxonomy: TaxonomyNcbi
    """
    for directory in os.listdir(REFERENCE_DIR_ROOT):  # for each species
        # is the directory an integer
        try:
            taxonId = int(directory)
        except ValueError:
            continue

        # does the directory represent species
        if taxonomy.getRank(taxonId) != 'species':
            print('The directory taxon id (%s) is not a species but "%s"!' % (taxonId, taxonomy.getRank(taxonId)))
            continue

        speciesName = taxonomy.getScientificName(taxonId)
        refGenomeDir = os.path.join(REFERENCE_DIR_ROOT, directory, 'bacteria_genomes')
        # refDraftGenomeDir = os.path.join(REFERENCE_DIR_ROOT, directory, 'bacteria_draft')
        genomeCount = len(os.listdir(refGenomeDir))
        # draftCount = len(os.listdir(refDraftGenomeDir))

        print("Species: %s" % speciesName)
        print("Genome count: %s" % genomeCount)
        # print("Draft genome count: %s" % draftCount)

        # check genomes directory
        gbkCompleteGenomeList = []
        for genomeDir in os.listdir(refGenomeDir):  # for each genome in a directory
            genomeDirPath = os.path.join(refGenomeDir, genomeDir)
            gbkCompleteGenome = None
            print('Processing: %s' % genomeDir)
            for f in os.listdir(genomeDirPath):

                if f.split('.')[-1] == 'gbk':  # for each gbk file
                    pass
                    # check whether contains "complete genome"
                    records = []
                    fPath = os.path.join(genomeDirPath, f)
                    for record in SeqIO.parse(fPath, "genbank"):
                        records.append(record)

                    if len(records) > 1:
                        print('File "%s" contains more than one record! (took first record)' % fPath)
                    assert len(records) > 0, 'No record contained in file "%s' % fPath
                    record = records[0]
                    if 'complete genome' in str(record.description).lower():
                        if gbkCompleteGenome is None:
                            gbkCompleteGenome = fPath
                        else:
                            print('There are at least two files describing a "complete genome": "%s" and "%s"' % (gbkCompleteGenome, fPath))


            if gbkCompleteGenome is None:
                print('No file describing a complete genome found in "%s"' % genomeDirPath)
            else:
                # gbkCompleteGenomeList.append(gbkCompleteGenome)

                # check whether there are corresponding files for this record
                # gbkCompleteGenome
                d = {}
                for suffix in ['ptt', 'rpt', 'fna', 'gbk']:

                    f = ".".join(gbkCompleteGenome.split('.')[0:-1]) + '.' + suffix
                    if not os.path.isfile(f):
                        print('File "%s" does not exists!')
                    d[suffix] = f
                gbkCompleteGenomeList.append(d)

            # write a file with the paths!!!
            # directory
            out = csv.OutFileBuffer(os.path.join(REFERENCE_DIR_ROOT, directory, 'bacteria_genomes_sources.csv'))
            for d in gbkCompleteGenomeList:
                out.writeText(d.get('gbk') + ',' + d.get('ptt') + ',' + d.get('fna') + ',' + d.get('rpt') + '\n')
            out.close()


def readFromGbkFile(gbkFile, listOfFields=['gene_dict']):
    """
        http://biopython.org/DIST/docs/api/Bio.SeqFeature.SeqFeature-class.html

        Functionalities:
        -description: record.description (DEFINITION)
        -accession: record.annotations["accessions"] (ACCESSION) take the first one if more available
        -gi: record.annotations["gi"] (GI)
        -acc_version: record.id (VERSION), accesion version
        -taxonId: feature source, feature.qualifiers['db_xref']
        -seq: DNA sequence, str(record.seq)
        -seq_rev: reverse complement DNA sequence, str(record.reverse_complement().seq)
        -gene_name_list: list o gene names, feature gene, feature.qualifiers['gene'], (not including unnamed!)
        -gene_dict: dict of dict of genes, feature CDS (not including unnamed)

        @return: list of a dictionary of entries from the input list of fields for each record
        @rtype: dict
    """
    retList = []
    for record in SeqIO.parse(gbkFile, "genbank"):
        retDict = {}

        for item in listOfFields:
            if item == 'description':
                retDict[item] = record.description
            elif item == 'accession':
                accessionList = record.annotations['accessions']
                retDict[item] = accessionList[0]
                if len(accessionList) > 0:
                    print('There are more than one accession: "%s", the first was taken!' % accessionList)
            elif item == 'gi':
                retDict[item] = record.annotations['gi']
            elif item == 'acc_version':
                retDict[item] = record.id
            elif item == 'taxonId':
                for feature in record.features:
                    if feature.type == 'source':
                        dbxref = feature.qualifiers['db_xref']
                        assert len(dbxref) == 1, 'TaxonId: this list should be of length 1: "%s"' % dbxref
                        taxon, taxonId = dbxref[0].split(':')
                        assert taxon == 'taxon', 'This string should be "taxon", not "%s"' % taxon
                        retDict[item] = int(taxonId)
                        break
                if item not in retDict:
                    print('TaxonId not found for "%s" in "%s"' % (record.id, gbkFile))
            elif item == 'seq':
                retDict[item] = str(record.seq)
            elif item == 'seq_rev':
                retDict[item] = str(record.reverse_complement().seq)
            elif item == 'gene_name_list':
                retDict[item] = []
                for feature in record.features:
                    if feature.type == 'gene' and 'gene' in feature.qualifiers:
                        retDict[item].append(feature.qualifiers['gene'])
            elif item == 'gene_dict':
                genesDict = {}
                for feature in record.features:
                    if feature.type == 'gene' and 'gene' in feature.qualifiers:
                        # genesDict[feature.qualifiers['gene']] = {}
                        print feature.location

                        # print genesDict


            elif item == 'gene_dict__':
                dd = {}
                for feature in record.features:
                    if feature.type == 'CDS' and 'gene' in feature.qualifiers:
                        d = {}
                        gene = feature.qualifiers['gene']
                        assert len(gene) == 1, 'This list should be of length 1: "%s"' % gene
                        gene = gene[0]
                        # d['gene'] = gene

                        loc = str(feature.location).replace(']', ':', 1).replace('(', '', 1).strip('[)').split(':')
                        try:
                            d['start'] = int(loc[0]) - 1
                            d['stop'] = int(loc[1]) - 1
                            d['strand'] = loc[2]
                        except:
                            print gene, feature.location, '-------------------'


                        d['prot'] = feature.qualifiers['translation']
                        translTable = feature.qualifiers['transl_table']
                        assert len(translTable) == 1, 'This list should be of length 1: "%s"' % translTable
                        d['translTable'] = int(translTable[0])


                        # assert gene not in dd # !!!!!
                        dd[gene] = d
                retDict[item] = dd

            else:
                print('Item "%s" not supported!' % item)
                continue





        # print 'record.name (ACCESSION, LOCUS):', record.name
        # print 'record.annotations["accessions"] (ACCESSION):', record.annotations['accessions']
        # print 'record.id (VERSION):', record.id
        # print 'record.annotations["gi"] (GI):', record.annotations['gi']
        # for feature in record.features:
        #     if feature.type == 'source':
        #         print 'taxonId:', int(feature.qualifiers['db_xref'][0].split(':')[1])
        #
        #
        #
        # print 'seqLen:', len(record.seq)
        # print 'reverse complement len:', len(record.reverse_complement().seq)
        # print 'dbxrefs:', record.dbxrefs
        # print 'features:', type(record.features[0])
        # idx = 0
        #
        # for feature in record.features:
        #
        #     if feature.type == 'CDS' and 'gene' in feature.qualifiers:
        #         idx += 1
        #
        #         print 'gene (%s): %s' % (idx, feature.qualifiers['gene'])
        #         print feature.qualifiers['translation']  # take this translate to DNA using appropriate translation table, use the coordinates to get the string! and compare!!!
        #         print feature.location  # take this translate to DNA using appropriate translation table, use the coordinates to get the string! and compare!!!
        #         break
        #


                # break
                # print features
                # if idx > 10:
                #     break

        # try:
        #     pass
        # except:
        #     pass
        # try:
        #     print 'id:', record.ID
        # except:
        #     pass
        # try:
        #     print 'annotation:', record.annotation
        # except:
        #     pass
        # try:
        #     print 'annotations:', record.annotations
        # except:
        #     pass
        # try:
        #     print 'subannotations:', record.subannotations
        # except:
        #     pass

        retList.append(retDict)
        # print retList
    return retList


# get genome meta data
def getGenomeMetadata():
    for directory in os.listdir(REFERENCE_DIR_ROOT):  # for each species
        sourcesFilePath = os.path.join(REFERENCE_DIR_ROOT, directory, 'bacteria_genomes_sources.csv')
        if os.path.isfile(sourcesFilePath):
            print('Processing: %s' % sourcesFilePath)

            out = csv.OutFileBuffer(os.path.join(REFERENCE_DIR_ROOT, directory, 'bacteria_genomes_metadata.csv'))
            out.writeText('# Accession, AccessionVersion, GI, TaxonId, TaxonomyVersion, DataSource')
            recordCount = 0
            for line in open(sourcesFilePath):
                recordCount += 1
                gbk, ptt, fna, rpt = line.split(',')

                # accession =

                print line


            out.close()
            print('Records stored: %s' % recordCount)



        # is the directory an integer
        # try:
        #     taxonId = int(directory)
        # except ValueError:
        #     continue

        # does the directory represent species
        # if taxonomy.getRank(taxonId) != 'species':
        #     print('The directory taxon id (%s) is not a species but "%s"!' % (taxonId, taxonomy.getRank(taxonId)))
        #     continue



if __name__ == "__main__":
    # Escherichia coli, 562, species

    assert os.path.isfile(NCBI_TAXONOMY_FILE), "The taxonomy file does not exists!"
    assert os.path.isdir(REFERENCE_DIR_ROOT), "The root of the reference does not exist!"
    taxonomy = taxonomy_ncbi.TaxonomyNcbi(NCBI_TAXONOMY_FILE, considerNoRank=True)

    if False:
        getBacteriaGenomesFileList(taxonomy)

    if False:
    # if True:
        getGenomeMetadata()

    if True:
        readFromGbkFile('/Volumes/VerbatimSSD/work/hsim01/562/bacteria_genomes/Escherichia_coli__BL21_Gold_DE3_pLysS_AG__uid59245/NC_012947.gbk')


    taxonomy.close()










# ---------------------------------------------------------

def downloadStrains(taxonId, taxonomy, directory):
    """
        Given a species ncbi taxon id, download all strain sequences of that species.

        @param taxonId:
        @param taxonomy:
        @type taxonomy: TaxonomyNcbi
        @param directory:

    """
    pass

    # is the taxon id a species id?
    assert taxonomy.getRank(taxonId) == 'species', "Argument taxon Id (%s) is not a species!" % taxonId
    assert os.path.isdir(directory), 'Path "%s" is not a directory!' % directory


    # get all taxon ids of the species, print out how many, stop if no found
    strainList = taxonomy.getChildrenNcbids(taxonId)
    if strainList is None:
        print('No strains found for species "%s"' % taxonId)
        return
    print('There are %s strains for species "%s" in the taxonomy' % (len(strainList), taxonId))


    for strainId in strainList:
        print strainId

        break
    # a different approach needed!!!





def testDownloadStrains():
    # Escherichia coli, 562, species
    taxonId = 562
    taxonomy = taxonomy_ncbi.TaxonomyNcbi('/Users/ivan/Documents/work/binning/taxonomy/20140916/ncbitax_sqlite.db', considerNoRank=True)
    directory = '/Users/ivan/Documents/work/binning/database/hsim'

    downloadStrains(taxonId, taxonomy, directory)

    taxonomy.close()


