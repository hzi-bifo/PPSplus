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
"""

import os
import sys
import signal
import string
import shutil
from Bio import SeqIO
from algbioi.com import csv
from algbioi.com import gbk
from algbioi.com import parallel
# from algbioi.com import taxonomy_ncbi
import multiprocessing as mp

from algbioi.com import fasta as fas


#from FastaFileFunctions import fastaFileToDict
#from FastaFileFunctions import getSequenceToBpDict
#from TabSepFileFunctions import getMapping
#from TabSepFileFunctions import getColumnAsList
#from Common import noNewLine


def scanForPlasmids():
    """
        Reads in a genbank file from stdin, store an accession number of a sequence that contain word
        'plasmid' (ignorecase) in the record.description Store also for each sequence a list of locations
        that contain 'plasmid' (ignorecase) in feature values

        USAGE:
        nohup time zcat /local/johdro/refdata/static/ncbi-genomes-bacteria_20121122/dna.gbff.gz \
        /local/johdro/refdata/static/ncbi-draftgenomes-bacteria_20121122/dna-contigs.gbff.gz \
        /local/johdro/refdata/static/ncbi-draftgenomes-bacteria_20121122/dna-scaffolds.gbff.gz \
        /local/johdro/refdata/static/ncbi-hmp_20121016/dna-contigs.gbff.gz \
        /local/johdro/refdata/static/ncbi-hmp_20121016/dna-scaffolds.gbff.gz \
        /local/johdro/refdata/static/ncbi-refseq-microbial_56/dna.gbff.gz \
        | python /net/metagenomics/projects/PPSmg/scripts/scripts25/plasmids.py &
    """

    #plasmidAccessionFile = '/Users/ivan/Documents/nobackup/refseq/bacterialGenomes/accession_test.txt'
    plasmidAccessionFile = '/local/igregor/ref_20121122/nobackup/plasmid_accessions.txt'
#    plasmidRegionsFile = ''
    seqCount = 0
    uniqueSeqCount = 0
    plasmidSeqCount = 0

    accessionSet = set([])
#    accessionToLocationList = dict([])

    outPlasmidAccessions = csv.OutFileBuffer(plasmidAccessionFile)
#    outPlasmidLocations = OutFileBuffer(plasmidRegionsFile)
    for record in SeqIO.parse(sys.stdin, "genbank"):
        seqCount += 1

        if str(record.id) in accessionSet:
            continue
        else:
            accessionSet.add(str(record.id))

        uniqueSeqCount += 1

        #is the whole sequence a plasmid
        if 'PLASMID' in str(record.description).upper():
            outPlasmidAccessions.writeText(str(str(record.id) + '\n'))
            plasmidSeqCount += 1
#        else:
#            #get plasmid regions
#            for feature in record.features:
#
#                for key in feature.qualifiers:
#                    if str(key) not in ['misc_feature']: #check this
#                        continue
#                    val = feature.qualifiers[key]
#                    if 'PLASMID' in string.upper(str(val)):
#                        if str(record.id) in accessionToLocationList:
#                            accessionToLocationList[str(record.id)].append(str(feature.location))
#                        else:
#                            accessionToLocationList[str(record.id)] = str(feature.location)

#    for key, val in accessionToLocationList.iteritems():
#        outPlasmidLocations.writeText(str(key + '\t' + str(val).replace("'",'').replace(']','').replace('[','')   + '\n'))

    outPlasmidAccessions.close()
#    outPlasmidLocations.close()

    print 'seqCount', seqCount
    print 'uniqueSeqCount', uniqueSeqCount
    print 'plasmidSeqCount', plasmidSeqCount


def scan():
    recordCount = 0
    recordPlasmid = 0
    breakCount = 50
    featureType = set([])
    keySet = set([])

    for record in SeqIO.parse(sys.stdin, "genbank"):
        recordCount += 1

        seqId = record.id

        if string.find(record.description, 'plasmid') != -1:
            #print record.description
            recordPlasmid += 1

        for feature in record.features:
            #print feature.location, 'location'

            for key in feature.qualifiers:
                val = feature.qualifiers[key]
                if 'PLASMID' in string.upper(str(val)):
                    print feature
                    #keySet.add((feature.type, key, str(val)))
                    #keySet.add((feature.type, key))
            #print 'location', feature.location

            #print 'type', feature.type
            #if feature.type == "CDS":
            #    print feature

            #    for x in feature.qualifiers:
            #        print x, feature.qualifiers[x]

            #for xrefentry in feature.qualifiers:
            #    print xrefentry
                #    ( key, val ) = xrefentry.split(":")
                #    if key == "taxon":
                #        taxonId = int(val)
                #        break

            #for qualifier in feature.qualifiers:
            #    print qualifier

            #featureType.add(feature.type)
            #print feature

            #for xrefentry in feature.qualifiers["db_xref"]:
            #    print featurexrefentry



            #if feature.type == "source":
                #for xrefentry in feature.qualifiers["db_xref"]:
                #    ( key, val ) = xrefentry.split(":")
                #    if key == "taxon":
                #        taxonId = int(val)
                #        break

        if recordCount == breakCount:
            break

    print 'recordCount', recordCount
    #print 'recordPlasmid', recordPlasmid
    #print 'featureType', featureType
    #for key in keySet:
    #    print key




       # if seqId in seqIdSet:
       #     print seqId, 'already in set', seqId

       # seq = str(record.seq)
       # cumulativeLen += len(seq)

      #  if len(string.replace(noNewLine(seq),'N','')) == 0:
 #           zeros += 1

       # taxonId = None

       # for feature in record.features:
       #     if feature.type == "source":
       #         for xrefentry in feature.qualifiers["db_xref"]:
       #             ( key, val ) = xrefentry.split(":")
       #             if key == "taxon":
       #                 taxonId = int(val)
       #                 break
       #     if taxonId != None:
       #         break


       # if taxonId == None:
       #     print 'could not find taxonId for', seqId
       # else:
       #     taxonSet.add(taxonId)

#    print 'record count', recordCount
#    print 'seq count', len(seqIdSet)
#    print 'taxon id count', len(taxonSet)
#    if len(seqIdSet) > 0:
#        print 'avg. seq. len', cumulativeLen/len(seqIdSet)
#    print 'zeros', zeros


def extractPlasmidAccessions(gbDir, outFile):
    """

        @param gbDir: directory containing gbk files, where the taxon id is encoded in the same way as PPS reference
        @param outFile: output file containing (taxonId, tab, accession) of records whose description contain 'plasmid'
    """
    start = 1
    stop = 1000000000
    i = 0
    out = csv.OutFileBuffer(outFile)
    for f in os.listdir(gbDir):
        if not os.path.isfile(os.path.join(gbDir, f)):
            continue
        i += 1
        if i < start:
            continue
        # if f in ['911244.1.gb','933263.1.gb', '1151366.1.gb']:
        #     continue

        # taxonId = int(f.split('.')[0])

        try:
            taxonId = int(f.split('_')[1].split('.')[0])

            for record in SeqIO.parse(os.path.join(gbDir, f), "genbank"):
                if "plasmid" in str(record.description).lower():
                    out.writeText(str(taxonId) + '\t' + record.id + '\n')
        except Exception as e:
            print('An exception occured when processing: %s' % f)

        if i % 1000 == 0:
            print i

        if i == stop:
            break

    out.close()
    print 'done'


    # for record in SeqIO.parse("/Volumes/VerbatimSSD/work/plasmids/tmp/gb/139.1.gb", "genbank"):
    #     print record.id, record.description


    # handle = Entrez.efetch(db="nucleotide", id=str(intervalGi), seq_start=str(intervalFrom),
    #                                seq_stop=str(intervalTo), strand=strandStr, rettype="gb", retmode="xml")
    # dnaXml = handle.read()
    # print(handle.read())

    # handle = Entrez.efetch(db="nucleotide", id="NZ_AYVY01000002.1", rettype="gb", retmode="text")
    # l = SeqIO.parse(handle, "genbank")
    # for i in l:
    #     print "i.id: " + i.id
    #     print "i.name: " + i.name
    #     print "i.description: " + i.description
    #     print "i.annotations: " + i.annotations
        #
        # print "i.letter_annotations: " + i.letter_annotations
        # print "i.features: " + i.features

    # id_list = ["NZ_AYVY01000049.1", "NZ_AYVY01000044.1", "NZ_AYVY01000045.1"]
    # search_results = Entrez.read(Entrez.epost("pubmed", id=",".join(id_list)))
    # webenv = search_results["WebEnv"]
    # query_key = search_results["QueryKey"]
    # print webenv
    # print query_key

    # print dnaXml


def test1():
    # handle broken pipes
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

    if len(sys.argv) != 3:
        print 'Needs two parameters: a directory containing dbk files, an output file'
    extractPlasmidAccessions(sys.argv[1], sys.argv[2])

    # scanForPlasmids()
    #scan()


def extractPlasmidAccession2Part(gbkFile, taxFile, outDir, taxonId=None):  # new !

    try:
        # tax = taxonomy_ncbi.TaxonomyNcbi(taxFile)
        last = ''  # TODO: add this as a function to algbioi.com.gbk
        for line in open(gbkFile):
            line = line.strip()
            if len(line) > 0:
                last = line
        if '//' != last:
            print('Unusual end of file: |%s| In file: |%s|, file skipped!' % (last, gbkFile))
            return
        # print('Proc: %s' % os.path.basename(gbkFile))
        tax = None
        recList = gbk.readFromGbkFile(gbkFile, tax, ('description', 'acc_version'))
        lineList = []
        for rec in recList:
            if "plasmid" in rec['description'].lower():
                accVersion = rec.get('acc_version')
                # taxonId2 = rec.get('taxonId')
                # if taxonId is None:
                #     taxonId = taxonId2

                # if taxonId is not None and taxonId2 is not None and int(taxonId) != int(taxonId2):
                #     print('%s %s %s' % (taxonId, taxonId2, accVersion))

                if accVersion is not None and taxonId is not None:

                    line = '%s\t%s\n' % (taxonId, accVersion)
                    lineList.append(line)
                else:
                    print('%s %s' % (taxonId, accVersion))

        if len(lineList) > 0:
            out = csv.OutFileBuffer(os.path.join(outDir, 'plasmids_' + os.path.basename(gbkFile)))
            out.writeText(''.join(lineList))
            out.close()

        # tax.close()
    except Exception as e:
        print 'Exception', type(e), e.message, gbkFile, taxFile, taxonId


def extractPlasmidAccession2(gbkDir, taxFile, outDir):  # new !

    assert os.path.isdir(gbkDir)
    assert os.path.isfile(taxFile)
    assert os.path.isdir(outDir)

    taskList = []
    for f in os.listdir(gbkDir):
        fPath = os.path.join(gbkDir, f)
        if os.path.isfile(fPath):
            try:
                taxonId = int(f.split('_')[1].split('.')[0])
            except Exception:
                taxonId = None

            taskList.append(parallel.TaskThread(extractPlasmidAccession2Part, (fPath, taxFile, outDir, taxonId)))
            # rList.append(extractPlasmidAccession2Part(fPath, taxFile, taxonId))

    parallel.runThreadParallel(taskList, min(15, mp.cpu_count()), keepRetValues=False)


def filterOutPlasmids(inDir, outDir, plasmidFile):  # new !
    assert os.path.isdir(inDir)
    assert os.path.isdir(outDir)
    assert os.path.isfile(plasmidFile)

    # map: taxonId -> set-of-acc
    taxonIdToAccSet = {}

    for line in open(plasmidFile):
        line = line.strip()
        tokens = line.split('\t')
        if len(tokens) >= 2:
            taxonId = int(tokens[0])
            acc = tokens[1]

            if taxonId not in taxonIdToAccSet:
                taxonIdToAccSet[taxonId] = set([acc])
            else:
                taxonIdToAccSet[taxonId].add(acc)
    filtered = 0
    for inFile in os.listdir(inDir):
        inFilePath = os.path.join(inDir, inFile)
        if os.path.isfile(inFilePath):
            taxonId = int(inFile.split('.', 1)[0])
            outFilePath = os.path.join(outDir, '%s.1.fna' % taxonId)
            if taxonId in taxonIdToAccSet:
                notAllowedSet = taxonIdToAccSet[taxonId]
                out = csv.OutFileBuffer(outFilePath)
                for seqId, seq in fas.fastaFileToDictWholeNames(inFilePath).iteritems():
                    if seqId not in notAllowedSet:
                        out.writeText('>%s\n%s\n' % (seqId, seq))
                    else:
                        filtered += 1
                out.close()
            else:
                shutil.copy2(inFilePath, outFilePath)
    print('Filtered seq: %s' % filtered)


def _extractPlasmids():
    # gbkDir = '/home/igregor/Documents/work/ppsp/bac_arch/gbk2'
    gbkDir = '/media/igregor/verbatim/work/ppsp_ref/ref_cami1/bac_arch/gbk2'
    taxFile = '/home/igregor/Documents/work/taxonomy/NCBI201502/ncbitax_sqlite.db'
    outDir = '/home/igregor/Documents/work/ppsp/bac_arch/plasmids'

    # taxFile = '/net/metagenomics/projects/PPSmg/database/NCBI201502/nobackup/mg5/taxonomy/ncbitax_sqlite.db'
    # gbkDir = '/local/igregor/ref_cami_1/bac_arch/gbk2'
    # outDir = '/local/igregor/ref_cami_1/bac_arch/plasmids'
    # outFile = '////local/igregor/ref_cami_1/bac_arch/plasmids2.csv' TODO: don't use

    # gbkDir = '/net/metagenomics/projects/PPSmg/database/NCBI201502/nobackup/gbk2'
    # outDir = '/net/metagenomics/projects/PPSmg/database/NCBI201502/nobackup/plasmids'

    extractPlasmidAccession2(gbkDir, taxFile, outDir)


def _filterOutPlasmids():

    inDir = '/local/igregor/ref_cami_1/bac_arch_centroids'
    outDir = '/local/igregor/ref_cami_1/bac_arch/bac_arch_no_plas'
    plasmidFile = '/local/igregor/ref_cami_1/bac_arch/plasmids.txt'
    filterOutPlasmids(inDir, outDir, plasmidFile)

# MAIN
if __name__ == "__main__":
    # test1()
    # _extractPlasmids()
    _filterOutPlasmids()
