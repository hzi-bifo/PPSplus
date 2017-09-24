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

    Miscellaneous functions to test the results from the dataset from Dimidrij (Sponge assembly)
"""

import os
from algbioi.com import csv
from algbioi.com import fasta
from algbioi.com import taxonomy_ncbi

from algbioi.com import common

# FILE="/Users/ivan/Documents/nobackup/tmp8/output"

def concatenate(directory, outputFile):
    out = csv.OutFileBuffer(outputFile)
    for f in os.listdir(directory):
        path = os.path.join(directory, f)
        name = f.split('.')[0]
        seqIdToSeq = fasta.fastaFileToDict(path)
        out.writeText('>' + str(name) + '\n')
        for seqId, seq in seqIdToSeq.iteritems():
            out.writeText(str(seq) + 200*'N' + '\n')
    out.close()


def _main():
    # concatenate('/Volumes/My_Passport_Mac/work/rel_1_3/cami/assemblies',
    #             '/Volumes/My_Passport_Mac/work/rel_1_3/cami/cami_genomes_concat.fna')
    pass


def outToCami(ppspOutFile):
    """
        Creates a cami output file, in format:

        #CAMI Format for Binning
        @Task:Binning
        @Version:1.0
        @ContestantID:CONTESTANTID
        @SampleID:SAMPLEID
        @Referencebased:T
        @Assemblybased:T
        @ReplicateInfo:T

        @@SEQUENCEID	TAXID	BINID

        read1201	123	123
        read1202	123	123
        read1203	131564	131564
        read1204	562	562.1
        read1205	562	562.2

    """
    out = csv.OutFileBuffer(ppspOutFile + '.cami')
    out.writeText("""#CAMI Format for Binning
@Task:Binning
@Version:1.0
@ContestantID:CONTESTANTID
@SampleID:SAMPLEID
@Referencebased:T
@Assemblybased:T
@ReplicateInfo:T

@@SEQUENCEID	TAXID	BINID

""")
    for line in open(ppspOutFile):
        name, taxonId = line.strip('\n').split('\t', 2)
        out.writeText("%s\t%s\t%s\n" % (name, taxonId, taxonId))
    out.close()


def createCamiOut(fileDir):
    """
        Directory where the output files are located
    """
    outList = getPPSPOutPathList(fileDir)
    if len(outList) != 1:
        print('Unusual number of the output "pOUT" files detected: %s' % len(outList))
    for path in outList:
        outToCami(path)


def getPPSPOutPathList(outDir):
    retList = []
    for f in os.listdir(outDir):
        fList = f.split('.')
        if fList[-1] == 'pOUT' and fList[-2] != 'PP':
            retList.append(os.path.join(outDir, f))
    return retList


def readAssignments(assignmentFile):
    """
        Reads an assingment file, either in the cami format or in the PPS output (out) format

        @rtype: dict
        @return: mapping(name->taxonId)
    """
    if os.path.basename(assignmentFile).split('.')[-1] == 'cami':
        return readCami(assignmentFile)
    else:
        return csv.predToDict(assignmentFile)


def readCami(camiAssignFile):
    """
        Reads a file in the cami format

        @rtype: dict
    """
    ret = {}
    for line in open(camiAssignFile):
        line = line.strip()
        if not (line.startswith('#') or line.startswith('@') or len(line) == 0):
            name, taxonId = line.split('\t')[0:2]
            ret[name] = taxonId
        else:
            print line
    return ret


# get ranks at which taxon ids are defined:
def testRanks():
    # ids = '/Users/ivan/Documents/work/tmp/cami/ids.txt'
    ids = '/Users/ivan/Documents/work/tmp/cami/taxon_ids.txt'
    db = '/Users/ivan/Documents/work/binning/taxonomy/20140916/ncbitax_sqlite.db'
    taxonomy = taxonomy_ncbi.TaxonomyNcbi(db)
    ranks = set()
    for line in open(ids):
        line = line.strip()
        rank = taxonomy.getRank(line)
        while rank.strip() == 'no rank':
            line = taxonomy.getParentNcbid(line)
            rank = taxonomy.getRank(line)
        ranks.add(rank)
    print ranks
    taxonomy.close()


def dimRef():
    base = '/Users/ivan/Documents/nobackup/tmp8/dim'
    src = os.path.join(base, 'bins')
    src = os.path.join(base, 'markers')
    # dst = os.path.join(base, 'dimitrij_ref.csv')
    dst = os.path.join(base, 'dimitrij_ref_markers.csv')
    out = csv.OutFileBuffer(dst)
    for f in os.listdir(src):
        # print f
        fp = os.path.join(src, f)
        if os.path.isfile(fp):
            # print fp
            id = int(f.split('.', 1)[0])
            contigs = csv.getColumnAsList(fp)
            # print contigs
            for c in contigs:
                # print c
                out.writeText('%s\t%s\n' % (c, id))
    out.close()

def dimFasta():
    fastaIn = '/Volumes/Macintosh HD/Users/ivan/Documents/nobackup/tmp8/dim/assembly.fa'
    fastaOut = '/Volumes/Macintosh HD/Users/ivan/Documents/nobackup/tmp8/dim/assembly_dim_markers.fa'
    contigIdsSet = set(csv.getColumnAsList('/Users/ivan/Documents/nobackup/tmp8/dim/dimitrij_ref_markers.csv'))

    out = csv.OutFileBuffer(fastaOut)
    for k, v in fasta.fastaFileToDictWholeNames(fastaIn).iteritems():
        if k in contigIdsSet:
            out.writeText('>%s\n%s\n' % (k, v))
    out.close()


def dimDiscrepancyAll():
    f1 = '/Users/ivan/Documents/nobackup/tmp8/dim/discrepancy/300/test_2_2_assembly_dim_markers.fa.ids_all.mP'
    f2 = '/Users/ivan/Documents/nobackup/tmp8/dim/discrepancy/300/test_2_3_assembly_dim_markers.fa.ids_all.mP'
    f1b = '/Users/ivan/Documents/nobackup/tmp8/dim/discrepancy/1000/test2B0_assembly_dim_markers.fa.ids_all.mP'
    f2b = '/Users/ivan/Documents/nobackup/tmp8/dim/discrepancy/1000/test2B1_assembly_dim_markers.fa.ids_all.mP'
    discrepancyMothur(f1, f2)
    discrepancyMothur(f1b, f2b)

def discrepancyMothur(f1, f2):
    eList = []
    for e in open(f1):
        eList.append(e)
    for e in open(f2):
        eList.append(e)
    d = {}
    for e in eList:
        a = e.split()
        # print e
        k = '_'.join(a[0:2]) + '_' + '_'.join(a[3:5])
        v = float(a[2])
        if k in d:
            d[k].append(v)
        else:
            d[k] = [v]

    diffList = []
    sum = 0.0
    count = 0
    notInList = 0
    for k, v in d.iteritems():
        if len(v) == 2:
            difference = abs(v[0] - v[1])
            diffList.append(difference)
            sum += difference
            count += 1
        else:
            notInList += 1

    print('difference avg: %s, max: %s, notInList: %s' % (sum / float(count), max(diffList), notInList))



# def checkNovelty():
#     sampleAccSet = set(csv.getColumnAsList('/Users/ivan/Documents/nobackup/cami/gold_standard_clades.csv'))
#     novelty = csv.getMapping('/Users/ivan/Documents/nobackup/cami/metadata1000_novelty.tsv', 0, 1)
#     for k, v in novelty.iteritems():
#         assert len(v) == 1
#         v = v[0]
#         print k, v
#     for acc in sampleAccSet:
#         print('acc: %s novelty: %s' % (acc, novelty.get(acc)))

if __name__ == "__main__":
    pass
    import sys
    print sys.version
    # checkNovelty()
    # dimDiscrepancyAll()
    # dimFasta()
    # dimRef()
    # _main()
    # testRanks()
    # createCamiOut(FILE)
    # print readCami('/Users/ivan/Documents/work/tmp/cami/gold_standard_assembly.fasta.pOUT.cami')
