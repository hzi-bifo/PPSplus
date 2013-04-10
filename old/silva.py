#!/usr/bin/env python

"""
To handle the Silva database, exported FASTA (http://www.arb-silva.de/download/arb-files/)
"""

import re
import os

from com.csv import OutFileBuffer
from com.csv import getMapping
from com.csv import getColumnAsList
from com.fasta import fastaFileToDictWholeNames
from com.taxonomy_ncbi import TaxonomyNcbi
from core.taxonomy import Taxonomy


def exportedSilvaFastaToMothurFormat(inFastaFile, outDir, taxonomyNcbi, taxonomy, accessionToNcbiFile):
    """
        input exported silva fasta file

        Takes an exported fasta file from the Silva download page and transforms it to files that can be used
        used for the Mothur Bayesian classifier
    """
    seqIdToSeq = fastaFileToDictWholeNames(inFastaFile)
    accessionToNcbi = getMapping(accessionToNcbiFile, 0, 1, sep='\t', comment = '#')

    newFasta = os.path.join(outDir, os.path.basename(inFastaFile))
    newFastaTax = str(os.path.join(outDir, os.path.basename(inFastaFile)) + '.tax')
    fastaOut = OutFileBuffer(newFasta)
    fastaTaxOut = OutFileBuffer(newFastaTax)

    counterEukaryota = 0
    counterProblem = 0
    counterOk = 0
    counterNonDef = 0
    counterNoPath = 0

    for seqId in seqIdToSeq:

        #print counter
        seq = seqIdToSeq[seqId]
        if re.match(r'.*Eukaryota.*',seqId):
            #print 'skipped', seqId
            counterEukaryota += 1
            continue
        if not re.match(r'.*Bacteria.*',seqId) and not re.match(r'.*Archaea.*',seqId):
            counterNonDef += 1
            continue
        accessionNumber = re.sub(r'^([^\.]+)\.[0-9]+\.[0-9]+ [^\ ].*$', r'\1', seqId)
        start = re.sub(r'^[^\.]+\.([0-9]+)\.[0-9]+ [^\ ].*$', r'\1', seqId)
        stop = re.sub(r'^[^\.]+\.[0-9]+\.([0-9]+) [^\ ].*$', r'\1', seqId)
        #lineage = re.sub(r'^[^\.]+\.[0-9]+\.[0-9]+ ([^\ ].*)$', r'\1', seqId)
        #name = lineage.split(';')[-1]
        if accessionNumber not in accessionToNcbi:
            counterProblem += 1
            if counterProblem < 10:
                print seqId
            continue
        else:
            assert len(accessionToNcbi[accessionNumber]) == 1
            ncbid = accessionToNcbi[accessionNumber][0]




        #ncbid = taxonomyNcbi.getNcbid2(name, checkRank=False)
        #if ncbid == None:
        #    ncbid = 'NA'
        #    i = -2
        #    ncbid2 = None
        #    while ncbid2 == None:
        #        ncbid2 = taxonomyNcbi.getNcbid2(lineage.split(';')[i], checkRank=False)
        #        i -= 1
        #    taxonomyPath = taxonomy.getPathToRootSemicolonSeparated(ncbid2)
        #    print seqId, taxonomyPath
        #    counterProblem += 1
        #else:
        #    taxonomyPath = taxonomy.getPathToRootSemicolonSeparated(ncbid)

        seq = seq.replace('U','T').replace(' ','').replace('\r','')
        name = str(accessionNumber + ':' + start + '-' + stop + '|ncbid:' + str(ncbid))
        try:
            taxonomyPath = taxonomy.getPathToRootSemicolonSeparated(ncbid)
        except:
            print 'no path for', ncbid
            counterProblem += 1
            counterNoPath += 1
            continue


        #print 'seqId:', seqId
        #print 'accessionNumber:', accessionNumber
        #print 'start:', start
        #print 'stop:', stop
        #print 'lineage:',lineage
        #print 'scientificName:',scientificName
        #print 'ncbid:', ncbid
        #print 'taxonomyPath:', taxonomyPath
        #print 'seq:', seq
        #print 'name:', name

        fastaOut.writeText(str('>' + name + '\n' + seq + '\n'))
        fastaTaxOut.writeText(str(name + '\t' + taxonomyPath + '\n'))
        counterOk += 1
        #if counterOk % 1000 == 0:
        #    print counterOk #, counterProblem, counterEukaryota, counterNonDef

    print 'ok:', counterOk
    print 'problem:', counterProblem
    print 'eukaryota:', counterEukaryota
    print 'counterNonDef:', counterNonDef
    print 'counterNoPath', counterNoPath
    fastaOut.close()
    fastaTaxOut.close()




def test():
    #silva 111
    #inFastaFile = '/Users/ivan/Documents/work/binning/database/silva111/downloaded/LSURef_111_tax_silva_trunc.fasta'
    #inFastaFile = '/Users/ivan/Documents/work/binning/database/silva111/downloaded/SSURef_111_NR_tax_silva_trunc.fasta'

    #silva108
    #inFastaFile = '/Users/ivan/Documents/work/binning/database/silva108/LSURef_108_tax_silva_trunc.fasta'
    inFastaFile = '/Users/ivan/Documents/work/binning/database/silva108/SSURef_108_NR_tax_silva_trunc_v2.fasta'

    accessionToNcbiFile = '/Users/ivan/Documents/work/binning/database/silva2Ncbi/silva2ncbi/silva2ncbi.map'
    outDir = '/Users/ivan/Documents/work/binning/database/silva108/mothurFormat'
    taxonomyNcbi = TaxonomyNcbi('/Users/ivan/Documents/work/binning/taxonomy/20120828/ncbitax_sqlite.db')
    taxonomy = Taxonomy('/Users/ivan/Documents/work/binning/taxonomy/20120828/ncbitax_sqlite.db',
                        taxonomicRanks = ['superkingdom','phylum','class','order','family','genus','species'])
    exportedSilvaFastaToMothurFormat(inFastaFile, outDir, taxonomyNcbi, taxonomy,accessionToNcbiFile)
    taxonomyNcbi.close()
    taxonomy.close()

def checkMapping():

    accessionToNcbiFile = '/Users/ivan/Documents/work/binning/database/silva2Ncbi/silva2ncbi/silva2ncbi.map'
    #acessionsListFile = '/Users/ivan/Documents/work/binning/database/silva111/downloaded/LSURef_111_accession_numbers.txt'
    acessionsListFile = '/Users/ivan/Documents/work/binning/database/silva111/downloaded/SSURef_111_NR_accession_numbers.txt'

    accessionToNcbiDict = getMapping(accessionToNcbiFile, 0, 1, sep='\t', comment = '#')
    accessionList = getColumnAsList(acessionsListFile, entryModifyFunction=None, colNum=0, sep=None, comment='#')

    notFoundCount = 0
    for i in accessionList:
        if i not in accessionToNcbiDict:
            if notFoundCount < 10:
                print i
            notFoundCount += 1

    print 'not found', notFoundCount, 'out of', len(accessionList)


if __name__ == "__main__":
    test()
    #checkMapping()

