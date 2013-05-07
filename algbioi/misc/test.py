import re
import os
import sys

from Bio.Seq import Seq

from algbioi.com import fasta as fas
from algbioi.com import csv
from algbioi.com import taxonomy_ncbi as tax
import traceback


def sayHello2(a, b, d=4, e=5, c=3):
    """
    General description.

    @param t:
    @param n:
    @param uu: uu description
        @type uu: int
    @param ss: ssu param
        @type ss: str
    @param bb: bb something
        type bb: float
    """

    s = set()

    s.add('a')
    print s
    s.add('b')
    s = Seq('ATGC')
    print a

    #b='jkgh'

    print("%s %s %s %s %s" % (a, b, c, d, e))

    return int(2345)



def test2():
    print 'new test'
    # s = 'lsuparc_silva106_ncbitax.bacteria+archaea.tax'
    # print s[(s.rindex('.', 0, s.rindex('.')) + 1):s.rindex('.')]
    # print 'done'
    sayHello2('a', 'b', 44, e=33, c=33)
    # print 'more'
    sayHello2('a','a','a')



def stat():
    """


    """
    #file = '/Volumes/hhu-hera/PPSmg/data/nobackup/mercier51Strains/contigs_soapdenovo-20121119.fna'
    #file = '/Volumes/hhu-hera/PPSmg/data/mercier50/nobackup/seed0/soap_.contig'
    #file = '/Volumes/hhu-hera/PPSmg/data/mercier50/nobackup/seed5/soap_.contig'
    #file = '/Volumes/hhu-hera/PPSmg/data/mercier50/nobackup/seed5_63mer/soap_.contig'
    #file = '/Volumes/hhu-hera/PPSmg/data/mercier50/nobackup/seed5_127/soap_seed5_127.contig'
    #fas.fastaFileToDict(file)
    file = '/Volumes/hhu-hera/PPSmg/data/mercier50/nobackup/seed0/ReadsR_seed0.fa'
    c = 0
    bp = 0
    minLen = 1000
    maxLen = 0
    totalBp = 0
    totalCount = 0
    for k, v in fas.fastaFileToDict(file).iteritems():
        l = len(v)
        totalCount += 1
        if l > 1000:
            c += 1
            bp += l
        totalBp += l
        if l < minLen:
            minLen = l
        elif l > maxLen:
            maxLen = l

    print('Bigger than 1000bp (contigs, bp):', c, bp)
    print('maxLen, minLen, avgLen:', minLen, maxLen, (totalBp / totalCount))
    print('total:', totalCount)





def toPercent(costList, costIdxHundredPercent=2):
    """ For PTree, relative cost comparison. """
    percent = costList[costIdxHundredPercent] / 100.0
    return map(lambda x: round(x / percent, 3), costList)


def toPercent2(timeStrList, timeIdxHundredPercent=2):
    """ For PTree, relative time comparison. """
    timeList = timeStrList.split(';')
    secList = []
    for t in timeList:
        m = 0
        s = 0
        h = 0
        for i in t.split(' '):
            if 's' in i:
                s = float(re.sub(r'([0-9\.]+)s', r'\1', i))
            elif 'm' in i:
                m = float(re.sub(r'([0-9]+)m', r'\1', i))
            elif 'h' in i:
                h = float(re.sub(r'([0-9]+)h', r'\1', i))
            else:
                assert False, 'Unknown entry'
        sec = s + (m * 60) + (h * 60 * 60)
        secList.append(sec)
    percent = secList[timeIdxHundredPercent] / 100.0
    return map(lambda x: round(x / percent, 3), secList)


def filterOutSequences(fastaFile, predFile, outFastaFile, outPredFile, minBp=1000):
    seqIdToSeq = fas.fastaFileToDict(fastaFile)
    seqIdToBp = fas.getSequenceToBpDict(fastaFile)
    seqIdToPred = csv.predToDict(predFile)
    outFasta = csv.OutFileBuffer(outFastaFile)
    outPred = csv.OutFileBuffer(outPredFile)

    totalBp = 0
    taken = 0
    takenBp = 0
    for seqId, seq in seqIdToSeq.iteritems():
        bp = seqIdToBp[seqId]
        totalBp += bp
        if bp >= minBp:
            taxonId = seqIdToPred[seqId]
            outFasta.writeText('>' + str(seqId) + '\n' + str(seq) + '\n')
            outPred.writeText(str(seqId) + '\t' + str(taxonId) + '\n')
            taken += 1
            takenBp += bp
    outFasta.close()
    outPred.close()

    print('Total sequences: ', len(seqIdToSeq))
    print('Total size: ', totalBp)
    print('Taken sequences: ', taken)
    print('Taken size:', takenBp)


def fastaBySeqNameList(inSeqIdList, inFastaFile, outFastaFile):
    """
        Generates outFastaFile that contains sequences that are in the inSeqIdList and inFastaFile.
    """
    seqIdList = csv.getColumnAsList(inSeqIdList)
    seqIdToSeq = fas.fastaFileToDict(inFastaFile)
    out = csv.OutFileBuffer(outFastaFile)
    for seqId in seqIdList:
        seq = seqIdToSeq.get(seqId, None)
        if seq is None:
            print("Can't find sequence for seqId: %s" % seqId)
        else:
            out.writeText('>' + str(seqId) + '\n' + str(seq) + '\n')
    out.close()


def testExeption():
    d = {}
    try:
        print (str(1/0))
    except Exception as ex:
        #print(ex.message)
        #print(ex.args)
        print traceback.print_exc(file=sys.stdout)
    print('done')



if __name__ == "__main__":
    fastaBySeqNameList('/Volumes/hhu-hera/data/CowRumen/chunked070513/chunks2000.scaffold_names',
                       '/Users/ivan/Documents/work/binning/data/CowRumen/assembly/cow_rumen_fragmented_velvet_assembly_scaffolds.fas',
                       '/Volumes/hhu-hera/data/CowRumen/chunked070513/chunks2000.scaffolds')

    #refToClades('/Volumes/hera - net/metagenomics/projects/PPSmg/data/nobackup/NCBI20121122/sequences',
    #          '/Users/ivan/Documents/nobackup/species_list.txt',
    #          '/Users/ivan/Documents/work/binning/taxonomy/20121122/ncbitax_sqlite.db',
    #          rank='species')
    #testExeption()
    #stat()
    #test2()
    #print toPercent([12,24,50,209,3], 2)
    #filterOutSequences('/Users/ivan/Documents/work/binning/data/mercier51Strains/contigs_soapdenovo-20121119.fna',
    #                   '/Users/ivan/Documents/work/binning/data/mercier51Strains/binning_soapdenovo-20121119.tax',
    #                   '/Users/ivan/Documents/work/binning/data/mercier51Strains/contigs_soapdenovo-20121119_1000bp.fna',
    #                   '/Users/ivan/Documents/work/binning/data/mercier51Strains/binning_soapdenovo-20121119_1000bp.tax')