import re

from Bio.Seq import Seq

from algbioi import com as fas
from algbioi.com import csv


def sayHi(uu, ss, bb):
    """
    General description.

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
    a = Seq('ATGC')
    print a

    b='jkgh'

    print(uu + ss + bb)

    return int(2345)


def test2():
    s = 'lsuparc_silva106_ncbitax.bacteria+archaea.tax'
    print s[(s.rindex('.', 0, s.rindex('.')) + 1):s.rindex('.')]
    print 'done'

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


if __name__ == "__main__":
    #stat()
    #test2()
    #print toPercent([12,24,50,209,3], 2)
    filterOutSequences('/Users/ivan/Documents/work/binning/data/mercier51Strains/contigs_soapdenovo-20121119.fna',
                       '/Users/ivan/Documents/work/binning/data/mercier51Strains/binning_soapdenovo-20121119.tax',
                       '/Users/ivan/Documents/work/binning/data/mercier51Strains/contigs_soapdenovo-20121119_1000bp.fna',
                       '/Users/ivan/Documents/work/binning/data/mercier51Strains/binning_soapdenovo-20121119_1000bp.tax')