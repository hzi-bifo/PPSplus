from Bio.Seq import Seq
from com.taxonomy_ncbid import TaxonomyNcbi
import com.fasta as fas
import re
import os

__author__ = 'ivan'

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
    b.replace



    print(uu + ss + bb)

    return int(2345)


def taxPlay(db):
    """
        @param db: taxonomy
            @type db: TaxonomyNcbi
    """


    pass

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






if __name__ == "__main__":
    #stat()
    #test2()
    print toPercent([12,24,50,209,3], 2)