from Bio.Seq import Seq
from com.taxonomy_ncbid import TaxonomyNcbi
import com.fasta as fas
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
    file = '/Volumes/hhu-hera/PPSmg/data/mercier50/nobackup/seed0/soap_.contig'
    #file = '/Volumes/hhu-hera/PPSmg/data/mercier50/nobackup/seed5/soap_.contig'
    #file = '/Volumes/hhu-hera/PPSmg/data/mercier50/nobackup/seed5_63mer/soap_.contig'
    #file = '/Volumes/hhu-hera/PPSmg/data/mercier50/nobackup/seed5_127/soap_seed5_127.contig'
    #fas.fastaFileToDict(file)
    c = 0
    bp = 0
    for k,v in fas.fastaFileToDict(file).iteritems():
        if len(v) > 1000:
            c += 1
            bp += len(v)

    print(c, bp)



if __name__ == "__main__":
    stat()
    #test2()