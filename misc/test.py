from Bio.Seq import Seq
from com.taxonomy_ncbid import TaxonomyNcbi
import com.fasta as fas

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
    b=[]

    print(uu + ss + bb)

    return int(2345)


def taxPlay(db):
    """
        @param db: taxonomy
            @type db: TaxonomyNcbi
    """
    #db.getScientificName()
    pass

def stat():
    file = '/Volumes/hera - net/metagenomics/projects/PPSmg/data/mercier50/nobackup/seed0/soap_.contig'
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