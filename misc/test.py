from Bio.Seq import Seq
from com.taxonomy_ncbid import TaxonomyNcbi

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

    print(uu + ss + bb)

    return int(2345)


def taxPlay(db):
    """
        @param db: taxonomy
            @type db: TaxonomyNcbi
    """
    #db.getScientificName()
    pass



if __name__ == "__main__":
    sayHi()