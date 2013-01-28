from Bio.Seq import Seq

__author__ = 'ivan'

from sets import Set

def sayHi():
    """
    Test comment.
    """
    s = Set()
    s.add('a')
    print s
    s.add('b')
    a = Seq('ATGC')
    print a



if __name__ == "__main__":
    sayHi()