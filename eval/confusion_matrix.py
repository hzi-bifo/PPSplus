#!/usr/bin/env python

"""
This module is planned to replacement for the java script!!!
"""

import argparse




class ConfusionMatrix():
    def __init__(self, dbFile):
        pass




#Main function
def _main():
    parser = argparse.ArgumentParser(description='''Computes the confusion matrix.''',
                                     epilog='''
                                     ''')

    parser.add_argument('-f', '--fasta', nargs=1, type=file, required=True,
                        help='Fasta file.', metavar='contigs.fna',
                        dest='f')

    parser.add_argument('-r', '--reference', nargs=1, type=file, required=True,
                        help='Tab separated reference prediction file (first column ~ sequence_id; last column ~ predicted ncbid).', metavar='ref.csv',
                        dest='r')

    parser.add_argument('-p', '--predictions', nargs=1, type=file, required=True,
                        help='Tab separated prediction file (first column ~ sequence_id; last column ~ predicted ncbid).', metavar='pred.csv',
                        dest='p')

    parser.add_argument('-d', '--database', nargs=1, type=file, required=True,
                        help='Database file in the sqlite3 format.', metavar='ncbitax_sqlite.db',
                        dest='d')

    parser.add_argument('-l', '--ranks', nargs=1, type=file, required=True,
                        help='Taxonomic ranks (levels) at which the comparison will be made.(default: superkingdom,phylum,class,order,family,genus,species)',
                        metavar='family,genus,species',
                        dest='l')

    parser.add_argument('-o', '--output', nargs=1, type=file, required=True,
                        help='Prefix of the output path.',
                        dest='o')

    args = parser.parse_args()

    #if args.c:
    #assert len(args.f) == 1 and len(args.p) == 1 and len(args.m) == 1 and len(args.d) == 1 #make this nicer!

    #cons = Consistency(args.f[0].name, args.p[0].name, args.m[0].name, args.d[0].name, minScaffContigCount, minScaffBpLen, cladesSet, args.a)



def _test():
    pass

if __name__ == "__main__":
    print('Not implemented yet!')
    pass
    #_test()
    #_main()