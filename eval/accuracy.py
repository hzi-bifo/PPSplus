#!/usr/bin/env python

import argparse

from com.csv import predToDict
from com.fasta import getSequenceToBpDict
from com.taxonomy_ncbid import TaxonomyNcbi


class Accuracy():
    """
        Implements computation of the "precision" and "recall" according to different definitions.
    """
    def __init__(self, fastaFilePath, predFilePath, trueFilePath, databaseFile, ranks):
        self._seqToBp = getSequenceToBpDict(fastaFilePath)
        self._seqToPred = predToDict(predFilePath)
        self._seqToTrue = predToDict(trueFilePath)
        self._taxonomy = _TaxonomyWrapperA(databaseFile, ranks)

    def getAccuracy(self, rank, minFracClade, minFracPred=None, asBp=False, weightAccordingBinSize=False):
        """
            Precision (specificity) and Recall (sensitivity) according to PhyloPythiaS and PhyloPythia papers.

            The number of classes correspond to the number of classes in the true reference and param "minFracClades".

            @param rank: on which taxonomic ranks the predictions are made
            @param minFracClade: a clade is considered only if the dataset (true predictions) contain at least this
                          fraction of sequences that belong to the clade.
            @param minFracPred: a clade is considered only if the corresponding predicted bins contain at least this
                         fraction of the overall sequences (None ~ this criteria is not considered and only
                         true "reference" bins are used for the comparison).
            @param asBp: count it according to the sequence lengths
            @param weightAccordingBinSize: weight individual bins according to their bin size

            @return: [precision, recall, classPrecisionNum, classRecallNum]
        """
        predAtRankDict = self._taxonomy.getPredDictAtRank(self._seqToPred, rank)
        trueAtRankDict = self._taxonomy.getPredDictAtRank(self._seqToTrue, rank)
        tp = dict([]) # class label -> count of sequences correctly assigned to clade i
        t = dict([])  # class label -> true count of sequences of clade i
        p = dict([])  # class label -> count of sequences assigned to clade i
        tpOther = 0   # count of sequences correctly unassigned
        tOther = 0    # true count of sequences that are unassigned at given rank

        for seq in self._seqToPred:
            if asBp:
                bp = self._seqToBp[seq]
            else:
                bp = 1

            # true
            i = None
            if seq in trueAtRankDict:
                i = trueAtRankDict[seq]
                if i not in t:
                    t[i] = bp
                else:
                    t[i] += bp
            else:
                tOther += bp
                if seq not in predAtRankDict:
                    tpOther += bp

            # pred
            j = None
            if seq in predAtRankDict:
                j = predAtRankDict[seq]
                if j not in p:
                    p[j] =  bp
                else:
                    p[j] += bp

            # match
            if i == j and i is not None:
                if i not in tp:
                    tp[i] =  bp
                else:
                    tp[i] += bp

        classesR = t.keys()
        classesP = p.keys()

        #filter out least abundant TRUE clades
        sum = 0
        for i in classesR:
            sum += t[i]
        rmList = []
        for i in classesR:
            if float(t[i])/float(sum) < minFracClade:
                rmList.append(i)
        for i in rmList:
            classesR.remove(i)

        #filter out least abundant PREDICTED clades
        if minFracPred is None:
            classesP = classesR
        else:
            sum = 0
            for i in classesP:
                sum += p[i]
            rmList = []
            for i in classesP:
                if float(p[i])/float(sum) < minFracPred:
                    rmList.append(i)
            for i in rmList:
                classesP.remove(i)

        #null null elements
        for i in classesR:
            if i not in tp:
                tp[i] = 0
            if i not in p:
                p[i] = 0
        for i in classesP:
            if i not in tp:
                tp[i] = 0
            if i not in t:
                t[i] = 0

        #compute weights of individual bins
        wp = dict([]) #weights for precision
        wr = dict([]) #weights for recall
        if weightAccordingBinSize:
            #weights correspond to the number of sequences/bp assigned to individual bins (differ for precision and recall!)
            sumP = 0.0
            sumR = 0.0
            for i in classesP:
                sumP += p[i]
            for i in classesR:
                sumR += t[i]
            sumR += tOther

            for i in classesP:
                wp[i] = float(p[i])/sumP
            for i in classesR:
                wr[i] = float(t[i])/sumR
            if tOther > 0:
                wrOther = float(tOther)/sumR
        else:
            #all bins are equally important!
            for i in classesP: # precision
                wp[i] = 1.0/float(len(classesP))

            for i in classesR: # recall
                if tOther > 0:
                    w = 1.0/float(len(classesR)+1)
                    wr[i] = w
                    wrOther = w
                else:
                    wr[i] = 1.0/float(len(classesR))

        #sch = 0.0
        #for i in classesP:
        #    sch += wp[i]
        #print 'sch P', str(sch)
        #sch = 0.0
        #for i in classesR:
        #    sch += wr[i]
        #if tOther > 0:
        #    sch += wrOther
        #print 'sch R', str(sch)

        #precision
        precision = 0.0
        for i in classesP:
            if p[i] > 0:
                precision += (float(tp[i])/float(p[i]))*wp[i]

        #recall
        recall = 0.0
        for i in classesR:
            recall += (float(tp[i])/float(t[i]))*wr[i]
        if tOther != 0:
            recall += (float(tpOther)/float(tOther))*wrOther
        #
        return [precision, recall, len(classesP), len(classesR)]


    def test(self):
        pass

    def close(self):
        self._taxonomy.close()


class _TaxonomyWrapperA():
    """ Wraps the functionality of the database. """
    def __init__(self, databaseFile, ranks):
        self._taxonomy = TaxonomyNcbi(databaseFile)
        self._rankToId = dict([])
        self._ncbidToRankId = dict([])
        self._parentBuff = dict([]) # rank -> ncbid -> parent ncbid at rank
        id = 0
        for rank in ranks:
            self._rankToId[rank] = id
            self._parentBuff[id] = dict([])
            id += 1

    def _getRankId(self, ncbid):
        if ncbid in self._ncbidToRankId:
            return self._ncbidToRankId[ncbid]
        else:
            rank = self._taxonomy.getRank(ncbid)
            if rank == None:
                return None
            else:
                self._ncbidToRankId[ncbid] = self._rankToId[rank]
                return self._ncbidToRankId[ncbid]

    def getPredDictAtRank(self, seqToDict, rank):
        """
            Gets predictions at the given rank as dict based on the argument dict.
        """
        rankId = self._rankToId[rank]
        outDict = dict([])
        parentBuff = self._parentBuff[rankId]
        for seq in seqToDict:
            ncbid = seqToDict[seq]
            if ncbid in parentBuff:
                assert seq not in outDict
                outDict[seq] = parentBuff[ncbid]
            else:
                ncbidRank = self._taxonomy.getRank(ncbid)
                if ncbidRank != None and self._rankToId[ncbidRank] >= rankId:
                    if ncbidRank == rank:
                        outDict[seq] = ncbid
                    else:
                        #get parent
                        current = self._taxonomy.getParentNcbid(ncbid)
                        while self._taxonomy.getRank(current) != rank and current is not None:
                            current = self._taxonomy.getParentNcbid(current)
                        if self._taxonomy.getRank(current) == rank:
                            outDict[seq] = current
                            parentBuff[ncbid] = current
        return outDict

    def close(self):
        self._taxonomy.close()


def test():
    fastaFilePath = '/Users/ivan/Documents/work/binning/data/simMC/AMGN_AMD.Arachne.contigs.fna'
    predFilePath = '/Users/ivan/Documents/work/binning/tests/simMC/AMD05/output/AMGN_AMD.Arachne.contigs.fna.pOUT' #just marker genes
    #predFilePath = '/Users/ivan/Documents/work/binning/tests/simMC/AMD06/output/AMGN_AMD.Arachne.contigs.fna.pOUT' #with taxator
    trueFilePath = '/Users/ivan/Documents/work/binning/data/simMC/AMD.Arachne.genus'
    databaseFile = '/Users/ivan/Documents/work/binning/taxonomy/ncbi_taxonomy_20110629/ncbitax_sqlite.db'
    ranks = ['root','superkingdom','phylum','class','order','family','genus','species']
    acc = Accuracy(fastaFilePath, predFilePath, trueFilePath, databaseFile, ranks)
    #acc.test()
    print 'precision, recall, #classes precision, #classes recall, comment'
    for rank in ['superkingdom','phylum','class','order','family','genus']:
        print rank, '--------------------------'
        #p, r, cp, cr = acc.getAccuracy(rank, minFracClade=0.01, minFracPred=None, asBp=False, weightAccordingBinSize=False)
        #print str(str(round(p*100.0, 1)) + '%, ' + str(round(r*100.0, 1)) + '%, ' + str(cp) + ', ' + str(cr) + ', ~PP paper')

        #p, r, cp, cr = acc.getAccuracy(rank, minFracClade=0.01, minFracPred=None, asBp=True, weightAccordingBinSize=False)
        #print str(str(round(p*100.0, 1)) + '%, ' + str(round(r*100.0, 1)) + '%, ' + str(cp) + ', ' + str(cr) + ', ~PP paper; in "bp"')

        p, r, cp, cr = acc.getAccuracy(rank, minFracClade=0.001, minFracPred=0.001, asBp=False, weightAccordingBinSize=False)
        print str(str(round(p*100.0, 1)) + '%, ' + str(round(r*100.0, 1)) + '%, ' + str(cp) + ', ' + str(cr) + ', ~PP paper; in "seq"; P~pred classes')

        p, r, cp, cr = acc.getAccuracy(rank, minFracClade=0.001, minFracPred=0.001, asBp=True, weightAccordingBinSize=False)
        print str(str(round(p*100.0, 1)) + '%, ' + str(round(r*100.0, 1)) + '%, ' + str(cp) + ', ' + str(cr) + ', ~PP paper; in "bp"; P~pred classes')

        p, r, cp, cr = acc.getAccuracy(rank, minFracClade=0.001, minFracPred=0.001, asBp=False, weightAccordingBinSize=True)
        print str(str(round(p*100.0, 1)) + '%, ' + str(round(r*100.0, 1)) + '%, ' + str(cp) + ', ' + str(cr) + ', ~PP paper; in "seq. count"; P~pred classes; weighted bins')

        p, r, cp, cr = acc.getAccuracy(rank, minFracClade=0.001, minFracPred=0.001, asBp=True, weightAccordingBinSize=True)
        print str(str(round(p*100.0, 1)) + '%, ' + str(round(r*100.0, 1)) + '%, ' + str(cp) + ', ' + str(cr) + ', ~PP paper; in "bp"; P~pred classes; weighted bins')

    acc.close()


def main():
    parser = argparse.ArgumentParser(description='''Computes precision and recall measures according to different definitions".''',
                                     epilog='''
                                     ''')

    parser.add_argument('-f', '--fasta', nargs=1, type=file, required=True,
                        help='Fasta file.', metavar='contigs.fna',
                        dest='f')

    parser.add_argument('-p', '--predictions', nargs=1, type=file, required=True,
                        help='Tab separated prediction file (first column sequence name, last column predicted ncbid.', metavar='pred.csv',
                        dest='p')

    parser.add_argument('-t', '--true-assignments', nargs=1, type=file, required=True,
                        help='Tab separated true assignments file (first column sequence name, last column predicted ncbid.', metavar='true_assignments.csv',
                        dest='t')

    parser.add_argument('-d', '--database', nargs=1, type=file, required=True,
                        help='Database file in the sqlite3 format.', metavar='ncbitax_sqlite.db',
                        dest='d')

    parser.add_argument('-r', '--ranks', nargs=1,
                        help='Compute the measures only for these ranks (given as comma separated strings) Default ~ consider all ranks.',
                        metavar='order,family,genus',
                        dest='r')

    parser.add_argument('-c', '--min-frac-clade', nargs=1,
                        help=str('A clade is considered only if the dataset (true predictions) contain at least this ' +
                        'fraction of sequences that belong to the corresponding clade at the corresponding rank. '  +
                        '(e.g. value 0.01 means that all clades that are considered in the statistics represent at ' +
                        'least 1%% of the overall dataset) Default ~ 0.01'),
                        metavar='0.01',
                        dest='c')

    parser.add_argument('-b', '--min-frac-bin', nargs=1,
                        help='In the computation of precision. A clade is considered only if the corresponding predicted bins contain at least this' +
                        ' fraction of the overall predicted sequences at the corresponding rank (Default ~ this condition is not considered and only' +
                        ' true "reference" bins are used for the computation). ' +
                        ' Default ~ None',
                        metavar='0.01',
                        dest='b')

    parser.add_argument('-s', '--consider-seq-len', action='store_true',
                        help='Compute measures based on the sequence lengths (in bp). Default ~ consider sequence counts.',
                        dest='s')

    parser.add_argument('-w', '--weight-bins', action='store_true',
                        help='The measures are computed using weighted averages over bin sizes.' +
                        ' Size of true bins is used to compute "Recall". Size of predicted bins (if min-frac-bin is used) is used to compute "Precision". ' +
                        '(Default ~ consider all bins as equal)',
                        dest='w')

    parser.add_argument('-o', '--overview', action='store_true',
                        help='Compute the measures according to several default settings. You can still set the (-c) and (-b) options.',
                        dest='o')

    args = parser.parse_args()

    ranksAll = ['superkingdom','phylum','class','order','family','genus','species']

    if args.r:
        ranks = str(args.r[0].name).split(',')
    else:
        ranks = ranksAll

    #print args.f[0].name
    #print args.p[0].name
    #print args.t[0].name
    #print args.d[0].name
    #print ranksAll

    #minFracClade
    if args.c:
        mfc = float(args.c[0])
    else:
        mfc = 0.01

    #minFracPred (bin)
    if args.b:
        mfp = float(args.b[0])
    else:
        mfp = 0.01

    acc = Accuracy(args.f[0].name, args.p[0].name, args.t[0].name, args.d[0].name, ['root','superkingdom','phylum','class','order','family','genus','species'])

    print 'precision, recall, #classes precision, #classes recall, comment'
    for rank in ranks:
        if args.o: # overview
            print rank, '--------------------------'
            p, r, cp, cr = acc.getAccuracy(rank, minFracClade=mfc, minFracPred=mfp, asBp=False, weightAccordingBinSize=False)
            print str(str(round(p*100.0, 1)) + '%, ' + str(round(r*100.0, 1)) + '%, ' + str(cp) + ', ' + str(cr) + ', "in seq count, not weighted"')

            p, r, cp, cr = acc.getAccuracy(rank, minFracClade=mfc, minFracPred=mfp, asBp=True, weightAccordingBinSize=False)
            print str(str(round(p*100.0, 1)) + '%, ' + str(round(r*100.0, 1)) + '%, ' + str(cp) + ', ' + str(cr) + ', "in bp, not weighted"')

            p, r, cp, cr = acc.getAccuracy(rank, minFracClade=mfc, minFracPred=mfp, asBp=False, weightAccordingBinSize=True)
            print str(str(round(p*100.0, 1)) + '%, ' + str(round(r*100.0, 1)) + '%, ' + str(cp) + ', ' + str(cr) + ', in sec count"; weighted bins')

            p, r, cp, cr = acc.getAccuracy(rank, minFracClade=mfc, minFracPred=mfp, asBp=True, weightAccordingBinSize=True)
            print str(str(round(p*100.0, 1)) + '%, ' + str(round(r*100.0, 1)) + '%, ' + str(cp) + ', ' + str(cr) + ', in "bp"; weighted bins')

        else: #custom
            p, r, cp, cr = acc.getAccuracy(rank, minFracClade=mfc, minFracPred=mfp, asBp=bool(args.s), weightAccordingBinSize=bool(args.w))
            print str(str(round(p*100.0, 1)) + '%, ' + str(round(r*100.0, 1)) + '%, ' + str(cp) + ', ' + str(cr))

    acc.close()


if __name__ == "__main__":
    #if os.getcwd() == '/Users/ivan/Documents/work/python/workspace/pPPS/src':
    #    test()
    #else:
        main()