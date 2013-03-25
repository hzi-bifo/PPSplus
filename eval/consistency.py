#!/usr/bin/env python

import os
import sys
import argparse

from com.csv import predToDict
from com.csv import getMapping
from com.fasta import getSequenceToBpDict
from com.taxonomy_ncbid import TaxonomyNcbi


class _TaxonomyWrapper():
    """
        Wraps the taxonomy to buffer (speed up) calls.
    """
    def __init__(self, databaseFile):
        self._taxonomy = TaxonomyNcbi(databaseFile)
        self._ncbidToNcbidParent = dict([])

    def getParent(self, ncbid):
        if ncbid in self._ncbidToNcbidParent:
            return self._ncbidToNcbidParent[ncbid]
        else:
            parent = self._taxonomy.getParentNcbid(ncbid)
            self._ncbidToNcbidParent[ncbid] = parent
            return parent

    def getDist(self, ncbid, ncbidSet):
        """
            Gets the distance between ncbid and the closest clade in the ncbidSet which represents a path from the root to a clade.
        """
        current = ncbid
        dist = 0
        while (current not in ncbidSet):
            if current == None:
                sys.stderr.write('Consistency:_TaxonomyWrapper:getDist: current is "None" ' + str(ncbid) + ' ' + str(ncbidSet) + '\n')
                break
            current = self.getParent(current)
            dist += 1
        return dist

    def getDistantParent(self, ncbid, intDist):
        """
            Gets parent of "ncbid" that is in distance intDist.
        """
        parentNcbid = ncbid
        for i in range(int(intDist)):
            parentNcbid = self._taxonomy.getParentNcbid(parentNcbid)
        return parentNcbid

    def getDistTowardsRoot(self, ncbid, parentNcbid):
        """
            Gets the distance from ncbid to parentNcbid that is its parent and lies on the same path to the root.
        """
        current = ncbid
        dist = 0
        while (current != parentNcbid) and (current != 1):
            if current == None:
                sys.stderr.write('Consistency:_TaxonomyWrapper:getDistTowardsRoot: current is "None" ' + str(ncbid) + ' ' + str(parentNcbid) + '\n')
                break
            current = self._taxonomy.getParentNcbid(current)
            dist += 1
        return float(dist)

    def getScientificName(self,ncbid):
        return self._taxonomy.getScientificName(ncbid)

    def close(self):
        """ To free resources. """
        self._taxonomy.close()

class ScScaffold():
    """
        Represents one scaffold.
    """
    def __init__(self, name, ncbid, pathSet, weightedDistToPath, contigsNameList, contigNameToNcbid,
                 ncbidToBp, ncbidToWeight, ncbidToDist, ncbidToLeafDist):
        self._name = name
        self._ncbid = ncbid
        self._pathSet = pathSet
        self._weightedDistToPath = weightedDistToPath
        self._contigsNameList = contigsNameList
        self._contigNameToNcbid = contigNameToNcbid
        self._ncbidToBp = ncbidToBp
        self._ncbidToWeight = ncbidToWeight
        self._ncbidToDist = ncbidToDist
        self._ncbidToLeafDist = ncbidToLeafDist


    def getName(self):
        return self._name

    def getNcbid(self):
        """
            Returns ncbid according to which the consistency of this Scaffold was computed.
        """
        return self._ncbid


    def getPathSet(self):
        """
            Get the path from the root to the ncbid according to which this scaffold was assigned,
            the path is represented as a set.
        """
        return self._pathSet


    def getToLeafDist(self, ncbid):
        """
            Get the distance from the argument ncbid to the ncbid to which the scaffold was assigned.
        """
        if ncbid in self._ncbidToLeafDist:
            return self._ncbidToLeafDist[ncbid]
        else:
            return None


    def getToPathDist(self, ncbid):
        """
            Get the distance from the argument ncbid to the path according to which this scaffold was assigned.
        """
        if ncbid in self._ncbidToDist:
            return self._ncbidToDist[ncbid]
        else:
            return None


    def getContigsNameList(self):
        """
            Gets all contigs` names of this scaffold as a list.
        """
        return self._contigsNameList


    def getCollectiveLength(self):
        """
            Get the sum of the lengths (in bp) of all constituent contigs.
        """
        length = 0.0
        for ncbid in self._ncbidToBp:
            length += self._ncbidToBp[ncbid]
        return length


    def getConsistencyTotal(self, asCount=False):
        """
            Consistency: def 1.
        """
        consistentCount = 0.0
        for contigName in self._contigsNameList:
            ncbid = self._contigNameToNcbid[contigName]
            if ncbid in self._pathSet:
                consistentCount += 1
        if asCount:
            return float(consistentCount)
        else:
            return float(consistentCount)/len(self._contigsNameList)


    def getConsistencyTotalBp(self, asBpCount=False):
        """
            Consistency: def 2.
        """
        consistentBp = 0.0
        sumBp = 0.0
        for ncbid in self._ncbidToBp:
            bp = self._ncbidToBp[ncbid]
            if ncbid in self._pathSet:
                consistentBp += bp
            sumBp += bp
        if asBpCount:
            return float(consistentBp)
        else:
            return float(consistentBp)/sumBp


    def getConsistencyAvgDist(self, asTotalCount=False):
        """
            Consistency: def3.
        """
        sumDist = 0.0
        for contigName in self._contigsNameList:
            ncbid = self._contigNameToNcbid[contigName]
            sumDist += self._ncbidToDist[ncbid]
        if asTotalCount:
            return float(sumDist)
        else:
            return float(sumDist)/len(self._contigsNameList)


    def getConsistencyWeightedAvgDist(self):
        """
            Consistency: def4.
        """
        return self._weightedDistToPath


    def getConsistencyAvgDistLeaf(self, asTotalCount=False):
        """
            Consistency: def5.
        """
        sumDist = 0.0
        for contigName in self._contigsNameList:
            ncbid = self._contigNameToNcbid[contigName]
            sumDist += self._ncbidToLeafDist[ncbid]
        if asTotalCount:
            return float(sumDist)
        else:
            return float(sumDist)/len(self._contigsNameList)


    def getConsistencyAvgWeightedDistLeaf(self):
        """
            Consistency: def6.
        """
        sumDist = 0.0
        for ncbid in self._ncbidToWeight:
            sumDist += float(self._ncbidToWeight[ncbid])*self._ncbidToLeafDist[ncbid]
        return float(sumDist)


class Consistency():
    """ Main class """

    def __init__(self, fastaFilePath, predFilePath, mappingFilePath, databaseFile,
                 minScaffContigCount=None, minScaffBpLen=None, cladesSet=None, considerContigWithNoScaff=True):
        """
            Initializes the Consistency class.

            @param minScaffContigCount: consider only scaffolds that contain at least this number of contigs
            @param minScaffBpLen: consider only scaffolds with at least this collective length (in bp)
            @param cladesSet: consider only scaffolds that contain at least one contig from this set
            @param considerContigWithNoScaff: consider also contigs that are not assigned to scaffolds (as artificial scaffolds)
        """
        #read data
        self._contigNameToBp = getSequenceToBpDict(fastaFilePath)
        self._contigToPred = predToDict(predFilePath)
        self._scaffToContigsList = getMapping(mappingFilePath, 0, 1, '\t')
        self._taxonomy = _TaxonomyWrapper(databaseFile)

        #for contigName in self._scaffToContigsList['TS28_scaffold10505']:
        #    if contigName == 'TS28_contig188275':
        #        print 'here: ', contigName

        #check the consistency of the data!
        #st = set([])

        #if a contig that is defined in the mappind doesn`t exist (in the fasta file) we remove it
        for scaff in self._scaffToContigsList:
            contigsList = self._scaffToContigsList[scaff]
            #if scaff == 'TS28_scaffold10505':
            #    print 'here.. TS28_scaffold10505'
            #    for c in contigsList:
            #        print c
            removeList = []
            for contig in contigsList:
                #if scaff == 'TS28_scaffold10505':
                #    print contig
                #st.add(contig)
                #if contig == "TS28_contig188275":
                #    print "contig TS28_contig188275 was considered"
                if contig not in self._contigNameToBp:
                    removeList.append(contig)

            for contig in removeList:
                contigsList.remove(contig)
                    #contigsList.remove(contig)
                    #if contig == "TS28_contig188275":
                    #    print "contig TS28_contig188275 was removed"
                #else:
                #    if contig == "TS28_contig188275":
                #        print "contig TS28_contig188275 was NOT removed"

        #print st


        #if a contig was predicted but there is no scaffold assigned to it then this contig is assigned to an "artificial scaffold"
        if considerContigWithNoScaff:
            scaffContigSet = set([])
            for s in self._scaffToContigsList:
                list = self._scaffToContigsList[s]
                for c in list:
                    scaffContigSet.add(c)
            aloneContigSet = set([])
            for c in self._contigToPred:
                if c not in scaffContigSet:
                    aloneContigSet.add(c)

            for c in aloneContigSet:
                scaffName = str('scaffold_' + c) #make up a scaffold name
                assert scaffName not in self._scaffToContigsList, 'The names of contigs are ambiguous!'
                self._scaffToContigsList[scaffName] = [c]

        #filter out scaffolds
        self._scaffolds = dict([])
        for scaffName in self._scaffToContigsList:
            contigsList = self._scaffToContigsList[scaffName]
            if minScaffContigCount != None:
                if len(contigsList) < minScaffContigCount:
                    continue

            if minScaffBpLen != None:
                sum = 0
                for contig in contigsList:
                    sum += self._contigNameToBp[contig]
                if sum < minScaffBpLen:
                    continue

            if cladesSet != None:
                passScaff = False
                for contig in contigsList:
                    if (contig in self._contigToPred) and (self._contigToPred[contig] in cladesSet):
                        passScaff = True
                        break
                if not passScaff:
                    continue

            #process the scaffold, but if everything in the scaffold was assigned to the root, then ignore it!
            s = self._processScaffold(scaffName)
            if s.getNcbid() != 1:
                self._scaffolds[scaffName] = s


    def _getPred(self, contigName):
        """
            Gets prediction for contig "contigName" or 1 if it wasn`t assigned.
        """
        if contigName in self._contigToPred:
            return self._contigToPred[contigName]
        else:
            return int(1)


    def _processScaffold(self, scaffName):
        """
            Gets "Scaffold" - an object that contain pre-processed info.
        """
        allNcbidSet = set([])
        for contigName in self._scaffToContigsList[scaffName]:
            allNcbidSet.add(self._getPred(contigName))
        parentNcbidSet = set([1])
        leafNcbidSet = set([])
        for ncbid in allNcbidSet:
            if ncbid == 1:
                continue
            current = self._taxonomy.getParent(ncbid)
            while current not in parentNcbidSet:
                if current == None:
                    sys.stderr.write('Consistency:Consistency:_processScaffold: current is "None" ' + str(ncbid) + ' '
                                     + str(current) + ' ' + str(parentNcbidSet) + '\n')
                    break
                parentNcbidSet.add(current)
                current = self._taxonomy.getParent(current)
        for ncbid in allNcbidSet:
            if ncbid not in parentNcbidSet:
                leafNcbidSet.add(ncbid)

        ncbidToBp = dict([])
        sumBp = 0
        for contigName in self._scaffToContigsList[scaffName]:

            #if contigName == 'TS28_contig188275':
            #    print contigName, scaffName

            ncbid = self._getPred(contigName)

            bp = self._contigNameToBp[contigName]
            sumBp += int(bp)
            if ncbid not in ncbidToBp:
                ncbidToBp[ncbid] = int(bp)
            else:
                ncbidToBp[ncbid] += int(bp)

        ncbidToWeight = dict([])
        for ncbid in allNcbidSet:
            ncbidToWeight[ncbid] = float(ncbidToBp[ncbid])/float(sumBp)

        #for all paths defined by leaf ncbids compute the weighted distance to all other ncbids
        minDistW = sys.float_info.max
        minDist = sys.maxint
        minDistNcbid = 1
        minPath = set([1])
        for ncbid in leafNcbidSet:
            path = set([1])
            current = ncbid
            while current != 1:
                if current == None:
                    sys.stderr.write('Consistency:Consistency:_processScaffold: current is "None" (while current != 1) ' + str(current) + '\n')
                    break
                path.add(current)
                current = self._taxonomy.getParent(current)

            distW = 0.0
            dist  = 0.0
            for ncbidA in allNcbidSet:
                d = float(self._taxonomy.getDist(ncbidA, path))
                distW += d * ncbidToWeight[ncbidA]
                dist  += d
            if distW < minDistW:
                minDistW = distW
                minDist = dist
                minDistNcbid = ncbid
                minPath = path

        #if everything is assigned to the root, then the distance is 0
        if len(leafNcbidSet) == 0:
            minDistW = 0.0
            minDist = 0.0
            assert minDistNcbid == 1

        #for each ncbid compute the distance to the path
        ncbidToDist = dict([])
        for ncbid in allNcbidSet:
            ncbidToDist[ncbid] = float(self._taxonomy.getDist(ncbid, minPath))

        contigNameToNcbid = dict([])
        for contigName in self._scaffToContigsList[scaffName]:
            contigNameToNcbid[contigName] = self._getPred(contigName)

        #for each ncbid compute the distance to the leaf (minDistNcbid)
        ncbidToLeafDist = dict([])
        for ncbid in allNcbidSet:
            d = int(ncbidToDist[ncbid])
            lcaNcbid = self._taxonomy.getDistantParent(ncbid, d)
            ncbidToLeafDist[ncbid] = float(d + self._taxonomy.getDistTowardsRoot(minDistNcbid, lcaNcbid))

        return ScScaffold(scaffName, minDistNcbid, minPath, minDistW, self._scaffToContigsList[scaffName],
                         contigNameToNcbid, ncbidToBp, ncbidToWeight, ncbidToDist, ncbidToLeafDist)


    def getScaffoldsDict(self):
        """
            Gets all scaffolds (as a dict) that were filtered out according to the parameters.
        """
        return self._scaffolds


    def getScaffoldsPrint(self):
        """
            Gets a list of scaffolds to be printed out.
        """
        buff = ''
        scaffList = []
        for scaffName in self._scaffolds:
            scaffList.append(scaffName)
        scaffList.sort()

        for scaffName in scaffList:
            scaff = self._scaffolds[scaffName]
            contigCount = len(scaff.getContigsNameList())
            pathSet = scaff.getPathSet()
            scaffNcbid = scaff.getNcbid()
            buff += str(scaff.getName() + ', ' + str(scaffNcbid)  + ', ' + str(round(float(scaff.getCollectiveLength())/1000.0, 3)) + 'kbp,  ('
                        + str(int(scaff.getConsistencyTotal(asCount=True))) + '/' + str(contigCount) + ')')
            if abs(scaff.getConsistencyTotal() - 1) > 0.0001:
                buff += str(', ' +  str(round(scaff.getConsistencyTotal()*100,0)) + '%, ' + str(round(scaff.getConsistencyTotalBp()*100,0)) + '%bp')
            if scaff.getConsistencyAvgDist() > 0.0001:
                buff += ',  pathD:, ' + str(round(scaff.getConsistencyAvgDist(),2)) + ', ' + str(round(scaff.getConsistencyWeightedAvgDist(),2)) + 'w'
            if scaff.getConsistencyAvgDistLeaf() > 0.0001:
                buff += str(',  leafD:,' + str(round(scaff.getConsistencyAvgDistLeaf(),2)) + ', ' + str(round(scaff.getConsistencyAvgWeightedDistLeaf(),2)) + 'w')

            buff += ',  ('
            i = 0
            contigList = scaff.getContigsNameList()
            contigList.sort()
            for contig in contigList:
                #contigNcbid = 1
                #if contig in self._contigToPred:
                #    contigNcbid = self._contigToPred[contig]
                contigNcbid = self._getPred(contig)
                bp = 0
                if contig in self._contigNameToBp:
                    bp = self._contigNameToBp[contig]
                buff += str(contig + ' ' + str(int(bp)) + 'bp ' + str(contigNcbid))
                if contigNcbid == scaffNcbid:
                    buff += '*'
                elif contigNcbid in pathSet:
                    buff += str('+' + str(int(self._taxonomy.getDistTowardsRoot(scaffNcbid, contigNcbid))))
                else:
                    buff += str('-' + str(int(scaff.getToLeafDist(contigNcbid))) + '-' + str(int(scaff.getToPathDist(contigNcbid))))
                if (i + 1) == contigCount:
                    buff += ')'
                else:
                    buff += '; '
                i += 1
            buff += '\n'

        return buff


    def getGroupedScaffoldsPrint(self):
        """
            Gets scaffolds grouped according to their ncbid.
        """
        #ncbid -> list of scaffolds
        ncbidToScaffList = dict([])
        for scaffName in self._scaffolds:
            scaff = self._scaffolds[scaffName]
            scaffNcbid = scaff.getNcbid()
            if scaffNcbid not in ncbidToScaffList:
                ncbidToScaffList[scaffNcbid] = [scaff]
            else:
                ncbidToScaffList[scaffNcbid].append(scaff)

        #ncbids
        ncbidSet = set(ncbidToScaffList.keys())

        #ncbid -> scientific name
        scientificNameToNcbid = dict()
        nameList = []
        for ncbid in ncbidSet:
            name = self._taxonomy.getScientificName(ncbid)
            assert name not in scientificNameToNcbid
            scientificNameToNcbid[name] = ncbid
            nameList.append(name)
        nameList.sort()

        buff = ''
        nameList.append('Summary')

        for name in nameList:
            if name != 'Summary':
                ncbid = scientificNameToNcbid[name]
                scaffolds = ncbidToScaffList[ncbid]
            else:
                ncbid = str('all clades (' + str(len(set(ncbidToScaffList.keys()))) + ')')
                scaffolds = []
                for s in self._scaffolds:
                    scaffolds.append(self._scaffolds[s])

            totalContigCount = 0.0 # count the number of all contigs in all scaffolds
            totalBpLen = 0.0 # sum up lengths of all scaffolds` lengths
            totalConsistentContigCount = 0.0 # number of contigs that are consistent
            totalConsistentBpLen = 0.0 # number of Bp that are consistent
            totalPathDist = 0.0 # distance from all contigs to the path
            totalPathDistWeighted = 0.0 # weighted distance from all contigs to the path
            totalLeafDist = 0.0 # distance from all contigs to the respective leafs (of the path)
            totalLeafDistWeighted = 0.0 # weighted distance from all contigs to the respective leaf (of the path)

            for scaff in scaffolds:
                contigCount = float(len(scaff.getContigsNameList()))
                collectiveLength = float(scaff.getCollectiveLength())
                totalContigCount += contigCount #
                totalConsistentContigCount += float(scaff.getConsistencyTotal(asCount=True)) #
                totalBpLen += collectiveLength #
                totalConsistentBpLen += float(scaff.getConsistencyTotalBp(asBpCount=True)) #
                totalPathDist += float(scaff.getConsistencyAvgDist(asTotalCount=True)) #
                totalPathDistWeighted += collectiveLength*float(scaff.getConsistencyWeightedAvgDist()) #
                totalLeafDist += float(scaff.getConsistencyAvgDistLeaf(asTotalCount=True)) #
                totalLeafDistWeighted += collectiveLength*scaff.getConsistencyAvgWeightedDistLeaf()

            buff += str(name + ', (' + str(ncbid) + '), scaffolds: ' + str(len(scaffolds)) + ', contigs: (' + str(int(totalConsistentContigCount))
                        + '/' + str(int(totalContigCount)) + '), '
                        + str(round(((totalConsistentContigCount/totalContigCount)*100.0),2)) + '%, ('
                        + str(round(totalConsistentBpLen/1000.0, 1)) + '/' + str(round(totalBpLen/1000.0, 1)) + ' kb), '
                        + str(round(((totalConsistentBpLen/totalBpLen)*100.0),2)) + '% bp,  pathDist:, '
                        + str(round(totalPathDist/totalContigCount,2)) + ', '
                        + str(round(totalPathDistWeighted/totalBpLen,2)) + 'w, leafDist:, '
                        + str(round(totalLeafDist/totalContigCount,2)) + ', '
                        + str(round(totalLeafDistWeighted/totalBpLen,2)) + 'w'
                       )
            #print totalPathDistWeighted

            buff += '\n'

        return buff


    def close(self):
        self._taxonomy.close()


def test():
    """
        For debugging purposes.
    """
    db = '/Users/ivan/Documents/work/binning/taxonomy/ncbi_taxonomy_20110629/ncbitax_sqlite.db'

    minScaffContigCount = 6
    minScaffBpLen = None
    cladesSet = None #set([133925])

    considerContigWithNoScaff = False

    #SRM
    #fasta = '/Users/ivan/Documents/work/binning/data/Reindeer/contigs.fna'
    #pred = '/Users/ivan/Documents/work/binning/tests/SRM1/07/output/contigs.fna.pOUT'
    #map = '/Users/ivan/Documents/work/binning/data/Reindeer/scaffolds-contigs.tab'

    #HG
    fasta = '/Users/ivan/Documents/work/binning/data/HumanGut/working/contigs.fna'
    #pred = '/Users/ivan/Documents/work/binning/tests/HumanGut/02/output/contigs.fna.pOUT'
    pred = '/Users/ivan/Documents/work/binning/data/HumanGut/working/PPS_contigs_no_minus_one.txt'
    map = '/Users/ivan/Documents/work/binning/data/HumanGut/working/scafftocontig.txt'

    #TW
    #fasta = '/Users/ivan/Documents/work/binning/data/TW/assembly/contigs_tw.fna'
    #pred = '/Volumes/Macintosh HD/Users/ivan/Documents/work/binning/tests/TW/TW13/output/contigs_tw.fna.pOUT'
    #map = '/Users/ivan/Documents/work/binning/data/TW/scaff_contig.tab'

    cons = Consistency(fasta, pred, map, db, minScaffContigCount, minScaffBpLen, cladesSet, considerContigWithNoScaff)

    print cons.getScaffoldsPrint()

    print cons.getGroupedScaffoldsPrint()

    cons.close()


def main():
    """
        Main function.
    """
    parser = argparse.ArgumentParser(description='''Computes the scaffold-contig consistency based on the "maximum support path".''',
                                     epilog='''

                                     ''')

    parser.add_argument('-f', '--fasta', nargs=1, type=file, required=True,
                        help='Fasta file.', metavar='contigs.fna',
                        dest='f')

    parser.add_argument('-p', '--predictions', nargs=1, type=file, required=True,
                        help='Tab separated prediction file (first column sequence name, last column predicted ncbid.', metavar='pred.csv',
                        dest='p')

    parser.add_argument('-m', '--mapping', nargs=1, type=file, required=True,
                        help='Tab separated scaffold-contig mapping file (first column scaffold name, second column contig name.',
                        metavar='scaffoldToContig.csv',
                        dest='m')

    parser.add_argument('-d', '--database', nargs=1, type=file, required=True,
                        help='Database file in the sqlite3 format.', metavar='ncbitax_sqlite.db',
                        dest='d')

    parser.add_argument('-c', '--min-scaff-contig-count', nargs=1,
                        help='Consider only scaffolds that contain at least this number of contigs.', metavar='X',
                        dest='c')

    parser.add_argument('-b', '--min-scaff-bp', nargs=1,
                        help='Consider only scaffolds that has at least this collective length (sum of the lengths of all constituent contigs).',
                        metavar='X',
                        dest='b')

    parser.add_argument('-n', '--allowed-ncbids-list', nargs=1,
                        help='Consider only scaffolds that contain at least one contig that was assigned to one of this clades' +
                        ' (given as comma separated NCBIDs)',
                        metavar='186803,83763,128827',
                        dest='n')

    parser.add_argument('-a', '--consider-artificial-scaff', action='store_true',
                        help='Consider contigs that are not assigned to scaffolds. (For each such contig an artificial scaffold is created.)',
                        dest='a')

    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Print scaffold-contig consistency info for each scaffold.',
                        dest='v')

    args = parser.parse_args()

    #print str(args.c[0])
    #print str(args.b[0])

    if args.c:
        assert len(args.c) == 1
        minScaffContigCount = int(args.c[0])
    else:
        minScaffContigCount = None

    if args.b:
        assert len(args.b) == 1
        minScaffBpLen = int(args.b[0])
    else:
        minScaffBpLen = None

    if args.n:
        cladesSet = set([])
        assert len(args.n) == 1
        for i in str(args.n[0]).split(','):
            cladesSet.add(int(i))
    else:
        cladesSet = None


    assert len(args.f) == 1 and len(args.p) == 1 and len(args.m) == 1 and len(args.d) == 1 #make this nicer!
    #print args.p
    #print args.m
    #print args.d
    #print args.a

    cons = Consistency(args.f[0].name, args.p[0].name, args.m[0].name, args.d[0].name, minScaffContigCount, minScaffBpLen, cladesSet, args.a)

    if args.v:
        print cons.getScaffoldsPrint()

    print cons.getGroupedScaffoldsPrint()

    cons.close()


if __name__ == "__main__":
    if os.getcwd() == '/Users/ivan/Documents/work/python/workspace/pPPS/src':
        test()
    else:
        main()