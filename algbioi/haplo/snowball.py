#!/usr/bin/env python

"""
    Copyright (C) 2015  Ivan Gregor

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Implements the main "Snowball" algorithm of the package
"""
from algbioi.com import fq, common

from algbioi.haplo import join_rec, read_rec


class JoinInfo(object):
    def __init__(self, score, overlap, seedRec, neighborRec):
        """
            Represents an info used to help to join two records.
        """
        # TODO: move to the record class?
        self.score = score
        self.overlap = overlap
        self.seedRec = seedRec
        self.neighborRec = neighborRec


def getBestOverlap(inspectList, tripletMap, maxMismatchQSAllowed, minScore, minAnnotOverlapScore, scoreStopSearch):
    """
        Compute overlaps for all the record pairs and return the best one or None.

        @param inspectList: a list of tuples of records (seed, neighbor)
        @param tripletMap: map: triplet -> aminoacid
        @param maxMismatchQSAllowed: allow mismatches at positions with QS lower or equal to this QS
        @param minScore: minimum score required to join records
        @param minAnnotOverlapScore: minimum annotation overlap score required
        @param scoreStopSearch: stop searching alternative overlaps after this score has been found

        @type inspectList: list[(algbioi.haplo.read_rec.ReadRec, algbioi.haplo.read_rec.ReadRec)]
        @type tripletMap: dict[str,str]
        @type maxMismatchQSAllowed: int
        @type minScore: int
        @type minAnnotOverlapScore: int
        @type scoreStopSearch: int
        @rtype: algbioi.haplo.snowball.JoinInfo

        @return: best matching records overlap or None
    """
    # list of read-record pairs that can be joined
    candidateJoinList = []

    # go over all read-record pairs and inspect its possible overlap scores
    for seedRec, neighborRec in inspectList:

        # find possible overlap positions
        overlapList = getPossibleOverlaps(seedRec, neighborRec)

        # there are possible overlap positions to evaluate
        if overlapList is not None:

            # inspect the overlap list, store possible overlaps
            possibleJoinList = []
            for overlap in overlapList:

                # ge the score of this overlap
                score, annotScore = inspectOverlap(seedRec, neighborRec, overlap, tripletMap, maxMismatchQSAllowed)[:2]

                # the scores are sufficiently good
                if score is not None and score >= minScore and annotScore >= minAnnotOverlapScore:
                    possibleJoinList.append((score, overlap))

                    # stop the search, this score is good enough
                    if score >= scoreStopSearch:
                        break

            # there was at least one sufficient overlap of these two read-records
            if len(possibleJoinList) > 0:
                # get the overlap with the best score
                possibleJoinList.sort(key=lambda x: x[0], reverse=True)
                scoreB, overlapB = possibleJoinList[0]
                candidateJoinList.append(JoinInfo(scoreB, overlapB, seedRec, neighborRec))

    # return the read-record pair with the highest score
    if len(candidateJoinList) > 0:
        candidateJoinList.sort(key=lambda x: x.score, reverse=True)
        return candidateJoinList[0]
    else:
        # no sufficient overlap has been found
        return None


def inspectOverlap(rec1, rec2, overlapIdx, tripletMap, maxMismatchQSAllowed):
    """
        Get a score (number of matching positions) of the records overlap, or None if the overlap contains a mismatch.

        Allowed mismatches at positions with low QS or at positions that result in the same aminoacid, in the case
        it's within the annotated pars of the sequences.
        In addition to the score, annotation score is computed, which is the number of matching positions from within
        the Hmm annotation.

        @param rec1: first read-record
        @param rec2: second read-record
        @param overlapIdx: an overlap to be inspected, start pos of rec1 mapped onto rec2 (it can be a negative value)
        @param tripletMap: map: triplet -> aminoacid
        @param maxMismatchQSAllowed: mismatches at positions with QS lower or equal to this QS  are allowed

        @type rec1: algbioi.haplo.read_rec.ReadRec
        @type rec2: algbioi.haplo.read_rec.ReadRec
        @type overlapIdx: int
        @type tripletMap: dict[str,str]
        @type maxMismatchQSAllowed: int
        @rtype: (int, int, int)

        @return: counts of (mis-)matches (overlap score, annot. score, low-QS-mismatch, aminoacid-match-nucl-mismatch)
    """
    score = 0
    scoreAnnot = 0
    lowQsMismatch = 0
    aminoMatchMismatch = 0

    dna1 = rec1.dnaSeq
    dna2 = rec2.dnaSeq

    # get start positions withing the dna sequences
    if overlapIdx >= 0:
        r1 = 0
        r2 = overlapIdx
    else:
        r1 = - overlapIdx
        r2 = 0

    # length of the dna overlap
    overlapLen = min(len(dna1) - r1, len(dna2) - r2)

    # go along the overlapping dna sequences, count matching positions, stop at the first mismatch
    i = 0
    while i < overlapLen:

        # the position is within the Hmm annotation
        if (rec1.annotStart <= r1 < rec1.annotStart + rec1.annotLen) \
                and (rec2.annotStart <= r2 < rec2.annotStart + rec2.annotLen):
            withinAnnot = True
        else:
            withinAnnot = False

        # position match
        if dna1[r1] == dna2[r2]:
            score += 1
            if withinAnnot:
                scoreAnnot += 1

        # mismatch at low QS
        elif ord(rec1.qsArray[r1]) - 33 <= maxMismatchQSAllowed or ord(rec2.qsArray[r2]) - 33 <= maxMismatchQSAllowed:
            lowQsMismatch += 1

        # the position is within the Hmm annotation (coding region)
        elif withinAnnot:
            # do the corresponding triplets encode the same aminoacid

            # get the frame offset, the triplet start-position, the triplet
            offset = (r1 - rec1.annotStart) % 3
            start = r1 - offset
            triplet1 = dna1[start:start+3]

            offset = (r2 - rec2.annotStart) % 3
            start = r2 - offset
            triplet2 = dna2[start:start+3]
            # TODO: check index out of range !!!
            # both triplets encode the same aminoacid
            if tripletMap[triplet1] == tripletMap[triplet2]:
                aminoMatchMismatch += 1
            else:
                # mismatch found
                return tuple(4 * [None])
        else:
            # mismatch found
            return tuple(4 * [None])

        # move to the next position
        i += 1
        r1 += 1
        r2 += 1

    # no mismatch found
    return score, scoreAnnot, lowQsMismatch, aminoMatchMismatch


def getPossibleOverlaps(rec1, rec2):
    """
        Get possible overlaps of the read-records.

        @type rec1: in_parse.ReadRec
        @type rec1: algbioi.haplo.read_rec.ReadRec

        @return: a list of possible overlaps as index of rec1 start onto rec2 (i.e. negative values possible) or None
        @rtype: list[int]
    """
    # get the overlap length of the Hmm annotated regions (startHmmCoord + local annotation)
    start = max(3 * rec1.hmmCoordStart, 3 * rec2.hmmCoordStart)
    end = min(3 * rec1.hmmCoordStart + rec1.annotLen - 1, 3 * rec2.hmmCoordStart + rec2.annotLen - 1)

    # there is "approximately" an overlap of the Hmm annotated regions
    if start <= end:
        overlapList = []

        # difference of the Hmm coordinates start points
        coordDiff = 3 * (rec1.hmmCoordStart - rec2.hmmCoordStart)

        # pos within rec2 matching start of rec1
        middlePoint = rec2.annotStart + coordDiff - rec1.annotStart

        # how far to go to the left and right from the middle point
        goLeftMax = rec1.annotLen + coordDiff  # may be reduced !!! for performance purposes !!!
        goRightMax = rec2.annotLen - coordDiff  # may be reduced too !!!

        overlapList.append(middlePoint)

        # go left and right and add all possible overlaps
        diff = 3
        goLeft = True
        goRight = True
        while goLeft or goRight:
            if diff < goLeftMax:
                overlapList.append(middlePoint - diff)
            else:
                goLeft = False
            if diff < goRightMax:
                overlapList.append(middlePoint + diff)
            else:
                goRight = False
            diff += 3

        return overlapList

    else:
        # there is no overlap
        return None


# TODO: move parameters to the main !!! 9 15 3 30 94
# TODO: clean the results ? take only sufficient (3) coverage regions ?
# def runSnowball(recList, recSeedList, translTable=11, maxMismatchQSAllowed=9, minScoreReqiured=15,
#                 minAnnotOverlapScore=3, scoreStopSearch=30, maxQS=94):
# def runSnowball(recList, recSeedList, translTable=11, maxMismatchQSAllowed=9, minScoreReqiured=100,
#                 minAnnotOverlapScore=50, scoreStopSearch=200, maxQS=94):
def runSnowball(recList, recSeedList, translTable=11, maxMismatchQSAllowed=9, minScoreReqiured=75,
                minAnnotOverlapScore=40, scoreStopSearch=100, maxQS=94):
    """
        Runs the snowball assembly algorithm. Given read records representing individual reads,
        outputs super-reads (contigs).

        @param recList: read-record list sorted according to the start alignment coordinates within a gene family
        @param recSeedList: read-record list sorted according to the biggest overlap with all the annotations

        @type recList: list[read_rec.ReadRec]
        @type recSeedList: list[read_rec.ReadRec]

        @return: set of read-records representing assembled super-reads (contigs)
        @rtype: set[read_rec.ReadRec]
    """
    if len(recList) == 0 or len(recSeedList) == 0:
        return set()

    # initialize the working set by the seed
    seed = recSeedList[0]
    workingSet = {seed}

    # map: triplet -> aminoacid
    tripletMap = join_rec.getTripletMap(translTable)
    # generator of consensus dna sequences and QS arrays when merging the sequences
    consSeqGen = fq.QsMultMatrix(maxQS)

    # next record to be considered in the snowball
    left = None
    right = None

    # find the position of the seed in the recList
    for i in common.binarySearch(recList, seed, fun=lambda x: x.hmmCoordStart):
        if recList[i].recordId == seed.recordId:
            if i - 1 >= 0:
                left = i - 1
            if i + 1 < len(recList):
                right = i + 1
            break

    # while there is a record to be considered for joining with one of the seeds in the working set
    while (left is not None) or (right is not None):

        # list of pairs (seed, left or right) to be inspected for joining
        inspectList = []

        for seed in workingSet:
            if left is not None:
                inspectList.append((seed, recList[left]))
            if right is not None:
                inspectList.append((seed, recList[right]))

        # get the best overlap, compute overlap of the left and right to all the seeds
        joinInfo = getBestOverlap(inspectList, tripletMap, maxMismatchQSAllowed, minScoreReqiured, minAnnotOverlapScore,
                                  scoreStopSearch)

        moveRight = False
        moveLeft = False
        # there is a sufficient overlap of the seed and a neighbor, join them
        if joinInfo is not None:

            # remove the seed from the working set
            workingSet.remove(joinInfo.seedRec)

            # merge the records (seed and a neighbor)
            mergedRec = read_rec.mergeRecords(joinInfo.seedRec, joinInfo.neighborRec, joinInfo.overlap, consSeqGen,
                                              joinInfo.score)

            # replace the seed in the working set by the merged record
            workingSet.add(mergedRec)

            # was the neighbor left or right
            if left is not None and joinInfo.neighborRec.recordId == recList[left].recordId:
                moveLeft = True
            else:
                assert right is not None and joinInfo.neighborRec.recordId == recList[right].recordId
                moveRight = True
        else:
            # there is no sufficient overlap, add neighbors to the working set
            if left is not None:
                workingSet.add(recList[left])
                moveLeft = True
            if right is not None:
                workingSet.add(recList[right])
                moveRight = True

        # shift the indices to get new neighbors to consider
        if moveRight:
            if right + 1 < len(recList):
                right += 1
            else:
                right = None
        if moveLeft:
            if left - 1 >= 0:
                left -= 1
            else:
                left = None

    return workingSet