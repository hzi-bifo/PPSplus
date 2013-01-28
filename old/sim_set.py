#!/usr/bin/env python

import random

from com.csv import getMapping
from com.csv import getColumnAsList
from com.csv import OutFileBuffer
from com.fasta import fastaFileToDict


class RandLognormLength():
    def __init__(self, mu, sigma, minLen):
        self.mu = mu
        self.sigma = sigma
        self.minLen = minLen

    def getLen(self):
        while True:
            d = random.lognormvariate(self.mu, self.sigma)
            if d >= self.minLen:
                return int(round(d))


def generateSimSet():

    genomesFilePath = '/Users/ivan/Documents/work/binning/data/V35/Mercier_New.fasta'
    genomeList = '/Users/ivan/Documents/work/binning/data/V35/genome_list2.txt'
    outFastaFilePath = '/Users/ivan/Documents/work/binning/data/V35/simSet/simSetTest01.fna'
    outLabelsFilePath = '/Users/ivan/Documents/work/binning/data/V35/simSet/simSetTest01.txt'
    random.seed(95321)
    mu = 7.622665475
    sigma = 0.609863959
    minLength = 1000
    bpPerClade = 1000000

    rand = RandLognormLength(mu, sigma, minLength)
    genomesDict = fastaFileToDict(genomesFilePath)
    ncbidToGenomeName = getMapping(genomeList, 2, 0, sep=';', comment = '#')
    ncbidList = getColumnAsList(genomeList, entryModifyFunction=None, colNum=2, sep=';', comment='#')
    outFna = OutFileBuffer(outFastaFilePath)
    outLabels = OutFileBuffer(outLabelsFilePath)

    for ncbid in ncbidList:
        bpCounter = 0
        count = 0
        genome = genomesDict[ncbidToGenomeName[ncbid][0]]

        while True:
            subLength = rand.getLen()

            start = random.randint(0, (len(genome) - subLength + 1))

            end = start + subLength
            subSeq = genome[start:end]
            name = str(str(ncbid) + '_' + str(count) + '_' + str(len(subSeq)) + '_' + str(start) + '_' + str(end))
            outFna.writeText(str('>' + name + '\n' + subSeq + '\n'))
            outLabels.writeText(str(name + '\t' + str(ncbid) + '\n'))
            count =+ 1
            bpCounter += len(subSeq)
            if bpCounter >= bpPerClade:
                break


    outFna.close()
    outLabels.close()


if __name__ == "__main__":
    generateSimSet()
