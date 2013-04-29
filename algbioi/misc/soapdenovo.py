"""
    Contains functionality to process the SOAPdenovo output files.
"""

import sys
import re
from algbioi.com import csv
from algbioi.com import fasta as fas
from algbioi.com import taxonomy_ncbi


class FileHolder():
    def __init__(self, readsF=None, readsR=None, community=None, contig=None,
                 readOnContig=None, newContigIndex=None, scaf=None, scafSeq=None):
        """
            Holds the paths to the SOAPdenovo files.

            @param readsF:
            @param readsR:
            @param community:
            @param contig:
            @param contigL:
            @param taxonomyL:
            @param readOnContig:
            @param newContigIndex:
            @param scaf:
            @param scafSeq:
        """
        self.readsF = readsF
        self.readsR = readsR
        self.community = community
        self.contig = contig
        self.contigL = contig + '1k'
        self.contigM = self.contigL + '.M.fna'
        self.contigMisassembled = self.contigL + '.Misassembled.fna'
        self.contigLTaxL = self.contigL + '.taxL'
        self.contigLTaxW = self.contigL + '.taxW'
        self.contigLTaxM = self.contigL + '.M.tax'
        self.profile = self.contig + '.profile.csv'
        self.readOnContig = readOnContig
        self.newContigIndex = newContigIndex
        self.scaf = scaf
        self.scafSeq = scafSeq


def getLenStat(fileName, minLength=1000):
    """
        Get basic statistics concerning the lengths of the sequence.

        @param fileName: fasta file
        @type fileName: str
    """
    buf = ""
    c = 0
    bp = 0
    minLen = sys.maxint
    maxLen = 0
    totalBp = 0
    totalCount = 0
    for k, l in fas.getSequenceToBpDict(fileName).iteritems():
        totalCount += 1
        totalBp += l
        if l >= minLength:
            c += 1
            bp += l
            if l < minLen:
                minLen = l
            elif l > maxLen:
                maxLen = l

    buf += 'Bigger than %sbp (sequences: %s, Mbp: %s)\n' % (minLength, c, round(float(bp) / 1000000.0, 3))
    buf += 'Bigger than %sbp (min: %s, max %s, avg %s bp)\n' % (minLength, minLen, maxLen, round((float(bp) / c)))
    buf += 'Total (sequences: %s, Mbp: %s\n' % (totalCount, round(float(totalBp) / 1000000.0, 3))
    return buf


def toLongSeq(inFastaFileName, outFastaFileName, minLength=1000):
    """
        Creates a fasta file that contains sequences that are at least minLength long.

        @param inFastaFileName:
        @param outFastaFileName:
        @param minLength:
    """
    out = csv.OutFileBuffer(outFastaFileName)
    for seqId, seq in fas.fastaFileToDictWholeNames(inFastaFileName).iteritems():
        if len(seq) >= minLength:
            out.writeText('>' + str(seqId) + '\n' + str(seq) + '\n')
    out.close()


def getCommunityId(readHeader):
    """
        Gets the communityId given a readHeader.
    """
    return readHeader.strip('>').split('_', 1)[0]


def getCoverage(contigHeader):
    """
        Gets the coverage given a contigHeader.
    """
    return float(contigHeader.split(' ')[-1].split('_')[1])


def getShortContigId(contigHeader):
    """
        Gets the short contig name out of the contig Header.
    """
    return contigHeader.strip('>').split(' ', 1)[0]


def getReadsTaxonIdList(readsFile, communityFile, readHeaderToCommunityId=getCommunityId):
    """
        Gets list of taxonIds in the same order as they are in the readOnContig file.
        The first taxonId is at index 1.

        @param readsFile:
        @param communityFile:
        @param readHeaderToCommunityId:
        @return:
    """
    communityIdToTaxonId = csv.predToDict(communityFile)
    d = [None]
    rowList = csv.getColumnAsList(readsFile, colNum=0, sep='\n')
    for line in rowList:
        if str(line).startswith('>'):
            taxonId = int(communityIdToTaxonId.get(readHeaderToCommunityId(line)))
            d.append(taxonId)
            d.append(taxonId)
    return d


def toContigsLabelList(inFastaFileName, readsF, readsR, readOnContig, community, outMappingFileName):
    """
        Gets mapping from contigIds to lists of taxonIds of individual reads of the contigs.

        @param inFastaFileName:
        @param readsF:
        @param readsR:
        @param readOnContig:
        @param community:
        @param outMappingFileName:
    """
    # contigIds
    contigIdToBp = fas.getSequenceToBpDict(inFastaFileName)

    # map: contigId -> list of readIds
    contigIdToReadList = csv.getMapping(readOnContig, 1, 0, sep='\t', comment='r')

    # taxonIds as a list for reads
    readFTaxonIdList = getReadsTaxonIdList(readsF, community)
    print 's1'
    readRTaxonIdList = getReadsTaxonIdList(readsR, community)
    print 's2'

    if len(readFTaxonIdList) != len(readRTaxonIdList):
        print('toContigsLabels: different number of reads in the reads files, exit')
        return

    for i in range(len(readFTaxonIdList))[1:]:
        if readFTaxonIdList[i] != readRTaxonIdList[i]:
            print('toContigsLabels: at index %s different taxon ids %s and %s' %
                  (i, readFTaxonIdList[i], readRTaxonIdList[i] ))
        if readFTaxonIdList[i] is None or readRTaxonIdList[i] is None:
            print('toContigsLabels: at index %s, one is None %s or %s' % (i, readFTaxonIdList[i], readRTaxonIdList[i]))
    print 's3'
    #
    out = csv.OutFileBuffer(outMappingFileName)
    for contigId in contigIdToBp:
        readList = contigIdToReadList[contigId]
        taxonIdList = []
        for readId in readList:
            taxonIdList.append(readFTaxonIdList[int(readId)])
        out.writeText(str(contigId) + '\t' + ','.join(map(str, taxonIdList)) + '\n')
    out.close()
    print 's4'


def toContigsLabels(inMapFile, outMapFile):
    """
        Creates the label of contigs from the label of reads.

        @param inMapFile: maps contigId to a list of read taxonIds
        @param outMapFile: maps contigId to weight and the most prevalent taxonId
    """
    out = csv.OutFileBuffer(outMapFile)

    for line in csv.getColumnAsList(inMapFile, sep='\n'):
        contigId, taxonIds = str(line).split('\t')
        taxonIdsList = map(int, str(taxonIds).split(','))
        idToCount = {}
        totalCount = 0.0
        for taxonId in taxonIdsList:
            totalCount += 1
            if taxonId in idToCount:
                idToCount[taxonId] += 1
            else:
                idToCount[taxonId] = 1
        pairList = []
        for taxonId, count in idToCount.iteritems():
            pairList.append((taxonId, count))
        pairList.sort(key=lambda x: x[1], reverse=True)
        weight = round(float(pairList[0][1]) / totalCount, 3)
        out.writeText(str(contigId) + '\t' + str(weight) + '\t' + str(pairList[0][0]) + '\n')

    out.close()


def toWellMappedContigs(inFastaFile, inTaxonomyWFile,
                        outFastaFile, outFastaMisAssembledFile, outTaxonomyFile, weightThreshold=0.99):
    """
        Creates the fasta and mapping files that contain well assembled contigs (filter out misassembled contigs).

        @param inFastaFile: input fasta file with contigs
        @param inTaxonomyWFile: input file that contains taxonomy with weights (seqId, weight, taxonId)
        @param outFastaFile: fasta file containing well assembled sequences
        @param outFastaMisAssembledFile: fasta file containing misassembled contigs
        @param outTaxonomyFile: resulting taxonomy of the well assembled sequences (seqId, taxonId)
        @param weightThreshold: only contigs the weight of which is at least this value will be taken
        @return: statistics
    """
    seqIdToTaxonId = csv.predToDict(inTaxonomyWFile)
    seqIdToWeight = csv.getMapping(inTaxonomyWFile, 0, 1, '\t')
    outFastaOk = csv.OutFileBuffer(outFastaFile)
    outFastaMis = csv.OutFileBuffer(outFastaMisAssembledFile)
    outTaxonomyOk = csv.OutFileBuffer(outTaxonomyFile)

    totalBp = 0.0
    totalCount = 0.0
    okBp = 0.0
    okCount = 0.0
    avgSumBp = 0.0

    for seqId, seq in fas.fastaFileToDictWholeNames(inFastaFile).iteritems():
        bp = len(seq)
        totalBp += bp
        totalCount += 1
        seqIdPrefix = str(seqId).split(' ')[0]
        weight = seqIdToWeight[seqIdPrefix][0]
        fastaEntry = '>' + str(seqIdPrefix) + '\n' + str(seq) + '\n'
        if float(weight) >= weightThreshold:
            outFastaOk.writeText(fastaEntry)
            outTaxonomyOk.writeText(str(seqIdPrefix) + '\t' + str(seqIdToTaxonId[seqIdPrefix]) + '\n')
            okBp += bp
            okCount += 1
            avgSumBp += getCoverage(seqId) * bp
        else:
            outFastaMis.writeText(fastaEntry)

    outFastaOk.close()
    outFastaMis.close()
    outTaxonomyOk.close()

    return 'Taken: %s/%sMB, %s/%sseq, %s%% bp %s%% seq, avg coverage %s' % (round(okBp / 1000000, 2),
                                                                            round(totalBp / 1000000, 2),
                                                                            okCount, totalCount,
                                                                            round((okBp / totalBp) * 100, 2),
                                                                            round((okCount / totalCount) * 100, 2),
                                                                            round(avgSumBp / okBp, 3))


def getProfile(readsFFastaFile,
               communityFile, contigMFastaFile, contigLFastaFile, taxonomyMFile, taxonomyDbFile, outProfileFile):
    """
        Gets the profile of the dataset.

        @param readsFFastaFile:
        @param communityFile:
        @param contigMFastaFile:
        @param contigLFastaFile:
        @param taxonomyMFile:
        @param taxonomyDbFile: taxonomy in the sqlite3 format
        @param outProfileFile: output file
    """
    # get map: taxonId -> read count
    taxonIdToReadCount = {}
    readTotalCount = 0
    for taxonId in getReadsTaxonIdList(readsFFastaFile, communityFile, readHeaderToCommunityId=getCommunityId)[1:]:
        if taxonId in taxonIdToReadCount:
            taxonIdToReadCount[taxonId] += 1
        else:
            taxonIdToReadCount[taxonId] = 1
        readTotalCount += 1

    # get map: taxonId -> contig count
    # get map: taxonId -> contig bp
    taxonIdToContigCount = {}
    taxonIdToContigBp = {}
    totalContigCount = 0
    seqIdToTaxonId = csv.predToDict(taxonomyMFile)
    seqIdToBp = fas.getSequenceToBpDict(contigMFastaFile)
    for seqId, bp in seqIdToBp.iteritems():
        totalContigCount += 1
        taxonId = seqIdToTaxonId[seqId]
        if taxonId in taxonIdToContigBp:
            taxonIdToContigBp[taxonId] += bp
        else:
            taxonIdToContigBp[taxonId] = bp
        if taxonId in taxonIdToContigCount:
            taxonIdToContigCount[taxonId] += 1
        else:
            taxonIdToContigCount[taxonId] = 1

    taxonIdToTotalBp = {}
    taxonIdToAvgSumCov = {}
    taxonIdToAvgCov = {}
    totalBp = 0.0
    for taxonId in taxonIdToContigBp:
        taxonIdToTotalBp[taxonId] = 0.0
        taxonIdToAvgSumCov[taxonId] = 0.0
        taxonIdToAvgCov[taxonId] = 0.0

    for seqId in fas.fastaFileToDictWholeNames(contigLFastaFile):
        shortSeqId = getShortContigId(seqId)
        if shortSeqId in seqIdToBp:
            coverage = getCoverage(seqId)
            bp = seqIdToBp[shortSeqId]
            taxonId = seqIdToTaxonId[shortSeqId]
            taxonIdToTotalBp[taxonId] += bp
            taxonIdToAvgSumCov[taxonId] += float(coverage) * float(bp)
            totalBp += bp

    for taxonId, bp in taxonIdToTotalBp.iteritems():
        if bp > 0:
            taxonIdToAvgCov[taxonId] = taxonIdToAvgSumCov[taxonId] / float(bp)

    tupleList = []
    taxonomy = taxonomy_ncbi.TaxonomyNcbi(taxonomyDbFile, considerNoRank=True)
    avgCoverage = 0.0
    for taxonId, readCount in taxonIdToReadCount.iteritems():
        tupleList.append((taxonId,
                          round(100 * (readCount / float(readTotalCount)), 1),
                          round(100 * (taxonIdToTotalBp.get(taxonId, 0) / float(totalBp)), 1),
                          round(taxonIdToAvgCov.get(taxonId, 0), 2),
                          round(taxonIdToTotalBp.get(taxonId, 0) / 1000000.0, 2),
                          taxonIdToContigCount.get(taxonId, 0),
                          taxonomy.getScientificName(taxonId)
        ))
        avgCoverage += taxonIdToAvgCov.get(taxonId, 0) * taxonIdToTotalBp.get(taxonId, 0)
    avgCoverage = avgCoverage / float(totalBp)
    tupleList.sort(key=lambda x: x[1], reverse=True)

    out = csv.OutFileBuffer(outProfileFile)
    out.writeText('#taxonId, % reads, % contigs, avg coverage, MB contigs, contigs count, scientific name\n')
    for entry in tupleList:
        out.writeText(','.join(map(str, entry)) + '\n')



    out.writeText('#Sum/Avg., -, -, ' + str(round(avgCoverage, 2)) + ', ' + str(round(totalBp / 1000000.0, 2)) +
                  ', ' + str(totalContigCount) + ', -\n')
    out.close()
    taxonomy.close()


def _main():
    taxonomyDbFile = '/Users/ivan/Documents/work/binning/taxonomy/20121122/ncbitax_sqlite.db'

    fileHolderLn = FileHolder(readsF='/Users/ivan/Documents/nobackup/assembly/lognorm/ReadsF_lognormal.fasta',
                            readsR='/Users/ivan/Documents/nobackup/assembly/lognorm/ReadsR_lognormal.fasta',
                            community='/Users/ivan/Documents/work/binning/data/mercier51Strains/syn-mercier51strains/generation/community_20121116.tax',
                            contig='/Users/ivan/Documents/nobackup/assembly/lognorm/soap_lognorm.contig',
                            readOnContig='/Users/ivan/Documents/nobackup/assembly/lognorm/soap_lognorm.readOnContig',
                            newContigIndex='/Users/ivan/Documents/nobackup/assembly/lognorm/soap_lognorm.newContigIndex',
                            scaf='/Users/ivan/Documents/nobackup/assembly/lognorm/soap_lognorm.scaf',
                            scafSeq='/Users/ivan/Documents/nobackup/assembly/lognorm/soap_lognorm.scafSeq')

    fileHolderUn = FileHolder(readsF='/Users/ivan/Documents/nobackup/assembly/uniform/ReadsF_uniform.fasta',
                        readsR='/Users/ivan/Documents/nobackup/assembly/uniform/ReadsR_uniform.fasta',
                        community='/Users/ivan/Documents/work/binning/data/mercier51Strains/syn-mercier51strains/generation/community_20121116.tax',
                        contig='/Users/ivan/Documents/nobackup/assembly/uniform/soap_uniform.contig',
                        readOnContig='/Users/ivan/Documents/nobackup/assembly/uniform/soap_uniform.readOnContig',
                        newContigIndex='/Users/ivan/Documents/nobackup/assembly/uniform/soap_uniform.newContigIndex',
                        scaf='/Users/ivan/Documents/nobackup/assembly/uniform/soap_uniform.scaf',
                        scafSeq='/Users/ivan/Documents/nobackup/assembly/uniform/soap_uniform.scafSeq')

    # choose input data
    fh = fileHolderLn
    #fh = fileHolderUn

    # lengths statistics
    if False:
        print('contigs:')
        print getLenStat(fh.contig)
        print('scaffolds:')
        print getLenStat(fh.scaf)

    # take contigs of appropriate length
    if False:
        toLongSeq(fh.contig, fh.contigL, minLength=1000)

    # get labels (map: contigId -> list of read taxonIds)
    if False:
        toContigsLabelList(fh.contigL, fh.readsF, fh.readsR, fh.readOnContig, fh.community, fh.contigLTaxL)

    # get labels (map: contigId -> weight, majority taxonId)
    if False:
        toContigsLabels(fh.contigLTaxL, fh.contigLTaxW)

    if False:
        print toWellMappedContigs(fh.contigL, fh.contigLTaxW, fh.contigM, fh.contigMisassembled, fh.contigLTaxM,
                                  weightThreshold=1.0)

    if True:
        getProfile(fh.readsF, fh.community, fh.contigM, fh.contigL, fh.contigLTaxM, taxonomyDbFile, fh.profile)


def filterOutReads():
    inFasta = ''
    outFasta = ''
    out = csv.OutFileBuffer(outFasta)
    notAllowedSet = set(['BA000019.2'])  # Nostoc sp. PCC 7120
    for seqId, seq in fas.fastaFileToDict(inFasta).iteritems():
        id = re.sub(r'([^_]+)_.*', r'\1', seqId)
        if id not in notAllowedSet:
            out.writeText('>' + str(seqId) + '\n' + str(seq) + '\n')
    out.close()


def filterOutContigs():
    inFasta = ''
    outFasta = ''
    predFile = ''
    out = csv.OutFileBuffer(outFasta)
    seqIdToTaxonId = csv.predToDict(predFile)
    notAllowedTaxonIdSet = set([103690])  # Nostoc sp. PCC 7120
    for seqId, seq in fas.fastaFileToDict(inFasta).iteritems():
        if int(seqIdToTaxonId[seqId]) not in notAllowedTaxonIdSet:
            out.writeText('>' + str(seqId) + '\n' + str(seq) + '\n')
    out.close()


def tmpCmp():
    communityList = csv.getColumnAsList('/Users/ivan/Documents/work/binning/data/mercier51Strains/syn-mercier51strains/'
                                        'generation/community_20121116.tax',
                        colNum=1, sep='\t')
    profileList = csv.getColumnAsList('/Users/ivan/Documents/nobackup/assembly/uniform/soap_uniform.contig.profile.csv',
                                      colNum=0, sep=',')

    cSet = set(map(int, communityList))
    pSet = set(map(int, profileList))
    for i in cSet:
        if i not in pSet:
            print("Ncbid %s from community is not in profile" % i)
    for i in pSet:
        if i not in cSet:
            print("Ncbid %s from profile is not in community" % i)


if __name__ == "__main__":
    #tmpCmp()
    _main()


    # 266265, Burkholderia xenovorans LB400, Bacteria; Proteobacteria; Betaproteobacteria; Burkholderiales; Burkholderiaceae; Burkholderia; Burkholderia xenovorans
    # 70601, Pyrococcus horikoshii OT3, Archaea; Euryarchaeota; Thermococci; Thermococcales; Thermococcaceae; Pyrococcus; Pyrococcus horikoshii
    # 269796, Rhodospirillum rubrum ATCC 11170, Bacteria; Proteobacteria; Alphaproteobacteria; Rhodospirillales; Rhodospirillaceae; Rhodospirillum; Rhodospirillum rubrum
    # 525146, Desulfovibrio desulfuricans subsp. desulfuricans str. ATCC 27774, Bacteria; Proteobacteria; delta/epsilon subdivisions; Deltaproteobacteria; Desulfovibrionales; Desulfovibrionaceae; Desulfovibrio; Desulfovibrio desulfuricans; Desulfovibrio desulfuricans subsp. desulfuricans
