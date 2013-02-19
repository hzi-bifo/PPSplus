#!/usr/bin/env python

import os
import re

from com import csv
from com import common
from com import taxonomy_ncbid
from core.taxonomy import Taxonomy


def toRealNames(config, sequences):
    """
        Transforms a PPS file fileName.fas.PP.out that names sequences according to their ids to their real names.
    """
    outIdsPPSFile = str(config.get('inputIdsFastaFile') + '.PP.out')
    outNamesPPSFile = outIdsPPSFile + '.n'
    #os.path.normpath
    print outNamesPPSFile

    try:
        fr = open(os.path.normpath(outIdsPPSFile),'r')
        fw = open(os.path.normpath(outNamesPPSFile),'w')
    except Exception:
        print "Cannot open one of the files:", outIdsPPSFile, "or", outNamesPPSFile
        raise
    else:
        for line in fr:
            if re.match(r'^[0-9]+_[0-9]+[^0-9].*$', line):
                id = re.sub(r'^[0-9]+_([0-9]+)[^0-9].*$',r'\1' , line)
                rest = re.sub(r'^[0-9]+_[0-9]+([^0-9].*)$',r'\1' , line)
                seq = sequences.getSequence(int(id))
                fw.write(seq.name + rest) # seq.scaffold.name
            else:
                fw.write(line)
    finally:
        fr.close()
        fw.close()


def readPPSOutput(sequences, taxonomy, inputFastaIdsPPSFile, overwriteAllPlacements=False):
    """
        Reads the output file of PPS and for each sequence decides:
        if overwriteAllPlacements=True is, then the sequence is placed according to the PPS file regardless of its
        previous placement
        if overwriteAllPlacements=False then if a sequence is placed to a less specific rank, than PPS suggests then
        the sequence is placed according to the PPS file
    """

    infile = str(inputFastaIdsPPSFile + '.out')
    try:
        f = open(os.path.normpath(infile),'r')
    except Exception:
            print "Cannot open file:", infile
            raise
    else:
        #i = 0
        for line in f:
            line = common.noNewLine(line)
            if re.match(r'^[0-9]+_[0-9]+.*[^0-9]+[0-9]+[^0-9]*$', line):
                scaffoldId = int(re.sub(r'^([0-9]+)_[0-9]+.*[^0-9]+[0-9]+[^0-9]*$',r'\1' ,line))
                contigId = int(re.sub(r'^[0-9]+_([0-9]+).*[^0-9]+[0-9]+[^0-9]*$',r'\1' ,line))
                ncbid = int(re.sub(r'^[0-9]+_[0-9]+.*[^0-9]+([0-9]+)[^0-9]*$',r'\1' ,line))
                weight = None # the weight is not yet defined !!!
                if ncbid != 1:
                    #print line, ":", scaffoldId, contigId, ncbid
                    taxPathDictPPS = taxonomy.getPathToRoot(ncbid)
                    if taxPathDictPPS.keys() >= 1:
                        taxPathDictCurrent = sequences.getSequence(contigId).getTaxonomyPath()
                        if taxPathDictCurrent == None:
                            sequences.setTaxonomyPath(contigId, scaffoldId, taxPathDictPPS, weight)#weight = None !!!
                            #i += 1
                        else:
                            if ((overwriteAllPlacements) or (taxPathDictPPS.keys() > taxPathDictCurrent.keys())):
                                sequences.setTaxonomyPathOverride(contigId, scaffoldId, taxPathDictPPS, weight)#weight = None !!!
                                #i += 1
        #print "placed seq by PPS:", i

    finally:
        f.close()


#!!!set the colNum accordingly (may be 0 and 3???
def ppsOut2ppOut(inFile, outFile, taxonomy, taxaRanks):
    outBuff = csv.OutFileBuffer(outFile)
    namesList = csv.getColumnAsList(inFile, entryModifyFunction=None, colNum=0, sep='\t', comment='#')
    ncbidsList = csv.getColumnAsList(inFile, entryModifyFunction=None, colNum=1, sep='\t', comment='#')

    header = str('#Script ppsOut2ppOut\n#ID' + '\t' + 'root')
    for rank in taxaRanks:
        header += str('\t' + rank)
    outBuff.writeText(str(header + '\n'))

    for i in range(len(namesList)):
        name = namesList[i]
        ncbid = ncbidsList[i]
        taxPathDict = taxonomy.getPathToRoot(int(ncbid))
        buff = str(name)
        if taxPathDict == None:
            buff += str('\t')
        else:
            buff += str('\t' + 'root')

        for rank in taxaRanks:
            if (taxPathDict != None) and (rank in taxPathDict) and (not taxPathDict[rank].isCopy()):
                buff += str('\t' + taxPathDict[rank].name)
            else:
                buff += '\t'
        outBuff.writeText(str(buff + '\n'))
    outBuff.close()

class PP2PPSoutParser():

    def __init__(self,taxonomyNcbi, outBuffer):
        self.taxonomy = taxonomyNcbi
        self.out = outBuffer
        self.nameToNcbidDict = dict([])
        self.nameToNcbidDict['Spirochaetes (class)'] = 203691 # was changed to Spirochaetes
        self.nameToNcbidDict['Actinobacteria'] = 201174
        self.nameToNcbidDict['Actinobacteria (class)'] = 1760
        self.nameToNcbidDict['Fusobacteria (class)'] = 203490 #is considered to be: Fusobacteriia

    def parse(self, line):
        tokens = line.split('\t')
        id = str(tokens[0])
        length = len(tokens)
        label = None

        for i in range(len(tokens)-1):
            t = tokens[-(i+1)]
            if t != '':
                if t in self.nameToNcbidDict:
                    label = self.nameToNcbidDict[t]
                else:
                    label = self.taxonomy.getNcbid(t)
                    if label != None:
                        self.nameToNcbidDict[t] = label
                if label == None:
                    print 'CANNOT parse line:', line
                break

        if label != None:
            self.out.writeText(str(id + '\t' + str(label) + '\n'))

    def finalize(self):
        pass

def ppOut2PPSout():
    inFile = '/Users/ivan/Documents/work/binning/data/HumanGut/PP/TS29_scaff.file.0.5.txt'
    outFile = '/Users/ivan/Documents/work/binning/data/HumanGut/PP/TS29_scaff.file.0.5.PPS.txt'
    dbFile = '/Users/ivan/Documents/work/binning/taxonomy/20120828/ncbitax_sqlite.db' #DB
    taxonomy = taxonomy_ncbid.TaxonomyNcbi(dbFile)

    out = csv.OutFileBuffer(outFile)

    csv.forEachLine(inFile, PP2PPSoutParser(taxonomy, out))

    out.close()


def main01():
    #config = Config(open(os.path.normpath('/Users/ivan/Documents/work/binning/tests/CowRumen/01/config.cfg')), 'pPPS')
    #config = Config(open(os.path.normpath('/net/metagenomics/projects/PPSmg/tests/V35/config.cfg')), 'pPPS')
    #configMl = Config2(config, 'MLTreeMap')
    #configPPS = Config2(config, 'PPS')

    #read sequences
    #sequences = Sequences(config)

    #write ids file
    #sequences.writeSequences(config.get('inputIdsFastaFile'))

    #taxonomy = Taxonomy(config.get('databaseFile'), config.get('taxonomicRanks').split(','))

    taxonomicRanks = 'superkingdom,phylum,class,order,family,genus,species'.split(',')
    taxonomy = Taxonomy('/Users/ivan/Documents/work/binning/taxonomy/20120828/ncbitax_sqlite.db', taxonomicRanks)

    #ppsOut2ppOut('D:\\VM\\tmp\\simMC_AMD\\AMD.Arachne.genus', 'D:\\VM\\tmp\\simMC_AMD\\AMD.Arachne.genus.PP.out', taxonomy, config.get('taxonomicRanks').split(','))

    #ppsOut2ppOut('/Users/ivan/Documents/work/binning/data/CowRumen/cowRumenOrderNcbids.txt',
    #             '/Users/ivan/Documents/work/binning/data/CowRumen/cowRumenOrderNcbids.PP.txt', taxonomy, config.get('taxonomicRanks').split(','))

    #ppsOut2ppOut('/net/metagenomics/projects/PPSmg/data/V35/contigsMappedBlast1000LabelsSpecies.txt',
    #             '/net/metagenomics/projects/PPSmg/data/V35/contigsMappedBlast1000LabelsSpecies.PP.txt', taxonomy, config.get('taxonomicRanks').split(','))

    ppsOut2ppOut('/Users/ivan/Documents/work/binning/data/simMC/fromJohannes/contigs.genus.tax',
                 '/Users/ivan/Documents/work/binning/data/simMC/fromJohannes/contigs.genus.PP.tax', taxonomy, taxonomicRanks)


    #readPPSOutput(sequences, taxonomy, config.get('inputIdsFastaFile'))

    #sequences.writePlacements(str(config.get('inputIdsFastaFile') + '.pOUT'), config.get('taxonomicRanks').split(','))

    #toRealNames(config, sequences)
    taxonomy.close()


def collectChildren(taxonomy, ncbid):
    """
        Used in genomesToMask.
    """
    list = taxonomy.childrenNcbids(ncbid)
    if list == None:
        return [ncbid]
    else:
        resultList = []
        for i in list:
            #print 'i', i
            li = collectChildren(taxonomy, i)
            resultList.extend(li)
        return resultList


def genomesToMask():
    rank = 'genus' #which rank will be masked
    fileName = '/Users/ivan/Documents/work/binning/data/simMC/fromJohannes/contigs_genus_ncbids.txt'
    outFile = '/Users/ivan/Documents/work/binning/data/simMC/fromJohannes/genome_genus_masked.txt'
    outFile2 = '/Users/ivan/Documents/work/binning/data/simMC/fromJohannes/genome_ncbids_genus.txt'
    #outFile = '/Users/ivan/Documents/work/binning/data/V35/genome_species_masked.txt' #output file
    #outFile2 = '/Users/ivan/Documents/work/binning/data/V35/genome_ncbids_species.txt' #output file
    #fileName='/Users/ivan/Documents/work/binning/data/V35/genome_ncbids.txt' #list of all genome ncbids
    dbFile = '/Users/ivan/Documents/work/binning/taxonomy/20120828/ncbitax_sqlite.db' #DB
    out = csv.OutFileBuffer(outFile)
    out2 = csv.OutFileBuffer(outFile2)

    genomeNcbids = csv.getColumnAsList(fileName, entryModifyFunction=None, colNum=0, sep=None, comment='#')
    taxonomy = taxonomy_ncbid.TaxonomyNcbi(dbFile)

    maskNcbids = []
    #print len(genomeNcbids), genomeNcbids
    for ncbid in genomeNcbids:
        while taxonomy.getRank(ncbid) != rank:
            ncbid = taxonomy.getParentNcbid(ncbid)
            if int(ncbid) == 1:
                print 'root reached!'
                break
        maskNcbids.append(int(ncbid))

    #print len(Set(maskNcbids)), maskNcbids

    maskSet = set(maskNcbids)
    for i in maskSet:
        out2.writeText(str(str(i) + '\n'))

    resultList = []
    for ncbid in maskSet:
        list = collectChildren(taxonomy, ncbid)
        for i in list:
            out.writeText(str(str(i) + '\n'))
        print ncbid, list

    #print taxonomy.childrenNcbids(818) #997888,818


    out.close()
    out2.close()
    taxonomy.close()



if __name__ == "__main__":
  genomesToMask()
  #ppOut2PPSout()
  #main01()