#!/usr/bin/env python

#/AM/home-0/shared/python/Python-2.7.1/python

import re
import os
import Common
from Config import Config
from Taxonomy import Taxonomy

#takes scaffold-contigs mapping and scaffold placement, outputs contigs placement
#
#@param scafContigFile tab sepparated scaffold-contigs mapping (scaffold \t contig)
#@param scafPPSOutFile scaffold predictions (PPS output file)
#@param contigPPSOutFile contigs predictions (as if it was a PPS output file)
def scafToContigOutput(scafContigFile, scafPPSOutFile, contigPPSOutFile):
    scafToContigs = dict([])
    try:
        f = open(os.path.normpath(scafContigFile),'r')
    except Exception:
        print "Cannot open file:", scafContigFile
        raise
    else:
        for line in f:
            line = Common.noNewLine(line)
            scaffold = re.sub(r'^[ ]*([^ \t]+)\t[^ \t]*',r'\1', line)
            contig = re.sub(r'^[ ]*[^ \t]+\t([^ \t]*)',r'\1', line)
            if scaffold in scafToContigs:
                scafToContigs[scaffold].append(contig)
            else:
                temp = []
                temp.append(contig)
                scafToContigs[scaffold] = temp

    try:
        fr = open(os.path.normpath(scafPPSOutFile),'r')
        fw = open(os.path.normpath(contigPPSOutFile),'w')
    except Exception:
        print "Cannot open one of the files:", scafPPSOutFile, contigPPSOutFile
        raise
    else:
        for line in fr:
            line = Common.noNewLine(line)
            if len(line) == 0 or re.match('#', line):
                fw.write(line + '\n')
            else:
                scaffold = re.sub(r'^[ ]*([^ \t]+)[ \t]*.*$',r'\1', line)
                assignment = re.sub(r'^[ ]*[^ \t]+[ \t]*(.*$)',r'\1', line)
                if scaffold in scafToContigs:
                    contigsList = scafToContigs[scaffold]
                    for contig in contigsList:
                        lineW = str(contig + '\t' + assignment + '\n')
                        fw.write(lineW)
                    #print str(lineW),
                else:
                    print 'there is not scaffold-contigs mapping for scaffold:', scaffold

#transforms the PPS out file to a compatible PPS PP.out file
def ppsOutToPPOut(ppsOutFile, outPPOutFile, taxaRanks, taxonomy):
    print ppsOutFile

    #contig file to an ncbid
    contigToNcbid = dict([])
    try:
        f = open(os.path.normpath(ppsOutFile),'r')
    except Exception:
        print "Cannot open file:", ppsOutFile
        raise
    else:
        for line in f:
            line = Common.noNewLine(line)
            contig = re.sub(r'^[ ]*([^\t]+)\t.*$',r'\1', line)
            try:
                ncbid = int(re.sub(r'^.*\t([0-9]+)[ \t]*$',r'\1', line))
            except Exception:
                print 'line skipped:', line
                continue
            contigToNcbid[contig] = ncbid
            #print str('|' + contig + '|' + str(ncbid) + '|')

    try:
        f = open(os.path.normpath(outPPOutFile), 'w')
        f.write('#Translate output to PP.out format from: ' + ppsOutFile + '\n#\n'),
        header = str('#ID' + '\t' + 'root')
        for rank in taxaRanks:
            header += str('\t' + rank)
        f.write(header)

        for contig in contigToNcbid:
            taxPathDict = taxonomy.getPathToRoot(contigToNcbid[contig])
            entry = str('\n' + contig)
            if taxPathDict == None:
                entry += str('\t')
            else:
                entry += str('\t' + 'root')
            for rank in taxaRanks:
                if (taxPathDict != None) and (rank in taxPathDict) and (not taxPathDict[rank].isCopy()):
                    entry += str('\t' + taxPathDict[rank].name)
                else:
                    entry += '\t'
            f.write(entry)
    except Exception:
        print "Cannot create a file or write to it:", outFile
        raise
    finally:
        f.close()


    #open file for write

    def writePlacementsPPOut(self, outFile, taxaRanks, outputFileContigSubPattern):

        try:
            f = open(os.path.normpath(outFile), 'w')

            f.write('#Output of pPPS\n#\n'),
            header = str('#ID' + '\t' + 'root')
            for rank in taxaRanks:
                header += str('\t' + rank)
            f.write(header)

            for seq in self.sequences:
                entry = str('\n' + re.sub(outputFileContigSubPattern, r'\1' , seq.name))
                taxPathDict = seq.getTaxonomyPath()
                if taxPathDict == None:
                    entry += str('\t')
                else:
                    entry += str('\t' + 'root')
                for rank in taxaRanks:
                    if (taxPathDict != None) and (rank in taxPathDict) and (not taxPathDict[rank].isCopy()):
                        entry += str('\t' + taxPathDict[rank].name)
                    else:
                        entry += '\t'
                f.write(entry)
        except Exception:
            print "Cannot create a file or write to it:", outFile
            raise
        finally:
            f.close()

if __name__ == "__main__":
    #test 2
    #ppsOutFile = 'D:\A_Phylo\A_Metagenomic\data\humanGut\PPS_contigs.txt'
    #outPPOutFile = 'D:\A_Phylo\A_Metagenomic\data\humanGut\PPS_PP_contigs.txt'
    #ppsOutFile = 'C:/Documents and Settings/Administrator/Desktop/temp/johdroPred/inputTW.fas.ids04.lP'
    ppsOutFile = 'C:/Documents and Settings/Administrator/Desktop/temp/johdroPred/inputTW.fas.ids05.lP'
    #outPPOutFile = 'C:/Documents and Settings/Administrator/Desktop/temp/johdroPred/inputTW.fas.ids04.lP.PP.out'
    outPPOutFile = 'C:/Documents and Settings/Administrator/Desktop/temp/johdroPred/inputTW.fas.ids05.lP.PP.out'


    config = Config(open(os.path.normpath('D:\\A_Phylo\\A_Metagenomic\\pPPS\\workspace\\pPPS\\config01.cfg')), 'pPPS')
    databaseFile = os.path.normpath(config.get('databaseFile'))
    taxonomicRanks = config.get('taxonomicRanks').split(',')
    taxonomy= Taxonomy(databaseFile, taxonomicRanks)
    ppsOutToPPOut(ppsOutFile, outPPOutFile, taxonomicRanks, taxonomy)

    #test 1
    #scafContigFile = 'D:/A_Phylo/A_Metagenomic/reindeer/data/scaffolds-contigs.tab'
    #scafPPSOutFile = 'D:/A_Phylo/A_Metagenomic/reindeer/predictions/pps04/scaffoldsOut/SRM_Scaffolds_namesOnly.fna.PP.out'
    #contigPPSOutFile = 'D:/A_Phylo/A_Metagenomic/reindeer/predictions/pps04/scaffoldsOut/SRM_Scaffolds_namesOnly.fna.PP.out_contigs'
    #scafToContigOutput(scafContigFile, scafPPSOutFile, contigPPSOutFile)