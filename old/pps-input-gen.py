#!/usr/bin/env python

import os
import sys
import shutil
import re
import argparse
import sqlite3
import glob


def main():
    """
        #-s summary.txt -p predictions.txt -f tw.fna -e expertSSDDir -t ncbitax_sqlite.db -o ssdDir -g 3 -b 100000
    """
    parser = argparse.ArgumentParser(description='''Generate the input data for PPS''',
                                     epilog='''clades and the sample specific data''')

    parser.add_argument('-s', '--pred-summary', action='store', nargs=1, required=True,
                        help='Summary file that summarizes predictions by the 16S analysis', metavar='summary.txt',
                        dest='summary')

    parser.add_argument('-p', '--predictions', action='store', nargs=1, required=True,
                        help='Prediction file from the 16S analysis',
                        dest='predictions')

    parser.add_argument('-f', '--fasta-file', action='store', nargs=1,
                        help='Sample fasta file.',
                        dest='fastaFile')

    parser.add_argument('-e', '--expert-data-dir', action='store', nargs=1,
                        help='Directory that contains expert sample specific data.',
                        dest='expertDataDir')

    parser.add_argument('-t', '--taxonomy', action='store', nargs=1, required=True,
                        help='NCBI taxonomy in sqlite3 format', metavar='ncbitax_sqlite.db',
                        dest='taxonomy')

    parser.add_argument('-o', '--output-dir', action='store', nargs=1, required=True,
                        help='Output directory', metavar='outputDir',
                        dest='outputDir')

    parser.add_argument('-g', '--min-gen-count', action='store', nargs=1,
                        help='Minimum number of genomes from distinct species to model a clade', metavar='3',
                        dest='minGenCount')

    parser.add_argument('-b', '--min-bp-count', action='store', nargs=1,
                        help='Minimum number of base pairs to model a clade', metavar='100000',
                        dest='minBpCount')


    args = parser.parse_args()

    if not args.summary:
        sys.stderr.write('Enter a summary file.')
        sys.stderr.write(parser.print_help())
    else:
        summaryFile = args.summary[0]

    if not args.predictions:
        sys.stderr.write('Enter correct prediction file.')
        sys.stderr.write(parser.print_help())
    else:
        predictionsFile = args.predictions[0]

    if not args.fastaFile:
        sys.stderr.write('Enter correct fasta file.')
        sys.stderr.write(parser.print_help())
    else:
        fastaFile = args.fastaFile[0]

    if not args.expertDataDir:
        sys.stderr.write('Enter correct directory with expert sample specific data.')
        sys.stderr.write(parser.print_help())
    else:
        expertDataDir = args.expertDataDir[0]

    if not args.taxonomy:
        sys.stderr.write('Enter correct taxonomy')
        sys.stderr.write(parser.print_help())
    else:
        databaseFile = args.taxonomy[0]

    if not args.outputDir:
        sys.stderr.write('Enter the output directory')
        sys.stderr.write(parser.print_help())
    else:
        outputDir = args.outputDir[0]

    if not args.minGenCount:
        minGenCount = 3
    else:
        minGenCount = int(args.minGenCount[0])

    if not args.minBpCount:
        minBpCount = 100000
    else:
        minBpCount = int(args.minBpCount[0])

    #read the directory with the expert sample specific data and get the ncbids as a "expertList"

    expertSet = set([])
    expertFiles = []
    path = os.path.normpath(expertDataDir)
    listing = os.listdir(path)
    for file in listing:
        try:
            ncbid = int(re.sub(r'^([0-9]+)\..*',r'\1', file))
        except Exception:
            continue
        expertSet.add(ncbid)
        expertFiles.append(file)

    #read the summary file, for each line, include a clade if (or):
    sList = []
    try:
        f = open(os.path.normpath(summaryFile),'r')
    except Exception:
        print "Cannot open file:", summaryFile
        raise
    else:
        i = 0
        for line in f:
            if i == 0:
                i += 1
                continue

            line = noNewLine(line)

            bpCount = int(re.sub(r'^([^\t]+)\t.*',r'\1', line))
            #seqCount =
            genomeCount = int(re.sub(r'^[^\t]+\t[^\t]+\t([^\t]+).*',r'\1', line))
            #lineage
            ncbid = int(re.sub(r'^[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t([^\t]+).*',r'\1', line))

            if (bpCount >= minBpCount) or (genomeCount >= minGenCount):
                sList.append(ncbid)

    finally:
        f.close()

    #parent clades of all ncbids in the expert list or in the 16S list
    parentSet = set([])

    for ncbid in expertSet:
        parentSet = parentSet.union(ncbidToParentNcbids(ncbid, databaseFile))

    for ncbid in sList:
        parentSet = parentSet.union(ncbidToParentNcbids(ncbid, databaseFile))

    cladesList = []
    cladesListS = []
    for clade in expertSet:
        cladesList.append(clade)
    for clade in sList:
        if clade not in parentSet:
            incl = True
            clPathToRoot = ncbidToParentNcbids(clade, databaseFile)
            for ncbid in expertSet:
                if ncbid in clPathToRoot:
                    incl = False
                    break
            if incl:
                cladesList.append(clade)
                cladesListS.append(clade)

    #store the clades to model to a file!!!
    i=0
    buffer = ''
    for c in cladesList:
        if i == 0:
            buffer += str(c)
            i += 1
        else:
            buffer += str('\n' + str(c))
    writeToFile(buffer, outputDir, 'ncbids.txt')

    #create sample specific fasta files

    #create destination directory
    dirPath = os.path.join(os.path.normpath(outputDir),'ssd')
    try:
        os.mkdir(dirPath)
    except OSError:
        for filePath in glob.glob(os.path.join(os.path.normpath(dirPath),r'*.fas')):
            os.remove(filePath)
        for filePath in glob.glob(os.path.join(os.path.normpath(dirPath),r'*.fna')):
            os.remove(filePath)


    #copy expert data to the destination directory
    for eFile in expertFiles:
        pathSrc = os.path.join(os.path.normpath(expertDataDir), eFile)
        pathDst = os.path.join(os.path.normpath(dirPath), eFile)
        shutil.copy (pathSrc, pathDst)

    #create sample specific data out of the 16S DATA

    #get map: ncbid -> set of seqIds
    ncbidToSeqIdList = dict([])
    cladesSetS = set([])
    for c in cladesListS:
        cladesSetS.add(c)

    #read the prediction file
    #print cladesListS
    ncbidToSetOfSeqIds = readPredictions(predictionsFile, cladesSetS)

    #for ncbid in ncbidToSetOfSeqIds:
    #    print ncbid
    #    print ncbidToSetOfSeqIds[ncbid]

    #prepare buffer: ncbid -> sequences in fasta file
    buffer = dict([])
    for ncbid in ncbidToSetOfSeqIds:
        buffer[ncbid] = ''

    try:
        f = open(os.path.normpath(fastaFile),'r')
    except Exception:
        print "Cannot open file:", fastaFile
        raise
    else:
        name = ''
        seq = ''
        for line in f:
            line = noNewLine(line)
            if re.match('>', line):
                if seq != '':
                    assert name != ''
                    ###
                    for ncbid in ncbidToSetOfSeqIds:
                        if name in ncbidToSetOfSeqIds[ncbid]:
                            buffer[ncbid] += str('>' + name + '\n' + seq + '\n')
                    seq = ''
                name = line.replace('>','')

            else:
                seq += line
        if seq != '':
            assert name != ''
            ###
            for ncbid in ncbidToSetOfSeqIds:
                if name in ncbidToSetOfSeqIds[ncbid]:
                    buffer[ncbid] += str('>' + name + '\n' + seq + '\n')
    finally:
        f.close()

    #for ncbid in buffer:
    #    print str(ncbid) + '---------------------------------'
    #    print buffer[ncbid]

    #store buffers to files
    for ncbid in buffer:
        writeToFile(buffer[ncbid], os.path.normpath(dirPath), str(str(ncbid) + '.1.fna'))



#ncbidToSetOfSeqIds
def readPredictions(predFile, cladesSetS):
    ncbidToSetOfSeqIds = dict([])
    try:
        f = open(os.path.normpath(predFile),'r')
    except Exception:
        print "Cannot open file:", predFile
        raise
    else:
        for line in f:
            line = noNewLine(line)
            seqId = re.sub(r'^([^\t]+)\t.*',r'\1', line)
            ncbids = re.sub(r'^[^\t]+\t(.*)',r'\1', line)
            ncbid = lastNcbid(ncbids)
            if ncbid not in cladesSetS:
                continue
            if ncbid not in ncbidToSetOfSeqIds:
                ncbidToSetOfSeqIds[ncbid] = set([])
            ncbidToSetOfSeqIds[ncbid].add(seqId)
    finally:
        f.close()
    return ncbidToSetOfSeqIds


def writeToFile(buffer, dir, fileName):
    """
        Writes a string to a file.
    """
    try:
        filePath = os.path.join(os.path.normpath(dir), fileName)
        f = open(filePath, 'w')
        f.write(buffer)
    except Exception:
        print "Cannot create a file or write to it:", filePath
        raise
    finally:
        f.close()


def ncbidToParentNcbids(ncbid, databaseFile):
    """
        @param ncbid
        @return: set of ncbids that lie
    """
    parentSet = set([])
    try:
        conn = sqlite3.connect(os.path.normpath(databaseFile))
        cursor = conn.cursor()

        if ncbid == 1:
            return parentSet

        taxonNcbid = ncbid

        while True:

            if taxonNcbid == 1: #the root of the taxonomy reached
                break
            cursor.execute('SELECT taxon_id FROM taxon T WHERE T.ncbi_taxon_id=?',(taxonNcbid,))
            result = cursor.fetchall()

            assert len(result) == 1, str('Cannot find taxon_id for ncbi: ' + str(taxonNcbid)) #temp
            taxonId = result[0][0]

            #get parent
            cursor.execute('SELECT parent_taxon_id FROM taxon T WHERE T.taxon_id=?', (taxonId,))
            result = cursor.fetchall()
            assert len(result) == 1, str('Cannot find parent for taxon_id', taxonId)
            taxonNcbid = result[0][0]
            parentSet.add(taxonNcbid)

        return parentSet

    except Exception:
        print "Failed to create connection to a database:", databaseFile
        raise
    finally:
        cursor.close()
        conn.close()


def noNewLine(str):
    return str.replace('\n', '').replace('\r','')


def lastNcbid(ncbids):
    list = ncbids.split(';')
    return int(list[len(list)-2])

if __name__ == "__main__":
  main()