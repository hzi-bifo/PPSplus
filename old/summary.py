#!/usr/bin/env python

import os
import sys
import re
import argparse
import glob
import sqlite3


def main():
    """
    usage: -t ncbitax_sqlite.db -g genomesDir -i predictions.txt -f tw.fna
    """
    parser = argparse.ArgumentParser(description='''List clades that can be modeled by PhyloPythiaS''',
                                     epilog='''Format:#seq\t#kb\tNCBID_name;...;NCBID_name;\tNCBID''')

    parser.add_argument('-t', '--taxonomy', action='store', nargs=1, required=True,
                        help='NCBI taxonomy in sqlite3 format', metavar='ncbitax_sqlite.db',
                        dest='taxonomy')

    parser.add_argument('-g', '--genomes', action='store', nargs=1, required=True,
                        help='Directory containing all genomes available for PPS',
                        dest='genomes')

    parser.add_argument('-i', '--inputFile', action='store', nargs=1,
                        help='If the input file is not given, the script takes stdin as the input',
                        dest='inFile')

    parser.add_argument('-f', '--fastaFile', action='store', nargs=1,
                        help='Fasta file',
                        dest='fasta')

    args = parser.parse_args()

    if not args.taxonomy:
        sys.stderr.write('Enter correct taxonomy')
        sys.stderr.write(parser.print_help())
    else:
        databaseFile = args.taxonomy[0]

    if not args.genomes:
        sys.stderr.write('Enter correct genomes')
        sys.stderr.write(parser.print_help())
    else:
        genomesDir = os.path.normpath(args.genomes[0])

    if not args.fasta:
        sys.stderr.write('Enter correct path to a fasta file')
        sys.stderr.write(parser.print_help())
    else:
        fastaFile = os.path.normpath(args.fasta[0])

    if args.inFile:
        inFile = os.path.normpath(args.inFile[0])
    else:
        inFile = sys.stdin

    #read the fasta file and get map: seq_id -> seq_len
    seqIdToLen = readFasta(fastaFile)

    #read input file/stdin and get map: ncbid -> list of seq_ids
    ncbidToListOfSeqIds = readInput(inFile)

    #compute #seq and get map: ncbid -> #seq
    ncbidToSeqCount = dict([])
    for ncbid in ncbidToListOfSeqIds:
        ncbidToSeqCount[ncbid] = len(ncbidToListOfSeqIds[ncbid])

    #compute #bp and get map: ncbid -> #kb
    ncbidToBp = dict([])
    for ncbid in ncbidToListOfSeqIds:
        list = ncbidToListOfSeqIds[ncbid]
        bp = 0
        for seq in list:
            bp += seqIdToLen[seq]

        ncbidToBp[ncbid] = bp

    #for each ncbid compute # of genomes
    taxonomicRanks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    ncbidToGenomeCount = dict([])
    for ncbid in ncbidToListOfSeqIds:
        ncbidToGenomeCount[ncbid] = getGenomeWgsCount(ncbid, 3, genomesDir, databaseFile, taxonomicRanks)


    #go over all ncbid and output the result
    print "bpCount\tseqCount\tgenomeCount\tlineage\tncbid"
    for ncbid in ncbidToListOfSeqIds:
        print str(str(ncbidToBp[ncbid]) + '\t' +  str(ncbidToSeqCount[ncbid]) +
                  '\t' + str(ncbidToGenomeCount[ncbid]) + '\t' +
                  str(ncbidToPath(ncbid, databaseFile, taxonomicRanks))
                  +
                  '\t'+ str(ncbid) )



def readFasta(filePath):
    seqIdToLen = dict([])
    try:
        f = open(os.path.normpath(filePath),'r')
    except Exception:
        print "Cannot open file:", filePath
        raise
    else:
        name = ''
        seq = ''
        for line in f:
            line = noNewLine(line)
            if re.match('>', line):
                if seq != '':
                    assert name != ''
                    seqIdToLen[name] = len(seq)
                    seq = ''
                name = line.replace('>','')
            else:
                seq += line
        if seq != '':
            assert name != ''
            seqIdToLen[name] = len(seq)
    finally:
        f.close()
    return seqIdToLen


def readInput(inFile):
    ncbidToListOfSeqIds = dict([])
    try:
        f = open(os.path.normpath(inFile),'r')
    except Exception:
        print "Cannot open file:", inFile
        raise
    else:

        for line in f:
            line = noNewLine(line)
            seqId = re.sub(r'^([^\t]+)\t.*',r'\1', line)
            ncbids = re.sub(r'^[^\t]+\t(.*)',r'\1', line)
            ncbid = lastNcbid(ncbids)
            if ncbid not in ncbidToListOfSeqIds:
                ncbidToListOfSeqIds[ncbid] = []
            ncbidToListOfSeqIds[ncbid].append(seqId)
    finally:
        f.close()
    return ncbidToListOfSeqIds

def lastNcbid(ncbids):
    list = ncbids.split(';')
    return int(list[len(list)-2])


def getGenomeWgsCount(ncbid, threshold, dir, databaseFile, taxonomicRanks):
    """
        Get the number of genomes/wgs available in the directory

        @param threshold: it returns max. threshold genomes/wgs (after this threshold is reached, it returns)
        @param dir: directory that contains genomes/wgs in the form: "ncbid.[0-9]*.f[an][sa]"

        @return: the number of genomes/wgs from different species that are subclades of the input ncbid
    """
    try:
        conn = sqlite3.connect(os.path.normpath(databaseFile))
        cursor = conn.cursor()

        speciesIdsList = []
        collectSpecies(speciesIdsList, cursor, ncbid, dir, threshold)
        return len(speciesIdsList)
    except Exception:
        print "Failed to create connection to a database:", databaseFile
        raise
    finally:
        cursor.close()
        conn.close()


def genomeExists(listOfNcbids, dir):
    """
        Gets the first ncbid of species/subspecies from the input list that is contained in the directory

        @param dir: directory that contain genomes/wgs files

        @return None or the first ncbid for which there is a genome or a draft genome in the directory
    """
    for ncbid in listOfNcbids:
        if glob.glob(os.path.join(os.path.normpath(dir), str(str(ncbid) + '.[0-9]*.f[an][sa]'))):
            return ncbid
    return None


def collectSpecies(speciesIds, cursor, root, dir, threshold):
    """
        @param speciesIds: output list of ncbids (species or subspecies) for which there are genomes/wgs in the directory
    """
    if len(speciesIds) >= threshold:
        return
    cursor.execute('SELECT node_rank FROM taxon T WHERE T.ncbi_taxon_id=?',(root,))
    result = cursor.fetchall()
    if len(result) == 0:
        return
    assert len(result) == 1
    if result[0][0] == 'species':
        list = []
        list.append(root)
        collectSubSpecies(list, cursor, root)
        ncbid = genomeExists(list, dir) #now check if at least one genome exists from the list
        if ncbid != None:
            speciesIds.append(ncbid)
    else:
        cursor.execute('SELECT ncbi_taxon_id FROM taxon T WHERE T.parent_taxon_id=?',(root,))
        result = cursor.fetchall()
        for item in result:
            collectSpecies(speciesIds, cursor, int(item[0]), dir, threshold)


def collectSubSpecies(list, cursor, root):
    """
        Gets ncbids of all subspecies of the root.

        @param list: output list with all ncbids of the subspecies
    """
    cursor.execute('SELECT ncbi_taxon_id FROM taxon T WHERE T.parent_taxon_id=?',(root,))
    result = cursor.fetchall()
    if len(result) > 0:
        for item in result:
            id = int(item[0])
            list.append(id)
            collectSubSpecies(list, cursor, id)


def ncbidToPath(ncbid, databaseFile, taxonomicRanks):
    try:
        conn = sqlite3.connect(os.path.normpath(databaseFile))
        cursor = conn.cursor()

        path = []
        if ncbid == 1:
            return path

        taxonNcbid = ncbid
        taxonRank = ''

        while True:

            if taxonNcbid == 1: #the root of the taxonomy reached
                break
            cursor.execute('SELECT taxon_id FROM taxon T WHERE T.ncbi_taxon_id=?',(taxonNcbid,))
            result = cursor.fetchall()

            assert len(result) == 1, str('Cannot find taxon_id for ncbi: ' + str(taxonNcbid)) #temp
            taxonId = result[0][0]

            taxonRank = result[0][0]

            cursor.execute('SELECT name FROM taxon_name TN WHERE TN.taxon_id=? AND name_class="scientific name"', (taxonId,))
            result = cursor.fetchall()
            assert len(result) == 1, str('Cannot find scientific name for taxon_id: ' + str(taxonId))
            taxonName = result[0][0]

            path.append(taxonName)

            #get parent
            cursor.execute('SELECT parent_taxon_id FROM taxon T WHERE T.taxon_id=?', (taxonId,))
            result = cursor.fetchall()
            assert len(result) == 1, str('Cannot find parent for taxon_id', taxonId)
            taxonNcbid = result[0][0]

        path.reverse()
        buffer = ''
        for s in path:
            buffer += s + ','
        return buffer

    except Exception:
        print "Failed to create connection to a database:", databaseFile
        raise
    finally:
        cursor.close()
        conn.close()



def noNewLine(str):
    return str.replace('\n', '').replace('\r','')


if __name__ == "__main__":
  main()