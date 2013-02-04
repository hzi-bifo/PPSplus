#!/usr/bin/env python

"""
    Creates the marker gene database for the 16S and 23S genes as well as for the Amphora marker genes.

    The input 16S and 23S marker genes are generated using the "arb" software from the files downloaded from the
    SILVA database. The database with the "Amphora" marker genes was generated using script "download_mg.py"
    that downloads sequences from NCBI.
"""

import os
import glob
import argparse

from com import csv
from com import fasta as fas
from com import taxonomy_ncbid as tax


class _TaxonomyWrap():
    def __init__(self, taxonomyFile,
                 allowedRanks=['root','superkingdom','phylum','class','order','family','genus','species']):
        self._taxonomy = tax.TaxonomyNcbi(taxonomyFile)
        self._allowedRanks = allowedRanks
        self._parentDict = dict([]) # map: taxonId -> parent taxonId
        self._rankDict = dict([]) # map: taxonId -> rank

    def _getParent(self, taxonId):
        if taxonId in self._parentDict:
            return self._parentDict[taxonId]
        else:
            parent = self._taxonomy.getParentNcbid(taxonId)
            self._parentDict[taxonId] = parent
            return parent

    def _getRank(self, taxonId):
        if taxonId in self._rankDict:
            return self._rankDict[taxonId]
        else:
            rank = self._taxonomy.getRank(taxonId)
            self._rankDict[taxonId] = rank
            return rank

    def getPathToRoot(self, taxonId):
        rankToTaxonId = dict([])
        current = taxonId
        while current != 1:
            rank = self._getRank(current)
            if rank in self._allowedRanks:
                rankToTaxonId[rank] = current
            current = self._getParent(current)
            if current is None:
                return None
        ranks  = self._allowedRanks[1:]
        ranks.reverse()
        pathList = []
        current = taxonId
        for rank in ranks:
            if rank in rankToTaxonId:
                current = rankToTaxonId[rank]
            pathList.append(current)

        pathList.reverse()
        return str(";".join(map(str, pathList)) + ';')


def _main():

    parser = argparse.ArgumentParser(description=__doc__, epilog="""""")

    parser.add_argument('-i', '--input-data-dir', action='store', nargs=1, required=True,
        help=str('Directory that contains fasta files and corresponding mapping files, for each "*.tax" (or "*.csv") ' +
                 'file there must be a "*.fna" file with the same name. All files with suffix "tax" (or "*.csv") ' +
                 'will be considered.'),
        metavar='input_dir',
        dest='inDir')

    parser.add_argument('-o', '--output-dir', action='store', nargs=1, required=True,
        help='Directory that contains the output files.',
        metavar='out_dir',
        dest='outDir')

    parser.add_argument('-s', '--source-type', required=True,
        help='To determine the source, use "s" for the Silva database and "a" the Amphora database.',
        choices=["s","a"],
        dest='srcType')

    parser.add_argument('-t', '--taxonomy-file', nargs=1, type=file, required=True,
        help='NCBI taxonomy database file in the sqlite3 format.', metavar='ncbitax_sqlite.db',
        dest='taxonomy')

    args = parser.parse_args()

    inDir = args.inDir[0]
    outDir =  args.outDir[0]
    srcType = args.srcType[0]
    taxonomy = _TaxonomyWrap(args.taxonomy[0].name)

    for dir in [inDir, outDir]:
        assert os.path.isdir(dir), 'Path: "' + dir + '" does not exists!'

    mapDict = dict([])  # map: seqId -> ncbid
    fastaDict = dict([]) # map: seqId -> seq
    for mapFilePath in glob.glob(os.path.join(os.path.normpath(inDir), r'*.[ct][sa][vx]')): # *.csv or *.tax
        assert mapFilePath.endswith('.csv') or mapFilePath.endswith('.tax'), str(
        'The mapping files can either end with .csv or .tax ' + mapFilePath)
        base = os.path.basename(mapFilePath).rsplit('.', 1)[0]

        fastaDict = fas.fastaFileToDict(os.path.join(os.path.dirname(mapFilePath), (base + '.fna')))
        print("Processing: " + base + " seq:" + str(len(fastaDict)))

        if 'a' in srcType: # Amphora
            mapDict = dict([])
            for k in csv.getColumnAsList(mapFilePath, colNum=0, sep='\t'):
                v =  k.rsplit('|', 1)[1].split(':')[1]
                assert ((k not in mapDict) or (mapDict[k] == v)), str(
                    'There are at least two different values for key: ' + str(k) + ' in ' + mapFilePath)
                mapDict[k] = v
        elif 's' in srcType: # Silva
            mapTmp = csv.getMapping(mapFilePath, 0, 2, '\t')
            mapDict = dict([])
            for k, v in mapTmp.iteritems():
                mapDict[k] = v[0]
        else:
            assert False, 'Unsupported source type!'

        if len(mapDict) != len(fastaDict):
            print(base + ': The mapping file and the corresponding fasta file have different number of entries: "' +
                  str(len(mapDict)) + '" "' + str(len(fastaDict)) + '" these files will be skipped!')
            continue

        count = len(csv.getColumnAsList(mapFilePath))
        if len(mapDict) != count:
            print(base + ': The mapping file contained duplicates! unique: ' +
                  str(len(mapDict)) + ' non-unique: ' + str(count))


        # store the data
        outDna = csv.OutFileBuffer(os.path.join(outDir,str(base + '.fna')))
        outTax = csv.OutFileBuffer(os.path.join(outDir,str(base + '.tax')))
        count = 0
        for seqId, taxonId in mapDict.iteritems():
            path = taxonomy.getPathToRoot(taxonId)
            if path is None:
                print('Could not find: ' + str(taxonId) + ' for seqId: ' + seqId)
                continue
            seq = fastaDict[seqId]
            if 'a' in srcType: # Amphora
                id = seqId
            elif 's' in srcType: # Silva
                id = str(seqId + '|ncbid:' + str(taxonId))

            outTax.writeText(str(id + '\t' + path + '\n'))
            outDna.writeText(str('>' + id + '\n' + seq + '\n'))
            count += 1

        outDna.close()
        outTax.close()
        print('Stored entries: ' + str(count))

        # Silva:
        #-i /Users/ivan/Documents/work/binning/database/silva111/arbGenerated -s s -t /Users/ivan/Documents/work/binning/taxonomy/20121122/ncbitax_sqlite.db
        # -o /Users/ivan/Documents/work/binning/database/silva111/db

        # Amphora
        # -i /Users/ivan/Documents/work/binning/database/markerGenes3/mGenesExtracted -s a -t /Users/ivan/Documents/work/binning/taxonomy/20121122/ncbitax_sqlite.db
        # -o /Users/ivan/Documents/work/binning/database/markerGenes3/db

    print 'done'

def _testTaxonomyWrap():
    taxonomy = _TaxonomyWrap('/Users/ivan/Documents/work/binning/taxonomy/20121122/ncbitax_sqlite.db')
    print(str(taxonomy.getPathToRoot(870603)))

if __name__ == "__main__":
    _main()
    #_testTaxonomyWrap()