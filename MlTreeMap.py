#!/usr/bin/env python

from __future__ import with_statement
import os
import glob
import re
from Sequences import Sequences
from Taxonomy import Taxonomy
from Config import Config
from Config import Config2

def placeSequences(sequences, taxonomy, final_RAxML_outputs):
    """
        Go through the final_RAxML_outputs directory and assign/place sequences accordingly.
    """
    for filePath in glob.glob(os.path.join(os.path.normpath(final_RAxML_outputs),r'*.txt')):

        file = re.findall(r'p_[0-9]+_[0-9]+[_]+RAxML_parsed.txt', filePath)
        assert len(file) == 1, str('Could not get a filename from the filePath: ' + filePath)
        pair = re.findall(r'[0-9]+',file[0])
        assert len(pair) == 2, str('Wrong filename format!')
        scaffoldId = int(pair[0])
        contigId = int(pair[1])

        with open(filePath, 'r') as f:
            content = f.read()
            f.close()

        placements = re.findall(r'[0-9]+%[^%]+\)\.', content)
        map = dict([]) # map: percent -> placement
        for p in placements:
            percent =  int(re.sub(r'^([0-9]+)%:\ .*\.$',r'\1' ,p))
            placement =  re.sub(r'^[0-9]+%:\ (.*\.)$',r'\1' ,p)
            if percent not in map:
                map[percent] = placement
        maxPlacement = map[max(map.keys())]
        ncbisBracket = re.findall(r'\([0-9]+\)', maxPlacement)
        ncbis = []
        for ncbiBracket in ncbisBracket:
            ncbi = re.sub(r'\(([0-9]+)\)', r'\1', ncbiBracket)
            ncbis.append(ncbi)

        taxPathDict = taxonomy.getPathFromLowestCommonAncestorToRoot(ncbis)
        if taxPathDict.keys() >= 1:
            sequences.setTaxonomyPath(contigId, scaffoldId, taxPathDict)
        else:
            raise Exception('No taxonomic path found for ncbids: ', ncbis)

        #test
        #for k in Config.taxonomicRanks:
        #    if k not in taxPathDict:
        #        break
        #    n = taxPathDict[k]
        #    print n.ncbid, n.rank, n.name


def test():
    config = Config(open(os.path.normpath('D:\\A_Phylo\\A_Metagenomic\\pPPS\\workspace\\pPPS\\config01.cfg')), 'pPPS')

    #read sequences
    sequences = Sequences(config)

    #write ids file
    sequences.writeSequences(config.get('inputIdsFastaFile'))

    #init taxonomy
    databaseFile = os.path.normpath(config.get('databaseFile'))
    taxonomicRanks = config.get('taxonomicRanks').split(',')
    taxonomy = Taxonomy(databaseFile, taxonomicRanks)

    config2 = Config2(config, 'MLTreeMap')

    #place sequences
    placeSequences(sequences, taxonomy, config2.get('final_RAxML_outputs'))
    print 'Placed sequences:', len(sequences.placedSeqSet)

    #assign not placed contigs of one scaffold to the lowest common ancestor of assigned contigs
    if eval(config.get('placeContigsFromTheSameScaffold')):
        sequences.placeContigsFromTheSameScaffold(taxonomy)
        print 'Not assigned contigs of one scaffold are placed to the lowest common ancestor of assigned contigs'
        print 'Placed sequences now:', len(sequences.placedSeqSet)
    else:
        print 'Not assigned contigs from the same scaffold are NOT placed'

    taxonomy.close()

    #print placements for each sequence that is assigned/placed
    print '--------------------------------'
    #for s in sequences.sequences:
    #    dict = sequences.getTaxonomyPath(s.id)
    #    if dict == None:
    #        continue
    #    print s.name, '...'
    #    for rank in taxonomicRanks:
    #        if rank not in dict:
    #            break
    #        print rank, dict[rank].ncbid, dict[rank].name
    #    print ''




if __name__ == "__main__":
  test()
