#!/usr/bin/env python

"""
    Copyright (C) 2014  Ivan Gregor

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

    Note that we could have written some parts of this code in a nicer way,
    but didn't have time. Be careful when reusing the source code.


    Prepare scripts for the Ruben's comparison charts.

    Input:
        an array of files containing taxonomic assignments.
        an input FASTA file
        scale for the circles

    Output:
        a list of leaf clades
        abundance file (ncbid TAB abundance1 TAB abundance2 TAB abundance3 ...)



"""

import algbioi.com.fasta as fas
import algbioi.com.csv as csv
import algbioi.com.taxonomy_ncbi as tax
import math

TAXONOMY_FILE = "/Users/ivan/Documents/work/binning/taxonomy/20121122/ncbitax_sqlite.db"

def _main():

    storeClades = False
    storeAbundancies = True

    # totalBp = 161343000.0

    # taxonomy
    taxonomy = tax.TaxonomyNcbi(TAXONOMY_FILE)

    # output
    # abundanceFile = "/Users/ivan/Documents/work/Doc/PPSplus/abundance_plots/cmp_expert/abundance.txt"
    # abundanceFile = "/Users/ivan/Documents/work/Doc/PPSplus/abundance_plots/cmp_sim/abundance.txt"
    abundanceFile = '/Users/ivan/Documents/work/Doc/PPSplus/abundance_plots/cmp_sim2/abundance.txt'

    # leaf nodes
    # leafNodesFile = "/Users/ivan/Documents/work/Doc/PPSplus/abundance_plots/cmp_expert/leaf_nodes.txt"
    # leafNodesFile = "/Users/ivan/Documents/work/Doc/PPSplus/abundance_plots/cmp_sim/leaf_nodes.txt"
    leafNodesFile = '/Users/ivan/Documents/work/Doc/PPSplus/abundance_plots/cmp_sim2/leaf_nodes.txt'


    # fasta file
    # fastaFile = "/Users/ivan/Documents/work/binning/data/HumanGut/working/scaffolds.fna"
    fastaFile = "/Users/ivan/Documents/work/binning/data/mercier050513/melanieAssembly/lognorm/VelvetNoscafMinimus2_ma_merge_lognorm_min_1000bp.fna"

    # a list of assignments
    # predListFile = ["/Users/ivan/Documents/work/Doc/PPSplus/evaluation_real/hg_cmp_to_expert/hg_scaffolds_mg.csv",
    #                 "/Users/ivan/Documents/work/Doc/PPSplus/evaluation_real/hg_cmp_to_expert/hg_scaffolds_pp.csv",
    #                 "/Users/ivan/Documents/work/Doc/PPSplus/evaluation_real/hg_cmp_to_expert/hg_scaffolds_ppsp.csv"]

    # predListFile = ["/Users/ivan/Documents/work/binning/data/mercier050513/melanieAssembly/lognorm/VelvetNoscafMinimus2_ma_merge_lognorm_min_1000bp.tax",  # TRUE
    #                  "/Users/ivan/Documents/work/Doc/PPSplus/abundance_plots/cmp_sim/ppsp_lognorm_mg_no_rs_species_nc.out",  # PPS+ exclude MG exclude nothing, RS exclude species
    #                  "/Users/ivan/Documents/work/Doc/PPSplus/abundance_plots/cmp_sim/megan_lognorm_species_nc.out",  # MEGAN exclude species
    #                  "/Users/ivan/Documents/work/Doc/PPSplus/abundance_plots/cmp_sim/taxator_lognorm_species_nc.tax",  # Taxator-tk exclude species
    #                  "/Users/ivan/Documents/work/Doc/PPSplus/abundance_plots/cmp_sim/pps_lognorm_genus_nc.out"]  # PPS exclude species

    predListFile = ["/Users/ivan/Documents/work/binning/data/mercier050513/melanieAssembly/lognorm/VelvetNoscafMinimus2_ma_merge_lognorm_min_1000bp.tax",  # TRUE
                     "/Users/ivan/Documents/work/Doc/PPSplus/abundance_plots/cmp_sim2/ppsp_lognorm_mg_species_rs_species_nc.out",  # PPS+ exclude species
                     "/Users/ivan/Documents/work/Doc/PPSplus/abundance_plots/cmp_sim/megan_lognorm_species_nc.out",  # MEGAN exclude species
                     "/Users/ivan/Documents/work/Doc/PPSplus/abundance_plots/cmp_sim/taxator_lognorm_species_nc.tax",  # Taxator-tk exclude species
                     "/Users/ivan/Documents/work/Doc/PPSplus/abundance_plots/cmp_sim/pps_lognorm_genus_nc.out"]  # PPS exclude species



    # taxonIdOrderFile = '/Users/ivan/Documents/work/Doc/PPSplus/abundance_plots/taxon_id_order.txt'
    # taxonIdOrderFile = "/Users/ivan/Documents/work/Doc/PPSplus/abundance_plots/cmp_sim/taxon_id_order2.txt"
    # taxonIdOrderFile = "/Users/ivan/Documents/work/Doc/PPSplus/abundance_plots/cmp_sim2/taxon_id_order.txt"

    seqIdToBp = fas.getSequenceToBpDict(fastaFile)

    toParent = {}



    # transform assignments to mapping: seqId -> ncbid
    predDictList = []
    for predFile in predListFile:
        predDictList.append(mapToAllowedRanks(taxonomy, csv.predToDict(predFile)))

    # get all taxIds emerging in files
    taxonIdSet = set()
    for d in predDictList:
        for taxonId in d.values():
            taxonIdSet.add(abs(taxonId))

    # taxon ids of all nodes
    allTaxonIdSet = set([1])
    parentTaxonIdSet = set()
    for taxonId in taxonIdSet:
        if taxonId not in allTaxonIdSet:
            id = taxonId
            while int(id) != 1:
                allTaxonIdSet.add(id)
                idChild = id
                id = taxonomy.getParentNcbid(id)
                toParent[idChild] = id
                parentTaxonIdSet.add(id)


    # taxon ids of leaf nodes
    leafTaxonIdSet = set()
    for taxonId in taxonIdSet:
        if taxonId not in parentTaxonIdSet:
            leafTaxonIdSet.add(taxonId)

    # map: taxonId -> list of parental nodes
    idToList = {}
    for taxonId in allTaxonIdSet:
        l = []
        id = taxonId
        while int(id) != 1:
            l.append(id)
            id = toParent[id]
        l.append(1)
        idToList[taxonId] = l

    # abundances
    outList = []
    errorCount = 0
    for predDict in predDictList:  # list of map: seqId -> ncbid
        ncbidToBp = {}
        for id in allTaxonIdSet:
            ncbidToBp[id] = 0
        for seqId, ncbid in predDict.iteritems():
            l = idToList[abs(ncbid)]
            try:
                bp = seqIdToBp[seqId]
            except:
                bp = 0
                print seqId
                errorCount += 1

            for id in l:
                ncbidToBp[id] += bp
        outList.append(ncbidToBp)
    print("Key errors: %s" % errorCount)


    if storeClades:
        out = csv.OutFileBuffer(leafNodesFile)
        for taxonId in leafTaxonIdSet:
            out.writeText(str(taxonId) + '\n')


    # count for each node how many children it has
    d = {}
    for taxonId in allTaxonIdSet:
        d[taxonId] = 0

    for taxonId in allTaxonIdSet:
        parentId = taxonomy.getParentNcbid(taxonId)
        if parentId in d:
            d[parentId] += 1

    degreeTwoSet = set()
    for taxonId, count in d.iteritems():
        if count == 1:
            degreeTwoSet.add(taxonId)


    #
    if storeAbundancies:
        out = csv.OutFileBuffer(abundanceFile)
        # out.writeText('A' + '\t' + 'B' + '\t' + 'C' + '\n')
        out.writeText('A' + '\t' + 'B' + '\t' + 'C' + '\t' + 'D' + '\t' + 'E' + '\n')

        #for taxonId in allTaxonIdSet:
        # taxonIdOrder = map(int, csv.getColumnAsList(taxonIdOrderFile))
        #for taxonId in taxonIdOrder:
        for taxonId in allTaxonIdSet:
            line = str(taxonId)
            for entry in outList:
                # line += '\t' + str(math.log(float(entry[taxonId]) / float(totalBp)))
                if (entry[taxonId] < 0.0000001): # or (taxonId in degreeTwoSet):  # nodes with degree 2 will have zero abundance
                    line += '\t' + '0'
                else:
                    line += '\t' + str(0.01* math.log(float(entry[taxonId]), 30))
                    # line += '\t' + str(float(entry[taxonId]))
            out.writeText(line + '\n')
        out.close()


    taxonomy.close()
    print("done")


def replNamesByIds():

    taxonomy = tax.TaxonomyNcbi(TAXONOMY_FILE)
    treeIds = '/Users/ivan/Documents/work/Doc/PPSplus/abundance_plots/cmp_sim2/leaf_nodes.newick'
    treeNames = '/Users/ivan/Documents/work/Doc/PPSplus/abundance_plots/cmp_sim2/tree_taxon_names.newick'
    leafNodes = '/Users/ivan/Documents/work/Doc/PPSplus/abundance_plots/cmp_sim2/leaf_nodes.txt'

    leafSet = set(map(int, csv.getColumnAsList(leafNodes)))


    line = open(treeIds).readline()
    lineList = list(line)
    print lineList

    out = ""
    a = lineList
    i = 0
    l = []
    while i < len(a):
        if '0' <= a[i] and a[i] <= '9':
            l.append(a[i])
        else:
            if len(l) > 0:
                ncbid = int("".join(l))
                if ncbid in leafSet:
                    out += taxonomy.getScientificName(ncbid)
                else:
                    out += str(ncbid)
                l = []
            out += a[i]
        i += 1

    print out

    taxonomy.close()
    outBuffer = csv.OutFileBuffer(treeNames)
    outBuffer.writeText(out)
    outBuffer.close()


def mapToAllowedRanks(taxonomy, inPredDict):
    """
        Maps the values of the input dictionary that represent taxon ids to the major taxonomic ranks and
        outputs the dictionary.
    """
    allowedRanks = set(tax.TAXONOMIC_RANKS)
    outPredDict = {}
    toParent = {}

    for seqId, taxonId in inPredDict.iteritems():
        id = taxonId
        while (id != 1) and (taxonomy.getRank(id) not in allowedRanks):

            if id in toParent:
                id = toParent[id]
            else:
                id0 = id
                id = taxonomy.getParentNcbid(id)
                if id is None:
                    id = 1
                toParent[id0] = id

        outPredDict[seqId] = id
    return outPredDict




if __name__ == "__main__":
    # replNamesByIds()
    _main()