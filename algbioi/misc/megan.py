"""
    This file contains methods necessary to classify a sample using megan. Suppose blast 2.2.27.



"""
import os

from algbioi.com import fasta as fas
from algbioi.com import csv
from algbioi.com import taxonomy_ncbi
from algbioi.eval import accuracy
from algbioi.eval import confusion_matrix
from algbioi.eval import consistency
from algbioi.misc import out_proc
from algbioi.core import run
from algbioi.core.sequences import Sequences

class TaxonomyNcbiWrap():
    def __init__(self, databaseFile):
        self._taxonomy = taxonomy_ncbi.TaxonomyNcbi(databaseFile)
        self._taxonIdToSpeciesTaxonId = {}

    def getTaxonIdAtSpeciesRank(self, taxonId):
        id = self._taxonIdToSpeciesTaxonId.get(taxonId, None)
        if id is None:
            id = taxonId
            rank = self._taxonomy.getRank(id)
            while ((rank != 'species') and (id != 1)):
                id = self._taxonomy.getParentNcbid(id)
                rank = self._taxonomy.getRank(id)
            self._taxonIdToSpeciesTaxonId[taxonId] = id
            if id == 1:
                print("There is no ncbi taxon id defined at rank species for taxon id: %s" % taxonId)
        return id


def createFastaFile(inFasta, outFasta):
    """
        Add the length info to each fasta header. (E.g. >seqId -> seqId|length|1234)

        @param inFasta: input fasta file
        @param outFasta: output fasta file with length info
    """
    idToSeq = fas.fastaFileToDict(inFasta)
    out = csv.OutFileBuffer(outFasta)

    for seqId, seq in idToSeq.iteritems():
        out.writeText('>' + str(seqId) + '|length|' + str(len(''.join(seq.split()))) + '\n' + str(seq) + '\n')
    out.close()


def addTaxonIds(inBlastTabFile, outBlastTabFile, mapFile, taxonomy):
    """
        Adds taxonIds to the blast tab file as its last column.

        @param inBlastTabFile: in blast tab file in format: outfmt 7
        @param outBlastTabFile: out blast tab file with last column representing taxonId
        @param mapFile: mapping reference seqId -> species ncbid
        @type taxonomy: TaxonomyNcbiWrap
    """
    out = csv.OutFileBuffer(outBlastTabFile)
    refIdToTaxonId = csv.getMapping(mapFile, 0, 1, '\t')
    for line in open(inBlastTabFile, 'r'):
        line = line.strip()
        if line.startswith('#'):
            out.writeText(line + '\n')
            continue

        fields = line.split()
        if len(fields) < 2:
            continue

        refId = fields[1].strip()
        taxonId = int(refIdToTaxonId[refId][0])
        taxonId = int(taxonomy.getTaxonIdAtSpeciesRank(taxonId))
        if taxonId == 1:
            taxonId = -1
        line = line + '\t' + str(taxonId)
        out.writeText(line + '\n')
    out.close()


# filter out alignments at different ranks for the uniform and lognorm dataset !!!



def _prepareFastaAndAlignments():
    taxonomy = TaxonomyNcbiWrap('/Users/ivan/Documents/work/binning/taxonomy/20121122/ncbitax_sqlite.db')
    mapFile = '/Users/ivan/Documents/nobackup/megan/ref_all.tax'

    blastTabDir = '/Users/ivan/Documents/nobackup/megan/blastTab'
    blastTabMeganDir = '/Users/ivan/Documents/nobackup/megan/blastTabMegan'

    fastaMeganDir = '/Users/ivan/Documents/nobackup/megan/fastaMegan'
    fastaOrigDir = '/Users/ivan/Documents/nobackup/megan/fastaOrig'

    for f in os.listdir(blastTabDir):
        print(f)
        addTaxonIds(os.path.join(blastTabDir, f), os.path.join(blastTabMeganDir, f), mapFile, taxonomy)

    for f in os.listdir(fastaOrigDir):
        print(f)
        createFastaFile(os.path.join(fastaOrigDir, f), os.path.join(fastaMeganDir, f))


def filterOutAlignments(inAlignments, outAlignments, cladesToExclude, mapFile, refIdCol=1):
    """

        @param inAlignments: input blast tab alignments
        @param outAlignments: out filtered blast tab alignments
        @param cladesToExclude: taxonIds that needs to be excluded
        @param mapFile: map refId to taxonId
    """
    # get all refIds that match to cladesToExclude according to the mapFile
    taxonIdsToExclude = csv.getColumnAsList(cladesToExclude)
    taxonIdToRefIds = csv.getMapping(mapFile, 1, 0, sep='\t')

    refIdSet = set()
    for taxonId in taxonIdsToExclude:
        refIdSet.update(set(taxonIdToRefIds[taxonId]))

    out = csv.OutFileBuffer(outAlignments)
    # go through the alignment file and filter out all entries matching the refIdSet

    excludeCount = 0
    for line in open(inAlignments):
        line = line.strip()
        if line.startswith('#'):
            out.writeText(line + '\n')
        else:
            fields = line.split('\t')
            if fields >= refIdCol + 1:
                entry = fields[refIdCol].strip()
                if entry not in refIdSet:
                    out.writeText(line + '\n')
                else:
                    excludeCount += 1
    out.close()
    print('Alignments excluded: %s %s %s' % (os.path.basename(inAlignments), os.path.basename(cladesToExclude),
                                             excludeCount))


def _filterAlignments():
    ranks = ['strain', 'species', 'genus', 'family']
    clades = {}
    clades['strain'] = '/Users/ivan/Documents/nobackup/megan/exclude_clades_strain.txt'
    clades['species'] = '/Users/ivan/Documents/nobackup/megan/exclude_clades_species.txt'
    clades['genus'] = '/Users/ivan/Documents/nobackup/megan/exclude_clades_genus.txt'
    clades['family'] = '/Users/ivan/Documents/nobackup/megan/exclude_clades_family.txt'

    mapFile = '/Users/ivan/Documents/nobackup/megan/ref_all.tax'

    # for MEGAN
    # inFileList = ['/Users/ivan/Documents/nobackup/megan/blastTabMegan/VelvetNoscafMinimus2_ma_merge_uniform_min_1000bp_m7.blast',
    #            '/Users/ivan/Documents/nobackup/megan/blastTabMegan/VelvetNoscafMinimus2_ma_merge_lognorm_min_1000bp_m7.blast']
    #
    # for TAXATOR
    inFileList = ['/Users/ivan/Documents/nobackup/taxator/raw_alignments/VelvetNoscafMinimus2_ma_merge_uniform_min_1000bp_alignments',
                  '/Users/ivan/Documents/nobackup/taxator/raw_alignments/VelvetNoscafMinimus2_ma_merge_lognorm_min_1000bp_alignments']

    for f in inFileList:
        for rank in ranks:
            # for MEGAN
            # filterOutAlignments(f, str('.'.join(f.split('.')[0:-1]) + '_' + rank + '.blast'), clades[rank], mapFile)
            # FOR TAXATOR
            filterOutAlignments(f, str('_'.join(f.split('_')[0:-1]) + '_' + rank + '_alignments'), clades[rank],
                                mapFile, refIdCol=4)

def _checkMagenPred(predDir):
    for f in os.listdir(predDir):
        print f
        _checkDuplicity(os.path.join(predDir, f))


def _checkDuplicity(inPredFile):
    idToTaxonIdList = csv.getMapping(inPredFile, 0, 1, '\t')
    for id, taxonIdList in idToTaxonIdList.iteritems():
        if len(taxonIdList) > 1:
            print id, taxonIdList
            print inPredFile
            break


def _evalSimDataset(predFile, labelFile, fastaFile, outDir, databaseFile,
                    recallMinFracClade, precisionMinFracPred, correctLabelThreshold):
    """
        For the simulated datasets, compute accuracy, accuracy corrected, and confusion matrices.

        @param predFile: megan prediction file
        @param trueLabelsFile: true labels
        @param fastaFile:
        @param outDir: directory containing the results
    """
    # create output directories
    rootDir = os.path.join(outDir, ''.join(os.path.basename(predFile).split('.')[0:-1]))

    if not os.path.isdir(rootDir):
        os.mkdir(rootDir)
    tablesDir = os.path.join(rootDir, 'tables')
    if not os.path.isdir(tablesDir):
        os.mkdir(tablesDir)

    # compute accuracy
    buff = csv.OutFileBuffer(os.path.join(rootDir, 'precision_recall.csv'))
    acc = accuracy.Accuracy(fastaFile, predFile, labelFile, databaseFile)
    buff.writeText(acc.getAccuracyPrint(taxonomy_ncbi.TAXONOMIC_RANKS[1:], minFracClade=recallMinFracClade,
                                        minFracPred=precisionMinFracPred, overview=True))
    buff.close()
    acc.close()

    # compute accuracy corrections
    buff = csv.OutFileBuffer(os.path.join(rootDir, 'precision_recall_corrections.csv'))
    acc = accuracy.Accuracy(fastaFile, predFile, labelFile, databaseFile, correctLabelThreshold)
    buff.writeText(acc.getAccuracyPrint(taxonomy_ncbi.TAXONOMIC_RANKS[1:], minFracClade=recallMinFracClade,
                                        minFracPred=precisionMinFracPred, overview=True))
    buff.close()
    acc.close()

    # compute confusion tables
    ranks = taxonomy_ncbi.TAXONOMIC_RANKS[1:]
    cm = confusion_matrix.ConfusionMatrix(fastaFile, predFile, labelFile, databaseFile, ranks)
    for rank in ranks:
        cm.generateConfusionMatrix(rank, os.path.join(tablesDir, os.path.basename(fastaFile)))
    cm.close()


def _evalSimDatasetBatch():
    # predFileArray = [['/Users/ivan/Documents/nobackup/megan/predictions/VelvetNoscafMinimus2_ma_merge_uniform_min_1000bp_m7-ex.txt',
    #                   '/Users/ivan/Documents/nobackup/megan/predictions/VelvetNoscafMinimus2_ma_merge_uniform_min_1000bp_m7_strain-ex.txt',
    #                   '/Users/ivan/Documents/nobackup/megan/predictions/VelvetNoscafMinimus2_ma_merge_uniform_min_1000bp_m7_species-ex.txt',
    #                   '/Users/ivan/Documents/nobackup/megan/predictions/VelvetNoscafMinimus2_ma_merge_uniform_min_1000bp_m7_genus-ex.txt',
    #                   '/Users/ivan/Documents/nobackup/megan/predictions/VelvetNoscafMinimus2_ma_merge_uniform_min_1000bp_m7_family-ex.txt'],
    #
    #                  ['/Users/ivan/Documents/nobackup/megan/predictions/VelvetNoscafMinimus2_ma_merge_lognorm_min_1000bp_m7-ex.txt',
    #                   '/Users/ivan/Documents/nobackup/megan/predictions/VelvetNoscafMinimus2_ma_merge_lognorm_min_1000bp_m7_strain-ex.txt',
    #                   '/Users/ivan/Documents/nobackup/megan/predictions/VelvetNoscafMinimus2_ma_merge_lognorm_min_1000bp_m7_species-ex.txt',
    #                   '/Users/ivan/Documents/nobackup/megan/predictions/VelvetNoscafMinimus2_ma_merge_lognorm_min_1000bp_m7_genus-ex.txt',
    #                   '/Users/ivan/Documents/nobackup/megan/predictions/VelvetNoscafMinimus2_ma_merge_lognorm_min_1000bp_m7_family-ex.txt']]
    predFileArray = [['/Users/ivan/Documents/nobackup/taxator/predictions/VelvetNoscafMinimus2_ma_merge_uniform_min_1000bp.tax',
                      '/Users/ivan/Documents/nobackup/taxator/predictions/VelvetNoscafMinimus2_ma_merge_uniform_min_1000bp_strain.tax',
                      '/Users/ivan/Documents/nobackup/taxator/predictions/VelvetNoscafMinimus2_ma_merge_uniform_min_1000bp_species.tax',
                      '/Users/ivan/Documents/nobackup/taxator/predictions/VelvetNoscafMinimus2_ma_merge_uniform_min_1000bp_genus.tax',
                      '/Users/ivan/Documents/nobackup/taxator/predictions/VelvetNoscafMinimus2_ma_merge_uniform_min_1000bp_family.tax'],

                    ['/Users/ivan/Documents/nobackup/taxator/predictions/VelvetNoscafMinimus2_ma_merge_lognorm_min_1000bp.tax',
                     '/Users/ivan/Documents/nobackup/taxator/predictions/VelvetNoscafMinimus2_ma_merge_lognorm_min_1000bp_strain.tax',
                     '/Users/ivan/Documents/nobackup/taxator/predictions/VelvetNoscafMinimus2_ma_merge_lognorm_min_1000bp_species.tax',
                     '/Users/ivan/Documents/nobackup/taxator/predictions/VelvetNoscafMinimus2_ma_merge_lognorm_min_1000bp_genus.tax',
                     '/Users/ivan/Documents/nobackup/taxator/predictions/VelvetNoscafMinimus2_ma_merge_lognorm_min_1000bp_family.tax']]


    labelFile = ['/Users/ivan/Documents/work/binning/data/mercier050513/melanieAssembly/uniform/VelvetNoscafMinimus2_ma_merge_uniform_min_1000bp.tax',
                 '/Users/ivan/Documents/work/binning/data/mercier050513/melanieAssembly/lognorm/VelvetNoscafMinimus2_ma_merge_lognorm.tax']

    fastaFile = ['/Users/ivan/Documents/work/binning/data/mercier050513/melanieAssembly/uniform/VelvetNoscafMinimus2_ma_merge_uniform_min_1000bp.fna',
                 '/Users/ivan/Documents/work/binning/data/mercier050513/melanieAssembly/lognorm/VelvetNoscafMinimus2_ma_merge_lognorm_min_1000bp.fna']

    # outDir = '/Users/ivan/Documents/nobackup/megan/evaluation'
    outDir = '/Users/ivan/Documents/nobackup/taxator/evaluation'
    databaseFile = '/Users/ivan/Documents/work/binning/taxonomy/20121122/ncbitax_sqlite.db'
    recallMinFracClade = 0.001
    precisionMinFracPred = 0.001
    correctLabelThreshold = 0.9

    for predFileList, label, fasta in zip(predFileArray, labelFile, fastaFile):
        for pred in predFileList:
            _evalSimDataset(pred, label, fasta, outDir, databaseFile, recallMinFracClade,
                            precisionMinFracPred, correctLabelThreshold)


def _evalRealDataset():
    contigPredArray = ['/Users/ivan/Documents/nobackup/taxator/predictions/cr_chunks2000_contigs_alignment.tax', '/Users/ivan/Documents/nobackup/taxator/predictions/hg_contigs_alignment.tax']
    scaffPredArray = ['/Users/ivan/Documents/nobackup/taxator/predictions/cr_chunks2000_scaffolds_alignment.tax', '/Users/ivan/Documents/nobackup/taxator/predictions/hg_scaffolds_alignment.tax']
    # contigPredArray = ['/Users/ivan/Documents/nobackup/megan/predictions/cr_chunks2000_contigs_m7-ex_1.txt', '/Users/ivan/Documents/nobackup/megan/predictions/hg_contigs_m7-ex_1.txt']
    # scaffPredArray = ['/Users/ivan/Documents/nobackup/megan/predictions/cr_chunks2000_scaffolds_m7-ex_1.txt', '/Users/ivan/Documents/nobackup/megan/predictions/hg_scaffolds_m7-ex_1.txt']
    contigFnaArray = ['/Users/ivan/Documents/nobackup/megan/fastaOrig/cr_chunks2000_contigs.fna', '/Users/ivan/Documents/nobackup/megan/fastaOrig/hg_contigs.fna']
    scaffFnaArray = ['/Users/ivan/Documents/nobackup/megan/fastaOrig/cr_chunks2000_scaffolds.fna', '/Users/ivan/Documents/nobackup/megan/fastaOrig/hg_scaffolds.fna']
    scaffContigMapArray = ['/Users/ivan/Documents/nobackup/megan/cowRumenChunked/chunks2000.groups', '/Users/ivan/Documents/nobackup/megan/hg/scafftocontig.txt']
    # outDir = '/Users/ivan/Documents/nobackup/megan/evaluation'
    outDir = '/Users/ivan/Documents/nobackup/taxator/evaluation'
    databaseFile = '/Users/ivan/Documents/work/binning/taxonomy/20121122/ncbitax_sqlite.db'
    taxonomicRanks = taxonomy_ncbi.TAXONOMIC_RANKS[1:]
    minSeqLen = 1000
    minScaffLen = 1000

    for i in range(2):

        d = os.path.join(outDir, os.path.basename('.'.join(contigFnaArray[i].split('.')[0:-1])))
        if not os.path.isdir(d):
            os.mkdir(d)

        cToPred = csv.predToDict(contigPredArray[i])

        # compute scaffold contig consistency
        if True:
            cons = consistency.Consistency(contigFnaArray[i], cToPred, scaffContigMapArray[i], databaseFile,
                                           minScaffContigCount=None, minScaffBpLen=None, cladesSet=None,
                                           considerContigWithNoScaff=True)

            out = csv.OutFileBuffer(os.path.join(d, 'consistency.txt'))
            out.writeText(cons.getScaffoldsPrint() + '\n')
            out.writeText(cons.getGroupedScaffoldsPrint() + '\n')
            out.close()
            cons.close()

        if True:
            # confusion tables: prediction of contigs vs predictions of scaffolds
            cmpDir = os.path.join(d, 'contigs_vs_scaff')
            if not os.path.isdir(cmpDir):
                os.mkdir(cmpDir)

            # create prediction file for contigs from a prediction file for scaffolds
            scaffPredAsContigs = os.path.join(d, 'scaffPredAsContigs.out')
            out_proc.scafToContigOutput(scaffContigMapArray[i], scaffPredArray[i], scaffPredAsContigs)

            ranks = taxonomy_ncbi.TAXONOMIC_RANKS[1:]
            cm = confusion_matrix.ConfusionMatrix(contigFnaArray[i], contigPredArray[i],
                                                  scaffPredAsContigs, databaseFile, ranks)
            for rank in ranks:
                cm.generateConfusionMatrix(rank, cmpDir)
            cm.close()

        # compute abundance profiles
        if True:
            abundanceContigDir = os.path.join(d, 'abundance_contig')
            abundanceScaffDir = os.path.join(d, 'abundance_scaff')
            if not os.path.isdir(abundanceContigDir):
                os.mkdir(abundanceContigDir)
            if not os.path.isdir(abundanceScaffDir):
                os.mkdir(abundanceScaffDir)

            # for contigs !!!
            run.createScaffoldAbundanceProfiles(contigFnaArray[i], contigPredArray[i],  abundanceContigDir, databaseFile,
                                                taxonomicRanks, minScaffLen)
            # for scaffolds !!!
            run.createScaffoldAbundanceProfiles(scaffFnaArray[i], scaffPredArray[i],  abundanceScaffDir, databaseFile,
                                                taxonomicRanks, minScaffLen)



if __name__ == "__main__":
    # _prepareFastaAndAlignments()

    # _filterAlignments()

    #_checkMagenPred('/Users/ivan/Documents/nobackup/megan/predictions')

    _evalSimDatasetBatch()

    _evalRealDataset()


# to get species to exclude:
# python /Users/ivan/Documents/work/python/PyCharm/PPSplus/algbioi/ref/mask_db.py
# -a cl -c clades.txt -d /Volumes/hera_net/metagenomics/projects/PPSmg/data/nobackup/NCBI20121122/sequences -r family
# -o . -t /Users/ivan/Documents/work/binning/taxonomy/20121122/ncbitax_sqlite.db