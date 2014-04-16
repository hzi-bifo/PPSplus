"""
    To generate small simulated datasets from a large dataset

    source simset:
/Users/ivan/Documents/work/binning/data/mercier050513/melanieAssembly/uniform/VelvetNoscafMinimus2_ma_merge_uniform_min_1000bp.fna
/Users/ivan/Documents/work/binning/data/mercier050513/melanieAssembly/uniform/VelvetNoscafMinimus2_ma_merge_uniform_min_1000bp.tax

/Users/ivan/Documents/work/binning/data/mercier050513/soap_uniform.contig.profile.csv
    build dataset from:
"""


from algbioi.com import fasta as fas
from algbioi.com import csv


# build datasets for these ids, they are sorted according to the prevelance
PREVELANT_TAX_IDS = [243090, 226186, 316274, 257310, 435590, 198628, 324602, 379066, 240015, 246200, 203119, 226185]

IN_FASTA = '/Users/ivan/Documents/work/binning/data/mercier050513/melanieAssembly/uniform/VelvetNoscafMinimus2_ma_merge_uniform_min_1000bp.fna'
IN_TAX = '/Users/ivan/Documents/work/binning/data/mercier050513/melanieAssembly/uniform/VelvetNoscafMinimus2_ma_merge_uniform_min_1000bp.tax'

OUT_FASTA = '/Users/ivan/Documents/work/binning/data/mercier050513/test_subset/test_3strains.fna'
OUT_TAX = '/Users/ivan/Documents/work/binning/data/mercier050513/test_subset/test_3strains.tax'


def genSubset(taxonIdSet, inFasta, inTax, outFasta, outTax):
    """
        From the input files generate a subset output files that contain taxon ids given by the taxonIdSet.

        @param taxonIdSet:
        @type taxonIdSet: frozenset
        @param inFasta:
        @param inTax:
        @param outFasta:
        @param outTax:
    """
    outF = csv.OutFileBuffer(outFasta)
    outT = csv.OutFileBuffer(outTax)

    seqIdToSeq = fas.fastaFileToDictWholeNames(inFasta)
    seqIdToTaxonId = csv.getMapping(inTax, 0, 1, sep='\t')

    for seqId, seq in seqIdToSeq.iteritems():
        taxonId = int(seqIdToTaxonId[seqId][0])
        if taxonId in taxonIdSet:
            outF.writeText(">" + seqId + "\n" + seq + "\n")
            outT.writeText(seqId + "\t" + str(taxonId) + "\n")

    outF.close()
    outT.close()


if __name__ == "__main__":
    genSubset(frozenset(PREVELANT_TAX_IDS[:3]), IN_FASTA, IN_TAX, OUT_FASTA, OUT_TAX)
