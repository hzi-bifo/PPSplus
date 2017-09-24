#!/usr/bin/env python

"""
    Copyright (C) 2015  Ivan Gregor

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

    File preprocession for the container integration.
"""

import os
from algbioi.com import fasta as fas
from algbioi.com import csv
import sys


def createConcatFasta(inFastaList, outFasta):

    out = csv.OutFileBuffer(outFasta)
    for line in open(inFastaList):
        line = line.strip()
        if line == '':
            continue
        else:
            for seqId, seq in fas.fastaFileToDictWholeNames(line).iteritems():
                out.writeText('>%s\n%s\n' % (seqId, seq))
    out.close()


def createConfig(configPath, outFasta, toolsRoot, refRoot, pipelineDir, configParam):

    # e.g. "processors=4;maxLeafClades=100;minPercentInLeaf=0.5"
    paramDict = dict(map(lambda x: tuple(x.split('=')), configParam.split(';')))

    out = csv.OutFileBuffer(configPath)

    out.writeText("""[PhyloPythiaS_Plus]
#
# Configuration file of the PhyloPythiaS Plus pipeline.
#
# Project pipeline directory for output and temporary data
""")
    out.writeText('pipelineDir=%s\n' % pipelineDir)
    out.writeText("""
#
#
# INPUT FILES
#
# Fasta file containing contigs/sequences (for the classification)
""")
    out.writeText('inputFastaFile=%s\n' % outFasta)
    out.writeText("""
#
# Fasta file containing scaffolds (optional)
inputFastaScaffoldsFile=
#
# Scaffold-contig mapping, tab separated file; map: scaffold -> contigs; (optional)
scaffoldsToContigsMapFile=
#
# Reference assignments file in PPS (*.out) format (optional)
referencePlacementFileOut=
#
#
# REFERENCE
#
# Directory that contains file named "ncbitax_sqlite.db" that represents the NCBI taxonomy in the sqlite3 format
""")
    out.writeText('databaseFile=%s/reference_NCBI201502/taxonomy\n' % refRoot)

    out.writeText("""
#
# Directory containing all reference sequences that can be used to train PhyloPythiaS
""")
    out.writeText('refSeq=%s/reference_NCBI201502/centroids_noplasmids\n' % refRoot)

    out.writeText("""
#
# Exclude reference sequences according to the reference assignments file at a given rank (optional)
# allowed ranks are: phylum, class, order, family, genus, species, strain
excludeRefSeqRank=
#
# Directory containing the "16S" marker gene reference databases (16S, 23S, 5S)
""")
    out.writeText('s16Database=%s/reference_NCBI201502/silva115\n' % refRoot)

    out.writeText("""
#
# Directory containing the "31" marker gene databases
""")
    out.writeText('mgDatabase=%s/reference_NCBI201502/mg5\n' % refRoot)

    out.writeText("""
#
# Exclude all reference marker gene sequences according to the reference assignments file at a given rank (optional)
# allowed ranks are: phylum, class, order, family, genus, species, strain
excludeRefMgRank=
#
#
# TOOLS
#
# Directory containing the PhyloPythiaS installation (supported version: 1.3)
""")
    out.writeText('ppsInstallDir=%s/PhyloPythiaS/1_3\n' % toolsRoot)
    out.writeText("""
#
# Directory containing the Mothur installation (supported version: Debian: v.1.33.3)
""")
    out.writeText('mothurInstallDir=%s/mothur/mothur1333\n' % toolsRoot)
    out.writeText("""
#
# Directory containing the HMMR 3.0 installation (supported version: 3.0)
""")
    out.writeText('hmmerBinDir=%s/hmmer-3.0/binaries\n' % toolsRoot)
    out.writeText("""
#
# Tool for the 16S analysis, from http://weizhong-lab.ucsd.edu/meta_rna/
""")
    out.writeText('rnaHmmInstallDir=%s/rna_hmm3\n' % toolsRoot)
    out.writeText("""
#
#
# BASIC SETTINGS
#
# Number of processors to be used
""")
    out.writeText('processors=%s\n' % paramDict.get('processors', '1'))
    out.writeText("""
#
# If False, use the new python scripts implementation (Default)
ruby=False
#
# Ncbi taxon ids only up to this rank, including, (as seen from the superkingdom) will be considered
# (e.g.: 1 ~ phylum, 2 ~ class, 3 ~ order, 4 ~ family, 5 ~ genus, 6 ~ species)
rankIdCut=6
#
# The maximum number of leaf clades (ncbi taxon ids) that will be modelled by PhyloPythiaS
# (max 100 recommended)
""")
    out.writeText('maxLeafClades=%s\n' % paramDict.get('maxLeafClades', '100'))
    out.writeText("""
#
# An ncbi taxon id will be considered if at least this percentage of all sample specific data (assigned to leafs
# of the taxonomy) was assigned to it.
""")
    out.writeText('minPercentInLeaf=%s\n' % paramDict.get('minPercentInLeaf', '1.0'))
    out.writeText("""
#
# Min. length of a sequence/contig/scaffold, shorter sequences won't be considered
minSeqLen=1000
#
# ADVANCED SETTINGS
#
# It is not recommended to change the advanced settings, however it is possible to give you more control.
#
# Assign not placed contigs of one scaffold according to the placed contigs of this scaffold,
# when sample specific data is generated based on the marker gene analysis (16S, 23S, 5S, and "31" genes).
# (see also: agThreshold, assignedPartThreshold)
placeContigsFromTheSameScaffold=True
#
# Agreement threshold
# Let us consider that N contigs of a scaffold are placed along a common path from the root R
# down to some low level clade LC. All contigs of the respective scaffold are placed to the lowest
# clade C (that lies on the path from R to LC), s.t., |agThreshold*N| of the contigs are placed
# on the path from C to LC. (default: 0.3-0.5?) (see also: placeContigsFromTheSameScaffold, assignedPartThreshold)
agThreshold=0.3
#
# Let us consider that the sum of all contigs' lengths of a scaffold is L and
# the sum of all assigned contigs' lengths of the scaffold is AL
# if AL/L >= assignedPartThreshold then all not assigned contigs of the scaffold are assigned according
# to the assigned contigs (see also: placeContigsFromTheSameScaffold, agThreshold)
assignedPartThreshold=0.5
#
# Up to this rank (including starting from 0) all ncbi taxon ids (found by the marker gene analysis) will be included
# if there is enough data to model them in the database (see: minGenomesWgs) or
# in the sample specific data (see: minBpToModel)
# (default: set it to 0~superkingdom or 1~phylum since there is often enough data in the sequence database to model them)
rankIdAll=0
#
#
weightStayAll=60.0
#
# Only ncbi taxon ids to which at least this number of bp of the sample specific data was
# assigned will be considered (except for ncbi taxon ids that belong to the ranks denoted by rankIdAll)
# (see rankIdCut)
rankIdCutMinBp=10000
#
# Override placements computed by the marker gene analysis by the output of PPS?
# (recommended: True)
overrideWithPPSPlacements=True
#
# A clade can be modeled if there is at least minGenomesWgs of genomes or draft genomes,
# from different species in the reference sequence database (see parameter: refSeq; see also: minBpToModel)
minBpToModel=100000
#
# Minimum size of a file with sample specific data for a clade (ncbi taxon id)
minSSDfileSize=5000
#
# Maximum size of a file with sample specific data
maxSSDfileSize=400000
#
# A clade can be modeled if there is at least minGenomesWgs of genomes or draft genomes
# from different species in the reference sequence database (see parameter: refSeq; see also: minBpToModel)
minGenomesWgs=1
#
# If there is at least this number of bp of reference sequences in the reference directory (see: refSeq)
# for one species, it is considered to be a sufficient amount of reference data for this species and we don't require
# that there exists a reference genome for this species
minBpPerSpecies=300000
#
#
# A sequence will be assigned according to the top percent candidate assignments
# a candidate assignment will be considered if its weight is >= highestWeight*(1 - threshold)
candidatePlTopPercentThreshold=0.1
#
# The contigs' names in the reference file must correspond to the contigs' names in the input fasta files
# if this is not the case, you can enter a regexp. where the string in brackets denotes the part of the contig's name
# that will be copied to the (*.PP.pOUT) and (*.pOUT) pipeline output files
# e.g. the contigs' names in the fasta input file are: scaffold[0-9]+_contig[0-9]+
# and the contigs' names in the reference file are: contig[0-9]+
# to mach the names, you can enter: ^scaffold[0-9]+_(contig[0-9]+)
# e.g. ^(Scaffold_[0-9]+_[0-9]+[^0-9].*)
outputFileContigSubPattern=^(.*)
#
# Comma separated parameters of the mothur classify command
mothurClassifyParamOther=method=bayesian, cutoff=80, iters=300
#
# PhyloPythiaS configuration file (optional)
# Default configuration file will be generated, if left empty (recommended)
# Please make sure that all paths in this configuration file are correct.
# You can, e.g first run the pipeline without defining this configuration file, then a default configuration file will
# be generated (you can find it in: pipelineDir\working\PPS_config_generated.txt). Then, you can rename and modify this
# PPS configuration file and then use it here.
configPPS=
#
# A noise removal parameter for the recall computation. A clade from the sample is considered only if
# the dataset contains at least this fraction of sequences of a respective clade.
recallMinFracClade=0.001
#
# A noise removal parameter for the precision computation. A predicted clade is considered only if at least
# this fraction of sequences, out of all sequences predicted at a particular rank, is assigned to this clade.
precisionMinFracPred=0.001
#
# If all sequences from one clade are consistently predicted to a wrong clade, the wrong labels are corrected to
# the right ones. If we wanted to reconstruct a draft genome, it's more important that the sequences are predicted
# consistently, even though they are predicted to a wrong clade.
# This is done only if (correctLabelThreshold * 100)% bp are predicted in such a wrong way.
correctlabelthreshold=0.9
""""")

    out.close()


def _main():
    if len(sys.argv) < 3:
        print('Parameters:\ninFastaList\n outFasta\n configPath\n toolsRoot\n refRoot\n pipelineDir\nconfigParam\n')

    inFastaList, outFasta, configPath, toolsRoot, refRoot, pipelineDir, configParam = sys.argv[1:8]

    assert os.path.isfile(inFastaList)
    assert os.path.isdir(os.path.dirname(outFasta))
    assert os.path.isdir(os.path.dirname(configPath))
    assert os.path.isdir(toolsRoot)
    assert os.path.isdir(refRoot)
    assert os.path.isdir(pipelineDir)

    if not os.path.isfile(outFasta):
        createConcatFasta(inFastaList, outFasta)

    if not os.path.isfile(configPath):
        createConfig(configPath, outFasta, toolsRoot, refRoot, pipelineDir, configParam)

if __name__ == "__main__":
    _main()