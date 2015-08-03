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

    To wrap the pps python scripts.
"""

import os
import sys
from algbioi.com import csv
from algbioi.com import parallel

import multiprocessing as mp


def runTrain(ppsInstallDir, ppsConfigFilePathPy):
    cwd = os.path.join(ppsInstallDir, 'python_scripts')
    scriptPath = os.path.join(cwd, 'train.py')
    cmd = 'python %s -c %s -y' % (scriptPath, ppsConfigFilePathPy)
    if parallel.reportFailedCmd(parallel.runCmdSerial([parallel.TaskCmd(cmd, cwd)])) is not None:
        sys.exit(-1)


def runPredict(ppsInstallDir, ppsConfigFilePathPy, fastaFile):
    cwd = os.path.join(ppsInstallDir, 'python_scripts')
    scriptPath = os.path.join(cwd, 'predict.py')
    cmd = 'python %s -c %s -fasta %s ' % (scriptPath, ppsConfigFilePathPy, fastaFile)

    baseName = os.path.basename(fastaFile)

    for f in os.listdir(os.path.dirname(fastaFile)):
        if (f.endswith('.out') or f.endswith('.sl')) and baseName in f:
            os.remove(os.path.join(os.path.dirname(fastaFile), f))

    if parallel.reportFailedCmd(parallel.runCmdSerial([parallel.TaskCmd(cmd, cwd)])) is not None:
        sys.exit(-1)

    for f in os.listdir(os.path.dirname(fastaFile)):
        if f.endswith('.sl') and baseName in f:
            os.remove(os.path.join(os.path.dirname(fastaFile), f))


def createPPSConfigPython(ppsConfigFilePath, keyDict):
    """
        Generates the up to date PPS configuration file

        @param keyDict: must specify: NCBI_PROCESSED_DIR, NCBI_TAX_DIR, PROJECT_DIR, TREE_FILE, SAMPLE_SPECIFIC_DIR,
            can be specified: PROCESSORS, GENOMES_EXCLUDE
        @type keyDict: dict
    """
    out = csv.OutFileBuffer(ppsConfigFilePath)
    out.writeText(
        """#######################################################################
#configuration file (PPS+ GENERATED !!!)
#please make sure that there is no space before or after ":"
#lines starting with character "#" are treated as comments
#please provide complete paths instead of only file or directory names
#######################################################################
#directory where processed NCBI data is stored, provide empty directory to create new
#REUSABLE\n""")
    out.writeText('NCBI_PROCESSED_DIR:%s\n' % keyDict.get('NCBI_PROCESSED_DIR', ''))
    out.writeText(
        """#Directory containing NCBI taxonomy in SQlite3 format with file name "ncbitax_sqlite.db"
#provide empty directory to create new database
#REUSABLE\n""")
    out.writeText('NCBI_TAX_DIR:%s\n' % keyDict.get('NCBI_TAX_DIR', ''))
    out.writeText('#project directory, the directory must be empty\n')
    out.writeText('PROJECT_DIR:%s\n' % keyDict.get('PROJECT_DIR', ''))
    out.writeText(
        """#############################
#!!!FOLLOWING ARE OPTIONAL!!!
#############################
#The name of the SQLite3 database (leave it empty to use default name "ncbitax_sqlite.db"
NCBI_SQL:
###### Output space options #####
#a file containing a tree in newick format (see restrictions in INSTALL.txt)
#OR a file with ncbi taxon ids (one id per line) to create a tree from\n""")
    out.writeText('TREE_FILE:%s\n' % keyDict.get('TREE_FILE', ''))
    out.writeText(
        """#Taxonomic ranks (comma separated, no space) starting at the lowest rank. \
Please make sure that "root" is there at the end.
TAXONOMY_RANKS:species,genus,family,order,class,phylum,superkingdom,root
#number of minimum genomes a clade must contain to be included in generic model
#effective only if tree file is not provided
N_MIN_GENOMES_GENERIC:3
#action on loss 0:disabled, 1:invert
LOSS_ACTION:0
###### Input space options #####
#a directory with sample specific fasta files (file names must start with appropriate organism/species \
ncbi taxonomic id)
#leave empty if you don't have any\n""")
    out.writeText('SAMPLE_SPECIFIC_DIR:%s\n' % keyDict.get('SAMPLE_SPECIFIC_DIR', ''))
    out.writeText(
        """#kmer feature space for multiple kmers use kmer_min-kmer_max
KMER:4-6
#Fragment lengths for different models (comma separated, no space)
FRAGMENT_LEN:1000,3000,5000,10000,15000,50000
#kmer feature
#use reverse complement for computing kmer features?
REV_COMPLEMENT:1
#remove reverse complement features?
RM_REV_COMPLEMENT:1
#0:disabled, 1:sequence length, 2:sequence_length-k+1, 3:embedded monomer frequency
KMER_NORMALIZATION:1
#Number of examples per training file
NUMBER_EXAMPLES:10000
#step size for sample specific data; either a single number (for all fragment lengths) or an array separated with ","
SAMPLE_SPECIFIC_STEP:1000,300,500,1000,1500,5000
###### Training options #####
#C values for SVM, if single value is given then models will be build with that value.
#If comma separated (no space) values are given then cross-validation will be performed.
#If a single value is provided, all models will be built with it. Our experience shows that in general
#values less than 1 (e.g. 0.01 and 0.1) do not provide good models.
C_GRID:1000
#clean-up the data (sampled_fasta and train_data directories) created after training? TRUE/FALSE
CLEAN_UP_TRAIN:FALSE
#kernel type 0:linear, 1:polynomial, 2:rbf (on-linear kernels are computationally expensive)
KERNEL:0
##polynomial kernel degree
KERNEL_POLYNOMIAL_DEGREE:2
##rbf kernel gamma
KERNEL_RBF_GAMMA:1
##polynomial kernel s
KERNEL_POLYNOMIAL_S:1
######  Predictions options #####
#number of classifiers to use, keep this odd to avoid ties
N_CLASSIFIERS:3
#Create output in PhyloPhitia format TRUE/FALSE (in prediction)
POSTPROCESSING:FALSE
###### Misc options #####
#should the models be built in parallel (please make sure that you have enough number of
processors and main memory), number of CPUs to use\n""")
    out.writeText('PROCESSORS:%s\n' % keyDict.get('PROCESSORS', mp.cpu_count()))
    out.writeText(
        """#allowed file extensions
EXTENSIONS:fna,fasta,fas,fa
#genomes to exclude: file containing one ncbi tax_id per line\n""")
    out.writeText('GENOMES_EXCLUDE:%s\n' % keyDict.get('GENOMES_EXCLUDE', ''))
    out.writeText(
        """#if the training data is already there then just build models (TRUE/FALSE)
ONLY_MODELS:FALSE\n""")
    out.writeText("BALANCE_CLASSES:FALSE\n")
    out.close()


if __name__ == "__main__":
    pass