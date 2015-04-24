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

    This module is meant for the CAMi evaluation container as a default task.
    It encapsulates calling of the following algbioi.eval scripts:
    accuracy.py
    consistency.py
    confusion_matrix.py
"""

import os
import sys
import argparse

from algbioi.com import fasta
from algbioi.com import csv
from algbioi.com import ncbitax2sqlite
from algbioi.eval import cami
from algbioi.eval import accuracy
from algbioi.eval import consistency
from algbioi.eval import confusion_matrix
from algbioi.com import taxonomy_ncbi

CORRECT_LABEL_THRESHOLD = 0.9
MIN_FRAC_CLADE = 0.01
MIN_FRAC_PRED = 0.01
RANKS = taxonomy_ncbi.TAXONOMIC_RANKS[1:]

class Args():
    def __init__(self, db, dmp):
        self.db = db
        self.dmp = dmp

def _main():
    # define arguments
    parser = argparse.ArgumentParser(description='Default task: PPS+ evaluation', epilog='')

    parser.add_argument('-b', '--cont-binning-file', nargs=1, type=file, required=True,
                        help='Binning file containing labels assigned to contigs.', metavar='assignments.csv',
                        dest='b')

    parser.add_argument('-t', '--cont-true-binning-file', nargs=1, type=file, required=True,
                        help='Binning file containing true labels for the contigs.', metavar='labels.csv', dest='t')

    parser.add_argument('-f', '--cont-contigs-file-listing', nargs=1, type=file, required=False,
                        help='A list of paths of FASTA contigs files.', metavar='fasta_listing.txt', dest='f')

    parser.add_argument('-m', '--cont-scaffold-contig-mapping', nargs=1, type=file, required=False,
                        help='Scaffold contig mapping, tab separated.', metavar='mapping.csv', dest='m')

    parser.add_argument('-n', '--cont-ncbi-taxonomy', nargs=1, required=False,
                        help='Directory containing the NCBI names.dmp and nodes.dmp files.', metavar='taxonomy_dir', dest='n')

    parser.add_argument('-o', '--cont-output-dir', nargs=1, required=True,
                        help='Output directory.', metavar='output_dir', dest='o')

    parser.add_argument('-j', '--default-job', nargs='+',
                        help='What task/job should be performed (p~precision/recall, s~scaff-contig consistency, '
                             'c~confusion tables, default if not spec compute all)', metavar='', dest='j')

    args = parser.parse_args()

    # read and check the arguments
    seqIdToBp = None
    scaffToContig = None
    binning = None
    trueBinning = None
    outputDir = None
    job = None

    if args.o and len(args.o) == 1 and os.path.isdir(args.o[0]):
        outputDir = args.o[0]

    if args.b and len(args.b) == 1 and os.path.isfile(args.b[0].name):
        binningFile = args.b[0].name
        binning = cami.readAssignments(binningFile)

    if args.t and len(args.t) == 1 and os.path.isfile(args.t[0].name):
        trueBinningFile = args.t[0].name
        trueBinning = cami.readAssignments(trueBinningFile)

    if args.f and len(args.f) == 1 and os.path.isfile(args.f[0].name) :
        contigsFileListing = args.f[0].name
        for line in open(contigsFileListing):
            if os.path.isfile(line.strip()):
                d = fasta.getSequenceToBpDict(line.strip())
                if seqIdToBp is None:
                    seqIdToBp = d
                else:
                    count = len(d) + len(seqIdToBp)
                    seqIdToBp.update(d)
                    if count > len(seqIdToBp):
                        sys.stderr.write('The fasta files contain duplicate entries!')

    if args.m and len(args.m) == 1 and os.path.isfile(args.m[0].name):
        scaffoldContigMapping = args.m[0].name
        scaffToContig = csv.getMapping(scaffoldContigMapping, 0, 1, '\t')

    taxonomyPath = os.path.join(outputDir, 'taxonomy_ncbi.db')
    if not os.path.isfile(taxonomyPath):
        if args.n and len(args.n) == 1 and os.path.isdir(args.n[0]):
            # build the ncbi taxonomy in the case it doesn't exist
            ncbitax2sqlite.build_database(Args(db=taxonomyPath, dmp=args.n[0]))
        else:
            taxonomyPath = None

    if args.j and len(args.j) > 0 and len(set(args.j).intersection(set(['p', 's', 'c']))) > 0:
        job = set(args.j)

    # print job
    # print args.j
    # print len(seqIdToBp)
    # print len(binning)
    # print len(trueBinning)
    # print taxonomyPath
    # print outputDir

    if (job is None or 'p' in args.j) and seqIdToBp and binning and trueBinning and taxonomyPath and outputDir:
        print('Computing precision/recall')
        # precision/recall - no correction
        acc = accuracy.Accuracy(seqIdToBp, binning, trueBinning, taxonomyPath)
        out = csv.OutFileBuffer(os.path.join(outputDir, 'precision_recall.csv'))
        out.writeText(acc.getAccuracyPrint(RANKS, MIN_FRAC_CLADE, MIN_FRAC_CLADE))
        out.close()
        acc.close()

        # precision/recall - with correction
        acc = accuracy.Accuracy(seqIdToBp, binning, trueBinning, taxonomyPath, CORRECT_LABEL_THRESHOLD)
        out = csv.OutFileBuffer(os.path.join(outputDir, 'precision_recall_correction.csv'))
        out.writeText(acc.getAccuracyPrint(RANKS, MIN_FRAC_CLADE, MIN_FRAC_CLADE))
        out.close()
        acc.close()

    # compute confusion matrices
    if (job is None or 'c' in args.j) and seqIdToBp and binning and trueBinning and taxonomyPath and outputDir:
        print('Computing confusion matrices')
        confusionMatrix = confusion_matrix.ConfusionMatrix(seqIdToBp, binning, trueBinning, taxonomyPath, RANKS)
        for rank in RANKS:
            confusionMatrix.generateConfusionMatrix(rank, os.path.join(outputDir, 'confusion_matrix'))
        confusionMatrix.close()

    # compute scaffold contig consistency
    if (job is None or 's' in args.j) and seqIdToBp and binning and scaffToContig and taxonomyPath \
            and outputDir:
        print('Computing scaffold-contig consistency')
        cons = consistency.Consistency(seqIdToBp, binning, scaffToContig, taxonomyPath)
        out = csv.OutFileBuffer(os.path.join(outputDir, 'consistency.txt'))
        out.writeText(cons.getGroupedScaffoldsPrint())
        cons.close()
        out.close()

    createEvalMetaFile(outputDir)


def createEvalMetaFile(outputDir):

    precisionRecallFile = os.path.join(outputDir, 'precision_recall.csv')
    precisionRecallCorrectionFile = os.path.join(outputDir, 'precision_recall_correction.csv')
    confusionMatrixDir = os.path.join(outputDir, 'confusion_matrix')
    consistencyFile = os.path.join(outputDir, 'consistency.txt')

    metaOut = csv.OutFileBuffer(os.path.join(outputDir, 'output.yaml'))
    # creates a metafile describing the results

    if os.path.isfile(precisionRecallFile):
        metaOut.writeText('''title: Precision and recall
type:
    name: csv
    options:
        header: "#"
        separator: ","
    data:
        reference: %s\n\n''' % precisionRecallFile)

    if os.path.isfile(precisionRecallCorrectionFile):
        metaOut.writeText('''title: Precision and recall with correction
type:
    name: csv
    options:
        header: "#"
        separator: ","
    data:
        reference: %s\n\n''' % precisionRecallCorrectionFile)

    if os.path.isfile(consistencyFile):
                metaOut.writeText('''title: Consistency
type:
    name: txt
    data:
        reference: %s\n\n''' % consistencyFile)

    if os.path.isdir(confusionMatrixDir):

        for f in os.listdir(confusionMatrixDir):

            filePath = os.path.join(confusionMatrixDir, f)
            rank = filePath.rsplit('.', 2)[1].split('_')[0]
            metaOut.writeText('''title: Confusion table for %s
type:
    name: csv
    options:
        header: ""
        separator: ","
    description:
        inline: Where rows correspond to the true assignments and columns correspond to the assignments by a binning method.
    data:
        reference: %s\n\n''' % (rank, filePath))

    metaOut.close()


if __name__ == "__main__":
    _main()
    # createEvalMetaFile('/Users/ivan/Documents/nobackup/tmp10')  # test yaml
