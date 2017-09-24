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

    Contains functionality that enables the use of the evaluation scripts for the CAMI challenge.
"""

import os
from algbioi.com import csv
from algbioi.com import fasta


def concatenate(directory, outputFile):
    out = csv.OutFileBuffer(outputFile)
    for f in os.listdir(directory):
        path = os.path.join(directory, f)
        name = f.split('.')[0]
        seqIdToSeq = fasta.fastaFileToDict(path)
        out.writeText('>' + str(name) + '\n')
        for seqId, seq in seqIdToSeq.iteritems():
            out.writeText(str(seq) + 200*'N' + '\n')
    out.close()


def outToCami(ppspOutFile):
    """
        Creates a cami output file, in format.
    """
    out = csv.OutFileBuffer(ppspOutFile + '.binning')
    out.writeText("""# The bioboxes.org binning output format: https://github.com/bioboxes/rfc/tree/master/data-format

@Version:0.9.1
@SampleID:
@TaxonomyID:

@@SEQUENCEID	TAXID

""")
    for line in open(ppspOutFile):
        name, taxonId = line.strip('\n').split('\t', 2)
        out.writeText("%s\t%s\n" % (name, taxonId))
    out.close()


def createCamiOut(fileDir):
    """
        Directory where the output files are located
    """
    outList = getPPSPOutPathList(fileDir)
    if len(outList) != 1:
        print('Unusual number of the output "pOUT" files detected: %s' % len(outList))
    for path in outList:
        outToCami(path)


def getPPSPOutPathList(outDir):
    retList = []
    for f in os.listdir(outDir):
        fList = f.split('.')
        if fList[-1] == 'pOUT' and fList[-2] != 'PP':
            retList.append(os.path.join(outDir, f))
    return retList


def readAssignments(assignmentFile):
    """
        Reads an assingment file, either in the cami format or in the PPS output (out) format

        @rtype: dict
        @return: mapping(name->taxonId)
    """
    if os.path.basename(assignmentFile).split('.')[-1] == 'binning':
        return readCami(assignmentFile)
    else:
        return csv.predToDict(assignmentFile)


def readCami(camiAssignFile):
    """
        Reads a file in the cami format

        @rtype: dict
    """
    ret = {}
    for line in open(camiAssignFile):
        line = line.strip()
        if not (line.startswith('#') or line.startswith('@') or len(line) == 0):
            name, taxonId = line.split('\t')[0:2]
            ret[name] = taxonId
        else:
            print line
    return ret
