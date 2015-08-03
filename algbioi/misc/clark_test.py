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


"""

import sys
import os

# from algbioi.com import taxonomy_ncbi
from algbioi.ref import create_mg_db
from algbioi.com import csv
from algbioi.com import fasta as fas


if sys.platform == 'darwin':
    FS_ROOT = '/Volumes/net'
else:
    FS_ROOT = '/net'

# for the uniform dataset:
TAXONOMY_FILE = os.path.join(FS_ROOT, 'metagenomics/projects/PPSmg/taxonomy/20121122/ncbitax_sqlite.db')
REF_MASK_GENUS = os.path.join(FS_ROOT, 'metagenomics/projects/PPSmg/tests/mercierVelvetMinimus050513/uniformSpeciesGeneral/config_rs_genus_mg_no/working/ref_seq_mask_genus')
REF_MASK_SPECIES = os.path.join(FS_ROOT, 'metagenomics/projects/PPSmg/tests/mercierVelvetMinimus050513/uniformSpeciesGeneral/config_rs_species_mg_no/working/ref_seq_mask_species')
REF_MASK_STRAIN = os.path.join(FS_ROOT, 'metagenomics/projects/PPSmg/tests/mercierVelvetMinimus050513/uniformSpeciesGeneral/config_rs_strain_mg_no/working/ref_seq_mask_strain')
REF_MASK_NO = os.path.join(FS_ROOT, 'metagenomics/projects/PPSmg/database/NCBI20121122/sequences')

WORKING_DIR = os.path.join(FS_ROOT, 'metagenomics/projects/PPSmg/tests/clark_test/working')

# /net/metagenomics/projects/PPSmg/data/mercier050513/velvet_minimus2/uniform/VelvetNoscafMinimus2_ma_merge_uniform_min_1000bp.fna



def toClark():
    # at which ranks do I predict using clark
    refDict = {}
    refDict['species'] = REF_MASK_STRAIN
    refDict['genus'] = REF_MASK_SPECIES
    refDict['family'] = REF_MASK_GENUS
    refDict['strain'] = REF_MASK_NO

    rankIdx = {}
    rankIdx['species'] = -2
    rankIdx['genus'] = -3
    rankIdx['family'] = -4
    rankIdx['strain'] = -2

    # ranks = ['species']
    # ranks = ['species', 'genus', 'family']
    ranks = ['strain']

    tax = create_mg_db.TaxonomyWrap(TAXONOMY_FILE)

    for rank in ranks:
        out = csv.OutFileBuffer(os.path.join(WORKING_DIR, 'targets_addresses_%s.txt' % rank))

        count = 0
        for f in os.listdir(refDict[rank]):
            fPath = os.path.join(refDict[rank], f)
            if os.path.isfile(fPath):
                taxonId = int(f.split('.', 1)[0])
                taxPathList = tax.getPathToRoot(taxonId).split(';')
                # print taxonId, taxPathList, taxPathList[rankIdx[rank]], fPath
                out.writeText('%s\t%s\n' % (fPath, taxPathList[rankIdx[rank]]))
                count += 1
                # if count > 2:
                #     break
        out.close()
        print('Entries for "%s" : %s' % (rank, count))

    tax.close()


def toKraken():
    targetsAddressDir = '/net/metagenomics/projects/PPSmg/tests/clark_test/working'
    dstOutDir = '/local/igregor/kraken_db_tmp'

    # ranks = ['genus', 'family']
    # ranks = ['species']
    ranks = ['strain']

    for rank in ranks:
        count = 0
        dstDir = os.path.join(dstOutDir, 'for_kraken_%s' % rank)
        if not os.path.isdir(dstDir):
            os.mkdir(dstDir)

        targetFile = os.path.join(targetsAddressDir, 'targets_addresses_%s.txt' % rank)

        for line in open(targetFile):
            line = line.strip()
            tokens = line.split('\t')
            if len(tokens) >= 2:
                fastaPath = tokens[0]
                # taxonId = int(tokens[1])
                taxonId = int(os.path.basename(fastaPath).split('.')[0])
                if not os.path.isfile(fastaPath):
                    print('File not exists: %s' % fastaPath)
                else:
                    out = csv.OutFileBuffer(os.path.join(dstDir, os.path.basename(fastaPath)))

                    for seqId, seq in fas.fastaFileToDictWholeNames(fastaPath).iteritems():
                        count += 1
                        out.writeText('>sequence%s|kraken:taxid|%s\t%s\n%s\n' % (count, taxonId, seqId, seq))

                    out.close()


if __name__ == "__main__":
    toClark()
    # toKraken()