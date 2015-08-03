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

    To create the mg databases from the gbk files.
"""

import os
# import sys
# import signal
# import string
# from Bio import SeqIO
from algbioi.com import csv
from algbioi.com import gbk
from algbioi.com import parallel
from algbioi.com import taxonomy_ncbi
import multiprocessing as mp
from algbioi.com import fasta as fas
from algbioi.ref import create_mg_db

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


def extractMgForTaxon(taxonId, fasPath, gbkFileList, taxFile, mgNameSet, outFile):
    try:
        # create taxonomy
        # tax = taxonomy_ncbi.TaxonomyNcbi(taxFile)
        tax = None

        # read in sequences
        seqIdToSeq = fas.fastaFileToDictWholeNames(fasPath)

        # output FASTA entries
        outLineList = []

        # for each gbk file
        for gbkFile in gbkFileList:

            # is gbk valid?
            last = ''
            for line in open(gbkFile):  # TODO: add this as a function to algbioi.com.gbk
                line = line.strip()
                if len(line) > 0:
                    last = line
            if '//' != last:
                print('Unusual end of file: |%s| In file: |%s|, file skipped!' % (last, gbkFile))
                continue

            # read in a gbk file
            recList = gbk.readFromGbkFile(gbkFile, tax, ('acc_version', 'genes_annotation', 'description'), considerDnaSeq=False)

            # for each record
            for rec in recList:

                # get gene annotations
                # taxonId2 = rec.get('taxonId')
                accVersion = rec.get('acc_version')
                geneDict = rec.get('genes_annotation')
                description = rec.get('description')

                if description is None:
                    print('No description for: %s %s' % (taxonId, accVersion))
                    continue
                else:
                    if "plasmid" in str(description).lower():
                        continue

                # check tax id
                # if taxonId2 is not None and int(taxonId) != int(taxonId2):
                #     print('Taxon id mismatch: %s %s %s' % (taxonId, taxonId2, accVersion))

                if geneDict is not None and len(geneDict) > 0 and accVersion is not None:

                    # get the reference sequence
                    seq = seqIdToSeq.get(accVersion)
                    if seq is None:
                        print('Seq not found: %s %s' % (accVersion, taxonId))
                    else:
                        # for all gene-records
                        for geneName, gr in geneDict.iteritems():
                            if geneName in mgNameSet:
                                seqDna = Seq(gbk.getRegions(seq, gr.locationList), generic_dna)
                                try:
                                    protSeq = seqDna.translate(table=11, stop_symbol='', cds=True)
                                except Exception as e:
                                    protSeq = None
                                    print('Error when translating: %s %s %s %s' % (taxonId, accVersion, geneName, type(e)))
                                if str(gr.seqPROT) == str(protSeq):
                                    # print taxonId, gr.geneName, gr.sequenceAccessionVersion
                                    outLineList.append('>g|%s|%s|ncbid:%s\n%s\n' % (gr.geneName, accVersion, taxonId, str(seqDna).lower()))

                                    if gr.sequenceAccessionVersion != accVersion:
                                        print('Acc versions do not match: %s %s %s' % (gr.sequenceAccessionVersion, accVersion, taxonId))
                                else:
                                    print('Prot seq. do not match: %s %s %s %s %s' % (accVersion, geneName, gr.seqPROT, protSeq, taxonId))
                elif accVersion is None:
                    print('Acc version is none! %s' % taxonId)

        if len(outLineList) > 0:
            out = csv.OutFileBuffer(outFile)
            out.writeText(''.join(outLineList))
            out.close()

        # tax.close()
    except Exception as e:
        print e.message, type(e), taxonId, fasPath, gbkFileList, taxFile, mgNameSet, outFile


def extractMg(gbkDir, fasDir, taxFile, mgNameFile, outDir):
    assert os.path.isdir(gbkDir)
    assert os.path.isdir(fasDir)
    assert os.path.isfile(taxFile)
    assert os.path.isfile(mgNameFile)
    assert os.path.isdir(outDir)

    # read in MG-names
    mgNameSet = set()
    for name in open(mgNameFile):
        mgNameSet.add(name.strip())

    # map: taxonId -> fasta-file-path
    fasDict = {}
    for f in os.listdir(fasDir):
        fPath = os.path.join(fasDir, f)
        if os.path.isfile(fPath):
            taxonId = int(f.split('.')[0])
            fasDict[taxonId] = fPath

    # map: taxonId -> list of gbk files
    gbkDict = {}
    for f in os.listdir(gbkDir):
        fPath = os.path.join(gbkDir, f)
        if os.path.isfile(fPath):
            try:
                taxonId = int(f.split('_')[1].split('.')[0])
            except Exception as e:
                print('Exception: %s %s' % (type(e), f))
                raise e
            if taxonId not in gbkDict:
                gbkDict[taxonId] = [fPath]
            else:
                gbkDict[taxonId].append(fPath)
    taskList = []
    for taxonId, fasPath in fasDict.iteritems():

        if taxonId in gbkDict:
            gbkFileList = gbkDict[taxonId]
            outFile = os.path.join(outDir, str(taxonId) + '_mg.fna')

            taskList.append(parallel.TaskThread(extractMgForTaxon, (taxonId, fasPath, gbkFileList, taxFile, set(mgNameSet), outFile)))

        else:
            print('No gbk file for taxonId: %s' % taxonId)

            # taskList.append(parallel.TaskThread(extractPlasmidAccession2Part, (fPath, taxFile, taxonId)))

    parallel.runThreadParallel(taskList, mp.cpu_count(), keepRetValues=False)


def createMgDb(mgNameFile, taxFile, srcDir, outDir):
    assert os.path.isfile(mgNameFile)
    assert os.path.isfile(taxFile)
    assert os.path.isdir(srcDir)
    assert os.path.isdir(outDir)

    tax = create_mg_db.TaxonomyWrap(taxFile)

    # read in MG-names
    mgNameSet = set()
    for name in open(mgNameFile):
        mgNameSet.add(name.strip())

    # map: geneName -> out-fasta-file
    outFas = {}

    # map: geneName -> out-tax-file
    outTax = {}

    # open files for writing
    for geneName in mgNameSet:
        outFas[geneName] = csv.OutFileBuffer(os.path.join(outDir, geneName + '_bact+arch_dna.fna'))
        outTax[geneName] = csv.OutFileBuffer(os.path.join(outDir, geneName + '_bact+arch_dna.tax'))

    seqNameSet = set()
    seqSet = set()

    for f in os.listdir(srcDir):
        fPath = os.path.join(srcDir, f)
        if os.path.isfile(fPath):
            for seqId, seq in fas.fastaFileToDictWholeNames(fPath).iteritems():
                tokens = seqId.split('|')
                geneName = tokens[1]
                taxonId = int(tokens[3].split(':')[1])
                taxonPath = tax.getPathToRoot(taxonId)

                if seqId in seqNameSet:
                    print('SeqId already exist! %s' % seqId)
                    continue

                if seq in seqSet:
                    print('Seq already exist for: %s' % geneName)
                    continue

                outFas[geneName].writeText('>%s\n%s\n' % (seqId, seq))
                outTax[geneName].writeText('%s\t%s\n' % (seqId, taxonPath))
                seqNameSet.add(seqId)
                seqSet.add(seq)

    for geneName in mgNameSet:
        outFas[geneName].close()
        outTax[geneName].close()
    tax.close()
    # print mgNameFile, taxFile, srcDir, outDir


def _main():

    # gbkDir = '/home/igregor/Documents/work/ppsp/bac_arch/gbk3'
    # fasDir = '/home/igregor/Documents/work/ppsp/bac_arch/bac_arch_centroids'
    # taxFile = '/home/igregor/Documents/work/taxonomy/NCBI201502/ncbitax_sqlite.db'
    # mgNameFile = '/home/igregor/Documents/work/ppsp/bac_arch/mg_31.txt'
    # outDir = '/home/igregor/Documents/work/ppsp/bac_arch/mg_gbk'
    #
    # outDir2 = '/home/igregor/Documents/work/ppsp/bac_arch/mg_n1/db'

    # gbkDir = '/local/igregor/ref_cami_1/bac_arch/gbk2'
    # fasDir = '/local/igregor/ref_cami_1/bac_arch_centroids'
    taxFile = '/net/metagenomics/projects/PPSmg/database/NCBI201502/nobackup/mg5/taxonomy/ncbitax_sqlite.db'
    mgNameFile = '/local/igregor/ref_cami_1/bac_arch/mg_31.txt'
    outDir = '/local/igregor/ref_cami_1/bac_arch/mg_gbk'
    outDir2 = '/local/igregor/ref_cami_1/bac_arch/mg_n1/db'

    # extractMg(gbkDir, fasDir, taxFile, mgNameFile, outDir)
    #
    createMgDb(mgNameFile, taxFile, outDir, outDir2)


def _test1():
    taxFile = '/home/igregor/Documents/work/taxonomy/NCBI201502/ncbitax_sqlite.db'
    mgNameFile = '/home/igregor/Documents/work/ppsp/bac_arch/mg_31.txt'
    mgOrig = '/home/igregor/Documents/work/ppsp/mg/mg4/db'
    mgGbk = '/home/igregor/Documents/work/ppsp/mg/mg_n1/db'
    mgDst = '/home/igregor/Documents/work/ppsp/mg/mg5/db'
    refDir = '/media/igregor/verbatim/work/ppsp_ref/ref_cami1/bac_arch/bac_arch_no_plas'

    # init taxonomies
    tax = taxonomy_ncbi.TaxonomyNcbi(taxFile)
    taxW = create_mg_db.TaxonomyWrap(taxFile)

    mapping = getMapping()

    # for all genes
    for line in open(mgNameFile):
        gene = line.strip()

        fnaFile = '%s_bact+arch_dna.fna' % gene
        taxFile = '%s_bact+arch_dna.tax' % gene

        outFna = csv.OutFileBuffer(os.path.join(mgDst, fnaFile))
        outTax = csv.OutFileBuffer(os.path.join(mgDst, taxFile))

        # buffer mgGbk
        seqToEntryGbk = {}
        seqAllSet = set()
        for seqId, seq in fas.fastaFileToDictWholeNames(os.path.join(mgGbk, fnaFile)).iteritems():

            taxonId = int(seqId.split('|')[-1].split(':')[1])
            taxPath = taxW.getPathToRoot(taxonId)
            taxPathL = map(lambda x: int(x), taxPath.split(';')[:-1])
            assert seq not in seqToEntryGbk
            seqToEntryGbk[seq] = (seqId, taxonId, taxPathL, taxPath)
            seqAllSet.add(seq)

        # buffer orig mapping
        seqIdToTaxOrig = {}
        for entry in open(os.path.join(mgOrig, taxFile)):
            entry = entry.strip().split('\t')
            seqId = entry[0]
            taxonId = int(seqId.split('|')[-1].split(':')[1])
            taxonIdN = taxonId
            taxPath = None
            taxPathL = None

            if tax.exists(taxonId):
                taxPath = taxW.getPathToRoot(taxonId)
                taxPathL = map(lambda x: int(x), taxPath.split(';')[:-1])

            else:
                taxonId = int(map(lambda x: int(x), entry[1].split(';')[:-1])[-1])
                if tax.exists(taxonId):
                    taxPath = taxW.getPathToRoot(taxonId)
                    taxPathL = map(lambda x: int(x), taxPath.split(';')[:-1])
                else:
                    taxonId = mapping.get(taxonIdN)
                    print 'mapping:%s->%s\t%s\t%s' % (taxonIdN, taxonId, seqId, tax.exists(taxonId))
                    if taxonId is not None and tax.exists(taxonId):
                        taxPath = taxW.getPathToRoot(taxonId)
                        taxPathL = map(lambda x: int(x), taxPath.split(';')[:-1])
                    else:
                        print 'mapping from higher rank taken', seqId
                        taxonId = int(map(lambda x: int(x), entry[1].split(';')[:-1])[-2])
                        if tax.exists(taxonId):
                            taxPath = taxW.getPathToRoot(taxonId)
                            taxPathL = None

            assert seqId not in seqIdToTaxOrig
            seqIdToTaxOrig[seqId] = (taxonId, taxPathL, taxPath)

        # buffer orig sequences incl. mapping
        seqToEntryOrig = {}
        for seqId, seq in fas.fastaFileToDictWholeNames(os.path.join(mgOrig, fnaFile)).iteritems():
            taxonId, taxPathL, taxPath = seqIdToTaxOrig[seqId]
            entry = (seqId, taxonId, taxPathL, taxPath)
            if seq in seqToEntryOrig:
                seqToEntryOrig[seq].append(entry)
            else:
                seqToEntryOrig[seq] = [entry]
            seqAllSet.add(seq)

        # store all sequences
        for seq in seqAllSet:

            entryList = []

            if seq in seqToEntryGbk:
                entryList.append(seqToEntryGbk[seq])

            if seq in seqToEntryOrig:
                entryList += seqToEntryOrig[seq]

            assert len(entryList) > 0

            if len(entryList) == 1:
                entry = entryList[0]
            else:
                entry = entryList[0]
                maxSize = 0
                for e in entryList:
                    seqId, taxonId, taxPathL, taxPath = e
                    if not taxPathL is None:
                        size = getTotalSize(tax, taxPathL[-1], refDir)
                        if size > maxSize:
                            maxSize = size
                            entry = e
                seqId, taxonId, taxPathL, taxPath = entry
                assert taxPathL is not None  # TODO:

            seqId, taxonId, taxPathL, taxPath = entry
            if seqId not in ['gi|556503834:443051-443641-|gid:945066|ncbid:511145']:  # TODO: skip ambiguous sequences
                outFna.writeText('>%s\n%s\n' % (seqId, seq))
                outTax.writeText('%s\t%s\n' % (seqId, taxPath))

        outFna.close()
        outTax.close()
    print('done...')

    tax.close()
    taxW.close()


def checkUnique():
    mgNameFile = '/home/igregor/Documents/work/ppsp/bac_arch/mg_31.txt'
    mgDst = '/home/igregor/Documents/work/ppsp/mg/mg5/db'
    nameToGene = {}
    seqToGene = {}
    c = 0
    error = 0
    error2 = 0
    for line in open(mgNameFile):
        gene = line.strip()

        fnaFile = '%s_bact+arch_dna.fna' % gene
        taxFile = '%s_bact+arch_dna.tax' % gene

        fnaPath = os.path.join(mgDst, fnaFile)
        taxPath = os.path.join(mgDst, taxFile)


        # for line in open(taxPath):
        #     line = line.strip()
        #     if len(line) > 0:
        #         c += 1
        #         if line in nameToGene:
        #             print gene, nameToGene[line], line
        #             error += 1
        #         else:
        #             nameToGene[line] = gene

        for seqId, seq in fas.fastaFileToDictWholeNames(fnaPath).iteritems():
            if seq in seqToGene:
                print gene, seqToGene[seq], seqId
                error2 += 1
            else:
                seqToGene[seq] = gene

    print c, error, error2


def getSuccessors(tax, taxonId):
    """
        Get the id of all successors including itself.

        @type tax: taxonomy_ncbi.TaxonomyNcbi
        @rtype: list[int]
    """
    rList = [taxonId]
    childList = tax.getChildrenNcbids(taxonId)

    if childList is not None:

        for child in childList:
            rList += getSuccessors(tax, child)

    return rList


def getTotalSize(tax, taxonId, refDir):

    taxIdList = getSuccessors(tax, taxonId)

    sumBp = 0

    for taxId in taxIdList:

        fPath = os.path.join(refDir, '%s.1.fna' % taxId)
        if os.path.isfile(fPath):
            sumBp += os.path.getsize(fPath)

    return round(float(sumBp) / 1000000., 2)


def getBp():

    tax = taxonomy_ncbi.TaxonomyNcbi('/home/igregor/Documents/work/taxonomy/NCBI201502/ncbitax_sqlite.db')

    taxonId = 29461
    # taxonId = 437687

    # print len(getBpPart(tax, taxonId))
    c = 0
    for i in getSuccessors(tax, taxonId):
        c += 1
        print tax.getScientificName(i)
    # print len(tax.getChildrenNcbids(taxonId))
    print c
    tax.close()


def getMapping():
    m = {}

    # tax = taxonomy_ncbi.TaxonomyNcbi('/home/igregor/Documents/work/taxonomy/NCBI201502/ncbitax_sqlite.db')

    m[1125778] = 936595  # Lachnoanaerobaculum sp. OBRC5-5 True
    m[596373] = 545619  # Paraoerskovia marina
    m[915168] = 915173  # Faecalicoccus acidiformans
    m[515898] = 515897  # Pontibaca methylaminivorans
    m[268392] = 189966  # Rhizomicrobium palustre
    m[914154] = 849756  # Mongoliitalea lutea
    m[429433] = 652451  # uncultured Pelosinus sp.
    m[522600] = 1095917  # uncultured Gemmobacter sp.
    m[1219527] = 1216886  # Nguyenibacter vanlangensis
    m[469546] = 469547  # Aquamicrobium segne
    m[469550] = 469551  # Aquamicrobium ahrensii
    m[469548] = 469547  # Aquamicrobium segne
    m[287484] = 354246  # Spongiibacter marinus
    m[553320] = 741655  # uncultured Propionispira sp.
    m[89785] = 1605449  # Alicyclobacillus genomosp. 1
    m[868132] = 877455  # Methanobacterium lacus
    m[204471] = 1605449  # Alicyclobacillus genomosp. 1
    m[157697] = 1605449  # Alicyclobacillus genomosp. 1
    m[645276] = 645275  # Murinocardiopsis flavida
    m[564202] = 1795  # Mycobacterium neoaurum
    m[564200] = 1795  # Mycobacterium neoaurum
    m[90958] = 180  # Leptospirillum ferrooxidans
    m[1092647] = 294  # Pseudomonas fluorescens
    m[103818] = 83526  # actobacillus paralimentarius
    m[1004270] = 1004279  # Mesorhizobium muleiense
    m[1141658] = 1141657  # Nocardia vulneris
    m[434269] = 407234  # Rhodobacter vinaykumarii
    m[1141656] = 1141657  # Nocardia vulneris
    m[717855] = 717859  # Agromyces subtropicus
    m[564199] = 37916  # Mycobacterium chlorophenolicum
    m[1141659] = 1141657  # Nocardia vulneris
    m[930403] = 294  # Pseudomonas fluorescens
    m[980274] = 445327  # Rhodopirellula lusitana
    m[980287] = 445327  # Rhodopirellula lusitana
    m[1226328] = 1226327  # Acinetobacter kookii
    m[759840] = 759839  # Rhodopseudomonas harwoodiae
    m[796358] = 701526  # Streptococcus porcorum
    m[1179768] = 1179767  # Bradyrhizobium ganzhouense
    m[1108596] = 1108595  # Chromobacterium vaccinii
    m[657020] = 1422  # Geobacillus stearothermophilus
    m[81855] = 89059  # Lactobacillus acidipiscis
    m[472172] = 28254  # Marinomonas communis
    m[520090] = 520089  # Nocardioides hungaricus
    m[566319] = 470  # Acinetobacter baumannii
    m[657021] = 33941  # Geobacillus thermoleovorans
    m[1128627] = 1406  # Paenibacillus polymyxa
    m[564204] = 1795  # Mycobacterium neoaurum
    m[335018] = 56731  # Rhizobium giardinii
    m[564201] = 1795  # Mycobacterium neoaurum
    m[348777] = 1437882  # Pseudomonas nitroreducens HBP1
    m[1141653] = 1141657  # Nocardia vulneris
    m[363856] = 170367  # Paenibacillus anaericanus
    m[980268] = 980271  # Rhodopirellula rubra
    m[435962] = 1218493  # Lactobacillus kullabergensis
    m[404999] = 379347  # Ruegeria mobilis
    m[980264] = 980271  # Rhodopirellula rubra
    m[1113630] = 1203393  # Komagataeibacter maltaceti
    m[520837] = 1214071  # Ensifer sesbaniae
    m[564195] = 188917  # Mycobacterium petroleophilum
    m[288838] = 240427  # Lactobacillus paracollinoides
    m[1166682] = 1166683  # Rhizobium cauense
    m[980273] = 445327  # Rhodopirellula lusitana
    m[1271773] = 1609977  # Sphingomonas sp. WHSC-8
    m[1141654] = 1141657  # Nocardia vulneris
    m[1273082] = 1402  # Bacillus licheniformis
    m[652786] = 354630  # Mucilaginibacter lappiensis
    m[650063] = 69670  # Brevundimonas halotolerans
    m[980270] = 980271  # Rhodopirellula rubra
    m[413968] = 2242  # Halobacterium salinarum
    m[564197] = 1795  # Mycobacterium neoaurum
    m[980275] = 445327  # Rhodopirellula lusitana
    m[928851] = 106590  # Cupriavidus necator
    m[593520] = 392408  # Deinococcus depolymerans
    m[1247515] = 1247514  # Luteimonas abyssi
    m[980277] = 445327  # Rhodopirellula lusitana
    m[334851] = 60552  # Burkholderia vietnamiensis
    m[650065] = 69670  # Brevundimonas halotolerans
    m[980281] = 445327  # Rhodopirellula lusitana
    m[531950] = 1437453  # Streptomyces leeuwenhoekii
    m[435963] = 1218494  # Lactobacillus helsingborgensis
    m[1141657] = 1141657  # Nocardia vulneris
    m[1141652] = 1141657  # Nocardia vulneris
    m[438742] = 1930  # Streptomyces scabiei
    m[230092] = 47856  # Micromonospora coerulea
    m[796357] = 701526  # Streptococcus porcorum
    m[1247516] = 1247514  # Luteimonas abyssi
    m[1151581] = 1562970  # Porphyromonadaceae bacterium ING2-E5B
    m[118061] = 118060  # Enterococcus rotai
    m[877501] = 877500  # Arcobacter anaerophilus
    m[133922] = 1236974  # Paenibacillus sp. JCM 10914
    m[564196] = 258505  # Mycobacterium fluoranthenivorans
    m[980278] = 445327  # Rhodopirellula lusitana
    m[230093] = 47856  # Micromonospora coerulea
    m[457675] = 554117  # Paenibacillus phoenicis
    m[1214926] = 1214925  # Streptomyces chumphonensis
    m[564203] = 1795  # Mycobacterium neoaurum
    m[1156391] = 1406  # Paenibacillus polymyxa
    m[980272] = 980271  # Rhodopirellula rubra
    m[520088] = 520089  # Nocardioides hungaricus
    m[931073] = 931074  # Oenococcus alcoholitolerans
    m[501841] = 502049  # Thalassospira profundimaris
    m[1226353] = 1518501  # Bradyrhizobium valentinum
    m[1075394] = 1118451  # Rhizobium naphthalenivorans
    m[1141655] = 1141657  # Nocardia vulneris
    m[1262586] = 1262585  # Acinetobacter qingfengensis
    m[1166681] = 1166683  #

    # taxonId =  1166683

    # print tax.getScientificName(taxonId)
    # print tax.exists(taxonId)
    return m


def modifySilvaMapping():
    taxFile = '/home/igregor/Documents/work/taxonomy/NCBI201502/ncbitax_sqlite.db'
    inDir = '/home/igregor/Documents/work/ppsp/silva115/db_old'
    outDir = '/home/igregor/Documents/work/ppsp/silva115/db'

    tax = taxonomy_ncbi.TaxonomyNcbi(taxFile)
    taxW = create_mg_db.TaxonomyWrap(taxFile)

    mapping = getMapping()

    # for all tax files
    for f in os.listdir(inDir):
        fPath = os.path.join(inDir, f)
        if f.endswith('.tax') and os.path.isfile(fPath):

            # create output file
            out = csv.OutFileBuffer(os.path.join(outDir, f))

            for line in open(fPath):
                line = line.strip()
                if len(line) > 0:
                    tokens = line.split('\t')
                    seqId = tokens[0]
                    taxonId = int(seqId.split('|')[-1].split(':')[1])
                    taxonId = mapping.get(taxonId, taxonId)
                    if not tax.exists(taxonId):
                        taxonId = int(map(lambda x: int(x), tokens[1].split(';')[:-1])[-1])

                    if tax.exists(taxonId):
                        taxPath = taxW.getPathToRoot(taxonId)
                        out.writeText('%s\t%s\n' % (seqId, taxPath))
                    else:
                        print('Problems to parse: %s\t taxon id not found' % line)

            print('Processed: %s' % f)
            out.close()

    tax.close()
    taxW.close()


# MAIN
if __name__ == "__main__":
    # _main()
    # _test1()
    # checkUnique()
    # getBp()
    modifySilvaMapping()