#!/usr/bin/env python

import os


if __name__ == "__main__":
    inDir = '/local/igregor/ref_cami_1/bac_arch/centroids_accessions'
    outDir = '/local/igregor/ref_cami_1/bac_arch/centroids_accessions2'

    maxPerFile = 100

    for f in os.listdir(inDir):
        fPath = os.path.join(inDir, f)
        if os.path.isfile(fPath):
            lineList = []
            for line in open(fPath):
                line = line.strip()
                lineList.append(line)

            outFileCount = len(lineList) / maxPerFile
            if len(lineList) % maxPerFile > 0:
                outFileCount += 1

            writeLineCount = 0
            for i in range(outFileCount):
                subset = lineList[i*maxPerFile: i*maxPerFile + maxPerFile]

                fw = open(os.path.join(outDir, '%s_%s' % (i, f)), 'w')
                for e in subset:
                    fw.write('%s\n' % e)
                    writeLineCount += 1
                fw.close()
            assert writeLineCount == len(lineList), '%s %s' % (writeLineCount, len(lineList))
