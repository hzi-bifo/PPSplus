

import os
from algbioi.com import csv
from algbioi.com import fasta
from algbioi.com import common

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

def _main():
    concatenate('/Volumes/My_Passport_Mac/work/rel_1_3/cami/assemblies',
                '/Volumes/My_Passport_Mac/work/rel_1_3/cami/cami_genomes_concat.fna')

if __name__ == "__main__":
    _main()