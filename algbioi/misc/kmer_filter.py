"""
    Filter out reads containing rare kmers.

    Requires to install: jellyfish: http://www.cbcb.umd.edu/software/jellyfish/

"""

import os
import subprocess
import zlib
import pickle
from Bio import SeqIO
from Bio.Seq import Seq

# Settings
WORKING_DIR = '/net/4christina/projects/metagenomics/nobackup/igregor/kmer_filter/A'
# WORKING_DIR = '/Users/ivan/Documents/nobackup/kmer/working'
READS_FASTA_FILE = '/net/metagenomics/projects/barley_metagenomics/samples/A/trimmed/563_A_paired_trimmed_reads_named.fa'
# READS_FASTA_FILE = '/Users/ivan/Documents/nobackup/kmer/reads.fna'
# JELLYFISH = '/usr/local/bin/jellyfish'
JELLYFISH = '/net/metagenomics/projects/PPSmg/tools/jellyfish-1.1.10/install/bin/jellyfish'
KMER = '23'
THREADS = '10'
THRESHOLD = '2'
# KMER_SET_MAX_SIZE = 2000
KMER_SET_MAX_SIZE = 100000000

def indexToDna(index, length, mapping={0: 'A', 1: 'T', 2: 'G', 3: 'C'}):
    s = str(bin(index))[2:]  # get binary number, trim starting 0b
    s = ((length * 2) - len(s)) * '0' + s  # add starting zeros
    dna = ''
    for i in range(length):  # for each dna character
        dna += mapping[int(s[i*2:(i+1)*2], 2)]  # map binary to dna
    return dna

def indexLineToDnaLine(indexLine, kmerLen):
    line = ''
    for entry in indexLine.split('\t'):
        index, count = entry.split(':')
        if line != '':
            line += '\t'
        line += indexToDna(int(index), kmerLen) + ':' + str(count)
    return line

def dnaToIndex(kmer, kmerLen, mapping={'A': 0, 'T': 1, 'G': 2, 'C': 3}):
    index = 0
    for i in range(kmerLen):
        index <<= 2
        index += mapping[kmer[i]]
    return index

def _test():
    for index in range(64):
        dna = indexToDna(index, 3)
        index2 = dnaToIndex(dna, 3)
        if index != index2:
            print 'error'
        print index, index2, dna

def dnaToIndexSet(dnaSeq, kmerLen):
    s = set()
    for i in range(len(dnaSeq) - kmerLen + 1):
        try:
            dna = dnaSeq[i:(i + kmerLen)]
            index = dnaToIndex(dna, kmerLen)
            # dnaReverse = str(Seq(dna).complement())
            # indexReverse = dnaToIndex(dnaReverse, kmerLen)

        except KeyError:
            continue
        s.add(index)
        # s.add(indexReverse)
    return frozenset(s)




def fastaFileToBinaryFasta(fastaFile, outFilePath, kmerLen=None):

    fr = open(os.path.normpath(fastaFile), 'r')
    fw = open(outFilePath, 'w')
    readBuffer = SeqIO.parse(fr, 'fasta')
    count = 0
    for record in readBuffer:
        if count % 2 == 0:
            id1 = str(record.id)
            seq1 = str(record.seq)
        else:
            id2 = str(record.id)
            seq2 = str(record.seq)
            if kmerLen is not None:
                s = frozenset(dnaToIndexSet(seq1, kmerLen) | dnaToIndexSet(seq2, kmerLen))
            else:
                s = None
            entry = zlib.compress('>' + id1 + '\n' + seq1 + '\n>' + id2 + '\n' + seq2)
            #entry = '>' + id1 + '\n' + seq1 + '\n>' + id2 + '\n' + seq2 + '\n'

            pickle.dump((entry, s), fw)
        count += 1
    fr.close()
    fw.close()


def binaryFastaToFasta(fastaFileBinary, outFilePath):
    fr = open(os.path.normpath(fastaFileBinary), 'r')
    fw = open(outFilePath, 'w')
    while True:
        try:
            entry = pickle.load(fr)
        except EOFError:
            break
        fw.write(zlib.decompress(entry[0]) + '\n')
    fr.close()
    fw.close()




def main():
    #
    # Count kmers
    outFile = 'kmer' + KMER + '_max' + THRESHOLD
    if False:
        cmd = JELLYFISH + ' count -m ' + KMER + ' -s 1000000000 -t ' + THREADS + ' -o ' + outFile + ' -C -U ' + THRESHOLD + \
               ' --invalid-char=ignore ' + READS_FASTA_FILE
        #cmd = JELLYFISH + ' count -m ' + KMER + ' -s 1000000000 -t ' + THREADS + ' -o ' + outFile + ' ' + READS_FASTA_FILE
        # TODO: choose the right command !!!
        print('Counting kmers: ' + cmd)
        proc = subprocess.Popen(cmd, shell=True, bufsize=-1, cwd=WORKING_DIR)
        proc.wait()
        if proc.returncode != 0:
            print('Counting kmers returned with non-zero status %s' % proc.returncode)
            return
        else:
            print('Counting kmers completed')

    # Merge output files
    if False:
        if len(os.listdir(WORKING_DIR)) > 1:
            cmd = JELLYFISH + ' merge -o ' + outFile + '_merged ' + outFile + '_*'
            print('Merging output files:' + cmd)
            proc = subprocess.Popen(cmd, shell=True, bufsize=-1, cwd=WORKING_DIR)
            proc.wait()
            if proc.returncode != 0:
                print('Merging kmers returned with non-zero status %s' % proc.returncode)
                return
            else:
                print('Merging kmers completed')
        else:
            os.rename(os.path.join(WORKING_DIR, str(outFile + '_0')), os.path.join(WORKING_DIR, str(outFile + '_merged')))
    outFile = str(outFile + '_merged')

    # Dump kmers
    if False:
        cmd = JELLYFISH + " dump -c " + outFile + " | awk '{ if ($2 <= " + THRESHOLD + ") {print $0} }'  | cat > " + \
              outFile + ".csv"
        print('Dumping kmers: ' + cmd)
        proc = subprocess.Popen(cmd, shell=True, bufsize=-1, cwd=WORKING_DIR)
        proc.wait()
        if proc.returncode != 0:
            print('Dumping kmers returned with non-zero status %s' % proc.returncode)
            return
        else:
            print('Dumping kmers completed')
    outFile = str(outFile + '.csv')

    # Read in kmers
    kmerSetsFileName = os.path.join(WORKING_DIR, 'kmerSets')
    if False:
        print('Reading kmers in binary')
        kmerSetsFile = open(kmerSetsFileName, 'w')
        kmer = int(KMER)
        s = set()
        count = 0
        for line in open(os.path.join(WORKING_DIR, outFile), 'r'):
            try:
                dna = line.split()[0]
                index = dnaToIndex(dna, kmer)
                dnaReverse = str(Seq(dna).complement())
                indexReverse = dnaToIndex(dnaReverse, kmer)
                # print line.split()[0]
                # print index
                # print kmer
            except KeyError:
                continue
            s.add(index)
            s.add(indexReverse)
            count += 2
            if count >= KMER_SET_MAX_SIZE:
                count = 0
                pickle.dump(s, kmerSetsFile)
                s = set()
        pickle.dump(s, kmerSetsFile)
        kmerSetsFile.close()
        print('Reading kmers in binary done.')

    # Store fasta file in binary format !
    if False:
        print('Storing fasta file as binary')
        fastaFileToBinaryFasta(READS_FASTA_FILE, os.path.join(WORKING_DIR, 'fasta0'), kmerLen=int(KMER))
        print('Storing fasta file as binary done')

    # Filter out kmers
    if True:
        print('filtering kmers')
        kmerSetsFile = open(kmerSetsFileName, 'r')
        count = 0
        while True:
            count += 1
            print('filter iter: %s' % count)
            try:
                s = pickle.load(kmerSetsFile)
            except EOFError:
                break
            fr = open(os.path.join(WORKING_DIR, 'fasta0'), 'r')
            fw = open(os.path.join(WORKING_DIR, 'fasta1'), 'w')
            while True:
                try:
                    entry = pickle.load(fr)
                except EOFError:
                    break
                if len(s.intersection(entry[1])) == 0:
                    pickle.dump(entry, fw)
                # print(str(entry[1]))
                # print(str(s))

            fr.close()
            fw.close()
            os.remove(os.path.join(WORKING_DIR, 'fasta0'))
            os.rename(os.path.join(WORKING_DIR, 'fasta1'), os.path.join(WORKING_DIR, 'fasta0'))
        kmerSetsFile.close()
        print('filtering kmers done')


    # Store filtered sequences
    if False:
        print('storing filtered sequences')
        binaryFastaToFasta(os.path.join(WORKING_DIR, 'fasta0'), os.path.join(WORKING_DIR, 'fasta_result.fna'))
        print('DONE')

if __name__ == "__main__":
    main()
    # _test()

# count kmers OLD
# jellyfish count -m 23 -s 10000000 -t 20 -o kmer23_max2 -C -U 2 --invalid-char=ignore /net/metagenomics/projects/barley_metagenomics/samples/A/trimmed/563_A_paired_trimmed_reads_named.fa
#
# merge kmers
# jellyfish merge -o merged.txt kmer23_max2_*
#
# dump kmers
# jellyfish dump -c -U 2 -o kmer23_max2.csv merged_kmer23_max2
#
# jellyfish dump -c -U 2 merged_kmer23_max2 | awk '{ if ($2 < 3) {print $0} }'  | cat > kmer23_max2__.csv

