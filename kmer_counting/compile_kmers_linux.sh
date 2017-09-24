#!/bin/sh

gcc -O3 -Wall -fmessage-length=0 -MMD -MP -o fasta2kmers2 ./kmers/char_list.c  ./kmers/kmer_counter.c ./kmers/kmer_fasta.c ./kmers/kmer_features.c  ./kmers/kmer_main.c ./kmers/kmer_seq.c ./kmers/test.c -I. -lm
