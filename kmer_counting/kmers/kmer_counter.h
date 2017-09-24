/*
 *  kmer_counter.h
 *
 *  Created on: Nov 28, 2011
 *
 *
 *  Copyright (C) 2014  Ivan Gregor
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef KMER_COUNTER_H_
#define KMER_COUNTER_H_

#include "kmer_features.h"
#include "kmer_seq.h"

/* The maximum length of the longest kmer, make sure that the sliding window and all features`
 * indices fit into the "unsigned long", i.e.
 * (4^1 + 4^2 + ... + 4^KMER_COUNTER_MAX_KMER_LEN < (max. unsigned long)))
 *  */
#define KMER_COUNTER_MAX_KMER_LEN 15

typedef void KmerCounter;

KmerCounter * kmerCounterCreate(int * kmerArray, int kmerArrayLen);

void kmerCounterDestroy(KmerCounter * kmerCounter);

int kmerCounterCountKmers(KmerCounter * kmerCounter, DNASeq * seq, KmerFeatureVector * featureVector);

unsigned long kmerCounterGetKmerIndex(KmerCounter * kmerCounter, char * kmer, int kmerLen);

int kmerCounterComputeReverseComplementMap(KmerCounter * kmerCounter, unsigned long *outMapArray,
		unsigned long * outMapArrayLen);

unsigned long * kmerCounterGetReverseComplementMap(KmerCounter * kmerCounter);

unsigned long kmerCounterGetFeatureArrayLenRequired(KmerCounter * kmerCounter);

unsigned long * kmerCounterGetCompatibilityMap(KmerCounter * kmerCounter);

#endif /* KMER_COUNTER_H_ */
