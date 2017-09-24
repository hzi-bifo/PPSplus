/*
 *  kmer_features.h
 *
 *  Created on: Jan 3, 2012
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

#ifndef KMER_FEATURES_H_
#define KMER_FEATURES_H_

typedef void KmerFeatureVector;

KmerFeatureVector * kmerFeatureVectorCreate(unsigned long size);

double kmerFeatureGetFeatureAt(KmerFeatureVector * kmerFeatureVector, unsigned long index);

double * kmerFeatureGetFeatureArray(KmerFeatureVector * kmerFeatureVector);

unsigned long kmerFeatureGetFeatureArrayLen(KmerFeatureVector * kmerFeatureVector);

/*unsigned long * kmerFeatureGetKmerCounterLong(KmerFeatureVector * kmerFeatureVector);*/

int kmerFeatureGetModifiedVector(KmerFeatureVector * kmerFeatureVector, unsigned long * reverseComplementMap,
		unsigned long seqLen, int reverseComplement, int removeRedundantFeatures, int normalization,
		int backwardsCompatibility, unsigned long * backwardsCompatibilityMap,
		unsigned long * backwardsCompatibilityMapInv, double * outFeatureVector,
		unsigned long * outFeatureVectorLen);

void kmerFeatureVectorClear(KmerFeatureVector * kmerFeatureVector);


void kmerFeatureVectorDestroy(KmerFeatureVector * kmerFeatureVector);

#endif /* KMER_FEATURES_H_ */

