/*
 * kmer_counter.c
 *
 *  Created on: Nov 28, 2011
 *
 * Kmer counting based on a string matching algorithm ~ Rabin Karp
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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "kmer_counter.h"
#include "kmer_features.h"
#include "kmer_seq.h"

typedef struct {
	int * kmerArray;
	int kmerArrayLen;
	unsigned long * kmerOffsetArray;
	unsigned long featureArrayRequiredLen;
	unsigned long * reverseComplementMapArray;
	unsigned long * compatibilityMap;//this is allocated if needed (can be NULL)
} KCounter;

/* commnet out to remove superfluous checks */
#define KMER_COUNTER_DEBUG


/* helper method declaration */
inline void initWindow(char * sequence, int kmer, unsigned long * window, unsigned long * windowNBmp);

inline void updateWindow(char nextChar, int kmerMinus1, int twoTimesKmerMinus1,
		unsigned long * window, unsigned long * windowNBmp);

inline void processWindow(KCounter * kc, unsigned long * window, unsigned long * windowNBmp,
		double * featureArray/*, unsigned long * kmerCounterLong*/);

inline void getKmerOffsetArray(unsigned long * kmerOffsetArray, int * kmerArray, int kmerArrayLen);

inline void old_seq2str(unsigned long sequence, int length, char *str);

inline unsigned long * getCompatibilityMap(KmerCounter * kmerCounter);

inline unsigned long getReverseComplementIndex(int kmer, unsigned long index){
	int i;
	unsigned long reverseIndex = 0;
	for (i=0; i<kmer; i++){
		reverseIndex = (reverseIndex << 2) + ((index & 3) ^ 1);
		index >>= 2;
	}
	return reverseIndex;
}

/* PUBLIC methods */


/**
 * Returns mapping: indices (Kaustubh`s implementation) -> indices (in this kmer_counter.c implementation)
 *
 * Note that you cannot free this memory!!!
 *
 * */
inline unsigned long * kmerCounterGetCompatibilityMap(KmerCounter * kmerCounter){
	KCounter * kc = (KCounter *)kmerCounter;
	if (kc == NULL){
		return NULL;
	}
	if (kc->compatibilityMap == NULL){
		kc->compatibilityMap = getCompatibilityMap((KmerCounter *) kc);
	}
	return kc->compatibilityMap;
}




/**
 * @return Don`t free this array!
 * */
inline unsigned long * kmerCounterGetReverseComplementMap(KmerCounter * kmerCounter){
	KCounter * kc = (KCounter *)kmerCounter;
	if (kc == NULL || kc->reverseComplementMapArray == NULL){
		return NULL;
	} else {
		return kc->reverseComplementMapArray;
	}
}


/**
 * For each index defines defines an index of the reverse complement feature.
 *
 * @param outMapArray map: (feature_index) -> (reverse complement feature index in outMapArray[feature_index])
 * @param outMapArrayLen length of the outMapArray must be long enough, is set to the required array length
 *
 * @return 0 if OK
 * */
inline int kmerCounterComputeReverseComplementMap(KmerCounter * kmerCounter, unsigned long * outMapArray,
		unsigned long * outMapArrayLen){
	KCounter * kc = (KCounter *)kmerCounter;
	if (kc == NULL || outMapArray == NULL || outMapArrayLen == NULL
			|| *outMapArrayLen < kc->featureArrayRequiredLen || kc->kmerArray == NULL ||
			kc->kmerOffsetArray == NULL){
		return 1;
	}
	*outMapArrayLen = kc->featureArrayRequiredLen;

	unsigned long i,j,size,kmer,offset;
	for (i=0; i<kc->kmerArrayLen; i++){
		offset = kc->kmerOffsetArray[i];
		kmer = kc->kmerArray[i];
		size = pow(4,kmer);
		for (j=0; j<size; j++){
			outMapArray[j+offset] = getReverseComplementIndex(kmer, j) + offset;
		}
	}
	return 0;
}


inline unsigned long kmerCounterGetFeatureArrayLenRequired(KmerCounter * kmerCounter){
	KCounter * kc = (KCounter *)kmerCounter;
	if (kc != NULL){
		return kc->featureArrayRequiredLen;
	} else {
		return 0;
	}
}

inline KmerCounter * kmerCounterCreate(int * kmerArray, int kmerArrayLen){
	if (kmerArray == NULL){
		printf("kmerCounterCreate: kmer_array == NULL\n");
		return NULL;
	}
	if (kmerArrayLen <= 0){
		printf("kmerCounterCreate: kmer_array_len = %d\n", kmerArrayLen);
		return NULL;
	}
	if (kmerArray[kmerArrayLen-1] > KMER_COUNTER_MAX_KMER_LEN){
		printf("kmerCounterCreate: the maximum kmer length is %d (length %d not supported)",
				KMER_COUNTER_MAX_KMER_LEN, kmerArray[kmerArrayLen-1]);
		return NULL;
	}
	unsigned short i;
	unsigned long featuresNeeded = 0;
	unsigned short kmer = 0;
	unsigned short temp = 0;
	for (i=0; i<kmerArrayLen; i++){
		temp = kmer;
		kmer = kmerArray[i];
		if (temp >= kmer){
			printf("kmerCounterCreate: kmers are not in increasing order or they are negative\n");
			return NULL;
		}
		featuresNeeded += pow(4, kmerArray[i]);
	}

	KCounter * kc = (KCounter *)malloc(sizeof(KCounter));
	if (kc == NULL){
		return NULL;
	}

	kc->featureArrayRequiredLen = featuresNeeded;
	kc->kmerArray = kmerArray;
	kc->kmerArrayLen = kmerArrayLen;

	kc->kmerOffsetArray = (unsigned long *)malloc(kmerArrayLen * sizeof(unsigned long *));
	if (kc->kmerOffsetArray == NULL){
		free(kc);
		return NULL;
	}
	//for each kmer length compute offset in the feature array
	getKmerOffsetArray(kc->kmerOffsetArray, kmerArray, kmerArrayLen);

	kc->reverseComplementMapArray = (unsigned long *)malloc(sizeof(unsigned long)*(kc->featureArrayRequiredLen));
	if (kc->reverseComplementMapArray == NULL){
		free(kc->kmerOffsetArray);
		free(kc);
		return NULL;
	}
	//init reverse complement map array
	unsigned long len = kc->featureArrayRequiredLen;
	if (kmerCounterComputeReverseComplementMap((KmerCounter *)kc, kc->reverseComplementMapArray, &len) != 0
			|| len != kc->featureArrayRequiredLen){
		free(kc->reverseComplementMapArray);
		free(kc->kmerOffsetArray);
		free(kc);
		return NULL;
	}
	kc->compatibilityMap = NULL;//this is allocated only if needed for the first time

	return (KmerCounter *)kc;
}


inline void kmerCounterDestroy(KmerCounter * kmerCounter){
	KCounter * kc = (KCounter *)kmerCounter;
	if (kc != NULL){
		free(kc->kmerOffsetArray);
		free(kc->reverseComplementMapArray);
		if (kc->compatibilityMap != NULL){//this is allocated if needed
			free(kc->compatibilityMap);
		}
		free(kc);
	}
}


/*extract feature count from the array
double getKmerCount(int * featureArray, int featureArrayLen, char * feature, int featureLen,
		int * kmerArray, int kmerArrayLen);
*/

/*set mapping feature -> index
  void setLookupTable(int * tableArray, int tableArrayLen, char * feature, int featureLen,
		int * kmerArray, int kmerArrayLen);
*/


/**
 *
 * */
//inline unsigned long getKmerIndex(unsigned long kmerOffset, char * kmer, int kmerLen){
inline unsigned long kmerCounterGetKmerIndex(KmerCounter * kmerCounter, char * kmer, int kmerLen){
	KCounter * kc = (KCounter *)kmerCounter;
	if (kc == NULL){
		printf("kmerCounterGetKmerIndex: KmerCounter is NULL\n");
		exit(1);
	}

	unsigned char i;
	unsigned long kmerOffset = 0;
	unsigned char found = 0;
	for (i=0;i<kc->kmerArrayLen; i++){
		if (kc->kmerArray[i] == kmerLen){
			kmerOffset = kc->kmerOffsetArray[i];
			found = 1;
			break;
		}
	}
	if (found == 0){
		printf("kmerCounterGetKmerIndex: Cannot find index for kmer: %s",kmer);
		exit(1);
	}

	unsigned long w = 0;
	unsigned char c;
	for (i=0; i<kmerLen; i++){
		switch (kmer[i]){
			case 'a':
			case 'A': c = 0; break;
			case 't':
			case 'T': c = 1; break;
			case 'g':
			case 'G': c = 2; break;
			case 'c':
			case 'C': c = 3; break;
			default: {
				printf("kmerCounterGetKmerIndex: Kmer contains an undefined character: %c\n",kmer[i]);
				exit(1);
			};
			break;
		}
		w = (w << 2) + c;
	}
	return w + kmerOffset;
}


/**
 * Counts all kmers contained in the "sequence" and store its frequences in the "feature array".
 * Kmers are computed only from the substrings that contain "AaTtGgCc", i.e. substrings that contain
 * any other character are not counted. Indices of kmers correspond to the number representation using base 4,
 * where A~0,T~1,G~2,C~3 (e.g. index of kmer TGAC will be: 1*4^3 + 2*4^2 + 0*4^1 + 3*4^0 = 99 in the feature
 * array) If kmers of different length are stored in the feature array, indices of short kmers are stored
 * first, i.e. to store counts of kmers of length 1,2, and 3: kmers of length 1 will occupy indices 0-3
 * of the feature array, kmers of length 2 will occupy indices 4-19; and kmers of length 3 will occupy
 * indices 20-83. Each entry of the feature array correspond to the number - how many times a particular
 * kmer is contained in the sequence. The maximum kmer length is defined in the header.
 * (Kmer is a string that contain [ACGTacgt]of a specific length)
 *
 * @param sequence: an array of characters that correspond to a DNA sequence
 * @param seqLen: length of the sequence (max. unsigned long)
 * @param kmerArray: array of kmer lengths in ascending order! (e.g. [4,5,6])
 * @param kmerArrayLen: length of the kmerArray (e.g. 3)
 * @param featureArray: output array to store feature counts (max. unsigned long)
 * @param featureArrayLen: the length of the featureArray (max. unsigned long)
 * (@param featureShuffleTable: store mapping: feature internal index -> feature required index TODO if needed)
 *  */
inline int kmerCounterCountKmers(KmerCounter * kmerCounter, DNASeq * seq, KmerFeatureVector * featureVector){

	kmerFeatureVectorClear(featureVector);
	double * featureArray = kmerFeatureGetFeatureArray(featureVector);
	unsigned long featureArrayLen = kmerFeatureGetFeatureArrayLen(featureVector);
	/*unsigned long * kmerCounterLong =  kmerFeatureGetKmerCounterLong(featureVector);*/
	char * sequence = seqGetSeq(seq);
	unsigned long seqLen = seqGetSeqLen(seq);


	KCounter * kc = (KCounter *)kmerCounter;
	if (kc == NULL){
		return 1;
	}
	if (sequence == NULL){
		printf("kmerCounterCountKmers: sequence = NULL\n");
		return 1;
	}
	if (seqLen <= 0){
		printf("kmerCounterCountKmers: seq_len = %lu\n", seqLen);
		return 1;
	}
	if (featureArray == NULL){
		printf("kmerCounterCountKmers: feature_array = NULL\n");
		return 1;
	}
	/*if (kmerCounterLong == NULL){
		printf("kmerCounterCountKmers: kmerCounterLong = NULL\n");
		return 1;
	}*/
	if ((kc->kmerArray[kc->kmerArrayLen-1]) > seqLen){
		printf("kmerCounterCountKmers: the longest kmer (%d) is longer than the sequence (%lu)\n",
				(kc->kmerArray[kc->kmerArrayLen-1]), seqLen);
		return 1;
	}
	if (kc->featureArrayRequiredLen != featureArrayLen){
		printf("kmerCounterCountKmers: the length of the feature_array provided"
				" (%lu) is different from the length of the array needed %lu\n",
				featureArrayLen, kc->featureArrayRequiredLen);
		return 1;
	}

	/* init feature array */
	//unsigned long i;
	//for (i=0; i<featureArrayLen; i++){
	//	featureArray[i] = 0.0;
	//}

	/* init param */
	unsigned long window; //sliding window that represents the biggest kmer
	unsigned long windowNBmp;//sliding window that represents bitmap of undefined characters
	int maxKmer = kc->kmerArray[kc->kmerArrayLen-1];//length of the maximum kmer
	int maxKmerMinus1 = maxKmer -1;
	int twoTimesMaxKmerMinus1 = maxKmerMinus1 << 1;
	int minKmer = kc->kmerArray[0];//length of the minimum kmer
	//unsigned long kmerOffsetArray[kc->kmerArrayLen];//for each kmer length compute offset in the feature array
	//getKmerOffsetArray(kmerOffsetArray, kmerArray, kmerArrayLen);

	/* init window with the first longest kmer */
	initWindow(sequence, maxKmer, &window, &windowNBmp);

	/* increment indices of all kmers in the window that start at index 0 */
	processWindow(kc, &window, &windowNBmp, featureArray/*, kmerCounterLong*/);
	//printf("%d %d\n", window, windowNBmp);

	/* go with the sliding window till the end of the sequence */
	unsigned long index;
	for (index = maxKmer; index < seqLen; index++){

		/* update the window by reading one more character on the right (removing one on the left) */
		updateWindow(sequence[index], maxKmerMinus1, twoTimesMaxKmerMinus1, &window, &windowNBmp);

		/* increment kmers` indices */
		processWindow(kc, &window, &windowNBmp, featureArray/*, kmerCounterLong*/);
		//printf("%d %d\n", window, windowNBmp);
	}

	/* go with the sliding window beyond the end of the sequence to count shorter kmers */
	for (index=0; index < (maxKmer - minKmer); index++){
		updateWindow('N', maxKmerMinus1, twoTimesMaxKmerMinus1, &window, &windowNBmp);
		processWindow(kc, &window, &windowNBmp, featureArray/*, kmerCounterLong*/);
		//printf("%d %d\n", window, windowNBmp);
	}

	return 0;
}

/* PRIVATE methods */

/**
 * Increment the kmers` indices in the feature array for all kmers that start at the first index of the window.
 * If a non-defined character is a part of the window, then the kmer won`t be counted.
 *
 * (TODO add the featureShuffleTable if needed)
 * */
inline void processWindow(KCounter * kc, unsigned long * window, unsigned long * windowNBmp,
		double * featureArray/*, unsigned long * kmerCounterLong*/){

	short i;
	unsigned long index;
	unsigned long b;
	unsigned short offset;

	/* store the index if the window that represents the biggest kmer doesn`t contain any non-defined character */
	if ((*windowNBmp) == 0){
		index = (*window) + kc->kmerOffsetArray[kc->kmerArrayLen-1];
		++(featureArray[index]);
		/*++(*kmerCounterLong);*/
		//TODO handle also shorter kmers here not later
	}

	/* for all shorter kmers (than the longest one) store the index if the smaller window that represent a
	 * respective kmer doesn`t contain any non-defined character */
	for (i=kc->kmerArrayLen-2; i>=0 ; i--){
		offset = kc->kmerArray[kc->kmerArrayLen-1] - kc->kmerArray[i];//max kmer length - current kmer length
		b = (*windowNBmp) >> offset;//get bitmap of the current smaller window
		if (b == 0){//no non-defined characters in b
			index = ((*window) >> (offset << 1)) + kc->kmerOffsetArray[i];
			++(featureArray[index]);
			/*++(*kmerCounterLong);*/
		}
	}
}


/**
 * Initialize the window of the biggest kmer that starts at index 0.
 *
 * @param kmer: the length of the longest kmer
 *
 * */
inline void initWindow(char * sequence, int kmer, unsigned long * window, unsigned long * windowNBmp){

	unsigned short i;
	unsigned short c;
	unsigned short n;
	//unsigned long w = 0;
	*window = 0;
	*windowNBmp = 0;
	//unsigned long b = 0;
	for (i=0; i<kmer; i++){
		n = 0;
		switch (sequence[i]){
		case 'a':
		case 'A': c = 0; break;
		case 't':
		case 'T': c = 1; break;
		case 'g':
		case 'G': c = 2; break;
		case 'c':
		case 'C': c = 3; break;
		default:  {c = 0; n = 1;}; break;
		}
		*window = ((*window) << 2) + c;
		*windowNBmp = ((*windowNBmp) << 1) + n;
		//printf("%d %d\n",w,b);
	}
	//*window = w;
	//*windowNBmp = b;
}


/**
 * Updates the window of the biggest kmer, i.e. removes one char on the left and adds one on the right.
 *
 * @param nextChar: a new char that is being added
 * @param kmerMinus1: the biggest kmer length minus one
 * @param twoTimesKmerMinus1: two times the length of the biggest kmer minus one
 * */
inline void updateWindow(char nextChar, int kmerMinus1, int twoTimesKmerMinus1,
		unsigned long * window, unsigned long * windowNBmp){

	unsigned short c;
	unsigned short n = 0;
	switch (nextChar){
		case 'a':
		case 'A': c = 0; break;
		case 't':
		case 'T': c = 1; break;
		case 'g':
		case 'G': c = 2; break;
		case 'c':
		case 'C': c = 3; break;
		default:  {c = 0; n = 1;}; break;
	}

	/* updates the window and the bitmap
	 *
	 * index_i+1 = ( index_i - a_i*4^(k-1) )*4 + a_(i+k)
	 * */
	*window = (((*window) - (((*window) >> twoTimesKmerMinus1) << twoTimesKmerMinus1)) << 2) + c;
	*windowNBmp = (((*windowNBmp) - (((*windowNBmp) >> kmerMinus1) << kmerMinus1)) << 1) + n;
}


/**
 *
 * */
inline void getKmerOffsetArray(unsigned long * kmerOffsetArray, int * kmerArray, int kmerArrayLen){
	int i;
	kmerOffsetArray[0] = 0;//the shortest kmer`s indices start from 0 in the feature array
	for (i=1; i<kmerArrayLen; i++){
		kmerOffsetArray[i] = kmerOffsetArray[i-1] + pow(4,kmerArray[i-1]);
	}
}


/**
 * Maps an index of a kmer of certain length, in the original Kaustubh`s implementation
 * (function "seq2str" in "features.c", to a string of the same length.
 * Note that this mapping is different in this implementation (kmer_counter.c)!
 *
 * @param outKmerAsString output parameter must have length at least "kmerLength"
 * */
inline void old_seq2str(unsigned long index, int kmerLength, char *outKmerAsStr){
	int i;
	static const char map[] = {'A','C','T','G'};
	/* build kmer string */
	for (i=0;i<kmerLength;i++) {
		outKmerAsStr[i] = map[index&0x3];
		index = index >> 2;
	}
}


/**
 * @see kmerCounterGetCompatibilityMap
 * */
inline unsigned long * getCompatibilityMap(KmerCounter * kmerCounter){

	if (kmerCounter == NULL){
		return NULL;
	}

	KCounter * kc = (KCounter *)kmerCounter;
	int * kmerArray = kc->kmerArray;
	int kmerArrayLen = kc->kmerArrayLen;
	if (kmerArray == NULL || kmerArrayLen < 1){
		return NULL;
	}

	unsigned long outArrayLen = kmerCounterGetFeatureArrayLenRequired(kmerCounter);
	unsigned long * outArray = (unsigned long*)malloc(sizeof(unsigned long)*outArrayLen);
	if (outArray == NULL){
		return NULL;
	}

	unsigned long i, j, kmerLen, featureArrayLen, outIndex;
	char buffer[KMER_COUNTER_MAX_KMER_LEN];
	outIndex = 0;
	for (i=0; i<kmerArrayLen; i++){
		kmerLen = kmerArray[i];
		featureArrayLen = pow(4,kmerLen);
		for (j=0; j<featureArrayLen; j++){
			old_seq2str(j, kmerLen, buffer);
			outArray[outIndex] = kmerCounterGetKmerIndex(kmerCounter, buffer, kmerLen);
			++outIndex;
		}
	}
	return outArray;
}


