/*
 *  kmer_features.c
 *
 *  Created on: Jan 3, 2012
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
#include <stdlib.h>

#include "kmer_features.h"

typedef struct {
	double * featureArray;
	unsigned long featureArrayLen;
	/*unsigned long kmerCounter;*/
} KFeatureVector;

inline KmerFeatureVector * kmerFeatureVectorCreate(unsigned long size){

	if (size < 1){
		return NULL;
	}
	KFeatureVector * v = (KFeatureVector *)malloc(sizeof(KFeatureVector));
	if (v == NULL){
		return NULL;
	}
	v->featureArrayLen = size;
	/*v->kmerCounter = 0;*/
	v->featureArray = (double*)malloc(sizeof(double)*(v->featureArrayLen));
	if (v->featureArray == NULL){
		free(v);
		return NULL;
	}

	unsigned long i, arrayLen = v->featureArrayLen;
	double * array = v->featureArray;
	for (i=0; i<arrayLen; i++){
		array[i] = 0.0;
	}

	return (KmerFeatureVector *)v;
}


inline double kmerFeatureGetFeatureAt(KmerFeatureVector * kmerFeatureVector, unsigned long index){
	KFeatureVector * v = (KFeatureVector *)kmerFeatureVector;
	if (v != NULL && v->featureArray != NULL && 0 <= index && index < v->featureArrayLen){
		return v->featureArray[index];
	} else {
		printf("kmerFeatureGetFeatureAt: wrong index %lu or kmerFeatureVector is corrupted\n",index);
		exit(1);
	}
}


inline double * kmerFeatureGetFeatureArray(KmerFeatureVector * kmerFeatureVector){
	KFeatureVector * v = (KFeatureVector *)kmerFeatureVector;
	if (v != NULL && v->featureArray != NULL){
		return v->featureArray;
	} else {
		printf("kmerFeatureGetFeatureArray: kmerFeatureVector is corrupted\n");
		exit(1);
	}
}


inline unsigned long kmerFeatureGetFeatureArrayLen(KmerFeatureVector * kmerFeatureVector){
	KFeatureVector * v = (KFeatureVector *)kmerFeatureVector;
	if (v != NULL && v->featureArray != NULL){
		return v->featureArrayLen;
	} else {
		printf("kmerFeatureGetFeatureArrayLen: kmerFeatureVector is corrupted\n");
		exit(1);
	}
}

/**
 * Gets a modified feature vector
 *
 * @return 0 if OK
 * */
inline int kmerFeatureGetModifiedVector(KmerFeatureVector * kmerFeatureVector,
		unsigned long * reverseComplementMap, unsigned long seqLen, int reverseComplement,
		int removeRedundantFeatures, int normalization, int backwardsCompatibility,
		unsigned long * backwardsCompatibilityMap, unsigned long * backwardsCompatibilityMapInv,
		double * outFeatureVector, unsigned long * outFeatureVectorLen){


	KFeatureVector * v = (KFeatureVector *)kmerFeatureVector;
	if (v == NULL || v->featureArray == NULL || outFeatureVector == NULL || reverseComplementMap == NULL ||
			(backwardsCompatibility == 1 && (backwardsCompatibilityMap == NULL ||
					backwardsCompatibilityMapInv == NULL))){
		//printf("*\n");
		return 1;
	}

	unsigned long i, j, i2, j2, resultFeatureArrayLen, arrayLen = v->featureArrayLen;
	double * featureArray = v->featureArray;
	double tmp;

	/* Correct counters as if we computed also reverse complement */
	if (reverseComplement == 1 && backwardsCompatibility == 0){
		for (i=0; i<arrayLen; i++){
			j = reverseComplementMap[i];
			if (i < j){
				tmp = featureArray[i] + featureArray[j];
				outFeatureVector[i] = tmp;
				outFeatureVector[j] = tmp;
			} else if (i == j){
				outFeatureVector[i] = featureArray[i] * 2.0;
			}
		}
	} else if (reverseComplement == 0 && backwardsCompatibility == 0){
		for (i=0; i<arrayLen; i++){
			outFeatureVector[i] = featureArray[i];
		}
	} else if (backwardsCompatibility == 1 && reverseComplement == 0){
		for (i=0; i<arrayLen; i++){
			outFeatureVector[i] = featureArray[backwardsCompatibilityMap[i]];
		}
	} else if (backwardsCompatibility == 1 && reverseComplement == 1){
		for (i=0; i<arrayLen; i++){
			j = reverseComplementMap[i];
			if (i < j){
				i2 = backwardsCompatibilityMapInv[i];
				j2 = backwardsCompatibilityMapInv[j];
				tmp = featureArray[i] + featureArray[j];
				outFeatureVector[i2] = tmp;
				outFeatureVector[j2] = tmp;
			} else if (i == j){
				i2 = backwardsCompatibilityMapInv[i];
				outFeatureVector[i2] = featureArray[i] * 2.0;
			}
		}
	} else {
		//printf("**\n");
		return 1;
	}

	/* remove redundant features */
	resultFeatureArrayLen = arrayLen;
	if (removeRedundantFeatures == 1){
		if (backwardsCompatibility == 0){
			resultFeatureArrayLen = 0;
			for (i=0; i<arrayLen; i++){
				j = reverseComplementMap[i];
				if (i<=j){
					outFeatureVector[resultFeatureArrayLen] = outFeatureVector[i];
					//printf("%lg ", outFeatureVector[resultFeatureArrayLen]);
					++resultFeatureArrayLen;
				}
			}
		} else if (backwardsCompatibility == 1){//check this again!!!
			resultFeatureArrayLen = 0;
			for (i2=0; i2<arrayLen; i2++){
				i = backwardsCompatibilityMap[i2];
				j = reverseComplementMap[i];
				j2 = backwardsCompatibilityMapInv[j];
				if (i2<=j2){
					outFeatureVector[resultFeatureArrayLen] = outFeatureVector[i2];
					//printf("%lg ", outFeatureVector[resultFeatureArrayLen]);
					++resultFeatureArrayLen;
				}
			}
		} else {
			//printf("***\n");
			return 1;
		}
	}

	/* Normalization */
	if (normalization == 1){
		for (i=0; i<resultFeatureArrayLen; i++){
			outFeatureVector[i] /= (double)seqLen;
			//printf("%lg", outFeatureVector[i]);
		}
	}

	(*outFeatureVectorLen) = resultFeatureArrayLen;
	//printf("\n%lu\n",(*outFeatureVectorLen));
	return 0;
}


/*inline unsigned long * kmerFeatureGetKmerCounterLong(KmerFeatureVector * kmerFeatureVector){
	KFeatureVector * v = (KFeatureVector *)kmerFeatureVector;
	if (v != NULL){
		return &(v->kmerCounter);
	} else {
		return NULL;
	}
}*/


inline void kmerFeatureVectorClear(KmerFeatureVector * kmerFeatureVector){
	KFeatureVector * v = (KFeatureVector *)kmerFeatureVector;
	if (v != NULL && v->featureArray != NULL){
		//v->kmerCounter = 0;
		unsigned long i, arrayLen = v->featureArrayLen;
		double * array = v->featureArray;
		for (i=0; i<arrayLen; i++){
			array[i] = 0.0;
		}
	}
}


inline void kmerFeatureVectorDestroy(KmerFeatureVector * kmerFeatureVector){
	KFeatureVector * v = (KFeatureVector *)kmerFeatureVector;
	if (v != NULL){
		if (v->featureArray != NULL){
			free(v->featureArray);
		}
		free(v);
	}
}
