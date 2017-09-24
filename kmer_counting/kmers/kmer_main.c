/*
 *  kmer_main.c
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "kmer_main.h"
#include "kmer_seq.h"
#include "kmer_features.h"
#include "kmer_counter.h"
#include "kmer_fasta.h"


inline int storeOutputFeatureVector(FILE * outputFile, double * outFeatureVector,
		unsigned long outFeatureVectorLen, int outputHeaderFormat, int outputFeatureVectorFormat,
		int outputFeatureVectorIndexOffset, int labelFormat, long defaultLabelNcbid,
		DNASeq * seq);


/**
 * @param inputFastaFilePath path to a fasta file (not NULL)
 * @param outputFeatureFilePath path to an output file (not NULL)
 * @param outputFopenMode w~for writing (original file is erased), a~append to an existing file or create it
 * @param kmerFrom less or equal to kmerTo
 * @param kmerTo (e.g. if kmerFrom~4,kmerTo~6 then kmers 4,5,6 are computed) max. ~15
 * @param reverseComplement 1 ~ consider also reverse complement (results in half of the features) else 0
 * @param removeRedundantFeatures 1~remove redundant features if reverse complement~1 else 0
 * @param normalization 0~ no normalization, 1 ~ sequence length
 * @param labelFormat 0~ no labels, 1~ label defined in seq name after "label:", 2~label equal to "1"
 * @param outputFeatureVectorFormat 0 ~ plain, 1 ~ sparse
 * @param outputFeatureVectorIndexOffset offset of indices in the output feature file
 * @param outputHeaderFormat 1 ~ #seq name appended at the end of the line, 0 ~ no header
 * @param backwardsCompatibility [1 ~ order of the features correspond to the order in the original
 * implementation,] 0 ~ order of the features is according to this implementation
 *
 * @return 0 ~ OK
 * */
int kmerMainComputeFeatures(char * inputFastaFilePath, char * outputFeatureFilePath, char outputFopenMode,
		int kmerFrom, int kmerTo, int reverseComplement, int removeRedundantFeatures, int normalization,
		int labelFormat, int outputFeatureVectorFormat, int outputFeatureVectorIndexOffset,
		int outputHeaderFormat, int backwardsCompatibility){

	/* check parameters */
	if (kmerFrom <= 0 || kmerTo < kmerFrom || kmerTo > KMER_COUNTER_MAX_KMER_LEN){
		printf("kmerMainComputeFeatures: wrong kmerFrom(%d) or kmerTo(%d) parameters,"
				" max. kmer is: (%d)!\n", kmerFrom, kmerTo, KMER_COUNTER_MAX_KMER_LEN);
		return 1;
	}

	if (outputFopenMode != 'w' && outputFopenMode != 'a'){
		printf("kmerMainComputeFeatures: wrong outputFopenMode :%c, allowed modes are (w or a)!\n",
				outputFopenMode);
		return 2;
	}

	if (reverseComplement != 0 && reverseComplement != 1){
		printf("kmerMainComputeFeatures: wrong value of reverseComplement param: %d!\n",reverseComplement);
		return 3;
	}
	if (((removeRedundantFeatures == 1) && (reverseComplement != 1))){
		printf("kmerMainComputeFeatures: wrong combination of parameters: removeRedundantFeatures=%d "
				"and reverseComplement=%d!\n",
				removeRedundantFeatures, reverseComplement);
		return 4;
	}
	if (normalization < 0 || normalization > 1){
		printf("kmerMainComputeFeatures: not supported normalization param: %d!\n",normalization);
		return 5;
	}
	if (labelFormat < 0 || 2 < labelFormat){
		printf("kmerMainComputeFeatures: not supported label format: %d!\n",labelFormat);
		return 6;
	}
	if (outputFeatureVectorFormat != 0 && outputFeatureVectorFormat != 1){
		printf("kmerMainComputeFeatures: not supported output feature vector format: %d!\n",
				outputFeatureVectorFormat);
		return 7;
	}
	if (outputHeaderFormat != 0 && outputHeaderFormat != 1){
		printf("kmerMainComputeFeatures: not supported output header format: %d!\n",outputHeaderFormat);
		return 8;
	}
	if (backwardsCompatibility != 0 && backwardsCompatibility != 1){
		printf("kmerMainComputeFeatures: not supported backwardsCompatibility param: %d!\n",
				backwardsCompatibility);
		return 9;
	}

	/* Init structures */
	unsigned long i;
	int kmerArrayLen = kmerTo - kmerFrom + 1;
	int kmerArray[kmerArrayLen];
	for (i=0; i<kmerArrayLen; i++){
		kmerArray[i] = kmerFrom + i;
	}

	long defaultLabelNcbid = 0;
	if (inputFastaFilePath != NULL){
		char * fileName = strrchr(inputFastaFilePath,'/');
		if (fileName == NULL){
			fileName = inputFastaFilePath;
		}
		//printf("fileName%s\n",fileName);
		if ((strlen(fileName) > 1) && (sscanf(&fileName[1],"%ld", &defaultLabelNcbid) != 1)){
			defaultLabelNcbid = 0;
		}
		//printf("ncbid:%ld\n",defaultLabelNcbid);
	}

	FastaHandler * fh = readBufferInit(inputFastaFilePath);
	if (fh == NULL){
		printf("kmerMainComputeFeatures: Cannot initialize read buffer from: %s\n", inputFastaFilePath);
		return 10;
	}

	if (outputFeatureFilePath == NULL){
		printf("kmerMainComputeFeatures: Output feature file is NULL!\n");
		return 11;
	}
	FILE * outputFile = NULL;
	if (outputFopenMode == 'a'){
		outputFile = fopen(outputFeatureFilePath,"a");
	} else if (outputFopenMode == 'w'){
		outputFile = fopen(outputFeatureFilePath,"w");
	}
	if (outputFile == NULL){
		printf("kmerMainComputeFeatures: Cannot initialize output file: %s!\n", outputFeatureFilePath);
		readBufferDestroy(fh);
		return 12;
	}

	DNASeq * seq =  seqCreate();
	if (seq == NULL){
		printf("kmerMainComputeFeatures: Cannot initialize sequence struct (DNASeq)!\n");
		readBufferDestroy(fh);
		fclose(outputFile);
		return 13;
	}

	KmerCounter * kc = kmerCounterCreate(kmerArray, kmerArrayLen);
	if (kc == NULL){
		printf("kmerMainComputeFeatures: Cannot create KmerCounter!\n");
		readBufferDestroy(fh);
		fclose(outputFile);
		seqDestroy(seq);
		return 14;
	}

	KmerFeatureVector * featureVector = kmerFeatureVectorCreate(kmerCounterGetFeatureArrayLenRequired(kc));
	if (featureVector == NULL){
		printf("kmerMainComputeFeatures: Cannot create featureVector!\n");
		readBufferDestroy(fh);
		fclose(outputFile);
		seqDestroy(seq);
		kmerCounterDestroy(kc);
		return 15;
	}

	int errno = 0;
	double * outFeatureVector = (double *)malloc(sizeof(double)*kmerCounterGetFeatureArrayLenRequired(kc));
	if (outFeatureVector == NULL){
		printf("kmerMainComputeFeatures: Cannot create out feature vector (double*)!\n");
		readBufferDestroy(fh);
		fclose(outputFile);
		seqDestroy(seq);
		kmerCounterDestroy(kc);
		kmerFeatureVectorDestroy(featureVector);
		return 16;
	}
	unsigned long outFeatureVectorLen = 0;
	unsigned long * reverseComplementMap = kmerCounterGetReverseComplementMap(kc);
	unsigned long * backwardsCompatibilityMap = NULL;
	unsigned long * backwardsCompatibilityMapInv = NULL;
	if (backwardsCompatibility == 1){
		backwardsCompatibilityMap = kmerCounterGetCompatibilityMap(kc);//old->new
		if (backwardsCompatibilityMap != NULL){
			backwardsCompatibilityMapInv = (unsigned long *)malloc(sizeof(unsigned long)
					*kmerCounterGetFeatureArrayLenRequired(kc));
			if (backwardsCompatibilityMapInv != NULL){
				for (i=0; i<kmerCounterGetFeatureArrayLenRequired(kc); i++){
					backwardsCompatibilityMapInv[backwardsCompatibilityMap[i]] = i;//check this !!!
				}
			}
		}
		if (backwardsCompatibilityMap == NULL || backwardsCompatibilityMapInv == NULL){
			printf("kmerMainComputeFeatures: Cannot create backwards compatibility mapping!\n");
			free(outFeatureVector);
			readBufferDestroy(fh);
			fclose(outputFile);
			seqDestroy(seq);
			kmerCounterDestroy(kc);
			kmerFeatureVectorDestroy(featureVector);
			return 17;
		}
	}

	/* read sequence one by one, compute features, modify the feature vector and store it to the output file */
	while (readBufferHasNext(fh)){
		if (readBufferGetNextSeq(fh, seq) != 0){
			printf("kmerMainComputeFeatures: Cannot read sequence!\n");
			errno = 18;
			break;
		}

		/* count kmers */
		if (kmerCounterCountKmers(kc, seq, featureVector) != 0){
			printf("kmerMainComputeFeatures: count kmers returned non-zero value!\n");
			errno = 19;
			break;
		}

		/* get modified feature vector for the output */
		if (kmerFeatureGetModifiedVector(featureVector, reverseComplementMap, seqGetSeqLen(seq),
				reverseComplement,removeRedundantFeatures,normalization, backwardsCompatibility,
				backwardsCompatibilityMap, backwardsCompatibilityMapInv,
				outFeatureVector, &outFeatureVectorLen) != 0){
			printf("kmerMainComputeFeatures: cannot modify feature vector!\n");
			errno = 20;
			break;
		}

		/* store feature vector */
		if (storeOutputFeatureVector(outputFile, outFeatureVector, outFeatureVectorLen,
				outputHeaderFormat, outputFeatureVectorFormat, outputFeatureVectorIndexOffset,
				labelFormat, defaultLabelNcbid, seq) != 0){
			printf("kmerMainComputeFeatures: cannot store feature vector to a file!\n");
			errno = 21;
			break;
		}
	}

	/* free resources */
	if (backwardsCompatibilityMapInv != NULL){
		free(backwardsCompatibilityMapInv);
	}
	free(outFeatureVector);
	kmerFeatureVectorDestroy(featureVector);
	int errnoFclose = fclose(outputFile);
	if (errnoFclose != 0){
		printf("kmerMainComputeFeatures: fclose returned with non-zero status: %d!\n", errnoFclose);
		errno = errnoFclose;
	}
	kmerCounterDestroy(kc);
	seqDestroy(seq);
	readBufferDestroy(fh);
	return errno;
}


/**
 * Store a feature vector to a file.
 * Store first "outFeatureVectorLen" features of vector "outFeatureVector".
 * */
inline int storeOutputFeatureVector(FILE * outputFile, double * outFeatureVector,
		unsigned long outFeatureVectorLen, int outputHeaderFormat, int outputFeatureVectorFormat,
		int outputFeatureVectorIndexOffset, int labelFormat, long defaultLabelNcbid, DNASeq * seq){

 	/* sequence name */
 	char * seqName = seqGetName(seq);
 	if (seqName == NULL){
 		printf("kmer_main:storeOutputFeatureVector: seqName is NULL!\n");
 		return 1;
 	}

 	/* store sequence label if defined in sequence name as label:ncbid or 0 if label not found */
 	if (labelFormat == 1){

 		long ncbid;
 		char * labelNcbidStr = strstr(seqName, "label:");
 		if(labelNcbidStr != NULL){
 			if (sscanf(labelNcbidStr,"label:%ld", &ncbid) != 1){
 				ncbid = defaultLabelNcbid;
 				if (ncbid == 0){
 					printf("kmer_main:storeOutputFeatureVector: zero ncbid for sequence name(1): %s\n",seqName);
 				}
 			}
 		} else {
 			ncbid = defaultLabelNcbid;
 			if (ncbid == 0){
 				printf("kmer_main:storeOutputFeatureVector: zero ncbid for sequence name(2): %s\n",seqName);
 			}
 		}
 		if (fprintf(outputFile, "%ld\t", ncbid) < 0){
 			printf("kmer_main:storeOutputFeatureVector: cannot write label %ld to a file!\n",ncbid);
 			return 1;
 		}
 	} else if (labelFormat == 2){
 		if (fprintf(outputFile, "1\t") < 0){
 			printf("kmer_main:storeOutputFeatureVector: cannot write label `1` to a file "
 					"(option: labelFormat=2)!\n");
 			return 1;
 		}
 	} else if (labelFormat != 0){
 		printf("kmer_main:storeOutputFeatureVector: wrong labelFormat: %d\n", labelFormat);
 		return 1;
 	}

 	/* store feature vector */
 	unsigned long i;
 	if(outputFeatureVectorFormat == 1){
 		/* sparse format */
 		for(i=0; i<outFeatureVectorLen; i++) {
 			if (fabs(outFeatureVector[i]) > 0.000001){
 				if (fprintf(outputFile, "%ld:%lg\t", (long)(i+outputFeatureVectorIndexOffset),
 						outFeatureVector[i]) < 0){
 					printf("kmer_main:storeOutputFeatureVector: Cannot store feature (sparse) to a file!\n");
 					return 1;
 				}
 			}
 		}
 	} else if (outputFeatureVectorFormat == 0){
 		/* plain format*/
 		for(i=0; i<outFeatureVectorLen;i++){
 			if (fprintf(outputFile, "%lg\t", outFeatureVector[i]) < 0){
 				printf("kmer_main:storeOutputFeatureVector: Cannot store feature (plain) to a file!\n");
 				return 1;
 			}
 		}
 	} else {
 		/* wrong format */
 		printf("kmer_main:storeOutputFeatureVector: wrong outputFeatureVectorFormat: %d\n",
				outputFeatureVectorFormat);
		return 1;
	}

 	/* comment */
 	if (outputHeaderFormat == 1){
 		if (fprintf(outputFile, "#%s", seqName) < 0){
 			printf("kmer_main:storeOutputFeatureVector: Cannot store header to a file!\n");
 			return 1;
 		}
 	} else if (outputHeaderFormat != 0){
 		printf("kmer_main:storeOutputFeatureVector: wrong output header format: %d\n", outputHeaderFormat);
 		return 1;
 	}

	if (fprintf(outputFile,"\n") < 0){
		printf("kmer_main:storeOutputFeatureVector: Cannot store new line to a file!\n");
		return 1;
	}

	return 0;
}

