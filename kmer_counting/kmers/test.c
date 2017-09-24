/**
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
 *
 * */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//#include <limits.h>

#include "kmer_counter.h"
#include "kmer_fasta.h"
#include "kmer_seq.h"
#include "kmer_fasta.h"
#include "char_list.h"
#include "kmer_features.h"
#include "kmer_main.h"

/* taken from "features.c" */
inline void seq2str(unsigned int sequence, int length, char *str){
	int i;
	static const char map[] = {'A','C','T','G'};
	/* Build query string */
	for (i=0;i<length;i++) {
		str[i] = map[sequence&0x3];
		sequence = sequence >> 2;
	}
}


/* return all possible kmers of length k from DNA alphabet from "features.c" */
char **dna_kmers(int k, int map){
	unsigned long n_kmers, i;
	if(map==0)
		n_kmers = pow(4, k);
	else if(map==1)
		n_kmers = pow(2, k);
	else
		{printf("\ndna_kmers: invalid map %d\n", map); exit(1);};
	char **kmers = (char **)malloc(n_kmers*sizeof(char*));
	if (kmers == NULL){
		printf("dna_kmers: cannot allocate kmers\n");
		exit(1);
	}
	for(i=0;i<n_kmers;i++){
		kmers[i] = malloc((k+1)*sizeof(char));
		if (kmers[i] == NULL){
			printf("dna_kmers: cannot allocate kmers[i]\n");
			exit(1);
		}
		if(map==0){
			seq2str(i, k, kmers[i]);
		} else if (map==1){
			//seq2str_RY(i, k, kmers[i]);
			printf("Map=1 not defined");
			exit(1);
		}
		kmers[i][k]='\0';
	}
	return(kmers);
}

void free_dna_kmers(char *** kmers, int k){
	unsigned long i;
	unsigned long count = pow(4,k);
	for (i=0; i<count; i++){
		free((*kmers)[i]);
	}
	free(*kmers);
	*kmers = NULL;
}

/**
 * Simple kmer counting function.
 * */
double countKmerStd(char * sequence, unsigned long seqLen, char * kmer, int kmerLen){
	double count = 0.0;
	unsigned long i;
	for (i=0; i<seqLen-kmerLen+1; i++){
		if (strncmp(&(sequence[i]), kmer, kmerLen) == 0){
			++count;
		}
	}
	return count;
}

/**
 * Unit test.
 * */
void kmerCounterUnitTest(){

	printf("kmerCounterUnitTest: start\n");

	//char * sequence = "ATTGCAGCCNNAGATTNCGCATCTATN";
	//unsigned long seqLen = 27;
	int kmerArrayLen = 6;
	int kmerArray[kmerArrayLen];
	kmerArray[0] = 1;
	kmerArray[1] = 2;
	kmerArray[2] = 3;
	kmerArray[3] = 4;
	kmerArray[4] = 5;
	kmerArray[5] = 6;
	//unsigned long featureArrayLen = 5460;
	unsigned long i,j,k;
	//unsigned long kmerOffsetArray[kmerArrayLen];
	//getKmerOffsetArray(kmerOffsetArray, kmerArray, kmerArrayLen);

	KmerCounter * kc = kmerCounterCreate(kmerArray, kmerArrayLen);
	if (kc == NULL){
		printf("kmerCounterUnitTest: Cannot create KmerCounter\n");
		exit(1);
	}

	KmerFeatureVector * featureVector = kmerFeatureVectorCreate(kmerCounterGetFeatureArrayLenRequired(kc));
	if (featureVector == NULL){
		printf("kmerCounterUnitTest: Cannot create featureVector\n");
		exit(1);
	}
	//double * featureArray = (double*)malloc(sizeof(double)*featureArrayLen);
	//if (featureArray == NULL){
	//	printf("kmerCounterUnitTest: Cannot allocate feature array!\n");
	//	exit(1);
	//}

	char *** kmerArrayQueryAll = (char***)malloc(sizeof(char**)*kmerArrayLen);
	if (kmerArrayQueryAll == NULL){
		printf("kmerCounterUnitTest: Cannot allocate kmerArrayQuery!\n");
		exit(1);
	}
	for (i=0; i<kmerArrayLen; i++){
		kmerArrayQueryAll[i] = dna_kmers(kmerArray[i], 0);
	}


	FastaHandler * fh = readBufferInit("./input/test_fasta3.fna");
	//FastaHandler * fh = readBufferInit("/Users/ivan/Documents/nobackup/fasta/inputTW2.fas");
	if (fh == NULL){
		printf("kmerCounterUnitTest: Cannot initialize read buffer\n");
		return;
	}
	DNASeq * seq =  seqCreate();
	while (readBufferHasNext(fh)){
		//printf(".");
		if (readBufferGetNextSeq(fh, seq) != 0){
			printf("kmerCounterUnitTest: Cannot read sequence!\n");
			break;
		}

		//count kmers
		if (kmerCounterCountKmers(kc, seq, featureVector) != 0){
			printf("kmerCounterUnitTest: count kmers returned non-zero value!\n");
		}


		int kmerArrayQueryLen;
		//char ** kmerArrayQuery;
		double u,v;
		for (i=0; i<kmerArrayLen; i++){

			//kmerArrayQuery = kmerArrayQueryAll[i];

			kmerArrayQueryLen = pow(4,kmerArray[i]);
			for (j=0; j<kmerArrayQueryLen; j++){
				char * kmerQuery = kmerArrayQueryAll[i][j];
				u = countKmerStd(seqGetSeq(seq), seqGetSeqLen(seq), kmerQuery, kmerArray[i]);
				v = kmerFeatureGetFeatureAt(featureVector, kmerCounterGetKmerIndex(kc, kmerQuery, kmerArray[i]));
				//v = featureArray[kmerCounterGetKmerIndex(kc, kmerQuery, kmerArray[i])];

				if (fabs(u - v) > 0.0001){
					printf("feature counts are not correct\n");
					printf("%s\n", seqGetSeq(seq));
					printf("%s: Std:%lf, %lf\n", kmerQuery, u,v);
					for (k=0; k<kmerFeatureGetFeatureArrayLen(featureVector); k++){
						printf("%d %f\n", (int)k, kmerFeatureGetFeatureAt(featureVector,k));
					}
					exit(1);
				}
			}
			//free_dna_kmers(&kmerArrayQuery, kmerArray[i]);
		}
		//printf("OK\n");
	}
	seqDestroy(seq);
	kmerFeatureVectorDestroy(featureVector);
	readBufferDestroy(fh);
	kmerCounterDestroy(kc);
	for (i=0; i<kmerArrayLen; i++){
		free_dna_kmers(&(kmerArrayQueryAll[i]), kmerArray[i]);
	}
	free(kmerArrayQueryAll);
	//free(featureArray);

	printf("kmerCounterUnitTest: finished\n");
}

/**
 * Speed Test.
 * */
void kmerCounterSpeedTest(){
	printf("kmerCounterSpeedTest: start\n");
	int kmerArrayLen = 3;
	//unsigned long featureArrayLen = 5376;
	int kmerArray[kmerArrayLen];
	kmerArray[0] = 4;
	kmerArray[1] = 5;
	kmerArray[2] = 6;

	//unsigned long kmerOffsetArray[kmerArrayLen];
	//getKmerOffsetArray(kmerOffsetArray, kmerArray, kmerArrayLen);

	KmerCounter * kc = kmerCounterCreate(kmerArray, kmerArrayLen);
	if (kc == NULL){
		printf("kmerCounterSpeedTest: Cannot create KmerCounter\n");
		exit(1);
	}
	KmerFeatureVector * featureVector = kmerFeatureVectorCreate(kmerCounterGetFeatureArrayLenRequired(kc));
	if (featureVector == NULL){
		printf("kmerCounterUnitTest: Cannot create featureVector\n");
		exit(1);
	}
	//double * featureArray = (double*)malloc(sizeof(double)*featureArrayLen);
	//if (featureArray == NULL){
	//	printf("kmerCounterSpeedTest: Cannot allocate feature array!\n");
	//	exit(1);
	//}

	//FastaHandler * fh = readBufferInit("test_fasta.fna");
	//FastaHandler * fh = readBufferInit("/Users/ivan/Documents/nobackup/fasta/inputTW2.fas");
	//FastaHandler * fh = readBufferInit("/Users/ivan/Documents/nobackup/fasta/TS29_contigs.fna");
	FastaHandler * fh = readBufferInit("/Users/ivan/Documents/nobackup/fasta/contigs.fna");
	if (fh == NULL){
		printf("kmerCounterSpeedTest: Cannot initialize read buffer\n");
		return;
	}
	DNASeq * seq =  seqCreate();
	while (readBufferHasNext(fh)){
		//printf(".");
		if (readBufferGetNextSeq(fh, seq) != 0){
			printf("kmerCounterSpeedTest: Cannot read sequence!\n");
			break;
		}
		//count kmers


		if (kmerCounterCountKmers(kc, seq, featureVector) != 0){
			printf("kmerCounterSpeedTest: count kmers returned non-zero value!\n");
		}
	}
	seqDestroy(seq);
	kmerFeatureVectorDestroy(featureVector);
	readBufferDestroy(fh);
	kmerCounterDestroy(kc);
	//free(featureArray);
	printf("kmerCounterSpeedTest: stop\n");
}



void writeSpeedPutcTest(){

	char * filePath = "/AM/metagenomic/work/projects/pPPS/tests/fasta2kmers/tmp.out";
	//char * filePath = "/Users/ivan/Documents/work/CC/workspace/kmers/output/test_fasta.txt"
	FILE * outputFile = fopen(filePath,"w");
	unsigned int i,j,k,size=1024;
	for (i=0; i<size; i++){
		for (j=0;j<size;j++){
			for (k=0; k<size; k++){
				if (putc ((int)((k % 20) + 48), outputFile) == EOF){
					printf("writeSpeedTest putc error!\n");
					return;
				}
			}
		}
	}
	fclose(outputFile);
}

void kmerCounterReverseComplementTest(){
	printf("kmerCounterReverseComplementTest: start\n");
	int kmerArrayLen = 2;
	int kmerArray[kmerArrayLen];
	kmerArray[0] = 1;
	kmerArray[1] = 2;
	unsigned long i, outMapArrayLen = 20;//kmerCounterGetFeatureArrayLenRequired(kc);
	unsigned long results[20] = {1,0,3,2,9,5,17,13,8,4,16,12,11,7,19,15,10,6,18,14};

	KmerCounter * kc = kmerCounterCreate(kmerArray, kmerArrayLen);
	unsigned long outMapArray[outMapArrayLen];
	if (kmerCounterComputeReverseComplementMap(kc, outMapArray, &outMapArrayLen) != 0){
		printf("kmerReverseComplementTest: (1) kmerCounterGetReverseComplementMap returned non-zero value!\n");
		return;
	}
	for (i=0; i<outMapArrayLen; i++){
		if (outMapArray[i] != results[i]){
			printf("kmerReverseComplementTest: wrong kmer at index:%lu reported:%lu correct%lu",
					i,outMapArray[i],results[i]);
			return;
		}
		//printf("%lu %lu\n",i,outMapArray[i]);
	}
	if (outMapArrayLen != kmerCounterGetFeatureArrayLenRequired(kc)){
		printf("kmerReverseComplementTest: inconsistent required features value(*) %lu %lu!\n",
				outMapArrayLen, kmerCounterGetFeatureArrayLenRequired(kc));
		return;
	}

	kmerCounterDestroy(kc);

	int kmerArrayLen2 = 3;
	int kmerArray2[kmerArrayLen2];
	kmerArray2[0] = 4;
	kmerArray2[1] = 5;
	kmerArray2[2] = 6;
	unsigned long outMapArrayLen2 = 5376;//kmerCounterGetFeatureArrayLenRequired(kc);

	KmerCounter * kc2 = kmerCounterCreate(kmerArray2, kmerArrayLen2);
	unsigned long outMapArray2[outMapArrayLen2];
	if (kmerCounterComputeReverseComplementMap(kc2, outMapArray2, &outMapArrayLen2) != 0){
		printf("kmerReverseComplementTest: (2) kmerCounterGetReverseComplementMap returned non-zero value!\n");
		return;
	}
	for (i=0; i<outMapArrayLen2; i++){
		if (outMapArray2[outMapArray2[i]] != i){
			printf("kmerReverseComplementTest: wrong mapping: %lu -> %lu", i, outMapArray2[i]);
			return;
		}
		//printf("%lu %lu\n",i,outMapArray2[i]);
	}
	kmerCounterDestroy(kc2);

	printf("kmerCounterReverseComplementTest: end\n");
}

void kmerMainTest(){
	/* main test */
	//char * inputFastaFilePath = "/Users/ivan/Documents/work/CC/workspace/kmers/test_fasta.fna";
	//char * inputFastaFilePath = "/Users/ivan/Documents/work/CC/workspace/kmers/input/test_fasta2.fna";
	//char * inputFastaFilePath = "/Users/ivan/Documents/work/CC/workspace/kmers/input/test_fasta3.fna";
	char * inputFastaFilePath = "/Users/ivan/Documents/work/CC/workspace/kmers/input/34034.1.fas";
	//char * inputFastaFilePath = "/Users/ivan/Documents/nobackup/fasta/contigs.fna";
	char * outputFeatureFilePath = "/Users/ivan/Documents/work/CC/workspace/kmers/output/test_fasta.txt";
	char outputFopenMode = 'a';
	int kmerFrom = 1;
	int kmerTo = 2;
	int reverseComplement = 1; //1
	int removeRedundantFeatures = 1;//1
	int normalization = 0; //...
	int labelFormat = 1; //0 or 1
	int outputFeatureVectorFormat = 1; //0 or 1
	int outputFeatureVectorIndexOffset = 1;
	int outputHeaderFormat = 1;//..
	int backwardsCompatibility = 1;//1

	int retVal = kmerMainComputeFeatures(inputFastaFilePath, outputFeatureFilePath, outputFopenMode, kmerFrom,
			kmerTo, reverseComplement, removeRedundantFeatures, normalization, labelFormat,
			outputFeatureVectorFormat, outputFeatureVectorIndexOffset, outputHeaderFormat,
			backwardsCompatibility);
	if (retVal != 0){
		printf("kmerMainComputeFeatures: returned non-zero value: %d\n", retVal);
	}

}


void printHelp(char * programName){
	printf("Converts a fasta file to kmer feature vectors (author: Ivan Gregor version: 11.1.2012)\n\n");
	printf("Disclaimer:\n\t This program is provided WITHOUT ANY WARRANTY.\n\n");

	printf("Usage:\n\t%s -i inputFastaFilePath -f outputFeatureFilePath [options]\n\n", programName);

	printf("All options:\n");
	printf("\t i: inputFastaFilePath: path to a fasta file\n");
	printf("\t f: outputFeatureFilePath: path to an output file\n");
	printf("\t a: outputFopenMode: w~for writing (original file is erased), "
			      "a~append to an existing file or create a new file (default: w)\n");
	printf("\t j: kmerFrom: is less or equal to kmerTo\n");
	printf("\t k: kmerTo: (e.g. if kmerFrom~4,kmerTo~6 then kmers 4,5,6 are computed) max. ~15\n");
	printf("\t r: reverseComplement: 1 ~ consider also reverse complement "
			      "(results in half of the features) else 0 (default: 1)\n");
	printf("\t R: removeRedundantFeatures: 1~remove redundant features if reverse complement~1 else 0 "
			      "(default 1)\n");
	printf("\t n: normalization: 0~ no normalization, 1 ~ sequence length (default: 0)\n");
	printf("\t l: labelFormat: 0~ no labels, 1~ label defined in seq name after label:, "
			      "2~ label equal to `1` (default: 0)\n");
	printf("\t s: outputFeatureVectorFormat: 0 ~ plain, 1 ~ sparse (default: 1)\n");
	printf("\t o: outputFeatureVectorIndexOffset: offset of indices in the output feature file "
			      "(default: 0)\n");
	printf("\t h: outputHeaderFormat: 1 ~ #seq name appended at the end of the line, "
			      "0 ~ no header (default: 0)\n");
	printf("\t b: backwardsCompatibility 1 ~ order of the features correspond to the order "
			      "in the original implementation, 0 ~ order of the features is according to "
			      "this implementation (default: 0)\n");
}


int main(int argc, char *argv[]){


	if (argc == 0){
		//char * file = "/AM/metagenomic/work/projects/pPPS/data/HumanGut/working/contigs.fna";
		//clist tests (char_list.h)
		//clistUnitTest();
		//clistSpeedTest();

		//readBuffer FASTA tests (kmer_fasta.h)
		//readBufferUnitTest();
		//readBufferSpeedTest();

		//kmerCounterUnitTest();
		//kmerCounterSpeedTest();



		//kmerCounterReverseComplementTest();

		kmerMainTest();
		//
		//writeSpeedPutcTest();
		return 0;
	} else {
		//parse parameters
		char * inputFastaFilePath = NULL;
		char * outputFeatureFilePath = NULL;
		char outputFopenMode = 'w';
		int kmerFrom = 0;
		int kmerTo = 0;
		int reverseComplement = 1;
		int removeRedundantFeatures = 1;
		int normalization = 0;
		int labelFormat = 0;
		int outputFeatureVectorFormat = 1;
		int outputFeatureVectorIndexOffset = 0;
		int outputHeaderFormat = 0;
		int backwardsCompatibility = 0;

		if (argc > 4){
			int i;
			for (i=1; i<argc-1;){
				if (argv[i][0] == '-' && strlen(argv[i]) > 1){
					switch (argv[i][1]){
					case 'i': ++i; inputFastaFilePath = argv[i]; break;
					case 'f': ++i; outputFeatureFilePath = argv[i]; break;
					case 'a': ++i; outputFopenMode = argv[i][0]; break;
					case 'j': ++i; kmerFrom = atoi(argv[i]); break;
					case 'k': ++i; kmerTo = atoi(argv[i]); break;
					case 'r': ++i; reverseComplement = atoi(argv[i]); break;
					case 'R': ++i; removeRedundantFeatures = atoi(argv[i]); break;
					case 'n': ++i; normalization = atoi(argv[i]); break;
					case 'l': ++i; labelFormat = atoi(argv[i]); break;
					case 's': ++i; outputFeatureVectorFormat = atoi(argv[i]); break;
					case 'o': ++i; outputFeatureVectorIndexOffset = atoi(argv[i]); break;
					case 'h': ++i; outputHeaderFormat = atoi(argv[i]); break;
					case 'b': ++i; backwardsCompatibility = atoi(argv[i]); break;
					}
				}
				++i;
			}

			int retVal =  kmerMainComputeFeatures(inputFastaFilePath, outputFeatureFilePath, outputFopenMode,
					kmerFrom, kmerTo, reverseComplement, removeRedundantFeatures, normalization, labelFormat,
					outputFeatureVectorFormat, outputFeatureVectorIndexOffset, outputHeaderFormat,
					backwardsCompatibility);
			if (retVal != 0){
				printHelp(argv[0]);
			}
			return retVal;
		} else {
			printHelp(argv[0]);
			return -1;
		}
	}
	//printf("main end\n");
	//return 0;
}
