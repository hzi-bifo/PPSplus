/*
 *  kmer_fasta.c
 *
 *  Created on: Dec 12, 2011
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

#include "kmer_fasta.h"
#include "kmer_seq.h"

typedef struct {
	FILE *fr;
	int lastChar;
} FileHandler;



/**
 * Initializes the fasta file handler.
 * */
inline FastaHandler * readBufferInit(char * inputFile){

	if (inputFile == NULL){
		printf("readBuffer: Input fasta file is not specified!\n");
		return NULL;
	}

	FileHandler * fh = (FileHandler *)malloc(sizeof(FileHandler));
	if (fh == NULL){
		printf("readBuffer: Can`t allocate FileHandler!\n");
		return NULL;
	}

	/* open fasta file for reading */
	fh->fr = fopen(inputFile,"r");
	if (fh->fr == NULL){
		printf("readBuffer: Can`t open fasta file: %s\n", inputFile);
		return NULL;
	}

	/* move to the first ">" */
	while ((fh->lastChar = getc(fh->fr)) != EOF){
		if (((char)(fh->lastChar)) == '>'){
			break;
		}
	}

	return (FastaHandler *)fh;
}

/**
 * @return 0 if empty
 * */
inline int readBufferHasNext(FastaHandler * fastaHandler){
	return (((FileHandler *)fastaHandler)->lastChar == EOF)?0:1;
}

/**
 * @param seq: DNA seq or NULL
 * @return: DNA seq or NULL
 * */
inline int readBufferGetNextSeq(FastaHandler * fastaHandler, DNASeq * seq){
	FileHandler * fh = (FileHandler *)fastaHandler;
	if (fh == NULL){
		printf("readBufferGetNextSeq: read buffer is NULL!\n");
		return 1;
	}
	if (fh->lastChar == EOF){
		printf("readBufferGetNextSeq: read buffer is already empty\n");
		return 1;
	}
	if (seq == NULL){
		printf("readBufferGetNextSeq: seq is NULL\n");
		return 1;
	}
	seqClear(seq);

	// read seq name
	char c;
	while ((fh->lastChar = getc(fh->fr)) != EOF){
		c = (char)fh->lastChar;
		if ((c != '\n') && (c != '\r')){//&& c != ' '&& c != '\t' //remove \r !!!
			if (seqExtendName(seq, c) != 0){
				printf("readBufferGetNextSeq: Cannot extend sequence name!\n");
				seqClear(seq);
				return 1;
			}
		} else {
			break;
		}
	}

	/*if (fh->lastChar == EOF){
		printf("readBufferGetNextSeq: EOF at the end of sequence name!\n");
		seqClear(seq);
		return 1;
	}

	// skip seq description
	if ((char)fh->lastChar == '\t'){//(char)fh->lastChar == ' ' ||
		while ((fh->lastChar = getc(fh->fr)) != EOF){
			if ((char)fh->lastChar == '\n' || (char)fh->lastChar == '\r'){
				break;
			}
		}
	}*/

	if (fh->lastChar == EOF){
		printf("readBufferGetNextSeq: EOF at the end of sequence name within comment!\n");
		seqClear(seq);
		return 1;
	}

	// skip possible comment
	short comment = 0;
	while ((fh->lastChar = getc(fh->fr)) != EOF){
		c = (char)fh->lastChar;
		if (c == '\n' || c == '\r'){//remove \r !!!
			comment = 0;
			continue;
		}
		if (c == ';'){
			comment = 1;
			continue;
		}
		if (comment == 1){
			continue;
		}
		break;
	}

	if ((char)fh->lastChar == EOF){
		printf("readBufferGetNextSeq: EOF appeared while comments were being skipped!\n");
		seqClear(seq);
		return 1;
	}

	char nucleotide;
	//unsigned long nonDNAChar = 0;
	//unsigned long DNAChar = 0;

	// read sequence
	do {
		c = (char)fh->lastChar;


		if (c == '>'){
			break;
		}


		if (c == '\n' || c == '\r' || c == ' ' || c == '\t'){
			continue;
		}


		if (c == 'A' || c == 'a'){
			nucleotide = 'A';
			//++DNAChar;
		} else if (c == 'T' || c == 't'){
			nucleotide = 'T';
			//++DNAChar;
		} else if (c == 'G' || c == 'g'){
			nucleotide = 'G';
			//++DNAChar;
		} else if (c == 'C' || c == 'c'){
			nucleotide = 'C';
			//++DNAChar;
		} else {
			nucleotide = c;
			//nucleotide = '\0';
			//++nonDNAChar;
		}


		if (seqExtendSeq(seq, nucleotide) != 0){
			printf("readBufferGetNextSeq: Can`t extend sequence!\n");
			seqClear(seq);
			return 1;
		}


		/*if (nucleotide != '\0'){
			if (nonDNAChar > 0){//replace a stretch of nonDNA chars by one "N"
				if (seqExtendSeq(*seq, 'N') != 0){
					seqDestroy(*seq);
					*seq = NULL;
					printf("kf_getNextSeq: can`t extend sequence seq with N\n");
					return 1;
				}
				nonDNAChar = 0;
			}
			if (seqExtendSeq(*seq, nucleotide) != 0){
				seqDestroy(*seq);
				*seq = NULL;
				printf("kf_getNextSeq: can`t extend sequence seq\n");
				return 1;
			}
		}*/
	} while ((fh->lastChar = getc(fh->fr)) != EOF);

	return 0;
}


/**
 * Destroys the fasta handler.
 * */
inline void readBufferDestroy(FastaHandler * fastaHandler){
	FileHandler * fileHandler = (FileHandler *)fastaHandler;
	if (fileHandler != NULL){
		if (fileHandler->fr != NULL){
			fclose(fileHandler->fr);
		}
		free(fileHandler);
	}
}


void readBufferUnitTest(){
	int i,size = 5;
	char * names[size];
	char * seqArray[size];
	names[0] = "seq_name1";
	names[1] = "seq_name2";
	names[2] = "seq_NAME3";
	names[3] = "seq_name_4";
	names[4] = "seqName5";
	seqArray[0] = "ATGCATGCCTCTTATTCTTTTCCCCCCCCCCCCCCAAAAAAAAA";
	seqArray[1] = "AAAATTTTTTTCCCCCGGGGGCCCCTAAAATT";
	seqArray[2] = "AAATTTGGGCGTCATGCTATGDCCCCDDNNXXXXATXXGTCMMMTCGA";
	seqArray[3] = "ATGCTGANTGCTATCTTTTAARETGNXCTGTTTTTAAAAAAGGGGGGCCCCCCGTCATGCCTG";
	seqArray[4] = "NCTGTCHJKGTGNMTJHATGCAATATCTGCTGA";

	FastaHandler * fh = readBufferInit("./input/test_fasta.fna");
	if (fh == NULL){
		printf("readBufferUnitTest: Cannot initialize the read buffer");
		return;
	}
	DNASeq * seq =  seqCreate();
	int success = 1;
	for (i=0; readBufferHasNext(fh); i++){
		if (readBufferGetNextSeq(fh, seq) != 0){
			printf("readBufferUnitTest: Cannot read sequence!\n");
			success = 0;
			break;
		}
		if (strcmp(names[i], seqGetName(seq)) != 0){
			printf("readBufferUnitTest: Sequence name is different!\n");
			printf("%s\n", names[i]);
			printf("%s\n", seqGetName(seq));
			success = 0;
			break;
		}
		if (strcmp(seqArray[i], seqGetSeq(seq)) != 0){
			printf("readBufferUnitTest: Sequences are different!\n");
			printf("%s\n", seqArray[i]);
			printf("%s\n", seqGetSeq(seq));
			success = 0;
			break;
		}
	}
	seqDestroy(seq);
	readBufferDestroy(fh);
	if (success){
		printf("readBufferUnitTest: success\n");
	}
}


void readBufferSpeedTest(){
	printf("readBufferSpeedTest: start\n");
	//FastaHandler * fh = readBufferInit("/Users/ivan/Documents/nobackup/fasta/inputTW2.fas");
	FastaHandler * fh = readBufferInit("/Users/ivan/Documents/nobackup/fasta/contigs.fna");
	if (fh == NULL){
		printf("readBufferSpeedTest: Cannot initialize read buffer!\n");
		return;
	}
	DNASeq * seq =  seqCreate();
	while (readBufferHasNext(fh)){
		if (readBufferGetNextSeq(fh, seq) != 0){
			printf("readBufferSpeedTest: Cannot read sequence!\n");
			break;
		}
	}
	seqDestroy(seq);
	readBufferDestroy(fh);
	printf("readBufferSpeedTest: end\n");
}

