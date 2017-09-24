
/*
 *  kmer_seq.c
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
//#include <string.h>

#import "kmer_seq.h"
#import "char_list.h"


typedef struct {
	CList * name;
	CList * seq;
} Seq;


/**
 * Initialize the sequence structure.
 * */
inline DNASeq * seqCreate(){

	Seq * s = (Seq *)malloc(sizeof(Seq));
	if (s == NULL){
		return NULL;
	}
	s->name = clistCreate(INIT_SEQ_NAME_LEN);
	if (s->name == NULL){
		free(s);
		return NULL;
	}
	s->seq = clistCreate(INIT_SEQ_LEN);
	if (s->seq == NULL){
		free(s->name);
		free(s);
		return NULL;
	}
	return (DNASeq *)s;
}


inline char * seqGetName(DNASeq * seq){
	if (seq == NULL){
		return NULL;
	}
	return clistGetString(((Seq *)seq)->name);
}

inline char * seqGetSeq(DNASeq * seq){
	if (seq == NULL){
		return NULL;
	}
	return clistGetString(((Seq *)seq)->seq);
}

inline unsigned long seqGetNameLen(DNASeq * seq){
	if (seq == NULL){
		return 0;
	}
	return clistSize(((Seq *)seq)->name);
}

inline unsigned long seqGetSeqLen(DNASeq * seq){
	if (seq == NULL){
		return 0;
	}
	return clistSize(((Seq *)seq)->seq);
}

/**
 * @return 0 if successful
 * */
inline int seqExtendName(DNASeq * seq, char c){
	if (seq == NULL){
		return 1;
	}
	return clistAdd(((Seq *)seq)->name, c);
}

/**
 * @return 0 if successful
 * */
inline int seqExtendSeq(DNASeq * seq, char c){
	if (seq == NULL){
		return 1;
	}
	return clistAdd(((Seq *)seq)->seq, c);
}

inline void seqClear(DNASeq * seq){
	if (seq != NULL){
		clistClear(((Seq *)seq)->name);
		clistClear(((Seq *)seq)->seq);
	}
}

inline void seqDestroy(DNASeq * seq){
	if (seq != NULL){
		clistDestroy(((Seq *)seq)->name);
		clistDestroy(((Seq *)seq)->seq);
		free(seq);
	}
}
