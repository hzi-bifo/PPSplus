/*
 *  kmer_seq.h
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

#ifndef KMER_SEQ_H_
#define KMER_SEQ_H_



#define INIT_SEQ_LEN 4096
#define INIT_SEQ_NAME_LEN 64

typedef void DNASeq;


DNASeq * seqCreate();

char * seqGetName(DNASeq * seq);

char * seqGetSeq(DNASeq * seq);

unsigned long seqGetNameLen(DNASeq * seq);

unsigned long seqGetSeqLen(DNASeq * seq);

int seqExtendName(DNASeq * seq, char c);

int seqExtendSeq(DNASeq * seq, char c);

void seqClear(DNASeq * seq);

void seqDestroy(DNASeq * seq);


#endif /* KMER_SEQ_H_ */
