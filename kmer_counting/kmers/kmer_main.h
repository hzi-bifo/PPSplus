/*
 *  kmer_main.h
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

#ifndef KMER_MAIN_H_
#define KMER_MAIN_H_


int kmerMainComputeFeatures(char * inputFastaFilePath, char * outputFeatureFilePath, char outputFopenMode,
		int kmerFrom, int kmerTo, int reverseComplement, int removeRedundantFeatures, int normalization,
		int labelFormat, int outputFeatureVectorFormat, int outputFeatureVectorIndexOffset,
		int outputHeaderFormat, int backwardsCompatibility);

#endif /* KMER_MAIN_H_ */
