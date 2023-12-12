/******************************************************************************
 *  The Janelia Farm Research Campus Software Copyright 1.1
 *
 *  Copyright (c) 2010, Howard Hughes Medical Institute, All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *     1. Redistributions of source code must retain the above copyright
 *        notice, this list of conditions and the following disclaimer.
 *     2. Redistributions in binary form must reproduce the above copyright
 *        notice, this list of conditions and the following disclaimer in the
 *        documentation and/or other materials provided with the distribution.
 *     3. Neither the name of the Howard Hughes Medical Institute nor the
 *        names of its contributors may be used to endorse or promote products
 *        derived from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 *  TO, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, NON-INFRINGEMENT, OR FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 *  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; REASONABLE ROYALTIES; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 *  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
 *  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 ******************************************************************************/

/*
 *  Written by Michael Farrar, 2010.
 *  Please send bug reports and/or suggestions to farrar.michael@gmail.com.
 */

#ifndef INCLUDE_DUMPNCBI_H
#define INCLUDE_DUMPNCBI_H
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdint.h>
using namespace std;

class SingleDBReader {
public:
	SingleDBReader();
	~SingleDBReader();
	bool open(const string& fileName);
	int getSeq(char* res, int maxLen);
	void close();
	void load();
	void unload();
	string getSeqName(int index);
	string getSeqString(int index);
	// returns SBN of the last processed sequence
	long long getSBN();
	int totalSeqs();
	bool isDnaDatabase();
	int getMaxSequenceLength();
	long long getTotalSequencesLength();
private:
	FILE *fpin;
	FILE *fpsq;
	FILE *fphr;
	uint32_t num_seq;
	int seqIndex;
	uint32_t alpha_type;
	uint32_t seq_start, seq_end;
	uint32_t amb_start;
	uint32_t* arr_seq;
	uint32_t* arr_amb;
	uint32_t* arr_hdr;

	uint32_t hdr_off;
	uint32_t seq_off;
	uint32_t amb_off;

	uint32_t max_seq;

	// total residues count
	uint64_t total_res;

	char* file;
	char* hdr_file;
};

class BlastDBReader {
public:
	BlastDBReader();
	~BlastDBReader();
	bool open(const string& fileName);
	void close();
	// gets next sequence from database and writes it in res
	// if its length exceeds maxLen, no characters will be read and (maxLen + 1) will be returned
	int getSeq(char* res, int maxLen);
	// returns global index of the next sequence (not processed so far)
	int getSeqIndex();
	// returns SBN of the last processed sequence
	long long getSBN();
	string getSeqName(int index);
	string getSeqString(int index);
	bool isDnaDatabase();
	int getBasesNumber();
	int getMaxSequenceLength();
	long long getTotalSequencesLength();
private:
	bool getBaseFiles(const string& fileName, vector<string>* result);
	vector<string> baseFiles;
	int openedBaseFileIndex;
	vector<SingleDBReader> singleDBReaders;
	int seqIndex;
	bool isDna;
	int basesNumber;
};



#endif
