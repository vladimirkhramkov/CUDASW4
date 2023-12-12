/******************************************************************************
 *	The Janelia Farm Research Campus Software Copyright 1.1
 *
 *	Copyright (c) 2010, Howard Hughes Medical Institute, All rights reserved.
 *
 *	Redistribution and use in source and binary forms, with or without
 *	modification, are permitted provided that the following conditions are met:
 *
 *		 1. Redistributions of source code must retain the above copyright
 *				notice, this list of conditions and the following disclaimer.
 *		 2. Redistributions in binary form must reproduce the above copyright
 *				notice, this list of conditions and the following disclaimer in the
 *				documentation and/or other materials provided with the distribution.
 *		 3. Neither the name of the Howard Hughes Medical Institute nor the
 *				names of its contributors may be used to endorse or promote products
 *				derived from this software without specific prior written permission.
 *
 *	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *	"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 *	TO, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, NON-INFRINGEMENT, OR FITNESS
 *	FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *	OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *	SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 *	TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *	PROFITS; REASONABLE ROYALTIES; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *	AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 *	OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
 *	USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 ******************************************************************************/

/*
 *	Written by Michael Farrar, 2010.
 *	Please send bug reports and/or suggestions to farrar.michael@gmail.com.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <errno.h>
#ifdef WIN32
#include "dumpmisc.h"
#else
#include <stdint.h>
#include <endian.h>
#endif

#include "dumpseq.h"
#include "dumphdr.h"
#include "dumpncbi.h"

#define DNA_DATABASE		 0
#define PROTEIN_DATABASE 1

#define DATABASE_VERSION 4

#include <string>
#include <sstream>

using namespace std;

int letcnt[300];

static void error(char *format, ...)
{
	va_list	 argptr;

	va_start(argptr, format);
	vfprintf(stderr, format, argptr);
	va_end(argptr);

	exit(255);
}

int BlastDBReader::getBasesNumber() {
	return basesNumber;
}

int BlastDBReader::getMaxSequenceLength() {
	int maxSequenceLength = 0;
	for (int i = 0; i < singleDBReaders.size(); ++i)
		maxSequenceLength = max(maxSequenceLength, singleDBReaders[i].getMaxSequenceLength());
	return maxSequenceLength;
}

long long BlastDBReader::getTotalSequencesLength() {
	long long totalSequencesLength = 0;
	for (int i = 0; i < singleDBReaders.size(); ++i)
		totalSequencesLength += singleDBReaders[i].getTotalSequencesLength();
	return totalSequencesLength;
}

long long BlastDBReader::getSBN() {
	return singleDBReaders[openedBaseFileIndex].getSBN();
}

bool BlastDBReader::getBaseFiles(const string& fileName, vector<string>* result) {
cout << "BF: " << fileName << endl;
        vector<string> res;
        string filePath = fileName;
        while (filePath.size() && filePath[filePath.size() - 1] != '/') filePath.erase(filePath.size()-1, 1);
        if (!filePath.size()) return false;
        ifstream sqFile((fileName + ".psq").c_str());
        if (!sqFile.is_open()) {
                sqFile.open((fileName + ".nsq").c_str(), ios_base::in);
        }
        if (sqFile.is_open()) {
                res.push_back(fileName);
                sqFile.close();
        } else {
                ifstream alFile((fileName + ".nal").c_str());
                if (!alFile.is_open()) {
                        alFile.open((fileName + ".pal").c_str(), ios_base::in);
                }
                if (alFile.is_open()) {
                        string s;
                        while (getline(alFile, s)) {
                                stringstream ss(s);
                                string tmp; ss >> tmp;
                                if (tmp == "DBLIST") {
                                        while (ss >> tmp) {
                                                if (tmp[0] == '"') {
                                                        tmp = tmp.substr(1, tmp.length() - 2);
                                                }
                                                vector<string> tempResult;
                                                if (!getBaseFiles(filePath + tmp, &tempResult)) return false;
                                                res.insert(res.end(), tempResult.begin(), tempResult.end());
                                        }
                                        break;
                                }

                        }
                        alFile.close();
                }
        }
        result->insert(result->end(), res.begin(), res.end());
        return true;
}

bool BlastDBReader::open(const string& fileName) {
	baseFiles.clear();
	seqIndex = 0;
	openedBaseFileIndex = 0;
	stringstream ss(fileName);

	string baseFileName;
	basesNumber = 0;
	while (ss >> baseFileName) {
		basesNumber++;
		if (!getBaseFiles(baseFileName, &baseFiles)) return false;
	}
	singleDBReaders = vector<SingleDBReader>(baseFiles.size());
	for (int i = 0; i < singleDBReaders.size(); ++i) {
		if (!singleDBReaders[i].open(baseFiles[i])) return false;
	}
	singleDBReaders[0].load();
	return true;
}

BlastDBReader::BlastDBReader() {}

BlastDBReader::~BlastDBReader() {}

int BlastDBReader::getSeqIndex() {
	return seqIndex;
}

void BlastDBReader::close() {
	for (int i = 0; i < singleDBReaders.size(); ++i) singleDBReaders[i].close();
	baseFiles.clear();
}

int BlastDBReader::getSeq(char* res, int maxLen) {
	if (openedBaseFileIndex == baseFiles.size()) return 0;
	++seqIndex;
	int length = singleDBReaders[openedBaseFileIndex].getSeq(res, maxLen);
	if (length != 0) return length;
	++openedBaseFileIndex;
	if (openedBaseFileIndex == baseFiles.size()) return 0;
	singleDBReaders[openedBaseFileIndex - 1].unload();
	singleDBReaders[openedBaseFileIndex].load();
	return singleDBReaders[openedBaseFileIndex].getSeq(res, maxLen);
}

SingleDBReader::SingleDBReader() {
	file = 0;
	hdr_file = 0;
}
SingleDBReader::~SingleDBReader() {}

string SingleDBReader::getSeqName(int index) {
	fseek(fpin, hdr_off + sizeof(uint32_t) * index, SEEK_SET);
	uint32_t hdr_start, hdr_end;
	fread(&hdr_start, sizeof(uint32_t), 1, fpin);
	hdr_start = htobe32(hdr_start);
	fread(&hdr_end, sizeof(uint32_t), 1, fpin);
	hdr_end = htobe32(hdr_end);
	StaticString::visibleString = "";
	parseBlastHeader(fphr, hdr_start, hdr_end);
	return StaticString::visibleString;
}

void SingleDBReader::load() {
	//cout << "start processing SEQUENCE FILE" << endl;
	if (file) error("file is already loaded\n");
	fseek(fpsq, 0, SEEK_END);
	int fileSize = ftell(fpsq);
	file = (char*)malloc(sizeof(char)*fileSize);
	fseek(fpsq, 0, SEEK_SET);
	int readedcnt = fread(file, sizeof(char), fileSize, fpsq);
	//cout << "end processing SEQUENCE FILE" << endl;

	//cout << "start processing HEADER FILE" << endl;
	if (hdr_file) error("header file is already loaded\n");
	fseek(fphr, 0, SEEK_END);
	fileSize = ftell(fphr);
	hdr_file = (char*)malloc(sizeof(char)*fileSize);
	fseek(fphr, 0, SEEK_SET);
	readedcnt = fread(hdr_file, sizeof(char), fileSize, fphr);
	//cout << "end processing HEADER FILE" << endl;
}

void SingleDBReader::unload() {
	if (!file) error("file is not loaded\n");
	free(file), file = 0;
	if (!hdr_file) error("header file is not loaded\n");
	free(hdr_file), hdr_file = 0;
}

void SingleDBReader::close() {
	if (file) free(file), file = 0;
	if (fpin != NULL) fclose(fpin);
	if (fpsq != NULL) fclose(fpsq);
	if (fphr != NULL) fclose(fphr);
	if (arr_seq) free(arr_seq);
	if (arr_amb) free(arr_amb);
	if (arr_hdr) free(arr_hdr);
}

bool SingleDBReader::open(const string& fileName) {
cout << "FN: " << fileName << "\n";
	fpin = fpsq = fphr = NULL;
	arr_seq = arr_amb = arr_hdr = 0;
	seqIndex = 0;



	fpin = fpsq = fphr = NULL;


	int len;

	char *time = NULL;
	char *title = NULL;

	uint32_t info[4];

	uint32_t version;

	if ((fpin = fopen((fileName + ".pin").c_str(), "rb")) != NULL) {

		/* open up the sequence file */
		if ((fpsq = fopen((fileName + ".psq").c_str(), "rb")) == NULL) {

cout << "TTT 1\n";
			return false;
		}

		/* open up the header file */
		if ((fphr = fopen((fileName + ".phr").c_str(), "rb")) == NULL) {
cout << "TTT 2\n";
			return false;
		}
	}

	/* try to open a dna database */
	if (fpin == NULL) {
		char* test = (char*)malloc(500);
		strcpy(test, (fileName + ".nin").c_str());
		if ((fpin = fopen(test, "rb")) == NULL) {
cout << "TTT 4\n";
			return false;
		}

		/* open up the sequence file */
		if ((fpsq = fopen((fileName + ".nsq").c_str(), "rb")) == NULL) {
cout << "TTT 5\n" << endl;
			return false;
		}

		/* open up the header file */
		if ((fphr = fopen((fileName + ".nhr").c_str(), "rb")) == NULL) {
cout << "TTT 6\n" << endl;
			return false;
		}
	}

	/* read the head of the index file */
	if (fread(&info[0], sizeof(uint32_t), 3, fpin) != 3) {
		error("Reading version, alpha, length error %d\n", ferror(fpin));
	}

	version = htobe32(info[0]);
	if (version != DATABASE_VERSION) {
		error("Unsupported version %d\n", htobe32(info[0]));
	}

	alpha_type = htobe32(info[1]);
	if (alpha_type != DNA_DATABASE && alpha_type != PROTEIN_DATABASE) {
		error("Unknown database type %d\n", alpha_type);
	}

	/* read the database title */
	len = htobe32(info[2]);
	title = (char*)malloc(len+1);
	if (title == NULL) {
		error("Malloc title %d error\n", len+1);
	}
	if (fread(title, sizeof(char), len, fpin) != len) {
		error("Read title %d error %d\n", len, ferror(fpin));
	}
	title[len] = 0;

	/* read the database time stamp */
	if (fread(&info[0], sizeof(uint32_t), 1, fpin) != 1) {
		error("Read length error %d\n", ferror(fpin));
	}

	len = htobe32(info[0]);
	time = (char*)malloc(len+1);
	if (time == NULL) {
		error("Malloc timestamp %d error\n", len+1);
	}
	if (fread(time, sizeof(char), len, fpin) != len) {
		error("Read timestamp %d error %d\n", len, ferror(fpin));
	}
	time[len] = 0;

	/* read in database stats */
	if (fread(&info[0], sizeof(uint32_t), 4, fpin) != 4) {
		error("Read stats error %d\n", ferror(fpin));
	}
	num_seq	 = htobe32(info[0]);
	total_res = *(uint64_t *)(info+1);
	max_seq	 = htobe32(info[3]);

	/* save the offsets of the header and sequence indexes */
	hdr_off = ftell(fpin);
	seq_off = hdr_off + sizeof(uint32_t) * (num_seq + 1);
	amb_off = seq_off + sizeof(uint32_t) * (num_seq + 1);

	/* read in the first header offset */
	if (fseek(fpin, hdr_off, SEEK_SET) != 0) {
		error("Seek header off 0x%08X io error %d\n", hdr_off, ferror(fpin));
	}
	arr_hdr = (uint32_t*)malloc(sizeof(uint32_t) * (num_seq + 1));
	fread(arr_hdr, sizeof(uint32_t), num_seq + 1, fpin);
	/* read in the first sequence offset */
	if (fseek(fpin, seq_off, SEEK_SET) != 0) {
		error("Seek seq off 0x%08X io error %d\n", seq_off, ferror(fpin));
	}
	if (fread(&seq_start, sizeof(uint32_t), 1, fpin) != 1) {
		error("Read seq offset error %d\n", ferror(fpin));
	}
	seq_start = htobe32(seq_start);
	seq_off += sizeof(uint32_t);

	if (alpha_type == DNA_DATABASE) {
		arr_amb = (uint32_t*)malloc(sizeof(uint32_t) * (num_seq + 1));
		if (fseek(fpin, amb_off, SEEK_SET) != 0) {
			error("Seek amb off 0x%08X io error %d\n", amb_off, ferror(fpin));
		}

		fseek(fpin, amb_off, SEEK_SET);
		fread(arr_amb, sizeof(uint32_t), num_seq, fpin);
	}
	arr_seq = (uint32_t*)malloc(sizeof(uint32_t) * (num_seq + 1));
	if (fseek(fpin, seq_off, SEEK_SET) != 0) {
		error("Seek next seq off 0x%08X io error %d\n", seq_off, ferror(fpin));
	}

	int count = fread(arr_seq, sizeof(uint32_t), num_seq + 1, fpin);

	init_LUT();


	if (time != NULL) free(time);
	if (title != NULL) free(title);

	return true;
}

int SingleDBReader::getSeq(char* res, int maxLen) {
	if (seqIndex == num_seq) return 0;
		/* read in the start of the ambiguity table if dna */
		if (alpha_type == DNA_DATABASE) {
			amb_start = arr_amb[seqIndex];
			amb_start = htobe32(amb_start);
		}

		seq_end = arr_seq[seqIndex];
		seq_end = htobe32(seq_end);


		/* print out the sequence */
		int length;
		if (alpha_type == DNA_DATABASE) {
			length = printDnaSequence (file, seq_start, seq_end, amb_start, res, maxLen);
		} else {
			length = printProteinSequence(file, seq_start, seq_end, res, maxLen);
		}

		seq_start = seq_end;
		seqIndex++;
		return length;
}

string BlastDBReader::getSeqName(int index) {
	for (int i = 0; i < singleDBReaders.size(); ++i) {
		if (index < singleDBReaders[i].totalSeqs()) return singleDBReaders[i].getSeqName(index);
		else index -= singleDBReaders[i].totalSeqs();
	}
	return "";
}

string SingleDBReader::getSeqString(int index) {
	fseek(fpin, seq_off + (index - 1) * sizeof(uint32_t), SEEK_SET);
	uint32_t seqStart, seqEnd;
	fread(&seqStart, sizeof(uint32_t), 1, fpin);
	fread(&seqEnd, sizeof(uint32_t), 1, fpin);
	seqStart = htobe32(seqStart);
	seqEnd = htobe32(seqEnd);
	if (alpha_type == DNA_DATABASE) {
		fseek(fpin, amb_off + (index) * sizeof(uint32_t), SEEK_SET);
		uint32_t ambStart;
		fread(&ambStart, sizeof(uint32_t), 1, fpin);
		ambStart = htobe32(ambStart);
		return getDnaSequence(fpsq, seqStart, seqEnd, ambStart);
	} else {
		return getProteinSequence(fpsq, seqStart, seqEnd);
	}
}

int SingleDBReader::totalSeqs() {
	return num_seq;
}

string BlastDBReader::getSeqString(int index) {
	for (int i = 0; i < singleDBReaders.size(); ++i) {
		if (index < singleDBReaders[i].totalSeqs()) return singleDBReaders[i].getSeqString(index);
		else index -= singleDBReaders[i].totalSeqs();
	}
	return "";
}

bool SingleDBReader::isDnaDatabase() {
	return alpha_type == DNA_DATABASE;
}

bool BlastDBReader::isDnaDatabase() {
	return singleDBReaders[0].isDnaDatabase();
}

long long SingleDBReader::getSBN() {
	uint32_t hdr_start, hdr_end;

	hdr_start = arr_hdr[seqIndex - 1];
	hdr_start = htobe32(hdr_start);

	hdr_end = arr_hdr[seqIndex];
	hdr_end = htobe32(hdr_end);

	for (int i = hdr_start; i < hdr_end - 4; ++i) {
		if (hdr_file[i] == 'S') {
			if (hdr_file[i + 1] == 'B' && hdr_file[i + 2] == 'N') {
				long long result = 0;
				int j = i + 4;
				//string stringResult = "";
				while (j < hdr_end && ((hdr_file[j] >= '0' && hdr_file[j] <= '9') || (hdr_file[j] >= 'A' && hdr_file[j] <= 'Z'))) {
					result = result * 36 + hdr_file[j] - (hdr_file[j] <= '9' ? '0' : 'A' - 10);
					//stringResult += hdr_file[j];
					++j;
				}
				//cout << stringResult << endl;
				return result;
			}
		}
	}
	printf("ERR: sid: %d, hstart: %d, hend: %d\n", seqIndex, hdr_start, hdr_end);
	return -1;
}

int SingleDBReader::getMaxSequenceLength() {
	return max_seq;
}

long long SingleDBReader::getTotalSequencesLength() {
	return total_res;
}
