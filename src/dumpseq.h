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

#ifndef INCLUDE_DUMPSEQ_H
#define INCLUDE_DUMPSEQ_H

#ifdef WIN32
#include "dumpmisc.h"
#else
#include <stdint.h>
#endif
#include <string>
using namespace std;



void init_LUT();
int printDnaSequence (char* file, uint32_t start, uint32_t end, uint32_t amb_start, char* res, int maxLen);
int printProteinSequence(char* file, uint32_t start, uint32_t end, char* res, int maxLen);
string getDnaSequence(FILE* fpsq, uint32_t start, uint32_t end, uint32_t amb_start);
string getProteinSequence(FILE* fpsq, uint32_t start, uint32_t end);
#endif
