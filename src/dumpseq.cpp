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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#ifdef WIN32
#include "dumpmisc.h"
#else
#include <stdint.h>
#include <endian.h>
#endif

#include "dumpseq.h"

#define MAX_INT  0xffffffff;
#define LINE_WIDTH 60

const char str_2bit[] = "ACGT";
const char str_4bit[] = "-ACMGRSVTWYHKDBN";

uint32_t LUT[256];

static void error(char *format, ...)
{
  va_list   argptr;

  va_start(argptr, format);
  vfprintf(stderr, format, argptr);
  va_end(argptr);

  exit(255);
}

void init_LUT() {
	int c, i;
	int cc;
	for (cc = 0; cc < 256; ++cc) {
	  c = cc;
	  uint32_t res;
	  for (i = 0; i < 4; ++i) {
		  res <<= 8;
		  res |= str_2bit[c & 3];
		  c >>= 2;
	  }
	  LUT[cc] = res;
	}
}

int printProteinSequence(char* file, uint32_t start, uint32_t end, char* res, int maxLen)
{
  int size;

  const char str[] = "-ABCDEFGHIKLMNPQRSTVWXYZU*OJ";


  size = end - start - 1;
  if (size > maxLen) {
	  return size;
  }

  for (int i = 0; i < size; ++i) res[i] = str[file[start + i]];
  return size;
}

int printDnaSequence (char* file, uint32_t start, uint32_t end, uint32_t amb_start, char* res, int maxLen)
{
  int res_cnt;
  int remainder;
  int amb_res;
  int size;

  uint64_t x;
  uint64_t soff, eoff;

  uint32_t amb_cnt = 0;
  uint32_t large_amb = 0;

  unsigned char *amb_mem = NULL;
  unsigned char *amb_ptr = NULL;

  unsigned char c;

  // if there is an ambiguity table, read it in
  if (amb_start != end) {
	amb_cnt = ((uint32_t*)(file+amb_start))[0];
    amb_cnt = htobe32(amb_cnt);
/*
     if the most significant bit is set on the count, then each
     * correction will take two entries in the table.  the layout
     * is described below.*/

    large_amb = amb_cnt >> 31;
    amb_cnt = amb_cnt & 0x7fffffff;

    // allocate enough space for the ambiguity table
    amb_mem = (unsigned char *) malloc(amb_cnt * sizeof(uint32_t));
    if (amb_mem == NULL) {
      error("Malloc amb table error size %d\n", amb_cnt * sizeof(uint32_t));
    }

    // read the table

    amb_ptr = (unsigned char*)(file+amb_start + sizeof(uint32_t));//amb_mem;

    size = amb_start - start;

  } else {
    size = end - start;
  }

  /* read the last byte of the sequence, so we can calculate the
   * number of residues in the last byte of the sequence (0-3).*/

  c = file[start+size-1];

/*   the least two significant bits indicate how many residues
   * are encoded in the last byte.*/

  remainder = c & 0x03;
  res_cnt = (size - 1) * 4 + remainder;
  if (res_cnt > maxLen) {
	  return res_cnt;
  }

  int i;
  for (x = 0; x < size-1; ++x) {
	  ((uint32_t*)res)[x] = LUT[(unsigned char)file[start+x]];
  }
  int ptr = (size-1)*4;
  c = file[start+size-1];
  for (i = 0; i < remainder; ++i) {
	  res[ptr++] = str_2bit[c >> 6];
	  c <<= 2;
  }


  while (amb_cnt > 0) {
	// get the residue symbol
	amb_res = (int) (*amb_ptr >> 4);
	 /*the layout of the ambiguity table differs if it is using
	 * large offsets, i.e. offsets > 16 million.
	 *
	 * for small offsets the layout is:
	 *    4 bits - necleotide
	 *    4 bits - repeat count
	 *   24 bits - offset
	 *
	 * for large offsets the layout is:
	 *    4 bits - necleotide
	 *   12 bits - repeat count
	 *   48 bits - offset*/

	if (large_amb) {
	  // get the repeat count
	  eoff  = (((uint64_t) (*amb_ptr & 0x0f)) << 8);
	  eoff += (((uint64_t) *(amb_ptr+1)) << 0);

	  // get the offset
	  soff  = (((uint64_t) *(amb_ptr+2)) << 40);
	  soff += (((uint64_t) *(amb_ptr+3)) << 32);
	  soff += (((uint64_t) *(amb_ptr+4)) << 24);
	  soff += (((uint64_t) *(amb_ptr+5)) << 16);
	  soff += (((uint64_t) *(amb_ptr+6)) << 8);
	  soff += (((uint64_t) *(amb_ptr+7)) << 0);

	  amb_ptr += 8;
	  amb_cnt -= 2;
	} else {
	  // get the repeat count
	  eoff  = (uint64_t) (*amb_ptr & 0x0f);

	  // get the offset
	  soff  = (((uint64_t) *(amb_ptr+1)) << 16);
	  soff += (((uint64_t) *(amb_ptr+2)) << 8);
	  soff += (((uint64_t) *(amb_ptr+3)) << 0);

	  amb_ptr += 4;
	  amb_cnt -= 1;
	}
	eoff += soff;
	for (x = soff; x <= eoff; ++x) res[x] = str_4bit[amb_res];
  }
  return res_cnt;
}

string getDnaSequence(FILE* fpsq, uint32_t start, uint32_t end,
		uint32_t amb_start) {
  int res_cnt;
  int remainder;
  int amb_res;
  uint64_t size;

  uint64_t x;
  uint64_t soff, eoff;

  uint32_t amb_cnt = 0;
  uint32_t large_amb = 0;

  unsigned char *amb_mem = NULL;
  unsigned char *amb_ptr = NULL;

  unsigned char c;

  if (fseek(fpsq, start, SEEK_SET) != 0) {
    error("Seek seq staart 0x%08x io error %d\n", start, ferror(fpsq));
  }
  char* file = (char*)malloc(sizeof(char) * (end - start));
  fread(file, sizeof(char), end - start, fpsq);

  // if there is an ambiguity table, read it in
  if (amb_start != end) {
	amb_cnt = ((uint32_t*)(file+(amb_start-start)))[0];
	amb_cnt = htobe32(amb_cnt);
/*
	 if the most significant bit is set on the count, then each
	 * correction will take two entries in the table.  the layout
	 * is described below.*/

	large_amb = amb_cnt >> 31;
	amb_cnt = amb_cnt & 0x7fffffff;

	// allocate enough space for the ambiguity table
	amb_mem = (unsigned char *) malloc(amb_cnt * sizeof(uint32_t));
	if (amb_mem == NULL) {
	  error("Malloc amb table error size %d\n", amb_cnt * sizeof(uint32_t));
	}

	// read the table

	amb_ptr = (unsigned char*)(file+(amb_start - start + sizeof(uint32_t)));//amb_mem;

	size = amb_start - start;

  } else {
	size = end - start;
  }

  /* read the last byte of the sequence, so we can calculate the
   * number of residues in the last byte of the sequence (0-3).*/

  c = file[size-1];

/*   the least two significant bits indicate how many residues
   * are encoded in the last byte.*/

  remainder = c & 0x03;
  res_cnt = (size - 1) * 4 + remainder;
  char* res = (char*)malloc(sizeof(char) * (res_cnt + 1));
  res[res_cnt] = 0;


  int i;
  for (x = 0; x < size-1; ++x) {
	  ((uint32_t*)res)[x] = LUT[(unsigned char)file[x]];
  }
  int ptr = (size-1)*4;
  c = file[size-1];
  for (i = 0; i < remainder; ++i) {
	  res[ptr++] = str_2bit[c >> 6];
	  c <<= 2;
  }


  while (amb_cnt > 0) {
	// get the residue symbol
	amb_res = (int) (*amb_ptr >> 4);
	 /*the layout of the ambiguity table differs if it is using
	 * large offsets, i.e. offsets > 16 million.
	 *
	 * for small offsets the layout is:
	 *    4 bits - necleotide
	 *    4 bits - repeat count
	 *   24 bits - offset
	 *
	 * for large offsets the layout is:
	 *    4 bits - necleotide
	 *   12 bits - repeat count
	 *   48 bits - offset*/

	if (large_amb) {
	  // get the repeat count
	  eoff  = (((uint64_t) (*amb_ptr & 0x0f)) << 8);
	  eoff += (((uint64_t) *(amb_ptr+1)) << 0);

	  // get the offset
	  soff  = (((uint64_t) *(amb_ptr+2)) << 40);
	  soff += (((uint64_t) *(amb_ptr+3)) << 32);
	  soff += (((uint64_t) *(amb_ptr+4)) << 24);
	  soff += (((uint64_t) *(amb_ptr+5)) << 16);
	  soff += (((uint64_t) *(amb_ptr+6)) << 8);
	  soff += (((uint64_t) *(amb_ptr+7)) << 0);

	  amb_ptr += 8;
	  amb_cnt -= 2;
	} else {
	  // get the repeat count
	  eoff  = (uint64_t) (*amb_ptr & 0x0f);

	  // get the offset
	  soff  = (((uint64_t) *(amb_ptr+1)) << 16);
	  soff += (((uint64_t) *(amb_ptr+2)) << 8);
	  soff += (((uint64_t) *(amb_ptr+3)) << 0);

	  amb_ptr += 4;
	  amb_cnt -= 1;
	}
	eoff += soff;
	for (x = soff; x <= eoff; ++x) res[x] = str_4bit[amb_res];
  }
  string ans(res);
  free(res);
  return ans;
}

string getProteinSequence(FILE* fpsq, uint32_t start, uint32_t end) {
  int size;

  const char str[] = "-ABCDEFGHIKLMNPQRSTVWXYZU*OJ";

  if (fseek(fpsq, start, SEEK_SET) != 0) {
    error("Seek seq staart 0x%08x io error %d\n", start, ferror(fpsq));
  }
  char* file = (char*)malloc(sizeof(char) * (end - start));
  fread(file, sizeof(char), end - start, fpsq);

  size = end - start - 1;

  string res = "";
  for (int i = 0; i < size; ++i) res += str[file[i]];
  return res;
}
