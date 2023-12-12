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
#include <ctype.h>
#ifdef WIN32
#include "dumpmisc.h"
#else
#include <stdint.h>
#endif

#include "dumphdr.h"


string StaticString::visibleString = "";

typedef struct _ASNINFO {
  int            indent;
  int            remaining;
  unsigned char *ptr;
  int            size;
  uint32_t       offset;
  unsigned char *buf;
  char          *name;
} ASNINFO;

typedef void (*asnParserFn)(ASNINFO *asn);

typedef struct _ASN_TABLE
{
  char        *name;
  asnParserFn  parser;
  int          optional;
} ASN_TABLE;

static void error(ASNINFO *asn, char *format, ...)
{
  int       i;
  char     *ptr;
  char      str[32];
  uint32_t  offset = 0;

  va_list   argptr;

  if (asn->name != NULL) {
    fprintf(stderr, "Error when parsing: %s\n", asn->name);
  }

  va_start(argptr, format);
  vfprintf(stderr, format, argptr);
  va_end(argptr);

  /* dump the header information */
  if (asn->size > 0 && asn->buf != NULL) {
    fprintf(stderr, "Header offset: 0x%08x\n", asn->offset);
    if (asn->size > asn->remaining) offset = asn->size - asn->remaining;
    fprintf(stderr, "Error offset: 0x%08x\n", offset);
    for (i = 0; i < asn->size; ++i) {
      if ((i % 16) == 0) {
        if (i != 0) fprintf(stderr, " |%s|\n", str);
        fprintf(stderr, "%08x ", i);
        ptr = str;
      }
      if ((i % 8)  == 0) fprintf(stderr, " ");
      fprintf(stderr, "%02x ", asn->buf[i]);
      *ptr = (isprint(asn->buf[i])) ? asn->buf[i] : '.';
      *(ptr+1) = '\0';
      ++ptr;
    }

    while ((i % 16) != 0) {
      if ((i % 8)  == 0) fprintf(stderr, " ");
      fprintf(stderr, "   ");
      ++i;
    };

    fprintf(stderr, " |%s|\n", str);
  }

  exit(255);
}


static void indent(ASNINFO *asn)
{
  int c = asn->indent;
}

static void block(ASNINFO *asn)
{
  asn->indent += 2;
}

static void unblock(ASNINFO *asn)
{
  asn->indent -= 2;
  indent(asn);
}

static unsigned char get(ASNINFO *asn)
{
  unsigned char t;

  if (asn->remaining <= 0) {
    error(asn, "End of buffer encountered in get()\n");
  }

  t = *asn->ptr++;
  --asn->remaining;

  return t;
}

static int accept(ASNINFO *asn, const unsigned char *str, int len)
{
  int i;
  unsigned char *t;

  /* if there are not enough buffer left, dont even try */
  if (asn->remaining < len) return 0;

  t = asn->ptr;
  for (i = 0; i < len; ++i) {
    if (*t++ != *str++) break;
  }

  /* check for a match */
  if (i == len) {
    asn->ptr += len;
    asn->remaining -= len;
  }

  return (i == len);
}

static int expect(ASNINFO *asn, const unsigned char *str, int len)
{
  int i;
  int j;
  int xl;
  uint32_t e1 = 0;
  uint32_t e2 = 0;

  for (i = 0; i < len; ++i) {
    if (i >= asn->remaining) {

      /* lets try to print the expected value and what was actually
       * found.  when the code was written, the longest len could
       * be was two, so this should work.
       */
      xl = (len > sizeof(e1)) ? sizeof(e1) : len;
      xl = xl * 2;
      for (j = 0; j < len; ++j) {
        e1 = (e1 << 8) + (uint32_t)str[i];
      }

      error(asn, "End of buffer at %d of 0x%0*x\n", i, xl, e1);
    }
    if (asn->ptr[i] != str[i]) {

      /* lets try to print the expected value and what was actually
       * found.  when the code was written, the longest len could
       * be was two, so this should work.
       */
      xl = (len > sizeof(e1)) ? sizeof(e1) : len;
      xl = xl * 2;
      for (j = 0; j < len; ++j) {
        e1 = (e1 << 8) + (uint32_t)asn->ptr[i];
        e2 = (e2 << 8) + (uint32_t)str[i];
      }

      error(asn, "Expected 0x%0*x found 0x%0*x after %d\n", xl, e2, xl, e1, i);
    }
  }

  asn->ptr += len;
  asn->remaining -= len;

  return 1;
}



static void parseVisibleStringOut(ASNINFO *asn)
{
  int i, s;
  int len;

  unsigned char c;

  /* make sure we are dealing with a visible string which start
   * with a 0x1a followed by the length of the string, then the
   * string.  NOTE: the strings are NOT zero terminated.
   */
  expect(asn, (unsigned char *) "\x1a", 1);

  /* parse the length.  if the most significant bit is not set, i.e. the
   * length is less than 128, then the length of the string is encoded in
   * the single byte.  otherwise, the 7 least significant bits encode the
   * number of bytes used to encode the length of the string.
   */
  c = get(asn);
  if (c < 128) {
    len = c;
  } else {
    /* figure out how many bytes the length takes up */
    s = c - 128;
    if (s > sizeof(len)) {
      error(asn, "String length %d is too large\n", s);
    }

    /* read in the length one byte at a time, most significant bytes
     * are first. */
    len = 0;
    for (i = 0; i < s; ++i) {
      c = get(asn);
      len = (len << 8) + (unsigned int)c;
    }
  }

  /* print the string */
  while (len > 0) {
	char curChar = get(asn);
	StaticString::visibleString += curChar;

    --len;
  }
}


static void parseVisibleString(ASNINFO *asn)
{
  int i, s;
  int len;

  unsigned char c;

  /* make sure we are dealing with a visible string which start
   * with a 0x1a followed by the length of the string, then the
   * string.  NOTE: the strings are NOT zero terminated.
   */
  expect(asn, (unsigned char *) "\x1a", 1);

  /* parse the length.  if the most significant bit is not set, i.e. the
   * length is less than 128, then the length of the string is encoded in
   * the single byte.  otherwise, the 7 least significant bits encode the
   * number of bytes used to encode the length of the string.
   */
  c = get(asn);
  if (c < 128) {
    len = c;
  } else {
    /* figure out how many bytes the length takes up */
    s = c - 128;
    if (s > sizeof(len)) {
      error(asn, "String length %d is too large\n", s);
    }

    /* read in the length one byte at a time, most significant bytes
     * are first. */
    len = 0;
    for (i = 0; i < s; ++i) {
      c = get(asn);
      len = (len << 8) + (unsigned int)c;
    }
  }

  /* print the string */
  while (len > 0) {
	char curChar = get(asn);
    --len;
  }
}

static void parseInteger(ASNINFO *asn)
{
  int i;
  int len;
  int value;

  unsigned char c;

  /* make sure we are dealing with an integer.  integers start
   * with a 0x02 followed by the integers length.
   */
  expect(asn, (unsigned char *) "\x02", 1);

  len = get(asn);
  if (len > sizeof(value)) {
    error(asn, "Integer length %d is too large\n", len);
  }

  /* most significat bytes first */
  value = 0;
  for (i = 0; i < len; ++i) {
    c = get(asn);
    value = (value << 8) + (unsigned int)c;
  }
}

static void parseChoice(ASNINFO *asn, const ASN_TABLE *table)
{
  int i;
  unsigned char hdr[2];

  char *name = asn->name;

  /* which choice is indicated by the first byte, 0xa0 for the
   * first item, 0xa1 for the second, 0xa2 for the third, etc.
   */
  hdr[0] = '\xa0';
  hdr[1] = '\x80';

  /* parse the entries of the choice.  go through the list of
   * choices until we find a match or we hit the end of the
   * list.
   */
  for (i = 0; table[i].parser != NULL; ++i) {
    hdr[0] = 0xa0 + i;
    if (accept(asn, hdr, 2)) {
      asn->name = table[i].name;
      table[i].parser(asn);
      asn->name = name;

      expect(asn, (unsigned char *) "\x00\x00", 2);
      break;
    }
  }

  /* make sure we found a choice */
  if (table[i].parser == NULL) {
    error(asn, "No choice was found\n");
  }
}

static void parseSequence(ASNINFO *asn, const ASN_TABLE *table)
{
  int i;
  unsigned char hdr[2];

  char *name = asn->name;

  /* which item of the sequence is indicated by the first byte,
   * 0xa0 for the first item, 0xa1 for the second, 0xa2 for the
   * third, etc.
   */
  hdr[0] = '\xa0';
  hdr[1] = '\x80';

  /* beginning of the sequence starts with a 0x3080 */
  expect(asn, (unsigned char *) "\x30\x80", 2);

  block(asn);

  /* parse the entries of the sequence going from item to item.
   * if the item is missing and it is NOT optional, through an
   * error reporting the missing item.
   */
  for (i = 0; table[i].parser != NULL; ++i) {
    hdr[0] = 0xa0 + i;
    if (accept(asn, hdr, 2)) {

      indent(asn);

      asn->name = table[i].name;
      table[i].parser(asn);
      asn->name = name;

      expect(asn, (unsigned char *) "\x00\x00", 2);
    } else if (!table[i].optional) {
      error(asn, "Expecting %s\n", table[i].name);
    }
  }

  unblock(asn);

  /* end if the sequence ends with 0x0000 */
  expect(asn, (unsigned char *) "\x00\x00", 2);
}

static void parseSequenceOf(ASNINFO *asn, const asnParserFn parser)
{
  /* beginning of the sequence */
  expect(asn, (unsigned char *) "\x30\x80", 2);

  block(asn);

  /* parse the entries of the sequence of.  the sequence of structure is
   * terminated by two NULL bytes.
   */
  while (!accept(asn, (unsigned char *) "\x00\x00", 2)) {

    indent(asn);

    /* while not at the end of the list, call the designated parser */
    parser(asn);
  }

  unblock(asn);
}

static void parseDateStd(ASNINFO *asn)
{
  const ASN_TABLE table[] = {
    { "year",               parseInteger,         0 },
    { "month",              parseInteger,         1 },
    { "day",                parseInteger,         1 },
    { "season",             parseVisibleString,   1 },
    { "hour",               parseInteger,         1 },
    { "minute",             parseInteger,         1 },
    { "second",             parseInteger,         1 },
    { 0,                    0,                    0 },
  };

  parseSequence(asn, table);
}

static void parseDate(ASNINFO *asn)
{
  const ASN_TABLE table[] = {
    { "str",                parseVisibleString,   0 },
    { "std",                parseDateStd,         0 },
    { 0,                    0,                    0 },
  };

  parseChoice(asn, table);
}

static void parseObjectId(ASNINFO *asn)
{
  const ASN_TABLE table[] = {
    { "id",                 parseInteger,         0 },
    { "str",                parseVisibleStringOut,   0 },
    { 0,                    0,                    0 },
  };

  parseChoice(asn, table);
}

static void parseGiimportId(ASNINFO *asn)
{
  const ASN_TABLE table[] = {
    { "id",                 parseInteger,         0 },
    { "db",                 parseVisibleString,   1 },
    { "release",            parseVisibleString,   1 },
    { 0,                    0,                    0 },
  };

  parseSequence(asn, table);
}

static void parseTextseqId(ASNINFO *asn)
{
  const ASN_TABLE table[] = {
    { "name",               parseVisibleString,   1 },
    { "accession",          parseVisibleString,   1 },
    { "release",            parseVisibleString,   1 },
    { "version",            parseInteger,         1 },
    { 0,                    0,                    0 },
  };

  parseSequence(asn, table);
}

static void parseIdPatNumber(ASNINFO *asn)
{
  const ASN_TABLE table[] = {
    { "number",             parseVisibleString,   0 },
    { "app-number",         parseVisibleString,   0 },
    { 0,                    0,                    0 },
  };

  parseChoice(asn, table);
}

static void parseIdPat(ASNINFO *asn)
{
  const ASN_TABLE table[] = {
    { "country",            parseVisibleString,   0 },
    { "id",                 parseIdPatNumber,     0 },
    { "doc-type",           parseVisibleString,   1 },
    { 0,                    0,                    0 },
  };

  parseSequence(asn, table);
}

static void parsePatentSeqId(ASNINFO *asn)
{
  const ASN_TABLE table[] = {
    { "seqid",              parseInteger,         0 },
    { "cit",                parseIdPat,           0 },
    { 0,                    0,                    0 },
  };

  parseSequence(asn, table);
}

static void parseDbtag(ASNINFO *asn)
{
  const ASN_TABLE table[] = {
    { "db",                 parseVisibleString,   0 },
    { "tag",                parseObjectId,        0 },
    { 0,                    0,                    0 },
  };

  parseSequence(asn, table);
}

static void parsePdbSeqId(ASNINFO *asn)
{
  const ASN_TABLE table[] = {
    { "mol",                parseVisibleString,   0 },
    { "chain",              parseInteger,         1 },
    { "rel",                parseDate,            1 },
    { 0,                    0,                    0 },
  };

  parseSequence(asn, table);
}

static void parseSeqId(ASNINFO *asn)
{
  const ASN_TABLE table[] = {
    { "local",              parseObjectId,        0 },
    { "gibbsq",             parseInteger,         0 },
    { "gibbmt",             parseInteger,         0 },
    { "giim",               parseGiimportId,      0 },
    { "genbank",            parseTextseqId,       0 },
    { "embl",               parseTextseqId,       0 },
    { "pir",                parseTextseqId,       0 },
    { "swissprot",          parseTextseqId,       0 },
    { "patent",             parsePatentSeqId,     0 },
    { "other",              parseTextseqId,       0 },
    { "general",            parseDbtag,           0 },
    { "gi",                 parseInteger,         0 },
    { "ddbj",               parseTextseqId,       0 },
    { "prf",                parseTextseqId,       0 },
    { "pdb",                parsePdbSeqId,        0 },
    { "tpg",                parseTextseqId,       0 },
    { "tpe",                parseTextseqId,       0 },
    { "tpd",                parseTextseqId,       0 },
    { "gpipe",              parseTextseqId,       0 },
    { "named-annot",        parseTextseqId,       0 },
    { 0,                    0,                    0 },
  };

  parseChoice(asn, table);
}

static void parseSeqOfSeqIds(ASNINFO *asn)
{
  parseSequenceOf(asn, parseSeqId);
}

static void parseMembership(ASNINFO *asn)
{
  parseSequenceOf(asn, parseInteger);
}

static void parseLinks(ASNINFO *asn)
{
  parseSequenceOf(asn, parseInteger);
}

static void parseOtherInfo(ASNINFO *asn)
{
  parseSequenceOf(asn, parseInteger);
}

static void parseBlastDefLine(ASNINFO *asn)
{
  const ASN_TABLE table[] = {
    { "title",              parseVisibleStringOut,   1 },
    { "seqid",              parseSeqOfSeqIds,     1 },
    { "taxid",              parseInteger,         1 },
    { "memberships",        parseMembership,      1 },
    { "links",              parseLinks,           1 },
    { "other-info",         parseOtherInfo,       1 },
    { 0,                    0,                    0 },
  };

  asn->name = "Blast-def-line";
  parseSequence(asn, table);
}

void parseBlastHeader(FILE *fphr, uint32_t start, uint32_t end)
{
  uint32_t size;
  ASNINFO  asn;

  unsigned char *ptr;

  ptr      = NULL;
  asn.buf  = NULL;
  asn.name = NULL;
  asn.size = 0;

  /* figure out how big the header is and read it in */
  size = end - start;
  ptr = (unsigned char*)malloc(size);
  if (ptr == NULL) {
    error(&asn, "Malloc error of size %d\n", size);
  }

  /* read in the header */
  if (fseek(fphr, start, SEEK_SET) != 0) {
    error(&asn, "seek io error: %d\n", ferror(fphr));
  }

  if (fread(ptr, sizeof(char), size, fphr) != size) {
    error(&asn, "Read error: %d\n", ferror(fphr));
  }

  /* set up the asn structure */
  asn.ptr       = ptr;
  asn.buf       = ptr;
  asn.indent    = 0;
  asn.offset    = start;
  asn.size      = size;
  asn.remaining = size;
  asn.name      = "Blast-def-line-set";

  parseSequenceOf(&asn, parseBlastDefLine);

  /* the number of bytes in the buffer should be zero */
  if (asn.remaining != 0) {
    error(&asn, "Unexpected characters remaining (%d)\n", asn.remaining);
  }
  if (ptr != NULL) free(ptr);
}
