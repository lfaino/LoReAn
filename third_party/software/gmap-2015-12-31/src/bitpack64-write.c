static char rcsid[] = "$Id: bitpack64-write.c 179363 2015-11-20 22:13:52Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "bitpack64-write.h"

#ifdef WORDS_BIGENDIAN
#include "bigendian.h"		/* For FWRITE_UINTS */
#else
#include "littleendian.h"	/* For FWRITE_UINTS */
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For memset */
#include "mem.h"
#include "assert.h"
#include "fopen.h"
#include "popcount.h"

#ifdef HAVE_SSE2
#include <emmintrin.h>
#endif


/* #define ALLOW_ODD_PACKSIZES 1 */

/* #define USE_ONE_FILE_FOR_FIXED 1 */

#define DIFFERENTIAL_METAINFO_SIZE 2
#define PAIRED_METAINFO_SIZE 3
#define RANK_METAINFO_SIZE 1	/* A variant of differential, where packsize is always 6 (lg 64) */
#define DIRECT_METAINFO_SIZE 1
#define BLOCKSIZE 64
#define POSITIONS_PAGE 4294967296 /* 2^32 */

#define BUFFER_SIZE 1000000


/* Note: For offset pointers, where we need fast cumulative sums, we
   use vertical format (where successive values are in different
   packed unsigned ints).  For lcp, we want raw values, and vertical
   format is still slightly more efficient than horizontal format. */

#ifdef HAVE_SSE2
static int
write_reg_buffered_vert (FILE *strm_fp, Positionsptr_T *strm_buffer,
			 int strm_buffer_size, int strm_buffer_i, __m128i OutReg) {

#if 0
  /* Type casting method (when we passed in pointer to OutReg).  Needs a memory fence. */
  UINT4 *buffer = (UINT4 *) OutReg;
  _mm_lfence();  /* Needed to avoid storing incorrect values into strm_buffer */
#else
  /* Storing method.  Safer.  */
  UINT4 buffer[4];
  _mm_store_si128((__m128i *) buffer,OutReg);
#endif

  /* printf("Writing %08X %08X %08X %08X\n",buffer[0],buffer[1],buffer[2],buffer[3]); */

  strm_buffer[strm_buffer_i++] = buffer[0];
  if (strm_buffer_i == strm_buffer_size) {
    FWRITE_UINTS(strm_buffer,strm_buffer_size,strm_fp);
    strm_buffer_i = 0;
  }

  strm_buffer[strm_buffer_i++] = buffer[1];
  if (strm_buffer_i == strm_buffer_size) {
    FWRITE_UINTS(strm_buffer,strm_buffer_size,strm_fp);
    strm_buffer_i = 0;
  }

  strm_buffer[strm_buffer_i++] = buffer[2];
  if (strm_buffer_i == strm_buffer_size) {
    FWRITE_UINTS(strm_buffer,strm_buffer_size,strm_fp);
    strm_buffer_i = 0;
  }

  strm_buffer[strm_buffer_i++] = buffer[3];
  if (strm_buffer_i == strm_buffer_size) {
    FWRITE_UINTS(strm_buffer,strm_buffer_size,strm_fp);
    strm_buffer_i = 0;
  }

  return strm_buffer_i;
}
#else
static int
write_reg_buffered_vert (FILE *strm_fp, Positionsptr_T *strm_buffer,
			 int strm_buffer_size, int strm_buffer_i,
			 UINT4 *horizontal, int nwritten) {
  UINT4 vertical[64];
  int nrows = nwritten/4, row, column, k;

  /* Convert to vertical */
  for (column = 0; column < 4; column++) {
    k = column;
    for (row = 0; row < nrows; row++) {
      vertical[k] = *horizontal++;
      k += 4;
    }
  }
    
  /* Send to output buffer */
  for (k = 0; k < nwritten; k++) {
    /* printf("Writing %08X\n",vertical[k]); */
    strm_buffer[strm_buffer_i++] = vertical[k];
    if (strm_buffer_i == strm_buffer_size) {
      FWRITE_UINTS(strm_buffer,strm_buffer_size,strm_fp);
      strm_buffer_i = 0;
    }
  }

  return strm_buffer_i;
}
#endif



static int
write_reg_buffered_horiz (FILE *strm_fp, Positionsptr_T *strm_buffer,
			  int strm_buffer_size, int strm_buffer_i,
			  UINT4 *values, int nwritten) {
  int k;

  /* Send to output buffer */
  for (k = 0; k < nwritten; k++) {
    /* printf("Writing %08X\n",values[k]); */
    strm_buffer[strm_buffer_i++] = values[k];
    if (strm_buffer_i == strm_buffer_size) {
      FWRITE_UINTS(strm_buffer,strm_buffer_size,strm_fp);
      strm_buffer_i = 0;
    }
  }

  return strm_buffer_i;
}




#ifdef HAVE_SSE2
static __m128i mask1, mask2, mask3, mask4, mask5, mask6, mask7, mask8,
  mask9, mask10, mask11, mask12, mask13, mask14, mask15, mask16,
  mask17, mask18, mask19, mask20, mask21, mask22, mask23, mask24,
  mask25, mask26, mask27, mask28, mask29, mask30, mask31;
#endif


static void
write_setup () {

#ifdef HAVE_SSE2
  mask1 = _mm_set1_epi32(1U);
  mask2 = _mm_set1_epi32(3U);
  mask3 =  _mm_set1_epi32(7U);
  mask4 =  _mm_set1_epi32(15U);
  mask5 =  _mm_set1_epi32(31U);
  mask6 =  _mm_set1_epi32(63U);
  mask7 =  _mm_set1_epi32(127U);
  mask8 =  _mm_set1_epi32(255U);
  mask9 =  _mm_set1_epi32(511U);
  mask10 =  _mm_set1_epi32(1023U);
  mask11 =  _mm_set1_epi32(2047U);
  mask12 =  _mm_set1_epi32(4095U);
  mask13 =  _mm_set1_epi32(8191U);
  mask14 =  _mm_set1_epi32(16383U);
  mask15 =  _mm_set1_epi32(32767U);
  mask16 =  _mm_set1_epi32(65535U);
  mask17 =  _mm_set1_epi32(131071U);
  mask18 =  _mm_set1_epi32(262143U);
  mask19 =  _mm_set1_epi32(524287U);
  mask20 =  _mm_set1_epi32(1048575U);
  mask21 =  _mm_set1_epi32(2097151U);
  mask22 =  _mm_set1_epi32(4194303U);
  mask23 =  _mm_set1_epi32(8388607U);
  mask24 =  _mm_set1_epi32(16777215U);
  mask25 =  _mm_set1_epi32(33554431U);
  mask26 =  _mm_set1_epi32(67108863U);
  mask27 =  _mm_set1_epi32(134217727U);
  mask28 =  _mm_set1_epi32(268435455U);
  mask29 =  _mm_set1_epi32(536870911U);
  mask30 =  _mm_set1_epi32(1073741823U);
  mask31 =  _mm_set1_epi32(2147483647U);
#endif

  return;
}

#ifdef ALLOW_ODD_PACKSIZES
/* nwritten = 1 * 4 = 4 unsigned ints */
static int
write_01_vert (FILE *strm_fp, Positionsptr_T *strm_buffer, int strm_buffer_size, int strm_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask1);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask1);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 1));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask1);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 2));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask1);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 3));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask1);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask1);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 5));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask1);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 6));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask1);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 7));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask1);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask1);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 9));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask1);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 10));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask1);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 11));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask1);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask1);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 13));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask1);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 14));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask1);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 15));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    return strm_buffer_i;
}
#endif


#ifdef HAVE_SSE2
/* nwritten = 1 * 4 = 4 unsigned ints */
static int
write_02_vert (FILE *strm_fp, Positionsptr_T *strm_buffer, int strm_buffer_size, int strm_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask2);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask2);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 2));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask2);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask2);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 6));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask2);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask2);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 10));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask2);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask2);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 14));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask2);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask2);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 18));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask2);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask2);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 22));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask2);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask2);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 26));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask2);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask2);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 30));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    return strm_buffer_i;
}
#endif

static int
pack_02_horiz (UINT4 *out, const UINT4 *in) {
  int column;

  for (column = 0; column < 4; column++) {
    *out |= (*in)   % (1U << 2 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 2 )  ) <<  2 ;
    ++in;
    *out |= ( (*in)   % (1U << 2 )  ) <<  4 ;
    ++in;
    *out |= ( (*in)   % (1U << 2 )  ) <<  6 ;
    ++in;
    *out |= ( (*in)   % (1U << 2 )  ) <<  8 ;
    ++in;
    *out |= ( (*in)   % (1U << 2 )  ) <<  10 ;
    ++in;
    *out |= ( (*in)   % (1U << 2 )  ) <<  12 ;
    ++in;
    *out |= ( (*in)   % (1U << 2 )  ) <<  14 ;
    ++in;
    *out |= ( (*in)   % (1U << 2 )  ) <<  16 ;
    ++in;
    *out |= ( (*in)   % (1U << 2 )  ) <<  18 ;
    ++in;
    *out |= ( (*in)   % (1U << 2 )  ) <<  20 ;
    ++in;
    *out |= ( (*in)   % (1U << 2 )  ) <<  22 ;
    ++in;
    *out |= ( (*in)   % (1U << 2 )  ) <<  24 ;
    ++in;
    *out |= ( (*in)   % (1U << 2 )  ) <<  26 ;
    ++in;
    *out |= ( (*in)   % (1U << 2 )  ) <<  28 ;
    ++in;
    *out |= ( (*in)   % (1U << 2 )  ) <<  30 ;
    ++out;
    ++in;
  }

  return 4;
}



#ifdef ALLOW_ODD_PACKSIZES
/* nwritten = 2 * 4 = 8 unsigned ints */
static int
write_03_vert (FILE *strm_fp, Positionsptr_T *strm_buffer, int strm_buffer_size, int strm_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask3);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask3);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 3));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask3);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 6));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask3);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 9));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask3);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask3);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 15));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask3);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 18));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask3);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 21));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask3);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask3);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 27));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask3);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 30));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 3 - 1);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask3);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 1));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask3);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask3);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 7));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask3);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 10));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask3);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 13));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    return strm_buffer_i;
}
#endif



#ifdef ALLOW_ODD_PACKSIZES
/* nwritten = 3 * 4 = 12 unsigned ints */
static int
write_05_vert (FILE *strm_fp, Positionsptr_T *strm_buffer, int strm_buffer_size, int strm_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask5);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask5);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 5));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask5);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 10));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask5);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 15));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask5);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask5);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 25));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask5);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 30));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 5 - 3);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask5);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 3));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask5);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask5);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 13));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask5);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 18));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask5);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 23));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask5);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 5 - 1);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask5);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 1));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask5);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 6));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask5);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 11));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    return strm_buffer_i;
}
#endif


#ifdef HAVE_SSE2
/* nwritten = 3 * 4 = 12 unsigned ints */
static int
write_06_vert (FILE *strm_fp, Positionsptr_T *strm_buffer, int strm_buffer_size, int strm_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask6);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask6);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 6));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask6);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask6);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 18));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask6);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask6);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 30));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 6 - 4);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask6);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask6);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 10));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask6);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask6);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 22));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask6);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 6 - 2);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask6);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 2));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask6);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask6);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 14));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask6);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask6);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 26));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    return strm_buffer_i;
}
#endif

static int
pack_06_horiz (UINT4 *out, const UINT4 *in) {
  int column;

  for (column = 0; column < 4; column++) {
    *out |= (*in)   % (1U << 6 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 6 )  ) <<  6 ;
    ++in;
    *out |= ( (*in)   % (1U << 6 )  ) <<  12 ;
    ++in;
    *out |= ( (*in)   % (1U << 6 )  ) <<  18 ;
    ++in;
    *out |= ( (*in)   % (1U << 6 )  ) <<  24 ;
    ++in;
    *out |= ( (*in)   % (1U << 6 )  ) <<  30 ;
    ++out;
    *out |=  ( (*in)   % (1U << 6 ) ) >> ( 6  -  4 );
    ++in;
    *out |= ( (*in)   % (1U << 6 )  ) <<  4 ;
    ++in;
    *out |= ( (*in)   % (1U << 6 )  ) <<  10 ;
    ++in;
    *out |= ( (*in)   % (1U << 6 )  ) <<  16 ;
    ++in;
    *out |= ( (*in)   % (1U << 6 )  ) <<  22 ;
    ++in;
    *out |= ( (*in)   % (1U << 6 )  ) <<  28 ;
    ++out;
    *out |=  ( (*in)   % (1U << 6 ) ) >> ( 6  -  2 );
    ++in;
    *out |= ( (*in)   % (1U << 6 )  ) <<  2 ;
    ++in;
    *out |= ( (*in)   % (1U << 6 )  ) <<  8 ;
    ++in;
    *out |= ( (*in)   % (1U << 6 )  ) <<  14 ;
    ++in;
    *out |= ( (*in)   % (1U << 6 )  ) <<  20 ;
    ++in;
    *out |= ( (*in)   % (1U << 6 )  ) <<  26 ;
    ++out;
    ++in;
  }

  return 12;
}



#ifdef ALLOW_ODD_PACKSIZES
/* nwritten = 4 * 4 = 16 unsigned ints */
static int
write_07_vert (FILE *strm_fp, Positionsptr_T *strm_buffer, int strm_buffer_size, int strm_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask7);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask7);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 7));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask7);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 14));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask7);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 21));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask7);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 7 - 3);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask7);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 3));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask7);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 10));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask7);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 17));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask7);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask7);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 31));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 7 - 6);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask7);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 6));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask7);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 13));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask7);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask7);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 27));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 7 - 2);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask7);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 2));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask7);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 9));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    return strm_buffer_i;
}
#endif


#ifdef ALLOW_ODD_PACKSIZES
/* nwritten = 5 * 4 = 20 unsigned ints */
static int
write_09_vert (FILE *strm_fp, Positionsptr_T *strm_buffer, int strm_buffer_size, int strm_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask9);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask9);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 9));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask9);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 18));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask9);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 27));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 9 - 4);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask9);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask9);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 13));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask9);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 22));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask9);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 31));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 9 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask9);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask9);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 17));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask9);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 26));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 9 - 3);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask9);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 3));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask9);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask9);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 21));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask9);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 30));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 9 - 7);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask9);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 7));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    return strm_buffer_i;
}
#endif


#ifdef HAVE_SSE2
/* nwritten = 5 * 4 = 20 unsigned ints */
static int
write_10_vert (FILE *strm_fp, Positionsptr_T *strm_buffer, int strm_buffer_size, int strm_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask10);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask10);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 10));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask10);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask10);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 30));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 10 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask10);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask10);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 18));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask10);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 10 - 6);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask10);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 6));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask10);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask10);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 26));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 10 - 4);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask10);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask10);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 14));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask10);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 10 - 2);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask10);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 2));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask10);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask10);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 22));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    return strm_buffer_i;
}
#endif

static int
pack_10_horiz (UINT4 *out, const UINT4 *in) {
  int column;

  for (column = 0; column < 4; column++) {
    *out |= (*in)   % (1U << 10 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 10 )  ) <<  10 ;
    ++in;
    *out |= ( (*in)   % (1U << 10 )  ) <<  20 ;
    ++in;
    *out |= ( (*in)   % (1U << 10 )  ) <<  30 ;
    ++out;
    *out |=  ( (*in)   % (1U << 10 ) ) >> ( 10  -  8 );
    ++in;
    *out |= ( (*in)   % (1U << 10 )  ) <<  8 ;
    ++in;
    *out |= ( (*in)   % (1U << 10 )  ) <<  18 ;
    ++in;
    *out |= ( (*in)   % (1U << 10 )  ) <<  28 ;
    ++out;
    *out |=  ( (*in)   % (1U << 10 ) ) >> ( 10  -  6 );
    ++in;
    *out |= ( (*in)   % (1U << 10 )  ) <<  6 ;
    ++in;
    *out |= ( (*in)   % (1U << 10 )  ) <<  16 ;
    ++in;
    *out |= ( (*in)   % (1U << 10 )  ) <<  26 ;
    ++out;
    *out |=  ( (*in)   % (1U << 10 ) ) >> ( 10  -  4 );
    ++in;
    *out |= ( (*in)   % (1U << 10 )  ) <<  4 ;
    ++in;
    *out |= ( (*in)   % (1U << 10 )  ) <<  14 ;
    ++in;
    *out |= ( (*in)   % (1U << 10 )  ) <<  24 ;
    ++out;
    *out |=  ( (*in)   % (1U << 10 ) ) >> ( 10  -  2 );
    ++in;
    *out |= ( (*in)   % (1U << 10 )  ) <<  2 ;
    ++in;
    *out |= ( (*in)   % (1U << 10 )  ) <<  12 ;
    ++in;
    *out |= ( (*in)   % (1U << 10 )  ) <<  22 ;
    ++out;
    ++in;
  }

  return 20;
}




#ifdef ALLOW_ODD_PACKSIZES
/* nwritten = 6 * 4 = 24 unsigned ints */
static int
write_11_vert (FILE *strm_fp, Positionsptr_T *strm_buffer, int strm_buffer_size, int strm_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask11);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask11);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 11));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask11);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 22));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 11 - 1);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask11);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 1));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask11);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask11);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 23));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 11 - 2);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask11);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 2));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask11);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 13));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask11);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 11 - 3);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask11);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 3));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask11);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 14));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask11);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 25));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 11 - 4);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask11);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask11);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 15));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask11);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 26));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 11 - 5);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask11);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 5));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    return strm_buffer_i;
}
#endif


#ifdef HAVE_SSE2
/* nwritten = 6 * 4 = 24 unsigned ints */
static int
write_12_vert (FILE *strm_fp, Positionsptr_T *strm_buffer, int strm_buffer_size, int strm_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask12);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask12);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask12);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 12 - 4);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask12);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask12);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask12);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 12 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask12);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask12);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    InReg = _mm_and_si128(_mm_load_si128(++in), mask12);

    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask12);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask12);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 12 - 4);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask12);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask12);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask12);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 12 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask12);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask12);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    return strm_buffer_i;
}
#endif

static int
pack_12_horiz (UINT4 *out, const UINT4 *in) {
  int column;

  for (column = 0; column < 4; column++) {

    *out |= (*in)   % (1U << 12 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 12 )  ) <<  12 ;
    ++in;
    *out |= ( (*in)   % (1U << 12 )  ) <<  24 ;
    ++out;
    *out |=  ( (*in)   % (1U << 12 ) ) >> ( 12  -  4 );
    ++in;
    *out |= ( (*in)   % (1U << 12 )  ) <<  4 ;
    ++in;
    *out |= ( (*in)   % (1U << 12 )  ) <<  16 ;
    ++in;
    *out |= ( (*in)   % (1U << 12 )  ) <<  28 ;
    ++out;
    *out |=  ( (*in)   % (1U << 12 ) ) >> ( 12  -  8 );
    ++in;
    *out |= ( (*in)   % (1U << 12 )  ) <<  8 ;
    ++in;
    *out |= ( (*in)   % (1U << 12 )  ) <<  20 ;
    ++out;
    ++in;
    *out |= (*in)   % (1U << 12 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 12 )  ) <<  12 ;
    ++in;
    *out |= ( (*in)   % (1U << 12 )  ) <<  24 ;
    ++out;
    *out |=  ( (*in)   % (1U << 12 ) ) >> ( 12  -  4 );
    ++in;
    *out |= ( (*in)   % (1U << 12 )  ) <<  4 ;
    ++in;
    *out |= ( (*in)   % (1U << 12 )  ) <<  16 ;
    ++in;
    *out |= ( (*in)   % (1U << 12 )  ) <<  28 ;
    ++out;
    *out |=  ( (*in)   % (1U << 12 ) ) >> ( 12  -  8 );
    ++in;
    *out |= ( (*in)   % (1U << 12 )  ) <<  8 ;
    ++in;
    *out |= ( (*in)   % (1U << 12 )  ) <<  20 ;
    ++out;
    ++in;
  }

  return 24;
}




#ifdef ALLOW_ODD_PACKSIZES
/* nwritten = 7 * 4 = 28 unsigned ints */
static int
write_13_vert (FILE *strm_fp, Positionsptr_T *strm_buffer, int strm_buffer_size, int strm_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask13);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask13);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 13));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask13);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 26));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 13 - 7);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask13);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 7));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask13);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 13 - 1);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask13);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 1));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask13);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 14));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask13);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 27));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 13 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask13);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask13);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 21));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 13 - 2);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask13);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 2));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask13);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 15));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask13);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 13 - 9);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask13);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 9));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask13);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 22));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 13 - 3);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask13);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 3));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    return strm_buffer_i;
}
#endif


#ifdef HAVE_SSE2
/* nwritten = 7 * 4 = 28 unsigned ints */
static int
write_14_vert (FILE *strm_fp, Positionsptr_T *strm_buffer, int strm_buffer_size, int strm_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask14);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask14);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 14));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask14);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 14 - 10);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask14);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 10));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask14);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 14 - 6);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask14);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 6));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask14);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 14 - 2);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask14);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 2));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask14);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask14);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 30));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 14 - 12);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask14);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask14);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 26));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 14 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask14);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask14);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 22));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 14 - 4);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask14);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask14);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 18));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    return strm_buffer_i;
}
#endif


static int
pack_14_horiz (UINT4 *out, const UINT4 *in) {
  int column;

  for (column = 0; column < 4; column++) {
    *out |= (*in)   % (1U << 14 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 14 )  ) <<  14 ;
    ++in;
    *out |= ( (*in)   % (1U << 14 )  ) <<  28 ;
    ++out;
    *out |=  ( (*in)   % (1U << 14 ) ) >> ( 14  -  10 );
    ++in;
    *out |= ( (*in)   % (1U << 14 )  ) <<  10 ;
    ++in;
    *out |= ( (*in)   % (1U << 14 )  ) <<  24 ;
    ++out;
    *out |=  ( (*in)   % (1U << 14 ) ) >> ( 14  -  6 );
    ++in;
    *out |= ( (*in)   % (1U << 14 )  ) <<  6 ;
    ++in;
    *out |= ( (*in)   % (1U << 14 )  ) <<  20 ;
    ++out;
    *out |=  ( (*in)   % (1U << 14 ) ) >> ( 14  -  2 );
    ++in;
    *out |= ( (*in)   % (1U << 14 )  ) <<  2 ;
    ++in;
    *out |= ( (*in)   % (1U << 14 )  ) <<  16 ;
    ++in;
    *out |= ( (*in)   % (1U << 14 )  ) <<  30 ;
    ++out;
    *out |=  ( (*in)   % (1U << 14 ) ) >> ( 14  -  12 );
    ++in;
    *out |= ( (*in)   % (1U << 14 )  ) <<  12 ;
    ++in;
    *out |= ( (*in)   % (1U << 14 )  ) <<  26 ;
    ++out;
    *out |=  ( (*in)   % (1U << 14 ) ) >> ( 14  -  8 );
    ++in;
    *out |= ( (*in)   % (1U << 14 )  ) <<  8 ;
    ++in;
    *out |= ( (*in)   % (1U << 14 )  ) <<  22 ;
    ++out;
    *out |=  ( (*in)   % (1U << 14 ) ) >> ( 14  -  4 );
    ++in;
    *out |= ( (*in)   % (1U << 14 )  ) <<  4 ;
    ++in;
    *out |= ( (*in)   % (1U << 14 )  ) <<  18 ;
    ++out;
    ++in;
  }

  return 28;
}


#ifdef ALLOW_ODD_PACKSIZES
/* nwritten = 8 * 4 = 32 unsigned ints */
static int
write_15_vert (FILE *strm_fp, Positionsptr_T *strm_buffer, int strm_buffer_size, int strm_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask15);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask15);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 15));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask15);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 30));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 15 - 13);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask15);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 13));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask15);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 15 - 11);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask15);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 11));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask15);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 26));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 15 - 9);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask15);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 9));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask15);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 15 - 7);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask15);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 7));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask15);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 22));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 15 - 5);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask15);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 5));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask15);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 15 - 3);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask15);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 3));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask15);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 18));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 15 - 1);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask15);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 1));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    return strm_buffer_i;
}
#endif



#ifdef ALLOW_ODD_PACKSIZES
/* nwritten = 9 * 4 = 36 unsigned ints */
static int
write_17_vert (FILE *strm_fp, Positionsptr_T *strm_buffer, int strm_buffer_size, int strm_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask17);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask17);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 17));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 17 - 2);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask17);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 2));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask17);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 19));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 17 - 4);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask17);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask17);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 21));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 17 - 6);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask17);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 6));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask17);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 23));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 17 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask17);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask17);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 25));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 17 - 10);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask17);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 10));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask17);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 27));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 17 - 12);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask17);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask17);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 29));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 17 - 14);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask17);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 14));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask17);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 31));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);


    OutReg = _mm_srli_epi32(InReg, 17 - 16);
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    return strm_buffer_i;
}
#endif



#ifdef HAVE_SSE2
/* nwritten = 9 * 4 = 36 unsigned ints */
static int
write_18_vert (FILE *strm_fp, Positionsptr_T *strm_buffer, int strm_buffer_size, int strm_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask18);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask18);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 18));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 18 - 4);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask18);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask18);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 22));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 18 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask18);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask18);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 26));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 18 - 12);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask18);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask18);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 30));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 18 - 16);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask18);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 18 - 2);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask18);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 2));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask18);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 18 - 6);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask18);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 6));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask18);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 18 - 10);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask18);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 10));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask18);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 18 - 14);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask18);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 14));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    return strm_buffer_i;
}
#endif

static int
pack_18_horiz (UINT4 *out, const UINT4 *in) {
  int column;

  for (column = 0; column < 4; column++) {
    *out |= (*in)   % (1U << 18 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 18 )  ) <<  18 ;
    ++out;
    *out |=  ( (*in)   % (1U << 18 ) ) >> ( 18  -  4 );
    ++in;
    *out |= ( (*in)   % (1U << 18 )  ) <<  4 ;
    ++in;
    *out |= ( (*in)   % (1U << 18 )  ) <<  22 ;
    ++out;
    *out |=  ( (*in)   % (1U << 18 ) ) >> ( 18  -  8 );
    ++in;
    *out |= ( (*in)   % (1U << 18 )  ) <<  8 ;
    ++in;
    *out |= ( (*in)   % (1U << 18 )  ) <<  26 ;
    ++out;
    *out |=  ( (*in)   % (1U << 18 ) ) >> ( 18  -  12 );
    ++in;
    *out |= ( (*in)   % (1U << 18 )  ) <<  12 ;
    ++in;
    *out |= ( (*in)   % (1U << 18 )  ) <<  30 ;
    ++out;
    *out |=  ( (*in)   % (1U << 18 ) ) >> ( 18  -  16 );
    ++in;
    *out |= ( (*in)   % (1U << 18 )  ) <<  16 ;
    ++out;
    *out |=  ( (*in)   % (1U << 18 ) ) >> ( 18  -  2 );
    ++in;
    *out |= ( (*in)   % (1U << 18 )  ) <<  2 ;
    ++in;
    *out |= ( (*in)   % (1U << 18 )  ) <<  20 ;
    ++out;
    *out |=  ( (*in)   % (1U << 18 ) ) >> ( 18  -  6 );
    ++in;
    *out |= ( (*in)   % (1U << 18 )  ) <<  6 ;
    ++in;
    *out |= ( (*in)   % (1U << 18 )  ) <<  24 ;
    ++out;
    *out |=  ( (*in)   % (1U << 18 ) ) >> ( 18  -  10 );
    ++in;
    *out |= ( (*in)   % (1U << 18 )  ) <<  10 ;
    ++in;
    *out |= ( (*in)   % (1U << 18 )  ) <<  28 ;
    ++out;
    *out |=  ( (*in)   % (1U << 18 ) ) >> ( 18  -  14 );
    ++in;
    *out |= ( (*in)   % (1U << 18 )  ) <<  14 ;
    ++out;
    ++in;
  }

  return 36;
}



#ifdef ALLOW_ODD_PACKSIZES
/* nwritten = 10 * 4 = 40 unsigned ints */
static int
write_19_vert (FILE *strm_fp, Positionsptr_T *strm_buffer, int strm_buffer_size, int strm_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask19);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask19);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 19));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 19 - 6);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask19);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 6));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask19);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 25));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 19 - 12);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask19);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask19);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 31));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 19 - 18);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask19);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 18));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 19 - 5);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask19);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 5));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask19);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 19 - 11);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask19);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 11));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask19);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 30));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 19 - 17);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask19);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 17));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 19 - 4);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask19);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask19);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 23));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 19 - 10);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask19);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 10));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask19);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 29));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);


    OutReg = _mm_srli_epi32(InReg, 19 - 16);
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    return strm_buffer_i;
}
#endif


#ifdef HAVE_SSE2
/* nwritten = 10 * 4 = 40 unsigned ints */
static int
write_20_vert (FILE *strm_fp, Positionsptr_T *strm_buffer, int strm_buffer_size, int strm_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask20);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask20);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 20 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask20);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask20);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 20 - 16);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask20);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 20 - 4);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask20);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask20);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 20 - 12);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask20);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    InReg = _mm_and_si128(_mm_load_si128(++in), mask20);

    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask20);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 20 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask20);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask20);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 20 - 16);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask20);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 20 - 4);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask20);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask20);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 20 - 12);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask20);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    return strm_buffer_i;
}
#endif

static int
pack_20_horiz (UINT4 *out, const UINT4 *in) {
  int column;

  for (column = 0; column < 4; column++) {
    *out |= (*in)   % (1U << 20 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 20 )  ) <<  20 ;
    ++out;
    *out |=  ( (*in)   % (1U << 20 ) ) >> ( 20  -  8 );
    ++in;
    *out |= ( (*in)   % (1U << 20 )  ) <<  8 ;
    ++in;
    *out |= ( (*in)   % (1U << 20 )  ) <<  28 ;
    ++out;
    *out |=  ( (*in)   % (1U << 20 ) ) >> ( 20  -  16 );
    ++in;
    *out |= ( (*in)   % (1U << 20 )  ) <<  16 ;
    ++out;
    *out |=  ( (*in)   % (1U << 20 ) ) >> ( 20  -  4 );
    ++in;
    *out |= ( (*in)   % (1U << 20 )  ) <<  4 ;
    ++in;
    *out |= ( (*in)   % (1U << 20 )  ) <<  24 ;
    ++out;
    *out |=  ( (*in)   % (1U << 20 ) ) >> ( 20  -  12 );
    ++in;
    *out |= ( (*in)   % (1U << 20 )  ) <<  12 ;
    ++out;
    ++in;
    *out |= (*in)   % (1U << 20 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 20 )  ) <<  20 ;
    ++out;
    *out |=  ( (*in)   % (1U << 20 ) ) >> ( 20  -  8 );
    ++in;
    *out |= ( (*in)   % (1U << 20 )  ) <<  8 ;
    ++in;
    *out |= ( (*in)   % (1U << 20 )  ) <<  28 ;
    ++out;
    *out |=  ( (*in)   % (1U << 20 ) ) >> ( 20  -  16 );
    ++in;
    *out |= ( (*in)   % (1U << 20 )  ) <<  16 ;
    ++out;
    *out |=  ( (*in)   % (1U << 20 ) ) >> ( 20  -  4 );
    ++in;
    *out |= ( (*in)   % (1U << 20 )  ) <<  4 ;
    ++in;
    *out |= ( (*in)   % (1U << 20 )  ) <<  24 ;
    ++out;
    *out |=  ( (*in)   % (1U << 20 ) ) >> ( 20  -  12 );
    ++in;
    *out |= ( (*in)   % (1U << 20 )  ) <<  12 ;
    ++out;
    ++in;
  }

  return 40;
}



#ifdef ALLOW_ODD_PACKSIZES
/* nwritten = 11 * 4 = 44 unsigned ints */
static int
write_21_vert (FILE *strm_fp, Positionsptr_T *strm_buffer, int strm_buffer_size, int strm_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask21);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask21);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 21));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 21 - 10);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask21);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 10));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask21);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 31));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 21 - 20);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask21);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 21 - 9);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask21);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 9));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask21);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 30));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 21 - 19);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask21);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 19));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 21 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask21);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask21);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 29));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 21 - 18);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask21);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 18));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 21 - 7);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask21);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 7));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask21);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 21 - 17);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask21);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 17));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 21 - 6);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask21);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 6));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask21);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 27));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);


    OutReg = _mm_srli_epi32(InReg, 21 - 16);
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    return strm_buffer_i;
}
#endif


#ifdef HAVE_SSE2
/* nwritten = 11 * 4 = 44 unsigned ints */
static int
write_22_vert (FILE *strm_fp, Positionsptr_T *strm_buffer, int strm_buffer_size, int strm_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask22);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask22);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 22));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 22 - 12);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask22);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 22 - 2);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask22);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 2));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask22);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 22 - 14);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask22);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 14));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 22 - 4);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask22);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask22);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 26));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 22 - 16);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask22);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 22 - 6);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask22);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 6));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask22);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 22 - 18);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask22);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 18));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 22 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask22);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask22);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 30));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 22 - 20);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask22);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 22 - 10);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask22);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 10));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    return strm_buffer_i;
}
#endif

static int
pack_22_horiz (UINT4 *out, const UINT4 *in) {
  int column;

  for (column = 0; column < 4; column++) {
    *out |= (*in)   % (1U << 22 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 22 )  ) <<  22 ;
    ++out;
    *out |=  ( (*in)   % (1U << 22 ) ) >> ( 22  -  12 );
    ++in;
    *out |= ( (*in)   % (1U << 22 )  ) <<  12 ;
    ++out;
    *out |=  ( (*in)   % (1U << 22 ) ) >> ( 22  -  2 );
    ++in;
    *out |= ( (*in)   % (1U << 22 )  ) <<  2 ;
    ++in;
    *out |= ( (*in)   % (1U << 22 )  ) <<  24 ;
    ++out;
    *out |=  ( (*in)   % (1U << 22 ) ) >> ( 22  -  14 );
    ++in;
    *out |= ( (*in)   % (1U << 22 )  ) <<  14 ;
    ++out;
    *out |=  ( (*in)   % (1U << 22 ) ) >> ( 22  -  4 );
    ++in;
    *out |= ( (*in)   % (1U << 22 )  ) <<  4 ;
    ++in;
    *out |= ( (*in)   % (1U << 22 )  ) <<  26 ;
    ++out;
    *out |=  ( (*in)   % (1U << 22 ) ) >> ( 22  -  16 );
    ++in;
    *out |= ( (*in)   % (1U << 22 )  ) <<  16 ;
    ++out;
    *out |=  ( (*in)   % (1U << 22 ) ) >> ( 22  -  6 );
    ++in;
    *out |= ( (*in)   % (1U << 22 )  ) <<  6 ;
    ++in;
    *out |= ( (*in)   % (1U << 22 )  ) <<  28 ;
    ++out;
    *out |=  ( (*in)   % (1U << 22 ) ) >> ( 22  -  18 );
    ++in;
    *out |= ( (*in)   % (1U << 22 )  ) <<  18 ;
    ++out;
    *out |=  ( (*in)   % (1U << 22 ) ) >> ( 22  -  8 );
    ++in;
    *out |= ( (*in)   % (1U << 22 )  ) <<  8 ;
    ++in;
    *out |= ( (*in)   % (1U << 22 )  ) <<  30 ;
    ++out;
    *out |=  ( (*in)   % (1U << 22 ) ) >> ( 22  -  20 );
    ++in;
    *out |= ( (*in)   % (1U << 22 )  ) <<  20 ;
    ++out;
    *out |=  ( (*in)   % (1U << 22 ) ) >> ( 22  -  10 );
    ++in;
    *out |= ( (*in)   % (1U << 22 )  ) <<  10 ;
    ++out;
    ++in;
  }

  return 44;
}


#ifdef ALLOW_ODD_PACKSIZES
/* nwritten = 12 * 4 = 48 unsigned ints */
static int
write_23_vert (FILE *strm_fp, Positionsptr_T *strm_buffer, int strm_buffer_size, int strm_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask23);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask23);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 23));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 23 - 14);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask23);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 14));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 23 - 5);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask23);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 5));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask23);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 23 - 19);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask23);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 19));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 23 - 10);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask23);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 10));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 23 - 1);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask23);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 1));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask23);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 23 - 15);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask23);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 15));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 23 - 6);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask23);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 6));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask23);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 29));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 23 - 20);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask23);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 23 - 11);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask23);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 11));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 23 - 2);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask23);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 2));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask23);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 25));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 23 - 16);
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    return strm_buffer_i;
}
#endif


#ifdef HAVE_SSE2
/* nwritten = 12 * 4 = 48 unsigned ints */
static int
write_24_vert (FILE *strm_fp, Positionsptr_T *strm_buffer, int strm_buffer_size, int strm_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask24);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask24);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 24 - 16);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask24);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 24 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask24);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    InReg = _mm_and_si128(_mm_load_si128(++in), mask24);

    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask24);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 24 - 16);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask24);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 24 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask24);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    InReg = _mm_and_si128(_mm_load_si128(++in), mask24);

    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask24);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 24 - 16);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask24);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 24 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask24);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    InReg = _mm_and_si128(_mm_load_si128(++in), mask24);

    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask24);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 24 - 16);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask24);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 24 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask24);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    return strm_buffer_i;
}
#endif

static int
pack_24_horiz (UINT4 *out, const UINT4 *in) {
  int column;

  for (column = 0; column < 4; column++) {
    *out |= (*in)   % (1U << 24 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 24 )  ) <<  24 ;
    ++out;
    *out |=  ( (*in)   % (1U << 24 ) ) >> ( 24  -  16 );
    ++in;
    *out |= ( (*in)   % (1U << 24 )  ) <<  16 ;
    ++out;
    *out |=  ( (*in)   % (1U << 24 ) ) >> ( 24  -  8 );
    ++in;
    *out |= ( (*in)   % (1U << 24 )  ) <<  8 ;
    ++out;
    ++in;
    *out |= (*in)   % (1U << 24 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 24 )  ) <<  24 ;
    ++out;
    *out |=  ( (*in)   % (1U << 24 ) ) >> ( 24  -  16 );
    ++in;
    *out |= ( (*in)   % (1U << 24 )  ) <<  16 ;
    ++out;
    *out |=  ( (*in)   % (1U << 24 ) ) >> ( 24  -  8 );
    ++in;
    *out |= ( (*in)   % (1U << 24 )  ) <<  8 ;
    ++out;
    ++in;
    *out |= (*in)   % (1U << 24 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 24 )  ) <<  24 ;
    ++out;
    *out |=  ( (*in)   % (1U << 24 ) ) >> ( 24  -  16 );
    ++in;
    *out |= ( (*in)   % (1U << 24 )  ) <<  16 ;
    ++out;
    *out |=  ( (*in)   % (1U << 24 ) ) >> ( 24  -  8 );
    ++in;
    *out |= ( (*in)   % (1U << 24 )  ) <<  8 ;
    ++out;
    ++in;
    *out |= (*in)   % (1U << 24 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 24 )  ) <<  24 ;
    ++out;
    *out |=  ( (*in)   % (1U << 24 ) ) >> ( 24  -  16 );
    ++in;
    *out |= ( (*in)   % (1U << 24 )  ) <<  16 ;
    ++out;
    *out |=  ( (*in)   % (1U << 24 ) ) >> ( 24  -  8 );
    ++in;
    *out |= ( (*in)   % (1U << 24 )  ) <<  8 ;
    ++out;
    ++in;
  }

  return 48;
}


#ifdef ALLOW_ODD_PACKSIZES
/* nwritten = 13 * 4 = 52 unsigned ints */
static int
write_25_vert (FILE *strm_fp, Positionsptr_T *strm_buffer, int strm_buffer_size, int strm_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask25);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask25);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 25));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 25 - 18);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask25);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 18));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 25 - 11);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask25);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 11));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 25 - 4);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask25);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask25);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 29));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 25 - 22);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask25);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 22));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 25 - 15);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask25);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 15));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 25 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask25);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 25 - 1);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask25);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 1));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask25);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 26));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 25 - 19);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask25);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 19));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 25 - 12);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask25);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 25 - 5);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask25);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 5));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask25);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 30));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 25 - 23);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask25);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 23));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);


    OutReg = _mm_srli_epi32(InReg, 25 - 16);
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    return strm_buffer_i;
}
#endif


#ifdef HAVE_SSE2
/* nwritten = 13 * 4 = 52 unsigned ints */
static int
write_26_vert (FILE *strm_fp, Positionsptr_T *strm_buffer, int strm_buffer_size, int strm_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask26);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask26);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 26));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 26 - 20);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask26);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 26 - 14);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask26);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 14));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 26 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask26);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 26 - 2);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask26);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 2));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask26);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 26 - 22);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask26);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 22));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 26 - 16);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask26);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 26 - 10);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask26);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 10));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 26 - 4);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask26);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask26);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 30));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 26 - 24);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask26);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 26 - 18);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask26);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 18));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 26 - 12);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask26);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 26 - 6);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask26);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 6));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    return strm_buffer_i;
}
#endif


static int
pack_26_horiz (UINT4 *out, const UINT4 *in) {
  int column;

  for (column = 0; column < 4; column++) {
    *out |= (*in)   % (1U << 26 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 26 )  ) <<  26 ;
    ++out;
    *out |=  ( (*in)   % (1U << 26 ) ) >> ( 26  -  20 );
    ++in;
    *out |= ( (*in)   % (1U << 26 )  ) <<  20 ;
    ++out;
    *out |=  ( (*in)   % (1U << 26 ) ) >> ( 26  -  14 );
    ++in;
    *out |= ( (*in)   % (1U << 26 )  ) <<  14 ;
    ++out;
    *out |=  ( (*in)   % (1U << 26 ) ) >> ( 26  -  8 );
    ++in;
    *out |= ( (*in)   % (1U << 26 )  ) <<  8 ;
    ++out;
    *out |=  ( (*in)   % (1U << 26 ) ) >> ( 26  -  2 );
    ++in;
    *out |= ( (*in)   % (1U << 26 )  ) <<  2 ;
    ++in;
    *out |= ( (*in)   % (1U << 26 )  ) <<  28 ;
    ++out;
    *out |=  ( (*in)   % (1U << 26 ) ) >> ( 26  -  22 );
    ++in;
    *out |= ( (*in)   % (1U << 26 )  ) <<  22 ;
    ++out;
    *out |=  ( (*in)   % (1U << 26 ) ) >> ( 26  -  16 );
    ++in;
    *out |= ( (*in)   % (1U << 26 )  ) <<  16 ;
    ++out;
    *out |=  ( (*in)   % (1U << 26 ) ) >> ( 26  -  10 );
    ++in;
    *out |= ( (*in)   % (1U << 26 )  ) <<  10 ;
    ++out;
    *out |=  ( (*in)   % (1U << 26 ) ) >> ( 26  -  4 );
    ++in;
    *out |= ( (*in)   % (1U << 26 )  ) <<  4 ;
    ++in;
    *out |= ( (*in)   % (1U << 26 )  ) <<  30 ;
    ++out;
    *out |=  ( (*in)   % (1U << 26 ) ) >> ( 26  -  24 );
    ++in;
    *out |= ( (*in)   % (1U << 26 )  ) <<  24 ;
    ++out;
    *out |=  ( (*in)   % (1U << 26 ) ) >> ( 26  -  18 );
    ++in;
    *out |= ( (*in)   % (1U << 26 )  ) <<  18 ;
    ++out;
    *out |=  ( (*in)   % (1U << 26 ) ) >> ( 26  -  12 );
    ++in;
    *out |= ( (*in)   % (1U << 26 )  ) <<  12 ;
    ++out;
    *out |=  ( (*in)   % (1U << 26 ) ) >> ( 26  -  6 );
    ++in;
    *out |= ( (*in)   % (1U << 26 )  ) <<  6 ;
    ++out;
    ++in;
  }

  return 52;
}



#ifdef ALLOW_ODD_PACKSIZES
/* nwritten = 14 * 4 = 56 unsigned ints */
static int
write_27_vert (FILE *strm_fp, Positionsptr_T *strm_buffer, int strm_buffer_size, int strm_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask27);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask27);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 27));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 27 - 22);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask27);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 22));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 27 - 17);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask27);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 17));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 27 - 12);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask27);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 27 - 7);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask27);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 7));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 27 - 2);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask27);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 2));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask27);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 29));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 27 - 24);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask27);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 27 - 19);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask27);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 19));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 27 - 14);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask27);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 14));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 27 - 9);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask27);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 9));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 27 - 4);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask27);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask27);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 31));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 27 - 26);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask27);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 26));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 27 - 21);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask27);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 21));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 27 - 16);
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    return strm_buffer_i;
}
#endif


#ifdef HAVE_SSE2
/* nwritten = 14 * 4 = 56 unsigned ints */
static int
write_28_vert (FILE *strm_fp, Positionsptr_T *strm_buffer, int strm_buffer_size, int strm_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask28);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask28);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 28 - 24);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask28);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 28 - 20);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask28);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 28 - 16);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask28);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 28 - 12);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask28);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 28 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask28);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 28 - 4);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask28);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    InReg = _mm_and_si128(_mm_load_si128(++in), mask28);

    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask28);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 28 - 24);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask28);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 28 - 20);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask28);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 28 - 16);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask28);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 28 - 12);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask28);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 28 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask28);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 28 - 4);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask28);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    return strm_buffer_i;
}
#endif


static int
pack_28_horiz (UINT4 *out, const UINT4 *in) {
  int column;

  for (column = 0; column < 4; column++) {
    *out |= (*in)   % (1U << 28 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 28 )  ) <<  28 ;
    ++out;
    *out |=  ( (*in)   % (1U << 28 ) ) >> ( 28  -  24 );
    ++in;
    *out |= ( (*in)   % (1U << 28 )  ) <<  24 ;
    ++out;
    *out |=  ( (*in)   % (1U << 28 ) ) >> ( 28  -  20 );
    ++in;
    *out |= ( (*in)   % (1U << 28 )  ) <<  20 ;
    ++out;
    *out |=  ( (*in)   % (1U << 28 ) ) >> ( 28  -  16 );
    ++in;
    *out |= ( (*in)   % (1U << 28 )  ) <<  16 ;
    ++out;
    *out |=  ( (*in)   % (1U << 28 ) ) >> ( 28  -  12 );
    ++in;
    *out |= ( (*in)   % (1U << 28 )  ) <<  12 ;
    ++out;
    *out |=  ( (*in)   % (1U << 28 ) ) >> ( 28  -  8 );
    ++in;
    *out |= ( (*in)   % (1U << 28 )  ) <<  8 ;
    ++out;
    *out |=  ( (*in)   % (1U << 28 ) ) >> ( 28  -  4 );
    ++in;
    *out |= ( (*in)   % (1U << 28 )  ) <<  4 ;
    ++out;
    ++in;
    *out |= (*in)   % (1U << 28 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 28 )  ) <<  28 ;
    ++out;
    *out |=  ( (*in)   % (1U << 28 ) ) >> ( 28  -  24 );
    ++in;
    *out |= ( (*in)   % (1U << 28 )  ) <<  24 ;
    ++out;
    *out |=  ( (*in)   % (1U << 28 ) ) >> ( 28  -  20 );
    ++in;
    *out |= ( (*in)   % (1U << 28 )  ) <<  20 ;
    ++out;
    *out |=  ( (*in)   % (1U << 28 ) ) >> ( 28  -  16 );
    ++in;
    *out |= ( (*in)   % (1U << 28 )  ) <<  16 ;
    ++out;
    *out |=  ( (*in)   % (1U << 28 ) ) >> ( 28  -  12 );
    ++in;
    *out |= ( (*in)   % (1U << 28 )  ) <<  12 ;
    ++out;
    *out |=  ( (*in)   % (1U << 28 ) ) >> ( 28  -  8 );
    ++in;
    *out |= ( (*in)   % (1U << 28 )  ) <<  8 ;
    ++out;
    *out |=  ( (*in)   % (1U << 28 ) ) >> ( 28  -  4 );
    ++in;
    *out |= ( (*in)   % (1U << 28 )  ) <<  4 ;
    ++out;
    ++in;
  }

  return 56;
}



#ifdef ALLOW_ODD_PACKSIZES
/* nwritten = 15 * 4 = 60 unsigned ints */
static int
write_29_vert (FILE *strm_fp, Positionsptr_T *strm_buffer, int strm_buffer_size, int strm_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask29);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask29);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 29));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 29 - 26);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask29);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 26));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 29 - 23);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask29);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 23));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 29 - 20);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask29);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 29 - 17);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask29);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 17));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 29 - 14);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask29);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 14));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 29 - 11);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask29);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 11));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 29 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask29);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 29 - 5);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask29);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 5));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 29 - 2);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask29);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 2));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask29);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 31));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 29 - 28);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask29);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 29 - 25);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask29);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 25));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 29 - 22);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask29);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 22));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 29 - 19);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask29);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 19));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);


    OutReg = _mm_srli_epi32(InReg, 29 - 16);
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    return strm_buffer_i;
}
#endif


#ifdef HAVE_SSE2
/* nwritten = 15 * 4 = 60 unsigned ints */
static int
write_30_vert (FILE *strm_fp, Positionsptr_T *strm_buffer, int strm_buffer_size, int strm_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask30);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask30);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 30));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 30 - 28);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask30);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 30 - 26);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask30);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 26));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 30 - 24);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask30);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 30 - 22);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask30);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 22));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 30 - 20);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask30);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 30 - 18);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask30);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 18));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 30 - 16);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask30);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 30 - 14);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask30);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 14));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 30 - 12);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask30);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 30 - 10);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask30);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 10));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 30 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask30);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 30 - 6);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask30);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 6));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 30 - 4);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask30);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 30 - 2);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask30);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 2));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    return strm_buffer_i;
}
#endif


static int
pack_30_horiz (UINT4 *out, const UINT4 *in) {
  int column;

  for (column = 0; column < 4; column++) {
    *out |= (*in)   % (1U << 30 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 30 )  ) <<  30 ;
    ++out;
    *out |=  ( (*in)   % (1U << 30 ) ) >> ( 30  -  28 );
    ++in;
    *out |= ( (*in)   % (1U << 30 )  ) <<  28 ;
    ++out;
    *out |=  ( (*in)   % (1U << 30 ) ) >> ( 30  -  26 );
    ++in;
    *out |= ( (*in)   % (1U << 30 )  ) <<  26 ;
    ++out;
    *out |=  ( (*in)   % (1U << 30 ) ) >> ( 30  -  24 );
    ++in;
    *out |= ( (*in)   % (1U << 30 )  ) <<  24 ;
    ++out;
    *out |=  ( (*in)   % (1U << 30 ) ) >> ( 30  -  22 );
    ++in;
    *out |= ( (*in)   % (1U << 30 )  ) <<  22 ;
    ++out;
    *out |=  ( (*in)   % (1U << 30 ) ) >> ( 30  -  20 );
    ++in;
    *out |= ( (*in)   % (1U << 30 )  ) <<  20 ;
    ++out;
    *out |=  ( (*in)   % (1U << 30 ) ) >> ( 30  -  18 );
    ++in;
    *out |= ( (*in)   % (1U << 30 )  ) <<  18 ;
    ++out;
    *out |=  ( (*in)   % (1U << 30 ) ) >> ( 30  -  16 );
    ++in;
    *out |= ( (*in)   % (1U << 30 )  ) <<  16 ;
    ++out;
    *out |=  ( (*in)   % (1U << 30 ) ) >> ( 30  -  14 );
    ++in;
    *out |= ( (*in)   % (1U << 30 )  ) <<  14 ;
    ++out;
    *out |=  ( (*in)   % (1U << 30 ) ) >> ( 30  -  12 );
    ++in;
    *out |= ( (*in)   % (1U << 30 )  ) <<  12 ;
    ++out;
    *out |=  ( (*in)   % (1U << 30 ) ) >> ( 30  -  10 );
    ++in;
    *out |= ( (*in)   % (1U << 30 )  ) <<  10 ;
    ++out;
    *out |=  ( (*in)   % (1U << 30 ) ) >> ( 30  -  8 );
    ++in;
    *out |= ( (*in)   % (1U << 30 )  ) <<  8 ;
    ++out;
    *out |=  ( (*in)   % (1U << 30 ) ) >> ( 30  -  6 );
    ++in;
    *out |= ( (*in)   % (1U << 30 )  ) <<  6 ;
    ++out;
    *out |=  ( (*in)   % (1U << 30 ) ) >> ( 30  -  4 );
    ++in;
    *out |= ( (*in)   % (1U << 30 )  ) <<  4 ;
    ++out;
    *out |=  ( (*in)   % (1U << 30 ) ) >> ( 30  -  2 );
    ++in;
    *out |= ( (*in)   % (1U << 30 )  ) <<  2 ;
    ++out;
    ++in;
  }

  return 60;
}



#ifdef ALLOW_ODD_PACKSIZES
/* nwritten = 16 * 4 = 64 unsigned ints */
static int
write_31_vert (FILE *strm_fp, Positionsptr_T *strm_buffer, int strm_buffer_size, int strm_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask31);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask31);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 31));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 31 - 30);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask31);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 30));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 31 - 29);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask31);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 29));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 31 - 28);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask31);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 31 - 27);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask31);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 27));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 31 - 26);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask31);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 26));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 31 - 25);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask31);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 25));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 31 - 24);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask31);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 31 - 23);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask31);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 23));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 31 - 22);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask31);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 22));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 31 - 21);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask31);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 21));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 31 - 20);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask31);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 31 - 19);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask31);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 19));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 31 - 18);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask31);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 18));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 31 - 17);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask31);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 17));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);


    OutReg = _mm_srli_epi32(InReg, 31 - 16);
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    return strm_buffer_i;
}
#endif


#ifdef HAVE_SSE2
/* nwritten = 16 * 4 = 64 unsigned ints */
static int
write_32_vert (FILE *strm_fp, Positionsptr_T *strm_buffer, int strm_buffer_size, int strm_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_load_si128(in);
    OutReg = InReg;
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    InReg = _mm_load_si128(++in);

    OutReg = InReg;
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    InReg = _mm_load_si128(++in);

    OutReg = InReg;
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    InReg = _mm_load_si128(++in);

    OutReg = InReg;
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    InReg = _mm_load_si128(++in);

    OutReg = InReg;
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    InReg = _mm_load_si128(++in);

    OutReg = InReg;
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    InReg = _mm_load_si128(++in);

    OutReg = InReg;
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    InReg = _mm_load_si128(++in);

    OutReg = InReg;
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    InReg = _mm_load_si128(++in);

    OutReg = InReg;
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    InReg = _mm_load_si128(++in);

    OutReg = InReg;
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    InReg = _mm_load_si128(++in);

    OutReg = InReg;
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    InReg = _mm_load_si128(++in);

    OutReg = InReg;
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    InReg = _mm_load_si128(++in);

    OutReg = InReg;
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    InReg = _mm_load_si128(++in);

    OutReg = InReg;
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    InReg = _mm_load_si128(++in);

    OutReg = InReg;
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    InReg = _mm_load_si128(++in);

    OutReg = InReg;
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    return strm_buffer_i;
}
#endif


static int
pack_32_horiz (UINT4 *out, const UINT4 *in) {
  int column;

  for (column = 0; column < 4; column++) {
    *out = *in;
    ++out;
    ++in;
    *out = *in;
    ++out;
    ++in;
    *out = *in;
    ++out;
    ++in;
    *out = *in;
    ++out;
    ++in;
    *out = *in;
    ++out;
    ++in;
    *out = *in;
    ++out;
    ++in;
    *out = *in;
    ++out;
    ++in;
    *out = *in;
    ++out;
    ++in;
    *out = *in;
    ++out;
    ++in;
    *out = *in;
    ++out;
    ++in;
    *out = *in;
    ++out;
    ++in;
    *out = *in;
    ++out;
    ++in;
    *out = *in;
    ++out;
    ++in;
    *out = *in;
    ++out;
    ++in;
    *out = *in;
    ++out;
    ++in;
    *out = *in;
    ++out;
    ++in;
  }

  return 64;
}



#ifdef HAVE_SSE2
/* nwritten = 2 * 4 = 8 unsigned ints */
static int
write_04_vert (FILE *strm_fp, Positionsptr_T *strm_buffer, int strm_buffer_size, int strm_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask4);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask4);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask4);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask4);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask4);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask4);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask4);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask4);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    InReg = _mm_and_si128(_mm_load_si128(++in), mask4);

    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask4);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask4);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask4);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask4);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask4);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask4);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask4);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    return strm_buffer_i;
}
#endif


static int
pack_04_horiz (UINT4 *out, const UINT4 *in) {
  int column;

  for (column = 0; column < 4; column++) {
    *out |= (*in)   % (1U << 4 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 4 )  ) <<  4 ;
    ++in;
    *out |= ( (*in)   % (1U << 4 )  ) <<  8 ;
    ++in;
    *out |= ( (*in)   % (1U << 4 )  ) <<  12 ;
    ++in;
    *out |= ( (*in)   % (1U << 4 )  ) <<  16 ;
    ++in;
    *out |= ( (*in)   % (1U << 4 )  ) <<  20 ;
    ++in;
    *out |= ( (*in)   % (1U << 4 )  ) <<  24 ;
    ++in;
    *out |= ( (*in)   % (1U << 4 )  ) <<  28 ;
    ++out;
    ++in;
    *out |= (*in)   % (1U << 4 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 4 )  ) <<  4 ;
    ++in;
    *out |= ( (*in)   % (1U << 4 )  ) <<  8 ;
    ++in;
    *out |= ( (*in)   % (1U << 4 )  ) <<  12 ;
    ++in;
    *out |= ( (*in)   % (1U << 4 )  ) <<  16 ;
    ++in;
    *out |= ( (*in)   % (1U << 4 )  ) <<  20 ;
    ++in;
    *out |= ( (*in)   % (1U << 4 )  ) <<  24 ;
    ++in;
    *out |= ( (*in)   % (1U << 4 )  ) <<  28 ;
    ++out;
    ++in;
  }

  return 8;
}


#ifdef HAVE_SSE2
/* nwritten = 4 * 4 = 16 unsigned ints */
static int
write_08_vert (FILE *strm_fp, Positionsptr_T *strm_buffer, int strm_buffer_size, int strm_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask8);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask8);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask8);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask8);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    InReg = _mm_and_si128(_mm_load_si128(++in), mask8);

    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask8);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask8);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask8);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    InReg = _mm_and_si128(_mm_load_si128(++in), mask8);

    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask8);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask8);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask8);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    InReg = _mm_and_si128(_mm_load_si128(++in), mask8);

    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask8);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask8);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask8);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    return strm_buffer_i;
}
#endif


static int
pack_08_horiz (UINT4 *out, const UINT4 *in) {
  int column;

  for (column = 0; column < 4; column++) {
    *out |= (*in)   % (1U << 8 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 8 )  ) <<  8 ;
    ++in;
    *out |= ( (*in)   % (1U << 8 )  ) <<  16 ;
    ++in;
    *out |= ( (*in)   % (1U << 8 )  ) <<  24 ;
    ++out;
    ++in;
    *out |= (*in)   % (1U << 8 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 8 )  ) <<  8 ;
    ++in;
    *out |= ( (*in)   % (1U << 8 )  ) <<  16 ;
    ++in;
    *out |= ( (*in)   % (1U << 8 )  ) <<  24 ;
    ++out;
    ++in;
    *out |= (*in)   % (1U << 8 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 8 )  ) <<  8 ;
    ++in;
    *out |= ( (*in)   % (1U << 8 )  ) <<  16 ;
    ++in;
    *out |= ( (*in)   % (1U << 8 )  ) <<  24 ;
    ++out;
    ++in;
    *out |= (*in)   % (1U << 8 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 8 )  ) <<  8 ;
    ++in;
    *out |= ( (*in)   % (1U << 8 )  ) <<  16 ;
    ++in;
    *out |= ( (*in)   % (1U << 8 )  ) <<  24 ;
    ++out;
    ++in;
  }

  return 16;
}


#ifdef HAVE_SSE2
/* nwritten = 8 * 4 = 32 unsigned ints */
static int
write_16_vert (FILE *strm_fp, Positionsptr_T *strm_buffer, int strm_buffer_size, int strm_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask16);

    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask16);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    InReg = _mm_and_si128(_mm_load_si128(++in), mask16);

    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask16);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    InReg = _mm_and_si128(_mm_load_si128(++in), mask16);

    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask16);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    InReg = _mm_and_si128(_mm_load_si128(++in), mask16);

    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask16);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    InReg = _mm_and_si128(_mm_load_si128(++in), mask16);

    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask16);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    InReg = _mm_and_si128(_mm_load_si128(++in), mask16);

    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask16);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    InReg = _mm_and_si128(_mm_load_si128(++in), mask16);

    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask16);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    InReg = _mm_and_si128(_mm_load_si128(++in), mask16);

    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask16);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					       OutReg);

    return strm_buffer_i;
}
#endif


static int
pack_16_horiz (UINT4 *out, const UINT4 *in) {
  int column;

  for (column = 0; column < 4; column++) {
    *out |= (*in)   % (1U << 16 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 16 )  ) <<  16 ;
    ++out;
    ++in;
    *out |= (*in)   % (1U << 16 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 16 )  ) <<  16 ;
    ++out;
    ++in;
    *out |= (*in)   % (1U << 16 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 16 )  ) <<  16 ;
    ++out;
    ++in;
    *out |= (*in)   % (1U << 16 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 16 )  ) <<  16 ;
    ++out;
    ++in;
    *out |= (*in)   % (1U << 16 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 16 )  ) <<  16 ;
    ++out;
    ++in;
    *out |= (*in)   % (1U << 16 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 16 )  ) <<  16 ;
    ++out;
    ++in;
    *out |= (*in)   % (1U << 16 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 16 )  ) <<  16 ;
    ++out;
    ++in;
    *out |= (*in)   % (1U << 16 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 16 )  ) <<  16 ;
    ++out;
    ++in;
  }

  return 32;
}


/* Vertical format requires all values in a block to be decoded */
#ifdef HAVE_SSE2
static int
write_vert (FILE *strm_fp, Positionsptr_T *strm_buffer, int strm_buffer_size, int strm_buffer_i,
	    const UINT4 *_in, int packsize) {

#if 0
  int i;

  printf("Entering with packsize %d\n",packsize);
  for (i = 0; i < BLOCKSIZE; i++) {
    printf("%d ",_in[i]);
  }
  printf("\n");
#endif

  switch (packsize) {
#ifdef ALLOW_ODD_PACKSIZES
  case 0: return strm_buffer_i;
  case 1: return write_01_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 2: return write_02_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 3: return write_03_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 4: return write_04_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 5: return write_05_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 6: return write_06_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 7: return write_07_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 8: return write_08_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 9: return write_09_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 10: return write_10_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 11: return write_11_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 12: return write_12_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 13: return write_13_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 14: return write_14_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 15: return write_15_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 16: return write_16_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 17: return write_17_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 18: return write_18_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 19: return write_19_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 20: return write_20_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 21: return write_21_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 22: return write_22_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 23: return write_23_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 24: return write_24_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 25: return write_25_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 26: return write_26_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 27: return write_27_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 28: return write_28_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 29: return write_29_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 30: return write_30_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 31: return write_31_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 32: return write_32_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
#else
  case 0: return strm_buffer_i;
  case 2: return write_02_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 4: return write_04_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 6: return write_06_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 8: return write_08_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 10: return write_10_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 12: return write_12_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 14: return write_14_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 16: return write_16_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 18: return write_18_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 20: return write_20_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 22: return write_22_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 24: return write_24_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 26: return write_26_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 28: return write_28_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 30: return write_30_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
  case 32: return write_32_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,_in);
#endif
  default: fprintf(stderr,"packsize of %d not allowed\n",packsize); abort();
  }
}

#else

static void
reorder_values_vertically (Positionsptr_T *vertical, const Positionsptr_T *horizontal) {
  int column, row, k = 0;
  Positionsptr_T *out;

  out = &(vertical[0]);
  for (column = 0; column < 4; column++) {
    k = column;
    for (row = 0; row < BLOCKSIZE/4; row++) {
      *out++ = horizontal[k];
      k += 4;
    }
  }

#if 0
  printf("horizontal\n");
  for (k = 0; k < BLOCKSIZE; k++) {
    if (k % 4 == 0) {
      printf("\n");
    }
    printf("%u ",horizontal[k]);
  }
  printf("\n");

  printf("vertical\n");
  for (k = 0; k < BLOCKSIZE; k++) {
    if (k % (BLOCKSIZE/4) == 0) {
      printf("\n");
    }
    printf("%u ",vertical[k]);
  }
  printf("\n");
#endif

  return;
}


/* Non-SIMD code cannot write vertical format easily, so using
   horizontal code and conversions */
static int
write_vert (FILE *strm_fp, Positionsptr_T *strm_buffer, int strm_buffer_size, int strm_buffer_i,
	    const UINT4 *horizontal, int packsize) {
  int nwritten;
  UINT4 buffer[BLOCKSIZE], vertical[BLOCKSIZE];

#if 0
  int i;

  printf("Entering with packsize %d\n",packsize);
  for (i = 0; i < BLOCKSIZE; i++) {
    printf("%d ",_in[i]);
  }
  printf("\n");
#endif

  reorder_values_vertically(vertical,horizontal);
  memset((void *) buffer,0,BLOCKSIZE*sizeof(UINT4));

  switch (packsize) {
  case 0: return strm_buffer_i;
  case 2: nwritten = pack_02_horiz(buffer,&(vertical[0])); break;
  case 4: nwritten = pack_04_horiz(buffer,&(vertical[0])); break;
  case 6: nwritten = pack_06_horiz(buffer,&(vertical[0])); break;
  case 8: nwritten = pack_08_horiz(buffer,&(vertical[0])); break;
  case 10: nwritten = pack_10_horiz(buffer,&(vertical[0])); break;
  case 12: nwritten = pack_12_horiz(buffer,&(vertical[0])); break;
  case 14: nwritten = pack_14_horiz(buffer,&(vertical[0])); break;
  case 16: nwritten = pack_16_horiz(buffer,&(vertical[0])); break;
  case 18: nwritten = pack_18_horiz(buffer,&(vertical[0])); break;
  case 20: nwritten = pack_20_horiz(buffer,&(vertical[0])); break;
  case 22: nwritten = pack_22_horiz(buffer,&(vertical[0])); break;
  case 24: nwritten = pack_24_horiz(buffer,&(vertical[0])); break;
  case 26: nwritten = pack_26_horiz(buffer,&(vertical[0])); break;
  case 28: nwritten = pack_28_horiz(buffer,&(vertical[0])); break;
  case 30: nwritten = pack_30_horiz(buffer,&(vertical[0])); break;
  case 32: nwritten = pack_32_horiz(buffer,&(vertical[0])); break;
  default: fprintf(stderr,"packsize of %d not allowed\n",packsize); abort();
  }

  return write_reg_buffered_vert(strm_fp,strm_buffer,
				 strm_buffer_size,strm_buffer_i,
				 buffer,nwritten);
}
#endif


/* Columnar order allows just the necessary values in a block to be decoded */
static void
columnar_order (UINT4 *columnar, const UINT4 *vertical) {

  columnar[0] = vertical[0];		/* remainder 1 */
  columnar[1] = vertical[4];		/* remainder 5 */
  columnar[2] = vertical[8];		/* remainder 9 */
  columnar[3] = vertical[12];		/* remainder 13 */
  columnar[4] = vertical[16];		/* remainder 17 */
  columnar[5] = vertical[20];		/* remainder 21 */
  columnar[6] = vertical[24];		/* remainder 25 */
  columnar[7] = vertical[28];		/* remainder 29 */

  columnar[8] = vertical[1];		/* remainder 2 */
  columnar[9] = vertical[5];		/* remainder 6 */
  columnar[10] = vertical[9];		/* remainder 10 */
  columnar[11] = vertical[13];		/* remainder 14 */
  columnar[12] = vertical[17];		/* remainder 18 */
  columnar[13] = vertical[21];		/* remainder 22 */
  columnar[14] = vertical[25];		/* remainder 26 */
  columnar[15] = vertical[29];		/* remainder 30 */

  columnar[16] = vertical[2];		/* remainder 3 */
  columnar[17] = vertical[6];		/* remainder 7 */
  columnar[18] = vertical[10];		/* remainder 11 */
  columnar[19] = vertical[14];		/* remainder 15 */
  columnar[20] = vertical[18];		/* remainder 19 */
  columnar[21] = vertical[22];		/* remainder 23 */
  columnar[22] = vertical[26];		/* remainder 27 */
  columnar[23] = vertical[30];		/* remainder 31 */

  columnar[24] = vertical[3];		/* remainder 4 */
  columnar[25] = vertical[7];		/* remainder 8 */
  columnar[26] = vertical[11];		/* remainder 12 */
  columnar[27] = vertical[15];		/* remainder 16 */
  columnar[28] = vertical[19];		/* remainder 20 */
  columnar[29] = vertical[23];		/* remainder 24 */
  columnar[30] = vertical[27];		/* remainder 28 */
  columnar[31] = vertical[31];		/* remainder 32 */

  columnar[32] = vertical[32];		/* remainder 63 */
  columnar[33] = vertical[36];		/* remainder 59 */
  columnar[34] = vertical[40];		/* remainder 55 */
  columnar[35] = vertical[44];		/* remainder 51 */
  columnar[36] = vertical[48];		/* remainder 47 */
  columnar[37] = vertical[52];		/* remainder 43 */
  columnar[38] = vertical[56];		/* remainder 39 */
  columnar[39] = vertical[60];		/* remainder 35 */

  columnar[40] = vertical[33];		/* remainder 62 */
  columnar[41] = vertical[37];		/* remainder 58 */
  columnar[42] = vertical[41];		/* remainder 54 */
  columnar[43] = vertical[45];		/* remainder 50 */
  columnar[44] = vertical[49];		/* remainder 46 */
  columnar[45] = vertical[53];		/* remainder 42 */
  columnar[46] = vertical[57];		/* remainder 38 */
  columnar[47] = vertical[61];		/* remainder 34 */

  columnar[48] = vertical[34];		/* remainder 61 */
  columnar[49] = vertical[38];		/* remainder 57 */
  columnar[50] = vertical[42];		/* remainder 53 */
  columnar[51] = vertical[46];		/* remainder 49 */
  columnar[52] = vertical[50];		/* remainder 45 */
  columnar[53] = vertical[54];		/* remainder 41 */
  columnar[54] = vertical[58];		/* remainder 37 */
  columnar[55] = vertical[62];		/* remainder 33 */

  columnar[56] = vertical[35];		/* remainder 60 */
  columnar[57] = vertical[39];		/* remainder 56 */
  columnar[58] = vertical[43];		/* remainder 52 */
  columnar[59] = vertical[47];		/* remainder 48 */
  columnar[60] = vertical[51];		/* remainder 44 */
  columnar[61] = vertical[55];		/* remainder 40 */
  columnar[62] = vertical[59];		/* remainder 36 */
  columnar[63] = vertical[63];		/* remainder 32 */

  return;
}


#ifdef HAVE_SSE2

static int
write_columnar (FILE *strm_fp, Positionsptr_T *strm_buffer, int strm_buffer_size, int strm_buffer_i,
		const UINT4 *_in, int packsize) {
  UINT4 columnar[BLOCKSIZE];

#if 0
  int i;

  printf("Entering with packsize %d\n",packsize);
  for (i = 0; i < BLOCKSIZE; i++) {
    printf("%d ",_in[i]);
  }
  printf("\n");
#endif

  columnar_order(columnar,_in);

  switch (packsize) {
  case 0: return strm_buffer_i;
  case 2: return write_02_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,columnar);
  case 4: return write_04_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,columnar);
  case 6: return write_06_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,columnar);
  case 8: return write_08_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,columnar);
  case 10: return write_10_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,columnar);
  case 12: return write_12_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,columnar);
  case 14: return write_14_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,columnar);
  case 16: return write_16_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,columnar);
  case 18: return write_18_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,columnar);
  case 20: return write_20_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,columnar);
  case 22: return write_22_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,columnar);
  case 24: return write_24_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,columnar);
  case 26: return write_26_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,columnar);
  case 28: return write_28_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,columnar);
  case 30: return write_30_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,columnar);
  case 32: return write_32_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,columnar);

  default: fprintf(stderr,"packsize of %d not allowed\n",packsize); abort();
  }
}

#else

/* Non-SIMD code cannot write vertical format easily, so using
   horizontal code and conversions */

static int
write_columnar (FILE *strm_fp, Positionsptr_T *strm_buffer, int strm_buffer_size, int strm_buffer_i,
		const UINT4 *horizontal, int packsize) {
  int nwritten;
  UINT4 buffer[BLOCKSIZE], vertical[BLOCKSIZE];
  UINT4 columnar[BLOCKSIZE];

#if 0
  int i;

  printf("Entering with packsize %d\n",packsize);
  for (i = 0; i < BLOCKSIZE; i++) {
    printf("%d ",_in[i]);
  }
  printf("\n");
#endif

  columnar_order(columnar,horizontal);
  reorder_values_vertically(vertical,columnar);
  memset((void *) buffer,0,BLOCKSIZE*sizeof(UINT4));

  switch (packsize) {
  case 0: return strm_buffer_i;
  case 2: nwritten = pack_02_horiz(buffer,&(vertical[0])); break;
  case 4: nwritten = pack_04_horiz(buffer,&(vertical[0])); break;
  case 6: nwritten = pack_06_horiz(buffer,&(vertical[0])); break;
  case 8: nwritten = pack_08_horiz(buffer,&(vertical[0])); break;
  case 10: nwritten = pack_10_horiz(buffer,&(vertical[0])); break;
  case 12: nwritten = pack_12_horiz(buffer,&(vertical[0])); break;
  case 14: nwritten = pack_14_horiz(buffer,&(vertical[0])); break;
  case 16: nwritten = pack_16_horiz(buffer,&(vertical[0])); break;
  case 18: nwritten = pack_18_horiz(buffer,&(vertical[0])); break;
  case 20: nwritten = pack_20_horiz(buffer,&(vertical[0])); break;
  case 22: nwritten = pack_22_horiz(buffer,&(vertical[0])); break;
  case 24: nwritten = pack_24_horiz(buffer,&(vertical[0])); break;
  case 26: nwritten = pack_26_horiz(buffer,&(vertical[0])); break;
  case 28: nwritten = pack_28_horiz(buffer,&(vertical[0])); break;
  case 30: nwritten = pack_30_horiz(buffer,&(vertical[0])); break;
  case 32: nwritten = pack_32_horiz(buffer,&(vertical[0])); break;
  default: fprintf(stderr,"packsize of %d not allowed\n",packsize); abort();
  }

  return write_reg_buffered_vert(strm_fp,strm_buffer,
				 strm_buffer_size,strm_buffer_i,
				 buffer,nwritten);
}

#endif



/* Horizontal format is slightly more complicated for random access of individual values */
int
Bitpack64_write_horiz (FILE *strm_fp, Positionsptr_T *strm_buffer, int strm_buffer_size, int strm_buffer_i,
		       const UINT4 *horizontal, int packsize) {
  int nwritten;
  UINT4 buffer[BLOCKSIZE];

  write_setup();

#if 0
  int i;

  printf("Entering with packsize %d\n",packsize);
  for (i = 0; i < BLOCKSIZE; i++) {
    printf("%d ",_in[i]);
  }
  printf("\n");
#endif

  memset((void *) buffer,0,BLOCKSIZE*sizeof(UINT4));

  switch (packsize) {
  case 0: return strm_buffer_i;
  case 2: nwritten = pack_02_horiz(buffer,&(horizontal[0])); break;
  case 4: nwritten = pack_04_horiz(buffer,&(horizontal[0])); break;
  case 6: nwritten = pack_06_horiz(buffer,&(horizontal[0])); break;
  case 8: nwritten = pack_08_horiz(buffer,&(horizontal[0])); break;
  case 10: nwritten = pack_10_horiz(buffer,&(horizontal[0])); break;
  case 12: nwritten = pack_12_horiz(buffer,&(horizontal[0])); break;
  case 14: nwritten = pack_14_horiz(buffer,&(horizontal[0])); break;
  case 16: nwritten = pack_16_horiz(buffer,&(horizontal[0])); break;
  case 18: nwritten = pack_18_horiz(buffer,&(horizontal[0])); break;
  case 20: nwritten = pack_20_horiz(buffer,&(horizontal[0])); break;
  case 22: nwritten = pack_22_horiz(buffer,&(horizontal[0])); break;
  case 24: nwritten = pack_24_horiz(buffer,&(horizontal[0])); break;
  case 26: nwritten = pack_26_horiz(buffer,&(horizontal[0])); break;
  case 28: nwritten = pack_28_horiz(buffer,&(horizontal[0])); break;
  case 30: nwritten = pack_30_horiz(buffer,&(horizontal[0])); break;
  case 32: nwritten = pack_32_horiz(buffer,&(horizontal[0])); break;
  default: fprintf(stderr,"packsize of %d not allowed\n",packsize); abort();
  }

  return write_reg_buffered_horiz(strm_fp,strm_buffer,
				  strm_buffer_size,strm_buffer_i,
				  buffer,nwritten);
}



/* Processes 64 values at a time.  Returns packsize. */
/* Handles first 32 values from the initial value, and the last 32
   values from the final value.  More efficient since we need to
   process only half as many inputs. */
static int
compute_q4_diffs_bidir (UINT4 *diffs, UINT4 *values) {
  UINT4 packsize;
  int i;
  UINT4 maxdiff = 0;
  int firstbit, msb;

#if 0
  for (i = 0; i < 64; i++) {
    assert(values[i+1] >= values[i]);
  }
#endif

  maxdiff |= (diffs[32] = values[64] - values[63]);
  maxdiff |= (diffs[33] = values[64] - values[62]);
  maxdiff |= (diffs[34] = values[64] - values[61]);
  maxdiff |= (diffs[35] = values[64] - values[60]);
  for (i = 36; i < 64; i++) {
    maxdiff |= (diffs[i] = values[64+32-(i+1-4)] - values[64+32-(i+1)]);
  }
  for (i = 31; i >= 4; i--) {
    maxdiff |= (diffs[i] = values[i+1] - values[i+1-4]);
  }
  maxdiff |= (diffs[3] = values[4] - values[0]);
  maxdiff |= (diffs[2] = values[3] - values[0]);
  maxdiff |= (diffs[1] = values[2] - values[0]);
  maxdiff |= (diffs[0] = values[1] - values[0]);

  if (maxdiff == 0) {
    /* __builtin_clz() behaves oddly on zero */
    return 0;

  } else {
#ifdef HAVE_BUILTIN_CLZ
    firstbit = __builtin_clz(maxdiff);
    packsize = 32 - firstbit;
#elif defined(HAVE_ASM_BSR)
    asm("bsr %1,%0" : "=r"(msb) : "r"(maxdiff));
    packsize = msb + 1;
#else
    firstbit = ((maxdiff >> 16) ? clz_table[maxdiff >> 16] : 16 + clz_table[maxdiff]);
    packsize = 32 - firstbit;
#endif

#ifdef ALLOW_ODD_PACKSIZES
    return packsize;
#else
    return (packsize + 1) & ~1;	/* Converts packsizes to the next multiple of 2 */
#endif
  }
}


static int
compute_q1_diffs (UINT4 *diffs, UINT4 *values) {
  UINT4 packsize;
  int i;
  UINT4 maxdiff = 0;
  int firstbit, msb;

#if 0
  for (i = 0; i < 64; i++) {
    assert(values[i+1] >= values[i]);
  }
#endif

  for (i = 63; i >= 0; i--) {
    maxdiff |= (diffs[i] = values[i+1] - values[i]);
  }

  if (maxdiff == 0) {
    /* __builtin_clz() behaves oddly on zero */
    return 0;

  } else {
#ifdef HAVE_BUILTIN_CLZ
    firstbit = __builtin_clz(maxdiff);
    packsize = 32 - firstbit;
#elif defined(HAVE_ASM_BSR)
    asm("bsr %1,%0" : "=r"(msb) : "r"(maxdiff));
    packsize = msb + 1;
#else
    firstbit = ((maxdiff >> 16) ? clz_table[maxdiff >> 16] : 16 + clz_table[maxdiff]);
    packsize = 32 - firstbit;
#endif

#ifdef ALLOW_ODD_PACKSIZES
    return packsize;
#else
    return (packsize + 1) & ~1;	/* Converts packsizes to the next multiple of 2 */
#endif
  }
}



#ifdef HAVE_64_BIT
static int
compute_q4_diffs_bidir_huge (UINT4 *diffs, UINT8 *values) {
  UINT4 packsize;
  int i;
  UINT4 maxdiff = 0;
  int firstbit, msb;


  maxdiff |= (diffs[32] = (UINT4) (values[64] - values[63]));
  maxdiff |= (diffs[33] = (UINT4) (values[64] - values[62]));
  maxdiff |= (diffs[34] = (UINT4) (values[64] - values[61]));
  maxdiff |= (diffs[35] = (UINT4) (values[64] - values[60]));
  for (i = 36; i < 64; i++) {
    maxdiff |= (diffs[i] = (UINT4) (values[64+32-(i+1-4)] - values[64+32-(i+1)]));
  }
  for (i = 31; i >= 4; i--) {
    maxdiff |= (diffs[i] = (UINT4) (values[i+1] - values[i+1-4]));
  }
  maxdiff |= (diffs[3] = (UINT4) (values[4] - values[0]));
  maxdiff |= (diffs[2] = (UINT4) (values[3] - values[0]));
  maxdiff |= (diffs[1] = (UINT4) (values[2] - values[0]));
  maxdiff |= (diffs[0] = (UINT4) (values[1] - values[0]));

  if (maxdiff == 0) {
    /* __builtin_clz() behaves oddly on zero */
    return 0;

  } else {
#ifdef HAVE_BUILTIN_CLZ
    firstbit = __builtin_clz(maxdiff);
    packsize = 32 - firstbit;
#elif defined(HAVE_ASM_BSR)
    asm("bsr %1,%0" : "=r"(msb) : "r"(maxdiff));
    packsize = msb + 1;
#else
    firstbit = ((maxdiff >> 16) ? clz_table[maxdiff >> 16] : 16 + clz_table[maxdiff]);
    packsize = 32 - firstbit;
#endif

#ifdef ALLOW_ODD_PACKSIZES
    return packsize;
#else
    return (packsize + 1) & ~1;	/* Converts packsizes to the next multiple of 2 */
#endif
  }
}
#endif



/* We want to store values 0..n, with final value at ascending[n]
   possibly stored as the final metainfo value */
/* Stored in columnar order */
void
Bitpack64_write_differential (char *ptrsfile, char *compfile, UINT4 *ascending, Oligospace_T n) {
  FILE *ptrs_fp, *comp_fp;
  UINT4 *ptrs;
  int ptri, i;
  Oligospace_T positioni;

  /* Buffer is used to avoid frequent writes to the file */
  UINT4 *buffer;
  int buffer_size = BUFFER_SIZE;
  int buffer_i;

  UINT4 diffs[BLOCKSIZE], last_block[BLOCKSIZE+1];

  UINT4 nwritten;
  int packsize;


  write_setup();

  printf("Entered Bitpack64_write_differential with n %llu\n",n);

  /* 2 metavalues: nwritten (pointer) and cumulative sum for block.
     Packsize can be computed from difference between successive
     pointers, if only even packsizes are allowed */
  ptrs = (UINT4 *) CALLOC(((n + BLOCKSIZE)/BLOCKSIZE + 1) * DIFFERENTIAL_METAINFO_SIZE,sizeof(UINT4));
  ptri = 0;

  if ((comp_fp = FOPEN_WRITE_BINARY(compfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",compfile);
    exit(9);
  }
  buffer = (UINT4 *) CALLOC(buffer_size,sizeof(UINT4));
  buffer_i = 0;

  nwritten = 0U;

  /* Last value of ascending is at ascending[n] */
  /* Use <= n instead of < n, because we want ascending[n] to be taken care of by unpack_00, not a check for remainder == 0 */
  for (positioni = 0; positioni + BLOCKSIZE <= n; positioni += BLOCKSIZE) {
    /* Pointer */
    ptrs[ptri++] = nwritten/4;	/* In 128-bit registers */

    /* Value for start of block */
    ptrs[ptri++] = ascending[positioni];
	
    /* Pack block of 64 diffs */
    packsize = compute_q4_diffs_bidir(diffs,&(ascending[positioni]));
    buffer_i = write_columnar(comp_fp,buffer,buffer_size,buffer_i,diffs,packsize);

#ifdef ALLOW_ODD_PACKSIZES
    nwritten += 2 * ((packsize + 1) & ~1);
#else
    nwritten += 2 * packsize;
#endif
  }

  /* Old: Check for positioni < n, because if positioni == n, ascending[n] will be taken care of as metainfo */
  /* Use <= n instead of < n, because we want ascending[n] to be taken care of by unpack_00, not a check for remainder == 0 */
  if (positioni <= n) {
    /* Finish last block of 64 */
    ptrs[ptri++] = nwritten/4;	/* In 128-bit registers */

    /* Value for start of block */
    ptrs[ptri++] = ascending[positioni];

    /* For differential, want <=.  For direct, want < */
    for (i = 0; i <= (int) (n - positioni); i++) {
      last_block[i] = ascending[positioni+i];
    }
    for ( ; i <= BLOCKSIZE; i++) {
      /* Copy last value for rest of block */
      last_block[i] = ascending[n];
    }

    /* Pack block of < 64 diffs */
    packsize = compute_q4_diffs_bidir(diffs,last_block);
    buffer_i = write_columnar(comp_fp,buffer,buffer_size,buffer_i,diffs,packsize);

#ifdef ALLOW_ODD_PACKSIZES
    nwritten += 2 * ((packsize + 1) & ~1);
#else
    nwritten += 2 * packsize;
#endif
  }


  /* Write the final pointer, which will point after the end of the file */
  ptrs[ptri++] = nwritten/4;	/* In 128-bit registers */

  /* Value for end of block */
  ptrs[ptri++] = ascending[n];

  if ((ptrs_fp = FOPEN_WRITE_BINARY(ptrsfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",ptrsfile);
    exit(9);
  } else {
    FWRITE_UINTS(ptrs,ptri,ptrs_fp);
    FREE(ptrs);
    fclose(ptrs_fp);
  }
    
  /* Empty buffer */
  if (buffer_i > 0) {
    FWRITE_UINTS(buffer,buffer_i,comp_fp);	
    buffer_i = 0;
  }
  FREE(buffer);
  fclose(comp_fp);

  return;
}


/* We want to store values 0..n, with final value at ascending[n]
   possibly stored as the final metainfo value */
/* D4 stored in columnar order, plus D1 stored as direct */
void
Bitpack64_write_differential_paired (char *ptrsfile, char *compfile, UINT4 *ascending, Oligospace_T n) {
  FILE *ptrs_fp, *comp_fp;
  UINT4 *ptrs;
  int ptri, i;
  Oligospace_T positioni;

  /* Buffer is used to avoid frequent writes to the file */
  UINT4 *buffer;
  int buffer_size = BUFFER_SIZE;
  int buffer_i;

  UINT4 diffs[BLOCKSIZE], last_block[BLOCKSIZE+1];

  UINT4 nwritten;
  int packsize;


  write_setup();

  /* 2 metavalues: nwritten (pointer) and cumulative sum for block.
     Packsize can be computed from difference between successive
     pointers, if only even packsizes are allowed */
  ptrs = (UINT4 *) CALLOC(((n + BLOCKSIZE)/BLOCKSIZE + 1) * PAIRED_METAINFO_SIZE,sizeof(UINT4));
  ptri = 0;

  if ((comp_fp = FOPEN_WRITE_BINARY(compfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",compfile);
    exit(9);
  }
  buffer = (UINT4 *) CALLOC(buffer_size,sizeof(UINT4));
  buffer_i = 0;

  nwritten = 0U;

  /* Last value of ascending is at ascending[n] */
  /* Use <= n instead of < n, because we want ascending[n] to be taken care of by unpack_00, not a check for remainder == 0 */
  for (positioni = 0; positioni + BLOCKSIZE <= n; positioni += BLOCKSIZE) {
    /* Pointer to D4 */
    ptrs[ptri++] = nwritten/4;	/* In 128-bit registers */

    /* Prefix sum for start of block */
    ptrs[ptri++] = ascending[positioni];
	
    /* D4: Pack block of 64 diffs */
    packsize = compute_q4_diffs_bidir(diffs,&(ascending[positioni]));
    buffer_i = write_columnar(comp_fp,buffer,buffer_size,buffer_i,diffs,packsize);

#ifdef ALLOW_ODD_PACKSIZES
    nwritten += 2 * ((packsize + 1) & ~1);
#else
    nwritten += 2 * packsize;
#endif

    /* Pointer to D1 */
    ptrs[ptri++] = nwritten/4;	/* In 128-bit registers */

    /* D1: Pack block of 64 diffs */
    packsize = compute_q1_diffs(diffs,&(ascending[positioni]));
    buffer_i = write_vert(comp_fp,buffer,buffer_size,buffer_i,diffs,packsize);

#ifdef ALLOW_ODD_PACKSIZES
    nwritten += 2 * ((packsize + 1) & ~1);
#else
    nwritten += 2 * packsize;
#endif
  }

  /* Old: Check for positioni < n, because if positioni == n, ascending[n] will be taken care of as metainfo */
  /* Use <= n instead of < n, because we want ascending[n] to be taken care of by unpack_00, not a check for remainder == 0 */
  if (positioni <= n) {
    /* Finish last block of 64 */
    /* Pointer to D4 */
    ptrs[ptri++] = nwritten/4;	/* In 128-bit registers */

    /* Prefix sum for start of block */
    ptrs[ptri++] = ascending[positioni];

    /* For differential, want <=.  For direct, want < */
    for (i = 0; i <= (int) (n - positioni); i++) {
      last_block[i] = ascending[positioni+i];
    }
    for ( ; i <= BLOCKSIZE; i++) {
      /* Copy last value for rest of block */
      last_block[i] = ascending[n];
    }

    /* D4: Pack block of < 64 diffs */
    packsize = compute_q4_diffs_bidir(diffs,last_block);
    buffer_i = write_columnar(comp_fp,buffer,buffer_size,buffer_i,diffs,packsize);

#ifdef ALLOW_ODD_PACKSIZES
    nwritten += 2 * ((packsize + 1) & ~1);
#else
    nwritten += 2 * packsize;
#endif

    /* Pointer to D1 */
    ptrs[ptri++] = nwritten/4;	/* In 128-bit registers */

    /* D1: Pack block of < 64 diffs */
    packsize = compute_q1_diffs(diffs,last_block);
    buffer_i = write_vert(comp_fp,buffer,buffer_size,buffer_i,diffs,packsize);
  }

  /* Write the final pointer, which will point after the end of the file */
  ptrs[ptri++] = nwritten/4;	/* In 128-bit registers */

  /* Prefix sum for end of block */
  ptrs[ptri++] = ascending[n];

  if ((ptrs_fp = FOPEN_WRITE_BINARY(ptrsfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",ptrsfile);
    exit(9);
  } else {
    FWRITE_UINTS(ptrs,ptri,ptrs_fp);
    FREE(ptrs);
    fclose(ptrs_fp);
  }
    
  /* Empty buffer */
  if (buffer_i > 0) {
    FWRITE_UINTS(buffer,buffer_i,comp_fp);	
    buffer_i = 0;
  }
  FREE(buffer);
  fclose(comp_fp);

  return;
}




/* Worst case:
   64 128 192 256
   256 256 256 256 */

#define FIXED10_PACKSIZE 10 	/* Enough to hold +/- 256 */

/* We want to store values 0..n, with final value at ascending[n]
   possibly stored as the final metainfo value */
/* Stored in columnar order */
void
Bitpack64_write_fixed10 (char *ptrsfile, char *compfile, UINT4 *ascending, Oligospace_T n) {
#ifndef USE_ONE_FILE_FOR_FIXED
  FILE *ptrs_fp;
#endif
  FILE *comp_fp;
  UINT4 *ptrs;
  int ptri, i;
  Oligospace_T positioni;

  /* Buffer is used to avoid frequent writes to the file */
  UINT4 *buffer;
  int buffer_size = BUFFER_SIZE;
  int buffer_i;

  UINT4 diffs[BLOCKSIZE], last_block[BLOCKSIZE+1];

  UINT4 nwritten;
  int packsize;

  write_setup();

  /* 2 metavalues: nwritten (pointer) and cumulative sum for block.
     Packsize can be computed from difference between successive
     pointers, if only even packsizes are allowed */
#ifdef USE_ONE_FILE_FOR_FIXED
  ptrs = (UINT4 *) CALLOC(4,sizeof(UINT4));
  ptri = 0;
#else
  ptrs = (UINT4 *) CALLOC(((n + BLOCKSIZE)/BLOCKSIZE + 1) * RANK_METAINFO_SIZE,sizeof(UINT4));
  ptri = 0;
#endif

  if ((comp_fp = FOPEN_WRITE_BINARY(compfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",compfile);
    exit(9);
  }
  buffer = (UINT4 *) CALLOC(buffer_size,sizeof(UINT4));
  buffer_i = 0;

  nwritten = 0U;

  /* Last value of ascending is at ascending[n] */
  /* Use <= n instead of < n, because we want ascending[n] to be taken care of by unpack_00, not a check for remainder == 0 */
  for (positioni = 0; positioni + BLOCKSIZE <= n; positioni += BLOCKSIZE) {
#if 0
    /* Pointer */
    ptrs[ptri++] = nwritten/4;	/* In 128-bit registers */
#endif

    /* Value for start of block */
    ptrs[ptri++] = ascending[positioni];
#ifdef USE_ONE_FILE_FOR_FIXED
    if (ptri == 4) {
      FWRITE_UINTS(ptrs,4,comp_fp);
      ptri = 0;
    }
#endif
	
    /* Pack block of 64 diffs */
    packsize = compute_q4_diffs_bidir(diffs,&(ascending[positioni]));
    assert(packsize <= FIXED10_PACKSIZE);
    buffer_i = write_columnar(comp_fp,buffer,buffer_size,buffer_i,diffs,FIXED10_PACKSIZE);

    nwritten += 2 * FIXED10_PACKSIZE;
  }

  /* Old: Check for positioni < n, because if positioni == n, ascending[n] will be taken care of as metainfo */
  /* Use <= n instead of < n, because we want ascending[n] to be taken care of by unpack_00, not a check for remainder == 0 */
  if (positioni <= n) {
#if 0
    /* Finish last block of 64 */
    ptrs[ptri++] = nwritten/4;	/* In 128-bit registers */
#endif

    /* Value for start of block */
    ptrs[ptri++] = ascending[positioni];
#ifdef USE_ONE_FILE_FOR_FIXED
    if (ptri == 4) {
      FWRITE_UINTS(ptrs,4,comp_fp);
      ptri = 0;
    }
#endif

    /* For differential, want <=.  For direct, want < */
    for (i = 0; i <= (int) (n - positioni); i++) {
      last_block[i] = ascending[positioni+i];
    }
    for ( ; i <= BLOCKSIZE; i++) {
      /* Copy last value for rest of block */
      last_block[i] = ascending[n];
    }

    /* Pack block of < 64 diffs */
    packsize = compute_q4_diffs_bidir(diffs,last_block);
    assert(packsize <= FIXED10_PACKSIZE);
    buffer_i = write_columnar(comp_fp,buffer,buffer_size,buffer_i,diffs,FIXED10_PACKSIZE);

    nwritten += 2 * FIXED10_PACKSIZE;
  }


#if 0
  /* Write the final pointer, which will point after the end of the file */
  ptrs[ptri++] = nwritten/4;	/* In 128-bit registers */
#endif

  /* Value for end of block */
  ptrs[ptri++] = ascending[n];
#ifdef USE_ONE_FILE_FOR_FIXED
  for (i = ptri; i < 4; i++) {
    ptrs[i] = 0U;
  }
  FWRITE_UINTS(ptrs,4,comp_fp);
#else
  if ((ptrs_fp = FOPEN_WRITE_BINARY(ptrsfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",ptrsfile);
    exit(9);
  } else {
    FWRITE_UINTS(ptrs,ptri,ptrs_fp);
    fclose(ptrs_fp);
  }
#endif
  FREE(ptrs);
    
  /* Empty buffer */
  if (buffer_i > 0) {
    FWRITE_UINTS(buffer,buffer_i,comp_fp);	
    buffer_i = 0;
  }
  FREE(buffer);
  fclose(comp_fp);

  return;
}


void
Bitpack64_write_differential_huge (char *pagesfile, char *ptrsfile, char *compfile,
				   UINT8 *ascending, Oligospace_T n) {
  UINT8 currpage, nextpage;
  FILE *pages_fp, *ptrs_fp, *comp_fp;
  UINT4 pages[25];	/* Allows us to handle up to 100 billion positions */
  UINT4 *ptrs;
  int ptri;
  Oligospace_T positioni;

  /* Buffer is used to avoid frequent writes to the file */
  UINT4 *buffer;
  int buffer_size = BUFFER_SIZE;
  int buffer_i;

  UINT4 diffs[BLOCKSIZE];
  UINT8 last_block[BLOCKSIZE+1];

  int pagei = 0, i;
  UINT4 nwritten;
  int packsize;


  write_setup();

  /* 2 metavalues: nwritten (pointer) and cumulative sum for block.
     Packsize can be computed from difference between successive
     pointers, if only even packsizes are allowed */
  ptrs = (UINT4 *) CALLOC(((n + BLOCKSIZE)/BLOCKSIZE + 1) * DIFFERENTIAL_METAINFO_SIZE,sizeof(UINT4));
  ptri = 0;

  if ((comp_fp = FOPEN_WRITE_BINARY(compfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",compfile);
    exit(9);
  }
  buffer = (UINT4 *) CALLOC(buffer_size,sizeof(UINT4));
  buffer_i = 0;

  currpage = 0;
  nextpage = POSITIONS_PAGE;
  nwritten = 0U;

  /* Last value of ascending is at ascending[n] */
  /* Use <= n instead of < n, because we want ascending[n] to be taken care of by unpack_00, not a check for remainder == 0 */
  for (positioni = 0; positioni + BLOCKSIZE <= n; positioni += BLOCKSIZE) {
    /* Pointer */
    ptrs[ptri++] = nwritten/4;	/* In 128-bit registers */

    /* Value for start of block */
    while (ascending[positioni] >= nextpage) {
      fprintf(stderr,"\nAt position %u (block %u), ascending %llu >= nextpage %llu",
	      positioni,positioni/BLOCKSIZE,ascending[positioni],nextpage);
      pages[pagei++] = positioni/BLOCKSIZE;
      currpage = nextpage;
      nextpage += POSITIONS_PAGE;
    }
    ptrs[ptri++] = ascending[positioni] - currpage;

	
    /* Pack block of 64 diffs */
    packsize = compute_q4_diffs_bidir_huge(diffs,&(ascending[positioni]));
    buffer_i = write_columnar(comp_fp,buffer,buffer_size,buffer_i,diffs,packsize);

#ifdef ALLOW_ODD_PACKSIZES
    nwritten += 2 * ((packsize + 1) & ~1);
#else
    nwritten += 2 * packsize;
#endif
  }

  /* Old: Check for positioni < n, because if positioni == n, ascending[n] will be taken care of as metainfo */
  /* Use <= n instead of < n, because we want ascending[n] to be taken care of by unpack_00, not a check for remainder == 0 */
  if (positioni <= n) {
    /* Finish last block of 64 */
    ptrs[ptri++] = nwritten/4;	/* In 128-bit registers */

    /* Value for start of block */
    while (ascending[positioni] >= nextpage) {
      fprintf(stderr,"\nAt position %u (block %u), ascending %llu >= nextpage %llu",
	      positioni,positioni/BLOCKSIZE,ascending[positioni],nextpage);
      pages[pagei++] = positioni/BLOCKSIZE;
      currpage = nextpage;
      nextpage += POSITIONS_PAGE;
    }
    ptrs[ptri++] = ascending[positioni] - currpage;

    /* For differential, want <=.  For direct, want < */
    for (i = 0; i <= (int) (n - positioni); i++) {
      last_block[i] = ascending[positioni+i] - currpage;
    }
    for ( ; i <= BLOCKSIZE; i++) {
      /* Copy last value for rest of block */
      last_block[i] = ascending[n] - currpage;
    }

    /* Pack block of < 64 diffs */
    packsize = compute_q4_diffs_bidir_huge(diffs,last_block);
    buffer_i = write_columnar(comp_fp,buffer,buffer_size,buffer_i,diffs,packsize);

#ifdef ALLOW_ODD_PACKSIZES
    nwritten += 2 * ((packsize + 1) & ~1);
#else
    nwritten += 2 * packsize;
#endif
  }


  /* Write the final pointer, which will point after the end of the file */
  ptrs[ptri++] = nwritten/4;	/* In 128-bit registers */

  /* Value for end of block */
  if (ascending[n] >= nextpage) {
    fprintf(stderr,"\nAt final oligo %u (block %u), ascending %llu >= nextpage %llu",
	    n,n/BLOCKSIZE,ascending[n],nextpage);
    pages[pagei++] = n/BLOCKSIZE;
    currpage = nextpage;
    /* nextpage += POSITIONS_PAGE; */
  }
  ptrs[ptri++] = ascending[n] - currpage;


  /* Write pages */
  if (pagei > 0) {
    pages[pagei++] = -1U; /* Final value */
    if ((pages_fp = FOPEN_WRITE_BINARY(pagesfile)) == NULL) {
      fprintf(stderr,"Can't write to file %s\n",pagesfile);
      exit(9);
    } else {
      fprintf(stderr,"\nHave %d pages:",pagei);
      for (i = 0; i < pagei; i++) {
	fprintf(stderr," %u",pages[i]);
      }
      fprintf(stderr,"\n");
      FWRITE_UINTS(pages,pagei,pages_fp);
      /* FREE(pages); */
      fclose(pages_fp);
    }
  }

  if ((ptrs_fp = FOPEN_WRITE_BINARY(ptrsfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",ptrsfile);
    exit(9);
  } else {
    FWRITE_UINTS(ptrs,ptri,ptrs_fp);
    FREE(ptrs);
    fclose(ptrs_fp);
  }
    
  /* Empty buffer */
  if (buffer_i > 0) {
    FWRITE_UINTS(buffer,buffer_i,comp_fp);	
    buffer_i = 0;
  }
  FREE(buffer);
  fclose(comp_fp);

  return;
}


void
Bitpack64_write_fixed10_huge (char *pagesfile, char *ptrsfile, char *compfile,
			      UINT8 *ascending, Oligospace_T n) {
#ifndef USE_ONE_FILE_FOR_FIXED
  FILE *ptrs_fp;
#endif
  UINT8 currpage, nextpage;
  FILE *pages_fp, *comp_fp;
  UINT4 pages[25];	/* Allows us to handle up to 100 billion positions */
  UINT4 *ptrs;
  int ptri;
  Oligospace_T positioni;

  /* Buffer is used to avoid frequent writes to the file */
  UINT4 *buffer;
  int buffer_size = BUFFER_SIZE;
  int buffer_i;

  UINT4 diffs[BLOCKSIZE];
  UINT8 last_block[BLOCKSIZE+1];

  int pagei = 0, i;
  UINT4 nwritten;
  int packsize;


  write_setup();

  /* 2 metavalues: nwritten (pointer) and cumulative sum for block.
     Packsize can be computed from difference between successive
     pointers, if only even packsizes are allowed */
#ifdef USE_ONE_FILE_FOR_FIXED
  ptrs = (UINT *) CALLOC(4,sizeof(UINT4));
  ptri = 0;
#else
  ptrs = (UINT4 *) CALLOC(((n + BLOCKSIZE)/BLOCKSIZE + 1) * RANK_METAINFO_SIZE,sizeof(UINT4));
  ptri = 0;
#endif

  if ((comp_fp = FOPEN_WRITE_BINARY(compfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",compfile);
    exit(9);
  }
  buffer = (UINT4 *) CALLOC(buffer_size,sizeof(UINT4));
  buffer_i = 0;

  currpage = 0;
  nextpage = POSITIONS_PAGE;
  nwritten = 0U;

  /* Last value of ascending is at ascending[n] */
  /* Use <= n instead of < n, because we want ascending[n] to be taken care of by unpack_00, not a check for remainder == 0 */
  for (positioni = 0; positioni + BLOCKSIZE <= n; positioni += BLOCKSIZE) {
#if 0
    /* Pointer */
    ptrs[ptri++] = nwritten/4;	/* In 128-bit registers */
#endif

    /* Value for start of block */
    while (ascending[positioni] >= nextpage) {
      fprintf(stderr,"\nAt position %u (block %u), ascending %llu >= nextpage %llu",
	      positioni,positioni/BLOCKSIZE,ascending[positioni],nextpage);
      pages[pagei++] = positioni/BLOCKSIZE;
      currpage = nextpage;
      nextpage += POSITIONS_PAGE;
    }
    ptrs[ptri++] = ascending[positioni] - currpage;
#ifdef USE_ONE_FILE_FOR_FIXED
    if (ptri == 4) {
      FWRITE_UINTS(ptrs,4,comp_fp);
      ptri = 0;
    }
#endif
	
    /* Pack block of 64 diffs */
    packsize = compute_q4_diffs_bidir_huge(diffs,&(ascending[positioni]));
    assert(packsize <= FIXED10_PACKSIZE);
    buffer_i = write_columnar(comp_fp,buffer,buffer_size,buffer_i,diffs,FIXED10_PACKSIZE);

    nwritten += 2 * packsize;
  }

  /* Old: Check for positioni < n, because if positioni == n, ascending[n] will be taken care of as metainfo */
  /* Use <= n instead of < n, because we want ascending[n] to be taken care of by unpack_00, not a check for remainder == 0 */
  if (positioni <= n) {
#if 0
    /* Finish last block of 64 */
    ptrs[ptri++] = nwritten/4;	/* In 128-bit registers */
#endif

    /* Value for start of block */
    while (ascending[positioni] >= nextpage) {
      fprintf(stderr,"\nAt position %u (block %u), ascending %llu >= nextpage %llu",
	      positioni,positioni/BLOCKSIZE,ascending[positioni],nextpage);
      pages[pagei++] = positioni/BLOCKSIZE;
      currpage = nextpage;
      nextpage += POSITIONS_PAGE;
    }
    ptrs[ptri++] = ascending[positioni] - currpage;
#ifdef USE_ONE_FILE_FOR_FIXED
    if (ptri == 4) {
      FWRITE_UINTS(ptrs,4,comp_fp);
      ptri = 0;
    }
#endif

    /* For differential, want <=.  For direct, want < */
    for (i = 0; i <= (int) (n - positioni); i++) {
      last_block[i] = ascending[positioni+i] - currpage;
    }
    for ( ; i <= BLOCKSIZE; i++) {
      /* Copy last value for rest of block */
      last_block[i] = ascending[n] - currpage;
    }

    /* Pack block of < 64 diffs */
    packsize = compute_q4_diffs_bidir_huge(diffs,last_block);
    assert(packsize <= FIXED10_PACKSIZE);
    buffer_i = write_columnar(comp_fp,buffer,buffer_size,buffer_i,diffs,FIXED10_PACKSIZE);

    nwritten += 2 * packsize;
  }


#if 0
  /* Write the final pointer, which will point after the end of the file */
  ptrs[ptri++] = nwritten/4;	/* In 128-bit registers */
#endif

  /* Value for end of block */
  if (ascending[n] >= nextpage) {
    fprintf(stderr,"\nAt final oligo %u (block %u), ascending %llu >= nextpage %llu",
	    n,n/BLOCKSIZE,ascending[n],nextpage);
    pages[pagei++] = n/BLOCKSIZE;
    currpage = nextpage;
    /* nextpage += POSITIONS_PAGE; */
  }
  ptrs[ptri++] = ascending[n] - currpage;
#ifdef USE_ONE_FILE_FOR_FIXED
  for (i = ptri; i < 4; i++) {
    ptrs[i] = 0U;
  }
  FWRITE_UINTS(ptrs,4,comp_fp);
#else
  if ((ptrs_fp = FOPEN_WRITE_BINARY(ptrsfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",ptrsfile);
    exit(9);
  } else {
    FWRITE_UINTS(ptrs,ptri,ptrs_fp);
    fclose(ptrs_fp);
  }
#endif
  FREE(ptrs);

  /* Write pages */
  if (pagei > 0) {
    pages[pagei++] = -1U; /* Final value */
    if ((pages_fp = FOPEN_WRITE_BINARY(pagesfile)) == NULL) {
      fprintf(stderr,"Can't write to file %s\n",pagesfile);
      exit(9);
    } else {
      fprintf(stderr,"\nHave %d pages:",pagei);
      for (i = 0; i < pagei; i++) {
	fprintf(stderr," %u",pages[i]);
      }
      fprintf(stderr,"\n");
      FWRITE_UINTS(pages,pagei,pages_fp);
      /* FREE(pages); */
      fclose(pages_fp);
    }
  }

  /* Empty buffer */
  if (buffer_i > 0) {
    FWRITE_UINTS(buffer,buffer_i,comp_fp);	
    buffer_i = 0;
  }
  FREE(buffer);
  fclose(comp_fp);

  return;
}




static int
compute_packsize (UINT4 *values) {
  UINT4 packsize;
  UINT4 maxvalue = 0;
  int i;
  int firstbit, msb;

  for (i = 0; i < 64; i++) {
    maxvalue |= values[i];
  }

  if (maxvalue == 0) {
    /* __builtin_clz() behaves oddly on zero */
    return 0;

  } else {
#ifdef HAVE_BUILTIN_CLZ
    firstbit = __builtin_clz(maxvalue);
    packsize = 32 - firstbit;
#elif defined(HAVE_ASM_BSR)
    asm("bsr %1,%0" : "=r"(msb) : "r"(maxvalue));
    packsize = msb + 1;
#else
    firstbit = ((maxvalue >> 16) ? clz_table[maxvalue >> 16] : 16 + clz_table[maxvalue]);
    packsize = 32 - firstbit;
#endif

#ifdef ALLOW_ODD_PACKSIZES
    return packsize;
#else
    return (packsize + 1) & ~1;	/* Converts packsizes to the next multiple of 2 */
#endif
  }
}


/* Want to store values 0..n-1.  The value direct[n] does not exist.  */
/* Stored in vertical order */
void
Bitpack64_write_direct (char *ptrsfile, char *compfile, UINT4 *direct, Oligospace_T n) {
  FILE *ptrs_fp, *comp_fp;
  UINT4 *ptrs;
  int ptri, i;
  Oligospace_T positioni;

  UINT4 *buffer;
  int buffer_size = BUFFER_SIZE;
  int buffer_i;

  UINT4 last_block[BLOCKSIZE];

  UINT4 nwritten;
  int packsize;


  write_setup();

  /* 1 metavalue: nwritten (pointer).  Packsize can be
     computed from difference between successive pointers, if only
     even packsizes are allowed */
  ptrs = (UINT4 *) CALLOC(((n + BLOCKSIZE - 1)/BLOCKSIZE + 1) * DIRECT_METAINFO_SIZE,sizeof(UINT4));
  ptri = 0;

  if ((comp_fp = FOPEN_WRITE_BINARY(compfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",compfile);
    exit(9);
  }
  buffer = (UINT4 *) CALLOC(buffer_size,sizeof(UINT4));
  buffer_i = 0;

  nwritten = 0U;

  for (positioni = 0; positioni + BLOCKSIZE < n; positioni += BLOCKSIZE) {
    /* Pointer */
    ptrs[ptri++] = nwritten/4;	/* In 128-bit registers */

    /* Pack block of 64 diffs */
    packsize = compute_packsize(&(direct[positioni]));
    buffer_i = write_vert(comp_fp,buffer,buffer_size,buffer_i,&(direct[positioni]),packsize);

#ifdef ALLOW_ODD_PACKSIZES
    nwritten += 2 * ((packsize + 1) & ~1);
#else
    nwritten += 2 * packsize;
#endif
  }

  if (positioni < n) {
    /* Finish last block of 64 */
    ptrs[ptri++] = nwritten/4;	/* In 128-bit registers */
    
    i = 0;
    while (positioni < n) {
      last_block[i++] = direct[positioni++];
    }
    while (i < BLOCKSIZE) {
      last_block[i++] = 0;
    }

    packsize = compute_packsize(last_block);
    buffer_i = write_vert(comp_fp,buffer,buffer_size,buffer_i,last_block,packsize);

#ifdef ALLOW_ODD_PACKSIZES
    nwritten += 2 * ((packsize + 1) & ~1);
#else
    nwritten += 2 * packsize;
#endif
  }

  /* Write the final pointer, which will point after the end of the
     file */
  ptrs[ptri++] = nwritten/4;	/* In 128-bit registers */

  if ((ptrs_fp = FOPEN_WRITE_BINARY(ptrsfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",ptrsfile);
    exit(9);
  } else {
    FWRITE_UINTS(ptrs,ptri,ptrs_fp);
    FREE(ptrs);
    fclose(ptrs_fp);
  }
    
  /* Empty buffer */
  if (buffer_i > 0) {
    FWRITE_UINTS(buffer,buffer_i,comp_fp);	
    buffer_i = 0;
  }
  FREE(buffer);
  fclose(comp_fp);

  return;
}


