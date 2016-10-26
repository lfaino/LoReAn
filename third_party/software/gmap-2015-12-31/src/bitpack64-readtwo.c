static char rcsid[] = "$Id: bitpack64-readtwo.c 168395 2015-06-26 17:13:13Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "bitpack64-readtwo.h"

#include <stdio.h>
#include <stdlib.h>

#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
#elif defined(HAVE_SSE2)
#include <emmintrin.h>
#endif

#define POSITIONS_PAGE 4294967296 /* 2^32 */


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


/* #define ALLOW_ODD_PACKSIZES 1 */

/* Two ideas for branch-free code:
   BRANCH_FREE_ROW_SUM simplifies the summation of the difference values to always add 4 values.
   BRANCH_FREE_QTR_BLOCK avoids having if statements based on the quarter-block.

   If BRANCH_FREE_QTR_BLOCK is selected, then BRANCH_FREE_ROW_SUM must also be selected, so

   Case 1: BRANCH_FREE_ROW_SUM 0, BRANCH_FREE_QTR_BLOCK 0
   Case 2: BRANCH_FREE_ROW_SUM 1, BRANCH_FREE_QTR_BLOCK 0
   Case 3: BRANCH_FREE_ROW_SUM 1, BRANCH_FREE_QTR_BLOCK 1

   Note that BRANCH_FREE_QTR_BLOCK can be tricky for 8-byte quantities, e.g.,
   in Bitpack64_read_one_huge.  Would therefore recommend it be turned off.
*/

/* #define BRANCH_FREE_ROW_SUM 1 -- Not supported here */
/* #define BRANCH_FREE_QTR_BLOCK 1 */

#ifdef DEBUG
#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
#else
/* For debugging */
static void
print_vector_hex (__m128i x) {
  UINT4 *s = (UINT4 *) &x;

  printf("%08X %08X %08X %08X\n",s[0],s[1],s[2],s[3]);
  return;
}

static void
print_vector (__m128i x) {
  UINT4 *s = (UINT4 *) &x;

  printf("%u %u %u %u\n",s[0],s[1],s[2],s[3]);
  return;
}
#endif
#endif


#define BLOCKSIZE 64

#if 0
void
Bitpack64_read_setup () {

#ifdef HAVE_SSE2
#ifdef ALLOW_ODD_PACKSIZES
  mask1 = _mm_set1_epi32(1U);
  mask3 =  _mm_set1_epi32(7U);
  mask5 =  _mm_set1_epi32(31U);
  mask7 =  _mm_set1_epi32(127U);
  mask9 =  _mm_set1_epi32(511U);
  mask11 =  _mm_set1_epi32(2047U);
  mask13 =  _mm_set1_epi32(8191U);
  mask15 =  _mm_set1_epi32(32767U);
  mask17 =  _mm_set1_epi32(131071U);
  mask19 =  _mm_set1_epi32(524287U);
  mask21 =  _mm_set1_epi32(2097151U);
  mask23 =  _mm_set1_epi32(8388607U);
  mask25 =  _mm_set1_epi32(33554431U);
  mask27 =  _mm_set1_epi32(134217727U);
  mask29 =  _mm_set1_epi32(536870911U);
  mask31 =  _mm_set1_epi32(2147483647U);
#endif
  mask2 = _mm_set1_epi32(3U);
  mask4 =  _mm_set1_epi32(15U);
  mask6 =  _mm_set1_epi32(63U);
  mask8 =  _mm_set1_epi32(255U);
  mask10 =  _mm_set1_epi32(1023U);
  mask12 =  _mm_set1_epi32(4095U);
  mask14 =  _mm_set1_epi32(16383U);
  mask16 =  _mm_set1_epi32(65535U);
  mask18 =  _mm_set1_epi32(262143U);
  mask20 =  _mm_set1_epi32(1048575U);
  mask22 =  _mm_set1_epi32(4194303U);
  mask24 =  _mm_set1_epi32(16777215U);
  mask26 =  _mm_set1_epi32(67108863U);
  mask28 =  _mm_set1_epi32(268435455U);
  mask30 =  _mm_set1_epi32(1073741823U);
#endif

  return;
}
#endif


#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
static void
unpack_00 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  int i;

  for (i = 0; i < BLOCKSIZE; i++) {
    *out++ = 0;
  }

  return;
}

#else
static void
unpack_00 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i total = _mm_set1_epi32(0U);

  _mm_store_si128(out++, total);
  _mm_store_si128(out++, total);
  _mm_store_si128(out++, total);
  _mm_store_si128(out++, total);
  _mm_store_si128(out++, total);
  _mm_store_si128(out++, total);
  _mm_store_si128(out++, total);
  _mm_store_si128(out++, total);

  return;
}

/* Handles the case where remainder == 64 => column 3, row -1 */
static void
unpack_00_0 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {

  /* _mm_store_si128(out++, zero); -- Not needed, since row == -1 */

  return;
}


static void
unpack_00_1_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i zero = _mm_set1_epi32(0U);

  /* 1 */
  _mm_store_si128(out++, zero);

  /* Skip row */
  out++;

  /* 3 */
  _mm_store_si128(out++, zero);

  return;
}

static void
unpack_00_2_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i zero = _mm_set1_epi32(0U);

  /* 2 */
  _mm_store_si128(out++, zero);
  _mm_store_si128(out++, zero);

  /* 4 */
  _mm_store_si128(out++, zero);
  _mm_store_si128(out++, zero);

  return;
}
#endif


#ifdef ALLOW_ODD_PACKSIZES
static void
unpack_01 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg;
    const __m128i mask1 = _mm_set1_epi32(1U);

    OutReg = _mm_and_si128( InReg , mask1);
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,1) , mask1);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask1);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,3) , mask1);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask1);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,5) , mask1);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask1);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,7) , mask1);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask1);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,9) , mask1);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,10) , mask1);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,11) , mask1);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask1);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,13) , mask1);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,14) , mask1);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,15) ;
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}
#endif


#ifdef WORDS_BIGENDIAN
static void
unpack_02 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    *out = ( Bigendian_convert_uint(*in) >>  0  )   % (1U << 2 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  2  )   % (1U << 2 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  4  )   % (1U << 2 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  6  )   % (1U << 2 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  8  )   % (1U << 2 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  10  )   % (1U << 2 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  12  )   % (1U << 2 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  14  )   % (1U << 2 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  16  )   % (1U << 2 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  18  )   % (1U << 2 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  20  )   % (1U << 2 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  22  )   % (1U << 2 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  24  )   % (1U << 2 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  26  )   % (1U << 2 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  28  )   % (1U << 2 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  30  )   % (1U << 2 ) ;
    out++;
  }

  return;
}

#elif !defined(HAVE_SSE2)
static void
unpack_02 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    *out = ( (*in) >>  0  )   % (1U << 2 ) ;
    out++;
    *out = ( (*in) >>  2  )   % (1U << 2 ) ;
    out++;
    *out = ( (*in) >>  4  )   % (1U << 2 ) ;
    out++;
    *out = ( (*in) >>  6  )   % (1U << 2 ) ;
    out++;
    *out = ( (*in) >>  8  )   % (1U << 2 ) ;
    out++;
    *out = ( (*in) >>  10  )   % (1U << 2 ) ;
    out++;
    *out = ( (*in) >>  12  )   % (1U << 2 ) ;
    out++;
    *out = ( (*in) >>  14  )   % (1U << 2 ) ;
    out++;
    *out = ( (*in) >>  16  )   % (1U << 2 ) ;
    out++;
    *out = ( (*in) >>  18  )   % (1U << 2 ) ;
    out++;
    *out = ( (*in) >>  20  )   % (1U << 2 ) ;
    out++;
    *out = ( (*in) >>  22  )   % (1U << 2 ) ;
    out++;
    *out = ( (*in) >>  24  )   % (1U << 2 ) ;
    out++;
    *out = ( (*in) >>  26  )   % (1U << 2 ) ;
    out++;
    *out = ( (*in) >>  28  )   % (1U << 2 ) ;
    out++;
    *out = ( (*in) >>  30  )   % (1U << 2 ) ;
    out++;
  }

  return;
}

#else
static void
unpack_02_fwd (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg;
    const __m128i mask2 = _mm_set1_epi32(3U);

    OutReg = _mm_and_si128( InReg , mask2);
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask2);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask2);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask2);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask2);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,10) , mask2);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask2);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,14) , mask2);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_02_fwd_1_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask2 = _mm_set1_epi32(3U);


    /* 1 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask2);
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 3 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask2);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_02_fwd_2_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask2 = _mm_set1_epi32(3U);

    /* 2 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask2);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask2);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);


    /* 4 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask2);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask2);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_02_fwd_3_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask2 = _mm_set1_epi32(3U);

    /* 3 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask2);
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 5 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask2);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_02_fwd_4_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask2 = _mm_set1_epi32(3U);

    /* 4 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask2);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask2);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    /* 6 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask2);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,10) , mask2);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}


static void
unpack_02_fwd_5_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask2 = _mm_set1_epi32(3U);

    /* 5 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask2);
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 7 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask2);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_02_fwd_6_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask2 = _mm_set1_epi32(3U);

    /* 6 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask2);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,10) , mask2);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    /* 8 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask2);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,14) , mask2);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_02_fwd_7_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask2 = _mm_set1_epi32(3U);

    /* 7 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask2);
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 1 */
    total = /* OutReg = */ _mm_and_si128( InReg , mask2);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_02_fwd_8_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask2 = _mm_set1_epi32(3U);

    /* 8 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask2);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,14) , mask2);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    /* 2 */
    total = /* OutReg = */ _mm_and_si128( InReg , mask2);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask2);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_02_rev (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg;
    const __m128i mask2 = _mm_set1_epi32(3U);

    InReg = _mm_load_si128(in);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask2);
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,18) , mask2);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,20) , mask2);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,22) , mask2);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,24) , mask2);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,26) , mask2);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,28) , mask2);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_02_rev_1_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask2 = _mm_set1_epi32(3U);

    /* 1 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask2);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 3 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,20) , mask2);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_02_rev_2_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask2 = _mm_set1_epi32(3U);

    /* 2 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask2);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,18) , mask2);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    /* 4 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,20) , mask2);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,22) , mask2);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_02_rev_3_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask2 = _mm_set1_epi32(3U);

    /* 3 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,20) , mask2);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 5 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,24) , mask2);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_02_rev_4_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask2 = _mm_set1_epi32(3U);

    /* 4 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,20) , mask2);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,22) , mask2);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    /* 6 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,24) , mask2);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,26) , mask2);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_02_rev_5_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask2 = _mm_set1_epi32(3U);

    /* 5 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,24) , mask2);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 7 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,28) , mask2);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_02_rev_6_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask2 = _mm_set1_epi32(3U);

    /* 6 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,24) , mask2);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,26) , mask2);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    /* 8 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,28) , mask2);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,30) ;
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_02_rev_7_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask2 = _mm_set1_epi32(3U);

    /* 7 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,28) , mask2);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 1 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask2);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_02_rev_8_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask2 = _mm_set1_epi32(3U);

    /* 8 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,28) , mask2);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,30) ;
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    /* 2 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask2);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,18) , mask2);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

#endif




#ifdef ALLOW_ODD_PACKSIZES
static void
unpack_03 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg;
    const __m128i mask3 =  _mm_set1_epi32(7U);

    OutReg = _mm_and_si128( InReg , mask3);
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,3) , mask3);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask3);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,9) , mask3);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask3);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,15) , mask3);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,18) , mask3);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,21) , mask3);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,24) , mask3);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,27) , mask3);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 3 - 1), mask3));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,1) , mask3);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask3);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,7) , mask3);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,10) , mask3);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,13) , mask3);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}
#endif



#ifdef WORDS_BIGENDIAN
static void
unpack_04 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  UINT4 outer, inwordpointer;
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    for (outer = 0; outer < 2 ; outer++) {
      for (inwordpointer = 0; inwordpointer < 32; inwordpointer +=  4) {
	*(out++) = ( Bigendian_convert_uint(*in) >> inwordpointer )   % (1U << 4 ) ;
      }
      in += 4;
    }
  }

  return;
}

#elif !defined(HAVE_SSE2)
static void
unpack_04 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  UINT4 outer, inwordpointer;
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    for (outer = 0; outer < 2 ; outer++) {
      for (inwordpointer = 0; inwordpointer < 32; inwordpointer +=  4) {
	*(out++) = ( (*in) >> inwordpointer )   % (1U << 4 ) ;
      }
      in += 4;
    }
  }

  return;
}

#else
static void
unpack_04_fwd (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg;
    const __m128i mask4 = _mm_set1_epi32(15U);

    OutReg = _mm_and_si128( InReg , mask4);
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask4);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask4);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask4);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask4);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,20) , mask4);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,24) , mask4);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_04_fwd_1_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask4 = _mm_set1_epi32(15U);

    /* 1 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask4);
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 3 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask4);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_04_fwd_2_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask4 = _mm_set1_epi32(15U);

    /* 2 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask4);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask4);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    /* 4 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask4);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask4);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_04_fwd_3_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask4 = _mm_set1_epi32(15U);

    /* 3 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask4);
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 5 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask4);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_04_fwd_4_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask4 = _mm_set1_epi32(15U);

    /* 4 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask4);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask4);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    /* 6 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask4);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,20) , mask4);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_04_fwd_5_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask4 = _mm_set1_epi32(15U);

    /* 5 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask4);
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 7 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,24) , mask4);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_04_fwd_6_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask4 = _mm_set1_epi32(15U);

    /* 6 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask4);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,20) , mask4);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    /* 8 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,24) , mask4);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_04_fwd_7_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask4 = _mm_set1_epi32(15U);

    /* 7 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,24) , mask4);
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    total = /* OutReg = */ _mm_and_si128( InReg , mask4);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_04_fwd_8_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask4 = _mm_set1_epi32(15U);

    /* 8 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,24) , mask4);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    /* 2 */
    total = /* OutReg = */ _mm_and_si128( InReg , mask4);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask4);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}


static void
unpack_04_rev (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg;
    const __m128i mask4 = _mm_set1_epi32(15U);

    InReg = _mm_load_si128(++in);

    OutReg = _mm_and_si128( InReg , mask4);
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask4);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask4);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask4);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask4);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,20) , mask4);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,24) , mask4);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_04_rev_1_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask4 = _mm_set1_epi32(15U);

    /* 1 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask4);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 3 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask4);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_04_rev_2_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask4 = _mm_set1_epi32(15U);

    /* 2 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask4);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask4);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    /* 4 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask4);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask4);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_04_rev_3_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask4 = _mm_set1_epi32(15U);

    /* 3 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask4);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 5 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask4);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_04_rev_4_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask4 = _mm_set1_epi32(15U);

    /* 4 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask4);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask4);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    /* 6 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask4);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,20) , mask4);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_04_rev_5_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask4 = _mm_set1_epi32(15U);

    /* 5 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask4);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 7 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,24) , mask4);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_04_rev_6_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask4 = _mm_set1_epi32(15U);

    /* 6 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask4);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,20) , mask4);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    /* 8 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,24) , mask4);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_04_rev_7_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask4 = _mm_set1_epi32(15U);

    /* 7 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,24) , mask4);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 1 */
    total = /* OutReg = */ _mm_and_si128( InReg , mask4);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);


    return;
}

static void
unpack_04_rev_8_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask4 = _mm_set1_epi32(15U);

    /* 8 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,24) , mask4);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    /* 2 */
    total = /* OutReg = */ _mm_and_si128( InReg , mask4);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask4);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}
#endif


#ifdef ALLOW_ODD_PACKSIZES
static void
unpack_05 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg;
    const __m128i  mask5 =  _mm_set1_epi32(31U);

    OutReg = _mm_and_si128( InReg , mask5);
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,5) , mask5);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,10) , mask5);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,15) , mask5);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,20) , mask5);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,25) , mask5);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 5 - 3), mask5));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,3) , mask5);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask5);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,13) , mask5);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,18) , mask5);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,23) , mask5);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 5 - 1), mask5));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,1) , mask5);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask5);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,11) , mask5);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}
#endif



#ifdef WORDS_BIGENDIAN
static void
unpack_06 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    *out = ( Bigendian_convert_uint(*in) >>  0  )   % (1U << 6 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  6  )   % (1U << 6 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  12  )   % (1U << 6 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  18  )   % (1U << 6 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  24  )   % (1U << 6 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  30  )   % (1U << 6 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 4 ))<<( 6 - 4 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  4  )   % (1U << 6 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  10  )   % (1U << 6 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  16  )   % (1U << 6 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  22  )   % (1U << 6 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  28  )   % (1U << 6 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 2 ))<<( 6 - 2 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  2  )   % (1U << 6 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  8  )   % (1U << 6 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  14  )   % (1U << 6 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  20  )   % (1U << 6 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  26  )   % (1U << 6 ) ;
    out++;
  }

  return;
}

#elif !defined(HAVE_SSE2)
static void
unpack_06 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    *out = ( (*in) >>  0  )   % (1U << 6 ) ;
    out++;
    *out = ( (*in) >>  6  )   % (1U << 6 ) ;
    out++;
    *out = ( (*in) >>  12  )   % (1U << 6 ) ;
    out++;
    *out = ( (*in) >>  18  )   % (1U << 6 ) ;
    out++;
    *out = ( (*in) >>  24  )   % (1U << 6 ) ;
    out++;
    *out = ( (*in) >>  30  )   % (1U << 6 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 4 ))<<( 6 - 4 );
    out++;
    *out = ( (*in) >>  4  )   % (1U << 6 ) ;
    out++;
    *out = ( (*in) >>  10  )   % (1U << 6 ) ;
    out++;
    *out = ( (*in) >>  16  )   % (1U << 6 ) ;
    out++;
    *out = ( (*in) >>  22  )   % (1U << 6 ) ;
    out++;
    *out = ( (*in) >>  28  )   % (1U << 6 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 2 ))<<( 6 - 2 );
    out++;
    *out = ( (*in) >>  2  )   % (1U << 6 ) ;
    out++;
    *out = ( (*in) >>  8  )   % (1U << 6 ) ;
    out++;
    *out = ( (*in) >>  14  )   % (1U << 6 ) ;
    out++;
    *out = ( (*in) >>  20  )   % (1U << 6 ) ;
    out++;
    *out = ( (*in) >>  26  )   % (1U << 6 ) ;
    out++;
  }

  return;
}

#else
static void
unpack_06_fwd (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg;
    const __m128i mask6 =  _mm_set1_epi32(63U);

    OutReg = _mm_and_si128( InReg , mask6);
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask6);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask6);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,18) , mask6);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,24) , mask6);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 6 - 4), mask6));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask6);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,10) , mask6);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_06_fwd_1_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i total;
    const __m128i mask6 =  _mm_set1_epi32(63U);

    /* 1 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask6);
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 3 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask6);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_06_fwd_2_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask6 =  _mm_set1_epi32(63U);

    /* 2 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask6);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask6);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    /* 4 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask6);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,18) , mask6);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_06_fwd_3_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask6 =  _mm_set1_epi32(63U);

    /* 3 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask6);
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 5 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,24) , mask6);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_06_fwd_4_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask6 =  _mm_set1_epi32(63U);

    /* 4 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask6);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,18) , mask6);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    /* 6 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,24) , mask6);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 6 - 4), mask6));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_06_fwd_5_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask6 =  _mm_set1_epi32(63U);

    /* 5 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,24) , mask6);
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 7 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask6);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_06_fwd_6_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask6 =  _mm_set1_epi32(63U);

    /* 6 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,24) , mask6);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 6 - 4), mask6));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    /* 8 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask6);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,10) , mask6);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_06_fwd_7_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask6 =  _mm_set1_epi32(63U);

    /* 1 first */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask6);
    _mm_store_si128(&(out[2]), total);

    /* 7 second */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask6);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_06_fwd_8_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask6 =  _mm_set1_epi32(63U);

    /* 2 first */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask6);
    _mm_store_si128(&(out[2]), total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask6);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(&(out[3]), total);

    /* 8 second */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask6);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,10) , mask6);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_06_rev (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg;
    const __m128i mask6 =  _mm_set1_epi32(63U);

    InReg = _mm_load_si128(++in);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask6);
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,22) , mask6);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 6 - 2), mask6));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask6);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask6);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,14) , mask6);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,20) , mask6);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_06_rev_1_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask6 =  _mm_set1_epi32(63U);

    /* 1 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask6);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 3 */
    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 6 - 2), mask6));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_06_rev_2_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask6 =  _mm_set1_epi32(63U);

    /* 2 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask6);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,22) , mask6);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    /* 4 */
    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 6 - 2), mask6));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask6);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_06_rev_3_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask6 =  _mm_set1_epi32(63U);

    /* 3 */
    InReg = _mm_load_si128(++in);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 6 - 2), mask6));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 5 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask6);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_06_rev_4_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask6 =  _mm_set1_epi32(63U);

    /* 4 */
    InReg = _mm_load_si128(++in);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 6 - 2), mask6));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask6);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    /* 6 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask6);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,14) , mask6);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_06_rev_5_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask6 =  _mm_set1_epi32(63U);

    /* 5 */
    in += 2;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask6);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 7 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,20) , mask6);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_06_rev_6_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask6 =  _mm_set1_epi32(63U);

    /* 6 */
    in += 2;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask6);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,14) , mask6);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    /* 8 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,20) , mask6);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,26) ;
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;

}

static void
unpack_06_rev_7_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask6 =  _mm_set1_epi32(63U);

    /* 1 first */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask6);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(&(out[2]), total);

    /* 7 second */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,20) , mask6);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_06_rev_8_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask6 =  _mm_set1_epi32(63U);

    /* 2 first */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask6);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(&(out[2]), total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,22) , mask6);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(&(out[3]), total);

    /* 8 second */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,20) , mask6);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,26) ;
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}
#endif


#ifdef ALLOW_ODD_PACKSIZES
static void
unpack_07 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg;
    const __m128i mask7 =  _mm_set1_epi32(127U);

    OutReg = _mm_and_si128( InReg , mask7);
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,7) , mask7);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,14) , mask7);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,21) , mask7);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 7 - 3), mask7));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,3) , mask7);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,10) , mask7);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,17) , mask7);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,24) , mask7);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,31) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 7 - 6), mask7));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask7);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,13) , mask7);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,20) , mask7);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,27) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 7 - 2), mask7));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask7);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,9) , mask7);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}
#endif



#ifdef WORDS_BIGENDIAN
static void
unpack_08 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  UINT4 outer, inwordpointer;
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    for (outer = 0; outer < 4; outer++) {
      for (inwordpointer = 0; inwordpointer < 32; inwordpointer += 8) {
	*(out++) = ( Bigendian_convert_uint(*in) >> inwordpointer )   % (1U << 8 ) ;
      }
      in += 4;
    }
  }

  return;
}

#elif !defined(HAVE_SSE2)
static void
unpack_08 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  UINT4 outer, inwordpointer;
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    for (outer = 0; outer < 4; outer++) {
      for (inwordpointer = 0; inwordpointer < 32; inwordpointer += 8) {
	*(out++) = ( (*in) >> inwordpointer )   % (1U << 8 ) ;
      }
      in += 4;
    }
  }

  return;
}

#else
static void
unpack_08_fwd (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg;
    const __m128i mask8 =  _mm_set1_epi32(255U);

    OutReg = _mm_and_si128( InReg , mask8);
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask8);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask8);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128( InReg , mask8);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask8);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask8);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_08_fwd_1_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask8 =  _mm_set1_epi32(255U);

    /* 1 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask8);
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 3 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask8);
    _mm_store_si128(out++, total);

    return;
}


static void
unpack_08_fwd_2_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask8 =  _mm_set1_epi32(255U);

    /* 2 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask8);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask8);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    /* 4 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask8);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    /* InReg = _mm_load_si128(++in); */

    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_08_fwd_3_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask8 =  _mm_set1_epi32(255U);

    /* 3 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask8);
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 5 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask8);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_08_fwd_4_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask8 =  _mm_set1_epi32(255U);

    /* 4 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask8);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    /* InReg = _mm_load_si128(++in); */

    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    /* 6 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask8);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask8);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_08_fwd_5_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask8 =  _mm_set1_epi32(255U);

    /* 5 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask8);
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 7 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask8);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_08_fwd_6_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask8 =  _mm_set1_epi32(255U);

    /* 6 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask8);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask8);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    /* 8 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask8);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_08_fwd_7_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask8 =  _mm_set1_epi32(255U);

    /* 1 first */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask8);
    _mm_store_si128(&(out[2]), total);

    /* 7 second */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask8);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_08_fwd_8_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask8 =  _mm_set1_epi32(255U);

    /* 2 first */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask8);
    _mm_store_si128(&(out[2]), total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask8);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(&(out[3]), total);

    /* 8 second */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask8);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}


static void
unpack_08_rev (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg;
    const __m128i mask8 =  _mm_set1_epi32(255U);

    in += 2;
    InReg = _mm_load_si128(in);


    OutReg = _mm_and_si128( InReg , mask8);
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask8);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask8);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128( InReg , mask8);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask8);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask8);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_08_rev_1_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask8 =  _mm_set1_epi32(255U);

    /* 1 */
    in += 2;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask8);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 3 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask8);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_08_rev_2_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask8 =  _mm_set1_epi32(255U);

    /* 2 */
    in += 2;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask8);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask8);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    /* 4 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask8);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    /* InReg = _mm_load_si128(++in); */

#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_08_rev_3_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask8 =  _mm_set1_epi32(255U);

    /* 3 */
    in += 2;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask8);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 5 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask8);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_08_rev_4_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask8 =  _mm_set1_epi32(255U);

    /* 4 */
    in += 2;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask8);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    /* InReg = _mm_load_si128(++in); */

#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);


    /* 6 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask8);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask8);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_08_rev_5_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask8 =  _mm_set1_epi32(255U);

    /* 5 */
    in += 3;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask8);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 7 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask8);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_08_rev_6_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask8 =  _mm_set1_epi32(255U);

    /* 6 */
    in += 3;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask8);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask8);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);


    /* 8 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask8);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_08_rev_7_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask8 =  _mm_set1_epi32(255U);

    /* 1 first */
    in += 2;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask8);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(&(out[2]), total);

    /* 7 second */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask8);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_08_rev_8_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask8 =  _mm_set1_epi32(255U);

    /* 2 first */
    in += 2;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask8);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(&(out[2]), total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask8);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(&(out[3]), total);


    /* 8 second */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask8);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}
#endif



#ifdef ALLOW_ODD_PACKSIZES
static void
unpack_09 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg;
    const __m128i mask9 =  _mm_set1_epi32(511U);

    OutReg = _mm_and_si128( InReg , mask9);
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,9) , mask9);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,18) , mask9);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,27) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 9 - 4), mask9));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask9);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,13) , mask9);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,22) , mask9);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,31) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 9 - 8), mask9));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask9);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,17) , mask9);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 9 - 3), mask9));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,3) , mask9);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask9);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,21) , mask9);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 9 - 7), mask9));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,7) , mask9);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}
#endif



#ifdef WORDS_BIGENDIAN
static void
unpack_10 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    *out = ( Bigendian_convert_uint(*in) >>  0  )   % (1U << 10 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  10  )   % (1U << 10 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  20  )   % (1U << 10 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  30  )   % (1U << 10 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 8 ))<<( 10 - 8 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  8  )   % (1U << 10 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  18  )   % (1U << 10 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  28  )   % (1U << 10 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 6 ))<<( 10 - 6 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  6  )   % (1U << 10 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  16  )   % (1U << 10 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  26  )   % (1U << 10 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 4 ))<<( 10 - 4 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  4  )   % (1U << 10 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  14  )   % (1U << 10 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  24  )   % (1U << 10 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 2 ))<<( 10 - 2 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  2  )   % (1U << 10 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  12  )   % (1U << 10 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  22  )   % (1U << 10 ) ;
    out++;

  }
  return;
}

#elif !defined(HAVE_SSE2)
static void
unpack_10 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    *out = ( (*in) >>  0  )   % (1U << 10 ) ;
    out++;
    *out = ( (*in) >>  10  )   % (1U << 10 ) ;
    out++;
    *out = ( (*in) >>  20  )   % (1U << 10 ) ;
    out++;
    *out = ( (*in) >>  30  )   % (1U << 10 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 8 ))<<( 10 - 8 );
    out++;
    *out = ( (*in) >>  8  )   % (1U << 10 ) ;
    out++;
    *out = ( (*in) >>  18  )   % (1U << 10 ) ;
    out++;
    *out = ( (*in) >>  28  )   % (1U << 10 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 6 ))<<( 10 - 6 );
    out++;
    *out = ( (*in) >>  6  )   % (1U << 10 ) ;
    out++;
    *out = ( (*in) >>  16  )   % (1U << 10 ) ;
    out++;
    *out = ( (*in) >>  26  )   % (1U << 10 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 4 ))<<( 10 - 4 );
    out++;
    *out = ( (*in) >>  4  )   % (1U << 10 ) ;
    out++;
    *out = ( (*in) >>  14  )   % (1U << 10 ) ;
    out++;
    *out = ( (*in) >>  24  )   % (1U << 10 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 2 ))<<( 10 - 2 );
    out++;
    *out = ( (*in) >>  2  )   % (1U << 10 ) ;
    out++;
    *out = ( (*in) >>  12  )   % (1U << 10 ) ;
    out++;
    *out = ( (*in) >>  22  )   % (1U << 10 ) ;
    out++;

  }
  return;
}

#else
static void
unpack_10_fwd (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg;
    const __m128i mask10 =  _mm_set1_epi32(1023U);

    OutReg = _mm_and_si128( InReg , mask10);
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,10) , mask10);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,20) , mask10);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 10 - 8), mask10));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask10);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,18) , mask10);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 10 - 6), mask10));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask10);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_10_fwd_1_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask10 =  _mm_set1_epi32(1023U);

    /* 1 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask10);
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 3 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,20) , mask10);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_10_fwd_2_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask10 =  _mm_set1_epi32(1023U);

    /* 2 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask10);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,10) , mask10);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);


    /* 4 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,20) , mask10);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 10 - 8), mask10));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_10_fwd_3_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask10 =  _mm_set1_epi32(1023U);

    /* 3 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,20) , mask10);
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 5 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask10);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_10_fwd_4_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask10 =  _mm_set1_epi32(1023U);

    /* 4 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,20) , mask10);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 10 - 8), mask10));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);


    /* 6 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask10);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,18) , mask10);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_10_fwd_5_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask10 =  _mm_set1_epi32(1023U);

    /* 5 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask10);
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 7 */
    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 10 - 6), mask10));
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_10_fwd_6_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask10 =  _mm_set1_epi32(1023U);

    /* 6 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask10);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,18) , mask10);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);


    /* 8 */
    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 10 - 6), mask10));
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask10);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_10_fwd_7_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask10 =  _mm_set1_epi32(1023U);

    /* 1 first */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask10);
    _mm_store_si128(&(out[2]), total);


    /* 7 second */
    InReg = _mm_load_si128(++in);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 10 - 6), mask10));
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_10_fwd_8_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask10 =  _mm_set1_epi32(1023U);

    /* 2 first */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask10);
    _mm_store_si128(&(out[2]), total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,10) , mask10);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(&(out[3]), total);


    /* 8 second */
    InReg = _mm_load_si128(++in);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 10 - 6), mask10));
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask10);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_10_rev (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg;
    const __m128i mask10 =  _mm_set1_epi32(1023U);

    in += 2;
    InReg = _mm_load_si128(in);


    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask10);
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 10 - 4), mask10));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask10);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,14) , mask10);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 10 - 2), mask10));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask10);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask10);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_10_rev_1_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask10 =  _mm_set1_epi32(1023U);

    /* 1 */
    in += 2;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask10);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 3 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask10);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_10_rev_2_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask10 =  _mm_set1_epi32(1023U);

    /* 2 */
    in += 2;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask10);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 10 - 4), mask10));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);


    /* 4 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask10);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,14) , mask10);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_10_rev_3_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask10 =  _mm_set1_epi32(1023U);

    /* 3 */
    in += 3;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask10);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 5 */
    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 10 - 2), mask10));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_10_rev_4_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask10 =  _mm_set1_epi32(1023U);

    /* 4 */
    in += 3;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask10);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,14) , mask10);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);


    /* 6 */
    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 10 - 2), mask10));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask10);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_10_rev_5_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask10 =  _mm_set1_epi32(1023U);

    /* 5 */
    in += 3;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 10 - 2), mask10));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 7 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask10);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_10_rev_6_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask10 =  _mm_set1_epi32(1023U);

    /* 6 */
    in += 3;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 10 - 2), mask10));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask10);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);


    /* 8 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask10);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,22) ;
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_10_rev_7_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask10 =  _mm_set1_epi32(1023U);

    /* 1 first */
    in += 2;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask10);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(&(out[2]), total);


    /* 7 second */
    in += 2;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask10);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_10_rev_8_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask10 =  _mm_set1_epi32(1023U);

    /* 2 first */
    in += 2;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask10);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(&(out[2]), total);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 10 - 4), mask10));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(&(out[3]), total);


    /* 8 second */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask10);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,22) ;
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}
#endif


#ifdef ALLOW_ODD_PACKSIZES
static void
unpack_11 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg;
    const __m128i mask11 =  _mm_set1_epi32(2047U);

    OutReg = _mm_and_si128( InReg , mask11);
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,11) , mask11);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 11 - 1), mask11));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,1) , mask11);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask11);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,23) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 11 - 2), mask11));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask11);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,13) , mask11);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 11 - 3), mask11));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,3) , mask11);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,14) , mask11);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,25) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 11 - 4), mask11));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask11);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,15) , mask11);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 11 - 5), mask11));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,5) , mask11);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}
#endif


#ifdef WORDS_BIGENDIAN
static void
unpack_12 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    *out = ( Bigendian_convert_uint(*in) >>  0  )   % (1U << 12 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  12  )   % (1U << 12 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  24  )   % (1U << 12 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 4 ))<<( 12 - 4 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  4  )   % (1U << 12 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  16  )   % (1U << 12 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  28  )   % (1U << 12 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 8 ))<<( 12 - 8 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  8  )   % (1U << 12 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  20  )   % (1U << 12 ) ;
    in += 4;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  0  )   % (1U << 12 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  12  )   % (1U << 12 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  24  )   % (1U << 12 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 4 ))<<( 12 - 4 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  4  )   % (1U << 12 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  16  )   % (1U << 12 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  28  )   % (1U << 12 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 8 ))<<( 12 - 8 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  8  )   % (1U << 12 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  20  )   % (1U << 12 ) ;
    out++;
  }

  return;
}

#elif !defined(HAVE_SSE2)
static void
unpack_12 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    *out = ( (*in) >>  0  )   % (1U << 12 ) ;
    out++;
    *out = ( (*in) >>  12  )   % (1U << 12 ) ;
    out++;
    *out = ( (*in) >>  24  )   % (1U << 12 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 4 ))<<( 12 - 4 );
    out++;
    *out = ( (*in) >>  4  )   % (1U << 12 ) ;
    out++;
    *out = ( (*in) >>  16  )   % (1U << 12 ) ;
    out++;
    *out = ( (*in) >>  28  )   % (1U << 12 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 8 ))<<( 12 - 8 );
    out++;
    *out = ( (*in) >>  8  )   % (1U << 12 ) ;
    out++;
    *out = ( (*in) >>  20  )   % (1U << 12 ) ;
    in += 4;
    out++;
    *out = ( (*in) >>  0  )   % (1U << 12 ) ;
    out++;
    *out = ( (*in) >>  12  )   % (1U << 12 ) ;
    out++;
    *out = ( (*in) >>  24  )   % (1U << 12 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 4 ))<<( 12 - 4 );
    out++;
    *out = ( (*in) >>  4  )   % (1U << 12 ) ;
    out++;
    *out = ( (*in) >>  16  )   % (1U << 12 ) ;
    out++;
    *out = ( (*in) >>  28  )   % (1U << 12 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 8 ))<<( 12 - 8 );
    out++;
    *out = ( (*in) >>  8  )   % (1U << 12 ) ;
    out++;
    *out = ( (*in) >>  20  )   % (1U << 12 ) ;
    out++;
  }

  return;
}

#else
static void
unpack_12_fwd (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg;
    const __m128i mask12 =  _mm_set1_epi32(4095U);

    OutReg = _mm_and_si128( InReg , mask12);
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask12);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 12 - 4), mask12));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask12);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask12);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 12 - 8), mask12));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask12);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_12_fwd_1_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask12 =  _mm_set1_epi32(4095U);

    /* 1 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask12);
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 3 */
    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 12 - 4), mask12));
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_12_fwd_2_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask12 =  _mm_set1_epi32(4095U);

    /* 2 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask12);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask12);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    /* 4 */
    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 12 - 4), mask12));
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask12);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_12_fwd_3_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask12 =  _mm_set1_epi32(4095U);

    /* 3 */
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 12 - 4), mask12));
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 5 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask12);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_12_fwd_4_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask12 =  _mm_set1_epi32(4095U);

    /* 4 */
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 12 - 4), mask12));
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask12);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);


    /* 6 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask12);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 12 - 8), mask12));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_12_fwd_5_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask12 =  _mm_set1_epi32(4095U);

    /* 5 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask12);
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 7 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask12);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_12_fwd_6_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask12 =  _mm_set1_epi32(4095U);

    /* 6 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask12);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 12 - 8), mask12));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);


    /* 8 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask12);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_12_fwd_7_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask12 =  _mm_set1_epi32(4095U);

    /* 1 first */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask12);
    _mm_store_si128(&(out[2]), total);


    /* 7 second */
    in += 2;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask12);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_12_fwd_8_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask12 =  _mm_set1_epi32(4095U);

    /* 2 first */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask12);
    _mm_store_si128(&(out[2]), total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask12);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(&(out[3]), total);


    /* 8 second */
    in += 2;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask12);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_12_rev (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg;
    const __m128i mask12 =  _mm_set1_epi32(4095U);

    in += 3;
    InReg = _mm_load_si128(in);


    OutReg = _mm_and_si128( InReg , mask12);
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask12);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 12 - 4), mask12));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask12);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask12);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 12 - 8), mask12));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask12);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_12_rev_1_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask12 =  _mm_set1_epi32(4095U);

    /* 1 */
    in += 3;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask12);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 3 */
    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 12 - 4), mask12));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_12_rev_2_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask12 =  _mm_set1_epi32(4095U);

    /* 2 */
    in += 3;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask12);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask12);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);


    /* 4 */
    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 12 - 4), mask12));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask12);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_12_rev_3_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask12 =  _mm_set1_epi32(4095U);

    /* 3 */
    in += 3;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 12 - 4), mask12));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 5 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask12);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_12_rev_4_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask12 =  _mm_set1_epi32(4095U);

    /* 4 */
    in += 3;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 12 - 4), mask12));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask12);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);


    /* 6 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask12);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 12 - 8), mask12));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_12_rev_5_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask12 =  _mm_set1_epi32(4095U);

    /* 5 */
    in += 4;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask12);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 7 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask12);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_12_rev_6_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask12 =  _mm_set1_epi32(4095U);

    /* 6 */
    in += 4;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask12);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 12 - 8), mask12));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);


    /* 8 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask12);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,20) ;
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_12_rev_7_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask12 =  _mm_set1_epi32(4095U);

    /* 1 first */
    in += 3;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask12);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(&(out[2]), total);


    /* 7 second */
    in += 2;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask12);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_12_rev_8_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask12 =  _mm_set1_epi32(4095U);

    /* 2 first */
    in += 3;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask12);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(&(out[2]), total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask12);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(&(out[3]), total);


    /* 8 second */
    in += 2;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask12);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,20) ;
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}
#endif


#ifdef ALLOW_ODD_PACKSIZES
static void
unpack_13 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg;
    const __m128i mask13 =  _mm_set1_epi32(8191U);

    OutReg = _mm_and_si128( InReg , mask13);
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,13) , mask13);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 13 - 7), mask13));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,7) , mask13);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 13 - 1), mask13));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,1) , mask13);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,14) , mask13);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,27) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 13 - 8), mask13));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask13);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,21) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 13 - 2), mask13));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask13);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,15) , mask13);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 13 - 9), mask13));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,9) , mask13);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 13 - 3), mask13));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,3) , mask13);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}
#endif


#ifdef WORDS_BIGENDIAN
static void
unpack_14 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    *out = ( Bigendian_convert_uint(*in) >>  0  )   % (1U << 14 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  14  )   % (1U << 14 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  28  )   % (1U << 14 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 10 ))<<( 14 - 10 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  10  )   % (1U << 14 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  24  )   % (1U << 14 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 6 ))<<( 14 - 6 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  6  )   % (1U << 14 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  20  )   % (1U << 14 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 2 ))<<( 14 - 2 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  2  )   % (1U << 14 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  16  )   % (1U << 14 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  30  )   % (1U << 14 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 12 ))<<( 14 - 12 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  12  )   % (1U << 14 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  26  )   % (1U << 14 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 8 ))<<( 14 - 8 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  8  )   % (1U << 14 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  22  )   % (1U << 14 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 4 ))<<( 14 - 4 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  4  )   % (1U << 14 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  18  )   % (1U << 14 ) ;
    out++;
  }

  return;
}

#elif !defined(HAVE_SSE2)
static void
unpack_14 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    *out = ( (*in) >>  0  )   % (1U << 14 ) ;
    out++;
    *out = ( (*in) >>  14  )   % (1U << 14 ) ;
    out++;
    *out = ( (*in) >>  28  )   % (1U << 14 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 10 ))<<( 14 - 10 );
    out++;
    *out = ( (*in) >>  10  )   % (1U << 14 ) ;
    out++;
    *out = ( (*in) >>  24  )   % (1U << 14 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 6 ))<<( 14 - 6 );
    out++;
    *out = ( (*in) >>  6  )   % (1U << 14 ) ;
    out++;
    *out = ( (*in) >>  20  )   % (1U << 14 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 2 ))<<( 14 - 2 );
    out++;
    *out = ( (*in) >>  2  )   % (1U << 14 ) ;
    out++;
    *out = ( (*in) >>  16  )   % (1U << 14 ) ;
    out++;
    *out = ( (*in) >>  30  )   % (1U << 14 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 12 ))<<( 14 - 12 );
    out++;
    *out = ( (*in) >>  12  )   % (1U << 14 ) ;
    out++;
    *out = ( (*in) >>  26  )   % (1U << 14 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 8 ))<<( 14 - 8 );
    out++;
    *out = ( (*in) >>  8  )   % (1U << 14 ) ;
    out++;
    *out = ( (*in) >>  22  )   % (1U << 14 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 4 ))<<( 14 - 4 );
    out++;
    *out = ( (*in) >>  4  )   % (1U << 14 ) ;
    out++;
    *out = ( (*in) >>  18  )   % (1U << 14 ) ;
    out++;
  }

  return;
}

#else
static void
unpack_14_fwd (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg;
    const __m128i mask14 =  _mm_set1_epi32(16383U);

    OutReg = _mm_and_si128( InReg , mask14);
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,14) , mask14);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 14 - 10), mask14));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,10) , mask14);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 14 - 6), mask14));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask14);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 14 - 2), mask14));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask14);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_14_fwd_1_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask14 =  _mm_set1_epi32(16383U);

    /* 1 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask14);
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 3 */
    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 14 - 10), mask14));
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_14_fwd_2_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask14 =  _mm_set1_epi32(16383U);

    /* 2 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask14);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,14) , mask14);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);


    /* 4 */
    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 14 - 10), mask14));
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,10) , mask14);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_14_fwd_3_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask14 =  _mm_set1_epi32(16383U);

    /* 3 */
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 14 - 10), mask14));
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 5 */
    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 14 - 6), mask14));
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_14_fwd_4_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask14 =  _mm_set1_epi32(16383U);

    /* 4 */
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 14 - 10), mask14));
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,10) , mask14);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);


    /* 6 */
    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 14 - 6), mask14));
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask14);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_14_fwd_5_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask14 =  _mm_set1_epi32(16383U);

    /* 5 */
    InReg = _mm_load_si128(++in);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 14 - 6), mask14));
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 7 */
    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 14 - 2), mask14));
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_14_fwd_6_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask14 =  _mm_set1_epi32(16383U);

    /* 6 */
    InReg = _mm_load_si128(++in);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 14 - 6), mask14));
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask14);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);


    /* 8 */
    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 14 - 2), mask14));
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask14);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_14_fwd_7_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask14 =  _mm_set1_epi32(16383U);

    /* 1 first */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask14);
    _mm_store_si128(&(out[2]), total);


    /* 7 second */
    in += 2;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 14 - 2), mask14));
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_14_fwd_8_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask14 =  _mm_set1_epi32(16383U);

    /* 2 first */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask14);
    _mm_store_si128(&(out[2]), total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,14) , mask14);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(&(out[3]), total);


    /* 8 second */
    in += 2;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 14 - 2), mask14));
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask14);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_14_rev (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg;
    const __m128i mask14 =  _mm_set1_epi32(16383U);

    in += 3;
    InReg = _mm_load_si128(in);


    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask14);
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 14 - 12), mask14));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask14);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 14 - 8), mask14));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask14);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 14 - 4), mask14));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask14);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,18) ;
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_14_rev_1_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask14 =  _mm_set1_epi32(16383U);

    /* 1 */
    in += 3;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask14);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 3 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask14);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_14_rev_2_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask14 =  _mm_set1_epi32(16383U);

    /* 2 */
    in += 3;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask14);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 14 - 12), mask14));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);


    /* 4 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask14);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 14 - 8), mask14));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_14_rev_3_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask14 =  _mm_set1_epi32(16383U);

    /* 3 */
    in += 4;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask14);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 5 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask14);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_14_rev_4_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask14 =  _mm_set1_epi32(16383U);

    /* 4 */
    in += 4;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask14);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 14 - 8), mask14));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);


    /* 6 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask14);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 14 - 4), mask14));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_14_rev_5_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask14 =  _mm_set1_epi32(16383U);

    /* 5 */
    in += 5;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask14);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 7 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask14);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_14_rev_6_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask14 =  _mm_set1_epi32(16383U);

    /* 6 */
    in += 5;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask14);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 14 - 4), mask14));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);


    /* 8 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask14);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,18) ;
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_14_rev_7_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask14 =  _mm_set1_epi32(16383U);

    /* 1 first */
    in += 3;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask14);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(&(out[2]), total);


    /* 7 second */
    in += 3;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask14);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_14_rev_8_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask14 =  _mm_set1_epi32(16383U);

    /* 2 first */
    in += 3;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask14);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(&(out[2]), total);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 14 - 12), mask14));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(&(out[3]), total);


    /* 8 second */
    in += 2;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask14);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,18) ;
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}
#endif



#ifdef ALLOW_ODD_PACKSIZES
static void
unpack_15 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg;
    const __m128i mask15 =  _mm_set1_epi32(32767U);

    OutReg = _mm_and_si128( InReg , mask15);
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,15) , mask15);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 15 - 13), mask15));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,13) , mask15);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 15 - 11), mask15));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,11) , mask15);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 15 - 9), mask15));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,9) , mask15);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 15 - 7), mask15));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,7) , mask15);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 15 - 5), mask15));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,5) , mask15);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 15 - 3), mask15));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,3) , mask15);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,18) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 15 - 1), mask15));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,1) , mask15);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}
#endif



#ifdef WORDS_BIGENDIAN
static void
unpack_16 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  UINT4 outer, inwordpointer;
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    for (outer = 0; outer < 8; outer++) {
      for(inwordpointer =  0; inwordpointer <32; inwordpointer += 16) {
	*(out++) = ( Bigendian_convert_uint(*in) >> inwordpointer )   % (1U << 16 ) ;
      }
      in += 4;
    }
  }

  return;
}

#elif !defined(HAVE_SSE2)
static void
unpack_16 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  UINT4 outer, inwordpointer;
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    for (outer = 0; outer < 8; outer++) {
      for(inwordpointer =  0; inwordpointer <32; inwordpointer += 16) {
	*(out++) = ( (*in) >> inwordpointer )   % (1U << 16 ) ;
      }
      in += 4;
    }
  }

  return;
}

#else
static void
unpack_16_fwd (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg;
    const __m128i mask16 =  _mm_set1_epi32(65535U);

    OutReg = _mm_and_si128( InReg , mask16);
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128( InReg , mask16);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128( InReg , mask16);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128( InReg , mask16);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_16_fwd_1_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask16 =  _mm_set1_epi32(65535U);

    /* 1 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 3 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_16_fwd_2_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask16 =  _mm_set1_epi32(65535U);

    /* 2 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    /* InReg = _mm_load_si128(++in); */

    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);


    /* 4 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    /* InReg = _mm_load_si128(++in); */

    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_16_fwd_3_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask16 =  _mm_set1_epi32(65535U);

    /* 3 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 5 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_16_fwd_4_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask16 =  _mm_set1_epi32(65535U);

    /* 4 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    /* InReg = _mm_load_si128(++in); */

    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);


    /* 6 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    /* InReg = _mm_load_si128(++in); */

    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_16_fwd_5_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask16 =  _mm_set1_epi32(65535U);

    /* 5 */
    in += 2;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 7 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_16_fwd_6_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask16 =  _mm_set1_epi32(65535U);

    /* 6 */
    in += 2;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    /* InReg = _mm_load_si128(++in); */

    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);


    /* 8 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_16_fwd_7_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask16 =  _mm_set1_epi32(65535U);

    /* 1 first */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
    _mm_store_si128(&(out[2]), total);


    /* 7 second */
    in += 3;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_16_fwd_8_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask16 =  _mm_set1_epi32(65535U);

    /* 2 first */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
    _mm_store_si128(&(out[2]), total);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    /* InReg = _mm_load_si128(++in); */

    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(&(out[3]), total);


    /* 8 second */
    in += 3;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_16_rev (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg;
    const __m128i mask16 =  _mm_set1_epi32(65535U);

    in += 4;
    InReg = _mm_load_si128(in);


    OutReg = _mm_and_si128( InReg , mask16);
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128( InReg , mask16);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128( InReg , mask16);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128( InReg , mask16);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_16_rev_1_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask16 =  _mm_set1_epi32(65535U);

    /* 1 */
    in += 4;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 3 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_16_rev_2_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask16 =  _mm_set1_epi32(65535U);

    /* 2 */
    in += 4;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    /* InReg = _mm_load_si128(++in); */

#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);


    /* 4 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    /* InReg = _mm_load_si128(++in); */

#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_16_rev_3_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask16 =  _mm_set1_epi32(65535U);

    /* 3 */
    in += 5;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 5 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_16_rev_4_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask16 =  _mm_set1_epi32(65535U);

    /* 4 */
    in += 5;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    /* InReg = _mm_load_si128(++in); */

#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);


    /* 6 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    /* InReg = _mm_load_si128(++in); */

#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_16_rev_5_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask16 =  _mm_set1_epi32(65535U);

    /* 5 */
    in += 6;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 7 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_16_rev_6_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask16 =  _mm_set1_epi32(65535U);

    /* 6 */
    in += 6;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    /* InReg = _mm_load_si128(++in); */

#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);


    /* 8 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,16) ;
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_16_rev_7_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask16 =  _mm_set1_epi32(65535U);

    /* 1 first */
    in += 4;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(&(out[2]), total);


    /* 7 second */
    in += 3;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_16_rev_8_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask16 =  _mm_set1_epi32(65535U);

    /* 2 first */
    in += 4;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(&(out[2]), total);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    /* InReg = _mm_load_si128(++in); */

#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(&(out[3]), total);


    /* 8 second */
    in += 3;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,16) ;
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}
#endif


#ifdef ALLOW_ODD_PACKSIZES
static void
unpack_17 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg;
    const __m128i mask17 =  _mm_set1_epi32(131071U);

    OutReg = _mm_and_si128( InReg , mask17);
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,17) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 17 - 2), mask17));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask17);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,19) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 17 - 4), mask17));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask17);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,21) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 17 - 6), mask17));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask17);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,23) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 17 - 8), mask17));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask17);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,25) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 17 - 10), mask17));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,10) , mask17);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,27) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 17 - 12), mask17));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask17);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,29) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 17 - 14), mask17));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,14) , mask17);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,31) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 17 - 16), mask17));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}
#endif



#ifdef WORDS_BIGENDIAN
static void
unpack_18 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    *out = ( Bigendian_convert_uint(*in) >>  0  )   % (1U << 18 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  18  )   % (1U << 18 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 4 ))<<( 18 - 4 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  4  )   % (1U << 18 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  22  )   % (1U << 18 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 8 ))<<( 18 - 8 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  8  )   % (1U << 18 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  26  )   % (1U << 18 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 12 ))<<( 18 - 12 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  12  )   % (1U << 18 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  30  )   % (1U << 18 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 16 ))<<( 18 - 16 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  16  )   % (1U << 18 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 2 ))<<( 18 - 2 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  2  )   % (1U << 18 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  20  )   % (1U << 18 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 6 ))<<( 18 - 6 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  6  )   % (1U << 18 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  24  )   % (1U << 18 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 10 ))<<( 18 - 10 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  10  )   % (1U << 18 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  28  )   % (1U << 18 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 14 ))<<( 18 - 14 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  14  )   % (1U << 18 ) ;
    out++;
  }

  return;
}

#elif !defined(HAVE_SSE2)
static void
unpack_18 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    *out = ( (*in) >>  0  )   % (1U << 18 ) ;
    out++;
    *out = ( (*in) >>  18  )   % (1U << 18 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 4 ))<<( 18 - 4 );
    out++;
    *out = ( (*in) >>  4  )   % (1U << 18 ) ;
    out++;
    *out = ( (*in) >>  22  )   % (1U << 18 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 8 ))<<( 18 - 8 );
    out++;
    *out = ( (*in) >>  8  )   % (1U << 18 ) ;
    out++;
    *out = ( (*in) >>  26  )   % (1U << 18 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 12 ))<<( 18 - 12 );
    out++;
    *out = ( (*in) >>  12  )   % (1U << 18 ) ;
    out++;
    *out = ( (*in) >>  30  )   % (1U << 18 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 16 ))<<( 18 - 16 );
    out++;
    *out = ( (*in) >>  16  )   % (1U << 18 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 2 ))<<( 18 - 2 );
    out++;
    *out = ( (*in) >>  2  )   % (1U << 18 ) ;
    out++;
    *out = ( (*in) >>  20  )   % (1U << 18 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 6 ))<<( 18 - 6 );
    out++;
    *out = ( (*in) >>  6  )   % (1U << 18 ) ;
    out++;
    *out = ( (*in) >>  24  )   % (1U << 18 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 10 ))<<( 18 - 10 );
    out++;
    *out = ( (*in) >>  10  )   % (1U << 18 ) ;
    out++;
    *out = ( (*in) >>  28  )   % (1U << 18 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 14 ))<<( 18 - 14 );
    out++;
    *out = ( (*in) >>  14  )   % (1U << 18 ) ;
    out++;
  }

  return;
}

#else
static void
unpack_18_fwd (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg;
    const __m128i mask18 =  _mm_set1_epi32(262143U);

    OutReg = _mm_and_si128( InReg , mask18);
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,18) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 4), mask18));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask18);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 8), mask18));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask18);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 12), mask18));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask18);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 16), mask18));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_18_fwd_1_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask18 =  _mm_set1_epi32(262143U);

    /* 1 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask18);
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 3 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask18);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_18_fwd_2_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask18 =  _mm_set1_epi32(262143U);

    /* 2 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask18);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,18) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 4), mask18));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);


    /* 4 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask18);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 8), mask18));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_18_fwd_3_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask18 =  _mm_set1_epi32(262143U);

    /* 3 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask18);
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 5 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask18);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_18_fwd_4_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask18 =  _mm_set1_epi32(262143U);

    /* 4 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask18);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 8), mask18));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);


    /* 6 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask18);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 12), mask18));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_18_fwd_5_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask18 =  _mm_set1_epi32(262143U);

    /* 5 */
    in += 2;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask18);
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 7 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask18);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_18_fwd_6_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask18 =  _mm_set1_epi32(262143U);

    /* 6 */
    in += 2;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask18);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 12), mask18));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);


    /* 8 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask18);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 16), mask18));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_18_fwd_7_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask18 =  _mm_set1_epi32(262143U);

    /* 1 first */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask18);
    _mm_store_si128(&(out[2]), total);


    /* 7 second */
    in += 3;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask18);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_18_fwd_8_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask18 =  _mm_set1_epi32(262143U);

    /* 2 first */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask18);
    _mm_store_si128(&(out[2]), total);

    OutReg =   _mm_srli_epi32(InReg,18) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 4), mask18));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(&(out[3]), total);



    /* 8 second */
    in += 2;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask18);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 16), mask18));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_18_rev (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg;
    const __m128i mask18 =  _mm_set1_epi32(262143U);

    in += 4;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 2), mask18));
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask18);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 6), mask18));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask18);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 10), mask18));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,10) , mask18);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 14), mask18));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,14) ;
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_18_rev_1_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask18 =  _mm_set1_epi32(262143U);

    /* 1 */
    in += 4;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 2), mask18));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 3 */
    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 6), mask18));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_18_rev_2_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask18 =  _mm_set1_epi32(262143U);

    /* 2 */
    in += 4;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 2), mask18));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask18);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);


    /* 4 */
    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 6), mask18));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask18);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);


    return;
}

static void
unpack_18_rev_3_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask18 =  _mm_set1_epi32(262143U);

    /* 3 */
    in += 5;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 6), mask18));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 5 */
    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 10), mask18));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_18_rev_4_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask18 =  _mm_set1_epi32(262143U);

    /* 4 */
    in += 5;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 6), mask18));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask18);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);


    /* 6 */
    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 10), mask18));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,10) , mask18);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_18_rev_5_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask18 =  _mm_set1_epi32(262143U);

    /* 5 */
    in += 6;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 10), mask18));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 7 */
    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 14), mask18));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_18_rev_6_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask18 =  _mm_set1_epi32(262143U);

    /* 6 */
    in += 6;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 10), mask18));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,10) , mask18);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);


    /* 8 */
    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 14), mask18));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,14) ;
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_18_rev_7_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask18 =  _mm_set1_epi32(262143U);

    /* 1 first */
    in += 4;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 2), mask18));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(&(out[2]), total);


    /* 7 second */
    in += (7 - 5);
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 14), mask18));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_18_rev_8_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask18 =  _mm_set1_epi32(262143U);

    /* 2 first */
    in += 4;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 2), mask18));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(&(out[2]), total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask18);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(&(out[3]), total);


    /* 8 second */
    in += (7 - 5);
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 14), mask18));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,14) ;
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}
#endif


#ifdef ALLOW_ODD_PACKSIZES
static void
unpack_19 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg;
    const __m128i mask19 =  _mm_set1_epi32(524287U);

    OutReg = _mm_and_si128( InReg , mask19);
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,19) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 19 - 6), mask19));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask19);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,25) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 19 - 12), mask19));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask19);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,31) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 19 - 18), mask19));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,18) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 19 - 5), mask19));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,5) , mask19);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 19 - 11), mask19));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,11) , mask19);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 19 - 17), mask19));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,17) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 19 - 4), mask19));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask19);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,23) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 19 - 10), mask19));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,10) , mask19);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,29) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 19 - 16), mask19));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}
#endif



#ifdef WORDS_BIGENDIAN
static void
unpack_20 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    *out = ( Bigendian_convert_uint(*in) >>  0  )   % (1U << 20 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  20  )   % (1U << 20 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 8 ))<<( 20 - 8 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  8  )   % (1U << 20 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  28  )   % (1U << 20 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 16 ))<<( 20 - 16 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  16  )   % (1U << 20 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 4 ))<<( 20 - 4 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  4  )   % (1U << 20 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  24  )   % (1U << 20 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 12 ))<<( 20 - 12 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  12  )   % (1U << 20 ) ;
    in += 4;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  0  )   % (1U << 20 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  20  )   % (1U << 20 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 8 ))<<( 20 - 8 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  8  )   % (1U << 20 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  28  )   % (1U << 20 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 16 ))<<( 20 - 16 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  16  )   % (1U << 20 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 4 ))<<( 20 - 4 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  4  )   % (1U << 20 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  24  )   % (1U << 20 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 12 ))<<( 20 - 12 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  12  )   % (1U << 20 ) ;
    out++;
  }

  return;
}

#elif !defined(HAVE_SSE2)
static void
unpack_20 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    *out = ( (*in) >>  0  )   % (1U << 20 ) ;
    out++;
    *out = ( (*in) >>  20  )   % (1U << 20 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 8 ))<<( 20 - 8 );
    out++;
    *out = ( (*in) >>  8  )   % (1U << 20 ) ;
    out++;
    *out = ( (*in) >>  28  )   % (1U << 20 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 16 ))<<( 20 - 16 );
    out++;
    *out = ( (*in) >>  16  )   % (1U << 20 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 4 ))<<( 20 - 4 );
    out++;
    *out = ( (*in) >>  4  )   % (1U << 20 ) ;
    out++;
    *out = ( (*in) >>  24  )   % (1U << 20 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 12 ))<<( 20 - 12 );
    out++;
    *out = ( (*in) >>  12  )   % (1U << 20 ) ;
    in += 4;
    out++;
    *out = ( (*in) >>  0  )   % (1U << 20 ) ;
    out++;
    *out = ( (*in) >>  20  )   % (1U << 20 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 8 ))<<( 20 - 8 );
    out++;
    *out = ( (*in) >>  8  )   % (1U << 20 ) ;
    out++;
    *out = ( (*in) >>  28  )   % (1U << 20 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 16 ))<<( 20 - 16 );
    out++;
    *out = ( (*in) >>  16  )   % (1U << 20 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 4 ))<<( 20 - 4 );
    out++;
    *out = ( (*in) >>  4  )   % (1U << 20 ) ;
    out++;
    *out = ( (*in) >>  24  )   % (1U << 20 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 12 ))<<( 20 - 12 );
    out++;
    *out = ( (*in) >>  12  )   % (1U << 20 ) ;
    out++;
  }

  return;
}

#else
static void
unpack_20_fwd (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg;
    const __m128i mask20 =  _mm_set1_epi32(1048575U);

    OutReg = _mm_and_si128( InReg , mask20);
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 8), mask20));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask20);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 16), mask20));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 4), mask20));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask20);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 12), mask20));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_20_fwd_1_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask20 =  _mm_set1_epi32(1048575U);

    /* 1 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask20);
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 3 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask20);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_20_fwd_2_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask20 =  _mm_set1_epi32(1048575U);

    /* 2 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask20);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 8), mask20));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);


    /* 4 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask20);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 16), mask20));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_20_fwd_3_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask20 =  _mm_set1_epi32(1048575U);

    /* 3 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask20);
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 5 */
    InReg = _mm_load_si128(++in);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 4), mask20));
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_20_fwd_4_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask20 =  _mm_set1_epi32(1048575U);

    /* 4 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask20);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 16), mask20));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);


    /* 6 */
    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 4), mask20));
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask20);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);


    return;

}

static void
unpack_20_fwd_5_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask20 =  _mm_set1_epi32(1048575U);

    /* 5 */
    in += 2;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 4), mask20));
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 7 */
    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 12), mask20));
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_20_fwd_6_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask20 =  _mm_set1_epi32(1048575U);

    /* 6 */
    in += 2;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 4), mask20));
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask20);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);


    /* 8 */
    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 12), mask20));
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_20_fwd_7_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask20 =  _mm_set1_epi32(1048575U);

    /* 1 first */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask20);
    _mm_store_si128(&(out[2]), total);


    /* 7 second */
    in += 3;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 12), mask20));
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_20_fwd_8_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask20 =  _mm_set1_epi32(1048575U);

    /* 2 first */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask20);
    _mm_store_si128(&(out[2]), total);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 8), mask20));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(&(out[3]), total);


    /* 8 second */
    in += (3 - 1);
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 12), mask20));
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_20_rev (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg;
    const __m128i mask20 =  _mm_set1_epi32(1048575U);

    in += 5;
    InReg = _mm_load_si128(in);

    OutReg = _mm_and_si128( InReg , mask20);
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 8), mask20));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask20);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 16), mask20));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 4), mask20));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask20);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 12), mask20));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_20_rev_1_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask20 =  _mm_set1_epi32(1048575U);

    /* 1 */
    in += 5;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask20);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 3 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask20);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_20_rev_2_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask20 =  _mm_set1_epi32(1048575U);

    /* 2 */
    in += 5;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask20);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 8), mask20));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);


    /* 4 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask20);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 16), mask20));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
} 

static void
unpack_20_rev_3_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask20 =  _mm_set1_epi32(1048575U);

    /* 3 */
    in += 6;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask20);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 5 */
    InReg = _mm_load_si128(++in);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 4), mask20));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_20_rev_4_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask20 =  _mm_set1_epi32(1048575U);

    /* 4 */
    in += 6;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask20);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 16), mask20));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);


    /* 6 */
    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 4), mask20));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask20);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;

} 

static void
unpack_20_rev_5_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask20 =  _mm_set1_epi32(1048575U);

    /* 5 */
    in += 7;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 4), mask20));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 7 */
    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 12), mask20));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_20_rev_6_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask20 =  _mm_set1_epi32(1048575U);

    /* 6 */
    in += 7;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 4), mask20));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask20);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);


    /* 8 */
    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 12), mask20));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,12) ;
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_20_rev_7_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask20 =  _mm_set1_epi32(1048575U);

    /* 1 first */
    in += 5;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask20);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(&(out[2]), total);


    /* 7 second */
    in += (8 - 5);
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 12), mask20));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_20_rev_8_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask20 =  _mm_set1_epi32(1048575U);

    /* 2 first */
    in += 5;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask20);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(&(out[2]), total);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 8), mask20));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(&(out[3]), total);


    /* 8 second */
    in += (8 - 6);
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 12), mask20));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,12) ;
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}
#endif


#ifdef ALLOW_ODD_PACKSIZES
static void
unpack_21 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg;
    const __m128i mask21 =  _mm_set1_epi32(2097151U);

    OutReg = _mm_and_si128( InReg , mask21);
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,21) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 21 - 10), mask21));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,10) , mask21);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,31) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 21 - 20), mask21));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 21 - 9), mask21));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,9) , mask21);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 21 - 19), mask21));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,19) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 21 - 8), mask21));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask21);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,29) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 21 - 18), mask21));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,18) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 21 - 7), mask21));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,7) , mask21);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 21 - 17), mask21));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,17) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 21 - 6), mask21));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask21);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,27) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 21 - 16), mask21));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}
#endif


#ifdef WORDS_BIGENDIAN
static void
unpack_22 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    *out = ( Bigendian_convert_uint(*in) >>  0  )   % (1U << 22 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  22  )   % (1U << 22 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 12 ))<<( 22 - 12 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  12  )   % (1U << 22 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 2 ))<<( 22 - 2 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  2  )   % (1U << 22 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  24  )   % (1U << 22 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 14 ))<<( 22 - 14 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  14  )   % (1U << 22 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 4 ))<<( 22 - 4 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  4  )   % (1U << 22 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  26  )   % (1U << 22 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 16 ))<<( 22 - 16 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  16  )   % (1U << 22 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 6 ))<<( 22 - 6 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  6  )   % (1U << 22 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  28  )   % (1U << 22 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 18 ))<<( 22 - 18 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  18  )   % (1U << 22 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 8 ))<<( 22 - 8 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  8  )   % (1U << 22 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  30  )   % (1U << 22 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 20 ))<<( 22 - 20 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  20  )   % (1U << 22 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 10 ))<<( 22 - 10 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  10  )   % (1U << 22 ) ;
    out++;
  }

  return;
}

#elif !defined(HAVE_SSE2)
static void
unpack_22 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    *out = ( (*in) >>  0  )   % (1U << 22 ) ;
    out++;
    *out = ( (*in) >>  22  )   % (1U << 22 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 12 ))<<( 22 - 12 );
    out++;
    *out = ( (*in) >>  12  )   % (1U << 22 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 2 ))<<( 22 - 2 );
    out++;
    *out = ( (*in) >>  2  )   % (1U << 22 ) ;
    out++;
    *out = ( (*in) >>  24  )   % (1U << 22 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 14 ))<<( 22 - 14 );
    out++;
    *out = ( (*in) >>  14  )   % (1U << 22 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 4 ))<<( 22 - 4 );
    out++;
    *out = ( (*in) >>  4  )   % (1U << 22 ) ;
    out++;
    *out = ( (*in) >>  26  )   % (1U << 22 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 16 ))<<( 22 - 16 );
    out++;
    *out = ( (*in) >>  16  )   % (1U << 22 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 6 ))<<( 22 - 6 );
    out++;
    *out = ( (*in) >>  6  )   % (1U << 22 ) ;
    out++;
    *out = ( (*in) >>  28  )   % (1U << 22 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 18 ))<<( 22 - 18 );
    out++;
    *out = ( (*in) >>  18  )   % (1U << 22 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 8 ))<<( 22 - 8 );
    out++;
    *out = ( (*in) >>  8  )   % (1U << 22 ) ;
    out++;
    *out = ( (*in) >>  30  )   % (1U << 22 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 20 ))<<( 22 - 20 );
    out++;
    *out = ( (*in) >>  20  )   % (1U << 22 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 10 ))<<( 22 - 10 );
    out++;
    *out = ( (*in) >>  10  )   % (1U << 22 ) ;
    out++;
  }

  return;
}

#else
static void
unpack_22_fwd (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg;
    const __m128i mask22 =  _mm_set1_epi32(4194303U);

    OutReg = _mm_and_si128( InReg , mask22);
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 12), mask22));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 2), mask22));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask22);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 14), mask22));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,14) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 4), mask22));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask22);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 16), mask22));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_22_fwd_1_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask22 =  _mm_set1_epi32(4194303U);

    /* 1 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask22);
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 3 */
    InReg = _mm_load_si128(++in);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 2), mask22));
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_22_fwd_2_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask22 =  _mm_set1_epi32(4194303U);

    /* 2 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask22);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 12), mask22));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);


    /* 4 */
    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 2), mask22));
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask22);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_22_fwd_3_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask22 =  _mm_set1_epi32(4194303U);

    /* 3 */
    InReg = _mm_load_si128(++in);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 2), mask22));
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 5 */
    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 14), mask22));
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_22_fwd_4_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask22 =  _mm_set1_epi32(4194303U);

    /* 4 */
    InReg = _mm_load_si128(++in);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 2), mask22));
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask22);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);


    /* 6 */
    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 14), mask22));
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,14) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 4), mask22));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_22_fwd_5_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask22 =  _mm_set1_epi32(4194303U);

    /* 5 */
    in += 2;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 14), mask22));
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 7 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask22);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_22_fwd_6_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask22 =  _mm_set1_epi32(4194303U);

    /* 6 */
    in += 2;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 14), mask22));
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,14) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 4), mask22));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);


    /* 8 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask22);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 16), mask22));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_22_fwd_7_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask22 =  _mm_set1_epi32(4194303U);

    /* 1 first */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask22);
    _mm_store_si128(&(out[2]), total);


    /* 7 second */
    in += (4 - 0);
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask22);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_22_fwd_8_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask22 =  _mm_set1_epi32(4194303U);

    /* 2 first */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask22);
    _mm_store_si128(&(out[2]), total);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 12), mask22));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(&(out[3]), total);


    /* 8 second */
    in += (4 - 1);
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask22);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 16), mask22));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_22_rev (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg;
    const __m128i mask22 =  _mm_set1_epi32(4194303U);

    in += 5;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 6), mask22));
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask22);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 18), mask22));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,18) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 8), mask22));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask22);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 20), mask22));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 10), mask22));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,10) ;
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_22_rev_1_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask22 =  _mm_set1_epi32(4194303U);

    /* 1 */
    in += 5;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 6), mask22));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 3 */
    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 18), mask22));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_22_rev_2_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask22 =  _mm_set1_epi32(4194303U);

    /* 2 */
    in += 5;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 6), mask22));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask22);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);


    /* 4 */
    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 18), mask22));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,18) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 8), mask22));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);


    return;
}

static void
unpack_22_rev_3_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask22 =  _mm_set1_epi32(4194303U);

    /* 3 */
    in += 6;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 18), mask22));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 5 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask22);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_22_rev_4_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask22 =  _mm_set1_epi32(4194303U);

    /* 4 */
    in += 6;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 18), mask22));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,18) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 8), mask22));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);


    /* 6 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask22);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 20), mask22));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_22_rev_5_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask22 =  _mm_set1_epi32(4194303U);

    /* 5 */
    in += 8;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask22);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 7 */
    InReg = _mm_load_si128(++in);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 10), mask22));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_22_rev_6_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask22 =  _mm_set1_epi32(4194303U);

    /* 6 */
    in += 8;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask22);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 20), mask22));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);


    /* 8 */
    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 10), mask22));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,10) ;
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_22_rev_7_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask22 =  _mm_set1_epi32(4194303U);

    /* 1 first */
    in += 5;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 6), mask22));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(&(out[2]), total);


    /* 7 second */
    in += (9 - 6);
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 10), mask22));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_22_rev_8_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask22 =  _mm_set1_epi32(4194303U);

    /* 2 first */
    in += 5;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 6), mask22));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(&(out[2]), total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask22);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(&(out[3]), total);


    /* 8 second */
    in += (9 - 6);
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 10), mask22));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,10) ;
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}
#endif


#ifdef ALLOW_ODD_PACKSIZES
static void
unpack_23 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg;
    const __m128i mask23 =  _mm_set1_epi32(8388607U);

    OutReg = _mm_and_si128( InReg , mask23);
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,23) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 23 - 14), mask23));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,14) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 23 - 5), mask23));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,5) , mask23);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 23 - 19), mask23));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,19) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 23 - 10), mask23));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,10) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 23 - 1), mask23));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,1) , mask23);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 23 - 15), mask23));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,15) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 23 - 6), mask23));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask23);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,29) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 23 - 20), mask23));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 23 - 11), mask23));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,11) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 23 - 2), mask23));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask23);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,25) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 23 - 16), mask23));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}
#endif



#ifdef WORDS_BIGENDIAN
static void
unpack_24 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    *out = ( Bigendian_convert_uint(*in) >>  0  )   % (1U << 24 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  24  )   % (1U << 24 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 16 ))<<( 24 - 16 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  16  )   % (1U << 24 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 8 ))<<( 24 - 8 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  8  )   % (1U << 24 ) ;
    in += 4;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  0  )   % (1U << 24 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  24  )   % (1U << 24 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 16 ))<<( 24 - 16 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  16  )   % (1U << 24 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 8 ))<<( 24 - 8 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  8  )   % (1U << 24 ) ;
    in += 4;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  0  )   % (1U << 24 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  24  )   % (1U << 24 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 16 ))<<( 24 - 16 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  16  )   % (1U << 24 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 8 ))<<( 24 - 8 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  8  )   % (1U << 24 ) ;
    in += 4;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  0  )   % (1U << 24 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  24  )   % (1U << 24 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 16 ))<<( 24 - 16 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  16  )   % (1U << 24 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 8 ))<<( 24 - 8 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  8  )   % (1U << 24 ) ;
    out++;
  }

  return;
}

#elif !defined(HAVE_SSE2)
static void
unpack_24 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    *out = ( (*in) >>  0  )   % (1U << 24 ) ;
    out++;
    *out = ( (*in) >>  24  )   % (1U << 24 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 16 ))<<( 24 - 16 );
    out++;
    *out = ( (*in) >>  16  )   % (1U << 24 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 8 ))<<( 24 - 8 );
    out++;
    *out = ( (*in) >>  8  )   % (1U << 24 ) ;
    in += 4;
    out++;
    *out = ( (*in) >>  0  )   % (1U << 24 ) ;
    out++;
    *out = ( (*in) >>  24  )   % (1U << 24 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 16 ))<<( 24 - 16 );
    out++;
    *out = ( (*in) >>  16  )   % (1U << 24 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 8 ))<<( 24 - 8 );
    out++;
    *out = ( (*in) >>  8  )   % (1U << 24 ) ;
    in += 4;
    out++;
    *out = ( (*in) >>  0  )   % (1U << 24 ) ;
    out++;
    *out = ( (*in) >>  24  )   % (1U << 24 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 16 ))<<( 24 - 16 );
    out++;
    *out = ( (*in) >>  16  )   % (1U << 24 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 8 ))<<( 24 - 8 );
    out++;
    *out = ( (*in) >>  8  )   % (1U << 24 ) ;
    in += 4;
    out++;
    *out = ( (*in) >>  0  )   % (1U << 24 ) ;
    out++;
    *out = ( (*in) >>  24  )   % (1U << 24 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 16 ))<<( 24 - 16 );
    out++;
    *out = ( (*in) >>  16  )   % (1U << 24 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 8 ))<<( 24 - 8 );
    out++;
    *out = ( (*in) >>  8  )   % (1U << 24 ) ;
    out++;
  }

  return;
}

#else
static void
unpack_24_fwd (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg;
    const __m128i mask24 =  _mm_set1_epi32(16777215U);

    OutReg = _mm_and_si128( InReg , mask24);
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 16), mask24));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 8), mask24));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128( InReg , mask24);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 16), mask24));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 8), mask24));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_24_fwd_1_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask24 =  _mm_set1_epi32(16777215U);

    /* 1 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask24);
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 3 */
    InReg = _mm_load_si128(++in);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 8), mask24));
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_24_fwd_2_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask24 =  _mm_set1_epi32(16777215U);

    /* 2 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask24);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 16), mask24));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);


    /* 4 */
    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 8), mask24));
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_24_fwd_3_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask24 =  _mm_set1_epi32(16777215U);

    /* 3 */
    InReg = _mm_load_si128(++in);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 8), mask24));
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 5 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask24);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_24_fwd_4_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask24 =  _mm_set1_epi32(16777215U);

    /* 4 */
    InReg = _mm_load_si128(++in);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 8), mask24));
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);


    /* 6 */
    total = /* OutReg = */ _mm_and_si128( InReg , mask24);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 16), mask24));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_24_fwd_5_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask24 =  _mm_set1_epi32(16777215U);

    /* 5 */
    in += 3;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask24);
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 7 */
    InReg = _mm_load_si128(++in);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 8), mask24));
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_24_fwd_6_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask24 =  _mm_set1_epi32(16777215U);

    /* 6 */
    in += 3;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask24);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 16), mask24));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);


    /* 8 */
    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 8), mask24));
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_24_fwd_7_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask24 =  _mm_set1_epi32(16777215U);

    /* 1 first */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask24);
    _mm_store_si128(&(out[2]), total);


    /* 7 second */
    in += 4;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 8), mask24));
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_24_fwd_8_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask24 =  _mm_set1_epi32(16777215U);

    /* 2 first */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask24);
    _mm_store_si128(&(out[2]), total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 16), mask24));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(&(out[3]), total);


    /* 8 second */
    in += (4 - 1);
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 8), mask24));
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_24_rev (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg;
    const __m128i mask24 =  _mm_set1_epi32(16777215U);

    in += 6;
    InReg = _mm_load_si128(in);


    OutReg = _mm_and_si128( InReg , mask24);
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 16), mask24));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 8), mask24));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128( InReg , mask24);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 16), mask24));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 8), mask24));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_24_rev_1_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask24 =  _mm_set1_epi32(16777215U);

    /* 1 */
    in += 6;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask24);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 3 */
    InReg = _mm_load_si128(++in);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 8), mask24));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_24_rev_2_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask24 =  _mm_set1_epi32(16777215U);

    /* 2 */
    in += 6;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask24);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 16), mask24));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);


    /* 4 */
    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 8), mask24));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_24_rev_3_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask24 =  _mm_set1_epi32(16777215U);

    /* 3 */
    in += 7;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 8), mask24));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 5 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask24);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_24_rev_4_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask24 =  _mm_set1_epi32(16777215U);

    /* 4 */
    in += 7;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 8), mask24));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);


    /* 6 */
    total = /* OutReg = */ _mm_and_si128( InReg , mask24);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 16), mask24));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_24_rev_5_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask24 =  _mm_set1_epi32(16777215U);

    /* 5 */
    in += 9;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask24);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 7 */
    InReg = _mm_load_si128(++in);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 8), mask24));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_24_rev_6_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask24 =  _mm_set1_epi32(16777215U);

    /* 6 */
    in += 9;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask24);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 16), mask24));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);


    /* 8 */
    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 8), mask24));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,8) ;
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_24_rev_7_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask24 =  _mm_set1_epi32(16777215U);

    /* 1 first */
    in += 6;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask24);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(&(out[2]), total);


    /* 7 second */
    in += (10 - 6);
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 8), mask24));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_24_rev_8_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask24 =  _mm_set1_epi32(16777215U);

    /* 2 first */
    in += 6;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask24);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(&(out[2]), total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 16), mask24));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(&(out[3]), total);


    /* 8 second */
    in += (10 - 7);
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 8), mask24));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,8) ;
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}
#endif


#ifdef ALLOW_ODD_PACKSIZES
static void
unpack_25 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg;
    const __m128i mask25 =  _mm_set1_epi32(33554431U);

    OutReg = _mm_and_si128( InReg , mask25);
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,25) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 25 - 18), mask25));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,18) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 25 - 11), mask25));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,11) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 25 - 4), mask25));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask25);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,29) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 25 - 22), mask25));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 25 - 15), mask25));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,15) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 25 - 8), mask25));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 25 - 1), mask25));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,1) , mask25);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 25 - 19), mask25));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,19) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 25 - 12), mask25));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 25 - 5), mask25));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,5) , mask25);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 25 - 23), mask25));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,23) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 25 - 16), mask25));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}
#endif



#ifdef WORDS_BIGENDIAN
static void
unpack_26 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    *out = ( Bigendian_convert_uint(*in) >>  0  )   % (1U << 26 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  26  )   % (1U << 26 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 20 ))<<( 26 - 20 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  20  )   % (1U << 26 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 14 ))<<( 26 - 14 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  14  )   % (1U << 26 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 8 ))<<( 26 - 8 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  8  )   % (1U << 26 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 2 ))<<( 26 - 2 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  2  )   % (1U << 26 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  28  )   % (1U << 26 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 22 ))<<( 26 - 22 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  22  )   % (1U << 26 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 16 ))<<( 26 - 16 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  16  )   % (1U << 26 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 10 ))<<( 26 - 10 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  10  )   % (1U << 26 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 4 ))<<( 26 - 4 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  4  )   % (1U << 26 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  30  )   % (1U << 26 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 24 ))<<( 26 - 24 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  24  )   % (1U << 26 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 18 ))<<( 26 - 18 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  18  )   % (1U << 26 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 12 ))<<( 26 - 12 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  12  )   % (1U << 26 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 6 ))<<( 26 - 6 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  6  )   % (1U << 26 ) ;
    out++;
  }
  
  return;
}

#elif !defined(HAVE_SSE2)
static void
unpack_26 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    *out = ( (*in) >>  0  )   % (1U << 26 ) ;
    out++;
    *out = ( (*in) >>  26  )   % (1U << 26 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 20 ))<<( 26 - 20 );
    out++;
    *out = ( (*in) >>  20  )   % (1U << 26 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 14 ))<<( 26 - 14 );
    out++;
    *out = ( (*in) >>  14  )   % (1U << 26 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 8 ))<<( 26 - 8 );
    out++;
    *out = ( (*in) >>  8  )   % (1U << 26 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 2 ))<<( 26 - 2 );
    out++;
    *out = ( (*in) >>  2  )   % (1U << 26 ) ;
    out++;
    *out = ( (*in) >>  28  )   % (1U << 26 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 22 ))<<( 26 - 22 );
    out++;
    *out = ( (*in) >>  22  )   % (1U << 26 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 16 ))<<( 26 - 16 );
    out++;
    *out = ( (*in) >>  16  )   % (1U << 26 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 10 ))<<( 26 - 10 );
    out++;
    *out = ( (*in) >>  10  )   % (1U << 26 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 4 ))<<( 26 - 4 );
    out++;
    *out = ( (*in) >>  4  )   % (1U << 26 ) ;
    out++;
    *out = ( (*in) >>  30  )   % (1U << 26 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 24 ))<<( 26 - 24 );
    out++;
    *out = ( (*in) >>  24  )   % (1U << 26 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 18 ))<<( 26 - 18 );
    out++;
    *out = ( (*in) >>  18  )   % (1U << 26 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 12 ))<<( 26 - 12 );
    out++;
    *out = ( (*in) >>  12  )   % (1U << 26 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 6 ))<<( 26 - 6 );
    out++;
    *out = ( (*in) >>  6  )   % (1U << 26 ) ;
    out++;
  }
  
  return;
}

#else
static void
unpack_26_fwd (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg;
    const __m128i mask26 =  _mm_set1_epi32(67108863U);

    OutReg = _mm_and_si128( InReg , mask26);
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 20), mask26));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 14), mask26));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,14) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 8), mask26));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 2), mask26));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask26);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 22), mask26));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 16), mask26));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_26_fwd_1_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask26 =  _mm_set1_epi32(67108863U);

    /* 1 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask26);
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 3 */
    InReg = _mm_load_si128(++in);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 14), mask26));
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_26_fwd_2_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask26 =  _mm_set1_epi32(67108863U);

    /* 2 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask26);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 20), mask26));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);


    /* 4 */
    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 14), mask26));
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,14) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 8), mask26));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_26_fwd_3_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask26 =  _mm_set1_epi32(67108863U);

    /* 3 */
    InReg = _mm_load_si128(++in);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 14), mask26));
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 5 */
    InReg = _mm_load_si128(++in);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 2), mask26));
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_26_fwd_4_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask26 =  _mm_set1_epi32(67108863U);

    /* 4 */
    InReg = _mm_load_si128(++in);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 14), mask26));
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,14) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 8), mask26));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);


    /* 6 */
    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 2), mask26));
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask26);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_26_fwd_5_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask26 =  _mm_set1_epi32(67108863U);

    /* 5 */
    in += 3;    
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 2), mask26));
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 7 */
    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 22), mask26));
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_26_fwd_6_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask26 =  _mm_set1_epi32(67108863U);

    /* 6 */
    in += 3;    
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 2), mask26));
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask26);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);


    /* 8 */
    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 22), mask26));
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 16), mask26));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);
    return;
}

static void
unpack_26_fwd_7_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask26 =  _mm_set1_epi32(67108863U);

    /* 1 first */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask26);
    _mm_store_si128(&(out[2]), total);


    /* 7 second */
    in += 4;    
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 22), mask26));
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_26_fwd_8_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask26 =  _mm_set1_epi32(67108863U);

    /* 2 first */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask26);
    _mm_store_si128(&(out[2]), total);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 20), mask26));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(&(out[3]), total);


    /* 8 second */
    in += (4 - 1);    
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 22), mask26));
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 16), mask26));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_26_rev (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg;
    const __m128i mask26 =  _mm_set1_epi32(67108863U);

    in += 6;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 10), mask26));
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,10) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 4), mask26));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask26);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 24), mask26));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 18), mask26));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,18) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 12), mask26));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 6), mask26));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,6) ;
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_26_rev_1_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask26 =  _mm_set1_epi32(67108863U);

    /* 1 */
    in += 6;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 10), mask26));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 3 */
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask26);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_26_rev_2_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask26 =  _mm_set1_epi32(67108863U);

    /* 2 */
    in += 6;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 10), mask26));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,10) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 4), mask26));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);


    /* 4 */
    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask26);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 24), mask26));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_26_rev_3_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask26 =  _mm_set1_epi32(67108863U);

    /* 3 */
    in += 8;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask26);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 5 */
    InReg = _mm_load_si128(++in);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 18), mask26));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_26_rev_4_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask26 =  _mm_set1_epi32(67108863U);

    /* 4 */
    in += 8;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask26);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 24), mask26));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);


    /* 6 */
    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 18), mask26));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,18) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 12), mask26));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_26_rev_5_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask26 =  _mm_set1_epi32(67108863U);

    /* 5 */
    in += 9;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 18), mask26));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 7 */
    InReg = _mm_load_si128(++in);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 6), mask26));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_26_rev_6_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask26 =  _mm_set1_epi32(67108863U);

    /* 6 */
    in += 9;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 18), mask26));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,18) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 12), mask26));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);


    /* 8 */
    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 6), mask26));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,6) ;
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_26_rev_7_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask26 =  _mm_set1_epi32(67108863U);

    /* 1 first */
    in += 6;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 10), mask26));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(&(out[2]), total);


    /* 7 second */
    in += (11 - 7);
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 6), mask26));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_26_rev_8_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask26 =  _mm_set1_epi32(67108863U);

    /* 2 first */
    in += 6;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 10), mask26));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(&(out[2]), total);

    OutReg =   _mm_srli_epi32(InReg,10) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 4), mask26));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(&(out[3]), total);


    /* 8 second */
    in += (11 - 8);
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 6), mask26));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,6) ;
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}
#endif


#ifdef ALLOW_ODD_PACKSIZES
static void
unpack_27 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg;
    const __m128i mask27 =  _mm_set1_epi32(134217727U);

    OutReg = _mm_and_si128( InReg , mask27);
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,27) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 27 - 22), mask27));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 27 - 17), mask27));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,17) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 27 - 12), mask27));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 27 - 7), mask27));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,7) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 27 - 2), mask27));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask27);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,29) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 27 - 24), mask27));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 27 - 19), mask27));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,19) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 27 - 14), mask27));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,14) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 27 - 9), mask27));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,9) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 27 - 4), mask27));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask27);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,31) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 27 - 26), mask27));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 27 - 21), mask27));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,21) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 27 - 16), mask27));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}
#endif



#ifdef WORDS_BIGENDIAN
static void
unpack_28 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    *out = ( Bigendian_convert_uint(*in) >>  0  )   % (1U << 28 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  28  )   % (1U << 28 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 24 ))<<( 28 - 24 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  24  )   % (1U << 28 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 20 ))<<( 28 - 20 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  20  )   % (1U << 28 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 16 ))<<( 28 - 16 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  16  )   % (1U << 28 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 12 ))<<( 28 - 12 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  12  )   % (1U << 28 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 8 ))<<( 28 - 8 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  8  )   % (1U << 28 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 4 ))<<( 28 - 4 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  4  )   % (1U << 28 ) ;
    in += 4;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  0  )   % (1U << 28 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  28  )   % (1U << 28 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 24 ))<<( 28 - 24 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  24  )   % (1U << 28 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 20 ))<<( 28 - 20 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  20  )   % (1U << 28 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 16 ))<<( 28 - 16 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  16  )   % (1U << 28 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 12 ))<<( 28 - 12 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  12  )   % (1U << 28 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 8 ))<<( 28 - 8 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  8  )   % (1U << 28 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 4 ))<<( 28 - 4 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  4  )   % (1U << 28 ) ;
    out++;
  }

  return;
}

#elif !defined(HAVE_SSE2)
static void
unpack_28 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    *out = ( (*in) >>  0  )   % (1U << 28 ) ;
    out++;
    *out = ( (*in) >>  28  )   % (1U << 28 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 24 ))<<( 28 - 24 );
    out++;
    *out = ( (*in) >>  24  )   % (1U << 28 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 20 ))<<( 28 - 20 );
    out++;
    *out = ( (*in) >>  20  )   % (1U << 28 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 16 ))<<( 28 - 16 );
    out++;
    *out = ( (*in) >>  16  )   % (1U << 28 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 12 ))<<( 28 - 12 );
    out++;
    *out = ( (*in) >>  12  )   % (1U << 28 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 8 ))<<( 28 - 8 );
    out++;
    *out = ( (*in) >>  8  )   % (1U << 28 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 4 ))<<( 28 - 4 );
    out++;
    *out = ( (*in) >>  4  )   % (1U << 28 ) ;
    in += 4;
    out++;
    *out = ( (*in) >>  0  )   % (1U << 28 ) ;
    out++;
    *out = ( (*in) >>  28  )   % (1U << 28 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 24 ))<<( 28 - 24 );
    out++;
    *out = ( (*in) >>  24  )   % (1U << 28 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 20 ))<<( 28 - 20 );
    out++;
    *out = ( (*in) >>  20  )   % (1U << 28 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 16 ))<<( 28 - 16 );
    out++;
    *out = ( (*in) >>  16  )   % (1U << 28 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 12 ))<<( 28 - 12 );
    out++;
    *out = ( (*in) >>  12  )   % (1U << 28 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 8 ))<<( 28 - 8 );
    out++;
    *out = ( (*in) >>  8  )   % (1U << 28 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 4 ))<<( 28 - 4 );
    out++;
    *out = ( (*in) >>  4  )   % (1U << 28 ) ;
    out++;
  }

  return;
}

#else
static void
unpack_28_fwd (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg;
    const __m128i mask28 =  _mm_set1_epi32(268435455U);

    OutReg = _mm_and_si128( InReg , mask28);
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 24), mask28));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 20), mask28));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 16), mask28));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 12), mask28));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 8), mask28));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 4), mask28));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,4) ;
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_28_fwd_1_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask28 =  _mm_set1_epi32(268435455U);

    /* 1 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask28);
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 3 */
    InReg = _mm_load_si128(++in);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 20), mask28));
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_28_fwd_2_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask28 =  _mm_set1_epi32(268435455U);

    /* 2 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask28);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 24), mask28));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);


    /* 4 */
    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 20), mask28));
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 16), mask28));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_28_fwd_3_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask28 =  _mm_set1_epi32(268435455U);

    /* 3 */
    InReg = _mm_load_si128(++in);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 20), mask28));
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 5 */
    InReg = _mm_load_si128(++in);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 12), mask28));
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_28_fwd_4_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask28 =  _mm_set1_epi32(268435455U);

    /* 4 */
    InReg = _mm_load_si128(++in);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 20), mask28));
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 16), mask28));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);


    /* 6 */
    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 12), mask28));
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 8), mask28));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_28_fwd_5_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask28 =  _mm_set1_epi32(268435455U);

    /* 5 */
    in += 3;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 12), mask28));
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 7 */
    InReg = _mm_load_si128(++in);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 4), mask28));
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_28_fwd_6_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask28 =  _mm_set1_epi32(268435455U);

    /* 6 */
    in += 3;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 12), mask28));
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 8), mask28));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);


    /* 8 */
    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 4), mask28));
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,4) ;
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_28_fwd_7_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask28 =  _mm_set1_epi32(268435455U);

    /* 1 first  */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask28);
    _mm_store_si128(&(out[2]), total);


    /* 7 second */
    in += 5;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 4), mask28));
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_28_fwd_8_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask28 =  _mm_set1_epi32(268435455U);

    /* 2 first */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask28);
    _mm_store_si128(&(out[2]), total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 24), mask28));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(&(out[3]), total);


    /* 8 second */
    in += (5 - 1);
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 4), mask28));
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,4) ;
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_28_rev (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg;
    const __m128i mask28 =  _mm_set1_epi32(268435455U);

    in += 7;
    InReg = _mm_load_si128(in);


    OutReg = _mm_and_si128( InReg , mask28);
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 24), mask28));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 20), mask28));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 16), mask28));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 12), mask28));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 8), mask28));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 4), mask28));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,4) ;
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_28_rev_1_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask28 =  _mm_set1_epi32(268435455U);

    /* 1 */
    in += 7;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask28);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 3 */
    InReg = _mm_load_si128(++in);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 20), mask28));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_28_rev_2_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask28 =  _mm_set1_epi32(268435455U);

    /* 2 */
    in += 7;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask28);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 24), mask28));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);


    /* 4 */
    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 20), mask28));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 16), mask28));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);
    return;
}

static void
unpack_28_rev_3_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask28 =  _mm_set1_epi32(268435455U);

    /* 3 */
    in += 8;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 20), mask28));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 5 */
    InReg = _mm_load_si128(++in);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 12), mask28));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_28_rev_4_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask28 =  _mm_set1_epi32(268435455U);

    /* 4 */
    in += 8;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 20), mask28));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 16), mask28));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);


    /* 6 */
    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 12), mask28));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 8), mask28));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_28_rev_5_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask28 =  _mm_set1_epi32(268435455U);

    /* 5 */
    in += 10;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 12), mask28));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 7 */
    InReg = _mm_load_si128(++in);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 4), mask28));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_28_rev_6_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask28 =  _mm_set1_epi32(268435455U);

    /* 6 */
    in += 10;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 12), mask28));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 8), mask28));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);


    /* 8 */
    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 4), mask28));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,4) ;
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_28_rev_7_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask28 =  _mm_set1_epi32(268435455U);

    /* 1 first */
    in += 7;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask28);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(&(out[2]), total);


    /* 7 second */
    in += (12 - 7);
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 4), mask28));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_28_rev_8_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask28 =  _mm_set1_epi32(268435455U);

    /* 2 first */
    in += 7;
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask28);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(&(out[2]), total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 24), mask28));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(&(out[3]), total);


    /* 8 second */
    in += (12 - 8);
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 4), mask28));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,4) ;
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}
#endif


#ifdef ALLOW_ODD_PACKSIZES
static void
unpack_29 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg;
    const __m128i mask29 =  _mm_set1_epi32(536870911U);

    OutReg = _mm_and_si128( InReg , mask29);
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,29) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 29 - 26), mask29));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 29 - 23), mask29));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,23) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 29 - 20), mask29));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 29 - 17), mask29));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,17) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 29 - 14), mask29));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,14) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 29 - 11), mask29));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,11) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 29 - 8), mask29));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 29 - 5), mask29));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,5) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 29 - 2), mask29));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask29);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,31) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 29 - 28), mask29));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 29 - 25), mask29));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,25) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 29 - 22), mask29));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 29 - 19), mask29));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,19) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 29 - 16), mask29));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}
#endif



#ifdef WORDS_BIGENDIAN
static void
unpack_30 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    *out = ( Bigendian_convert_uint(*in) >>  0  )   % (1U << 30 ) ;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  30  )   % (1U << 30 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 28 ))<<( 30 - 28 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  28  )   % (1U << 30 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 26 ))<<( 30 - 26 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  26  )   % (1U << 30 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 24 ))<<( 30 - 24 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  24  )   % (1U << 30 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 22 ))<<( 30 - 22 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  22  )   % (1U << 30 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 20 ))<<( 30 - 20 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  20  )   % (1U << 30 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 18 ))<<( 30 - 18 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  18  )   % (1U << 30 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 16 ))<<( 30 - 16 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  16  )   % (1U << 30 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 14 ))<<( 30 - 14 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  14  )   % (1U << 30 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 12 ))<<( 30 - 12 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  12  )   % (1U << 30 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 10 ))<<( 30 - 10 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  10  )   % (1U << 30 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 8 ))<<( 30 - 8 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  8  )   % (1U << 30 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 6 ))<<( 30 - 6 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  6  )   % (1U << 30 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 4 ))<<( 30 - 4 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  4  )   % (1U << 30 ) ;
    in += 4;
    *out |= (Bigendian_convert_uint(*in) % (1U<< 2 ))<<( 30 - 2 );
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  2  )   % (1U << 30 ) ;
    out++;
  }

  return;
}

#elif !defined(HAVE_SSE2)
static void
unpack_30 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    *out = ( (*in) >>  0  )   % (1U << 30 ) ;
    out++;
    *out = ( (*in) >>  30  )   % (1U << 30 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 28 ))<<( 30 - 28 );
    out++;
    *out = ( (*in) >>  28  )   % (1U << 30 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 26 ))<<( 30 - 26 );
    out++;
    *out = ( (*in) >>  26  )   % (1U << 30 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 24 ))<<( 30 - 24 );
    out++;
    *out = ( (*in) >>  24  )   % (1U << 30 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 22 ))<<( 30 - 22 );
    out++;
    *out = ( (*in) >>  22  )   % (1U << 30 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 20 ))<<( 30 - 20 );
    out++;
    *out = ( (*in) >>  20  )   % (1U << 30 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 18 ))<<( 30 - 18 );
    out++;
    *out = ( (*in) >>  18  )   % (1U << 30 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 16 ))<<( 30 - 16 );
    out++;
    *out = ( (*in) >>  16  )   % (1U << 30 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 14 ))<<( 30 - 14 );
    out++;
    *out = ( (*in) >>  14  )   % (1U << 30 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 12 ))<<( 30 - 12 );
    out++;
    *out = ( (*in) >>  12  )   % (1U << 30 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 10 ))<<( 30 - 10 );
    out++;
    *out = ( (*in) >>  10  )   % (1U << 30 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 8 ))<<( 30 - 8 );
    out++;
    *out = ( (*in) >>  8  )   % (1U << 30 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 6 ))<<( 30 - 6 );
    out++;
    *out = ( (*in) >>  6  )   % (1U << 30 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 4 ))<<( 30 - 4 );
    out++;
    *out = ( (*in) >>  4  )   % (1U << 30 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 2 ))<<( 30 - 2 );
    out++;
    *out = ( (*in) >>  2  )   % (1U << 30 ) ;
    out++;
  }

  return;
}

#else
static void
unpack_30_fwd (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg;
    const __m128i mask30 =  _mm_set1_epi32(1073741823U);

    OutReg = _mm_and_si128( InReg , mask30);
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 28), mask30));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 26), mask30));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 24), mask30));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 22), mask30));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 20), mask30));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 18), mask30));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,18) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 16), mask30));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_30_fwd_1_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask30 =  _mm_set1_epi32(1073741823U);

    /* 1 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask30);
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 3 */
    InReg = _mm_load_si128(++in);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 26), mask30));
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_30_fwd_2_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask30 =  _mm_set1_epi32(1073741823U);

    /* 2 */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask30);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 28), mask30));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);


    /* 4 */
    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 26), mask30));
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 24), mask30));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_30_fwd_3_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask30 =  _mm_set1_epi32(1073741823U);

    /* 3 */
    InReg = _mm_load_si128(++in);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 26), mask30));
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 5 */
    InReg = _mm_load_si128(++in);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 22), mask30));
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_30_fwd_4_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask30 =  _mm_set1_epi32(1073741823U);

    /* 4 */
    InReg = _mm_load_si128(++in);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 26), mask30));
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 24), mask30));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);


    /* 6 */
    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 22), mask30));
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 20), mask30));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_30_fwd_5_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask30 =  _mm_set1_epi32(1073741823U);

    /* 5 */
    in += 3;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 22), mask30));
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 7 */
    InReg = _mm_load_si128(++in);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 18), mask30));
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_30_fwd_6_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask30 =  _mm_set1_epi32(1073741823U);

    /* 6 */
    in += 3;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 22), mask30));
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 20), mask30));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);


    /* 8 */
    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 18), mask30));
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,18) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 16), mask30));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_30_fwd_7_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask30 =  _mm_set1_epi32(1073741823U);

    /* 1 first */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask30);
    _mm_store_si128(&(out[2]), total);


    /* 7 second */
    in += 5;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 18), mask30));
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_30_fwd_8_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask30 =  _mm_set1_epi32(1073741823U);

    /* 2 first */
    InReg = _mm_load_si128(in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask30);
    _mm_store_si128(&(out[2]), total);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 28), mask30));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(&(out[3]), total);


    /* 8 second */
    in += (5 - 1);
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 18), mask30));
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,18) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 16), mask30));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_30_rev (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg;
    const __m128i mask30 =  _mm_set1_epi32(1073741823U);

    in += 7;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 14), mask30));
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,14) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 12), mask30));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 10), mask30));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,10) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 8), mask30));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 6), mask30));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,6) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 4), mask30));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,4) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 2), mask30));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,2) ;
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_30_rev_1_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask30 =  _mm_set1_epi32(1073741823U);

    /* 1 */
    in += 7;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 14), mask30));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 3 */
    InReg = _mm_load_si128(++in);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 10), mask30));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_30_rev_2_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask30 =  _mm_set1_epi32(1073741823U);

    /* 2 */
    in += 7;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 14), mask30));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,14) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 12), mask30));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);


    /* 4 */
    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 10), mask30));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,10) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 8), mask30));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);
    return;
}

static void
unpack_30_rev_3_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask30 =  _mm_set1_epi32(1073741823U);

    /* 3 */
    in += 9;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 10), mask30));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 5 */
    InReg = _mm_load_si128(++in);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 6), mask30));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_30_rev_4_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask30 =  _mm_set1_epi32(1073741823U);

    /* 4 */
    in += 9;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 10), mask30));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,10) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 8), mask30));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);


    /* 6 */
    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 6), mask30));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,6) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 4), mask30));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_30_rev_5_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask30 =  _mm_set1_epi32(1073741823U);

    /* 5 */
    in += 11;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 6), mask30));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 7 */
    InReg = _mm_load_si128(++in);

    OutReg =   _mm_srli_epi32(InReg,4) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 2), mask30));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_30_rev_6_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask30 =  _mm_set1_epi32(1073741823U);

    /* 6 */
    in += 11;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 6), mask30));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,6) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 4), mask30));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);


    /* 8 */
    OutReg =   _mm_srli_epi32(InReg,4) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 2), mask30));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,2) ;
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_30_rev_7_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask30 =  _mm_set1_epi32(1073741823U);

    /* 1 first */
    in += 7;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 14), mask30));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(&(out[2]), total);


    /* 7 second */
    in += (13 - 8);
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,4) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 2), mask30));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_30_rev_8_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask30 =  _mm_set1_epi32(1073741823U);

    /* 2 first */
    in += 7;
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 14), mask30));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(&(out[2]), total);

    OutReg =   _mm_srli_epi32(InReg,14) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 12), mask30));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(&(out[3]), total);


    /* 8 second */
    in += (13 - 9);
    InReg = _mm_load_si128(in);

    OutReg =   _mm_srli_epi32(InReg,4) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 2), mask30));
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,2) ;
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}
#endif


#ifdef ALLOW_ODD_PACKSIZES
static void
unpack_31 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg;
    const __m128i mask31 =  _mm_set1_epi32(2147483647U);

    OutReg = _mm_and_si128( InReg , mask31);
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,31) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 31 - 30), mask31));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 31 - 29), mask31));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,29) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 31 - 28), mask31));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 31 - 27), mask31));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,27) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 31 - 26), mask31));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 31 - 25), mask31));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,25) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 31 - 24), mask31));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 31 - 23), mask31));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,23) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 31 - 22), mask31));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 31 - 21), mask31));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,21) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 31 - 20), mask31));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 31 - 19), mask31));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,19) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 31 - 18), mask31));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,18) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 31 - 17), mask31));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,17) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 31 - 16), mask31));
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}
#endif



#ifdef WORDS_BIGENDIAN
static void
unpack_32 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    *out = ( Bigendian_convert_uint(*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( Bigendian_convert_uint(*in) >>  0  )   ;
    out++;
  }

  return;
}

#elif !defined(HAVE_SSE2)
static void
unpack_32 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    *out = ( (*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( (*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( (*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( (*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( (*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( (*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( (*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( (*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( (*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( (*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( (*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( (*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( (*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( (*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( (*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( (*in) >>  0  )   ;
    out++;
  }

  return;
}

#else
static void
unpack_32_fwd (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i OutReg;

    OutReg = _mm_load_si128(in++);
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_load_si128(in++);
    /* total = _mm_add_epi32(total, _mm_load_si128(in++)); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_load_si128(in++);
    /* total = _mm_add_epi32(total, _mm_load_si128(in++)); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_load_si128(in++);
    /* total = _mm_add_epi32(total, _mm_load_si128(in++)); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_load_si128(in++);
    /* total = _mm_add_epi32(total, _mm_load_si128(in++)); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_load_si128(in++);
    /* total = _mm_add_epi32(total, _mm_load_si128(in++)); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_load_si128(in++);
    /* total = _mm_add_epi32(total, _mm_load_si128(in++)); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_load_si128(in++);
    /* total = _mm_add_epi32(total, _mm_load_si128(in++)); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_32_fwd_1_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i total;

    /* 1 */
    total = _mm_load_si128(in);
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 3 */
    in += 2;
    total = _mm_load_si128(in);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_32_fwd_2_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i total;

    /* 2 */
    total = _mm_load_si128(in);
    _mm_store_si128(out++, total);

    total = _mm_add_epi32(total, _mm_load_si128(++in));
    _mm_store_si128(out++, total);

    /* 4 */
    total = _mm_load_si128(++in);
    _mm_store_si128(out++, total);

    total = _mm_add_epi32(total, _mm_load_si128(++in));
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_32_fwd_3_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i total;

    /* 3 */
    in += 2;
    total = _mm_load_si128(in);
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 5 */
    in += 2;
    total = _mm_load_si128(in);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_32_fwd_4_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i total;

    /* 4 */
    in += 2;
    total = _mm_load_si128(in);
    _mm_store_si128(out++, total);

    total = _mm_add_epi32(total, _mm_load_si128(++in));
    _mm_store_si128(out++, total);

    /* 6 */
    total = _mm_load_si128(++in);
    _mm_store_si128(out++, total);

    total = _mm_add_epi32(total, _mm_load_si128(++in));
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_32_fwd_5_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i total;

    /* 5 */
    in += 4;
    total = _mm_load_si128(in);
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 7 */
    in += 2;
    total = _mm_load_si128(in);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_32_fwd_6_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i total;

    /* 6 */
    in += 4;
    total = _mm_load_si128(in);
    _mm_store_si128(out++, total);

    total = _mm_add_epi32(total, _mm_load_si128(++in));
    _mm_store_si128(out++, total);

    /* 8 */
    total = _mm_load_si128(++in);
    _mm_store_si128(out++, total);

    total = _mm_add_epi32(total, _mm_load_si128(++in));
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_32_fwd_7_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i total;

    /* 1 first */
    total = _mm_load_si128(in);
    _mm_store_si128(&(out[2]), total);


    /* 7 second */
    in += 6;
    total = _mm_load_si128(in);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_32_fwd_8_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i total;

    /* 2 first */
    total = _mm_load_si128(in);
    _mm_store_si128(&(out[2]), total);

    total = _mm_add_epi32(total, _mm_load_si128(++in));
    _mm_store_si128(&(out[3]), total);


    /* 8 second */
    in += (6 - 1);
    total = _mm_load_si128(in);
    _mm_store_si128(out++, total);

    total = _mm_add_epi32(total, _mm_load_si128(++in));
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_32_rev (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i OutReg;

    in += 8;

    OutReg = _mm_load_si128(in++);
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_load_si128(in++);
    /* total = _mm_add_epi32(total, _mm_load_si128(in++)); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_load_si128(in++);
    /* total = _mm_add_epi32(total, _mm_load_si128(in++)); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_load_si128(in++);
    /* total = _mm_add_epi32(total, _mm_load_si128(in++)); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_load_si128(in++);
    /* total = _mm_add_epi32(total, _mm_load_si128(in++)); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_load_si128(in++);
    /* total = _mm_add_epi32(total, _mm_load_si128(in++)); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_load_si128(in++);
    /* total = _mm_add_epi32(total, _mm_load_si128(in++)); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_load_si128(in++);
    /* total = _mm_add_epi32(total, _mm_load_si128(in++)); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_32_rev_1_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i total;

    /* 1 */
    in += 8;
    total = _mm_load_si128(in);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 3 */
    in += 2;
    total = _mm_load_si128(in);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_32_rev_2_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i total;

    /* 2 */
    in += 8;
    total = _mm_load_si128(in);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, _mm_load_si128(++in));
#else
    total = _mm_add_epi32(total, _mm_load_si128(++in));
#endif
    _mm_store_si128(out++, total);


    /* 4 */
    total = _mm_load_si128(++in);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, _mm_load_si128(++in));
#else
    total = _mm_add_epi32(total, _mm_load_si128(++in));
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_32_rev_3_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i total;

    /* 3 */
    in += 10;
    total = _mm_load_si128(in);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 5 */
    in += 2;
    total = _mm_load_si128(in);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_32_rev_4_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i total;

    /* 4 */
    in += 10;
    total = _mm_load_si128(in);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, _mm_load_si128(++in));
#else
    total = _mm_add_epi32(total, _mm_load_si128(++in));
#endif
    _mm_store_si128(out++, total);


    /* 6 */
    total = _mm_load_si128(++in);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, _mm_load_si128(++in));
#else
    total = _mm_add_epi32(total, _mm_load_si128(++in));
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_32_rev_5_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i total;

    /* 5 */
    in += 12;
    total = _mm_load_si128(in);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    /* Skip row */
    out++;

    /* 7 */
    in += 2;
    total = _mm_load_si128(in);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_32_rev_6_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i total;

    /* 6 */
    in += 12;
    total = _mm_load_si128(in);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, _mm_load_si128(++in));
#else
    total = _mm_add_epi32(total, _mm_load_si128(++in));
#endif
    _mm_store_si128(out++, total);


    /* 8 */
    total = _mm_load_si128(++in);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, _mm_load_si128(++in));
#else
    total = _mm_add_epi32(total, _mm_load_si128(++in));
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_32_rev_7_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i total;

    /* 1 first */
    in += 8;
    total = _mm_load_si128(in);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(&(out[2]), total);


    /* 7 second */
    in += (14 - 8);
    total = _mm_load_si128(in);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_32_rev_8_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i total;

    /* 2 first */
    in += 8;
    total = _mm_load_si128(in);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(&(out[2]), total);

#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, _mm_load_si128(++in));
#else
    total = _mm_add_epi32(total, _mm_load_si128(++in));
#endif
    _mm_store_si128(&(out[3]), total);


    /* 8 second */
    in += (14 - 9);
    total = _mm_load_si128(in);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, _mm_load_si128(++in));
#else
    total = _mm_add_epi32(total, _mm_load_si128(++in));
#endif
    _mm_store_si128(out++, total);

    return;
}
#endif


#ifdef PARTIAL_BRANCH_ROW_SUM

#define assign_sum_fwd(x,offset0,diffs,row) switch (row) {	\
  case -1: x = offset0; break;					\
  case 0: x = offset0 + diffs[0]; break;			\
  case 1: x = offset0 + diffs[0] + diffs[1]; break;			\
  case 2: x = offset0 + diffs[0] + diffs[1] + diffs[2]; break;		\
  default: x = offset0 + diffs[row-3] + diffs[row-2] + diffs[row-1] + diffs[row]; \
  }

#define assign_sum_rev(x,offset1,diffs,row) switch (row) {	\
  case -1: x = offset1; break;					\
  case 0: x = offset1 - diffs[0]; break;			\
  case 1: x = offset1 - diffs[0] - diffs[1]; break;			\
  case 2: x = offset1 - diffs[0] - diffs[1] - diffs[2]; break;		\
  default: x = offset1 - diffs[row-3] - diffs[row-2] - diffs[row-1] - diffs[row]; \
  }

#define return_sum_fwd(offset0,diffs,row) switch (row) {	\
  case -1: return offset0;					\
  case 0: return offset0 + diffs[0];				\
  case 1: return offset0 + diffs[0] + diffs[1];			\
  case 2: return offset0 + diffs[0] + diffs[1] + diffs[2];		\
  default: return offset0 + diffs[row-3] + diffs[row-2] + diffs[row-1] + diffs[row]; \
  }

#define return_sum_rev(offset1,diffs,row) switch (row) {	\
  case -1: return offset1;					\
  case 0: return offset1 - diffs[0];				\
  case 1: return offset1 - diffs[0] - diffs[1];			\
  case 2: return offset1 - diffs[0] - diffs[1] - diffs[2];		\
  default: return offset1 - diffs[row-3] - diffs[row-2] - diffs[row-1] - diffs[row]; \
  }

#else

#define assign_sum_fwd(x,offset0,diffs,row) switch (row) {	\
  case -1: x = offset0; break;						\
  case 0: x = offset0 + diffs[0]; break;				\
  case 1: x = offset0 + diffs[0] + diffs[1]; break;			\
  case 2: x = offset0 + diffs[0] + diffs[1] + diffs[2]; break;		\
  case 3: x = offset0 + diffs[0] + diffs[1] + diffs[2] + diffs[3]; break; \
  case 4: x = offset0 + diffs[1] + diffs[2] + diffs[3] + diffs[4]; break; \
  case 5: x = offset0 + diffs[2] + diffs[3] + diffs[4] + diffs[5]; break; \
  case 6: x = offset0 + diffs[3] + diffs[4] + diffs[5] + diffs[6]; break; \
  case 7: x = offset0 + diffs[4] + diffs[5] + diffs[6] + diffs[7]; break; \
  default: abort();							\
  }

#define assign_sum_rev(x,offset1,diffs,row) switch (row) {	\
  case -1: x = offset1; break;						\
  case 0: x = offset1 - diffs[0]; break;				\
  case 1: x = offset1 - diffs[0] - diffs[1]; break;			\
  case 2: x = offset1 - diffs[0] - diffs[1] - diffs[2]; break;		\
  case 3: x = offset1 - diffs[0] - diffs[1] - diffs[2] - diffs[3]; break; \
  case 4: x = offset1 - diffs[1] - diffs[2] - diffs[3] - diffs[4]; break; \
  case 5: x = offset1 - diffs[2] - diffs[3] - diffs[4] - diffs[5]; break; \
  case 6: x = offset1 - diffs[3] - diffs[4] - diffs[5] - diffs[6]; break; \
  case 7: x = offset1 - diffs[4] - diffs[5] - diffs[6] - diffs[7]; break; \
  default: abort();							\
  }

#define return_sum_fwd(offset0,diffs,row) switch (row) {	\
  case -1: return offset0; break;					\
  case 0: return offset0 + diffs[0]; break;				\
  case 1: return offset0 + diffs[0] + diffs[1]; break;			\
  case 2: return offset0 + diffs[0] + diffs[1] + diffs[2]; break;	\
  case 3: return offset0 + diffs[0] + diffs[1] + diffs[2] + diffs[3]; break; \
  case 4: return offset0 + diffs[1] + diffs[2] + diffs[3] + diffs[4]; break; \
  case 5: return offset0 + diffs[2] + diffs[3] + diffs[4] + diffs[5]; break; \
  case 6: return offset0 + diffs[3] + diffs[4] + diffs[5] + diffs[6]; break; \
  case 7: return offset0 + diffs[4] + diffs[5] + diffs[6] + diffs[7]; break; \
  default: abort();							\
  }

#define return_sum_rev(offset1,diffs,row) switch (row) {	\
  case -1: return offset1; break;					\
  case 0: return offset1 - diffs[0]; break;				\
  case 1: return offset1 - diffs[0] - diffs[1]; break;			\
  case 2: return offset1 - diffs[0] - diffs[1] - diffs[2]; break;	\
  case 3: return offset1 - diffs[0] - diffs[1] - diffs[2] - diffs[3]; break; \
  case 4: return offset1 - diffs[1] - diffs[2] - diffs[3] - diffs[4]; break; \
  case 5: return offset1 - diffs[2] - diffs[3] - diffs[4] - diffs[5]; break; \
  case 6: return offset1 - diffs[3] - diffs[4] - diffs[5] - diffs[6]; break; \
  case 7: return offset1 - diffs[4] - diffs[5] - diffs[6] - diffs[7]; break; \
  default: abort();							\
  }

#endif


#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
typedef void (*Unpacker_T) (UINT4* __restrict__, const UINT4* __restrict__);
#else
typedef void (*Unpacker_T) (__m128i* __restrict__, const __m128i* __restrict__);
#endif


#ifdef ALLOW_ODD_PACKSIZES
static Unpacker_T unpacker_table[33] =
  {unpack_00,
   unpack_01, unpack_02, unpack_03, unpack_04,
   unpack_05, unpack_06, unpack_07, unpack_08,
   unpack_09, unpack_10, unpack_11, unpack_12,
   unpack_13, unpack_14, unpack_15, unpack_16,
   unpack_17, unpack_18, unpack_19, unpack_20,
   unpack_21, unpack_22, unpack_23, unpack_24,
   unpack_25, unpack_26, unpack_27, unpack_28,
   unpack_29, unpack_30, unpack_31, unpack_32};

#elif defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
static Unpacker_T unpacker_all_table[33] =
  {unpack_00,
   unpack_00, unpack_02, unpack_00, unpack_04,
   unpack_00, unpack_06, unpack_00, unpack_08,
   unpack_00, unpack_10, unpack_00, unpack_12,
   unpack_00, unpack_14, unpack_00, unpack_16,
   unpack_00, unpack_18, unpack_00, unpack_20,
   unpack_00, unpack_22, unpack_00, unpack_24,
   unpack_00, unpack_26, unpack_00, unpack_28,
   unpack_00, unpack_30, unpack_00, unpack_32};

#else
static Unpacker_T unpacker_all_table[34] =
  {unpack_00, unpack_00,
   unpack_02_fwd, unpack_02_rev, unpack_04_fwd, unpack_04_rev,
   unpack_06_fwd, unpack_06_rev, unpack_08_fwd, unpack_08_rev,
   unpack_10_fwd, unpack_10_rev, unpack_12_fwd, unpack_12_rev,
   unpack_14_fwd, unpack_14_rev, unpack_16_fwd, unpack_16_rev,
   unpack_18_fwd, unpack_18_rev, unpack_20_fwd, unpack_20_rev,
   unpack_22_fwd, unpack_22_rev, unpack_24_fwd, unpack_24_rev,
   unpack_26_fwd, unpack_26_rev, unpack_28_fwd, unpack_28_rev,
   unpack_30_fwd, unpack_30_rev, unpack_32_fwd, unpack_32_rev};

/* Entry 16 in each packsize handles remainder == 64 => quarter_block == 4, column 3, row -1 */
static Unpacker_T unpacker_table[17][17] = 
  {{unpack_00_1_3, unpack_00_2_4, unpack_00_2_4, unpack_00_1_3,
    unpack_00_1_3, unpack_00_2_4, unpack_00_2_4, unpack_00_1_3,
    unpack_00_1_3, unpack_00_2_4, unpack_00_2_4, unpack_00_1_3,
    unpack_00_1_3, unpack_00_2_4, unpack_00_2_4, unpack_00_1_3,
    unpack_00_0},

   {unpack_02_fwd_1_3, unpack_02_fwd_2_4, unpack_02_rev_2_4, unpack_02_rev_1_3,
    unpack_02_fwd_3_5, unpack_02_fwd_4_6, unpack_02_rev_4_6, unpack_02_rev_3_5,
    unpack_02_fwd_5_7, unpack_02_fwd_6_8, unpack_02_rev_6_8, unpack_02_rev_5_7, 
    unpack_02_fwd_7_1, unpack_02_fwd_8_2, unpack_02_rev_8_2, unpack_02_rev_7_1,
    unpack_00_0},

   {unpack_04_fwd_1_3, unpack_04_fwd_2_4, unpack_04_rev_2_4, unpack_04_rev_1_3,
    unpack_04_fwd_3_5, unpack_04_fwd_4_6, unpack_04_rev_4_6, unpack_04_rev_3_5,
    unpack_04_fwd_5_7, unpack_04_fwd_6_8, unpack_04_rev_6_8, unpack_04_rev_5_7, 
    unpack_04_fwd_7_1, unpack_04_fwd_8_2, unpack_04_rev_8_2, unpack_04_rev_7_1,
    unpack_00_0},

   {unpack_06_fwd_1_3, unpack_06_fwd_2_4, unpack_06_rev_2_4, unpack_06_rev_1_3,
    unpack_06_fwd_3_5, unpack_06_fwd_4_6, unpack_06_rev_4_6, unpack_06_rev_3_5,
    unpack_06_fwd_5_7, unpack_06_fwd_6_8, unpack_06_rev_6_8, unpack_06_rev_5_7, 
    unpack_06_fwd_7_1, unpack_06_fwd_8_2, unpack_06_rev_8_2, unpack_06_rev_7_1,
    unpack_00_0},

   {unpack_08_fwd_1_3, unpack_08_fwd_2_4, unpack_08_rev_2_4, unpack_08_rev_1_3,
    unpack_08_fwd_3_5, unpack_08_fwd_4_6, unpack_08_rev_4_6, unpack_08_rev_3_5,
    unpack_08_fwd_5_7, unpack_08_fwd_6_8, unpack_08_rev_6_8, unpack_08_rev_5_7, 
    unpack_08_fwd_7_1, unpack_08_fwd_8_2, unpack_08_rev_8_2, unpack_08_rev_7_1,
    unpack_00_0},

   {unpack_10_fwd_1_3, unpack_10_fwd_2_4, unpack_10_rev_2_4, unpack_10_rev_1_3,
    unpack_10_fwd_3_5, unpack_10_fwd_4_6, unpack_10_rev_4_6, unpack_10_rev_3_5,
    unpack_10_fwd_5_7, unpack_10_fwd_6_8, unpack_10_rev_6_8, unpack_10_rev_5_7, 
    unpack_10_fwd_7_1, unpack_10_fwd_8_2, unpack_10_rev_8_2, unpack_10_rev_7_1,
    unpack_00_0},

   {unpack_12_fwd_1_3, unpack_12_fwd_2_4, unpack_12_rev_2_4, unpack_12_rev_1_3,
    unpack_12_fwd_3_5, unpack_12_fwd_4_6, unpack_12_rev_4_6, unpack_12_rev_3_5,
    unpack_12_fwd_5_7, unpack_12_fwd_6_8, unpack_12_rev_6_8, unpack_12_rev_5_7, 
    unpack_12_fwd_7_1, unpack_12_fwd_8_2, unpack_12_rev_8_2, unpack_12_rev_7_1,
    unpack_00_0},

   {unpack_14_fwd_1_3, unpack_14_fwd_2_4, unpack_14_rev_2_4, unpack_14_rev_1_3,
    unpack_14_fwd_3_5, unpack_14_fwd_4_6, unpack_14_rev_4_6, unpack_14_rev_3_5,
    unpack_14_fwd_5_7, unpack_14_fwd_6_8, unpack_14_rev_6_8, unpack_14_rev_5_7, 
    unpack_14_fwd_7_1, unpack_14_fwd_8_2, unpack_14_rev_8_2, unpack_14_rev_7_1,
    unpack_00_0},

   {unpack_16_fwd_1_3, unpack_16_fwd_2_4, unpack_16_rev_2_4, unpack_16_rev_1_3,
    unpack_16_fwd_3_5, unpack_16_fwd_4_6, unpack_16_rev_4_6, unpack_16_rev_3_5,
    unpack_16_fwd_5_7, unpack_16_fwd_6_8, unpack_16_rev_6_8, unpack_16_rev_5_7, 
    unpack_16_fwd_7_1, unpack_16_fwd_8_2, unpack_16_rev_8_2, unpack_16_rev_7_1,
    unpack_00_0},

   {unpack_18_fwd_1_3, unpack_18_fwd_2_4, unpack_18_rev_2_4, unpack_18_rev_1_3,
    unpack_18_fwd_3_5, unpack_18_fwd_4_6, unpack_18_rev_4_6, unpack_18_rev_3_5,
    unpack_18_fwd_5_7, unpack_18_fwd_6_8, unpack_18_rev_6_8, unpack_18_rev_5_7, 
    unpack_18_fwd_7_1, unpack_18_fwd_8_2, unpack_18_rev_8_2, unpack_18_rev_7_1,
    unpack_00_0},

   {unpack_20_fwd_1_3, unpack_20_fwd_2_4, unpack_20_rev_2_4, unpack_20_rev_1_3,
    unpack_20_fwd_3_5, unpack_20_fwd_4_6, unpack_20_rev_4_6, unpack_20_rev_3_5,
    unpack_20_fwd_5_7, unpack_20_fwd_6_8, unpack_20_rev_6_8, unpack_20_rev_5_7, 
    unpack_20_fwd_7_1, unpack_20_fwd_8_2, unpack_20_rev_8_2, unpack_20_rev_7_1,
    unpack_00_0},

   {unpack_22_fwd_1_3, unpack_22_fwd_2_4, unpack_22_rev_2_4, unpack_22_rev_1_3,
    unpack_22_fwd_3_5, unpack_22_fwd_4_6, unpack_22_rev_4_6, unpack_22_rev_3_5,
    unpack_22_fwd_5_7, unpack_22_fwd_6_8, unpack_22_rev_6_8, unpack_22_rev_5_7, 
    unpack_22_fwd_7_1, unpack_22_fwd_8_2, unpack_22_rev_8_2, unpack_22_rev_7_1,
    unpack_00_0},

   {unpack_24_fwd_1_3, unpack_24_fwd_2_4, unpack_24_rev_2_4, unpack_24_rev_1_3,
    unpack_24_fwd_3_5, unpack_24_fwd_4_6, unpack_24_rev_4_6, unpack_24_rev_3_5,
    unpack_24_fwd_5_7, unpack_24_fwd_6_8, unpack_24_rev_6_8, unpack_24_rev_5_7, 
    unpack_24_fwd_7_1, unpack_24_fwd_8_2, unpack_24_rev_8_2, unpack_24_rev_7_1,
    unpack_00_0},

   {unpack_26_fwd_1_3, unpack_26_fwd_2_4, unpack_26_rev_2_4, unpack_26_rev_1_3,
    unpack_26_fwd_3_5, unpack_26_fwd_4_6, unpack_26_rev_4_6, unpack_26_rev_3_5,
    unpack_26_fwd_5_7, unpack_26_fwd_6_8, unpack_26_rev_6_8, unpack_26_rev_5_7, 
    unpack_26_fwd_7_1, unpack_26_fwd_8_2, unpack_26_rev_8_2, unpack_26_rev_7_1,
    unpack_00_0},

   {unpack_28_fwd_1_3, unpack_28_fwd_2_4, unpack_28_rev_2_4, unpack_28_rev_1_3,
    unpack_28_fwd_3_5, unpack_28_fwd_4_6, unpack_28_rev_4_6, unpack_28_rev_3_5,
    unpack_28_fwd_5_7, unpack_28_fwd_6_8, unpack_28_rev_6_8, unpack_28_rev_5_7, 
    unpack_28_fwd_7_1, unpack_28_fwd_8_2, unpack_28_rev_8_2, unpack_28_rev_7_1,
    unpack_00_0},

   {unpack_30_fwd_1_3, unpack_30_fwd_2_4, unpack_30_rev_2_4, unpack_30_rev_1_3,
    unpack_30_fwd_3_5, unpack_30_fwd_4_6, unpack_30_rev_4_6, unpack_30_rev_3_5,
    unpack_30_fwd_5_7, unpack_30_fwd_6_8, unpack_30_rev_6_8, unpack_30_rev_5_7, 
    unpack_30_fwd_7_1, unpack_30_fwd_8_2, unpack_30_rev_8_2, unpack_30_rev_7_1,
    unpack_00_0},

   {unpack_32_fwd_1_3, unpack_32_fwd_2_4, unpack_32_rev_2_4, unpack_32_rev_1_3,
    unpack_32_fwd_3_5, unpack_32_fwd_4_6, unpack_32_rev_4_6, unpack_32_rev_3_5,
    unpack_32_fwd_5_7, unpack_32_fwd_6_8, unpack_32_rev_6_8, unpack_32_rev_5_7, 
    unpack_32_fwd_7_1, unpack_32_fwd_8_2, unpack_32_rev_8_2, unpack_32_rev_7_1,
    unpack_00_0},

};
#endif


#define METAINFO_SIZE 2

#define get_column(s) (s) & 3 /* Not s % 4, which fails on negative values */
#define get_row(s) (s) >> 2 /* Not s / 4, which fails on negative values */


#ifndef LARGE_GENOMES
/* bitpackpages: A list of b-mers (12-mers by default), ending with -1U */
UINT4
Bitpack64_read_two (UINT4 *end0, Storedoligomer_T oligo, UINT4 *bitpackptrs, UINT4 *bitpackcomp) {
  Storedoligomer_T bmer;
  UINT4 *info, nwritten, packsize_div2;
  int remainder0, remainder1, column;
#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
  UINT4 offset0, offset1;
  UINT4 ptr;
  int remainder, row, k, i;
  UINT4 diffs[BLOCKSIZE+1], *bitpack;
#else
  __m128i diffs[4];  /* Need to provide space for 8 rows (or 2 128-bit registers) for ptr and for end0 */
  int delta, row0, row1;
#ifdef BRANCH_FREE_QTR_BLOCK
  UINT4 psums[5];		/* Need 5 to handle case where remainder == 64 */
#endif
  __m128i *bitpack;
  UINT4 *_diffs;
#endif
#ifdef DEBUG
  UINT4 offsets[BLOCKSIZE+1];
#endif


  bmer = oligo/BLOCKSIZE;
  info = &(bitpackptrs[bmer * METAINFO_SIZE]);

  debug(printf("Entered Bitpack64_read_two with oligo %u => bmer %u\n",oligo,bmer));

#ifdef WORDS_BIGENDIAN
  nwritten = Bigendian_convert_uint(info[0]);		/* In 128-bit registers */
  bitpack = (UINT4 *) &(bitpackcomp[nwritten*4]);
  packsize_div2 = (Bigendian_convert_uint(info[METAINFO_SIZE]) - nwritten);

#elif !defined(HAVE_SSE2)
  nwritten = info[0];		/* In 128-bit registers */
  bitpack = (UINT4 *) &(bitpackcomp[nwritten*4]);
  packsize_div2 = (info[METAINFO_SIZE] - nwritten);

#else
  nwritten = info[0];		/* In 128-bit registers */
  bitpack = (__m128i *) &(bitpackcomp[nwritten*4]);
  /* packsize = (info[METAINFO_SIZE] - nwritten)*2; */
  packsize_div2 = (info[METAINFO_SIZE] - nwritten);
#endif

  remainder0 = oligo % BLOCKSIZE;
  remainder1 = remainder0 + 1;

  debug(printf("nwritten %u, packsize %d\n",nwritten,packsize_div2 * 2));
  debug(Bitpack64_block_offsets(offsets,oligo,bitpackptrs,bitpackcomp));

#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
#ifdef WORDS_BIGENDIAN
  offset0 = Bigendian_convert_uint(info[1]);
  offset1 = Bigendian_convert_uint(info[METAINFO_SIZE+1]);
#else
  offset0 = info[1];
  offset1 = info[METAINFO_SIZE+1];
#endif

  /* Unpack all 64 diffs for non-SIMD */
  (unpacker_all_table[packsize_div2*2])(&(diffs[1]),bitpack);

#ifdef DEBUG
#ifdef WORDS_BIGENDIAN
  printf("oligo: %08X, remainder %d, offset0 %u, offset1 %u\n",
         oligo,oligo % BLOCKSIZE,Bigendian_convert_uint(info[1]),Bigendian_convert_uint(info[METAINFO_SIZE+1]));
#else
  printf("oligo: %08X, remainder %d, offset0 %u, offset1 %u\n",
	 oligo,oligo % BLOCKSIZE,info[1],info[METAINFO_SIZE+1]);
#endif
  printf("bitpack:\n");


  for (i = 1; i <= BLOCKSIZE; i++) {
    printf("%d ",diffs[i]);
    if (i % (BLOCKSIZE/4) == 0) {
      printf("\n");
    } else if (i % (BLOCKSIZE/8) == 0) {
      printf("| ");
    }
  }
  printf("\n");
  printf("end of diffs\n");
#endif  

  if ((remainder = oligo % BLOCKSIZE) == 0) {
    ptr = offset0;

  } else if (remainder <= 16) {
    ptr = offset0;

    column = (remainder - 1) % 4; /* Goes from 0 to 3 */
    row = (remainder - 1) / 4;
    debug(printf("column %d, row %d\n",column,row));
    
    for (k = column*2 + 1, i = 0; i <= row; k += BLOCKSIZE/4, i++) {
      debug(printf("Adding diffs[%d] = %u\n",k,diffs[k]));
      ptr += diffs[k];
    }

  } else if (remainder <= 32) {
    ptr = offset0;

    column = (remainder - 1) % 4; /* Goes from 0 to 3 */
    row = (remainder - 1) / 4;
    debug(printf("column %d, row %d\n",column,row));
    
    for (k = column*2 + 1, i = 0; i < 4; k += BLOCKSIZE/4, i++) {
      debug(printf("Adding diffs[%d] = %u\n",k,diffs[k]));
      ptr += diffs[k];
    }

    for (k = column*2 + 2; i <= row; k += BLOCKSIZE/4, i++) {
      debug(printf("Adding diffs[%d] = %u\n",k,diffs[k]));
      ptr += diffs[k];
    }

  } else if (remainder <= 48) {
    ptr = offset1;

    column = (63 - remainder) % 4; /* Goes from 0 to 3.  Assert remainder < 64 */
    row = (63 - remainder) / 4;
    debug(printf("column %d, row %d\n",column,row));

    for (k = column*2 + 9, i = 0; i < 4; k += BLOCKSIZE/4, i++) {
      debug(printf("Subtracting diffs[%d] = %u\n",k,diffs[k]));
      ptr -= diffs[k];
    }

    for (k = column*2 + 10; i <= row; k += BLOCKSIZE/4, i++) {
      debug(printf("Subtracting diffs[%d] = %u\n",k,diffs[k]));
      ptr -= diffs[k];
    }

  } else {
    ptr = offset1;

    column = (63 - remainder) % 4; /* Goes from 0 to 3.  Assert remainder < 64 */
    row = (63 - remainder) / 4;
    debug(printf("column %d, row %d\n",column,row));

    for (k = column*2 + 9, i = 0; i <= row; k += BLOCKSIZE/4, i++) {
      debug(printf("Subtracting diffs[%d] = %u\n",k,diffs[k]));
      ptr -= diffs[k];
    }
  }

  remainder++;
  if (remainder == 64) {
    *end0 = offset1;

  } else if (remainder <= 16) {
    /* Compute necessary cumulative sums */
    *end0 = offset0;

    /* Add 1 for start at diffs[1], and 1 to leave the first element intact */
    diffs[0] = 0;
    column = (remainder - 1) % 4; /* Goes from 0 to 3 */
    row = (remainder - 1) / 4;
    debug(printf("column %d, row %d\n",column,row));
    
    for (k = column*2 + 1, i = 0; i <= row; k += BLOCKSIZE/4, i++) {
      debug(printf("Adding diffs[%d] = %u\n",k,diffs[k]));
      *end0 += diffs[k];
    }

  } else if (remainder <= 32) {
    /* Compute necessary cumulative sums */
    *end0 = offset0;

    /* Add 1 for start at diffs[1], and 1 to leave the first element intact */
    diffs[0] = 0;
    column = (remainder - 1) % 4; /* Goes from 0 to 3 */
    row = (remainder - 1) / 4;
    debug(printf("column %d, row %d\n",column,row));
    
    for (k = column*2 + 1, i = 0; i < 4; k += BLOCKSIZE/4, i++) {
      debug(printf("Adding diffs[%d] = %u\n",k,diffs[k]));
      *end0 += diffs[k];
    }

    for (k = column*2 + 2; i <= row; k += BLOCKSIZE/4, i++) {
      debug(printf("Adding diffs[%d] = %u\n",k,diffs[k]));
      *end0 += diffs[k];
    }

  } else if (remainder <= 48) {
    *end0 = offset1;

    column = (63 - remainder) % 4; /* Goes from 0 to 3.  Assert remainder < 64 */
    row = (63 - remainder) / 4;
    debug(printf("column %d, row %d\n",column,row));

    for (k = column*2 + 9, i = 0; i < 4; k += BLOCKSIZE/4, i++) {
      debug(printf("Subtracting diffs[%d] = %u\n",k,diffs[k]));
      *end0 -= diffs[k];
    }

    for (k = column*2 + 10; i <= row; k += BLOCKSIZE/4, i++) {
      debug(printf("Subtracting diffs[%d] = %u\n",k,diffs[k]));
      *end0 -= diffs[k];
    }

  } else {
    *end0 = offset1;

    column = (63 - remainder) % 4; /* Goes from 0 to 3.  Assert remainder < 64 */
    row = (63 - remainder) / 4;
    debug(printf("column %d, row %d\n",column,row));

    for (k = column*2 + 9, i = 0; i <= row; k += BLOCKSIZE/4, i++) {
      debug(printf("Subtracting diffs[%d] = %u\n",k,diffs[k]));
      *end0 -= diffs[k];
    }
  }

  return ptr;

#else			    /* littleendian and SSE2 */
  _diffs = (UINT4 *) diffs;	/* Assumes a dummy register in diffs[0] */

#ifdef BRANCH_FREE_QTR_BLOCK
  psums[0] = psums[1] = info[1];
  psums[2] = psums[3] = psums[4] = info[METAINFO_SIZE+1];

  delta = 31 - abs(remainder1 - 32);
  column = get_column(delta);
  row = get_row(delta);
  debug(printf("quarter-block %d, delta %d, column %d, row %d\n",quarter_block_1,delta,column,row));
  
  (unpacker_table[packsize_div2][column*4 + quarter_block_1])(diffs,bitpack);
  *end0 = psums[quarter_block_1] + _diffs[row+1] + _diffs[row+2] + _diffs[row+3] + _diffs[row+4];


  delta = 31 - abs(remainder0 - 32);
  column = get_column(delta);
  row = get_row(delta);
  debug(printf("quarter-block %d, delta %d, column %d, row %d\n",quarter_block_0,delta,column,row));

  (unpacker_table[packsize_div2][column*4 + quarter_block_0])(diffs,bitpack);
  return psums[quarter_block_0] + _diffs[row+1] + _diffs[row+2] + _diffs[row+3] + _diffs[row+4];

#else

  if (remainder0 < 16) {
    /* Quarter-block 0 */
    delta = remainder0 - 1;
    column = get_column(delta);
    row0 = get_row(delta);
    row1 = get_row(delta + 1);
    (unpacker_table[packsize_div2][column*4 + 0])(diffs,bitpack);

    _diffs = (UINT4 *) &(diffs[2]);
    assign_sum_fwd(*end0,info[1],_diffs,row1);

    _diffs = (UINT4 *) &(diffs[0]);
    return_sum_fwd(info[1],_diffs,row0);

  } else if (remainder0 < 32) {
    /* Quarter-block 1 */
    delta = remainder0 - 1;
    column = get_column(delta);
    row0 = get_row(delta);
    row1 = get_row(delta + 1);
    (unpacker_table[packsize_div2][column*4 + 1])(diffs,bitpack);

    _diffs = (UINT4 *) &(diffs[2]);
    assign_sum_fwd(*end0,info[1],_diffs,row1);

    _diffs = (UINT4 *) &(diffs[0]);
    return_sum_fwd(info[1],_diffs,row0);

  } else if (remainder0 < 48) {
    /* Quarter-block 2 */
    delta = 63 - remainder1;
    column = get_column(delta);
    row1 = get_row(delta);
    row0 = get_row(delta + 1);
    (unpacker_table[packsize_div2][column*4 + 2])(diffs,bitpack);

    _diffs = (UINT4 *) &(diffs[0]);
    assign_sum_rev(*end0,info[METAINFO_SIZE+1],_diffs,row1);

    _diffs = (UINT4 *) &(diffs[2]);
    return_sum_rev(info[METAINFO_SIZE+1],_diffs,row0);

  } else {
    /* Quarter-block 3 */
    delta = 63 - remainder1;
    column = get_column(delta);
    row1 = get_row(delta);
    row0 = get_row(delta + 1);
    (unpacker_table[packsize_div2][column*4 + 3])(diffs,bitpack);

    _diffs = (UINT4 *) &(diffs[0]);
    assign_sum_rev(*end0,info[METAINFO_SIZE+1],_diffs,row1);

    _diffs = (UINT4 *) &(diffs[2]);
    return_sum_rev(info[METAINFO_SIZE+1],_diffs,row0);
  }

#endif	/* BRANCH_FREE_QTR_BLOCK */
#endif	/* HAVE_SSE2 */

}
#endif


#ifdef LARGE_GENOMES
/* bitpackpages: A list of b-mers (12-mers by default), ending with -1U */
UINT8
Bitpack64_read_two_huge (UINT8 *end0, Storedoligomer_T oligo,
			 UINT4 *bitpackpages, UINT4 *bitpackptrs, UINT4 *bitpackcomp) {
  Storedoligomer_T bmer;
  UINT4 *info, nwritten;
  UINT8 offset0, offset1;
  UINT4 packsize_div2;
  int remainder0, remainder1, column;
#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
  UINT4 ptr;
  int remainder, row, k, i;
  UINT4 diffs[BLOCKSIZE+1], *bitpack;
#else
  int delta, row0, row1;
#ifdef BRANCH_FREE_ROW_SUM
  __m128i diffs[3];
#else
  __m128i diffs[4];  /* Need to provide space for 8 rows (or 2 128-bit registers) for ptr and for end0 */
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
  UINT8 psums[5];		/* Need 5 to handle case where remainder == 64 */
#endif
  __m128i *bitpack;
  UINT4 *_diffs;
#endif
  UINT4 *pageptr;
#ifdef DEBUG
  UINT4 offsets[BLOCKSIZE+1];
#endif


  bmer = oligo/BLOCKSIZE;
  info = &(bitpackptrs[bmer * METAINFO_SIZE]);

  debug(printf("Entered Bitpack64_read_two_huge with oligo %u => bmer %u\n",oligo,bmer));

#ifdef WORDS_BIGENDIAN
  nwritten = Bigendian_convert_uint(info[0]);
  bitpack = (UINT4 *) &(bitpackcomp[nwritten*4]);
  offset0 = (UINT8) Bigendian_convert_uint(info[1]);
  offset1 = (UINT8) Bigendian_convert_uint(info[METAINFO_SIZE+1]);

#elif !defined(HAVE_SSE2)
  nwritten = info[0];
  bitpack = (UINT4 *) &(bitpackcomp[nwritten*4]);
  offset0 = (UINT8) info[1];
  offset1 = (UINT8) info[METAINFO_SIZE+1];

#else
  nwritten = info[0];
  bitpack = (__m128i *) &(bitpackcomp[nwritten*4]);
  offset0 = (UINT8) info[1];
  offset1 = (UINT8) info[METAINFO_SIZE+1];
#endif

  debug(printf("offsets are %llu, %llu\n",offset0,offset1));

#ifdef WORDS_BIGENDIAN
  if (bitpackpages != NULL) {
    pageptr = bitpackpages;
    debug(printf("  compare bmer %u with pageptr %u\n",bmer,*pageptr));
    while (bmer >= Bigendian_convert_uint(*pageptr)) {
      offset0 += POSITIONS_PAGE;
      offset1 += POSITIONS_PAGE;
      pageptr++;
    }

    if (bmer + 1 >= Bigendian_convert_uint(*pageptr)) {
      offset1 += POSITIONS_PAGE;
      /* pageptr++; */
    }
  }
  debug(printf("offsets are %llu, %llu\n",offset0,offset1));
  packsize_div2 = (Bigendian_convert_uint(info[METAINFO_SIZE]) - nwritten);

#else
  if (bitpackpages != NULL) {
    pageptr = bitpackpages;
    debug(printf("  compare bmer %u with pageptr %u\n",bmer,*pageptr));
    while (bmer >= *pageptr) {
      offset0 += POSITIONS_PAGE;
      offset1 += POSITIONS_PAGE;
      pageptr++;
    }

    if (bmer + 1 >= *pageptr) {
      offset1 += POSITIONS_PAGE;
      /* pageptr++; */
    }
  }
  debug(printf("offsets are %llu, %llu\n",offset0,offset1));
  packsize_div2 = (info[METAINFO_SIZE] - nwritten);
#endif

  remainder0 = oligo % BLOCKSIZE;
  remainder1 = remainder0 + 1;

  /* debug(Bitpack64_block_offsets_huge(offsets,oligo,bitpackpages,bitpackptrs,bitpackcomp)); */

#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)

  /* Unpack all 64 diffs for non-SIMD */
  (unpacker_all_table[packsize_div2*2])(&(diffs[1]),bitpack);

#ifdef DEBUG
  printf("oligo: %08X, remainder %d, offset0 %u, offset1 %u\n",
	 oligo,oligo % BLOCKSIZE,info[1],info[METAINFO_SIZE+1]);
  printf("bitpack:\n");

  for (i = 1; i <= BLOCKSIZE; i++) {
    printf("%d ",diffs[i]);
    if (i % (BLOCKSIZE/4) == 0) {
      printf("\n");
    } else if (i % (BLOCKSIZE/8) == 0) {
      printf("| ");
    }
  }
  printf("\n");
  printf("end of diffs\n");
#endif  

  if ((remainder = oligo % BLOCKSIZE) == 0) {
    ptr = offset0;

  } else if (remainder <= 16) {
    ptr = offset0;

    /* Add 1 for start at diffs[1], and 1 to leave the first element intact */
    column = (remainder - 1) % 4; /* Goes from 0 to 3 */
    row = (remainder - 1) / 4;
    debug(printf("column %d, row %d\n",column,row));
    
    for (k = column*2 + 1, i = 0; i <= row; k += BLOCKSIZE/4, i++) {
      debug(printf("Adding diffs[%d] = %u\n",k,diffs[k]));
      ptr += diffs[k];
    }

  } else if (remainder <= 32) {
    ptr = offset0;

    column = (remainder - 1) % 4; /* Goes from 0 to 3 */
    row = (remainder - 1) / 4;
    debug(printf("column %d, row %d\n",column,row));
    
    for (k = column*2 + 1, i = 0; i < 4; k += BLOCKSIZE/4, i++) {
      debug(printf("Adding diffs[%d] = %u\n",k,diffs[k]));
      ptr += diffs[k];
    }

    for (k = column*2 + 2; i <= row; k += BLOCKSIZE/4, i++) {
      debug(printf("Adding diffs[%d] = %u\n",k,diffs[k]));
      ptr += diffs[k];
    }

  } else if (remainder <= 48) {
    ptr = offset1;

    column = (63 - remainder) % 4; /* Goes from 0 to 3.  Assert remainder < 64 */
    row = (63 - remainder) / 4;
    debug(printf("column %d, row %d\n",column,row));

    for (k = column*2 + 9, i = 0; i < 4; k += BLOCKSIZE/4, i++) {
      debug(printf("Subtracting diffs[%d] = %u\n",k,diffs[k]));
      ptr -= diffs[k];
    }

    for (k = column*2 + 10; i <= row; k += BLOCKSIZE/4, i++) {
      debug(printf("Subtracting diffs[%d] = %u\n",k,diffs[k]));
      ptr -= diffs[k];
    }

  } else {
    ptr = offset1;

    column = (63 - remainder) % 4; /* Goes from 0 to 3.  Assert remainder < 64 */
    row = (63 - remainder) / 4;
    debug(printf("column %d, row %d\n",column,row));

    for (k = column*2 + 9, i = 0; i <= row; k += BLOCKSIZE/4, i++) {
      debug(printf("Subtracting diffs[%d] = %u\n",k,diffs[k]));
      ptr -= diffs[k];
    }
  }

  remainder++;
  if (remainder == 64) {
    *end0 = offset1;

  } else if (remainder <= 16) {
    /* Compute necessary cumulative sums */
    *end0 = offset0;

    /* Add 1 for start at diffs[1], and 1 to leave the first element intact */
    diffs[0] = 0;
    column = (remainder - 1) % 4; /* Goes from 0 to 3 */
    row = (remainder - 1) / 4;
    debug(printf("column %d, row %d\n",column,row));
    
    for (k = column*2 + 1, i = 0; i <= row; k += BLOCKSIZE/4, i++) {
      debug(printf("Adding diffs[%d] = %u\n",k,diffs[k]));
      *end0 += diffs[k];
    }

  } else if (remainder <= 32) {
    /* Compute necessary cumulative sums */
    *end0 = offset0;

    /* Add 1 for start at diffs[1], and 1 to leave the first element intact */
    diffs[0] = 0;
    column = (remainder - 1) % 4; /* Goes from 0 to 3 */
    row = (remainder - 1) / 4;
    debug(printf("column %d, row %d\n",column,row));
    
    for (k = column*2 + 1, i = 0; i < 4; k += BLOCKSIZE/4, i++) {
      debug(printf("Adding diffs[%d] = %u\n",k,diffs[k]));
      *end0 += diffs[k];
    }

    for (k = column*2 + 2; i <= row; k += BLOCKSIZE/4, i++) {
      debug(printf("Adding diffs[%d] = %u\n",k,diffs[k]));
      *end0 += diffs[k];
    }

  } else if (remainder <= 48) {
    *end0 = offset1;

    column = (63 - remainder) % 4; /* Goes from 0 to 3.  Assert remainder < 64 */
    row = (63 - remainder) / 4;
    debug(printf("column %d, row %d\n",column,row));

    for (k = column*2 + 9, i = 0; i < 4; k += BLOCKSIZE/4, i++) {
      debug(printf("Subtracting diffs[%d] = %u\n",k,diffs[k]));
      *end0 -= diffs[k];
    }

    for (k = column*2 + 10; i <= row; k += BLOCKSIZE/4, i++) {
      debug(printf("Subtracting diffs[%d] = %u\n",k,diffs[k]));
      *end0 -= diffs[k];
    }

  } else {
    *end0 = offset1;

    column = (63 - remainder) % 4; /* Goes from 0 to 3.  Assert remainder < 64 */
    row = (63 - remainder) / 4;
    debug(printf("column %d, row %d\n",column,row));

    for (k = column*2 + 9, i = 0; i <= row; k += BLOCKSIZE/4, i++) {
      debug(printf("Subtracting diffs[%d] = %u\n",k,diffs[k]));
      *end0 -= diffs[k];
    }
  }

  return ptr;


#else			    /* littleendian and SSE2 */
  _diffs = (UINT4 *) diffs;	/* Assumes a dummy register in diffs[0] */

#ifdef BRANCH_FREE_QTR_BLOCK
  psums[0] = psums[1] = offset0;
  psums[2] = psums[3] = psums[4] = offset1;

  delta = 31 - abs(remainder1 - 32);
  column = get_column(delta);
  row = get_row(delta);
  debug(printf("quarter-block %d, delta %d, column %d, row %d\n",quarter_block_1,delta,column,row));
  
  (unpacker_table[packsize_div2][column*4 + quarter_block_1])(diffs,bitpack);
  *end0 = psums[quarter_block_1] + (INT4) (_diffs[row+1] + _diffs[row+2] + _diffs[row+3] + _diffs[row+4]);


  delta = 31 - abs(remainder0 - 32);
  column = get_column(delta);
  row = get_row(delta);
  debug(printf("quarter-block %d, delta %d, column %d, row %d\n",quarter_block_0,delta,column,row));

  (unpacker_table[packsize_div2][column*4 + quarter_block_0])(diffs,bitpack);
  return psums[quarter_block_0] + (INT4) (_diffs[row+1] + _diffs[row+2] + _diffs[row+3] + _diffs[row+4]);

#else

  if (remainder0 < 16) {
    /* Quarter-block 0 */
    delta = remainder0 - 1;
    column = get_column(delta);
    row0 = get_row(delta);
    row1 = get_row(delta + 1);
    debug(printf("quarter_block 0, remainder %d, delta %d, column %d, row %d\n",remainder0,delta,column,row0));
    (unpacker_table[packsize_div2][column*4 + 0])(diffs,bitpack);

    _diffs = (UINT4 *) &(diffs[2]);
    assign_sum_fwd(*end0,offset0,_diffs,row1);

    _diffs = (UINT4 *) &(diffs[0]);
    return_sum_fwd(offset0,_diffs,row0);

  } else if (remainder0 < 32) {
    /* Quarter-block 1 */
    delta = remainder0 - 1;
    column = get_column(delta);
    row0 = get_row(delta);
    row1 = get_row(delta + 1);
    debug(printf("quarter_block 1, remainder %d, delta %d, column %d, row %d\n",remainder0,delta,column,row0));
    (unpacker_table[packsize_div2][column*4 + 1])(diffs,bitpack);

    _diffs = (UINT4 *) &(diffs[2]);
    assign_sum_fwd(*end0,offset0,_diffs,row1);

    _diffs = (UINT4 *) &(diffs[0]);
    return_sum_fwd(offset0,_diffs,row0);

  } else if (remainder0 < 48) {
    /* Quarter-block 2 */
    delta = 63 - remainder1;
    column = get_column(delta);
    row1 = get_row(delta);
    row0 = get_row(delta + 1);
    debug(printf("quarter_block 2, remainder %d, delta %d, column %d, row %d\n",remainder0,delta,column,row0));
    (unpacker_table[packsize_div2][column*4 + 2])(diffs,bitpack);

    _diffs = (UINT4 *) &(diffs[0]);
    assign_sum_rev(*end0,offset1,_diffs,row1);

    _diffs = (UINT4 *) &(diffs[2]);
    return_sum_rev(offset1,_diffs,row0);

  } else {
    /* Quarter-block 3 */
    delta = 63 - remainder1;
    column = get_column(delta);
    row1 = get_row(delta);
    row0 = get_row(delta + 1);
    debug(printf("quarter_block 3, remainder %d, delta %d, column %d, row %d\n",remainder0,delta,column,row0));
    (unpacker_table[packsize_div2][column*4 + 3])(diffs,bitpack);

    _diffs = (UINT4 *) &(diffs[0]);
    assign_sum_rev(*end0,offset1,_diffs,row1);

    _diffs = (UINT4 *) &(diffs[2]);
    return_sum_rev(offset1,_diffs,row0);
  }

#endif	/* BRANCH_FREE_QTR_BLOCK */
#endif  /* HAVE_SSE2 */
}
#endif
