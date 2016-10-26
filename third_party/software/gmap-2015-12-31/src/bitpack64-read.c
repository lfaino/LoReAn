static char rcsid[] = "$Id: bitpack64-read.c 168395 2015-06-26 17:13:13Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "bitpack64-read.h"

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
   in Bitpack64_offsetptr_huge.  Would therefore recommend it be turned off.
*/


/* #define BRANCH_FREE_ROW_SUM 1 */
/* #define BRANCH_FREE_QTR_BLOCK 1 */

#ifdef HAVE_SSE2
#ifdef DEBUG
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


#if 0
#ifdef HAVE_SSE2
#ifdef ALLOW_ODD_PACKSIZES
static __m128i mask1, mask2, mask3, mask4, mask5, mask6, mask7, mask8,
  mask9, mask10, mask11, mask12, mask13, mask14, mask15, mask16,
  mask17, mask18, mask19, mask20, mask21, mask22, mask23, mask24,
  mask25, mask26, mask27, mask28, mask29, mask30, mask31;
#else
static __m128i mask2, mask4, mask6, mask8, mask10, mask12, mask14, mask16,
  mask18, mask20, mask22, mask24, mask26, mask28, mask30;
#endif
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
#ifdef BRANCH_FREE_ROW_SUM
  __m128i zero = _mm_set1_epi32(0U);
  _mm_store_si128(out++, zero); /* dummy */
#endif

  /* _mm_store_si128(out++, zero); -- Not needed, since row == -1 */

  return;
}


static void
unpack_00_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i zero = _mm_set1_epi32(0U);

#ifdef BRANCH_FREE_ROW_SUM
  _mm_store_si128(out++, zero); /* dummy */
#endif

  _mm_store_si128(out++, zero);

  return;
}

static void
unpack_00_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i zero = _mm_set1_epi32(0U);

#ifdef BRANCH_FREE_ROW_SUM
  out++;			/* dummy */
#endif

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
unpack_02_fwd_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask2 = _mm_set1_epi32(3U);


#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask2);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_02_fwd_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask2 = _mm_set1_epi32(3U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask2);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask2);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_02_fwd_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask2 = _mm_set1_epi32(3U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask2);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_02_fwd_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask2 = _mm_set1_epi32(3U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask2);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask2);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}


static void
unpack_02_fwd_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask2 = _mm_set1_epi32(3U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask2);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_02_fwd_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask2 = _mm_set1_epi32(3U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask2);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,10) , mask2);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_02_fwd_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask2 = _mm_set1_epi32(3U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask2);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_02_fwd_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask2 = _mm_set1_epi32(3U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask2);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,14) , mask2);
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
unpack_02_rev_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask2 = _mm_set1_epi32(3U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask2);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_02_rev_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask2 = _mm_set1_epi32(3U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

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

    return;
}

static void
unpack_02_rev_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask2 = _mm_set1_epi32(3U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,20) , mask2);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_02_rev_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask2 = _mm_set1_epi32(3U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

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

    return;
}

static void
unpack_02_rev_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask2 = _mm_set1_epi32(3U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,24) , mask2);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_02_rev_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask2 = _mm_set1_epi32(3U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

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

    return;
}

static void
unpack_02_rev_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask2 = _mm_set1_epi32(3U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,28) , mask2);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_02_rev_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask2 = _mm_set1_epi32(3U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask1), 3-1));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 3 - 1), mask3));
#endif
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
unpack_04_fwd_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask4 = _mm_set1_epi32(15U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask4);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_04_fwd_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask4 = _mm_set1_epi32(15U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask4);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask4);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_04_fwd_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask4 = _mm_set1_epi32(15U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask4);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_04_fwd_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask4 = _mm_set1_epi32(15U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask4);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask4);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_04_fwd_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask4 = _mm_set1_epi32(15U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask4);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_04_fwd_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask4 = _mm_set1_epi32(15U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask4);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,20) , mask4);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_04_fwd_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask4 = _mm_set1_epi32(15U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,24) , mask4);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_04_fwd_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask4 = _mm_set1_epi32(15U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,24) , mask4);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
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
unpack_04_rev_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask4 = _mm_set1_epi32(15U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(++in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask4);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_04_rev_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask4 = _mm_set1_epi32(15U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

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

    return;
}

static void
unpack_04_rev_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask4 = _mm_set1_epi32(15U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(++in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask4);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_04_rev_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask4 = _mm_set1_epi32(15U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

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

    return;
}

static void
unpack_04_rev_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask4 = _mm_set1_epi32(15U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(++in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask4);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_04_rev_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask4 = _mm_set1_epi32(15U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

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

    return;
}

static void
unpack_04_rev_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask4 = _mm_set1_epi32(15U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(++in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,24) , mask4);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_04_rev_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask4 = _mm_set1_epi32(15U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask3), 5-3));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 5 - 3), mask5));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask1), 5-1));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 5 - 1), mask5));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 6-4));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 6 - 4), mask6));
#endif
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
unpack_06_fwd_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i total;
    const __m128i mask6 =  _mm_set1_epi32(63U);

#ifdef BRANCH_FREE_ROW_SUM
  __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask6);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_06_fwd_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask6 =  _mm_set1_epi32(63U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask6);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask6);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_06_fwd_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask6 =  _mm_set1_epi32(63U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask6);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_06_fwd_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask6 =  _mm_set1_epi32(63U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask6);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,18) , mask6);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_06_fwd_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask6 =  _mm_set1_epi32(63U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,24) , mask6);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_06_fwd_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask6 =  _mm_set1_epi32(63U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,24) , mask6);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 6-4));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 6 - 4), mask6));
#endif
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_06_fwd_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask6 =  _mm_set1_epi32(63U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(++in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask6);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_06_fwd_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask6 =  _mm_set1_epi32(63U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask2), 6-2));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 6 - 2), mask6));
#endif
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
unpack_06_rev_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask6 =  _mm_set1_epi32(63U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(++in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask6);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_06_rev_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask6 =  _mm_set1_epi32(63U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

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

    return;
}

static void
unpack_06_rev_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask6 =  _mm_set1_epi32(63U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(++in);


    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask2), 6-2));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 6 - 2), mask6));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_06_rev_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask6 =  _mm_set1_epi32(63U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    InReg = _mm_load_si128(++in);


    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask2), 6-2));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 6 - 2), mask6));
#endif
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
unpack_06_rev_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask6 =  _mm_set1_epi32(63U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 2;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask6);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_06_rev_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask6 =  _mm_set1_epi32(63U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

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

    return;

}

static void
unpack_06_rev_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask6 =  _mm_set1_epi32(63U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 2;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,20) , mask6);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_06_rev_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask6 =  _mm_set1_epi32(63U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 2;
    InReg = _mm_load_si128(in);


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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask3), 7-3));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 7 - 3), mask7));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask6), 7-6));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 7 - 6), mask7));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask2), 7-2));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 7 - 2), mask7));
#endif
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
unpack_08_fwd_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask8 =  _mm_set1_epi32(255U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask8);
    _mm_store_si128(out++, total);

    return;
}


static void
unpack_08_fwd_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask8 =  _mm_set1_epi32(255U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask8);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask8);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_08_fwd_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask8 =  _mm_set1_epi32(255U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask8);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_08_fwd_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask8 =  _mm_set1_epi32(255U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask8);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    /* InReg = _mm_load_si128(++in); */

    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_08_fwd_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask8 =  _mm_set1_epi32(255U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(++in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask8);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_08_fwd_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask8 =  _mm_set1_epi32(255U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    InReg = _mm_load_si128(++in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask8);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask8);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_08_fwd_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask8 =  _mm_set1_epi32(255U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(++in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask8);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_08_fwd_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask8 =  _mm_set1_epi32(255U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

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
unpack_08_rev_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask8 =  _mm_set1_epi32(255U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 2;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask8);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_08_rev_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask8 =  _mm_set1_epi32(255U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

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

    return;
}

static void
unpack_08_rev_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask8 =  _mm_set1_epi32(255U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 2;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask8);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_08_rev_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask8 =  _mm_set1_epi32(255U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

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

    return;
}

static void
unpack_08_rev_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask8 =  _mm_set1_epi32(255U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 3;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask8);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_08_rev_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask8 =  _mm_set1_epi32(255U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

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

    return;
}

static void
unpack_08_rev_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask8 =  _mm_set1_epi32(255U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 3;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask8);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_08_rev_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask8 =  _mm_set1_epi32(255U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 3;
    InReg = _mm_load_si128(in);


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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 9-4));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 9 - 4), mask9));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 9-8));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 9 - 8), mask9));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask3), 9-3));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 9 - 3), mask9));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask7), 9-7));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 9 - 7), mask9));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 10-8));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 10 - 8), mask10));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask6), 10-6));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 10 - 6), mask10));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask10);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_10_fwd_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask10 =  _mm_set1_epi32(1023U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask10);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_10_fwd_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask10 =  _mm_set1_epi32(1023U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask10);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,10) , mask10);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_10_fwd_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask10 =  _mm_set1_epi32(1023U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,20) , mask10);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_10_fwd_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask10 =  _mm_set1_epi32(1023U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,20) , mask10);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 10-8));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 10 - 8), mask10));
#endif
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_10_fwd_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask10 =  _mm_set1_epi32(1023U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(++in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask10);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_10_fwd_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask10 =  _mm_set1_epi32(1023U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    InReg = _mm_load_si128(++in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask10);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,18) , mask10);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_10_fwd_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask10 =  _mm_set1_epi32(1023U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(++in);


    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask6), 10-6));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 10 - 6), mask10));
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_10_fwd_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask10 =  _mm_set1_epi32(1023U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    InReg = _mm_load_si128(++in);


    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask6), 10-6));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 10 - 6), mask10));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 10-4));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 10 - 4), mask10));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask2), 10-2));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 10 - 2), mask10));
#endif
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
unpack_10_rev_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask10 =  _mm_set1_epi32(1023U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 2;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask10);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_10_rev_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask10 =  _mm_set1_epi32(1023U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 2;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask10);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 10-4));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 10 - 4), mask10));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_10_rev_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask10 =  _mm_set1_epi32(1023U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 3;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask10);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_10_rev_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask10 =  _mm_set1_epi32(1023U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

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

    return;
}

static void
unpack_10_rev_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask10 =  _mm_set1_epi32(1023U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 3;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask2), 10-2));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 10 - 2), mask10));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_10_rev_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask10 =  _mm_set1_epi32(1023U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 3;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask2), 10-2));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 10 - 2), mask10));
#endif
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
unpack_10_rev_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask10 =  _mm_set1_epi32(1023U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 4;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask10);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_10_rev_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask10 =  _mm_set1_epi32(1023U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 4;
    InReg = _mm_load_si128(in);


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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask1), 11-1));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 11 - 1), mask11));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask2), 11-2));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 11 - 2), mask11));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask3), 11-3));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 11 - 3), mask11));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 11-4));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 11 - 4), mask11));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask5), 11-5));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 11 - 5), mask11));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 12-4));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 12 - 4), mask12));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 12-8));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 12 - 8), mask12));
#endif
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
unpack_12_fwd_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask12 =  _mm_set1_epi32(4095U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask12);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_12_fwd_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask12 =  _mm_set1_epi32(4095U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask12);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask12);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_12_fwd_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask12 =  _mm_set1_epi32(4095U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 12-4));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 12 - 4), mask12));
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_12_fwd_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask12 =  _mm_set1_epi32(4095U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 12-4));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 12 - 4), mask12));
#endif
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask12);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_12_fwd_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask12 =  _mm_set1_epi32(4095U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(++in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask12);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_12_fwd_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask12 =  _mm_set1_epi32(4095U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    InReg = _mm_load_si128(++in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask12);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 12-8));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 12 - 8), mask12));
#endif
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_12_fwd_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask12 =  _mm_set1_epi32(4095U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 2;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask12);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_12_fwd_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask12 =  _mm_set1_epi32(4095U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 12-4));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 12 - 4), mask12));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 12-8));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 12 - 8), mask12));
#endif
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
unpack_12_rev_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask12 =  _mm_set1_epi32(4095U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 3;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask12);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_12_rev_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask12 =  _mm_set1_epi32(4095U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

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

    return;
}

static void
unpack_12_rev_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask12 =  _mm_set1_epi32(4095U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 3;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 12-4));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 12 - 4), mask12));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_12_rev_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask12 =  _mm_set1_epi32(4095U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 3;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 12-4));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 12 - 4), mask12));
#endif
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
unpack_12_rev_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask12 =  _mm_set1_epi32(4095U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 4;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask12);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_12_rev_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask12 =  _mm_set1_epi32(4095U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 4;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask12);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 12-8));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 12 - 8), mask12));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_12_rev_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask12 =  _mm_set1_epi32(4095U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 5;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask12);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_12_rev_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask12 =  _mm_set1_epi32(4095U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 5;
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask7), 13-7));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 13 - 7), mask13));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,7) , mask13);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask1), 13-1));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 13 - 1), mask13));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 13-8));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 13 - 8), mask13));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask13);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,21) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask2), 13-2));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 13 - 2), mask13));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask9), 13-9));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 13 - 9), mask13));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,9) , mask13);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask3), 13-3));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 13 - 3), mask13));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask10), 14-10));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 14 - 10), mask14));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,10) , mask14);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask6), 14-6));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 14 - 6), mask14));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask14);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask2), 14-2));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 14 - 2), mask14));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask14);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_14_fwd_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask14 =  _mm_set1_epi32(16383U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask14);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_14_fwd_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask14 =  _mm_set1_epi32(16383U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask14);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,14) , mask14);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_14_fwd_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask14 =  _mm_set1_epi32(16383U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask10), 14-10));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 14 - 10), mask14));
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_14_fwd_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask14 =  _mm_set1_epi32(16383U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask10), 14-10));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 14 - 10), mask14));
#endif
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,10) , mask14);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_14_fwd_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask14 =  _mm_set1_epi32(16383U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(++in);


    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask6), 14-6));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 14 - 6), mask14));
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_14_fwd_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask14 =  _mm_set1_epi32(16383U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    InReg = _mm_load_si128(++in);


    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask6), 14-6));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 14 - 6), mask14));
#endif
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask14);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_14_fwd_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask14 =  _mm_set1_epi32(16383U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 2;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask2), 14-2));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 14 - 2), mask14));
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_14_fwd_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask14 =  _mm_set1_epi32(16383U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 2;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask2), 14-2));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 14 - 2), mask14));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask12), 14-12));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 14 - 12), mask14));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask14);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 14-8));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 14 - 8), mask14));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask14);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 14-4));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 14 - 4), mask14));
#endif
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
unpack_14_rev_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask14 =  _mm_set1_epi32(16383U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 3;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask14);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_14_rev_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask14 =  _mm_set1_epi32(16383U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 3;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask14);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask12), 14-12));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 14 - 12), mask14));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_14_rev_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask14 =  _mm_set1_epi32(16383U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 4;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask14);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_14_rev_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask14 =  _mm_set1_epi32(16383U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 4;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask14);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 14-8));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 14 - 8), mask14));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_14_rev_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask14 =  _mm_set1_epi32(16383U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 5;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask14);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_14_rev_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask14 =  _mm_set1_epi32(16383U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 5;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask14);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 14-4));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 14 - 4), mask14));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_14_rev_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask14 =  _mm_set1_epi32(16383U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 6;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask14);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_14_rev_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask14 =  _mm_set1_epi32(16383U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 6;
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask13), 15-13));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 15 - 13), mask15));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,13) , mask15);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask11), 15-11));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 15 - 11), mask15));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,11) , mask15);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask9), 15-9));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 15 - 9), mask15));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,9) , mask15);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask7), 15-7));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 15 - 7), mask15));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,7) , mask15);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask5), 15-5));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 15 - 5), mask15));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,5) , mask15);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask3), 15-3));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 15 - 3), mask15));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,3) , mask15);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,18) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask1), 15-1));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 15 - 1), mask15));
#endif
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
unpack_16_fwd_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask16 =  _mm_set1_epi32(65535U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_16_fwd_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask16 =  _mm_set1_epi32(65535U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    /* InReg = _mm_load_si128(++in); */

    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_16_fwd_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask16 =  _mm_set1_epi32(65535U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(++in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_16_fwd_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask16 =  _mm_set1_epi32(65535U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

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
unpack_16_fwd_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask16 =  _mm_set1_epi32(65535U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 2;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_16_fwd_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask16 =  _mm_set1_epi32(65535U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 2;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    /* InReg = _mm_load_si128(++in); */

    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_16_fwd_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask16 =  _mm_set1_epi32(65535U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 3;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_16_fwd_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask16 =  _mm_set1_epi32(65535U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

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
unpack_16_rev_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask16 =  _mm_set1_epi32(65535U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 4;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_16_rev_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask16 =  _mm_set1_epi32(65535U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

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

    return;
}

static void
unpack_16_rev_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask16 =  _mm_set1_epi32(65535U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 5;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_16_rev_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask16 =  _mm_set1_epi32(65535U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

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

    return;
}

static void
unpack_16_rev_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask16 =  _mm_set1_epi32(65535U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 6;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_16_rev_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask16 =  _mm_set1_epi32(65535U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

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

    return;
}

static void
unpack_16_rev_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask16 =  _mm_set1_epi32(65535U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 7;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_16_rev_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask16 =  _mm_set1_epi32(65535U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 7;
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask2), 17-2));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 17 - 2), mask17));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask17);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,19) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 17-4));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 17 - 4), mask17));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask17);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,21) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask6), 17-6));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 17 - 6), mask17));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask17);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,23) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 17-8));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 17 - 8), mask17));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask17);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,25) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask10), 17-10));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 17 - 10), mask17));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,10) , mask17);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,27) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask12), 17-12));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 17 - 12), mask17));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask17);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,29) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask14), 17-14));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 17 - 14), mask17));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,14) , mask17);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,31) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 17-16));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 17 - 16), mask17));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 18-4));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 4), mask18));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask18);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 18-8));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 8), mask18));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask18);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask12), 18-12));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 12), mask18));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask18);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 18-16));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 16), mask18));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_18_fwd_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask18 =  _mm_set1_epi32(262143U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask18);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_18_fwd_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask18 =  _mm_set1_epi32(262143U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask18);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,18) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 18-4));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 4), mask18));
#endif
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_18_fwd_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask18 =  _mm_set1_epi32(262143U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(++in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask18);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_18_fwd_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask18 =  _mm_set1_epi32(262143U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    InReg = _mm_load_si128(++in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask18);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 18-8));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 8), mask18));
#endif
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_18_fwd_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask18 =  _mm_set1_epi32(262143U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 2;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask18);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_18_fwd_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask18 =  _mm_set1_epi32(262143U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 2;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask18);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask12), 18-12));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 12), mask18));
#endif
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_18_fwd_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask18 =  _mm_set1_epi32(262143U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 3;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask18);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_18_fwd_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask18 =  _mm_set1_epi32(262143U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 3;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask18);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 18-16));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 16), mask18));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask2), 18-2));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 2), mask18));
#endif
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask18);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask6), 18-6));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 6), mask18));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask18);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask10), 18-10));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 10), mask18));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,10) , mask18);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask14), 18-14));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 14), mask18));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,14) ;
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_18_rev_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask18 =  _mm_set1_epi32(262143U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 4;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask2), 18-2));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 2), mask18));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_18_rev_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask18 =  _mm_set1_epi32(262143U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 4;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask2), 18-2));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 2), mask18));
#endif
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

    return;
}

static void
unpack_18_rev_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask18 =  _mm_set1_epi32(262143U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 5;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask6), 18-6));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 6), mask18));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_18_rev_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask18 =  _mm_set1_epi32(262143U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 5;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask6), 18-6));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 6), mask18));
#endif
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
unpack_18_rev_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask18 =  _mm_set1_epi32(262143U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 6;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask10), 18-10));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 10), mask18));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_18_rev_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask18 =  _mm_set1_epi32(262143U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 6;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask10), 18-10));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 10), mask18));
#endif
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
unpack_18_rev_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask18 =  _mm_set1_epi32(262143U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 7;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask14), 18-14));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 14), mask18));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_18_rev_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask18 =  _mm_set1_epi32(262143U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 7;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask14), 18-14));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 18 - 14), mask18));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask6), 19-6));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 19 - 6), mask19));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask19);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,25) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask12), 19-12));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 19 - 12), mask19));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask19);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,31) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask18), 19-18));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 19 - 18), mask19));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,18) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask5), 19-5));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 19 - 5), mask19));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,5) , mask19);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask11), 19-11));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 19 - 11), mask19));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,11) , mask19);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask17), 19-17));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 19 - 17), mask19));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,17) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 19-4));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 19 - 4), mask19));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask19);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,23) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask10), 19-10));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 19 - 10), mask19));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,10) , mask19);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,29) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 19-16));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 19 - 16), mask19));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 20-8));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 8), mask20));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask20);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 20-16));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 16), mask20));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 20-4));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 4), mask20));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask20);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask12), 20-12));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 12), mask20));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_20_fwd_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask20 =  _mm_set1_epi32(1048575U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask20);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_20_fwd_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask20 =  _mm_set1_epi32(1048575U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask20);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 20-8));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 8), mask20));
#endif
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_20_fwd_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask20 =  _mm_set1_epi32(1048575U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(++in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask20);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_20_fwd_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask20 =  _mm_set1_epi32(1048575U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    InReg = _mm_load_si128(++in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask20);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 20-16));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 16), mask20));
#endif
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;

}

static void
unpack_20_fwd_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask20 =  _mm_set1_epi32(1048575U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 2;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 20-4));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 4), mask20));
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_20_fwd_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask20 =  _mm_set1_epi32(1048575U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 2;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 20-4));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 4), mask20));
#endif
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask20);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_20_fwd_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask20 =  _mm_set1_epi32(1048575U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 3;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask12), 20-12));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 12), mask20));
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_20_fwd_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask20 =  _mm_set1_epi32(1048575U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 3;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask12), 20-12));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 12), mask20));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 20-8));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 8), mask20));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask20);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 20-16));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 16), mask20));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 20-4));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 4), mask20));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask20);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask12), 20-12));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 12), mask20));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_20_rev_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask20 =  _mm_set1_epi32(1048575U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 5;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask20);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_20_rev_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask20 =  _mm_set1_epi32(1048575U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 5;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask20);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 20-8));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 8), mask20));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
} 

static void
unpack_20_rev_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask20 =  _mm_set1_epi32(1048575U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 6;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask20);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_20_rev_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask20 =  _mm_set1_epi32(1048575U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 6;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask20);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 20-16));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 16), mask20));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;

} 

static void
unpack_20_rev_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask20 =  _mm_set1_epi32(1048575U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 7;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 20-4));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 4), mask20));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_20_rev_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask20 =  _mm_set1_epi32(1048575U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 7;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 20-4));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 4), mask20));
#endif
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
unpack_20_rev_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask20 =  _mm_set1_epi32(1048575U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 8;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask12), 20-12));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 12), mask20));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_20_rev_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask20 =  _mm_set1_epi32(1048575U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 8;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask12), 20-12));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 20 - 12), mask20));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask10), 21-10));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 21 - 10), mask21));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,10) , mask21);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,31) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask20), 21-20));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 21 - 20), mask21));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask9), 21-9));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 21 - 9), mask21));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,9) , mask21);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask19), 21-19));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 21 - 19), mask21));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,19) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 21-8));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 21 - 8), mask21));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask21);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,29) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask18), 21-18));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 21 - 18), mask21));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,18) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask7), 21-7));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 21 - 7), mask21));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,7) , mask21);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask17), 21-17));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 21 - 17), mask21));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,17) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask6), 21-6));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 21 - 6), mask21));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask21);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,27) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 21-16));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 21 - 16), mask21));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask12), 22-12));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 12), mask22));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask2), 22-2));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 2), mask22));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask22);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask14), 22-14));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 14), mask22));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,14) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 22-4));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 4), mask22));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask22);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 22-16));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 16), mask22));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_22_fwd_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask22 =  _mm_set1_epi32(4194303U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask22);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_22_fwd_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask22 =  _mm_set1_epi32(4194303U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask22);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask12), 22-12));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 12), mask22));
#endif
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_22_fwd_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask22 =  _mm_set1_epi32(4194303U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(++in);


    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask2), 22-2));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 2), mask22));
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_22_fwd_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask22 =  _mm_set1_epi32(4194303U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    InReg = _mm_load_si128(++in);


    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask2), 22-2));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 2), mask22));
#endif
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask22);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_22_fwd_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask22 =  _mm_set1_epi32(4194303U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 2;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask14), 22-14));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 14), mask22));
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_22_fwd_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask22 =  _mm_set1_epi32(4194303U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 2;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask14), 22-14));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 14), mask22));
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,14) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 22-4));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 4), mask22));
#endif
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_22_fwd_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask22 =  _mm_set1_epi32(4194303U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 4;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask22);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_22_fwd_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask22 =  _mm_set1_epi32(4194303U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 4;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask22);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 22-16));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 16), mask22));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask6), 22-6));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 6), mask22));
#endif
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask22);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask18), 22-18));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 18), mask22));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,18) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 22-8));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 8), mask22));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask22);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask20), 22-20));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 20), mask22));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask10), 22-10));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 10), mask22));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,10) ;
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_22_rev_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask22 =  _mm_set1_epi32(4194303U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 5;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask6), 22-6));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 6), mask22));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_22_rev_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask22 =  _mm_set1_epi32(4194303U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 5;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask6), 22-6));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 6), mask22));
#endif
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

    return;
}

static void
unpack_22_rev_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask22 =  _mm_set1_epi32(4194303U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 6;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask18), 22-18));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 18), mask22));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_22_rev_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask22 =  _mm_set1_epi32(4194303U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 6;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask18), 22-18));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 18), mask22));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,18) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 22-8));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 8), mask22));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_22_rev_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask22 =  _mm_set1_epi32(4194303U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 8;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask22);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_22_rev_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask22 =  _mm_set1_epi32(4194303U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 8;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask22);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask20), 22-20));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 20), mask22));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_22_rev_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask22 =  _mm_set1_epi32(4194303U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 9;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask10), 22-10));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 10), mask22));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_22_rev_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask22 =  _mm_set1_epi32(4194303U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 9;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask10), 22-10));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 22 - 10), mask22));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask14), 23-14));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 23 - 14), mask23));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,14) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask5), 23-5));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 23 - 5), mask23));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,5) , mask23);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask19), 23-19));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 23 - 19), mask23));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,19) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask10), 23-10));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 23 - 10), mask23));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,10) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask1), 23-1));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 23 - 1), mask23));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,1) , mask23);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask15), 23-15));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 23 - 15), mask23));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,15) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask6), 23-6));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 23 - 6), mask23));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask23);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,29) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask20), 23-20));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 23 - 20), mask23));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask11), 23-11));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 23 - 11), mask23));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,11) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask2), 23-2));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 23 - 2), mask23));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask23);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,25) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 23-16));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 23 - 16), mask23));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 24-16));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 16), mask24));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 24-8));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 8), mask24));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 24-16));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 16), mask24));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 24-8));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 8), mask24));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_24_fwd_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask24 =  _mm_set1_epi32(16777215U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask24);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_24_fwd_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask24 =  _mm_set1_epi32(16777215U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask24);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 24-16));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 16), mask24));
#endif
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_24_fwd_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask24 =  _mm_set1_epi32(16777215U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(++in);


    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 24-8));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 8), mask24));
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_24_fwd_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask24 =  _mm_set1_epi32(16777215U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    InReg = _mm_load_si128(++in);


    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 24-8));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 8), mask24));
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_24_fwd_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask24 =  _mm_set1_epi32(16777215U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 3;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask24);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_24_fwd_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask24 =  _mm_set1_epi32(16777215U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 3;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask24);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 24-16));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 16), mask24));
#endif
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_24_fwd_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask24 =  _mm_set1_epi32(16777215U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 4;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 24-8));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 8), mask24));
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_24_fwd_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask24 =  _mm_set1_epi32(16777215U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 4;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 24-8));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 8), mask24));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 24-16));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 16), mask24));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 24-8));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 8), mask24));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 24-16));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 16), mask24));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 24-8));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 8), mask24));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_24_rev_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask24 =  _mm_set1_epi32(16777215U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 6;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask24);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_24_rev_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask24 =  _mm_set1_epi32(16777215U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 6;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask24);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 24-16));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 16), mask24));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_24_rev_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask24 =  _mm_set1_epi32(16777215U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 7;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 24-8));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 8), mask24));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_24_rev_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask24 =  _mm_set1_epi32(16777215U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 7;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 24-8));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 8), mask24));
#endif
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
unpack_24_rev_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask24 =  _mm_set1_epi32(16777215U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 9;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask24);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_24_rev_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask24 =  _mm_set1_epi32(16777215U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 9;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask24);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 24-16));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 16), mask24));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_24_rev_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask24 =  _mm_set1_epi32(16777215U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 10;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 24-8));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 8), mask24));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_24_rev_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask24 =  _mm_set1_epi32(16777215U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 10;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 24-8));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 24 - 8), mask24));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask18), 25-18));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 25 - 18), mask25));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,18) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask11), 25-11));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 25 - 11), mask25));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,11) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 25-4));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 25 - 4), mask25));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask25);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,29) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask22), 25-22));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 25 - 22), mask25));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask15), 25-15));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 25 - 15), mask25));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,15) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 25-8));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 25 - 8), mask25));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask1), 25-1));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 25 - 1), mask25));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,1) , mask25);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask19), 25-19));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 25 - 19), mask25));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,19) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask12), 25-12));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 25 - 12), mask25));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask5), 25-5));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 25 - 5), mask25));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,5) , mask25);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask23), 25-23));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 25 - 23), mask25));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,23) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 25-16));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 25 - 16), mask25));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask20), 26-20));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 20), mask26));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask14), 26-14));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 14), mask26));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,14) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 26-8));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 8), mask26));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask2), 26-2));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 2), mask26));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask26);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask22), 26-22));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 22), mask26));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 26-16));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 16), mask26));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_26_fwd_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask26 =  _mm_set1_epi32(67108863U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask26);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_26_fwd_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask26 =  _mm_set1_epi32(67108863U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask26);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask20), 26-20));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 20), mask26));
#endif
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_26_fwd_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask26 =  _mm_set1_epi32(67108863U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(++in);


    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask14), 26-14));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 14), mask26));
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_26_fwd_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask26 =  _mm_set1_epi32(67108863U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    InReg = _mm_load_si128(++in);


    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask14), 26-14));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 14), mask26));
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,14) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 26-8));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 8), mask26));
#endif
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_26_fwd_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask26 =  _mm_set1_epi32(67108863U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 3;    
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask2), 26-2));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 2), mask26));
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_26_fwd_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask26 =  _mm_set1_epi32(67108863U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 3;    
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask2), 26-2));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 2), mask26));
#endif
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask26);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_26_fwd_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask26 =  _mm_set1_epi32(67108863U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 4;    
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask22), 26-22));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 22), mask26));
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_26_fwd_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask26 =  _mm_set1_epi32(67108863U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 4;    
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask22), 26-22));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 22), mask26));
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 26-16));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 16), mask26));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask10), 26-10));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 10), mask26));
#endif
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,10) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 26-4));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 4), mask26));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask26);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask24), 26-24));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 24), mask26));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask18), 26-18));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 18), mask26));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,18) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask12), 26-12));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 12), mask26));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask6), 26-6));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 6), mask26));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,6) ;
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_26_rev_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask26 =  _mm_set1_epi32(67108863U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 6;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask10), 26-10));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 10), mask26));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_26_rev_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask26 =  _mm_set1_epi32(67108863U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 6;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask10), 26-10));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 10), mask26));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,10) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 26-4));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 4), mask26));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_26_rev_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask26 =  _mm_set1_epi32(67108863U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 8;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask26);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_26_rev_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask26 =  _mm_set1_epi32(67108863U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 8;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask26);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask24), 26-24));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 24), mask26));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_26_rev_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask26 =  _mm_set1_epi32(67108863U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 9;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask18), 26-18));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 18), mask26));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_26_rev_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask26 =  _mm_set1_epi32(67108863U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 9;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask18), 26-18));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 18), mask26));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,18) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask12), 26-12));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 12), mask26));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_26_rev_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask26 =  _mm_set1_epi32(67108863U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 11;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask6), 26-6));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 6), mask26));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_26_rev_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask26 =  _mm_set1_epi32(67108863U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 11;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask6), 26-6));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 26 - 6), mask26));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask22), 27-22));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 27 - 22), mask27));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask17), 27-17));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 27 - 17), mask27));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,17) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask12), 27-12));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 27 - 12), mask27));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask7), 27-7));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 27 - 7), mask27));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,7) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask2), 27-2));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 27 - 2), mask27));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask27);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,29) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask24), 27-24));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 27 - 24), mask27));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask19), 27-19));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 27 - 19), mask27));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,19) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask14), 27-14));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 27 - 14), mask27));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,14) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask9), 27-9));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 27 - 9), mask27));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,9) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 27-4));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 27 - 4), mask27));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask27);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,31) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask26), 27-26));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 27 - 26), mask27));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask21), 27-21));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 27 - 21), mask27));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,21) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 27-16));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 27 - 16), mask27));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask24), 28-24));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 24), mask28));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask20), 28-20));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 20), mask28));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 28-16));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 16), mask28));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask12), 28-12));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 12), mask28));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 28-8));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 8), mask28));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 28-4));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 4), mask28));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,4) ;
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_28_fwd_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask28 =  _mm_set1_epi32(268435455U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask28);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_28_fwd_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask28 =  _mm_set1_epi32(268435455U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask28);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask24), 28-24));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 24), mask28));
#endif
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_28_fwd_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask28 =  _mm_set1_epi32(268435455U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(++in);


    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask20), 28-20));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 20), mask28));
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_28_fwd_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask28 =  _mm_set1_epi32(268435455U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    InReg = _mm_load_si128(++in);


    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask20), 28-20));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 20), mask28));
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 28-16));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 16), mask28));
#endif
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_28_fwd_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask28 =  _mm_set1_epi32(268435455U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 3;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask12), 28-12));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 12), mask28));
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_28_fwd_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask28 =  _mm_set1_epi32(268435455U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 3;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask12), 28-12));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 12), mask28));
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 28-8));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 8), mask28));
#endif
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_28_fwd_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask28 =  _mm_set1_epi32(268435455U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 5;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 28-4));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 4), mask28));
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_28_fwd_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask28 =  _mm_set1_epi32(268435455U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 5;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 28-4));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 4), mask28));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask24), 28-24));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 24), mask28));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask20), 28-20));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 20), mask28));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 28-16));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 16), mask28));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask12), 28-12));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 12), mask28));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 28-8));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 8), mask28));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 28-4));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 4), mask28));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,4) ;
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_28_rev_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask28 =  _mm_set1_epi32(268435455U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 7;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask28);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_28_rev_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask28 =  _mm_set1_epi32(268435455U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 7;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask28);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask24), 28-24));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 24), mask28));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_28_rev_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask28 =  _mm_set1_epi32(268435455U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 8;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask20), 28-20));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 20), mask28));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_28_rev_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask28 =  _mm_set1_epi32(268435455U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 8;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask20), 28-20));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 20), mask28));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 28-16));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 16), mask28));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_28_rev_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask28 =  _mm_set1_epi32(268435455U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 10;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask12), 28-12));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 12), mask28));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_28_rev_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask28 =  _mm_set1_epi32(268435455U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 10;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask12), 28-12));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 12), mask28));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 28-8));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 8), mask28));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_28_rev_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask28 =  _mm_set1_epi32(268435455U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 12;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 28-4));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 4), mask28));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_28_rev_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask28 =  _mm_set1_epi32(268435455U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 12;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 28-4));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 28 - 4), mask28));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask26), 29-26));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 29 - 26), mask29));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask23), 29-23));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 29 - 23), mask29));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,23) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask20), 29-20));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 29 - 20), mask29));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask17), 29-17));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 29 - 17), mask29));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,17) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask14), 29-14));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 29 - 14), mask29));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,14) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask11), 29-11));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 29 - 11), mask29));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,11) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 29-8));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 29 - 8), mask29));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask5), 29-5));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 29 - 5), mask29));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,5) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask2), 29-2));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 29 - 2), mask29));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask29);
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,31) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask28), 29-28));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 29 - 28), mask29));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask25), 29-25));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 29 - 25), mask29));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,25) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask22), 29-22));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 29 - 22), mask29));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask19), 29-19));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 29 - 19), mask29));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,19) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 29-16));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 29 - 16), mask29));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask28), 30-28));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 28), mask30));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask26), 30-26));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 26), mask30));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask24), 30-24));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 24), mask30));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask22), 30-22));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 22), mask30));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask20), 30-20));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 20), mask30));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask18), 30-18));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 18), mask30));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,18) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 30-16));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 16), mask30));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_30_fwd_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i total;
    const __m128i mask30 =  _mm_set1_epi32(1073741823U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask30);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_30_fwd_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask30 =  _mm_set1_epi32(1073741823U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask30);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask28), 30-28));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 28), mask30));
#endif
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_30_fwd_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask30 =  _mm_set1_epi32(1073741823U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    InReg = _mm_load_si128(++in);


    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask26), 30-26));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 26), mask30));
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_30_fwd_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask30 =  _mm_set1_epi32(1073741823U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    InReg = _mm_load_si128(++in);


    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask26), 30-26));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 26), mask30));
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask24), 30-24));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 24), mask30));
#endif
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_30_fwd_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask30 =  _mm_set1_epi32(1073741823U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 3;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask22), 30-22));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 22), mask30));
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_30_fwd_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask30 =  _mm_set1_epi32(1073741823U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 3;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask22), 30-22));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 22), mask30));
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask20), 30-20));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 20), mask30));
#endif
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_30_fwd_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask30 =  _mm_set1_epi32(1073741823U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 5;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask18), 30-18));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 18), mask30));
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_30_fwd_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask30 =  _mm_set1_epi32(1073741823U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 5;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask18), 30-18));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 18), mask30));
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,18) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 30-16));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 16), mask30));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask14), 30-14));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 14), mask30));
#endif
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,14) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask12), 30-12));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 12), mask30));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask10), 30-10));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 10), mask30));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,10) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 30-8));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 8), mask30));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask6), 30-6));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 6), mask30));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,6) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 30-4));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 4), mask30));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,4) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask2), 30-2));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 2), mask30));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,2) ;
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    return;
}

static void
unpack_30_rev_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask30 =  _mm_set1_epi32(1073741823U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 7;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask14), 30-14));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 14), mask30));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_30_rev_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask30 =  _mm_set1_epi32(1073741823U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 7;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask14), 30-14));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 14), mask30));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,14) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask12), 30-12));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 12), mask30));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_30_rev_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask30 =  _mm_set1_epi32(1073741823U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 9;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask10), 30-10));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 10), mask30));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_30_rev_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask30 =  _mm_set1_epi32(1073741823U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 9;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask10), 30-10));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 10), mask30));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,10) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 30-8));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 8), mask30));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_30_rev_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask30 =  _mm_set1_epi32(1073741823U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 11;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask6), 30-6));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 6), mask30));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_30_rev_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask30 =  _mm_set1_epi32(1073741823U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 11;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask6), 30-6));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 6), mask30));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,6) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 30-4));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 4), mask30));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, OutReg);
#else
    total = _mm_add_epi32(total, OutReg);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_30_rev_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask30 =  _mm_set1_epi32(1073741823U);

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 13;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,4) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask2), 30-2));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 2), mask30));
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_30_rev_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;
    const __m128i mask30 =  _mm_set1_epi32(1073741823U);

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 13;
    InReg = _mm_load_si128(in);


    OutReg =   _mm_srli_epi32(InReg,4) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask2), 30-2));
#else
    total = /* OutReg = */ _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 30 - 2), mask30));
#endif
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

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask30), 31-30));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 31 - 30), mask31));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask29), 31-29));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 31 - 29), mask31));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,29) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask28), 31-28));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 31 - 28), mask31));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask27), 31-27));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 31 - 27), mask31));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,27) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask26), 31-26));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 31 - 26), mask31));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask25), 31-25));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 31 - 25), mask31));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,25) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask24), 31-24));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 31 - 24), mask31));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask23), 31-23));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 31 - 23), mask31));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,23) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask22), 31-22));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 31 - 22), mask31));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask21), 31-21));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 31 - 21), mask31));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,21) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask20), 31-20));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 31 - 20), mask31));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask19), 31-19));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 31 - 19), mask31));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,19) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask18), 31-18));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 31 - 18), mask31));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,18) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask17), 31-17));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 31 - 17), mask31));
#endif
    /* total = _mm_add_epi32(total, OutReg); */
    _mm_store_si128(out++, OutReg);

    OutReg =   _mm_srli_epi32(InReg,17) ;
    InReg = _mm_load_si128(++in);

#ifdef MULTIMASK
    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 31-16));
#else
    OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi32(InReg, 31 - 16), mask31));
#endif
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
unpack_32_fwd_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i total;

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    total = _mm_load_si128(in++);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_32_fwd_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i total;

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    total = _mm_load_si128(in++);
    _mm_store_si128(out++, total);

    total = _mm_add_epi32(total, _mm_load_si128(in++));
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_32_fwd_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i total;

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 2;

    total = _mm_load_si128(in++);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_32_fwd_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i total;

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 2;

    total = _mm_load_si128(in++);
    _mm_store_si128(out++, total);

    total = _mm_add_epi32(total, _mm_load_si128(in++));
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_32_fwd_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i total;

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 4;

    total = _mm_load_si128(in++);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_32_fwd_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i total;

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 4;

    total = _mm_load_si128(in++);
    _mm_store_si128(out++, total);

    total = _mm_add_epi32(total, _mm_load_si128(in++));
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_32_fwd_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i total;

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 6;

    total = _mm_load_si128(in++);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_32_fwd_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i total;

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 6;

    total = _mm_load_si128(in++);
    _mm_store_si128(out++, total);

    total = _mm_add_epi32(total, _mm_load_si128(in++));
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
unpack_32_rev_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i total;

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 8;

    total = _mm_load_si128(in++);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_32_rev_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i total;

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 8;

    total = _mm_load_si128(in++);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, _mm_load_si128(in++));
#else
    total = _mm_add_epi32(total, _mm_load_si128(in++));
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_32_rev_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i total;

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 10;

    total = _mm_load_si128(in++);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_32_rev_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i total;

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 10;

    total = _mm_load_si128(in++);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, _mm_load_si128(in++));
#else
    total = _mm_add_epi32(total, _mm_load_si128(in++));
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_32_rev_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i total;

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 12;

    total = _mm_load_si128(in++);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_32_rev_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i total;

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 12;

    total = _mm_load_si128(in++);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, _mm_load_si128(in++));
#else
    total = _mm_add_epi32(total, _mm_load_si128(in++));
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_32_rev_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i total;

#ifdef BRANCH_FREE_ROW_SUM
    __m128i dummy = _mm_set1_epi32(0U);
    _mm_store_si128(out++, dummy);
#endif

    in += 14;

    total = _mm_load_si128(in++);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(dummy,total);
#endif
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_32_rev_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i total;

#ifdef BRANCH_FREE_ROW_SUM
    out++;			/* dummy */
#endif

    in += 14;

    total = _mm_load_si128(in++);
#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(_mm_set1_epi32(0U),total);
#endif
    _mm_store_si128(out++, total);

#ifdef BRANCH_FREE_QTR_BLOCK
    total = _mm_sub_epi32(total, _mm_load_si128(in++));
#else
    total = _mm_add_epi32(total, _mm_load_si128(in++));
#endif
    _mm_store_si128(out++, total);

    return;
}



#endif


#ifndef BRANCH_FREE_ROW_SUM

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

#endif


#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)

static void
vertical_order (UINT4 *vertical, UINT4 *columnar) {

  vertical[0] = columnar[0];		/* remainder 1 */
  vertical[4] = columnar[1];		/* remainder 5 */
  vertical[8] = columnar[2];		/* remainder 9 */
  vertical[12] = columnar[3];		/* remainder 13 */
  vertical[16] = columnar[4];		/* remainder 17 */
  vertical[20] = columnar[5];		/* remainder 21 */
  vertical[24] = columnar[6];		/* remainder 25 */
  vertical[28] = columnar[7];		/* remainder 29 */

  vertical[1] = columnar[8];		/* remainder 2 */
  vertical[5] = columnar[9];		/* remainder 6 */
  vertical[9] = columnar[10];		/* remainder 10 */
  vertical[13] = columnar[11];		/* remainder 14 */
  vertical[17] = columnar[12];		/* remainder 18 */
  vertical[21] = columnar[13];		/* remainder 22 */
  vertical[25] = columnar[14];		/* remainder 26 */
  vertical[29] = columnar[15];		/* remainder 30 */

  vertical[2] = columnar[16];		/* remainder 3 */
  vertical[6] = columnar[17];		/* remainder 7 */
  vertical[10] = columnar[18];		/* remainder 11 */
  vertical[14] = columnar[19];		/* remainder 15 */
  vertical[18] = columnar[20];		/* remainder 19 */
  vertical[22] = columnar[21];		/* remainder 23 */
  vertical[26] = columnar[22];		/* remainder 27 */
  vertical[30] = columnar[23];		/* remainder 31 */

  vertical[3] = columnar[24];		/* remainder 4 */
  vertical[7] = columnar[25];		/* remainder 8 */
  vertical[11] = columnar[26];		/* remainder 12 */
  vertical[15] = columnar[27];		/* remainder 16 */
  vertical[19] = columnar[28];		/* remainder 20 */
  vertical[23] = columnar[29];		/* remainder 24 */
  vertical[27] = columnar[30];		/* remainder 28 */
  vertical[31] = columnar[31];		/* remainder 32 */

  vertical[32] = columnar[32];		/* remainder 63 */
  vertical[36] = columnar[33];		/* remainder 59 */
  vertical[40] = columnar[34];		/* remainder 55 */
  vertical[44] = columnar[35];		/* remainder 51 */
  vertical[48] = columnar[36];		/* remainder 47 */
  vertical[52] = columnar[37];		/* remainder 43 */
  vertical[56] = columnar[38];		/* remainder 39 */
  vertical[60] = columnar[39];		/* remainder 35 */

  vertical[33] = columnar[40];		/* remainder 62 */
  vertical[37] = columnar[41];		/* remainder 58 */
  vertical[41] = columnar[42];		/* remainder 54 */
  vertical[45] = columnar[43];		/* remainder 50 */
  vertical[49] = columnar[44];		/* remainder 46 */
  vertical[53] = columnar[45];		/* remainder 42 */
  vertical[57] = columnar[46];		/* remainder 38 */
  vertical[61] = columnar[47];		/* remainder 34 */

  vertical[34] = columnar[48];		/* remainder 61 */
  vertical[38] = columnar[49];		/* remainder 57 */
  vertical[42] = columnar[50];		/* remainder 53 */
  vertical[46] = columnar[51];		/* remainder 49 */
  vertical[50] = columnar[52];		/* remainder 45 */
  vertical[54] = columnar[53];		/* remainder 41 */
  vertical[58] = columnar[54];		/* remainder 37 */
  vertical[62] = columnar[55];		/* remainder 33 */

  vertical[35] = columnar[56];		/* remainder 60 */
  vertical[39] = columnar[57];		/* remainder 56 */
  vertical[43] = columnar[58];		/* remainder 52 */
  vertical[47] = columnar[59];		/* remainder 48 */
  vertical[51] = columnar[60];		/* remainder 44 */
  vertical[55] = columnar[61];		/* remainder 40 */
  vertical[59] = columnar[62];		/* remainder 36 */
  vertical[63] = columnar[63];		/* remainder 32 */

  return;
}

static void
vertical_order_huge (UINT8 *vertical, UINT4 *columnar) {

  vertical[0] = (UINT8) columnar[0];		/* remainder 1 */
  vertical[4] = (UINT8) columnar[1];		/* remainder 5 */
  vertical[8] = (UINT8) columnar[2];		/* remainder 9 */
  vertical[12] = (UINT8) columnar[3];		/* remainder 13 */
  vertical[16] = (UINT8) columnar[4];		/* remainder 17 */
  vertical[20] = (UINT8) columnar[5];		/* remainder 21 */
  vertical[24] = (UINT8) columnar[6];		/* remainder 25 */
  vertical[28] = (UINT8) columnar[7];		/* remainder 29 */

  vertical[1] = (UINT8) columnar[8];		/* remainder 2 */
  vertical[5] = (UINT8) columnar[9];		/* remainder 6 */
  vertical[9] = (UINT8) columnar[10];		/* remainder 10 */
  vertical[13] = (UINT8) columnar[11];		/* remainder 14 */
  vertical[17] = (UINT8) columnar[12];		/* remainder 18 */
  vertical[21] = (UINT8) columnar[13];		/* remainder 22 */
  vertical[25] = (UINT8) columnar[14];		/* remainder 26 */
  vertical[29] = (UINT8) columnar[15];		/* remainder 30 */

  vertical[2] = (UINT8) columnar[16];		/* remainder 3 */
  vertical[6] = (UINT8) columnar[17];		/* remainder 7 */
  vertical[10] = (UINT8) columnar[18];		/* remainder 11 */
  vertical[14] = (UINT8) columnar[19];		/* remainder 15 */
  vertical[18] = (UINT8) columnar[20];		/* remainder 19 */
  vertical[22] = (UINT8) columnar[21];		/* remainder 23 */
  vertical[26] = (UINT8) columnar[22];		/* remainder 27 */
  vertical[30] = (UINT8) columnar[23];		/* remainder 31 */

  vertical[3] = (UINT8) columnar[24];		/* remainder 4 */
  vertical[7] = (UINT8) columnar[25];		/* remainder 8 */
  vertical[11] = (UINT8) columnar[26];		/* remainder 12 */
  vertical[15] = (UINT8) columnar[27];		/* remainder 16 */
  vertical[19] = (UINT8) columnar[28];		/* remainder 20 */
  vertical[23] = (UINT8) columnar[29];		/* remainder 24 */
  vertical[27] = (UINT8) columnar[30];		/* remainder 28 */
  vertical[31] = (UINT8) columnar[31];		/* remainder 32 */

  vertical[32] = (UINT8) columnar[32];		/* remainder 63 */
  vertical[36] = (UINT8) columnar[33];		/* remainder 59 */
  vertical[40] = (UINT8) columnar[34];		/* remainder 55 */
  vertical[44] = (UINT8) columnar[35];		/* remainder 51 */
  vertical[48] = (UINT8) columnar[36];		/* remainder 47 */
  vertical[52] = (UINT8) columnar[37];		/* remainder 43 */
  vertical[56] = (UINT8) columnar[38];		/* remainder 39 */
  vertical[60] = (UINT8) columnar[39];		/* remainder 35 */

  vertical[33] = (UINT8) columnar[40];		/* remainder 62 */
  vertical[37] = (UINT8) columnar[41];		/* remainder 58 */
  vertical[41] = (UINT8) columnar[42];		/* remainder 54 */
  vertical[45] = (UINT8) columnar[43];		/* remainder 50 */
  vertical[49] = (UINT8) columnar[44];		/* remainder 46 */
  vertical[53] = (UINT8) columnar[45];		/* remainder 42 */
  vertical[57] = (UINT8) columnar[46];		/* remainder 38 */
  vertical[61] = (UINT8) columnar[47];		/* remainder 34 */

  vertical[34] = (UINT8) columnar[48];		/* remainder 61 */
  vertical[38] = (UINT8) columnar[49];		/* remainder 57 */
  vertical[42] = (UINT8) columnar[50];		/* remainder 53 */
  vertical[46] = (UINT8) columnar[51];		/* remainder 49 */
  vertical[50] = (UINT8) columnar[52];		/* remainder 45 */
  vertical[54] = (UINT8) columnar[53];		/* remainder 41 */
  vertical[58] = (UINT8) columnar[54];		/* remainder 37 */
  vertical[62] = (UINT8) columnar[55];		/* remainder 33 */

  vertical[35] = (UINT8) columnar[56];		/* remainder 60 */
  vertical[39] = (UINT8) columnar[57];		/* remainder 56 */
  vertical[43] = (UINT8) columnar[58];		/* remainder 52 */
  vertical[47] = (UINT8) columnar[59];		/* remainder 48 */
  vertical[51] = (UINT8) columnar[60];		/* remainder 44 */
  vertical[55] = (UINT8) columnar[61];		/* remainder 40 */
  vertical[59] = (UINT8) columnar[62];		/* remainder 36 */
  vertical[63] = (UINT8) columnar[63];		/* remainder 32 */

  return;
}


#else
static void
vertical_order_fwd (UINT4 *vertical, UINT4 *columnar) {

  vertical[0] = columnar[0];		/* remainder 1 */
  vertical[4] = columnar[1];		/* remainder 5 */
  vertical[8] = columnar[2];		/* remainder 9 */
  vertical[12] = columnar[3];		/* remainder 13 */
  vertical[16] = columnar[4];		/* remainder 17 */
  vertical[20] = columnar[5];		/* remainder 21 */
  vertical[24] = columnar[6];		/* remainder 25 */
  vertical[28] = columnar[7];		/* remainder 29 */

  vertical[1] = columnar[8];		/* remainder 2 */
  vertical[5] = columnar[9];		/* remainder 6 */
  vertical[9] = columnar[10];		/* remainder 10 */
  vertical[13] = columnar[11];		/* remainder 14 */
  vertical[17] = columnar[12];		/* remainder 18 */
  vertical[21] = columnar[13];		/* remainder 22 */
  vertical[25] = columnar[14];		/* remainder 26 */
  vertical[29] = columnar[15];		/* remainder 30 */

  vertical[2] = columnar[16];		/* remainder 3 */
  vertical[6] = columnar[17];		/* remainder 7 */
  vertical[10] = columnar[18];		/* remainder 11 */
  vertical[14] = columnar[19];		/* remainder 15 */
  vertical[18] = columnar[20];		/* remainder 19 */
  vertical[22] = columnar[21];		/* remainder 23 */
  vertical[26] = columnar[22];		/* remainder 27 */
  vertical[30] = columnar[23];		/* remainder 31 */

  vertical[3] = columnar[24];		/* remainder 4 */
  vertical[7] = columnar[25];		/* remainder 8 */
  vertical[11] = columnar[26];		/* remainder 12 */
  vertical[15] = columnar[27];		/* remainder 16 */
  vertical[19] = columnar[28];		/* remainder 20 */
  vertical[23] = columnar[29];		/* remainder 24 */
  vertical[27] = columnar[30];		/* remainder 28 */
  vertical[31] = columnar[31];		/* remainder 32 */

  return;
}

static void
vertical_order_rev (UINT4 *vertical, UINT4 *columnar) {

  vertical[0] = columnar[0];		/* remainder 63 */
  vertical[4] = columnar[1];		/* remainder 59 */
  vertical[8] = columnar[2];		/* remainder 55 */
  vertical[12] = columnar[3];		/* remainder 51 */
  vertical[16] = columnar[4];		/* remainder 47 */
  vertical[20] = columnar[5];		/* remainder 43 */
  vertical[24] = columnar[6];		/* remainder 39 */
  vertical[28] = columnar[7];		/* remainder 35 */

  vertical[1] = columnar[8];		/* remainder 62 */
  vertical[5] = columnar[9];		/* remainder 58 */
  vertical[9] = columnar[10];		/* remainder 54 */
  vertical[13] = columnar[11];		/* remainder 50 */
  vertical[17] = columnar[12];		/* remainder 46 */
  vertical[21] = columnar[13];		/* remainder 42 */
  vertical[25] = columnar[14];		/* remainder 38 */
  vertical[29] = columnar[15];		/* remainder 34 */

  vertical[2] = columnar[16];		/* remainder 61 */
  vertical[6] = columnar[17];		/* remainder 57 */
  vertical[10] = columnar[18];		/* remainder 53 */
  vertical[14] = columnar[19];		/* remainder 49 */
  vertical[18] = columnar[20];		/* remainder 45 */
  vertical[22] = columnar[21];		/* remainder 41 */
  vertical[26] = columnar[22];		/* remainder 37 */
  vertical[30] = columnar[23];		/* remainder 33 */

  vertical[3] = columnar[24];		/* remainder 60 */
  vertical[7] = columnar[25];		/* remainder 56 */
  vertical[11] = columnar[26];		/* remainder 52 */
  vertical[15] = columnar[27];		/* remainder 48 */
  vertical[19] = columnar[28];		/* remainder 44 */
  vertical[23] = columnar[29];		/* remainder 40 */
  vertical[27] = columnar[30];		/* remainder 36 */
  vertical[31] = columnar[31];		/* remainder 32 */

  return;
}

#if defined(HAVE_64_BIT) && (defined(UTILITYP) || defined(LARGE_GENOMES))
static void
vertical_order_huge_fwd (UINT8 *vertical, UINT4 *columnar) {

  vertical[0] = (UINT8) columnar[0];		/* remainder 1 */
  vertical[4] = (UINT8) columnar[1];		/* remainder 5 */
  vertical[8] = (UINT8) columnar[2];		/* remainder 9 */
  vertical[12] = (UINT8) columnar[3];		/* remainder 13 */
  vertical[16] = (UINT8) columnar[4];		/* remainder 17 */
  vertical[20] = (UINT8) columnar[5];		/* remainder 21 */
  vertical[24] = (UINT8) columnar[6];		/* remainder 25 */
  vertical[28] = (UINT8) columnar[7];		/* remainder 29 */

  vertical[1] = (UINT8) columnar[8];		/* remainder 2 */
  vertical[5] = (UINT8) columnar[9];		/* remainder 6 */
  vertical[9] = (UINT8) columnar[10];		/* remainder 10 */
  vertical[13] = (UINT8) columnar[11];		/* remainder 14 */
  vertical[17] = (UINT8) columnar[12];		/* remainder 18 */
  vertical[21] = (UINT8) columnar[13];		/* remainder 22 */
  vertical[25] = (UINT8) columnar[14];		/* remainder 26 */
  vertical[29] = (UINT8) columnar[15];		/* remainder 30 */

  vertical[2] = (UINT8) columnar[16];		/* remainder 3 */
  vertical[6] = (UINT8) columnar[17];		/* remainder 7 */
  vertical[10] = (UINT8) columnar[18];		/* remainder 11 */
  vertical[14] = (UINT8) columnar[19];		/* remainder 15 */
  vertical[18] = (UINT8) columnar[20];		/* remainder 19 */
  vertical[22] = (UINT8) columnar[21];		/* remainder 23 */
  vertical[26] = (UINT8) columnar[22];		/* remainder 27 */
  vertical[30] = (UINT8) columnar[23];		/* remainder 31 */

  vertical[3] = (UINT8) columnar[24];		/* remainder 4 */
  vertical[7] = (UINT8) columnar[25];		/* remainder 8 */
  vertical[11] = (UINT8) columnar[26];		/* remainder 12 */
  vertical[15] = (UINT8) columnar[27];		/* remainder 16 */
  vertical[19] = (UINT8) columnar[28];		/* remainder 20 */
  vertical[23] = (UINT8) columnar[29];		/* remainder 24 */
  vertical[27] = (UINT8) columnar[30];		/* remainder 28 */
  vertical[31] = (UINT8) columnar[31];		/* remainder 32 */

  return;
}
#endif

#if defined(HAVE_64_BIT) && (defined(UTILITYP) || defined(LARGE_GENOMES))
static void
vertical_order_huge_rev (UINT8 *vertical, UINT4 *columnar) {

  vertical[0] = (UINT8) columnar[0];		/* remainder 63 */
  vertical[4] = (UINT8) columnar[1];		/* remainder 59 */
  vertical[8] = (UINT8) columnar[2];		/* remainder 55 */
  vertical[12] = (UINT8) columnar[3];		/* remainder 51 */
  vertical[16] = (UINT8) columnar[4];		/* remainder 47 */
  vertical[20] = (UINT8) columnar[5];		/* remainder 43 */
  vertical[24] = (UINT8) columnar[6];		/* remainder 39 */
  vertical[28] = (UINT8) columnar[7];		/* remainder 35 */

  vertical[1] = (UINT8) columnar[8];		/* remainder 62 */
  vertical[5] = (UINT8) columnar[9];		/* remainder 58 */
  vertical[9] = (UINT8) columnar[10];		/* remainder 54 */
  vertical[13] = (UINT8) columnar[11];		/* remainder 50 */
  vertical[17] = (UINT8) columnar[12];		/* remainder 46 */
  vertical[21] = (UINT8) columnar[13];		/* remainder 42 */
  vertical[25] = (UINT8) columnar[14];		/* remainder 38 */
  vertical[29] = (UINT8) columnar[15];		/* remainder 34 */

  vertical[2] = (UINT8) columnar[16];		/* remainder 61 */
  vertical[6] = (UINT8) columnar[17];		/* remainder 57 */
  vertical[10] = (UINT8) columnar[18];		/* remainder 53 */
  vertical[14] = (UINT8) columnar[19];		/* remainder 49 */
  vertical[18] = (UINT8) columnar[20];		/* remainder 45 */
  vertical[22] = (UINT8) columnar[21];		/* remainder 41 */
  vertical[26] = (UINT8) columnar[22];		/* remainder 37 */
  vertical[30] = (UINT8) columnar[23];		/* remainder 33 */

  vertical[3] = (UINT8) columnar[24];		/* remainder 60 */
  vertical[7] = (UINT8) columnar[25];		/* remainder 56 */
  vertical[11] = (UINT8) columnar[26];		/* remainder 52 */
  vertical[15] = (UINT8) columnar[27];		/* remainder 48 */
  vertical[19] = (UINT8) columnar[28];		/* remainder 44 */
  vertical[23] = (UINT8) columnar[29];		/* remainder 40 */
  vertical[27] = (UINT8) columnar[30];		/* remainder 36 */
  vertical[31] = (UINT8) columnar[31];		/* remainder 32 */

  return;
}
#endif

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
  {{unpack_00_1, unpack_00_2, unpack_00_2, unpack_00_1,
    unpack_00_1, unpack_00_2, unpack_00_2, unpack_00_1,
    unpack_00_1, unpack_00_2, unpack_00_2, unpack_00_1,
    unpack_00_1, unpack_00_2, unpack_00_2, unpack_00_1,
    unpack_00_0},

   {unpack_02_fwd_1, unpack_02_fwd_2, unpack_02_rev_2, unpack_02_rev_1,
    unpack_02_fwd_3, unpack_02_fwd_4, unpack_02_rev_4, unpack_02_rev_3,
    unpack_02_fwd_5, unpack_02_fwd_6, unpack_02_rev_6, unpack_02_rev_5, 
    unpack_02_fwd_7, unpack_02_fwd_8, unpack_02_rev_8, unpack_02_rev_7,
    unpack_00_0},

   {unpack_04_fwd_1, unpack_04_fwd_2, unpack_04_rev_2, unpack_04_rev_1,
    unpack_04_fwd_3, unpack_04_fwd_4, unpack_04_rev_4, unpack_04_rev_3,
    unpack_04_fwd_5, unpack_04_fwd_6, unpack_04_rev_6, unpack_04_rev_5, 
    unpack_04_fwd_7, unpack_04_fwd_8, unpack_04_rev_8, unpack_04_rev_7,
    unpack_00_0},

   {unpack_06_fwd_1, unpack_06_fwd_2, unpack_06_rev_2, unpack_06_rev_1,
    unpack_06_fwd_3, unpack_06_fwd_4, unpack_06_rev_4, unpack_06_rev_3,
    unpack_06_fwd_5, unpack_06_fwd_6, unpack_06_rev_6, unpack_06_rev_5, 
    unpack_06_fwd_7, unpack_06_fwd_8, unpack_06_rev_8, unpack_06_rev_7,
    unpack_00_0},

   {unpack_08_fwd_1, unpack_08_fwd_2, unpack_08_rev_2, unpack_08_rev_1,
    unpack_08_fwd_3, unpack_08_fwd_4, unpack_08_rev_4, unpack_08_rev_3,
    unpack_08_fwd_5, unpack_08_fwd_6, unpack_08_rev_6, unpack_08_rev_5, 
    unpack_08_fwd_7, unpack_08_fwd_8, unpack_08_rev_8, unpack_08_rev_7,
    unpack_00_0},

   {unpack_10_fwd_1, unpack_10_fwd_2, unpack_10_rev_2, unpack_10_rev_1,
    unpack_10_fwd_3, unpack_10_fwd_4, unpack_10_rev_4, unpack_10_rev_3,
    unpack_10_fwd_5, unpack_10_fwd_6, unpack_10_rev_6, unpack_10_rev_5, 
    unpack_10_fwd_7, unpack_10_fwd_8, unpack_10_rev_8, unpack_10_rev_7,
    unpack_00_0},

   {unpack_12_fwd_1, unpack_12_fwd_2, unpack_12_rev_2, unpack_12_rev_1,
    unpack_12_fwd_3, unpack_12_fwd_4, unpack_12_rev_4, unpack_12_rev_3,
    unpack_12_fwd_5, unpack_12_fwd_6, unpack_12_rev_6, unpack_12_rev_5, 
    unpack_12_fwd_7, unpack_12_fwd_8, unpack_12_rev_8, unpack_12_rev_7,
    unpack_00_0},

   {unpack_14_fwd_1, unpack_14_fwd_2, unpack_14_rev_2, unpack_14_rev_1,
    unpack_14_fwd_3, unpack_14_fwd_4, unpack_14_rev_4, unpack_14_rev_3,
    unpack_14_fwd_5, unpack_14_fwd_6, unpack_14_rev_6, unpack_14_rev_5, 
    unpack_14_fwd_7, unpack_14_fwd_8, unpack_14_rev_8, unpack_14_rev_7,
    unpack_00_0},

   {unpack_16_fwd_1, unpack_16_fwd_2, unpack_16_rev_2, unpack_16_rev_1,
    unpack_16_fwd_3, unpack_16_fwd_4, unpack_16_rev_4, unpack_16_rev_3,
    unpack_16_fwd_5, unpack_16_fwd_6, unpack_16_rev_6, unpack_16_rev_5, 
    unpack_16_fwd_7, unpack_16_fwd_8, unpack_16_rev_8, unpack_16_rev_7,
    unpack_00_0},

   {unpack_18_fwd_1, unpack_18_fwd_2, unpack_18_rev_2, unpack_18_rev_1,
    unpack_18_fwd_3, unpack_18_fwd_4, unpack_18_rev_4, unpack_18_rev_3,
    unpack_18_fwd_5, unpack_18_fwd_6, unpack_18_rev_6, unpack_18_rev_5, 
    unpack_18_fwd_7, unpack_18_fwd_8, unpack_18_rev_8, unpack_18_rev_7,
    unpack_00_0},

   {unpack_20_fwd_1, unpack_20_fwd_2, unpack_20_rev_2, unpack_20_rev_1,
    unpack_20_fwd_3, unpack_20_fwd_4, unpack_20_rev_4, unpack_20_rev_3,
    unpack_20_fwd_5, unpack_20_fwd_6, unpack_20_rev_6, unpack_20_rev_5, 
    unpack_20_fwd_7, unpack_20_fwd_8, unpack_20_rev_8, unpack_20_rev_7,
    unpack_00_0},

   {unpack_22_fwd_1, unpack_22_fwd_2, unpack_22_rev_2, unpack_22_rev_1,
    unpack_22_fwd_3, unpack_22_fwd_4, unpack_22_rev_4, unpack_22_rev_3,
    unpack_22_fwd_5, unpack_22_fwd_6, unpack_22_rev_6, unpack_22_rev_5, 
    unpack_22_fwd_7, unpack_22_fwd_8, unpack_22_rev_8, unpack_22_rev_7,
    unpack_00_0},

   {unpack_24_fwd_1, unpack_24_fwd_2, unpack_24_rev_2, unpack_24_rev_1,
    unpack_24_fwd_3, unpack_24_fwd_4, unpack_24_rev_4, unpack_24_rev_3,
    unpack_24_fwd_5, unpack_24_fwd_6, unpack_24_rev_6, unpack_24_rev_5, 
    unpack_24_fwd_7, unpack_24_fwd_8, unpack_24_rev_8, unpack_24_rev_7,
    unpack_00_0},

   {unpack_26_fwd_1, unpack_26_fwd_2, unpack_26_rev_2, unpack_26_rev_1,
    unpack_26_fwd_3, unpack_26_fwd_4, unpack_26_rev_4, unpack_26_rev_3,
    unpack_26_fwd_5, unpack_26_fwd_6, unpack_26_rev_6, unpack_26_rev_5, 
    unpack_26_fwd_7, unpack_26_fwd_8, unpack_26_rev_8, unpack_26_rev_7,
    unpack_00_0},

   {unpack_28_fwd_1, unpack_28_fwd_2, unpack_28_rev_2, unpack_28_rev_1,
    unpack_28_fwd_3, unpack_28_fwd_4, unpack_28_rev_4, unpack_28_rev_3,
    unpack_28_fwd_5, unpack_28_fwd_6, unpack_28_rev_6, unpack_28_rev_5, 
    unpack_28_fwd_7, unpack_28_fwd_8, unpack_28_rev_8, unpack_28_rev_7,
    unpack_00_0},

   {unpack_30_fwd_1, unpack_30_fwd_2, unpack_30_rev_2, unpack_30_rev_1,
    unpack_30_fwd_3, unpack_30_fwd_4, unpack_30_rev_4, unpack_30_rev_3,
    unpack_30_fwd_5, unpack_30_fwd_6, unpack_30_rev_6, unpack_30_rev_5, 
    unpack_30_fwd_7, unpack_30_fwd_8, unpack_30_rev_8, unpack_30_rev_7,
    unpack_00_0},

   {unpack_32_fwd_1, unpack_32_fwd_2, unpack_32_rev_2, unpack_32_rev_1,
    unpack_32_fwd_3, unpack_32_fwd_4, unpack_32_rev_4, unpack_32_rev_3,
    unpack_32_fwd_5, unpack_32_fwd_6, unpack_32_rev_6, unpack_32_rev_5, 
    unpack_32_fwd_7, unpack_32_fwd_8, unpack_32_rev_8, unpack_32_rev_7,
    unpack_00_0},

};
   
#endif


#define DIFFERENTIAL_METAINFO_SIZE 2
#define PAIRED_METAINFO_SIZE 3

#define get_column(s) (s) & 3 /* Not s % 4, which fails on negative values */
#define get_row(s) (s) >> 2 /* Not s / 4, which fails on negative values */


#if 0
/* Use Bitpack64_read_two instead */

/* bitpackpages: A list of b-mers (12-mers by default), ending with -1U */
UINT4
Bitpack64_offsetptr (UINT4 *end0, Storedoligomer_T oligo, UINT4 *bitpackptrs, UINT4 *bitpackcomp) {
  Storedoligomer_T bmer;
  UINT4 *info, nwritten, packsize_div2;
  int delta, remainder0, remainder1, quarter_block_0, quarter_block_1, column, row;
#ifdef HAVE_SSE2
#ifdef BRANCH_FREE_ROW_SUM
  __m128i diffs[3];
#else
  __m128i diffs[2];
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
  UINT4 psums[5];		/* Need 5 to handle case where remainder == 64 */
#else
#endif
  __m128i *bitpack;
  UINT4 *_diffs;
#else
  UINT4 ptr;
  int k, i;
  UINT4 diffs[BLOCKSIZE+1], *bitpack;
#endif
#ifdef DEBUG
  UINT4 offsets[BLOCKSIZE+1];
#endif


  bmer = oligo/BLOCKSIZE;
  info = &(bitpackptrs[bmer * DIFFERENTIAL_METAINFO_SIZE]);

  debug(printf("Entered Bitpack64_offsetptr with oligo %u => bmer %u\n",oligo,bmer));

  nwritten = info[0];		/* In 128-bit registers */
#ifdef HAVE_SSE2  
  bitpack = (__m128i *) &(bitpackcomp[nwritten*4]);
#else
  bitpack = (UINT4 *) &(bitpackcomp[nwritten*4]);
#endif

  /* packsize = (info[DIFFERENTIAL_METAINFO_SIZE] - nwritten)*2; */
  packsize_div2 = (info[DIFFERENTIAL_METAINFO_SIZE] - nwritten);

  remainder0 = oligo % BLOCKSIZE;
  remainder1 = remainder0 + 1;
  if (remainder1 == 32) {
    /* For the case (31,32), treat 32 as being in the first half-block.  Otherwise, 32 is in the second half-block */
    quarter_block_0 = quarter_block_1 = 1;
  } else {
    quarter_block_0 = remainder0 / 16;
    quarter_block_1 = remainder1 / 16;
  }

  debug(Bitpack64_block_offsets(offsets,oligo,bitpackptrs,bitpackcomp));

#ifdef HAVE_SSE2
  _diffs = (UINT4 *) diffs;	/* Assumes a dummy register in diffs[0] */

#ifdef BRANCH_FREE_QTR_BLOCK
  psums[0] = psums[1] = info[1];
  psums[2] = psums[3] = psums[4] = info[DIFFERENTIAL_METAINFO_SIZE+1];

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

  /* Don't have to handle remainder == 0 as a special case for obtaining
     ascending[n], if we use " <= n" in Bitpack64_write_differential */

  if (quarter_block_1 <= 1) {
    delta = remainder1 - 1;
    column = get_column(delta);
    row = get_row(delta);
    debug(printf("quarter-block %d, delta %d, column %d, row %d\n",quarter_block_1,delta,column,row));

    (unpacker_table[packsize_div2][column*4 + quarter_block_1])(diffs,bitpack);
#ifdef BRANCH_FREE_ROW_SUM
    *end0 = info[1] + _diffs[row+1] + _diffs[row+2] + _diffs[row+3] + _diffs[row+4];
#else
    assign_sum_fwd(*end0,info[1],_diffs,row);
#endif

  } else {
    delta = 63 - remainder1;
    column = get_column(delta);
    row = get_row(delta);
    debug(printf("quarter-block %d, delta %d, column %d, row %d\n",quarter_block_1,delta,column,row));

    (unpacker_table[packsize_div2][column*4 + quarter_block_1])(diffs,bitpack);
#ifdef BRANCH_FREE_ROW_SUM
    *end0 = info[DIFFERENTIAL_METAINFO_SIZE+1] - _diffs[row+1] - _diffs[row+2] - _diffs[row+3] - _diffs[row+4];
#else
    assign_sum_rev(*end0,info[DIFFERENTIAL_METAINFO_SIZE+1],_diffs,row);
#endif
  }

  if (quarter_block_0 <= 1) {
    delta = remainder0 - 1;
    column = get_column(delta);
    row = get_row(delta);
    debug(printf("quarter-block %d, delta %d, column %d, row %d\n",quarter_block_0,delta,column,row));

    (unpacker_table[packsize_div2][column*4 + quarter_block_0])(diffs,bitpack);
#ifdef BRANCH_FREE_ROW_SUM
    return info[1] + _diffs[row+1] + _diffs[row+2] + _diffs[row+3] + _diffs[row+4];
#else
    return_sum_fwd(info[1],_diffs,row);
#endif

  } else {
    delta = 63 - remainder0;
    column = get_column(delta);
    row = get_row(delta);
    debug(printf("quarter-block %d, delta %d, column %d, row %d\n",quarter_block_0,delta,column,row));

    (unpacker_table[packsize_div2][column*4 + quarter_block_0])(diffs,bitpack);
#ifdef BRANCH_FREE_ROW_SUM
    return info[DIFFERENTIAL_METAINFO_SIZE+1] - _diffs[row+1] - _diffs[row+2] - _diffs[row+3] - _diffs[row+4];
#else
    return_sum_rev(info[DIFFERENTIAL_METAINFO_SIZE+1],_diffs,row);
#endif
  }

#endif

#else  /* HAVE_SSE2 */

  offset0 = info[1];
  offset1 = info[DIFFERENTIAL_METAINFO_SIZE+1];

  /* Unpack all 64 diffs for non-SIMD */
  (unpacker_all_table[packsize_div2*2])(&(diffs[1]),bitpack);

#ifdef DEBUG
  printf("oligo: %08X, remainder %d, offset0 %u, offset1 %u\n",
	 oligo,oligo % BLOCKSIZE,info[1],info[DIFFERENTIAL_METAINFO_SIZE+1]);
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

#endif	/* HAVE_SSE2 */

}
#endif


#if 0
/* Use Bitpack64_read_two_huge instead */

/* bitpackpages: A list of b-mers (12-mers by default), ending with -1U */
UINT8
Bitpack64_offsetptr_huge (UINT8 *end0, Storedoligomer_T oligo,
			  UINT4 *bitpackpages, UINT4 *bitpackptrs, UINT4 *bitpackcomp) {
  Storedoligomer_T bmer;
  UINT4 *info, nwritten;
  UINT8 offset0, offset1;
  UINT4 packsize_div2;
  int delta, remainder0, remainder1, quarter_block_0, quarter_block_1, column, row;
#ifdef HAVE_SSE2
#ifdef BRANCH_FREE_ROW_SUM
  __m128i diffs[3];
#else
  __m128i diffs[2];
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
  UINT8 psums[5];		/* Need 5 to handle case where remainder == 64 */
#endif
  __m128i *bitpack;
  UINT4 *_diffs;
#else
  UINT4 ptr;
  int column, row, k, i;
  UINT4 diffs[BLOCKSIZE+1], *bitpack;
#endif
  UINT4 *pageptr;
#ifdef DEBUG
  UINT4 offsets[BLOCKSIZE+1];
#endif


  bmer = oligo/BLOCKSIZE;
  info = &(bitpackptrs[bmer * DIFFERENTIAL_METAINFO_SIZE]);

  debug(printf("Entered Bitpack64_offsetptr_huge with oligo %u => bmer %u\n",oligo,bmer));

  nwritten = info[0];
#ifdef HAVE_SSE2  
  bitpack = (__m128i *) &(bitpackcomp[nwritten*4]);
#else
  bitpack = (UINT4 *) &(bitpackcomp[nwritten*4]);
#endif

  if (bitpackpages == NULL) {
    offset0 = info[1];
    offset1 = info[DIFFERENTIAL_METAINFO_SIZE+1];
  } else {
    offset0 = 0UL;
    pageptr = bitpackpages;
    while (bmer >= *pageptr) {
      offset0 += POSITIONS_PAGE;
      pageptr++;
    }

    offset1 = offset0;
    if (bmer + 1 >= *pageptr) {
      offset1 += POSITIONS_PAGE;
      /* pageptr++; */
    }

    offset0 += info[1];
    offset1 += info[DIFFERENTIAL_METAINFO_SIZE+1];
  }

  /* packsize = (info[DIFFERENTIAL_METAINFO_SIZE] - nwritten)*2; */
  packsize_div2 = (info[DIFFERENTIAL_METAINFO_SIZE] - nwritten);

  remainder0 = oligo % BLOCKSIZE;
  remainder1 = remainder0 + 1;
  if (remainder1 == 32) {
    /* For the case (31,32), treat 32 as being in the first half-block.  Otherwise, 32 is in the second half-block */
    quarter_block_0 = quarter_block_1 = 1;
  } else {
    quarter_block_0 = remainder0 / 16;
    quarter_block_1 = remainder1 / 16;
  }

  debug(Bitpack64_block_offsets(offsets,oligo,bitpackptrs,bitpackcomp));

#ifdef HAVE_SSE2
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

  if (quarter_block_1 <= 1) {
    delta = remainder1 - 1;
    column = get_column(delta);
    row = get_row(delta);
    debug(printf("quarter-block %d, delta %d, column %d, row %d\n",quarter_block_1,delta,column,row));

    (unpacker_table[packsize_div2][column*4 + quarter_block_1])(diffs,bitpack);
#ifdef BRANCH_FREE_ROW_SUM
    *end0 = offset0 + _diffs[row+1] + _diffs[row+2] + _diffs[row+3] + _diffs[row+4];
#else
    assign_sum_fwd(*end0,offset0,_diffs,row);
#endif

  } else {
    delta = 63 - remainder1;
    column = get_column(delta);
    row = get_row(delta);
    debug(printf("quarter-block %d, delta %d, column %d, row %d\n",quarter_block_1,delta,column,row));

    (unpacker_table[packsize_div2][column*4 + quarter_block_1])(diffs,bitpack);
#ifdef BRANCH_FREE_ROW_SUM
    *end0 = offset1 - _diffs[row+1] - _diffs[row+2] - _diffs[row+3] - _diffs[row+4];
#else
    assign_sum_rev(*end0,offset1,_diffs,row);
#endif
  }

  if (quarter_block_0 <= 1) {
    delta = remainder0 - 1;
    column = get_column(delta);
    row = get_row(delta);
    debug(printf("quarter-block %d, delta %d, column %d, row %d\n",quarter_block_0,delta,column,row));

    (unpacker_table[packsize_div2][column*4 + quarter_block_0])(diffs,bitpack);
#ifdef BRANCH_FREE_ROW_SUM
    return offset0 + _diffs[row+1] + _diffs[row+2] + _diffs[row+3] + _diffs[row+4];
#else
    return_sum_fwd(offset0,_diffs,row);
#endif

  } else {
    delta = 63 - remainder0;
    column = get_column(delta);
    row = get_row(delta);
    debug(printf("quarter-block %d, delta %d, column %d, row %d\n",quarter_block_0,delta,column,row));

    (unpacker_table[packsize_div2][column*4 + quarter_block_0])(diffs,bitpack);
#ifdef BRANCH_FREE_ROW_SUM
    return offset1 - _diffs[row+1] - _diffs[row+2] - _diffs[row+3] - _diffs[row+4];
#else
    return_sum_rev(offset1,_diffs,row);
#endif
  }

#endif


#else  /* HAVE_SSE2 */

  /* Unpack all 64 diffs for non-SIMD */
  (unpacker_all_table[packsize_div2*2])(&(diffs[1]),bitpack);

#ifdef DEBUG
  printf("oligo: %08X, remainder %d, offset0 %u, offset1 %u\n",
	 oligo,oligo % BLOCKSIZE,info[1],info[DIFFERENTIAL_METAINFO_SIZE+1]);
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
    /* Compute necessary cumulative sums */
    ptr = offset0 = info[1];

    /* Add 1 for start at diffs[1], and 1 to leave the first element intact */
    diffs[0] = 0;
    column = (remainder - 1) % 4; /* Goes from 0 to 3 */
    row = (remainder - 1) / 4;
    debug(printf("column %d, row %d\n",column,row));
    
    for (k = column*2 + 1, i = 0; i <= row; k += BLOCKSIZE/4, i++) {
      debug(printf("Adding diffs[%d] = %u\n",k,diffs[k]));
      ptr += diffs[k];
    }

  } else if (remainder <= 32) {
    /* Compute necessary cumulative sums */
    ptr = offset0;

    /* Add 1 for start at diffs[1], and 1 to leave the first element intact */
    diffs[0] = 0;
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

#endif	/* HAVE_SSE2 */
}
#endif


#if 0
/* bitpackpages: A list of b-mers (12-mers by default), ending with -1U */
UINT4
Bitpack64_offsetptr_paired (UINT4 *end0, Storedoligomer_T oligo, UINT4 *bitpackptrs, UINT4 *bitpackcomp) {
  UINT4 ptr0;
  Storedoligomer_T bmer;
  UINT4 *info, nwritten, packsize_div2;
  int delta, remainder, quarter_block, column, row;
#ifdef HAVE_SSE2
#ifdef BRANCH_FREE_ROW_SUM
  __m128i diffs[3];
#else
  __m128i diffs[2];
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
  UINT4 psums[5];		/* Need 5 to handle case where remainder == 64 */
#else
#endif
  __m128i *bitpack;
  UINT4 *_diffs;
#else
  UINT4 ptr;
  int k, i;
  UINT4 diffs[BLOCKSIZE+1], *bitpack;
#endif
#ifdef DEBUG
  UINT4 offsets[BLOCKSIZE+1];
#endif


  bmer = oligo/BLOCKSIZE;
  info = &(bitpackptrs[bmer * PAIRED_METAINFO_SIZE]);

  debug(printf("Entered Bitpack64_offsetptr_paired with oligo %u => bmer %u\n",oligo,bmer));

  nwritten = info[0];		/* In 128-bit registers */
#ifdef HAVE_SSE2  
  bitpack = (__m128i *) &(bitpackcomp[nwritten*4]);
#else
  bitpack = (UINT4 *) &(bitpackcomp[nwritten*4]);
#endif

  /* packsize = (info[DIFFERENTIAL_METAINFO_SIZE] - nwritten)*2; */
  packsize_div2 = (info[2] - nwritten);

  remainder = oligo % BLOCKSIZE;
  if (remainder == 31) {
    /* For the case (31,32), treat 32 as being in the first half-block.  Otherwise, 32 is in the second half-block */
    quarter_block = 1;
  } else {
    quarter_block = remainder / 16;
  }

  debug(Bitpack64_block_offsets(offsets,oligo,bitpackptrs,bitpackcomp));

#ifdef HAVE_SSE2
  _diffs = (UINT4 *) diffs;	/* Assumes a dummy register in diffs[0] */

  if (quarter_block <= 1) {
    delta = remainder - 1;
    column = get_column(delta);
    row = get_row(delta);
    debug(printf("quarter-block %d, delta %d, column %d, row %d\n",quarter_block,delta,column,row));

    (unpacker_table[packsize_div2][column*4 + quarter_block])(diffs,bitpack);
#ifdef BRANCH_FREE_ROW_SUM
    ptr0 = info[1] + _diffs[row+1] + _diffs[row+2] + _diffs[row+3] + _diffs[row+4];
#else
    assign_sum_fwd(ptr0,info[1],_diffs,row);
#endif

  } else {
    delta = 63 - remainder;
    column = get_column(delta);
    row = get_row(delta);
    debug(printf("quarter-block %d, delta %d, column %d, row %d\n",quarter_block,delta,column,row));

    (unpacker_table[packsize_div2][column*4 + quarter_block])(diffs,bitpack);
#ifdef BRANCH_FREE_ROW_SUM
    ptr0 = info[PAIRED_METAINFO_SIZE+1] - _diffs[row+1] - _diffs[row+2] - _diffs[row+3] - _diffs[row+4];
#else
    assign_sum_rev(ptr0,info[PAIRED_METAINFO_SIZE+1],_diffs,row);
#endif
  }

  *end0 = ptr0 + Bitpack64_access_dispatch(bitpack,/*nwritten*/info[PAIRED_METAINFO_SIZE]-info[2],
					   remainder);

  return ptr0;

#else  /* HAVE_SSE2 */

  offset0 = info[1];
  offset1 = info[DIFFERENTIAL_METAINFO_SIZE+1];

  /* Unpack all 64 diffs for non-SIMD */
  (unpacker_all_table[packsize_div2*2])(&(diffs[1]),bitpack);

#ifdef DEBUG
  printf("oligo: %08X, remainder %d, offset0 %u, offset1 %u\n",
	 oligo,oligo % BLOCKSIZE,info[1],info[DIFFERENTIAL_METAINFO_SIZE+1]);
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

#endif	/* HAVE_SSE2 */

}
#endif



/* Needed for poly-T to avoid computing on metablock after the last
   one to find end0 */
UINT4
Bitpack64_read_one (Storedoligomer_T oligo, UINT4 *bitpackptrs, UINT4 *bitpackcomp) {
  Storedoligomer_T bmer;
  UINT4 *info, nwritten, packsize_div2;
  int delta, remainder, quarter_block, column, row;
#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
  UINT4 ptr;
  UINT4 diffs[BLOCKSIZE+1], *bitpack;
  int k, i;
#else
#ifdef BRANCH_FREE_ROW_SUM
  __m128i diffs[3];
#else
  __m128i diffs[2];
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
  UINT4 psums[4];
#endif
  __m128i *bitpack;
  UINT4 *_diffs;
#endif
#ifdef DEBUG
  UINT4 offsets[BLOCKSIZE+1];
#endif


  bmer = oligo/BLOCKSIZE;
  info = &(bitpackptrs[bmer * DIFFERENTIAL_METAINFO_SIZE]);

  debug(printf("Entered Bitpack64_read_one with oligo %u => bmer %u\n",oligo,bmer));

#if defined(WORDS_BIGENDIAN)
  nwritten = Bigendian_convert_uint(info[0]);		/* In 128-bit registers */
  bitpack = (UINT4 *) &(bitpackcomp[nwritten*4]);
  packsize_div2 = (Bigendian_convert_uint(info[DIFFERENTIAL_METAINFO_SIZE]) - nwritten);

#elif !defined(HAVE_SSE2)
  nwritten = info[0];		/* In 128-bit registers */
  bitpack = (UINT4 *) &(bitpackcomp[nwritten*4]);
  packsize_div2 = (info[DIFFERENTIAL_METAINFO_SIZE] - nwritten);

#else
  nwritten = info[0];		/* In 128-bit registers */
  bitpack = (__m128i *) &(bitpackcomp[nwritten*4]);
  /* packsize = (info[DIFFERENTIAL_METAINFO_SIZE] - nwritten)*2; */
  packsize_div2 = (info[DIFFERENTIAL_METAINFO_SIZE] - nwritten);
#endif

  remainder = oligo % BLOCKSIZE;
  quarter_block = remainder / 16;

#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)

  /* Unpack all 64 diffs for non-SIMD */
  (unpacker_all_table[packsize_div2*2])(&(diffs[1]),bitpack);

  if (remainder <= 16) {
#ifdef WORDS_BIGENDIAN
    ptr = Bigendian_convert_uint(/*offset0*/info[1]);
#else
    ptr = /*offset0*/info[1];
#endif

    delta = remainder - 1;
    column = get_column(delta);
    row = get_row(delta);
    debug(printf("column %d, row %d\n",column,row));
    
    for (k = column*2 + 1, i = 0; i <= row; k += BLOCKSIZE/4, i++) {
      debug(printf("Adding diffs[%d] = %u\n",k,diffs[k]));
      ptr += diffs[k];
    }

  } else if (remainder <= 32) {
#ifdef WORDS_BIGENDIAN
    ptr = Bigendian_convert_uint(/*offset0*/info[1]);
#else
    ptr = /*offset0*/info[1];
#endif

    delta = remainder - 1;
    column = get_column(delta);
    row = get_row(delta);
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
#ifdef WORDS_BIGENDIAN
    ptr = Bigendian_convert_uint(/*offset1*/info[DIFFERENTIAL_METAINFO_SIZE+1]);
#else
    ptr = /*offset1*/info[DIFFERENTIAL_METAINFO_SIZE+1];
#endif

    delta = 63 - remainder;
    column = get_column(delta);
    row = get_row(delta);
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
#ifdef WORDS_BIGENDIAN
    ptr = Bigendian_convert_uint(/*offset1*/info[DIFFERENTIAL_METAINFO_SIZE+1]);
#else
    ptr = /*offset1*/info[DIFFERENTIAL_METAINFO_SIZE+1];
#endif

    delta = 63 - remainder;
    column = get_column(delta);
    row = get_row(delta);
    debug(printf("column %d, row %d\n",column,row));

    for (k = column*2 + 9, i = 0; i <= row; k += BLOCKSIZE/4, i++) {
      debug(printf("Subtracting diffs[%d] = %u\n",k,diffs[k]));
      ptr -= diffs[k];
    }
  }

  return ptr;

#else  /* littleendian and SSE2 */
  _diffs = (UINT4 *) diffs;	/* Assumes a dummy register in diffs[0] */

#ifdef BRANCH_FREE_QTR_BLOCK
  psums[0] = psums[1] = info[1];
  psums[2] = psums[3] = info[DIFFERENTIAL_METAINFO_SIZE+1];

  delta = 31 - abs(remainder - 32);
  column = get_column(delta);
  row = get_row(delta);
  debug(printf("quarter-block %d, delta %d, column %d, row %d\n",quarter_block,delta,column,row));

  (unpacker_table[packsize_div2][column*4 + quarter_block])(diffs,bitpack);
  return psums[quarter_block] + _diffs[row+1] + _diffs[row+2] + _diffs[row+3] + _diffs[row+4];

#else

  if (quarter_block <= 1) {
    delta = remainder - 1;
    column = get_column(delta);
    row = get_row(delta);
    debug(printf("quarter-block %d, delta %d, column %d, row %d\n",quarter_block,delta,column,row));

    (unpacker_table[packsize_div2][column*4 + quarter_block])(diffs,bitpack);
#ifdef BRANCH_FREE_ROW_SUM
    return info[1] + _diffs[row+1] + _diffs[row+2] + _diffs[row+3] + _diffs[row+4];
#else
    return_sum_fwd(info[1],_diffs,row);
#endif

  } else {
    delta = 63 - remainder;
    column = get_column(delta);
    row = get_row(delta);
    debug(printf("quarter-block %d, delta %d, column %d, row %d\n",quarter_block,delta,column,row));

    (unpacker_table[packsize_div2][column*4 + quarter_block])(diffs,bitpack);
#ifdef BRANCH_FREE_ROW_SUM
    return info[DIFFERENTIAL_METAINFO_SIZE+1] - _diffs[row+1] - _diffs[row+2] - _diffs[row+3] - _diffs[row+4];
#else
    return_sum_rev(info[DIFFERENTIAL_METAINFO_SIZE+1],_diffs,row);
#endif
  }

#endif	/* BRANCH_FREE_QTR_BLOCK */
#endif	/* littleendian and SSE2 */
}


/* Needed for poly-T to avoid computing on metablock after the last
   one to find end0 */
UINT8
Bitpack64_read_one_huge (Storedoligomer_T oligo, UINT4 *bitpackpages,
			 UINT4 *bitpackptrs, UINT4 *bitpackcomp) {
  UINT4 *pageptr;
  Storedoligomer_T bmer;
  UINT4 *info, nwritten, packsize_div2;
  UINT8 offset0, offset1;
  int delta, remainder, quarter_block, column, row;
#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
  UINT8 ptr;
  UINT4 diffs[BLOCKSIZE+1], *bitpack;
  int k;
#else
#ifdef BRANCH_FREE_ROW_SUM
  __m128i diffs[3];
#else
  __m128i diffs[2];
#endif
#ifdef BRANCH_FREE_QTR_BLOCK
  UINT8 psums[4];
#endif
  __m128i *bitpack;
  UINT4 *_diffs;
#endif
  int i;


  bmer = oligo/BLOCKSIZE;
  info = &(bitpackptrs[bmer * DIFFERENTIAL_METAINFO_SIZE]);

  debug(printf("Entered Bitpack64_read_one_huge with oligo %u => bmer %u\n",oligo,bmer));

#ifdef WORDS_BIGENDIAN
  nwritten = Bigendian_convert_uint(info[0]);		/* In 128-bit registers */
  bitpack = (UINT4 *) &(bitpackcomp[nwritten*4]);
  packsize_div2 = (Bigendian_convert_uint(info[DIFFERENTIAL_METAINFO_SIZE]) - nwritten);

#elif !defined(HAVE_SSE2)
  nwritten = info[0];		/* In 128-bit registers */
  bitpack = (UINT4 *) &(bitpackcomp[nwritten*4]);
  packsize_div2 = (info[DIFFERENTIAL_METAINFO_SIZE] - nwritten);

#else
  nwritten = info[0];		/* In 128-bit registers */
  bitpack = (__m128i *) &(bitpackcomp[nwritten*4]);
  /* packsize = (info[DIFFERENTIAL_METAINFO_SIZE] - nwritten)*2; */
  packsize_div2 = (info[DIFFERENTIAL_METAINFO_SIZE] - nwritten);
#endif


#ifdef DEBUG
  printf("bitpack (for packsize %d):\n",packsize_div2*2);
  for (i = 0; i < packsize_div2; i++) {
    print_vector_hex(bitpack[i]);
  }
  printf("\n");
#endif  

  remainder = oligo % BLOCKSIZE;
  quarter_block = remainder / 16;

#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)

  /* Unpack all 64 diffs for non-SIMD */
  (unpacker_all_table[packsize_div2*2])(&(diffs[1]),bitpack);

  if ((remainder = oligo % BLOCKSIZE) == 0) {
#ifdef WORDS_BIGENDIAN
    ptr = Bigendian_convert_uint(/*offset0*/info[1]);
#else
    ptr = /*offset0*/info[1];
#endif

  } else if (remainder <= 16) {
#ifdef WORDS_BIGENDIAN
    ptr = Bigendian_convert_uint(/*offset0*/info[1]);
    if (bitpackpages != NULL) {
      pageptr = bitpackpages;
      while (bmer+1 >= Bigendian_convert_uint(*pageptr)) {
	ptr += POSITIONS_PAGE;
	pageptr++;
      }
    }
#else
    ptr = /*offset0*/info[1];
    if (bitpackpages != NULL) {
      pageptr = bitpackpages;
      while (bmer+1 >= *pageptr) {
	ptr += POSITIONS_PAGE;
	pageptr++;
      }
    }
#endif

    column = (remainder - 1) % 4; /* Goes from 0 to 3 */
    row = (remainder - 1) / 4;
    debug(printf("column %d, row %d\n",column,row));
    
    for (k = column*2 + 1, i = 0; i <= row; k += BLOCKSIZE/4, i++) {
      debug(printf("Adding diffs[%d] = %u\n",k,diffs[k]));
      ptr += diffs[k];
    }

  } else if (remainder <= 32) {
#ifdef WORDS_BIGENDIAN
    ptr = Bigendian_convert_uint(/*offset0*/info[1]);
    if (bitpackpages != NULL) {
      pageptr = bitpackpages;
      while (bmer+1 >= Bigendian_convert_uint(*pageptr)) {
	ptr += POSITIONS_PAGE;
	pageptr++;
      }
    }
#else
    ptr = /*offset0*/info[1];
    if (bitpackpages != NULL) {
      pageptr = bitpackpages;
      while (bmer+1 >= *pageptr) {
	ptr += POSITIONS_PAGE;
	pageptr++;
      }
    }
#endif

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
#ifdef WORDS_BIGENDIAN
    ptr = Bigendian_convert_uint(/*offset1*/info[DIFFERENTIAL_METAINFO_SIZE+1]);
    if (bitpackpages != NULL) {
      pageptr = bitpackpages;
      while (bmer+1 >= Bigendian_convert_uint(*pageptr)) {
	ptr += POSITIONS_PAGE;
	pageptr++;
      }
    }
#else
    ptr = /*offset1*/info[DIFFERENTIAL_METAINFO_SIZE+1];
    if (bitpackpages != NULL) {
      pageptr = bitpackpages;
      while (bmer+1 >= *pageptr) {
	ptr += POSITIONS_PAGE;
	pageptr++;
      }
    }
#endif

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
#ifdef WORDS_BIGENDIAN
    ptr = Bigendian_convert_uint(/*offset1*/info[DIFFERENTIAL_METAINFO_SIZE+1]);
    if (bitpackpages != NULL) {
      pageptr = bitpackpages;
      while (bmer+1 >= Bigendian_convert_uint(*pageptr)) {
	ptr += POSITIONS_PAGE;
	pageptr++;
      }
    }
#else
    ptr = /*offset1*/info[DIFFERENTIAL_METAINFO_SIZE+1];
    if (bitpackpages != NULL) {
      pageptr = bitpackpages;
      while (bmer+1 >= *pageptr) {
	ptr += POSITIONS_PAGE;
	pageptr++;
      }
    }
#endif

    column = (63 - remainder) % 4; /* Goes from 0 to 3.  Assert remainder < 64 */
    row = (63 - remainder) / 4;
    debug(printf("column %d, row %d\n",column,row));

    for (k = column*2 + 9, i = 0; i <= row; k += BLOCKSIZE/4, i++) {
      debug(printf("Subtracting diffs[%d] = %u\n",k,diffs[k]));
      ptr -= diffs[k];
    }
  }

  return ptr;

#else			    /* littleendian and SSE2 */
  _diffs = (UINT4 *) diffs;	/* Assumes a dummy register in diffs[0] */

#ifdef BRANCH_FREE_QTR_BLOCK
  if (bitpackpages == NULL) {
    offset0 = info[1];
    offset1 = info[DIFFERENTIAL_METAINFO_SIZE+1];
  } else {
    offset0 = 0UL;
    pageptr = bitpackpages;
    while (bmer >= *pageptr) {
      offset0 += POSITIONS_PAGE;
      pageptr++;
    }

    offset1 = offset0;
    if (bmer + 1 >= *pageptr) {
      offset1 += POSITIONS_PAGE;
      pageptr++;
    }

    offset0 += info[1];
    offset1 += info[DIFFERENTIAL_METAINFO_SIZE+1];
  }
  debug(printf("offset0 = %u, offset1 = %u\n",offset0,offset1));

  psums[0] = psums[1] = offset0;
  psums[2] = psums[3] = offset1;

  delta = 31 - abs(remainder - 32);
  column = get_column(delta);
  row = get_row(delta);
  debug(printf("quarter-block %d, delta %d, column %d, row %d\n",quarter_block,delta,column,row));

  (unpacker_table[packsize_div2][column*4 + quarter_block])(diffs,bitpack);
#ifdef DEBUG
  printf("%d %d %d %d\n",_diffs[0],_diffs[1],_diffs[2],_diffs[3]);
  printf("%d %d %d %d\n",_diffs[4],_diffs[5],_diffs[6],_diffs[7]);
  printf("%d %d %d %d\n",_diffs[8],_diffs[9],_diffs[10],_diffs[11]);
  printf("%d %d %d %d\n",_diffs[12],_diffs[13],_diffs[14],_diffs[15]);
#endif

  debug(printf("Returning %u + %d + %d + %d + %d\n",
	       psums[quarter_block],_diffs[row+1],_diffs[row+2],_diffs[row+3],_diffs[row+4]));
  return psums[quarter_block] + (INT4) (_diffs[row+1] + _diffs[row+2] + _diffs[row+3] + _diffs[row+4]);

#else

  if (quarter_block <= 1) {
    offset0 = (UINT8) info[1];
    if (bitpackpages != NULL) {
      pageptr = bitpackpages;
      while (bmer >= *pageptr) {
	offset0 += POSITIONS_PAGE;
	pageptr++;
      }
    }

    delta = remainder - 1;
    column = get_column(delta);
    row = get_row(delta);
    debug(printf("quarter-block %d, delta %d, column %d, row %d\n",quarter_block,delta,column,row));

    (unpacker_table[packsize_div2][column*4 + quarter_block])(diffs,bitpack);
#ifdef DEBUG
    printf("%u %u %u %u\n",_diffs[0],_diffs[1],_diffs[2],_diffs[3]);
    printf("%u %u %u %u\n",_diffs[4],_diffs[5],_diffs[6],_diffs[7]);
    printf("%u %u %u %u\n",_diffs[8],_diffs[9],_diffs[10],_diffs[11]);
    printf("%u %u %u %u\n",_diffs[12],_diffs[13],_diffs[14],_diffs[15]);
#endif

#ifdef BRANCH_FREE_ROW_SUM
    return offset0 + _diffs[row+1] + _diffs[row+2] + _diffs[row+3] + _diffs[row+4];
#else
    return_sum_fwd(offset0,_diffs,row);
#endif

  } else {
    offset1 = (UINT8) info[DIFFERENTIAL_METAINFO_SIZE+1];
    if (bitpackpages != NULL) {
      pageptr = bitpackpages;
      while (bmer + 1 >= *pageptr) {
	offset1 += POSITIONS_PAGE;
	pageptr++;
      }
    }

    delta = 63 - remainder;
    column = get_column(delta);
    row = get_row(delta);
    debug(printf("quarter-block %d, delta %d, column %d, row %d\n",quarter_block,delta,column,row));

    (unpacker_table[packsize_div2][column*4 + quarter_block])(diffs,bitpack);
#ifdef DEBUG
    printf("%u %u %u %u\n",_diffs[0],_diffs[1],_diffs[2],_diffs[3]);
    printf("%u %u %u %u\n",_diffs[4],_diffs[5],_diffs[6],_diffs[7]);
    printf("%u %u %u %u\n",_diffs[8],_diffs[9],_diffs[10],_diffs[11]);
    printf("%u %u %u %u\n",_diffs[12],_diffs[13],_diffs[14],_diffs[15]);
#endif

#ifdef BRANCH_FREE_ROW_SUM
    return offset1 - _diffs[row+1] - _diffs[row+2] - _diffs[row+3] - _diffs[row+4];
#else
    return_sum_rev(offset1,_diffs,row);
#endif
  }

#endif	/* BRANCH_FREE_QTR_BLOCK */
#endif	/* littleendian and SSE2 */
}


#ifndef PMAP
/* Unpack all offsets.  Can treat offset0 as a special case */
void
Bitpack64_block_offsets (UINT4 *offsets, Storedoligomer_T oligo,
			 UINT4 *bitpackptrs, UINT4 *bitpackcomp) {
  UINT4 *info, nwritten;
  UINT4 offset0, offset1, temp;
  int packsize, k;
#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
  int column, row;
  UINT4 diffs[BLOCKSIZE], columnar[BLOCKSIZE], *bitpack, *vertical;
#else
  __m128i diffs[8], *bitpack;
  UINT4 *_diffs;
#endif
#ifdef DEBUG
  int i;
#endif


  info = &(bitpackptrs[oligo/BLOCKSIZE * DIFFERENTIAL_METAINFO_SIZE]);
#ifdef WORDS_BIGENDIAN
  nwritten = Bigendian_convert_uint(info[0]);
  bitpack = (UINT4 *) &(bitpackcomp[nwritten*4]);
  offset0 = Bigendian_convert_uint(info[1]);
  offset1 = Bigendian_convert_uint(info[DIFFERENTIAL_METAINFO_SIZE+1]);
  packsize = (Bigendian_convert_uint(info[DIFFERENTIAL_METAINFO_SIZE]) - nwritten)*2;

#elif !defined(HAVE_SSE2)
  nwritten = info[0];		/* In 128-bit registers */
  bitpack = (UINT4 *) &(bitpackcomp[nwritten*4]);
  offset0 = info[1];
  offset1 = info[DIFFERENTIAL_METAINFO_SIZE+1];
  packsize = (info[DIFFERENTIAL_METAINFO_SIZE] - nwritten)*2;

#else
  nwritten = info[0];		/* In 128-bit registers */
  bitpack = (__m128i *) &(bitpackcomp[nwritten*4]);
  offset0 = info[1];
  offset1 = info[DIFFERENTIAL_METAINFO_SIZE+1];
  packsize = (info[DIFFERENTIAL_METAINFO_SIZE] - nwritten)*2;
#endif


#ifdef DEBUG
  printf("oligo: %08X, nwritten %u, offset0 %u, offset1 %u, packsize %d\n",
	 oligo,nwritten,offset0,offset1,packsize);
#endif


#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
  /* Unpack all 64 diffs for non-SIMD */
  (unpacker_all_table[packsize])(&(diffs[0]),bitpack);

#ifdef DEBUG
  for (i = 0; i < BLOCKSIZE i += 16) {
    printf("%u %u %u %u ",diffs[i],diffs[i+1],diffs[i+2],diffs[i+3]);
    printf("%u %u %u %u ",diffs[i+4],diffs[i+5],diffs[i+6],diffs[i+7]);
    printf("%u %u %u %u ",diffs[i+8],diffs[i+9],diffs[i+10],diffs[i+11]);
    printf("%u %u %u %u\n",diffs[i+12],diffs[i+13],diffs[i+14],diffs[i+15]);
  }
  printf("end of diffs horizontal (because non-SIMD unpackers are horizontal)\n");
#endif

  /* Convert to columnar */
  vertical = &(diffs[0]);
  for (column = 0; column < 4; column++) {
    k = column;
    for (row = 0; row < BLOCKSIZE/4; row++) {
      columnar[k] = *vertical++;
      k += 4;
    }
  }

  /* Convert to vertical, shifting by 1 */
  vertical_order(&(offsets[1]),columnar);

#ifdef DEBUG
  printf("%u\n",offset0);
  for (i = 1; i <= BLOCKSIZE; i += 4) {
    printf("%u %u %u %u\n",offsets[i],offsets[i+1],offsets[i+2],offsets[i+3]);
  }
  printf("end of diffs vertical\n");
#endif

#else  /* littleendian and SSE2 */

#ifdef DEBUG
  printf("bitpack:\n");
  for (i = 0; i < packsize/2; i++) {
    print_vector_hex(bitpack[i]);
  }
  printf("\n");
#endif  

  _diffs = (UINT4 *) &(diffs[0]);

  /* Unpack fwd 32 cumulative sums under SIMD */
  (unpacker_all_table[packsize])(&(diffs[0]),bitpack);
  vertical_order_fwd(&(offsets[1]),_diffs);

  /* Unpack rev 32 cumulative sums under SIMD */
  (unpacker_all_table[packsize+1])(&(diffs[0]),bitpack);
  vertical_order_rev(&(offsets[33]),_diffs);

#ifdef DEBUG
  printf("%u\n",offsets[i]);
  for (i = 1; i <= BLOCKSIZE; i += 4) {
    printf("%u %u %u %u\n",offsets[i],offsets[i+1],offsets[i+2],offsets[i+3]);
  }
  printf("end of diffs vertical\n");
#endif

#endif	/* littleendian and SSE2 */

  /* Perform cumulative sum */
  offsets[0] = offset0;
  offsets[1] += offset0;
  offsets[2] += offset0;
  offsets[3] += offset0;
  offsets[4] += offset0;
  for (k = 5; k <= 32; k++) {
    offsets[k] += offsets[k-4];
  }

  /* Skip offsets[33] through offsets[36] */
  for (k = 37; k <= BLOCKSIZE; k++) {
    offsets[k] += offsets[k-4];
  }

  /* Now swap offsets */
  for (k = 33; k <= 48; k++) {
    temp = offsets[96-k];
    offsets[96-k] = offset1 - offsets[k];
    offsets[k] = offset1 - temp;
  }
  offsets[64] = offset1;


#ifdef DEBUG
  printf("%u\n",offsets[0]);
  for (i = 1; i <= 32; i += 4) {
    printf("%u %u %u %u\n",offsets[i],offsets[i+1],offsets[i+2],offsets[i+3]);
  }
  printf("\n");
  for (i = 33; i <= 64; i += 4) {
    printf("%u %u %u %u\n",offsets[i],offsets[i+1],offsets[i+2],offsets[i+3]);
  }
  printf("end of offsets\n");
#endif

  return;
}
#endif


#ifndef PMAP
#if defined(HAVE_64_BIT) && (defined(UTILITYP) || defined(LARGE_GENOMES))
/* Unpack all offsets.  Can treat offset0 as a special case */
void
Bitpack64_block_offsets_huge (UINT8 *offsets, Storedoligomer_T oligo,
			      UINT4 *bitpackpages, UINT4 *bitpackptrs, UINT4 *bitpackcomp) {
  UINT4 *pageptr;
  UINT4 *info, nwritten;
  Storedoligomer_T bmer;
  UINT8 offset0, offset1, temp;
  int packsize, k;
#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
  int column, row;
  UINT4 diffs[BLOCKSIZE], columnar[BLOCKSIZE], *bitpack, *vertical;
#else
  __m128i diffs[8], *bitpack;
  UINT4 *_diffs;
#endif
#ifdef DEBUG
  int i;
#endif

  bmer = oligo/BLOCKSIZE;

  info = &(bitpackptrs[bmer * DIFFERENTIAL_METAINFO_SIZE]);

#ifdef WORDS_BIGENDIAN
  nwritten = Bigendian_convert_uint(info[0]); /* In 128-bit registers */
  bitpack = (UINT4 *) &(bitpackcomp[nwritten*4]);

#elif !defined(HAVE_SSE2)
  nwritten = info[0];		/* In 128-bit registers */
  bitpack = (UINT4 *) &(bitpackcomp[nwritten*4]);

#else
  nwritten = info[0];		/* In 128-bit registers */
  bitpack = (__m128i *) &(bitpackcomp[nwritten*4]);
#endif

#ifdef WORDS_BIGENDIAN
  offset0 = offset1 = 0UL;
  pageptr = bitpackpages;
  while (bmer >= Bigendian_convert_uint(*pageptr)) {
    offset0 += POSITIONS_PAGE;
    pageptr++;
  }

  offset1 = offset0;
  if (bmer+1 >= Bigendian_convert_uint(*pageptr)) {
    offset1 += POSITIONS_PAGE;
  }
#else
  offset0 = offset1 = 0UL;
  pageptr = bitpackpages;
  while (bmer >= *pageptr) {
    offset0 += POSITIONS_PAGE;
    pageptr++;
  }

  offset1 = offset0;
  if (bmer+1 >= *pageptr) {
    offset1 += POSITIONS_PAGE;
  }
#endif


#ifdef WORDS_BIGENDIAN
  offset0 += Bigendian_convert_uint(info[1]);
  offset1 += Bigendian_convert_uint(info[DIFFERENTIAL_METAINFO_SIZE+1]);
  packsize = (Bigendian_convert_uint(info[DIFFERENTIAL_METAINFO_SIZE]) - nwritten)*2;
#else
  offset0 += info[1];
  offset1 += info[DIFFERENTIAL_METAINFO_SIZE+1];
  packsize = (info[DIFFERENTIAL_METAINFO_SIZE] - nwritten)*2;
#endif


#ifdef DEBUG
  printf("oligo: %08X, nwritten %u, offset0 %u, offset1 %u, packsize %d\n",
	 oligo,nwritten,offset0,offset1,packsize);
#endif


#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
  /* Unpack all 64 diffs for non-SIMD */
  (unpacker_all_table[packsize])(&(diffs[0]),bitpack);

#ifdef DEBUG
  for (i = 0; i < 64; i += 16) {
    printf("%u %u %u %u ",diffs[i],diffs[i+1],diffs[i+2],diffs[i+3]);
    printf("%u %u %u %u ",diffs[i+4],diffs[i+5],diffs[i+6],diffs[i+7]);
    printf("%u %u %u %u ",diffs[i+8],diffs[i+9],diffs[i+10],diffs[i+11]);
    printf("%u %u %u %u\n",diffs[i+12],diffs[i+13],diffs[i+14],diffs[i+15]);
  }
  printf("end of diffs horizontal (because non-SIMD unpackers are horizontal)\n");
#endif

  /* Convert to columnar */
  vertical = &(diffs[0]);
  for (column = 0; column < 4; column++) {
    k = column;
    for (row = 0; row < 16; row++) {
      columnar[k] = (UINT8) *vertical++;
      k += 4;
    }
  }

  /* Convert to vertical, shifting by 1 */
  vertical_order_huge(&(offsets[1]),columnar);

#ifdef DEBUG
  printf("%u\n",offset0);
  for (i = 1; i <= 64; i += 4) {
    printf("%u %u %u %u\n",offsets[i],offsets[i+1],offsets[i+2],offsets[i+3]);
  }
  printf("end of diffs vertical\n");
#endif


#else
#ifdef DEBUG
  printf("bitpack:\n");
  for (i = 0; i < packsize/2; i++) {
    print_vector_hex(bitpack[i]);
  }
  printf("\n");
#endif

  _diffs = (UINT4 *) &(diffs[0]);

  /* Unpack fwd 32 cumulative sums under SIMD */
  (unpacker_all_table[packsize])(&(diffs[0]),bitpack);
  vertical_order_huge_fwd(&(offsets[1]),_diffs);

  /* Unpack rev 32 cumulative sums under SIMD */
  (unpacker_all_table[packsize+1])(&(diffs[0]),bitpack);
  vertical_order_huge_rev(&(offsets[33]),_diffs);

#ifdef DEBUG
  printf("%u\n",offsets[i]);
  for (i = 1; i <= 64; i += 4) {
    printf("%u %u %u %u\n",offsets[i],offsets[i+1],offsets[i+2],offsets[i+3]);
  }
  printf("end of diffs vertical\n");
#endif

#endif	/* HAVE_SSE2 */

  /* Perform cumulative sum */
  offsets[0] = offset0;
  offsets[1] += offset0;
  offsets[2] += offset0;
  offsets[3] += offset0;
  offsets[4] += offset0;
  for (k = 5; k <= 32; k++) {
    offsets[k] += offsets[k-4];
  }

  /* Skip offsets[33] through offsets[36] */
  for (k = 37; k <= 64; k++) {
    offsets[k] += offsets[k-4];
  }

  /* Now swap offsets */
  for (k = 33; k <= 48; k++) {
    temp = offsets[96-k];
    offsets[96-k] = offset1 - offsets[k];
    offsets[k] = offset1 - temp;
  }
  offsets[64] = offset1;


#ifdef DEBUG
  printf("%u\n",offsets[0]);
  for (i = 1; i <= 32; i += 4) {
    printf("%u %u %u %u\n",offsets[i],offsets[i+1],offsets[i+2],offsets[i+3]);
  }
  printf("\n");
  for (i = 33; i <= 64; i += 4) {
    printf("%u %u %u %u\n",offsets[i],offsets[i+1],offsets[i+2],offsets[i+3]);
  }
  printf("end of offsets\n");
#endif

  return;
}
#endif
#endif
