static char rcsid[] = "$Id: compress.c 168395 2015-06-26 17:13:13Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
# define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
# define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#include "compress.h"
#include "compress-write.h"


#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <ctype.h>		/* For isalpha, toupper */
#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
#else
#include "littleendian.h"
#endif
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>		/* For off_t */
#endif
#include "complement.h"
#include "assert.h"
#include "mem.h"		/* For Compress_new */
#include "assert.h"

#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
/* Skip */
#else
#include <emmintrin.h>
#endif
#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSSE3)
/* Skip */
#else
#include <tmmintrin.h>
#endif
#ifdef HAVE_SSE4_1
#include <smmintrin.h>
#endif


#ifdef DEBUG0
#define debug0(x) x
#else
#define debug0(x)
#endif


#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* fragment_left and fragment_right */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif


/* SIMD */
#ifdef DEBUG9
#define debug9(x) x
#else
#define debug9(x)
#endif

/* Checking SSSE3 procedure */
#ifdef DEBUG14
#define debug14(x) x
#else
#define debug14(x)
#endif


#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
#define STEP_SIZE 32
#else
#define STEP_SIZE 128
#endif


#define T Compress_T
struct T {
  Genomecomp_T *blocks;
  int nblocks;
  Genomecomp_T **shift_array;
  bool availp[STEP_SIZE];
#ifdef DEBUG14
  int querylength;
#endif
};


void
Compress_free (T *old) {
  if (*old) {
#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
    FREE((*old)->shift_array[0]);
#else
    _mm_free((*old)->shift_array[0]);
#endif
    FREE((*old)->shift_array);
#if 0
#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
    FREE((*old)->blocks);
#else
    _mm_free((*old)->blocks);
#endif
#endif
    FREE(*old);
  }
  return;
}
void
Compress_print (T this) {
  int ptr = 0;

  while (ptr < this->nblocks*COMPRESS_BLOCKSIZE) {
    printf("high: %08X  low: %08X  flags: %08X\n",
	   this->blocks[ptr],this->blocks[ptr+1],this->blocks[ptr+2]);
    ptr += COMPRESS_BLOCKSIZE;
  }
  printf("\n");
  return;
}


int
Compress_nblocks (T this) {
  return this->nblocks;
}


static void
write_chars (Genomecomp_T high, Genomecomp_T low, Genomecomp_T flags) {
  char Buffer[33];
  int i;

  Buffer[32] = '\0';
  /* printf("%08X %08X %08X => ",high,low,flags); */

  for (i = 0; i < 32; i++) {
    switch (((high & 0x01) << 1) | (low & 0x01)) {
    case 0U: Buffer[i] = 'A'; break;
    case 1U: Buffer[i] = 'C'; break;
    case 2U: Buffer[i] = 'G'; break;
    case 3U: Buffer[i] = 'T'; break;
    default: abort();
    }
    high >>= 1;
    low >>= 1;
  }

  if (flags != 0U) {
    for (i = 0; i < 32; i++) {
      if (flags & 0x01) {
	Buffer[i] = 'N';
      }
      flags >>= 1;
    }
  }

  printf("%s",Buffer);
  return;
}


#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
void
Compress_print_blocks (Genomecomp_T *blocks, int nshift, int pos5, int pos3) {
  int ptr, endptr;

  endptr = (nshift + pos3)/32U*3;   /* /STEP_SIZE*COMPRESS_BLOCKSIZE */
  ptr = (nshift + pos5)/32U*3;

  while (ptr <= endptr) {
    printf("high: %08X  low: %08X  flags: %08X\t",
	   blocks[ptr],blocks[ptr+1],blocks[ptr+2]);
    write_chars(blocks[ptr],blocks[ptr+1],blocks[ptr+2]);
    printf("\n");
    ptr += COMPRESS_BLOCKSIZE;
  }
  printf("\n");
  return;
}

#else
void
Compress_print_blocks (Genomecomp_T *blocks, int nshift, int pos5, int pos3) {
  int ptr, endptr;
  int startcolumni, endcolumni;


  startcolumni = ((nshift + pos5) % 128) / 32;
  endcolumni = ((nshift + pos3) % 128) / 32;

  endptr = (nshift + pos3)/128U*12;
  ptr = (nshift + pos5)/128U*12;

  while (ptr < endptr) {
    printf("high: %08X  low: %08X  flags: %08X\t",
	   blocks[ptr],blocks[ptr+4],blocks[ptr+8]);
    write_chars(blocks[ptr],blocks[ptr+4],blocks[ptr+8]);
    printf("\n");

    printf("high: %08X  low: %08X  flags: %08X\t",
	   blocks[ptr+1],blocks[ptr+5],blocks[ptr+9]);
    write_chars(blocks[ptr+1],blocks[ptr+5],blocks[ptr+9]);
    printf("\n");

    printf("high: %08X  low: %08X  flags: %08X\t",
	   blocks[ptr+2],blocks[ptr+6],blocks[ptr+10]);
    write_chars(blocks[ptr+2],blocks[ptr+6],blocks[ptr+10]);
    printf("\n");

    printf("high: %08X  low: %08X  flags: %08X\t",
	   blocks[ptr+3],blocks[ptr+7],blocks[ptr+11]);
    write_chars(blocks[ptr+3],blocks[ptr+7],blocks[ptr+11]);
    printf("\n");

    printf("\n");

    ptr += COMPRESS_BLOCKSIZE;
  }

  /* Last block */
  printf("high: %08X  low: %08X  flags: %08X\t",
	 blocks[ptr],blocks[ptr+4],blocks[ptr+8]);
  write_chars(blocks[ptr],blocks[ptr+4],blocks[ptr+8]);
  printf("\n");

  if (0 && endcolumni == 0) {
    printf("\n");
    return;
  }

  printf("high: %08X  low: %08X  flags: %08X\t",
	 blocks[ptr+1],blocks[ptr+5],blocks[ptr+9]);
  write_chars(blocks[ptr+1],blocks[ptr+5],blocks[ptr+9]);
  printf("\n");

  if (0 && endcolumni == 1) {
    printf("\n");
    return;
  }

  printf("high: %08X  low: %08X  flags: %08X\t",
	 blocks[ptr+2],blocks[ptr+6],blocks[ptr+10]);
  write_chars(blocks[ptr+2],blocks[ptr+6],blocks[ptr+10]);
  printf("\n");

  if (0 && endcolumni == 2) {
    printf("\n");
    return;
  }

  printf("high: %08X  low: %08X  flags: %08X\t",
	 blocks[ptr+3],blocks[ptr+7],blocks[ptr+11]);
  write_chars(blocks[ptr+3],blocks[ptr+7],blocks[ptr+11]);
  printf("\n");

  printf("\n");
  return;
}


void
Compress_print_one_block (Genomecomp_T *blocks) {
  int ptr = 0;

  printf("high: %08X  low: %08X  flags: %08X\t",
	 blocks[ptr],blocks[ptr+4],blocks[ptr+8]);
  write_chars(blocks[ptr],blocks[ptr+4],blocks[ptr+8]);
  printf("\n");
  
  printf("high: %08X  low: %08X  flags: %08X\t",
	 blocks[ptr+1],blocks[ptr+5],blocks[ptr+9]);
  write_chars(blocks[ptr+1],blocks[ptr+5],blocks[ptr+9]);
  printf("\n");

  printf("high: %08X  low: %08X  flags: %08X\t",
	 blocks[ptr+2],blocks[ptr+6],blocks[ptr+10]);
  write_chars(blocks[ptr+2],blocks[ptr+6],blocks[ptr+10]);
  printf("\n");

  printf("high: %08X  low: %08X  flags: %08X\t",
	 blocks[ptr+3],blocks[ptr+7],blocks[ptr+11]);
  write_chars(blocks[ptr+3],blocks[ptr+7],blocks[ptr+11]);
  printf("\n");

  printf("\n");
  return;
}

#endif



/*                   87654321 */
#define LEFT_SET   0x80000000
#define LEFT_CLEAR 0x00000000


T
Compress_new_fwd (char *gbuffer, Chrpos_T length) {
  T new = (T) MALLOC(sizeof(*new));
  Genomecomp_T high, low, flags;
  Chrpos_T ptr;
  Chrpos_T position;
  int c, i;
  int in_counter = 0;

#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
  new->nblocks = (length+31)/32U;
  new->shift_array = (Genomecomp_T **) MALLOC(STEP_SIZE * sizeof(Genomecomp_T *));
  new->shift_array[0] = (Genomecomp_T *) MALLOC(STEP_SIZE*(new->nblocks+1)*COMPRESS_BLOCKSIZE * sizeof(Genomecomp_T));
#else
  new->nblocks = (length+127)/128U;
  new->shift_array = (Genomecomp_T **) MALLOC(STEP_SIZE * sizeof(Genomecomp_T *));
  new->shift_array[0] = (Genomecomp_T *) _mm_malloc(STEP_SIZE*(new->nblocks+1)*COMPRESS_BLOCKSIZE * sizeof(Genomecomp_T),16);
#endif
#ifdef DEBUG14
  new->querylength = length;
#endif

  /* Note that elements of shift_array do not have extra block at beginning */
  new->blocks = new->shift_array[0];
  new->availp[0] = true;	/* Same as new->blocks */
  for (i = 1; i < STEP_SIZE; i++) {
    new->shift_array[i] = &(new->shift_array[i-1][(new->nblocks+1)*COMPRESS_BLOCKSIZE]);
    new->availp[i] = false;
  }

  ptr = 0;

  position = 0U;
  while (position < length) {

#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
    high = low = flags = 0U;
    in_counter = 0;
    while (position < length && in_counter < 32) {
      c = gbuffer[position++];
      high >>= 1;
      low >>= 1;
      flags >>= 1;

      /* Assume that gbuffer is upper case */
      switch /*(uppercaseCode[c])*/ (c) {
      case 'A': /* high |= LEFT_CLEAR; */ /* low |= LEFT_CLEAR; */ /* flags |= LEFT_CLEAR; */ break;
      case 'C': /* high |= LEFT_CLEAR; */    low |= LEFT_SET;      /* flags |= LEFT_CLEAR; */ break;
      case 'G':    high |= LEFT_SET;      /* low |= LEFT_CLEAR; */ /* flags |= LEFT_CLEAR; */ break;
      case 'T':    high |= LEFT_SET;         low |= LEFT_SET;      /* flags |= LEFT_CLEAR; */ break;
      default:  /* high |= LEFT_CLEAR; */ /* low |= LEFT_CLEAR; */    flags |= LEFT_SET;
      }
      in_counter++;
    }
      
    while (in_counter < 32) {
      high >>= 1;
      low >>= 1;
      flags >>= 1;
      in_counter++;
    }

    /* Use old storage method */
    new->blocks[ptr] = high;
    new->blocks[ptr+1] = low;
    new->blocks[ptr+2] = flags;

#else
    for (i = 0; i < 4; i++) {
      /* Word i */
      high = low = flags = 0U;
      in_counter = 0;
      while (position < length && in_counter < 32) {
	c = gbuffer[position++];
	high >>= 1;
	low >>= 1;
	flags >>= 1;

	/* Assume that gbuffer is upper case */
	switch /*(uppercaseCode[c])*/ (c) {
	case 'A': /* high |= LEFT_CLEAR; */ /* low |= LEFT_CLEAR; */ /* flags |= LEFT_CLEAR; */ break;
	case 'C': /* high |= LEFT_CLEAR; */    low |= LEFT_SET;      /* flags |= LEFT_CLEAR; */ break;
	case 'G':    high |= LEFT_SET;      /* low |= LEFT_CLEAR; */ /* flags |= LEFT_CLEAR; */ break;
	case 'T':    high |= LEFT_SET;         low |= LEFT_SET;      /* flags |= LEFT_CLEAR; */ break;
	default:  /* high |= LEFT_CLEAR; */ /* low |= LEFT_CLEAR; */    flags |= LEFT_SET;
	}
	in_counter++;
      }
      
      while (in_counter < 32) {
	high >>= 1;
	low >>= 1;
	flags >>= 1;
	in_counter++;
      }

      new->blocks[ptr + i] = high;
      new->blocks[ptr + i + 4] = low;
      new->blocks[ptr + i + 8] = flags;
    }
#endif

    ptr += COMPRESS_BLOCKSIZE;
  }

#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
  /* Compress_shift will access these values */
  new->blocks[ptr] = 0U;
  new->blocks[ptr+1] = 0U;
  new->blocks[ptr+2] = 0U;
#else
  /* Compress_shift will access these values */
  new->blocks[ptr] = new->blocks[ptr+1] = new->blocks[ptr+2] = new->blocks[ptr+3] = 0U;
  new->blocks[ptr+4] = new->blocks[ptr+5] = new->blocks[ptr+6] = new->blocks[ptr+7] = 0U;
  new->blocks[ptr+8] = new->blocks[ptr+9] = new->blocks[ptr+10] = new->blocks[ptr+11] = 0U;
#endif

  debug0(printf("Compress_new_fwd\n"));
  debug0(Compress_print_blocks(new->blocks,new->nblocks));
  debug0(printf("\n"));

#if 0
  /* Pre-fill nshift of 0, so we don't have to check for that as a special case */
  new->availp[0] = true;
  memcpy(&(new->shift_array[0]),new->blocks,new->nblocks*COMPRESS_BLOCKSIZE*sizeof(Genomecomp_T));
#endif

  return new;
}

T
Compress_new_rev (char *gbuffer, Chrpos_T length) {
  T new = (T) MALLOC(sizeof(*new));
  Genomecomp_T high, low, flags;
  Chrpos_T ptr;
  Chrpos_T position;
  int c, i;
  int in_counter = 0;

#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
  new->nblocks = (length+31)/32U;
  new->shift_array = (Genomecomp_T **) MALLOC(STEP_SIZE * sizeof(Genomecomp_T *));
  new->shift_array[0] = (Genomecomp_T *) MALLOC(STEP_SIZE*(new->nblocks+1)*COMPRESS_BLOCKSIZE * sizeof(Genomecomp_T));
#else
  new->nblocks = (length+127)/128U;
  new->shift_array = (Genomecomp_T **) MALLOC(STEP_SIZE * sizeof(Genomecomp_T *));
  new->shift_array[0] = (Genomecomp_T *) _mm_malloc(STEP_SIZE*(new->nblocks+1)*COMPRESS_BLOCKSIZE * sizeof(Genomecomp_T),16);
#endif
#ifdef DEBUG14
  new->querylength = length;
#endif

  /* Note that elements of shift_array do not have extra block at beginning */
  new->blocks = new->shift_array[0];
  new->availp[0] = true;	/* Same as new->blocks */
  for (i = 1; i < STEP_SIZE; i++) {
    new->shift_array[i] = &(new->shift_array[i-1][(new->nblocks+1)*COMPRESS_BLOCKSIZE]);
    new->availp[i] = false;
  }

  ptr = 0;

  position = length;
  while (position > 0) {

#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
    high = low = flags = 0U;
    in_counter = 0;
    while (position > 0 && in_counter < 32) {
      c = gbuffer[--position];
      high >>= 1;
      low >>= 1;
      flags >>= 1;

      /* Assume that gbuffer is upper case */
      switch /*(uppercaseCode[c])*/ (c) {
      case 'T': /* high |= LEFT_CLEAR; */ /* low |= LEFT_CLEAR; */ /* flags |= LEFT_CLEAR; */ break;
      case 'G': /* high |= LEFT_CLEAR; */    low |= LEFT_SET;      /* flags |= LEFT_CLEAR; */ break;
      case 'C':    high |= LEFT_SET;      /* low |= LEFT_CLEAR; */ /* flags |= LEFT_CLEAR; */ break;
      case 'A':    high |= LEFT_SET;         low |= LEFT_SET;      /* flags |= LEFT_CLEAR; */ break;
      default:  /* high |= LEFT_CLEAR; */ /* low |= LEFT_CLEAR; */    flags |= LEFT_SET;
      }
      in_counter++;
    }

    while (in_counter < 32) {
      high >>= 1;
      low >>= 1;
      flags >>= 1;
      in_counter++;
    }

    new->blocks[ptr] = high;
    new->blocks[ptr+1] = low;
    new->blocks[ptr+2] = flags;

#else
    for (i = 0; i < 4; i++) {
      /* Word i */
      high = low = flags = 0U;
      in_counter = 0;
      while (position > 0 && in_counter < 32) {
	c = gbuffer[--position];
	high >>= 1;
	low >>= 1;
	flags >>= 1;

	/* Assume that gbuffer is upper case */
	switch /*(uppercaseCode[c])*/ (c) {
	case 'T': /* high |= LEFT_CLEAR; */ /* low |= LEFT_CLEAR; */ /* flags |= LEFT_CLEAR; */ break;
	case 'G': /* high |= LEFT_CLEAR; */    low |= LEFT_SET;      /* flags |= LEFT_CLEAR; */ break;
	case 'C':    high |= LEFT_SET;      /* low |= LEFT_CLEAR; */ /* flags |= LEFT_CLEAR; */ break;
	case 'A':    high |= LEFT_SET;         low |= LEFT_SET;      /* flags |= LEFT_CLEAR; */ break;
	default:  /* high |= LEFT_CLEAR; */ /* low |= LEFT_CLEAR; */    flags |= LEFT_SET;
	}
	in_counter++;
      }

      while (in_counter < 32) {
	high >>= 1;
	low >>= 1;
	flags >>= 1;
	in_counter++;
      }

      new->blocks[ptr + i] = high;
      new->blocks[ptr + i + 4] = low;
      new->blocks[ptr + i + 8] = flags;
    }
#endif
    
    ptr += COMPRESS_BLOCKSIZE;
  }

#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
  /* Compress_shift will access these values */
  new->blocks[ptr] = 0U;
  new->blocks[ptr+1] = 0U;
  new->blocks[ptr+2] = 0U;
#else
  /* Compress_shift will access these values */
  new->blocks[ptr] = new->blocks[ptr+1] = new->blocks[ptr+2] = new->blocks[ptr+3] = 0U;
  new->blocks[ptr+4] = new->blocks[ptr+5] = new->blocks[ptr+6] = new->blocks[ptr+7] = 0U;
  new->blocks[ptr+8] = new->blocks[ptr+9] = new->blocks[ptr+10] = new->blocks[ptr+11] = 0U;
#endif

  debug0(printf("Compress_new_rev\n"));
  debug0(Compress_print_blocks(new->blocks,new->nblocks));
  debug0(printf("\n"));

#if 0
  /* Pre-fill nshift of 0, so we don't have to check for that as a special case */
  new->availp[0] = true;
  memcpy(&(new->shift_array[0]),new->blocks,new->nblocks*COMPRESS_BLOCKSIZE*sizeof(Genomecomp_T));
#endif

  return new;
}


#ifdef DEBUG14
static void
print_vector_hex (__m128i x) {
  printf("%08X %08X %08X %08X\n",
	 _mm_extract_epi32(x,0),_mm_extract_epi32(x,1),_mm_extract_epi32(x,2),_mm_extract_epi32(x,3));
  return;
}

static bool
Compress_shift_check (T this, int nshift) {
  Genomecomp_T *shifted;
  int rightshift;
  int ptr;

  /* printf("Entered Compress_shift_check with nshift = %d, nblocks = %d\n",nshift,this->nblocks); */
  shifted = this->shift_array[nshift];

  /* Shift */
  ptr = this->nblocks*12;
  if (nshift == 0) {
    while (ptr >= 0) {
      assert(shifted[ptr] == this->blocks[ptr]);
      ptr--;
    }
    
  } else if (nshift < 32) {
    rightshift = 32 - nshift;
    
    while (ptr > 0) {
      assert(shifted[ptr+11] == ((this->blocks[ptr+11] << nshift) | (this->blocks[ptr+10] >> rightshift)));
      assert(shifted[ptr+7] == ((this->blocks[ptr+7] << nshift) | (this->blocks[ptr+6] >> rightshift)));
      assert(shifted[ptr+3] == ((this->blocks[ptr+3] << nshift) | (this->blocks[ptr+2] >> rightshift)));

      assert(shifted[ptr+10] == ((this->blocks[ptr+10] << nshift) | (this->blocks[ptr+9] >> rightshift)));
      assert(shifted[ptr+6] == ((this->blocks[ptr+6] << nshift) | (this->blocks[ptr+5] >> rightshift)));
      assert(shifted[ptr+2] == ((this->blocks[ptr+2] << nshift) | (this->blocks[ptr+1] >> rightshift)));

      assert(shifted[ptr+9] == ((this->blocks[ptr+9] << nshift) | (this->blocks[ptr+8] >> rightshift)));
      assert(shifted[ptr+5] == ((this->blocks[ptr+5] << nshift) | (this->blocks[ptr+4] >> rightshift)));
      assert(shifted[ptr+1] == ((this->blocks[ptr+1] << nshift) | (this->blocks[ptr] >> rightshift)));

      assert(shifted[ptr+8] == ((this->blocks[ptr+8] << nshift) | (this->blocks[ptr-1] >> rightshift)));
      assert(shifted[ptr+4] == ((this->blocks[ptr+4] << nshift) | (this->blocks[ptr-5] >> rightshift)));
      assert(shifted[ptr] == ((this->blocks[ptr] << nshift) | (this->blocks[ptr-9] >> rightshift)));
      
      ptr -= 12;
    }

    assert(shifted[11] == ((this->blocks[11] << nshift) | (this->blocks[10] >> rightshift)));
    assert(shifted[7] == ((this->blocks[7] << nshift) | (this->blocks[6] >> rightshift)));
    assert(shifted[3] == ((this->blocks[3] << nshift) | (this->blocks[2] >> rightshift)));

    assert(shifted[10] == ((this->blocks[10] << nshift) | (this->blocks[9] >> rightshift)));
    assert(shifted[6] == ((this->blocks[6] << nshift) | (this->blocks[5] >> rightshift)));
    assert(shifted[2] == ((this->blocks[2] << nshift) | (this->blocks[1] >> rightshift)));

    assert(shifted[9] == ((this->blocks[9] << nshift) | (this->blocks[8] >> rightshift)));
    assert(shifted[5] == ((this->blocks[5] << nshift) | (this->blocks[4] >> rightshift)));
    assert(shifted[1] == ((this->blocks[1] << nshift) | (this->blocks[0] >> rightshift)));

    assert(shifted[8] == (this->blocks[8] << nshift));
    assert(shifted[4] == (this->blocks[4] << nshift));
    assert(shifted[0] == (this->blocks[0] << nshift));

  } else if (nshift == 32) {

    while (ptr > 0) {
      assert(shifted[ptr+11] == this->blocks[ptr+10]);
      assert(shifted[ptr+7] == this->blocks[ptr+6]);
      assert(shifted[ptr+3] == this->blocks[ptr+2]);

      assert(shifted[ptr+10] == this->blocks[ptr+9]);
      assert(shifted[ptr+6] == this->blocks[ptr+5]);
      assert(shifted[ptr+2] == this->blocks[ptr+1]);

      assert(shifted[ptr+9] == this->blocks[ptr+8]);
      assert(shifted[ptr+5] == this->blocks[ptr+4]);
      assert(shifted[ptr+1] == this->blocks[ptr]);

      assert(shifted[ptr+8] == this->blocks[ptr-1]);
      assert(shifted[ptr+4] == this->blocks[ptr-5]);
      assert(shifted[ptr] == this->blocks[ptr-9]);

      ptr -= 12;
    }

    assert(shifted[11] == this->blocks[10]);
    assert(shifted[7] == this->blocks[6]);
    assert(shifted[3] == this->blocks[2]);

    assert(shifted[10] == this->blocks[9]);
    assert(shifted[6] == this->blocks[5]);
    assert(shifted[2] == this->blocks[1]);

    assert(shifted[9] == this->blocks[8]);
    assert(shifted[5] == this->blocks[4]);
    assert(shifted[1] == this->blocks[0]);

    assert(shifted[8] == 0U);
    assert(shifted[4] == 0U);
    assert(shifted[0] == 0U);

  } else if (nshift < 64) {
    nshift -= 32;
    rightshift = 32 - nshift;

    while (ptr > 0) {
      assert(shifted[ptr+11] == ((this->blocks[ptr+10] << nshift) | (this->blocks[ptr+9] >> rightshift)));
      assert(shifted[ptr+7] == ((this->blocks[ptr+6] << nshift) | (this->blocks[ptr+5] >> rightshift)));
      assert(shifted[ptr+3] == ((this->blocks[ptr+2] << nshift) | (this->blocks[ptr+1] >> rightshift)));

      assert(shifted[ptr+10] == ((this->blocks[ptr+9] << nshift) | (this->blocks[ptr+8] >> rightshift)));
      assert(shifted[ptr+6] == ((this->blocks[ptr+5] << nshift) | (this->blocks[ptr+4] >> rightshift)));
      assert(shifted[ptr+2] == ((this->blocks[ptr+1] << nshift) | (this->blocks[ptr] >> rightshift)));

      assert(shifted[ptr+9] == ((this->blocks[ptr+8] << nshift) | (this->blocks[ptr-1] >> rightshift)));
      assert(shifted[ptr+5] == ((this->blocks[ptr+4] << nshift) | (this->blocks[ptr-5] >> rightshift)));
      assert(shifted[ptr+1] == ((this->blocks[ptr] << nshift) | (this->blocks[ptr-9] >> rightshift)));

      assert(shifted[ptr+8] == ((this->blocks[ptr-1] << nshift) | (this->blocks[ptr-2] >> rightshift)));
      assert(shifted[ptr+4] == ((this->blocks[ptr-5] << nshift) | (this->blocks[ptr-6] >> rightshift)));
      assert(shifted[ptr] == ((this->blocks[ptr-9] << nshift) | (this->blocks[ptr-10] >> rightshift)));

      ptr -= 12;
    }

    assert(shifted[11] == ((this->blocks[10] << nshift) | (this->blocks[9] >> rightshift)));
    assert(shifted[7] == ((this->blocks[6] << nshift) | (this->blocks[5] >> rightshift)));
    assert(shifted[3] == ((this->blocks[2] << nshift) | (this->blocks[1] >> rightshift)));

    assert(shifted[10] == ((this->blocks[9] << nshift) | (this->blocks[8] >> rightshift)));
    assert(shifted[6] == ((this->blocks[5] << nshift) | (this->blocks[4] >> rightshift)));
    assert(shifted[2] == ((this->blocks[1] << nshift) | (this->blocks[0] >> rightshift)));

    assert(shifted[9] == (this->blocks[8] << nshift));
    assert(shifted[5] == (this->blocks[4] << nshift));
    assert(shifted[1] == (this->blocks[0] << nshift));

    assert(shifted[8] == 0U);
    assert(shifted[4] == 0U);
    assert(shifted[0] == 0U);


  } else if (nshift == 64) {

    while (ptr > 0) {
      assert(shifted[ptr+11] == this->blocks[ptr+9]);
      assert(shifted[ptr+7] == this->blocks[ptr+5]);
      assert(shifted[ptr+3] == this->blocks[ptr+1]);

      assert(shifted[ptr+10] == this->blocks[ptr+8]);
      assert(shifted[ptr+6] == this->blocks[ptr+4]);
      assert(shifted[ptr+2] == this->blocks[ptr]);

      assert(shifted[ptr+9] == this->blocks[ptr-1]);
      assert(shifted[ptr+5] == this->blocks[ptr-5]);
      assert(shifted[ptr+1] == this->blocks[ptr-9]);

      assert(shifted[ptr+8] == this->blocks[ptr-2]);
      assert(shifted[ptr+4] == this->blocks[ptr-6]);
      assert(shifted[ptr] == this->blocks[ptr-10]);

      ptr -= 12;
    }

    assert(shifted[11] == this->blocks[9]);
    assert(shifted[7] == this->blocks[5]);
    assert(shifted[3] == this->blocks[1]);

    assert(shifted[10] == this->blocks[8]);
    assert(shifted[6] == this->blocks[4]);
    assert(shifted[2] == this->blocks[0]);

    assert(shifted[9] == 0U);
    assert(shifted[5] == 0U);
    assert(shifted[1] == 0U);

    assert(shifted[8] == 0U);
    assert(shifted[4] == 0U);
    assert(shifted[0] == 0U);

  } else if (nshift < 96) {
    nshift -= 64;
    rightshift = 32 - nshift;

    while (ptr > 0) {
      assert(shifted[ptr+11] == ((this->blocks[ptr+9] << nshift) | (this->blocks[ptr+8] >> rightshift)));
      assert(shifted[ptr+7] == ((this->blocks[ptr+5] << nshift) | (this->blocks[ptr+4] >> rightshift)));
      assert(shifted[ptr+3] == ((this->blocks[ptr+1] << nshift) | (this->blocks[ptr] >> rightshift)));

      assert(shifted[ptr+10] == ((this->blocks[ptr+8] << nshift) | (this->blocks[ptr-1] >> rightshift)));
      assert(shifted[ptr+6] == ((this->blocks[ptr+4] << nshift) | (this->blocks[ptr-5] >> rightshift)));
      assert(shifted[ptr+2] == ((this->blocks[ptr] << nshift) | (this->blocks[ptr-9] >> rightshift)));

      assert(shifted[ptr+9] == ((this->blocks[ptr-1] << nshift) | (this->blocks[ptr-2] >> rightshift)));
      assert(shifted[ptr+5] == ((this->blocks[ptr-5] << nshift) | (this->blocks[ptr-6] >> rightshift)));
      assert(shifted[ptr+1] == ((this->blocks[ptr-9] << nshift) | (this->blocks[ptr-10] >> rightshift)));

      assert(shifted[ptr+8] == ((this->blocks[ptr-2] << nshift) | (this->blocks[ptr-3] >> rightshift)));
      assert(shifted[ptr+4] == ((this->blocks[ptr-6] << nshift) | (this->blocks[ptr-7] >> rightshift)));
      assert(shifted[ptr] == ((this->blocks[ptr-10] << nshift) | (this->blocks[ptr-11] >> rightshift)));

      ptr -= 12;
    }

    assert(shifted[11] == ((this->blocks[9] << nshift) | (this->blocks[8] >> rightshift)));
    assert(shifted[7] == ((this->blocks[5] << nshift) | (this->blocks[4] >> rightshift)));
    assert(shifted[3] == ((this->blocks[1] << nshift) | (this->blocks[0] >> rightshift)));

    assert(shifted[10] == (this->blocks[8] << nshift));
    assert(shifted[6] == (this->blocks[4] << nshift));
    assert(shifted[2] == (this->blocks[0] << nshift));

    assert(shifted[9] == 0U);
    assert(shifted[5] == 0U);
    assert(shifted[1] == 0U);

    assert(shifted[8] == 0U);
    assert(shifted[4] == 0U);
    assert(shifted[0] == 0U);

  } else if (nshift == 96) {

    while (ptr > 0) {
      assert(shifted[ptr+11] == this->blocks[ptr+8]);
      assert(shifted[ptr+7] == this->blocks[ptr+4]);
      assert(shifted[ptr+3] == this->blocks[ptr]);

      assert(shifted[ptr+10] == this->blocks[ptr-1]);
      assert(shifted[ptr+6] == this->blocks[ptr-5]);
      assert(shifted[ptr+2] == this->blocks[ptr-9]);

      assert(shifted[ptr+9] == this->blocks[ptr-2]);
      assert(shifted[ptr+5] == this->blocks[ptr-6]);
      assert(shifted[ptr+1] == this->blocks[ptr-10]);

      assert(shifted[ptr+8] == this->blocks[ptr-3]);
      assert(shifted[ptr+4] == this->blocks[ptr-7]);
      assert(shifted[ptr] == this->blocks[ptr-11]);

      ptr -= 12;
    }

    assert(shifted[11] == this->blocks[8]);
    assert(shifted[7] == this->blocks[4]);
    assert(shifted[3] == this->blocks[0]);

    assert(shifted[10] == 0U);
    assert(shifted[6] == 0U);
    assert(shifted[2] == 0U);

    assert(shifted[9] == 0U);
    assert(shifted[5] == 0U);
    assert(shifted[1] == 0U);

    assert(shifted[8] == 0U);
    assert(shifted[4] == 0U);
    assert(shifted[0] == 0U);

  } else {
    nshift -= 96;
    rightshift = 32 - nshift;

    while (ptr > 0) {
      assert(shifted[ptr+11] == ((this->blocks[ptr+8] << nshift) | (this->blocks[ptr-1] >> rightshift)));
      assert(shifted[ptr+7] == ((this->blocks[ptr+4] << nshift) | (this->blocks[ptr-5] >> rightshift)));
      assert(shifted[ptr+3] == ((this->blocks[ptr] << nshift) | (this->blocks[ptr-9] >> rightshift)));

      assert(shifted[ptr+10] == ((this->blocks[ptr-1] << nshift) | (this->blocks[ptr-2] >> rightshift)));
      assert(shifted[ptr+6] == ((this->blocks[ptr-5] << nshift) | (this->blocks[ptr-6] >> rightshift)));
      assert(shifted[ptr+2] == ((this->blocks[ptr-9] << nshift) | (this->blocks[ptr-10] >> rightshift)));

      assert(shifted[ptr+9] == ((this->blocks[ptr-2] << nshift) | (this->blocks[ptr-3] >> rightshift)));
      assert(shifted[ptr+5] == ((this->blocks[ptr-6] << nshift) | (this->blocks[ptr-7] >> rightshift)));
      assert(shifted[ptr+1] == ((this->blocks[ptr-10] << nshift) | (this->blocks[ptr-11] >> rightshift)));

      assert(shifted[ptr+8] == ((this->blocks[ptr-3] << nshift) | (this->blocks[ptr-4] >> rightshift)));
      assert(shifted[ptr+4] == ((this->blocks[ptr-7] << nshift) | (this->blocks[ptr-8] >> rightshift)));
      assert(shifted[ptr] == ((this->blocks[ptr-11] << nshift) | (this->blocks[ptr-12] >> rightshift)));

      ptr -= 12;
    }

    assert(shifted[11] == (this->blocks[8] << nshift));
    assert(shifted[7] == (this->blocks[4] << nshift));
    assert(shifted[3] == (this->blocks[0] << nshift));

    assert(shifted[10] == 0U);
    assert(shifted[6] == 0U);
    assert(shifted[2] == 0U);

    assert(shifted[9] == 0U);
    assert(shifted[5] == 0U);
    assert(shifted[1] == 0U);

    assert(shifted[8] == 0U);
    assert(shifted[4] == 0U);
    assert(shifted[0] == 0U);
  }

  return true;
}
#endif


#ifdef DEBUG14
static void
shift_sse2 (T this, int nshift) {
  Genomecomp_T shifted[1000];
  int rightshift;
  int ptr;

  /* Shift */
  ptr = this->nblocks*12;
  if (nshift == 0) {
    while (ptr >= 0) {
      memcpy(&(shifted[ptr]),&(this->blocks[ptr]),12*sizeof(Genomecomp_T));
      ptr -= 12;
    }

  } else if (nshift < 32) {
    rightshift = 32 - nshift;

    while (ptr > 0) {
      shifted[ptr+11] = (this->blocks[ptr+11] << nshift) | (this->blocks[ptr+10] >> rightshift);
      shifted[ptr+7] = (this->blocks[ptr+7] << nshift) | (this->blocks[ptr+6] >> rightshift);
      shifted[ptr+3] = (this->blocks[ptr+3] << nshift) | (this->blocks[ptr+2] >> rightshift);

      shifted[ptr+10] = (this->blocks[ptr+10] << nshift) | (this->blocks[ptr+9] >> rightshift);
      shifted[ptr+6] = (this->blocks[ptr+6] << nshift) | (this->blocks[ptr+5] >> rightshift);
      shifted[ptr+2] = (this->blocks[ptr+2] << nshift) | (this->blocks[ptr+1] >> rightshift);

      shifted[ptr+9] = (this->blocks[ptr+9] << nshift) | (this->blocks[ptr+8] >> rightshift);
      shifted[ptr+5] = (this->blocks[ptr+5] << nshift) | (this->blocks[ptr+4] >> rightshift);
      shifted[ptr+1] = (this->blocks[ptr+1] << nshift) | (this->blocks[ptr] >> rightshift);

      shifted[ptr+8] = (this->blocks[ptr+8] << nshift) | (this->blocks[ptr-1] >> rightshift);
      shifted[ptr+4] = (this->blocks[ptr+4] << nshift) | (this->blocks[ptr-5] >> rightshift);
      shifted[ptr] = (this->blocks[ptr] << nshift) | (this->blocks[ptr-9] >> rightshift);

      ptr -= 12;
    }

    shifted[11] = (this->blocks[11] << nshift) | (this->blocks[10] >> rightshift);
    shifted[7] = (this->blocks[7] << nshift) | (this->blocks[6] >> rightshift);
    shifted[3] = (this->blocks[3] << nshift) | (this->blocks[2] >> rightshift);

    shifted[10] = (this->blocks[10] << nshift) | (this->blocks[9] >> rightshift);
    shifted[6] = (this->blocks[6] << nshift) | (this->blocks[5] >> rightshift);
    shifted[2] = (this->blocks[2] << nshift) | (this->blocks[1] >> rightshift);

    shifted[9] = (this->blocks[9] << nshift) | (this->blocks[8] >> rightshift);
    shifted[5] = (this->blocks[5] << nshift) | (this->blocks[4] >> rightshift);
    shifted[1] = (this->blocks[1] << nshift) | (this->blocks[0] >> rightshift);

    shifted[8] = this->blocks[8] << nshift;
    shifted[4] = this->blocks[4] << nshift;
    shifted[0] = this->blocks[0] << nshift;

  } else if (nshift == 32) {

    while (ptr > 0) {
      shifted[ptr+11] = this->blocks[ptr+10];
      shifted[ptr+7] = this->blocks[ptr+6];
      shifted[ptr+3] = this->blocks[ptr+2];

      shifted[ptr+10] = this->blocks[ptr+9];
      shifted[ptr+6] = this->blocks[ptr+5];
      shifted[ptr+2] = this->blocks[ptr+1];

      shifted[ptr+9] = this->blocks[ptr+8];
      shifted[ptr+5] = this->blocks[ptr+4];
      shifted[ptr+1] = this->blocks[ptr];

      shifted[ptr+8] = this->blocks[ptr-1];
      shifted[ptr+4] = this->blocks[ptr-5];
      shifted[ptr] = this->blocks[ptr-9];

      ptr -= 12;
    }

    shifted[11] = this->blocks[10];
    shifted[7] = this->blocks[6];
    shifted[3] = this->blocks[2];

    shifted[10] = this->blocks[9];
    shifted[6] = this->blocks[5];
    shifted[2] = this->blocks[1];

    shifted[9] = this->blocks[8];
    shifted[5] = this->blocks[4];
    shifted[1] = this->blocks[0];

    shifted[8] = 0U;
    shifted[4] = 0U;
    shifted[0] = 0U;

  } else if (nshift < 64) {
    nshift -= 32;
    rightshift = 32 - nshift;

    while (ptr > 0) {
      shifted[ptr+11] = (this->blocks[ptr+10] << nshift) | (this->blocks[ptr+9] >> rightshift);
      shifted[ptr+7] = (this->blocks[ptr+6] << nshift) | (this->blocks[ptr+5] >> rightshift);
      shifted[ptr+3] = (this->blocks[ptr+2] << nshift) | (this->blocks[ptr+1] >> rightshift);

      shifted[ptr+10] = (this->blocks[ptr+9] << nshift) | (this->blocks[ptr+8] >> rightshift);
      shifted[ptr+6] = (this->blocks[ptr+5] << nshift) | (this->blocks[ptr+4] >> rightshift);
      shifted[ptr+2] = (this->blocks[ptr+1] << nshift) | (this->blocks[ptr] >> rightshift);

      shifted[ptr+9] = (this->blocks[ptr+8] << nshift) | (this->blocks[ptr-1] >> rightshift);
      shifted[ptr+5] = (this->blocks[ptr+4] << nshift) | (this->blocks[ptr-5] >> rightshift);
      shifted[ptr+1] = (this->blocks[ptr] << nshift) | (this->blocks[ptr-9] >> rightshift);

      shifted[ptr+8] = (this->blocks[ptr-1] << nshift) | (this->blocks[ptr-2] >> rightshift);
      shifted[ptr+4] = (this->blocks[ptr-5] << nshift) | (this->blocks[ptr-6] >> rightshift);
      shifted[ptr] = (this->blocks[ptr-9] << nshift) | (this->blocks[ptr-10] >> rightshift);

      ptr -= 12;
    }

    shifted[11] = (this->blocks[10] << nshift) | (this->blocks[9] >> rightshift);
    shifted[7] = (this->blocks[6] << nshift) | (this->blocks[5] >> rightshift);
    shifted[3] = (this->blocks[2] << nshift) | (this->blocks[1] >> rightshift);

    shifted[10] = (this->blocks[9] << nshift) | (this->blocks[8] >> rightshift);
    shifted[6] = (this->blocks[5] << nshift) | (this->blocks[4] >> rightshift);
    shifted[2] = (this->blocks[1] << nshift) | (this->blocks[0] >> rightshift);

    shifted[9] = this->blocks[8] << nshift;
    shifted[5] = this->blocks[4] << nshift;
    shifted[1] = this->blocks[0] << nshift;

    shifted[8] = 0U;
    shifted[4] = 0U;
    shifted[0] = 0U;


  } else if (nshift == 64) {

    while (ptr > 0) {
      shifted[ptr+11] = this->blocks[ptr+9];
      shifted[ptr+7] = this->blocks[ptr+5];
      shifted[ptr+3] = this->blocks[ptr+1];

      shifted[ptr+10] = this->blocks[ptr+8];
      shifted[ptr+6] = this->blocks[ptr+4];
      shifted[ptr+2] = this->blocks[ptr];

      shifted[ptr+9] = this->blocks[ptr-1];
      shifted[ptr+5] = this->blocks[ptr-5];
      shifted[ptr+1] = this->blocks[ptr-9];

      shifted[ptr+8] = this->blocks[ptr-2];
      shifted[ptr+4] = this->blocks[ptr-6];
      shifted[ptr] = this->blocks[ptr-10];

      ptr -= 12;
    }

    shifted[11] = this->blocks[9];
    shifted[7] = this->blocks[5];
    shifted[3] = this->blocks[1];

    shifted[10] = this->blocks[8];
    shifted[6] = this->blocks[4];
    shifted[2] = this->blocks[0];

    shifted[9] = 0U;
    shifted[5] = 0U;
    shifted[1] = 0U;

    shifted[8] = 0U;
    shifted[4] = 0U;
    shifted[0] = 0U;

  } else if (nshift < 96) {
    nshift -= 64;
    rightshift = 32 - nshift;

    while (ptr > 0) {
      shifted[ptr+11] = (this->blocks[ptr+9] << nshift) | (this->blocks[ptr+8] >> rightshift);
      shifted[ptr+7] = (this->blocks[ptr+5] << nshift) | (this->blocks[ptr+4] >> rightshift);
      shifted[ptr+3] = (this->blocks[ptr+1] << nshift) | (this->blocks[ptr] >> rightshift);

      shifted[ptr+10] = (this->blocks[ptr+8] << nshift) | (this->blocks[ptr-1] >> rightshift);
      shifted[ptr+6] = (this->blocks[ptr+4] << nshift) | (this->blocks[ptr-5] >> rightshift);
      shifted[ptr+2] = (this->blocks[ptr] << nshift) | (this->blocks[ptr-9] >> rightshift);

      shifted[ptr+9] = (this->blocks[ptr-1] << nshift) | (this->blocks[ptr-2] >> rightshift);
      shifted[ptr+5] = (this->blocks[ptr-5] << nshift) | (this->blocks[ptr-6] >> rightshift);
      shifted[ptr+1] = (this->blocks[ptr-9] << nshift) | (this->blocks[ptr-10] >> rightshift);

      shifted[ptr+8] = (this->blocks[ptr-2] << nshift) | (this->blocks[ptr-3] >> rightshift);
      shifted[ptr+4] = (this->blocks[ptr-6] << nshift) | (this->blocks[ptr-7] >> rightshift);
      shifted[ptr] = (this->blocks[ptr-10] << nshift) | (this->blocks[ptr-11] >> rightshift);

      ptr -= 12;
    }

    shifted[11] = (this->blocks[9] << nshift) | (this->blocks[8] >> rightshift);
    shifted[7] = (this->blocks[5] << nshift) | (this->blocks[4] >> rightshift);
    shifted[3] = (this->blocks[1] << nshift) | (this->blocks[0] >> rightshift);

    shifted[10] = this->blocks[8] << nshift;
    shifted[6] = this->blocks[4] << nshift;
    shifted[2] = this->blocks[0] << nshift;

    shifted[9] = 0U;
    shifted[5] = 0U;
    shifted[1] = 0U;

    shifted[8] = 0U;
    shifted[4] = 0U;
    shifted[0] = 0U;

  } else if (nshift == 96) {

    while (ptr > 0) {
      shifted[ptr+11] = this->blocks[ptr+8];
      shifted[ptr+7] = this->blocks[ptr+4];
      shifted[ptr+3] = this->blocks[ptr];

      shifted[ptr+10] = this->blocks[ptr-1];
      shifted[ptr+6] = this->blocks[ptr-5];
      shifted[ptr+2] = this->blocks[ptr-9];

      shifted[ptr+9] = this->blocks[ptr-2];
      shifted[ptr+5] = this->blocks[ptr-6];
      shifted[ptr+1] = this->blocks[ptr-10];

      shifted[ptr+8] = this->blocks[ptr-3];
      shifted[ptr+4] = this->blocks[ptr-7];
      shifted[ptr] = this->blocks[ptr-11];

      ptr -= 12;
    }

    shifted[11] = this->blocks[8];
    shifted[7] = this->blocks[4];
    shifted[3] = this->blocks[0];

    shifted[10] = 0U;
    shifted[6] = 0U;
    shifted[2] = 0U;

    shifted[9] = 0U;
    shifted[5] = 0U;
    shifted[1] = 0U;

    shifted[8] = 0U;
    shifted[4] = 0U;
    shifted[0] = 0U;

  } else {
    nshift -= 96;
    rightshift = 32 - nshift;

    while (ptr > 0) {
      shifted[ptr+11] = (this->blocks[ptr+8] << nshift) | (this->blocks[ptr-1] >> rightshift);
      shifted[ptr+7] = (this->blocks[ptr+4] << nshift) | (this->blocks[ptr-5] >> rightshift);
      shifted[ptr+3] = (this->blocks[ptr] << nshift) | (this->blocks[ptr-9] >> rightshift);

      shifted[ptr+10] = (this->blocks[ptr-1] << nshift) | (this->blocks[ptr-2] >> rightshift);
      shifted[ptr+6] = (this->blocks[ptr-5] << nshift) | (this->blocks[ptr-6] >> rightshift);
      shifted[ptr+2] = (this->blocks[ptr-9] << nshift) | (this->blocks[ptr-10] >> rightshift);

      shifted[ptr+9] = (this->blocks[ptr-2] << nshift) | (this->blocks[ptr-3] >> rightshift);
      shifted[ptr+5] = (this->blocks[ptr-6] << nshift) | (this->blocks[ptr-7] >> rightshift);
      shifted[ptr+1] = (this->blocks[ptr-10] << nshift) | (this->blocks[ptr-11] >> rightshift);

      shifted[ptr+8] = (this->blocks[ptr-3] << nshift) | (this->blocks[ptr-4] >> rightshift);
      shifted[ptr+4] = (this->blocks[ptr-7] << nshift) | (this->blocks[ptr-8] >> rightshift);
      shifted[ptr] = (this->blocks[ptr-11] << nshift) | (this->blocks[ptr-12] >> rightshift);

      ptr -= 12;
    }

    shifted[11] = this->blocks[8] << nshift;
    shifted[7] = this->blocks[4] << nshift;
    shifted[3] = this->blocks[0] << nshift;

    shifted[10] = 0U;
    shifted[6] = 0U;
    shifted[2] = 0U;

    shifted[9] = 0U;
    shifted[5] = 0U;
    shifted[1] = 0U;

    shifted[8] = 0U;
    shifted[4] = 0U;
    shifted[0] = 0U;
  }

  Compress_print_blocks(shifted,nshift,0,this->querylength);

  return;
}
#endif



#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
Genomecomp_T *
Compress_shift (T this, int nshift) {
  Genomecomp_T *shifted;
  int rightshift;
  int ptr;
#ifdef DEBUG9
  Genomecomp_T high, low, flags;
#endif

  if (this->availp[nshift] == true) {
    return this->shift_array[nshift];

  } else {
    shifted = this->shift_array[nshift];

    /* Shift */
    ptr = this->nblocks*COMPRESS_BLOCKSIZE;
    rightshift = 32 - nshift;

    while (ptr > 0) {
      shifted[ptr+2] = (this->blocks[ptr+2] << nshift) | (this->blocks[ptr-1] >> rightshift);
      shifted[ptr+1] = (this->blocks[ptr+1] << nshift) | (this->blocks[ptr-2] >> rightshift);
      shifted[ptr] = (this->blocks[ptr] << nshift) | (this->blocks[ptr-3] >> rightshift);
      ptr -= COMPRESS_BLOCKSIZE;
    }

    shifted[2] = this->blocks[2] << nshift;
    shifted[1] = this->blocks[1] << nshift;
    shifted[0] = this->blocks[0] << nshift;
  }

  this->availp[nshift] = true;
  
  debug1(Compress_print_blocks(shifted,this->nblocks+1));
  return shifted;
}

#elif defined(HAVE_SSSE3) && !defined(DEFECTIVE_SSE2_COMPILER)
Genomecomp_T *
Compress_shift (T this, int nshift) {
  Genomecomp_T *shifted;
  int leftshift_words, leftshift_bits, rightshift_bits;
  int ptr, startptr, blocki;
  __m128i out, current, prev, leftpart, rightpart, zeroes;
#ifdef DEBUG9
  Genomecomp_T high, low, flags;
#endif

  if (this->availp[nshift] == true) {
    return this->shift_array[nshift];

  } else {
    /* printf("Original blocks\n"); */
    /* Compress_print_blocks(this->shift_array[0],0,0,this->querylength); */

    shifted = this->shift_array[nshift];
    leftshift_words = nshift / 32;
    leftshift_bits = nshift % 32;
    rightshift_bits = 32 - leftshift_bits;
    zeroes = _mm_set1_epi32(0U);
    /* printf("Entered Compress_shift with nshift %d => leftshift_words %d, leftshift_bits %d, rightshift_bits %d\n",
       nshift,leftshift_words,leftshift_bits,rightshift_bits); */

    /* Take care of high (startptr 0), then low (startptr 4), then flags (startptr 8) */
    for (startptr = 0; startptr < 12; startptr += 4) {
      /* printf("startptr = %d\n",startptr); */
      ptr = startptr;
      prev = zeroes;

      for (blocki = 0; blocki <= this->nblocks; blocki++) {
	current = _mm_load_si128((__m128i *) &(this->blocks[ptr]));

	/* printf("prev:      "); */
	/* print_vector_hex(prev); */
	/* printf("current:   "); */
	/* print_vector_hex(current); */

	switch (leftshift_words) {
	case 0:
	  leftpart = _mm_alignr_epi8(current,prev,12);
	  rightpart = current;	/* _mm_alignr_epi(current,prev,16) */
	  break;
	case 1:
	  leftpart = _mm_alignr_epi8(current,prev,8);
	  rightpart = _mm_alignr_epi8(current,prev,12);
	  break;
	case 2:
	  leftpart = _mm_alignr_epi8(current,prev,4);
	  rightpart = _mm_alignr_epi8(current,prev,8);
	  break;
	case 3:
	  leftpart = prev;	/* _mm_alignr_epi(current,prev,0) */
	  rightpart = _mm_alignr_epi8(current,prev,4);
	  break;
	}

	/* printf("leftpart:  "); */
	/* print_vector_hex(leftpart); */
	/* printf("rightpart: "); */
	/* print_vector_hex(rightpart); */
	/* printf("\n"); */

	/* Note: This line needs a non-defective SSE2 compiler */
	out = _mm_or_si128(_mm_srli_epi32(leftpart,rightshift_bits),_mm_slli_epi32(rightpart,leftshift_bits));
	_mm_store_si128((__m128i *) &(shifted[ptr]),out);
	prev = current;
	ptr += 12;
      }
    }
  }

  this->availp[nshift] = true;

  /* printf("Correct answer\n"); */
  /* shift_sse2(this,nshift); */

  /* printf("New answer\n"); */
  /* Compress_print_blocks(shifted,nshift,0,this->querylength); */
  debug14(Compress_shift_check(this,nshift));
  debug1(Compress_print_blocks(shifted,this->nblocks+1));
  return shifted;
}

#else
/* HAVE_SSE2 but not SSSE3 */
Genomecomp_T *
Compress_shift (T this, int nshift) {
  Genomecomp_T *shifted;
  int rightshift;
  int ptr;

  if (this->availp[nshift] == true) {
    return this->shift_array[nshift];

  } else {
    this->availp[nshift] = true; /* Need to set here, before nshift changes below */
    shifted = this->shift_array[nshift];

    /* Shift */
    ptr = this->nblocks*12;
    if (nshift == 0) {
      while (ptr >= 0) {
	memcpy(&(shifted[ptr]),&(this->blocks[ptr]),12*sizeof(Genomecomp_T));
	ptr -= 12;
      }

    } else if (nshift < 32) {
      rightshift = 32 - nshift;

      while (ptr > 0) {
	shifted[ptr+11] = (this->blocks[ptr+11] << nshift) | (this->blocks[ptr+10] >> rightshift);
	shifted[ptr+7] = (this->blocks[ptr+7] << nshift) | (this->blocks[ptr+6] >> rightshift);
	shifted[ptr+3] = (this->blocks[ptr+3] << nshift) | (this->blocks[ptr+2] >> rightshift);

	shifted[ptr+10] = (this->blocks[ptr+10] << nshift) | (this->blocks[ptr+9] >> rightshift);
	shifted[ptr+6] = (this->blocks[ptr+6] << nshift) | (this->blocks[ptr+5] >> rightshift);
	shifted[ptr+2] = (this->blocks[ptr+2] << nshift) | (this->blocks[ptr+1] >> rightshift);

	shifted[ptr+9] = (this->blocks[ptr+9] << nshift) | (this->blocks[ptr+8] >> rightshift);
	shifted[ptr+5] = (this->blocks[ptr+5] << nshift) | (this->blocks[ptr+4] >> rightshift);
	shifted[ptr+1] = (this->blocks[ptr+1] << nshift) | (this->blocks[ptr] >> rightshift);

	shifted[ptr+8] = (this->blocks[ptr+8] << nshift) | (this->blocks[ptr-1] >> rightshift);
	shifted[ptr+4] = (this->blocks[ptr+4] << nshift) | (this->blocks[ptr-5] >> rightshift);
	shifted[ptr] = (this->blocks[ptr] << nshift) | (this->blocks[ptr-9] >> rightshift);

	ptr -= 12;
      }

      shifted[11] = (this->blocks[11] << nshift) | (this->blocks[10] >> rightshift);
      shifted[7] = (this->blocks[7] << nshift) | (this->blocks[6] >> rightshift);
      shifted[3] = (this->blocks[3] << nshift) | (this->blocks[2] >> rightshift);

      shifted[10] = (this->blocks[10] << nshift) | (this->blocks[9] >> rightshift);
      shifted[6] = (this->blocks[6] << nshift) | (this->blocks[5] >> rightshift);
      shifted[2] = (this->blocks[2] << nshift) | (this->blocks[1] >> rightshift);

      shifted[9] = (this->blocks[9] << nshift) | (this->blocks[8] >> rightshift);
      shifted[5] = (this->blocks[5] << nshift) | (this->blocks[4] >> rightshift);
      shifted[1] = (this->blocks[1] << nshift) | (this->blocks[0] >> rightshift);

      shifted[8] = this->blocks[8] << nshift;
      shifted[4] = this->blocks[4] << nshift;
      shifted[0] = this->blocks[0] << nshift;

    } else if (nshift == 32) {

      while (ptr > 0) {
	shifted[ptr+11] = this->blocks[ptr+10];
	shifted[ptr+7] = this->blocks[ptr+6];
	shifted[ptr+3] = this->blocks[ptr+2];

	shifted[ptr+10] = this->blocks[ptr+9];
	shifted[ptr+6] = this->blocks[ptr+5];
	shifted[ptr+2] = this->blocks[ptr+1];

	shifted[ptr+9] = this->blocks[ptr+8];
	shifted[ptr+5] = this->blocks[ptr+4];
	shifted[ptr+1] = this->blocks[ptr];

	shifted[ptr+8] = this->blocks[ptr-1];
	shifted[ptr+4] = this->blocks[ptr-5];
	shifted[ptr] = this->blocks[ptr-9];

	ptr -= 12;
      }

      shifted[11] = this->blocks[10];
      shifted[7] = this->blocks[6];
      shifted[3] = this->blocks[2];

      shifted[10] = this->blocks[9];
      shifted[6] = this->blocks[5];
      shifted[2] = this->blocks[1];

      shifted[9] = this->blocks[8];
      shifted[5] = this->blocks[4];
      shifted[1] = this->blocks[0];

      shifted[8] = 0U;
      shifted[4] = 0U;
      shifted[0] = 0U;

    } else if (nshift < 64) {
      nshift -= 32;
      rightshift = 32 - nshift;

      while (ptr > 0) {
	shifted[ptr+11] = (this->blocks[ptr+10] << nshift) | (this->blocks[ptr+9] >> rightshift);
	shifted[ptr+7] = (this->blocks[ptr+6] << nshift) | (this->blocks[ptr+5] >> rightshift);
	shifted[ptr+3] = (this->blocks[ptr+2] << nshift) | (this->blocks[ptr+1] >> rightshift);

	shifted[ptr+10] = (this->blocks[ptr+9] << nshift) | (this->blocks[ptr+8] >> rightshift);
	shifted[ptr+6] = (this->blocks[ptr+5] << nshift) | (this->blocks[ptr+4] >> rightshift);
	shifted[ptr+2] = (this->blocks[ptr+1] << nshift) | (this->blocks[ptr] >> rightshift);

	shifted[ptr+9] = (this->blocks[ptr+8] << nshift) | (this->blocks[ptr-1] >> rightshift);
	shifted[ptr+5] = (this->blocks[ptr+4] << nshift) | (this->blocks[ptr-5] >> rightshift);
	shifted[ptr+1] = (this->blocks[ptr] << nshift) | (this->blocks[ptr-9] >> rightshift);

	shifted[ptr+8] = (this->blocks[ptr-1] << nshift) | (this->blocks[ptr-2] >> rightshift);
	shifted[ptr+4] = (this->blocks[ptr-5] << nshift) | (this->blocks[ptr-6] >> rightshift);
	shifted[ptr] = (this->blocks[ptr-9] << nshift) | (this->blocks[ptr-10] >> rightshift);

	ptr -= 12;
      }

      shifted[11] = (this->blocks[10] << nshift) | (this->blocks[9] >> rightshift);
      shifted[7] = (this->blocks[6] << nshift) | (this->blocks[5] >> rightshift);
      shifted[3] = (this->blocks[2] << nshift) | (this->blocks[1] >> rightshift);

      shifted[10] = (this->blocks[9] << nshift) | (this->blocks[8] >> rightshift);
      shifted[6] = (this->blocks[5] << nshift) | (this->blocks[4] >> rightshift);
      shifted[2] = (this->blocks[1] << nshift) | (this->blocks[0] >> rightshift);

      shifted[9] = this->blocks[8] << nshift;
      shifted[5] = this->blocks[4] << nshift;
      shifted[1] = this->blocks[0] << nshift;

      shifted[8] = 0U;
      shifted[4] = 0U;
      shifted[0] = 0U;


    } else if (nshift == 64) {

      while (ptr > 0) {
	shifted[ptr+11] = this->blocks[ptr+9];
	shifted[ptr+7] = this->blocks[ptr+5];
	shifted[ptr+3] = this->blocks[ptr+1];

	shifted[ptr+10] = this->blocks[ptr+8];
	shifted[ptr+6] = this->blocks[ptr+4];
	shifted[ptr+2] = this->blocks[ptr];

	shifted[ptr+9] = this->blocks[ptr-1];
	shifted[ptr+5] = this->blocks[ptr-5];
	shifted[ptr+1] = this->blocks[ptr-9];

	shifted[ptr+8] = this->blocks[ptr-2];
	shifted[ptr+4] = this->blocks[ptr-6];
	shifted[ptr] = this->blocks[ptr-10];

	ptr -= 12;
      }

      shifted[11] = this->blocks[9];
      shifted[7] = this->blocks[5];
      shifted[3] = this->blocks[1];

      shifted[10] = this->blocks[8];
      shifted[6] = this->blocks[4];
      shifted[2] = this->blocks[0];

      shifted[9] = 0U;
      shifted[5] = 0U;
      shifted[1] = 0U;

      shifted[8] = 0U;
      shifted[4] = 0U;
      shifted[0] = 0U;

    } else if (nshift < 96) {
      nshift -= 64;
      rightshift = 32 - nshift;

      while (ptr > 0) {
	shifted[ptr+11] = (this->blocks[ptr+9] << nshift) | (this->blocks[ptr+8] >> rightshift);
	shifted[ptr+7] = (this->blocks[ptr+5] << nshift) | (this->blocks[ptr+4] >> rightshift);
	shifted[ptr+3] = (this->blocks[ptr+1] << nshift) | (this->blocks[ptr] >> rightshift);

	shifted[ptr+10] = (this->blocks[ptr+8] << nshift) | (this->blocks[ptr-1] >> rightshift);
	shifted[ptr+6] = (this->blocks[ptr+4] << nshift) | (this->blocks[ptr-5] >> rightshift);
	shifted[ptr+2] = (this->blocks[ptr] << nshift) | (this->blocks[ptr-9] >> rightshift);

	shifted[ptr+9] = (this->blocks[ptr-1] << nshift) | (this->blocks[ptr-2] >> rightshift);
	shifted[ptr+5] = (this->blocks[ptr-5] << nshift) | (this->blocks[ptr-6] >> rightshift);
	shifted[ptr+1] = (this->blocks[ptr-9] << nshift) | (this->blocks[ptr-10] >> rightshift);

	shifted[ptr+8] = (this->blocks[ptr-2] << nshift) | (this->blocks[ptr-3] >> rightshift);
	shifted[ptr+4] = (this->blocks[ptr-6] << nshift) | (this->blocks[ptr-7] >> rightshift);
	shifted[ptr] = (this->blocks[ptr-10] << nshift) | (this->blocks[ptr-11] >> rightshift);

	ptr -= 12;
      }

      shifted[11] = (this->blocks[9] << nshift) | (this->blocks[8] >> rightshift);
      shifted[7] = (this->blocks[5] << nshift) | (this->blocks[4] >> rightshift);
      shifted[3] = (this->blocks[1] << nshift) | (this->blocks[0] >> rightshift);

      shifted[10] = this->blocks[8] << nshift;
      shifted[6] = this->blocks[4] << nshift;
      shifted[2] = this->blocks[0] << nshift;

      shifted[9] = 0U;
      shifted[5] = 0U;
      shifted[1] = 0U;

      shifted[8] = 0U;
      shifted[4] = 0U;
      shifted[0] = 0U;

    } else if (nshift == 96) {

      while (ptr > 0) {
	shifted[ptr+11] = this->blocks[ptr+8];
	shifted[ptr+7] = this->blocks[ptr+4];
	shifted[ptr+3] = this->blocks[ptr];

	shifted[ptr+10] = this->blocks[ptr-1];
	shifted[ptr+6] = this->blocks[ptr-5];
	shifted[ptr+2] = this->blocks[ptr-9];

	shifted[ptr+9] = this->blocks[ptr-2];
	shifted[ptr+5] = this->blocks[ptr-6];
	shifted[ptr+1] = this->blocks[ptr-10];

	shifted[ptr+8] = this->blocks[ptr-3];
	shifted[ptr+4] = this->blocks[ptr-7];
	shifted[ptr] = this->blocks[ptr-11];

	ptr -= 12;
      }

      shifted[11] = this->blocks[8];
      shifted[7] = this->blocks[4];
      shifted[3] = this->blocks[0];

      shifted[10] = 0U;
      shifted[6] = 0U;
      shifted[2] = 0U;

      shifted[9] = 0U;
      shifted[5] = 0U;
      shifted[1] = 0U;

      shifted[8] = 0U;
      shifted[4] = 0U;
      shifted[0] = 0U;

    } else {
      nshift -= 96;
      rightshift = 32 - nshift;

      while (ptr > 0) {
	shifted[ptr+11] = (this->blocks[ptr+8] << nshift) | (this->blocks[ptr-1] >> rightshift);
	shifted[ptr+7] = (this->blocks[ptr+4] << nshift) | (this->blocks[ptr-5] >> rightshift);
	shifted[ptr+3] = (this->blocks[ptr] << nshift) | (this->blocks[ptr-9] >> rightshift);

	shifted[ptr+10] = (this->blocks[ptr-1] << nshift) | (this->blocks[ptr-2] >> rightshift);
	shifted[ptr+6] = (this->blocks[ptr-5] << nshift) | (this->blocks[ptr-6] >> rightshift);
	shifted[ptr+2] = (this->blocks[ptr-9] << nshift) | (this->blocks[ptr-10] >> rightshift);

	shifted[ptr+9] = (this->blocks[ptr-2] << nshift) | (this->blocks[ptr-3] >> rightshift);
	shifted[ptr+5] = (this->blocks[ptr-6] << nshift) | (this->blocks[ptr-7] >> rightshift);
	shifted[ptr+1] = (this->blocks[ptr-10] << nshift) | (this->blocks[ptr-11] >> rightshift);

	shifted[ptr+8] = (this->blocks[ptr-3] << nshift) | (this->blocks[ptr-4] >> rightshift);
	shifted[ptr+4] = (this->blocks[ptr-7] << nshift) | (this->blocks[ptr-8] >> rightshift);
	shifted[ptr] = (this->blocks[ptr-11] << nshift) | (this->blocks[ptr-12] >> rightshift);

	ptr -= 12;
      }

      shifted[11] = this->blocks[8] << nshift;
      shifted[7] = this->blocks[4] << nshift;
      shifted[3] = this->blocks[0] << nshift;

      shifted[10] = 0U;
      shifted[6] = 0U;
      shifted[2] = 0U;

      shifted[9] = 0U;
      shifted[5] = 0U;
      shifted[1] = 0U;

      shifted[8] = 0U;
      shifted[4] = 0U;
      shifted[0] = 0U;
    }

    debug1(Compress_print_blocks(shifted,this->nblocks+1));
    return shifted;
  }
}
#endif


#if 0
/* Needed for Genome_count_mismatches_fragment_left/right */
Genomecomp_T *
Compress32_shift (T this, int nshift) {
  Genomecomp_T *shifted;
  int rightshift;
  int ptr;
#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
#else
  __m128i out, current, next;
#endif
#ifdef DEBUG9
  Genomecomp_T high, low, flags;
#endif

  if (this->availp[nshift] == true) {
    return this->shift_array[nshift];

  } else {
    shifted = this->shift_array[nshift];

    /* Shift */
    ptr = this->nblocks*COMPRESS_BLOCKSIZE;
    if (nshift == 0) {
      while (ptr >= 0) {
	memcpy(&(shifted[ptr]),&(this->blocks[ptr]),COMPRESS_BLOCKSIZE*sizeof(Genomecomp_T));
	ptr -= COMPRESS_BLOCKSIZE;
      }

    } else {
      rightshift = 32 - nshift;

#ifdef DEFECTIVE_SSE2_COMPILER
      while (ptr > 0) {
	shifted[ptr+2] = (this->blocks[ptr+2] << nshift) | (this->blocks[ptr-2] >> rightshift);
	shifted[ptr+1] = (this->blocks[ptr+1] << nshift) | (this->blocks[ptr-3] >> rightshift);
	shifted[ptr] = (this->blocks[ptr] << nshift) | (this->blocks[ptr-4] >> rightshift);
	ptr -= COMPRESS_BLOCKSIZE;
      }

      shifted[2] = this->blocks[2] << nshift;
      shifted[1] = this->blocks[1] << nshift;
      shifted[0] = this->blocks[0] << nshift;

#elif defined(HAVE_SSE2) && !defined(WORDS_BIGENDIAN)
      next = _mm_load_si128((__m128i *) &(this->blocks[ptr]));
      while (ptr > 0) {
	current = next;
	next = _mm_load_si128((__m128i *) &(this->blocks[ptr-4]));
	out = _mm_or_si128(_mm_slli_epi32(current,nshift),_mm_srli_epi32(next,rightshift));
	_mm_store_si128((__m128i *) &(shifted[ptr]),out);

#ifdef DEBUG9
	flags = (this->blocks[ptr+2] << nshift) | (this->blocks[ptr-2] >> rightshift);
	low = (this->blocks[ptr+1] << nshift) | (this->blocks[ptr-3] >> rightshift);
	high = (this->blocks[ptr] << nshift) | (this->blocks[ptr-4] >> rightshift);
	/* printf("Comparing %08X with %08X, %08X with %08X, %08X with %08X\n",
	   shifted[ptr],high,shifted[ptr+1],low,shifted[ptr+2],flags); */
	if (shifted[ptr+2] != flags) abort();
	if (shifted[ptr+1] != low) abort();
	if (shifted[ptr] != high) abort();
#endif

	ptr -= COMPRESS_BLOCKSIZE;
      }

      out = _mm_slli_epi32(next,nshift);
      _mm_store_si128((__m128i *) &(shifted[0]),out);

#ifdef DEBUG9
      flags = this->blocks[2] << nshift;
      low = this->blocks[1] << nshift;
      high = this->blocks[0] << nshift;
      /* printf("Comparing %08X with %08X, %08X with %08X, %08X with %08X\n",
	 shifted[0],high,shifted[1],low,shifted[2],flags); */
      if (shifted[2] != flags) abort();
      if (shifted[1] != low) abort();
      if (shifted[0] != high) abort();
#endif

#else
      while (ptr > 0) {
	shifted[ptr+2] = (this->blocks[ptr+2] << nshift) | (this->blocks[ptr-2] >> rightshift);
	shifted[ptr+1] = (this->blocks[ptr+1] << nshift) | (this->blocks[ptr-3] >> rightshift);
	shifted[ptr] = (this->blocks[ptr] << nshift) | (this->blocks[ptr-4] >> rightshift);
	ptr -= COMPRESS_BLOCKSIZE;
      }

      shifted[2] = this->blocks[2] << nshift;
      shifted[1] = this->blocks[1] << nshift;
      shifted[0] = this->blocks[0] << nshift;
#endif

    }

    this->availp[nshift] = true;

    debug1(Compress_print_blocks(shifted,this->nblocks+1));
    return shifted;
  }
}
#endif


/* Fragment from pos5 to pos3 is in low 16 bits of each UINT4, with
   pos3 sitting at bit 15. */
void
Compress_get_16mer_left (UINT4 *high, UINT4 *low, UINT4 *flags, T this, int pos3) {
  int leftshift, rightshift;
  int columni, blocki;
  Genomecomp_T *ptr, curr_high, curr_low, curr_flags, prev_high, prev_low, prev_flags;

#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
  /* query is stored as 3 x 32-bit words */
  blocki = pos3/32U*3;

  ptr = &(this->blocks[blocki]);
  curr_high = ptr[0];
  curr_low = ptr[1];
  curr_flags = ptr[2];

  if (blocki == 0) {
    prev_high = prev_low = prev_flags = 0U;
  } else {
    ptr -= 3;
    prev_high = ptr[0];
    prev_low = ptr[1];
    prev_flags = ptr[2];
  }

#else
  /* query is stored as 3 x 128-bit words */
  columni = (pos3 % 128) / 32;
  blocki = pos3/128U*12 + columni;

  ptr = &(this->blocks[blocki]);
  curr_high = ptr[0];
  curr_low = ptr[4];
  curr_flags = ptr[8];

  ptr -= 1;
  if (columni != 0) {
    prev_high = ptr[0];
    prev_low = ptr[4];
    prev_flags = ptr[8];
  } else if (blocki == 0) {
    prev_high = prev_low = prev_flags = 0U;
  } else {
    ptr -= 8;
    prev_high = ptr[0];
    prev_low = ptr[4];
    prev_flags = ptr[8];
  }
#endif


  debug2(printf("high:  %08X %08X\n",prev_high,curr_high));
  debug2(printf("low:   %08X %08X\n",prev_low,curr_low));
  debug2(printf("flags: %08X %08X\n",prev_flags,curr_flags));

  rightshift = pos3 % 32;
  leftshift = 32 - rightshift;
  *high = (curr_high << leftshift) | (prev_high >> rightshift);
  *low = (curr_low << leftshift) | (prev_low >> rightshift);
  *flags = (curr_flags << leftshift) | (prev_flags >> rightshift);

  *high >>= 16;
  *low >>= 16;
  *flags >>= 16;

  debug2(printf("high:  %08X\n",*high));
  debug2(printf("low:   %08X\n",*low));
  debug2(printf("flags: %08X\n",*flags));

  return;
}



/* Fragment from pos5 to pos3 is in low 16 bits of each UINT4, with
   pos5 sitting at bit 0. */
void
Compress_get_16mer_right (UINT4 *high, UINT4 *low, UINT4 *flags, T this, int pos5) {
  int leftshift, rightshift;
  int columni, blocki;
  Genomecomp_T *ptr, curr_high, curr_low, curr_flags, next_high, next_low, next_flags;

#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
  /* query is stored as 3 x 32-bit words */
  blocki = pos5/32U*3;

  ptr = &(this->blocks[blocki]);
  curr_high = ptr[0];
  curr_low = ptr[1];
  curr_flags = ptr[2];

  ptr += 3;
  next_high = ptr[0];
  next_low = ptr[1];
  next_flags = ptr[2];

#else
  /* query is stored as 3 x 128-bit words */
  columni = (pos5 % 128) / 32;
  blocki = pos5/128U*12 + columni;

  ptr = &(this->blocks[blocki]);
  curr_high = ptr[0];
  curr_low = ptr[4];
  curr_flags = ptr[8];

  ptr += 1;
  if (columni != 3) {
    next_high = ptr[0];
    next_low = ptr[4];
    next_flags = ptr[8];
  } else {
    ptr += 8;
    next_high = ptr[0];
    next_low = ptr[4];
    next_flags = ptr[8];
  }
#endif

  debug2(printf("high:  %08X %08X\n",curr_high,next_high));
  debug2(printf("low:   %08X %08X\n",curr_low,next_low));
  debug2(printf("flags: %08X %08X\n",curr_flags,next_flags));

  rightshift = pos5 % 32;
  leftshift = 32 - rightshift;
  *high = (next_high << leftshift) | (curr_high >> rightshift);
  *low = (next_low << leftshift) | (curr_low >> rightshift);
  *flags = (next_flags << leftshift) | (curr_flags >> rightshift);

  debug2(printf("high:  %08X\n",*high));
  debug2(printf("low:   %08X\n",*low));
  debug2(printf("flags: %08X\n",*flags));

  return;
}



