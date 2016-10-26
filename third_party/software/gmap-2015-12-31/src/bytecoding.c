static char rcsid[] = "$Id: bytecoding.c 170515 2015-07-23 23:03:24Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "bytecoding.h"

#include <stdio.h>
#include <stdlib.h>
#include "mem.h"
#include "fopen.h"

#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
#else
#include "littleendian.h"
#endif


#ifdef DEBUG10
#define debug10(x) x
#else
#define debug10(x)
#endif


void
Bytecoding_write (char *bytesfile, char *excfile, char *guidefile, UINT4 *values,
		  UINT4 genomelength, int guide_interval) {
  FILE *fp_bytes, *fp_guide, *fp_exceptions;
  unsigned char *bytes;
  UINT4 nexceptions = 0;

  UINT4 n = genomelength, i;
  UINT4 guide_value = 0;


  bytes = (unsigned char *) MALLOC((n+1)*sizeof(unsigned char));
  fp_exceptions = FOPEN_WRITE_BINARY(excfile);
  fp_guide = FOPEN_WRITE_BINARY(guidefile);

  for (i = 0; i <= n; i++) {
    if (values[i] < 255) {
      bytes[i] = (unsigned char) values[i];
    } else {
      bytes[i] = (unsigned char) 255; /* Indicates an exception */

      while (i >= guide_value) {
	FWRITE_UINT(nexceptions,fp_guide);
	guide_value += guide_interval;
      }

      FWRITE_UINT(i,fp_exceptions);
      FWRITE_UINT(values[i],fp_exceptions);
      nexceptions++;
    }
  }

#if 0
  /* Overkill */
  while (i >= guide_value) {
    FWRITE_UINT(nexceptions,fp_guide);
    guide_value += guide_interval;
  }
#else
  /* Just need this */
  FWRITE_UINT(nexceptions,fp_guide);
#endif


  fclose(fp_exceptions);
  fclose(fp_guide);

  fprintf(stderr,"Byte-coding: %u values < 255, %u exceptions >= 255 (%.1f%%)\n",
	  (n+1)-nexceptions,nexceptions,100*(double) nexceptions/(double) (n+1));
  fprintf(stderr,"Writing bytes file...");
  fp_bytes = FOPEN_WRITE_BINARY(bytesfile);
  fwrite(bytes,sizeof(unsigned char),n+1,fp_bytes);
  fclose(fp_bytes);
  fprintf(stderr,"done\n");

  FREE(bytes);

  return;
}


unsigned char *
Bytecoding_write_exceptions_only (char *excfile, char *guidefile, UINT4 *values,
				  UINT4 genomelength, int guide_interval) {
  unsigned char *bytes;
  FILE *fp_guide, *fp_exceptions;
  UINT4 nexceptions = 0;

  UINT4 n = genomelength, i;
  UINT4 guide_value = 0;


  bytes = (unsigned char *) MALLOC((n+1)*sizeof(unsigned char));
  fp_exceptions = FOPEN_WRITE_BINARY(excfile);
  fp_guide = FOPEN_WRITE_BINARY(guidefile);

  for (i = 0; i <= n; i++) {
    if (values[i] < 255) {
      bytes[i] = (unsigned char) values[i];
    } else {
      bytes[i] = (unsigned char) 255; /* Indicates an exception */

      while (i >= guide_value) {
	FWRITE_UINT(nexceptions,fp_guide);
	guide_value += guide_interval;
      }

      FWRITE_UINT(i,fp_exceptions);
      FWRITE_UINT(values[i],fp_exceptions);
      nexceptions++;
    }
  }

#if 0
  /* Overkill */
  while (i >= guide_value) {
    FWRITE_UINT(nexceptions,fp_guide);
    guide_value += guide_interval;
  }
#else
  /* Just need this */
  FWRITE_UINT(nexceptions,fp_guide);
#endif


  fclose(fp_exceptions);
  fclose(fp_guide);

  fprintf(stderr,"Byte-coding: %u values < 255, %u exceptions >= 255 (%.1f%%)\n",
	  (n+1)-nexceptions,nexceptions,100*(double) nexceptions/(double) (n+1));

  return bytes;
}


#define LCPCHILDDC_BLOCKSIZE 5
#define BUFFER_NBLOCKS 10000000

/* Interleaved byte array with lcp info, child info, and
   discriminating chars.  Each lcp and child element takes 1 byte, and
   the discriminating chars takes 1 nibble.  Therefore, can store in
   blocks of 5 bytes: lcp0, lcp1, discrim (with nibble1 in bits 4-7
   and nibble0 in bits 0-3), child0, and child1, which is the
   approximate order of information needed by Sarray_get_child.
   Assumes that lcpbytes, lcp_exceptions, and lcp_guide have already
   been computed and written, although this format obviates the
   separate lcpbytes file. */

void
Bytecoding_write_lcpchilddc (char *bytesfile, char *excfile, char *guidefile, UINT4 *child,
			     unsigned char *discrim_chars, unsigned char *lcpbytes,
			     UINT4 genomelength, int guide_interval) {
  FILE *fp_bytes, *fp_guide, *fp_exceptions;
  unsigned char *bytes_buffer, *bytes_ptr;
  UINT4 nexceptions = 0;

  UINT4 n = genomelength, i;
  /* size_t nblocks; */
  int b;
  
  UINT4 guide_value = 0;

  /* nblocks = ((n + 1) + 1)/2; */

  bytes_buffer = (unsigned char *) MALLOC(BUFFER_NBLOCKS * LCPCHILDDC_BLOCKSIZE * sizeof(unsigned char));

  fp_exceptions = FOPEN_WRITE_BINARY(excfile);
  fp_guide = FOPEN_WRITE_BINARY(guidefile);
  fp_bytes = FOPEN_WRITE_BINARY(bytesfile);

  i = 0;
  bytes_ptr = &(bytes_buffer[0]);
  b = 0;
  fprintf(stderr,"Writing file %s",bytesfile);
  while (i + 1 <= n) {
    *bytes_ptr++ = lcpbytes[i];	/* Byte 0 */
    *bytes_ptr++ = lcpbytes[i+1]; /* Byte 1 */

    *bytes_ptr++ = *discrim_chars++; /* Byte 2 */

    if (child[i] < 255) {
      *bytes_ptr++ = (unsigned char) child[i]; /* Byte 3 */
    } else {
      *bytes_ptr++ = (unsigned char) 255; /* Byte 3.  Indicates an exception */

      while (i >= guide_value) {
	FWRITE_UINT(nexceptions,fp_guide);
	guide_value += guide_interval;
      }

      FWRITE_UINT(i,fp_exceptions);
      FWRITE_UINT(child[i],fp_exceptions);
      nexceptions++;
    }
    i++;

    if (child[i] < 255) {
      *bytes_ptr++ = (unsigned char) child[i]; /* Byte 4 */
    } else {
      *bytes_ptr++ = (unsigned char) 255; /* Byte 4.  Indicates an exception */

      while (i >= guide_value) {
	FWRITE_UINT(nexceptions,fp_guide);
	guide_value += guide_interval;
      }

      FWRITE_UINT(i,fp_exceptions);
      FWRITE_UINT(child[i],fp_exceptions);
      nexceptions++;
    }
    i++;

    if (++b >= BUFFER_NBLOCKS) {
      fwrite(bytes_buffer,sizeof(unsigned char),BUFFER_NBLOCKS*LCPCHILDDC_BLOCKSIZE,fp_bytes);
      bytes_ptr = &(bytes_buffer[0]);
      b = 0;
      fprintf(stderr,".");
    }
  }

  if (i <= n) {
    *bytes_ptr++ = lcpbytes[i];	/* Byte 0 */
    *bytes_ptr++ = 0;		/* Byte 1 */

    *bytes_ptr++ = *discrim_chars++; /* Byte 2 */

    if (child[i] < 255) {
      *bytes_ptr++ = (unsigned char) child[i]; /* Byte 3 */
    } else {
      *bytes_ptr++ = (unsigned char) 255; /* Byte 3.  Indicates an exception */

      while (i >= guide_value) {
	FWRITE_UINT(nexceptions,fp_guide);
	guide_value += guide_interval;
      }

      FWRITE_UINT(i,fp_exceptions);
      FWRITE_UINT(child[i],fp_exceptions);
      nexceptions++;
    }

    *bytes_ptr++ = 0x00;	/* Byte 4 */

    if (++b >= BUFFER_NBLOCKS) {
      fwrite(bytes_buffer,sizeof(unsigned char),BUFFER_NBLOCKS*LCPCHILDDC_BLOCKSIZE,fp_bytes);
      bytes_ptr = &(bytes_buffer[0]);
      b = 0;
      fprintf(stderr,".");
    }
  }


#if 0
  /* Overkill */
  while (i >= guide_value) {
    FWRITE_UINT(nexceptions,fp_guide);
    guide_value += guide_interval;
  }
#else
  /* Just need this */
  FWRITE_UINT(nexceptions,fp_guide);
#endif


  if (b > 0) {
    fwrite(bytes_buffer,sizeof(unsigned char),b*LCPCHILDDC_BLOCKSIZE,fp_bytes);
  }
  fprintf(stderr,"done\n");

  fclose(fp_bytes);
  fclose(fp_exceptions);
  fclose(fp_guide);

  fprintf(stderr,"Byte-coding: %u values < 255, %u exceptions >= 255 (%.1f%%)\n",
	  (n+1)-nexceptions,nexceptions,100*(double) nexceptions/(double) (n+1));

  FREE(bytes_buffer);

  return;
}


#define get_bit(i,bitvector) ((bitvector)[(i)/64] & (1UL << ((i)%64)))


#if 0
void
Bytecoding_write_lcpchilddcn (char *bytesfile, char *excfile, char *guidefile, UINT4 *child,
			      unsigned char *discrim_chars, unsigned char *lcpbytes, UINT8 *predictive_nextp,
			      UINT4 genomelength, int guide_interval) {
  FILE *fp_bytes, *fp_guide, *fp_exceptions;
  unsigned char *bytes, *bytes_orig;
  UINT4 nexceptions = 0;

  UINT4 n = genomelength, i;
  size_t nblocks;
  
  UINT4 guide_value = 0;

  nblocks = ((n + 1) + 1)/2;

  bytes_orig = bytes = (unsigned char *) MALLOC(nblocks * LCPCHILDDC_BLOCKSIZE * sizeof(unsigned char));
  fp_exceptions = FOPEN_WRITE_BINARY(excfile);
  fp_guide = FOPEN_WRITE_BINARY(guidefile);

  i = 0;
  while (i + 1 <= n) {
    *bytes++ = lcpbytes[i];
    *bytes++ = lcpbytes[i+1];

    *bytes++ = *discrim_chars++;

    if (child[i] < 127) {
      if (get_bit(i,predictive_nextp)) {
	*bytes++ = (unsigned char) child[i] | 0x80;
      } else {
	*bytes++ = (unsigned char) child[i];
      }
    } else {
      if (get_bit(i,predictive_nextp)) {
	*bytes++ = (unsigned char) 127 | 0x80; /* Indicates an exception */
      } else {
	*bytes++ = (unsigned char) 127; /* Indicates an exception */
      }

      while (i >= guide_value) {
	FWRITE_UINT(nexceptions,fp_guide);
	guide_value += guide_interval;
      }

      FWRITE_UINT(i,fp_exceptions);
      FWRITE_UINT(child[i],fp_exceptions);
      nexceptions++;
    }
    i++;

    if (child[i] < 127) {
      if (get_bit(i,predictive_nextp)) {
	*bytes++ = (unsigned char) child[i] | 0x80;
      } else {
	*bytes++ = (unsigned char) child[i];
      }
    } else {
      if (get_bit(i,predictive_nextp)) {
	*bytes++ = (unsigned char) 127 | 0x80; /* Indicates an exception */
      } else {
	*bytes++ = (unsigned char) 127; /* Indicates an exception */
      }

      while (i >= guide_value) {
	FWRITE_UINT(nexceptions,fp_guide);
	guide_value += guide_interval;
      }

      FWRITE_UINT(i,fp_exceptions);
      FWRITE_UINT(child[i],fp_exceptions);
      nexceptions++;
    }
    i++;

  }

  if (i <= n) {
    *bytes++ = lcpbytes[i];
    *bytes++ = 0;

    *bytes++ = *discrim_chars++;

    if (child[i] < 127) {
      if (get_bit(i,predictive_nextp)) {
	*bytes++ = (unsigned char) child[i] | 0x80;
      } else {
	*bytes++ = (unsigned char) child[i];
      }
    } else {
      if (get_bit(i,predictive_nextp)) {
	*bytes++ = (unsigned char) 127 | 0x80; /* Indicates an exception */
      } else {
	*bytes++ = (unsigned char) 127; /* Indicates an exception */
      }

      while (i >= guide_value) {
	FWRITE_UINT(nexceptions,fp_guide);
	guide_value += guide_interval;
      }

      FWRITE_UINT(i,fp_exceptions);
      FWRITE_UINT(child[i],fp_exceptions);
      nexceptions++;
    }

    *bytes++ = 0x00;
  }


#if 0
  /* Overkill */
  while (i >= guide_value) {
    FWRITE_UINT(nexceptions,fp_guide);
    guide_value += guide_interval;
  }
#else
  /* Just need this */
  FWRITE_UINT(nexceptions,fp_guide);
#endif


  fclose(fp_exceptions);
  fclose(fp_guide);

  fprintf(stderr,"Byte-coding: %u values < 127, %u exceptions >= 127 (%.1f%%)\n",
	  (n+1)-nexceptions,nexceptions,100*(double) nexceptions/(double) (n+1));
  fprintf(stderr,"Writing bytes file...");
  fp_bytes = FOPEN_WRITE_BINARY(bytesfile);
  fwrite(bytes_orig,sizeof(unsigned char),nblocks*LCPCHILDDC_BLOCKSIZE,fp_bytes);
  fclose(fp_bytes);
  fprintf(stderr,"done\n");

  FREE(bytes_orig);

  return;
}
#endif


UINT4
Bytecoding_read (UINT4 key, unsigned char *bytes, UINT4 *exceptions, int nexceptions) {
  unsigned char byte;
  UINT4 lowi, middlei, highi;

  if ((byte = bytes[key]) < 255) {
    debug10(printf("value %d < 255\n",byte));
    return (UINT4) byte;
  } else {

    lowi = 0;
    highi = nexceptions;
    debug10(printf("entered binary search with lowi=%d, highi=%d, goal=%u\n",lowi,highi,key));
    
    while (lowi < highi) {
      middlei = lowi + ((highi - lowi) / 2);
#ifdef WORDS_BIGENDIAN
      debug10(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		     lowi,exceptions[2*lowi],middlei,exceptions[2*middlei],
		     highi,exceptions[2*highi],key));
      if (key < Bigendian_convert_uint(exceptions[2*middlei])) {
	highi = middlei;
      } else if (key > Bigendian_convert_uint(exceptions[2*middlei])) {
	lowi = middlei + 1;
      } else {
	debug10(printf("binary search returns %d => %u\n",middlei,exceptions[2*middlei+1]));
	return Bigendian_convert_uint(exceptions[2*middlei+1]);
      }
#else
      debug10(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		     lowi,exceptions[2*lowi],middlei,exceptions[2*middlei],
		     highi,exceptions[2*highi],key));
      if (key < exceptions[2*middlei]) {
	highi = middlei;
      } else if (key > exceptions[2*middlei]) {
	lowi = middlei + 1;
      } else {
	debug10(printf("binary search returns %d => %u\n",middlei,exceptions[2*middlei+1]));
	return exceptions[2*middlei+1];
      }
#endif
    }

    /* debug10(printf("binary search returns %d => %u\n",highi,exceptions[highi+1])); */
    /* return exceptions[highi + 1]; */

    fprintf(stderr,"Bytecoding_read should have found index %u as an exception, but failed\n",key);
    abort();
  }
}

UINT4
Bytecoding_read_wguide (UINT4 key, unsigned char *bytes, UINT4 *guide, UINT4 *exceptions, int guide_interval) {
  unsigned char byte;
  UINT4 lowi, middlei, highi;
  UINT4 guidei;

  if ((byte = bytes[key]) < 255) {
    debug10(printf("value %d < 255\n",byte));
    return (UINT4) byte;

  } else {
    guidei = key/guide_interval;
#ifdef WORDS_BIGENDIAN
    lowi = Bigendian_convert_uint(guide[guidei]);
    highi = Bigendian_convert_uint(guide[guidei+1]);
#else
    lowi = guide[guidei];
    highi = guide[guidei+1];
#endif

    debug10(printf("entered binary search with lowi=%d, highi=%d, goal=%u\n",lowi,highi,key));
    
    while (lowi < highi) {
      middlei = lowi + ((highi - lowi) / 2);
#ifdef WORDS_BIGENDIAN
      debug10(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		     lowi,exceptions[2*lowi],middlei,exceptions[2*middlei],
		     highi,exceptions[2*highi],key));
      if (key < Bigendian_convert_uint(exceptions[2*middlei])) {
	highi = middlei;
      } else if (key > Bigendian_convert_uint(exceptions[2*middlei])) {
	lowi = middlei + 1;
      } else {
	debug10(printf("binary search returns %d => %u\n",middlei,exceptions[2*middlei+1]));
	return Bigendian_convert_uint(exceptions[2*middlei+1]);
      }
#else
      debug10(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		     lowi,exceptions[2*lowi],middlei,exceptions[2*middlei],
		     highi,exceptions[2*highi],key));
      if (key < exceptions[2*middlei]) {
	highi = middlei;
      } else if (key > exceptions[2*middlei]) {
	lowi = middlei + 1;
      } else {
	debug10(printf("binary search returns %d => %u\n",middlei,exceptions[2*middlei+1]));
	return exceptions[2*middlei+1];
      }
#endif
    }

    /* debug10(printf("binary search returns %d => %u\n",highi,exceptions[highi+1])); */
    /* return exceptions[highi + 1]; */

    fprintf(stderr,"Bytecoding_read_wguide should have found index %u as an exception, but failed\n",key);
    abort();
  }
}
  

UINT4
Bytecoding_lcpchilddc_lcp (UINT4 key, unsigned char *bytes, UINT4 *exceptions, int nexceptions) {
  UINT8 blocki = key/2;		/* Needs to be UINT8, because 5 * 2^32 will overflow UINT4 */
  unsigned char *block = &(bytes[blocki * LCPCHILDDC_BLOCKSIZE]);
  unsigned char byte;
  UINT4 lowi, middlei, highi;

  if ((byte = block[0 + (key % 2)]) < 255) {
    debug10(printf("value %d < 255\n",byte));
    return (UINT4) byte;

  } else {
    lowi = 0;
    highi = nexceptions;
    debug10(printf("entered binary search with lowi=%d, highi=%d, goal=%u\n",lowi,highi,key));
    
    while (lowi < highi) {
      middlei = lowi + ((highi - lowi) / 2);
#ifdef WORDS_BIGENDIAN
      debug10(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		     lowi,exceptions[2*lowi],middlei,exceptions[2*middlei],
		     highi,exceptions[2*highi],key));
      if (key < Bigendian_convert_uint(exceptions[2*middlei])) {
	highi = middlei;
      } else if (key > Bigendian_convert_uint(exceptions[2*middlei])) {
	lowi = middlei + 1;
      } else {
	debug10(printf("binary search returns %d => %u\n",middlei,exceptions[2*middlei+1]));
	return Bigendian_convert_uint(exceptions[2*middlei+1]);
      }
#else
      debug10(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		     lowi,exceptions[2*lowi],middlei,exceptions[2*middlei],
		     highi,exceptions[2*highi],key));
      if (key < exceptions[2*middlei]) {
	highi = middlei;
      } else if (key > exceptions[2*middlei]) {
	lowi = middlei + 1;
      } else {
	debug10(printf("binary search returns %d => %u\n",middlei,exceptions[2*middlei+1]));
	return exceptions[2*middlei+1];
      }
#endif
    }

    /* debug10(printf("binary search returns %d => %u\n",highi,exceptions[highi+1])); */
    /* return exceptions[highi + 1]; */

    fprintf(stderr,"Bytecoding_lcp should have found index %u as an exception, but failed\n",key);
    abort();
  }
}


/*                                      0   1   2   3   4   5   6   7   8   9   A   B   C   D   E   F */
static char discrim_char_before[16] = {'?','$','$','$','$','$','A','A','A','A','C','C','C','G','G','T'};
static char discrim_char_after[16]  = {'?','A','C','G','T','X','C','G','T','X','G','T','X','T','X','X'};


char
Bytecoding_lcpchilddc_dc (char *c1, UINT4 key, unsigned char *bytes) {
  UINT8 blocki = key/2;		/* Needs to be UINT8, because 5 * 2^32 will overflow UINT4 */
  unsigned char *block = &(bytes[blocki * LCPCHILDDC_BLOCKSIZE]);
  /* int pos = key % 2; */

  unsigned char nibble;

  nibble = (block[2] >> (4 * (key % 2))) & 0x0F;
  *c1 = discrim_char_before[nibble];
  return discrim_char_after[nibble]; /* c2 */
}


UINT4
Bytecoding_lcpchilddc_child_up (UINT4 key, unsigned char *bytes, UINT4 *guide, UINT4 *exceptions, int guide_interval) {
  UINT8 blocki = key/2;		/* Needs to be UINT8, because 5 * 2^32 will overflow UINT4 */
  unsigned char *block = &(bytes[blocki * LCPCHILDDC_BLOCKSIZE]);
  unsigned char byte;
  UINT4 lowi, middlei, highi;
  UINT4 guidei;

  if ((byte = block[3 + (key % 2)]) < 255) {
    debug10(printf("value %d < 255\n",byte));
    return key - (UINT4) byte;

  } else {
    guidei = key/guide_interval;
#ifdef WORDS_BIGENDIAN
    lowi = Bigendian_convert_uint(guide[guidei]);
    highi = Bigendian_convert_uint(guide[guidei+1]);
#else
    lowi = guide[guidei];
    highi = guide[guidei+1];
#endif

    debug10(printf("entered binary search with lowi=%d, highi=%d, goal=%u\n",lowi,highi,key));
    
    while (lowi < highi) {
      middlei = lowi + ((highi - lowi) / 2);
#ifdef WORDS_BIGENDIAN
      debug10(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		     lowi,exceptions[2*lowi],middlei,exceptions[2*middlei],
		     highi,exceptions[2*highi],key));
      if (key < Bigendian_convert_uint(exceptions[2*middlei])) {
	highi = middlei;
      } else if (key > Bigendian_convert_uint(exceptions[2*middlei])) {
	lowi = middlei + 1;
      } else {
	debug10(printf("binary search returns %d => %u\n",middlei,exceptions[2*middlei+1]));
	return key - Bigendian_convert_uint(exceptions[2*middlei+1]);
      }
#else
      debug10(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		     lowi,exceptions[2*lowi],middlei,exceptions[2*middlei],
		     highi,exceptions[2*highi],key));
      if (key < exceptions[2*middlei]) {
	highi = middlei;
      } else if (key > exceptions[2*middlei]) {
	lowi = middlei + 1;
      } else {
	debug10(printf("binary search returns %d => %u\n",middlei,exceptions[2*middlei+1]));
	return key - exceptions[2*middlei+1];
      }
#endif
    }

    /* debug10(printf("binary search returns %d => %u\n",highi,exceptions[highi+1])); */
    /* return exceptions[highi + 1]; */

    fprintf(stderr,"Bytecoding_lcpchilddc_child_up should have found index %u as an exception, but failed\n",key);
    abort();
  }
}

UINT4
Bytecoding_lcpchilddc_child_next (UINT4 key, unsigned char *bytes, UINT4 *guide, UINT4 *exceptions, int guide_interval) {
  UINT8 blocki = key/2;		/* Needs to be UINT8, because 5 * 2^32 will overflow UINT4 */
  unsigned char *block = &(bytes[blocki * LCPCHILDDC_BLOCKSIZE]);
  unsigned char byte;
  UINT4 lowi, middlei, highi;
  UINT4 guidei;

  if ((byte = block[3 + (key % 2)]) < 255) {
    debug10(printf("value %d < 255\n",byte));
    return (UINT4) byte + key + 1;

  } else {
    guidei = key/guide_interval;
#ifdef WORDS_BIGENDIAN
    lowi = Bigendian_convert_uint(guide[guidei]);
    highi = Bigendian_convert_uint(guide[guidei+1]);
#else
    lowi = guide[guidei];
    highi = guide[guidei+1];
#endif

    debug10(printf("entered binary search with lowi=%d, highi=%d, goal=%u\n",lowi,highi,key));
    
    while (lowi < highi) {
      middlei = lowi + ((highi - lowi) / 2);
#ifdef WORDS_BIGENDIAN
      debug10(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		     lowi,exceptions[2*lowi],middlei,exceptions[2*middlei],
		     highi,exceptions[2*highi],key));
      if (key < Bigendian_convert_uint(exceptions[2*middlei])) {
	highi = middlei;
      } else if (key > Bigendian_convert_uint(exceptions[2*middlei])) {
	lowi = middlei + 1;
      } else {
	debug10(printf("binary search returns %d => %u\n",middlei,exceptions[2*middlei+1]));
	return Bigendian_convert_uint(exceptions[2*middlei+1]) + key + 1;
      }
#else
      debug10(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		     lowi,exceptions[2*lowi],middlei,exceptions[2*middlei],
		     highi,exceptions[2*highi],key));
      if (key < exceptions[2*middlei]) {
	highi = middlei;
      } else if (key > exceptions[2*middlei]) {
	lowi = middlei + 1;
      } else {
	debug10(printf("binary search returns %d => %u\n",middlei,exceptions[2*middlei+1]));
	return exceptions[2*middlei+1] + key + 1;
      }
#endif
    }

    /* debug10(printf("binary search returns %d => %u\n",highi,exceptions[highi+1])); */
    /* return exceptions[highi + 1]; */

    fprintf(stderr,"Bytecoding_lcpchilddc_child_next should have found index %u as an exception, but failed\n",key);
    abort();
  }
}


UINT4
Bytecoding_lcpchilddc_lcp_next (UINT4 *child_next, UINT4 key, unsigned char *bytes, UINT4 *child_guide,
				UINT4 *child_exceptions, int child_guide_interval,
				UINT4 *lcp_exceptions, int n_lcp_exceptions) {
  UINT8 blocki = key/2;		/* Needs to be UINT8, because 5 * 2^32 will overflow UINT4 */
  unsigned char *block = &(bytes[blocki * LCPCHILDDC_BLOCKSIZE]);
  unsigned char byte;
  UINT4 lowi, middlei, highi;
  UINT4 guidei;

  if ((byte = block[3 + (key % 2)]) < 255) {
    debug10(printf("value %d < 255\n",byte));
    *child_next = (UINT4) byte + key + 1;
    return Bytecoding_lcpchilddc_lcp(*child_next,bytes,lcp_exceptions,n_lcp_exceptions);

  } else {
    guidei = key/child_guide_interval;
#ifdef WORDS_BIGENDIAN
    lowi = Bigendian_convert_uint(child_guide[guidei]);
    highi = Bigendian_convert_uint(child_guide[guidei+1]);
#else
    lowi = child_guide[guidei];
    highi = child_guide[guidei+1];
#endif

    debug10(printf("entered binary search with lowi=%d, highi=%d, goal=%u\n",lowi,highi,key));
    
    while (lowi < highi) {
      middlei = lowi + ((highi - lowi) / 2);
#ifdef WORDS_BIGENDIAN
      debug10(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		     lowi,child_exceptions[2*lowi],middlei,child_exceptions[2*middlei],
		     highi,child_exceptions[2*highi],key));
      if (key < Bigendian_convert_uint(child_exceptions[2*middlei])) {
	highi = middlei;
  } else if (key > Bigendian_convert_uint(child_exceptions[2*middlei])) {
	lowi = middlei + 1;
      } else {
	debug10(printf("binary search returns %d => %u\n",middlei,child_exceptions[2*middlei+1]));
	*child_next = Bigendian_convert_uint(child_exceptions[2*middlei+1]) + key + 1;
	return Bytecoding_lcpchilddc_lcp(*child_next,bytes,lcp_exceptions,n_lcp_exceptions);
      }
#else
      debug10(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		     lowi,child_exceptions[2*lowi],middlei,child_exceptions[2*middlei],
		     highi,child_exceptions[2*highi],key));
      if (key < child_exceptions[2*middlei]) {
	highi = middlei;
      } else if (key > child_exceptions[2*middlei]) {
	lowi = middlei + 1;
      } else {
	debug10(printf("binary search returns %d => %u\n",middlei,child_exceptions[2*middlei+1]));
	*child_next = child_exceptions[2*middlei+1] + key + 1;
	return Bytecoding_lcpchilddc_lcp(*child_next,bytes,lcp_exceptions,n_lcp_exceptions);
      }
#endif
    }

    /* debug10(printf("binary search returns %d => %u\n",highi,exceptions[highi+1])); */
    /* return exceptions[highi + 1]; */

    fprintf(stderr,"Bytecoding_lcpchilddc_lcp_next should have found index %u as an exception, but failed\n",key);
    abort();
  }
}


#if 0
UINT4
Bytecoding_lcpchilddcn_child_up (bool *nextp, UINT4 key, unsigned char *bytes, UINT4 *guide, UINT4 *exceptions, int guide_interval) {
  UINT8 blocki = key/2;		/* Needs to be UINT8, because 5 * 2^32 will overflow UINT4 */
  unsigned char *block = &(bytes[blocki * LCPCHILDDC_BLOCKSIZE]);
  unsigned char byte;
  UINT4 lowi, middlei, highi;
  UINT4 guidei;

  byte = block[3 + (key % 2)];
  if (byte & 0x80) {
    *nextp = true;
  } else {
    *nextp = false;
  }
  byte &= 0x7F;

  if (byte < 127) {
    debug10(printf("value %d < 127\n",byte));
    return key - (UINT4) byte;

  } else {
    guidei = key/guide_interval;
#ifdef WORDS_BIGENDIAN
    lowi = Bigendian_convert_uint(guide[guidei]);
    highi = Bigendian_convert_uint(guide[guidei+1]);
#else
    lowi = guide[guidei];
    highi = guide[guidei+1];
#endif

    debug10(printf("entered binary search with lowi=%d, highi=%d, goal=%u\n",lowi,highi,key));
    
    while (lowi < highi) {
      middlei = lowi + ((highi - lowi) / 2);
#ifdef WORDS_BIGENDIAN
      debug10(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		     lowi,exceptions[2*lowi],middlei,exceptions[2*middlei],
		     highi,exceptions[2*highi],key));
      if (key < Bigendian_convert_uint(exceptions[2*middlei])) {
	highi = middlei;
     } else if (key > Bigendian_convert_uint(exceptions[2*middlei])) {
	lowi = middlei + 1;
      } else {
	debug10(printf("binary search returns %d => %u\n",middlei,exceptions[2*middlei+1]));
	return key - Bigendian_convert_uint(exceptions[2*middlei+1]);
      }
#else
      debug10(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		     lowi,exceptions[2*lowi],middlei,exceptions[2*middlei],
		     highi,exceptions[2*highi],key));
      if (key < exceptions[2*middlei]) {
	highi = middlei;
      } else if (key > exceptions[2*middlei]) {
	lowi = middlei + 1;
      } else {
	debug10(printf("binary search returns %d => %u\n",middlei,exceptions[2*middlei+1]));
	return key - exceptions[2*middlei+1];
      }
#endif
    }

    /* debug10(printf("binary search returns %d => %u\n",highi,exceptions[highi+1])); */
    /* return exceptions[highi + 1]; */

    fprintf(stderr,"Bytecoding_lcpchilddcn_child_up should have found index %u as an exception, but failed\n",key);
    abort();
  }
}
#endif


#if 0
UINT4
Bytecoding_lcpchilddcn_child_next (bool *nextp, UINT4 key, unsigned char *bytes, UINT4 *guide, UINT4 *exceptions, int guide_interval) {
  UINT8 blocki = key/2;		/* Needs to be UINT8, because 5 * 2^32 will overflow UINT4 */
  unsigned char *block = &(bytes[blocki * LCPCHILDDC_BLOCKSIZE]);
  unsigned char byte;
  UINT4 lowi, middlei, highi;
  UINT4 guidei;

  byte = block[3 + (key % 2)];
  if (byte & 0x80) {
    *nextp = true;
  } else {
    *nextp = false;
  }
  byte &= 0x7F;

  if (byte < 127) {
    debug10(printf("value %d < 127\n",byte));
    return (UINT4) byte + key + 1;

  } else {
    guidei = key/guide_interval;
#ifdef WORDS_BIGENDIAN
    lowi = Bigendian_convert_uint(guide[guidei]);
    highi = Bigendian_convert_uint(guide[guidei+1]);
#else
    lowi = guide[guidei];
    highi = guide[guidei+1];
#endif

    debug10(printf("entered binary search with lowi=%d, highi=%d, goal=%u\n",lowi,highi,key));
    
    while (lowi < highi) {
      middlei = lowi + ((highi - lowi) / 2);
#ifdef WORDS_BIGENDIAN
      debug10(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		     lowi,exceptions[2*lowi],middlei,exceptions[2*middlei],
		     highi,exceptions[2*highi],key));
      if (key < Bigendian_convert_uint(exceptions[2*middlei])) {
	highi = middlei;
      } else if (key > Bigendian_convert_uint(exceptions[2*middlei])) {
	lowi = middlei + 1;
      } else {
	debug10(printf("binary search returns %d => %u\n",middlei,exceptions[2*middlei+1]));
	return Bigendian_convert_uint(exceptions[2*middlei+1]) + key + 1;
      }
#else
      debug10(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		     lowi,exceptions[2*lowi],middlei,exceptions[2*middlei],
		     highi,exceptions[2*highi],key));
      if (key < exceptions[2*middlei]) {
	highi = middlei;
      } else if (key > exceptions[2*middlei]) {
	lowi = middlei + 1;
      } else {
	debug10(printf("binary search returns %d => %u\n",middlei,exceptions[2*middlei+1]));
	return exceptions[2*middlei+1] + key + 1;
      }
#endif
    }

    /* debug10(printf("binary search returns %d => %u\n",highi,exceptions[highi+1])); */
    /* return exceptions[highi + 1]; */

    fprintf(stderr,"Bytecoding_lcpchilddcn_child_next should have found index %u as an exception, but failed\n",key);
    abort();
  }
}
#endif



