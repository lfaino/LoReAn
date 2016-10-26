/* $Id: types.h 168395 2015-06-26 17:13:13Z twu $ */
#ifndef TYPES_INCLUDED
#define TYPES_INCLUDED
#ifdef HAVE_CONFIG_H
#include <config.h>		/* For SIZEOF_UNSIGNED_LONG_LONG, SIZEOF_UNSIGNED_LONG, HAVE_64_BIT */
#endif

/* Number of bits, such as index1part or basesize.  Need to allow for negative values. */
typedef int Width_T;

/* Number of entries, such as offsetscomp_blocksize */
typedef unsigned int Blocksize_T;


/* A 2-byte word */
typedef unsigned short UINT2;

/* A 4-byte word */
typedef unsigned int UINT4;
typedef int INT4;


/* Compressed representation of genome (high, low, flags).  Always
   UINT4.  Can think of as a genome block unit.  */
typedef UINT4 Genomecomp_T;


/* An 8-byte word */
/* Oligospace_T needs to hold 1 more than maximum Storedoligomer_T.
   If 8-byte words are not available, then maximum k-mer is 15 */
/* Prefer to use unsigned long long, whic should be 8 bytes on all systems */
#if (SIZEOF_UNSIGNED_LONG_LONG == 8)
#define HAVE_64_BIT 1
#define MAXIMUM_KMER 16
typedef unsigned long long UINT8;
typedef unsigned long long Oligospace_T;

#elif (SIZEOF_UNSIGNED_LONG == 8)
#define HAVE_64_BIT 1
#define MAXIMUM_KMER 16
typedef unsigned long UINT8;
typedef unsigned long Oligospace_T;

#else
#define MAXIMUM_KMER 15
#define OLIGOSPACE_NOT_LONG
typedef unsigned int Oligospace_T;

#endif

/* Contents of compressed offsets file.  Storing as UINT4, even for
   large genomes, to reduce zero-padding of bitstreams.  For large
   genomes, need to store 64-bit Positionsptr_T quantity in 2 UINT4
   words. */
typedef UINT4 Offsetscomp_T;

/* Holds a k-mer.  Can be UINT4 as long as k <= 16. */
/* Some procedures use Shortoligomer_T, which should be the same */
typedef UINT4 Storedoligomer_T;


/* Definitions */
/* Large genome: Genomic length > 2^32, needing 8-byte Univcoord_T */
/* Huge genome: Entries in positions file > 2^32, needing 8-byte Positionsptr_T */

/* An offset into the positions file of an IndexDB.  For small genomes
   < 2^32 bp such as human, need 3 billion divided by sampling
   interval (default 3), requiring a maximum of 32 bits or 4 bytes
   (Positionsptr_T).  For huge genomes or more frequent sampling,
   need 8 bytes, or Hugepositionsptr_T. */
#ifdef HAVE_64_BIT

#ifdef UTILITYP
typedef UINT8 Hugepositionsptr_T;
typedef UINT4 Positionsptr_T;
#elif defined(LARGE_GENOMES)
/* Don't really need offsets to be 8-byte unless we have a huge
   genome, but this simplifies the code */
typedef UINT8 Hugepositionsptr_T;
typedef UINT8 Positionsptr_T;
#else
typedef UINT4 Positionsptr_T;
#endif

#else
typedef UINT4 Positionsptr_T;
#endif


/* For definition of Chrpos_T, see genomicpos.h */

/* For intervals and IIT files */
#ifdef HAVE_64_BIT

#ifdef UTILITYP
#include "uint8list.h"
typedef UINT8 Univcoord_T;
typedef Uint8list_T Univcoordlist_T;

#elif defined(LARGE_GENOMES)
#include "uint8list.h"
typedef UINT8 Univcoord_T;
typedef Uint8list_T Univcoordlist_T;

#else
#include "uintlist.h"
typedef UINT4 Univcoord_T;
typedef Uintlist_T Univcoordlist_T;

#endif

#else
#include "uintlist.h"
typedef UINT4 Univcoord_T;
typedef Uintlist_T Univcoordlist_T;

#endif

/* For univintervals and Univ_IIT (chromosome_iit) files.  Use the largest word size allowable on the machine.  */
#ifdef HAVE_64_BIT
typedef UINT8 Univ_IIT_coord_T;
#else
typedef UINT4 Univ_IIT_coord_T;
#endif

/* For splicetrie */
typedef UINT4 Trieoffset_T;
typedef UINT4 Triecontent_T;

/* For suffix array */
typedef UINT4 Sarrayptr_T;

#endif

