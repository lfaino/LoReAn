static char rcsid[] = "$Id: indexdb_hr.c 168395 2015-06-26 17:13:13Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
# define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
# define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif

#include "indexdb_hr.h"
#include "indexdbdef.h"
#include "genome128_hr.h"
#include "bitpack64-read.h"
#include "bitpack64-readtwo.h"


#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
/* #define CONVERT_TO_LITTLEENDIAN 1 */
/* Because only call is to generate intersection_diagonals, which are assumed to be in native form */
#else
#include "littleendian.h"
#endif

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>		/* For memcpy */
#include "mem.h"
#include "listdef.h"


/* ALLOW_DUPLICATES is possible only if we permit alternative strains */


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* results of merge_batches */
#ifdef DEBUG0
#define debug0(x) x
#else
#define debug0(x)
#endif

/* merge_batches */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* binary_search, identify_doubles */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* heapify */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif

/* Compoundpos_intersect */
#ifdef DEBUG4
#define debug4(x) x
#else
#define debug4(x)
#endif

/* Compoundpos_find */
#ifdef DEBUG6
#define debug6(x) x
#else
#define debug6(x)
#endif

/* straddling at beginning of genome.  May want to turn on DEBUG11 in stage1hr.c */
#ifdef DEBUG11
#define debug11(x) x
#else
#define debug11(x)
#endif


#define T Indexdb_T

static int index1part;		/* For debugging */
static unsigned int kmer_mask;	/* Was LOW12MER     0x00FFFFFF */
#define right_subst  0x00000001
static unsigned int left_subst; /* Was LEFT_SUBST   0x00100000 */
static unsigned int top_subst;  /* Was TOP_SUBST    0x00400000 */

void
Indexdb_hr_setup (int index1part_in) {
  index1part = index1part_in;
  kmer_mask = ~(~0UL << 2*index1part);

  top_subst = (1 << 2*(index1part-1));
  left_subst = (1 << 2*(index1part-2));

  return;
}



typedef struct Batch_T *Batch_T;
struct Batch_T {
  int nentries;
  Univcoord_T position;
#ifdef LARGE_GENOMES
  unsigned char *positionptr_high;
  UINT4 *positionptr_low;
#else
  Univcoord_T *positionptr;
#endif
};

typedef struct Header_T *Header_T;
struct Header_T {
  int heapsize;
  int delta;
};


/* We want to handle 16 nodes.  In the typical heap structure, we need
a node 8 with one left child, node 16, and then sentinels to handle
nodes 17-31.  Instead, we use a different heap structure, where node 1
has only one child, node 2.  Then, each parent node i has a right node
(2*i) and a left node (2*i-1). */

#define PARENT2(i) ((i+1) >> 1)
#define LEFTSIBLING2(i) (i-1)

static const unsigned int PARENT_EVEN[17] =
/* 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 */
  {0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8};

static void
heap_insert_even (Batch_T *heap, int *heapsize, Batch_T elt, Univcoord_T key) {
  unsigned int i, parenti;

  i = ++(*heapsize);
  while (i > 1 && heap[parenti = PARENT_EVEN[i]]->position > key) {
    heap[i] = heap[parenti];
    i = parenti;
  }
  heap[i] = elt;

  return;
}


#define NCASES 48


#if 0
static void
check_heap_even (Batch_T *heap, int heapsize) {
  int i, j;

  for (i = 1; i <= heapsize; i++) {
    if (heap[i]->position > heap[2*i-1]->position) {
      fprintf(stderr,"Failed because position %llu at heap %d is > position %llu at heap %d\n",
	      (unsigned long long) heap[i]->position,i,(unsigned long long) heap[2*i-1]->position,2*i-1);
      for (j = 1; j <= heapsize*2; j++) {
	fprintf(stderr,"%02d %u\n",j,heap[j]->position);
      }
      abort();
    }
    if (heap[i]->position > heap[2*i]->position) {
      fprintf(stderr,"Failed because position %llu at heap %d is > position %llu at heap %d\n",
	      (unsigned long long) heap[i]->position,i,(unsigned long long) heap[2*i]->position,2*i);
      for (j = 1; j <= heapsize*2; j++) {
	fprintf(stderr,"%02d %u\n",j,heap[j]->position);
      }
      abort();
    }
  }
}
#endif

#define READ_THEN_WRITE 1

#ifdef LARGE_GENOMES
static unsigned char sentinel_position_high = (unsigned char) -1;
static UINT4 sentinel_position_low = (UINT4) -1;
#endif


static Univcoord_T *
merge_batches_one_heap_16_existing (int *nmerged, struct Batch_T *batchpool, int nentries, int diagterm) {
  Univcoord_T *positions, *ptr, position, last_position, this_position;
  struct Batch_T sentinel_struct;
  Batch_T batch, sentinel, heap[17];
  int heapsize;
  unsigned int i;
#ifdef READ_THEN_WRITE
  unsigned int smallesti_1, smallesti_2, smallesti;
#else
  unsigned int parenti, smallesti;
#endif

  debug3(printf("starting merge_batches_one_heap_16_existing\n"));

  debug0(int nentries_save = nentries);

  ptr = positions = (Univcoord_T *) CALLOC(nentries,sizeof(Univcoord_T));

  /* Set up heap */
  heapsize = 0;
  for (i = 0; i < 16; i++) {
    batch = &(batchpool[i]);
    if (batch->nentries > 0) {
#ifdef LARGE_GENOMES
      batch->position = (((Univcoord_T) *batch->positionptr_high++) << 32) + (*batch->positionptr_low++);
#elif defined(WORDS_BIGENDIAN)
      batch->position = Bigendian_convert_univcoord(*batch->positionptr++);
#else
      batch->position = *batch->positionptr++;
#endif
      heap_insert_even(heap,&heapsize,batch,batch->position);
    }
  }

  sentinel_struct.position = (Univcoord_T) -1; /* infinity */
#ifdef LARGE_GENOMES
  sentinel_struct.positionptr_high = &sentinel_position_high;
  sentinel_struct.positionptr_low = &sentinel_position_low;
#else
  sentinel_struct.positionptr = &(sentinel_struct.position);
#endif
  sentinel = &sentinel_struct;

  for (i = heapsize+1; i <= 16; i++) {
    heap[i] = sentinel;
  }

  last_position = 0U;
  while (--nentries >= 1) {
    debug3(printf("nentries = %d, top of heap is %u (%d)\n",
		  nentries+1,heap[1]->position,heapsize));

    /* Get minimum */
    batch = heap[1];
#ifdef CONVERT_TO_LITTLEENDIAN
    this_position = Bigendian_convert_univcoord(batch->position) + diagterm;
#else
    this_position = batch->position + diagterm;
#endif
    if (this_position != last_position) {
      *ptr++ = this_position;
    }
    last_position = this_position;

    if (--batch->nentries <= 0) {
      /* Use last batch (or sentinel) in heap for insertion */
      heap[1] = batch = (heapsize == 1) ? sentinel : heap[heapsize];
      heap[heapsize--] = sentinel;

    } else {
      /* Advance heap, and use this batch for insertion */
#ifdef LARGE_GENOMES
      batch->position = (((Univcoord_T) *batch->positionptr_high++) << 32) + (*batch->positionptr_low++);
#elif defined(WORDS_BIGENDIAN)
      batch->position = Bigendian_convert_univcoord(*batch->positionptr++);
#else
      batch->position = *batch->positionptr++;
#endif
    }

    position = batch->position;
    debug3(printf("starting heapify with %u\n",position));

#ifdef READ_THEN_WRITE
    /* Comparison 0/3 */
    debug3(printf("Comparing right %d: %u\n",2,heap[2]->position));
    if (position <= heap[2]->position) {
      debug3(printf("Inserting at 1\n"));
      /* heap[1] = batch; -- not necessary because batch is already at heap[1] */
    } else {
      /* Comparison 1/3 */
      debug3(printf("Comparing left %d/right %d: %u and %u\n",
		    3,4,heap[3]->position,heap[4]->position));
      smallesti = 4 - (heap[3]->position < heap[4]->position);
      if (position <= heap[smallesti]->position) {
	debug3(printf("Inserting at 2\n"));
	heap[1] = heap[2];
	heap[2] = batch;
      } else {
	smallesti_1 = smallesti;
	smallesti <<= 1;
	/* Comparison 2/3 */
	debug3(printf("Comparing left %d/right %d: %u and %u\n",
		      smallesti-1,smallesti,heap[smallesti-1]->position,heap[smallesti]->position));
	smallesti -= (heap[LEFTSIBLING2(smallesti)]->position < heap[smallesti]->position);
	if (position <= heap[smallesti]->position) {
	  debug3(printf("Inserting at %d\n",smallesti_1));
	  heap[1] = heap[2];
	  heap[2] = heap[smallesti_1];
	  heap[smallesti_1] = batch;
	} else {
	  smallesti_2 = smallesti;
	  smallesti <<= 1;
	  /* Comparison 3/3 */
	  debug3(printf("Comparing left %d/right %d: %u and %u\n",
			smallesti-1,smallesti,heap[smallesti-1]->position,heap[smallesti]->position));
	  smallesti -= (heap[LEFTSIBLING2(smallesti)]->position < heap[smallesti]->position);
	  if (position <= heap[smallesti]->position) {
	    debug3(printf("Inserting at %d\n",smallesti_2));
	    heap[1] = heap[2];
	    heap[2] = heap[smallesti_1];
	    heap[smallesti_1] = heap[smallesti_2];
	    heap[smallesti_2] = batch;
	  } else {
	    debug3(printf("Inserting at %d\n",smallesti));
	    heap[1] = heap[2];
	    heap[2] = heap[smallesti_1];
	    heap[smallesti_1] = heap[smallesti_2];
	    heap[smallesti_2] = heap[smallesti];
	    heap[smallesti] = batch;
	  }
	}
      }
    }
#else
    /* Comparison 0/3 */
    debug3(printf("Comparing right %d: %u\n",2,heap[2]->position));
    if (position <= heap[2]->position) {
      debug3(printf("Inserting at 1\n"));
      /* heap[1] = batch; -- not necessary because batch is already at heap[1] */
    } else {
      heap[1] = heap[2];
      /* Comparison 1/3 */
      debug3(printf("Comparing left %d/right %d: %u and %u\n",
		    3,4,heap[3]->position,heap[4]->position));
      smallesti = 4 - (heap[3]->position < heap[4]->position);
      if (position <= heap[smallesti]->position) {
	debug3(printf("Inserting at 2\n"));
	heap[2] = batch;
      } else {
	heap[2] = heap[smallesti];
	parenti = smallesti;
	smallesti <<= 1;
	/* Comparison 2/3 */
	debug3(printf("Comparing left %d/right %d: %u and %u\n",
		      smallesti-1,smallesti,heap[smallesti-1]->position,heap[smallesti]->position));
	smallesti -= (heap[LEFTSIBLING2(smallesti)]->position < heap[smallesti]->position);
	if (position <= heap[smallesti]->position) {
	  debug3(printf("Inserting at %d\n",parenti));
	  heap[parenti] = batch;
	} else {
	  heap[parenti] = heap[smallesti];
	  parenti = smallesti;
	  smallesti <<= 1;
	  /* Comparison 3/3 */
	  debug3(printf("Comparing left %d/right %d: %u and %u\n",
			smallesti-1,smallesti,heap[smallesti-1]->position,heap[smallesti]->position));
	  smallesti -= (heap[LEFTSIBLING2(smallesti)]->position < heap[smallesti]->position);
	  if (position <= heap[smallesti]->position) {
	    debug3(printf("Inserting at %d\n",parenti));
	    heap[parenti] = batch;
	  } else {
	    heap[parenti] = heap[smallesti];
	    debug3(printf("Inserting at %d\n",smallesti));
	    heap[smallesti] = batch;
	  }
	}
      }
    }
#endif
  }

#ifdef CONVERT_TO_LITTLEENDIAN
  this_position = Bigendian_convert_univcoord(heap[1]->position) + diagterm;
#else
  this_position = heap[1]->position + diagterm;
#endif
  if (this_position != last_position) {
    *ptr++ = this_position;
  }

  *nmerged = (ptr - positions);

#if 0
  position = positions[0];
  for (i = 1; i < nentries_save; i++) {
    if (positions[i] <= position) {
      abort();
    }
    position = positions[i];
  }
#endif

  debug0(
	 for (i = 0; i < nentries_save; i++) {
	   printf("%u\n",positions[i]);
	 }
	 printf("\n");
	 )

  return positions;
}


static Univcoord_T *
merge_batches_one_heap_4_existing (int *nmerged, struct Batch_T *batchpool, int nentries, int diagterm) {
  Univcoord_T *positions, *ptr, position, last_position, this_position;
  struct Batch_T sentinel_struct;
  Batch_T batch, sentinel, heap[5];
  int heapsize;
  unsigned int i;
#ifdef READ_THEN_WRITE
  unsigned int smallesti;
#else
  unsigned int parenti, smallesti;
#endif

  debug3(printf("starting merge_batches_one_heap_4_existing\n"));

  debug0(int nentries_save = nentries);

  ptr = positions = (Univcoord_T *) CALLOC(nentries,sizeof(Univcoord_T));

  /* Set up heap */
  heapsize = 0;
  for (i = 0; i < 4; i++) {
    batch = &(batchpool[i]);
    if (batch->nentries > 0) {
#ifdef LARGE_GENOMES
      batch->position = (((Univcoord_T) *batch->positionptr_high++) << 32) + (*batch->positionptr_low++);
#elif defined(WORDS_BIGENDIAN)
      batch->position = Bigendian_convert_univcoord(*batch->positionptr++);
#else
      batch->position = *batch->positionptr++;
#endif
      heap_insert_even(heap,&heapsize,batch,batch->position);
    }
  }

  sentinel_struct.position = (Univcoord_T) -1; /* infinity */
#ifdef LARGE_GENOMES
  sentinel_struct.positionptr_high = &sentinel_position_high;
  sentinel_struct.positionptr_low = &sentinel_position_low;
#else
  sentinel_struct.positionptr = &(sentinel_struct.position);
#endif
  sentinel = &sentinel_struct;

  for (i = heapsize+1; i <= 4; i++) {
    heap[i] = sentinel;
  }

  last_position = 0U;
  while (--nentries >= 1) {
    debug3(printf("nentries = %d, top of heap is %u (%d)\n",
		  nentries+1,heap[1]->position,heapsize));

    /* Get minimum */
    batch = heap[1];
#ifdef CONVERT_TO_LITTLEENDIAN
    this_position = Bigendian_convert_univcoord(batch->position) + diagterm;
#else
    this_position = batch->position + diagterm;
#endif
    if (this_position != last_position) {
      *ptr++ = this_position;
    }
    last_position = this_position;


    if (--batch->nentries <= 0) {
      /* Use last batch (or sentinel) in heap for insertion */
      heap[1] = batch = (heapsize == 1) ? sentinel : heap[heapsize];
      heap[heapsize--] = sentinel;

    } else {
      /* Advance heap, and use this batch for insertion */
#ifdef LARGE_GENOMES
      batch->position = (((Univcoord_T) *batch->positionptr_high++) << 32) + (*batch->positionptr_low++);
#elif defined(WORDS_BIGENDIAN)
      batch->position = Bigendian_convert_univcoord(*batch->positionptr++);
#else
      batch->position = *batch->positionptr++;
#endif
    }

    position = batch->position;
    debug3(printf("starting heapify with %u\n",position));

#ifdef READ_THEN_WRITE
    /* Comparison 0/3 */
    debug3(printf("Comparing right %d: %u\n",2,heap[2]->position));
    if (position <= heap[2]->position) {
      debug3(printf("Inserting at 1\n"));
      /* heap[1] = batch; -- not necessary because batch is already at heap[1] */
    } else {
      /* Comparison 1/3 */
      debug3(printf("Comparing left %d/right %d: %u and %u\n",
		    3,4,heap[3]->position,heap[4]->position));
      smallesti = 4 - (heap[3]->position < heap[4]->position);
      if (position <= heap[smallesti]->position) {
	debug3(printf("Inserting at 2\n"));
	heap[1] = heap[2];
	heap[2] = batch;
      } else {
	debug3(printf("Inserting at %d\n",smallesti));
	heap[1] = heap[2];
	heap[2] = heap[smallesti];
	heap[smallesti] = batch;
      }
    }

#else
    /* Comparison 0/3 */
    debug3(printf("Comparing right %d: %u\n",2,heap[2]->position));
    if (position <= heap[2]->position) {
      debug3(printf("Inserting at 1\n"));
      /* heap[1] = batch; -- not necessary because batch is already at heap[1] */
    } else {
      heap[1] = heap[2];
      /* Comparison 1/3 */
      debug3(printf("Comparing left %d/right %d: %u and %u\n",
		    3,4,heap[3]->position,heap[4]->position));
      smallesti = 4 - (heap[3]->position < heap[4]->position);
      if (position <= heap[smallesti]->position) {
	debug3(printf("Inserting at 2\n"));
	heap[2] = batch;
      } else {
	heap[2] = heap[smallesti];
	heap[smallesti] = batch;
      }
    }

#endif
  }

#ifdef CONVERT_TO_LITTLEENDIAN
  this_position = Bigendian_convert_univcoord(heap[1]->position) + diagterm;
#else
  this_position = heap[1]->position + diagterm;
#endif
  if (this_position != last_position) {
    *ptr++ = this_position;
  }

  *nmerged = (ptr - positions);

#if 0
  position = positions[0];
  for (i = 1; i < nentries_save; i++) {
    if (positions[i] <= position) {
      abort();
    }
    position = positions[i];
  }
#endif

  debug0(
	 for (i = 0; i < nentries_save; i++) {
	   printf("%u\n",positions[i]);
	 }
	 printf("\n");
	 )


  return positions;
}


/************************************************************************
 *  The following positions functions are take from indexdb.c
 ************************************************************************/

#ifndef LARGE_GENOMES
static void
positions_move_absolute (int positions_fd, Positionsptr_T ptr) {
  off_t offset = ptr*((off_t) sizeof(Univcoord_T));

  if (lseek(positions_fd,offset,SEEK_SET) < 0) {
    fprintf(stderr,"Attempted to do lseek on offset %zd*%d=%zd\n",
	    ptr,(int) sizeof(Univcoord_T),offset);
    perror("Error in indexdb.c, positions_move_absolute_4");
    exit(9);
  }
  return;
}

static void
positions_read_multiple (int positions_fd, Univcoord_T *values, int n) {
  int i;
  Univcoord_T value;
  unsigned char buffer[4];

#ifdef WORDS_BIGENDIAN
  /* Need to keep in bigendian format */
  for (i = 0; i < n; i++) {
    read(positions_fd,buffer,4);

    value = (buffer[0] & 0xff);
    value <<= 8;
    value |= (buffer[1] & 0xff);
    value <<= 8;
    value |= (buffer[2] & 0xff);
    value <<= 8;
    value |= (buffer[3] & 0xff);

    values[i] = value;
  }
#else
  for (i = 0; i < n; i++) {
    read(positions_fd,buffer,4);

    value = (buffer[3] & 0xff);
    value <<= 8;
    value |= (buffer[2] & 0xff);
    value <<= 8;
    value |= (buffer[1] & 0xff);
    value <<= 8;
    value |= (buffer[0] & 0xff);

    values[i] = value;
  }
#endif

  return;
}
#endif



#ifdef LARGE_GENOMES
static UINT4 *
point_one_shift (int *nentries, unsigned char **positions_high, T this, Storedoligomer_T subst) {
  UINT4 *positions_low;
  Positionsptr_T ptr0, end0;
#ifdef DEBUG
  int i;
#endif

  if (this->compression_type == NO_COMPRESSION) {
#ifdef WORDS_BIGENDIAN
    abort();
#else
    ptr0 = this->offsetsstrm[subst];
    end0 = this->offsetsstrm[subst+1];
#endif

  } else if (this->compression_type == BITPACK64_COMPRESSION) {
    ptr0 = Bitpack64_read_two_huge(&end0,subst,this->offsetspages,this->offsetsmeta,this->offsetsstrm);
  }


  debug(printf("point_one_shift: %08X %u %u\n",subst,ptr0,end0));

  if ((*nentries = end0 - ptr0) == 0) {
    *positions_high = (unsigned char *) NULL;
    return (UINT4 *) NULL;
  } else {
    if (this->positions_access == FILEIO) {
      abort();

    } else {
      /* ALLOCATED or MMAPPED */
      *positions_high = &(this->positions_high[ptr0]);
      positions_low = &(this->positions_low[ptr0]);
    }
  }
      
  debug(
	printf("%d entries:",*nentries);
	for (i = 0; i < *nentries; i++) {
	  printf(" %u",(Univcoord_T) positions_high[i] << 32 + positions_low[i]);
	}
	printf("\n");
	);
  
  return positions_low;
}

#else

static Univcoord_T *
point_one_shift (int *nentries, T this, Storedoligomer_T subst) {
  Univcoord_T *positions;
  Positionsptr_T ptr0, end0;
#ifdef DEBUG
  int i;
#endif

  if (this->compression_type == NO_COMPRESSION) {
#ifdef WORDS_BIGENDIAN
#if 0
    if (this->offsetsstrm_access == ALLOCATED) {
      ptr0 = this->offsetsstrm[subst];
      end0 = this->offsetsstrm[subst+1];
    } else {
      ptr0 = Bigendian_convert_uint(this->offsetsstrm[subst]);
      end0 = Bigendian_convert_uint(this->offsetsstrm[subst+1]);
    }
#else
    abort();
#endif
#else
    ptr0 = this->offsetsstrm[subst];
    end0 = this->offsetsstrm[subst+1];
#endif

  } else if (this->compression_type == BITPACK64_COMPRESSION) {
    ptr0 = Bitpack64_read_two(&end0,subst,this->offsetsmeta,this->offsetsstrm);
  }


  debug(printf("point_one_shift: %08X %u %u\n",subst,ptr0,end0));

  if ((*nentries = end0 - ptr0) == 0) {
    return (Univcoord_T *) NULL;
  } else {
    if (this->positions_access == FILEIO) {
      positions = (Univcoord_T *) CALLOC(*nentries,sizeof(Univcoord_T));
#ifdef HAVE_PTHREAD
      pthread_mutex_lock(&this->positions_read_mutex);
#endif
      positions_move_absolute(this->positions_fd,ptr0);
      positions_read_multiple(this->positions_fd,positions,*nentries);
#ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&this->positions_read_mutex);
#endif

    } else {
      /* ALLOCATED or MMAPPED */
      positions = &(this->positions[ptr0]);
    }
  }
      
#ifdef WORDS_BIGENDIAN
  debug(
	printf("%d entries:",*nentries);
	for (i = 0; i < *nentries; i++) {
	  printf(" %u",Bigendian_convert_univcoord(positions[i]));
	}
	printf("\n");
	);
#else
  debug(
	printf("%d entries:",*nentries);
	for (i = 0; i < *nentries; i++) {
	  printf(" %u",positions[i]);
	}
	printf("\n");
	);
#endif
  
  return positions;
}

#endif



/*                87654321 */
#define RIGHT_A 0x00000000
#define RIGHT_C 0x00000001
#define RIGHT_G 0x00000002
#define RIGHT_T 0x00000003

/*                      87654321 */
#define LOW_TWO_BITS  0x00000003

#ifdef DEBUG
static char *
shortoligo_nt (Storedoligomer_T oligo, int oligosize) {
  char *nt;
  int i, j;
  Storedoligomer_T lowbits;

  nt = (char *) CALLOC(oligosize+1,sizeof(char));
  j = oligosize-1;
  for (i = 0; i < oligosize; i++) {
    lowbits = oligo & LOW_TWO_BITS;
    switch (lowbits) {
    case RIGHT_A: nt[j] = 'A'; break;
    case RIGHT_C: nt[j] = 'C'; break;
    case RIGHT_G: nt[j] = 'G'; break;
    case RIGHT_T: nt[j] = 'T'; break;
    }
    oligo >>= 2;
    j--;
  }

  return nt;
}
#endif


#ifdef LARGE_GENOMES
static int
count_one_shift (T this, Storedoligomer_T subst, int nadjacent) {
  Positionsptr_T ptr0, end0;

  if (this->compression_type == NO_COMPRESSION) {
#ifdef WORDS_BIGENDIAN
#if 0
    if (this->offsetsstrm_access == ALLOCATED) {
      ptr0 = this->offsetsstrm[subst];
      end0 = this->offsetsstrm[subst+nadjacent];
    } else {
      ptr0 = Bigendian_convert_uint(this->offsetsstrm[subst]);
      end0 = Bigendian_convert_uint(this->offsetsstrm[subst+nadjacent]);
    }
#else
    abort();
#endif
#else
    ptr0 = this->offsetsstrm[subst];
    end0 = this->offsetsstrm[subst+nadjacent];
#endif

  } else if (this->compression_type == BITPACK64_COMPRESSION) {
    ptr0 = Bitpack64_read_one_huge(subst,this->offsetspages,this->offsetsmeta,this->offsetsstrm);
    end0 = Bitpack64_read_one_huge(subst+nadjacent,this->offsetspages,this->offsetsmeta,this->offsetsstrm);

  } else {
    abort();
  }

  debug(printf("count_one_shift: oligo = %06X (%s), %u - %u = %u\n",
	       subst,shortoligo_nt(subst,index1part),end0,ptr0,end0-ptr0));
  return (end0 - ptr0);

}

#else
static int
count_one_shift (T this, Storedoligomer_T subst, int nadjacent) {
  Positionsptr_T ptr0, end0;

  if (this->compression_type == NO_COMPRESSION) {
#ifdef WORDS_BIGENDIAN
#if 0
    if (this->offsetsstrm_access == ALLOCATED) {
      ptr0 = this->offsetsstrm[subst];
      end0 = this->offsetsstrm[subst+nadjacent];
    } else {
      ptr0 = Bigendian_convert_uint(this->offsetsstrm[subst]);
      end0 = Bigendian_convert_uint(this->offsetsstrm[subst+nadjacent]);
    }
#else
    abort();
#endif
#else
    ptr0 = this->offsetsstrm[subst];
    end0 = this->offsetsstrm[subst+nadjacent];
#endif

  } else if (this->compression_type == BITPACK64_COMPRESSION) {
    ptr0 = Bitpack64_read_one(subst,this->offsetsmeta,this->offsetsstrm);
    end0 = Bitpack64_read_one(subst+nadjacent,this->offsetsmeta,this->offsetsstrm);

  } else {
    abort();
  }

  debug(printf("count_one_shift: oligo = %06X (%s), %u - %u = %u\n",
	       subst,shortoligo_nt(subst,index1part),end0,ptr0,end0-ptr0));
  return (end0 - ptr0);

}

#endif


/************************************************************************
 *   Counting procedures
 ************************************************************************/

/* Don't mask out leftmost nucleotides with LOWXXMER */
int
Indexdb_count_left_subst_2 (T this, Storedoligomer_T oligo) {
  int nentries = 0;
  Storedoligomer_T base;
  int i;

  debug(printf("count_left_subst_2: oligo = %06X (%s)\n",oligo,shortoligo_nt(oligo,index1part)));

#ifdef ALLOW_DUPLICATES
  /* Right shift */
  base = (oligo >> 4);
  for (i = 0; i < 16; i++, base += left_subst) {
    nentries += count_one_shift(this,base);
  }
#else
  /* Right shift */
  base = (oligo >> 4);
  debug(printf("shift right => %06X (%s)\n",base,shortoligo_nt(base,index1part)));
  for (i = 0; i < 16; i++, base += left_subst) {
#if 0
    nentries += count_one_shift(this,base,/*nadjacent*/1);
#else
    nentries += Indexdb_count_no_subst(this,base);
#endif
  }
#endif
      
  return nentries;
}


/* Don't mask out leftmost nucleotides with LOWXXMER */
int
Indexdb_count_left_subst_1 (T this, Storedoligomer_T oligo) {
  int nentries = 0;
  Storedoligomer_T base;
  int i;

  debug(printf("count_left_subst_1: oligo = %06X (%s)\n",oligo,shortoligo_nt(oligo,index1part)));

#ifdef ALLOW_DUPLICATES
  /* Zero shift. */
  base = (oligo >> 2);
  for (i = 0; i < 4; i++, base += top_subst) {
    nentries += count_one_shift(this,base);
  }
#else
  /* Zero shift. */
  base = (oligo >> 2);
  for (i = 0; i < 4; i++, base += top_subst) {
#if 0
    nentries += count_one_shift(this,base,/*nadjacent*/1);
#else
    nentries += Indexdb_count_no_subst(this,base);
#endif
  }
#endif
      
  return nentries;
}


int
Indexdb_count_right_subst_2 (T this, Storedoligomer_T oligo) {
  int nentries;
  Storedoligomer_T base;
#ifdef ALLOW_DUPLICATES
  int i;
#endif
#ifdef DEBUG
  int i;
#endif

  debug(printf("count_right_subst_2: oligo = %06X (%s)\n",oligo,shortoligo_nt(oligo,index1part)));

#ifdef ALLOW_DUPLICATES
  /* Left shift */
  base = (oligo << 4) & kmer_mask;
  nentries = 0;
  for (i = 0; i < 16; i++, base += right_subst) {
    nentries += count_one_shift(this,base);
  }
#else
  /* Left shift */
  base = (oligo << 4) & kmer_mask;
  nentries = count_one_shift(this,base,/*nadjacent*/16);

  debug(
	printf("Details\n");
	nentries = 0;
	for (i = 0; i < 16; i++, base += right_subst) {
	  nentries += count_one_shift(this,base,/*nadjacent*/1);
	}
	);
#endif
      
  return nentries;
}


int
Indexdb_count_right_subst_1 (T this, Storedoligomer_T oligo) {
  int nentries;
  Storedoligomer_T base;
#ifdef ALLOW_DUPLICATES
  int i;
#endif
#ifdef DEBUG
  int i;
#endif

  debug(printf("count_right_subst_1: oligo = %06X (%s)\n",oligo,shortoligo_nt(oligo,index1part)));

#ifdef ALLOW_DUPLICATES
  /* Zero shift */
  base = (oligo << 2) & kmer_mask;
  nentries = 0;
  for (i = 0; i < 4; i++, base += right_subst) {
    nentries += count_one_shift(this,base);
  }
#else
  /* Zero shift */
  base = (oligo << 2) & kmer_mask;
  nentries = count_one_shift(this,base,/*nadjacent*/4);

  debug(
	printf("Details\n");
	nentries = 0;
	for (i = 0; i < 4; i++, base += right_subst) {
	  nentries += count_one_shift(this,base,/*nadjacent*/1);
	}
	);
#endif
      
  return nentries;
}


/************************************************************************/


static bool free_positions_p;	/* Needs to be true if Indexdb positions are FILEIO */

void
Compoundpos_init_positions_free (bool positions_fileio_p) {
  if (positions_fileio_p == true) {
    free_positions_p = true;
  } else {
    free_positions_p = false;
  }
  return;
}



struct Compoundpos_T {
  int n;

#ifdef LARGE_GENOMES
  unsigned char *positions_high[16];
  UINT4 *positions_low[16];
#else
  Univcoord_T *positions[16];
#endif
  int npositions[16];

  struct Batch_T batchpool[16];
  Batch_T heap[17];
  int heapsize;
  struct Batch_T sentinel_struct;
  Batch_T sentinel;

#ifdef LARGE_GENOMES
  unsigned char *positions_high_reset[16]; /* altered by find_nomiss_aux and find_onemiss_aux */
  UINT4 *positions_low_reset[16]; /* altered by find_nomiss_aux and find_onemiss_aux */
#else
  Univcoord_T *positions_reset[16]; /* altered by find_nomiss_aux and find_onemiss_aux */
#endif
  int npositions_reset[16]; /* altered by find_nomiss_aux and find_onemiss_aux */
};


void
Compoundpos_set (Compoundpos_T compoundpos) {
  int i;

  for (i = 0; i < compoundpos->n; i++) {
#ifdef LARGE_GENOMES
    compoundpos->positions_high_reset[i] = compoundpos->positions_high[i];
    compoundpos->positions_low_reset[i] = compoundpos->positions_low[i];
#else
    compoundpos->positions_reset[i] = compoundpos->positions[i];
#endif
    compoundpos->npositions_reset[i] = compoundpos->npositions[i];
  }
  return;
}

void
Compoundpos_reset (Compoundpos_T compoundpos) {
  int i;

  for (i = 0; i < compoundpos->n; i++) {
#ifdef LARGE_GENOMES
    compoundpos->positions_high[i] = compoundpos->positions_high_reset[i];
    compoundpos->positions_low[i] = compoundpos->positions_low_reset[i];
#else
    compoundpos->positions[i] = compoundpos->positions_reset[i];
#endif
    compoundpos->npositions[i] = compoundpos->npositions_reset[i];
  }
  return;
}


void
Compoundpos_print_sizes (Compoundpos_T compoundpos) {
  int i;

  for (i = 0; i < compoundpos->n; i++) {
    printf(" %d",compoundpos->npositions[i]);
  }

  return;
}


void
Compoundpos_dump (Compoundpos_T compoundpos, int diagterm) {
  int i, j;

  printf("%d diagonals: ",compoundpos->n);
  for (i = 0; i < compoundpos->n; i++) {
    printf(" %d",compoundpos->npositions[i]);
  }
  printf("\n");

  for (i = 0; i < compoundpos->n; i++) {
    for (j = 0; j < compoundpos->npositions[i]; j++) {
#ifdef LARGE_GENOMES
      printf(" compound%d.%d:%llu+%d\n",
	     i,j,((Univcoord_T) compoundpos->positions_high[i][j] << 32) + compoundpos->positions_low[i][j],diagterm);
#elif defined(WORDS_BIGENDIAN)
      printf(" compound%d.%d:%u+%d\n",
	     i,j,Bigendian_convert_univcoord(compoundpos->positions[i][j]),diagterm);
#else
      printf(" compound%d.%d:%u+%d\n",i,j,compoundpos->positions[i][j],diagterm);
#endif
    }
  }
  return;
}


void
Compoundpos_free (Compoundpos_T *old) {
  int i;

  if (*old) {
    if (free_positions_p == true) {
      for (i = 0; i < (*old)->n; i++) {
#ifdef LARGE_GENOMES
	FREE((*old)->positions_high[i]);
	FREE((*old)->positions_low[i]);
#else
	FREE((*old)->positions[i]);
#endif
      }
    }

    /* No need, since allocated statically.  FREE((*old)->npositions); */
    /* No need, since allocated statically.  FREE((*old)->positions); */
  
    FREE(*old);
  }
  return;
}


Compoundpos_T
Indexdb_compoundpos_left_subst_2 (T this, Storedoligomer_T oligo) {
  Compoundpos_T compoundpos = (Compoundpos_T) MALLOC(sizeof(*compoundpos));
  Storedoligomer_T base;
  int i;

  debug(printf("compoundpos_left_subst_2: %06X (%s)\n",oligo,shortoligo_nt(oligo,index1part)));

  compoundpos->n = 16;
  /* compoundpos->npositions = (int *) CALLOC(16,sizeof(int)); */
  /* compoundpos->positions = (Univcoord_T **) CALLOC(16,sizeof(Univcoord_T *)); */

  /* Right shift */
  base = (oligo >> 4);
  for (i = 0; i < 16; i++, base += left_subst) {
#ifdef LARGE_GENOMES
    compoundpos->positions_low[i] =
      point_one_shift(&(compoundpos->npositions[i]),&(compoundpos->positions_high[i]),this,base);
#else
    compoundpos->positions[i] = point_one_shift(&(compoundpos->npositions[i]),this,base);
#endif
  }

  return compoundpos;
}

Compoundpos_T
Indexdb_compoundpos_left_subst_1 (T this, Storedoligomer_T oligo) {
  Compoundpos_T compoundpos = (Compoundpos_T) MALLOC(sizeof(*compoundpos));
  Storedoligomer_T base;
  int i;

  debug(printf("compoundpos_left_subst_1: %06X (%s)\n",oligo,shortoligo_nt(oligo,index1part)));

  compoundpos->n = 4;
  /* compoundpos->npositions = (int *) CALLOC(4,sizeof(int)); */
  /* compoundpos->positions = (Univcoord_T **) CALLOC(4,sizeof(Univcoord_T *)); */

  /* Zero shift */
  base = (oligo >> 2);
  for (i = 0; i < 4; i++, base += top_subst) {
#ifdef LARGE_GENOMES
    compoundpos->positions_low[i] =
      point_one_shift(&(compoundpos->npositions[i]),&(compoundpos->positions_high[i]),this,base);
#else
    compoundpos->positions[i] = point_one_shift(&(compoundpos->npositions[i]),this,base);
#endif
  }

  return compoundpos;
}

Compoundpos_T
Indexdb_compoundpos_right_subst_2 (T this, Storedoligomer_T oligo) {
  Compoundpos_T compoundpos = (Compoundpos_T) MALLOC(sizeof(*compoundpos));
  Storedoligomer_T base;
  int i;

  debug(printf("compoundpos_right_subst_2: %06X (%s)\n",oligo,shortoligo_nt(oligo,index1part)));

  compoundpos->n = 16;
  /* compoundpos->npositions = (int *) CALLOC(16,sizeof(int)); */
  /* compoundpos->positions = (Univcoord_T **) CALLOC(16,sizeof(Univcoord_T *)); */

  /* Left shift */
  base = (oligo << 4) & kmer_mask;
  for (i = 0; i < 16; i++, base += right_subst) {
#ifdef LARGE_GENOMES
    compoundpos->positions_low[i] =
      point_one_shift(&(compoundpos->npositions[i]),&(compoundpos->positions_high[i]),this,base);
#else
    compoundpos->positions[i] = point_one_shift(&(compoundpos->npositions[i]),this,base);
#endif
  }

  return compoundpos;
}

Compoundpos_T
Indexdb_compoundpos_right_subst_1 (T this, Storedoligomer_T oligo) {
  Compoundpos_T compoundpos = (Compoundpos_T) MALLOC(sizeof(*compoundpos));
  Storedoligomer_T base;
  int i;

  debug(printf("compoundpos_right_subst_1: %06X (%s)\n",oligo,shortoligo_nt(oligo,index1part)));

  compoundpos->n = 4;
  /* compoundpos->npositions = (int *) CALLOC(4,sizeof(int)); */
  /* compoundpos->positions = (Univcoord_T **) CALLOC(4,sizeof(Univcoord_T *)); */

  /* Zero shift */
  base = (oligo << 2) & kmer_mask;
  for (i = 0; i < 4; i++, base += right_subst) {
#ifdef LARGE_GENOMES
    compoundpos->positions_low[i] =
      point_one_shift(&(compoundpos->npositions[i]),&(compoundpos->positions_high[i]),this,base);
#else
    compoundpos->positions[i] = point_one_shift(&(compoundpos->npositions[i]),this,base);
#endif
  }

  return compoundpos;
}



/************************************************************************/

#ifdef LARGE_GENOMES
static int
binary_search (int lowi, int highi, unsigned char *positions_high, UINT4 *positions_low, Univcoord_T goal) {
  bool foundp = false;
  int middlei;
  Univcoord_T position;

#ifdef NOBINARY
  return lowi;
#endif

  if (goal == 0U) {
    return lowi;
  }

  while (!foundp && lowi < highi) {
    middlei = lowi + ((highi - lowi) / 2);
    position = ((Univcoord_T) positions_high[middlei] << 32) + positions_low[middlei];
    debug2(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		  lowi,(positions_high[lowi] << 32) + positions_low[lowi],
		  middlei,position,
		  highi,(positions_high[highi] << 32) + positions_low[highi],goal));
    if (goal < position) {
      highi = middlei;
    } else if (goal > position) {
      lowi = middlei + 1;
    } else {
      foundp = true;
    }
  }

  if (foundp == true) {
    return middlei;
  } else {
    return highi;
  }
}

#else

static int
binary_search (int lowi, int highi, Univcoord_T *positions, Univcoord_T goal) {
  bool foundp = false;
  int middlei;

#ifdef NOBINARY
  return lowi;
#endif

  if (goal == 0U) {
    return lowi;
  }

  while (!foundp && lowi < highi) {
    middlei = lowi + ((highi - lowi) / 2);
#ifdef WORDS_BIGENDIAN
    debug2(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		  lowi,Bigendian_convert_univcoord(positions[lowi]),
		  middlei,Bigendian_convert_univcoord(positions[middlei]),
		  highi,Bigendian_convert_univcoord(positions[highi]),goal));
    if (goal < Bigendian_convert_univcoord(positions[middlei])) {
      highi = middlei;
    } else if (goal > Bigendian_convert_univcoord(positions[middlei])) {
      lowi = middlei + 1;
    } else {
      foundp = true;
    }
#else
    debug2(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		  lowi,positions[lowi],middlei,positions[middlei],
		  highi,positions[highi],goal));
    if (goal < positions[middlei]) {
      highi = middlei;
    } else if (goal > positions[middlei]) {
      lowi = middlei + 1;
    } else {
      foundp = true;
    }
#endif
  }

  if (foundp == true) {
    return middlei;
  } else {
    return highi;
  }
}

#endif


void
Compoundpos_heap_init (Compoundpos_T compoundpos, int querylength, int diagterm) {
  Batch_T batch;
  int startbound, i;

  compoundpos->heapsize = 0;
  for (i = 0; i < compoundpos->n; i++) {
    batch = &(compoundpos->batchpool[i]);
#ifdef LARGE_GENOMES
    batch->positionptr_high = compoundpos->positions_high[i];
    batch->positionptr_low = compoundpos->positions_low[i];
#else
    batch->positionptr = compoundpos->positions[i];
#endif
    batch->nentries = compoundpos->npositions[i];
    if (diagterm < querylength) {
      startbound = querylength - diagterm;
#ifdef LARGE_GENOMES
      while (batch->nentries > 0 && (((Univcoord_T) *batch->positionptr_high) << 32) + (*batch->positionptr_low) < (unsigned int) startbound) {
	debug11(printf("Eliminating diagonal %u as straddling beginning of genome (Compoundpos_heap_init)\n",
		       ((Univcoord_T) *batch->positionptr_high << 32) + *batch->positionptr_low));
	++batch->positionptr_high;
	++batch->positionptr_low;
	--batch->nentries;
      }
#elif defined(WORDS_BIGENDIAN)
      while (batch->nentries > 0 && Bigendian_convert_univcoord(*batch->positionptr) < (unsigned int) startbound) {
	debug11(printf("Eliminating diagonal %u as straddling beginning of genome (Compoundpos_heap_init)\n",
		       Bigendian_convert_univcoord(*batch->positionptr)));
	++batch->positionptr;
	--batch->nentries;
      }
#else
      while (batch->nentries > 0 && *batch->positionptr < (unsigned int) startbound) {
	debug11(printf("Eliminating diagonal %u as straddling beginning of genome (Compoundpos_heap_init)\n",
		       *batch->positionptr));
	++batch->positionptr;
	--batch->nentries;
      }
#endif
    }
    if (batch->nentries > 0) {
#ifdef LARGE_GENOMES
      batch->position = (((Univcoord_T) *batch->positionptr_high) << 32) + (*batch->positionptr_low);
#elif defined(WORDS_BIGENDIAN)
      batch->position = Bigendian_convert_univcoord(*batch->positionptr);
#else
      batch->position = *batch->positionptr;
#endif
      heap_insert_even(compoundpos->heap,&compoundpos->heapsize,batch,batch->position);
    }
  }

  compoundpos->sentinel_struct.position = (Univcoord_T) -1; /* infinity */
#ifdef LARGE_GENOMES
  compoundpos->sentinel_struct.positionptr_high = &sentinel_position_high;
  compoundpos->sentinel_struct.positionptr_low = &sentinel_position_low;
#else
  compoundpos->sentinel_struct.positionptr = &(compoundpos->sentinel_struct.position);
#endif
  compoundpos->sentinel = &compoundpos->sentinel_struct;

  for (i = compoundpos->heapsize+1; i <= compoundpos->n; i++) {
    compoundpos->heap[i] = compoundpos->sentinel;
  }

  return;
}


/* Used by DEBUG3 and DEBUG6 */
static void
heap_even_dump (Batch_T *heap, int heapsize) {
  int i;
  Batch_T batch;

  for (i = 1; i <= heapsize; i++) {
    batch = heap[i];
    printf("#%d--%d:%u  ",i,batch->nentries,batch->position);
  }
  printf("\n");
}



/* Returns true if found.  emptyp is true only if every batch is
   empty.  If procedure returns true, empty is guaranteed to be
   false. */
bool
Compoundpos_find (bool *emptyp, Compoundpos_T compoundpos, Univcoord_T local_goal) {
  Batch_T *heap = compoundpos->heap, batch;
  int i, j;

  debug6(printf("\nEntering Compoundpos_find with local_goal %u\n",local_goal));

  *emptyp = true;
  i = 1;
  while (i <= compoundpos->heapsize) {
    debug6(printf("Compoundpos_find iteration, heapsize %d:\n",compoundpos->heapsize));
    debug6(heap_even_dump(heap,compoundpos->heapsize));

    batch = heap[i];
#ifdef LARGE_GENOMES
    if (batch->nentries > 0 && (((Univcoord_T) *batch->positionptr_high) << 32) + (*batch->positionptr_low) < local_goal) {
      j = 1;
      while (j < batch->nentries &&
	     ((Univcoord_T) batch->positionptr_high[j] << 32) + batch->positionptr_low[j] < local_goal) {
	j <<= 1;		/* gallop by 2 */
      }
      if (j >= batch->nentries) {
	j = binary_search(j >> 1,batch->nentries,batch->positionptr_high,batch->positionptr_low,local_goal);
      } else {
	j = binary_search(j >> 1,j,batch->positionptr_high,batch->positionptr_low,local_goal);
      }
      batch->positionptr_high += j;
      batch->positionptr_low += j;
      batch->nentries -= j;
      debug6(printf("binary search jump %d positions to %d:%u\n",
		    j,batch->nentries,(((Univcoord_T) *batch->positionptr_high) << 32) + (*batch->positionptr_low)));
    }
#elif defined(WORDS_BIGENDIAN)
    if (batch->nentries > 0 && Bigendian_convert_univcoord(*batch->positionptr) < local_goal) {
      j = 1;
      while (j < batch->nentries && Bigendian_convert_univcoord(batch->positionptr[j]) < local_goal) {
	j <<= 1;		/* gallop by 2 */
      }
      if (j >= batch->nentries) {
	j = binary_search(j >> 1,batch->nentries,batch->positionptr,local_goal);
      } else {
	j = binary_search(j >> 1,j,batch->positionptr,local_goal);
      }
      batch->positionptr += j;
      batch->nentries -= j;
      debug6(printf("binary search jump %d positions to %d:%u\n",
		    j,batch->nentries,Bigendian_convert_univcoord(*batch->positionptr)));
    }
#else
    if (batch->nentries > 0 && *batch->positionptr < local_goal) {
      j = 1;
      while (j < batch->nentries && batch->positionptr[j] < local_goal) {
	j <<= 1;		/* gallop by 2 */
      }
      if (j >= batch->nentries) {
	j = binary_search(j >> 1,batch->nentries,batch->positionptr,local_goal);
      } else {
	j = binary_search(j >> 1,j,batch->positionptr,local_goal);
      }
      batch->positionptr += j;
      batch->nentries -= j;
      debug6(printf("binary search jump %d positions to %d:%u\n",
		    j,batch->nentries,*batch->positionptr));
    }
#endif

    if (batch->nentries <= 0) {
      /* Empty, so continue with loop */
      /* Move last heap to this one, and reduce heapsize */
      compoundpos->heap[i] = compoundpos->heap[compoundpos->heapsize];
      --compoundpos->heapsize;

#ifdef LARGE_GENOMES
    } else if (((Univcoord_T) *batch->positionptr_high << 32) + (*batch->positionptr_low) > local_goal) {
      /* Already advanced past goal, so continue with loop */
      debug6(printf("Setting emptyp to be false\n"));
      *emptyp = false;
      i++;
#elif defined(WORDS_BIGENDIAN)
    } else if (Bigendian_convert_univcoord(*batch->positionptr) > local_goal) {
      /* Already advanced past goal, so continue with loop */
      debug6(printf("Setting emptyp to be false\n"));
      *emptyp = false;
      i++;
#else
    } else if (*batch->positionptr > local_goal) {
      /* Already advanced past goal, so continue with loop */
      debug6(printf("Setting emptyp to be false\n"));
      *emptyp = false;
      i++;
#endif
    } else {
      /* Found goal, so return */
      debug6(printf("Setting emptyp to be false\n"));
      *emptyp = false;
#ifdef LARGE_GENOMES
      debug6(printf("Found! Returning position %llu\n",(((Univcoord_T) *batch->positionptr_high) << 32) + (*batch->positionptr_low)));
#elif defined(WORDS_BIGENDIAN)
      debug6(printf("Found! Returning position %u\n",Bigendian_convert_univcoord(*batch->positionptr)));
#else
      debug6(printf("Found! Returning position %u\n",*batch->positionptr));
#endif
#ifdef LARGE_GENOMES
      ++batch->positionptr_high;
      ++batch->positionptr_low;
#else
      ++batch->positionptr;
#endif
      --batch->nentries;
      return true;
    }
  }

  /* Done with loop: Fail. */
  debug6(printf("Returning emptyp %d\n",*emptyp));
  return false;
}



/* Returns 0 if heapsize is 0, else 1, and returns smallest value >= local_goal */
int
Compoundpos_search (Univcoord_T *value, Compoundpos_T compoundpos, Univcoord_T local_goal) {
  int parenti, smallesti, j;
  Batch_T batch, *heap = compoundpos->heap;
  Univcoord_T position;

  debug3(printf("\nEntering Compoundpos_search with local_goal %u\n",local_goal));
  if (compoundpos->heapsize <= 0) {
    debug3(printf("Returning because heapsize is %d\n",compoundpos->heapsize));
    return 0;
  }

  if (compoundpos->n == 4) {
    while (compoundpos->heapsize > 0 && (batch = heap[1])->position < local_goal) {
      debug3(printf("Compoundpos_search iteration, heapsize %d:\n",compoundpos->heapsize));
      debug3(heap_even_dump(heap,compoundpos->heapsize));
#ifdef LARGE_GENOMES
      if (batch->nentries > 0 && (((Univcoord_T) *batch->positionptr_high) << 32) + (*batch->positionptr_low) < local_goal) {
	j = 1;
	while (j < batch->nentries &&
	       ((Univcoord_T) batch->positionptr_high[j] << 32) + batch->positionptr_low[j] < local_goal) {
	  j <<= 1;		/* gallop by 2 */
	}
	if (j >= batch->nentries) {
	  j = binary_search(j >> 1,batch->nentries,batch->positionptr_high,batch->positionptr_low,local_goal);
	} else {
	  j = binary_search(j >> 1,j,batch->positionptr_high,batch->positionptr_low,local_goal);
	}
	batch->positionptr_high += j;
	batch->positionptr_low += j;
	batch->nentries -= j;
	debug3(printf("binary search jump %d positions to %d:%u\n",
		      j,batch->nentries,(((Univcoord_T) *batch->positionptr_high) << 32) + (*batch->positionptr_low)));
      }
      batch->position = (((Univcoord_T) *batch->positionptr_high) << 32) + (*batch->positionptr_low);
#elif defined(WORDS_BIGENDIAN)
      if (batch->nentries > 0 && Bigendian_convert_univcoord(*batch->positionptr) < local_goal) {
	j = 1;
	while (j < batch->nentries && Bigendian_convert_univcoord(batch->positionptr[j]) < local_goal) {
	  j <<= 1;		/* gallop by 2 */
	}
	if (j >= batch->nentries) {
	  j = binary_search(j >> 1,batch->nentries,batch->positionptr,local_goal);
	} else {
	  j = binary_search(j >> 1,j,batch->positionptr,local_goal);
	}
	batch->positionptr += j;
	batch->nentries -= j;
	debug3(printf("binary search jump %d positions to %d:%u\n",
		      j,batch->nentries,Bigendian_convert_univcoord(*batch->positionptr)));
      }
      batch->position = Bigendian_convert_univcoord(*batch->positionptr);
#else
      if (batch->nentries > 0 && *batch->positionptr < local_goal) {
	j = 1;
	while (j < batch->nentries && batch->positionptr[j] < local_goal) {
	  j <<= 1;		/* gallop by 2 */
	}
	if (j >= batch->nentries) {
	  j = binary_search(j >> 1,batch->nentries,batch->positionptr,local_goal);
	} else {
	  j = binary_search(j >> 1,j,batch->positionptr,local_goal);
	}
	batch->positionptr += j;
	batch->nentries -= j;
	debug3(printf("binary search jump %d positions to %d:%u\n",
		      j,batch->nentries,*batch->positionptr));
      }
      batch->position = *batch->positionptr;
#endif

      if (batch->nentries <= 0) {
	debug3(printf("top of heap found to be empty\n"));
	heap[1] = batch = (compoundpos->heapsize == 1) ? 
	  compoundpos->sentinel : heap[compoundpos->heapsize];
	heap[compoundpos->heapsize--] = compoundpos->sentinel;
      }
      
      position = batch->position;
      debug3(printf("heapify downward on %u\n",position));
      debug3(printf("Comparing right %d: %u\n",2,heap[2]->position));
      if (position <= heap[2]->position) {
	debug3(printf("Inserting at 1\n"));
	/* heap[1] = batch; -- not necessary because batch is already at heap[1] */
      } else {
	heap[1] = heap[2];
	debug3(printf("Comparing left %d/right %d: %u and %u\n",
		      3,4,heap[3]->position,heap[4]->position));
	smallesti = 4 - (heap[3]->position < heap[4]->position);
	if (position <= heap[smallesti]->position) {
	  debug3(printf("Inserting at 2\n"));
	  heap[2] = batch;
	} else {
	  debug3(printf("Inserting at %d\n",smallesti));
	  heap[2] = heap[smallesti];
	  heap[smallesti] = batch;
	}
      }
    }
    if (batch->position == local_goal) {
      *value = batch->position;
      debug3(printf("Found! Returning position %llu\n",(unsigned long long) *value));
      return 1;
    }

  } else {
    /* 16 batches */
    while (compoundpos->heapsize > 0 && (batch = heap[1])->position < local_goal) {
      debug3(printf("Compoundpos_search iteration, heapsize %d:\n",compoundpos->heapsize));
      debug3(heap_even_dump(heap,compoundpos->heapsize));
#ifdef LARGE_GENOMES
      if (batch->nentries > 0 && (((Univcoord_T) *batch->positionptr_high) << 32) + (*batch->positionptr_low) < local_goal) {
	j = 1;
	while (j < batch->nentries &&
	       ((Univcoord_T) batch->positionptr_high[j] << 32) + batch->positionptr_low[j] < local_goal) {
	  j <<= 1;		/* gallop by 2 */
	}
	if (j >= batch->nentries) {
	  j = binary_search(j >> 1,batch->nentries,batch->positionptr_high,batch->positionptr_low,local_goal);
	} else {
	  j = binary_search(j >> 1,j,batch->positionptr_high,batch->positionptr_low,local_goal);
	}
	batch->positionptr_high += j;
	batch->positionptr_low += j;
	batch->nentries -= j;
	debug3(printf("binary search jump %d positions to %d:%u\n",
		      j,batch->nentries,(((Univcoord_T) *batch->positionptr_high) << 32 + (*batch->positionptr_low))));
      }
      batch->position = (((Univcoord_T) *batch->positionptr_high) << 32) + (*batch->positionptr_low);
#elif defined(WORDS_BIGENDIAN)
      if (batch->nentries > 0 && Bigendian_convert_univcoord(*batch->positionptr) < local_goal) {
	j = 1;
	while (j < batch->nentries && Bigendian_convert_univcoord(batch->positionptr[j]) < local_goal) {
	  j <<= 1;		/* gallop by 2 */
	}
	if (j >= batch->nentries) {
	  j = binary_search(j >> 1,batch->nentries,batch->positionptr,local_goal);
	} else {
	  j = binary_search(j >> 1,j,batch->positionptr,local_goal);
	}
	batch->positionptr += j;
	batch->nentries -= j;
	debug3(printf("binary search jump %d positions to %d:%u\n",
		      j,batch->nentries,Bigendian_convert_univcoord(*batch->positionptr)));
      }
      batch->position = Bigendian_convert_univcoord(*batch->positionptr);
#else
      if (batch->nentries > 0 && *batch->positionptr < local_goal) {
	j = 1;
	while (j < batch->nentries && batch->positionptr[j] < local_goal) {
	  j <<= 1;		/* gallop by 2 */
	}
	if (j >= batch->nentries) {
	  j = binary_search(j >> 1,batch->nentries,batch->positionptr,local_goal);
	} else {
	  j = binary_search(j >> 1,j,batch->positionptr,local_goal);
	}
	batch->positionptr += j;
	batch->nentries -= j;
	debug3(printf("binary search jump %d positions to %d:%u\n",
		      j,batch->nentries,*batch->positionptr));
      }
      batch->position = *batch->positionptr;
#endif

      if (batch->nentries <= 0) {
	debug3(printf("top of heap found to be empty\n"));
	heap[1] = batch = (compoundpos->heapsize == 1) ? 
	  compoundpos->sentinel : heap[compoundpos->heapsize];
	heap[compoundpos->heapsize--] = compoundpos->sentinel;
      }
      
      position = batch->position;
      debug3(printf("heapify downward on %u\n",position));
      /* Comparison 0/3 */
      debug3(printf("Comparing right %d: %u\n",2,heap[2]->position));
      if (position <= heap[2]->position) {
	debug3(printf("Inserting at 1\n"));
	/* heap[1] = batch; -- not necessary because batch is already at heap[1] */
      } else {
	heap[1] = heap[2];
	/* Comparison 1/3 */
	debug3(printf("Comparing left %d/right %d: %u and %u\n",
		      3,4,heap[3]->position,heap[4]->position));
	smallesti = 4 - (heap[3]->position < heap[4]->position);
	if (position <= heap[smallesti]->position) {
	  debug3(printf("Inserting at 2\n"));
	  heap[2] = batch;
	} else {
	  heap[2] = heap[smallesti];
	  parenti = smallesti;
	  smallesti <<= 1;
	  /* Comparison 2/3 */
	  debug3(printf("Comparing left %d/right %d: %u and %u\n",
			smallesti-1,smallesti,heap[smallesti-1]->position,heap[smallesti]->position));
	  smallesti -= (heap[LEFTSIBLING2(smallesti)]->position < heap[smallesti]->position);
	  if (position <= heap[smallesti]->position) {
	    debug3(printf("Inserting at %d\n",parenti));
	    heap[parenti] = batch;
	  } else {
	    heap[parenti] = heap[smallesti];
	    parenti = smallesti;
	    smallesti <<= 1;
	    /* Comparison 3/3 */
	    debug3(printf("Comparing left %d/right %d: %u and %u\n",
			  smallesti-1,smallesti,heap[smallesti-1]->position,heap[smallesti]->position));
	    smallesti -= (heap[LEFTSIBLING2(smallesti)]->position < heap[smallesti]->position);
	    if (position <= heap[smallesti]->position) {
	      debug3(printf("Inserting at %d\n",parenti));
	      heap[parenti] = batch;
	    } else {
	      heap[parenti] = heap[smallesti];
	      debug3(printf("Inserting at %d\n",smallesti));
	      heap[smallesti] = batch;
	    }
	  }
	}
      }
    }
    if (batch->position == local_goal) {
      *value = batch->position;
      debug3(printf("Found! Returning position %llu\n",(unsigned long long) *value));
      return 1;
    }
  }

  *value = batch->position;
  debug3(printf("Returning position %llu\n",(unsigned long long) *value));
  return 1;
}


Univcoord_T *
Indexdb_merge_compoundpos (int *nmerged, Compoundpos_T compoundpos, int diagterm) {
  int i;
  Batch_T batch;
  struct Batch_T batchpool[16];
  int nentries = 0;

  debug(printf("merge_compoundpos, sizes:"));

  for (i = 0; i < compoundpos->n; i++) {
    batch = &(batchpool[i]);
#ifdef LARGE_GENOMES
    batch->positionptr_high = compoundpos->positions_high[i];
    batch->positionptr_low = compoundpos->positions_low[i];
#else
    batch->positionptr = compoundpos->positions[i];
#endif
    batch->nentries = compoundpos->npositions[i];
    debug(printf(" %d",batch->nentries));
    nentries += batch->nentries;
  }
  debug(printf("\n"));

  if (nentries == 0) {
    *nmerged = 0;
    return (Univcoord_T *) NULL;
  } else if (compoundpos->n == 4) {
    return merge_batches_one_heap_4_existing(&(*nmerged),batchpool,nentries,diagterm);
  } else {
    return merge_batches_one_heap_16_existing(&(*nmerged),batchpool,nentries,diagterm);
  }
}



/* Should be the same as count_one_shift(this,oligo,1) */
int
Indexdb_count_no_subst (T this, Storedoligomer_T oligo) {
  Positionsptr_T ptr0, end0;

  if (this->compression_type == NO_COMPRESSION) {
#ifdef WORDS_BIGENDIAN
#if 0
    if (this->offsetsstrm_access == ALLOCATED) {
      ptr0 = this->offsetsstrm[oligo];
      end0 = this->offsetsstrm[oligo+1];
    } else {
      ptr0 = Bigendian_convert_uint(this->offsetsstrm[oligo]);
      end0 = Bigendian_convert_uint(this->offsetsstrm[oligo+1]);
    }
#else
    abort();
#endif
#else
    ptr0 = this->offsetsstrm[oligo];
    end0 = this->offsetsstrm[oligo+1];
#endif

  } else if (this->compression_type == BITPACK64_COMPRESSION) {
#ifdef LARGE_GENOMES
    ptr0 = Bitpack64_read_two_huge(&end0,oligo,this->offsetspages,this->offsetsmeta,this->offsetsstrm);
#else
    ptr0 = Bitpack64_read_two(&end0,oligo,this->offsetsmeta,this->offsetsstrm);
#endif

  } else {
    abort();
  }


  debug(printf("count_one_shift: oligo = %06X (%s), %u - %u = %u\n",
	       oligo,shortoligo_nt(oligo,index1part),end0,ptr0,end0-ptr0));
  return (end0 - ptr0);
}


#if 0
int
Indexdb_gsnapbase (T this) {
  if (this->index1interval == 1) {
    return index1part;
  } else if (this->index1interval == 3) {
    return index1part - 2;
  } else {
    abort();
  }
}
#endif
