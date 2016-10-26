static char rcsid[] = "$Id: iit-read-univ.c 168395 2015-06-26 17:13:13Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "iit-read-univ.h"

#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
#else
#include "littleendian.h"
#endif

#include <stdlib.h>		/* For qsort */
#include <string.h>		/* For memset */
#include <strings.h>
#include <ctype.h>		/* For isspace */
#ifdef HAVE_UNISTD_H
#include <unistd.h>		/* For mmap on Linux */
#endif
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>		/* For open, fstat, and mmap */
#endif
/* Not sure why this was included
#include <sys/param.h>
*/
#ifdef HAVE_FCNTL_H
#include <fcntl.h>		/* For open */
#endif
#ifdef HAVE_SYS_STAT_H
#include <sys/stat.h>		/* For open and fstat */
#endif
#include <sys/mman.h>		/* For mmap and madvise */
#include <math.h>		/* For qsort */
#include <errno.h>		/* For perror */
#include "assert.h"
#include "except.h"
#include "mem.h"
#include "access.h"
#include "fopen.h"
#include "uintlist.h"
#include "intlist.h"


/* Note: if sizeof(int) or sizeof(unsigned int) are not 4, then the below code is faulty */


/* Integer interval tree. */

/*
 * n intervals;
 *   specified by their indices e[1..n]
 *   and endpoint-access function:
 *                low  (e[i])
 *                high (e[i])
 *        is_contained (x, e[i])
 *   eg:
 *        interval e[i]          ... "[" low (e[i]) "," high (e[i]) ")"
 *        is_contained (x, e[i]) ... (    (low (e[i]) <= x
 *                                    and (x < high (e[i]))
 */

/*--------------------------------------------------------------------------*/

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Timing */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* Flanking */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Binary search */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif


/* Cannot use version in iitdef.h, because that uses Chrpos_T for the
   value, and utility functions need to have Univcoord_T (or UINT8 if
   available). */
typedef struct Univ_FNode_T *Univ_FNode_T;
struct Univ_FNode_T {
  Univ_IIT_coord_T value;
  int a;
  int b;
  int leftindex;
  int rightindex;
};

#define T Univ_IIT_T
struct T {
  bool coord_values_8p;

  int fd;
  Access_T access;		/* access type */

#ifdef HAVE_PTHREAD
  pthread_mutex_t read_mutex;
#endif

  int ntypes;			/* Always >= 1 */

  int total_nintervals;
  int nnodes;

  int *sigmas;			/* Ordering for IIT */
  int *omegas;			/* Ordering for IIT */

  struct Univ_FNode_T *nodes;
  struct Univinterval_T *intervals;

  UINT4 *typepointers;
  char *typestrings;

  off_t labelorder_offset;
  size_t labelorder_length; /* mmap length (mmap uses size_t, not off_t) */
  char *labelorder_mmap;

  off_t labelpointers_offset;
  size_t labelpointers_length; /* mmap length (mmap uses size_t, not off_t) */
  char *labelpointers_mmap;

  off_t label_offset;
  size_t label_length;		/* mmap length (mmap uses size_t, not off_t) */
  char *label_mmap;

  off_t annotpointers_offset;
  size_t annotpointers_length; /* mmap length (mmap uses size_t, not off_t) */
  char *annotpointers_mmap;

  off_t annot_offset;
  size_t annot_length;		/* mmap length (mmap uses size_t, not off_t) */
  char *annot_mmap;

  int *labelorder;
  UINT4 *labelpointers;
  char *labels;

  UINT4 *annotpointers;
  char *annotations;

  void **datapointers;
};



static void
file_move_absolute (int fd, off_t offset, off_t objsize, Univcoord_T n) {
  off_t position = offset + n*objsize;

  if (lseek(fd,position,SEEK_SET) < 0) {
    perror("Error in gmap, file_move_label");
    exit(9);
  }
  return;
}


bool
Univ_IIT_coord_values_8p (T this) {
  return this->coord_values_8p;
}

int
Univ_IIT_total_nintervals (T this) {
  return this->total_nintervals;
}

int
Univ_IIT_ntypes (T this) {
  return this->ntypes;
}

Univcoord_T
Univ_IIT_length (T this, int index) {
  Univinterval_T interval;

  interval = &(this->intervals[index-1]);
  return Univinterval_length(interval);
}


/* Assumes intervals are stored using universal coordinates */
Univcoord_T
Univ_IIT_genomelength (T chromosome_iit, bool with_circular_alias_p) {
  Univcoord_T max = 0U, this_max;
  Univinterval_T interval;
  int i;
  int circular_typeint;

  circular_typeint = Univ_IIT_typeint(chromosome_iit,"circular");

  for (i = 0; i < chromosome_iit->total_nintervals; i++) {
    interval = &(chromosome_iit->intervals[i]);
    if (with_circular_alias_p == true && Univinterval_type(interval) == circular_typeint) {
      this_max = Univinterval_high(interval) + Univinterval_length(interval);
    } else {
      this_max = Univinterval_high(interval);
    }
    if (this_max > max) {
      max = this_max;
    }
  }

  /* Convert from zero-based coordinate */
  return max+1U;
}


bool *
Univ_IIT_circularp (bool *any_circular_p, T chromosome_iit) {
  bool *circularp;
  Univinterval_T interval;
  int chrnum, nchromosomes;
  int circular_typeint;

  nchromosomes = chromosome_iit->total_nintervals;
  circularp = (bool *) CALLOC(nchromosomes+1,sizeof(bool));

  *any_circular_p = false;
  circularp[0] = false;		/* chrnum of 0 indicates translocation */
  if ((circular_typeint = Univ_IIT_typeint(chromosome_iit,"circular")) >= 0) {
    for (chrnum = 0; chrnum < nchromosomes; chrnum++) {
      interval = &(chromosome_iit->intervals[chrnum]);
      if (Univinterval_type(interval) == circular_typeint) {
	*any_circular_p = true;
	circularp[chrnum+1] = true;
      }
    }
  }

  return circularp;
}


Univinterval_T
Univ_IIT_interval (T this, int index) {
  assert(index <= this->total_nintervals);
  return &(this->intervals[index-1]); /* Convert to 0-based */
}

Univcoord_T
Univ_IIT_interval_low (T this, int index) {
  Univinterval_T interval;

  assert(index <= this->total_nintervals);
  interval = &(this->intervals[index-1]);
  return Univinterval_low(interval);
}

Univcoord_T
Univ_IIT_interval_high (T this, int index) {
  Univinterval_T interval;

  assert(index <= this->total_nintervals);
  interval = &(this->intervals[index-1]);
  return Univinterval_high(interval);
}

Univcoord_T
Univ_IIT_interval_length (T this, int index) {
  Univinterval_T interval;

  assert(index <= this->total_nintervals);
  interval = &(this->intervals[index-1]);
  return Univinterval_length(interval);
}

int
Univ_IIT_interval_type (T this, int index) {
  Univinterval_T interval;

  assert(index <= this->total_nintervals);
  interval = &(this->intervals[index-1]);
  return Univinterval_type(interval);
}


Univcoord_T
Univ_IIT_next_chrbound (T this, int index, int circular_typeint) {
  Univinterval_T interval;

  assert(index <= this->total_nintervals);
  interval = &(this->intervals[index-1]);
  if (Univinterval_type(interval) == circular_typeint) {
    return Univinterval_high(interval) + Univinterval_length(interval);
  } else {
    return Univinterval_high(interval);
  }
}

/* chrhigh is one past the highest position in the chromosome */
void
Univ_IIT_interval_bounds (Univcoord_T *low, Univcoord_T *high, Chrpos_T *length, T this,
			  int index, int circular_typeint) {
  Univinterval_T interval;

  assert(index > 0);
  assert(index <= this->total_nintervals);

  interval = &(this->intervals[index-1]);
  *low = Univinterval_low(interval);
  *length = Univinterval_length(interval);
  if (Univinterval_type(interval) == circular_typeint) {
    *high = Univinterval_high(interval) + 1 + (*length);
  } else {
    *high = Univinterval_high(interval) + 1;
  }
  return;
}

void
Univ_IIT_intervals_setup (Univcoord_T **chroffsets, Univcoord_T **chrhighs, Chrpos_T **chrlengths,
			  T this, int nchromosomes, int circular_typeint) {
  Univinterval_T interval;
  Univcoord_T *chroffsets_ptr, *chrhighs_ptr;
  Chrpos_T *chrlengths_ptr, length;
  int index;

  chroffsets_ptr = *chroffsets = (Univcoord_T *) CALLOC(nchromosomes,sizeof(Univcoord_T));
  chrhighs_ptr = *chrhighs = (Univcoord_T *) CALLOC(nchromosomes,sizeof(Univcoord_T));
  chrlengths_ptr = *chrlengths = (Chrpos_T *) CALLOC(nchromosomes,sizeof(Chrpos_T));

  for (index = 1; index <= nchromosomes; index++) {
    interval = &(this->intervals[index-1]);
    *chroffsets_ptr++ = Univinterval_low(interval);
    length = *chrlengths_ptr++ = Univinterval_length(interval);
    if (Univinterval_type(interval) == circular_typeint) {
      *chrhighs_ptr++ = Univinterval_high(interval) + 1 + length;
    } else {
      *chrhighs_ptr++ = Univinterval_high(interval) + 1;
    }
  }

  return;
}


/* Maps from chromosome_iit chrnums to iit divints */
int *
Univ_IIT_divint_crosstable (T chromosome_iit, IIT_T iit) {
  int *crosstable;
  int chrnum, nchromosomes;
  char *chr;
  UINT4 start;

  nchromosomes = chromosome_iit->total_nintervals;
  crosstable = (int *) CALLOC(nchromosomes+1,sizeof(int));

  for (chrnum = 0; chrnum < nchromosomes; chrnum++) {
#ifdef WORDS_BIGENDIAN
    /* chromosome_iit should be version 1 */
    start = Bigendian_convert_uint(chromosome_iit->labelpointers[chrnum]);
#else
    start = chromosome_iit->labelpointers[chrnum];
#endif
    chr = &(chromosome_iit->labels[start]);

    /* upon lookup, chrnum from Univ_IIT_get_one(chromosome_iit)
       is 1-based, so we need to store in chrnum+1 */
    crosstable[chrnum+1] = IIT_divint(iit,chr);
#if 0
    if (crosstable[chrnum+1] < 0) {
      fprintf(stderr,"Note: No splicesites are provided in chr %s\n",chr);
    } else {
      fprintf(stderr,"chrnum %d maps to splicesite divint %d\n",chrnum,crosstable[chrnum]);
    }
#endif
  }

  return crosstable;
}


/* The iit file has a '\0' after each string, so functions know where
   it ends */
char *
Univ_IIT_typestring (T this, int type) {
  UINT4 start;

  start = this->typepointers[type];
  return &(this->typestrings[start]);
}

int
Univ_IIT_typeint (T this, char *typestring) {
  int i = 0;
  UINT4 start;

  while (i < this->ntypes) {
    start = this->typepointers[i];
    if (!strcmp(typestring,&(this->typestrings[start]))) {
      return i;
    }
    i++;
  }

  return -1;
}

char *
Univ_IIT_label (T this, int index, bool *allocp) {
  int recno;
  UINT4 start;

  recno = index - 1; /* Convert to 0-based */

#ifdef WORDS_BIGENDIAN
  start = Bigendian_convert_uint(this->labelpointers[recno]);
#else
  start = this->labelpointers[recno];
#endif
  *allocp = false;
  return &(this->labels[start]);
}


static char EMPTY_STRING[1] = {'\0'};

/* The iit file has a '\0' after each string, so functions know where
   it ends */
/* Note: annotation itself is never allocated */
char *
Univ_IIT_annotation (char **restofheader, T this, int index, bool *alloc_header_p) {
  int recno;
  UINT4 start;
  char *annotation, *p;
  int len;

  recno = index - 1; /* Convert to 0-based */
#ifdef WORDS_BIGENDIAN
  start = Bigendian_convert_uint(this->annotpointers[recno]);
#else
  start = this->annotpointers[recno];
#endif

  annotation =  &(this->annotations[start]);
  if (annotation[0] == '\0') {
    *restofheader = annotation; /* Both are empty strings */

    *alloc_header_p = false;
    return annotation;

  } else if (annotation[0] == '\n') {
    *restofheader = EMPTY_STRING;

    *alloc_header_p = false;
    return &(annotation[1]);

  } else {
    p = annotation;
    while (*p != '\0' && *p != '\n') p++;
    len = (p - annotation)/sizeof(char);
    *restofheader = (char *) CALLOC(1+len+1,+sizeof(char));
    *restofheader[0] = ' ';
    strncpy(&((*restofheader)[1]),annotation,len);

    if (*p == '\n') p++;

    *alloc_header_p = true;
    return p;
  }
}


void
Univ_IIT_dump_typestrings (FILE *fp, T this) {
  int type;
  UINT4 start;

  for (type = 0; type < this->ntypes; type++) {
    start = this->typepointers[type];
    fprintf(fp,"%d\t%s\n",type,&(this->typestrings[start]));
  }
  return;
}

void
Univ_IIT_dump (T this) {
  int index = 0, i;
  Univinterval_T interval;
  Univcoord_T startpos, endpos;
  char *label, *annotation, *restofheader;
  bool allocp;

  for (i = 0; i < this->total_nintervals; i++) {
    interval = &(this->intervals[i]);
    label = Univ_IIT_label(this,index+1,&allocp);
    printf(">%s",label);
    if (allocp == true) {
      FREE(label);
    }
    startpos = Univinterval_low(interval);
    endpos = startpos + Univinterval_length(interval) - 1U;

    printf(" %llu..%llu",(unsigned long long) startpos,(unsigned long long) endpos);

    if (Univinterval_type(interval) > 0) {
      printf(" %s",Univ_IIT_typestring(this,Univinterval_type(interval)));
    }

    annotation = Univ_IIT_annotation(&restofheader,this,index+1,&allocp);
    printf("%s\n",restofheader);
    printf("%s",annotation);
    if (allocp) {
      FREE(restofheader);
    }

    index++;
  }

  return;
}

void
Univ_IIT_dump_table (T this, bool zerobasedp) {
  int index = 0, i;
  Univinterval_T interval;
  Univcoord_T startpos, endpos;
  char *label;
  bool allocp;

  for (i = 0; i < this->total_nintervals; i++) {
    interval = &(this->intervals[i]);
    label = Univ_IIT_label(this,index+1,&allocp);
    printf("%s\t",label);
    if (allocp == true) {
      FREE(label);
    }
    startpos = Univinterval_low(interval);
    endpos = startpos + Univinterval_length(interval) - 1U;

    if (zerobasedp) {
      printf("%llu..%llu\t",(unsigned long long) startpos,(unsigned long long) endpos);
    } else {
      printf("%llu..%llu\t",(unsigned long long) startpos+1,(unsigned long long) endpos+1);
    }

    printf("%u",Univinterval_length(interval));
    if (Univinterval_type(interval) > 0) {
      printf("\t%s",Univ_IIT_typestring(this,Univinterval_type(interval)));
    }
    printf("\n");

    index++;
  }

  return;
}

void
Univ_IIT_dump_fai (T this) {
  int index = 0, i;
  Univinterval_T interval;
  char *label;
  bool allocp;

  for (i = 0; i < this->total_nintervals; i++) {
    interval = &(this->intervals[i]);
    label = Univ_IIT_label(this,index+1,&allocp);
    printf("%s %u\n",label,Univinterval_length(interval));
    if (allocp == true) {
      FREE(label);
    }
    index++;
  }

  return;
}


#ifdef USE_MPI
/* For chromosome.iit file, which is stored in version 1 */
void
Univ_IIT_dump_sam (MPI_File fp, T this, char *sam_read_group_id, char *sam_read_group_name,
		   char *sam_read_group_library, char *sam_read_group_platform) {
  int index = 0, i;
  Univinterval_T interval;
  Chrpos_T interval_length;
  char *label, buffer[20];
  bool allocp;
  int circular_typeint;

  if (this == NULL) {
    return;
  } else {
    circular_typeint = Univ_IIT_typeint(this,"circular");
  }

  for (i = 0; i < this->total_nintervals; i++) {
    interval = &(this->intervals[i]);
    label = Univ_IIT_label(this,index+1,&allocp);
    MPI_File_write_shared(fp,"@SQ\tSN:",strlen("@SQ\tSN:"),MPI_CHAR,MPI_STATUS_IGNORE);
    MPI_File_write_shared(fp,label,strlen(label),MPI_CHAR,MPI_STATUS_IGNORE);
    if (allocp == true) {
      FREE(label);
    }
    /* startpos = Univinterval_low(interval); */
    /* endpos = startpos + Univinterval_length(interval) - 1U; */

    interval_length = Univinterval_length(interval);
    sprintf(buffer,"%u",interval_length);
    MPI_File_write_shared(fp,"\tLN:%s",strlen("\tLN:")+strlen(buffer),MPI_CHAR,MPI_STATUS_IGNORE);
    if (Univinterval_type(interval) == circular_typeint) {
      MPI_File_write_shared(fp,"\ttp:circular",strlen("\ttp:circular"),MPI_CHAR,MPI_STATUS_IGNORE);
    }
    MPI_File_write_shared(fp,"\n",1,MPI_CHAR,MPI_STATUS_IGNORE);

    index++;
  }

  if (sam_read_group_id != NULL) {
    MPI_File_write_shared(fp,"@RG\tID:",strlen("@RG\tID:"),MPI_CHAR,MPI_STATUS_IGNORE);
    MPI_File_write_shared(fp,sam_read_group_id,strlen(sam_read_group_id),MPI_CHAR,MPI_STATUS_IGNORE);

    if (sam_read_group_platform != NULL) {
      MPI_File_write_shared(fp,"\tPL:",strlen("\tPL:"),MPI_CHAR,MPI_STATUS_IGNORE);
      MPI_File_write_shared(fp,sam_read_group_platform,strlen(sam_read_group_platform),MPI_CHAR,MPI_STATUS_IGNORE);
    }
    if (sam_read_group_library != NULL) {
      MPI_File_write_shared(fp,"\tLB:",strlen("\tLB:"),MPI_CHAR,MPI_STATUS_IGNORE);
      MPI_File_write_shared(fp,sam_read_group_library,strlen(sam_read_group_library),MPI_CHAR,MPI_STATUS_IGNORE);
    }
    MPI_File_write_shared(fp,"\tSM:",strlen("\tSM:"),MPI_CHAR,MPI_STATUS_IGNORE);
    MPI_File_write_shared(fp,sam_read_group_name,strlen(sam_read_group_name),MPI_CHAR,MPI_STATUS_IGNORE);
    MPI_File_write_shared(fp,"\n",1,MPI_CHAR,MPI_STATUS_IGNORE);
  }

  return;
}


int
Univ_IIT_reserve_sam (T this, char *sam_read_group_id, char *sam_read_group_name,
		      char *sam_read_group_library, char *sam_read_group_platform) {
  int nchars = 0;
  int index = 0, i;
  Univinterval_T interval;
  Chrpos_T interval_length;
  char *label, buffer[20];
  bool allocp;
  int circular_typeint;

  if (this == NULL) {
    return 0;
  } else {
    circular_typeint = Univ_IIT_typeint(this,"circular");
  }

  for (i = 0; i < this->total_nintervals; i++) {
    interval = &(this->intervals[i]);
    label = Univ_IIT_label(this,index+1,&allocp);
    nchars += strlen("@SQ\tSN:");
    nchars += strlen(label);
    if (allocp == true) {
      FREE(label);
    }
    /* startpos = Univinterval_low(interval); */
    /* endpos = startpos + Univinterval_length(interval) - 1U; */

    interval_length = Univinterval_length(interval);
    sprintf(buffer,"%u",interval_length);
    nchars += strlen("\tLN:")+strlen(buffer);
    if (Univinterval_type(interval) == circular_typeint) {
      nchars += strlen("\ttp:circular");
    }
    nchars += strlen("\n");

    index++;
  }

  if (sam_read_group_id != NULL) {
    nchars += strlen("@RG\tID:");
    nchars += strlen(sam_read_group_id);

    if (sam_read_group_platform != NULL) {
      nchars += strlen("\tPL:");
      nchars += strlen(sam_read_group_platform);
    }
    if (sam_read_group_library != NULL) {
      nchars += strlen("\tLB:");
      nchars += strlen(sam_read_group_library);
    }
    nchars += strlen("\tSM:");
    nchars += strlen(sam_read_group_name);
    nchars += strlen("\n");
  }

  return nchars;
}


#else
/* For chromosome.iit file, which is stored in version 1 */
void
Univ_IIT_dump_sam (FILE *fp, T this, char *sam_read_group_id, char *sam_read_group_name,
		   char *sam_read_group_library, char *sam_read_group_platform) {
  int index = 0, i;
  Univinterval_T interval;
  char *label;
  bool allocp;
  int circular_typeint;

  if (this == NULL) {
    return;
  } else {
    circular_typeint = Univ_IIT_typeint(this,"circular");
  }

  for (i = 0; i < this->total_nintervals; i++) {
    interval = &(this->intervals[i]);
    label = Univ_IIT_label(this,index+1,&allocp);
    fprintf(fp,"@SQ\tSN:%s",label);
    if (allocp == true) {
      FREE(label);
    }
    /* startpos = Univinterval_low(interval); */
    /* endpos = startpos + Univinterval_length(interval) - 1U; */

    fprintf(fp,"\tLN:%u",Univinterval_length(interval));
    if (Univinterval_type(interval) == circular_typeint) {
      fprintf(fp,"\ttp:circular");
    }
    fprintf(fp,"\n");

    index++;
  }

  if (sam_read_group_id != NULL) {
    fprintf(fp,"@RG\tID:%s",sam_read_group_id);
    if (sam_read_group_platform != NULL) {
      fprintf(fp,"\tPL:%s",sam_read_group_platform);
    }
    if (sam_read_group_library != NULL) {
      fprintf(fp,"\tLB:%s",sam_read_group_library);
    }
    fprintf(fp,"\tSM:%s",sam_read_group_name);
    fprintf(fp,"\n");
  }

  return;
}
#endif



Chrpos_T *
Univ_IIT_chrlengths (T this) {
  Chrpos_T *chrlengths;
  int i;
  Univinterval_T interval;

  chrlengths = (Chrpos_T *) MALLOC(this->total_nintervals * sizeof(Chrpos_T));
  for (i = 0; i < this->total_nintervals; i++) {
    interval = &(this->intervals[i]);
    chrlengths[i] = Univinterval_length(interval);
  }

  return chrlengths;
}


void
Univ_IIT_dump_labels (FILE *fp, T this) {
  int i;
  UINT4 start;
  char *label;

  for (i = 0; i < this->total_nintervals; i++) {
#ifdef WORDS_BIGENDIAN
    start = Bigendian_convert_uint(this->labelpointers[i]);
#else
    start = this->labelpointers[i];
#endif
    label = &(this->labels[start]);
    fprintf(fp,"%s ",label);
  }
  fprintf(fp,"\n");
  return;
}



/* The iit file has a '\0' after each string, so functions know where
   it ends */
char
Univ_IIT_annotation_firstchar (T this, int index) {
  int recno;
  UINT4 start;

  recno = index - 1; /* Convert to 0-based */

#ifdef WORDS_BIGENDIAN
  start = Bigendian_convert_uint(this->annotpointers[recno]);
#else
  start = this->annotpointers[recno];
#endif

  return this->annotations[start];
}




/* For contig.iit file, which is stored in version 1 */
void
Univ_IIT_dump_contigs (T this, T chromosome_iit, bool directionalp) {
  int index = 0, i, chromosome_index;
  Univinterval_T interval;
  Univcoord_T startpos, endpos, chroffset;
  Chrpos_T chrstart, chrend;
  char *label, firstchar, *chrstring;
  bool allocp;

  for (i = 0; i < this->total_nintervals; i++) {
    interval = &(this->intervals[i]);
    label = Univ_IIT_label(this,index+1,&allocp);
    printf("%s\t",label);
    if (allocp == true) {
      FREE(label);
    }
    startpos = Univinterval_low(interval);
    endpos = startpos + Univinterval_length(interval) - 1U;

    chromosome_index = Univ_IIT_get_one(chromosome_iit,startpos,startpos);
    chroffset = Univinterval_low(Univ_IIT_interval(chromosome_iit,chromosome_index));
    chrstart = startpos - chroffset;
    chrend = endpos - chroffset;

    chrstring = Univ_IIT_label(chromosome_iit,chromosome_index,&allocp); 

    if (directionalp == false) {
      printf("%llu..%llu\t",(unsigned long long) startpos+1U,(unsigned long long) endpos+1U);
      printf("%s:%u..%u\t",chrstring,chrstart+1U,chrend+1U);

    } else {
      firstchar = Univ_IIT_annotation_firstchar(this,index+1);
      if (firstchar == '-') {
	printf("%llu..%llu\t",(unsigned long long) endpos+1U,(unsigned long long) startpos+1U);
	printf("%s:%u..%u\t",chrstring,chrend+1U,chrstart+1U);
      } else {
	printf("%llu..%llu\t",(unsigned long long) startpos+1U,(unsigned long long) endpos+1U);
	printf("%s:%u..%u\t",chrstring,chrstart+1U,chrend+1U);
      }
    }
    if (allocp == true) {
      FREE(chrstring);
    }
    
    printf("%u",Univinterval_length(interval));
    if (Univinterval_type(interval) > 0) {
      printf("\t%s",Univ_IIT_typestring(this,Univinterval_type(interval)));
    }
    printf("\n");
    
    index++;
  }

  return;
}


/************************************************************************
 * For file format, see iit-write-univ.c
 ************************************************************************/

void
Univ_IIT_free (T *old) {

  if (*old != NULL) {
    if ((*old)->access == MMAPPED) {
#ifdef HAVE_MMAP
      munmap((void *) (*old)->annot_mmap,(*old)->annot_length);
      munmap((void *) (*old)->annotpointers_mmap,(*old)->annotpointers_length);
      munmap((void *) (*old)->label_mmap,(*old)->label_length);
      munmap((void *) (*old)->labelpointers_mmap,(*old)->labelpointers_length);
      munmap((void *) (*old)->labelorder_mmap,(*old)->labelorder_length);
#endif
      close((*old)->fd);

    } else if ((*old)->access == FILEIO) {
      FREE((*old)->annotations);
      FREE((*old)->annotpointers);
      FREE((*old)->labels);
      FREE((*old)->labelpointers);
      FREE((*old)->labelorder);
      /* close((*old)->fd); -- closed in read_annotations */

    } else if ((*old)->access == ALLOCATED_PRIVATE) {
      /* Nothing to close.  IIT must have been created by Univ_IIT_new. */

    } else if ((*old)->access == ALLOCATED_SHARED) {
      /* Nothing to close.  IIT must have been created by Univ_IIT_new. */

    } else {
      abort();
    }

    FREE((*old)->typestrings);
    FREE((*old)->typepointers);

    FREE((*old)->intervals);

    /* Note: we are depending on Mem_free() to check that these are non-NULL */
    FREE((*old)->nodes);
    FREE((*old)->omegas);
    FREE((*old)->sigmas);

    FREE(*old);

  }

  return;
}



static void
move_relative (FILE *fp, off_t offset) {

#ifdef HAVE_FSEEKO
  if (fseeko(fp,offset,SEEK_CUR) < 0) {
    fprintf(stderr,"Error in move_relative, seek\n");
    abort();
  }
#else
  if (fseek(fp,(long) offset,SEEK_CUR) < 0) {
    fprintf(stderr,"Error in move_relative, seek\n");
    abort();
  }
#endif

  return;
}


static off_t
read_tree_univ (off_t offset, off_t filesize, FILE *fp, char *filename, T new) {
  size_t items_read;
  int i;
  UINT4 uint4;

  if ((offset += sizeof(int)*(new->total_nintervals+1)) > filesize) {
    fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after sigmas %lld, filesize %lld).  Did you generate it using iit_store?\n",
	    filename,(long long int) offset,(long long int) filesize);
    exit(9);
  } else {
    new->sigmas = (int *) CALLOC(new->total_nintervals+1,sizeof(int));
    if ((items_read = FREAD_INTS(new->sigmas,new->total_nintervals+1,fp)) != (unsigned int) new->total_nintervals + 1) {
      fprintf(stderr,"IIT file %s appears to be truncated\n",filename);
      exit(9);
    }
  }

  if ((offset += sizeof(int)*(new->total_nintervals+1)) > filesize) {
    fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after omegas %lld, filesize %lld).  Did you generate it using iit_store?\n",
	    filename,(long long int) offset,(long long int) filesize);
    exit(9);
  } else {
    new->omegas = (int *) CALLOC(new->total_nintervals+1,sizeof(int));
    if ((items_read = FREAD_INTS(new->omegas,new->total_nintervals+1,fp)) != (unsigned int) new->total_nintervals + 1) {
      fprintf(stderr,"IIT file %s appears to be truncated\n",filename);
      exit(9);
    }
  }

  debug(printf("nnodes: %d\n",new->nnodes));
  if (new->nnodes == 0) {
    new->nodes = (struct Univ_FNode_T *) NULL;
  } else {
    new->nodes = (struct Univ_FNode_T *) CALLOC(new->nnodes,sizeof(struct Univ_FNode_T));
#ifdef WORDS_BIGENDIAN
    if (new->coord_values_8p == true) {
#ifdef HAVE_64_BIT
      for (i = 0; i < new->nnodes; i++) {
	Bigendian_fread_uint8(&(new->nodes[i].value),fp);
	Bigendian_fread_int(&(new->nodes[i].a),fp);
	Bigendian_fread_int(&(new->nodes[i].b),fp);
	Bigendian_fread_int(&(new->nodes[i].leftindex),fp);
	Bigendian_fread_int(&(new->nodes[i].rightindex),fp);
      }
      offset += (sizeof(UINT8)+sizeof(int)+sizeof(int)+sizeof(int)+sizeof(int))*new->nnodes;
#else
      fprintf(stderr,"IIT file contains 64-bit coordinates, but this computer is only 32-bit.  Cannot continue.\n");
      exit(9);
#endif
    } else {
      for (i = 0; i < new->nnodes; i++) {
	Bigendian_fread_uint(&uint4,fp);
	new->nodes[i].value = (UINT8) uint4;
	Bigendian_fread_int(&(new->nodes[i].a),fp);
	Bigendian_fread_int(&(new->nodes[i].b),fp);
	Bigendian_fread_int(&(new->nodes[i].leftindex),fp);
	Bigendian_fread_int(&(new->nodes[i].rightindex),fp);
      }
      offset += (sizeof(UINT4)+sizeof(int)+sizeof(int)+sizeof(int)+sizeof(int))*new->nnodes;
    }
#else
    if (new->coord_values_8p == true) {
#ifdef HAVE_64_BIT
#if 1
      offset += sizeof(struct Univ_FNode_T)*fread(new->nodes,sizeof(struct Univ_FNode_T),new->nnodes,fp);
#else
      for (i = 0; i < new->nnodes; i++) {
	FREAD_UINT8(&(new->nodes[i].value),fp);
	FREAD_INT(&(new->nodes[i].a),fp);
	FREAD_INT(&(new->nodes[i].b),fp);
	FREAD_INT(&(new->nodes[i].leftindex),fp);
	FREAD_INT(&(new->nodes[i].rightindex),fp);
	printf("i %d, node value %llu\n",i,(unsigned long long) new->nodes[i].value);
      }
      offset += (sizeof(UINT8)+sizeof(int)+sizeof(int)+sizeof(int)+sizeof(int))*new->nnodes;
#endif
#else
      fprintf(stderr,"IIT file contains 64-bit coordinates, but this computer is only 32-bit.  Cannot continue.\n");
      exit(9);
#endif
    } else {
      for (i = 0; i < new->nnodes; i++) {
	FREAD_UINT(&uint4,fp);
	new->nodes[i].value = (UINT8) uint4;
	FREAD_INT(&(new->nodes[i].a),fp);
	FREAD_INT(&(new->nodes[i].b),fp);
	FREAD_INT(&(new->nodes[i].leftindex),fp);
	FREAD_INT(&(new->nodes[i].rightindex),fp);
      }
      offset += (sizeof(UINT4)+sizeof(int)+sizeof(int)+sizeof(int)+sizeof(int))*new->nnodes;
    }
#endif
    if (offset > filesize) {
      fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after nodes %lld, filesize %lld).  Did you generate it using iit_store?\n",
	      filename,(long long int) offset,(long long int) filesize);
      exit(9);
    }
  }
  debug(printf("\n"));

  return offset;
}


static off_t
read_intervals_univ (off_t offset, off_t filesize, FILE *fp, char *filename, T new) {
  int i;
  UINT4 uint4;

#ifdef WORDS_BIGENDIAN
  if (new->coord_values_8p == true) {
#ifdef HAVE_64_BIT
    for (i = 0; i < new->total_nintervals; i++) {
      Bigendian_fread_uint8(&(new->intervals[i].low),fp);
      Bigendian_fread_uint8(&(new->intervals[i].high),fp);
      Bigendian_fread_int(&(new->intervals[i].type),fp);
    }
#else
    fprintf(stderr,"IIT file contains 64-bit coordinates, but this computer is only 32-bit.  Cannot continue.\n");
    exit(9);
#endif
  } else {
    for (i = 0; i < new->total_nintervals; i++) {
      Bigendian_fread_uint(&uint4,fp);
      new->intervals[i].low = (Univcoord_T) uint4;
      Bigendian_fread_uint(&uint4,fp);
      new->intervals[i].high = (Univcoord_T) uint4;
      Bigendian_fread_int(&(new->intervals[i].type),fp);
    }
    offset += (sizeof(UINT4)+sizeof(UINT4)+sizeof(int))*new->total_nintervals;
  }
#else
  if (new->coord_values_8p == true) {
#ifdef HAVE_64_BIT
    for (i = 0; i < new->total_nintervals; i++) {
      FREAD_UINT8(&(new->intervals[i].low),fp);
      FREAD_UINT8(&(new->intervals[i].high),fp);
      FREAD_INT(&(new->intervals[i].type),fp);
    }
    offset += (sizeof(UINT8)+sizeof(UINT8)+sizeof(int))*new->total_nintervals;
#else
    fprintf(stderr,"IIT file contains 64-bit coordinates, but this computer is only 32-bit.  Cannot continue.\n");
    exit(9);
#endif
  } else {
    for (i = 0; i < new->total_nintervals; i++) {
      FREAD_UINT(&uint4,fp);
      new->intervals[i].low = (Univcoord_T) uint4;
      FREAD_UINT(&uint4,fp);
      new->intervals[i].high = (Univcoord_T) uint4;
      FREAD_INT(&(new->intervals[i].type),fp);
    }
    offset += (sizeof(UINT4)+sizeof(UINT4)+sizeof(int))*new->total_nintervals;
  }
#endif
  if (offset > filesize) {
    fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after intervals %lld, filesize %lld).  Did you generate it using iit_store?\n",
	    filename,(long long int) offset,(long long int) filesize);
    exit(9);
  }

  return offset;
}


static void
read_words (off_t offset, off_t filesize, FILE *fp, T new) {
  off_t stringlen;
  UINT4 length;
#ifdef DEBUG
  int i;
#endif

  new->typepointers = (UINT4 *) CALLOC(new->ntypes+1,sizeof(UINT4));
  offset += sizeof(int)*FREAD_UINTS(new->typepointers,new->ntypes+1,fp);
  debug(
	printf("typepointers:");
	for (i = 0; i < new->ntypes+1; i++) {
	  printf(" %u",new->typepointers[i]);
	}
	printf("\n");
	);

  stringlen = new->typepointers[new->ntypes];
  if (stringlen == 0) {
    new->typestrings = (char *) NULL;
  } else {
    new->typestrings = (char *) CALLOC(stringlen,sizeof(char));
    offset += sizeof(char)*FREAD_CHARS(new->typestrings,stringlen,fp);
  }
  debug(
	printf("typestrings:\n");
	for (i = 0; i < stringlen; i++) {
	  printf("%c",new->typestrings[i]);
	}
	printf("\n");
	);

  debug1(printf("Starting read of labelorder offset/length\n"));
  new->labelorder_offset = offset;
  new->labelorder_length = (size_t) (new->total_nintervals*sizeof(int));
  /* fprintf(stderr,"Doing a move_relative for labelorder_length %zu\n",new->labelorder_length); */
  move_relative(fp,new->labelorder_length);
  offset += new->labelorder_length;

  debug1(printf("Starting read of labelpointer offset/length\n"));
  new->labelpointers_offset = offset;
  new->labelpointers_length = (size_t) ((new->total_nintervals+1)*sizeof(UINT4));
  /* fprintf(stderr,"Doing a move_relative for labelpointer %zu\n",new->total_nintervals * sizeof(UINT4)); */
  move_relative(fp,new->total_nintervals * sizeof(UINT4));
  FREAD_UINT(&length,fp);
  new->label_length = (size_t) length;
  offset += new->labelpointers_length;

  debug1(printf("Starting read of label offset/length\n"));
  new->label_offset = offset;
  /* new->label_length computed above */
  /* fprintf(stderr,"Doing a move_relative for label_length %zu\n",new->label_length); */
  move_relative(fp,new->label_length);
  offset += new->label_length;

  debug1(printf("Starting read of annotpointers offset/length\n"));
  new->annotpointers_offset = offset;
  new->annotpointers_length = (size_t) ((new->total_nintervals+1)*sizeof(UINT4));
  offset += new->annotpointers_length;

  new->annot_offset = offset;

#ifdef BAD_32BIT
  /* This fails if length > 4 GB */
  move_relative(fp,new->total_nintervals * sizeof(unsigned int));
  FREAD_UINT(&length,fp);
  new->annot_length = (size_t) length;
  fprintf(stderr,"Incorrect length: %u\n",length);
#else
  new->annot_length = filesize - new->annot_offset;
  /* fprintf(stderr,"annot_length: %zu\n",new->annot_length); */
#endif

#if 0
  /* To do this check, we need to get stringlen for annotation similarly to that for labels */
  last_offset = offset + sizeof(char)*stringlen;
  if (last_offset != filesize) {
    fprintf(stderr,"Problem with last_offset (%lld) not equal to filesize = (%lld)\n",
	    (long long int) last_offset,(long long int) filesize);
    exit(9);
  }
#endif

  return;
}

static void
read_words_debug (off_t offset, off_t filesize, FILE *fp, T new) {
  off_t stringlen;
  UINT4 length;
  int i;
#if 0
  off_t last_offset;
#endif

  new->typepointers = (unsigned int *) CALLOC(new->ntypes+1,sizeof(unsigned int));
  offset += sizeof(int)*FREAD_UINTS(new->typepointers,new->ntypes+1,fp);
  printf("typepointers:");
  for (i = 0; i < new->ntypes+1; i++) {
    printf(" %u",new->typepointers[i]);
  }
  printf("\n");

  stringlen = new->typepointers[new->ntypes];
  if (stringlen == 0) {
    new->typestrings = (char *) NULL;
  } else {
    new->typestrings = (char *) CALLOC(stringlen,sizeof(char));
    offset += sizeof(char)*FREAD_CHARS(new->typestrings,stringlen,fp);
  }
  printf("typestrings:\n");
  for (i = 0; i < stringlen; i++) {
    printf("%c",new->typestrings[i]);
  }
  printf("\n");

  debug1(printf("Starting read of labelorder offset/length\n"));
  new->labelorder_offset = offset;
  new->labelorder_length = (size_t) (new->total_nintervals*sizeof(int));
  move_relative(fp,new->labelorder_length);
  offset += new->labelorder_length;

  debug1(printf("Starting read of labelpointers offset/length\n"));
  new->labelpointers_offset = offset;
  new->labelpointers_length = (size_t) ((new->total_nintervals+1)*sizeof(UINT4));
  move_relative(fp,new->total_nintervals * sizeof(UINT4));
  FREAD_UINT(&length,fp);
  new->label_length = (size_t) length;
  offset += new->labelpointers_length;

  fprintf(stderr,"label_length: %zu\n",new->label_length);
  debug1(printf("Starting read of label offset/length\n"));
  new->label_offset = offset;
  /* new->label_length computed above */
  move_relative(fp,new->label_length);
  offset += new->label_length;

  debug1(printf("Starting read of annotpointers offset/length\n"));
  new->annotpointers_offset = offset;
  new->annotpointers_length = (size_t) ((new->total_nintervals+1)*sizeof(UINT4));
  offset += new->annotpointers_length;

  new->annot_offset = offset;

#ifdef BAD_32BIT
  /* This fails if length > 4 GB */
  move_relative(fp,new->total_nintervals * sizeof(unsigned int));
  FREAD_UINT(&length,fp);
  new->annot_length = (size_t) length;
  fprintf(stderr,"Incorrect length: %u\n",length);
#else
  new->annot_length = filesize - new->annot_offset;
  fprintf(stderr,"annot_length: %zu\n",new->annot_length);
#endif

#if 0
  /* To do this check, we need to get stringlen for annotation similarly to that for labels */
  last_offset = offset + sizeof(char)*stringlen;
  if (last_offset != filesize) {
    fprintf(stderr,"Problem with last_offset (%lld) not equal to filesize = (%lld)\n",
	    (long long int) last_offset,(long long int) filesize);
    exit(9);
  }
#endif

  return;
}

/* This function only assigns pointers.  Subsequent accesses to
   memory, other than char *, still need to be read correctly
   by bigendian machines */
#ifdef HAVE_MMAP
static bool
mmap_annotations (char *filename, T new, bool readonlyp) {
  int remainder;

  if (readonlyp == true) {
    if ((new->fd = open(filename,O_RDONLY,0764)) < 0) {
      fprintf(stderr,"Error: can't open file %s with open for reading\n",filename);
      exit(9);
    }

    new->labelorder_mmap = Access_mmap_offset(&remainder,new->fd,new->labelorder_offset,new->labelorder_length,
					      sizeof(char),/*randomp*/true);
    debug(fprintf(stderr,"labelorder_mmap is %p\n",new->labelorder_mmap));
    new->labelorder = (int *) &(new->labelorder_mmap[remainder]);
    new->labelorder_length += (size_t) remainder;

    new->labelpointers_mmap = Access_mmap_offset(&remainder,new->fd,new->labelpointers_offset,new->labelpointers_length,
						 sizeof(char),/*randomp*/true);
    debug(fprintf(stderr,"labelpointers_mmap is %p\n",new->labelpointers_mmap));
    new->labelpointers = (UINT4 *) &(new->labelpointers_mmap[remainder]);
    new->labelpointers_length += (size_t) remainder;

    new->label_mmap = Access_mmap_offset(&remainder,new->fd,new->label_offset,new->label_length,
					  sizeof(char),/*randomp*/true);
    debug(fprintf(stderr,"labels_mmap is %p\n",new->label_mmap));
    new->labels = (char *) &(new->label_mmap[remainder]);
    new->label_length += (size_t) remainder;

    new->annotpointers_mmap = Access_mmap_offset(&remainder,new->fd,new->annotpointers_offset,new->annotpointers_length,
						 sizeof(char),/*randomp*/true);
    debug(fprintf(stderr,"annotpointers_mmap is %p\n",new->annotpointers_mmap));
    new->annotpointers = (UINT4 *) &(new->annotpointers_mmap[remainder]);
    new->annotpointers_length += (size_t) remainder;

    new->annot_mmap = Access_mmap_offset(&remainder,new->fd,new->annot_offset,new->annot_length,
					 sizeof(char),/*randomp*/true);
    debug(fprintf(stderr,"annots_mmap is %p\n",new->annot_mmap));
    new->annotations = (char *) &(new->annot_mmap[remainder]);
    new->annot_length += (size_t) remainder;

  } else {
    if ((new->fd = open(filename,O_RDWR,0764)) < 0) {
      fprintf(stderr,"Error: can't open file %s with open for reading/writing\n",filename);
      exit(9);
    }

    new->labelorder_mmap = Access_mmap_offset_rw(&remainder,new->fd,new->labelorder_offset,new->labelorder_length,
						 sizeof(char),/*randomp*/true);
    new->labelorder = (int *) &(new->labelorder_mmap[remainder]);
    new->labelorder_length += (size_t) remainder;

    new->labelpointers_mmap = Access_mmap_offset_rw(&remainder,new->fd,new->labelpointers_offset,new->labelpointers_length,
						    sizeof(char),/*randomp*/true);
    new->labelpointers = (UINT4 *) &(new->labelpointers_mmap[remainder]);
    new->labelpointers_length += (size_t) remainder;

    new->label_mmap = Access_mmap_offset_rw(&remainder,new->fd,new->label_offset,new->label_length,
					    sizeof(char),/*randomp*/true);
    new->labels = (char *) &(new->label_mmap[remainder]);
    new->label_length += (size_t) remainder;

    new->annotpointers_mmap = Access_mmap_offset_rw(&remainder,new->fd,new->annotpointers_offset,new->annotpointers_length,
						    sizeof(char),/*randomp*/true);
    new->annotpointers = (UINT4 *) &(new->annotpointers_mmap[remainder]);
    new->annotpointers_length += (size_t) remainder;

    new->annot_mmap = Access_mmap_offset_rw(&remainder,new->fd,new->annot_offset,new->annot_length,
					    sizeof(char),/*randomp*/true);
    new->annotations = (char *) &(new->annot_mmap[remainder]);
    new->annot_length += (size_t) remainder;
  }

  if (new->labelorder == NULL || new->labelpointers == NULL || new->labels == NULL) {
    fprintf(stderr,"Memory mapping failed in reading IIT file %s.  Using slow file IO instead.\n",filename);
    return false;
  }

  if (new->annotpointers == NULL || new->annotations == NULL) {
    fprintf(stderr,"Memory mapping failed in reading IIT file %s.  Using slow file IO instead.\n",filename);
    return false;
  }

  return true;
}
#endif


/* Used if access is FILEIO.  Subsequent accesses by bigendian
   machines to anything but (char *) will still need to convert. */
static void
read_annotations (T new) {

  file_move_absolute(new->fd,new->labelorder_offset,sizeof(int),/*n*/0);
  new->labelorder = (int *) CALLOC(new->total_nintervals,sizeof(int));
  read(new->fd,new->labelorder,new->total_nintervals*sizeof(int));

  file_move_absolute(new->fd,new->labelpointers_offset,sizeof(UINT4),/*n*/0);
  new->labelpointers = (UINT4 *) CALLOC(new->total_nintervals+1,sizeof(UINT4));
  read(new->fd,new->labelpointers,(new->total_nintervals+1)*sizeof(UINT4));

  file_move_absolute(new->fd,new->label_offset,sizeof(char),/*n*/0);
  new->labels = (char *) CALLOC(new->label_length,sizeof(char));
  read(new->fd,new->labels,new->label_length*sizeof(char));

  file_move_absolute(new->fd,new->annotpointers_offset,sizeof(UINT4),/*n*/0);
  new->annotpointers = (UINT4 *) CALLOC(new->total_nintervals+1,sizeof(UINT4));
  read(new->fd,new->annotpointers,(new->total_nintervals+1)*sizeof(UINT4));

  file_move_absolute(new->fd,new->annot_offset,sizeof(char),/*n*/0);
  new->annotations = (char *) CALLOC(new->annot_length,sizeof(char));
  read(new->fd,new->annotations,new->annot_length*sizeof(char));

  return;
}



T
Univ_IIT_read (char *filename, bool readonlyp, bool add_iit_p) {
  T new;
  FILE *fp;
  char *newfile = NULL;
  off_t offset = 0, filesize;

  /* printf("Reading IIT file %s\n",filename); */
  if (add_iit_p == true) {
    newfile = (char *) CALLOC(strlen(filename)+strlen(".iit")+1,sizeof(char));
    sprintf(newfile,"%s.iit",filename);
    if ((fp = FOPEN_READ_BINARY(newfile)) != NULL) {
      filename = newfile;
    } else if ((fp = FOPEN_READ_BINARY(filename)) == NULL) {
      /* fprintf(stderr,"Cannot open IIT file %s or %s\n",filename,newfile); */
      FREE(newfile);
      return NULL;
    }
  } else if ((fp = FOPEN_READ_BINARY(filename)) == NULL) {
    fprintf(stderr,"Cannot open IIT file %s\n",filename);
    return NULL;
  }

  new = (T) MALLOC(sizeof(*new));

  filesize = Access_filesize(filename);

  if (FREAD_INT(&new->total_nintervals,fp) < 1) {
    fprintf(stderr,"IIT file %s appears to be empty\n",filename);
    return NULL;
  } else if ((offset += sizeof(int)) > filesize) {
    fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after first byte %lld, filesize %lld).  Did you generate it using iit_store?\n",
	    filename,(long long int) offset,(long long int) filesize);
    return NULL;
  }

  if (new->total_nintervals < 0) {
#ifdef LARGE_GENOMES
    new->coord_values_8p = true;
    new->total_nintervals = -new->total_nintervals;
#elif defined(UTILITYP)
    new->coord_values_8p = true;
    new->total_nintervals = -new->total_nintervals;
#else
    fprintf(stderr,"This is a large genome of more than 2^32 (4 billion) bp.\n");
#ifdef GSNAP
    fprintf(stderr,"You should run gsnapl instead.\n");
#else
    fprintf(stderr,"You should run gmapl instead.\n");
#endif
    exit(9);
#endif

  } else if (new->total_nintervals > 0) {
    new->coord_values_8p = false;

  } else {
    abort();
  }

  debug(printf("version: 1\n"));
  debug(printf("total_nintervals: %d\n",new->total_nintervals));
  debug(printf("coord values 8p: %d\n",new->coord_values_8p));


  if (FREAD_INT(&new->ntypes,fp) < 1) {
    fprintf(stderr,"IIT file %s appears to be truncated\n",filename);
    return NULL;
  } else if (new->ntypes < 0) {
    fprintf(stderr,"IIT file %s appears to have a negative number of types\n",filename);
    return NULL;
  } else if ((offset += sizeof(int)) > filesize) {
    fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after ntypes %lld, filesize %lld).  Did you generate it using iit_store?\n",
	    filename,(long long int) offset,(long long int) filesize);
    return NULL;
  }
  debug(printf("ntypes: %d\n",new->ntypes));

  if (FREAD_INT(&new->nnodes,fp) < 1) {
    fprintf(stderr,"IIT file %s appears to be truncated\n",filename);
    return NULL;
  } else if (new->nnodes < 0) {
    fprintf(stderr,"IIT file %s appears to have a negative number of nodes\n",filename);
    return NULL;
  } else if ((offset += sizeof(int)) > filesize) {
    fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after nnodes %lld, filesize %lld).  Did you generate it using iit_store?\n",
	    filename,(long long int) offset,(long long int) filesize);
    return NULL;
  }

  /* new->divsort = NO_SORT; */
  offset = read_tree_univ(offset,filesize,fp,filename,new);

  new->intervals = (struct Univinterval_T *) CALLOC(new->total_nintervals,sizeof(struct Univinterval_T));
  offset = read_intervals_univ(offset,filesize,fp,filename,new);

  read_words(offset,filesize,fp,new);
  fclose(fp);

#ifndef HAVE_MMAP
  debug1(printf("No mmap available.  Reading annotations\n"));
  new->access = FILEIO;
  new->fd = Access_fileio(filename);
  read_annotations(new);
  close(new->fd);
  /* pthread_mutex_init(&new->read_mutex,NULL); */
#else
  debug1(printf("mmap available.  Setting up pointers to annotations\n"));
  new->access = MMAPPED;
  if (mmap_annotations(filename,new,readonlyp) == false) {
    debug1(printf("  Failed.  Reading annotations\n"));
    new->access = FILEIO;
    new->fd = Access_fileio(filename);
    read_annotations(new);
    close(new->fd);
    /* pthread_mutex_init(&new->read_mutex,NULL); */
  }
#endif
    
  if (newfile != NULL) {
    FREE(newfile);
  }

  return new;
}


void
Univ_IIT_debug (char *filename) {
  T new;
  FILE *fp;
  char *newfile = NULL;
  off_t offset = 0, filesize;
  bool add_iit_p = false;

  if (add_iit_p == true) {
    newfile = (char *) CALLOC(strlen(filename)+strlen(".iit")+1,sizeof(char));
    sprintf(newfile,"%s.iit",filename);
    if ((fp = FOPEN_READ_BINARY(newfile)) != NULL) {
      filename = newfile;
    } else if ((fp = FOPEN_READ_BINARY(filename)) == NULL) {
      /* fprintf(stderr,"Cannot open IIT file %s or %s\n",filename,newfile); */
      FREE(newfile);
      return;
    }
  } else if ((fp = FOPEN_READ_BINARY(filename)) == NULL) {
    /* fprintf(stderr,"Cannot open IIT file %s\n",filename); */
    return;
  }

  new = (T) MALLOC(sizeof(*new));

  filesize = Access_filesize(filename);

  if (FREAD_INT(&new->total_nintervals,fp) < 1) {
    fprintf(stderr,"IIT file %s appears to be empty\n",filename);
    return;
  } else if ((offset += sizeof(int)) > filesize) {
    fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after first byte %lld, filesize %lld).  Did you generate it using iit_store?\n",
	    filename,(long long int) offset,(long long int) filesize);
    return;
  }

  if (new->total_nintervals < 0) {
#ifdef HAVE_64_BIT
    new->coord_values_8p = true;
    new->total_nintervals = -new->total_nintervals;
#else
    fprintf(stderr,"IIT file has 64-bit coordinates, but this machine is only 32-bit, so it cannot process it.\n");
    exit(9);
#endif

  } else if (new->total_nintervals > 0) {
    new->coord_values_8p = false;

  } else {
    abort();
  }

  printf("version: 1\n");
  printf("total_nintervals: %d\n",new->total_nintervals);
  printf("coord values 8p: %d\n",new->coord_values_8p);


  if (FREAD_INT(&new->ntypes,fp) < 1) {
    fprintf(stderr,"IIT file %s appears to be truncated\n",filename);
    return;
  } else if (new->ntypes < 0) {
    fprintf(stderr,"IIT file %s appears to have a negative number of types\n",filename);
    return;
  } else if ((offset += sizeof(int)) > filesize) {
    fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after ntypes %lld, filesize %lld).  Did you generate it using iit_store?\n",
	    filename,(long long int) offset,(long long int) filesize);
    return;
  }
  printf("ntypes: %d\n",new->ntypes);

  if (FREAD_INT(&new->nnodes,fp) < 1) {
    fprintf(stderr,"IIT file %s appears to be truncated\n",filename);
    return;
  } else if (new->nnodes < 0) {
    fprintf(stderr,"IIT file %s appears to have a negative number of nodes\n",filename);
    return;
  } else if ((offset += sizeof(int)) > filesize) {
    fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after nnodes %lld, filesize %lld).  Did you generate it using iit_store?\n",
	    filename,(long long int) offset,(long long int) filesize);
    return;
  }

  /* new->divsort = NO_SORT; */
  offset = read_tree_univ(offset,filesize,fp,filename,new);

  new->intervals = (struct Univinterval_T *) CALLOC(new->total_nintervals,sizeof(struct Univinterval_T));
  offset = read_intervals_univ(offset,filesize,fp,filename,new);

  read_words_debug(offset,filesize,fp,new);
  fclose(fp);

#ifndef HAVE_MMAP
  debug1(printf("No mmap available.  Reading annotations\n"));
  new->access = FILEIO;
  new->fd = Access_fileio(filename);
  read_annotations(new);
  close(new->fd);
  /* pthread_mutex_init(&new->read_mutex,NULL); */
#else
  debug1(printf("mmap available.  Setting up pointers to annotations\n"));
  new->access = MMAPPED;
  if (mmap_annotations(filename,new,/*readonlyp*/true) == false) {
    debug1(printf("  Failed.  Reading annotations\n"));
    new->access = FILEIO;
    new->fd = Access_fileio(filename);
    read_annotations(new);
    close(new->fd);
    /* pthread_mutex_init(&new->read_mutex,NULL); */
  }
#endif

  if (newfile != NULL) {
    FREE(newfile);
  }

  Univ_IIT_free(&new);

  return;
}


/************************************************************************/

static void 
fnode_query_aux (int *min, int *max, T this, int nodeindex, Univcoord_T x) {
  int lambda;
  Univ_FNode_T node;

  if (nodeindex == -1) {
    return;
  }

  node = &(this->nodes[nodeindex]);
  debug(printf("Entered fnode_query_aux with nodeindex %d: a %d, b %d, leftindex %d, rightindex %d, value %llu\n",
	       nodeindex,node->a,node->b,node->leftindex,node->rightindex,(unsigned long long) node->value));

  if (x == node->value) {
    debug(printf("%lluD:\n",(unsigned long long) node->value));
    if (node->a < *min) {
      *min = node->a;
    }
    if (node->b > *max) {
      *max = node->b;
    }
    return;
  } else if (x < node->value) {
    debug(printf("x %llu < node->value %llu\n",(unsigned long long) x,(unsigned long long) node->value));
    fnode_query_aux(&(*min),&(*max),this,node->leftindex,x);
    debug(printf("%lluL:\n",(unsigned long long) node->value));
    if (node->a < *min) {
      *min = node->a;
    }
    for (lambda = node->a; lambda <= node->b; lambda++) {
      debug(printf("Looking at lambda %d, segment %d\n",
		   lambda,this->sigmas[lambda]));
      if (Univinterval_is_contained(x,this->intervals,this->sigmas[lambda]) == true) {
	if (lambda > *max) {
	  *max = lambda;
	}
      } else {
	return;
      }
    }
    return;
  } else { 
    /* (node->value < x) */
    debug(printf("x %llu > node->value %llu\n",(unsigned long long) x,(unsigned long long) node->value));
    fnode_query_aux(&(*min),&(*max),this,node->rightindex,x);
    debug(printf("%lluR:\n",(unsigned long long) node->value));
    if (node->b > *max) {
      *max = node->b;
    }
    for (lambda = node->b; lambda >= node->a; lambda--) {
      debug(printf("Looking at lambda %d, segment %d\n",
		   lambda,this->omegas[lambda]));
      if (Univinterval_is_contained(x,this->intervals,this->omegas[lambda]) == true) {
	if (lambda < *min) {
	  *min = lambda;
	}
      } else {
	return;
      }
    }
    return;
  }
}

/************************************************************************/

int *
Univ_IIT_find (int *nmatches, T this, char *label) {
  int *matches = NULL, j;
  int low, middle, high, recno;
  bool foundp = false;
  int cmp;

  low = 0;
  high = this->total_nintervals;
  *nmatches = 0;

  while (!foundp && low < high) {
    middle = low + (high - low)/2;

#ifdef WORDS_BIGENDIAN
    cmp = strcmp(label,&(this->labels[Bigendian_convert_uint(this->labelpointers[Bigendian_convert_int(this->labelorder[middle])])]));
#else
    cmp = strcmp(label,&(this->labels[this->labelpointers[this->labelorder[middle]]]));
#endif

    if (cmp < 0) {
      high = middle;
    } else if (cmp > 0) {
      low = middle + 1;
    } else {
      foundp = true;
    }
  }

  if (foundp == true) {
    low = middle;
#ifdef WORDS_BIGENDIAN
    while (low-1 >= 0 && 
	   !strcmp(label,&(this->labels[Bigendian_convert_uint(this->labelpointers[Bigendian_convert_int(this->labelorder[low-1])])]))) {
      low--;
    }
#else
    while (low-1 >= 0 && 
	   !strcmp(label,&(this->labels[this->labelpointers[this->labelorder[low-1]]]))) {
      low--;
    }
#endif
    
    high = middle;
#ifdef WORDS_BIGENDIAN
    while (high+1 < this->total_nintervals && 
	   !strcmp(label,&(this->labels[Bigendian_convert_uint(this->labelpointers[Bigendian_convert_int(this->labelorder[high+1])])]))) {
      high++;
    }
#else
    while (high+1 < this->total_nintervals && 
	   !strcmp(label,&(this->labels[this->labelpointers[this->labelorder[high+1]]]))) {
      high++;
    }
#endif

    
    *nmatches = high - low + 1;
    if (*nmatches > 0) {
      matches = (int *) CALLOC(*nmatches,sizeof(int));
      j = 0;
      for (recno = low; recno <= high; recno++) {
#ifdef WORDS_BIGENDIAN
#ifdef DEBUG
	printf("Pushing %d:%d\n",recno,Bigendian_convert_int(this->labelorder[recno]));
#endif
	matches[j++] = Bigendian_convert_int(this->labelorder[recno])+1;
	
#else
#ifdef DEBUG
	printf("Pushing %d:%d\n",recno,this->labelorder[recno]);
#endif
	matches[j++] = this->labelorder[recno]+1;
#endif
      }
    }
  }

  return matches;
}

/* Slow.  Used before binary search method above. */
int
Univ_IIT_find_linear (T this, char *label) {
  int i;
  char *p;

  for (i = 0; i < this->total_nintervals; i++) {
#ifdef WORDS_BIGENDIAN
    p = &(this->labels[Bigendian_convert_uint(this->labelpointers[i])]);
#else
    p = &(this->labels[this->labelpointers[i]]);
#endif
    while (isspace((int) *p)) {
      p++;
    }
    if (!strcmp(label,p)) {
      return i + 1;
    }
  }

  return -1;
}

/* Returns 1-based index */
int
Univ_IIT_find_one (T this, char *label) {
  int index;
  int *matches, nmatches;

  matches = Univ_IIT_find(&nmatches,this,label);
  if (nmatches == 0) {
    /*
    fprintf(stderr,"Expected one match for %s, but got 0\n",
	    label);
    */
    index = -1;
  } else {
    if (nmatches > 1) {
      fprintf(stderr,"Expected one match for %s, but got %d\n",
	      label,nmatches);
    }
    index = matches[0];
    FREE(matches);
  }

  return index;
}    


/************************************************************************/


static const Except_T iit_error = { "IIT problem" };

static int
int_compare (const void *a, const void *b) {
  int x = * (int *) a;
  int y = * (int *) b;

  if (x < y) {
    return -1;
  } else if (y < x) {
    return +1;
  } else {
    return 0;
  }
}

int *
Univ_IIT_get (int *nmatches, T this, Univcoord_T x, Univcoord_T y) {
  int *sorted, *matches = NULL, *uniq, neval, nuniq, i;
  int lambda, prev;
  int min1, max1 = 0, min2, max2 = 0;
  int nintervals;

  if ((nintervals = this->total_nintervals) == 0) {
    *nmatches = 0;
    return (int *) NULL;
  } else {
    min1 = min2 = nintervals + 1;
  }

  debug(printf("Entering Univ_IIT_get with query %u %u\n",x,y));
  fnode_query_aux(&min1,&max1,this,0,x);
  fnode_query_aux(&min2,&max2,this,0,y);
  debug(printf("min1=%d max1=%d  min2=%d max2=%d\n",min1,max1,min2,max2));

  *nmatches = 0;
  if (max2 >= min1) {
    neval = (max2 - min1 + 1) + (max2 - min1 + 1);
    matches = (int *) CALLOC(neval,sizeof(int));
    uniq = (int *) CALLOC(neval,sizeof(int));

    i = 0;
    for (lambda = min1; lambda <= max2; lambda++) {
      matches[i++] = this->sigmas[lambda];
      matches[i++] = this->omegas[lambda];
    }

    /* Eliminate duplicates */
    qsort(matches,neval,sizeof(int),int_compare);
    nuniq = 0;
    prev = 0;
    debug(printf("unique segments in lambda %d to %d:",min1,max2));
    for (i = 0; i < neval; i++) {
      if (matches[i] != prev) {
	debug(printf(" %d",matches[i]));
	uniq[nuniq++] = matches[i];
	prev = matches[i];
      }
    }
    debug(printf("\n"));

    for (i = 0; i < nuniq; i++) {
      if (Univinterval_overlap_p(x,y,this->intervals,uniq[i]) == true) {
	matches[(*nmatches)++] = uniq[i];
	debug(printf("Pushing overlapping segment %d (%u..%u)\n",uniq[i],
		     Univinterval_low(&(this->intervals[uniq[i]-1])),
		     Univinterval_high(&(this->intervals[uniq[i]-1]))));
      } else {
	debug(printf("Not pushing non-overlapping segment %d (%u..%u)\n",uniq[i],
		     Univinterval_low(&(this->intervals[uniq[i]-1])),
		     Univinterval_high(&(this->intervals[uniq[i]-1]))));
      }
    }

    FREE(uniq);
  }

  return matches;

#if 0
  /* For some reason, sort_matches_by_position is annotated as being for version 3 and later */
  if (sortp == false) {
    return matches;
  } else if (this->version <= 2) {
    sorted = sort_matches_by_type(this,matches,*nmatches,/*alphabetizep*/true);
    FREE(matches);
    return sorted;
  } else {
    sorted = sort_matches_by_position(this,matches,*nmatches);
    FREE(matches);
    return sorted;
  }
#endif
}


/* Guaranteed to return one result, even if coordinate is out of bounds */
int
Univ_IIT_get_one (T this, Univcoord_T x, Univcoord_T y) {
  int lambda;
  int min1, max1 = 0, min2, max2 = 0;
  bool stopp;
  Univinterval_T interval;

  min1 = min2 = this->total_nintervals + 1;

  debug(printf("Entering Univ_IIT_get_one with query %llu %llu\n",(unsigned long long) x,(unsigned long long) y));
  fnode_query_aux(&min1,&max1,this,0,x);
  fnode_query_aux(&min2,&max2,this,0,y);
  debug(printf("min1=%d max1=%d  min2=%d max2=%d\n",min1,max1,min2,max2));

  if (max2 >= min1) {
    for (lambda = min1; lambda <= max2; lambda++) {
      if (Univinterval_overlap_p(x,y,this->intervals,this->sigmas[lambda]) == true) {
	return this->sigmas[lambda];
      }
    }
    for (lambda = min1; lambda <= max2; lambda++) {
      if (Univinterval_overlap_p(x,y,this->intervals,this->omegas[lambda]) == true) {
	return this->omegas[lambda];
      }
    }
  }

  /* fprintf(stderr,"Expected one match for %u--%u, but got none\n",x,y); */
  /* If we miss (e.g., for circular chromosome), then report the chromosome below */
  /* Look at betas or omegas for left flank */
  lambda = min1 - 1;
  stopp = false;
  while (lambda >= 1 && stopp == false) {
    interval = &(this->intervals[this->omegas[lambda]-1]);
    if (Univinterval_high(interval) >= x) {
      lambda--;
    } else {
      return this->omegas[lambda];
    }
  }

  return this->omegas[/*lambda*/1];
}

/* Generally called where intervals don't overlap, like chromosomes,
   and where x == y. */
/*
int
Univ_IIT_get_one_safe (T this, Univcoord_T x, Univcoord_T y) {
  int index;
  int *matches, nmatches;

  matches = Univ_IIT_get(&nmatches,this,x,y,sortp);
  if (nmatches != 1) {
    fprintf(stderr,"Expected one match for %u--%u, but got %d\n",
	    x,y,nmatches);
    abort();
  }
  index = matches[0];
  FREE(matches);
  return index;
}
*/


/* Note: Procedure call from get-genome.c needed to subtract 1 from
   position and then add 1 to chrpos */
char *
Univ_IIT_string_from_position (Chrpos_T *chrpos, Univcoord_T position, T chromosome_iit) {
  char *string, *chrstring;
  int index;
  bool allocp;

  index = Univ_IIT_get_one(chromosome_iit,position,position);
  *chrpos = position - Univinterval_low(Univ_IIT_interval(chromosome_iit,index));
  chrstring = Univ_IIT_label(chromosome_iit,index,&allocp); 
  if (allocp == true) {
    return chrstring;
  } else {
    string = (char *) CALLOC(strlen(chrstring)+1,sizeof(char));
    strcpy(string,chrstring);
    return string;
  }
}


