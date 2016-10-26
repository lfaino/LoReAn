static char rcsid[] = "$Id: iit-read.c 164702 2015-05-01 20:22:25Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "iit-read.h"
#include "iitdef.h"

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



#define T IIT_T

static void
file_move_absolute (int fd, off_t offset, off_t objsize, Chrpos_T n) {
  off_t position = offset + n*objsize;

  if (lseek(fd,position,SEEK_SET) < 0) {
    perror("Error in gmap, file_move_label");
    exit(9);
  }
  return;
}


bool
IIT_universalp (char *filename, bool add_iit_p) {
  char *newfile;
  FILE *fp;
  int total_nintervals;

  if (add_iit_p == true) {
    newfile = (char *) CALLOC(strlen(filename)+strlen(".iit")+1,sizeof(char));
    sprintf(newfile,"%s.iit",filename);
    if ((fp = FOPEN_READ_BINARY(newfile)) != NULL) {
      filename = newfile;
    } else if ((fp = FOPEN_READ_BINARY(filename)) == NULL) {
      /* fprintf(stderr,"Cannot open IIT file %s or %s\n",filename,newfile); */
      FREE(newfile);
      return false;
    }
  } else if ((fp = FOPEN_READ_BINARY(filename)) == NULL) {
    /* fprintf(stderr,"Cannot open IIT file %s\n",filename); */
    return false;
  }

  if (FREAD_INT(&total_nintervals,fp) < 1) {
    fprintf(stderr,"IIT file %s appears to be empty\n",filename);
    fclose(fp);
    if (add_iit_p == true) {
      FREE(newfile);
    }
    return false;
  } else if (total_nintervals == 0) {
    /* Need to use Univ_IIT_read instead */
    fclose(fp);
    if (add_iit_p == true) {
      FREE(newfile);
    }
    return false;
  } else {
    fclose(fp);
    if (add_iit_p == true) {
      FREE(newfile);
    }
    return true;
  }
}


bool
IIT_valuep (T this) {
  return this->valuep;
}


char *
IIT_name (T this) {
  return this->name;
}

int
IIT_version (T this) {
  return this->version;
}

int
IIT_total_nintervals (T this) {
  return this->total_nintervals;
}

int
IIT_nintervals (T this, int divno) {
  return this->nintervals[divno];
}


int
IIT_ntypes (T this) {
  return this->ntypes;
}

int
IIT_nfields (T this) {
  return this->nfields;
}


Chrpos_T
IIT_length (T this, int index) {
  Interval_T interval;

  interval = &(this->intervals[0][index-1]);
  return Interval_length(interval);
}


Chrpos_T
IIT_divlength (T this, char *divstring) {
  Chrpos_T max = 0U;
  Interval_T interval;
  int divno, i;

  divno = IIT_divint(this,divstring);
  for (i = 0; i < this->nintervals[divno]; i++) {
    interval = &(this->intervals[divno][i]);
    if (Interval_high(interval) > max) {
      max = Interval_high(interval);
    }
  }
  /* Convert from zero-based coordinate */
  return max+1U;
}


/* Assumes intervals are stored using universal coordinates */
Chrpos_T
IIT_totallength (T this) {
  Chrpos_T max = 0U;
  Interval_T interval;
  int divno, i;

  for (divno = 0; divno < this->ndivs; divno++) {
    for (i = 0; i < this->nintervals[divno]; i++) {
      interval = &(this->intervals[divno][i]);
      if (Interval_high(interval) > max) {
	max = Interval_high(interval);
      }
    }
  }
  /* Convert from zero-based coordinate */
  return max+1U;
}


Interval_T
IIT_interval (T this, int index) {
  assert(index <= this->total_nintervals);
  return &(this->intervals[0][index-1]); /* Convert to 0-based */
}

Chrpos_T
IIT_interval_low (T this, int index) {
  Interval_T interval;

  assert(index <= this->total_nintervals);
  interval = &(this->intervals[0][index-1]);
  return Interval_low(interval);
}

Chrpos_T
IIT_interval_high (T this, int index) {
  Interval_T interval;

  assert(index <= this->total_nintervals);
  interval = &(this->intervals[0][index-1]);
  return Interval_high(interval);
}

Chrpos_T
IIT_interval_length (T this, int index) {
  Interval_T interval;

  assert(index <= this->total_nintervals);
  interval = &(this->intervals[0][index-1]);
  return Interval_length(interval);
}

int
IIT_interval_type (T this, int index) {
  Interval_T interval;

  assert(index <= this->total_nintervals);
  interval = &(this->intervals[0][index-1]);
  return Interval_type(interval);
}


int
IIT_interval_sign (T this, int index) {
  Interval_T interval;

  assert(index <= this->total_nintervals);
  interval = &(this->intervals[0][index-1]);
  return Interval_sign(interval);
}


/* chrhigh is one past the highest position in the chromosome */
void
IIT_interval_bounds (Chrpos_T *low, Chrpos_T *high, Chrpos_T *length, T this,
		     int index, int circular_typeint) {
  Interval_T interval;

  assert(index > 0);
  assert(index <= this->total_nintervals);

  interval = &(this->intervals[0][index-1]);
  *low = Interval_low(interval);
  *length = Interval_length(interval);
  if (Interval_type(interval) == circular_typeint) {
    *high = Interval_high(interval) + 1 + (*length);
  } else {
    *high = Interval_high(interval) + 1;
  }
  return;
}

int
IIT_index (T this, int divno, int i) {
  return this->cum_nintervals[divno] + i + 1; /* 1-based */
}



int
IIT_ndivs (T this) {
  return this->ndivs;
}

/* The iit file has a '\0' after each string, so functions know where
   it ends */
char *
IIT_divstring (T this, int divno) {
  UINT4 start;

  start = this->divpointers[divno];
  return &(this->divstrings[start]);
}

int
IIT_divint (T this, char *divstring) {
  int i = 0;			/* Actually divstring for divno 0 is NULL */
  UINT4 start;

  if (divstring == NULL) {
    return 0;
  } else if (divstring[0] == '\0') {
    return 0;
  } else {
    while (i < this->ndivs) {
      start = this->divpointers[i];
      if (!strcmp(divstring,&(this->divstrings[start]))) {
	return i;
      }
      i++;
    }

    return -1;
  }
}

char *
IIT_divstring_from_index (T this, int index) {
  int divno = 1;
  UINT4 start;

  while (divno <= this->ndivs) {
    /* Checked on existing iit file to confirm we need >= and not > */
    if (this->cum_nintervals[divno] >= index) {
      start = this->divpointers[divno-1];
      return &(this->divstrings[start]);
    }
    divno++;
  }
  return (char *) NULL;
}


/* The iit file has a '\0' after each string, so functions know where
   it ends */
char *
IIT_typestring (T this, int type) {
  UINT4 start;

  start = this->typepointers[type];
  return &(this->typestrings[start]);
}

int
IIT_typeint (T this, char *typestring) {
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
IIT_fieldstring (T this, int fieldint) {
  UINT4 start;

  start = this->fieldpointers[fieldint];
  return &(this->fieldstrings[start]);
}

int
IIT_fieldint (T this, char *fieldstring) {
  int i = 0;
  UINT4 start;

  while (i < this->nfields) {
    start = this->fieldpointers[i];
    if (!strcmp(fieldstring,&(this->fieldstrings[start]))) {
      return i;
    }
    i++;
  }

  return -1;
}


char *
IIT_label (T this, int index, bool *allocp) {
  int recno;
#ifdef HAVE_64_BIT
  UINT8 start;
#else
  UINT4 start;
#endif

  recno = index - 1; /* Convert to 0-based */

#ifdef WORDS_BIGENDIAN
#ifdef HAVE_64_BIT
  if (this->label_pointers_8p == true) {
    start = Bigendian_convert_uint8(this->labelpointers8[recno]);
  } else {
    start = (UINT8) Bigendian_convert_uint(this->labelpointers[recno]);
  }
#else
  start = Bigendian_convert_uint(this->labelpointers[recno]);
#endif
#else
#ifdef HAVE_64_BIT
  if (this->label_pointers_8p == true) {
    start = this->labelpointers8[recno];
  } else {
    start = (UINT8) this->labelpointers[recno];
  }
#else
  start = this->labelpointers[recno];
#endif
#endif
  *allocp = false;
  return &(this->labels[start]);
}


static char EMPTY_STRING[1] = {'\0'};

/* The iit file has a '\0' after each string, so functions know where
   it ends */
/* Note: annotation itself is never allocated */
char *
IIT_annotation (char **restofheader, T this, int index, bool *alloc_header_p) {
  int recno;
  char *annotation, *p;
  int len;
#ifdef HAVE_64_BIT
  UINT8 start;
#else
  UINT4 start;
#endif

  recno = index - 1; /* Convert to 0-based */
#ifdef WORDS_BIGENDIAN
#ifdef HAVE_64_BIT
  if (this->annot_pointers_8p == true) {
    start = Bigendian_convert_uint8(this->annotpointers8[recno]);
  } else {
    start = (UINT8) Bigendian_convert_uint(this->annotpointers[recno]);
  }
#else
  start = Bigendian_convert_uint(this->annotpointers[recno]);
#endif
#else
#ifdef HAVE_64_BIT
  if (this->annot_pointers_8p == true) {
    start = this->annotpointers8[recno];
  } else {
    start = (UINT8) this->annotpointers[recno];
  }
#else
  start = this->annotpointers[recno];
#endif
#endif

  if (this->version <= 4) {
    *restofheader = EMPTY_STRING;

    *alloc_header_p = false;
    return &(this->annotations[start]);
  } else {
    /* Versions 5 and higher include rest of header with
       annotation.  Don't return initial '\n', unless annotation is empty */
    annotation = &(this->annotations[start]);
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
}

/* The iit file has a '\0' after each string, so functions know where
   it ends */
char
IIT_annotation_firstchar (T this, int index) {
  int recno;
#ifdef HAVE_64_BIT
  UINT8 start;
#else
  UINT4 start;
#endif

  recno = index - 1; /* Convert to 0-based */

#ifdef WORDS_BIGENDIAN
#ifdef HAVE_64_BIT
  if (this->annot_pointers_8p == true) {
    start = Bigendian_convert_uint8(this->annotpointers8[recno]);
  } else {
    start = (UINT8) Bigendian_convert_uint(this->annotpointers[recno]);
  }
#else
  start = Bigendian_convert_uint(this->annotpointers[recno]);
#endif
#else
#ifdef HAVE_64_BIT
  if (this->annot_pointers_8p == true) {
    start = this->annotpointers8[recno];
  } else {
    start = (UINT8) this->annotpointers[recno];
  }
#else
  start = this->annotpointers[recno];
#endif
#endif

  return this->annotations[start];
}

#ifdef HAVE_64_BIT
UINT8
#else
UINT4
#endif
IIT_annotation_strlen (T this, int index) {
  int recno;
#ifdef HAVE_64_BIT
  UINT8 start, end;
#else
  UINT4 start, end;
#endif

  recno = index - 1; /* Convert to 0-based */

#ifdef WORDS_BIGENDIAN
#ifdef HAVE_64_BIT
  if (this->annot_pointers_8p == true) {
    start = Bigendian_convert_uint8(this->annotpointers8[recno]);
    end = Bigendian_convert_uint8(this->annotpointers8[recno+1]);
  } else {
    start = (UINT8) Bigendian_convert_uint(this->annotpointers[recno]);
    end = (UINT8) Bigendian_convert_uint(this->annotpointers[recno+1]);
  }
#else
  start = Bigendian_convert_uint(this->annotpointers[recno]);
  end = Bigendian_convert_uint(this->annotpointers[recno+1]);
#endif
#else
#ifdef HAVE_64_BIT
  if (this->annot_pointers_8p == true) {
    start = this->annotpointers8[recno];
    end = this->annotpointers8[recno+1];
  } else {
    start = (UINT8) this->annotpointers[recno];
    end = (UINT8) this->annotpointers[recno+1];
  }
#else
  start = this->annotpointers[recno];
  end = this->annotpointers[recno+1];
#endif
#endif

  /*
  if (strlen(&(this->annotations[start])) != (end - start - 1)) {
    printf("Problem with %s: %d != %u\n",
    &(this->labels[this->labelpointers[recno]]),strlen(&(this->annotations[start])),end-start-1);
    abort();
  } else {
    printf("Okay %s: %d == %u\n",
    &(this->labels[this->labelpointers[recno]]),strlen(&(this->annotations[start])),end-start-1);
  }
  */

  return (end - start - 1);	/* Subtract terminal '\0' */
}

/* Always allocated */
char *
IIT_fieldvalue (T this, int index, int fieldint) {
  char *fieldvalue, *annotation, *p, *q;
  int recno, fieldno = 0, fieldlen;
#ifdef HAVE_64_BIT
  UINT8 start;
#else
  UINT4 start;
#endif
  bool allocp;

  recno = index - 1; /* Convert to 0-based */
#ifdef WORDS_BIGENDIAN
#ifdef HAVE_64_BIT
  if (this->annot_pointers_8p == true) {
    start = Bigendian_convert_uint8(this->annotpointers8[recno]);
  } else {
    start = (UINT8) Bigendian_convert_uint(this->annotpointers[recno]);
  }
#else
  start = Bigendian_convert_uint(this->annotpointers[recno]);
#endif
#else
#ifdef HAVE_64_BIT
  if (this->annot_pointers_8p == true) {
    start = this->annotpointers8[recno];
  } else {
    start = (UINT8) this->annotpointers[recno];
  }
#else
  start = this->annotpointers[recno];
#endif
#endif
  annotation = &(this->annotations[start]);
  allocp = false;

  p = annotation;

  /* Starting with version 5, annotation should have '\n' from the header line.  */
  while (*p != '\0' && *p != '\n') p++;
  if (*p == '\n') p++;

  while (*p != '\0' && fieldno < fieldint) {
    if (*p == '\n') {
      fieldno++;
    }
    p++;
  }

  if (*p == '\0') {
    fieldvalue = (char *) CALLOC(1,sizeof(char));
    fieldvalue[0] = '\0';
  } else {
    q = p;
    while (*q != '\0' && *q != '\n') {
      q++;
    }
    fieldlen = (q - p)/sizeof(char);
    fieldvalue = (char *) CALLOC(fieldlen+1,sizeof(char));
    strncpy(fieldvalue,p,fieldlen);
  }

  if (allocp == true) {
    FREE(annotation);
  }

  return fieldvalue;
}


void
IIT_dump_divstrings (FILE *fp, T this) {
  int divno;
  UINT4 start;

  /* Start with 1, because first divno has no name */
  for (divno = 1; divno < this->ndivs; divno++) {
    start = this->divpointers[divno];
    fprintf(fp,"%s ",&(this->divstrings[start]));
  }
  fprintf(fp,"\n");

  return;
}


void
IIT_dump_typestrings (FILE *fp, T this) {
  int type;
  UINT4 start;

  for (type = 0; type < this->ntypes; type++) {
    start = this->typepointers[type];
    fprintf(fp,"%d\t%s\n",type,&(this->typestrings[start]));
  }
  return;
}

void
IIT_dump_fieldstrings (FILE *fp, T this) {
  int field;
  UINT4 start;

  for (field = 0; field < this->nfields; field++) {
    start = this->fieldpointers[field];
    fprintf(fp,"%d\t%s\n",field,&(this->fieldstrings[start]));
  }
  return;
}

void
IIT_dump_labels (FILE *fp, T this) {
  int i;
#ifdef HAVE_64_BIT
  UINT8 start;
#else
  UINT4 start;
#endif
  char *label;

  for (i = 0; i < this->total_nintervals; i++) {
#ifdef WORDS_BIGENDIAN
#ifdef HAVE_64_BIT
    if (this->label_pointers_8p == true) {
      start = Bigendian_convert_uint8(this->labelpointers8[i]);
    } else {
      start = (UINT8) Bigendian_convert_uint(this->labelpointers[i]);
    }
#else
    start = Bigendian_convert_uint(this->labelpointers[i]);
#endif
#else
#ifdef HAVE_64_BIT
    if (this->label_pointers_8p == true) {
      start = this->labelpointers8[i];
    } else {
      start = (UINT8) this->labelpointers[i];
    }
#else
    start = this->labelpointers[i];
#endif
#endif
    label = &(this->labels[start]);
    fprintf(fp,"%s ",label);
  }
  fprintf(fp,"\n");
  return;
}


void
IIT_dump (T this, bool annotationonlyp, bool sortp) {
  int divno, i;
  Interval_T interval;
  char *divstring;
  char *labelptr, *annotptr, c;
  int *matches, nmatches, index;
  char *label, *annotation, *restofheader;
  bool allocp;

  if (sortp == false) {
    labelptr = this->labels;
    annotptr = this->annotations;
  }

  for (divno = 0; divno < this->ndivs; divno++) {
    divstring = IIT_divstring(this,divno);

    if (sortp == true) {
      if (this->nintervals[divno] > 0) {
	matches = IIT_get(&nmatches,this,divstring,/*x*/0,/*y*/-1U,/*sortp*/true);
	for (i = 0; i < nmatches; i++) {
	  index = matches[i];
	  label = IIT_label(this,index,&allocp);
	  printf(">%s ",label);
	  if (allocp == true) {
	    FREE(label);
	  }

	  if (divno > 0) {
	    /* zeroth divno has empty string */
	    printf("%s:",divstring);
	  }

	  interval = IIT_interval(this,index);
	  if (Interval_sign(interval) < 0) {
	    printf("%u..%u",Interval_high(interval),Interval_low(interval));
	  } else {
	    printf("%u..%u",Interval_low(interval),Interval_high(interval));
	  }
	  if (Interval_type(interval) > 0) {
	    printf(" %s",IIT_typestring(this,Interval_type(interval)));
	  }

	  annotation = IIT_annotation(&restofheader,this,index,&allocp);
	  printf("%s\n",restofheader);
	  printf("%s",annotation);

	  if (allocp == true) {
	    FREE(restofheader);
	  }
	}

	FREE(matches);
      }

    } else {
      for (i = 0; i < this->nintervals[divno]; i++) {
	printf(">");
	while ((c = *labelptr++) != '\0') {
	  printf("%c",c);
	}
	printf(" ");

	if (divno > 0) {
	  /* zeroth divno has empty string */
	  printf("%s:",divstring);
	}

	interval = &(this->intervals[divno][i]);
	if (Interval_sign(interval) < 0) {
	  printf("%u..%u",Interval_high(interval),Interval_low(interval));
	} else {
	  printf("%u..%u",Interval_low(interval),Interval_high(interval));
	}
	if (Interval_type(interval) > 0) {
	  printf(" %s",IIT_typestring(this,Interval_type(interval)));
	}
	if (this->version <= 4) {
	  printf("\n");
	  while ((c = *annotptr++) != '\0') {
	    printf("%c",c);
	  }
	} else {
	  /* Versions 5 and higher include rest of header with
	     annotation.  Don't print initial '\n', unless annotation is empty */
	  if (*annotptr == '\0') {
	    printf("\n");
	    annotptr++;
	  } else if (*annotptr == '\n') {
	    /* No rest of header */
	    while ((c = *annotptr++) != '\0') {
	      printf("%c",c);
	    }
	  } else {
	    printf(" ");
	    while ((c = *annotptr++) != '\0') {
	      printf("%c",c);
	    }
	  }
	}
      }
    }
  }

  return;
}


/* For chromosome.iit file, which is stored in version 1 */
void
IIT_dump_simple (T this) {
  int index = 0, i;
  Interval_T interval;
  Chrpos_T startpos, endpos;
  char *label;
  bool allocp;

  for (i = 0; i < this->nintervals[0]; i++) {
    interval = &(this->intervals[0][i]);
    label = IIT_label(this,index+1,&allocp);
    printf("%s\t",label);
    if (allocp == true) {
      FREE(label);
    }
    startpos = Interval_low(interval);
    endpos = startpos + Interval_length(interval) - 1U;

    printf("%u..%u\t",startpos+1U,endpos+1U);

    printf("%u",Interval_length(interval));
    if (Interval_type(interval) > 0) {
      printf("\t%s",IIT_typestring(this,Interval_type(interval)));
    }
    printf("\n");

    index++;
  }

  return;
}


#if 0
/* For higher version files, which are divided into divs */
void
IIT_dump_formatted (T this, bool directionalp) {
  int divno, index = 0, i;
  Interval_T interval;
  Chrpos_T startpos, endpos;
  char *label, *divstring, firstchar;
  bool allocp;

  for (divno = 0; divno < this->ndivs; divno++) {
    divstring = IIT_divstring(this,divno);
    for (i = 0; i < this->nintervals[divno]; i++) {
      interval = &(this->intervals[divno][i]);
      label = IIT_label(this,index+1,&allocp);
      printf("%s\t",label);
      if (allocp == true) {
	FREE(label);
      }
      startpos = Interval_low(interval);
      endpos = startpos + Interval_length(interval) - 1U;

      if (divno > 0) {
	printf("%s:",divstring);
      }
      if (directionalp == false) {
	printf("%u..%u\t",startpos+1U,endpos+1U);
      } else if (this->version <= 1) {
	firstchar = IIT_annotation_firstchar(this,index+1);
	if (firstchar == '-') {
	  printf("%u..%u\t",endpos+1U,startpos+1U);
	} else {
	  printf("%u..%u\t",startpos+1U,endpos+1U);
	}
      } else {
	if (Interval_sign(interval) < 0) {
	  printf("%u..%u\t",endpos+1U,startpos+1U);
	} else {
	  printf("%u..%u\t",startpos+1U,endpos+1U);
	}
      }

      printf("%u",Interval_length(interval));
      if (Interval_type(interval) > 0) {
	printf("\t%s",IIT_typestring(this,Interval_type(interval)));
      }
      printf("\n");

      index++;
    }
  }

  return;
}
#endif


#if 0
static int
uint_cmp (const void *x, const void *y) {
  unsigned int a = * (unsigned int *) x;
  unsigned int b = * (unsigned int *) y;

  if (a < b) {
    return -1;
  } else if (a > b) {
    return +1;
  } else {
    return 0;
  }
}

/* Need to work on */
UINT4 *
IIT_transitions (int **signs, int *nedges, T this) { 
  UINT4 *edges, *starts, *ends;
  int nintervals, i, j, k;
  Interval_T interval;
  Uintlist_T startlist = NULL, endlist = NULL;

  for (i = 0; i < this->nintervals; i++) {
    interval = &(this->intervals[i]);
    startlist = Uintlist_push(startlist,Interval_low(interval));
    endlist = Uintlist_push(endlist,Interval_high(interval));
  }

  if (Uintlist_length(startlist) == 0) {
    edges = (unsigned int *) NULL;
    *signs = (int *) NULL;
    *nedges = 0;
  } else {
    starts = Uintlist_to_array(&nintervals,startlist);
    ends = Uintlist_to_array(&nintervals,endlist);
    qsort(starts,nintervals,sizeof(unsigned int),uint_cmp);
    qsort(ends,nintervals,sizeof(unsigned int),uint_cmp);

    *nedges = nintervals+nintervals;
    *signs = (int *) CALLOC(*nedges,sizeof(int));
    edges = (unsigned int *) CALLOC(*nedges,sizeof(unsigned int));
    i = j = k = 0;
    while (i < nintervals && j < nintervals) {
      if (starts[i] <= ends[j]) {
	(*signs)[k] = +1;
	edges[k++] = starts[i++];
      } else {
	(*signs)[k] = -1;
	edges[k++] = ends[j++];
      }
    }
    while (i < nintervals) {
      (*signs)[k] = +1;
      edges[k++] = starts[i++];
    }
    while (j < nintervals) {
      (*signs)[k] = -1;
      edges[k++] = ends[j++];
    }

    FREE(ends);
    FREE(starts);
  }

  Uintlist_free(&endlist);
  Uintlist_free(&startlist);

  return edges;
}

UINT4 *
IIT_transitions_subset (int **signs, int *nedges, T this, int *indices, int nindices) { 
  UINT4 *edges, *starts, *ends;
  int nintervals, i, j, k;
  Interval_T interval;
  Uintlist_T startlist = NULL, endlist = NULL;

  for (k = 0; k < nindices; k++) {
    i = indices[k] - 1;
    interval = &(this->intervals[i]);
    startlist = Uintlist_push(startlist,Interval_low(interval));
    endlist = Uintlist_push(endlist,Interval_high(interval));
  }

  if (Uintlist_length(startlist) == 0) {
    edges = (unsigned int *) NULL;
    *signs = (int *) NULL;
    *nedges = 0;
  } else {
    starts = Uintlist_to_array(&nintervals,startlist);
    ends = Uintlist_to_array(&nintervals,endlist);
    qsort(starts,nintervals,sizeof(unsigned int),uint_cmp);
    qsort(ends,nintervals,sizeof(unsigned int),uint_cmp);

    *nedges = nintervals+nintervals;
    *signs = (int *) CALLOC(*nedges,sizeof(int));
    edges = (unsigned int *) CALLOC(*nedges,sizeof(unsigned int));
    i = j = k = 0;
    while (i < nintervals && j < nintervals) {
      if (starts[i] <= ends[j]) {
	(*signs)[k] = +1;
	edges[k++] = starts[i++];
      } else {
	(*signs)[k] = -1;
	edges[k++] = ends[j++];
      }
    }
    while (i < nintervals) {
      (*signs)[k] = +1;
      edges[k++] = starts[i++];
    }
    while (j < nintervals) {
      (*signs)[k] = -1;
      edges[k++] = ends[j++];
    }

    FREE(ends);
    FREE(starts);
  }

  Uintlist_free(&endlist);
  Uintlist_free(&startlist);

  return edges;
}
#endif


/* For IIT versions <= 2.  Previously sorted by Chrom_compare, but now
   we assume that chromosomes are represented by divs, which are
   pre-sorted by iit_store. */
#if 0
static int
string_compare (const void *x, const void *y) {
  char *a = (char *) x;
  char *b = (char *) y;

  return strcmp(a,b);
}

static int *
sort_matches_by_type (T this, int *matches, int nmatches, bool alphabetizep) {
  int *sorted;
  int type, index, i, j, k = 0, t;
  List_T *intervallists;
  Interval_T *intervals, interval;
  int *matches1, nmatches1, nintervals;
  char *typestring;
  char **strings;

  if (nmatches == 0) {
    return (int *) NULL;
  } else {
    sorted = (int *) CALLOC(nmatches,sizeof(int));
  }
  
  intervallists = (List_T *) CALLOC(this->ntypes,sizeof(List_T));
  for (i = 0; i < nmatches; i++) {
    index = matches[i];
    interval = &(this->intervals[0][index-1]);
    type = Interval_type(interval);
    intervallists[type] = List_push(intervallists[type],(void *) interval);
  }

  if (alphabetizep == true) {
    strings = (char **) CALLOC(this->ntypes,sizeof(char *));

    for (type = 0; type < this->ntypes; type++) {
      typestring = IIT_typestring(this,type);
      strings[type] = (char *) CALLOC(strlen(typestring)+1,sizeof(char));
      strcpy(strings[type],typestring);
    }
    qsort(strings,this->ntypes,sizeof(char *),string_compare);
  }

  for (t = 0; t < this->ntypes; t++) {
    if (alphabetizep == false) {
      type = t;
      typestring = IIT_typestring(this,type);
    } else {
      typestring = strings[t];
      type = IIT_typeint(this,typestring);
    }

    if ((nintervals = List_length(intervallists[type])) > 0) {
      intervals = (Interval_T *) List_to_array(intervallists[type],/*end*/NULL);
      qsort(intervals,nintervals,sizeof(Interval_T),Interval_cmp);

      i = 0;
      while (i < nintervals) {
	interval = intervals[i];
	matches1 = IIT_get_exact_multiple(&nmatches1,this,/*divstring*/NULL,Interval_low(interval),Interval_high(interval),type);
	if (matches1 != NULL) {
	  for (j = 0; j < nmatches1; j++) {
	    sorted[k++] = matches1[j];
	  }
	  i += nmatches1;
	  FREE(matches1);
	}
      }

      FREE(intervals);
      List_free(&(intervallists[type]));
    }

  }

  if (alphabetizep == true) {
    for (t = 0; t < this->ntypes; t++) {
      FREE(strings[t]);
    }
    FREE(strings);
  }

  FREE(intervallists);
  return sorted;
}
#endif


/* For IIT versions >= 3.  Assumes that matches are all in the same
   div */
static int *
sort_matches_by_position (T this, int *matches, int nmatches) {
  int *sorted, index, i;
  struct Interval_windex_T *intervals;

  if (nmatches == 0) {
    return (int *) NULL;
  } else {
    intervals = (struct Interval_windex_T *) CALLOC(nmatches,sizeof(struct Interval_windex_T));
    for (i = 0; i < nmatches; i++) {
      index = intervals[i].index = matches[i];
      intervals[i].interval = &(this->intervals[0][index-1]); /* Ignore divno here, because we have offset index */
    }
    qsort(intervals,nmatches,sizeof(struct Interval_windex_T),Interval_windex_cmp);

    sorted = (int *) CALLOC(nmatches,sizeof(int));
    for (i = 0; i < nmatches; i++) {
      sorted[i] = intervals[i].index;
    }

    FREE(intervals);
    return sorted;
  }
}




#if 0
/* Need to work on */
void
IIT_dump_counts (T this, bool alphabetizep) { 
  int type, divno, index, i, j, k, t;
  Interval_T interval;
  Uintlist_T *startlists, *endlists;
  int *matches, nmatches, nintervals;
  unsigned int *starts, *ends, edge;
  char *typestring;
  Chrom_T *chroms;

  startlists = (Uintlist_T *) CALLOC(this->ntypes,sizeof(Uintlist_T));
  endlists = (Uintlist_T *) CALLOC(this->ntypes,sizeof(Uintlist_T));
  for (i = 0; i < this->nintervals; i++) {
    interval = &(this->intervals[i]);
    type = Interval_type(interval);
    startlists[type] = Uintlist_push(startlists[type],Interval_low(interval));
    endlists[type] = Uintlist_push(endlists[type],Interval_high(interval));
  }

  if (alphabetizep == true) {
    chroms = (Chrom_T *) CALLOC(this->ntypes,sizeof(Chrom_T));

    for (type = 0; type < this->ntypes; type++) {
      typestring = IIT_typestring(this,type);
      chroms[type] = Chrom_from_string(typestring,/*mitochondrial_string*/NULL,/*order*/0U,/*circularp*/false);
    }
    qsort(chroms,this->ntypes,sizeof(Chrom_T),Chrom_compare);
  }

  for (t = 0; t < this->ntypes; t++) {
    if (alphabetizep == false) {
      type = t;
      typestring = IIT_typestring(this,type);
    } else {
      typestring = Chrom_string(chroms[t]); /* Not allocated; do not free */
      type = IIT_typeint(this,typestring);
    }

    if (Uintlist_length(startlists[type]) > 0) {
      starts = Uintlist_to_array(&nintervals,startlists[type]);
      ends = Uintlist_to_array(&nintervals,endlists[type]);
      qsort(starts,nintervals,sizeof(unsigned int),uint_cmp);
      qsort(ends,nintervals,sizeof(unsigned int),uint_cmp);

      i = j = 0;
      while (i < nintervals || j < nintervals) {
	if (i >= nintervals && j >= nintervals) {
	  /* done */
	  matches = (int *) NULL;
	} else if (i >= nintervals) {
	  /* work on remaining ends */
	  edge = ends[j++];
	  matches = IIT_get_typed(&nmatches,this,edge,edge,type,/*sortp*/false);
	  printf("%s\t%u\tend\t%d",typestring,edge,nmatches);
	  while (j < nintervals && ends[j] == edge) {
	    j++;
	  }
	} else if (j >= nintervals) {
	  /* work on remaining starts */
	  edge = starts[i++];
	  matches = IIT_get_typed(&nmatches,this,edge,edge,type,/*sortp*/false);
	  printf("%s\t%u\tstart\t%d",typestring,edge,nmatches);
	  while (i < nintervals && starts[i] == edge) {
	    i++;
	  }
	} else if (starts[i] <= ends[j]) {
	  edge = starts[i++];
	  matches = IIT_get_typed(&nmatches,this,edge,edge,type,/*sortp*/false);
	  printf("%s\t%u\tstart\t%d",typestring,edge,nmatches);
	  while (i < nintervals && starts[i] == edge) {
	    i++;
	  }
	} else {
	  edge = ends[j++];
	  matches = IIT_get_typed(&nmatches,this,edge,edge,type,/*sortp*/false);
	  printf("%s\t%u\tend\t%d",typestring,edge,nmatches);
	  while (j < nintervals && ends[j] == edge) {
	    j++;
	  }
	}

	if (matches != NULL) {
	  index = matches[0];
	  label = IIT_label(this,index,&allocp);
	  printf("\t%s",label);
	  if (allocp == true) {
	    FREE(label);
	  }

	  for (k = 1; k < nmatches; k++) {
	    index = matches[k];
	    label = IIT_label(this,index,&allocp);
	    printf(",%s",label);
	    if (allocp == true) {
	      FREE(label);
	    }
	  }
	  printf("\n");
	  FREE(matches);
	}
      }

      Uintlist_free(&(endlists[type]));
      Uintlist_free(&(startlists[type]));
      FREE(ends);
      FREE(starts);
    }

  }

  if (alphabetizep == true) {
    for (t = 0; t < this->ntypes; t++) {
      Chrom_free(&(chroms[t]));
    }
    FREE(chroms);
  }

  FREE(endlists);
  FREE(startlists);

  return;
}
#endif


/************************************************************************
 * For file format, see iit-write.c
 ************************************************************************/

void
IIT_free (T *old) {
  int divno;

  if (*old != NULL) {
    if ((*old)->name != NULL) {
      FREE((*old)->name);
    }

    if ((*old)->access == MMAPPED) {
#ifdef HAVE_MMAP
      munmap((void *) (*old)->annot_mmap,(*old)->annot_length);
      munmap((void *) (*old)->annotpointers_mmap,(*old)->annotpointers_length);
      munmap((void *) (*old)->label_mmap,(*old)->label_length);
      munmap((void *) (*old)->labelpointers_mmap,(*old)->labelpointers_length);
      munmap((void *) (*old)->labelorder_mmap,(*old)->labelorder_length);
      if ((*old)->valuep == true) {
	munmap((void *) (*old)->value_mmap,(*old)->value_length);
	munmap((void *) (*old)->valueorder_mmap,(*old)->valueorder_length);
      }
#endif
      close((*old)->fd);

    } else if ((*old)->access == FILEIO) {
      FREE((*old)->annotations);
#ifdef HAVE_64_BIT
      if ((*old)->annot_pointers_8p == true) {
	FREE((*old)->annotpointers8);
      } else {
	FREE((*old)->annotpointers);
      }
#else
      FREE((*old)->annotpointers);
#endif
      FREE((*old)->labels);
#ifdef HAVE_64_BIT
      if ((*old)->label_pointers_8p == true) {
	FREE((*old)->labelpointers8);
      } else {
	FREE((*old)->labelpointers);
      }
#else
      FREE((*old)->labelpointers);
#endif
      FREE((*old)->labelorder);
      /* close((*old)->fd); -- closed in read_annotations */

      if ((*old)->valuep == true) {
	FREE((*old)->values);
	FREE((*old)->valueorder);
      }

    } else if ((*old)->access == ALLOCATED_PRIVATE) {
      /* Nothing to close.  IIT must have been created by IIT_new. */

    } else if ((*old)->access == ALLOCATED_SHARED) {
      /* Nothing to close.  IIT must have been created by IIT_new. */

    } else {
      abort();
    }

    if ((*old)->fieldstrings != NULL) {
      FREE((*old)->fieldstrings);
    }
    FREE((*old)->fieldpointers);
    FREE((*old)->typestrings);
    FREE((*old)->typepointers);

    FREE((*old)->intervals[0]);
    FREE((*old)->intervals);

    for (divno = 0; divno < (*old)->ndivs; divno++) {
      /* Note: we are depending on Mem_free() to check that these are non-NULL */
      FREE((*old)->nodes[divno]);
      FREE((*old)->omegas[divno]);
      FREE((*old)->sigmas[divno]);
      if ((*old)->alphas != NULL) {
	FREE((*old)->betas[divno]);
	FREE((*old)->alphas[divno]);
      }
    }

    FREE((*old)->nodes);
    FREE((*old)->omegas);
    FREE((*old)->sigmas);
    if ((*old)->alphas != NULL) {
      FREE((*old)->betas);
      FREE((*old)->alphas);
    }

    FREE((*old)->divstrings);
    FREE((*old)->divpointers);
    FREE((*old)->cum_nnodes);
    FREE((*old)->nnodes);
    FREE((*old)->cum_nintervals);
    FREE((*old)->nintervals);

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
skip_trees (off_t offset, off_t filesize, FILE *fp, char *filename,
	    int skip_ndivs, int skip_nintervals, int skip_nnodes) {

  off_t skipsize;

  /* 4 is for alphas, betas, sigmas, and omegas */
  skipsize = (skip_nintervals + skip_ndivs) * 4 * sizeof(int);
  skipsize += skip_nnodes * sizeof(struct FNode_T);

  if ((offset += skipsize) > filesize) {
    fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after skip_trees %jd, filesize %jd).  Did you generate it using iit_store?\n",
	    filename,offset,filesize);
    exit(9);
  } else {
    move_relative(fp,skipsize);
  }

  return offset;
}



static off_t
read_tree (off_t offset, off_t filesize, FILE *fp, char *filename, T new, int divno) {
  size_t items_read;
  int i;

  if (new->version < 2) {
#if 0
    /* Computing only if needed */
    compute_flanking(new);
#else
    new->alphas[divno] = new->betas[divno] = (int *) NULL;
#endif

  } else {
    if ((offset += sizeof(int)*(new->nintervals[divno]+1)) > filesize) {
      fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after alphas %lld, filesize %lld).  Did you generate it using iit_store?\n",
	      filename,(long long int) offset,(long long int) filesize);
      exit(9);
    } else {
      new->alphas[divno] = (int *) CALLOC(new->nintervals[divno]+1,sizeof(int));
      if ((items_read = FREAD_INTS(new->alphas[divno],new->nintervals[divno]+1,fp)) != (unsigned int) new->nintervals[divno] + 1) {
	fprintf(stderr,"IIT file %s appears to be truncated.  items_read = %lld\n",
		filename,(long long int) items_read);
	exit(9);
      }
    }

    if ((offset += sizeof(int)*(new->nintervals[divno]+1)) > filesize) {
      fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after betas %lld, filesize %lld).  Did you generate it using iit_store?\n",
	      filename,(long long int) offset,(long long int) filesize);
      exit(9);
    } else {
      new->betas[divno] = (int *) CALLOC(new->nintervals[divno]+1,sizeof(int));
      if ((items_read = FREAD_INTS(new->betas[divno],new->nintervals[divno]+1,fp)) != (unsigned int) new->nintervals[divno] + 1) {
	fprintf(stderr,"IIT file %s appears to be truncated.  items_read = %zu\n",filename,items_read);
	exit(9);
      }
#if 0
      debug(
	    printf("betas[%d]:",divno);
	    for (i = 0; i < new->nintervals[divno]+1; i++) {
	      printf(" %d",new->betas[divno][i]);
	    }
	    printf("\n");
	    );
#endif
    }
  }

  if ((offset += sizeof(int)*(new->nintervals[divno]+1)) > filesize) {
    fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after sigmas %lld, filesize %lld).  Did you generate it using iit_store?\n",
	    filename,(long long int) offset,(long long int) filesize);
    exit(9);
  } else {
    new->sigmas[divno] = (int *) CALLOC(new->nintervals[divno]+1,sizeof(int));
    if ((items_read = FREAD_INTS(new->sigmas[divno],new->nintervals[divno]+1,fp)) != (unsigned int) new->nintervals[divno] + 1) {
      fprintf(stderr,"IIT file %s appears to be truncated\n",filename);
      exit(9);
    }
#if 0
    debug(
	  printf("sigmas[%d]:",divno);
	  for (i = 0; i < new->nintervals[divno]+1; i++) {
	    printf(" %d",new->sigmas[divno][i]);
	  }
	  printf("\n");
	  );
#endif
  }

  if ((offset += sizeof(int)*(new->nintervals[divno]+1)) > filesize) {
    fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after omegas %lld, filesize %lld).  Did you generate it using iit_store?\n",
	    filename,(long long int) offset,(long long int) filesize);
    exit(9);
  } else {
    new->omegas[divno] = (int *) CALLOC(new->nintervals[divno]+1,sizeof(int));
    if ((items_read = FREAD_INTS(new->omegas[divno],new->nintervals[divno]+1,fp)) != (unsigned int) new->nintervals[divno] + 1) {
      fprintf(stderr,"IIT file %s appears to be truncated\n",filename);
      exit(9);
    }
#if 0
    debug(
	  printf("omegas[%d]:",divno);
	  for (i = 0; i < new->nintervals[divno]+1; i++) {
	    printf(" %d",new->omegas[divno][i]);
	  }
	  printf("\n");
	  );
#endif
  }

  debug(printf("nnodes[%d]: %d\n",divno,new->nnodes[divno]));
  if (new->nnodes[divno] == 0) {
    new->nodes[divno] = (struct FNode_T *) NULL;
  } else {
    new->nodes[divno] = (struct FNode_T *) CALLOC(new->nnodes[divno],sizeof(struct FNode_T));
#ifdef WORDS_BIGENDIAN
    for (i = 0; i < new->nnodes[divno]; i++) {
      Bigendian_fread_uint(&(new->nodes[divno][i].value),fp);
      Bigendian_fread_int(&(new->nodes[divno][i].a),fp);
      Bigendian_fread_int(&(new->nodes[divno][i].b),fp);
      Bigendian_fread_int(&(new->nodes[divno][i].leftindex),fp);
      Bigendian_fread_int(&(new->nodes[divno][i].rightindex),fp);
    }
    offset += (sizeof(unsigned int)+sizeof(int)+sizeof(int)+sizeof(int)+sizeof(int))*new->nnodes[divno];
#else
    if (sizeof(struct FNode_T) == sizeof(unsigned int)+sizeof(int)+sizeof(int)+sizeof(int)+sizeof(int)) {
      offset += sizeof(struct FNode_T)*fread(new->nodes[divno],sizeof(struct FNode_T),new->nnodes[divno],fp);
    } else {
      for (i = 0; i < new->nnodes[divno]; i++) {
	fread(&(new->nodes[divno][i].value),sizeof(unsigned int),1,fp);
	fread(&(new->nodes[divno][i].a),sizeof(int),1,fp);
	fread(&(new->nodes[divno][i].b),sizeof(int),1,fp);
	fread(&(new->nodes[divno][i].leftindex),sizeof(int),1,fp);
	fread(&(new->nodes[divno][i].rightindex),sizeof(int),1,fp);
      }
      offset += (sizeof(unsigned int)+sizeof(int)+sizeof(int)+sizeof(int)+sizeof(int))*new->nnodes[divno];
    }
#endif
    if (offset > filesize) {
      fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after nodes %lld, filesize %lld).  Did you generate it using iit_store?\n",
	      filename,(long long int) offset,(long long int) filesize);
      exit(9);
    }

#if 1
    debug(
	  for (i = 0; i < new->nnodes[divno]; i++) {
	    printf("Read node %d %d %d\n",new->nodes[divno][i].value,new->nodes[divno][i].a,new->nodes[divno][i].b);
	  }
	  );
#endif

  }
  debug(printf("\n"));

  return offset;
}


static off_t
skip_intervals (int *skip_nintervals, off_t offset, off_t filesize, FILE *fp, char *filename, T new, 
		int divstart, int divend) {
  int divno;
  off_t skipsize = 0;

  *skip_nintervals = 0;
  for (divno = divstart; divno <= divend; divno++) {
    *skip_nintervals += new->nintervals[divno];
  }
  if (new->version >= 2) {
    skipsize += (sizeof(unsigned int)+sizeof(unsigned int)+sizeof(int)+sizeof(int))*(*skip_nintervals);
  } else {
    skipsize += (sizeof(unsigned int)+sizeof(unsigned int)+sizeof(int))*(*skip_nintervals);
  }

  if ((offset += skipsize) > filesize) {
    fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after skip_intervals %lld, filesize %lld).  Did you generate it using iit_store?\n",
	    filename,(long long int) offset,(long long int) filesize);
    exit(9);
  } else {
    move_relative(fp,skipsize);
  }

  return offset;
}


static off_t
read_intervals (off_t offset, off_t filesize, FILE *fp, char *filename, T new, int divno) {
  int i;

#ifdef WORDS_BIGENDIAN
  for (i = 0; i < new->nintervals[divno]; i++) {
    Bigendian_fread_uint(&(new->intervals[divno][i].low),fp);
    Bigendian_fread_uint(&(new->intervals[divno][i].high),fp);
    if (new->version >= 2) {
      Bigendian_fread_int(&(new->intervals[divno][i].sign),fp);
    } else {
      new->intervals[divno][i].sign = +1;
    }
    Bigendian_fread_int(&(new->intervals[divno][i].type),fp);
  }
  if (new->version >= 2) {
    offset += (sizeof(unsigned int)+sizeof(unsigned int)+sizeof(int)+sizeof(int))*new->nintervals[divno];
  } else {
    offset += (sizeof(unsigned int)+sizeof(unsigned int)+sizeof(int))*new->nintervals[divno];
  }
#else
  if (new->version >= 2 && sizeof(struct Interval_T) == sizeof(unsigned int)+sizeof(unsigned int)+sizeof(int)+sizeof(int)) {
    offset += sizeof(struct Interval_T)*fread(new->intervals[divno],sizeof(struct Interval_T),new->nintervals[divno],fp);
  } else if (new->version <= 1 && sizeof(struct Interval_T) == sizeof(unsigned int)+sizeof(unsigned int)+sizeof(int)) {
    offset += sizeof(struct Interval_T)*fread(new->intervals[divno],sizeof(struct Interval_T),new->nintervals[divno],fp);
  } else {
    for (i = 0; i < new->nintervals[divno]; i++) {
      fread(&(new->intervals[divno][i].low),sizeof(unsigned int),1,fp);
      fread(&(new->intervals[divno][i].high),sizeof(unsigned int),1,fp);
      if (new->version >= 2) {
	fread(&(new->intervals[divno][i].sign),sizeof(int),1,fp);
      } else {
	new->intervals[divno][i].sign = +1;
      }
      fread(&(new->intervals[divno][i].type),sizeof(int),1,fp);
    }
    if (new->version >= 2) {
      offset += (sizeof(unsigned int)+sizeof(unsigned int)+sizeof(int)+sizeof(int))*new->nintervals[divno];
    } else {
      offset += (sizeof(unsigned int)+sizeof(unsigned int)+sizeof(int))*new->nintervals[divno];
    }
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
#ifdef HAVE_64_BIT
  UINT8 length8;
#endif
  UINT4 length;
#ifdef DEBUG
  int i;
#endif

  new->typepointers = (unsigned int *) CALLOC(new->ntypes+1,sizeof(unsigned int));
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

  new->fieldpointers = (unsigned int *) CALLOC(new->nfields+1,sizeof(unsigned int));
  if (new->version < 2) {
    new->fieldpointers[0] = '\0';
  } else {
    offset += sizeof(int)*FREAD_UINTS(new->fieldpointers,new->nfields+1,fp);
  }
  stringlen = new->fieldpointers[new->nfields];
  if (stringlen == 0) {
    new->fieldstrings = (char *) NULL;
  } else {
    new->fieldstrings = (char *) CALLOC(stringlen,sizeof(char));
    offset += sizeof(char)*FREAD_CHARS(new->fieldstrings,stringlen,fp);
  }
  debug(
	printf("fieldstrings:\n");
	for (i = 0; i < stringlen; i++) {
	  printf("%c",new->fieldstrings[i]);
	}
	printf("\n");
	);

  if (new->valuep == true) {
    debug1(printf("Starting read of valueorder offset/length\n"));
    new->valueorder_offset = offset;
    new->valueorder_length = (size_t) (new->total_nintervals*sizeof(int));
    /* fprintf(stderr,"Doing a move_relative for valueorder_length %zu\n",new->valueorder_length); */
    move_relative(fp,new->valueorder_length);
    offset += new->valueorder_length;

    debug1(printf("Starting read of value offset/length\n"));
    new->value_offset = offset;
    new->value_length = (size_t) (new->total_nintervals*sizeof(double));
    /* fprintf(stderr,"Doing a move_relative for value_length %zu\n",new->value_length); */
    move_relative(fp,new->value_length);
    offset += new->value_length;
  }

  debug1(printf("Starting read of labelorder offset/length\n"));
  new->labelorder_offset = offset;
  new->labelorder_length = (size_t) (new->total_nintervals*sizeof(int));
  /* fprintf(stderr,"Doing a move_relative for labelorder_length %zu\n",new->labelorder_length); */
  move_relative(fp,new->labelorder_length);
  offset += new->labelorder_length;

  debug1(printf("Starting read of labelpointer offset/length\n"));
  new->labelpointers_offset = offset;
#ifdef HAVE_64_BIT
  if (new->label_pointers_8p == true) {
    new->labelpointers_length = (size_t) ((new->total_nintervals+1)*sizeof(UINT8));
    move_relative(fp,new->total_nintervals * sizeof(UINT8));
    FREAD_UINT8(&length8,fp);
    new->label_length = (size_t) length8;
  } else {
    new->labelpointers_length = (size_t) ((new->total_nintervals+1)*sizeof(UINT4));
    /* fprintf(stderr,"Doing a move_relative for labelpointer %zu\n",new->total_nintervals * sizeof(UINT4)); */
    move_relative(fp,new->total_nintervals * sizeof(UINT4));
    FREAD_UINT(&length,fp);
    new->label_length = (size_t) length;
  }
#else
  new->labelpointers_length = (size_t) ((new->total_nintervals+1)*sizeof(UINT4));
  /* fprintf(stderr,"Doing a move_relative for labelpointer %zu\n",new->total_nintervals * sizeof(UINT4)); */
  move_relative(fp,new->total_nintervals * sizeof(UINT4));
  FREAD_UINT(&length,fp);
  new->label_length = (size_t) length;
#endif
  offset += new->labelpointers_length;

  debug1(printf("Starting read of label offset/length\n"));
  new->label_offset = offset;
  /* new->label_length computed above */
  /* fprintf(stderr,"Doing a move_relative for label_length %zu\n",new->label_length); */
  move_relative(fp,new->label_length);
  offset += new->label_length;

  debug1(printf("Starting read of annotpointers offset/length\n"));
  new->annotpointers_offset = offset;
#ifdef HAVE_64_BIT
  if (new->annot_pointers_8p == true) {
    new->annotpointers_length = (size_t) ((new->total_nintervals+1)*sizeof(UINT8));
  } else {
    new->annotpointers_length = (size_t) ((new->total_nintervals+1)*sizeof(UINT4));
  }
#else
  new->annotpointers_length = (size_t) ((new->total_nintervals+1)*sizeof(UINT4));
#endif
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
#ifdef HAVE_64_BIT
  UINT8 length8;
#endif
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

  new->fieldpointers = (unsigned int *) CALLOC(new->nfields+1,sizeof(unsigned int));
  if (new->version < 2) {
    new->fieldpointers[0] = '\0';
  } else {
    offset += sizeof(int)*FREAD_UINTS(new->fieldpointers,new->nfields+1,fp);
  }
  stringlen = new->fieldpointers[new->nfields];
  if (stringlen == 0) {
    new->fieldstrings = (char *) NULL;
  } else {
    new->fieldstrings = (char *) CALLOC(stringlen,sizeof(char));
    offset += sizeof(char)*FREAD_CHARS(new->fieldstrings,stringlen,fp);
  }
  printf("fieldstrings:\n");
  for (i = 0; i < stringlen; i++) {
    printf("%c",new->fieldstrings[i]);
  }
  printf("\n");

  if (new->valuep == true) {
    debug1(printf("Starting read of valueorder offset/length\n"));
    new->valueorder_offset = offset;
    new->valueorder_length = (size_t) (new->total_nintervals*sizeof(int));
    /* fprintf(stderr,"Doing a move_relative for valueorder_length %zu\n",new->valueorder_length); */
    move_relative(fp,new->valueorder_length);
    offset += new->valueorder_length;

    debug1(printf("Starting read of value offset/length\n"));
    new->value_offset = offset;
    new->value_length = (size_t) (new->total_nintervals*sizeof(double));
    /* fprintf(stderr,"Doing a move_relative for value_length %zu\n",new->value_length); */
    move_relative(fp,new->value_length);
    offset += new->value_length;
  }

  debug1(printf("Starting read of labelorder offset/length\n"));
  new->labelorder_offset = offset;
  new->labelorder_length = (size_t) (new->total_nintervals*sizeof(int));
  move_relative(fp,new->labelorder_length);
  offset += new->labelorder_length;

  debug1(printf("Starting read of labelpointers offset/length\n"));
  new->labelpointers_offset = offset;
#ifdef HAVE_64_BIT
  if (new->label_pointers_8p == true) {
    new->labelpointers_length = (size_t) ((new->total_nintervals+1)*sizeof(UINT8));
    move_relative(fp,new->total_nintervals * sizeof(UINT8));
    FREAD_UINT8(&length8,fp);
    new->label_length = (size_t) length8;
  } else {
    new->labelpointers_length = (size_t) ((new->total_nintervals+1)*sizeof(UINT4));
    move_relative(fp,new->total_nintervals * sizeof(UINT4));
    FREAD_UINT(&length,fp);
    new->label_length = (size_t) length;
  }
#else
  new->labelpointers_length = (size_t) ((new->total_nintervals+1)*sizeof(UINT4));
  move_relative(fp,new->total_nintervals * sizeof(UINT4));
  FREAD_UINT(&length,fp);
  new->label_length = (size_t) length;
#endif
  offset += new->labelpointers_length;

  fprintf(stderr,"label_length: %zu\n",new->label_length);
  debug1(printf("Starting read of label offset/length\n"));
  new->label_offset = offset;
  /* new->label_length computed above */
  move_relative(fp,new->label_length);
  offset += new->label_length;

  debug1(printf("Starting read of annotpointers offset/length\n"));
  new->annotpointers_offset = offset;
#ifdef HAVE_64_BIT
  if (new->annot_pointers_8p == true) {
    new->annotpointers_length = (size_t) ((new->total_nintervals+1)*sizeof(UINT8));
  } else {
    new->annotpointers_length = (size_t) ((new->total_nintervals+1)*sizeof(UINT4));
  }
#else
  new->annotpointers_length = (size_t) ((new->total_nintervals+1)*sizeof(UINT4));
#endif
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

    if (new->valuep == true) {
      new->valueorder_mmap = Access_mmap_offset(&remainder,new->fd,new->valueorder_offset,new->valueorder_length,
						sizeof(char),/*randomp*/true);
      debug(fprintf(stderr,"valueorder_mmap is %p\n",new->valueorder_mmap));
      new->valueorder = (int *) &(new->valueorder_mmap[remainder]);
      new->valueorder_length += (size_t) remainder;
      
      new->value_mmap = Access_mmap_offset(&remainder,new->fd,new->value_offset,new->value_length,
					   sizeof(char),/*randomp*/true);
      debug(fprintf(stderr,"values_mmap is %p\n",new->value_mmap));
      new->values = (double *) &(new->value_mmap[remainder]);
      new->value_length += (size_t) remainder;
    }

    new->labelorder_mmap = Access_mmap_offset(&remainder,new->fd,new->labelorder_offset,new->labelorder_length,
					      sizeof(char),/*randomp*/true);
    debug(fprintf(stderr,"labelorder_mmap is %p\n",new->labelorder_mmap));
    new->labelorder = (int *) &(new->labelorder_mmap[remainder]);
    new->labelorder_length += (size_t) remainder;

    new->labelpointers_mmap = Access_mmap_offset(&remainder,new->fd,new->labelpointers_offset,new->labelpointers_length,
						 sizeof(char),/*randomp*/true);
    debug(fprintf(stderr,"labelpointers_mmap is %p\n",new->labelpointers_mmap));
#ifdef HAVE_64_BIT
    if (new->label_pointers_8p == true) {
      new->labelpointers8 = (UINT8 *) &(new->labelpointers_mmap[remainder]);
      new->labelpointers = (UINT4 *) NULL;
    } else {
      new->labelpointers8 = (UINT8 *) NULL;
      new->labelpointers = (UINT4 *) &(new->labelpointers_mmap[remainder]);
    }
#else
    new->labelpointers = (UINT4 *) &(new->labelpointers_mmap[remainder]);
#endif
    new->labelpointers_length += (size_t) remainder;

    new->label_mmap = Access_mmap_offset(&remainder,new->fd,new->label_offset,new->label_length,
					  sizeof(char),/*randomp*/true);
    debug(fprintf(stderr,"labels_mmap is %p\n",new->label_mmap));
    new->labels = (char *) &(new->label_mmap[remainder]);
    new->label_length += (size_t) remainder;

    new->annotpointers_mmap = Access_mmap_offset(&remainder,new->fd,new->annotpointers_offset,new->annotpointers_length,
						 sizeof(char),/*randomp*/true);
    debug(fprintf(stderr,"annotpointers_mmap is %p\n",new->annotpointers_mmap));
#ifdef HAVE_64_BIT
    if (new->annot_pointers_8p == true) {
      new->annotpointers8 = (UINT8 *) &(new->annotpointers_mmap[remainder]);
      new->annotpointers = (UINT4 *) NULL;
    } else {
      new->annotpointers8 = (UINT8 *) NULL;
      new->annotpointers = (UINT4 *) &(new->annotpointers_mmap[remainder]);
    }
#else
    new->annotpointers = (UINT4 *) &(new->annotpointers_mmap[remainder]);
#endif
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
#ifdef HAVE_64_BIT
    if (new->label_pointers_8p == true) {
      new->labelpointers8 = (UINT8 *) &(new->labelpointers_mmap[remainder]);
      new->labelpointers = (UINT4 *) NULL;
    } else {
      new->labelpointers8 = (UINT8 *) NULL;
      new->labelpointers = (UINT4 *) &(new->labelpointers_mmap[remainder]);
    }
#else
    new->labelpointers = (UINT4 *) &(new->labelpointers_mmap[remainder]);
#endif
    new->labelpointers_length += (size_t) remainder;

    new->label_mmap = Access_mmap_offset_rw(&remainder,new->fd,new->label_offset,new->label_length,
					    sizeof(char),/*randomp*/true);
    new->labels = (char *) &(new->label_mmap[remainder]);
    new->label_length += (size_t) remainder;

    new->annotpointers_mmap = Access_mmap_offset_rw(&remainder,new->fd,new->annotpointers_offset,new->annotpointers_length,
						    sizeof(char),/*randomp*/true);
#ifdef HAVE_64_BIT
    if (new->annot_pointers_8p == true) {
      new->annotpointers8 = (UINT8 *) &(new->annotpointers_mmap[remainder]);
      new->annotpointers = (UINT4 *) NULL;
    } else {
      new->annotpointers8 = (UINT8 *) NULL;
      new->annotpointers = (UINT4 *) &(new->annotpointers_mmap[remainder]);
    }
#else
    new->annotpointers = (UINT4 *) &(new->annotpointers_mmap[remainder]);
#endif
    new->annotpointers_length += (size_t) remainder;

    new->annot_mmap = Access_mmap_offset_rw(&remainder,new->fd,new->annot_offset,new->annot_length,
					    sizeof(char),/*randomp*/true);
    new->annotations = (char *) &(new->annot_mmap[remainder]);
    new->annot_length += (size_t) remainder;
  }

#ifdef HAVE_64_BIT
  if (new->label_pointers_8p == true) {
    if (new->labelorder == NULL || new->labelpointers8 == NULL || new->labels == NULL) {
      fprintf(stderr,"Memory mapping failed in reading IIT file %s.  Using slow file IO instead.\n",filename);
      return false;
    }
  } else {
    if (new->labelorder == NULL || new->labelpointers == NULL || new->labels == NULL) {
      fprintf(stderr,"Memory mapping failed in reading IIT file %s.  Using slow file IO instead.\n",filename);
      return false;
    }
  }
#else
  if (new->labelorder == NULL || new->labelpointers == NULL || new->labels == NULL) {
    fprintf(stderr,"Memory mapping failed in reading IIT file %s.  Using slow file IO instead.\n",filename);
    return false;
  }
#endif

#ifdef HAVE_64_BIT
  if (new->annot_pointers_8p == true) {
    if (new->annotpointers8 == NULL || new->annotations == NULL) {
      fprintf(stderr,"Memory mapping failed in reading IIT file %s.  Using slow file IO instead.\n",filename);
      return false;
    }
  } else {
    if (new->annotpointers == NULL || new->annotations == NULL) {
      fprintf(stderr,"Memory mapping failed in reading IIT file %s.  Using slow file IO instead.\n",filename);
      return false;
    }
  }
#else
  if (new->annotpointers == NULL || new->annotations == NULL) {
    fprintf(stderr,"Memory mapping failed in reading IIT file %s.  Using slow file IO instead.\n",filename);
    return false;
  }
#endif

  return true;
}
#endif


/* Used if access is FILEIO.  Subsequent accesses by bigendian
   machines to anything but (char *) will still need to convert. */
static void
read_annotations (T new) {

  if (new->valuep == true) {
    file_move_absolute(new->fd,new->valueorder_offset,sizeof(int),/*n*/0);
    new->valueorder = (int *) CALLOC(new->total_nintervals,sizeof(int));
    read(new->fd,new->valueorder,new->total_nintervals*sizeof(int));

    file_move_absolute(new->fd,new->value_offset,sizeof(char),/*n*/0);
    new->values = (double *) CALLOC(new->value_length,sizeof(char));
    read(new->fd,new->values,new->value_length*sizeof(char));
  }

  file_move_absolute(new->fd,new->labelorder_offset,sizeof(int),/*n*/0);
  new->labelorder = (int *) CALLOC(new->total_nintervals,sizeof(int));
  read(new->fd,new->labelorder,new->total_nintervals*sizeof(int));

#ifdef HAVE_64_BIT
  if (new->label_pointers_8p == true) {
    file_move_absolute(new->fd,new->labelpointers_offset,sizeof(UINT8),/*n*/0);
    new->labelpointers8 = (UINT8 *) CALLOC(new->total_nintervals+1,sizeof(UINT8));
    read(new->fd,new->labelpointers8,(new->total_nintervals+1)*sizeof(UINT8));
    new->labelpointers = (UINT4 *) NULL;
  } else {
    file_move_absolute(new->fd,new->labelpointers_offset,sizeof(UINT4),/*n*/0);
    new->labelpointers = (UINT4 *) CALLOC(new->total_nintervals+1,sizeof(UINT4));
    read(new->fd,new->labelpointers,(new->total_nintervals+1)*sizeof(UINT4));
    new->labelpointers8 = (UINT8 *) NULL;
  }
#else
  file_move_absolute(new->fd,new->labelpointers_offset,sizeof(UINT4),/*n*/0);
  new->labelpointers = (UINT4 *) CALLOC(new->total_nintervals+1,sizeof(UINT4));
  read(new->fd,new->labelpointers,(new->total_nintervals+1)*sizeof(UINT4));
#endif

  file_move_absolute(new->fd,new->label_offset,sizeof(char),/*n*/0);
  new->labels = (char *) CALLOC(new->label_length,sizeof(char));
  read(new->fd,new->labels,new->label_length*sizeof(char));

#ifdef HAVE_64_BIT
  if (new->annot_pointers_8p == true) {
    file_move_absolute(new->fd,new->annotpointers_offset,sizeof(UINT8),/*n*/0);
    new->annotpointers8 = (UINT8 *) CALLOC(new->total_nintervals+1,sizeof(UINT8));
    read(new->fd,new->annotpointers8,(new->total_nintervals+1)*sizeof(UINT8));
    new->annotpointers = (UINT4 *) NULL;
  } else {
    file_move_absolute(new->fd,new->annotpointers_offset,sizeof(UINT4),/*n*/0);
    new->annotpointers = (UINT4 *) CALLOC(new->total_nintervals+1,sizeof(UINT4));
    read(new->fd,new->annotpointers,(new->total_nintervals+1)*sizeof(UINT4));
    new->annotpointers8 = (UINT8 *) NULL;
  }
#else
  file_move_absolute(new->fd,new->annotpointers_offset,sizeof(UINT4),/*n*/0);
  new->annotpointers = (UINT4 *) CALLOC(new->total_nintervals+1,sizeof(UINT4));
  read(new->fd,new->annotpointers,(new->total_nintervals+1)*sizeof(UINT4));
#endif

  file_move_absolute(new->fd,new->annot_offset,sizeof(char),/*n*/0);
  new->annotations = (char *) CALLOC(new->annot_length,sizeof(char));
  read(new->fd,new->annotations,new->annot_length*sizeof(char));

  return;
}


int
IIT_read_divint (char *filename, char *divstring, bool add_iit_p) {
  char *newfile = NULL;
  FILE *fp;
  int version;
  off_t offset, filesize, skipsize;
  int total_nintervals, ntypes, nfields, divsort;
  int label_pointer_size, annot_pointer_size;

  int i, ndivs;
  UINT4 *divpointers, stringlen, start;
  char *divstrings;

  if (add_iit_p == true) {
    newfile = (char *) CALLOC(strlen(filename)+strlen(".iit")+1,sizeof(char));
    sprintf(newfile,"%s.iit",filename);
    if ((fp = FOPEN_READ_BINARY(newfile)) != NULL) {
      filename = newfile;
    } else if ((fp = FOPEN_READ_BINARY(filename)) == NULL) {
      /* fprintf(stderr,"Cannot open IIT file %s or %s\n",filename,newfile); */
      FREE(newfile);
      return -1;
    }
  } else if ((fp = FOPEN_READ_BINARY(filename)) == NULL) {
    /* fprintf(stderr,"Cannot open IIT file %s\n",filename); */
    return -1;
  }

  filesize = Access_filesize(filename);
  offset = 0U;

  if (FREAD_INT(&total_nintervals,fp) < 1) {
    fprintf(stderr,"IIT file %s appears to be empty\n",filename);
    fclose(fp);
    return -1;
  } else if ((offset += sizeof(int)) > filesize) {
    fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after first byte %lld, filesize %lld).  Did you generate it using iit_store?\n",
	    filename,(long long int) offset,(long long int) filesize);
    return -1;
  }

  if (total_nintervals > 0) {
    version = 1;

  } else {
    /* New format to indicate version > 1 */
    FREAD_INT(&version,fp);
    if (version > IIT_LATEST_VERSION_NOVALUES && version > IIT_LATEST_VERSION_VALUES) {
      fprintf(stderr,"This file is version %d, but this software can only read up to versions %d and %d\n",
	      version,IIT_LATEST_VERSION_NOVALUES,IIT_LATEST_VERSION_VALUES);
      return -1;
    } else if ((offset += sizeof(int)) > filesize) {
      fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after version %lld, filesize %lld).  Did you generate it using iit_store?\n",
	      filename,(long long int) offset,(long long int) filesize);
      return -1;
    }

    if (version < 5) {
    } else {
      /* Read new variables indicating sizes of label and annot pointers */
      if (FREAD_INT(&label_pointer_size,fp) < 1) {
	fprintf(stderr,"IIT file %s appears to be truncated\n",filename);
	return -1;
      } else if ((offset += sizeof(int)) > filesize) {
	fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after nintervals %lld, filesize %lld).  Did you generate it using iit_store?\n",
		filename,(long long int) offset,(long long int) filesize);
	return -1;
      }

      if (FREAD_INT(&annot_pointer_size,fp) < 1) {
	fprintf(stderr,"IIT file %s appears to be truncated\n",filename);
	return -1;
      } else if ((offset += sizeof(int)) > filesize) {
	fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after nintervals %lld, filesize %lld).  Did you generate it using iit_store?\n",
		filename,(long long int) offset,(long long int) filesize);
	return -1;
      }

      if (label_pointer_size == 4) {
      } else if (label_pointer_size == 8) {
      } else {
	fprintf(stderr,"IIT file %s has a problem with label_pointer_size being %d, expecting 4 or 8\n",
		filename,label_pointer_size);
	return -1;
      }

      if (annot_pointer_size == 4) {
      } else if (annot_pointer_size == 8) {
      } else {
	fprintf(stderr,"IIT file %s has a problem with annot_pointer_size being %d, expecting 4 or 8\n",
		filename,annot_pointer_size);
	return -1;
      }
    }

    /* Re-read total_nintervals */
    if (FREAD_INT(&total_nintervals,fp) < 1) {
      fprintf(stderr,"IIT file %s appears to be truncated\n",filename);
      return -1;
    } else if ((offset += sizeof(int)) > filesize) {
      fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after nintervals %lld, filesize %lld).  Did you generate it using iit_store?\n",
	      filename,(long long int) offset,(long long int) filesize);
      return -1;
    }
  }

  debug(printf("version: %d\n",version));
  debug(printf("total_nintervals: %d\n",total_nintervals));


  if (FREAD_INT(&ntypes,fp) < 1) {
    fprintf(stderr,"IIT file %s appears to be truncated\n",filename);
    return -1;
  } else if (ntypes < 0) {
    fprintf(stderr,"IIT file %s appears to have a negative number of types\n",filename);
    return -1;
  } else if ((offset += sizeof(int)) > filesize) {
    fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after ntypes %lld, filesize %lld).  Did you generate it using iit_store?\n",
	    filename,(long long int) offset,(long long int) filesize);
    return -1;
  }
  debug(printf("ntypes: %d\n",ntypes));


  if (version < 2) {
    nfields = 0;
  } else {
    if (FREAD_INT(&nfields,fp) < 1) {
      fprintf(stderr,"IIT file %s appears to be truncated\n",filename);
      return -1;
    } else if (nfields < 0) {
      fprintf(stderr,"IIT file %s appears to have a negative number of fields\n",filename);
      return -1;
    } else if ((offset += sizeof(int)) > filesize) {
      fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after nfields %lld, filesize %lld).  Did you generate it using iit_store?\n",
	      filename,(long long int) offset,(long long int) filesize);
      return -1;
    }
  }
  debug(printf("nfields: %d\n",nfields));


  if (version <= 2) {
    return -1;

  } else {

    if (FREAD_INT(&ndivs,fp) < 1) {
      fprintf(stderr,"IIT file %s appears to be truncated\n",filename);
      return -1;
    } else if (ndivs < 0) {
      fprintf(stderr,"IIT file %s appears to have a negative number of divs\n",filename);
      return -1;
    } else if ((offset += sizeof(int)) > filesize) {
      fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after ndivs %lld, filesize %lld).  Did you generate it using iit_store?\n",
	      filename,(long long int) offset,(long long int) filesize);
      return -1;
    }
    debug(printf("ndivs: %d\n",ndivs));

    /* Skip nintervals */
    offset += skipsize = sizeof(int)*ndivs;
    move_relative(fp,skipsize);

    /* Skip cum_nintervals */
    offset += skipsize = sizeof(int)*(ndivs+1);
    move_relative(fp,skipsize);

    /* Skip nnodes */
    offset += skipsize = sizeof(int)*ndivs;
    move_relative(fp,skipsize);

    /* Skip cum_nnodes */
    offset += skipsize = sizeof(int)*(ndivs+1);
    move_relative(fp,skipsize);

    if (FREAD_INT(&divsort,fp) < 1) {
      fprintf(stderr,"IIT file %s appears to be truncated\n",filename);
      return -1;
    } else if (divsort < 0) {
      fprintf(stderr,"IIT file %s appears to have a negative value for divsort\n",filename);
      return -1;
    } else if ((offset += sizeof(int)) > filesize) {
      fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after divsort %lld, filesize %lld).  Did you generate it using iit_store?\n",
	      filename,(long long int) offset,(long long int) filesize);
      return -1;
    }
    debug(printf("divsort: %d\n",divsort));

    divpointers = (UINT4 *) CALLOC(ndivs+1,sizeof(UINT4));
    offset += sizeof(int)*FREAD_UINTS(divpointers,ndivs+1,fp);
    debug(
	  printf("divpointers:");
	  for (i = 0; i < ndivs+1; i++) {
	    printf(" %u",divpointers[i]);
	  }
	  printf("\n");
	  );

    stringlen = divpointers[ndivs];
    if (stringlen == 0) {
      fprintf(stderr,"Problem with divstring stringlen being 0\n");
      exit(9);
    } else {
      divstrings = (char *) CALLOC(stringlen,sizeof(char));
    }
    offset += sizeof(char)*FREAD_CHARS(divstrings,stringlen,fp);
    debug(
	  printf("divstrings:\n");
	  for (i = 0; i < stringlen; i++) {
	    if (divstrings[i] == '\0') {
	      printf("\n");
	    } else {
	      printf("%c",divstrings[i]);
	    }
	  }
	  printf("(end of divstrings)\n");
	  );

    i = 0;
    while (i < ndivs) {
      start = divpointers[i];
      if (!strcmp(divstring,&(divstrings[start]))) {
	fclose(fp);
	FREE(divstrings);
	FREE(divpointers);
	if (newfile != NULL) {
	  FREE(newfile);
	}
	return i;
      }
      i++;
    }
    
    fclose(fp);
    FREE(divstrings);
    FREE(divpointers);
    if (newfile != NULL) {
      FREE(newfile);
    }
    return -1;
  }
}



T
IIT_read (char *filename, char *name, bool readonlyp, Divread_T divread, char *divstring,
	  bool add_iit_p, bool labels_read_p) {
  T new;
  FILE *fp;
  char *newfile = NULL;
  off_t offset = 0, filesize, stringlen;
  int skip_nintervals, desired_divno, divno;
  int label_pointer_size, annot_pointer_size;
#ifdef DEBUG
  int i;
  Interval_T interval;
#endif

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
    /* fprintf(stderr,"Cannot open IIT file %s\n",filename); */
    return NULL;
  }

  new = (T) MALLOC(sizeof(*new));

  filesize = Access_filesize(filename);

  if (name == NULL) {
    new->name = NULL;
  } else {
    new->name = (char *) CALLOC(strlen(name)+1,sizeof(char));
    strcpy(new->name,name);
  }

  if (FREAD_INT(&new->total_nintervals,fp) < 1) {
    fprintf(stderr,"IIT file %s appears to be empty\n",filename);
    fclose(fp);
    return NULL;
  } else if ((offset += sizeof(int)) > filesize) {
    fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after first byte %lld, filesize %lld).  Did you generate it using iit_store?\n",
	    filename,(long long int) offset,(long long int) filesize);
    return NULL;
  }

  if (new->total_nintervals != 0) {
    /* Need to use Univ_IIT_read instead */
    fprintf(stderr,"Unexpected error.  Using IIT_read code on a version 1 IIT\n");
    abort();

  } else {
    /* New format to indicate version > 1 */
    FREAD_INT(&new->version,fp);
    if (new->version > IIT_LATEST_VERSION_NOVALUES && new->version > IIT_LATEST_VERSION_VALUES) {
      fprintf(stderr,"This file is version %d, but this software can only read up to versions %d and %d\n",
	      new->version,IIT_LATEST_VERSION_NOVALUES,IIT_LATEST_VERSION_VALUES);
      return NULL;
    } else if ((offset += sizeof(int)) > filesize) {
      fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after version %lld, filesize %lld).  Did you generate it using iit_store?\n",
	      filename,(long long int) offset,(long long int) filesize);
      return NULL;
    }

    if (new->version == IIT_LATEST_VERSION_VALUES) {
      /* If IIT_LATEST_VERSION_VALUES increases, need to revise this code to handle version 6 */
      new->valuep = true;
    } else {
      new->valuep = false;
    }

    if (new->version <= 3) {
      new->label_pointers_8p = false;
      new->annot_pointers_8p = false;
    } else if (new->version == 4) {
      new->label_pointers_8p = true;
      new->annot_pointers_8p = true;
    } else {
      /* Read new variables indicating sizes of label and annot pointers */
      if (FREAD_INT(&label_pointer_size,fp) < 1) {
	fprintf(stderr,"IIT file %s appears to be truncated\n",filename);
	return NULL;
      } else if ((offset += sizeof(int)) > filesize) {
	fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after nintervals %lld, filesize %lld).  Did you generate it using iit_store?\n",
		filename,(long long int) offset,(long long int) filesize);
	return NULL;
      }

      if (FREAD_INT(&annot_pointer_size,fp) < 1) {
	fprintf(stderr,"IIT file %s appears to be truncated\n",filename);
	return NULL;
      } else if ((offset += sizeof(int)) > filesize) {
	fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after nintervals %lld, filesize %lld).  Did you generate it using iit_store?\n",
		filename,(long long int) offset,(long long int) filesize);
	return NULL;
      }

      if (label_pointer_size == 4) {
	new->label_pointers_8p = false;
      } else if (label_pointer_size == 8) {
	new->label_pointers_8p = true;
      } else {
	fprintf(stderr,"IIT file %s has a problem with label_pointer_size being %d, expecting 4 or 8\n",
		filename,label_pointer_size);
      }

      if (annot_pointer_size == 4) {
	new->annot_pointers_8p = false;
      } else if (annot_pointer_size == 8) {
	new->annot_pointers_8p = true;
      } else {
	fprintf(stderr,"IIT file %s has a problem with annot_pointer_size being %d, expecting 4 or 8\n",
		filename,annot_pointer_size);
      }
    }

    /* Re-read total_nintervals */
    if (FREAD_INT(&new->total_nintervals,fp) < 1) {
      fprintf(stderr,"IIT file %s appears to be truncated\n",filename);
      return NULL;
    } else if ((offset += sizeof(int)) > filesize) {
      fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after nintervals %lld, filesize %lld).  Did you generate it using iit_store?\n",
	      filename,(long long int) offset,(long long int) filesize);
      return NULL;
    }
  }

  debug(printf("version: %d\n",new->version));
  debug(printf("total_nintervals: %d\n",new->total_nintervals));


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


  if (new->version < 2) {
    new->nfields = 0;
  } else {
    if (FREAD_INT(&new->nfields,fp) < 1) {
      fprintf(stderr,"IIT file %s appears to be truncated\n",filename);
      return NULL;
    } else if (new->nfields < 0) {
      fprintf(stderr,"IIT file %s appears to have a negative number of fields\n",filename);
      return NULL;
    } else if ((offset += sizeof(int)) > filesize) {
      fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after nfields %lld, filesize %lld).  Did you generate it using iit_store?\n",
	      filename,(long long int) offset,(long long int) filesize);
      return NULL;
    }
  }
  debug(printf("nfields: %d\n",new->nfields));


  if (new->version <= 2) {
    new->ndivs = 1;

    new->nintervals = (int *) CALLOC(new->ndivs,sizeof(int));
    new->nintervals[0] = new->total_nintervals;
    new->cum_nintervals = (int *) CALLOC(new->ndivs+1,sizeof(int));
    new->cum_nintervals[0] = 0;
    new->cum_nintervals[1] = new->total_nintervals;

    new->nnodes = (int *) CALLOC(new->ndivs,sizeof(int));
    if (FREAD_INT(&(new->nnodes[0]),fp) < 1) {
      fprintf(stderr,"IIT file %s appears to be truncated\n",filename);
      return NULL;
    } else if (new->nnodes[0] < 0) {
      fprintf(stderr,"IIT file %s appears to have a negative number of nodes\n",filename);
      return NULL;
    } else if ((offset += sizeof(int)) > filesize) {
      fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after nnodes %lld, filesize %lld).  Did you generate it using iit_store?\n",
	      filename,(long long int) offset,(long long int) filesize);
      return NULL;
    }
    new->cum_nnodes = (int *) CALLOC(new->ndivs+1,sizeof(int));
    new->cum_nnodes[0] = 0;
    new->cum_nnodes[1] = new->nnodes[0];

    new->divsort = NO_SORT;

    new->divpointers = (UINT4 *) CALLOC(new->ndivs+1,sizeof(UINT4));
    new->divpointers[0] = 0;
    new->divpointers[1] = 1;

    new->divstrings = (char *) CALLOC(1,sizeof(char));
    new->divstrings[0] = '\0';

  } else {

    if (FREAD_INT(&new->ndivs,fp) < 1) {
      fprintf(stderr,"IIT file %s appears to be truncated\n",filename);
      return NULL;
    } else if (new->ndivs < 0) {
      fprintf(stderr,"IIT file %s appears to have a negative number of divs\n",filename);
      return NULL;
    } else if ((offset += sizeof(int)) > filesize) {
      fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after ndivs %lld, filesize %lld).  Did you generate it using iit_store?\n",
	      filename,(long long int) offset,(long long int) filesize);
      return NULL;
    }
    debug(printf("ndivs: %d\n",new->ndivs));

    new->nintervals = (int *) CALLOC(new->ndivs,sizeof(int));
    offset += sizeof(int)*FREAD_INTS(new->nintervals,new->ndivs,fp);
    debug(
	  printf("nintervals:");
	  for (i = 0; i < new->ndivs; i++) {
	    printf(" %d",new->nintervals[i]);
	  }
	  printf("\n");
	  );

    new->cum_nintervals = (int *) CALLOC(new->ndivs+1,sizeof(int));
    offset += sizeof(int)*FREAD_INTS(new->cum_nintervals,new->ndivs+1,fp);
    debug(
	  printf("cum_nintervals:");
	  for (i = 0; i <= new->ndivs; i++) {
	    printf(" %d",new->cum_nintervals[i]);
	  }
	  printf("\n");
	  );

    new->nnodes = (int *) CALLOC(new->ndivs,sizeof(int));
    offset += sizeof(int)*FREAD_INTS(new->nnodes,new->ndivs,fp);
    debug(
	  printf("nnodes:");
	  for (i = 0; i < new->ndivs; i++) {
	    printf(" %d",new->nnodes[i]);
	  }
	  printf("\n");
	  );

    new->cum_nnodes = (int *) CALLOC(new->ndivs+1,sizeof(int));
    offset += sizeof(int)*FREAD_INTS(new->cum_nnodes,new->ndivs+1,fp);
    debug(
	  printf("cum_nnodes:");
	  for (i = 0; i <= new->ndivs; i++) {
	    printf(" %d",new->cum_nnodes[i]);
	  }
	  printf("\n");
	  );

    if (FREAD_INT(&new->divsort,fp) < 1) {
      fprintf(stderr,"IIT file %s appears to be truncated\n",filename);
      return NULL;
    } else if (new->divsort < 0) {
      fprintf(stderr,"IIT file %s appears to have a negative value for divsort\n",filename);
      return NULL;
    } else if ((offset += sizeof(int)) > filesize) {
      fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after divsort %lld, filesize %lld).  Did you generate it using iit_store?\n",
	      filename,(long long int) offset,(long long int) filesize);
      return NULL;
    }
    debug(printf("divsort: %d\n",new->divsort));

    new->divpointers = (UINT4 *) CALLOC(new->ndivs+1,sizeof(UINT4));
    offset += sizeof(int)*FREAD_UINTS(new->divpointers,new->ndivs+1,fp);
    debug(
	  printf("divpointers:");
	  for (i = 0; i < new->ndivs+1; i++) {
	    printf(" %u",new->divpointers[i]);
	  }
	  printf("\n");
	  );

    stringlen = new->divpointers[new->ndivs];
    if (stringlen == 0) {
      new->divstrings = (char *) NULL;
    } else {
      new->divstrings = (char *) CALLOC(stringlen,sizeof(char));
      offset += sizeof(char)*FREAD_CHARS(new->divstrings,stringlen,fp);
    }
    debug(
	  printf("divstrings:\n");
	  for (i = 0; i < stringlen; i++) {
	    if (new->divstrings[i] == '\0') {
	      printf("\n");
	    } else {
	      printf("%c",new->divstrings[i]);
	    }
	  }
	  printf("(end of divstrings)\n");
	  );
  }

  new->alphas = (int **) CALLOC(new->ndivs,sizeof(int *));
  new->betas = (int **) CALLOC(new->ndivs,sizeof(int *));
  new->sigmas = (int **) CALLOC(new->ndivs,sizeof(int *));
  new->omegas = (int **) CALLOC(new->ndivs,sizeof(int *));
  new->nodes = (struct FNode_T **) CALLOC(new->ndivs,sizeof(struct FNode_T *));

  if (new->version == 1) {
    abort();
  }

  new->intervals = (struct Interval_T **) CALLOC(new->ndivs,sizeof(struct Interval_T *));

  if (divread == READ_ALL) {
    /* Read all divs */
    debug(printf("Reading all divs\n"));
    for (divno = 0; divno < new->ndivs; divno++) {
      debug(fprintf(stderr,"Starting read of div\n"));
      offset = read_tree(offset,filesize,fp,filename,new,divno);
      debug(fprintf(stderr,"Ending read of div\n"));
    }

    new->intervals[0] = (struct Interval_T *) CALLOC(new->total_nintervals,sizeof(struct Interval_T));
    offset = read_intervals(offset,filesize,fp,filename,new,/*divno*/0);
    for (divno = 1; divno < new->ndivs; divno++) {
      new->intervals[divno] = &(new->intervals[divno-1][new->nintervals[divno-1]]);
      offset = read_intervals(offset,filesize,fp,filename,new,divno);
    }

  } else if (divread == READ_NONE) {
    debug(printf("Reading no divs\n"));
    offset = skip_trees(offset,filesize,fp,filename,new->ndivs,
			new->cum_nintervals[new->ndivs],new->cum_nnodes[new->ndivs]);

    new->intervals[0] = (struct Interval_T *) CALLOC(new->total_nintervals,sizeof(struct Interval_T));
    offset = read_intervals(offset,filesize,fp,filename,new,/*divno*/0);
    for (divno = 1; divno < new->ndivs; divno++) {
      new->intervals[divno] = &(new->intervals[divno-1][new->nintervals[divno-1]]);
      offset = read_intervals(offset,filesize,fp,filename,new,divno);
    }

  } else if (divread == READ_ONE) {
    debug(printf("Reading only div %s\n",divstring));
    if ((desired_divno = IIT_divint(new,divstring)) < 0) {
      fprintf(stderr,"Cannot find div %s in IIT_read.  Ignoring div.\n",divstring);
      desired_divno = 0;
    }
    offset = skip_trees(offset,filesize,fp,filename,desired_divno,
			new->cum_nintervals[desired_divno],new->cum_nnodes[desired_divno]);
    debug1(fprintf(stderr,"Starting read of div\n"));
    offset = read_tree(offset,filesize,fp,filename,new,desired_divno);
    debug1(fprintf(stderr,"Ending read of div\n"));
    offset = skip_trees(offset,filesize,fp,filename,new->ndivs - (desired_divno + 1),
			new->cum_nintervals[new->ndivs] - new->cum_nintervals[desired_divno+1],
			new->cum_nnodes[new->ndivs] - new->cum_nnodes[desired_divno+1]);

    new->intervals[0] = (struct Interval_T *) CALLOC(new->total_nintervals,sizeof(struct Interval_T));
    offset = skip_intervals(&skip_nintervals,offset,filesize,fp,filename,new,0,desired_divno-1);
    debug1(fprintf(stderr,"Starting read of intervals\n"));
    new->intervals[desired_divno] = &(new->intervals[0][skip_nintervals]);
    offset = read_intervals(offset,filesize,fp,filename,new,desired_divno);
    debug1(fprintf(stderr,"Ending read of intervals\n"));
    offset = skip_intervals(&skip_nintervals,offset,filesize,fp,filename,new,desired_divno+1,new->ndivs-1);

    debug(
	  /*
	    printf("sigmas[%d]:\n",desired_divno);
	    for (i = 0; i < new->nintervals[desired_divno]+1; i++) {
	    interval = &(new->intervals[desired_divno][new->sigmas[desired_divno][i]]);
	    printf("%d %u..%u\n",new->sigmas[desired_divno][i],Interval_low(interval),Interval_high(interval));
	    }
	    printf("\n");
	  */

	  printf("alphas[%d]:\n",desired_divno);
	  for (i = 0; i < new->nintervals[desired_divno]+1; i++) {
	    interval = &(new->intervals[desired_divno][new->alphas[desired_divno][i]]);
	    printf("%d %u..%u\n",new->alphas[desired_divno][i],Interval_low(interval),Interval_high(interval));
	  }
	  printf("\n");
	  );


  } else {
    abort();
  }

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
IIT_debug (char *filename) {
  T new;
  FILE *fp;
  char *newfile = NULL;
  off_t offset = 0, filesize, stringlen;
  int skip_nintervals, desired_divno, divno, i;
  int label_pointer_size, annot_pointer_size;
  Divread_T divread = READ_ALL;
  char *divstring = NULL;
  bool add_iit_p = false;
#ifdef DEBUG
  Interval_T interval;
#endif

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

  new->name = NULL;

  if (FREAD_INT(&new->total_nintervals,fp) < 1) {
    fprintf(stderr,"IIT file %s appears to be empty\n",filename);
    fclose(fp);
    return;
  } else if ((offset += sizeof(int)) > filesize) {
    fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after first byte %lld, filesize %lld).  Did you generate it using iit_store?\n",
	    filename,(long long int) offset,(long long int) filesize);
    return;
  }

  if (new->total_nintervals > 0) {
    new->version = 1;
    new->valuep = false;
    
  } else {
    /* New format to indicate version > 1 */
    FREAD_INT(&new->version,fp);
    if (new->version > IIT_LATEST_VERSION_NOVALUES && new->version > IIT_LATEST_VERSION_VALUES) {
      fprintf(stderr,"This file is version %d, but this software can only read up to versions %d and %d\n",
	      new->version,IIT_LATEST_VERSION_NOVALUES,IIT_LATEST_VERSION_VALUES);
      return;
    } else if ((offset += sizeof(int)) > filesize) {
      fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after version %lld, filesize %lld).  Did you generate it using iit_store?\n",
	      filename,(long long int) offset,(long long int) filesize);
      return;
    }

    if (new->version == IIT_LATEST_VERSION_VALUES) {
      /* If IIT_LATEST_VERSION_VALUES increases, need to revise this code to handle version 6 */
      new->valuep = true;
    } else {
      new->valuep = false;
    }

    if (new->version <= 3) {
      new->label_pointers_8p = false;
      new->annot_pointers_8p = false;
    } else if (new->version == 4) {
      new->label_pointers_8p = true;
      new->annot_pointers_8p = true;
    } else {
      /* Read new variables indicating sizes of label and annot pointers */
      if (FREAD_INT(&label_pointer_size,fp) < 1) {
	fprintf(stderr,"IIT file %s appears to be truncated\n",filename);
	return;
      } else if ((offset += sizeof(int)) > filesize) {
	fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after nintervals %lld, filesize %lld).  Did you generate it using iit_store?\n",
		filename,(long long int) offset,(long long int) filesize);
	return;
      }

      if (FREAD_INT(&annot_pointer_size,fp) < 1) {
	fprintf(stderr,"IIT file %s appears to be truncated\n",filename);
	return;
      } else if ((offset += sizeof(int)) > filesize) {
	fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after nintervals %lld, filesize %lld).  Did you generate it using iit_store?\n",
		filename,(long long int) offset,(long long int) filesize);
	return;
      }

      if (label_pointer_size == 4) {
	new->label_pointers_8p = false;
      } else if (label_pointer_size == 8) {
	new->label_pointers_8p = true;
      } else {
	fprintf(stderr,"IIT file %s has a problem with label_pointer_size being %d, expecting 4 or 8\n",
		filename,label_pointer_size);
      }

      if (annot_pointer_size == 4) {
	new->annot_pointers_8p = false;
      } else if (annot_pointer_size == 8) {
	new->annot_pointers_8p = true;
      } else {
	fprintf(stderr,"IIT file %s has a problem with annot_pointer_size being %d, expecting 4 or 8\n",
		filename,annot_pointer_size);
      }
    }

    /* Re-read total_nintervals */
    if (FREAD_INT(&new->total_nintervals,fp) < 1) {
      fprintf(stderr,"IIT file %s appears to be truncated\n",filename);
      return;
    } else if ((offset += sizeof(int)) > filesize) {
      fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after nintervals %lld, filesize %lld).  Did you generate it using iit_store?\n",
	      filename,(long long int) offset,(long long int) filesize);
      return;
    }
  }
  if (new->total_nintervals < 0) {
    fprintf(stderr,"IIT file %s appears to have a negative number of intervals\n",filename);
    return;
  }

  printf("version: %d\n",new->version);
  printf("total_nintervals: %d\n",new->total_nintervals);

  if (new->version >= 5) {
    printf("label_pointer_size: %d\n",label_pointer_size);
    printf("annot_pointer_size: %d\n",annot_pointer_size);
  }


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


  if (new->version < 2) {
    new->nfields = 0;
  } else {
    if (FREAD_INT(&new->nfields,fp) < 1) {
      fprintf(stderr,"IIT file %s appears to be truncated\n",filename);
      return;
    } else if (new->nfields < 0) {
      fprintf(stderr,"IIT file %s appears to have a negative number of fields\n",filename);
      return;
    } else if ((offset += sizeof(int)) > filesize) {
      fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after nfields %lld, filesize %lld).  Did you generate it using iit_store?\n",
	      filename,(long long int) offset,(long long int) filesize);
      return;
    }
  }
  printf("nfields: %d\n",new->nfields);


  if (new->version <= 2) {
    new->ndivs = 1;

    new->nintervals = (int *) CALLOC(new->ndivs,sizeof(int));
    new->nintervals[0] = new->total_nintervals;
    new->cum_nintervals = (int *) CALLOC(new->ndivs+1,sizeof(int));
    new->cum_nintervals[0] = 0;
    new->cum_nintervals[1] = new->total_nintervals;

    new->nnodes = (int *) CALLOC(new->ndivs,sizeof(int));
    if (FREAD_INT(&(new->nnodes[0]),fp) < 1) {
      fprintf(stderr,"IIT file %s appears to be truncated\n",filename);
      return;
    } else if (new->nnodes[0] < 0) {
      fprintf(stderr,"IIT file %s appears to have a negative number of nodes\n",filename);
      return;
    } else if ((offset += sizeof(int)) > filesize) {
      fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after nnodes %lld, filesize %lld).  Did you generate it using iit_store?\n",
	      filename,(long long int) offset,(long long int) filesize);
      return;
    }
    new->cum_nnodes = (int *) CALLOC(new->ndivs+1,sizeof(int));
    new->cum_nnodes[0] = 0;
    new->cum_nnodes[1] = new->nnodes[0];

    new->divsort = NO_SORT;

    new->divpointers = (UINT4 *) CALLOC(new->ndivs+1,sizeof(UINT4));
    new->divpointers[0] = 0;
    new->divpointers[1] = 1;

    new->divstrings = (char *) CALLOC(1,sizeof(char));
    new->divstrings[0] = '\0';

  } else {

    if (FREAD_INT(&new->ndivs,fp) < 1) {
      fprintf(stderr,"IIT file %s appears to be truncated\n",filename);
      return;
    } else if (new->ndivs < 0) {
      fprintf(stderr,"IIT file %s appears to have a negative number of divs\n",filename);
      return;
    } else if ((offset += sizeof(int)) > filesize) {
      fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after ndivs %lld, filesize %lld).  Did you generate it using iit_store?\n",
	      filename,(long long int) offset,(long long int) filesize);
      return;
    }
    printf("ndivs: %d\n",new->ndivs);

    new->nintervals = (int *) CALLOC(new->ndivs,sizeof(int));
    offset += sizeof(int)*FREAD_INTS(new->nintervals,new->ndivs,fp);
    printf("nintervals:");
    for (i = 0; i < new->ndivs; i++) {
      printf(" %d",new->nintervals[i]);
    }
    printf("\n");

    new->cum_nintervals = (int *) CALLOC(new->ndivs+1,sizeof(int));
    offset += sizeof(int)*FREAD_INTS(new->cum_nintervals,new->ndivs+1,fp);
    printf("cum_nintervals:");
    for (i = 0; i <= new->ndivs; i++) {
      printf(" %d",new->cum_nintervals[i]);
    }
    printf("\n");

    new->nnodes = (int *) CALLOC(new->ndivs,sizeof(int));
    offset += sizeof(int)*FREAD_INTS(new->nnodes,new->ndivs,fp);
    printf("nnodes:");
    for (i = 0; i < new->ndivs; i++) {
      printf(" %d",new->nnodes[i]);
    }
    printf("\n");

    new->cum_nnodes = (int *) CALLOC(new->ndivs+1,sizeof(int));
    offset += sizeof(int)*FREAD_INTS(new->cum_nnodes,new->ndivs+1,fp);
    printf("cum_nnodes:");
    for (i = 0; i <= new->ndivs; i++) {
      printf(" %d",new->cum_nnodes[i]);
    }
    printf("\n");

    if (FREAD_INT(&new->divsort,fp) < 1) {
      fprintf(stderr,"IIT file %s appears to be truncated\n",filename);
      return;
    } else if (new->divsort < 0) {
      fprintf(stderr,"IIT file %s appears to have a negative value for divsort\n",filename);
      return;
    } else if ((offset += sizeof(int)) > filesize) {
      fprintf(stderr,"IIT file %s has an invalid binary format -- offset is too large (offset after divsort %lld, filesize %lld).  Did you generate it using iit_store?\n",
	      filename,(long long int) offset,(long long int) filesize);
      return;
    }
    printf("divsort: %d\n",new->divsort);

    new->divpointers = (UINT4 *) CALLOC(new->ndivs+1,sizeof(UINT4));
    offset += sizeof(int)*FREAD_UINTS(new->divpointers,new->ndivs+1,fp);
    printf("divpointers:");
    for (i = 0; i < new->ndivs+1; i++) {
      printf(" %u",new->divpointers[i]);
    }
    printf("\n");

    stringlen = new->divpointers[new->ndivs];
    if (stringlen == 0) {
      new->divstrings = (char *) NULL;
    } else {
      new->divstrings = (char *) CALLOC(stringlen,sizeof(char));
      offset += sizeof(char)*FREAD_CHARS(new->divstrings,stringlen,fp);
    }
    printf("divstrings:\n");
    for (i = 0; i < stringlen; i++) {
      if (new->divstrings[i] == '\0') {
	printf("\n");
      } else {
	printf("%c",new->divstrings[i]);
      }
    }
  }

  new->alphas = (int **) CALLOC(new->ndivs,sizeof(int *));
  new->betas = (int **) CALLOC(new->ndivs,sizeof(int *));
  new->sigmas = (int **) CALLOC(new->ndivs,sizeof(int *));
  new->omegas = (int **) CALLOC(new->ndivs,sizeof(int *));
  new->nodes = (struct FNode_T **) CALLOC(new->ndivs,sizeof(struct FNode_T *));
  new->intervals = (struct Interval_T **) CALLOC(new->ndivs,sizeof(struct Interval_T *));

  if (divread == READ_ALL) {
    /* Read all divs */
    debug(printf("Reading all divs\n"));
    for (divno = 0; divno < new->ndivs; divno++) {
      debug(printf("Div %d tree\n",divno));
      offset = read_tree(offset,filesize,fp,filename,new,divno);
    }

    debug(printf("Div 0 intervals\n"));
    new->intervals[0] = (struct Interval_T *) CALLOC(new->total_nintervals,sizeof(struct Interval_T));
    offset = read_intervals(offset,filesize,fp,filename,new,/*divno*/0);
    for (divno = 1; divno < new->ndivs; divno++) {
      debug(printf("Div %d intervals\n",divno));
      new->intervals[divno] = &(new->intervals[divno-1][new->nintervals[divno-1]]);
      offset = read_intervals(offset,filesize,fp,filename,new,divno);
    }

  } else if (divread == READ_NONE) {
    debug(printf("Reading no divs\n"));
    offset = skip_trees(offset,filesize,fp,filename,new->ndivs,
			new->cum_nintervals[new->ndivs],new->cum_nnodes[new->ndivs]);

    new->intervals[0] = (struct Interval_T *) CALLOC(new->total_nintervals,sizeof(struct Interval_T));
    offset = read_intervals(offset,filesize,fp,filename,new,/*divno*/0);
    for (divno = 1; divno < new->ndivs; divno++) {
      new->intervals[divno] = &(new->intervals[divno-1][new->nintervals[divno-1]]);
      offset = read_intervals(offset,filesize,fp,filename,new,divno);
    }

  } else if (divread == READ_ONE) {
    debug(printf("Reading only div %s\n",divstring));
    if ((desired_divno = IIT_divint(new,divstring)) < 0) {
      fprintf(stderr,"Cannot find div %s in IIT_read.  Ignoring div.\n",divstring);
      desired_divno = 0;
    }
    offset = skip_trees(offset,filesize,fp,filename,desired_divno,
			new->cum_nintervals[desired_divno],new->cum_nnodes[desired_divno]);
    debug1(fprintf(stderr,"Starting read of div\n"));
    offset = read_tree(offset,filesize,fp,filename,new,desired_divno);
    debug1(fprintf(stderr,"Ending read of div\n"));
    offset = skip_trees(offset,filesize,fp,filename,new->ndivs - (desired_divno + 1),
			new->cum_nintervals[new->ndivs] - new->cum_nintervals[desired_divno+1],
			new->cum_nnodes[new->ndivs] - new->cum_nnodes[desired_divno+1]);

    new->intervals[0] = (struct Interval_T *) CALLOC(new->total_nintervals,sizeof(struct Interval_T));
    offset = skip_intervals(&skip_nintervals,offset,filesize,fp,filename,new,0,desired_divno-1);
    debug1(fprintf(stderr,"Starting read of intervals\n"));
    new->intervals[desired_divno] = &(new->intervals[0][skip_nintervals]);
    offset = read_intervals(offset,filesize,fp,filename,new,desired_divno);
    debug1(fprintf(stderr,"Ending read of intervals\n"));
    offset = skip_intervals(&skip_nintervals,offset,filesize,fp,filename,new,desired_divno+1,new->ndivs-1);

    debug(
	  /*
	  printf("sigmas[%d]:\n",desired_divno);
	  for (i = 0; i < new->nintervals[desired_divno]+1; i++) {
	    interval = &(new->intervals[desired_divno][new->sigmas[desired_divno][i]]);
	    printf("%d %u..%u\n",new->sigmas[desired_divno][i],Interval_low(interval),Interval_high(interval));
	  }
	  printf("\n");
	  */

	  printf("alphas[%d]:\n",desired_divno);
	  for (i = 0; i < new->nintervals[desired_divno]+1; i++) {
	    interval = &(new->intervals[desired_divno][new->alphas[desired_divno][i]]);
	    printf("%d %u..%u\n",new->alphas[desired_divno][i],Interval_low(interval),Interval_high(interval));
	  }
	  printf("\n");
	  );


  } else {
    abort();
  }

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

  return;
}


/************************************************************************/

static void 
fnode_query_aux (int *min, int *max, T this, int divno, int nodeindex, Chrpos_T x) {
  int lambda;
  FNode_T node;

  if (nodeindex == -1) {
    return;
  }

  node = &(this->nodes[divno][nodeindex]);
  if (x == node->value) {
    debug(printf("%uD:\n",node->value));
    if (node->a < *min) {
      *min = node->a;
    }
    if (node->b > *max) {
      *max = node->b;
    }
    return;
  } else if (x < node->value) {
    fnode_query_aux(&(*min),&(*max),this,divno,node->leftindex,x);
    debug(printf("%uL:\n",node->value));
    if (node->a < *min) {
      *min = node->a;
    }
    for (lambda = node->a; lambda <= node->b; lambda++) {
      debug(printf("Looking at lambda %d, segment %d\n",
		   lambda,this->sigmas[divno][lambda]));
      if (Interval_is_contained(x,this->intervals[divno],this->sigmas[divno][lambda]) == true) {
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
    fnode_query_aux(&(*min),&(*max),this,divno,node->rightindex,x);
    debug(printf("%uR:\n", node->value));
    if (node->b > *max) {
      *max = node->b;
    }
    for (lambda = node->b; lambda >= node->a; lambda--) {
      debug(printf("Looking at lambda %d, segment %d\n",
		   lambda,this->omegas[divno][lambda]));
      if (Interval_is_contained(x,this->intervals[divno],this->omegas[divno][lambda]) == true) {
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
IIT_find (int *nmatches, T this, char *label) {
  int *matches = NULL, j;
  int low, middle, high, recno;
  bool foundp = false;
  int cmp;

  low = 0;
  high = this->total_nintervals;
  *nmatches = 0;

#ifdef DEBUG
#ifndef WORDS_BIGENDIAN
#ifdef HAVE_64_BIT
  for (middle = low; middle < high; middle++) {
    printf("%d:%d:%s, ",middle,this->labelorder[middle],
	   &(this->labels[this->labelpointers8[this->labelorder[middle]]]));
  }
  printf("\n");
#endif
#endif
#endif

  while (!foundp && low < high) {
    middle = (low+high)/2;
#ifdef DEBUG
#ifndef WORDS_BIGENDIAN
#ifdef HAVE_64_BIT
    printf("low %d:%d:%s. middle %d:%d:%s high %d:%d:%s\n",
	   low,this->labelorder[low],
	   &(this->labels[this->labelpointers8[this->labelorder[low]]]),
	   middle,this->labelorder[middle],
	   &(this->labels[this->labelpointers8[this->labelorder[middle]]]),
	   high,this->labelorder[high],
	   &(this->labels[this->labelpointers8[this->labelorder[high]]]));
#endif
#endif
#endif

#ifdef WORDS_BIGENDIAN
#ifdef HAVE_64_BIT
    if (this->label_pointers_8p == true) {
      cmp = strcmp(label,&(this->labels[Bigendian_convert_uint8(this->labelpointers8[Bigendian_convert_int(this->labelorder[middle])])]));
    } else {
      cmp = strcmp(label,&(this->labels[Bigendian_convert_uint(this->labelpointers[Bigendian_convert_int(this->labelorder[middle])])]));
    }
#else
    cmp = strcmp(label,&(this->labels[Bigendian_convert_uint(this->labelpointers[Bigendian_convert_int(this->labelorder[middle])])]));
#endif
#else
#ifdef HAVE_64_BIT
    if (this->label_pointers_8p == true) {
      cmp = strcmp(label,&(this->labels[this->labelpointers8[this->labelorder[middle]]]));
    } else {
      cmp = strcmp(label,&(this->labels[this->labelpointers[this->labelorder[middle]]]));
    }
#else
    cmp = strcmp(label,&(this->labels[this->labelpointers[this->labelorder[middle]]]));
#endif
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
#ifdef HAVE_64_BIT
    if (this->label_pointers_8p == true) {
      while (low-1 >= 0 && 
	     !strcmp(label,&(this->labels[Bigendian_convert_uint8(this->labelpointers8[Bigendian_convert_int(this->labelorder[low-1])])]))) {
	low--;
      }
    } else {
      while (low-1 >= 0 && 
	     !strcmp(label,&(this->labels[Bigendian_convert_uint(this->labelpointers[Bigendian_convert_int(this->labelorder[low-1])])]))) {
	low--;
      }
    }
#else
    while (low-1 >= 0 && 
	   !strcmp(label,&(this->labels[Bigendian_convert_uint(this->labelpointers[Bigendian_convert_int(this->labelorder[low-1])])]))) {
      low--;
    }
#endif
#else
#ifdef HAVE_64_BIT
    if (this->label_pointers_8p == true) {
      while (low-1 >= 0 && 
	     !strcmp(label,&(this->labels[this->labelpointers8[this->labelorder[low-1]]]))) {
	low--;
      }
    } else {
      while (low-1 >= 0 && 
	     !strcmp(label,&(this->labels[this->labelpointers[this->labelorder[low-1]]]))) {
	low--;
      }
    }
#else
    while (low-1 >= 0 && 
	   !strcmp(label,&(this->labels[this->labelpointers[this->labelorder[low-1]]]))) {
      low--;
    }
#endif
#endif
    
    high = middle;
#ifdef WORDS_BIGENDIAN
#ifdef HAVE_64_BIT
    if (this->label_pointers_8p == true) {
      while (high+1 < this->total_nintervals && 
	     !strcmp(label,&(this->labels[Bigendian_convert_uint8(this->labelpointers8[Bigendian_convert_int(this->labelorder[high+1])])]))) {
	high++;
      }
    } else {
      while (high+1 < this->total_nintervals && 
	     !strcmp(label,&(this->labels[Bigendian_convert_uint(this->labelpointers[Bigendian_convert_int(this->labelorder[high+1])])]))) {
	high++;
      }
    }
#else
    while (high+1 < this->total_nintervals && 
	   !strcmp(label,&(this->labels[Bigendian_convert_uint(this->labelpointers[Bigendian_convert_int(this->labelorder[high+1])])]))) {
      high++;
    }
#endif
#else
#ifdef HAVE_64_BIT
    if (this->label_pointers_8p == true) {
      while (high+1 < this->total_nintervals && 
	     !strcmp(label,&(this->labels[this->labelpointers8[this->labelorder[high+1]]]))) {
	high++;
      }
    } else {
      while (high+1 < this->total_nintervals && 
	     !strcmp(label,&(this->labels[this->labelpointers[this->labelorder[high+1]]]))) {
	high++;
      }
    }
#else
    while (high+1 < this->total_nintervals && 
	   !strcmp(label,&(this->labels[this->labelpointers[this->labelorder[high+1]]]))) {
      high++;
    }
#endif
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
IIT_find_linear (T this, char *label) {
  int i;
  char *p;

  for (i = 0; i < this->total_nintervals; i++) {
#ifdef WORDS_BIGENDIAN
#ifdef HAVE_64_BIT
    if (this->label_pointers_8p == true) {
      p = &(this->labels[Bigendian_convert_uint8(this->labelpointers8[i])]);
    } else {
      p = &(this->labels[Bigendian_convert_uint(this->labelpointers[i])]);
    }
#else
    p = &(this->labels[Bigendian_convert_uint(this->labelpointers[i])]);
#endif
#else
#ifdef HAVE_64_BIT
    if (this->label_pointers_8p == true) {
      p = &(this->labels[this->labelpointers8[i]]);
    } else {
      p = &(this->labels[this->labelpointers[i]]);
    }
#else
    p = &(this->labels[this->labelpointers[i]]);
#endif
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

int
IIT_find_one (T this, char *label) {
  int index;
  int *matches, nmatches;

  matches = IIT_find(&nmatches,this,label);
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


static int
uint_compare_ascending (const void *a, const void *b) {
  unsigned int x = * (unsigned int *) a;
  unsigned int y = * (unsigned int *) b;

  if (x < y) {
    return -1;
  } else if (y < x) {
    return +1;
  } else {
    return 0;
  }
}


static int
uint_compare_descending (const void *a, const void *b) {
  unsigned int x = * (unsigned int *) a;
  unsigned int y = * (unsigned int *) b;

  if (x > y) {
    return -1;
  } else if (y > x) {
    return +1;
  } else {
    return 0;
  }
}


Chrpos_T *
IIT_get_highs_for_low (int *nuniq, T this, int divno, Chrpos_T x) {
  Chrpos_T *uniq = NULL, *coords = NULL, prev;
  int neval, ncoords, i;
  int match, lambda, min1, max1 = 0;
  struct Interval_T interval;

  if (divno < 0) {
    /* fprintf(stderr,"No div %s found in iit file\n",divstring); */
    *nuniq = 0;
    return (Chrpos_T *) NULL;
  }
  min1 = this->nintervals[divno] + 1;

  debug(printf("Entering IIT_get_highs_for_low with divno %d and query %u\n",divno,x));
  fnode_query_aux(&min1,&max1,this,divno,0,x);
  debug(printf("min1=%d max1=%d\n",min1,max1));

  if (max1 < min1) {
    *nuniq = 0;
    return (Chrpos_T *) NULL;
  } else {
    neval = (max1 - min1 + 1) + (max1 - min1 + 1);
    coords = (Chrpos_T *) CALLOC(neval,sizeof(Chrpos_T));
    ncoords = 0;

    for (lambda = min1; lambda <= max1; lambda++) {
      match = this->sigmas[divno][lambda];
      /* Have to subtract 1 because intervals array is zero-based */
      interval = this->intervals[divno][match - 1];
      if (interval.low == x) {
	coords[ncoords++] = interval.high;
      }

      match = this->omegas[divno][lambda];
      /* Have to subtract 1 because intervals array is zero-based */
      interval = this->intervals[divno][match - 1];
      if (interval.low == x) {
	coords[ncoords++] = interval.high;
      }
    }

    if (ncoords == 0) {
      *nuniq = 0;
      FREE(coords);
      return (Chrpos_T *) NULL;

    } else {
      /* Eliminate duplicates */
      qsort(coords,ncoords,sizeof(Chrpos_T),uint_compare_ascending);

      uniq = (Chrpos_T *) CALLOC(ncoords,sizeof(Chrpos_T));
      *nuniq = 0;
      prev = 0;
      for (i = 0; i < ncoords; i++) {
	if (coords[i] != prev) {
	  uniq[(*nuniq)++] = coords[i];
	  prev = coords[i];
	}
      }
      
      FREE(coords);
      return uniq;
    }
  }
}


Chrpos_T *
IIT_get_lows_for_high (int *nuniq, T this, int divno, Chrpos_T x) {
  Chrpos_T *uniq = NULL, *coords = NULL, prev;
  int neval, ncoords, i;
  int match, lambda, min1, max1 = 0;
  struct Interval_T interval;

  if (divno < 0) {
    /* fprintf(stderr,"No div %s found in iit file\n",divstring); */
    *nuniq = 0;
    return (Chrpos_T *) NULL;
  }
  min1 = this->nintervals[divno] + 1;

  debug(printf("Entering IIT_get_lows_for_high with divno %d and query %u\n",divno,x));
  fnode_query_aux(&min1,&max1,this,divno,0,x);
  debug(printf("min1=%d max1=%d\n",min1,max1));

  if (max1 < min1) {
    *nuniq = 0;
    return (Chrpos_T *) NULL;
  } else {
    neval = (max1 - min1 + 1) + (max1 - min1 + 1);
    coords = (Chrpos_T *) CALLOC(neval,sizeof(Chrpos_T));
    ncoords = 0;

    for (lambda = min1; lambda <= max1; lambda++) {
      match = this->sigmas[divno][lambda];
      /* Have to subtract 1 because intervals array is zero-based */
      interval = this->intervals[divno][match - 1];
      if (interval.high == x) {
	coords[ncoords++] = interval.low;
      }

      match = this->omegas[divno][lambda];
      /* Have to subtract 1 because intervals array is zero-based */
      interval = this->intervals[divno][match - 1];
      if (interval.high == x) {
	coords[ncoords++] = interval.low;
      }
    }

    if (ncoords == 0) {
      *nuniq = 0;
      FREE(coords);
      return (Chrpos_T *) NULL;

    } else {
      /* Eliminate duplicates */
      qsort(coords,ncoords,sizeof(Chrpos_T),uint_compare_descending);

      uniq = (Chrpos_T *) CALLOC(ncoords,sizeof(Chrpos_T));
      *nuniq = 0;
      prev = 0;
      for (i = 0; i < ncoords; i++) {
	if (coords[i] != prev) {
	  uniq[(*nuniq)++] = coords[i];
	  prev = coords[i];
	}
      }
      
      FREE(coords);
      return uniq;
    }
  }
}


bool
IIT_low_exists_signed_p (T this, int divno, Chrpos_T x, int sign) {
  int match, lambda, min1, max1 = 0;
  struct Interval_T interval;

  if (divno < 0) {
    /* fprintf(stderr,"No div %s found in iit file\n",divstring); */
    return false;
  }
  min1 = this->nintervals[divno] + 1;

  debug(printf("Entering IIT_get_highs_for_low with divno %d and query %u\n",divno,x));
  fnode_query_aux(&min1,&max1,this,divno,0,x);
  debug(printf("min1=%d max1=%d\n",min1,max1));

  if (max1 < min1) {
    return false;
  } else {
    for (lambda = min1; lambda <= max1; lambda++) {
      match = this->sigmas[divno][lambda];
      /* Have to subtract 1 because intervals array is zero-based */
      interval = this->intervals[divno][match - 1];
      if (interval.low == x && (sign == 0 || interval.sign == sign)) {
	return true;
      }

      match = this->omegas[divno][lambda];
      /* Have to subtract 1 because intervals array is zero-based */
      interval = this->intervals[divno][match - 1];
      if (interval.low == x && (sign == 0 || interval.sign == sign)) {
	return true;
      }
    }

    return false;
  }
}

bool
IIT_high_exists_signed_p (T this, int divno, Chrpos_T x, int sign) {
  int match, lambda, min1, max1 = 0;
  struct Interval_T interval;

  if (divno < 0) {
    /* fprintf(stderr,"No div %s found in iit file\n",divstring); */
    return false;
  }
  min1 = this->nintervals[divno] + 1;

  debug(printf("Entering IIT_get_lows_for_high with divno %d and query %u\n",divno,x));
  fnode_query_aux(&min1,&max1,this,divno,0,x);
  debug(printf("min1=%d max1=%d\n",min1,max1));

  if (max1 < min1) {
    return false;
  } else {
    for (lambda = min1; lambda <= max1; lambda++) {
      match = this->sigmas[divno][lambda];
      /* Have to subtract 1 because intervals array is zero-based */
      interval = this->intervals[divno][match - 1];
      if (interval.high == x && (sign == 0 || interval.sign == sign)) {
	return true;
      }

      match = this->omegas[divno][lambda];
      /* Have to subtract 1 because intervals array is zero-based */
      interval = this->intervals[divno][match - 1];
      if (interval.high == x && (sign == 0 || interval.sign == sign)) {
	return true;
      }
    }

    return false;
  }
}


int *
IIT_get_lows_signed (int *nmatches, T this, int divno, Chrpos_T x, Chrpos_T y, int sign) {
  int *uniq = NULL, *matches, matchstart, neval, nfound, i;
  int match, lambda, prev;
  int min1, max1 = 0, min2, max2 = 0;
  struct Interval_T interval;

  if (divno < 0) {
    /* fprintf(stderr,"No div %s found in iit file\n",divstring); */
    *nmatches = 0;
    return (int *) NULL;
  } else {
    min1 = min2 = this->nintervals[divno] + 1;
  }

  debug(printf("Entering IIT_low_signed_p with divno %d and query %u..%u\n",divno,x,y));
  fnode_query_aux(&min1,&max1,this,divno,0,x);
  fnode_query_aux(&min2,&max2,this,divno,0,y);
  debug(printf("min1=%d max1=%d  min2=%d max2=%d\n",min1,max1,min2,max2));

  *nmatches = 0;
  if (max2 >= min1) {
    neval = (max2 - min1 + 1) + (max2 - min1 + 1);
    matches = (int *) CALLOC(neval,sizeof(int));

    nfound = 0;
    for (lambda = min1; lambda <= max2; lambda++) {
      match = this->sigmas[divno][lambda];
      /* Have to subtract 1 because intervals array is zero-based */
      interval = this->intervals[divno][match - 1];
      if (interval.low >= x && interval.low <= y && (sign == 0 || interval.sign == sign)) {
	matches[nfound++] = match;
      }

      match = this->omegas[divno][lambda];
      /* Have to subtract 1 because intervals array is zero-based */
      interval = this->intervals[divno][match - 1];
      if (interval.low >= x && interval.low <= y && (sign == 0 || interval.sign == sign)) {
	matches[nfound++] = match;
      }
    }

    if (nfound == 0) {
      return (int *) NULL;
    } else {
      /* Eliminate duplicates */
      uniq = (int *) CALLOC(nfound,sizeof(int));
      qsort(matches,nfound,sizeof(int),int_compare);
      prev = 0;
      debug(printf("unique segments in lambda %d to %d:",min1,max2));
      for (i = 0; i < nfound; i++) {
	if (matches[i] != prev) {
	  debug(printf(" %d",matches[i]));
	  uniq[(*nmatches)++] = matches[i];
	  prev = matches[i];
	}
      }
      debug(printf("\n"));

      /* No need to check for interval overlap */
    }
  }

  matchstart = this->cum_nintervals[divno];
  for (i = 0; i < *nmatches; i++) {
    uniq[i] += matchstart;
  }

  return uniq;
}


int *
IIT_get_highs_signed (int *nmatches, T this, int divno, Chrpos_T x, Chrpos_T y, int sign) {
  int *uniq = NULL, *matches, matchstart, neval, nfound, i;
  int match, lambda, prev;
  int min1, max1 = 0, min2, max2 = 0;
  struct Interval_T interval;

  if (divno < 0) {
    /* fprintf(stderr,"No div %s found in iit file\n",divstring); */
    *nmatches = 0;
    return (int *) NULL;
  } else {
    min1 = min2 = this->nintervals[divno] + 1;
  }

  debug(printf("Entering IIT_low_signed_p with divno %d and query %u..%u\n",divno,x,y));
  fnode_query_aux(&min1,&max1,this,divno,0,x);
  fnode_query_aux(&min2,&max2,this,divno,0,y);
  debug(printf("min1=%d max1=%d  min2=%d max2=%d\n",min1,max1,min2,max2));

  *nmatches = 0;
  if (max2 >= min1) {
    neval = (max2 - min1 + 1) + (max2 - min1 + 1);
    matches = (int *) CALLOC(neval,sizeof(int));

    nfound = 0;
    for (lambda = min1; lambda <= max2; lambda++) {
      match = this->sigmas[divno][lambda];
      /* Have to subtract 1 because intervals array is zero-based */
      interval = this->intervals[divno][match - 1];
      if (interval.high >= x && interval.high <= y && (sign == 0 || interval.sign == sign)) {
	matches[nfound++] = match;
      }

      match = this->omegas[divno][lambda];
      /* Have to subtract 1 because intervals array is zero-based */
      interval = this->intervals[divno][match - 1];
      if (interval.high >= x && interval.high <= y && (sign == 0 || interval.sign == sign)) {
	matches[nfound++] = match;
      }
    }

    if (nfound == 0) {
      return (int *) NULL;
    } else {
      /* Eliminate duplicates */
      uniq = (int *) CALLOC(nfound,sizeof(int));
      qsort(matches,nfound,sizeof(int),int_compare);
      prev = 0;
      debug(printf("unique segments in lambda %d to %d:",min1,max2));
      for (i = 0; i < nfound; i++) {
	if (matches[i] != prev) {
	  debug(printf(" %d",matches[i]));
	  uniq[(*nmatches)++] = matches[i];
	  prev = matches[i];
	}
      }
      debug(printf("\n"));

      /* No need to check for interval overlap */
    }
  }

  matchstart = this->cum_nintervals[divno];
  for (i = 0; i < *nmatches; i++) {
    uniq[i] += matchstart;
  }

  return uniq;
}



int *
IIT_get (int *nmatches, T this, char *divstring, Chrpos_T x, Chrpos_T y, bool sortp) {
  int *sorted, *matches = NULL, matchstart, *uniq, neval, nuniq, i;
  int lambda, prev;
  int divno;
  int min1, max1 = 0, min2, max2 = 0;
  int nintervals;

  divno = IIT_divint(this,divstring);
#if 1
  /* Usually don't need to check, unless crossing between iits,
     because divstring comes from same iit */
  if (divno < 0) {
    /* fprintf(stderr,"No div %s found in iit file\n",divstring); */
    *nmatches = 0;
    return (int *) NULL;
  }
#endif
  if ((nintervals = this->nintervals[divno]) == 0) {
    *nmatches = 0;
    return (int *) NULL;
  } else {
    min1 = min2 = nintervals + 1;
  }

  debug(printf("Entering IIT_get with query %u %u\n",x,y));
  fnode_query_aux(&min1,&max1,this,divno,0,x);
  fnode_query_aux(&min2,&max2,this,divno,0,y);
  debug(printf("min1=%d max1=%d  min2=%d max2=%d\n",min1,max1,min2,max2));

  *nmatches = 0;
  if (max2 >= min1) {
    neval = (max2 - min1 + 1) + (max2 - min1 + 1);
    matches = (int *) CALLOC(neval,sizeof(int));
    uniq = (int *) CALLOC(neval,sizeof(int));

    i = 0;
    for (lambda = min1; lambda <= max2; lambda++) {
      matches[i++] = this->sigmas[divno][lambda];
      matches[i++] = this->omegas[divno][lambda];
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
      if (Interval_overlap_p(x,y,this->intervals[divno],uniq[i]) == true) {
	matches[(*nmatches)++] = uniq[i];
	debug(printf("Pushing overlapping segment %d (%u..%u)\n",uniq[i],
		     Interval_low(&(this->intervals[divno][uniq[i]-1])),
		     Interval_high(&(this->intervals[divno][uniq[i]-1]))));
      } else {
	debug(printf("Not pushing non-overlapping segment %d (%u..%u)\n",uniq[i],
		     Interval_low(&(this->intervals[divno][uniq[i]-1])),
		     Interval_high(&(this->intervals[divno][uniq[i]-1]))));
      }
    }

    FREE(uniq);
  }

  /* Convert to universal indices */
  matchstart = this->cum_nintervals[divno];
  for (i = 0; i < *nmatches; i++) {
    matches[i] += matchstart;
  }

  if (sortp == false) {
    return matches;
#if 0
  } else if (this->version <= 2) {
    sorted = sort_matches_by_type(this,matches,*nmatches,/*alphabetizep*/true);
    FREE(matches);
    return sorted;
#endif
  } else {
    sorted = sort_matches_by_position(this,matches,*nmatches);
    FREE(matches);
    return sorted;
  }
}


int *
IIT_get_signed (int *nmatches, T this, char *divstring, Chrpos_T x, Chrpos_T y, int sign, bool sortp) {
  int *sorted, *matches = NULL, matchstart, *uniq, neval, nuniq, i;
  int lambda, prev;
  int divno;
  int min1, max1 = 0, min2, max2 = 0;
  int nintervals;
  int index;

  divno = IIT_divint(this,divstring);
#if 1
  /* Usually don't need to check, unless crossing between iits,
     because divstring comes from same iit */
  if (divno < 0) {
    /* fprintf(stderr,"No div %s found in iit file\n",divstring); */
    *nmatches = 0;
    return (int *) NULL;
  }
#endif
  if ((nintervals = this->nintervals[divno]) == 0) {
    *nmatches = 0;
    return (int *) NULL;
  } else {
    min1 = min2 = nintervals + 1;
  }

  debug(printf("Entering IIT_get with query %u %u\n",x,y));
  fnode_query_aux(&min1,&max1,this,divno,0,x);
  fnode_query_aux(&min2,&max2,this,divno,0,y);
  debug(printf("min1=%d max1=%d  min2=%d max2=%d\n",min1,max1,min2,max2));

  *nmatches = 0;
  if (max2 >= min1) {
    neval = (max2 - min1 + 1) + (max2 - min1 + 1);
    matches = (int *) CALLOC(neval,sizeof(int));
    uniq = (int *) CALLOC(neval,sizeof(int));

    i = 0;
    for (lambda = min1; lambda <= max2; lambda++) {
      index = this->sigmas[divno][lambda];
      if (sign == 0 || Interval_sign(&(this->intervals[divno][index-1])) == sign) {
	matches[i++] = index;
      }
      index = this->omegas[divno][lambda];
      if (sign == 0 || Interval_sign(&(this->intervals[divno][index-1])) == sign) {
	matches[i++] = index;
      }
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
      if (Interval_overlap_p(x,y,this->intervals[divno],uniq[i]) == true) {
	matches[(*nmatches)++] = uniq[i];
	debug(printf("Pushing overlapping segment %d (%u..%u)\n",uniq[i],
		     Interval_low(&(this->intervals[divno][uniq[i]-1])),
		     Interval_high(&(this->intervals[divno][uniq[i]-1]))));
      } else {
	debug(printf("Not pushing non-overlapping segment %d (%u..%u)\n",uniq[i],
		     Interval_low(&(this->intervals[divno][uniq[i]-1])),
		     Interval_high(&(this->intervals[divno][uniq[i]-1]))));
      }
    }

    FREE(uniq);
  }

  /* Convert to universal indices */
  matchstart = this->cum_nintervals[divno];
  for (i = 0; i < *nmatches; i++) {
    matches[i] += matchstart;
  }

  if (sortp == false) {
    return matches;
#if 0
  } else if (this->version <= 2) {
    sorted = sort_matches_by_type(this,matches,*nmatches,/*alphabetizep*/true);
    FREE(matches);
    return sorted;
#endif
  } else {
    sorted = sort_matches_by_position(this,matches,*nmatches);
    FREE(matches);
    return sorted;
  }
}


bool
IIT_exists_with_divno (T this, int divno, Chrpos_T x, Chrpos_T y) {
  int match;
  int lambda;
  int min1, max1 = 0, min2, max2 = 0;

  if (divno < 0) {
    /* fprintf(stderr,"No div %s found in iit file\n",divstring); */
    return false;
  }
  min1 = min2 = this->nintervals[divno] + 1;

  debug(printf("Entering IIT_get_with_divno with divno %d and query %u %u\n",divno,x,y));
  fnode_query_aux(&min1,&max1,this,divno,0,x);
  fnode_query_aux(&min2,&max2,this,divno,0,y);
  debug(printf("min1=%d max1=%d  min2=%d max2=%d\n",min1,max1,min2,max2));

  for (lambda = min1; lambda <= max2; lambda++) {
    match = this->sigmas[divno][lambda];
    if (Interval_overlap_p(x,y,this->intervals[divno],match) == true) {
      return true;
    }
    match = this->omegas[divno][lambda];
    if (Interval_overlap_p(x,y,this->intervals[divno],match) == true) {
      return true;
    }
  }

  return false;
}


bool
IIT_exists_with_divno_signed (T this, int divno, Chrpos_T x, Chrpos_T y, int sign) {
  int match;
  int lambda;
  int min1, max1 = 0, min2, max2 = 0;
  Interval_T interval;

  if (divno < 0) {
    /* fprintf(stderr,"No div %s found in iit file\n",divstring); */
    return false;
  }
  min1 = min2 = this->nintervals[divno] + 1;

  debug(printf("Entering IIT_exists_with_divno_signed with divno %d and query %u %u\n",divno,x,y));
  fnode_query_aux(&min1,&max1,this,divno,0,x);
  fnode_query_aux(&min2,&max2,this,divno,0,y);
  debug(printf("min1=%d max1=%d  min2=%d max2=%d\n",min1,max1,min2,max2));

  for (lambda = min1; lambda <= max2; lambda++) {
    match = this->sigmas[divno][lambda];
    interval = &(this->intervals[divno][match - 1]);
    if (Interval_low(interval) == x && Interval_high(interval) == y &&
	(sign == 0 || Interval_sign(interval) == sign)) {
      return true;
    }

    match = this->omegas[divno][lambda];
    interval = &(this->intervals[divno][match - 1]);
    if (Interval_low(interval) == x && Interval_high(interval) == y &&
	(sign == 0 || Interval_sign(interval) == sign)) {
      return true;
    }
  }

  return false;
}


bool
IIT_exists_with_divno_typed_signed (T this, int divno, Chrpos_T x, Chrpos_T y, int type, int sign) {
  int match;
  int lambda;
  int min1, max1 = 0, min2, max2 = 0;
  Interval_T interval;

  if (divno < 0) {
    /* fprintf(stderr,"No div %s found in iit file\n",divstring); */
    return false;
  }
  min1 = min2 = this->nintervals[divno] + 1;

  debug(printf("Entering IIT_exists_with_divno_typed_signed with divno %d and query %u %u\n",divno,x,y));
  fnode_query_aux(&min1,&max1,this,divno,0,x);
  fnode_query_aux(&min2,&max2,this,divno,0,y);
  debug(printf("min1=%d max1=%d  min2=%d max2=%d\n",min1,max1,min2,max2));

  for (lambda = min1; lambda <= max2; lambda++) {
    match = this->sigmas[divno][lambda];
    interval = &(this->intervals[divno][match - 1]);
    if (Interval_low(interval) == x && Interval_high(interval) == y &&
	Interval_type(interval) == type && (sign == 0 || Interval_sign(interval) == sign)) {
      return true;
    }

    match = this->omegas[divno][lambda];
    interval = &(this->intervals[divno][match - 1]);
    if (Interval_low(interval) == x && Interval_high(interval) == y &&
	Interval_type(interval) == type && (sign == 0 || Interval_sign(interval) == sign)) {
      return true;
    }
  }

  return false;
}


#if 0
bool
IIT_exists_with_divno_typed_signed (T this, int divno, Chrpos_T x, Chrpos_T y, int type, int sign) {
  int match;
  int lambda;
  int min1, max1 = 0, min2, max2 = 0;
  Interval_T interval;

  if (divno < 0) {
    /* fprintf(stderr,"No div %s found in iit file\n",divstring); */
    return false;
  }
  min1 = min2 = this->nintervals[divno] + 1;

  debug(printf("Entering IIT_get_with_divno with divno %d and query %u %u\n",divno,x,y));
  fnode_query_aux(&min1,&max1,this,divno,0,x);
  fnode_query_aux(&min2,&max2,this,divno,0,y);
  debug(printf("min1=%d max1=%d  min2=%d max2=%d\n",min1,max1,min2,max2));

  for (lambda = min1; lambda <= max2; lambda++) {
    match = this->sigmas[divno][lambda];
    interval = &(this->intervals[divno][match - 1]);
    if (Interval_overlap_p(x,y,this->intervals[divno],match) == true &&
	Interval_type(interval) == type && (sign == 0 || Interval_sign(interval) == sign)) {
      return true;
    }
    match = this->omegas[divno][lambda];
    interval = &(this->intervals[divno][match - 1]);
    if (Interval_overlap_p(x,y,this->intervals[divno],match) == true &&
	Interval_type(interval) == type && (sign == 0 || Interval_sign(interval) == sign)) {
      return true;
    }
  }

  return false;
}
#endif



int *
IIT_get_with_divno (int *nmatches, T this, int divno, Chrpos_T x, Chrpos_T y, bool sortp) {
  int *sorted, *matches = NULL, matchstart, *uniq, neval, nuniq, i;
  int lambda, prev;
  int min1, max1 = 0, min2, max2 = 0;

  if (divno < 0) {
    /* fprintf(stderr,"No div %s found in iit file\n",divstring); */
    *nmatches = 0;
    return (int *) NULL;
  }
  min1 = min2 = this->nintervals[divno] + 1;

  debug(printf("Entering IIT_get_with_divno with divno %d and query %u %u\n",divno,x,y));
  fnode_query_aux(&min1,&max1,this,divno,0,x);
  fnode_query_aux(&min2,&max2,this,divno,0,y);
  debug(printf("min1=%d max1=%d  min2=%d max2=%d\n",min1,max1,min2,max2));

  *nmatches = 0;
  if (max2 >= min1) {
    neval = (max2 - min1 + 1) + (max2 - min1 + 1);
    matches = (int *) CALLOC(neval,sizeof(int));
    uniq = (int *) CALLOC(neval,sizeof(int));

    i = 0;
    for (lambda = min1; lambda <= max2; lambda++) {
      matches[i++] = this->sigmas[divno][lambda];
      matches[i++] = this->omegas[divno][lambda];
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
      if (Interval_overlap_p(x,y,this->intervals[divno],uniq[i]) == true) {
	matches[(*nmatches)++] = uniq[i];
	debug(printf("Pushing overlapping segment %d (%u..%u)\n",uniq[i],
		     Interval_low(&(this->intervals[divno][uniq[i]-1])),
		     Interval_high(&(this->intervals[divno][uniq[i]-1]))));
      } else {
	debug(printf("Not pushing non-overlapping segment %d (%u..%u)\n",uniq[i],
		     Interval_low(&(this->intervals[divno][uniq[i]-1])),
		     Interval_high(&(this->intervals[divno][uniq[i]-1]))));
      }
    }

    FREE(uniq);
  }

  /* Convert to universal indices */
  matchstart = this->cum_nintervals[divno];
  for (i = 0; i < *nmatches; i++) {
    matches[i] += matchstart;
  }

  if (sortp == false) {
    return matches;
#if 0
  } else if (this->version <= 2) {
    sorted = sort_matches_by_type(this,matches,*nmatches,/*alphabetizep*/true);
    FREE(matches);
    return sorted;
#endif
  } else {
    sorted = sort_matches_by_position(this,matches,*nmatches);
    FREE(matches);
    return sorted;
  }
}



int *
IIT_get_signed_with_divno (int *nmatches, T this, int divno, Chrpos_T x, Chrpos_T y, bool sortp,
			   int sign) {
  int *sorted, *matches = NULL, matchstart, *uniq, neval, nuniq, i;
  int lambda, prev;
  int min1, max1 = 0, min2, max2 = 0;
  int index;

  if (divno < 0) {
    /* fprintf(stderr,"No div %s found in iit file\n",divstring); */
    *nmatches = 0;
    return (int *) NULL;
  }
  min1 = min2 = this->nintervals[divno] + 1;

  debug(printf("Entering IIT_get_with_divno with divno %d and query %u %u\n",divno,x,y));
  fnode_query_aux(&min1,&max1,this,divno,0,x);
  fnode_query_aux(&min2,&max2,this,divno,0,y);
  debug(printf("min1=%d max1=%d  min2=%d max2=%d\n",min1,max1,min2,max2));

  *nmatches = 0;
  if (max2 >= min1) {
    neval = (max2 - min1 + 1) + (max2 - min1 + 1);
    matches = (int *) CALLOC(neval,sizeof(int));
    uniq = (int *) CALLOC(neval,sizeof(int));

    i = 0;
    for (lambda = min1; lambda <= max2; lambda++) {
      index = this->sigmas[divno][lambda];
      if (sign == 0 || Interval_sign(&(this->intervals[divno][index-1])) == sign) {
	matches[i++] = index;
      }
      index = this->omegas[divno][lambda];
      if (sign == 0 || Interval_sign(&(this->intervals[divno][index-1])) == sign) {
	matches[i++] = index;
      }
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
      if (Interval_overlap_p(x,y,this->intervals[divno],uniq[i]) == true) {
	matches[(*nmatches)++] = uniq[i];
	debug(printf("Pushing overlapping segment %d (%u..%u)\n",uniq[i],
		     Interval_low(&(this->intervals[divno][uniq[i]-1])),
		     Interval_high(&(this->intervals[divno][uniq[i]-1]))));
      } else {
	debug(printf("Not pushing non-overlapping segment %d (%u..%u)\n",uniq[i],
		     Interval_low(&(this->intervals[divno][uniq[i]-1])),
		     Interval_high(&(this->intervals[divno][uniq[i]-1]))));
      }
    }

    FREE(uniq);
  }

  /* Convert to universal indices */
  matchstart = this->cum_nintervals[divno];
  for (i = 0; i < *nmatches; i++) {
    matches[i] += matchstart;
  }

  if (sortp == false) {
    return matches;
#if 0
  } else if (this->version <= 2) {
    sorted = sort_matches_by_type(this,matches,*nmatches,/*alphabetizep*/true);
    FREE(matches);
    return sorted;
#endif
  } else {
    sorted = sort_matches_by_position(this,matches,*nmatches);
    FREE(matches);
    return sorted;
  }
}


static int
coord_search_low (T this, int divno, Chrpos_T x) {
  int low, middle, high;
  bool foundp = false;
  Chrpos_T middlevalue;
  int index;

  low = 1;			/* not 0, because alphas[divno][0] not used */
  high = this->nintervals[divno];

  debug3(printf("low = %d, high = %d\n",low,high));
  while (!foundp && low < high) {
    middle = (low+high)/2;
    index = this->alphas[divno][middle];
    middlevalue = Interval_low(&(this->intervals[divno][index-1]));

    debug3(printf("  compare x %u with middlevalue %u (for interval %d)\n",x,middlevalue,this->alphas[divno][middle]-1));
    if (x < middlevalue) {
      high = middle;
    } else if (x > middlevalue) {
      low = middle + 1;
    } else {
      foundp = true;
    }
    debug3(printf("low = %d, high = %d, middle = %d\n",low,high,middle));
  }

  if (foundp == true) {
    debug3(printf("found\n"));
    return middle;
  } else {
    debug3(printf("not found\n"));
    return low;
  }
}

static int
coord_search_high (T this, int divno, Chrpos_T x) {
  int low, middle, high;
  bool foundp = false;
  Chrpos_T middlevalue;
  int index;

  low = 1;			/* not 0, because betas[divno][0] not used */
  high = this->nintervals[divno];

  while (!foundp && low < high) {
    middle = (low+high)/2;
    index = this->betas[divno][middle];
    middlevalue = Interval_high(&(this->intervals[divno][index-1]));

    if (x < middlevalue) {
      high = middle;
    } else if (x > middlevalue) {
      low = middle + 1;
    } else {
      foundp = true;
    }
  }

  if (foundp == true) {
    return middle;
  } else {
    return high;
  }
}


void
IIT_get_flanking (int **leftflanks, int *nleftflanks, int **rightflanks, int *nrightflanks,
		  T this, char *divstring, Chrpos_T x, Chrpos_T y, int nflanking, int sign) {
  int lambda, matchstart, i;
  Interval_T interval;
  bool stopp;
  int divno;

  divno = IIT_divint(this,divstring);

  debug2(printf("Entering IIT_get_flanking with divno %d, query %u %u, nflanking = %d, sign %d\n",divno,x,y,nflanking,sign));

  if (this->alphas[divno] == NULL) {
#if 0
    compute_flanking(this);
#else
    fprintf(stderr,"Flanking hits not supported on version %d of iit files.  Please use iit_update to update your file\n",
	    this->version);
    exit(9);
#endif
  }

  /* Look at alphas for right flank */
  lambda = coord_search_low(this,divno,y);
  debug2(printf("coord_search_low lambda = %d\n",lambda));

  *rightflanks = (int *) CALLOC(nflanking,sizeof(int));
  *nrightflanks = 0;
  stopp = false;
  while (lambda <= this->nintervals[divno] && stopp == false) {
    interval = &(this->intervals[divno][this->alphas[divno][lambda]-1]);
    if (Interval_low(interval) <= y) {
      debug2(printf("Advancing because interval_low %u <= %u\n",Interval_low(interval),y));
      lambda++;
    } else if (sign != 0 && Interval_sign(interval) != sign) {
      debug2(printf("Advancing because sign != 0 && interval_sign %d != %d\n",Interval_sign(interval),sign));
      lambda++;
    } else {
      (*rightflanks)[(*nrightflanks)++] = this->alphas[divno][lambda];
      debug2(printf("Storing right flank %d\n",this->alphas[divno][lambda]));
      if (*nrightflanks < nflanking) {
	debug2(printf("Advancing because need more\n"));
	lambda++;
      } else {
	stopp = true;
      }
    }
  }

  /* Look at betas for left flank */
  lambda = coord_search_high(this,divno,x);

  *leftflanks = (int *) CALLOC(nflanking,sizeof(int));
  *nleftflanks = 0;
  stopp = false;
  while (lambda >= 1 && stopp == false) {
    interval = &(this->intervals[divno][this->betas[divno][lambda]-1]);
    if (Interval_high(interval) >= x) {
      lambda--;
    } else if (sign != 0 && Interval_sign(interval) != sign) {
      lambda--;
    } else {
      (*leftflanks)[(*nleftflanks)++] = this->betas[divno][lambda];
      if (*nleftflanks < nflanking) {
	lambda--;
      } else {
	stopp = true;
      }
    }
  }

  /* Convert to universal indices */
  matchstart = this->cum_nintervals[divno];
  for (i = 0; i < *nrightflanks; i++) {
    (*rightflanks)[i] += matchstart;
  }
  for (i = 0; i < *nleftflanks; i++) {
    (*leftflanks)[i] += matchstart;
  }

  return;
}

void
IIT_get_flanking_with_divno (int **leftflanks, int *nleftflanks, int **rightflanks, int *nrightflanks,
			     T this, int divno, Chrpos_T x, Chrpos_T y, int nflanking, int sign) {
  int lambda, matchstart, i;
  Interval_T interval;
  bool stopp;

  debug2(printf("Entering IIT_get_flanking_with_divno with divno %d, query %u %u, nflanking = %d, sign %d\n",divno,x,y,nflanking,sign));

  if (this->alphas[divno] == NULL) {
#if 0
    compute_flanking(this);
#else
    fprintf(stderr,"Flanking hits not supported on version %d of iit files.  Please use iit_update to update your file\n",
	    this->version);
    exit(9);
#endif
  }

  /* Look at alphas for right flank */
  lambda = coord_search_low(this,divno,y);
  debug2(printf("coord_search_low lambda = %d\n",lambda));

  *rightflanks = (int *) CALLOC(nflanking,sizeof(int));
  *nrightflanks = 0;
  stopp = false;
  while (lambda <= this->nintervals[divno] && stopp == false) {
    interval = &(this->intervals[divno][this->alphas[divno][lambda]-1]);
    if (Interval_low(interval) <= y) {
      debug2(printf("Advancing because interval_low %u <= %u\n",Interval_low(interval),y));
      lambda++;
    } else if (sign != 0 && Interval_sign(interval) != sign) {
      debug2(printf("Advancing because sign != 0 && interval_sign %d != %d\n",Interval_sign(interval),sign));
      lambda++;
    } else {
      (*rightflanks)[(*nrightflanks)++] = this->alphas[divno][lambda];
      debug2(printf("Storing right flank %d\n",this->alphas[divno][lambda]));
      if (*nrightflanks < nflanking) {
	debug2(printf("Advancing because need more\n"));
	lambda++;
      } else {
	stopp = true;
      }
    }
  }

  /* Look at betas for left flank */
  lambda = coord_search_high(this,divno,x);

  *leftflanks = (int *) CALLOC(nflanking,sizeof(int));
  *nleftflanks = 0;
  stopp = false;
  while (lambda >= 1 && stopp == false) {
    interval = &(this->intervals[divno][this->betas[divno][lambda]-1]);
    if (Interval_high(interval) >= x) {
      lambda--;
    } else if (sign != 0 && Interval_sign(interval) != sign) {
      lambda--;
    } else {
      (*leftflanks)[(*nleftflanks)++] = this->betas[divno][lambda];
      if (*nleftflanks < nflanking) {
	lambda--;
      } else {
	stopp = true;
      }
    }
  }

  /* Convert to universal indices */
  matchstart = this->cum_nintervals[divno];
  for (i = 0; i < *nrightflanks; i++) {
    (*rightflanks)[i] += matchstart;
  }
  for (i = 0; i < *nleftflanks; i++) {
    (*leftflanks)[i] += matchstart;
  }

  return;
}

void
IIT_get_flanking_typed (int **leftflanks, int *nleftflanks, int **rightflanks, int *nrightflanks,
			T this, char *divstring, Chrpos_T x, Chrpos_T y, int nflanking, int type,
			int sign) {
  int lambda, matchstart, i;
  Interval_T interval;
  bool stopp;
  int divno;

  divno = IIT_divint(this,divstring);

  debug2(printf("Entering IIT_get_flanking_typed with query %u %u => divno is %d\n",x,y,divno));

  if (this->alphas[divno] == NULL) {
#if 0
    IIT_compute_flanking(this);
#else
    fprintf(stderr,"Flanking hits not supported on version %d of iit files.  Please use iit_update to update your file\n",
	    this->version);
    exit(9);
#endif
  }

  /* Look at alphas for right flank */
  lambda = coord_search_low(this,divno,y);
  debug2(printf("coord_search_low yields lambda %d\n",lambda));

  *rightflanks = (int *) CALLOC(nflanking,sizeof(int));
  *nrightflanks = 0;
  stopp = false;
  while (lambda <= this->nintervals[divno] && stopp == false) {
    interval = &(this->intervals[divno][this->alphas[divno][lambda]-1]);
    if (sign != 0 && Interval_sign(interval) != sign) {
      debug2(printf("Advancing because sign != 0 && interval_sign %d != %d\n",Interval_sign(interval),sign));
      lambda++;
    } else if (Interval_low(interval) <= y) {
      debug2(printf("Advancing because interval_low %u <= %u\n",Interval_low(interval),y));
      lambda++;
    } else if (Interval_type(interval) != type) {
      debug2(printf("Advancing because interval_type %d != %d\n",Interval_type(interval),type));
      lambda++;
    } else {
      (*rightflanks)[(*nrightflanks)++] = this->alphas[divno][lambda];
      debug2(printf("Storing right flank %d\n",this->alphas[divno][lambda]));
      if (*nrightflanks < nflanking) {
	debug2(printf("Advancing because need more\n"));
	lambda++;
      } else {
	stopp = true;
      }
    }
  }

  /* Look at betas for left flank */
  lambda = coord_search_high(this,divno,x);

  *leftflanks = (int *) CALLOC(nflanking,sizeof(int));
  *nleftflanks = 0;
  stopp = false;
  while (lambda >= 1 && stopp == false) {
    interval = &(this->intervals[divno][this->betas[divno][lambda]-1]);
    if (sign != 0 && Interval_sign(interval) != sign) {
      lambda--;
    } else if (Interval_high(interval) >= x) {
      lambda--;
    } else if (Interval_type(interval) != type) {
      lambda--;
    } else {
      (*leftflanks)[(*nleftflanks)++] = this->betas[divno][lambda];
      if (*nleftflanks < nflanking) {
	lambda--;
      } else {
	stopp = true;
      }
    }
  }

  /* Convert to universal indices */
  matchstart = this->cum_nintervals[divno];
  for (i = 0; i < *nrightflanks; i++) {
    (*rightflanks)[i] += matchstart;
  }
  for (i = 0; i < *nleftflanks; i++) {
    (*leftflanks)[i] += matchstart;
  }

  return;
}

void
IIT_get_flanking_multiple_typed (int **leftflanks, int *nleftflanks, int **rightflanks, int *nrightflanks,
				 T this, char *divstring, Chrpos_T x, Chrpos_T y, int nflanking, int *types, int ntypes) {
  int k, i;
  int lambda, matchstart;
  Interval_T interval;
  bool stopp;
  int divno;

  divno = IIT_divint(this,divstring);

  debug(printf("Entering IIT_get_flanking_multiple_typed with query %u %u\n",x,y));

  if (this->alphas[divno] == NULL) {
#if 0
    IIT_compute_flanking(this);
#else
    fprintf(stderr,"Flanking hits not supported on version %d of iit files.  Please use iit_update to update your file\n",
	    this->version);
    exit(9);
#endif
  }

  /* Look at alphas for right flank */
  lambda = coord_search_low(this,divno,y);

  *rightflanks = (int *) CALLOC(nflanking,sizeof(int));
  *nrightflanks = 0;
  stopp = false;
  while (lambda <= this->nintervals[divno] && stopp == false) {
    interval = &(this->intervals[divno][this->alphas[divno][lambda]-1]);
    if (Interval_low(interval) <= y) {
      lambda++;
    } else {
      k = 0;
      while (k < ntypes && Interval_type(interval) != types[k]) {
	k++;
      }
      if (k >= ntypes) {
	lambda++;
      } else {
	(*rightflanks)[(*nrightflanks)++] = this->alphas[divno][lambda];
	if (*nrightflanks < nflanking) {
	  lambda++;
	} else {
	  stopp = true;
	}
      }
    }
  }


  /* Look at betas for left flank */
  lambda = coord_search_high(this,divno,x);

  *leftflanks = (int *) CALLOC(nflanking,sizeof(int));
  *nleftflanks = 0;
  stopp = false;
  while (lambda >= 1 && stopp == false) {
    interval = &(this->intervals[divno][this->betas[divno][lambda]-1]);
    if (Interval_high(interval) >= x) {
      lambda--;
    } else {
      k = 0;
      while (k < ntypes && Interval_type(interval) != types[k]) {
	k++;
      }
      if (k >= ntypes) {
	lambda--;
      } else {
	(*leftflanks)[(*nleftflanks)++] = this->betas[divno][lambda];
	if (*nleftflanks < nflanking) {
	  lambda--;
	} else {
	  stopp = true;
	}
      }
    }
  }

  /* Convert to universal indices */
  matchstart = this->cum_nintervals[divno];
  for (i = 0; i < *nrightflanks; i++) {
    (*rightflanks)[i] += matchstart;
  }
  for (i = 0; i < *nleftflanks; i++) {
    (*leftflanks)[i] += matchstart;
  }

  return;
}


static const Except_T iit_error = { "IIT problem" };

int
IIT_get_one (T this, char *divstring, Chrpos_T x, Chrpos_T y) {
  int lambda;
  int min1, max1 = 0, min2, max2 = 0;
  int divno;
  bool stopp;
  Interval_T interval;

  divno = IIT_divint(this,divstring);
  min1 = min2 = this->nintervals[divno] + 1;

  debug(printf("Entering IIT_get_one with query %u %u\n",x,y));
  fnode_query_aux(&min1,&max1,this,divno,0,x);
  fnode_query_aux(&min2,&max2,this,divno,0,y);
  debug(printf("min1=%d max1=%d  min2=%d max2=%d\n",min1,max1,min2,max2));

  if (max2 >= min1) {
    for (lambda = min1; lambda <= max2; lambda++) {
      if (Interval_overlap_p(x,y,this->intervals[divno],this->sigmas[divno][lambda]) == true) {
	return this->sigmas[divno][lambda];
      }
    }
    for (lambda = min1; lambda <= max2; lambda++) {
      if (Interval_overlap_p(x,y,this->intervals[divno],this->omegas[divno][lambda]) == true) {
	return this->omegas[divno][lambda];
      }
    }
  }

  /* fprintf(stderr,"Expected one match for %u--%u, but got none\n",x,y); */
  /* If we miss (e.g., for circular chromosome), then report the chromosome below */
  /* Look at betas or omegas for left flank */
  lambda = min1 - 1;
  stopp = false;
  while (lambda >= 1 && stopp == false) {
    interval = &(this->intervals[divno][this->omegas[divno][lambda]-1]);
    if (Interval_high(interval) >= x) {
      lambda--;
    } else {
      return this->omegas[divno][lambda];
    }
  }

  return this->omegas[divno][/*lambda*/1];
}

/* Generally called where intervals don't overlap, like chromosomes,
   and where x == y. */
/*
int
IIT_get_one_safe (T this, Chrpos_T x, Chrpos_T y) {
  int index;
  int *matches, nmatches;

  matches = IIT_get(&nmatches,this,x,y,sortp);
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

int *
IIT_get_typed (int *ntypematches, T this, char *divstring, Chrpos_T x, Chrpos_T y, int type, bool sortp) {
  int *sorted;
  int index, divno;
  int *typematches = NULL, *matches, nmatches, i, j;
  Interval_T interval;

  *ntypematches = 0;
  matches = IIT_get(&nmatches,this,divstring,x,y,/*sortp*/false);
  for (i = 0; i < nmatches; i++) {
    index = matches[i];
    interval = &(this->intervals[0][index-1]);
    if (Interval_type(interval) == type) {
      (*ntypematches)++;
    }
  }

  if (*ntypematches > 0) {
    typematches = (int *) CALLOC(*ntypematches,sizeof(int));
    j = 0;
    for (i = 0; i < nmatches; i++) {
      index = matches[i];
      interval = &(this->intervals[0][index-1]);
      if (Interval_type(interval) == type) {
	typematches[j++] = index;
      }
    }
  }
  
  if (matches != NULL) {
    FREE(matches);
  }

  if (sortp == false) {
    return typematches;
#if 0
  } else if (this->version <= 2) {
    sorted = sort_matches_by_type(this,typematches,*ntypematches,/*alphabetizep*/false);
    FREE(typematches);
    return sorted;
#endif
  } else {
    divno = IIT_divint(this,divstring);
    sorted = sort_matches_by_position(this,typematches,*ntypematches);
    FREE(typematches);
    return sorted;
  }
}

int *
IIT_get_typed_with_divno (int *ntypematches, T this, int divno, Chrpos_T x, Chrpos_T y,
			  int type, bool sortp) {
  int *sorted;
  int index;
  int *typematches = NULL, *matches, nmatches, i, j;
  Interval_T interval;

  if (divno < 0) {
    /* fprintf(stderr,"No div %s found in iit file\n",divstring); */
    *ntypematches = 0;
    return (int *) NULL;
  }

  *ntypematches = 0;
  matches = IIT_get_with_divno(&nmatches,this,divno,x,y,/*sortp*/false);
  for (i = 0; i < nmatches; i++) {
    index = matches[i];
    interval = &(this->intervals[0][index-1]);
    if (Interval_type(interval) == type) {
      (*ntypematches)++;
    }
  }

  if (*ntypematches > 0) {
    typematches = (int *) CALLOC(*ntypematches,sizeof(int));
    j = 0;
    for (i = 0; i < nmatches; i++) {
      index = matches[i];
      interval = &(this->intervals[0][index-1]);
      if (Interval_type(interval) == type) {
	typematches[j++] = index;
      }
    }
  }
  
  if (matches != NULL) {
    FREE(matches);
  }

  if (sortp == false) {
    return typematches;
#if 0
  } else if (this->version <= 2) {
    sorted = sort_matches_by_type(this,typematches,*ntypematches,/*alphabetizep*/false);
    FREE(typematches);
    return sorted;
#endif
  } else {
    sorted = sort_matches_by_position(this,typematches,*ntypematches);
    FREE(typematches);
    return sorted;
  }
}


int *
IIT_get_typed_signed (int *ntypematches, T this, char *divstring, Chrpos_T x, Chrpos_T y,
		      int type, int sign, bool sortp) {
  int *sorted;
  int index, divno;
  int *typematches = NULL, *matches, nmatches, i, j;
  Interval_T interval;

  *ntypematches = 0;
  matches = IIT_get(&nmatches,this,divstring,x,y,/*sortp*/false);
  for (i = 0; i < nmatches; i++) {
    index = matches[i];
    interval = &(this->intervals[0][index-1]);
    if (Interval_type(interval) == type && (sign == 0 || Interval_sign(interval) == sign)) {
      (*ntypematches)++;
    }
  }

  if (*ntypematches > 0) {
    typematches = (int *) CALLOC(*ntypematches,sizeof(int));
    j = 0;
    for (i = 0; i < nmatches; i++) {
      index = matches[i];
      interval = &(this->intervals[0][index-1]);
      if (Interval_type(interval) == type && (sign == 0 || Interval_sign(interval) == sign)) {
	typematches[j++] = index;
      }
    }
  }
  
  if (matches != NULL) {
    FREE(matches);
  }

  if (sortp == false) {
    return typematches;
#if 0
  } else if (this->version <= 2) {
    sorted = sort_matches_by_type(this,typematches,*ntypematches,/*alphabetizep*/false);
    FREE(typematches);
    return sorted;
#endif
  } else {
    divno = IIT_divint(this,divstring);
    sorted = sort_matches_by_position(this,typematches,*ntypematches);
    FREE(typematches);
    return sorted;
  }
}


int *
IIT_get_typed_signed_with_divno (int *ntypematches, T this, int divno, Chrpos_T x, Chrpos_T y, 
				 int type, int sign, bool sortp) {
  int *sorted;
  int index;
  int *typematches = NULL, *matches, nmatches, i, j;
  Interval_T interval;

  if (divno < 0) {
    /* fprintf(stderr,"No div %s found in iit file\n",divstring); */
    *ntypematches = 0;
    return (int *) NULL;
  }

  *ntypematches = 0;
  matches = IIT_get_with_divno(&nmatches,this,divno,x,y,/*sortp*/false);
  for (i = 0; i < nmatches; i++) {
    index = matches[i];
    interval = &(this->intervals[0][index-1]);
    if (Interval_type(interval) == type && (sign == 0 || Interval_sign(interval) == sign)) {
      (*ntypematches)++;
    }
  }

  if (*ntypematches > 0) {
    typematches = (int *) CALLOC(*ntypematches,sizeof(int));
    j = 0;
    for (i = 0; i < nmatches; i++) {
      index = matches[i];
      interval = &(this->intervals[0][index-1]);
      if (Interval_type(interval) == type && (sign == 0 || Interval_sign(interval) == sign)) {
	typematches[j++] = index;
      }
    }
  }
  
  if (matches != NULL) {
    FREE(matches);
  }

  if (sortp == false) {
    return typematches;
#if 0
  } else if (this->version <= 2) {
    sorted = sort_matches_by_type(this,typematches,*ntypematches,/*alphabetizep*/false);
    FREE(typematches);
    return sorted;
#endif
  } else {
    sorted = sort_matches_by_position(this,typematches,*ntypematches);
    FREE(typematches);
    return sorted;
  }
}


int *
IIT_get_multiple_typed (int *ntypematches, T this, char *divstring, Chrpos_T x, Chrpos_T y, 
			int *types, int ntypes, bool sortp) {
  int *sorted;
  int index, divno;
  int *typematches = NULL, *matches, nmatches, i, j, k;
  Interval_T interval;

  *ntypematches = 0;
  matches = IIT_get(&nmatches,this,divstring,x,y,/*sortp*/false);
  for (i = 0; i < nmatches; i++) {
    index = matches[i];
    interval = &(this->intervals[0][index-1]);
    k = 0;
    while (k < ntypes && Interval_type(interval) != types[k]) {
      k++;
    }
    if (k < ntypes) {
      (*ntypematches)++;
    }
  }

  if (*ntypematches > 0) {
    typematches = (int *) CALLOC(*ntypematches,sizeof(int));
    j = 0;
    for (i = 0; i < nmatches; i++) {
      index = matches[i];
      interval = &(this->intervals[0][index-1]);
      k = 0;
      while (k < ntypes && Interval_type(interval) != types[k]) {
	k++;
      }
      if (k < ntypes) {
	typematches[j++] = index;
      }
    }
  }
  
  if (matches != NULL) {
    FREE(matches);
  }

  if (sortp == false || this->version >= 3) {
    return typematches;
#if 0
  } else if (this->version <= 2) {
    sorted = sort_matches_by_type(this,typematches,*ntypematches,/*alphabetizep*/true);
    FREE(typematches);
    return sorted;
#endif
  } else {
    divno = IIT_divint(this,divstring);
    sorted = sort_matches_by_position(this,typematches,*ntypematches);
    FREE(typematches);
    return sorted;
  }
}

int
IIT_get_exact (T this, char *divstring, Chrpos_T x, Chrpos_T y, int type) {
  int index;
  int *matches, nmatches, i;
  Interval_T interval;

  matches = IIT_get(&nmatches,this,divstring,x,y,/*sortp*/false);
  for (i = 0; i < nmatches; i++) {
    index = matches[i];
    interval = &(this->intervals[0][index-1]);
    if (Interval_low(interval) == x && Interval_high(interval) == y &&
	Interval_type(interval) == type) {
      FREE(matches);
      return index;
    }
  }

  FREE(matches);
  return -1;
}

bool
IIT_exact_p (T this, char *divstring, Chrpos_T x, Chrpos_T y, int type) {
  int index;
  int *matches, nmatches, i;
  Interval_T interval;

  if (x == y) {
    matches = IIT_get(&nmatches,this,divstring,x,y,/*sortp*/false);
    for (i = 0; i < nmatches; i++) {
      index = matches[i];
      interval = &(this->intervals[0][index-1]);
      if (Interval_low(interval) == x && Interval_high(interval) == y &&
	  Interval_sign(interval) == 0 && Interval_type(interval) == type) {
	FREE(matches);
	return true;
      }
    }

  } else if (x < y) {
    matches = IIT_get(&nmatches,this,divstring,x,y,/*sortp*/false);
    for (i = 0; i < nmatches; i++) {
      index = matches[i];
      interval = &(this->intervals[0][index-1]);
      if (Interval_low(interval) == x && Interval_high(interval) == y &&
	  Interval_sign(interval) > 0 && Interval_type(interval) == type) {
	FREE(matches);
	return true;
      }
    }

  } else {
    matches = IIT_get(&nmatches,this,divstring,y,x,/*sortp*/false);
    for (i = 0; i < nmatches; i++) {
      index = matches[i];
      interval = &(this->intervals[0][index-1]);
      if (Interval_low(interval) == x && Interval_high(interval) == y &&
	  Interval_sign(interval) < 0 && Interval_type(interval) == type) {
	FREE(matches);
	return true;
      }
    }
  }

  FREE(matches);
  return false;
}


int *
IIT_get_exact_multiple (int *nexactmatches, T this, char *divstring, Chrpos_T x, Chrpos_T y, int type) {
  int *exactmatches;
  int index;
  int *matches, nmatches, i, j;
  Interval_T interval;

  *nexactmatches = 0;
  matches = IIT_get(&nmatches,this,divstring,x,y,/*sortp*/false);
  for (i = 0; i < nmatches; i++) {
    index = matches[i];
    interval = &(this->intervals[0][index-1]);
    if (Interval_low(interval) == x && Interval_high(interval) == y &&
	Interval_type(interval) == type) {
      (*nexactmatches)++;
    }
  }

  if (*nexactmatches == 0) {
    FREE(matches);
    return (int *) NULL;
  } else {
    exactmatches = (int *) CALLOC(*nexactmatches,sizeof(int));
    j = 0;
    for (i = 0; i < nmatches; i++) {
      index = matches[i];
      interval = &(this->intervals[0][index-1]);
      if (Interval_low(interval) == x && Interval_high(interval) == y &&
	  Interval_type(interval) == type) {
	exactmatches[j++] = index;
      }
    }
    FREE(matches);
    return exactmatches;
  }
}

int *
IIT_get_exact_multiple_with_divno (int *nexactmatches, T this, int divno, Chrpos_T x, Chrpos_T y, int type) {
  int *exactmatches;
  int index;
  int *matches, nmatches, i, j;
  Interval_T interval;

  *nexactmatches = 0;
  matches = IIT_get_with_divno(&nmatches,this,divno,x,y,/*sortp*/false);
  for (i = 0; i < nmatches; i++) {
    index = matches[i];
    interval = &(this->intervals[0][index-1]);
    if (Interval_low(interval) == x && Interval_high(interval) == y &&
	Interval_type(interval) == type) {
      (*nexactmatches)++;
    }
  }

  if (*nexactmatches == 0) {
    FREE(matches);
    return (int *) NULL;
  } else {
    exactmatches = (int *) CALLOC(*nexactmatches,sizeof(int));
    j = 0;
    for (i = 0; i < nmatches; i++) {
      index = matches[i];
      interval = &(this->intervals[0][index-1]);
      if (Interval_low(interval) == x && Interval_high(interval) == y &&
	  Interval_type(interval) == type) {
	exactmatches[j++] = index;
      }
    }
    FREE(matches);
    return exactmatches;
  }
}


/************************************************************************/

/* Modified from IIT_find */
int *
IIT_get_values_between (int *nmatches, T this, double lowval, double highval, bool sortp) {
  int *matches = NULL, j;
  double val;
  int start, end;
  int low, middle, high, recno;
  bool foundp;

  debug(printf("Entering IIT_get_values_between with %f to %f\n",lowval,highval));

  /* Find start */
  foundp = false;
  low = 0;
  high = this->total_nintervals;

#ifdef DEBUG
#ifndef WORDS_BIGENDIAN
  for (middle = low; middle < high; middle++) {
    printf("%d:%d:%f\n",middle,this->valueorder[middle],
	   this->values[this->valueorder[middle]]);
  }
  printf("\n");
#endif
#endif

  while (!foundp && low < high) {
    middle = (low+high)/2;

#ifdef DEBUG
#ifndef WORDS_BIGENDIAN
    printf("low %d middle %d:%d:%f high %d\n",
	   low,middle,this->valueorder[middle],
	   this->values[this->valueorder[middle]],high);
#endif
#endif

#ifdef WORDS_BIGENDIAN
    val = Bigendian_convert_double(this->values[Bigendian_convert_int(this->valueorder[middle])]);
#else
    val = this->values[this->valueorder[middle]];
#endif

    if (val > lowval) {
      high = middle;
      debug(printf("Decreasing high to %d\n",high));
    } else if (val < lowval) {
      low = middle + 1;
      debug(printf("Increasing low to %d\n",low));
    } else {
      foundp = true;
    }
  }

  if (foundp == true) {
    start = middle;
    debug(printf("start is middle = %d\n\n",start));

#ifdef WORDS_BIGENDIAN
    while (start-1 >= 0 && 
	   lowval == Bigendian_convert_double(this->values[Bigendian_convert_int(this->valueorder[start-1])])) {
      start--;
    }
#else
    while (start-1 >= 0 && 
	   lowval == this->values[this->valueorder[start-1]]) {
      start--;
      debug(printf("Regressing start to %d\n",start));
    }
#endif

  } else if ((start = low) >= this->total_nintervals) {
    *nmatches = 0;
    return (int *) NULL;

  } else {
    debug(printf("start is low = %d\n\n",start));
#ifdef WORDS_BIGENDIAN
    val = Bigendian_convert_double(this->values[Bigendian_convert_int(this->valueorder[start])]);
#else
    val = this->values[this->valueorder[start]];
#endif
    debug(printf("Final value for low bound = %f\n",val));
    if (val < lowval) {
      *nmatches = 0;
      return (int *) NULL;
    }
  }
    

  /* Find end */
  foundp = false;
  low = 0;
  high = this->total_nintervals;
  while (!foundp && low < high) {
    middle = (low+high)/2;

#ifdef DEBUG
#ifndef WORDS_BIGENDIAN
    printf("low %d middle %d:%d:%f high %d\n",
	   low,middle,this->valueorder[middle],
	   this->values[this->valueorder[middle]],high);
#endif
#endif

#ifdef WORDS_BIGENDIAN
    val = Bigendian_convert_double(this->values[Bigendian_convert_int(this->valueorder[middle])]);
#else
    val = this->values[this->valueorder[middle]];
#endif

    if (val > highval) {
      high = middle;
      debug(printf("Decreasing high to %d\n",high));
    } else if (val < highval) {
      low = middle + 1;
      debug(printf("Increasing low to %d\n",low));
    } else {
      foundp = true;
    }
  }

  if (foundp == true) {
    end = middle;
    debug(printf("end is middle = %d\n\n",end));

#ifdef WORDS_BIGENDIAN
    while (end+1 < this->total_nintervals && 
	   highval == Bigendian_convert_double(this->values[Bigendian_convert_int(this->valueorder[end+1])])) {
      end++;
    }
#else
    while (end+1 < this->total_nintervals && 
	   highval == this->values[this->valueorder[end+1]]) {
      end++;
      debug(printf("Advancing end to %d\n",end));
    }
#endif

  } else if ((end = high - 1) < 0) {
    *nmatches = 0;
    return (int *) NULL;

  } else {
    debug(printf("end is high - 1 = %d\n\n",end));

#ifdef WORDS_BIGENDIAN
    val = Bigendian_convert_double(this->values[Bigendian_convert_int(this->valueorder[end])]);
#else
    val = this->values[this->valueorder[end]];
#endif
    debug(printf("Final value for high bound = %f\n",val));
  
    if (val > highval) {
      *nmatches = 0;
      return (int *) NULL;
    }
  }
    
  *nmatches = end - start + 1;
  if (*nmatches <= 0) {
    *nmatches = 0;
    return (int *) NULL;
  } else {
    matches = (int *) CALLOC(*nmatches,sizeof(int));
    j = 0;
    for (recno = start; recno <= end; recno++) {
#ifdef WORDS_BIGENDIAN
#ifdef DEBUG
      printf("Pushing %d:%d\n",recno,Bigendian_convert_int(this->valueorder[recno]));
#endif
      matches[j++] = Bigendian_convert_int(this->valueorder[recno])+1;
	
#else
#ifdef DEBUG
      printf("Pushing %d:%d\n",recno,this->valueorder[recno]);
#endif
      matches[j++] = this->valueorder[recno]+1;
#endif
    }
    
    return matches;
  }
}


int *
IIT_get_values_below (int *nmatches, T this, double highval, bool sortp) {
  int *matches = NULL, j;
  double val;
  int start = 0, end;
  int low, middle, high, recno;
  bool foundp;

  debug(printf("Entering IIT_get_values_below with %f\n",highval));

  /* Find end */
  foundp = false;
  low = 0;
  high = this->total_nintervals;
  while (!foundp && low < high) {
    middle = (low+high)/2;

#ifdef DEBUG
#ifndef WORDS_BIGENDIAN
    printf("low %d middle %d:%d:%f high %d\n",
	   low,middle,this->valueorder[middle],
	   this->values[this->valueorder[middle]],high);
#endif
#endif

#ifdef WORDS_BIGENDIAN
    val = Bigendian_convert_double(this->values[Bigendian_convert_int(this->valueorder[middle])]);
#else
    val = this->values[this->valueorder[middle]];
#endif

    if (val > highval) {
      high = middle;
      debug(printf("Decreasing high to %d\n",high));
    } else if (val < highval) {
      low = middle + 1;
      debug(printf("Increasing low to %d\n",low));
    } else {
      foundp = true;
    }
  }

  if (foundp == true) {
    end = middle;
    debug(printf("end is middle = %d\n\n",end));

#ifdef WORDS_BIGENDIAN
    while (end+1 < this->total_nintervals && 
	   highval == Bigendian_convert_double(this->values[Bigendian_convert_int(this->valueorder[end+1])])) {
      end++;
    }
#else
    while (end+1 < this->total_nintervals && 
	   highval == this->values[this->valueorder[end+1]]) {
      end++;
      debug(printf("Advancing end to %d\n",end));
    }
#endif

  } else if ((end = high - 1) < 0) {
    *nmatches = 0;
    return (int *) NULL;

  } else {
    debug(printf("end is high - 1 = %d\n\n",end));

#ifdef WORDS_BIGENDIAN
    val = Bigendian_convert_double(this->values[Bigendian_convert_int(this->valueorder[end])]);
#else
    val = this->values[this->valueorder[end]];
#endif
    debug(printf("Final value for high bound = %f\n",val));
  
    if (val > highval) {
      *nmatches = 0;
      return (int *) NULL;
    }
  }

    
  *nmatches = end - start + 1;
  if (*nmatches <= 0) {
    *matches = 0;
    return (int *) NULL;
  } else {
    matches = (int *) CALLOC(*nmatches,sizeof(int));
    j = 0;
    for (recno = start; recno <= end; recno++) {
#ifdef WORDS_BIGENDIAN
#ifdef DEBUG
      printf("Pushing %d:%d\n",recno,Bigendian_convert_int(this->valueorder[recno]));
#endif
      matches[j++] = Bigendian_convert_int(this->valueorder[recno])+1;
	
#else
#ifdef DEBUG
      printf("Pushing %d:%d\n",recno,this->valueorder[recno]);
#endif
      matches[j++] = this->valueorder[recno]+1;
#endif
    }
    
    return matches;
  }
}


int *
IIT_get_values_above (int *nmatches, T this, double lowval, bool sortp) {
  int *matches = NULL, j;
  double val;
  int start, end = this->total_nintervals - 1;
  int low, middle, high, recno;
  bool foundp;

  debug(printf("Entering IIT_get_values_above with %f\n",lowval));

  /* Find start */
  foundp = false;
  low = 0;
  high = this->total_nintervals;

  while (!foundp && low < high) {
    middle = (low+high)/2;

#ifdef DEBUG
#ifndef WORDS_BIGENDIAN
    printf("low %d middle %d:%d:%f high %d\n",
	   low,middle,this->valueorder[middle],
	   this->values[this->valueorder[middle]],high);
#endif
#endif

#ifdef WORDS_BIGENDIAN
    val = Bigendian_convert_double(this->values[Bigendian_convert_int(this->valueorder[middle])]);
#else
    val = this->values[this->valueorder[middle]];
#endif

    if (val > lowval) {
      high = middle;
      debug(printf("Decreasing high to %d\n",high));
    } else if (val < lowval) {
      low = middle + 1;
      debug(printf("Increasing low to %d\n",low));
    } else {
      foundp = true;
    }
  }

  if (foundp == true) {
    start = middle;
    debug(printf("start is middle = %d\n\n",start));

#ifdef WORDS_BIGENDIAN
    while (start-1 >= 0 && 
	   lowval == Bigendian_convert_double(this->values[Bigendian_convert_int(this->valueorder[start-1])])) {
      start--;
    }
#else
    while (start-1 >= 0 && 
	   lowval == this->values[this->valueorder[start-1]]) {
      start--;
      debug(printf("Regressing start to %d\n",start));
    }
#endif

  } else if ((start = low) >= this->total_nintervals) {
    *nmatches = 0;
    return (int *) NULL;

  } else {
    debug(printf("start is low = %d\n\n",start));
#ifdef WORDS_BIGENDIAN
    val = Bigendian_convert_double(this->values[Bigendian_convert_int(this->valueorder[start])]);
#else
    val = this->values[this->valueorder[start]];
#endif
    debug(printf("Final value for low bound = %f\n",val));
    if (val < lowval) {
      *nmatches = 0;
      return (int *) NULL;
    }
  }
    

  *nmatches = end - start + 1;
  if (*nmatches <= 0) {
    *matches = 0;
    return (int *) NULL;
  } else {
    matches = (int *) CALLOC(*nmatches,sizeof(int));
    j = 0;
    for (recno = start; recno <= end; recno++) {
#ifdef WORDS_BIGENDIAN
#ifdef DEBUG
      printf("Pushing %d:%d\n",recno,Bigendian_convert_int(this->valueorder[recno]));
#endif
      matches[j++] = Bigendian_convert_int(this->valueorder[recno])+1;
	
#else
#ifdef DEBUG
      printf("Pushing %d:%d\n",recno,this->valueorder[recno]);
#endif
      matches[j++] = this->valueorder[recno]+1;
#endif
    }
    
    return matches;
  }
}



/************************************************************************/

#if 0
/* Need to work on */
/* Retrieves intervals from an IIT where type > 0.  Used by gmapindex to 
   construct altstrain_iit.  Here, the iit is a contig_iit.  */
List_T
IIT_intervallist_typed (List_T *labellist, Uintlist_T *seglength_list, T this) {
  List_T intervallist = NULL;
  Interval_T interval;
  char *label, *annotation, *restofheader, firstchar;
  bool allocp;
  int i;
  Chrpos_T seglength;

  *labellist = NULL;
  *seglength_list = NULL;
  for (i = 0; i < this->nintervals; i++) {
    interval = &(this->intervals[i]);
    if (Interval_type(interval) > 0) {
      intervallist = List_push(intervallist,Interval_copy(interval));
      label = IIT_label(this,i+1,&allocp);
      *labellist = List_push(*labellist,label);

      if (this->version <= 1) {
	/* Annotation may be negative to indicate contig is reverse complement */
	annotation = IIT_annotation(&restofheader,this,i+1,&allocp);
	firstchar = annotation[0];
	if (firstchar == '-') {
	  seglength = (Chrpos_T) strtoul(&(annotation[1]),NULL,10);
	} else {
	  seglength = (Chrpos_T) strtoul(annotation,NULL,10);
	  *seglength_list = Uintlist_push(*seglength_list,seglength);
	}
	if (allocp == true) {
	  FREE(restofheader);
	}
      } else {
	seglength = (Chrpos_T) strtoul(annotation,NULL,10);
	*seglength_list = Uintlist_push(*seglength_list,seglength);
      }
    }
  }
  *labellist = List_reverse(*labellist);
  *seglength_list = Uintlist_reverse(*seglength_list);
  return List_reverse(intervallist);
}
#endif

List_T
IIT_typelist (T this) {
  List_T typelist = NULL;
  int i;
  char *typestring, *copy;

  for (i = 0; i < this->ntypes; i++) {
    typestring = IIT_typestring(this,i);
    copy = (char *) CALLOC(strlen(typestring)+1,sizeof(char));
    strcpy(copy,typestring);
    typelist = List_push(typelist,copy);
  }
  return List_reverse(typelist);
}
					  

/************************************************************************/

/* Assume 0-based index */
static void
print_header (Filestring_T fp, T this, int recno, char *chr, bool map_bothstrands_p,
	      bool relativep, Chrpos_T left, bool print_comment_p) {
  char *string, *restofheader, *p;
  Interval_T interval;
  bool allocp;
#if 0
  int typeint;
#endif

  string = IIT_label(this,recno+1,&allocp);

  FPRINTF(fp,"\t%s",this->name);

  interval = &(this->intervals[0][recno]);
  if (relativep == true) {
    if (Interval_sign(interval) >= 0) {
      FPRINTF(fp,"\t%u..%u",Interval_low(interval)-left,Interval_high(interval)-left);
    } else {
      FPRINTF(fp,"\t%u..%u",Interval_high(interval)-left,Interval_low(interval)-left);
    }
  } else {
    if (Interval_sign(interval) >= 0) {
      FPRINTF(fp,"\t%s:%u..%u",chr,Interval_low(interval),Interval_high(interval));
    } else {
      FPRINTF(fp,"\t%s:%u..%u",chr,Interval_high(interval),Interval_low(interval));
    }
  }

#if 0
  if (map_bothstrands_p == true) {
    if ((typeint = Interval_type(interval)) <= 0) {
      FPRINTF(fp,"\t\t%s",string);
    } else {
      FPRINTF(fp,"\t%s\t%s",IIT_typestring(this,typeint),string);
    }
  } else {
#endif
    FPRINTF(fp,"\t");
    p = string;
    while (*p != '\0' && *p != '\n') {
      PUTC(*p,fp);
      p++;
    }

#if 0
  }
#endif

  if (allocp == true) {
    FREE(string);
  }

  if (print_comment_p == true) {
    p = IIT_annotation(&restofheader,this,recno+1,&allocp);
    FPRINTF(fp,"\t");
    while (*p != '\0' && *p != '\n') {
      PUTC(*p,fp);
      p++;
    }

    if (allocp == true) {
      FREE(restofheader);
    }
  }

  FPRINTF(fp,"\n");

  return;
}


void
IIT_print_header (Filestring_T fp, T this, int *matches, int nmatches, bool map_bothstrands_p,
		  char *chr, bool reversep, bool relativep, Chrpos_T left,
		  bool print_comment_p) {
  int recno, i;

  if (reversep == true) {
    for (i = nmatches-1; i >= 0; i--) {
      recno = matches[i] - 1;	/* Convert to 0-based */
      print_header(fp,this,recno,chr,map_bothstrands_p,relativep,left,print_comment_p);
    }
  } else {
    for (i = 0; i < nmatches; i++) {
      recno = matches[i] - 1;	/* Convert to 0-based */
      print_header(fp,this,recno,chr,map_bothstrands_p,relativep,left,print_comment_p);
    }
  }

  return;
}


Overlap_T
IIT_gene_overlap (T map_iit, int divno, Chrpos_T x, Chrpos_T y, bool favor_multiexon_p) {
  int *matches, index;
  int nmatches, i;
  Chrpos_T exonstart, exonend;
  int observed_genestrand;
  char *annot, *restofheader, *p;
  bool allocp = false;
  bool multiexon_p;
  bool foundp = false;

  matches = IIT_get_with_divno(&nmatches,map_iit,divno,x,y,/*sortp*/false);

  for (i = 0; i < nmatches; i++) {
    index = matches[i];
    observed_genestrand = IIT_interval_sign(map_iit,index);
#if 0
    if (observed_genestrand > 0 && desired_genestrand < 0) {
      /* Inconsistent */
    } else if (observed_genestrand < 0 && desired_genestrand > 0) {
      /* Inconsistent */
    } else {
#endif
      annot = IIT_annotation(&restofheader,map_iit,index,&allocp);

      /* Skip header */
      p = annot;
      while (*p != '\0' && *p != '\n') {
	p++;
      }
      if (*p == '\n') p++;

      if (observed_genestrand > 0) {
	multiexon_p = false;
	while (*p != '\0') {
	  if (sscanf(p,"%u %u",&exonstart,&exonend) != 2) {
	    fprintf(stderr,"Can't parse exon coordinates in %s\n",p);
	    abort();
	  } else {
	    /* Advance to next exon */
	    while (*p != '\0' && *p != '\n') p++;
	    if (*p == '\n') p++;
	    if (*p != '\0') {
	      multiexon_p = true;
	    }
	    
	    if (exonend < x) {
	      /* No overlap */
	    } else if (exonstart > y) {
	      /* No overlap */
	    } else if (favor_multiexon_p == true) {
	      if (multiexon_p == true) {
		FREE(matches);
		if (allocp) FREE(annot);
		return KNOWN_GENE_MULTIEXON;
	      } else {
		/* Keep searching for a multi-exon gene */
		foundp = true;
	      }
	    } else {
	      FREE(matches);
	      if (allocp) FREE(annot);
	      return KNOWN_GENE;
	    }
	  }

	}

      } else {
	multiexon_p = false;
	while (*p != '\0') {
	  if (sscanf(p,"%u %u",&exonstart,&exonend) != 2) {
	    fprintf(stderr,"Can't parse exon coordinates in %s\n",p);
	    abort();
	  } else {
	    /* Advance to next exon */
	    while (*p != '\0' && *p != '\n') p++;
	    if (*p == '\n') p++;
	    if (*p != '\0') {
	      multiexon_p = true;
	    }

	    if (exonstart < x) {
	      /* No overlap */
	    } else if (exonend > y) {
	      /* No overlap */
	    } else if (favor_multiexon_p == true) {
	      if (multiexon_p == true) {
		FREE(matches);
		if (allocp) FREE(annot);
		return KNOWN_GENE_MULTIEXON;
	      } else {
		/* Keep searching for a multi-exon gene */
		foundp = true;
	      }
	    } else {
	      FREE(matches);
	      if (allocp) FREE(annot);
	      return KNOWN_GENE;
	    }
	  }
	}
      }
#if 0
    }
#endif
  }
  
  FREE(matches);
  if (allocp) FREE(annot);
  if (foundp == true) {
    return KNOWN_GENE;
  } else {
    return NO_KNOWN_GENE;
  }
}


