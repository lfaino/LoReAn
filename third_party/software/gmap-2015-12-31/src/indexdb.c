static char rcsid[] = "$Id: indexdb.c 168395 2015-06-26 17:13:13Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
# define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
# define memmove(d,s,n) bcopy((s),(d),(n))
#endif
#ifdef HAVE_SSE2
#include <emmintrin.h>
#endif

#include "indexdb.h"
#include "indexdbdef.h"


#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
#else
#include "littleendian.h"
#endif

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>		/* For memset */
#include <ctype.h>		/* For toupper */
#include <sys/mman.h>		/* For munmap */

#ifdef HAVE_UNISTD_H
#include <unistd.h>		/* For lseek and close */
#endif
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>		/* For off_t */
#endif
#if HAVE_DIRENT_H
# include <dirent.h>
# define NAMLEN(dirent) strlen((dirent)->d_name)
#else
# define dirent direct
# define NAMLEN(dirent) (dirent)->d_namlen
# if HAVE_SYS_NDIR_H
#  include <sys/ndir.h>
# endif
# if HAVE_SYS_DIR_H
#  include <sys/dir.h>
# endif
# if HAVE_NDIR_H
#  include <ndir.h>
# endif
#endif

#include "mem.h"
#include "fopen.h"
#include "types.h"		/* For Oligospace_T */

#include "compress.h"
#include "interval.h"
#include "complement.h"
#include "bitpack64-read.h"
#include "bitpack64-readtwo.h"


#ifdef HAVE_PTHREAD
#include <pthread.h>		/* sys/types.h already included above */
#endif

#define MAXENTRIES 20

/* Note: NONMODULAR is the old behavior.  Now we store only when
   startposition % index1interval == 0 */


/* Low-level codon hacking */
#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Calls to Indexdb_read */
#ifdef DEBUG0
#define debug0(x) x
#else
#define debug0(x)
#endif

/* Shifting of high word to low word for PMAP */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif


#ifdef PMAP

#if (defined(DEBUG) || defined(DEBUG0))
static Width_T index1part_aa;
#endif

void
Indexdb_setup (Width_T index1part_aa_in) {
#if (defined(DEBUG) || defined(DEBUG0))
  index1part_aa = index1part_aa_in;
#endif
  return;
}

#else

#define poly_A 0U
static Storedoligomer_T poly_T;  /* Was LOW12MER 0x00FFFFFF */

#if (defined(DEBUG) || defined(DEBUG0))
static Width_T index1part;
#endif

void
Indexdb_setup (Width_T index1part_in) {
#if (defined(DEBUG) || defined(DEBUG0))
  index1part = index1part_in;
#endif

  poly_T = ~(~0UL << 2*index1part_in);
  return;
}
#endif


#define T Indexdb_T


void
Indexdb_free (T *old) {
  if (*old) {
    if ((*old)->positions_access == ALLOCATED_PRIVATE) {
#ifdef LARGE_GENOMES
      FREE((*old)->positions_high);
      FREE((*old)->positions_low);
#else
      FREE((*old)->positions);
#endif

    } else if ((*old)->positions_access == ALLOCATED_SHARED) {
#ifdef LARGE_GENOMES
      Access_deallocate((*old)->positions_high,(*old)->positions_high_shmid);
      Access_deallocate((*old)->positions_low,(*old)->positions_low_shmid);
#else
      Access_deallocate((*old)->positions,(*old)->positions_shmid);
#endif

#ifdef HAVE_MMAP
    } else if ((*old)->positions_access == MMAPPED) {
#ifdef LARGE_GENOMES
      munmap((void *) (*old)->positions_high,(*old)->positions_high_len);
      close((*old)->positions_high_fd);
      munmap((void *) (*old)->positions_low,(*old)->positions_low_len);
      close((*old)->positions_low_fd);
#else
      munmap((void *) (*old)->positions,(*old)->positions_len);
      close((*old)->positions_fd);
#endif
#endif
    } else if ((*old)->positions_access == FILEIO) {
#ifdef HAVE_PTHREAD
      pthread_mutex_destroy(&(*old)->positions_read_mutex);
#endif
#ifdef LARGE_GENOMES
      close((*old)->positions_high_fd);
      close((*old)->positions_low_fd);
#else
      close((*old)->positions_fd);
#endif
    }

    if ((*old)->offsetsstrm_access == ALLOCATED_PRIVATE) {
      FREE((*old)->offsetsstrm);

    } else if ((*old)->offsetsstrm_access == ALLOCATED_SHARED) {
      Access_deallocate((*old)->offsetsstrm,(*old)->offsetsstrm_shmid);

#ifdef HAVE_MMAP
    } else if ((*old)->offsetsstrm_access == MMAPPED) {
      munmap((void *) (*old)->offsetsstrm,(*old)->offsetsstrm_len);
      close((*old)->offsetsstrm_fd);
#endif
    }
      
    if ((*old)->offsetsmeta_access == ALLOCATED_PRIVATE) {
      FREE((*old)->offsetsmeta);
    } else if ((*old)->offsetsmeta_access == ALLOCATED_SHARED) {
      Access_deallocate((*old)->offsetsmeta,(*old)->offsetsmeta_shmid);
    } else {
      /* Always ALLOCATED */
      abort();
    }

#ifdef LARGE_GENOMES
    if ((*old)->offsetspages_access == ALLOCATED_PRIVATE) {
      FREE((*old)->offsetspages);
    } else if ((*old)->offsetspages_access == ALLOCATED_SHARED) {
      Access_deallocate((*old)->offsetspages,(*old)->offsetspages_shmid);
    } else {
      /* Always ALLOCATED */
      abort();
    }
#endif

    FREE(*old);
  }
  return;
}


Width_T
Indexdb_interval (T this) {
  return this->index1interval;
}


bool
Indexdb_positions_fileio_p (T this) {
  if (this->positions_access == FILEIO) {
    return true;
  } else {
    return false;
  }
}

static Oligospace_T
power (int base, Width_T exponent) {
#ifdef OLIGOSPACE_NOT_LONG
  Oligospace_T result = 1U;
#else
  Oligospace_T result = 1UL;
#endif
  int i;

  for (i = 0; i < exponent; i++) {
    result *= base;
  }
  return result;
}

double
Indexdb_mean_size (T this, Mode_T mode, Width_T index1part) {
  Oligospace_T oligospace, n;

#ifdef PMAP
  /* index1part should be in aa */
  n = oligospace = power(this->alphabet_size,index1part);
#else
  n = oligospace = power(4,index1part);
  if (mode != STANDARD) {
    n = power(3,index1part);
  }
#endif

#ifdef WORDS_BIGENDIAN
  /* Also holds for ALLOCATED_PRIVATE and ALLOCATED_SHARED */
  return (double) Bigendian_convert_uint(this->offsetsstrm[Bigendian_convert_uint(this->offsetsmeta[oligospace/this->blocksize])])/(double) n;
#else
  return (double) this->offsetsstrm[this->offsetsmeta[oligospace/this->blocksize]]/(double) n;
#endif
}



static Filenames_T
Filenames_new (char *pages_filename, char *pointers_filename, char *offsets_filename,
	       char *positions_high_filename, char *positions_low_filename,
	       char *pointers_basename_ptr, char *offsets_basename_ptr,
	       char *positions_high_basename_ptr, char *positions_low_basename_ptr,
	       char *pointers_index1info_ptr, char *offsets_index1info_ptr,
	       char *positions_high_index1info_ptr, char *positions_low_index1info_ptr) {
  Filenames_T new = (Filenames_T) MALLOC(sizeof(*new));

  new->pages_filename = pages_filename;
  new->pointers_filename = pointers_filename;
  new->offsets_filename = offsets_filename;
  new->positions_high_filename = positions_high_filename;
  new->positions_low_filename = positions_low_filename;

  new->pointers_basename_ptr = pointers_basename_ptr;
  new->offsets_basename_ptr = offsets_basename_ptr;
  new->positions_high_basename_ptr = positions_high_basename_ptr;
  new->positions_low_basename_ptr = positions_low_basename_ptr;

  new->pointers_index1info_ptr = pointers_index1info_ptr;
  new->offsets_index1info_ptr = offsets_index1info_ptr;
  new->positions_high_index1info_ptr = positions_high_index1info_ptr;
  new->positions_low_index1info_ptr = positions_low_index1info_ptr;

  return new;
}

void
Filenames_free (Filenames_T *old) {

  FREE((*old)->pages_filename);
  FREE((*old)->pointers_filename);
  FREE((*old)->offsets_filename);
  FREE((*old)->positions_high_filename);
  FREE((*old)->positions_low_filename);

  FREE(*old);

  return;
}


Filenames_T
Indexdb_get_filenames_no_compression (Width_T *index1part, Width_T *index1interval,
				      char *genomesubdir, char *fileroot, char *idx_filesuffix, char *snps_root,
				      Width_T required_interval, bool offsets_only_p) {
  char *offsets_filename, *positions_high_filename, *positions_low_filename,
    *offsets_basename_ptr, *positions_high_basename_ptr, *positions_low_basename_ptr,
    *offsets_index1info_ptr, *positions_high_index1info_ptr, *positions_low_index1info_ptr;

  char *base_filename, *filename;
  char *pattern, interval_char, digit_string[2], *p, *q;
  char tens, ones;
  Width_T found_index1part, found_interval;
  int rootlength, patternlength;

  char *offsets_suffix, *positions_high_suffix, *positions_low_suffix;
  struct dirent *entry;
  DIR *dp;


  if (snps_root == NULL) {
    offsets_suffix = "offsets";
    positions_high_suffix = POSITIONS_HIGH_FILESUFFIX;
    positions_low_suffix = POSITIONS_LOW_FILESUFFIX;
  } else {
    offsets_suffix = (char *) CALLOC(strlen("offsets")+strlen(".")+strlen(snps_root)+1,sizeof(char));
    sprintf(offsets_suffix,"%s.%s","offsets",snps_root);
    positions_high_suffix = (char *) CALLOC(strlen(POSITIONS_HIGH_FILESUFFIX)+strlen(".")+strlen(snps_root)+1,sizeof(char));
    sprintf(positions_high_suffix,"%s.%s",POSITIONS_HIGH_FILESUFFIX,snps_root);
    positions_low_suffix = (char *) CALLOC(strlen(POSITIONS_LOW_FILESUFFIX)+strlen(".")+strlen(snps_root)+1,sizeof(char));
    sprintf(positions_low_suffix,"%s.%s",POSITIONS_LOW_FILESUFFIX,snps_root);
  }

  *index1part = 0;
  *index1interval = 1000;
  base_filename = (char *) NULL;

  if ((dp = opendir(genomesubdir)) == NULL) {
    fprintf(stderr,"Unable to open directory %s\n",genomesubdir);
    exit(9);
  }

  pattern = (char *) CALLOC(strlen(fileroot)+strlen(".")+strlen(idx_filesuffix)+1,sizeof(char));
  sprintf(pattern,"%s.%s",fileroot,idx_filesuffix);
  patternlength = strlen(pattern);

  digit_string[1] = '\0';	/* Needed for atoi */
  while ((entry = readdir(dp)) != NULL) {
    filename = entry->d_name;
    if (!strncmp(filename,pattern,patternlength)) {
      p = &(filename[strlen(pattern)]); /* Points after idx_filesuffix, e.g., "ref" */
      if ((q = strstr(p,offsets_suffix)) != NULL && !strcmp(q,offsets_suffix)) {
	if (q - p == 3) {
#ifdef PMAP
	  /* e.g., pf677 */
	  if (sscanf(p,"%c%c%c",&ones0,&ones,&interval_char) == 3) {
	    /* digit_string[0] = ones0; */
	    /* found_basesize = atoi(digit_string); */

	    digit_string[0] = ones;
	    found_index1part = atoi(digit_string);
	    
	    digit_string[0] = interval_char;
	    found_interval = atoi(digit_string);
	  } else {
	    abort();
	  }
#else
	  /* e.g., ref123offsets */
	  if (sscanf(p,"%c%c%c",&tens,&ones,&interval_char) == 3) {
	    digit_string[0] = tens;
	    found_index1part = 10*atoi(digit_string);
	    digit_string[0] = ones;
	    found_index1part += atoi(digit_string);

	    digit_string[0] = interval_char;
	    found_interval = atoi(digit_string);
	  } else {
	    abort();
	  }
#endif

	} else if (q - p == 1) {
	  /* Old style, e.g, idx or ref3 */
	  if (sscanf(p,"%c",&interval_char) == 1) {
	    if (interval_char == 'x') {
	      found_interval = 6;
	    } else {
	      digit_string[0] = interval_char;
	      found_interval = atoi(digit_string);
	    }
	  } else {
	    abort();
	  }
#ifdef PMAP
	  found_index1part = found_interval;
#else
	  found_index1part = 12;
#endif

	} else {
	  fprintf(stderr,"Cannot parse part between %s and offsets in filename %s\n",idx_filesuffix,filename);
	  if (snps_root != NULL) {
	    FREE(offsets_suffix);
	    FREE(positions_high_suffix);
	    FREE(positions_low_suffix);
	  }
	  return (Filenames_T) NULL;
	}

	if (required_interval != 0) {
	  if (found_interval == required_interval) {
	    *index1part = found_index1part;
	    *index1interval = found_interval;
	    FREE(base_filename);
	    base_filename = (char *) CALLOC(strlen(filename)+1,sizeof(char));
	    strcpy(base_filename,filename);
	  }
	} else {
	  if (found_interval < *index1interval) {
	    *index1part = found_index1part;
	    *index1interval = found_interval;
	    FREE(base_filename);
	    base_filename = (char *) CALLOC(strlen(filename)+1,sizeof(char));
	    strcpy(base_filename,filename);
	  }
	}
      }
    }
  }

  FREE(pattern);

  if (closedir(dp) < 0) {
    fprintf(stderr,"Unable to close directory %s\n",genomesubdir);
  }

  /* Construct full filenames */
  if (base_filename == NULL) {
#if 0
    fprintf(stderr,"Cannot find offsets file containing %s and %s",idx_filesuffix,offsets_suffix);
    if (required_interval > 0) {
      fprintf(stderr," and having sampling interval of %d",required_interval);
    }
    fprintf(stderr,"\n");
#endif

    /* offsets_filename = (char *) NULL; */
    /* positions_high_filename = (char *) NULL; */
    /* positions_low_filename = (char *) NULL; */
    if (snps_root != NULL) {
      FREE(offsets_suffix);
      FREE(positions_high_suffix);
      FREE(positions_low_suffix);
    }
    return (Filenames_T) NULL;

  } else {
    offsets_filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(base_filename)+1,sizeof(char));
    offsets_basename_ptr = &(offsets_filename[strlen(genomesubdir)+strlen("/")]);
    offsets_index1info_ptr = &(offsets_basename_ptr[patternlength]);

    sprintf(offsets_filename,"%s/%s",genomesubdir,base_filename);
    if (Access_file_exists_p(offsets_filename) == false) {
      fprintf(stderr,"Offsets filename %s does not exist\n",offsets_filename);
      FREE(offsets_filename);
      /* offsets_filename = (char *) NULL; */
      /* positions_high_filename = (char *) NULL; */
      /* positions_low_filename = (char *) NULL; */
      FREE(base_filename);
      if (snps_root != NULL) {
	FREE(offsets_suffix);
	FREE(positions_high_suffix);
	FREE(positions_low_suffix);
      }
      return (Filenames_T) NULL;
    }


    if ((q = strstr(base_filename,offsets_suffix)) == NULL) {
      abort();
    } else {
      rootlength = q - base_filename;
    }

    positions_high_filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+rootlength+strlen(positions_high_suffix)+1,sizeof(char));
    positions_high_basename_ptr = &(positions_high_filename[strlen(genomesubdir)+strlen("/")]);
    positions_high_index1info_ptr = &(positions_high_basename_ptr[patternlength]);

    positions_low_filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+rootlength+strlen(positions_low_suffix)+1,sizeof(char));
    positions_low_basename_ptr = &(positions_low_filename[strlen(genomesubdir)+strlen("/")]);
    positions_low_index1info_ptr = &(positions_low_basename_ptr[patternlength]);

    sprintf(positions_high_filename,"%s/",genomesubdir);
    strncpy(positions_high_basename_ptr,base_filename,rootlength);
    strcpy(&(positions_high_basename_ptr[rootlength]),positions_high_suffix);

    sprintf(positions_low_filename,"%s/",genomesubdir);
    strncpy(positions_low_basename_ptr,base_filename,rootlength);
    strcpy(&(positions_low_basename_ptr[rootlength]),positions_low_suffix);

    if (offsets_only_p == true) {
      /* Do not look for a positions file */
    } else if (Access_file_exists_p(positions_low_filename) == false) {
      fprintf(stderr,"Positions filename %s does not exist\n",positions_low_filename);
      FREE(offsets_filename);
      FREE(positions_high_filename);
      FREE(positions_low_filename);
      /* offsets_filename = (char *) NULL; */
      /* positions_high_filename = (char *) NULL; */
      /* positions_low_filename = (char *) NULL; */
      FREE(base_filename);
      if (snps_root != NULL) {
	FREE(offsets_suffix);
	FREE(positions_high_suffix);
	FREE(positions_low_suffix);
      }
      return (Filenames_T) NULL;
    }

    if (Access_file_exists_p(positions_high_filename) == false) {
      /* Not a large genome */
      FREE(positions_high_filename);
      positions_high_filename = (char *) NULL;
    }

    if (snps_root != NULL) {
      FREE(offsets_suffix);
      FREE(positions_high_suffix);
      FREE(positions_low_suffix);
    }

    FREE(base_filename);
    fprintf(stderr,"Looking for index files in directory %s (offsets not compressed)\n",genomesubdir);
    fprintf(stderr,"  Offsets file is %s\n",offsets_basename_ptr);
    if (positions_high_filename == NULL) {
      fprintf(stderr,"  Positions file is %s\n",positions_low_basename_ptr);
    } else {
      fprintf(stderr,"  Positions files are %s and %s\n",positions_high_basename_ptr,positions_low_basename_ptr);
    }
    return Filenames_new(/*pages_filename*/NULL,/*pointers_filename*/NULL,offsets_filename,
			 positions_high_filename,positions_low_filename,
			 /*pointers_basename_ptr*/NULL,offsets_basename_ptr,
			 positions_high_basename_ptr,positions_low_basename_ptr,
			 /*pointers_index1info_ptr*/NULL,offsets_index1info_ptr,
			 positions_high_index1info_ptr,positions_low_index1info_ptr);
  }
}



#ifdef PMAP
#define BASE_KMER_SAMPLING 3   /* e.g., 677 */
#define KMER_SAMPLING 2   /* e.g., 77 */
#else
#define BASE_KMER_SAMPLING 5   /* e.g., 12153 */
#define KMER_SAMPLING 3   /* e.g., 153 */
#endif


Filenames_T
Indexdb_get_filenames_bitpack (Width_T *index1part, Width_T *index1interval,
			       char *genomesubdir, char *fileroot, char *idx_filesuffix, char *snps_root,
			       Width_T required_index1part, Width_T required_interval,
			       Blocksize_T blocksize, bool offsets_only_p) {
  char *offsetspages_filename, *offsetsmeta_filename, *offsetsstrm_filename,
    *positions_high_filename, *positions_low_filename,
    *offsetspages_basename_ptr, *offsetsmeta_basename_ptr, *offsetsstrm_basename_ptr,
    *positions_high_basename_ptr, *positions_low_basename_ptr,
    *offsetspages_index1info_ptr, *offsetsmeta_index1info_ptr, *offsetsstrm_index1info_ptr,
    *positions_high_index1info_ptr, *positions_low_index1info_ptr;

  char *base_filename, *filename;
#ifdef PMAP
  char *pattern1, *pattern2, *a;
  int patternlength1, patternlength2, alphabet_strlen;
  Alphabet_T found_alphabet;
#else
  char *pattern;
  char tens;
#endif
  char interval_char, digit_string[2], *p, *q;
  Width_T found_index1part = 0, found_interval = 0;
  int rootlength, patternlength;

  char ones;
  char *offsetspages_suffix, *offsetsmeta_suffix, *offsetsstrm_suffix,
    *positions_high_suffix, *positions_low_suffix;
  struct dirent *entry;
  DIR *dp;


  if (snps_root == NULL) {
    if (blocksize == 32) {
      offsetspages_suffix = "offsets32pages";
      offsetsmeta_suffix = "offsets32meta";
      offsetsstrm_suffix = "offsets32strm";
      positions_high_suffix = POSITIONS_HIGH_FILESUFFIX;
      positions_low_suffix = POSITIONS_LOW_FILESUFFIX;
    } else if (blocksize == 64) {
      offsetspages_suffix = "offsets64pages";
      offsetsmeta_suffix = "offsets64meta";
      offsetsstrm_suffix = "offsets64strm";
      positions_high_suffix = POSITIONS_HIGH_FILESUFFIX;
      positions_low_suffix = POSITIONS_LOW_FILESUFFIX;
    } else {
      fprintf(stderr,"Unexpected blocksize %d\n",blocksize);
      abort();
    }
  } else {
    if (blocksize == 32) {
      offsetspages_suffix = (char *) CALLOC(strlen("offsets32pages.")+strlen(snps_root)+1,sizeof(char));
      offsetsmeta_suffix = (char *) CALLOC(strlen("offsets32meta.")+strlen(snps_root)+1,sizeof(char));
      offsetsstrm_suffix = (char *) CALLOC(strlen("offsets32strm.")+strlen(snps_root)+1,sizeof(char));
      positions_high_suffix = (char *) CALLOC(strlen(POSITIONS_HIGH_FILESUFFIX)+strlen(".")+strlen(snps_root)+1,sizeof(char));
      positions_low_suffix = (char *) CALLOC(strlen(POSITIONS_LOW_FILESUFFIX)+strlen(".")+strlen(snps_root)+1,sizeof(char));

      sprintf(offsetspages_suffix,"offsets32pages.%s",snps_root);
      sprintf(offsetsmeta_suffix,"offsets32meta.%s",snps_root);
      sprintf(offsetsstrm_suffix,"offsets32strm.%s",snps_root);
      sprintf(positions_high_suffix,"%s.%s",POSITIONS_HIGH_FILESUFFIX,snps_root);
      sprintf(positions_low_suffix,"%s.%s",POSITIONS_LOW_FILESUFFIX,snps_root);

    } else if (blocksize == 64) {
      offsetspages_suffix = (char *) CALLOC(strlen("offsets64pages.")+strlen(snps_root)+1,sizeof(char));
      offsetsmeta_suffix = (char *) CALLOC(strlen("offsets64meta.")+strlen(snps_root)+1,sizeof(char));
      offsetsstrm_suffix = (char *) CALLOC(strlen("offsets64strm.")+strlen(snps_root)+1,sizeof(char));
      positions_high_suffix = (char *) CALLOC(strlen(POSITIONS_HIGH_FILESUFFIX)+strlen(".")+strlen(snps_root)+1,sizeof(char));
      positions_low_suffix = (char *) CALLOC(strlen(POSITIONS_LOW_FILESUFFIX)+strlen(".")+strlen(snps_root)+1,sizeof(char));

      sprintf(offsetspages_suffix,"offsets64pages.%s",snps_root);
      sprintf(offsetsmeta_suffix,"offsets64meta.%s",snps_root);
      sprintf(offsetsstrm_suffix,"offsets64strm.%s",snps_root);
      sprintf(positions_high_suffix,"%s.%s",POSITIONS_HIGH_FILESUFFIX,snps_root);
      sprintf(positions_low_suffix,"%s.%s",POSITIONS_LOW_FILESUFFIX,snps_root);

    } else {
      fprintf(stderr,"Unexpected blocksize %d\n",blocksize);
      abort();
    }
  }


#ifdef PMAP
  *alphabet = NALPHABETS + 1;
#endif
  *index1part = 0;
  *index1interval = 1000;
  base_filename = (char *) NULL;

  if ((dp = opendir(genomesubdir)) == NULL) {
    fprintf(stderr,"Unable to open directory %s\n",genomesubdir);
    exit(9);
  }

#ifdef PMAP
  pattern1 = (char *) CALLOC(strlen(fileroot)+strlen(".")+1,sizeof(char)); /* e.g., "hg19." */
  sprintf(pattern1,"%s.",fileroot);
  patternlength1 = strlen(pattern1);

  pattern2 = (char *) CALLOC(strlen(".")+strlen(idx_filesuffix)+1,sizeof(char)); /* e.g., ".pr" */
  sprintf(pattern2,".%s",idx_filesuffix);
  patternlength2 = strlen(pattern2);

  digit_string[1] = '\0';	/* Needed for atoi */
  while ((entry = readdir(dp)) != NULL) {
    filename = entry->d_name;
    if (!strncmp(filename,pattern1,patternlength1)) {
      a = &(filename[strlen(pattern1)]); /* Points after fileroot, e.g., "hg19." */
      if ((p = strstr(a,pattern2)) != NULL && (q = strstr(p,offsetsstrm_suffix)) != NULL && !strcmp(q,offsetsstrm_suffix)) {
	if ((found_alphabet = Alphabet_find(a)) != AA0) {
	  alphabet_strlen = p - a;
	  p += patternlength2;

	  if (q - p == KMER_SAMPLING) {
	    /* Latest style, e.g., pf77 */
	    if (sscanf(p,"%c%c",&ones,&interval_char) == 2) {
	      digit_string[0] = ones;
	      found_index1part = atoi(digit_string);

	      digit_string[0] = interval_char;
	      found_interval = atoi(digit_string);
	    } else {
	      abort();
	    }

	  } else if (q - p == BASE_KMER_SAMPLING) {
	    /* Previous style, e.g., pf677 */
	    if (sscanf(p,"%c%c%c",&ones0,&ones,&interval_char) == 3) {
	      digit_string[0] = ones0;
	      found_basesize = atoi(digit_string);

	      digit_string[0] = ones;
	      found_index1part = atoi(digit_string);

	      digit_string[0] = interval_char;
	      found_interval = atoi(digit_string);
	    } else {
	      abort();
	    }

	  } else {
	    /* fprintf(stderr,"Cannot parse part between %s and offsets in filename %s\n",idx_filesuffix,filename); */
	    if (snps_root != NULL) {
	      FREE(offsetsstrm_suffix);
	      FREE(offsetsmeta_suffix);
	      FREE(offsetspages_suffix);
	      FREE(positions_high_suffix);
	      FREE(positions_low_suffix);
	    }
	    return (Filenames_T) NULL;
	  }

	  if ((required_alphabet == AA0 || found_alphabet == required_alphabet) &&
	      (required_index1part == 0 || found_index1part == required_index1part) &&
	      (required_interval == 0 || found_interval == required_interval)) {
	    if (required_alphabet == AA0 && found_alphabet > *alphabet) {
	      /* Skip, since we have already found an earlier alphabet */
	    } else if (required_index1part == 0 && found_index1part < *index1part) {
	      /* Skip, since we have already found a larger index1part */
	    } else if (required_interval == 0 && found_interval > *index1interval) {
	      /* Skip, since we have already found a smaller interval */
	    } else {
	      patternlength = patternlength1 + alphabet_strlen + patternlength2;
	      *basesize = found_basesize;
	      *index1part = found_index1part;
	      *index1interval = found_interval;
	      *alphabet = found_alphabet;
	      FREE(base_filename);
	      base_filename = (char *) CALLOC(strlen(filename)+1,sizeof(char));
	      strcpy(base_filename,filename);
	    }
	  }
	}
      }
    }
  }

  FREE(pattern2);
  FREE(pattern1);

#else

  pattern = (char *) CALLOC(strlen(fileroot)+strlen(".")+strlen(idx_filesuffix)+1,sizeof(char));
  sprintf(pattern,"%s.%s",fileroot,idx_filesuffix);
  patternlength = strlen(pattern);

  digit_string[1] = '\0';	/* Needed for atoi */
  while ((entry = readdir(dp)) != NULL) {
    filename = entry->d_name;
    if (!strncmp(filename,pattern,patternlength)) {
      p = &(filename[strlen(pattern)]); /* Points after idx_filesuffix, e.g., "ref" */
      if ((q = strstr(p,offsetsstrm_suffix)) != NULL && !strcmp(q,offsetsstrm_suffix)) {

	if (q - p == KMER_SAMPLING) {
	  /* New style, e.g., ref153 */
	  if (sscanf(p,"%c%c%c",&tens,&ones,&interval_char) == 3) {
	    digit_string[0] = tens;
	    found_index1part = 10*atoi(digit_string);
	    digit_string[0] = ones;
	    found_index1part += atoi(digit_string);

	    digit_string[0] = interval_char;
	    found_interval = atoi(digit_string);
	  } else {
	    abort();
	  }

	} else {
	  fprintf(stderr,"Cannot parse part between %s and offsets in filename %s: found %ld characters, expecting %d\n",
		  idx_filesuffix,filename,q-p,BASE_KMER_SAMPLING);
	  if (snps_root != NULL) {
	    FREE(offsetsstrm_suffix);
	    FREE(offsetsmeta_suffix);
	    FREE(offsetspages_suffix);
	    FREE(positions_high_suffix);
	    FREE(positions_low_suffix);
	  }
	  return (Filenames_T) NULL;
	}

	if ((required_index1part == 0 || found_index1part == required_index1part) &&
	    (required_interval == 0 || found_interval == required_interval)) {
	  if (required_index1part == 0 && found_index1part < *index1part) {
	    /* Skip, since we have already found a larger index1part */
	  } else if (required_interval == 0 && found_interval > *index1interval) {
	    /* Skip, since we have already found a smaller interval */
	  } else {
	    *index1part = found_index1part;
	    *index1interval = found_interval;
	    FREE(base_filename);
	    base_filename = (char *) CALLOC(strlen(filename)+1,sizeof(char));
	    strcpy(base_filename,filename);
	  }
	}
      }
    }
  }

  FREE(pattern);
#endif


  if (closedir(dp) < 0) {
    fprintf(stderr,"Unable to close directory %s\n",genomesubdir);
  }

  /* Construct full filenames */
  if (base_filename == NULL) {
    /* offsetspages_filename = (char *) NULL; */
    /* offsetsmeta_filename = (char *) NULL; */
    /* offsetsstrm_filename = (char *) NULL; */
    /* positions_high_filename = (char *) NULL; */
    /* positions_low_filename = (char *) NULL; */
    if (snps_root != NULL) {
      FREE(offsetsstrm_suffix);
      FREE(offsetsmeta_suffix);
      FREE(offsetspages_suffix);
      FREE(positions_high_suffix);
      FREE(positions_low_suffix);
    }
    return (Filenames_T) NULL;

  } else {
    offsetsstrm_filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(base_filename)+1,sizeof(char));
    offsetsstrm_basename_ptr = &(offsetsstrm_filename[strlen(genomesubdir)+strlen("/")]);
    offsetsstrm_index1info_ptr = &(offsetsstrm_basename_ptr[patternlength]);

    sprintf(offsetsstrm_filename,"%s/%s",genomesubdir,base_filename);
    if (Access_file_exists_p(offsetsstrm_filename) == false) {
      fprintf(stderr,"Offsets filename %s does not exist\n",offsetsstrm_filename);
      FREE(offsetsstrm_filename);
      /* offsetsstrm_filename = (char *) NULL; */
      /* positions_high_filename = (char *) NULL; */
      /* positions_low_filename = (char *) NULL; */
      FREE(base_filename);
      if (snps_root != NULL) {
	FREE(offsetsstrm_suffix);
	FREE(offsetsmeta_suffix);
	FREE(offsetspages_suffix);
	FREE(positions_high_suffix);
	FREE(positions_low_suffix);
      }
      return (Filenames_T) NULL;
    }


    if ((q = strstr(base_filename,offsetsstrm_suffix)) == NULL) {
      abort();
    } else {
      rootlength = q - base_filename;
    }

    offsetsmeta_filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+rootlength+strlen(offsetsmeta_suffix)+1,sizeof(char));
    offsetsmeta_basename_ptr = &(offsetsmeta_filename[strlen(genomesubdir)+strlen("/")]);
    offsetsmeta_index1info_ptr = &(offsetsmeta_basename_ptr[patternlength]);

    sprintf(offsetsmeta_filename,"%s/",genomesubdir);
    strncpy(offsetsmeta_basename_ptr,base_filename,rootlength);
    strcpy(&(offsetsmeta_basename_ptr[rootlength]),offsetsmeta_suffix);

    if (Access_file_exists_p(offsetsmeta_filename) == false) {
      fprintf(stderr,"Offsetsmeta filename %s does not exist\n",offsetsmeta_filename);
      FREE(offsetsstrm_filename);
      /* offsetsstrm_filename = (char *) NULL; */
      FREE(base_filename);
      if (snps_root != NULL) {
	FREE(offsetsstrm_suffix);
	FREE(offsetsmeta_suffix);
	FREE(positions_high_suffix);
	FREE(positions_low_suffix);
      }
      return (Filenames_T) NULL;
    }


    offsetspages_filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+rootlength+strlen(offsetspages_suffix)+1,sizeof(char));
    offsetspages_basename_ptr = &(offsetspages_filename[strlen(genomesubdir)+strlen("/")]);
    offsetspages_index1info_ptr = &(offsetspages_basename_ptr[patternlength]);

    sprintf(offsetspages_filename,"%s/",genomesubdir);
    strncpy(offsetspages_basename_ptr,base_filename,rootlength);
    strcpy(&(offsetspages_basename_ptr[rootlength]),offsetspages_suffix);
    if (Access_file_exists_p(offsetspages_filename) == false) {
      /* Not a huge genome */
      FREE(offsetspages_filename);
      offsetspages_filename = (char *) NULL;
    }

    positions_high_filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+rootlength+strlen(positions_high_suffix)+1,sizeof(char));
    positions_high_basename_ptr = &(positions_high_filename[strlen(genomesubdir)+strlen("/")]);
    positions_high_index1info_ptr = &(positions_high_basename_ptr[patternlength]);

    positions_low_filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+rootlength+strlen(positions_low_suffix)+1,sizeof(char));
    positions_low_basename_ptr = &(positions_low_filename[strlen(genomesubdir)+strlen("/")]);
    positions_low_index1info_ptr = &(positions_low_basename_ptr[patternlength]);

    sprintf(positions_high_filename,"%s/",genomesubdir);
    strncpy(positions_high_basename_ptr,base_filename,rootlength);
    strcpy(&(positions_high_basename_ptr[rootlength]),positions_high_suffix);

    sprintf(positions_low_filename,"%s/",genomesubdir);
    strncpy(positions_low_basename_ptr,base_filename,rootlength);
    strcpy(&(positions_low_basename_ptr[rootlength]),positions_low_suffix);

    if (offsets_only_p == true) {
      /* Do not look for a positions file */
    } else if (Access_file_exists_p(positions_low_filename) == false) {
      /* Try newer naming scheme: ref153positions instead of ref12153positions */
      sprintf(positions_high_filename,"%s/",genomesubdir);
      strncpy(positions_high_basename_ptr,base_filename,rootlength-BASE_KMER_SAMPLING); /* e.g., skip "12153" */
      strncpy(&(positions_high_basename_ptr[rootlength-BASE_KMER_SAMPLING]),&(base_filename[rootlength-KMER_SAMPLING]),KMER_SAMPLING);
      strcpy(&(positions_high_basename_ptr[rootlength+KMER_SAMPLING-BASE_KMER_SAMPLING]),positions_high_suffix);

      sprintf(positions_low_filename,"%s/",genomesubdir);
      strncpy(positions_low_basename_ptr,base_filename,rootlength-BASE_KMER_SAMPLING); /* e.g., skip "12153" */
      strncpy(&(positions_low_basename_ptr[rootlength-BASE_KMER_SAMPLING]),&(base_filename[rootlength-KMER_SAMPLING]),KMER_SAMPLING);
      strcpy(&(positions_low_basename_ptr[rootlength+KMER_SAMPLING-BASE_KMER_SAMPLING]),positions_low_suffix);

      if (Access_file_exists_p(positions_low_filename) == false) {
	fprintf(stderr,"Positions filename %s does not exist\n",positions_low_filename);
	FREE(offsetspages_filename);
	FREE(offsetsmeta_filename);
	FREE(offsetsstrm_filename);
	FREE(positions_high_filename);
	FREE(positions_low_filename);
	/* offsetsmeta_filename = (char *) NULL; */
	/* offsetsstrm_filename = (char *) NULL; */
	/* positions_high_filename = (char *) NULL; */
	/* positions_low_filename = (char *) NULL; */
	FREE(base_filename);
	if (snps_root != NULL) {
	  FREE(offsetsstrm_suffix);
	  FREE(offsetsmeta_suffix);
	  FREE(offsetspages_suffix);
	  FREE(positions_high_suffix);
	  FREE(positions_low_suffix);
	}
	return (Filenames_T) NULL;
      }
    }

    if (Access_file_exists_p(positions_high_filename) == false) {
      /* Not a large genome */
      FREE(positions_high_filename);
      positions_high_filename = (char *) NULL;
    }

    if (snps_root != NULL) {
      FREE(offsetsstrm_suffix);
      FREE(offsetsmeta_suffix);
      FREE(offsetspages_suffix);
      FREE(positions_high_suffix);
      FREE(positions_low_suffix);
    }

    FREE(base_filename);

    fprintf(stderr,"Looking for index files in directory %s\n",genomesubdir);
    if (offsetspages_filename != NULL)  {
      fprintf(stderr,"  Pages file is %s\n",offsetspages_basename_ptr);
    }
    fprintf(stderr,"  Pointers file is %s\n",offsetsmeta_basename_ptr);
    fprintf(stderr,"  Offsets file is %s\n",offsetsstrm_basename_ptr);
    if (positions_high_filename == NULL) {
      fprintf(stderr,"  Positions file is %s\n",positions_low_basename_ptr);
    } else {
      fprintf(stderr,"  Positions files are %s and %s\n",positions_high_basename_ptr,positions_low_basename_ptr);
    }
    return Filenames_new(offsetspages_filename,offsetsmeta_filename,offsetsstrm_filename,
			 positions_high_filename,positions_low_filename,
			 offsetsmeta_basename_ptr,offsetsstrm_basename_ptr,
			 positions_high_basename_ptr,positions_low_basename_ptr,
			 offsetsmeta_index1info_ptr,offsetsstrm_index1info_ptr,
			 positions_high_index1info_ptr,positions_low_index1info_ptr);

  }
}


Filenames_T
Indexdb_get_filenames (int *compression_type,
#ifdef PMAP
		       Alphabet_T *alphabet, Alphabet_T required_alphabet,
#endif
		       Width_T *index1part, Width_T *index1interval, char *genomesubdir,
		       char *fileroot, char *idx_filesuffix, char *snps_root,
		       Width_T required_index1part, Width_T required_interval,
		       bool offsets_only_p) {
  Filenames_T filenames;

  if ((filenames = Indexdb_get_filenames_no_compression(&(*index1part),&(*index1interval),
							genomesubdir,fileroot,idx_filesuffix,snps_root,
							required_interval,offsets_only_p)) != NULL) {
    *compression_type = NO_COMPRESSION;
    return filenames;
    

  } else if ((filenames = Indexdb_get_filenames_bitpack(
#ifdef PMAP
							&(*alphabet),required_alphabet,
#endif
							&(*index1part),&(*index1interval),
							genomesubdir,fileroot,idx_filesuffix,snps_root,
							required_index1part,required_interval,
							/*blocksize*/64,offsets_only_p)) != NULL) {
    *compression_type = BITPACK64_COMPRESSION;
    return filenames;
    
  } else {
    return (Filenames_T) NULL;
  }
}


void
Indexdb_shmem_remove (char *genomesubdir, char *fileroot, char *idx_filesuffix, char *snps_root,
#ifdef PMAP
		      Alphabet_T *alphabet, int *alphabet_size, Alphabet_T required_alphabet,
#endif
		      Width_T required_index1part, Width_T required_interval, bool expand_offsets_p) {
  Filenames_T filenames;
  int index1part, index1interval;

  if ((filenames = Indexdb_get_filenames_no_compression(&index1part,&index1interval,
							genomesubdir,fileroot,idx_filesuffix,snps_root,
							required_interval,/*offsets_only_p*/false)) != NULL) {
    /* Try non-compressed files */
    Access_shmem_remove(filenames->offsets_filename);

  } else if ((filenames = Indexdb_get_filenames_bitpack(
#ifdef PMAP
							&(*alphabet),required_alphabet,
#endif
							&index1part,&index1interval,
							genomesubdir,fileroot,idx_filesuffix,snps_root,
							required_index1part,required_interval,
							/*blocksize*/64,/*offsets_only_p*/false)) != NULL) {
    if (expand_offsets_p == true) {
      /* ALLOCATED_PRIVATE */

    } else {
      Access_shmem_remove(filenames->pointers_filename);
      Access_shmem_remove(filenames->offsets_filename);
#ifdef LARGE_GENOMES
      if (filenames->pages_filename != NULL) {
	Access_shmem_remove(filenames->pages_filename);
      }
#endif
    }
  }

#ifdef LARGE_GENOMES
  Access_shmem_remove(filenames->positions_high_filename);
  Access_shmem_remove(filenames->positions_low_filename);
#else
  Access_shmem_remove(filenames->positions_low_filename);
#endif

  Filenames_free(&filenames);

  return;
}


T
Indexdb_new_genome (Width_T *index1part, Width_T *index1interval,
		    char *genomesubdir, char *fileroot, char *idx_filesuffix, char *snps_root,
#ifdef PMAP
		    Alphabet_T *alphabet, int *alphabet_size, Alphabet_T required_alphabet,
#endif
		    Width_T required_index1part, Width_T required_interval, bool expand_offsets_p,
		    Access_mode_T offsetsstrm_access, Access_mode_T positions_access, bool sharedp) {
  T new = (T) MALLOC(sizeof(*new));
  Filenames_T filenames;
  Oligospace_T basespace, base;

  unsigned int poly_T;
  Positionsptr_T ptr0, end0;	/* UINT8 or UINT4 */
  off_t filesize;

#ifdef LARGE_GENOMES
  size_t offsetspages_len;
#endif

  char *comma;
  double seconds;
#ifdef HAVE_MMAP
  int npages;
#endif

  if ((filenames = Indexdb_get_filenames_no_compression(&new->index1part,&new->index1interval,
							genomesubdir,fileroot,idx_filesuffix,snps_root,
							required_interval,/*offsets_only_p*/false)) != NULL) {
    /* Try non-compressed files */
    fprintf(stderr,"Offsets compression type: none\n");
    new->compression_type = NO_COMPRESSION;
    new->blocksize = 1;

    *index1part = new->index1part;
    *index1interval = new->index1interval;

#ifdef PMAP
    basespace = power(*alphabet_size,new->index1part);
#else
    basespace = power(4,new->index1part);
#endif
    new->offsetsmeta = (UINT4 *) CALLOC(basespace+1,sizeof(UINT4));
    for (base = 0; base <= basespace; base++) {
      new->offsetsmeta[base] = base;
    }
    new->offsetsmeta_access = ALLOCATED_PRIVATE;

    if (offsetsstrm_access == USE_ALLOCATE) {
      if (snps_root) {
	fprintf(stderr,"Allocating memory for %s (%s) offsets, kmer %d, interval %d...",
		idx_filesuffix,snps_root,new->index1part,new->index1interval);
      } else {
	fprintf(stderr,"Allocating memory for %s offsets, kmer %d, interval %d...",
		idx_filesuffix,new->index1part,new->index1interval);
      }
      new->offsetsstrm = (UINT4 *) Access_allocate(&new->offsetsstrm_shmid,&new->offsetsstrm_len,&seconds,
						   filenames->offsets_filename,sizeof(UINT4),sharedp);
      if (new->offsetsstrm == NULL) {
	fprintf(stderr,"insufficient memory (need to use a lower batch mode (-B))\n");
	exit(9);
      } else {
	comma = Genomicpos_commafmt(new->offsetsstrm_len);
	fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma,seconds);
	FREE(comma);
	if (sharedp == true) {
	  new->offsetsstrm_access = ALLOCATED_SHARED;
	} else {
	  new->offsetsstrm_access = ALLOCATED_PRIVATE;
	}
      }

#ifdef HAVE_MMAP
    } else if (offsetsstrm_access == USE_MMAP_PRELOAD) {
      if (snps_root) {
	fprintf(stderr,"Pre-loading %s (%s) offsets, kmer %d, interval %d...",
		idx_filesuffix,snps_root,new->index1part,new->index1interval);
      } else {
	fprintf(stderr,"Pre-loading %s offsets, kmer %d, interval %d...",
		idx_filesuffix,new->index1part,new->index1interval);
      }
      new->offsetsstrm = (UINT4 *) Access_mmap_and_preload(&new->offsetsstrm_fd,&new->offsetsstrm_len,&npages,&seconds,
								   filenames->offsets_filename,sizeof(UINT4));
      if (new->offsetsstrm == NULL) {
	fprintf(stderr,"insufficient memory (will use disk file instead, but program may not run)\n");
#ifdef PMAP
	new->offsetsstrm_access = FILEIO;
#else
	exit(9);
#endif
      } else {
	comma = Genomicpos_commafmt(new->offsetsstrm_len);
	fprintf(stderr,"done (%s bytes, %d pages, %.2f sec)\n",comma,npages,seconds);
	FREE(comma);
	new->offsetsstrm_access = MMAPPED;
      }

    } else if (offsetsstrm_access == USE_MMAP_ONLY) {
      new->offsetsstrm = (UINT4 *) Access_mmap(&new->offsetsstrm_fd,&new->offsetsstrm_len,
						       filenames->offsets_filename,sizeof(UINT4),/*randomp*/false);
      if (new->offsetsstrm == NULL) {
	fprintf(stderr,"Insufficient memory for mmap of %s (will use disk file instead, but program may not run)\n",
		filenames->offsets_filename);
#ifdef PMAP
	new->offsetsstrm_access = FILEIO;
#else
	exit(9);
#endif
      } else {
	new->offsetsstrm_access = MMAPPED;
      }
#endif

    } else if (offsetsstrm_access == USE_FILEIO) {
      fprintf(stderr,"Offsets file I/O access of %s not allowed\n",filenames->offsets_filename);
      exit(9);

    } else {
      fprintf(stderr,"Don't recognize offsetsstrm_access type %d\n",offsetsstrm_access);
      abort();
    }

 } else if ((filenames = Indexdb_get_filenames_bitpack(
#ifdef PMAP
						       &(*alphabet),required_alphabet,
#endif
						       &new->index1part,&new->index1interval,
						       genomesubdir,fileroot,idx_filesuffix,snps_root,
						       required_index1part,required_interval,
						       /*blocksize*/64,/*offsets_only_p*/false)) != NULL) {
    /* Try bitpack compression  */
    *index1part = new->index1part;
    *index1interval = new->index1interval;

    if (expand_offsets_p == true) {
#ifdef LARGE_GENOMES
      fprintf(stderr,"Expansion of bitpack offsets not supported for large genomes\n");
      abort();
#else      
      fprintf(stderr,"Offsets compression type: none (bitpack expanded)\n");
      new->compression_type = NO_COMPRESSION;

#ifdef PMAP
      *alphabet_size = Alphabet_get_size(*alphabet);
      new->blocksize = 1;
      basespace = power(*alphabet_size,new->index1part);
#else
      new->blocksize = 1;
      basespace = power(4,new->index1part);
#endif
      new->offsetsmeta = (UINT4 *) CALLOC(basespace+1,sizeof(UINT4));
      for (base = 0; base <= basespace; base++) {
	new->offsetsmeta[base] = base;
      }
      new->offsetsmeta_access = ALLOCATED_PRIVATE;

#ifdef PMAP
      new->offsetsstrm = Indexdb_offsets_from_bitpack(filenames->pointers_filename,filenames->offsets_filename,
						      *alphabet_size,new->index1part);
#else
      new->offsetsstrm = Indexdb_offsets_from_bitpack(filenames->pointers_filename,filenames->offsets_filename,
						      new->index1part);
#endif
      new->offsetsstrm_access = ALLOCATED_PRIVATE;

#endif

#ifdef PMAP
#else
      /* Sanity check on positions filesize */
      poly_T = ~(~0UL << 2*new->index1part);
      end0 = new->offsetsstrm[poly_T+1];
#ifdef LARGE_GENOMES
      if ((filesize = Access_filesize(filenames->positions_high_filename)) != end0 * (off_t) sizeof(unsigned char)) {
	fprintf(stderr,"Something is wrong with the genomic index: expected file size for %s is %zu, but observed %zu.\n",
		filenames->positions_high_filename,end0*sizeof(unsigned char),filesize);
	abort();
      }
#endif
      if ((filesize = Access_filesize(filenames->positions_low_filename)) != end0 * (off_t) sizeof(UINT4)) {
	fprintf(stderr,"Something is wrong with the genomic index: expected file size for %s is %zu, but observed %zu.\n",
		filenames->positions_low_filename,end0*sizeof(UINT4),filesize);
	abort();
      }
#endif	/* PMAP */

    } else {
      fprintf(stderr,"Offsets compression type: bitpack64\n");
      new->compression_type = BITPACK64_COMPRESSION;

#ifdef PMAP
      *alphabet_size = Alphabet_get_size(*alphabet);
      /* new->blocksize = power(*alphabet_size,(*index1part) - new->offsetsstrm_basesize); */
#else
      new->blocksize = 64;  /* power(4,(*index1part) - new->offsetsstrm_basesize); */
#endif

      /* offsetsmeta and offsetsmeta always ALLOCATED */
      if (snps_root) {
	fprintf(stderr,"Allocating memory for %s (%s) offset pointers, kmer %d, interval %d...",
		idx_filesuffix,snps_root,new->index1part,new->index1interval);
      } else {
	fprintf(stderr,"Allocating memory for %s offset pointers, kmer %d, interval %d...",
		idx_filesuffix,new->index1part,new->index1interval);
      }
      new->offsetsmeta = (UINT4 *) Access_allocate(&new->offsetsmeta_shmid,&new->offsetsmeta_len,&seconds,
						   filenames->pointers_filename,sizeof(UINT4),sharedp);
      if (sharedp == true) {
	new->offsetsmeta_access = ALLOCATED_SHARED;
      } else {
	new->offsetsmeta_access = ALLOCATED_PRIVATE;
      }

      comma = Genomicpos_commafmt(new->offsetsmeta_len);
      fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma,seconds);
      FREE(comma);

      /* offsetsstrm could be ALLOCATED or MMAPPED +/- PRELOAD */
      if (offsetsstrm_access == USE_ALLOCATE) {
	if (snps_root) {
	  fprintf(stderr,"Allocating memory for %s (%s) offsets, kmer %d, interval %d...",
		  idx_filesuffix,snps_root,new->index1part,new->index1interval);
	} else {
	  fprintf(stderr,"Allocating memory for %s offsets, kmer %d, interval %d...",
		  idx_filesuffix,new->index1part,new->index1interval);
	}
	new->offsetsstrm = (UINT4 *) Access_allocate(&new->offsetsstrm_shmid,&new->offsetsstrm_len,&seconds,
						     filenames->offsets_filename,sizeof(UINT4),sharedp);
	if (new->offsetsstrm == NULL) {
	  fprintf(stderr,"insufficient memory (need to use a lower batch mode (-B))\n");
	  exit(9);
	} else {
	  comma = Genomicpos_commafmt(new->offsetsstrm_len);
	  fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma,seconds);
	  FREE(comma);
	  if (sharedp == true) {
	    new->offsetsstrm_access = ALLOCATED_SHARED;
	  } else {
	    new->offsetsstrm_access = ALLOCATED_PRIVATE;
	  }
	}

#ifdef HAVE_MMAP
      } else if (offsetsstrm_access == USE_MMAP_PRELOAD) {
	if (snps_root) {
	  fprintf(stderr,"Pre-loading %s (%s) offsets, kmer %d, interval %d...",
		  idx_filesuffix,snps_root,new->index1part,new->index1interval);
	} else {
	  fprintf(stderr,"Pre-loading %s offsets, kmer %d, interval %d...",
		  idx_filesuffix,new->index1part,new->index1interval);
	}
	new->offsetsstrm = (UINT4 *) Access_mmap_and_preload(&new->offsetsstrm_fd,&new->offsetsstrm_len,&npages,&seconds,
							     filenames->offsets_filename,sizeof(UINT4));
	if (new->offsetsstrm == NULL) {
	  fprintf(stderr,"insufficient memory (will use disk file instead, but program may not run)\n");
#ifdef PMAP
	  new->offsetsstrm_access = FILEIO;
#else
	  exit(9);
#endif
	} else {
	  comma = Genomicpos_commafmt(new->offsetsstrm_len);
	  fprintf(stderr,"done (%s bytes, %d pages, %.2f sec)\n",comma,npages,seconds);
	  FREE(comma);
	  new->offsetsstrm_access = MMAPPED;
	}

      } else if (offsetsstrm_access == USE_MMAP_ONLY) {
	new->offsetsstrm = (UINT4 *) Access_mmap(&new->offsetsstrm_fd,&new->offsetsstrm_len,
						 filenames->offsets_filename,sizeof(UINT4),/*randomp*/false);
	if (new->offsetsstrm == NULL) {
	  fprintf(stderr,"Insufficient memory for mmap of %s (will use disk file instead, but program may not run)\n",
		  filenames->offsets_filename);
#ifdef PMAP
	  new->offsetsstrm_access = FILEIO;
#else
	  exit(9);
#endif
	} else {
	  new->offsetsstrm_access = MMAPPED;
	}
#endif

      } else if (offsetsstrm_access == USE_FILEIO) {
	fprintf(stderr,"Offsetsstrm file I/O access of %s not allowed\n",filenames->offsets_filename);
	exit(9);

      } else {
	fprintf(stderr,"Don't recognize offsetsstrm_access type %d\n",offsetsstrm_access);
	abort();
      }

#ifdef PMAP
#else
      /* Sanity check on positions filesize */
#ifdef LARGE_GENOMES
      if (filenames->pages_filename != NULL) {
	new->offsetspages = (UINT4 *) Access_allocate(&new->offsetspages_shmid,&offsetspages_len,&seconds,filenames->pages_filename,sizeof(UINT4),sharedp);
	if (sharedp == true) {
	  new->offsetspages_access = ALLOCATED_SHARED;
	} else {
	  new->offsetspages_access = ALLOCATED_PRIVATE;
	}
      } else {
	new->offsetspages_access = ALLOCATED_PRIVATE;
	new->offsetspages = (UINT4 *) MALLOC(1*sizeof(UINT4));
	new->offsetspages[0] = -1U;
      }
#endif
      poly_T = ~(~0UL << 2*new->index1part);
#ifdef LARGE_GENOMES
      ptr0 = Bitpack64_read_two_huge(&end0,poly_T,new->offsetspages,new->offsetsmeta,new->offsetsstrm);
#else
      ptr0 = Bitpack64_read_two(&end0,poly_T,new->offsetsmeta,new->offsetsstrm);
#endif

#ifdef LARGE_GENOMES
      ptr0 = Bitpack64_read_two_huge(&end0,poly_T,new->offsetspages,new->offsetsmeta,new->offsetsstrm);
      if ((filesize = Access_filesize(filenames->positions_high_filename)) != end0 * (off_t) sizeof(unsigned char)) {
	fprintf(stderr,"Something is wrong with the genomic index: expected file size for %s is %zu, but observed %zu.\n",
		filenames->positions_high_filename,end0*sizeof(unsigned char),filesize);
	abort();
      }
#endif
      if ((filesize = Access_filesize(filenames->positions_low_filename)) != end0 * (off_t) sizeof(UINT4)) {
	fprintf(stderr,"Something is wrong with the genomic index: expected file size for %s is %zu, but observed %zu.\n",
		filenames->positions_low_filename,end0*sizeof(UINT4),filesize);
	abort();
      }
#endif	/* PMAP */

    }

  } else {
    fprintf(stderr,"Cannot find genomic index files in either current or old format.  Looking for files containing %s\n",idx_filesuffix);
    exit(9);
  }


  /* Read or memory map positions file */
  if (positions_access == USE_ALLOCATE) {
    if (snps_root) {
      fprintf(stderr,"Allocating memory for %s (%s) positions, kmer %d, interval %d...",
	      idx_filesuffix,snps_root,new->index1part,new->index1interval);
    } else {
      fprintf(stderr,"Allocating memory for %s positions, kmer %d, interval %d...",
	      idx_filesuffix,new->index1part,new->index1interval);
    }
#ifdef LARGE_GENOMES
    new->positions_high = (unsigned char *) Access_allocate(&new->positions_high_shmid,&new->positions_high_len,&seconds,
							    filenames->positions_high_filename,sizeof(unsigned char),sharedp);
    if (new->positions_high == NULL) {
      fprintf(stderr,"insufficient memory (need to use a lower batch mode (-B)\n");
      exit(9);
    } else {
      comma = Genomicpos_commafmt(new->positions_high_len);
      fprintf(stderr,"done (%s bytes, %.2f sec), ",comma,seconds);
      FREE(comma);

      new->positions_low = (UINT4 *) Access_allocate(&new->positions_low_shmid,&new->positions_low_len,&seconds,
						     filenames->positions_low_filename,sizeof(UINT4),sharedp);
      if (new->positions_low == NULL) {
	fprintf(stderr,"insufficient memory (need to use a lower batch mode (-B)\n");
	exit(9);
      } else {
	comma = Genomicpos_commafmt(new->positions_low_len);
	fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma,seconds);
	FREE(comma);

	if (sharedp == true) {
	  new->positions_access = ALLOCATED_SHARED;
	} else {
	  new->positions_access = ALLOCATED_PRIVATE;
	}
      }
    }
#else
    new->positions = (UINT4 *) Access_allocate(&new->positions_shmid,&new->positions_len,&seconds,
					       filenames->positions_low_filename,sizeof(UINT4),sharedp);
    if (new->positions == NULL) {
      fprintf(stderr,"insufficient memory (need to use a lower batch mode (-B)\n");
      exit(9);
    } else {
      comma = Genomicpos_commafmt(new->positions_len);
      fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma,seconds);
      FREE(comma);
      if (sharedp == true) {
	new->positions_access = ALLOCATED_SHARED;
      } else {
	new->positions_access = ALLOCATED_PRIVATE;
      }
    }
#endif
    

#ifdef HAVE_MMAP
  } else if (positions_access == USE_MMAP_PRELOAD) {
    if (snps_root) {
      fprintf(stderr,"Pre-loading %s (%s) positions, kmer %d, interval %d...",
	      idx_filesuffix,snps_root,new->index1part,new->index1interval);
    } else {
      fprintf(stderr,"Pre-loading %s positions, kmer %d, interval %d...",
	      idx_filesuffix,new->index1part,new->index1interval);
    }
#ifdef LARGE_GENOMES
    new->positions_high = (unsigned char *) Access_mmap_and_preload(&new->positions_high_fd,&new->positions_high_len,&npages,&seconds,
								  filenames->positions_high_filename,sizeof(unsigned char));
    if (new->positions_high == NULL) {
      fprintf(stderr,"insufficient memory (will use disk file instead, but program will be slow)\n");
      new->positions_access = FILEIO;
    } else {
      comma = Genomicpos_commafmt(new->positions_high_len);
      fprintf(stderr,"done (%s bytes, %d pages, %.2f sec), ",comma,npages,seconds);
      FREE(comma);

      new->positions_low = (UINT4 *) Access_mmap_and_preload(&new->positions_low_fd,&new->positions_low_len,&npages,&seconds,
							     filenames->positions_low_filename,sizeof(UINT4));
      if (new->positions_low == NULL) {
	fprintf(stderr,"insufficient memory (will use disk file instead, but program will be slow)\n");
	new->positions_access = FILEIO;
      } else {
	comma = Genomicpos_commafmt(new->positions_low_len);
	fprintf(stderr,"done (%s bytes, %d pages, %.2f sec)\n",comma,npages,seconds);
	FREE(comma);

	new->positions_access = MMAPPED;
      }
    }
#else
    new->positions = (UINT4 *) Access_mmap_and_preload(&new->positions_fd,&new->positions_len,&npages,&seconds,
						       filenames->positions_low_filename,sizeof(UINT4));
    if (new->positions == NULL) {
      fprintf(stderr,"insufficient memory (will use disk file instead, but program will be slow)\n");
      new->positions_access = FILEIO;
    } else {
      comma = Genomicpos_commafmt(new->positions_len);
      fprintf(stderr,"done (%s bytes, %d pages, %.2f sec)\n",comma,npages,seconds);
      FREE(comma);
      new->positions_access = MMAPPED;
    }
#endif

  } else if (positions_access == USE_MMAP_ONLY) {
#ifdef LARGE_GENOMES
    new->positions_high = (unsigned char *) Access_mmap(&new->positions_high_fd,&new->positions_high_len,
							filenames->positions_high_filename,sizeof(unsigned char),/*randomp*/true);
    if (new->positions_high == NULL) {
      fprintf(stderr,"Insufficient memory for mmap of %s (will use disk file instead, but program will be slow)\n",
	      filenames->positions_high_filename);
      new->positions_access = FILEIO;
    } else {
      new->positions_low = (UINT4 *) Access_mmap(&new->positions_low_fd,&new->positions_low_len,
						 filenames->positions_low_filename,sizeof(UINT4),/*randomp*/true);
      if (new->positions_low == NULL) {
	fprintf(stderr,"Insufficient memory for mmap of %s (will use disk file instead, but program will be slow)\n",
		filenames->positions_low_filename);
	new->positions_access = FILEIO;
      } else {
	new->positions_access = MMAPPED;
      }
    }
#else
    new->positions = (UINT4 *) Access_mmap(&new->positions_fd,&new->positions_len,
					   filenames->positions_low_filename,sizeof(UINT4),/*randomp*/true);

    if (new->positions == NULL) {
      fprintf(stderr,"Insufficient memory for mmap of %s (will use disk file instead, but program will be slow)\n",
	      filenames->positions_low_filename);
      new->positions_access = FILEIO;
    } else {
      new->positions_access = MMAPPED;
    }
#endif

#endif

  } else if (positions_access == USE_FILEIO) {
    new->positions_access = FILEIO;
  } else {
    fprintf(stderr,"Don't recognize positions_access %d\n",positions_access);
    abort();
  }

#ifdef HAVE_PTHREAD
  if (new->positions_access == FILEIO) {
    pthread_mutex_init(&new->positions_read_mutex,NULL);
  }
#endif

  Filenames_free(&filenames);

  return new;
}


/************************************************************************
 *   Debugging procedures
 ************************************************************************/

#ifndef PMAP

/*                87654321 */
#define RIGHT_A 0x00000000
#define RIGHT_C 0x00000001
#define RIGHT_G 0x00000002
#define RIGHT_T 0x00000003

/*                      87654321 */
#define LOW_TWO_BITS  0x00000003

#if (defined(DEBUG0) || defined(DEBUG1) || defined(DEBUG2))
static char *
shortoligo_nt (Storedoligomer_T oligo, Width_T oligosize) {
  char *nt;
  Width_T i, j;
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

#endif


/************************************************************************
 *   Read procedures
 ************************************************************************/


static void
positions_move_absolute_1 (int positions_fd, off_t ptr) {
  off_t offset = ptr*((off_t) sizeof(unsigned char));

  if (lseek(positions_fd,offset,SEEK_SET) < 0) {
    fprintf(stderr,"Attempted to do lseek on offset %jd*%d=%jd\n",
	    ptr,(int) sizeof(unsigned char),offset);
    perror("Error in indexdb.c, positions_move_absolute");
    exit(9);
  }
  return;
}

static void
positions_move_absolute_4 (int positions_fd, off_t ptr) {
  off_t offset = ptr*((off_t) sizeof(UINT4));

  if (lseek(positions_fd,offset,SEEK_SET) < 0) {
    fprintf(stderr,"Attempted to do lseek on offset %jd*%d=%jd\n",
	    ptr,(int) sizeof(UINT4),offset);
    perror("Error in indexdb.c, positions_move_absolute");
    exit(9);
  }
  return;
}

#if 0
static Univcoord_T
positions_read_forward (int positions_fd) {
  Univcoord_T value;
  char buffer[4];

  read(positions_fd,buffer,4);

  value = (buffer[3] & 0xff);
  value <<= 8;
  value |= (buffer[2] & 0xff);
  value <<= 8;
  value |= (buffer[1] & 0xff);
  value <<= 8;
  value |= (buffer[0] & 0xff);

  return value;
}
#endif

#ifdef LARGE_GENOMES
static void
positions_read_multiple_large (int positions_high_fd, int positions_low_fd, Univcoord_T *values, int n) {
  int i;
  Univcoord_T value;
  unsigned char buffer[4];

#ifdef WORDS_BIGENDIAN
  /* Need to keep in bigendian format */
  for (i = 0; i < n; i++) {
    read(positions_high_fd,buffer,1);
    value = (buffer[0] & 0xff);
    value <<= 8;

    read(positions_low_fd,buffer,4);
    value |= (buffer[0] & 0xff);
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
    read(positions_high_fd,buffer,1);
    value = (buffer[0] & 0xff);
    value <<= 8;

    read(positions_low_fd,buffer,4);
    value |= (buffer[3] & 0xff);
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

static void
positions_copy_multiple_large (Univcoord_T *positions, unsigned char *positions_high, UINT4 *positions_low, int n) {
  int i;

  for (i = 0; i < n; i++) {
    positions[i] = ((Univcoord_T) positions_high[i] << 32) + positions_low[i];
  }

  return;
}

#else

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




#if 0
static Univcoord_T
positions_read_backward (int positions_fd) {
  Univcoord_T value;
  char buffer[4];
  off_t reloffset = -2*((off_t) sizeof(Univcoord_T)); /* 1 to undo the effect of read */

  read(positions_fd,buffer,4);

  value = (buffer[3] & 0xff);
  value <<= 8;
  value |= (buffer[2] & 0xff);
  value <<= 8;
  value |= (buffer[1] & 0xff);
  value <<= 8;
  value |= (buffer[0] & 0xff);

  if (lseek(positions_fd,reloffset,SEEK_CUR) < 0) {
    fprintf(stderr,"Attempted to do lseek on relative offset %ld\n",(long int) reloffset);
    perror("Error in indexdb.c, positions_read_backward");
    exit(9);
  }
    
  return value;
}
#endif


/* Used by non-utility programs */
Positionsptr_T *
Indexdb_offsets_from_bitpack (char *offsetsmetafile, char *offsetsstrmfile,
#ifdef PMAP
			      int alphabet_size, Width_T index1part_aa
#else
			      Width_T index1part
#endif
			      ) {
#ifndef HAVE_MMAP
  int shmid;
#endif
  UINT4 *offsetsmeta;
  UINT4 *offsetsstrm;
  int offsetsmeta_fd, offsetsstrm_fd;
  size_t offsetsmeta_len, offsetsstrm_len;
  Positionsptr_T *offsets = NULL;
  Oligospace_T oligospace, oligoi;
#ifndef PMAP
  Blocksize_T blocksize;
#endif
  double seconds;

#ifdef PMAP
  oligospace = power(alphabet_size,index1part_aa);
#else
  oligospace = power(4,index1part);
  blocksize = 64;
#endif


#ifdef HAVE_MMAP
  offsetsmeta = (UINT4 *) Access_mmap(&offsetsmeta_fd,&offsetsmeta_len,offsetsmetafile,sizeof(UINT4),/*randomp*/false);
  offsetsstrm = (UINT4 *) Access_mmap(&offsetsstrm_fd,&offsetsstrm_len,offsetsstrmfile,sizeof(UINT4),/*randomp*/false);
#else
  offsetsmeta = (UINT4 *) Access_allocate(&shmid,&offsetsmeta_len,&seconds,offsetsmetafile,sizeof(UINT4),/*sharedp*/false);
  offsetsstrm = (UINT4 *) Access_allocate(&shmid,&offsetsstrm_len,&seconds,offsetsstrmfile,sizeof(UINT4),/*sharedp*/false);
#endif

#ifdef OLIGOSPACE_NOT_LONG
  fprintf(stderr,"Allocating memory (%u 4-byte words) for offsets, kmer %d...",oligospace+1U,
#ifdef PMAP
	  index1part_aa
#else
	  index1part
#endif
	  );
#else
  fprintf(stderr,"Allocating memory (%llu 4-byte words) for offsets, kmer %d...",(unsigned long long) oligospace+1UL,
#ifdef PMAP
	  index1part_aa
#else
	  index1part
#endif
	  );
#endif

  /* Bitpack procedures start from offsets[1], so we need to print offsets[0] as a special case */
  offsets = (Positionsptr_T *) MALLOC((oligospace+1) * sizeof(Positionsptr_T));

  if (offsets == NULL) {
    fprintf(stderr,"cannot allocated requested memory.  Cannot run expand offsets on this machine.\n");
    exit(9);
  } else {
    fprintf(stderr,"done\n");
  }

  fprintf(stderr,"Expanding offsetsstrm into offsets...");
#ifdef PMAP
  for (oligoi = 0UL; oligoi <= oligospace; oligoi += 1) {
    offsets[oligoi] = Bitpack64_read_one(oligoi,offsetsmeta,offsetsstrm);
  }
#elif defined(LARGE_GENOMES)
  fprintf(stderr,"Cannot do expand offsets on large genomes\n");
  exit(9);
#else
  for (oligoi = 0UL; oligoi < oligospace; oligoi += blocksize) {
    Bitpack64_block_offsets(&(offsets[oligoi]),oligoi,offsetsmeta,offsetsstrm);
  }
#endif

  fprintf(stderr,"done\n");

#ifdef HAVE_MMAP
  munmap((void *) offsetsstrm,offsetsstrm_len);
  munmap((void *) offsetsmeta,offsetsmeta_len);
#else
  FREE(offsetsstrm);
  FREE(offsetsmeta);
#endif

  return offsets;
}


#if defined(HAVE_64_BIT) && defined(UTILITYP)
/* Used by utility programs for writing indexdb */
Hugepositionsptr_T *
Indexdb_offsets_from_bitpack_huge (char *offsetspagesfile, char *offsetsmetafile, char *offsetsstrmfile,
#ifdef PMAP
				   int alphabet_size, Width_T index1part_aa
#else
				   Width_T index1part
#endif
				   ) {
  UINT4 *offsetspages;
  UINT4 *offsetsmeta;
  UINT4 *offsetsstrm;

  int shmid;
  int offsetsmeta_fd, offsetsstrm_fd;
  size_t offsetspages_len, offsetsmeta_len, offsetsstrm_len;
  Hugepositionsptr_T *offsets = NULL;
  Oligospace_T oligospace, oligoi;
  Blocksize_T blocksize;
  double seconds;

#ifdef PMAP
  oligospace = power(alphabet_size,index1part_aa);
  /* blocksize = power(alphabet_size,index1part_aa - offsetsstrm_basesize); */
#else
  oligospace = power(4,index1part);
  blocksize = 64;
#endif

  if (blocksize == 1) {
    return (Hugepositionsptr_T *) Access_allocate(&shmid,&offsetsstrm_len,&seconds,offsetsstrmfile,sizeof(Hugepositionsptr_T),/*sharedp*/false);

  } else {
    if (offsetspagesfile == NULL) {
      offsetspages = (UINT4 *) MALLOC(1*sizeof(UINT4));
      offsetspages[0] = -1U;
    } else {
      offsetspages = (UINT4 *) Access_allocate(&shmid,&offsetspages_len,&seconds,offsetspagesfile,sizeof(UINT4),/*sharedp*/false);
    }
#ifdef HAVE_MMAP
    offsetsmeta = (UINT4 *) Access_mmap(&offsetsmeta_fd,&offsetsmeta_len,offsetsmetafile,sizeof(UINT4),/*randomp*/false);
    offsetsstrm = (UINT4 *) Access_mmap(&offsetsstrm_fd,&offsetsstrm_len,offsetsstrmfile,sizeof(UINT4),/*randomp*/false);
#else
    offsetsmeta = (UINT4 *) Access_allocate(&shmid,&offsetsmeta_len,&seconds,offsetsmetafile,sizeof(UINT4),/*sharedp*/false);
    offsetsstrm = (UINT4 *) Access_allocate(&shmid,&offsetsstrm_len,&seconds,offsetsstrmfile,sizeof(UINT4),/*sharedp*/false);
#endif

#ifdef OLIGOSPACE_NOT_LONG
    fprintf(stderr,"Allocating memory (%u 8-byte words) for offsets, kmer %d...",oligospace+1U,
#ifdef PMAP
	    index1part_aa
#else
	    index1part
#endif
	    );
#else
    fprintf(stderr,"Allocating memory (%llu 8-byte words) for offsets, kmer %d...",(unsigned long long) oligospace+1UL,
#ifdef PMAP
	    index1part_aa
#else
	    index1part
#endif
	    );
#endif

    /* Bitpack procedures start from offsets[1], so we need to print offsets[0] as a special case */
    offsets = (Hugepositionsptr_T *) MALLOC((oligospace+1) * sizeof(Hugepositionsptr_T));

    if (offsets == NULL) {
      fprintf(stderr,"cannot allocated requested memory.  Cannot run expand offsets on this machine.\n");
      exit(9);
    } else {
      fprintf(stderr,"done\n");
    }

    fprintf(stderr,"Expanding offsetsstrm into offsets...");
#ifdef PMAP
    for (oligoi = 0UL; oligoi <= oligospace; oligoi += 1) {
      offsets[oligoi] = Bitpack64_read_one_huge(oligoi,offsetspages,offsetsmeta,offsetsstrm);
    }
#else
    for (oligoi = 0UL; oligoi < oligospace; oligoi += blocksize) {
      Bitpack64_block_offsets_huge(&(offsets[oligoi]),oligoi,offsetspages,offsetsmeta,offsetsstrm);
    }
#endif

    fprintf(stderr,"done\n");

#ifdef HAVE_MMAP
    munmap((void *) offsetsstrm,offsetsstrm_len);
    munmap((void *) offsetsmeta,offsetsmeta_len);
#else
    FREE(offsetsstrm);
    FREE(offsetsmeta);
#endif
    FREE(offsetspages);

    return offsets;
  }
}
#endif


#ifdef PMAP

/* PMAP version.  Doesn't mask bottom 12 nt. */
Univcoord_T *
Indexdb_read (int *nentries, T this, Storedoligomer_T aaindex) {
  Univcoord_T *positions, bigendian, littleendian;
  Positionsptr_T ptr, ptr0, end0;
  int i;
  char byte1, byte2, byte3;

  debug0(printf("%u (%s)\n",aaindex,Alphabet_aaindex_aa(aaindex,this->alphabet)));

  if (this->compression_type == NO_COMPRESSION) {
#ifdef WORDS_BIGENDIAN
    /* Also holds for ALLOCATED_PRIVATE and ALLOCATED_SHARED */
    ptr0 = Bigendian_convert_uint(this->offsetsstrm[aaindex]);
    end0 = Bigendian_convert_uint(this->offsetsstrm[aaindex+1]);
#else
    ptr0 = this->offsetsstrm[aaindex];
    end0 = this->offsetsstrm[aaindex+1];
#endif

  } else if (this->compression_type == BITPACK64_COMPRESSION) {
    ptr0 = Bitpack64_read_two(&end0,aaindex,this->offsetsmeta,this->offsetsstrm);

  } else {
    abort();
  }


  debug0(printf("offset pointers are %u and %u\n",ptr0,end0));

#ifdef ALLOW_DUPLICATES
  /* Not used */
  /* Skip backward over bad values, due to duplicates */
  if (this->positions_access == FILEIO) {
#ifdef HAVE_PTHREAD
    pthread_mutex_lock(&this->positions_read_mutex);
#endif
    positions_move_absolute(this->positions_fd,end0-1);
    while (end0 > ptr0 && positions_read_backward(this->positions_fd) == BADVAL) {
      end0--;
    }
#ifdef HAVE_PTHREAD
    pthread_mutex_unlock(&this->positions_read_mutex);
#endif
  } else {
    while (end0 > ptr0 && this->positions[end0-1] == BADVAL) {
      end0--;
    }
  }
#endif	/* ALLOW_DUPLICATES */

  if ((*nentries = end0 - ptr0) == 0) {
    return NULL;
  } else {
    positions = (Univcoord_T *) CALLOC(*nentries,sizeof(Univcoord_T));
    if (this->positions_access == FILEIO) {
#ifdef HAVE_PTHREAD
      pthread_mutex_lock(&this->positions_read_mutex);
#endif
#ifdef LARGE_GENOMES
      positions_move_absolute_1(this->positions_high_fd,ptr0);
      positions_move_absolute_4(this->positions_low_fd,ptr0);
      positions_read_multiple_large(this->positions_high_fd,this->positions_low_fd,positions,*nentries);
#else
      positions_move_absolute_4(this->positions_fd,ptr0);
      positions_read_multiple(this->positions_fd,positions,*nentries);
#endif
#ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&this->positions_read_mutex);
#endif

    } else if (this->positions_access == ALLOCATED_PRIVATE || this->positions_access == ALLOCATED_SHARED) {
#ifdef LARGE_GENOMES
      positions_copy_multiple_large(positions,&(this->positions_high[ptr0]),&(this->positions_low[ptr0]),*nentries);
#else
      memcpy(positions,&(this->positions[ptr0]),(*nentries)*sizeof(Univcoord_T));
#endif

    } else {
#ifdef WORDS_BIGENDIAN
#ifdef LARGE_GENOMES
      for (ptr = ptr0, i = 0; ptr < end0; ptr++, i++) {
	bigendian = (Univcoord_T) this->positions_high[ptr];
	bigendian <<= 8;

	littleendian = this->positions_low[ptr];
	bigendian |= littleendian & 0xff; /* 0 */
	bigendian <<= 8;
	bigendian |= ((littleendian >>= 8) & 0xff); /* 1 */
	bigendian <<= 8;
	bigendian |= ((littleendian >>= 8) & 0xff); /* 2 */
	bigendian <<= 8;
	bigendian |= ((littleendian >>= 8) & 0xff); /* 3 */
	positions[i] = bigendian;
      }
#else
      for (ptr = ptr0, i = 0; ptr < end0; ptr++, i++) {
	littleendian = this->positions[ptr];
	bigendian = littleendian & 0xff; /* 0 */
	bigendian <<= 8;
	bigendian |= ((littleendian >>= 8) & 0xff); /* 1 */
	bigendian <<= 8;
	bigendian |= ((littleendian >>= 8) & 0xff); /* 2 */
	bigendian <<= 8;
	bigendian |= ((littleendian >>= 8) & 0xff); /* 3 */
	positions[i] = bigendian;
      }
#endif

#else  /* littleendian */

#ifdef LARGE_GENOMES
      positions_copy_multiple_large(positions,&(this->positions_high[ptr0]),&(this->positions_low[ptr0]),*nentries);
#else
      memcpy(positions,&(this->positions[ptr0]),(*nentries)*sizeof(Univcoord_T));
#endif
#endif
    }
    debug0(
	   printf("%d entries:",*nentries);
	   for (i = 0; i < *nentries; i++) {
	     printf(" %u",positions[i]);
	   }
	   printf("\n");
	   );

    return positions;
  }
}

#else

/* GMAP version: Allocates memory */
Univcoord_T *
Indexdb_read (int *nentries, T this, Storedoligomer_T oligo) {
  Univcoord_T *positions;
  Positionsptr_T ptr0, end0;
  Storedoligomer_T part0;
#ifdef WORDS_BIGENDIAN
  int i;
  Positionsptr_T ptr;
  Univcoord_T bigendian, littleendian;
  char byte1, byte2, byte3;
#endif
#ifdef DEBUG0
  int j;
#endif


#if 0
  debug0(printf("%06X (%s)\n",oligo,shortoligo_nt(oligo,index1part)));
#endif
  part0 = oligo & poly_T;

  /* Ignore poly A and poly T on stage 1 */
  /* Was commented out */
  if (part0 == poly_A || part0 == poly_T) {
    *nentries = 0;
    return NULL;
  }

  if (this->compression_type == NO_COMPRESSION) {
#ifdef WORDS_BIGENDIAN
    /* Also holds for ALLOCATED_PRIVATE and ALLOCATED_SHARED */
    ptr0 = Bigendian_convert_uint(this->offsetsstrm[part0]);
    end0 = Bigendian_convert_uint(this->offsetsstrm[part0+1]);
#else
    ptr0 = this->offsetsstrm[part0];
    end0 = this->offsetsstrm[part0+1];
#endif

  } else if (this->compression_type == BITPACK64_COMPRESSION) {
#ifdef LARGE_GENOMES
    ptr0 = Bitpack64_read_two_huge(&end0,part0,this->offsetspages,this->offsetsmeta,this->offsetsstrm);
#else
    ptr0 = Bitpack64_read_two(&end0,part0,this->offsetsmeta,this->offsetsstrm);
#endif

  } else {
    abort();
  }


#ifdef ALLOW_DUPLICATES
  /* Not used */
  /* Skip backward over bad values, due to duplicates */
  if (this->positions_access == FILEIO) {
#ifdef HAVE_PTHREAD
    pthread_mutex_lock(&this->positions_read_mutex);
#endif
    positions_move_absolute(this->positions_fd,end0-1);
    while (end0 > ptr0 && positions_read_backward(this->positions_fd) == BADVAL) {
      end0--;
    }
#ifdef HAVE_PTHREAD
    pthread_mutex_unlock(&this->positions_read_mutex);
#endif
  } else {
    while (end0 > ptr0 && this->positions[end0-1] == BADVAL) {
      end0--;
    }
  }
#endif

  if ((*nentries = end0 - ptr0) == 0) {
    return NULL;
  } else {
    debug0(printf("Indexdb_read: offset pointers are %u and %u\n",ptr0,end0));
  
    positions = (Univcoord_T *) CALLOC(*nentries,sizeof(Univcoord_T));
    if (this->positions_access == FILEIO) {
#ifdef HAVE_PTHREAD
      pthread_mutex_lock(&this->positions_read_mutex);
#endif
#ifdef LARGE_GENOMES
      positions_move_absolute_1(this->positions_high_fd,ptr0);
      positions_move_absolute_4(this->positions_low_fd,ptr0);
      positions_read_multiple_large(this->positions_high_fd,this->positions_low_fd,positions,*nentries);
#else
      positions_move_absolute_4(this->positions_fd,ptr0);
      positions_read_multiple(this->positions_fd,positions,*nentries);
#endif
#ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&this->positions_read_mutex);
#endif
    } else if (this->positions_access == ALLOCATED_PRIVATE || this->positions_access == ALLOCATED_SHARED) {
#ifdef LARGE_GENOMES
      positions_copy_multiple_large(positions,&(this->positions_high[ptr0]),&(this->positions_low[ptr0]),*nentries);
#else
      memcpy(positions,&(this->positions[ptr0]),(*nentries)*sizeof(Univcoord_T));
#endif

    } else {
#ifdef WORDS_BIGENDIAN
#ifdef LARGE_GENOMES
      for (ptr = ptr0, i = 0; ptr < end0; ptr++, i++) {
	bigendian = (Univcoord_T) this->positions_high[ptr];
	bigendian <<= 8;

	littleendian = this->positions_low[ptr];
	bigendian |= littleendian & 0xff; /* 0 */
	bigendian <<= 8;
	bigendian |= ((littleendian >>= 8) & 0xff); /* 1 */
	bigendian <<= 8;
	bigendian |= ((littleendian >>= 8) & 0xff); /* 2 */
	bigendian <<= 8;
	bigendian |= ((littleendian >>= 8) & 0xff); /* 3 */
	positions[i] = bigendian;
      }
#else
      for (ptr = ptr0, i = 0; ptr < end0; ptr++, i++) {
	littleendian = this->positions[ptr];
	bigendian = littleendian & 0xff; /* 0 */
	bigendian <<= 8;
	bigendian |= ((littleendian >>= 8) & 0xff); /* 1 */
	bigendian <<= 8;
	bigendian |= ((littleendian >>= 8) & 0xff); /* 2 */
	bigendian <<= 8;
	bigendian |= ((littleendian >>= 8) & 0xff); /* 3 */
	positions[i] = bigendian;
      }
#endif

#else  /* littleendian */

#ifdef LARGE_GENOMES
      positions_copy_multiple_large(positions,&(this->positions_high[ptr0]),&(this->positions_low[ptr0]),*nentries);
#else
      memcpy(positions,&(this->positions[ptr0]),(*nentries)*sizeof(Univcoord_T));
#endif
#endif
    }

    debug0(
	   printf("%d entries:",*nentries);
	   for (j = 0; j < *nentries; j++) {
	     printf(" %u",positions[j]);
	   }
	   printf("\n");
	   );

    return positions;
  }
}


/* GSNAP version.  Expects calling procedure to handle bigendian conversion. */
UINT4 *
Indexdb_read_inplace (int *nentries,
#ifdef LARGE_GENOMES
		      unsigned char **positions_high,
#endif
		      T this, Storedoligomer_T oligo) {
  Positionsptr_T ptr0, end0;
  Storedoligomer_T part0;
#ifdef DEBUG0
  Positionsptr_T ptr;
#endif

  debug0(printf("%08X (%s)\n",oligo,shortoligo_nt(oligo,index1part)));
  part0 = oligo & poly_T;	/* Probably unnecessary, since stage1 procedure already masks oligo */

#if 0
  /* Needed to avoid overflow on 15-mers */
  if (part0 == poly_A || part0 == poly_T) {
    *nentries = 0;
    return NULL;
  }
#endif

  if (this->compression_type == NO_COMPRESSION) {
#ifdef WORDS_BIGENDIAN
    /* Also holds for ALLOCATED_PRIVATE and ALLOCATED_SHARED */
    ptr0 = Bigendian_convert_uint(this->offsetsstrm[part0]);
    end0 = Bigendian_convert_uint(this->offsetsstrm[part0+1]);
#else
    ptr0 = this->offsetsstrm[part0];
    end0 = this->offsetsstrm[part0+1];
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

#ifdef LARGE_GENOMES
  debug0(printf("Indexdb_read_inplace: offset pointers are %llu and %llu\n",ptr0,end0));
#else
  debug0(printf("Indexdb_read_inplace: offset pointers are %u and %u\n",ptr0,end0));
#endif

  *nentries = end0 - ptr0;

  if (*nentries == 0) {
    return NULL;
  } else if (this->positions_access == FILEIO) {
    abort();
  } else {
#ifdef LARGE_GENOMES
    debug0(
	   printf("%d entries:",*nentries);
	   for (ptr = ptr0; ptr < end0; ptr++) {
	     printf(" %d %u",this->positions_high[ptr],this->positions_low[ptr]);
	   }
	   printf("\n");
	   );

    *positions_high = &(this->positions_high[ptr0]);
    return &(this->positions_low[ptr0]);
#else
    debug0(
	   printf("%d entries:",*nentries);
	   for (ptr = ptr0; ptr < end0; ptr++) {
	     printf(" %u",this->positions[ptr]);
	   }
	   printf("\n");
	   );

    return &(this->positions[ptr0]);
#endif
  }
}

#endif	/* ifdef PMAP */


/* Analogous to Indexdb_read, except this includes diagterm.  Always allocates memory. */
Univcoord_T *
Indexdb_read_with_diagterm (int *nentries, T this, Storedoligomer_T oligo, int diagterm) {
  Univcoord_T *positions;
  Positionsptr_T ptr0, end0, ptr;
  int i;

  if (this->compression_type == NO_COMPRESSION) {
#ifdef WORDS_BIGENDIAN
    /* Also holds for ALLOCATED_PRIVATE and ALLOCATED_SHARED */
    ptr0 = Bigendian_convert_uint(this->offsetsstrm[oligo]);
    end0 = Bigendian_convert_uint(this->offsetsstrm[oligo+1]);
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

  debug0(printf("read_zero_shift: oligo = %06X, offset pointers are %u and %u\n",oligo,ptr0,end0));

  if ((*nentries = end0 - ptr0) == 0) {
    return (Univcoord_T *) NULL;
  } else {
    positions = (Univcoord_T *) CALLOC(*nentries,sizeof(Univcoord_T));
    if (this->positions_access == FILEIO) {
#ifdef HAVE_PTHREAD
      pthread_mutex_lock(&this->positions_read_mutex);
#endif
#ifdef LARGE_GENOMES
      positions_move_absolute_1(this->positions_high_fd,ptr0);
      positions_move_absolute_4(this->positions_low_fd,ptr0);
      positions_read_multiple_large(this->positions_high_fd,this->positions_low_fd,positions,*nentries);
#else
      positions_move_absolute_4(this->positions_fd,ptr0);
      positions_read_multiple(this->positions_fd,positions,*nentries);
#endif
#ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&this->positions_read_mutex);
#endif

    } else if (this->positions_access == ALLOCATED_PRIVATE || this->positions_access == ALLOCATED_SHARED) {
#ifdef LARGE_GENOMES
      for (ptr = ptr0, i = 0; ptr < end0; ptr++) {
	positions[i++] = ((Univcoord_T) this->positions_high[ptr] << 32) + this->positions_low[ptr] + diagterm;
      }
#else
      for (ptr = ptr0, i = 0; ptr < end0; ptr++) {
	positions[i++] = this->positions[ptr] + diagterm;
      }
#endif

    } else {

#ifdef WORDS_BIGENDIAN
#ifdef LARGE_GENOMES
      for (ptr = ptr0, i = 0; ptr < end0; ptr++) {
	positions[i++] = ((Univcoord_T) this->positions_high[ptr] << 32) + Bigendian_convert_uint(this->positions_low[ptr]) + diagterm;
      }
#else
      for (ptr = ptr0, i = 0; ptr < end0; ptr++) {
	positions[i++] = Bigendian_convert_univcoord(this->positions[ptr]) + diagterm;
      }
#endif
#else
#ifdef LARGE_GENOMES
      for (ptr = ptr0, i = 0; ptr < end0; ptr++) {
	positions[i++] = ((Univcoord_T) this->positions_high[ptr] << 32) + this->positions_low[ptr] + diagterm;
      }
#else
      for (ptr = ptr0, i = 0; ptr < end0; ptr++) {
	positions[i++] = this->positions[ptr] + diagterm;
      }
#endif
#endif
    }
  }
      
  debug0(
	printf("%d entries:",*nentries);
	for (i = 0; i < *nentries; i++) {
	  printf(" %u",positions[i]);
	}
	printf("\n");
	);
  
  return positions;
}


/* Analogous to Indexdb_read, except this includes diagterm.  Always allocates memory. */
Univcoord_T *
Indexdb_read_with_diagterm_sizelimit (int *nentries, T this, Storedoligomer_T oligo, int diagterm,
				      int size_threshold) {
  Univcoord_T *positions;
  Positionsptr_T ptr0, end0, ptr;
  int i;

  if (this->compression_type == NO_COMPRESSION) {
#ifdef WORDS_BIGENDIAN
    /* Also holds for ALLOCATED_PRIVATE and ALLOCATED_SHARED */
    ptr0 = Bigendian_convert_uint(this->offsetsstrm[oligo]);
    end0 = Bigendian_convert_uint(this->offsetsstrm[oligo+1]);
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


  debug0(printf("read_zero_shift: oligo = %06X, offset pointers are %u and %u\n",oligo,ptr0,end0));

  if ((*nentries = end0 - ptr0) == 0) {
    return (Univcoord_T *) NULL;

  } else if (*nentries > size_threshold) {
    *nentries = 0;
    return (Univcoord_T *) NULL;

  } else {
    positions = (Univcoord_T *) CALLOC(*nentries,sizeof(Univcoord_T));
    if (this->positions_access == FILEIO) {
#ifdef HAVE_PTHREAD
      pthread_mutex_lock(&this->positions_read_mutex);
#endif
#ifdef LARGE_GENOMES
      positions_move_absolute_1(this->positions_high_fd,ptr0);
      positions_move_absolute_4(this->positions_low_fd,ptr0);
      positions_read_multiple_large(this->positions_high_fd,this->positions_low_fd,positions,*nentries);
#else
      positions_move_absolute_4(this->positions_fd,ptr0);
      positions_read_multiple(this->positions_fd,positions,*nentries);
#endif

#ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&this->positions_read_mutex);
#endif

    } else if (this->positions_access == ALLOCATED_PRIVATE || this->positions_access == ALLOCATED_SHARED) {
#ifdef LARGE_GENOMES
      for (ptr = ptr0, i = 0; ptr < end0; ptr++) {
	positions[i++] = ((Univcoord_T) this->positions_high[ptr] << 32) + this->positions_low[ptr] + diagterm;
      }
#else
      for (ptr = ptr0, i = 0; ptr < end0; ptr++) {
	positions[i++] = this->positions[ptr] + diagterm;
      }
#endif

    } else {

#ifdef WORDS_BIGENDIAN
#ifdef LARGE_GENOMES
      for (ptr = ptr0, i = 0; ptr < end0; ptr++) {
	positions[i++] = ((Univcoord_T) this->positions_high[ptr] << 32) + Bigendian_convert_uint(this->positions_low[ptr]) + diagterm;
      }
#else
      for (ptr = ptr0, i = 0; ptr < end0; ptr++) {
	positions[i++] = Bigendian_convert_univcoord(this->positions[ptr]) + diagterm;
      }
#endif
#else
#ifdef LARGE_GENOMES
      for (ptr = ptr0, i = 0; ptr < end0; ptr++) {
	positions[i++] = ((Univcoord_T) this->positions_high[ptr] << 32) + this->positions_low[ptr] + diagterm;
      }
#else
      for (ptr = ptr0, i = 0; ptr < end0; ptr++) {
	positions[i++] = this->positions[ptr] + diagterm;
      }
#endif
#endif
    }
  }
      
  debug0(
	printf("%d entries:",*nentries);
	for (i = 0; i < *nentries; i++) {
	  printf(" %u",positions[i]);
	}
	printf("\n");
	);
  
  return positions;
}


/************************************************************************
 *   Create procedure -- for user-provided genomic segment
 ************************************************************************/

#if defined(UTILITYP) || defined(LARGE_GENOMES)
#else
T
Indexdb_new_segment (char *genomicseg,
#ifdef PMAP
		     int alphabet_size, Width_T index1part_aa, bool watsonp,
#else
		     Width_T index1part,
#endif
		     Width_T index1interval) {
  T new = (T) MALLOC(sizeof(*new));
  char *uppercaseCode;
  Positionsptr_T *work_offsets;	/* Working set for use in calculating positions */
  int totalcounts = 0;
  int c;
  Oligospace_T oligospace, oligoi;
  char *p;
  Univcoord_T position = 0UL;
#ifdef PMAP
  int frame = -1, between_counter[3], in_counter[3];
  Storedoligomer_T high = 0U, low = 0U, carry;
  Storedoligomer_T aaindex;
  int index1part_nt = 3*index1part_aa;
#else
  int between_counter = 0, in_counter = 0;
  Storedoligomer_T oligo = 0U, masked, mask;
#endif

  uppercaseCode = UPPERCASE_U2T;

#ifdef PMAP
  oligospace = power(alphabet_size,index1part_aa);
  between_counter[0] = between_counter[1] = between_counter[2] = 0;
  in_counter[0] = in_counter[1] = in_counter[2] = 0;
  new->index1part = index1part_aa;
#else
  mask = ~(~0UL << 2*index1part);
  oligospace = power(4,index1part);
  new->index1part = index1part;
#endif
  new->index1interval = 1;

  new->compression_type = NO_COMPRESSION;

  new->offsetsmeta = (UINT4 *) CALLOC(oligospace+1,sizeof(UINT4));
  for (oligoi = 0; oligoi <= oligospace; oligoi++) {
    new->offsetsmeta[oligoi] = oligoi;
  }
  new->offsetsmeta_access = ALLOCATED_PRIVATE;

  new->offsetsstrm = (UINT4 *) CALLOC(oligospace+1,sizeof(UINT4));
  new->offsetsstrm_access = ALLOCATED_PRIVATE;

  p = genomicseg;
  while ((c = *(p++)) != '\0') {
#ifdef PMAP
    if (++frame == 3) {
      frame = 0;
    }
    between_counter[frame] += 1;
    in_counter[frame] += 1;
#else
    between_counter++;
    in_counter++;
#endif

#ifdef PMAP
    carry = (low >> 30);
    switch (uppercaseCode[c]) {
    case 'A': low = (low << 2); break;
    case 'C': low = (low << 2) | 1U; break;
    case 'G': low = (low << 2) | 2U; break;
    case 'T': low = (low << 2) | 3U; break;
    default:
      high = low = carry = 0U;
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
      break;
    }
    high = (high << 2) | carry; 
#else
    switch (uppercaseCode[c]) {
    case 'A': oligo = (oligo << 2); break;
    case 'C': oligo = (oligo << 2) | 1U; break;
    case 'G': oligo = (oligo << 2) | 2U; break;
    case 'T': oligo = (oligo << 2) | 3U; break;
    default: oligo = 0U; in_counter = 0; break;
    }
#endif

    /*
    debug(printf("char=%c bc=%d ic=%d oligo=%016lX\n",
		 c,between_counter,in_counter,oligo));
    */

#ifdef PMAP
    if (in_counter[frame] > 0) {
      if (watsonp == true) {
	if (Alphabet_get_codon_fwd(low) == AA_STOP) {
	  in_counter[frame] = 0; 
	}
      } else {
	if (Alphabet_get_codon_rev(low) == AA_STOP) {
	  in_counter[frame] = 0; 
	}
      }
    }
    if (in_counter[frame] == index1part_aa + 1) {
      if (between_counter[frame] >= index1interval) {
	aaindex = Alphabet_get_aa_index(high,low,watsonp,index1part_nt);
#ifdef OLIGOSPACE_NOT_LONG
	oligoi = (Oligospace_T) aaindex + 1U;
#else
	oligoi = (Oligospace_T) aaindex + 1UL;
#endif
	new->offsetsstrm[oligoi] += 1;
	between_counter[frame] = 0;
      }
      in_counter[frame] -= 1;
    }
#else
    if (in_counter == index1part) {
      if (
#ifdef NONMODULAR
	  between_counter >= index1interval
#else
	  /* Actually, modular condition not needed for user-supplied genomic segment */
	  (position-index1part+1U) % index1interval == 0
#endif
	  ) {
	masked = oligo & mask;
#ifdef OLIGOSPACE_NOT_LONG
	oligoi = (Oligospace_T) masked + 1U;
#else
	oligoi = (Oligospace_T) masked + 1UL;
#endif
	new->offsetsstrm[oligoi] += 1;
	debug(printf("Found oligo %06X.  Incremented offsets for %llu to be %u\n",
		     masked,(unsigned long long) oligoi,new->offsetsstrm[oligoi]));
	between_counter = 0;
      }
      in_counter--;
    }
#endif

    position++;
  }

#ifdef ADDSENTINEL
  for (oligoi = 1; oligoi <= oligospace; oligoi++) {
    new->offsetsstrm[oligoi] = new->offsetsstrm[oligoi] + new->offsetsstrm[oligoi-1] + 1;
    debug(printf("Offset for %06X: %u\n",oligoi,new->offsetsstrm[oligoi]));
  }
#else
  for (oligoi = 1; oligoi <= oligospace; oligoi++) {
    new->offsetsstrm[oligoi] = new->offsetsstrm[oligoi] + new->offsetsstrm[oligoi-1];
    debug(printf("Offset for %06X: %u\n",oligoi,new->offsetsstrm[oligoi]));
  }
#endif


  /* Create positions */
  position = 0U;
#ifdef PMAP
  frame = -1;
  between_counter[0] = between_counter[1] = between_counter[2] = 0;
  in_counter[0] = in_counter[1] = in_counter[2] = 0;
  high = low = 0U;
#else
  between_counter = in_counter = 0;
  oligo = 0UL;
#endif

  work_offsets = (Positionsptr_T *) CALLOC(oligospace+1,sizeof(Positionsptr_T));
  for (oligoi = 0; oligoi <= oligospace; oligoi++) {
    work_offsets[oligoi] = new->offsetsstrm[oligoi];
  }

  totalcounts = new->offsetsstrm[oligospace];
  if (totalcounts == 0) {
#ifdef PMAP
    fprintf(stderr,"Error: user-provided genomic segment has no valid oligomers of size %d\n",index1part_nt);
#else
    fprintf(stderr,"Error: user-provided genomic segment has no valid oligomers of size %d\n",index1part);
#endif
    exit(9);
  }
  new->positions = (Univcoord_T *) CALLOC(totalcounts,sizeof(Univcoord_T));
  new->positions_access = ALLOCATED_PRIVATE;

  p = genomicseg;
  while ((c = *(p++)) != '\0') {
#ifdef PMAP
    if (++frame == 3) {
      frame = 0;
    }
    between_counter[frame] += 1;
    in_counter[frame] += 1;
#else
    between_counter++;
    in_counter++;
#endif

#ifdef PMAP
    carry = (low >> 30);
    switch (uppercaseCode[c]) {
    case 'A': low = (low << 2); break;
    case 'C': low = (low << 2) | 1; break;
    case 'G': low = (low << 2) | 2; break;
    case 'T': low = (low << 2) | 3; break;
    default:
      high = low = carry = 0U;
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
      break;
    }
    high = (high << 2) | carry; 
#else
    switch (uppercaseCode[c]) {
    case 'A': oligo = (oligo << 2); break;
    case 'C': oligo = (oligo << 2) | 1; break;
    case 'G': oligo = (oligo << 2) | 2; break;
    case 'T': oligo = (oligo << 2) | 3; break;
    default: oligo = 0U; in_counter = 0; break;
    }
#endif

    /*
    debug(printf("char=%c bc=%d ic=%d oligo=%06X\n",
		 c,between_counter,in_counter,oligo));
    */
    
#ifdef PMAP
    if (in_counter[frame] > 0) {
      if (watsonp == true) {
	if (Alphabet_get_codon_fwd(low) == AA_STOP) {
	  in_counter[frame] = 0; 
	}
      } else {
	if (Alphabet_get_codon_rev(low) == AA_STOP) {
	  in_counter[frame] = 0; 
	}
      }
    }

    if (in_counter[frame] == index1part_aa + 1) {
      if (between_counter[frame] >= index1interval) {
	aaindex = Alphabet_get_aa_index(high,low,watsonp,index1part_nt);
	if (watsonp == true) {
	  new->positions[work_offsets[aaindex]++] = position-index1part_nt+1;
	} else {
	  new->positions[work_offsets[aaindex]++] = position;
	}
	between_counter[frame] = 0;
      }
      in_counter[frame] -= 1;
    }
#else
    if (in_counter == index1part) {
      if (
#ifdef NONMODULAR
	  between_counter >= index1interval
#else
	  /* Actually, modular condition not needed for user-supplied genomic segment */
	  (position-index1part+1) % index1interval == 0
#endif
	  ) {
	masked = oligo & mask;
	new->positions[work_offsets[masked]++] = position-index1part+1;
	between_counter = 0;
      }
      in_counter--;
    }
#endif

    position++;
  }

#ifdef ADDSENTINEL
  for (oligoi = 0; oligoi < oligospace; oligoi++) {
    new->positions[work_offsets[oligoi]] = (Univcoord_T) -1;
  }
#endif

  FREE(work_offsets);

  return new;
}
#endif


int
Storedoligomer_compare (const void *a, const void *b) {
  Storedoligomer_T x = * (Storedoligomer_T *) a;
  Storedoligomer_T y = * (Storedoligomer_T *) b;

  if (x < y) {
    return -1;
  } else if (y < x) {
    return 1;
  } else {
    return 0;
  }
}

