static char rcsid[] = "$Id: gmapindex.c 167265 2015-06-11 00:04:50Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif
#include <ctype.h>
#include <string.h>
#include <strings.h>		/* For rindex */
#include <sys/mman.h>		/* For munmap */
#include "types.h"
#include "bool.h"
#include "assert.h"
#include "mem.h"
#include "fopen.h"

#include "table.h"
#ifdef HAVE_64_BIT
#include "tableuint8.h"
typedef Tableuint8_T Table_chrpos_T;
#else
typedef Tableuint_T Table_chrpos_T;
#endif
#include "tableuint.h"
#include "compress.h"
#include "chrom.h"
#include "segmentpos.h"
#include "univinterval.h"
#include "iit-write-univ.h"
#include "iit-read-univ.h"
#include "genome.h"
#include "genome128_hr.h"
#include "genome-write.h"
#include "indexdb-write.h"
#include "compress-write.h"
#include "intlist.h"
#include "indexdbdef.h"		/* For compression types */
#include "indexdb.h"		/* For Filenames_T */
#include "sarray-write.h"
#include "bytecoding.h"
#include "bitpack64-write.h"

#ifdef HAVE_SSE2
#include <emmintrin.h>
#endif
#ifdef HAVE_SSE4_1
#include <smmintrin.h>
#endif
#ifdef HAVE_POPCNT
#include <immintrin.h>
#endif
#if defined(HAVE_MM_POPCNT)
#include <nmmintrin.h>
#endif
#if defined(HAVE_LZCNT) || defined(HAVE_BMI1)
#include <immintrin.h>
#endif


#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
#else
#include "littleendian.h"
#endif

#define BUFFERSIZE 8192

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


/* Program variables */
typedef enum {NONE, AUXFILES, GENOME, UNSHUFFLE, COUNT, OFFSETS, POSITIONS, SUFFIX_ARRAY, LCP_CHILD,
	      COMPRESSED_SUFFIX_ARRAY, ARRAY_UNCOMPRESS, CHILD_UNCOMPRESS} Action_T;
static Action_T action = NONE;
static char *sourcedir = ".";
static char *destdir = ".";
static char *fileroot = NULL;
static int compression_types = BITPACK64_COMPRESSION;
static int compression_type;
static int index1part = 15;
static int index1interval = 3;	/* Interval for storing 12-mers */
static bool genome_lc_p = false;
static bool rawp = false;
static bool writefilep = false;
/* static bool sortchrp = true;	? Sorting now based on order in .coords file */
static int wraplength = 0;
static bool mask_lowercase_p = false;
static int nmessages = 50;

static char *mitochondrial_string = NULL;
static Sorttype_T divsort = CHROM_SORT;
static char *sortfilename = NULL;
static bool huge_offsets_p = false;


static void
check_compiler_assumptions () {
  unsigned int x = rand(), y = rand();
#ifdef HAVE_SSE2
  int z;
  __m128i a;
#ifdef HAVE_SSE4_1
  char negx, negy;
#endif
#endif


  fprintf(stderr,"Checking compiler assumptions for popcnt: ");
  fprintf(stderr,"%08X ",x);
#ifdef HAVE_LZCNT
  fprintf(stderr,"_lzcnt_u32=%d ",_lzcnt_u32(x));
#endif
#ifdef HAVE_BUILTIN_CLZ
  fprintf(stderr,"__builtin_clz=%d ",__builtin_clz(x));
#endif
#ifdef HAVE_BMI1
  fprintf(stderr,"_tzcnt_u32=%d ",_tzcnt_u32(x));
#endif
#ifdef HAVE_BUILTIN_CTZ
  fprintf(stderr,"__builtin_ctz=%d ",__builtin_ctz(x));
#endif

#ifdef HAVE_POPCNT
  fprintf(stderr,"_popcnt32=%d ",_popcnt32(x));
#endif
#if defined(HAVE_MM_POPCNT)
  fprintf(stderr,"_mm_popcnt_u32=%d ",_mm_popcnt_u32(x));
#endif
#if defined(HAVE_BUILTIN_POPCOUNT)
  fprintf(stderr,"__builtin_popcount=%d ",__builtin_popcount(x));
#endif

  fprintf(stderr,"\n");

#ifdef HAVE_SSE2
  fprintf(stderr,"Checking compiler assumptions for SSE2: ");
  fprintf(stderr,"%08X %08X",x,y);
  a = _mm_xor_si128(_mm_set1_epi32(x),_mm_set1_epi32(y));
  z = _mm_cvtsi128_si32(a);
  fprintf(stderr," xor=%08X\n",z);

#ifdef HAVE_SSE4_1
  if ((negx = (char) x) > 0) {
    negx = -negx;
  }
  if ((negy = (char) y) > 0) {
    negy = -negy;
  }

  fprintf(stderr,"Checking compiler assumptions for SSE4.1: ");
  fprintf(stderr,"%d %d",negx,negy);
  a = _mm_max_epi8(_mm_set1_epi8(negx),_mm_set1_epi8(negy));
  z = _mm_extract_epi8(a,0);
  fprintf(stderr," max=%d => ",z);
  if (negx > negy) {
    if (z == (int) negx) {
      fprintf(stderr,"compiler sign extends\n"); /* technically incorrect, but SIMD procedures behave properly */
    } else {
      fprintf(stderr,"compiler zero extends\n");
    }
  } else {
    if (z == (int) negy) {
      fprintf(stderr,"compiler sign extends\n"); /* technically incorrect, but SIMD procedures behave properly */
    } else {
      fprintf(stderr,"compiler zero extends\n");
    }
  }

#endif

#endif

  fprintf(stderr,"Finished checking compiler assumptions\n");

  return;
}


#if 0
/************************************************************************
 *   Reading strain from file
 ************************************************************************/

static char *
read_strain_from_strainfile (char *strainfile) {
  FILE *fp;
  char *refstrain = NULL, Buffer[1024], strain[1024], straintype[1024];

  if (strainfile != NULL) {
    fp = fopen(strainfile,"r");
    if (fp == NULL) {
      fprintf(stderr,"Cannot open strain file %s\n",strainfile);
    } else {
      while (fgets(Buffer,1024,fp) != NULL) {
	if (Buffer[0] == '#') {
	  /* Skip */
	} else if (sscanf(Buffer,"%s %s",strain,straintype) == 2) {
	  if (!strcmp(straintype,"reference") || !strcmp(straintype,"Reference") || 
	      !strcmp(straintype,"REFERENCE")) {
	    if (refstrain != NULL) {
	      fprintf(stderr,"More than one reference strain seen in %s\n",strainfile);
	      exit(9);
	    }
	    refstrain = (char *) CALLOC(strlen(strain)+1,sizeof(char));
	    strcpy(refstrain,strain);
	  }
	}
      }

      fclose(fp);
    }
  }

  if (refstrain != NULL) {
    return refstrain;
  } else {
    refstrain = (char *) CALLOC(strlen("reference")+1,sizeof(char));
    strcpy(refstrain,"reference");
    return refstrain;
  }
}

static char *
read_strain_from_coordsfile (char *coordsfile) {
  FILE *fp;
  char *refstrain = NULL, Buffer[1024], strain[1024], *ptr;

  if (coordsfile != NULL) {
    fp = fopen(coordsfile,"r");
    if (fp == NULL) {
      fprintf(stderr,"Cannot open coords file %s\n",coordsfile);
    } else {
      while (fgets(Buffer,1024,fp) != NULL) {
	if (Buffer[0] == '#') {
	  if ((ptr = strstr(Buffer,"Reference strain:")) != NULL) {
	    if (sscanf(ptr,"Reference strain: %s",strain) == 1) {
	      if (refstrain != NULL) {
		fprintf(stderr,"More than one reference strain seen in %s\n",coordsfile);
		exit(9);
	      }
	      refstrain = (char *) CALLOC(strlen(strain)+1,sizeof(char));
	      strcpy(refstrain,strain);
	    }
	  }
	}	    
      }

      fclose(fp);
    }
  }

  if (refstrain != NULL) {
    return refstrain;
  } else {
    refstrain = (char *) CALLOC(strlen("reference")+1,sizeof(char));
    strcpy(refstrain,"reference");
    return refstrain;
  }
}
#endif




/************************************************************************
 *   Creating aux file
 ************************************************************************/

/* accsegmentpos_table: char *accession -> Segmentpos_T segmentpos
   chrlength_table:     Chrom_T chrom -> Chrpos_T chrlength
*/

static void
chrlength_update (Table_chrpos_T chrlength_table, Chrom_T chrom, Univcoord_T segend) {
  Univcoord_T oldsegend;

#ifdef HAVE_64_BIT
  if ((oldsegend = (Univcoord_T) Tableuint8_get(chrlength_table,chrom)) == 0) {
    /* Initial entry for this chromosome */
    Tableuint8_put(chrlength_table,chrom,segend);

  } else if (segend > oldsegend) {
    /* Revise */
    Tableuint8_put(chrlength_table,chrom,segend);
  }
#else
  if ((oldsegend = (Univcoord_T) Tableuint_get(chrlength_table,chrom)) == 0) {
    /* Initial entry for this chromosome */
    Tableuint_put(chrlength_table,chrom,segend);

  } else if (segend > oldsegend) {
    /* Revise */
    Tableuint_put(chrlength_table,chrom,segend);
  }
#endif

  return;
}

static void
store_accession (Table_T accsegmentpos_table, Table_chrpos_T chrlength_table, Tableuint_T chrorder_table,
		 char *accession, char *chr_string, Chrpos_T chrpos1, 
		 Chrpos_T chrpos2, bool revcompp, Chrpos_T seglength, 
		 int contigtype, unsigned int universal_coord, bool circularp) {
  Chrom_T chrom;
  Segmentpos_T segmentpos;
  unsigned int order;

  if (chrorder_table != NULL) {
    order = Tableuint_get(chrorder_table,chr_string);
    chrom = Chrom_from_string(chr_string,mitochondrial_string,order,circularp);
  } else {
    chrom = Chrom_from_string(chr_string,mitochondrial_string,/*order*/universal_coord,circularp);
  }

  segmentpos = Segmentpos_new(chrom,chrpos1,chrpos2,revcompp,seglength,contigtype);
  Table_put(accsegmentpos_table,(void *) accession,(void *) segmentpos);

  /* Update chrlength */
  if (chrpos2 > chrpos1 + seglength) {
    chrlength_update(chrlength_table,chrom,chrpos2);
  } else {
    chrlength_update(chrlength_table,chrom,chrpos1+seglength);
  }

  return;
}


/* We assume that header has already been read.  We need to check each
   new line for a new header */
static Chrpos_T
count_sequence () {
  Chrpos_T seglength = 0U;
  int c;
  char Buffer[BUFFERSIZE], *p;
  bool newline = true;

  while (1) {
    /* Start of new line */
    if (newline == true) {
      if ((c = getc(stdin)) == EOF || c == '>') {
	return seglength;
      } else {
	seglength += 1U;
      }
    }

    if (fgets(Buffer,BUFFERSIZE,stdin) == NULL) {
      return seglength;
    } else {
      if ((p = rindex(Buffer,'\n')) != NULL) {
	*p = '\0';
	newline = true;
      } else {
	newline = false;
      }
      seglength += (Chrpos_T) strlen(Buffer);
    }
  }
}

static void
skip_sequence (Chrpos_T seglength) {
  int c;
  char Buffer[BUFFERSIZE];

  while (seglength > BUFFERSIZE) {
    if (fread(Buffer,sizeof(char),BUFFERSIZE,stdin) < BUFFERSIZE) {
      fprintf(stderr,"End of file reached.  Expecting %u more characters\n",seglength);
      exit(9);
    }
    seglength -= BUFFERSIZE;
  }

  if (seglength > 0U) {
    if (fread(Buffer,sizeof(char),seglength,stdin) < seglength) {
      fprintf(stderr,"End of file reached.  Expecting %u more characters\n",seglength);
      exit(9);
    }
  }

  if ((c = getchar()) != EOF && c != '\n') {
    fprintf(stderr,"Expecting linefeed at end of sequence.  Saw %d (%c) instead\n",c,c);
    exit(9);
  }

  if ((c = getchar()) != EOF && c != '>') {
    fprintf(stderr,"Expecting new FASTA line.  Saw %d (%c) instead\n",c,c);
    exit(9);
  }

  return;
}

static bool
process_sequence_aux (Chrpos_T *seglength, Table_T accsegmentpos_table, Table_chrpos_T chrlength_table,
		      Tableuint_T chrorder_table, char *fileroot, int ncontigs) {
  char Buffer[BUFFERSIZE], accession_p[BUFFERSIZE], *accession, 
    chrpos_string[BUFFERSIZE], *chr_string, *coords, *ptr, *p;
  Chrpos_T chrpos1, chrpos2, lower, upper;
  Univcoord_T universal_coord = 0U;
  bool revcompp, circularp;
  int nitems;

  /* Store sequence info */
  if (fgets(Buffer,BUFFERSIZE,stdin) == NULL) {
    return false;
  }

  nitems = sscanf(Buffer,"%s %s %llu",accession_p,chrpos_string,&universal_coord);
  if (nitems < 2) {
    fprintf(stderr,"Can't parse line %s\n",Buffer);
    exit(1);
  } else {
    if (ncontigs < nmessages) {
      fprintf(stderr,"Logging contig %s at %s in genome %s\n",accession_p,chrpos_string,fileroot);
    } else if (ncontigs == nmessages) {
      fprintf(stderr,"More than %d contigs.  Will stop printing messages\n",nmessages);
    }

    if (!index(chrpos_string,':')) {
      fprintf(stderr,"Can't parse chromosomal coordinates %s\n",chrpos_string);
      exit(1);
    } else {
      chr_string = strtok(chrpos_string,":");
      coords = strtok(NULL,":");
      if (sscanf(coords,"%u..%u",&chrpos1,&chrpos2) == 2) {
	/* 1:3..5, one-based, inclusive => (2,5), zero-based, boundaries */
	if (chrpos1 <= chrpos2) {
	  chrpos1--;
	  revcompp = false;
	  lower = chrpos1;
	  upper = chrpos2;
	} else {
	  chrpos2--;
	  revcompp = true;
	  lower = chrpos2;
	  upper = chrpos1;
	}
      } else if (sscanf(coords,"%u",&chrpos1) == 1) {
	/* 1:3, one-based, inclusive => (3,3), zero-based, boundaries */
	revcompp = false;
	lower = upper = chrpos1;
      } else {
	fprintf(stderr,"Can't parse chromosomal coordinates %s\n",coords);
	exit(1);
      }
    }

#if 0
    /* No longer supporting strains/types */
    p = Buffer;
    while (*p != '\0' && !isspace((int) *p)) { p++; } /* Skip to first space */
    while (*p != '\0' && isspace((int) *p)) { p++; } /* Skip past first space */
    while (*p != '\0' && !isspace((int) *p)) { p++; } /* Skip to second space */
    while (*p != '\0' && isspace((int) *p)) { p++; } /* Skip past second space */

    if (*p == '\0') {
      contigtype = 0;		/* Empty type string */
    } else {
      if ((ptr = rindex(p,'\n')) != NULL) {
	while (isspace((int) *ptr)) { ptr--; } /* Erase empty space */
	ptr++;
	*ptr = '\0';
      }
      if ((contigtype = Tableuint_get(contigtype_table,(void *) p)) == 0) {
	debug(printf("Storing type %s.\n",p));
	/* Store types as 1-based */
	contigtype = Tableuint_length(contigtype_table) + 1;
	typestring = (char *) CALLOC(strlen(p)+1,sizeof(char));
	strcpy(typestring,p);
	Tableuint_put(contigtype_table,(void *) typestring,contigtype);
	*contigtypelist = List_push(*contigtypelist,typestring);
      }
    }
#endif

    /* Look for circular.  Code modeled after parsing for strain above. */
    p = Buffer;
    while (*p != '\0' && !isspace((int) *p)) { p++; } /* Skip to first space */
    while (*p != '\0' && isspace((int) *p)) { p++; } /* Skip past first space */
    while (*p != '\0' && !isspace((int) *p)) { p++; } /* Skip to second space */
    while (*p != '\0' && isspace((int) *p)) { p++; } /* Skip past second space */
    while (*p != '\0' && !isspace((int) *p)) { p++; } /* Skip to third space */
    while (*p != '\0' && isspace((int) *p)) { p++; } /* Skip past third space */

    circularp = false;
    if (*p == '\0') {
      /* circularp = false; */
    } else {
      if ((ptr = rindex(p,'\n')) != NULL) {
	while (isspace((int) *ptr)) { ptr--; } /* Erase empty space */
	ptr++;
	*ptr = '\0';
      }
      if (!strcmp(p,"circular")) {
	fprintf(stderr,"%s is a circular chromosome\n",chr_string);
	circularp = true;
      } else {
	fprintf(stderr,"%s has an unrecognized tag %s\n",chr_string,p);
      }
    }

    /* The '>' character was already stripped off by the last call to count_sequence() */
    accession = (char *) CALLOC(strlen(accession_p)+1,sizeof(char));
    strcpy(accession,accession_p);
  }

  if (rawp == true) {
    *seglength = upper - lower;
    fprintf(stderr,"Skipping %u characters\n",*seglength);
    skip_sequence(*seglength);
  } else {
    *seglength = count_sequence();
    if (*seglength != upper - lower) {
      fprintf(stderr,"%s has expected sequence length %u-%u=%u but actual length %u\n",
	      accession,upper,lower,upper-lower,*seglength);
    }
  }

  if (nitems < 3) {
    universal_coord = 0U;
  }

  store_accession(accsegmentpos_table,chrlength_table,chrorder_table,
		  accession,chr_string,lower,upper,revcompp,
		  *seglength,/*contigtype*/0,universal_coord,circularp);

  return true;
}


/************************************************************************
 *   Creating genome and related files
 ************************************************************************/

/* Modifies chrlength_table to store offsets, rather than chrlengths */
static void
write_chromosome_file (char *genomesubdir, char *fileroot, Table_chrpos_T chrlength_table,
		       bool coord_values_8p) {
  FILE *textfp, *chrsubsetfp;
  char *divstring, *textfile, *chrsubsetfile, *iitfile, *chr_string, emptystring[1];
  int n, i;
  int typeint;
  Chrom_T *chroms;
  Univcoord_T chroffset = 0;
  Chrpos_T chrlength;
  List_T divlist = NULL;
  List_T intervallist = NULL, chrtypelist = NULL, labellist = NULL, annotlist = NULL, p;
  Table_T intervaltable, labeltable, annottable;
  Univinterval_T interval;

  emptystring[0] = '\0';

  if (divsort == NO_SORT) {
    fprintf(stderr,"divsort == NO_SORT\n");
#ifdef HAVE_64_BIT
    chroms = (Chrom_T *) Tableuint8_keys_by_timeindex(chrlength_table,0U);
    n = Tableuint8_length(chrlength_table);
#else
    chroms = (Chrom_T *) Tableuint_keys_by_timeindex(chrlength_table,0U);
    n = Tableuint_length(chrlength_table);
#endif

  } else {
    /* Get chromosomes in order */
#ifdef HAVE_64_BIT
    chroms = (Chrom_T *) Tableuint8_keys(chrlength_table,0U);
    n = Tableuint8_length(chrlength_table);
#else
    chroms = (Chrom_T *) Tableuint_keys(chrlength_table,0U);
    n = Tableuint_length(chrlength_table);
#endif
    switch (divsort) {
    case ALPHA_SORT: qsort(chroms,n,sizeof(Chrom_T),Chrom_compare_alpha); break;
    case NUMERIC_ALPHA_SORT: qsort(chroms,n,sizeof(Chrom_T),Chrom_compare_numeric_alpha); break;
    case CHROM_SORT: qsort(chroms,n,sizeof(Chrom_T),Chrom_compare_chrom); break;
    case FILENAME_SORT: qsort(chroms,n,sizeof(Chrom_T),Chrom_compare_order); break;
    default: abort();
    }
  }
  fprintf(stderr,"Have a total of %d chromosomes\n",n);
    

  /* Write chromosome text file and chrsubset file */
  textfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			     strlen(fileroot)+strlen(".chromosome")+1,sizeof(char));
  sprintf(textfile,"%s/%s.chromosome",genomesubdir,fileroot);
  /* Use binary, not text, so files are Unix-compatible */
  if ((textfp = FOPEN_WRITE_BINARY(textfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",textfile);
    exit(9);
  }
  FREE(textfile);

  chrsubsetfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
				  strlen(fileroot)+strlen(".chrsubset")+1,sizeof(char));
  sprintf(chrsubsetfile,"%s/%s.chrsubset",genomesubdir,fileroot);
  /* Use binary, not text, so files are Unix-compatible */
  if ((chrsubsetfp = FOPEN_WRITE_BINARY(chrsubsetfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",chrsubsetfile);
    exit(9);
  }
  FREE(chrsubsetfile);
  fprintf(chrsubsetfp,">all\n");
  fprintf(chrsubsetfp,"\n");

  chrtypelist = List_push(chrtypelist,"circular");
  chrtypelist = List_push(chrtypelist,"");
  for (i = 0; i < n; i++) {
#ifdef HAVE_64_BIT
    chrlength = (Chrpos_T) Tableuint8_get(chrlength_table,chroms[i]);
#else
    chrlength = (Chrpos_T) Tableuint_get(chrlength_table,chroms[i]);
#endif
    assert(chroffset <= chroffset + (Univcoord_T) chrlength - 1);
    chr_string = Chrom_string(chroms[i]);
    if (i < nmessages) {
      fprintf(stderr,"Chromosome %s has universal coordinates %llu..%llu\n",
	      chr_string,(unsigned long long) chroffset+1,
	      (unsigned long long) chroffset + 1 + (Univcoord_T) chrlength - 1);
    } else if (i == nmessages) {
      fprintf(stderr,"More than %d contigs.  Will stop printing messages\n",nmessages);
    }
      
    if (n <= 100) {
      fprintf(chrsubsetfp,">%s\n",chr_string);
      fprintf(chrsubsetfp,"+%s\n",chr_string);
    }

    fprintf(textfp,"%s\t%llu..%llu\t%u",
	    chr_string,(unsigned long long) chroffset+1,
	    (unsigned long long) chroffset + (Univcoord_T) chrlength,chrlength);
    if (Chrom_circularp(chroms[i]) == true) {
      fprintf(textfp,"\tcircular");
      typeint = 1;
    } else {
      typeint = 0;
    }
    fprintf(textfp,"\n");

    intervallist = List_push(intervallist,(void *) Univinterval_new(chroffset,chroffset + (Univcoord_T) chrlength - 1U,typeint));
    labellist = List_push(labellist,(void *) chr_string);
    annotlist = List_push(annotlist,(void *) emptystring); /* No annotations */

    /* Now chrlength_table holds chroffsets, not chrlengths */
#ifdef HAVE_64_BIT
    Tableuint8_put(chrlength_table,chroms[i],chroffset);
#else
    Tableuint_put(chrlength_table,chroms[i],chroffset);
#endif
    if (Chrom_circularp(chroms[i]) == true) {
      chroffset += (Univcoord_T) chrlength;
      chroffset += (Univcoord_T) chrlength;
    } else {
      chroffset += (Univcoord_T) chrlength;
    }
  }
  FREE(chroms);
  intervallist = List_reverse(intervallist);
  labellist = List_reverse(labellist);

  fclose(chrsubsetfp);
  fclose(textfp);

  /* Write chromosome IIT file */
  divstring = (char *) CALLOC(1,sizeof(char));
  divstring[0] = '\0';
  divlist = List_push(NULL,(void *) divstring);

  intervaltable = Table_new(65522,Table_string_compare,Table_string_hash);
  labeltable = Table_new(65522,Table_string_compare,Table_string_hash);
  annottable = Table_new(65522,Table_string_compare,Table_string_hash);

  Table_put(intervaltable,(void *) divstring,intervallist);
  Table_put(labeltable,(void *) divstring,labellist);
  Table_put(annottable,(void *) divstring,annotlist);

  iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			    strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
  sprintf(iitfile,"%s/%s.chromosome.iit",genomesubdir,fileroot);
  IIT_write_univ(iitfile,divlist,chrtypelist,
		 intervaltable,labeltable,annottable,
		 coord_values_8p,/*label_pointers_8p*/false,/*annot_pointers_8p*/false);
  FREE(iitfile);

  List_free(&divlist);
  FREE(divstring);

  Table_free(&annottable);
  Table_free(&labeltable);
  Table_free(&intervaltable);

  List_free(&annotlist);

  /* Do not free strings in labellist, since they are not allocated */
  List_free(&labellist);

  /* chrtypelist has no dynamically allocated strings */
  List_free(&chrtypelist);

  for (p = intervallist; p != NULL; p = List_next(p)) {
    interval = (Univinterval_T) List_head(p);
    Univinterval_free(&interval);
  }
  List_free(&intervallist);

  return;
}

static Table_T current_accsegmentpos_table;

/* x is really a char ** */
static int
bysegmentpos_compare_alpha (const void *x, const void *y) {
  char *acc1 = * (char **) x;
  char *acc2 = * (char **) y;
  Segmentpos_T a = (Segmentpos_T) Table_get(current_accsegmentpos_table,(void *) acc1);
  Segmentpos_T b = (Segmentpos_T) Table_get(current_accsegmentpos_table,(void *) acc2);

  return Segmentpos_compare_alpha(&a,&b);
}

static int
bysegmentpos_compare (const void *x, const void *y) {
  char *acc1 = * (char **) x;
  char *acc2 = * (char **) y;
  Segmentpos_T a = (Segmentpos_T) Table_get(current_accsegmentpos_table,(void *) acc1);
  Segmentpos_T b = (Segmentpos_T) Table_get(current_accsegmentpos_table,(void *) acc2);

  if (divsort == ALPHA_SORT) {
    return Segmentpos_compare_alpha(&a,&b);
  } else if (divsort == NUMERIC_ALPHA_SORT) {
    return Segmentpos_compare_numeric_alpha(&a,&b);
  } else if (divsort == CHROM_SORT) {
    return Segmentpos_compare_chrom(&a,&b);
  } else if (divsort == FILENAME_SORT) {
    return Segmentpos_compare_order(&a,&b);
  } else {
    abort();
  }
}

static void
write_contig_file (char *genomesubdir, char *fileroot, Table_T accsegmentpos_table, 
		   Table_chrpos_T chrlength_table, List_T contigtypelist, bool coord_values_8p) {
  FILE *textfp;
  char *textfile, *iitfile, *annot;
  int naccessions, i;
  char **accessions, *divstring, seglength[12]; /* 2^32 = 4*10^9 */
  Segmentpos_T segmentpos;
  Chrom_T chrom;
  Univcoord_T chroffset, universalpos1, universalpos2;
  List_T divlist = NULL, intervallist = NULL, labellist = NULL, annotlist = NULL, p;
  Table_T intervaltable, labeltable, annottable;
  Univinterval_T interval;
#if 0
  void **keys;
  int *values, ntypes;
#endif
  
  if (divsort == NO_SORT) {
    accessions = (char **) Table_keys_by_timeindex(accsegmentpos_table,NULL);
    naccessions = Table_length(accsegmentpos_table);
  } else {
    /* Get accessions in order */
    accessions = (char **) Table_keys(accsegmentpos_table,NULL);
    naccessions = Table_length(accsegmentpos_table);
    current_accsegmentpos_table = accsegmentpos_table;
    qsort(accessions,naccessions,sizeof(char *),bysegmentpos_compare);
  }

#if 0
  /* Get types in order */
  keys = Tableuint_keys(contigtype_table,NULL);
  values = Tableuint_values(contigtype_table,0);
  ntypes = Tableuint_length(contigtype_table);
  contigtypes = (char **) CALLOC(ntypes+1,sizeof(char *)); /* Add 1 for type 0 */
  contigtypes[0] = "";
  for (j = 0; j < ntypes; j++) {
    contigtypes[values[j]] = keys[j];
  }
  FREE(values);
  FREE(keys);
#endif

  /* Write contig text file */
  textfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			     strlen(fileroot)+strlen(".contig")+1,sizeof(char));
  sprintf(textfile,"%s/%s.contig",genomesubdir,fileroot);
  /* Use binary, not text, so files are Unix-compatible */
  if ((textfp = FOPEN_WRITE_BINARY(textfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",textfile);
    exit(9);
  }
  FREE(textfile);

  for (i = 0; i < naccessions; i++) {
    segmentpos = (Segmentpos_T) Table_get(accsegmentpos_table,(void *) accessions[i]);
    chrom = Segmentpos_chrom(segmentpos);
#ifdef HAVE_64_BIT
    chroffset = (Univcoord_T) Tableuint8_get(chrlength_table,chrom);
#else
    chroffset = (Univcoord_T) Tableuint_get(chrlength_table,chrom);
#endif
    universalpos1 = chroffset + (Univcoord_T) Segmentpos_chrpos1(segmentpos);
    universalpos2 = chroffset + (Univcoord_T) Segmentpos_chrpos2(segmentpos);

    /* Print as 1-based, inclusive [a,b] */
    if (Segmentpos_revcompp(segmentpos) == true) {
      fprintf(textfp,"%s\t%llu..%llu\t%s:%u..%u\t%u",
	      accessions[i],(unsigned long long) universalpos2+1U,(unsigned long long) universalpos1,
	      Chrom_string(chrom),Segmentpos_chrpos2(segmentpos)+1U,Segmentpos_chrpos1(segmentpos),
	      Segmentpos_length(segmentpos));
    } else {
      fprintf(textfp,"%s\t%llu..%llu\t%s:%u..%u\t%u",
	      accessions[i],(unsigned long long) universalpos1+1U,(unsigned long long) universalpos2,
	      Chrom_string(chrom),Segmentpos_chrpos1(segmentpos)+1U,Segmentpos_chrpos2(segmentpos),
	      Segmentpos_length(segmentpos));
    }


#if 0
    if (Segmentpos_type(segmentpos) > 0) {
      fprintf(textfp,"\t%s",contigtypes[Segmentpos_type(segmentpos)]);
    }
#endif

    fprintf(textfp,"\n");

    /* Store as 0-based, inclusive [a,b] */
    labellist = List_push(labellist,(void *) accessions[i]);
    if (Segmentpos_revcompp(segmentpos) == true) {
      /* The negative sign in the interval is the indication that the
	 contig was reverse complement */
      intervallist = List_push(intervallist, 
			       (void *) Univinterval_new(universalpos2-1U,universalpos1,
							 Segmentpos_type(segmentpos)));
    } else {
      intervallist = List_push(intervallist, 
			       (void *) Univinterval_new(universalpos1,universalpos2-1U,
							 Segmentpos_type(segmentpos)));
    }

#if 1
    /* IIT version 1.  Need to indicate revcomp with '-' as first character. */
    sprintf(seglength,"%u",Segmentpos_length(segmentpos));
    annot = (char *) CALLOC(strlen(seglength)+1,sizeof(char));
    strcpy(annot,seglength);
#else
    /* IIT versions >= 2.  Not used anymore, since universal IITs are restricted to version 1. */
    annot = (char *) CALLOC(1,sizeof(char));
    annot[0] = '\0';
#endif
    annotlist = List_push(annotlist,(void *) annot);

  }
  fclose(textfp);

#if 0
  FREE(contigtypes);
#endif
  FREE(accessions);
  intervallist = List_reverse(intervallist);
  /* contigtypelist = List_reverse(contigtypelist); -- Done by caller */ 
  labellist = List_reverse(labellist);
  annotlist = List_reverse(annotlist);

  /* Write contig IIT file */
  divstring = (char *) CALLOC(1,sizeof(char));
  divstring[0] = '\0';
  divlist = List_push(NULL,divstring);

  intervaltable = Table_new(65522,Table_string_compare,Table_string_hash);
  labeltable = Table_new(65522,Table_string_compare,Table_string_hash);
  annottable = Table_new(65522,Table_string_compare,Table_string_hash);

  Table_put(intervaltable,(void *) divstring,intervallist);
  Table_put(labeltable,(void *) divstring,labellist);
  Table_put(annottable,(void *) divstring,annotlist);

  iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			    strlen(fileroot)+strlen(".contig.iit")+1,sizeof(char));
  sprintf(iitfile,"%s/%s.contig.iit",genomesubdir,fileroot);

#if 0
  debug(
	for (p = contigtypelist; p != NULL; p = List_next(p)) {
	  printf("Type %s\n",(char *) List_head(p));
	}
	);
#endif
  IIT_write_univ(iitfile,divlist,contigtypelist,
		 intervaltable,labeltable,annottable,
		 coord_values_8p,/*label_pointers_8p*/false,/*annot_pointers_8p*/false);
  FREE(iitfile);

  List_free(&divlist);
  FREE(divstring);

  Table_free(&annottable);
  Table_free(&labeltable);
  Table_free(&intervaltable);

  for (p = annotlist; p != NULL; p = List_next(p)) {
    annot = (char *) List_head(p);
    FREE(annot);
  }
  List_free(&annotlist);

  /* Labels (accessions) are freed by accsegmentpos_table_gc */
  List_free(&labellist);

  for (p = intervallist; p != NULL; p = List_next(p)) {
    interval = (Univinterval_T) List_head(p);
    Univinterval_free(&interval);
  }
  List_free(&intervallist);

  return;
}

/************************************************************************/


#if 0
static void
stringlist_gc (List_T *list) {
  List_T p;
  char *string;

  for (p = *list; p != NULL; p = List_next(p)) {
    string = (char *) List_head(p);
    FREE(string);
  }
  List_free(&(*list));
  return;
}
#endif


static void
chrlength_table_gc (Table_chrpos_T *chrlength_table) {
  /* Don't free chrom entries in table, because they are freed by Segmentpos_free */
#ifdef HAVE_64_BIT
  Tableuint8_free(&(*chrlength_table));
#else
  Tableuint_free(&(*chrlength_table));
#endif

  return;
}

static void
accsegmentpos_table_gc (Table_T *accsegmentpos_table) {
  int n, i = 0;
  char *accession;
  Segmentpos_T segmentpos;
  void **keys, **values;

  /* For some reason, this was failing on some computers, perhaps
     because we weren't checking for n > 0 */
  if ((n = Table_length(*accsegmentpos_table)) > 0) {
    keys = Table_keys(*accsegmentpos_table,NULL);
    values = Table_values(*accsegmentpos_table,NULL);
    for (i = 0; i < n; i++) {
      accession = (char *) keys[i];
      FREE(accession);
    }
    for (i = 0; i < n; i++) {
      segmentpos = (Segmentpos_T) values[i];
      Segmentpos_free(&segmentpos);
    }
    FREE(values);
    FREE(keys);
  }

  Table_free(&(*accsegmentpos_table));
  return;
}

#if 0
static char *
remove_slashes (char *buffer) {
  char *copy, *p;

  copy = (char *) CALLOC(strlen(buffer)+1,sizeof(char));
  strcpy(copy,buffer);
  p = copy;
  while (*p != '\0') {
    if (*p == '/') {
      *p = '_';
    }
    p++;
  }
  return copy;
}
#endif


static int
add_compression_type (char *string) {
  if (!strcmp(string,"none")) {
    compression_types = NO_COMPRESSION;
    return 0;
  } else if (!strcmp(string,"all")) {
    compression_types = BITPACK64_COMPRESSION;
    return 1;
  } else {
    if (!strcmp(string,"bitpack64")) {
      compression_types |= BITPACK64_COMPRESSION;
    } else {
      fprintf(stderr,"Don't recognize compression type %s\n",string);
      fprintf(stderr,"Allowed values are: none, all, bitpack64\n");
      exit(9);
    }
    return 1;
  }
}


static char CHARTABLE[4] = {'A','C','G','T'};
static char *mode_prefix = ".";
/* static char *mode_prefix = ".metct."; */


#ifdef __STRICT_ANSI__
int getopt (int argc, char *const argv[], const char *optstring);
#endif

int 
main (int argc, char *argv[]) {
  int ncontigs;
  Table_T accsegmentpos_table;
  int shmid;

  FILE *fp;
  char *key, **keys, chrname[1024], chrname_alt[1024], Buffer[1024];
  Tableuint_T chrorder_table = NULL;
  unsigned int order;
  int i;

  Table_chrpos_T chrlength_table;
  List_T contigtypelist = NULL, p;
  Genome_T genomecomp, genomebits;
  Univ_IIT_T chromosome_iit, contig_iit;
  char *typestring;
  Univcoord_T n, genomelength, totalnts;
  char *chromosomefile, *iitfile, *positionsfile_high, *positionsfile_low, interval_char;
  char *sarrayfile, *lcpexcfile, *lcpguidefile;
  char *rankfile, *permuted_sarray_file; /* temporary files */

  /* For compressed suffix array */
  char *csaptrfiles[5], *csacompfiles[5], *sasampleqfile, *sasamplesfile, *saindex0file;

  char *childexcfile, *childguidefile;
  char *lcpchilddcfile;
#ifdef USE_SEPARATE_BUCKETS
  char *indexiptrsfile, *indexicompfile, *indexjptrsfile, *indexjcompfile;
#else
  char *indexijptrsfile, *indexijcompfile;
#endif
  UINT4 start, end;
  UINT4 *SA, *lcp, *child;
  UINT4 nbytes;
  Filenames_T filenames;
  Chrpos_T seglength;
  bool coord_values_8p;

  unsigned char *discrim_chars;
  unsigned char *lcp_bytes;
  UINT4 *lcp_guide, *lcp_exceptions;
  int n_lcp_exceptions;

  unsigned char *lcpchilddc;
  UINT4 *child_guide, *child_exceptions;
  int n_child_exceptions;

  int sa_fd, lcpchilddc_fd;
  size_t sa_len, lcpchilddc_len, lcpguide_len, lcpexc_len, childguide_len, childexc_len;
  double seconds;

#ifdef HAVE_64_BIT
  UINT8 noffsets;
#endif

  int c;
  extern int optind;
  extern char *optarg;
  char *string;

  while ((c = getopt(argc,argv,"F:D:d:z:k:q:ArlGUNHOPSLCXYWw:e:Ss:n:m9")) != -1) {
    switch (c) {
    case 'F': sourcedir = optarg; break;
    case 'D': destdir = optarg; break;
    case 'd': fileroot = optarg; break;

    case 'z':
      compression_types = NO_COMPRESSION; /* Initialize */
      string = strtok(optarg,",");
      if (add_compression_type(string) != 0) {
	while ((string = strtok(NULL,",")) != NULL && add_compression_type(string) != 0) {
	}
      }
      break;

    case 'k': index1part = atoi(optarg);
      if (index1part > MAXIMUM_KMER) {
	fprintf(stderr,"The choice of k-mer size must be %d or less\n",MAXIMUM_KMER);
	exit(9);
      }
      break;
    case 'q': index1interval = atoi(optarg); break;

    case 'A': action = AUXFILES; break;
    case 'r': rawp = true; break;
    case 'l': genome_lc_p = true; break;
    case 'G': action = GENOME; break;
    case 'U': action = UNSHUFFLE; break;
    case 'N': action = COUNT; break;
    case 'H': huge_offsets_p = true; break;
    case 'O': action = OFFSETS; break;
    case 'P': action = POSITIONS; break;
    case 'S': action = SUFFIX_ARRAY; break;
    case 'L': action = LCP_CHILD; break;
    case 'C': action = COMPRESSED_SUFFIX_ARRAY; break;
    case 'X': action = ARRAY_UNCOMPRESS; break;
    case 'Y': action = CHILD_UNCOMPRESS; break;
    case 'W': writefilep = true; break;
    case 'w': wraplength = atoi(optarg); break;
    case 'e': nmessages = atoi(optarg); break;

    case 's': 
      if (!strcmp(optarg,"none")) {
	divsort = NO_SORT;
      } else if (!strcmp(optarg,"alpha")) {
	divsort = ALPHA_SORT;
      } else if (!strcmp(optarg,"numeric-alpha")) {
	divsort = NUMERIC_ALPHA_SORT;
      } else if (!strcmp(optarg,"chrom")) {
	divsort = CHROM_SORT;
      } else if (!strcmp(optarg,"names")) {
	divsort = FILENAME_SORT;
      } else {
	fprintf(stderr,"Don't recognize sort type %s.  Allowed values are none, alpha, or chrom.",optarg);
	exit(9);
      }
      break;

    case 'n': sortfilename = optarg; break;

    case 'm': mask_lowercase_p = true; break;

    case '9': check_compiler_assumptions(); return 0; break;

    default: fprintf(stderr,"Unknown flag %c\n",c); exit(9);

    }
  }
  argc -= optind;
  argv += optind;

  if (index1interval == 3) {
    interval_char = '3';
  } else if (index1interval == 2) {
    interval_char = '2';
  } else if (index1interval == 1) {
    interval_char = '1';
  } else {
    fprintf(stderr,"Selected indexing interval %d is not allowed.  Only values allowed are 3, 2, or 1\n",index1interval);
    exit(9);
  }

  if (action != UNSHUFFLE) {
    if (fileroot == NULL) {
      fprintf(stderr,"Missing name of genome database.  Must specify with -d flag.\n");
      exit(9);
    }
  }

  if (action == AUXFILES) {
    /* Usage: cat <fastafile> | gmapindex [-F <sourcedir>] [-D <destdir>] -d <dbname> -A
       Requires <fastafile> in appropriate format
       Writes <destdir>/<dbname>.chromosome and <destdir>/<dbname>.contig files 
       and corresponding .iit files */

    if (divsort == FILENAME_SORT) {
      if (sortfilename == NULL) {
	fprintf(stderr,"For sorting by names file, need to provide file to -n flag");
	exit(9);
      } else if ((fp = fopen(sortfilename,"r")) == NULL) {
	fprintf(stderr,"Unable to open file %s provided to -n flag",sortfilename);
	exit(9);
      } else {
	chrorder_table = Tableuint_new(65522,Table_string_compare,Table_string_hash);
	order = 1;
	while (fgets(Buffer,1024,fp) != NULL) {
   	  if (sscanf(Buffer,"%s %s",chrname,chrname_alt) == 2) {	  
	    key = (char *) CALLOC(strlen(chrname_alt)+1,sizeof(char));
	    strcpy(key,chrname_alt);
	    Tableuint_put(chrorder_table,(void *) key,order);
          } else if (sscanf(Buffer,"%s",chrname) == 1) {
	    key = (char *) CALLOC(strlen(chrname)+1,sizeof(char));
	    strcpy(key,chrname);
	    Tableuint_put(chrorder_table,(void *) key,order);
          } else {
	    fprintf(stderr,"Unable to parse line %s\n",Buffer);
	  }
          order += 1;
	}
	fclose(fp);
      }
    }

    if (getc(stdin) != '>') {
      fprintf(stderr,"Expected file to start with '>'\n");
      exit(9);
    }

    /* Holds contigs.  keys are strings; values are structs. */
    accsegmentpos_table = Table_new(65522,Table_string_compare,Table_string_hash);
    /* Hold chromosomes.  keys are Chrom_Ts; values are uints. */
#ifdef HAVE_64_BIT
    chrlength_table = Tableuint8_new(65522,Chrom_compare_table,Chrom_hash_table);
#else
    chrlength_table = Tableuint_new(65522,Chrom_compare_table,Chrom_hash_table);
#endif

#if 0
    /* No longer supporting strains */
    /* keys are strings; values are ints */
    contigtype_table = Tableuint_new(100,Table_string_compare,Table_string_hash);

    refstrain = read_strain_from_coordsfile(coordsfile);
    fprintf(stderr,"Reference strain is %s\n",refstrain);
    contigtypelist = List_push(NULL,refstrain);
#endif
    /* The zeroth type is empty */
    typestring = (char *) CALLOC(1,sizeof(char));
    typestring[0] = '\0';
    contigtypelist = List_push(NULL,typestring);

    ncontigs = 0;
    totalnts = 0;
    while (process_sequence_aux(&seglength,accsegmentpos_table,chrlength_table,chrorder_table,fileroot,ncontigs) == true) {
      if (totalnts + seglength < totalnts) {
	/* Exceeds 32 bits */
	fprintf(stderr,"The total length of genomic sequence exceeds 2^32 = 4,294,967,296 bp, which the GMAP index format cannot handle\n");
	exit(9);
      } else {
	totalnts = totalnts + seglength;
      }
      ncontigs++;
    }
    fprintf(stderr,"Total genomic length = %llu bp\n",(unsigned long long) totalnts);

    if (ncontigs == 0) {
      fprintf(stderr,"No contig information was provided to gmapindex\n");
      exit(9);
    }

#ifdef HAVE_64_BIT
    if (totalnts > 4294967295) {
      coord_values_8p = true;
    } else {
      coord_values_8p = false;
    }
#else
    coord_values_8p = false;
#endif

    write_chromosome_file(destdir,fileroot,chrlength_table,coord_values_8p);
    write_contig_file(destdir,fileroot,accsegmentpos_table,chrlength_table,contigtypelist,coord_values_8p);

    for (p = contigtypelist; p != NULL; p = List_next(p)) {
      typestring = (char *) List_head(p);
      FREE(typestring);
    }
    List_free(&contigtypelist);

    chrlength_table_gc(&chrlength_table);
    accsegmentpos_table_gc(&accsegmentpos_table);

    if (chrorder_table != NULL) {
      keys = (char **) Tableuint_keys(chrorder_table,NULL);
      for (i = 0; i < Tableuint_length(chrorder_table); i++) {
	FREE(keys[i]);
      }
      FREE(keys);
      Tableuint_free(&chrorder_table);
    }

  } else if (action == GENOME) {
    /* Usage: cat <fastafile> | gmapindex [-F <sourcedir>] [-D <destdir>] -d <dbname> -G
       Requires <fastafile> in appropriate format and <sourcedir>/<dbname>.chromosome.iit 
       and <sourcedir>/<dbname>.contig.iit files.
       Creates <destdir>/<dbname>.genomecomp (in blocks of 32) */

    chromosomefile = (char *) CALLOC(strlen(sourcedir)+strlen("/")+
				     strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
    sprintf(chromosomefile,"%s/%s.chromosome.iit",sourcedir,fileroot);
    if ((chromosome_iit = Univ_IIT_read(chromosomefile,/*readonlyp*/true,/*add_iit_p*/false)) == NULL) {
      fprintf(stderr,"IIT file %s is not valid\n",chromosomefile);
      exit(9);
    }
    genomelength = Univ_IIT_genomelength(chromosome_iit,/*with_circular_alias_p*/true);
    FREE(chromosomefile);

    iitfile = (char *) CALLOC(strlen(sourcedir)+strlen("/")+
			      strlen(fileroot)+strlen(".contig.iit")+1,sizeof(char));
    sprintf(iitfile,"%s/%s.contig.iit",sourcedir,fileroot);
    if ((contig_iit = Univ_IIT_read(iitfile,/*readonlyp*/true,/*add_iit_p*/false)) == NULL) {
      fprintf(stderr,"IIT file %s is not valid\n",iitfile);
      exit(9);
    }
    FREE(iitfile);
    
    if (Univ_IIT_ntypes(contig_iit) == 1) {
      /* index1part needed only if writing an uncompressed genome using a file */
      Genome_write_comp32(destdir,fileroot,stdin,contig_iit,/*altstrain_iit*/NULL,
			  chromosome_iit,genome_lc_p,rawp,writefilep,
			  genomelength,index1part,nmessages);
    } else if (Univ_IIT_ntypes(contig_iit) > 1) {
      fprintf(stderr,"GMAPINDEX no longer supports alternate strains\n");
      abort();
    }

    Univ_IIT_free(&chromosome_iit);
    Univ_IIT_free(&contig_iit);

  } else if (action == UNSHUFFLE) {
    /* Usage: cat <fastafile> | gmapindex [-F <sourcedir>] [-D <destdir>] -d <dbname> -G
       Requires <fastafile> in appropriate format and <sourcedir>/<dbname>.chromosome.iit 
       and <sourcedir>/<dbname>.contig.iit files.
       Creates <destdir>/<dbname>.genomebits128 (in blocks of 128) */

    Compress_unshuffle_bits128(stdout,stdin);

  } else if (action == COUNT) {
    /* Usage: cat <genomefile> | gmapindex [-F <sourcedir>] [-D <destdir>] -d <dbname> -C */

#ifndef HAVE_64_BIT
    printf("0\n");
#else
    chromosomefile = (char *) CALLOC(strlen(sourcedir)+strlen("/")+
				     strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
    sprintf(chromosomefile,"%s/%s.chromosome.iit",sourcedir,fileroot);
    if ((chromosome_iit = Univ_IIT_read(chromosomefile,/*readonlyp*/true,/*add_iit_p*/false)) == NULL) {
      fprintf(stderr,"IIT file %s is not valid\n",chromosomefile);
      exit(9);
    }
    FREE(chromosomefile);

    noffsets = Indexdb_count_offsets(stdin,chromosome_iit,index1part,index1interval,
				     genome_lc_p,fileroot,mask_lowercase_p);
    printf("%llu\n",(unsigned long long) noffsets);

    Univ_IIT_free(&chromosome_iit);
#endif

  } else if (action == OFFSETS) {
    /* Usage: gmapindex [-F <sourcedir>] [-D <destdir>] -d <dbname> -O <genomefile>
       Creates <destdir>/<dbname>.idxoffsets */

    if (argc == 0) {
      fp = stdin;
    } else if ((fp = fopen(argv[0],"rb")) == NULL) {
      fprintf(stderr,"Could not open file %s\n",argv[0]);
      exit(9);
    }

    chromosomefile = (char *) CALLOC(strlen(sourcedir)+strlen("/")+
				     strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
    sprintf(chromosomefile,"%s/%s.chromosome.iit",sourcedir,fileroot);
    if ((chromosome_iit = Univ_IIT_read(chromosomefile,/*readonlyp*/true,/*add_iit_p*/false)) == NULL) {
      fprintf(stderr,"IIT file %s is not valid\n",chromosomefile);
      exit(9);
    }
    FREE(chromosomefile);

    if (huge_offsets_p == false) {
      fprintf(stderr,"Offset compression types:");
      if ((compression_types & BITPACK64_COMPRESSION) != 0) {
	fprintf(stderr," bitpack64");
      }
      fprintf(stderr,"\n");
      
      Indexdb_write_offsets(destdir,interval_char,fp,chromosome_iit,
			    index1part,index1interval,
			    genome_lc_p,fileroot,mask_lowercase_p,compression_types);
    } else {
      fprintf(stderr,"Offset compression types:");
      if ((compression_types & BITPACK64_COMPRESSION) != 0) {
	fprintf(stderr," bitpack64");
      }
      fprintf(stderr,"\n");
      
      Indexdb_write_offsets_huge(destdir,interval_char,fp,chromosome_iit,
				 index1part,index1interval,
				 genome_lc_p,fileroot,mask_lowercase_p,compression_types);
    }

    if (argc > 0) {
      fclose(fp);
    }

    Univ_IIT_free(&chromosome_iit);

  } else if (action == POSITIONS) {
    /* Usage: gmapindex [-F <sourcedir>] [-D <destdir>] -d <dbname> -P <genomefile>
       Requires <sourcedir>/<dbname>.idxoffsets.
       Creates <destdir>/<dbname>.idxpositions */

    if (argc == 0) {
      fp = stdin;
    } else if ((fp = fopen(argv[0],"rb")) == NULL) {
      fprintf(stderr,"Could not open file %s\n",argv[0]);
      exit(9);
    }

    chromosomefile = (char *) CALLOC(strlen(sourcedir)+strlen("/")+
				     strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
    sprintf(chromosomefile,"%s/%s.chromosome.iit",sourcedir,fileroot);
    if ((chromosome_iit = Univ_IIT_read(chromosomefile,/*readonlyp*/true,/*add_iit_p*/false)) == NULL) {
      fprintf(stderr,"IIT file %s is not valid\n",chromosomefile);
      exit(9);
    }
    genomelength = Univ_IIT_genomelength(chromosome_iit,/*with_circular_alias_p*/true);
    FREE(chromosomefile);

    filenames = Indexdb_get_filenames(&compression_type,&index1part,&index1interval,
				      sourcedir,fileroot,IDX_FILESUFFIX,/*snps_root*/NULL,
				      /*required_index1part*/index1part,
				      /*required_interval*/index1interval,/*offsets_only_p*/true);

    if (Univ_IIT_coord_values_8p(chromosome_iit) == true) {
      coord_values_8p = true;

      positionsfile_high = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
					   strlen(".")+strlen(IDX_FILESUFFIX)+
					   /*for kmer*/2+/*for interval char*/1+
					   strlen(POSITIONS_HIGH_FILESUFFIX)+1,sizeof(char));
      sprintf(positionsfile_high,"%s/%s.%s%02d%c%s",
	    destdir,fileroot,IDX_FILESUFFIX,index1part,interval_char,POSITIONS_HIGH_FILESUFFIX);

      positionsfile_low = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
					   strlen(".")+strlen(IDX_FILESUFFIX)+
					   /*for kmer*/2+/*for interval char*/1+
					   strlen(POSITIONS_LOW_FILESUFFIX)+1,sizeof(char));
      sprintf(positionsfile_low,"%s/%s.%s%02d%c%s",
	    destdir,fileroot,IDX_FILESUFFIX,index1part,interval_char,POSITIONS_LOW_FILESUFFIX);

    } else {
      coord_values_8p = false;
      
      positionsfile_high = (char *) NULL;

      positionsfile_low = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
					   strlen(".")+strlen(IDX_FILESUFFIX)+
					   /*for kmer*/2+/*for interval char*/1+
					   strlen(POSITIONS_LOW_FILESUFFIX)+1,sizeof(char));
      sprintf(positionsfile_low,"%s/%s.%s%02d%c%s",
	    destdir,fileroot,IDX_FILESUFFIX,index1part,interval_char,POSITIONS_LOW_FILESUFFIX);
    }

    if (huge_offsets_p == false) {
      Indexdb_write_positions(positionsfile_high,positionsfile_low,filenames->pointers_filename,
			      filenames->offsets_filename,fp,chromosome_iit,
			      index1part,index1interval,genomelength,
			      genome_lc_p,writefilep,fileroot,mask_lowercase_p,
			      compression_type,coord_values_8p);
    } else {
      Indexdb_write_positions_huge(positionsfile_high,positionsfile_low,filenames->pages_filename,filenames->pointers_filename,
				   filenames->offsets_filename,fp,chromosome_iit,
				   index1part,index1interval,genomelength,
				   genome_lc_p,writefilep,fileroot,mask_lowercase_p,
				   compression_type,coord_values_8p);
    }

    if (argc > 0) {
      fclose(fp);
    }

    Filenames_free(&filenames);

    FREE(positionsfile_high);
    FREE(positionsfile_low);
    Univ_IIT_free(&chromosome_iit);

  } else if (action == SUFFIX_ARRAY) {
    /* Usage: gmapindex [-F <sourcedir>] [-D <destdir>] -d <dbname> -S
       Creates <destdir>/<dbname>.sarray, .lcp, and .saindex */
    chromosomefile = (char *) CALLOC(strlen(sourcedir)+strlen("/")+
				     strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
    sprintf(chromosomefile,"%s/%s.chromosome.iit",sourcedir,fileroot);
    if ((chromosome_iit = Univ_IIT_read(chromosomefile,/*readonlyp*/true,/*add_iit_p*/false)) == NULL) {
      fprintf(stderr,"IIT file %s is not valid\n",chromosomefile);
      exit(9);
    }
    FREE(chromosomefile);

    /* Suffix array */
    genomelength = Univ_IIT_genomelength(chromosome_iit,/*with_circular_alias_p*/true);
    Univ_IIT_free(&chromosome_iit);

    fprintf(stderr,"Genome length is %llu\n",(unsigned long long) genomelength);
    if (genomelength > 4294967295) {
      fprintf(stderr,"Suffix arrays not yet supported for large genomes with more than 2^32 bp.  Will use hash table only.\n");

    } else {
      fprintf(stderr,"Building suffix array\n");
      sarrayfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".sarray")+1,sizeof(char));
      sprintf(sarrayfile,"%s/%s.sarray",destdir,fileroot);

      genomecomp = Genome_new(sourcedir,fileroot,/*snps_root*/NULL,/*genometype*/GENOME_OLIGOS,
			      /*uncompressedp*/false,/*access*/USE_MMAP_ONLY,/*sharedp*/false);
      Sarray_write_array(sarrayfile,genomecomp,genomelength);
      
      /* Bucket array */
#ifdef USE_SEPARATE_BUCKETS
      indexiptrsfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".saindexi64meta")+1,sizeof(char));
      sprintf(indexiptrsfile,"%s/%s.saindexi64meta",destdir,fileroot);
      indexicompfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".saindexi64strm")+1,sizeof(char));
      sprintf(indexicompfile,"%s/%s.saindexi64strm",destdir,fileroot);
      indexjptrsfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".saindexj64meta")+1,sizeof(char));
      sprintf(indexjptrsfile,"%s/%s.saindexj64meta",destdir,fileroot);
      indexjcompfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".saindexj64strm")+1,sizeof(char));
      sprintf(indexjcompfile,"%s/%s.saindexj64strm",destdir,fileroot);
      
      Sarray_write_index_separate(indexiptrsfile,indexicompfile,indexjptrsfile,indexjcompfile,
				  sarrayfile,genomecomp,genomelength,/*compressp*/true);
      FREE(indexjcompfile);
      FREE(indexjptrsfile);
      FREE(indexicompfile);
      FREE(indexiptrsfile);
#else
      indexijptrsfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".saindex64meta")+1,sizeof(char));
      sprintf(indexijptrsfile,"%s/%s.saindex64meta",destdir,fileroot);
      indexijcompfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".saindex64strm")+1,sizeof(char));
      sprintf(indexijcompfile,"%s/%s.saindex64strm",destdir,fileroot);
      Sarray_write_index_interleaved(indexijptrsfile,indexijcompfile,
				     sarrayfile,genomecomp,genomelength,/*compressp*/true,
				     CHARTABLE);
      FREE(indexijcompfile);
      FREE(indexijptrsfile);
#endif

      Genome_free(&genomecomp);
      FREE(sarrayfile);
    }

  } else if (action == LCP_CHILD) {
    /* Usage: gmapindex [-F <sourcedir>] [-D <destdir>] -d <dbname> -L
       Creates <destdir>/<dbname>.lcp, .saindex, and .rank (needed by COMPRESSED_SUFFIX_ARRAY) */

    chromosomefile = (char *) CALLOC(strlen(sourcedir)+strlen("/")+
				     strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
    sprintf(chromosomefile,"%s/%s.chromosome.iit",sourcedir,fileroot);
    if ((chromosome_iit = Univ_IIT_read(chromosomefile,/*readonlyp*/true,/*add_iit_p*/false)) == NULL) {
      fprintf(stderr,"IIT file %s is not valid\n",chromosomefile);
      exit(9);
    }
    FREE(chromosomefile);

    genomelength = Univ_IIT_genomelength(chromosome_iit,/*with_circular_alias_p*/true);
    Univ_IIT_free(&chromosome_iit);

    if (genomelength > 4294967295) {
      /* Warning message already printed for SUFFIX_ARRAY */
      /* fprintf(stderr,"Suffix arrays not yet supported for large genomes with more than 2^32 bp.  Will use hash table only.\n"); */
    } else {
      fprintf(stderr,"Building LCP array\n");

      n = genomelength;

      /* No need to mmap SA anymore */
      sarrayfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".sarray")+1,sizeof(char));
      sprintf(sarrayfile,"%s/%s.sarray",destdir,fileroot);

      /* Required for computing LCP, but uses non-SIMD instructions */
      genomebits = Genome_new(sourcedir,fileroot,/*snps_root*/NULL,/*genometype*/GENOME_BITS,
			      /*uncompressedp*/false,/*access*/USE_MMAP_ONLY,/*sharedp*/false);
      Genome_hr_setup(Genome_blocks(genomebits),/*snp_blocks*/NULL,
		      /*query_unk_mismatch_p*/false,/*genome_unk_mismatch_p*/false,
		      /*mode*/STANDARD);


      rankfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".rank")+1,sizeof(char));
      sprintf(rankfile,"%s/%s.rank",destdir,fileroot);
      permuted_sarray_file = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".permuted_sarray")+1,sizeof(char));
      sprintf(permuted_sarray_file,"%s/%s.permuted_sarray",destdir,fileroot);
      lcp = Sarray_compute_lcp(rankfile,permuted_sarray_file,sarrayfile,n);
      FREE(permuted_sarray_file);
      FREE(rankfile);

      Genome_free(&genomebits);

      /* Write lcp exceptions/guide, but return lcp_bytes */
      lcpexcfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".salcpexc")+1,sizeof(char));
      sprintf(lcpexcfile,"%s/%s.salcpexc",destdir,fileroot);
      lcpguidefile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".salcpguide1024")+1,sizeof(char));
      sprintf(lcpguidefile,"%s/%s.salcpguide1024",destdir,fileroot);

      lcp_bytes = Bytecoding_write_exceptions_only(lcpexcfile,lcpguidefile,lcp,genomelength,/*guide_interval*/1024);

      FREE(lcpguidefile);
      FREE(lcpexcfile);

      FREE(lcp);			/* Use lcp_bytes, which are more memory-efficient than lcp */


      fprintf(stderr,"Building DC array\n");

      /* Assume we have lcp_bytes already in memory.  Don't need to use guide for speed. */
      lcpguidefile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".salcpguide1024")+1,sizeof(char));
      sprintf(lcpguidefile,"%s/%s.salcpguide1024",destdir,fileroot);
      lcp_guide = (UINT4 *) Access_allocate(&shmid,&lcpguide_len,&seconds,lcpguidefile,sizeof(UINT4),/*sharedp*/false);
      FREE(lcpguidefile);

      lcpexcfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".salcpexc")+1,sizeof(char));
      sprintf(lcpexcfile,"%s/%s.salcpexc",destdir,fileroot);
      lcp_exceptions = (UINT4 *) Access_allocate(&shmid,&lcpexc_len,&seconds,lcpexcfile,sizeof(UINT4),/*sharedp*/false);
      n_lcp_exceptions = lcpexc_len/(sizeof(UINT4) + sizeof(UINT4));
      FREE(lcpexcfile);

      genomecomp = Genome_new(sourcedir,fileroot,/*snps_root*/NULL,/*genometype*/GENOME_OLIGOS,
			      /*uncompressedp*/false,/*access*/USE_MMAP_ONLY,/*sharedp*/false);

      /* Compute discriminating chars (DC) array */
      discrim_chars = Sarray_discriminating_chars(&nbytes,sarrayfile,genomecomp,lcp_bytes,lcp_guide,
						  lcp_exceptions,/*guide_interval*/1024,n,CHARTABLE);
      FREE(sarrayfile);
      Genome_free(&genomecomp);
      /* No need to munmap SA anymore */


      fprintf(stderr,"Building child array\n");

      /* Compute child array (relative values) */
      child = Sarray_compute_child(lcp_bytes,lcp_guide,lcp_exceptions,n);
      FREE(lcp_exceptions);
      FREE(lcp_guide);

      /* Write combined lcpchilddc file */
      lcpchilddcfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".salcpchilddc")+1,sizeof(char));
      sprintf(lcpchilddcfile,"%s/%s.salcpchilddc",destdir,fileroot);
      childexcfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".sachildexc")+1,sizeof(char));
      sprintf(childexcfile,"%s/%s.sachildexc",destdir,fileroot);
      childguidefile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".sachildguide1024")+1,sizeof(char));
      sprintf(childguidefile,"%s/%s.sachildguide1024",destdir,fileroot);
      Bytecoding_write_lcpchilddc(lcpchilddcfile,childexcfile,childguidefile,child,
				  discrim_chars,lcp_bytes,genomelength,/*guide_interval*/1024);
      FREE(childguidefile);
      FREE(childexcfile);
      FREE(lcpchilddcfile);

      FREE(child);
      FREE(discrim_chars);
      FREE(lcp_bytes);
    }

  } else if (action == COMPRESSED_SUFFIX_ARRAY) {
    /* Usage: gmapindex [-F <sourcedir>] [-D <destdir>] -d <dbname> -C
       Creates <destdir>/<dbname>.lcp and .csa.  Removes .sarray and .inverse_sarray */

    chromosomefile = (char *) CALLOC(strlen(sourcedir)+strlen("/")+
				     strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
    sprintf(chromosomefile,"%s/%s.chromosome.iit",sourcedir,fileroot);
    if ((chromosome_iit = Univ_IIT_read(chromosomefile,/*readonlyp*/true,/*add_iit_p*/false)) == NULL) {
      fprintf(stderr,"IIT file %s is not valid\n",chromosomefile);
      exit(9);
    }
    FREE(chromosomefile);

    genomelength = Univ_IIT_genomelength(chromosome_iit,/*with_circular_alias_p*/true);
    Univ_IIT_free(&chromosome_iit);

    if (genomelength > 4294967295) {
      /* Warning message already printed for SUFFIX_ARRAY */
      /* fprintf(stderr,"Suffix arrays not yet supported for large genomes with more than 2^32 bp.  Will use hash table only.\n"); */
    } else {
      fprintf(stderr,"Building compressed suffix array\n");

      /* No need to mmap SA anymore */
      csaptrfiles[0] = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".csaAmeta")+1,sizeof(char));
      sprintf(csaptrfiles[0],"%s/%s.csaAmeta",destdir,fileroot);
      csaptrfiles[1] = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".csaCmeta")+1,sizeof(char));
      sprintf(csaptrfiles[1],"%s/%s.csaCmeta",destdir,fileroot);
      csaptrfiles[2] = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".csaGmeta")+1,sizeof(char));
      sprintf(csaptrfiles[2],"%s/%s.csaGmeta",destdir,fileroot);
      csaptrfiles[3] = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".csaTmeta")+1,sizeof(char));
      sprintf(csaptrfiles[3],"%s/%s.csaTmeta",destdir,fileroot);
      csaptrfiles[4] = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".csaXmeta")+1,sizeof(char));
      sprintf(csaptrfiles[4],"%s/%s.csaXmeta",destdir,fileroot);

      csacompfiles[0] = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".csaAstrm")+1,sizeof(char));
      sprintf(csacompfiles[0],"%s/%s.csaAstrm",destdir,fileroot);
      csacompfiles[1] = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".csaCstrm")+1,sizeof(char));
      sprintf(csacompfiles[1],"%s/%s.csaCstrm",destdir,fileroot);
      csacompfiles[2] = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".csaGstrm")+1,sizeof(char));
      sprintf(csacompfiles[2],"%s/%s.csaGstrm",destdir,fileroot);
      csacompfiles[3] = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".csaTstrm")+1,sizeof(char));
      sprintf(csacompfiles[3],"%s/%s.csaTstrm",destdir,fileroot);
      csacompfiles[4] = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".csaXstrm")+1,sizeof(char));
      sprintf(csacompfiles[4],"%s/%s.csaXstrm",destdir,fileroot);

      sasampleqfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".sasampleq")+1,sizeof(char));
      sprintf(sasampleqfile,"%s/%s.sasampleq",destdir,fileroot);
      sasamplesfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".sasamples")+1,sizeof(char));
      sprintf(sasamplesfile,"%s/%s.sasamples",destdir,fileroot);
      saindex0file = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".saindex0")+1,sizeof(char));
      sprintf(saindex0file,"%s/%s.saindex0",destdir,fileroot);


#if 0
      csafile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".csa")+1,sizeof(char));
      sprintf(csafile,"%s/%s.csa",destdir,fileroot);
#endif


      sarrayfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".sarray")+1,sizeof(char));
      sprintf(sarrayfile,"%s/%s.sarray",destdir,fileroot);
      rankfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".rank")+1,sizeof(char));
      sprintf(rankfile,"%s/%s.rank",destdir,fileroot);

      genomecomp = Genome_new(sourcedir,fileroot,/*snps_root*/NULL,/*genometype*/GENOME_OLIGOS,
			      /*uncompressedp*/false,/*access*/USE_MMAP_ONLY,/*sharedp*/false);
      Sarray_write_csa(csaptrfiles,csacompfiles,sasampleqfile,sasamplesfile,saindex0file,
		       sarrayfile,rankfile,genomecomp,genomelength,CHARTABLE);
      FREE(genomecomp);

      remove(rankfile);		/* Need to delete remove(rankfile) from Sarray_compute_lcp */
      remove(sarrayfile);

      FREE(rankfile);
      FREE(sarrayfile);
      FREE(saindex0file);
      FREE(sasamplesfile);
      FREE(sasampleqfile);
      FREE(csacompfiles[3]); FREE(csaptrfiles[3]);
      FREE(csacompfiles[2]); FREE(csaptrfiles[2]);
      FREE(csacompfiles[1]); FREE(csaptrfiles[1]);
      FREE(csacompfiles[0]); FREE(csaptrfiles[0]);
    }

  } else if (action == ARRAY_UNCOMPRESS) {
    if (argc <= 2) {
      start = end = 0;
    } else {
      start = strtoul(argv[1],NULL,10);
      end = strtoul(argv[2],NULL,10);
    }


  } else if (action == CHILD_UNCOMPRESS) {
    if (argc <= 2) {
      start = end = 0;
    } else {
      start = strtoul(argv[1],NULL,10);
      end = strtoul(argv[2],NULL,10);
    }

    chromosomefile = (char *) CALLOC(strlen(sourcedir)+strlen("/")+
				     strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
    sprintf(chromosomefile,"%s/%s.chromosome.iit",sourcedir,fileroot);
    if ((chromosome_iit = Univ_IIT_read(chromosomefile,/*readonlyp*/true,/*add_iit_p*/false)) == NULL) {
      fprintf(stderr,"IIT file %s is not valid\n",chromosomefile);
      exit(9);
    }
    FREE(chromosomefile);

    genomelength = Univ_IIT_genomelength(chromosome_iit,/*with_circular_alias_p*/true);
    Univ_IIT_free(&chromosome_iit);

    genomecomp = Genome_new(sourcedir,fileroot,/*snps_root*/NULL,/*genometype*/GENOME_OLIGOS,
			    /*uncompressedp*/false,/*access*/USE_MMAP_ONLY,/*sharedp*/false);

    sarrayfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("sarray")+1,sizeof(char));
    sprintf(sarrayfile,"%s/%s%ssarray",destdir,fileroot,mode_prefix);
    SA = (UINT4 *) Access_mmap(&sa_fd,&sa_len,sarrayfile,sizeof(UINT4),/*randomp*/true);
    FREE(sarrayfile);

    lcpchilddcfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("salcpchilddc")+1,sizeof(char));
    sprintf(lcpchilddcfile,"%s/%s%ssalcpchilddc",destdir,fileroot,mode_prefix);
    lcpchilddc = (unsigned char *) Access_mmap(&lcpchilddc_fd,&lcpchilddc_len,lcpchilddcfile,
					       sizeof(unsigned char),/*randomp*/true);
    FREE(lcpchilddcfile);

    lcpguidefile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("salcpguide1024")+1,sizeof(char));
    sprintf(lcpguidefile,"%s/%s%ssalcpguide1024",destdir,fileroot,mode_prefix);
    lcp_guide = (UINT4 *) Access_allocate(&shmid,&lcpguide_len,&seconds,lcpguidefile,sizeof(UINT4),/*sharedp*/false);
    FREE(lcpguidefile);

    lcpexcfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("salcpexc")+1,sizeof(char));
    sprintf(lcpexcfile,"%s/%s%ssalcpexc",destdir,fileroot,mode_prefix);
    lcp_exceptions = (UINT4 *) Access_allocate(&shmid,&lcpexc_len,&seconds,lcpexcfile,sizeof(UINT4),/*sharedp*/false);
    n_lcp_exceptions = lcpexc_len/(sizeof(UINT4) + sizeof(UINT4));
    FREE(lcpexcfile);

    childguidefile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("sachildguide1024")+1,sizeof(char));
    sprintf(childguidefile,"%s/%s%ssachildguide1024",destdir,fileroot,mode_prefix);
    child_guide = (UINT4 *) Access_allocate(&shmid,&childguide_len,&seconds,childguidefile,sizeof(UINT4),/*sharedp*/false);
    FREE(childguidefile);

    childexcfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("sachildexc")+1,sizeof(char));
    sprintf(childexcfile,"%s/%s%ssachildexc",destdir,fileroot,mode_prefix);
    child_exceptions = (UINT4 *) Access_allocate(&shmid,&childexc_len,&seconds,childexcfile,sizeof(UINT4),/*sharedp*/false);
    n_child_exceptions = childexc_len/(sizeof(UINT4) + sizeof(UINT4));
    FREE(childexcfile);

    Sarray_child_uncompress(genomecomp,lcpchilddc,lcp_guide,lcp_exceptions,n_lcp_exceptions,
			    child_guide,child_exceptions,n_child_exceptions,SA,genomelength,start,end);

    Genome_free(&genomecomp);

    munmap((void *) lcpchilddc,lcpchilddc_len);
    close(lcpchilddc_fd);

    munmap((void *) SA,sa_len);
    close(sa_fd);
  }

  return 0;
}

