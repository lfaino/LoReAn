static char rcsid[] = "$Id: sam_sort.c 155408 2014-12-16 07:00:44Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>		/* For off_t */
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bool.h"
#include "mem.h"
#include "access.h"
#include "complement.h"
#include "littleendian.h"
#include "genomicpos.h"
#include "chrnum.h"
#include "samheader.h"
#include "samread.h"
#include "samflags.h"
#include "stopwatch.h"
#include "datadir.h"
#include "filestring.h"
#include "getopt.h"


/************************************************************************
 *
 *  Check for correctness:
 *
 *    Run ./sam_sort -d <genome> --mark-dups --mark-first --dups-only --no-sam-headers on a SAM file
 *    Do "cut -f 1 | sort | uniq" to find all duplicate accessions
 *    
 *    Run the following Perl program on the original FASTA file
 *    Do "grep '^>' | sed -e 's/>//' | sort | uniq" to use as a gold standard
 * 
 * use IO::File;
 * $fastafile = $ARGV[0];
 * 
 * $FP = new IO::File($fastafile);
 * while (defined($header = <$FP>) && defined($queryseq5 = <$FP>) && defined($queryseq3 = <$FP>)) {
 *     chop $queryseq5; chop $queryseq3;
 * 
 *     if ($queryseq5 lt $queryseq3) {
 *       $nseen{$queryseq5 . $queryseq3} += 1;
 *     } else {
 *       $nseen{$queryseq3 . $queryseq5} += 1;
 *     }
 * }
 * close($FP);
 * 
 * $FP = new IO::File($fastafile);
 * while (defined($header = <$FP>) && defined($queryseq5 = <$FP>) && defined($queryseq3 = <$FP>)) {
 *     chop $queryseq5; chop $queryseq3;
 * 
 *     if ($queryseq5 lt $queryseq3) {
 * 	$queryseq = $queryseq5 . $queryseq3;
 *     } else {
 * 	$queryseq = $queryseq3 . $queryseq5;
 *     }
 * 
 *     if (($n = $nseen{$queryseq}) > 1) {
 *  	 print $header; print $queryseq5 . "\n"; print $queryseq3 . "\n";
 *     }
 * }
 * close($FP);
 * 
 * exit;
 * 
 ************************************************************************/


#define CHUNK 1024
typedef enum {NO_SECONDARY_SORT, ORIG_SECONDARY_SORT, ACC_SECONDARY_SORT, MATEFWD_SECONDARY_SORT, MATEREV_SECONDARY_SORT} Secondary_sort_T;

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Details of getting queryseqs */
#ifdef DEBUG9
#define debug9(x) x
#else
#define debug9(x)
#endif

/* Cell_binary_search */
#ifdef DEBUG10
#define debug10(x) x
#else
#define debug10(x)
#endif

/* Testing of linelen algorithm */
#ifdef DEBUG14
#define debug14(x) x
#else
#define debug14(x)
#endif



#ifdef HAVE_FSEEKO
#define moveto(fp,offset) fseeko(fp,offset,SEEK_SET)
#else
#define moveto(fp,offset) fseek(fp,offset,SEEK_SET)
#endif



/* Program Options */
static char *genomesubdir = NULL;
static char *dbroot = NULL;
static char *dbversion = NULL;
static char *user_genomedir = NULL;

static bool sam_headers_p = true;
static bool mark_duplicates_p = false;
static bool mark_first_p = false;
static bool print_unique_p = true;
static bool print_duplicates_p = true;
static bool secondary_sort_method = NO_SECONDARY_SORT;
static bool multiple_primaries_p = false;

static Stopwatch_T stopwatch = NULL;

static char *split_output_root = NULL;
static bool appendp = false;
static FILE **outputs = NULL;


static struct option long_options[] = {
  /* Input options */
  {"dir", required_argument, 0, 'D'}, /* user_genomedir */
  {"db", required_argument, 0, 'd'}, /* dbroot */
  {"split-output", required_argument, 0, 0}, /* outputs */
  {"append-output", no_argument, 0, 0},	     /* appendp */

  {"sort2", required_argument, 0, 0}, /* secondary_sort_method */

  {"mark-dups", no_argument, 0, 0}, /* mark_duplicates_p, print_unique_p, print_duplicates_p */
  {"mark-first", no_argument, 0, 0}, /* mark_first_p */
  {"dups-only", no_argument, 0, 0}, /* mark_duplicates_p, print_unique_p, print_duplicates_p */
  {"uniq-only", no_argument, 0, 0}, /* mark_duplicates_p, print_unique_p, print_duplicates_p */
  {"multiple-primaries", no_argument, 0, 0}, /* multiple_primaries_p */
  {"no-sam-headers", no_argument, 0, 0},     /* sam_headers_p */

  /* Help options */
  {"version", no_argument, 0, '^'}, /* print_program_version */
  {"help", no_argument, 0, '?'}, /* print_program_usage */
  {0, 0, 0, 0}
};


static void
print_program_version () {
  fprintf(stdout,"\n");
  fprintf(stdout,"sam_sort: sorting and duplication-marking utility for SAM files\n");
  fprintf(stdout,"Part of GMAP package, version %s\n",PACKAGE_VERSION);
  fprintf(stdout,"Thomas D. Wu, Genentech, Inc.\n");
  fprintf(stdout,"Contact: twu@gene.com\n");
  fprintf(stdout,"\n");
  return;
}

static void
print_program_usage () {
  fprintf(stdout,"\
Usage: sam_sort [OPTIONS...] -d genome <sam file>\n\
\n\
Input options\n\
  -D, --dir=STRING        Genome directory\n\
  -d, --db=STRING         Genome database.  If argument is '?' (with\n\
                            the quotes), this command lists available databases.\n\
Output file options\n\
  --split-output=STRING   Basename for multiple-file output, separately for nomapping,\n\
                            halfmapping_uniq, halfmapping_mult, unpaired_uniq, unpaired_mult,\n\
                            paired_uniq, paired_mult, concordant_uniq, and concordant_mult results\n\
  --append-output         When --split-output is given, this flag will append output\n\
                            to the existing files.  Otherwise, the default is to create new files.\n\
\n\
Other options\n\
  --sort2=STRING          For positions with the same genomic position, sort secondarily by\n\
                             none: no guarantee about the secondary sort (default)\n\
                             orig: original order in the SAM output file (what samtools sort does)\n\
                             accession: alphabetically by accession name\n\
                             mate-fwd: by genomic position of the mate, in ascending order\n\
                             mate-rev: by genomic position of the mate, in descending order\n\
                          For sorting by mate genomic position, a nomapping mate is treated as genomic position 0\n\
  --mark-dups             Mark duplicate reads by altering the flag accordingly\n\
  --mark-first            Also mark the first occurrence of a read that has a subsequent duplicate\n\
\n\
  --dups-only             Print only duplicate reads\n\
  --uniq-only             Print only unique reads\n\
  --multiple-primaries    Specify if GSNAP or GMAP was run with the --multiple-primaries flag\n\
  --no-sam-headers        Do not print SAM header lines\n\
");
  return;
}


static char complCode[128] = COMPLEMENT_LC;

static void
make_complement_inplace (char *sequence, unsigned int length) {
  char temp;
  unsigned int i, j;

  for (i = 0, j = length-1; i < length/2; i++, j--) {
    temp = complCode[(int) sequence[i]];
    sequence[i] = complCode[(int) sequence[j]];
    sequence[j] = temp;
  }
  if (i == j) {
    sequence[i] = complCode[(int) sequence[i]];
  }

  return;
}


#define T Cell_T
typedef struct T *T;
struct T {
  unsigned int flag;
  SAM_split_output_type split_output;
  bool low_read_p;

  Chrnum_T chrnum;
  Univcoord_T genomicpos;
  Univcoord_T genomicpos_extend_softclip;
  Univcoord_T mate_genomicpos;
  int filei;
  off_t linestart;
  int linelen;

  char *acc;			/* Needed for ACC_SECONDARY_SORT */

  int readindex;		/* inputi or outputi.  Needed for marking duplicates to find the other queryseq */
  char *queryseq5;
  char *queryseq3;

  char *queryseq_alpha1;	/* Pointers to queryseq5 and queryseq3.  Needed for handling 5/3 swaps. */
  char *queryseq_alpha2;

};


static void
Cell_standardize_queryseqs (T this) {

  if (this->queryseq5 == NULL && this->queryseq3 == NULL) {
    this->queryseq_alpha1 = this->queryseq_alpha2 = NULL;
  } else if (this->queryseq5 == NULL) {
    this->queryseq_alpha1 = this->queryseq3;
    this->queryseq_alpha2 = NULL;
  } else if (this->queryseq3 == NULL) {
    this->queryseq_alpha1 = this->queryseq5;
    this->queryseq_alpha2 = NULL;
  } else if (strcmp(this->queryseq5,this->queryseq3) < 0) {
    this->queryseq_alpha1 = this->queryseq5;
    this->queryseq_alpha2 = this->queryseq3;
  } else {
    this->queryseq_alpha1 = this->queryseq3;
    this->queryseq_alpha2 = this->queryseq5;
  }

  return;
}
    

/* initial_softclip needs to be determined only if we are marking duplicates */
static void
Cell_fill (struct T *this, int readindex, unsigned int flag, SAM_split_output_type split_output,
	   bool query_lowp, int initial_softclip, Univcoord_T genomicpos,
	   int filei, off_t fileposition, int linelen) {

  this->readindex = readindex;

  this->flag = flag;
  this->split_output = split_output;
  this->low_read_p = query_lowp;

  this->genomicpos = genomicpos;
  this->genomicpos_extend_softclip = this->genomicpos - initial_softclip;

  this->filei = filei;
  this->linestart = fileposition;
  this->linelen = linelen;

  this->queryseq5 = (char *) NULL;
  this->queryseq3 = (char *) NULL;

  return;
}

/* initial_softclip needs to be determined only if we are marking duplicates */
static void
Cell_fill_nodups (struct T *this, unsigned int flag, SAM_split_output_type split_output,
		  Univcoord_T genomicpos, int filei, off_t fileposition, int linelen) {

  this->readindex = 0;

  this->flag = flag;
  this->split_output = split_output;
  this->low_read_p = true;

  this->genomicpos = genomicpos;
  this->genomicpos_extend_softclip = genomicpos;

  this->filei = filei;
  this->linestart = fileposition;
  this->linelen = linelen;

  this->queryseq5 = (char *) NULL;
  this->queryseq3 = (char *) NULL;

  return;
}

#if 0
static void
print_fromfile (FILE *fp, off_t fileposition, int linelength) {
  char buffer[CHUNK];

  moveto(fp,fileposition);

  while (linelength > CHUNK) {
    fread(buffer,sizeof(char),CHUNK,fp);
    fwrite(buffer,sizeof(char),CHUNK,stdout);
    linelength -= CHUNK;
  }
  if (linelength > 0) {
    fread(buffer,sizeof(char),linelength,fp);
    fwrite(buffer,sizeof(char),linelength,stdout);
  }

  return;
}
#endif



static void
Cell_print_fromfile (FILE *fp_input, T this, Filestring_T headers) {
  char buffer[CHUNK];
  int linelength = this->linelen;
  FILE *fp_output;

#if 0
  if (nofailsp == true && this->split_output == OUTPUT_NM) {
    /* Skip */
    return;

  } else if (failsonlyp == true && this->split_output != OUTPUT_NM &&
	     this->split_output != OUTPUT_HX && this->split_output != OUTPUT_UX &&
	     this->split_output != OUTPUT_PX && this->split_output != OUTPUT_CX) {
    return;
  }

  if (failedinput_root != NULL && primaryp(this->flag) == true) {
    /* Convert SAM line to FASTA or FASTQ and write to a failedinput file */
  }
#endif

  if (split_output_root == NULL) {
    if ((fp_output = outputs[0]) == NULL) {
      fp_output = outputs[0] = stdout;
      Filestring_print(fp_output,headers);
    }

  } else {
    if ((fp_output = outputs[this->split_output]) == NULL) {
      fp_output = outputs[this->split_output] = SAM_header_open_file(this->split_output,split_output_root,/*appendp*/false);
      Filestring_print(fp_output,headers);
    }
  }

  moveto(fp_input,this->linestart);

#ifdef DEBUG
  printf("readindex %d: ",this->readindex);
#endif

  while (linelength > CHUNK) {
    fread(buffer,sizeof(char),CHUNK,fp_input);
    fwrite(buffer,sizeof(char),CHUNK,fp_output);
    linelength -= CHUNK;
  }
  if (linelength > 0) {
    fread(buffer,sizeof(char),linelength,fp_input);
    fwrite(buffer,sizeof(char),linelength,fp_output);
  }

  return;
}



static int
Cell_genomicpos_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->genomicpos != 0 && y->genomicpos == 0) {
    return -1;
  } else if (y->genomicpos != 0 && x->genomicpos == 0) {
    return +1;
  } else if (x->genomicpos < y->genomicpos) {
    return -1;
  } else if (y->genomicpos < x->genomicpos) {
    return +1;
  } else {
    return 0;
  }
}

static int
Cell_genomicpos_linestart_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->genomicpos != 0 && y->genomicpos == 0) {
    return -1;
  } else if (y->genomicpos != 0 && x->genomicpos == 0) {
    return +1;
  } else if (x->genomicpos < y->genomicpos) {
    return -1;
  } else if (y->genomicpos < x->genomicpos) {
    return +1;
  } else if (x->linestart < y->linestart) {
    return -1;
  } else if (y->linestart < x->linestart) {
    return +1;
  } else {
    return 0;
  }
}

static int
Cell_genomicpos_extend_softclip_lowhigh_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->genomicpos_extend_softclip != 0 && y->genomicpos_extend_softclip == 0) {
    return -1;
  } else if (y->genomicpos_extend_softclip != 0 && x->genomicpos_extend_softclip == 0) {
    return +1;
  } else if (x->genomicpos_extend_softclip < y->genomicpos_extend_softclip) {
    return -1;
  } else if (y->genomicpos_extend_softclip < x->genomicpos_extend_softclip) {
    return +1;

  } else if (x->low_read_p == true && y->low_read_p == false) {
    return -1;
  } else if (y->low_read_p == true && x->low_read_p == false) {
    return +1;
  } else {
    return 0;
  }
}

static int
Cell_accession_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  return strcmp(x->acc,y->acc);
}

static int
Cell_matefwd_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->mate_genomicpos < y->mate_genomicpos) {
    return -1;
  } else if (y->mate_genomicpos < x->mate_genomicpos) {
    return +1;
  } else {
    return 0;
  }
}

static int
Cell_materev_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->mate_genomicpos > y->mate_genomicpos) {
    return -1;
  } else if (y->mate_genomicpos > x->mate_genomicpos) {
    return +1;
  } else {
    return 0;
  }
}


static int
Cell_queryseq_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;
  int cmp;

  if (x->queryseq_alpha1 != NULL && y->queryseq_alpha1 == NULL) {
    return -1;
  } else if (y->queryseq_alpha1 == NULL && x->queryseq_alpha1 != NULL) {
    return +1;
  } else if (x->queryseq_alpha1 == NULL && y->queryseq_alpha1 == NULL) {
    return 0;
  } else if ((cmp = strcmp(x->queryseq_alpha1,y->queryseq_alpha1)) != 0) {
    return cmp;
  } else if (x->queryseq_alpha2 != NULL && y->queryseq_alpha2 == NULL) {
    return -1;
  } else if (y->queryseq_alpha2 == NULL && x->queryseq_alpha2 != NULL) {
    return +1;
  } else if (x->queryseq_alpha2 == NULL && y->queryseq_alpha2 == NULL) {
    return 0;
  } else {
    return strcmp(x->queryseq_alpha2,y->queryseq_alpha2);
  }
}


#if 0
static int
Cell_binary_search (int lowi, int highi, T *cells, Univcoord_T goal) {
  int middlei;

  debug10(printf("entered binary search with lowi=%d, highi=%d, goal=%u\n",lowi,highi,goal));

  while (lowi < highi) {
    middlei = lowi + ((highi - lowi) / 2);
    debug10(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		   lowi,cells[lowi]->genomicpos,middlei,cells[middlei]->genomicpos,
		   highi,cells[highi]->genomicpos,goal));
    if (goal < cells[middlei]->genomicpos) {
      highi = middlei;
    } else if (goal > cells[middlei]->genomicpos) {
      lowi = middlei + 1;
    } else {
      debug10(printf("binary search returns %d\n",middlei));
      /* Rewind to first cell having the goal */
      while (middlei - 1 >= lowi && cells[middlei - 1]->genomicpos == goal) {
	middlei--;
      }
      return middlei;
    }
  }

  debug10(printf("binary search returns %d\n",highi));
  return highi;
}
#endif

#if 0
static int
Cell_find (int lowi, int highi, T *cells, Univcoord_T goal, int readindex) {
  int i;

  i = Cell_binary_search(lowi,highi,cells,goal);
  while (i < highi && cells[i]->genomicpos == goal) {
    if (cells[i]->readindex == readindex) {
      return i;
    } else {
      i++;
    }
  }

  fprintf(stderr,"Cannot find cell in %d..%d with genomicpos %u and readindex %d\n",
	  lowi,highi,goal,readindex);
  return -1;
}
#endif


static void
process_without_dups (FILE **sam_inputs, int *headerlengths, int *ncells, int ninputs,
		      Intlist_T linelengths, int ncells_total, Univ_IIT_T chromosome_iit,
		      Univcoord_T *chroffsets, Filestring_T headers) {
  T *cells;
  FILE *fp_sam;
  int filei, linei;
  int n_mappers = 0, n_nomappers = 0;
  Intlist_T l;
  struct T *cells_allocated, *ptr;
  int i, j, k;

  off_t fileposition;
  int linelen;
  unsigned int flag;
  SAM_split_output_type split_output;
  Univcoord_T genomicpos;

  int acclength;

  ptr = cells_allocated = (struct T *) MALLOC(ncells_total * sizeof(struct T));
  cells = (T *) MALLOC(ncells_total * sizeof(T));
  for (i = 0; i < ncells_total; i++) {
    cells[i] = &(ptr[i]);
  }

  fprintf(stderr,"Reading SAM files...\n");

  k = 0;
  l = linelengths;
  for (filei = 0; filei < ninputs; filei++) {
    fprintf(stderr,"  Reading file %d...",filei+1);
    fp_sam = sam_inputs[filei];
    fileposition = headerlengths[filei];
    for (linei = 0; linei < ncells[filei]; linei++) {
      linelen = Intlist_head(l);
      moveto(fp_sam,fileposition);
      genomicpos = Samread_parse_genomicpos_fromfile(fp_sam,&flag,&split_output,
						     chromosome_iit,chroffsets,linelen);
      Cell_fill_nodups(cells[k++],flag,split_output,genomicpos,filei,fileposition,linelen);
      if (flag & QUERY_UNMAPPED) {
	n_nomappers++;
      } else {
	n_mappers++;
      }
      fileposition += linelen;
      l = Intlist_next(l);
    }

    fprintf(stderr,"done\n");
  }

  /* Sort and print */
  if (secondary_sort_method == NO_SECONDARY_SORT) {
    Stopwatch_start(stopwatch);
    fprintf(stderr,"Sorting entries by genomicpos...");
    qsort(cells,ncells_total,sizeof(T),Cell_genomicpos_cmp);
    fprintf(stderr,"done (%.1f seconds)\n",Stopwatch_stop(stopwatch));

    Stopwatch_start(stopwatch);
    fprintf(stderr,"Printing entries...");
    for (k = 0; k < ncells_total; k++) {
      debug(printf("%u\t%u\t%d\n",cells[k]->genomicpos,cells[k]->linestart,cells[k]->linelen));
      fp_sam = sam_inputs[cells[k]->filei];
      Cell_print_fromfile(fp_sam,cells[k],headers);
    }
    fprintf(stderr,"done (%.1f seconds)\n",Stopwatch_stop(stopwatch));

  } else if (secondary_sort_method == ORIG_SECONDARY_SORT) {
    Stopwatch_start(stopwatch);
    fprintf(stderr,"Sorting entries by genomicpos and original file position...");
    qsort(cells,ncells_total,sizeof(T),Cell_genomicpos_linestart_cmp);
    fprintf(stderr,"done (%.1f seconds)\n",Stopwatch_stop(stopwatch));

    Stopwatch_start(stopwatch);
    fprintf(stderr,"Printing entries...");
    for (k = 0; k < ncells_total; k++) {
      fp_sam = sam_inputs[cells[k]->filei];
      Cell_print_fromfile(fp_sam,cells[k],headers);
    }
    fprintf(stderr,"done (%.1f seconds)\n",Stopwatch_stop(stopwatch));

  } else if (secondary_sort_method == ACC_SECONDARY_SORT) {
    Stopwatch_start(stopwatch);
    fprintf(stderr,"Sorting entries by genomicpos...");
    qsort(cells,ncells_total,sizeof(T),Cell_genomicpos_cmp);
    fprintf(stderr,"done (%.1f seconds)\n",Stopwatch_stop(stopwatch));

    Stopwatch_start(stopwatch);
    fprintf(stderr,"Subsorting by accession and printing entries...");
    i = 0;
    while (i < n_mappers) {
      j = i + 1;
      while (j < n_mappers && cells[j]->genomicpos == cells[i]->genomicpos) {
	j++;
      }
      
      if (j > i + 1) {
	for (k = i; k < j; k++) {
	  fp_sam = sam_inputs[cells[k]->filei];
	  moveto(fp_sam,cells[k]->linestart);
	  cells[k]->acc = Samread_get_acc_fromfile(&acclength,fp_sam,cells[k]->linelen);
	}
	
	qsort(&(cells[i]),j - i,sizeof(T),Cell_accession_cmp);
	for (k = i; k < j; k++) {
	  FREE(cells[k]->acc);
	}
      }
      
      for (k = i; k < j; k++) {
	fp_sam = sam_inputs[cells[k]->filei];
	Cell_print_fromfile(fp_sam,cells[k],headers);
      }
      
      i = j;
    }
    
    if (ncells_total > n_mappers + 1) {
      for (k = n_mappers; k < ncells_total; k++) {
	fp_sam = sam_inputs[cells[k]->filei];
	moveto(fp_sam,cells[k]->linestart);
	cells[k]->acc = Samread_get_acc_fromfile(&acclength,fp_sam,cells[k]->linelen);
      }
	
      qsort(&(cells[n_mappers]),n_nomappers,sizeof(T),Cell_accession_cmp);
      for (k = n_mappers; k < ncells_total; k++) {
	FREE(cells[k]->acc);
      }
    }
      
    for (k = n_mappers; k < ncells_total; k++) {
      fp_sam = sam_inputs[cells[k]->filei];
      Cell_print_fromfile(fp_sam,cells[k],headers);
    }
    fprintf(stderr,"done (%.1f seconds)\n",Stopwatch_stop(stopwatch));

  } else if (secondary_sort_method == MATEFWD_SECONDARY_SORT || secondary_sort_method == MATEREV_SECONDARY_SORT) {
    Stopwatch_start(stopwatch);
    fprintf(stderr,"Sorting entries by genomicpos...");
    qsort(cells,ncells_total,sizeof(T),Cell_genomicpos_cmp);
    fprintf(stderr,"done (%.1f seconds)\n",Stopwatch_stop(stopwatch));

    Stopwatch_start(stopwatch);
    fprintf(stderr,"Subsorting by mate position and printing entries...");
    i = 0;
    while (i < n_mappers) {
      j = i + 1;
      while (j < n_mappers && cells[j]->genomicpos == cells[i]->genomicpos) {
	j++;
      }
      
      if (j > i + 1) {
	for (k = i; k < j; k++) {
	  moveto(fp_sam,cells[k]->linestart);
	  cells[k]->mate_genomicpos = Samread_parse_mate_genomicpos_fromfile(fp_sam,chromosome_iit,chroffsets,cells[k]->linelen);
	}
	
	if (secondary_sort_method == MATEFWD_SECONDARY_SORT) {
	  qsort(&(cells[i]),j - i,sizeof(T),Cell_matefwd_cmp);
	} else {
	  qsort(&(cells[i]),j - i,sizeof(T),Cell_materev_cmp);
	}
      }
      
      for (k = i; k < j; k++) {
	fp_sam = sam_inputs[cells[k]->filei];
	Cell_print_fromfile(fp_sam,cells[k],headers);
      }

      i = j;
    }
      
    if (ncells_total > n_mappers + 1) {
      for (k = n_mappers; k < ncells_total; k++) {
	moveto(fp_sam,cells[k]->linestart);
	cells[k]->mate_genomicpos = Samread_parse_mate_genomicpos_fromfile(fp_sam,chromosome_iit,chroffsets,cells[k]->linelen);
      }
      
      if (secondary_sort_method == MATEFWD_SECONDARY_SORT) {
	qsort(&(cells[n_mappers]),n_nomappers,sizeof(T),Cell_matefwd_cmp);
      } else {
	qsort(&(cells[n_mappers]),n_nomappers,sizeof(T),Cell_materev_cmp);
      }
    }
      
    for (k = n_mappers; k < ncells_total; k++) {
      fp_sam = sam_inputs[cells[k]->filei];
      Cell_print_fromfile(fp_sam,cells[k],headers);
    }
    fprintf(stderr,"done (%.1f seconds)\n",Stopwatch_stop(stopwatch));

  } else {
    fprintf(stderr,"Secondary sort method not recognized\n");
    abort();
  }

  FREE(cells);
  FREE(cells_allocated);

  return;
}


static int
process_with_dups (FILE **sam_inputs, int *headerlengths, int *ncells, int ninputs,
		   Intlist_T linelengths, int ncells_total, Univ_IIT_T chromosome_iit,
		   Univcoord_T *chroffsets, Filestring_T headers) {
  FILE *fp_sam;
  int filei, linei;
  int nmarked = 0;
  int n_mappers = 0, n_nomappers = 0;
  T *cells, mate;
  struct T *cells_allocated, *ptr;
  int *queryseq5_index, *queryseq3_index;
  int i, k;
  int j, j_low, j_high, mate_allocated;

  off_t fileposition;
  int linelen;
  Univcoord_T genomicpos;
  char *hiti;

  Intlist_T l;
  unsigned int flag;
  SAM_split_output_type split_output;
  int initial_softclip;
  char *acc, *last_acc, *read;
  int readindex, nreads;
  int acclength, last_acclength, readlength;
  bool query_lowp;
  bool *duplicatep = NULL;
  bool all_duplicates_p;



  /* Actually, array lengths should be nreads, but we don't know that yet */
  queryseq5_index = (int *) CALLOC(ncells_total,sizeof(int));
  queryseq3_index = (int *) CALLOC(ncells_total,sizeof(int));

  ptr = cells_allocated = (struct T *) MALLOC(ncells_total * sizeof(struct T));
  cells = (T *) MALLOC(ncells_total * sizeof(T));
  for (i = 0; i < ncells_total; i++) {
    cells[i] = &(ptr[i]);
  }
    

  last_acc = MALLOC(sizeof(char));
  last_acc[0] = '\0';
  last_acclength = 0;
  readindex = -1;		/* readindex is 0-based */

  fprintf(stderr,"Reading SAM files...\n");

  k = 0;
  l = linelengths;
  for (filei = 0; filei < ninputs; filei++) {
    fprintf(stderr,"  Reading file %d...",filei+1);
    fp_sam = sam_inputs[filei];
    fileposition = headerlengths[filei];
    for (linei = 0; linei < ncells[filei]; linei++) {
      linelen = Intlist_head(l);
      moveto(fp_sam,fileposition);
      acc = Samread_parse_acc_and_softclip_fromfile(&acclength,&flag,&split_output,&hiti,
						    &genomicpos,&initial_softclip,&query_lowp,
						    fp_sam,chromosome_iit,chroffsets,linelen);
      if (acclength != last_acclength) {
	readindex++;
      } else if (strcmp(acc,last_acc)) {
	readindex++;
      }
      FREE(last_acc);
      last_acc = acc;
      last_acclength = acclength;

      if (flag & QUERY_UNMAPPED) {
	n_nomappers++;
      } else {
	n_mappers++;
      }

      /* debug(printf("Read readindex %d, chrnum %d, chrpos %u, linelen %d\n",readindex,chrnum,chrpos,linelen)); */
      if (flag & NOT_PRIMARY) {
      /* Don't use secondary hit for accessing reads */

      } else if (multiple_primaries_p == true) {
#if 0
	/* Now always parsed */
	hiti = Samread_parse_aux_fromfile(fp_sam,/*auxfield*/"HI",linelen);
#endif
	if (strcmp(hiti,"1")) {
	  /* Don't use second or later primary hit for accessing reads */
	} else if (flag & FIRST_READ_P) {
	  queryseq5_index[readindex] = k;
	} else {
	  queryseq3_index[readindex] = k;
	}
	
      } else {
	if (flag & FIRST_READ_P) {
	  queryseq5_index[readindex] = k;
	} else {
	  queryseq3_index[readindex] = k;
	}
      }
      
      FREE(hiti);
      Cell_fill(cells[k++],readindex,flag,split_output,query_lowp,initial_softclip,genomicpos,filei,fileposition,linelen);

      fileposition += linelen;
      l = Intlist_next(l);
    }

    fprintf(stderr,"done\n");
  }
  FREE(last_acc);

  nreads = readindex + 1;
  duplicatep = (bool *) CALLOC(nreads,sizeof(bool));
  
  /* Sort entries, based on genomicpos_extend_softclip */
  Stopwatch_start(stopwatch);
  fprintf(stderr,"Sorting SAM lines...");
  qsort(cells,ncells_total,sizeof(T),Cell_genomicpos_extend_softclip_lowhigh_cmp);
  fprintf(stderr,"done (%.1f seconds)\n",Stopwatch_stop(stopwatch));

  /* Mark all duplicates within mappers, based on genomicpos_extend_softclip */
  Stopwatch_start(stopwatch);
  fprintf(stderr,"Finding duplicates...");

  i = 0;
  while (i < n_mappers) {
    j_low = i + 1;
    while (j_low < n_mappers && cells[j_low]->genomicpos_extend_softclip == cells[i]->genomicpos_extend_softclip && 
	   cells[j_low]->low_read_p == true) {
      j_low++;
    }

    j_high = j_low;
    while (j_high < n_mappers && cells[j_high]->genomicpos_extend_softclip == cells[i]->genomicpos_extend_softclip) {
      j_high++;
    }

    if (j_low > i + 1) {
      debug(printf("\nFound multiple low hits from %d to (%d - 1) at genomicpos_extend_softclip %u\n",
		   i,j_low,cells[i]->genomicpos_extend_softclip));
	  
      /* Multiple low hits with same chrpos, so opportunity to mark duplicatep */
      /* Find queryseqs for each */
      for (k = i; k < j_low; k++) {
	fp_sam = sam_inputs[cells[k]->filei];
	debug9(printf("Looking for queryseqs for "));
	debug9(Cell_print_fromfile(fp_sam,cells[k],headers));

	if (cells[k]->flag & FIRST_READ_P) {
	  debug9(printf("Flag for entry %d is %u, indicating a first read\n",k,cells[k]->flag));
	  moveto(fp_sam,cells[k]->linestart);
	  Samread_parse_read_fromfile(fp_sam,&flag,&readlength,&read,cells[k]->linelen);
	  if (cells[k]->flag & QUERY_MINUSP) {
	    debug(printf("complementing queryseq5\n"));
	    make_complement_inplace(read,readlength);
	  }
	  cells[k]->queryseq5 = read;
	  debug9(printf("queryseq5 is %s\n",read));
	    
	  mate_allocated = queryseq3_index[cells[k]->readindex];
	  mate = &(cells_allocated[mate_allocated]);

	  debug9(printf("Mate is "));
	  debug9(Cell_print_fromfile(fp_sam,mate,headers));
	  moveto(fp_sam,mate->linestart);
	  Samread_parse_read_fromfile(fp_sam,&flag,&readlength,&read,mate->linelen);
	  if (mate->flag & QUERY_MINUSP) {
	    debug(printf("complementing queryseq3\n"));
	    make_complement_inplace(read,readlength);
	  }
	  cells[k]->queryseq3 = read;
	  debug9(printf("queryseq3 is %s\n",read));

	} else {
	  debug9(printf("Flag for entry %d is %u, indicating a second read\n",k,cells[k]->flag));
	  moveto(fp_sam,cells[k]->linestart);
	  Samread_parse_read_fromfile(fp_sam,&flag,&readlength,&read,cells[k]->linelen);
	  if (cells[k]->flag & QUERY_MINUSP) {
	    debug(printf("complementing queryseq3\n"));
	    make_complement_inplace(read,readlength);
	  }
	  cells[k]->queryseq3 = read;
	  debug9(printf("queryseq3 is %s\n",read));

	  mate_allocated = queryseq5_index[cells[k]->readindex];
	  mate = &(cells_allocated[mate_allocated]);
	  debug9(printf("Mate is "));
	  debug9(Cell_print_fromfile(fp_sam,mate,headers));
	  moveto(fp_sam,mate->linestart);
	  Samread_parse_read_fromfile(fp_sam,&flag,&readlength,&read,mate->linelen);
	  if (mate->flag & QUERY_MINUSP) {
	    debug(printf("complementing queryseq5\n"));
	    make_complement_inplace(read,readlength);
	  }
	  cells[k]->queryseq5 = read;
	  debug9(printf("queryseq5 is %s\n",read));
	}
	  
	Cell_standardize_queryseqs(cells[k]);
      }

      qsort(&(cells[i]),j_low - i,sizeof(T),Cell_queryseq_cmp);

      for (k = i + 1; k < j_low; k++) {
	debug(printf("Comparing cell %d with %d => cmp %d\n",k,k-1,Cell_queryseq_cmp(&(cells[k]),&(cells[k-1]))));
	debug(printf("  %s %s\n",cells[k-1]->queryseq_alpha1,cells[k-1]->queryseq_alpha2));
	debug(printf("  %s %s\n",cells[k]->queryseq_alpha1,cells[k]->queryseq_alpha2));
	if (Cell_queryseq_cmp(&(cells[k]),&(cells[k-1])) == 0) {
	  readindex = cells[k]->readindex;
	  duplicatep[readindex] = true;
	  if (mark_first_p == true) {
	    readindex = cells[k-1]->readindex;
	    duplicatep[readindex] = true;
	  }
	}
      }

      for (k = i; k < j_low; k++) {
	FREE(cells[k]->queryseq5);
	FREE(cells[k]->queryseq3);
      }
    }


    /* Also analyze high ends, to avoid a false negative when the initial_softclip extension is wrong
       due to an end indel */
    if (j_high > j_low + 1) {
      all_duplicates_p = true;
      for (k = j_low; k < j_high; k++) {
	readindex = cells[k]->readindex;
	if (duplicatep[readindex] == false) {
	  all_duplicates_p = false;
	}
      }

      if (all_duplicates_p == false) {
	debug(printf("\nFound multiple high hits from %d to (%d - 1) at genomicpos_extend_softclip %u\n",
		     j_low,j_high,cells[j_low]->genomicpos_extend_softclip));
	
	/* Multiple high hits with same chrpos, so opportunity to mark duplicatep */
        /* Find queryseqs for each */
	for (k = j_low; k < j_high; k++) {
	  fp_sam = sam_inputs[cells[k]->filei];
	  debug9(printf("Looking for queryseqs for "));
	  debug9(Cell_print_fromfile(fp_sam,cells[k],headers));

	  if (cells[k]->flag & FIRST_READ_P) {
	    debug9(printf("Flag for entry %d is %u, indicating a first read\n",k,cells[k]->flag));
	    moveto(fp_sam,cells[k]->linestart);
	    Samread_parse_read_fromfile(fp_sam,&flag,&readlength,&read,cells[k]->linelen);
	    if (cells[k]->flag & QUERY_MINUSP) {
	      debug(printf("complementing queryseq5\n"));
	      make_complement_inplace(read,readlength);
	    }
	    cells[k]->queryseq5 = read;
	    debug9(printf("queryseq5 is %s\n",read));
	    
	    mate_allocated = queryseq3_index[cells[k]->readindex];
	    mate = &(cells_allocated[mate_allocated]);
	    debug9(printf("Mate is "));
	    debug9(Cell_print_fromfile(fp_sam,mate,headers));
	    moveto(fp_sam,mate->linestart);
	    Samread_parse_read_fromfile(fp_sam,&flag,&readlength,&read,mate->linelen);
	    if (mate->flag & QUERY_MINUSP) {
	      debug(printf("complementing queryseq3\n"));
	      make_complement_inplace(read,readlength);
	    }
	    cells[k]->queryseq3 = read;
	    debug9(printf("queryseq3 is %s\n",read));

	  } else {
	    debug9(printf("Flag for entry %d is %u, indicating a second read\n",k,cells[k]->flag));
	    moveto(fp_sam,cells[k]->linestart);
	    Samread_parse_read_fromfile(fp_sam,&flag,&readlength,&read,cells[k]->linelen);
	    if (cells[k]->flag & QUERY_MINUSP) {
	      debug(printf("complementing queryseq3\n"));
	      make_complement_inplace(read,readlength);
	    }
	    cells[k]->queryseq3 = read;
	    debug9(printf("queryseq3 is %s\n",read));

	    mate_allocated = queryseq5_index[cells[k]->readindex];
	    mate = &(cells_allocated[mate_allocated]);
	    debug9(printf("Mate is "));
	    debug9(Cell_print_fromfile(fp_sam,mate,headers));
	    moveto(fp_sam,mate->linestart);
	    Samread_parse_read_fromfile(fp_sam,&flag,&readlength,&read,mate->linelen);
	    if (mate->flag & QUERY_MINUSP) {
	      debug(printf("complementing queryseq5\n"));
	      make_complement_inplace(read,readlength);
	    }
	    cells[k]->queryseq5 = read;
	    debug9(printf("queryseq5 is %s\n",read));
	  }
	  
	  Cell_standardize_queryseqs(cells[k]);
	}

	qsort(&(cells[j_low]),j_high - j_low,sizeof(T),Cell_queryseq_cmp);

	for (k = j_low + 1; k < j_high; k++) {
	  debug(printf("Comparing cell %d with %d => cmp %d\n",k,k-1,Cell_queryseq_cmp(&(cells[k]),&(cells[k-1]))));
	  debug(printf("  %s %s\n",cells[k-1]->queryseq_alpha1,cells[k-1]->queryseq_alpha2));
	  debug(printf("  %s %s\n",cells[k]->queryseq_alpha1,cells[k]->queryseq_alpha2));
	  if (Cell_queryseq_cmp(&(cells[k]),&(cells[k-1])) == 0) {
	    readindex = cells[k]->readindex;
	    duplicatep[readindex] = true;
	    if (mark_first_p == true) {
	      readindex = cells[k-1]->readindex;
	      duplicatep[readindex] = true;
	    }
	  }
	}

	for (k = j_low; k < j_high; k++) {
	  FREE(cells[k]->queryseq5);
	  FREE(cells[k]->queryseq3);
	}
      }
    }

    i = j_high;
  }

  /* Mark all duplicates within nomappers, based on queryseq */
  for (k = n_mappers; k < ncells_total; k++) {
    if (duplicatep[cells[k]->readindex] == true) {
      cells[k]->queryseq5 = cells[k]->queryseq3 = NULL; /* Will be sorted to end of list */

    } else if (cells[k]->flag & FIRST_READ_P) {
      fp_sam = sam_inputs[cells[k]->filei];
      moveto(fp_sam,cells[k]->linestart);
      Samread_parse_read_fromfile(fp_sam,&flag,&readlength,&read,cells[k]->linelen);
      if (cells[k]->flag & QUERY_MINUSP) {
	make_complement_inplace(read,readlength);
      }
      cells[k]->queryseq5 = read;
	    
      mate_allocated = queryseq3_index[cells[k]->readindex];
      mate = &(cells_allocated[mate_allocated]);
      moveto(fp_sam,mate->linestart);
      Samread_parse_read_fromfile(fp_sam,&flag,&readlength,&read,mate->linelen);
      if (mate->flag & QUERY_MINUSP) {
	make_complement_inplace(read,readlength);
      }
      cells[k]->queryseq3 = read;

      Cell_standardize_queryseqs(cells[k]);

    } else {
      fp_sam = sam_inputs[cells[k]->filei];
      moveto(fp_sam,cells[k]->linestart);
      Samread_parse_read_fromfile(fp_sam,&flag,&readlength,&read,cells[k]->linelen);
      if (cells[k]->flag & QUERY_MINUSP) {
	make_complement_inplace(read,readlength);
      }
      cells[k]->queryseq3 = read;
      
      mate_allocated = queryseq5_index[cells[k]->readindex];
      mate = &(cells_allocated[mate_allocated]);
      moveto(fp_sam,mate->linestart);
      Samread_parse_read_fromfile(fp_sam,&flag,&readlength,&read,mate->linelen);
      if (mate->flag & QUERY_MINUSP) {
	make_complement_inplace(read,readlength);
      }
      cells[k]->queryseq5 = read;
	  
      Cell_standardize_queryseqs(cells[k]);
    }
  }

  FREE(queryseq3_index);
  FREE(queryseq5_index);


  /* Sort non-mapping entries based on queryseqs */
  qsort(&(cells[n_mappers]),n_nomappers,sizeof(T),Cell_queryseq_cmp);

  for (k = n_mappers + 1; k < ncells_total; k++) {
    debug(printf("Comparing cell %d with %d => cmp %d\n",k,k-1,Cell_queryseq_cmp(&(cells[k]),&(cells[k-1]))));
    if (Cell_queryseq_cmp(&(cells[k]),&(cells[k-1])) == 0) {
      readindex = cells[k]->readindex;
      duplicatep[readindex] = true;
      if (mark_first_p == true) {
	readindex = cells[k-1]->readindex;
	duplicatep[readindex] = true;
      }
    }
  }

  for (k = n_mappers; k < ncells_total; k++) {
    FREE(cells[k]->queryseq5);
    FREE(cells[k]->queryseq3);
  }
  fprintf(stderr,"done (%.1f seconds)\n",Stopwatch_stop(stopwatch));


  /* Re-sort entries, based on genomicpos (not extended by initial_softclip), and secondary criterion */
  Stopwatch_start(stopwatch);
  fprintf(stderr,"Re-sorting entries...");
  if (secondary_sort_method == NO_SECONDARY_SORT) {
    qsort(cells,ncells_total,sizeof(T),Cell_genomicpos_cmp);

  } else if (secondary_sort_method == ORIG_SECONDARY_SORT) {
    qsort(cells,ncells_total,sizeof(T),Cell_genomicpos_linestart_cmp);

  } else if (secondary_sort_method == ACC_SECONDARY_SORT) {
    qsort(cells,ncells_total,sizeof(T),Cell_genomicpos_cmp);

    i = 0;
    while (i < ncells_total) {
      j = i + 1;
      while (j < ncells_total && cells[j]->genomicpos == cells[i]->genomicpos) {
	j++;
      }

      if (j > i + 1) {
	for (k = i; k < j; k++) {
	  fp_sam = sam_inputs[cells[k]->filei];
	  moveto(fp_sam,cells[k]->linestart);
	  cells[k]->acc = Samread_get_acc_fromfile(&acclength,fp_sam,cells[k]->linelen);
	}
	qsort(&(cells[i]),j - i,sizeof(T),Cell_accession_cmp);
	for (k = i; k < j; k++) {
	  FREE(cells[k]->acc);
	}
      }

      i = j;
    }

  } else if (secondary_sort_method == MATEFWD_SECONDARY_SORT ||
	     secondary_sort_method == MATEREV_SECONDARY_SORT) {
    qsort(cells,ncells_total,sizeof(T),Cell_genomicpos_cmp);

    i = 0;
    while (i < ncells_total) {
      j = i + 1;
      while (j < ncells_total && cells[j]->genomicpos == cells[i]->genomicpos) {
	j++;
      }

      if (j > i + 1) {
	for (k = i; k < j; k++) {
	  fp_sam = sam_inputs[cells[k]->filei];
	  moveto(fp_sam,cells[k]->linestart);
	  cells[k]->mate_genomicpos = Samread_parse_mate_genomicpos_fromfile(fp_sam,chromosome_iit,chroffsets,cells[k]->linelen);
	}
	if (secondary_sort_method == MATEFWD_SECONDARY_SORT) {
	  qsort(&(cells[i]),j - i,sizeof(T),Cell_matefwd_cmp);
	} else {
	  qsort(&(cells[i]),j - i,sizeof(T),Cell_materev_cmp);
	}
      }

      i = j;
    }
  }
  fprintf(stderr,"done (%.1f seconds)\n",Stopwatch_stop(stopwatch));

  /* Print results */
  Stopwatch_start(stopwatch);
  fprintf(stderr,"Printing results...");

  for (k = 0; k < ncells_total; k++) {
    if (duplicatep[cells[k]->readindex] == true) {
      if (print_duplicates_p == true) {
	fp_sam = sam_inputs[cells[k]->filei];
	moveto(fp_sam,cells[k]->linestart);
	Samread_print_as_duplicate_fromfile(fp_sam,cells[k]->linelen);
      }
      nmarked++;
    } else {
      if (print_unique_p == true) {
	/* Non-duplicate */
	fp_sam = sam_inputs[cells[k]->filei];
	Cell_print_fromfile(fp_sam,cells[k],headers);
      }
    }
  }
  fprintf(stderr,"done (%.1f seconds)\n",Stopwatch_stop(stopwatch));

  FREE(duplicatep);

  FREE(cells);
  FREE(cells_allocated);

  return nmarked;
}



#define BUFFERLEN 1024

int
main (int argc, char *argv[]) {
  FILE **sam_inputs, *fp_sam;
  int ninputs, filei;
  int nchromosomes, i;
  Univcoord_T *chroffsets;
  Chrpos_T *chrlengths;
  off_t fileposition;
  int lastchar;

  char buffer[BUFFERLEN], *lastp, *p;
  Intlist_T linelengths;
  int *headerlengths, linelen;
  Filestring_T headers = NULL;
#ifdef DEBUG14
  Intlist_T linelengths_goldstd;
  int linelen_goldstd;
#endif

  char *fileroot = NULL, *iitfile;
  Univ_IIT_T chromosome_iit;
  int *ncells, ncells_total, nmarked;

  int opt;
  extern int optind;
  extern char *optarg;
  int long_option_index = 0;
  const char *long_name;

  while ((opt = getopt_long(argc,argv,"D:d:^?",
			    long_options,&long_option_index)) != -1) {
    switch (opt) {
    case 0:
      long_name = long_options[long_option_index].name;
      if (!strcmp(long_name,"version")) {
	print_program_version();
	exit(0);
      } else if (!strcmp(long_name,"help")) {
	print_program_usage();
	exit(0);

      } else if (!strcmp(long_name,"split-output")) {
	split_output_root = optarg;
      } else if (!strcmp(long_name,"append-output")) {
	appendp = true;

      } else if (!strcmp(long_name,"sort2")) {
	if (!strcmp(optarg,"none")) {
	  secondary_sort_method = NO_SECONDARY_SORT;
	} else if (!strcmp(optarg,"orig")) {
	  secondary_sort_method = ORIG_SECONDARY_SORT;
	} else if (!strcmp(optarg,"accession")) {
	  secondary_sort_method = ACC_SECONDARY_SORT;
	} else if (!strcmp(optarg,"mate-fwd")) {
	  secondary_sort_method = MATEFWD_SECONDARY_SORT;
	} else if (!strcmp(optarg,"mate-rev")) {
	  secondary_sort_method = MATEREV_SECONDARY_SORT;
	} else {
	  fprintf(stderr,"--sort2 must be none, orig, accession, mate-fwd, or mate-rev\n");
	  exit(9);
	}

      } else if (!strcmp(long_name,"mark-dups")) {
	mark_duplicates_p = true;
	print_unique_p = true;
	print_duplicates_p = true;

      } else if (!strcmp(long_name,"mark-first")) {
	mark_first_p = true;

      } else if (!strcmp(long_name,"dups-only")) {
	mark_duplicates_p = true;
	print_unique_p = false;
	print_duplicates_p = true;

      } else if (!strcmp(long_name,"uniq-only")) {
	mark_duplicates_p = true;
	print_unique_p = true;
	print_duplicates_p = false;

      } else if (!strcmp(long_name,"multiple-primaries")) {
	multiple_primaries_p = true;

      } else if (!strcmp(long_name,"no-sam-headers")) {
	sam_headers_p = false;

      } else {
	/* Shouldn't reach here */
	fprintf(stderr,"Don't recognize option %s.  For usage, run 'get-genome --help'",long_name);
	exit(9);
      }
      break;

    case 'D': user_genomedir = optarg; break;
    case 'd': 
      dbroot = (char *) CALLOC(strlen(optarg)+1,sizeof(char));
      strcpy(dbroot,optarg);
      break;

    case '^': print_program_version(); exit(0);
    case '?': print_program_usage(); exit(0);
    default: exit(9);
    }
  }
  argc -= optind;
  argv += optind;

  if (dbroot == NULL) {
    print_program_usage();
    exit(9);
  } else if (!strcmp(dbroot,"?")) {
    Datadir_avail_gmap_databases(stdout,user_genomedir);
    exit(0);
  } else {
    genomesubdir = Datadir_find_genomesubdir(&fileroot,&dbversion,user_genomedir,dbroot);
    iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			      strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
    sprintf(iitfile,"%s/%s.chromosome.iit",genomesubdir,fileroot);
    chromosome_iit = Univ_IIT_read(iitfile,/*readonlyp*/true,/*add_iit_p*/false);
    FREE(iitfile);

    FREE(dbversion);
    FREE(genomesubdir);
    FREE(fileroot);
    FREE(dbroot);

    nchromosomes = Univ_IIT_total_nintervals(chromosome_iit);
    chrlengths = Univ_IIT_chrlengths(chromosome_iit);
    chroffsets = MALLOC(nchromosomes * sizeof(Univcoord_T));
    chroffsets[0] = 0;
    for (i = 1; i < nchromosomes; i++) {
      chroffsets[i] = chroffsets[i-1] + chrlengths[i-1];
    }
    FREE(chrlengths);
  }
    
  /* Open all outputs, even if --split-output is not used */
  outputs = (FILE **) CALLOC((1+N_SPLIT_OUTPUTS),sizeof(FILE *));


  /* Inputs */
  ninputs = argc;
  sam_inputs = (FILE **) CALLOC(ninputs,sizeof(FILE *));
  headerlengths = (int *) CALLOC(ninputs,sizeof(int));
  ncells = (int *) CALLOC(ninputs,sizeof(int));
  for (filei = 0; filei < ninputs; filei++) {
    if ((sam_inputs[filei] = fopen(argv[filei],"r")) == NULL) {
      fprintf(stderr,"Cannot open SAM file %s\n",argv[i]);
      exit(9);
    }
  }

  stopwatch = Stopwatch_new();
  Stopwatch_start(stopwatch);
  fprintf(stderr,"Analyzing %d SAM files...\n",ninputs);

  linelengths = (Intlist_T) NULL;
  ncells_total = 0;
  for (filei = 0; filei < ninputs; filei++) {
    fp_sam = sam_inputs[filei];
    fileposition = headerlengths[filei] = SAM_header_length(&lastchar,fp_sam); /* Ignore lastchar */

    /* Take care of char read by SAM_header_length */
#ifdef HAVE_FSEEKO
    fseeko(fp_sam,-1,SEEK_CUR);
#else
    fseek(fp_sam,-1,SEEK_CUR);
#endif

    linelen = 0;
    ncells[filei] = 0;
    while (fgets(buffer,BUFFERLEN,fp_sam) != NULL) {
      /* printf("Read %s\n",buffer); */
      lastp = buffer;
      while ((p = index(lastp,'\n')) != NULL) {
	linelen += (p - lastp)/sizeof(char) + 1;

	linelengths = Intlist_push(linelengths,linelen);
	fileposition += linelen;
	ncells[filei] += 1;

	linelen = 0;
	lastp = p + 1;
      }
      linelen += strlen(lastp);
      /* printf("Adding %d to get linelen %d\n",strlen(buffer),linelen); */
    }

    ncells_total += ncells[filei];

    if (fileposition != Access_filesize(argv[filei])) {
      fprintf(stderr,"Something is wrong with parsing of SAM file %s\n",argv[filei]);
      fprintf(stderr,"Final file position using sortinfo: %llu\n",(unsigned long long) fileposition);
      fprintf(stderr,"File size of SAM output file:       %llu\n",(unsigned long long) Access_filesize(argv[0]));
      exit(9);
    } else {
      fprintf(stderr,"  File %d has %d SAM lines.\n",filei+1,ncells[filei]);
    }

  }

  fprintf(stderr,"Done with analysis (%.1f seconds).  Found %d SAM lines total.\n",
	  Stopwatch_stop(stopwatch),ncells_total);

  if (ncells_total == 0) {
    /* Exit without printing header */

  } else if (sam_headers_p == false) {
    /* Don't print SAM headers */

  } else {
    moveto(sam_inputs[0],0);
    headers = SAM_header_change_HD_tosorted(sam_inputs[0],headerlengths[0]);
  }

  linelengths = Intlist_reverse(linelengths);

  if (mark_duplicates_p == false) {
    process_without_dups(sam_inputs,headerlengths,ncells,ninputs,linelengths,ncells_total,
			 chromosome_iit,chroffsets,headers);
  } else {
    nmarked = process_with_dups(sam_inputs,headerlengths,ncells,ninputs,linelengths,ncells_total,
				chromosome_iit,chroffsets,headers);
    fprintf(stderr,"Marked %d out of %d SAM lines as duplicates (%.1f%%)\n",
	    nmarked,ncells_total,100.0*(double) nmarked/(double) (ncells_total));
  }

  for (filei = 0; filei < ninputs; filei++) {
    fclose(sam_inputs[filei]);
  }
  FREE(sam_inputs);
  FREE(headerlengths);
  FREE(ncells);

  if (headers != NULL) {
    Filestring_free(&headers);
  }

  /* SAM_header_touch(outputs,split_output_root,appendp); -- Don't want to destroy other SAM files */
  for (i = 1; i <= N_SPLIT_OUTPUTS; i++) {
    if (outputs[i] != NULL) {
      fclose(outputs[i]);
    }
  }
  FREE(outputs);


  Intlist_free(&linelengths);

  FREE(chroffsets);
  Univ_IIT_free(&chromosome_iit);

  return 0;
}
