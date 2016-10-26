static char rcsid[] = "$Id: iit_get.c 153955 2014-11-24 17:54:45Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>		/* For getopt */
#endif
#include <string.h>
#include <strings.h>		/* For index, rindex */
#include <ctype.h>
#include <math.h>		/* For log */
#include "bool.h"
#include "mem.h"
#include "access.h"
#include "intlist.h"
#include "univinterval.h"
#include "interval.h"
#include "iit-read-univ.h"
#include "iit-read.h"
#include "complement.h"
#include "fopen.h"
#include "parserange.h"
#include "getopt.h"

#define BUFLEN 1024


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/************************************************************************
 *   Program options
 ************************************************************************/

static bool value_matches_p = false;
static bool user_lowvalue_p = false;
static bool user_highvalue_p = false;
static double lowval;
static double highval;

static char *fieldstring = NULL;
static int fieldint = -1;
static bool annotationonlyp = false;
static bool exactp = false;
static bool sortp = false;
static bool signedp = true;
static int nflanking = 0;
static bool centerp = false;
static int centerlength = 0;
static bool centeruc = false;
static bool runlengthp = false;
static bool tallyp = false;
static bool zeroesp = false;
static bool statsp = false;
static bool force_label_p = false;
static bool force_coords_p = false;
static bool overall_total_p = false;

static struct option long_options[] = {
  /* Input options */
  {"lowval", required_argument, 0, 'a'}, /* value_matches_p, user_lowvalue_p, lowval */
  {"highval", required_argument, 0, 'b'},  /* value_matches_p, user_highvalue_p, highval */
  {"field", required_argument, 0, 'f'},	/* fieldstring */
  {"annotonly", no_argument, 0, 'A'},	/* annotationonlyp */
  {"exact", no_argument, 0, 0},		/* exactp */
  {"sort", no_argument, 0, 'S'},	/* sortp */
  {"unsigned", no_argument, 0, 'U'},	/* signedp */
  {"flanking", required_argument, 0, 'u'},	/* nflanking */
  {"center", required_argument, 0, 'c'}, /* centerp, centerlength */
  {"centeruc", no_argument, 0, 'H'}, /* centeruc */
  {"runlength", no_argument, 0, 'R'}, /* runlengthp */
  {"tally", no_argument, 0, 'T'}, /* tallyp */
  {"zeroes", no_argument, 0, 'Z'}, /* zeroesp */
  {"stats", no_argument, 0, 'N'}, /* statsp */
  {"label", no_argument, 0, 'L'}, /* force_label_p */
  {"coords", no_argument, 0, 'C'}, /* force_coords_p */

  /* Help options */
  {"version", no_argument, 0, 'V'}, /* print_program_version */
  {"help", no_argument, 0, '?'}, /* print_program_usage */
  {0, 0, 0, 0}
};

static void
print_program_version () {
  fprintf(stdout,"\n");
  fprintf(stdout,"iit_get: retrieval utility for Interval Index Trees\n");
  fprintf(stdout,"Part of GMAP package, version %s\n",PACKAGE_VERSION);
  fprintf(stdout,"Thomas D. Wu, Genentech, Inc.\n");
  fprintf(stdout,"Contact: twu@gene.com\n");
  fprintf(stdout,"\n");
  return;
}

static void
print_program_usage () {
  fprintf(stdout,"\
Usage: iit_get [OPTIONS...] iitfile query\n\
       iit_get [OPTIONS...] iitfile query types...\n\
       iit_get [OPTIONS...] iitfile\n\
\n\
where query is one of the following forms:\n\
\n\
   chr:start..end\n\
   chr:start\n\
   chr:\n\
   start..end\n\
   start\n\
   label\n\
\n\
Options\n\
  -f, --field=STRING      Show given field part of the annotation\n\
  -L, --label             Interpret query or queries as a label, even if it is numeric\n\
  -C, --coords            Interpret query or queries as coords\n\
  -A, --annotonly         Show annotation lines only (no headers)\n\
  -S, --sort              Sort results by coordinates\n\
  -U, --unsigned          Print all intervals as low..high, even those entered as reverse (high < low)\n\
  -u, --flanking=INT      Show flanking segments on left and right\n\
\n\
Options for specific IIT formats\n\
  -a, --lowval=DOUBLE     Low bound on a values IIT (default -Inf)\n\
  -b, --highval=DOUBLE    High bound on a values IIT (default +Inf)\n\
  -c, --center=INT        Align reads so given position is centered at given column\n\
  -H, --centeruc          Report only reads with upper-case letter at given position\n\
  -R, --runlength         Report runlength IIT file in tally format\n\
  -T, --tally             Report tally IIT file in tally format\n\
  -Z, --zeroes            Include zeroes in tally format\n\
  -N, --stats             Statistics (count, npositions, mean) of tally format\n\
\n\
  -V, --version           Show version\n\
  -?, --help              Show this help message\n\
\n\
\n\
The iit_get program retrieves segments from an iit file that overlap a\n\
given coordinate or pair of coordinates.  Retrieval is done in\n\
logarithmic time (for more details on IIT files, see Wu and Watanabe,\n\
Bioinformatics 21:1859-1875, 2005).  The start coordinate should be\n\
less than or equal to the end coordinate.  If they are not, the program\n\
will reverse them for you.  If only a single coordinate is provided,\n\
this is equivalent to providing the same number for the start and end\n\
coordinate.\n\
\n\
The given iit file may contain types (which can be displayed by using\n\
the -T flag of the iit_dump program).  These types may be used to\n\
filter the output of the iit file.  Multiple types may be specified,\n\
which indicates a disjunctive query, such that iit_get returns entries\n\
that match any one of the given tags.\n\
\n\
If no query is provided on the command line, the program will expect\n\
one or more queries from stdin, one per line.\n\
\n\
See also: iit_store, iit_dump\n\
");
  return;
}


#ifdef __STRICT_ANSI__
int getopt (int argc, char *const argv[], const char *optstring);
#endif

static char complCode[128] = COMPLEMENT_LC;

static void
print_forward (char *sequence, int start, int end) {
  int i;

  for (i = start; i <= end; i++) {
    printf("%c",sequence[i]);
  }
  return;
}

static void
print_complement (char *sequence, int end, int start) {
  int i;

  for (i = end; i >= start; i--) {
    printf("%c",complCode[(int) sequence[i]]);
  }
  return;
}

static void
print_spaces (int n) {
  while (--n >= 0) {
    printf(" ");
  }
  return;
}


/* Need to store just the part of the query specified (e.g., 1..10) */
static void
print_interval_centered (char *divstring, Chrpos_T coordstart, int index, IIT_T iit, int fieldint) {
  Interval_T interval;
  char *label, *annotation, *restofheader, centerchar;
  bool allocp;
  int annotlength, left, centerpos;

  if (fieldint < 0) {
    annotation = IIT_annotation(&restofheader,iit,index,&allocp);
    if (allocp == true) {
      FREE(restofheader);
    }
  } else {
    annotation = IIT_fieldvalue(iit,index,fieldint);
    allocp = true;
  }
  annotlength = strlen(annotation);
  if (annotation[annotlength-1] == '\n') {
    annotlength--;
  }

  interval = IIT_interval(iit,index);
  left = coordstart - Interval_low(interval); /* + length(query) - queryend */
  if (Interval_sign(interval) < 0) {
    centerpos = annotlength-left-1;
  } else {
    centerpos = left;
  }
  centerchar = annotation[centerpos];

  if (centeruc == true && islower(centerchar)) {
    if (fieldint >= 0 && allocp == true) {
      FREE(annotation);
    }
  } else {
    print_spaces(centerlength-left);
    if (Interval_sign(interval) < 0) {
      print_complement(annotation,annotlength-1,centerpos+1);
      printf("[%c]",complCode[(int) centerchar]);
      print_complement(annotation,centerpos-1,0);
    } else {
      print_forward(annotation,0,centerpos-1);
      printf("[%c]",centerchar);
      print_forward(annotation,centerpos+1,annotlength-1);
    }
    print_spaces(centerlength+left-annotlength);
    if (fieldint >= 0 && allocp == true) {
      FREE(annotation);
    }
  
    printf("\t");
    if (Interval_type(interval) > 0) {
      printf("%s\t",IIT_typestring(iit,Interval_type(interval)));
    }

    if (divstring != NULL) {
      if (Interval_sign(interval) < 0) {
	printf("-%s:",divstring);
      } else {
	printf("+%s:",divstring);
      }
    }

    if (signedp == false) {
      printf("%u..%u",Interval_low(interval),Interval_high(interval));
    } else if (Interval_sign(interval) < 0) {
      printf("%u..%u",Interval_high(interval),Interval_low(interval));
    } else {
      printf("%u..%u",Interval_low(interval),Interval_high(interval));
    }
    printf("\t");

    label = IIT_label(iit,index,&allocp);
    printf("%s",label);
    if (allocp == true) {
      FREE(label);
    }
    printf("\n");
  }

  return;
}

static char *
get_total_tally (long int *tally, char *ptr) {
  int n;
  char *end;

  if ((end = index(ptr,'\n')) == NULL) {
    fprintf(stderr,"Premature end of line %s\n",ptr);
    return 0;
  }
  /* fprintf(stderr,"Getting tally for %.*s\n",end-ptr,ptr); */

  while (ptr < end) {
    while (ptr < end && !isdigit((int) *ptr)) {
      ptr++;
    }
    if (ptr < end) {
      sscanf(ptr,"%d",&n);
      (*tally) += n;
      while (ptr < end && !isspace(*ptr)) {
	ptr++;
      }
      while (ptr < end && isspace(*ptr)) {
	ptr++;
      }
    }
  }

  return ptr;
}


static void
print_line (char *ptr) {
  while (*ptr != '\0' && *ptr != '\n') {
    printf("%c",*ptr);
    ptr++;
  }
  return;
}



/* Need to store just the part of the query specified (e.g., 1..10) */
static long int
print_interval_tally (Chrpos_T *lastcoord, char *divstring, Chrpos_T coordstart, Chrpos_T coordend,
		      int indexi, IIT_T iit, bool zeroesp) {
  long int total = 0, subtotal;
  Interval_T interval;
  char *annotation, *restofheader, *ptr, *nextptr;
  bool allocp;
  Chrpos_T chrpos, intervalend;

  annotation = IIT_annotation(&restofheader,iit,indexi,&allocp);

  interval = IIT_interval(iit,indexi);
  chrpos = Interval_low(interval);
  intervalend = Interval_high(interval);

  ptr = annotation;

  if (zeroesp == true) {
    while (*lastcoord < chrpos) {
      if (statsp == false) {
	printf("%s\t%u\t%d\t\n",divstring,*lastcoord,0);
      }
      (*lastcoord)++;
    }
  }

  while (chrpos < coordstart) {
    if ((ptr = index(ptr,'\n')) == NULL) {
      fprintf(stderr,"Premature end of tally from %u to %u\n",
	      Interval_low(interval),Interval_high(interval));
      return total;
    } else {
      ptr++;
    }
    chrpos++;
  }

  while (chrpos <= intervalend && chrpos <= coordend) {
    subtotal = 0;
    nextptr = get_total_tally(&subtotal,ptr);
    if (subtotal > 0 || zeroesp == true) {
      if (statsp == false) {
	printf("%s\t%u\t%ld\t",divstring,chrpos,total);
	print_line(ptr);
	printf("\n");
      }
    }
    total += subtotal;
    ptr = nextptr;
    if ((ptr = index(ptr,'\n')) == NULL) {
      fprintf(stderr,"Premature end of tally from %u to %u\n",
	      Interval_low(interval),Interval_high(interval));
      return total;
    } else {
      ptr++;
    }
    chrpos++;
  }

  *lastcoord = chrpos;

  if (allocp == true) {
    FREE(restofheader);
  }
  
  return total;
}


/* Need to store just the part of the query specified (e.g., 1..10) */
static void
print_interval_runlength (Chrpos_T *lastcoord, char *divstring, Chrpos_T coordstart, Chrpos_T coordend,
			  int indexi, IIT_T iit, bool zeroesp) {
  Interval_T interval;
  char *label;
  bool allocp;
  Chrpos_T chrpos, intervalend;
  int value;

  label = IIT_label(iit,indexi,&allocp);
  value = atoi(label);

  interval = IIT_interval(iit,indexi);
  chrpos = Interval_low(interval);
  intervalend = Interval_high(interval);

  if (zeroesp == true) {
    while (*lastcoord < chrpos) {
      printf("%s\t%u\t%d\n",divstring,*lastcoord,0);
      (*lastcoord)++;
    }
  }

  while (chrpos < coordstart) {
    chrpos++;
  }

  while (chrpos <= intervalend && chrpos <= coordend) {
    printf("%s\t%u\t%d\n",divstring,chrpos,value);
    chrpos++;
  }

  *lastcoord = chrpos;

  if (allocp == true) {
    FREE(label);
  }
  
  return;
}


/************************************************************************
 *   Totals
 ************************************************************************/

/* Need to store just the part of the query specified (e.g., 1..10) */
static void
compute_totals_tally (long int *total, int *n, Chrpos_T coordstart, Chrpos_T coordend,
		      int indexi, IIT_T iit) {
  Interval_T interval;
  char *annotation, *restofheader, *ptr;
  bool allocp;
  Chrpos_T chrpos, intervalend;

  annotation = IIT_annotation(&restofheader,iit,indexi,&allocp);

  interval = IIT_interval(iit,indexi);
  chrpos = Interval_low(interval);
  intervalend = Interval_high(interval);

  ptr = annotation;

  while (chrpos < coordstart) {
    if ((ptr = index(ptr,'\n')) == NULL) {
      fprintf(stderr,"Premature end of tally from %u to %u\n",
	      Interval_low(interval),Interval_high(interval));
      return;
    } else {
      ptr++;
    }
    chrpos++;
  }

  while (chrpos <= intervalend && chrpos <= coordend) {
    ptr = get_total_tally(&(*total),ptr);
    *n += 1;
    ptr++;
    chrpos++;
  }

  if (allocp == true) {
    FREE(restofheader);
  }
  
  return;
}


/* Need to store just the part of the query specified (e.g., 1..10) */
static double
compute_logtotal_tally (long int *total, int *n, Chrpos_T coordstart, Chrpos_T coordend,
			int indexi, IIT_T iit) {
  double logtotal = 0.0;
  Interval_T interval;
  char *annotation, *restofheader, *ptr;
  bool allocp;
  Chrpos_T chrpos, intervalend;
  long int count;

  annotation = IIT_annotation(&restofheader,iit,indexi,&allocp);

  interval = IIT_interval(iit,indexi);
  chrpos = Interval_low(interval);
  intervalend = Interval_high(interval);

  ptr = annotation;

  while (chrpos < coordstart) {
    if ((ptr = index(ptr,'\n')) == NULL) {
      fprintf(stderr,"Premature end of tally from %u to %u\n",
	      Interval_low(interval),Interval_high(interval));
      return logtotal;
    } else {
      ptr++;
    }
    chrpos++;
  }

  while (chrpos <= intervalend && chrpos <= coordend) {
    count = 0;
    ptr = get_total_tally(&count,ptr);

    logtotal += log((double) count + 1.0);
    *total += count;
    *n += 1;
    ptr++;
    chrpos++;
  }

  if (allocp == true) {
    FREE(restofheader);
  }
  
  return logtotal;
}



/* coordstart used only if centerp or tallyp is true */
static long int
print_interval (Chrpos_T *lastcoord, long int total,
		char *divstring, Chrpos_T coordstart, Chrpos_T coordend, 
		int index, IIT_T iit, int ndivs, int fieldint) {
  Interval_T interval;
  char *label, *annotation, *restofheader;
  bool allocp;

  if (centerp == true) {
    print_interval_centered(divstring,coordstart,index,iit,fieldint);
    return 0;
  } else if (tallyp == true) {
    total += print_interval_tally(&(*lastcoord),divstring,coordstart,coordend,index,iit,zeroesp);
    return total;
  } else if (runlengthp == true) {
    print_interval_runlength(&(*lastcoord),divstring,coordstart,coordend,index,iit,zeroesp);
    return 0;
  }

  if (annotationonlyp == false) {
    label = IIT_label(iit,index,&allocp);
    printf(">%s ",label);
    if (allocp == true) {
      FREE(label);
    }
      
    if (ndivs > 1) {
      if (divstring == NULL) {
	/* For example, if interval was retrieved by label */
	divstring = IIT_divstring_from_index(iit,index);
      }
      printf("%s:",divstring);
    }

    debug(printf("index is %d\n",index));
    interval = IIT_interval(iit,index);
    if (signedp == false) {
      printf("%u..%u",Interval_low(interval),Interval_high(interval));
    } else if (Interval_sign(interval) < 0) {
      printf("%u..%u",Interval_high(interval),Interval_low(interval));
    } else {
      printf("%u..%u",Interval_low(interval),Interval_high(interval));
    }
    if (Interval_type(interval) > 0) {
      printf(" %s",IIT_typestring(iit,Interval_type(interval)));
    }
#if 0
    /* Unnecessary because of "\n" after restofheader below */
    if (IIT_version(iit) < 5) {
      printf("\n");
    }
#endif
  }


  if (fieldint < 0) {
    annotation = IIT_annotation(&restofheader,iit,index,&allocp);
    printf("%s\n",restofheader);
    printf("%s",annotation);
    if (allocp == true) {
      FREE(restofheader);
    }
  } else {
    annotation = IIT_annotation(&restofheader,iit,index,&allocp);
    printf("%s\n",restofheader);
    if (allocp == true) {
      FREE(restofheader);
    }
    annotation = IIT_fieldvalue(iit,index,fieldint);
    printf("%s\n",annotation);
    FREE(annotation);
  }

  return 0;
}


static void
print_interval_univ (Univcoord_T coordstart, Univcoord_T coordend, int index,
		     Univ_IIT_T chromosome_iit) {
  Univinterval_T interval;
  char *label, *annotation, *restofheader;
  bool allocp;

  if (annotationonlyp == true) {
    annotation = Univ_IIT_annotation(&restofheader,chromosome_iit,index,&allocp);
    printf("%s",annotation);
    if (allocp) {
      FREE(restofheader);
    }

  } else {
    label = Univ_IIT_label(chromosome_iit,index,&allocp);
    printf(">%s ",label);
    if (allocp == true) {
      FREE(label);
    }
      
    debug(printf("index is %d\n",index));
    interval = Univ_IIT_interval(chromosome_iit,index);
    printf("%llu..%llu",
	   (unsigned long long) Univinterval_low(interval),
	   (unsigned long long) Univinterval_high(interval));
    if (Univinterval_type(interval) > 0) {
      printf(" %s",Univ_IIT_typestring(chromosome_iit,Univinterval_type(interval)));
    }

    annotation = Univ_IIT_annotation(&restofheader,chromosome_iit,index,&allocp);
    printf("%s\n",restofheader);
    printf("%s",annotation);
    if (allocp) {
      FREE(restofheader);
    }

  }

  return;
}


static int *
get_matches (int *nmatches, char **divstring, Univcoord_T *coordstart, Univcoord_T *coordend,
	     int **leftflanks, int *nleftflanks, int **rightflanks, int *nrightflanks,
	     char *query, char *typestring, IIT_T *iit, char *filename) {
  int *matches;
  bool revcomp;
  int typeint;

  debug(printf("Entering get_matches with query %s.\n",query));

  if (force_label_p == true || (force_coords_p == false && Parserange_query(&(*divstring),&(*coordstart),&(*coordend),&revcomp,query,filename) == false)) {
    /* Treat query as a label */
    *divstring = (char *) NULL;
    if (*iit == NULL) {
      /* Read no divs */
      if ((*iit = IIT_read(filename,/*name*/NULL,true,/*divread*/READ_NONE,/*divstring*/NULL,
			   /*add_iit_p*/true,/*labels_read_p*/true)) == NULL) {
	if (Access_file_exists_p(filename) == false) {
	  fprintf(stderr,"Cannot read file %s\n",filename);
	} else {
	  fprintf(stderr,"File %s appears to be an invalid IIT file\n",filename);
	}
	exit(9);

      } else if (fieldstring != NULL) {
	if ((fieldint = IIT_fieldint(*iit,fieldstring)) < 0) {
	  fprintf(stderr,"No field %s defined in iit file.\n",fieldstring);
	  exit(9);
	}
      }
    }
    *nleftflanks = *nrightflanks = 0;
    matches = IIT_find(&(*nmatches),*iit,query);

  } else {
    if (*iit == NULL) {
      if ((*iit = IIT_read(filename,/*name*/NULL,true,/*divread*/READ_ONE,*divstring,
			   /*add_iit_p*/true,/*labels_read_p*/false)) == NULL) {
	if (Access_file_exists_p(filename) == false) {
	  fprintf(stderr,"Cannot read file %s\n",filename);
	} else {
	  fprintf(stderr,"File %s appears to be an invalid IIT file\n",filename);
	}
	exit(9);

      } else if (fieldstring != NULL) {
	if ((fieldint = IIT_fieldint(*iit,fieldstring)) < 0) {
	  fprintf(stderr,"No field %s defined in iit file.\n",fieldstring);
	  exit(9);
	}
      }
    }

    if (typestring == NULL) {
      /* Treat query as coordinates, without a typestring */
      if (exactp == true) {
	matches = IIT_get_exact_multiple(&(*nmatches),*iit,*divstring,*coordstart,*coordend,/*type*/0);
      } else {
	matches = IIT_get(&(*nmatches),*iit,*divstring,*coordstart,*coordend,sortp);
      }
      if (nflanking > 0) {
	IIT_get_flanking(&(*leftflanks),&(*nleftflanks),&(*rightflanks),&(*nrightflanks),*iit,*divstring,
			 *coordstart,*coordend,nflanking,/*sign*/0);
      }

    } else if ((typeint = IIT_typeint(*iit,typestring)) < 0) {
      fprintf(stderr,"No such type as %s.\n",typestring);
#if 0
      /* Treat query as coordinates, without a typestring */
      matches = IIT_get(&(*nmatches),*iit,*divstring,*coordstart,*coordend,sortp);
      if (nflanking > 0) {
	IIT_get_flanking(&(*leftflanks),&(*nleftflanks),&(*rightflanks),&(*nrightflanks),*iit,*divstring,
			 *coordstart,*coordend,nflanking,/*sign*/0);
      }
#else
      matches = (int *) NULL;
      *nleftflanks = *nrightflanks = 0;
      nmatches = 0;
#endif      

    } else {
      /* Treat query as coordinates, with a typestring */
      if (exactp == true) {
	matches = IIT_get_exact_multiple(&(*nmatches),*iit,*divstring,*coordstart,*coordend,typeint);
      } else {
	matches = IIT_get_typed(&(*nmatches),*iit,*divstring,*coordstart,*coordend,typeint,sortp);
      }
      if (nflanking > 0) {
	debug(printf("Running IIT_get_flanking_typed\n"));
	IIT_get_flanking_typed(&(*leftflanks),&(*nleftflanks),&(*rightflanks),&(*nrightflanks),*iit,*divstring,
			       *coordstart,*coordend,nflanking,typeint,/*sign*/0);
      }
    }

  }

  return matches;
}

static int *
get_matches_univ (int *nmatches, char **divstring, Univcoord_T *coordstart, Univcoord_T *coordend,
		  int **leftflanks, int *nleftflanks, int **rightflanks, int *nrightflanks,
		  char *query, char *typestring, Univ_IIT_T *chromosome_iit, char *filename) {
  int *matches;
  bool revcomp;
  int typeint;

  debug(printf("Entering get_matches_univ with query %s.\n",query));

  if (force_label_p == true || (force_coords_p == false && Parserange_query(&(*divstring),&(*coordstart),&(*coordend),&revcomp,query,filename) == false)) {
    /* Treat query as a label */
    *divstring = (char *) NULL;
    if (*chromosome_iit == NULL) {
      /* Read no divs */
      if ((*chromosome_iit = Univ_IIT_read(filename,/*readonlyp*/true,/*add_iit_p*/true)) == NULL) {
	if (Access_file_exists_p(filename) == false) {
	  fprintf(stderr,"Cannot read file %s\n",filename);
	} else {
	  fprintf(stderr,"File %s appears to be an invalid IIT file\n",filename);
	}
	exit(9);
      }
    }
    *nleftflanks = *nrightflanks = 0;
    matches = Univ_IIT_find(&(*nmatches),*chromosome_iit,query);

  } else {
    if (*chromosome_iit == NULL) {
      if ((*chromosome_iit = Univ_IIT_read(filename,/*readonlyp*/true,/*add_iit_p*/true)) == NULL) {
	if (Access_file_exists_p(filename) == false) {
	  fprintf(stderr,"Cannot read file %s\n",filename);
	} else {
	  fprintf(stderr,"File %s appears to be an invalid IIT file\n",filename);
	}
	exit(9);
      }
    }

    if (typestring == NULL) {
      /* Treat query as coordinates, without a typestring */
      matches = Univ_IIT_get(&(*nmatches),*chromosome_iit,*coordstart,*coordend);
      if (nflanking > 0) {
	fprintf(stderr,"Flanking not supported on chromosome IIT\n");
      }
      *nleftflanks = *nrightflanks = 0;

    } else if ((typeint = Univ_IIT_typeint(*chromosome_iit,typestring)) < 0) {
      fprintf(stderr,"No such type as %s.\n",typestring);
#if 0
      /* Treat query as coordinates, without a typestring */
      matches = Univ_IIT_get(&(*nmatches),*chromosome_iit,*coordstart,*coordend);
      if (nflanking > 0) {
	fprintf(stderr,"Flanking not supported on chromosome IIT\n");
      }
      *nleftflanks = *nrightflanks = 0;
#else
      matches = (int *) NULL;
      *nleftflanks = *nrightflanks = 0;
      *nmatches = 0;
#endif      

    } else {
      /* Treat query as coordinates, with a typestring */
      fprintf(stderr,"Ignoring types on chromosome IIT\n");
      matches = Univ_IIT_get(&(*nmatches),*chromosome_iit,*coordstart,*coordend);
      *nleftflanks = *nrightflanks = 0;
    }

  }

  return matches;
}


static int *
get_matches_multiple_typed (int *nmatches, char **divstring, Univcoord_T *coordstart, Univcoord_T *coordend,
			    int **leftflanks, int *nleftflanks, int **rightflanks, int *nrightflanks,
			    char *query, int *types, int ntypes, IIT_T *iit, char *filename) {
  int *matches;
  bool revcomp;

  if (force_label_p == true || (force_coords_p == false && Parserange_query(&(*divstring),&(*coordstart),&(*coordend),&revcomp,query,filename) == false)) {
    /* Not expecting a label */
    abort();
  }

  if ((*iit = IIT_read(filename,/*name*/NULL,true,/*divread*/READ_ONE,*divstring,/*add_iit_p*/true,
		       /*labels_read_p*/false)) == NULL) {
    if (Access_file_exists_p(filename) == false) {
      fprintf(stderr,"Cannot read file %s\n",filename);
    } else {
      fprintf(stderr,"File %s appears to be an invalid IIT file\n",filename);
    }
    exit(9);

  } else if (fieldstring != NULL) {
    if ((fieldint = IIT_fieldint(*iit,fieldstring)) < 0) {
      fprintf(stderr,"No field %s defined in iit file.\n",fieldstring);
      exit(9);
    }
  }

  matches = IIT_get_multiple_typed(&(*nmatches),*iit,*divstring,*coordstart,*coordend,types,ntypes,sortp);
  if (nflanking > 0) {
    IIT_get_flanking_multiple_typed(&(*leftflanks),&(*nleftflanks),&(*rightflanks),&(*nrightflanks),*iit,*divstring,
				    *coordstart,*coordend,nflanking,types,ntypes);
  }

  return matches;
}


static int
int_cmp (const void *x, const void *y) {
  int a = * (int *) x;
  int b = * (int *) y;

  if (a < b) {
    return -1;
  } else if (b < a) {
    return +1;
  } else {
    return 0;
  }
}


static int *
match_intersection (int *nmatches, int *matches1, int nmatches1, int *matches2, int nmatches2) {
  int *intersection;
  int i, j, k;

  qsort(matches1,nmatches1,sizeof(int),int_cmp);
  qsort(matches2,nmatches2,sizeof(int),int_cmp);
  
  i = j = 0;
  k = 0;
  while (i < nmatches1 && j < nmatches2) {
    if (matches1[i] < matches2[j]) {
      i++;
    } else if (matches1[i] > matches2[j]) {
      j++;
    } else {
      i++;
      j++;
      k++;
    }
  }

  if ((*nmatches = k) == 0) {
    FREE(matches1);
    return (int *) NULL;

  } else {
    intersection = (int *) CALLOC(*nmatches,sizeof(int));

    i = j = 0;
    k = 0;
    while (i < nmatches1 && j < nmatches2) {
      if (matches1[i] < matches2[j]) {
	i++;
      } else if (matches1[i] > matches2[j]) {
	j++;
      } else {
	i++;
	j++;
	intersection[k++] = matches1[i];
      }
    }

    FREE(matches1);
    return intersection;
  }
}


int 
main (int argc, char *argv[]) {
  char *filename;
  char *divstring = NULL, *lasttypestring, *ptr;
  Univcoord_T univ_coordstart, univ_coordend;
  Chrpos_T coordstart, coordend, lastcoord = 0U;
  char Buffer[BUFLEN], nocomment[BUFLEN], query[BUFLEN], typestring[BUFLEN];
  int typeint, *types, c;
  int nargs, ntypes, ndivs;
  int *value_matches = NULL, *matches = NULL;
  int n_value_matches = 0, nmatches = 0, i;
  int *leftflanks, *rightflanks, nleftflanks = 0, nrightflanks = 0;
  long int total;
  int n;
  IIT_T iit = NULL;
  Univ_IIT_T chromosome_iit = NULL;
  bool skipp, universalp;
  
  int opt;
  extern int optind;
  extern char *optarg;
  const char *long_name;
  int long_option_index = 0;

  while ((opt = getopt_long(argc,argv,"a:b:f:LCASUu:c:HRTZN",long_options,&long_option_index)) != -1) {
    switch (opt) {
    case 0:
      long_name = long_options[long_option_index].name;
      if (!strcmp(long_name,"version")) {
	print_program_version();
	exit(0);
      } else if (!strcmp(long_name,"help")) {
	print_program_usage();
	exit(0);
      } else if (!strcmp(long_name,"exact")) {
	exactp = true;
      } else {
	/* Shouldn't reach here */
	fprintf(stderr,"Don't recognize option %s.  For usage, run 'gsnap --help'",long_name);
	exit(9);
      }
      break;

    case 'a': lowval = atof(optarg); user_lowvalue_p = true; value_matches_p = true; break;
    case 'b': highval = atof(optarg); user_highvalue_p = true; value_matches_p = true; break;

    case 'f': fieldstring = optarg; break;
    case 'L': force_label_p = true; break;
    case 'C': force_coords_p = true; break;
    case 'A': annotationonlyp = true; break;
    case 'S': sortp = true; break;
    case 'U': signedp = false; break;
    case 'u': nflanking = atoi(optarg); break;

    case 'c': centerp = true; centerlength = atoi(optarg); break;
    case 'H': centeruc = true; break;
    case 'R': runlengthp = true; break;
    case 'T': tallyp = true; break;
    case 'Z': zeroesp = true; break;
    case 'N': statsp = true; break;

    case 'V': print_program_version(); exit(0);
    case '?': print_program_usage(); exit(0);
    default: exit(9);
    }
  }

  argc -= (optind - 1);
  argv += (optind - 1);

  if (argc <= 1) {
    fprintf(stderr,"Need to specify an iit file.  Type \"iit_get --help\" for help.\n");
    exit(9);
  } else {
    filename = argv[1];
  }

  if (value_matches_p == true) {
    if ((iit = IIT_read(filename,/*name*/NULL,true,/*divread*/READ_ALL,/*divstring*/NULL,
			/*add_iit_p*/true,/*labels_read_p*/true)) == NULL) {
      if (Access_file_exists_p(filename) == false) {
	fprintf(stderr,"Cannot read file %s\n",filename);
      } else {
	fprintf(stderr,"File %s appears to be an invalid IIT file\n",filename);
      }
      exit(9);

    } else if (fieldstring != NULL) {
      if ((fieldint = IIT_fieldint(iit,fieldstring)) < 0) {
	fprintf(stderr,"No field %s defined in iit file.\n",fieldstring);
	exit(9);
      }
    }

    if (IIT_valuep(iit) == false) {
      fprintf(stderr,"Error: This IIT file does not have values stored\n");
      exit(9);
    }

    if (user_lowvalue_p == true && user_highvalue_p == true) {
      if (lowval > highval) {
	fprintf(stderr,"Cannot have lowval %f > highval %f\n",lowval,highval);
	exit(9);
      } else {
	value_matches = IIT_get_values_between(&n_value_matches,iit,lowval,highval,/*sortp*/false);
      }
    } else if (user_lowvalue_p == true) {
      value_matches = IIT_get_values_above(&n_value_matches,iit,lowval,/*sortp*/false);
      
    } else { /* user_highvalue_p == true */
      value_matches = IIT_get_values_below(&n_value_matches,iit,highval,/*sortp*/false);
    }
  }

  if (0 && statsp == true && argc == 2) {
    /* Want total over entire IIT */
    if ((iit = IIT_read(filename,NULL,true,/*divread*/READ_ALL,/*divstring*/NULL,/*add_iit_p*/true,
			/*labels_read_p*/false)) == NULL) {
      if (Access_file_exists_p(filename) == false) {
	fprintf(stderr,"Cannot read file %s\n",filename);
      } else {
	fprintf(stderr,"File %s appears to be an invalid IIT file\n",filename);
      }
      exit(9);

    } else if (fieldstring != NULL) {
      if ((fieldint = IIT_fieldint(iit,fieldstring)) < 0) {
	fprintf(stderr,"No field %s defined in iit file.\n",fieldstring);
	exit(9);
      }
    } 

    total = 0;
    n = 0;
    for (i = 0; i < IIT_total_nintervals(iit); i++) {
      debug(printf("index = %d\n",matches[i]));
      compute_totals_tally(&total,&n,/*coordstart*/0,/*coordend*/-1U,i,iit);
    }
    printf("counts:%ld non-zero-positions:%u mean-over-nonzero:%.3f\n",total,n,(double) total/(double) n);
      
    IIT_free(&iit);
    return 0;

  } else if ((universalp = IIT_universalp(filename,/*add_iit_p*/true)) == true) {
    chromosome_iit = Univ_IIT_read(filename,/*readonlyp*/true,/*add_iit_p*/true);
    if (argc != 3) {
      fprintf(stderr,"For chromosome IIT file, need to specify a query on the command line\n");
      exit(9);
    } else {
      /* Try as 0:<iitfile> 1:<query> */
      matches = get_matches_univ(&nmatches,&divstring,&univ_coordstart,&univ_coordend,
				 &leftflanks,&nleftflanks,&rightflanks,&nrightflanks,
				 argv[2],/*typestring*/NULL,&chromosome_iit,filename);

      for (i = 0; i < nmatches; i++) {
	debug(printf("\nindex = %d\n",matches[i]));
	print_interval_univ(univ_coordstart,univ_coordend,matches[i],chromosome_iit);
      }
    }

    if (divstring != NULL) {
      FREE(divstring);
    }
    Univ_IIT_free(&chromosome_iit);
    return 0;

  } else if (argc == 2 && value_matches_p == true) {
    /* Note: Could potentially handle input from stdin, but currently just deal with value_matches */

    ndivs = IIT_ndivs(iit);
    for (i = 0; i < n_value_matches; i++) {
      debug(printf("\nindex = %d\n",matches[i]));
      print_interval(&lastcoord,/*total*/0,/*divstring*/NULL,/*coordstart*/0,/*coordend*/0,
		     value_matches[i],iit,ndivs,fieldint);
    }
    
    FREE(value_matches);
    IIT_free(&iit);
    return 0;

  } else if (argc == 2) {
    debug(printf("Running under argc 2\n"));
    /* Expecting input from stdin */
      
    if ((iit = IIT_read(filename,NULL,true,/*divread*/READ_ALL,/*divstring*/NULL,/*add_iit_p*/true,
			/*labels_read_p*/true)) == NULL) {
      if (Access_file_exists_p(filename) == false) {
	fprintf(stderr,"Cannot read file %s\n",filename);
      } else {
	fprintf(stderr,"File %s appears to be an invalid IIT file\n",filename);
      }
      exit(9);
    } else if (fieldstring != NULL) {
      if ((fieldint = IIT_fieldint(iit,fieldstring)) < 0) {
	fprintf(stderr,"No field %s defined in iit file.\n",fieldstring);
	exit(9);
      }
    }
	
    while (fgets(Buffer,BUFLEN,stdin) != NULL) {
      if ((ptr = rindex(Buffer,'\n')) != NULL) {
	*ptr = '\0';
      }
      strcpy(nocomment,Buffer);
      if ((ptr = rindex(nocomment,'#')) != NULL) {
	*ptr = '\0';
      }

      skipp = false;

      if ((nargs = sscanf(nocomment,"%s %s",query,typestring)) == 2) {
	debug(printf("typestring is %s\n",typestring));
	matches = get_matches(&nmatches,&divstring,&univ_coordstart,&univ_coordend,
			      &leftflanks,&nleftflanks,&rightflanks,&nrightflanks,
			      query,typestring,&iit,filename);
	coordstart = (Chrpos_T) univ_coordstart;
	coordend = (Chrpos_T) univ_coordend;

      } else if (nargs == 1) {
	debug(printf("typestring is NULL\n"));
	matches = get_matches(&nmatches,&divstring,&univ_coordstart,&univ_coordend,
			      &leftflanks,&nleftflanks,&rightflanks,&nrightflanks,
			      query,/*typestring*/NULL,&iit,filename);
	coordstart = (Chrpos_T) univ_coordstart;
	coordend = (Chrpos_T) univ_coordend;

      } else {
	fprintf(stderr,"Can't parse line %s.  Ignoring.\n",nocomment);
	skipp = true;
      }
	
      total = 0;
      if (skipp == false) {
	fprintf(stdout,"# Query: %s\n",Buffer);
	ndivs = IIT_ndivs(iit);
	if (nflanking > 0) {
	  for (i = nleftflanks-1; i >= 0; i--) {
	    debug(printf("\nleft index = %d\n",leftflanks[i]));
	    print_interval(&lastcoord,/*total*/0,divstring,coordstart,coordend,leftflanks[i],iit,ndivs,fieldint);
	  }
	  printf("====================\n");
	  FREE(leftflanks);
	}

	lastcoord = coordstart;
	for (i = 0; i < nmatches; i++) {
	  debug(printf("\nindex = %d\n",matches[i]));
	  total = print_interval(&lastcoord,total,divstring,coordstart,coordend,matches[i],iit,ndivs,fieldint);
	}

	if (nflanking > 0) {
	  printf("====================\n");
	  for (i = 0; i < nrightflanks; i++) {
	    debug(printf("\nright index = %d\n",rightflanks[i]));
	    print_interval(&lastcoord,/*total*/0,divstring,coordstart,coordend,rightflanks[i],iit,ndivs,fieldint);
	  }
	  FREE(rightflanks);
	}

	if (zeroesp == true) {
	  while (lastcoord <= coordend) {
	    printf("%s\t%u\t%d\n",divstring,lastcoord,0);
	    lastcoord++;
	  }
	}
      }

      if (divstring != NULL) {
	FREE(divstring);
      }
      FREE(matches);
      printf("%ld\n",total);
      fprintf(stdout,"# End\n");
      fflush(stdout);
    }

    IIT_free(&iit);
    return 0;

  } else {
    /* Get coordinates/type from command line */
    if (argc == 3) {
      /* Try as 0:<iitfile> 1:<query> */
      matches = get_matches(&nmatches,&divstring,&univ_coordstart,&univ_coordend,
			    &leftflanks,&nleftflanks,&rightflanks,&nrightflanks,
			    argv[2],/*typestring*/NULL,&iit,filename);
      coordstart = (Chrpos_T) univ_coordstart;
      coordend = (Chrpos_T) univ_coordend;
    
    } else if (argc == 4) {
      /* Try as 0:<iitfile> 1:<query> 2:<type> */
      debug(printf("Running under argc 4\n"));
      matches = get_matches(&nmatches,&divstring,&univ_coordstart,&univ_coordend,
			    &leftflanks,&nleftflanks,&rightflanks,&nrightflanks,
			    argv[2],argv[3],&iit,filename);
      coordstart = (Chrpos_T) univ_coordstart;
      coordend = (Chrpos_T) univ_coordend;

    } else {
      types = (int *) CALLOC(argc-3,sizeof(int));
      for (c = 3, ntypes = 0; c < argc; c++) {
	if ((typeint = IIT_typeint(iit,argv[c])) < 0) {
	  fprintf(stderr,"No such type as %s.  Ignoring the type.\n",argv[c]);
	} else {
	  types[ntypes++] = typeint;
	  lasttypestring = argv[c];
	}
      }
      if (ntypes == 0) {
	matches = get_matches(&nmatches,&divstring,&univ_coordstart,&univ_coordend,
			      &leftflanks,&nleftflanks,&rightflanks,&nrightflanks,
			      argv[2],/*typestring*/NULL,&iit,filename);
	coordstart = (Chrpos_T) univ_coordstart;
	coordend = (Chrpos_T) univ_coordend;

      } else if (ntypes == 1) {
	matches = get_matches(&nmatches,&divstring,&univ_coordstart,&univ_coordend,
			      &leftflanks,&nleftflanks,&rightflanks,&nrightflanks,
			      argv[2],lasttypestring,&iit,filename);
	coordstart = (Chrpos_T) univ_coordstart;
	coordend = (Chrpos_T) univ_coordend;

      } else {
	matches = get_matches_multiple_typed(&nmatches,&divstring,&univ_coordstart,&univ_coordend,
					     &leftflanks,&nleftflanks,&rightflanks,&nrightflanks,
					     argv[2],types,ntypes,&iit,filename);
	coordstart = (Chrpos_T) univ_coordstart;
	coordend = (Chrpos_T) univ_coordend;
      }
    }

    if (value_matches_p == true) {
      matches = match_intersection(&nmatches,/*matches1*/matches,/*nmatches1*/nmatches,
				   /*matches2*/value_matches,/*nmatches2*/n_value_matches);
      FREE(value_matches);
    }

#if 0
    if (centerp == true) {
      print_spaces(centerlength);
      printf("*");
      print_spaces(centerlength-1);
      printf("\n");
    }
#endif

    if (statsp == true) {
      total = 0;
      n = 0;
      for (i = 0; i < nmatches; i++) {
	debug(printf("index = %d\n",matches[i]));
	compute_totals_tally(&total,&n,coordstart,coordend,matches[i],iit);
      }
      n = coordend - coordstart + 1;
      printf("counts:%ld width:%u mean:%.3f\n",total,n,(double)total/(double) n);

#if 0
    } else if (geomeanp == true) {
      logtotal = 0.0;
      total = 0;
      n = 0;
      for (i = 0; i < nmatches; i++) {
	debug(printf("index = %d\n",matches[i]));
	logtotal = compute_logtotal_tally(&total,&n,coordstart,coordend,matches[i],iit);
      }
      printf("geomean:%f totalcounts:%ld posrange:%d\n",
	     exp(logtotal/(double) (coordend - coordstart + 1)) - 1.0,total,n);
#endif

    } else {
      ndivs = IIT_ndivs(iit);
      if (nflanking > 0) {
	for (i = nleftflanks-1; i >= 0; i--) {
	  debug(printf("\nleft index = %d\n",leftflanks[i]));
	  print_interval(&lastcoord,/*total*/0,divstring,coordstart,coordend,leftflanks[i],iit,ndivs,fieldint);
	}
	printf("====================\n");
	FREE(leftflanks);
      }

      lastcoord = coordstart;
      for (i = 0; i < nmatches; i++) {
	debug(printf("\nindex = %d\n",matches[i]));
	print_interval(&lastcoord,/*total*/0,divstring,coordstart,coordend,matches[i],iit,ndivs,fieldint);
      }
      
      if (nflanking > 0) {
	printf("====================\n");
	for (i = 0; i < nrightflanks; i++) {
	  debug(printf("\nright index = %d\n",rightflanks[i]));
	  print_interval(&lastcoord,/*total*/0,divstring,coordstart,coordend,rightflanks[i],iit,ndivs,fieldint);
	}
	FREE(rightflanks);
      }
    }

    if (divstring != NULL) {
      FREE(divstring);
    }
    FREE(matches);
    IIT_free(&iit);
    return 0;
  }
}

