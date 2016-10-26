static char rcsid[] = "$Id: iit_dump.c 99737 2013-06-27 19:33:03Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>		/* For getopt */
#endif
#include <ctype.h>		/* For isdigit */

#include "bool.h"
#include "mem.h"
#include "iit-read-univ.h"
#include "iit-read.h"
#include "list.h"
#include "getopt.h"

/************************************************************************
 *   Program options
 ************************************************************************/

static bool debugp = false;
static bool sortp = false;
static bool tagsp = false;
static bool countsp = false;
static bool integratep = false;
static bool annotationonlyp = false;

static struct option long_options[] = {
  /* Input options */
  {"debug", no_argument, 0, '9'}, /* debugp */
  {"sort", no_argument, 0, 'S'},	/* sortp */
  {"tags", no_argument, 0, 'T'}, /* tagsp */
  {"counts", no_argument, 0, 'C'}, /* countsp */
  {"integrate", no_argument, 0, 'I'}, /* integratep */
  {"annotonly", no_argument, 0, 'A'},	/* annotationonlyp */

  /* Help options */
  {"version", no_argument, 0, 'V'}, /* print_program_version */
  {"help", no_argument, 0, '?'}, /* print_program_usage */
  {0, 0, 0, 0}
};

static void
print_program_version () {
  fprintf(stdout,"\n");
  fprintf(stdout,"iit_dump: debugging utility for Interval Index Trees\n");
  fprintf(stdout,"Part of GMAP package, version %s\n",PACKAGE_VERSION);
  fprintf(stdout,"Thomas D. Wu, Genentech, Inc.\n");
  fprintf(stdout,"Contact: twu@gene.com\n");
  fprintf(stdout,"\n");
  return;
}

static void
print_program_usage () {
  fprintf(stdout,"\
Usage: iit_dump [OPTIONS...] iitfile\n\
\n\
Options\n\
  -S, --sort              Sort results by coordinates\n\
  -T, --tags              Show tags present in iit file\n\
  -C, --counts            Show counts for every boundary in iit file\n\
  -I, --integrate         Print intervals as integral output\n\
  -9, --debug             Provide debugging information\n\
  -A, --annotonly         Dump annotation lines only (no headers)\n\
\n\
  -V, --version           Show version\n\
  -?, --help              Show this help message\n\
\n\
The iit_dump program shows the entire contents of a given iit file.\n\
The default behavior is generate FASTA-type output, with both headers\n\
and annotations.  If only the annotations are desired, the -A flag\n\
may be used.  This flag may be useful for iit files created using the -G\n\
flag to iit_store, which stores the original gff3-formatted lines as\n\
the annotation.\n\
\n\
See also: iit_store, iit_get\n\
");

  return;
}

/*
static void
show_types (IIT_T iit) {
  List_T typelist, p;

  typelist = IIT_typelist(iit);
  for (p = typelist; p != NULL; p = List_next(p)) {
    printf("%s.\n",(char *) List_head(p));
  }
  return;
}
*/


/* Note that isnumber is a function in ctype.h on some systems */
static bool
isnumberp (long int *result, char *string) {
  char *p = string;

  *result = 0U;
  if (*p == '\0') {
    /* Empty string */
    return false;
  } else {
    while (*p != '\0') {
      if (*p == ',') {
	/* Skip commas */
      } else if (!isdigit((int) *p)) {
	return false;
      } else {
	*result = (*result) * 10 + (*p - '0');
      }
      p++;
    }
    return true;
  }
}


static void
print_runlengths (IIT_T iit, char *divstring) {
  long int *cum, level, n;
  unsigned int divlength, chrlength, low, high, lastpos, pos;
  int index, *matches;
  int nmatches, i;
  char *label;
  bool allocp;

  divlength = IIT_divlength(iit,divstring);
  cum = (long int *) CALLOC(divlength+1,sizeof(long int));

  matches = IIT_get(&nmatches,iit,divstring,/*x*/0,/*y*/-1U,/*sortp*/false);
  for (i = 0; i < nmatches; i++) {
    index = matches[i];
    label = IIT_label(iit,index,&allocp);
    if (isnumberp(&n,label) == false) {
      n = 1;
    }
    if (allocp == true) {
      FREE(label);
    }
    IIT_interval_bounds(&low,&high,&chrlength,iit,index,/*circular_typeint*/-1);
    cum[low] += n;
    cum[high] -= n;
    
  }
  FREE(matches);

  /* Print runlengths */
  lastpos = 1U;
  level = 0;
  for (pos = 1; pos <= divlength; pos++) {
    if (cum[pos] != 0) {
      /* printf("cum at pos %u is %d\n",pos,cum[pos]); */
      if (divstring[0] == '\0') {
	printf(">%ld %u..%u\n",level,lastpos,pos-1);
      } else {
	printf(">%ld %s:%u..%u\n",level,divstring,lastpos,pos-1);
      }
      level += cum[pos];
      lastpos = pos;
    }
  }

  /* Should print the final level of zero */
  if (level != 0) {
    fprintf(stderr,"Ended with a non-zero level\n");
    abort();
  }

  if (divstring[0] == '\0') {
    printf(">%ld %u..%u\n",level,lastpos,pos-1);
  } else {
    printf(">%ld %s:%u..%u\n",level,divstring,lastpos,pos-1);
  }

  FREE(cum);

  return;
}

int
main (int argc, char *argv[]) {
  char *iitfile;
  Univ_IIT_T chromosome_iit;
  IIT_T iit;
  int type;
  int divno;
  char *divstring;
  
  int opt;
  extern int optind;
  extern char *optarg;
  int long_option_index = 0;

  while ((opt = getopt_long(argc,argv,"9STCIA",long_options,&long_option_index)) != -1) {
    switch (opt) {
    case '9': debugp = true; break;
    case 'S': sortp = true; break;
    case 'T': tagsp = true; break;
    case 'C': countsp = true; break;
    case 'I': integratep = true; break;
    case 'A': annotationonlyp = true; break;

    case 'V': print_program_version(); exit(0);
    case '?': print_program_usage(); exit(0);
    default: exit(9);
    }
  }

  argc -= optind;
  argv += optind;

  iitfile = argv[0];

  if (IIT_universalp(iitfile,/*add_iit_p*/true) == true) {
    if (debugp == true) {
      Univ_IIT_debug(iitfile);
      return 0;

    } else if ((chromosome_iit = Univ_IIT_read(iitfile,/*readonlyp*/true,/*add_iit_p*/true)) == NULL) {
      fprintf(stderr,"Unable to open or parse IIT file %s\n",iitfile);
      exit(9);
    } else {
      Univ_IIT_dump(chromosome_iit);
    }
    Univ_IIT_free(&chromosome_iit);

  } else {
    if (debugp == true) {
      IIT_debug(iitfile);
      return 0;

    } else if ((iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
			       /*divstring*/NULL,/*add_iit_p*/true,/*labels_read_p*/true)) == NULL) {
      fprintf(stderr,"Unable to open or parse IIT file %s\n",iitfile);
      exit(9);
    }

    if (tagsp == true) {
      for (type = 1; type < IIT_ntypes(iit); type++) {
	printf("%s\n",IIT_typestring(iit,type));
      }
    } else if (countsp == true) {
      fprintf(stderr,"Flag -C not implemented\n");
#if 0
      IIT_dump_counts(iit,/*alphabetizep*/true);
#endif

    } else if (integratep == true) {
      for (divno = 1; divno < IIT_ndivs(iit); divno++) {
	divstring = IIT_divstring(iit,divno);
	print_runlengths(iit,divstring);
      }

    } else {
      IIT_dump(iit,annotationonlyp,sortp);
    }

    IIT_free(&iit);
  }

  return 0;
}
