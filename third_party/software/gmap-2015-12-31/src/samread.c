static char rcsid[] = "$Id: samread.c 154570 2014-12-03 22:06:33Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For strcpy */
#include <strings.h>		/* For rindex */
#include <ctype.h>

#include "except.h"
#include "mem.h"
#include "assert.h"
#include "bool.h"

#include "samread.h"
#include "samflags.h"
#include "chrnum.h"


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


static int
cigar_string_readlength (int *hardclip_low, int *hardclip_high, char *cigar) {
  int readlength;
  unsigned int npos;
  char *p, type;
  bool firstp = true;

  *hardclip_low = *hardclip_high = 0;
  if (cigar[0] == '*') {
    return 0;
  } else {
    readlength = 0;

    p = cigar;
    while (*p != '\0') {
      if (sscanf(p,"%u",&npos) != 1) {
	fprintf(stderr,"Unable to parse cigar %s.  No number in %s\n",cigar,p);
	exit(9);
      }

      while (*p != '\0' && isdigit(*p)) {
	p++;
      }
      if (*p == '\0') {
	fprintf(stderr,"Unable to parse cigar %s.  No letter after number %u\n",cigar,npos);
	exit(9);
      } else {
	type = *p++;
      }

      if (type == 'S' || type == 'M' || type == 'I') {
	readlength += (int) npos;

      } else if (type == 'H') {
	if (firstp == true) {
	  *hardclip_low = npos;
	} else {
	  *hardclip_high = npos;
	}
	readlength += (int) npos;

      } else if (type == 'D' || type == 'N') {
	/* Ignore */
      } else if (type == 'P') {
	/* Ignore */
      } else {
	fprintf(stderr,"samread.c cannot parse type %c\n",type);
	exit(9);
      }

      firstp = false;
    }
  }

  return readlength;
}


static int
cigar_string_initial_softclip (char *cigar) {
  unsigned int npos;
  char *p, type;

  if (cigar[0] == '*') {
    return 0;
  } else {
    p = cigar;
    while (*p != '\0') {
      if (sscanf(p,"%u",&npos) != 1) {
	fprintf(stderr,"Unable to parse cigar %s.  No number in %s\n",cigar,p);
	exit(9);
      }

      while (*p != '\0' && isdigit(*p)) {
	p++;
      }
      if (*p == '\0') {
	fprintf(stderr,"Unable to parse cigar %s.  No letter after number %u\n",cigar,npos);
	exit(9);
      } else {
	type = *p++;
      }

      if (type != 'S') {
	return 0;
      } else {
	return (int) npos;
      }
    }

    return 0;
  }
}


#if 0
char *
Samread_get_acc (unsigned int *flag, char *line) {
  char *acc, *p;
  int length;

  p = line;
  while (*p != '\0' && *p != '\t') p++;
  length = (p - line)/sizeof(char);
  acc = (char *) CALLOC(length+1,sizeof(char));
  strncpy(acc,line,length);

  if (*p == '\0') {
    fprintf(stderr,"Can't parse flag part of %s\n",line);
    abort();
  } else {
    p++;			/* Skip over tab */
  }
  *flag = strtoul(p,NULL,10);

  return acc;
}
#endif


/* Takes advantage of the fact that sam_sort knows the linelength */
char *
Samread_get_acc_fromfile (int *acclength, FILE *fp, int linelength) {
  char *acc, *p;

  p = acc = MALLOC((linelength + 1) * sizeof(char));
  while ((*p++ = fgetc(fp)) != '\t') ;
  *--p = '\0';			/* terminating char */
  *acclength = (p - acc)/sizeof(char);

  return acc;
}


/* Called just after we read in '\t', so should start at a field */
static SAM_split_output_type
parse_XO_fromfile (FILE *fp) {
  char c = 1, c0, c1;
  char abbrev0, abbrev1;

  c0 = fgetc(fp);
  c1 = fgetc(fp);
  if (c0 == 'X' && c1 == 'O') {
    fgetc(fp);		/* : */
    fgetc(fp);		/* type */
    fgetc(fp);		/* : */
    abbrev0 = fgetc(fp);
    abbrev1 = fgetc(fp);
    switch (abbrev0) {
    case 'N':
      if (abbrev1 == 'M') {
	return OUTPUT_NM;
      } else {
	fprintf(stderr,"Unexpected output type %c%c\n",abbrev0,abbrev1);
	return OUTPUT_NONE;
      }
    case 'C':
      switch (abbrev1) {
      case 'U': return OUTPUT_CU;
      case 'C': return OUTPUT_CC;
      case 'T': return OUTPUT_CT;
      case 'M': return OUTPUT_CM;
      case 'X': return OUTPUT_CX;
      default: fprintf(stderr,"Unexpected output type %c%c\n",abbrev0,abbrev1); return OUTPUT_NONE;
      }
    case 'H':
      switch (abbrev1) {
      case 'U': return OUTPUT_HU;
      case 'C': return OUTPUT_HC;
      case 'T': return OUTPUT_HT;
      case 'M': return OUTPUT_HM;
      case 'X': return OUTPUT_HX;
      default: fprintf(stderr,"Unexpected output type %c%c\n",abbrev0,abbrev1); return OUTPUT_NONE;
      }
    case 'U':
      switch (abbrev1) {
      case 'U': return OUTPUT_UU;
      case 'C': return OUTPUT_UC;
      case 'T': return OUTPUT_UT;
      case 'M': return OUTPUT_UM;
      case 'X': return OUTPUT_UX;
      default: fprintf(stderr,"Unexpected output type %c%c\n",abbrev0,abbrev1); return OUTPUT_NONE;
      }
    case 'P':
      switch (abbrev1) {
      case 'C': return OUTPUT_PC;
      case 'I': return OUTPUT_PI;
      case 'S': return OUTPUT_PS;
      case 'L': return OUTPUT_PL;
      case 'M': return OUTPUT_PM;
      case 'X': return OUTPUT_PX;
      default: fprintf(stderr,"Unexpected output type %c%c\n",abbrev0,abbrev1); return OUTPUT_NONE;
      }
    default: fprintf(stderr,"Unexpected output type %c%c\n",abbrev0,abbrev1); return OUTPUT_NONE;
    }
  }

  while (c != '\0') {
    while ((c = fgetc(fp)) != '\0' && c != '\t') ;
    if (c == '\0') {
      return OUTPUT_NONE;
    } else {
      c0 = fgetc(fp);
      c1 = fgetc(fp);
      if (c0 == 'X' && c1 == 'O') {
	fgetc(fp);		/* : */
	fgetc(fp);		/* type */
	fgetc(fp);		/* : */
	abbrev0 = fgetc(fp);
	abbrev1 = fgetc(fp);
	switch (abbrev0) {
	case 'N':
	  if (abbrev1 == 'M') {
	    return OUTPUT_NM;
	  } else {
	    fprintf(stderr,"Unexpected output type %c%c\n",abbrev0,abbrev1);
	    return OUTPUT_NONE;
	  }
	case 'C':
	  switch (abbrev1) {
	  case 'U': return OUTPUT_CU;
	  case 'C': return OUTPUT_CC;
	  case 'T': return OUTPUT_CT;
	  case 'M': return OUTPUT_CM;
	  case 'X': return OUTPUT_CX;
	  default: fprintf(stderr,"Unexpected output type %c%c\n",abbrev0,abbrev1); return OUTPUT_NONE;
	  }
	case 'H':
	  switch (abbrev1) {
	  case 'U': return OUTPUT_HU;
	  case 'C': return OUTPUT_HC;
	  case 'T': return OUTPUT_HT;
	  case 'M': return OUTPUT_HM;
	  case 'X': return OUTPUT_HX;
	  default: fprintf(stderr,"Unexpected output type %c%c\n",abbrev0,abbrev1); return OUTPUT_NONE;
	  }
	case 'U':
	  switch (abbrev1) {
	  case 'U': return OUTPUT_UU;
	  case 'C': return OUTPUT_UC;
	  case 'T': return OUTPUT_UT;
	  case 'M': return OUTPUT_UM;
	  case 'X': return OUTPUT_UX;
	  default: fprintf(stderr,"Unexpected output type %c%c\n",abbrev0,abbrev1); return OUTPUT_NONE;
	  }
	case 'P':
	  switch (abbrev1) {
	  case 'C': return OUTPUT_PC;
	  case 'I': return OUTPUT_PI;
	  case 'S': return OUTPUT_PS;
	  case 'L': return OUTPUT_PL;
	  case 'M': return OUTPUT_PM;
	  case 'X': return OUTPUT_PX;
	  default: fprintf(stderr,"Unexpected output type %c%c\n",abbrev0,abbrev1); return OUTPUT_NONE;
	  }
	default: fprintf(stderr,"Unexpected output type %c%c\n",abbrev0,abbrev1); return OUTPUT_NONE;
	}
      }
    }
  }

  return OUTPUT_NONE;
}


#define HITI_MAXDIGITS 10

/* Called just after we read in '\t', so should start at a field */
static SAM_split_output_type
parse_XO_and_HI_fromfile (char **hiti, FILE *fp) {
  SAM_split_output_type split_output = OUTPUT_NONE;
  char *p, c = 1, c0, c1;
  char abbrev0, abbrev1;

  *hiti = MALLOC((HITI_MAXDIGITS + 1) * sizeof(char));

  c0 = fgetc(fp);
  c1 = fgetc(fp);
  if (c0 == 'H' && c1 == 'I') {
    fgetc(fp);		/* : */
    fgetc(fp);		/* type */
    fgetc(fp);		/* : */

    p = *hiti;
    while ((c = *p++ = fgetc(fp)) != '\0' && c != '\t') ;
    *--p = '\0';			/* terminating char */

  } else if (c0 == 'X' && c1 == 'O') {
    fgetc(fp);		/* : */
    fgetc(fp);		/* type */
    fgetc(fp);		/* : */
    abbrev0 = fgetc(fp);
    abbrev1 = fgetc(fp);
    switch (abbrev0) {
    case 'N':
      if (abbrev1 == 'M') {
	split_output = OUTPUT_NM;
      } else {
	fprintf(stderr,"Unexpected output type %c%c\n",abbrev0,abbrev1);
	split_output = OUTPUT_NONE;
      }
    case 'C':
      switch (abbrev1) {
      case 'U': split_output = OUTPUT_CU; break;
      case 'C': split_output = OUTPUT_CC; break;
      case 'T': split_output = OUTPUT_CT; break;
      case 'M': split_output = OUTPUT_CM; break;
      case 'X': split_output = OUTPUT_CX; break;
      default:
	fprintf(stderr,"Unexpected output type %c%c\n",abbrev0,abbrev1);
	split_output = OUTPUT_NONE;
      }
    case 'H':
      switch (abbrev1) {
      case 'U': split_output = OUTPUT_HU; break;
      case 'C': split_output = OUTPUT_HC; break;
      case 'T': split_output = OUTPUT_HT; break;
      case 'M': split_output = OUTPUT_HM; break;
      case 'X': split_output = OUTPUT_HX; break;
      default:
	fprintf(stderr,"Unexpected output type %c%c\n",abbrev0,abbrev1);
	split_output = OUTPUT_NONE;
      }
    case 'U':
      switch (abbrev1) {
      case 'U': split_output = OUTPUT_UU; break;
      case 'C': split_output = OUTPUT_UC; break;
      case 'T': split_output = OUTPUT_UT; break;
      case 'M': split_output = OUTPUT_UM; break;
      case 'X': split_output = OUTPUT_UX; break;
      default:
	fprintf(stderr,"Unexpected output type %c%c\n",abbrev0,abbrev1);
	split_output = OUTPUT_NONE;
      }
    case 'P':
      switch (abbrev1) {
      case 'C': split_output = OUTPUT_PC; break;
      case 'I': split_output = OUTPUT_PI; break;
      case 'S': split_output = OUTPUT_PS; break;
      case 'L': split_output = OUTPUT_PL; break;
      case 'M': split_output = OUTPUT_PM; break;
      case 'X': split_output = OUTPUT_PX; break;
      default:
	fprintf(stderr,"Unexpected output type %c%c\n",abbrev0,abbrev1);
	split_output = OUTPUT_NONE;
      }
    default:
      fprintf(stderr,"Unexpected output type %c%c\n",abbrev0,abbrev1);
      split_output = OUTPUT_NONE;
    }
  }

  while (c != '\0') {
    while ((c = fgetc(fp)) != '\0' && c != '\t') ;
    if (c == '\0') {
      return split_output;
    } else {
      c0 = fgetc(fp);
      c1 = fgetc(fp);
      if (c0 == 'H' && c1 == 'I') {
	fgetc(fp);		/* : */
	fgetc(fp);		/* type */
	fgetc(fp);		/* : */

	p = *hiti;
	while ((c = *p++ = fgetc(fp)) != '\0' && c != '\t') ;
	*--p = '\0';			/* terminating char */

      } else if (c0 == 'X' && c1 == 'O') {
	fgetc(fp);		/* : */
	fgetc(fp);		/* type */
	fgetc(fp);		/* : */
	abbrev0 = fgetc(fp);
	abbrev1 = fgetc(fp);
	switch (abbrev0) {
	case 'N':
	  if (abbrev1 == 'M') {
	    split_output = OUTPUT_NM;
	  } else {
	    fprintf(stderr,"Unexpected output type %c%c\n",abbrev0,abbrev1);
	    split_output = OUTPUT_NONE;
	  }
	case 'C':
	  switch (abbrev1) {
	  case 'U': split_output = OUTPUT_CU; break;
	  case 'C': split_output = OUTPUT_CC; break;
	  case 'T': split_output = OUTPUT_CT; break;
	  case 'M': split_output = OUTPUT_CM; break;
	  case 'X': split_output = OUTPUT_CX; break;
	  default:
	    fprintf(stderr,"Unexpected output type %c%c\n",abbrev0,abbrev1);
	    split_output = OUTPUT_NONE;
	  }
	case 'H':
	  switch (abbrev1) {
	  case 'U': split_output = OUTPUT_HU; break;
	  case 'C': split_output = OUTPUT_HC; break;
	  case 'T': split_output = OUTPUT_HT; break;
	  case 'M': split_output = OUTPUT_HM; break;
	  case 'X': split_output = OUTPUT_HX; break;
	  default:
	    fprintf(stderr,"Unexpected output type %c%c\n",abbrev0,abbrev1);
	    split_output = OUTPUT_NONE;
	  }
	case 'U':
	  switch (abbrev1) {
	  case 'U': split_output = OUTPUT_UU; break;
	  case 'C': split_output = OUTPUT_UC; break;
	  case 'T': split_output = OUTPUT_UT; break;
	  case 'M': split_output = OUTPUT_UM; break;
	  case 'X': split_output = OUTPUT_UX; break;
	  default:
	    fprintf(stderr,"Unexpected output type %c%c\n",abbrev0,abbrev1);
	    split_output = OUTPUT_NONE;
	  }
	case 'P':
	  switch (abbrev1) {
	  case 'C': split_output = OUTPUT_PC; break;
	  case 'I': split_output = OUTPUT_PI; break;
	  case 'S': split_output = OUTPUT_PS; break;
	  case 'L': split_output = OUTPUT_PL; break;
	  case 'M': split_output = OUTPUT_PM; break;
	  case 'X': split_output = OUTPUT_PX; break;
	  default:
	    fprintf(stderr,"Unexpected output type %c%c\n",abbrev0,abbrev1);
	    split_output = OUTPUT_NONE;
	  }
	default:
	  fprintf(stderr,"Unexpected output type %c%c\n",abbrev0,abbrev1);
	  split_output = OUTPUT_NONE;
	}
      }
    }
  }

  return split_output;
}



/* Main parser for processing with dups */
char *
Samread_parse_acc_and_softclip_fromfile (int *acclength, unsigned int *flag, SAM_split_output_type *split_output,
					 char **hiti, Univcoord_T *genomicpos, int *initial_softclip, bool *query_lowp,
					 FILE *fp, Univ_IIT_T chromosome_iit, Univcoord_T *chroffsets, int linelength) {
  char *acc, *p;
  char *substring;
  Chrnum_T chrnum, mate_chrnum;
  Chrpos_T chrpos, mate_chrpos;
  Univcoord_T mate_genomicpos;

  substring = MALLOCA((linelength + 1) * sizeof(char));

  /* 1. QNAME */
  p = acc = MALLOC((linelength + 1) * sizeof(char));
  while ((*p++ = fgetc(fp)) != '\t') ;
  *--p = '\0';			/* terminating char */
  *acclength = (p - acc)/sizeof(char);

  /* 2. FLAG.  Skip */
  p = substring;
  while ((*p++ = fgetc(fp)) != '\t') ;
  *--p = '\0';

  if (sscanf(substring,"%u",&(*flag)) != 1) {
    fprintf(stderr,"Unable to find flag in %s\n",substring);
    abort();
  } else {
    debug(printf("  flag = %u\n",*flag));
  }

  /* 3. RNAME: chr.  Skip */
  p = substring;
  while ((*p++ = fgetc(fp)) != '\t') ;
  *--p = '\0';

  if (!strcmp(substring,"*")) {
    *genomicpos = 0;
  } else if ((chrnum = Univ_IIT_find_one(chromosome_iit,substring)) < 0) {
    fprintf(stderr,"Cannot find chromosome %s in chromosome IIT file\n",substring);
    exit(9);
  } else {
    *genomicpos = chroffsets[chrnum - 1];
  }

  /* 4. POS: chrpos */
  p = substring;
  while ((*p++ = fgetc(fp)) != '\t') ;
  *--p = '\0';

  if (sscanf(substring,"%u",&chrpos) != 1) {
    fprintf(stderr,"Unable to find chrpos in %s\n",substring);
    abort();
  } else {
    *genomicpos += chrpos;
  }

  /* 5. MAPQ: Mapping quality.  Skip */
  while (fgetc(fp) != '\t') ;

  /* 6. CIGAR.  Parse for initial_softclip */
  p = substring;
  while ((*p++ = fgetc(fp)) != '\t') ;
  *--p = '\0';

  *initial_softclip = cigar_string_initial_softclip(substring);


  /* 7. MRNM: Mate chr */
  p = substring;
  while ((*p++ = fgetc(fp)) != '\t') ;
  *--p = '\0';
  
  if (!strcmp(substring,"*")) {
    mate_genomicpos = 0;
  } else if (!strcmp(substring,"=")) {
    mate_chrnum = chrnum;
    mate_genomicpos = chroffsets[mate_chrnum - 1];
  } else if ((mate_chrnum = Univ_IIT_find_one(chromosome_iit,substring)) < 0) {
    fprintf(stderr,"Cannot find chromosome %s in chromosome IIT file\n",substring);
    exit(9);
  } else {
    mate_genomicpos = chroffsets[mate_chrnum - 1];
  }

  /* 8. MPOS: Mate chrpos */
  p = substring;
  while ((*p++ = fgetc(fp)) != '\t') ;
  *--p = '\0';

  if (sscanf(substring,"%d",&mate_chrpos) != 1) {
    fprintf(stderr,"Unable to find mate_chrpos_low in %s\n",substring);
    abort();
  } else {
    mate_genomicpos += mate_chrpos;
  }


  /* Determine if the query is low */
  if (*genomicpos == 0) {
    if (*flag & PAIRED_READ) {
      if (mate_genomicpos != 0) {
	/* Mate will be mapped, so this is the high end */
	*query_lowp = false;
      } else if (*flag & FIRST_READ_P) {
	*query_lowp = true;
      } else {
	*query_lowp = false;
      }
    } else {
      /* Single-end, so this is the low end */
      *query_lowp = true;
    }

  } else if (mate_chrpos == 0) {
    *query_lowp = true;

  } else if (*genomicpos < mate_genomicpos) {
    *query_lowp = true;

  } else {
    *query_lowp = false;

  }
    
  FREEA(substring);

  /* 9. ISIZE: Insert size.  Skip. */
  while (fgetc(fp) != '\t') ;

  /* 10. SEQ: queryseq.  Skip. */
  while (fgetc(fp) != '\t') ;

  /* 11. QUAL: quality scores.  Skip. */
  while (fgetc(fp) != '\t') ;

  *split_output = parse_XO_and_HI_fromfile(&(*hiti),fp);

  return acc;
}



/* ILLUMINA-A1CCE9_0004:1:1:1103:6310#0	0	20	33639850	255	55M21S	*	0	0	AAAAATTGTATACCGCAGATTCAGGCATGGATTCCGTGAAGGAACAACACCTAAANCCAAAGNTCGGAAGANCGGN	CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCDCCCCBCCBDBDCBCCDCDC@CCC&AAAA################	NM:i:2 */

#if 0
char *
Samread_parse_line (char **acc, unsigned int *flag, int *mapq, char **chr, Chrpos_T *chrpos, char **cigar,
		    char **mate_chr, Chrpos_T *mate_chrpos_low,
		    int *readlength, char **read, char **quality_string, char *line) {
  char *p, *q;
  int length, i;

  debug(printf("Entering Samread_parse_line with %s\n",line));

  p = line;
  while (!isspace(*p)) p++;
  length = (p - line)/sizeof(char);
  *acc = (char *) CALLOC(length+1,sizeof(char));
  strncpy(*acc,line,length);

  if (*p != '\0') {		/* Skip over tab */
    p++;
  }

  if (sscanf(p,"%u",&(*flag)) != 1) {
    fprintf(stderr,"Unable to find flag in %s\n",p);
    abort();
  } else {
    debug(printf("  flag = %u\n",*flag));
  }

  while (!isspace(*p)) p++;	/* Skip over flag */
  if (*p == '\0') {
    fprintf(stderr,"Can't parse chr part of %s\n",line);
    abort();
  } else {
    p++;			/* Skip over tab */
  }
  q = p;
  while (!isspace(*q)) q++;
  length = (q - p)/sizeof(char);
  *chr = (char *) CALLOC(length+1,sizeof(char));
  strncpy(*chr,p,length);
  debug(printf("  chr = %s\n",*chr));
  if (*q != '\0') {
    q++;
  }


  p = q;
  if (sscanf(p,"%u",&(*chrpos)) != 1) {
    fprintf(stderr,"Unable to find chrpos in %s\n",p);
    abort();
  } else {
    debug(printf("  chrpos = %u\n",*chrpos));
  }


  while (!isspace(*p)) p++;	/* Skip over chrpos */
  if (*p == '\0') {
    fprintf(stderr,"Can't parse chrpos part of %s\n",line);
    abort();
  } else {
    p++;			/* Skip over tab */
  }

  /* Read mapping quality */
  if (sscanf(p,"%d",&(*mapq)) != 1) {
    fprintf(stderr,"Unable to find mapq in %s\n",p);
    abort();
  } else {
    debug(printf("  mapq = %d\n",*mapq));
  }

  /* Skip past mapping quality */
  while (!isspace(*p)) p++;


  if (*p == '\0') {
    fprintf(stderr,"Can't parse cigar part of %s\n",line);
    abort();
  } else {
    p++;			/* Skip over tab */
  }
  q = p;
  while (!isspace(*q)) q++;
  length = (q - p)/sizeof(char);
  *cigar = (char *) CALLOC(length+1,sizeof(char));
  strncpy(*cigar,p,length);
  debug(printf("  cigar = %s\n",*cigar));
  

  /* mate chr */
  p = q;
  if (*p != '\0') {
    p++;			/* Should be a tab */
  }
  q = p;
  while (!isspace(*q)) q++;
  length = (q - p)/sizeof(char);
  *mate_chr = (char *) CALLOC(length+1,sizeof(char));
  strncpy(*mate_chr,p,length);
  debug(printf("  mate_chr = %s\n",*mate_chr));
  if (*q == '\0') {
    fprintf(stderr,"Can't parse mate chr part of %s\n",line);
    abort();
  } else {
    q++;
  }

  /* mate chrpos low */
  p = q;
  if (sscanf(p,"%u",&(*mate_chrpos_low)) != 1) {
    fprintf(stderr,"Unable to find mate_chrpos_low in %s\n",p);
    abort();
  } else {
    debug(printf("  mate_chrpos_low = %u\n",*mate_chrpos_low));
  }

  while (!isspace(*p)) p++;	/* Skip over mate_chrpos */
  if (*p == '\0') {
    fprintf(stderr,"Can't parse mate chrpos part of %s\n",line);
    abort();
  } else {
    p++;			/* Skip over tab */
  }


  /* Skip over insert size */
  while (!isspace(*p)) p++;
  if (*p == '\0') {
    fprintf(stderr,"Can't parse mate chrpos part of %s\n",line);
    abort();
  } else {
    p++;
  }


  q = p;
  while (!isspace(*q)) q++;
  *readlength = (q - p)/sizeof(char);
  if (*q == '\t') q++;
  debug(printf("  readlength = %d\n",*readlength));

  *read = (char *) CALLOC((*readlength)+1,sizeof(char));
  strncpy(*read,p,*readlength);
  debug(printf("  read = %s\n",*read));

  p = q;
  while (!isspace(*q)) q++;
  length = (q - p)/sizeof(char);
  *quality_string = (char *) CALLOC((*readlength)+1,sizeof(char));
  if (length == *readlength) {
    strncpy(*quality_string,p,length);
  } else {
    for (i = 0; i < *readlength; i++) {
      (*quality_string)[i] = ' ';
    }
  }

  if (*q == '\t') q++;

  return q;
}
#endif

#if 0
/* Leaves fp at start of auxinfo */
void
Samread_parse_line_fromfile (FILE *fp, char **acc, unsigned int *flag, int *mapq, char **chr, Chrpos_T *chrpos, char **cigar,
			     char **mate_chr, Chrpos_T *mate_chrpos_low,
			     int *readlength, char **read, char **quality_string, int linelength) {
  char *p;
  int length, i;
  char *substring;

  debug(printf("Entering Samread_parse_line_fromfile\n"));

  substring = MALLOCA((linelength + 1) * sizeof(char));

  /* 1. QNAME */
  p = *acc = (char *) MALLOC((linelength + 1) * sizeof(char));
  while ((*p++ = fgetc(fp)) != '\t') ;
  *--p = '\0';

  /* 2. FLAG */
  p = substring;
  while ((*p++ = fgetc(fp)) != '\t') ;
  *--p = '\0';

  if (sscanf(substring,"%u",&(*flag)) != 1) {
    fprintf(stderr,"Unable to find flag in %s\n",substring);
    abort();
  } else {
    debug(printf("  flag = %u\n",*flag));
  }

  /* 3. RNAME: chr */
  p = *chr = (char *) MALLOC((linelength + 1) * sizeof(char));
  while ((*p++ = fgetc(fp)) != '\t') ;
  *--p = '\0';

  /* 4. POS: chrpos */
  p = substring;
  while ((*p++ = fgetc(fp)) != '\t') ;
  *--p = '\0';

  if (sscanf(substring,"%u",&(*chrpos)) != 1) {
    fprintf(stderr,"Unable to find chrpos in %s\n",substring);
    abort();
  } else {
    debug(printf("  chrpos = %u\n",*chrpos));
  }

  /* 5. MAPQ: Mapping quality */
  p = substring;
  while ((*p++ = fgetc(fp)) != '\t') ;
  *--p = '\0';

  if (sscanf(substring,"%d",&(*mapq)) != 1) {
    fprintf(stderr,"Unable to find mapq in %s\n",substring);
    abort();
  } else {
    debug(printf("  mapq = %d\n",*mapq));
  }

  /* 5. CIGAR */
  p = *cigar = (char *) MALLOC((linelength + 1) * sizeof(char));
  while ((*p++ = fgetc(fp)) != '\t') ;
  *--p = '\0';

  /* 7. MRNM: Mate chr */
  p = *mate_chr = (char *) MALLOC((linelength + 1) * sizeof(char));
  while ((*p++ = fgetc(fp)) != '\t') ;
  *--p = '\0';

  /* 8. MPOS: Mate chrpos */
  p = substring;
  while ((*p++ = fgetc(fp)) != '\t') ;
  *--p = '\0';

  if (sscanf(substring,"%d",&(*mate_chrpos_low)) != 1) {
    fprintf(stderr,"Unable to find mate_chrpos_low in %s\n",substring);
    abort();
  } else {
    debug(printf("  mate_chrpos_low = %u\n",*mate_chrpos_low));
  }

  /* 9. ISIZE: Insert size.  Skip. */
  while (fgetc(fp) != '\t') ;

  /* 10. SEQ: queryseq */
  p = *read = (char *) MALLOC((linelength + 1) * sizeof(char));
  while ((*p++ = fgetc(fp)) != '\t') ;
  *--p = '\0';
  *readlength = (p - *read)/sizeof(char);

  /* 11. QUAL: quality scores */
  p = *quality_string = (char *) MALLOC((*readlength + 1) * sizeof(char));
  while ((*p++ = fgetc(fp)) != '\t') ;
  *--p = '\0';
  length = (p - *quality_string)/sizeof(char);

  if (length != *readlength) {
    for (i = 0; i < *readlength; i++) {
      (*quality_string)[i] = ' ';
    }
    *quality_string[i] = '\0';
  }

  FREEA(substring);

  return;
}
#endif


int
Samread_parse_linelen_fromfile (FILE *fp) {
  int linelen = 0;
  
  /* 1. QNAME: Skip */
  while (!feof(fp) && fgetc(fp) != '\t') {
    linelen++;
  }
  linelen++;

  if (feof(fp)) {
    return -1;
  }

  while (!feof(fp) && fgetc(fp) != '\n') {
    linelen++;
  }
  linelen++;

  return linelen;
}


/* Main parser for processing without dups */
Univcoord_T
Samread_parse_genomicpos_fromfile (FILE *fp, unsigned int *flag, SAM_split_output_type *split_output,
				   Univ_IIT_T chromosome_iit, Univcoord_T *chroffsets, int linelength) {
  Univcoord_T genomicpos;
  Chrnum_T chrnum;
  Chrpos_T chrpos;
  char *substring, *p;
  
  substring = MALLOCA((linelength + 1) * sizeof(char));

  /* 1. QNAME: Skip */
  while (fgetc(fp) != '\t') ;

  /* 2. FLAG */
  p = substring;
  while ((*p++ = fgetc(fp)) != '\t') ;
  *--p = '\0';

  if (sscanf(substring,"%u",&(*flag)) != 1) {
    fprintf(stderr,"Unable to find flag in %s\n",substring);
    abort();
  } else {
    debug(printf("  flag = %u\n",*flag));
  }

  /* 3. RNAME: chr */
  p = substring;
  while ((*p++ = fgetc(fp)) != '\t') ;
  *--p = '\0';

  if (!strcmp(substring,"*")) {
    genomicpos = 0;
  } else if ((chrnum = Univ_IIT_find_one(chromosome_iit,substring)) < 0) {
    fprintf(stderr,"Cannot find chromosome %s in chromosome IIT file\n",substring);
    exit(9);
  } else {
    genomicpos = chroffsets[chrnum - 1];
  }

  /* 4. POS: chrpos */
  p = substring;
  while ((*p++ = fgetc(fp)) != '\t') ;
  *--p = '\0';

  if (sscanf(substring,"%u",&chrpos) != 1) {
    fprintf(stderr,"Unable to find chrpos in %s\n",substring);
    abort();
  } else {
    genomicpos += chrpos;
  }
  
  FREEA(substring);

  /* 5. MAPQ: Mapping quality.  Skip */
  while (fgetc(fp) != '\t') ;

  /* 6. CIGAR.  Skip */
  while (fgetc(fp) != '\t') ;

  /* 7. MRNM: Mate chr.  Skip */
  while (fgetc(fp) != '\t') ;

  /* 8. MPOS: Mate chrpos.  Skip */
  while (fgetc(fp) != '\t') ;

  /* 9. ISIZE: Insert size.  Skip. */
  while (fgetc(fp) != '\t') ;

  /* 10. SEQ: queryseq.  Skip. */
  while (fgetc(fp) != '\t') ;

  /* 11. QUAL: quality scores.  Skip. */
  while (fgetc(fp) != '\t') ;

  *split_output = parse_XO_fromfile(fp);

  return genomicpos;
}


char *
Samread_parse_aux_fromfile (FILE *fp, char *auxfield, int linelength) {
  char *value, *p, c = 1, c0, c1;

  while (c != '\0') {
    while ((c = fgetc(fp)) != '\0' && c != '\t') ;
    if (c == '\0') {
      return (char *) NULL;
    } else {
      c0 = fgetc(fp);
      c1 = fgetc(fp);
      if (c0 == auxfield[0] && c1 == auxfield[1]) {
	fgetc(fp);		/* : */
	fgetc(fp);		/* type */
	fgetc(fp);		/* : */
	p = value = MALLOC((linelength+1) * sizeof(char));
	while ((c = *p++ = fgetc(fp)) != '\0' && c != '\t') ;
	*--p = '\0';			/* terminating char */
	return value;
      }
    }
  }

  return (char *) NULL;
}



#if 0
/* Returns only fields needed by sam_sort */
void
Samread_parse_read_and_mateinfo_fromfile (FILE *fp, unsigned int *flag, char **mate_chr, Chrpos_T *mate_chrpos,
					  int *readlength, char **read, int linelength) {
  char *p, *q, c;
  char *substring, *clipped;
  int hardclip_low, hardclip_high, subseq_length;

  substring = MALLOCA((linelength + 1) * sizeof(char));

  /* 1. QNAME.  Skip */
  while (fgetc(fp) != '\t') ;

  /* 2. FLAG */
  p = substring;
  while ((*p++ = fgetc(fp)) != '\t') ;
  *--p = '\0';

  if (sscanf(substring,"%u",&(*flag)) != 1) {
    fprintf(stderr,"Unable to find flag in %s\n",substring);
    abort();
  } else {
    debug(printf("  flag = %u\n",*flag));
  }

  /* 3. RNAME: chr.  Skip */
  while (fgetc(fp) != '\t') ;

  /* 4. POS: chrpos.  Skip */
  while (fgetc(fp) != '\t') ;

  /* 5. MAPQ: Mapping quality.  Skip */
  while (fgetc(fp) != '\t') ;

  /* 6. CIGAR.  Parse for cigar_readlength. */
  p = substring;
  while ((*p++ = fgetc(fp)) != '\t') ;
  *--p = '\0';

  /* For a nomapper, this readlength is incorrect */
  *readlength = cigar_string_readlength(&hardclip_low,&hardclip_high,substring);


  /* 7. MRNM: Mate chr */
  p = *mate_chr = (char *) MALLOC((linelength + 1) * sizeof(char));
  while ((*p++ = fgetc(fp)) != '\t') ;
  *--p = '\0';

  /* 8. MPOS: Mate chrpos */
w  p = substring;
  while ((*p++ = fgetc(fp)) != '\t') ;
  *--p = '\0';

  if (sscanf(substring,"%d",&(*mate_chrpos)) != 1) {
    fprintf(stderr,"Unable to find mate_chrpos_low in %s\n",substring);
    abort();
  } else {
    debug(printf("  mate_chrpos_low = %u\n",*mate_chrpos));
  }

  /* 9. ISIZE: Insert size.  Skip. */
  while (fgetc(fp) != '\t') ;

  /* 10. SEQ: queryseq */
  *read = (char *) MALLOC((linelength + 1) * sizeof(char));

  p = &((*read)[hardclip_low]);
  subseq_length = 0;
  while ((*p++ = fgetc(fp)) != '\t') {
    subseq_length++;
  }

  if (*readlength == 0) {
    *readlength = subseq_length;
  }

  if (subseq_length + hardclip_low + hardclip_high != *readlength) {
    fprintf(stderr,"Cigar readlength %d is not consistent with subsequence length %d + hardclip_low %d + hardclip_high %d\n",
	    *readlength,subseq_length,hardclip_low,hardclip_high);
  }

  if (hardclip_low > 0 || hardclip_high > 0) {
    if ((clipped = Samread_parse_aux_fromfile(fp,/*auxfield*/"XH",linelength)) == NULL) {
      fprintf(stderr,"Hard-clipped read needs XH field from a recent GSNAP version\n");
      exit(9);

    } else {
      if (hardclip_low > 0) {
	p = &((*read)[0]);
	q = clipped;
	while ((c = *q++) != '\0') {
	  *p++ = c;
	}
      } else {
	p = &((*read)[(*readlength) - hardclip_high]);
	q = clipped;
	while ((c = *q++) != '\0') {
	  *p++ = c;
	}
      }
      FREE(clipped);

    }
  }
  (*read)[*readlength] = '\0';

  FREEA(substring);

  return;
}
#endif


void
Samread_parse_read_fromfile (FILE *fp, unsigned int *flag, int *readlength, char **read, int linelength) {
  char *p, *q, c;
  char *substring, *clipped;
  int hardclip_low, hardclip_high, subseq_length;

  substring = MALLOCA((linelength + 1) * sizeof(char));

  /* 1. QNAME.  Skip */
  while (fgetc(fp) != '\t') ;

  /* 2. FLAG.  Skip */
  p = substring;
  while ((*p++ = fgetc(fp)) != '\t') ;
  *--p = '\0';

  if (sscanf(substring,"%u",&(*flag)) != 1) {
    fprintf(stderr,"Unable to find flag in %s\n",substring);
    abort();
  } else {
    debug(printf("  flag = %u\n",*flag));
  }

  /* 3. RNAME: chr.  Skip */
  while (fgetc(fp) != '\t') ;

  /* 4. POS: chrpos.  Skip */
  while (fgetc(fp) != '\t') ;

  /* 5. MAPQ: Mapping quality.  Skip */
  while (fgetc(fp) != '\t') ;

  /* 6. CIGAR.  Parse for cigar_readlength. */
  p = substring;
  while ((*p++ = fgetc(fp)) != '\t') ;
  *--p = '\0';

  /* For a nomapper, this readlength is incorrect */
  *readlength = cigar_string_readlength(&hardclip_low,&hardclip_high,substring);


  /* 7. MRNM: Mate chr.  Skip */
  while (fgetc(fp) != '\t') ;

  /* 8. MPOS: Mate chrpos.  Skip */
  while (fgetc(fp) != '\t') ;

  /* 9. ISIZE: Insert size.  Skip. */
  while (fgetc(fp) != '\t') ;

  /* 10. SEQ: queryseq */
  *read = (char *) MALLOC((linelength + 1) * sizeof(char));

  p = &((*read)[hardclip_low]);
  subseq_length = 0;
  while ((*p++ = fgetc(fp)) != '\t') {
    subseq_length++;
  }

  if (*readlength == 0) {
    *readlength = subseq_length;
  }

  if (subseq_length + hardclip_low + hardclip_high != *readlength) {
    fprintf(stderr,"Cigar readlength %d is not consistent with subsequence length %d + hardclip_low %d + hardclip_high %d\n",
	    *readlength,subseq_length,hardclip_low,hardclip_high);
  }

  if (hardclip_low > 0 || hardclip_high > 0) {
    if ((clipped = Samread_parse_aux_fromfile(fp,/*auxfield*/"XH",linelength)) == NULL) {
      fprintf(stderr,"Hard-clipped read needs XH field from a recent GSNAP version\n");
      exit(9);

    } else {
      if (hardclip_low > 0) {
	p = &((*read)[0]);
	q = clipped;
	while ((c = *q++) != '\0') {
	  *p++ = c;
	}
      } else {
	p = &((*read)[(*readlength) - hardclip_high]);
	q = clipped;
	while ((c = *q++) != '\0') {
	  *p++ = c;
	}
      }
      FREE(clipped);

    }
  }

  (*read)[*readlength] = '\0';

  FREEA(substring);

  return;
}


Univcoord_T
Samread_parse_mate_genomicpos_fromfile (FILE *fp, Univ_IIT_T chromosome_iit, Univcoord_T *chroffsets, int linelength) {
  Univcoord_T mate_genomicpos;
  Chrpos_T mate_chrpos;
  Chrnum_T chrnum, mate_chrnum;
  char *p;
  char *substring;

  substring = MALLOCA((linelength + 1) * sizeof(char));

  /* 1. QNAME.  Skip */
  while (fgetc(fp) != '\t') ;

  /* 2. FLAG.  Skip */
  while (fgetc(fp) != '\t') ;

  /* 3. RNAME: chr */
  p = substring;
  while ((*p++ = fgetc(fp)) != '\t') ;
  *--p = '\0';

  if (!strcmp(substring,"*")) {
    chrnum = -1;
  } else if ((chrnum = Univ_IIT_find_one(chromosome_iit,substring)) < 0) {
    fprintf(stderr,"Cannot find chromosome %s in chromosome IIT file\n",substring);
  }

  /* 4. POS: chrpos.  Skip */
  while (fgetc(fp) != '\t') ;

  /* 5. MAPQ: Mapping quality.  Skip */
  while (fgetc(fp) != '\t') ;

  /* 6. CIGAR.  Parse for cigar_readlength.  Skip. */
  while (fgetc(fp) != '\t') ;

  /* 7. MRNM: Mate chr */
  p = substring;
  while ((*p++ = fgetc(fp)) != '\t') ;
  *--p = '\0';
  
  if (!strcmp(substring,"*")) {
    mate_genomicpos = 0;
  } else if (!strcmp(substring,"=")) {
    mate_chrnum = chrnum;
    mate_genomicpos = chroffsets[mate_chrnum - 1];
  } else if ((mate_chrnum = Univ_IIT_find_one(chromosome_iit,substring)) < 0) {
    fprintf(stderr,"Cannot find chromosome %s in chromosome IIT file\n",substring);
    exit(9);
  } else {
    mate_genomicpos = chroffsets[mate_chrnum - 1];
  }

  /* 8. MPOS: Mate chrpos */
  p = substring;
  while ((*p++ = fgetc(fp)) != '\t') ;
  *--p = '\0';

  if (sscanf(substring,"%d",&mate_chrpos) != 1) {
    fprintf(stderr,"Unable to find mate_chrpos_low in %s\n",substring);
    abort();
  } else {
    mate_genomicpos += mate_chrpos;
  }

  FREEA(substring);

  return mate_genomicpos;
}


#if 0
char *
Samread_chrinfo (Chrpos_T *chrpos, char **cigar, char *line) {
  char *chr;
  unsigned int flag;
  int mapq;

  char *p, *q;
  int length;

  debug(printf("Entering Samread_chrinfo with %s\n",line));

  p = line;
  while (!isspace(*p)) p++;
  length = (p - line)/sizeof(char);
#if 0
  *acc = (char *) CALLOC(length+1,sizeof(char));
  strncpy(*acc,line,length);
#endif

  if (*p != '\0') {		/* Skip over tab */
    p++;
  }

  if (sscanf(p,"%u",&flag) != 1) {
    fprintf(stderr,"Unable to find flag in %s\n",p);
    abort();
  } else {
    debug(printf("  flag = %u\n",*flag));
  }

  while (!isspace(*p)) p++;	/* Skip over flag */
  if (*p == '\0') {
    fprintf(stderr,"Can't parse chr part of %s\n",line);
    abort();
  } else {
    p++;			/* Skip over tab */
  }
  q = p;
  while (!isspace(*q)) q++;
  length = (q - p)/sizeof(char);
  chr = (char *) CALLOC(length+1,sizeof(char));
  strncpy(chr,p,length);
  debug(printf("  chr = %s\n",chr));
  if (*q != '\0') {
    q++;
  }


  p = q;
  if (sscanf(p,"%u",&(*chrpos)) != 1) {
    fprintf(stderr,"Unable to find chrpos in %s\n",p);
    abort();
  } else {
    debug(printf("  chrpos = %u\n",*chrpos));
  }

  while (!isspace(*p)) p++;	/* Skip over chrpos */
  if (*p == '\0') {
    fprintf(stderr,"Can't parse chrpos part of %s\n",line);
    abort();
  } else {
    p++;			/* Skip over tab */
  }

  /* Read mapping quality */
  if (sscanf(p,"%d",&mapq) != 1) {
    fprintf(stderr,"Unable to find mapq in %s\n",p);
    abort();
  } else {
    debug(printf("  mapq = %d\n",mapq));
  }

  /* Skip past mapping quality */
  while (!isspace(*p)) p++;


  if (*p == '\0') {
    fprintf(stderr,"Can't parse cigar part of %s\n",line);
    abort();
  } else {
    p++;			/* Skip over tab */
  }
  q = p;
  while (!isspace(*q)) q++;
  length = (q - p)/sizeof(char);
  *cigar = (char *) CALLOC(length+1,sizeof(char));
  strncpy(*cigar,p,length);
  debug(printf("  cigar = %s\n",*cigar));

  return chr;
}
#endif



#define CHUNK 1024

void
Samread_print_as_duplicate_fromfile (FILE *fp, int linelength) {
  int nread = 0;
  char buffer[CHUNK], *p, c;
  unsigned int flag;

  /* 1. QNAME */
  while ((c = fgetc(fp)) != '\t') {
    putchar(c);
    nread++;
  }
  putchar('\t');
  nread++;

  /* 2. FLAG */
  p = buffer;
  while ((*p++ = fgetc(fp)) != '\t') {
    nread++;
  }
  nread++;
  *--p = '\0';

  if (sscanf(buffer,"%u",&flag) != 1) {
    fprintf(stderr,"Unable to find flag in %s\n",buffer);
    abort();
  } else {
    printf("%u\t",flag | DUPLICATE_READ);
  }

  /* 3... Rest */
  linelength -= nread;

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


#if 0
char
Samread_splice_strand (char *auxinfo) {
  char *p;
  char tag1, tag2;

  debug(printf("Entering Samread_splice_strand with %s\n",auxinfo));

  p = auxinfo;
  while (*p != '\0' && *p != '\n') {
    tag1 = p[0];
    tag2 = p[1];

    if (tag1 == 'X' && tag2 == 'S') {
      debug(printf("Found tag XS\n"));
      /* XS:A: */
      p += 5;

      if (*p == '+') {
	return '+';
      } else if (*p == '-') {
	return '-';
      } else if (*p == '?') {
	return '?';
      } else {
	fprintf(stderr,"Cannot parse strand %c after XS tag\n",*p);
	return ' ';
      }
    } else {
      while (*p != '\0' && *p != '\t') {
	p++;
      }
      if (*p == '\t') {
	p++;
      }
    }
  }

  return ' ';
}
#endif


#if 0
Intlist_T
Samread_parse_cigar (Uintlist_T *npositions, int *readlength, char *cigar) {
  Intlist_T types = NULL;
  unsigned int npos;
  char *p, type;

  *npositions = (Uintlist_T) NULL;
  *readlength = 0;

  if (cigar[0] == '*') {
    return (Intlist_T) NULL;
  }

  p = cigar;
  while (*p != '\0') {
    if (sscanf(p,"%u",&npos) != 1) {
      fprintf(stderr,"Unable to parse cigar %s.  No number in %s\n",cigar,p);
      abort();
    } else {
      *npositions = Uintlist_push(*npositions,npos);
    }

    while (*p != '\0' && isdigit(*p)) {
      p++;
    }
    if (*p == '\0') {
      fprintf(stderr,"Unable to parse cigar %s.  No letter after number %u\n",cigar,npos);
      exit(9);
    } else {
      type = *p++;
      types = Intlist_push(types,(int) type);
    }

    if (type == 'S' || type == 'M' || type == 'I') {
      *readlength += npos;
    } else if (type == 'H') {
      *readlength += npos;
    } else if (type == 'D' || type == 'N') {
      /* Ignore */
    } else {
      fprintf(stderr,"Unable to parse cigar %s.  Do not recognize letter %c\n",cigar,type);
      exit(9);
    }
  }

  *npositions = Uintlist_reverse(*npositions);
  return Intlist_reverse(types);
}
#endif


#if 0
void
Samread_print_cigar (Intlist_T types, Uintlist_T npositions) {
  Intlist_T p;
  Uintlist_T q;

  for (p = types, q = npositions; p != NULL; p = Intlist_next(p), q = Uintlist_next(q)) {
    printf("%u%c",Uintlist_head(q),Intlist_head(p));
  }
  return;
}
#endif


#if 0
Chrpos_T
Samread_chrpos_high (Intlist_T types, Uintlist_T npositions, Chrpos_T chrpos_low) {
  Intlist_T p;
  Uintlist_T q;
  Chrpos_T chrpos_high;
  int type;

  chrpos_high = chrpos_low;
  for (p = types, q = npositions; p != NULL; p = Intlist_next(p), q = Uintlist_next(q)) {
    if ((type = Intlist_head(p)) == 'S') {
      /* Ignore */

    } else if (type == 'H') {
      /* Ignore */

    } else if (type == 'M') {
      chrpos_high += Uintlist_head(q);

    } else if (type == 'N') {
      chrpos_high += Uintlist_head(q);

    } else if (type == 'I') {
      /* Do nothing */

    } else if (type == 'D') {
      /* CHECK */
      chrpos_high += Uintlist_head(q);

    } else {
      fprintf(stderr,"Cannot parse type %c\n",type);
      exit(9);
    }
    debug(printf("  type = %c, chrpos = %u\n",type,chrpos_high));
  }

  return chrpos_high - 1U;
}
#endif


#if 0
int
Samread_get_query_coordinates (int *query5, int *query3, Intlist_T types, Uintlist_T npositions,
			       int readlength, char *cigar) {
  int validlength;
  Intlist_T p;
  Uintlist_T q;
  int type;

  *query5 = 1;			/* 1-based */
  *query3 = readlength;
  validlength = 0;

  p = types;
  q = npositions;
  while (p != NULL) {
    if ((type = Intlist_head(p)) == 'S') {
      if (p == types) {
	*query5 = Uintlist_head(q) + 1; /* 1-based */
      } else if (Intlist_next(p) == NULL) {
	*query3 = readlength - Uintlist_head(q);
      } else {
	fprintf(stderr,"Cannot parse cigar %s.  Type S occurs in middle\n",cigar);
	exit(9);
      }
    } else if (type == 'H') {
      /* Do nothing */
    } else if (type == 'M') {
      validlength += Uintlist_head(q);
    } else if (type == 'N') {
      /* Do nothing */
    } else if (type == 'I') {
      validlength += Uintlist_head(q);
    } else if (type == 'D') {
      /* Do nothing */
    }
    p = Intlist_next(p);
    q = Uintlist_next(q);
  }

  debug(printf("Got query %d to %d, with length %d\n",*query5,*query3,validlength));
  if (validlength != (*query3) - (*query5) + 1) {
    fprintf(stderr,"Validlength %d from cigar != %d - %d + 1\n",validlength,*query3,*query5);
    abort();
  }

  return validlength;
}
#endif



#if 0
int
get_substrings (int *querylength, int **query_starts, Chrpos_T **genomic_starts, Chrpos_T **genomic_ends,
		char *cigar, Chrpos_T chrpos_low) {
  int nsubstrings = 0;
  unsigned int npos;
  char *p, type;

  int querypos = 0;
  Chrpos_T genomicpos = chrpos_low;
  Intlist_T query_starts_list = NULL;
  Uintlist_T genomic_starts_list = NULL, genomic_ends_list = NULL;

  if (cigar[0] == '*') {
    *querylength = 0;
    *query_starts = (int *) NULL;
    *genomic_starts = (Chrpos_T *) NULL;
    *genomic_ends = (Chrpos_T *) NULL;
    return 0;
  }

  query_starts_list = Intlist_push(NULL,querypos);
  genomic_starts_list = Uintlist_push(NULL,genomicpos);

  p = cigar;
  while (*p != '\0') {
    if (sscanf(p,"%u",&npos) != 1) {
      fprintf(stderr,"Unable to parse cigar %s in get_substrings.  No number in %s\n",cigar,p);
      abort();
    }

    while (*p != '\0' && isdigit(*p)) {
      p++;
    }
    if (*p == '\0') {
      fprintf(stderr,"Unable to parse cigar %s.  No letter after number %u\n",cigar,npos);
      exit(9);
    } else {
      type = *p++;
    }

    if (type == 'S') {
      querypos += npos;

    } else if (type == 'M') {
      querypos += npos;
      genomicpos += npos;

    } else if (type == 'I') {
      querypos += npos;

    } else if (type == 'H') {
      /* ? querypos += npos; */

    } else if (type == 'D') {
      genomicpos += npos;

    } else if (type == 'N') {
      genomic_ends_list = Uintlist_push(genomic_ends_list,genomicpos);
      /* nsubstrings++; */

      genomicpos += npos;

      query_starts_list = Intlist_push(query_starts_list,querypos);
      genomic_starts_list = Uintlist_push(genomic_starts_list,genomicpos);

    } else {
      fprintf(stderr,"Unable to parse cigar %s.  Do not recognize letter %c\n",cigar,type);
      exit(9);
    }
  }

  *querylength = querypos;
  genomic_ends_list = Uintlist_push(genomic_ends_list,genomicpos);
  /* nsubstrings++; */


  /* Convert lists to arrays */
  query_starts_list = Intlist_reverse(query_starts_list);
  *query_starts = Intlist_to_array(&nsubstrings,query_starts_list);
  Intlist_free(&query_starts_list);

  genomic_starts_list = Uintlist_reverse(genomic_starts_list);
  *genomic_starts = Uintlist_to_array(&nsubstrings,genomic_starts_list);
  Uintlist_free(&genomic_starts_list);

  genomic_ends_list = Uintlist_reverse(genomic_ends_list);
  *genomic_ends = Uintlist_to_array(&nsubstrings,genomic_ends_list);
  Uintlist_free(&genomic_ends_list);

  return nsubstrings;
}
#endif



#if 0
int
Samread_compute_insert_length (int *querylength5, int *querylength3,
			       char *cigar5, Chrpos_T chrpos_low_5, char *cigar3, Chrpos_T chrpos_low_3) {
  int insert_length;
  int nsubstrings5, nsubstrings3, i, j;
  int *query_starts_5, *query_starts_3;
  Chrpos_T *genomic_starts_5, *genomic_ends_5, *genomic_starts_3, *genomic_ends_3;
  Chrpos_T pos5, pos3;

  if (cigar5[0] == '*' || cigar3[0] == '*') {
    return 0;
  }

  nsubstrings5 = get_substrings(&(*querylength5),&query_starts_5,&genomic_starts_5,&genomic_ends_5,cigar5,chrpos_low_5);
  nsubstrings3 = get_substrings(&(*querylength3),&query_starts_3,&genomic_starts_3,&genomic_ends_3,cigar3,chrpos_low_3);

  for (i = 0; i < nsubstrings5; i++) {
    for (j = 0; j < nsubstrings3; j++) {
      if (genomic_ends_5[i] < genomic_starts_3[j]) {
	/* No overlap */
      } else if (genomic_starts_5[i] > genomic_ends_3[j]) {
	/* No overlap */
      } else {
	pos5 = genomic_starts_5[i] - query_starts_5[i];
	pos3 = genomic_starts_3[j] - query_starts_3[j];

	FREE(query_starts_5);
	FREE(genomic_starts_5);
	FREE(genomic_ends_5);
	FREE(query_starts_3);
	FREE(genomic_starts_3);
	FREE(genomic_ends_3);

	if (pos5 > pos3) {
	  return (int) (pos5 - pos3);
	} else {
	  return (int) (pos3 - pos5);
	}
      }
    }
  }

  if (genomic_ends_5[nsubstrings5-1] < genomic_starts_3[0]) {
    insert_length = genomic_starts_3[0] - genomic_ends_5[nsubstrings5-1] + (*querylength5) + (*querylength3);
  } else if (genomic_ends_3[nsubstrings3-1] < genomic_starts_5[0]) {
    insert_length = genomic_starts_5[0] - genomic_ends_3[nsubstrings3-1] + (*querylength5) + (*querylength3);
  } else {
    insert_length = 0;
  }

  FREE(query_starts_5);
  FREE(genomic_starts_5);
  FREE(genomic_ends_5);

  FREE(query_starts_3);
  FREE(genomic_starts_3);
  FREE(genomic_ends_3);

  return insert_length;
}
#endif


