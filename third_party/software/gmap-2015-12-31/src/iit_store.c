static char rcsid[] = "$Id: iit_store.c 153955 2014-11-24 17:54:45Z twu $";
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
#include <string.h>		/* For strlen */
#include <strings.h>		/* For rindex */
#include <ctype.h>
#include <math.h>		/* For qsort and NAN */
#include "bool.h"
#include "types.h"
#include "assert.h"
#include "mem.h"
#include "fopen.h"

#include "list.h"
#include "doublelist.h"
#include "univinterval.h"
#include "interval.h"
#include "table.h"
#include "tableint.h"
#include "chrom.h"
#include "iit-write-univ.h"
#include "iit-write.h"
#include "getopt.h"

#ifndef NAN
#define NAN nan("")
#endif

#ifndef NAN
static double NAN = nan("")
#endif


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

#define LINELENGTH 8192
#define MONITOR_INTERVAL 100000 /* 100 thousand entries */

/************************************************************************
 *   Program options
 ************************************************************************/

static char *outputfile = NULL;
static bool old_format_p = false;
static bool gff3_format_p = false;
static char *labelid = "ID";
static bool fieldsp = false;
static char iit_version = 0;

static Sorttype_T divsort = CHROM_SORT;
static char *mitochondrial_string = NULL;


static struct option long_options[] = {
  /* Input options */
  {"output", required_argument, 0, 'o'}, /* outputfile */
  {"oldformat", no_argument, 0, '1'}, /* old_format_p */
  {"fields", no_argument, 0, 'F'}, /* fieldsp */
  {"gff", no_argument, 0, 'G'}, /* gff3_format_p */
  {"label", required_argument, 0, 'l'}, /* labelid */
  {"iitversion", required_argument, 0, 'v'}, /* iit_version */
  {"sort", required_argument, 0, 's'}, /* sorttype */

  /* Help options */
  {"version", no_argument, 0, 'V'}, /* print_program_version */
  {"help", no_argument, 0, '?'}, /* print_program_usage */
  {0, 0, 0, 0}
};

static void
print_program_version () {
  fprintf(stdout,"\n");
  fprintf(stdout,"iit_store: indexing utility for Interval Index Trees\n");
  fprintf(stdout,"Part of GMAP package, version %s\n",PACKAGE_VERSION);
  fprintf(stdout,"Thomas D. Wu, Genentech, Inc.\n");
  fprintf(stdout,"Contact: twu@gene.com\n");
  fprintf(stdout,"\n");
  return;
}

static void
print_program_usage () {
  fprintf(stdout,"\
Usage: iit_store [OPTIONS...] -o outputfile inputfile, or\n\
       cat inputfile | iit_store [OPTIONS...] -o outputfile\n\
where\n\
   outputfile is the desired filename for the iit file\n\
       (.iit will be added as a suffix if necessary), and\n\
   inputfile is in either FASTA or GFF3 format, as described below.\n\
\n\
Options\n\
  -o, --output=STRING       Name of output iit file\n\
  -1, --oldformat           Old format for intervals:\n\
                             <start> <optional end> <optional div> <optional type>\n\
  -F, --fields              Annotation consists of separate fields\n\
  -G, --gff                 Parse input file in gff3 format\n\
  -l, --label=STRING        For gff input, the feature attribute to use (default is ID)\n\
\n\
  -s, --sort=STRING         Sorting of divisions: none, alpha, numeric-alpha, or chrom (default)\n\
                numeric-alpha: chr1 chr1_random chr2 chr10 chr10_random chrM chrUn chrX chrY\n\
                        chrom: chr1 chr2 chr10 chrX chrY chrM chr1_random chr10_random chrUn\n\
\n\
                               Note 1: For sorting purposes, any initial 'chr' will be ignored\n\
                               Note 2: For chrom, X, Y, M, MT (or chrX, chrY, and so on) are special\n\
\n\
  -v, --iitversion=STRING   Desired iit version for output iit\n\
                            (default = 0, which means latest version)\n\
\n\
  -V, --version             Show version\n\
  -?, --help                Show this help message\n\
\n\
\n\
Description of input format:\n\
\n\
The FASTA format for input files should be\n\
\n\
    >label interval [type] [/value=value]\n\
    optional_annotation (which may be zero, one, or multiple lines)\n\
\n\
where intervals have one of the following forms:\n\
   div:start..end\n\
   div:start\n\
   start..end\n\
   start\n\
and a given type, numeric value, or both is optional.  A numeric value\n\
allows intervals to be searched by a range of values using iit_get.\n\
\n\
Intervals may have directions.  To indicate a forward direction,\n\
the start coordinate should be less than the end coordinate.\n\
To indicate a reverse direction, the start coordinate should be\n\
greater than the end coordinate. If they are the same, then no\n\
direction is implied.  If no end coordinate is given, the end\n\
coordinate is assumed to be the same as the start coordinate.\n\
\n\
For example, the label may be a sequence accession, with the div representing\n\
a chromosome, and the type representing an additional piece of information\n\
A header might therefore look like\n\
\n\
    >NM_004448 17:35138441..35109780 refseq\n\
\n\
which indicates an interval on chromosome 17 in the reverse direction,\n\
and of type refseq.\n\
\n\
If the -F flag is provided, IIT files may store annotation for each interval\n\
as separate fields.  The input must contain the names of the fields, one per\n\
line, before the first interval header.  Each interval then contains annotation\n\
corresponding to each field, one value per line.\n\
\n\
The GFF3 format requires the -G flag and optionally the -l flag.\n\
The iit_store program will parse the chromosome from column 1, the start\n\
coordinate from column 4, the end coordinate from column 5, the strand\n\
from column 7, an if possible, the label from column 9.  The -l flag\n\
will indicate which feature from column 9 to retrieve, such as ID, Name,\n\
or Parent.  Appropriate choice of label may be helpful later on, because\n\
the iit_get program can retrieve information by label, as well as by\n\
coordinates.\n\
\n\
Limitations: Start and end coordinates must be non-negative integers, and are\n\
limited to the domain of a 64-bit quantity, which means coordinates must be\n\
less than 2^64.  If your machine is a 32-bit machine, coordinates must be less\n\
than 2^32 = 4294967295.\n\
\n\
See also: iit_get, iit_dump\n\
");
  return;
}

/* Empties contents of lines */
static char *
concatenate_lines (List_T lines, int content_size) {
  char *string, *temp;
  List_T l;

  string = (char *) CALLOC(content_size+1,sizeof(char));
  for (l = lines; l; l = List_next(l)) {
    temp = (char *) List_head(l);
    strcat(string,temp);
    FREE(temp);
  }
  
  /* Keep last return
  if (string[content_size-1] == '\n') {
    string[content_size-1] = '\0';
  }
  */

  return string;
}



/* Note that isnumber is a function in ctype.h on some systems */
static bool
isnumberp (Univcoord_T *result, char *string) {
  char *p = string;

  *result = 0U;
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

static bool
isrange (Univcoord_T *start, Univcoord_T *end, char *string) {
  bool result;
  Univcoord_T length;
  char *copy, *startstring, *endstring;

  copy = (char *) CALLOC(strlen(string)+1,sizeof(char));
  strcpy(copy,string);

  if (index(copy,'.')) {
    startstring = strtok(copy,"..");
    endstring = strtok(NULL,"..");
    return (isnumberp(&(*start),startstring) && isnumberp(&(*end),endstring));

  } else if (index(copy,'+')) {
    startstring = strtok(copy,"+");
    endstring = strtok(NULL,"+");
    if (!isnumberp(&(*start),startstring)) {
      result = false;
    } else if (endstring[0] == '-' && isnumberp(&length,&(endstring[1]))) {
      *end = (*start) - length;
      result = true;
    } else if (!isnumberp(&length,endstring)) {
      result = false;
    } else {
      *end = (*start) + length;
      result = true;
    }

  } else if (index(copy,'-')) {
    /* Old notation */
    startstring = strtok(copy,"--");
    endstring = strtok(NULL,"--");
    return (isnumberp(&(*start),startstring) && isnumberp(&(*end),endstring));

  } else {
    result = false;
  }

  FREE(copy);
  return result;
}


/* Example: >A X:1..10 red.  Here, A is a label, 1 and 10 are start and end, X is a div, and red is a type. */
/* Other variants: >A 1..10 red, or >A 1..10 */
static char *
scan_header_div (int *labellength, bool *seenp, List_T *divlist, List_T *typelist, Tableint_T div_seenp, Tableint_T typetable, 
		 bool *valuep, double *value, char **label, Univcoord_T *start, Univcoord_T *end, int *type,
		 char **restofheader, char *header) {
  char *divstring = NULL, *coords, *copy, Buffer[1024], query[1024], tag[1024], *typestring, *p;
  char *valueptr;

  *seenp = false;
  if (sscanf(header,">%s %s\n",Buffer,query) < 2) {
    fprintf(stderr,"Error parsing %s.  Expecting a FASTA type header with a label, coords (as <div>:<number>..<number>), and optional tag.\n",header);
    exit(9);
  }

  *labellength = strlen(Buffer);
  *label = (char *) CALLOC(*labellength+1,sizeof(char));
  strcpy(*label,Buffer);

  if (!index(query,':')) {
    debug(printf("Query %s has no div\n",query));
    divstring = (char *) CALLOC(1,sizeof(char));
    divstring[0] = '\0';
    coords = query;

  } else {
    debug(printf("Parsed query %s into ",query));
    p = strtok(query,":");

    divstring = (char *) CALLOC(strlen(p)+1,sizeof(char));
    strcpy(divstring,p);

    if (Tableint_get(div_seenp,(void *) divstring) == 0) {
      debug(printf("Entering new div %s.\n",divstring));
      Tableint_put(div_seenp,(void *) divstring,(int) true);
      copy = (char *) CALLOC(strlen(divstring)+1,sizeof(char));
      strcpy(copy,divstring);
      *divlist = List_push(*divlist,copy);
      *seenp = false;
    } else {
      *seenp = true;
    }

    coords = strtok(NULL,":");
    debug(printf("div %s and coords %s\n",divstring,coords));
  }
      
  if (coords == NULL) {
    fprintf(stderr,"Error parsing %s.  Expecting coords (as <div>:<number>..<number>)\n",query);
    fprintf(stderr,"Problematic line was: %s\n",header);
    exit(9);
  } else if (isnumberp(&(*start),coords)) {
    debug(printf("  and coords %s as a number\n",coords));
    *end = *start;
  } else if (isrange(&(*start),&(*end),coords)) {
    debug(printf("  and coords %s as a range starting at %llu and ending at %llu\n",
		 coords,(unsigned long long) *start,(unsigned long long) *end));
  } else {
    fprintf(stderr,"Error parsing %s:%s.  Expecting coords (as <div>:<number>..<number>)\n",query,coords);
    fprintf(stderr,"Problematic line was: %s\n",header);
    exit(9);
  }

  if ((valueptr = strstr(header,"/value=")) == NULL) {
    *value = NAN;
  } else {
    /* Note: Not checking for any errors */
    *valuep = true;
    valueptr += strlen("/value=");
    *value = atof(valueptr);
  }

  if (sscanf(header,">%s %s %s",Buffer,query,tag) < 3) {
    *type = 0;
    *restofheader = (char *) NULL;

  } else if (!strncmp(tag,"/value=",strlen("/value="))) {
    *type = 0;

    /* Get rest of header */
    p = header;
    while (!isspace(*p)) p++;	/* accession */
    while (isspace(*p)) p++;

    while (!isspace(*p)) p++;	/* coords */
    while (isspace(*p)) p++;

    if (*p == '\0') {
      *restofheader = (char *) NULL;
    } else {
      *restofheader = (char *) CALLOC(strlen(p)+1,sizeof(char));
      strcpy(*restofheader,p);
    }

  } else {
    if ((*type = Tableint_get(typetable,(void *) tag)) == 0) {
      /* Store types as 1-based */
      *type = Tableint_length(typetable) + 1;
      typestring = (char *) CALLOC(strlen(tag)+1,sizeof(char));
      strcpy(typestring,tag);
      Tableint_put(typetable,typestring,*type);
      *typelist = List_push(*typelist,typestring);
      /* debug(printf("Entering new type %s.\n",typestring)); */
    }

    /* Get rest of header */
    p = header;
    while (!isspace(*p)) p++;	/* accession */
    while (isspace(*p)) p++;

    while (!isspace(*p)) p++;	/* coords */
    while (isspace(*p)) p++;

    while (*p != '\0' && !isspace(*p)) p++;	/* tag */
    while (*p != '\0' && isspace(*p)) p++;

    if (*p == '\0') {
      *restofheader = (char *) NULL;
    } else {
      *restofheader = (char *) CALLOC(strlen(p)+1,sizeof(char));
      strcpy(*restofheader,p);
    }
  }

  return divstring;
}



/* Example: >A 1 10 X red.  Here, A is a label, 1 and 10 are start and end, X is a div, and red is a type. */
static char *
scan_header_spaces (int *labellength, bool *seenp, List_T *divlist, List_T *typelist, Tableint_T div_seenp, Tableint_T typetable, 
		    bool *valuep, double *value, char **label, Univcoord_T *start, Univcoord_T *end, int *type, char *header) {
  char *divstring, *copy, Buffer[1024], *typestring, *p, *ptr, *divstart;
  char *valueptr;
  int divlength;
  int nscanned;

  /* Example: >A 1 10 X red.  Here, A is a label, 1 and 10 are start and end, X is a div, and red is a type. */

  *seenp = false;
  nscanned = sscanf(header,">%s %llu %llu\n",Buffer,&(*start),&(*end));
  if (nscanned < 3) {
    fprintf(stderr,"Error parsing %s.  Expecting a FASTA type header with a label, two coordinates, and optional tag.\n",header);
    exit(9);
  } else {
    if ((valueptr = strstr(header,"/value=")) == NULL) {
      *value = NAN;
    } else {
      /* Note: Not checking for any errors */
      *valuep = true;
      valueptr += strlen("/value=");
      *value = atof(valueptr);
    }

    *labellength = strlen(Buffer);
    *label = (char *) CALLOC(*labellength+1,sizeof(char));
    strcpy(*label,Buffer);

    p = header;
    while (!isspace((int) *p)) { p++; } /* First word (label) */
    while (isspace((int) *p)) { p++; } /* First space */
    while (!isspace((int) *p)) { p++; } /* Second word (start coord) */
    while (isspace((int) *p)) { p++; } /* Second space */
    while (!isspace((int) *p)) { p++; } /* Third word (end coord) */
    while (*p != '\0' && isspace((int) *p)) { p++; } /* Third space */
    
    if (*p == '\0') {
      divstring = (char *) CALLOC(1,sizeof(char));
      divstring[0] = '\0';
      *type = 0;		/* Empty type string */
    } else {
      divstart = p;
      divlength = 0;
      while (!isspace((int) *p)) { p++; divlength++; } /* Fourth word (div) */

      divstring = (char *) CALLOC(divlength+1,sizeof(char));
      strncpy(divstring,divstart,divlength);
      if (Tableint_get(div_seenp,(void *) divstring) == 0) {
	debug(printf("Entering new div %s.\n",divstring));
	Tableint_put(div_seenp,(void *) divstring,(int) true);
	copy = (char *) CALLOC(divlength+1,sizeof(char));
	strcpy(copy,divstring);
	*divlist = List_push(*divlist,copy);
	*seenp = false;
      } else {
	*seenp = true;
      }

      while (*p != '\0' && isspace((int) *p)) { p++; } /* Fourth space */
      if (*p == '\0') {
	*type = 0;
      } else if (!strncmp(p,"/value=",strlen("/value="))) {
	*type = 0;
      } else {
	if ((ptr = rindex(p,'\n')) != NULL) {
	  while (isspace((int) *ptr)) { ptr--; } /* Erase empty space */
	  ptr++;
	  *ptr = '\0';
	}
	
	if ((*type = Tableint_get(typetable,(void *) p)) == 0) {
	  /* Store types as 1-based */
	  *type = Tableint_length(typetable) + 1;
	  typestring = (char *) CALLOC(strlen(p)+1,sizeof(char));
	  strcpy(typestring,p);
	  Tableint_put(typetable,typestring,*type);
	  *typelist = List_push(*typelist,typestring);
	  /* debug(printf("Entering new type %s.\n",typestring)); */
	}
      }
    }
  }

  return divstring;
}

static List_T
parse_fieldlist (char *firstchar, FILE *fp) {
  List_T fieldlist = NULL;
  char Buffer[LINELENGTH], *fieldname, *p;

  while (!feof(fp) && (*firstchar = fgetc(fp)) != '>') {
    if (*firstchar != EOF) {
      Buffer[0] = *firstchar;
      fgets(&(Buffer[1]),LINELENGTH-1,fp);
      if ((p = rindex(Buffer,'\n')) == NULL) {
	fprintf(stderr,"Line %s exceeds maximum length of %d\n",Buffer,LINELENGTH);
	exit(9);
      } else {
	*p = '\0';
      }
      fieldname = (char *) CALLOC(strlen(Buffer)+1,sizeof(char));
      strcpy(fieldname,Buffer);
      fieldlist = List_push(fieldlist,fieldname);
    }
  }

  return List_reverse(fieldlist);
}


static void
parse_fasta (bool *valuep, Univcoord_T *max_coordinate, Univcoord_T *label_totallength, Univcoord_T *annot_totallength,
	     List_T *divlist, List_T *typelist, Table_T intervaltable, Table_T valuetable, Table_T labeltable, Table_T annottable,
	     FILE *fp, Tableint_T div_seenp, Tableint_T typetable, char firstchar) {
  char Buffer[LINELENGTH], *divstring, *label, *restofheader = NULL, *tempstring;
  double value;
  Univcoord_T start, end;
  List_T lines, d;
  int labellength, content_size, type, nentries;
  bool seenp;

  /* *max_coordinate = 0; */
  *label_totallength = 0;
  *annot_totallength = 0;

  if (feof(fp)) {
    return;

  } else if (firstchar == '\0') {
    fgets(Buffer,LINELENGTH,fp);
  } else {
    Buffer[0] = firstchar;
    fgets(&(Buffer[1]),LINELENGTH-1,fp);
  }
  if (old_format_p == true) {
    divstring = scan_header_spaces(&labellength,&seenp,&(*divlist),&(*typelist),
				   div_seenp,typetable,&(*valuep),&value,&label,&start,&end,&type,Buffer);
    restofheader = (char *) NULL;
  } else {
    divstring = scan_header_div(&labellength,&seenp,&(*divlist),&(*typelist),
				div_seenp,typetable,&(*valuep),&value,&label,&start,&end,&type,&restofheader,Buffer);
  }
  *max_coordinate = start;
  if (end > *max_coordinate) {
    *max_coordinate = end;
  }
  
  Table_put(valuetable,(void *) divstring,
	    Doublelist_push(Table_get(valuetable,(void *) divstring),value));

  *label_totallength = labellength;
  Table_put(labeltable,(void *) divstring,
	    List_push(Table_get(labeltable,(void *) divstring),label));

  lines = NULL;
  content_size = 0;
  if (restofheader != NULL) {
    lines = List_push(lines,(void *) restofheader);
    content_size += strlen(restofheader);
  }

  nentries = 1;			/* Because we already processed the first entry above */
  while (fgets(Buffer,LINELENGTH,fp) != NULL) {
    if (Buffer[0] == '>') {
      if (++nentries % MONITOR_INTERVAL == 0) {
	fprintf(stderr,"Read %d entries in FASTA file...\n",nentries);
      }

      /* Store as Univinterval_T now, but may need to change to Interval_T later */
      Table_put(intervaltable,(void *) divstring,
		List_push(Table_get(intervaltable,(void *) divstring),
			  (void *) Univinterval_new(start,end,type)));

      lines = List_reverse(lines);
      if (restofheader == NULL && content_size > 0) {
	tempstring = (char *) CALLOC(2,sizeof(char));
	tempstring[0] = '\n';
	tempstring[1] = '\0';
	lines = List_push(lines,tempstring);
	content_size += 1;
      }
      *annot_totallength += content_size;

      Table_put(annottable,(void *) divstring,
		List_push(Table_get(annottable,(void *) divstring),
			  (void *) concatenate_lines(lines,content_size)));
      List_free(&lines);

      if (seenp == true) {
	FREE(divstring);
      }
      if (old_format_p == true) {
	divstring = scan_header_spaces(&labellength,&seenp,&(*divlist),&(*typelist),div_seenp,typetable,
				       &(*valuep),&value,&label,&start,&end,&type,Buffer);
	restofheader = (char *) NULL;
      } else {
	divstring = scan_header_div(&labellength,&seenp,&(*divlist),&(*typelist),div_seenp,typetable,
				    &(*valuep),&value,&label,&start,&end,&type,&restofheader,Buffer);
      }
      if (start > *max_coordinate) {
	*max_coordinate = start;
      }
      if (end > *max_coordinate) {
	*max_coordinate = end;
      }

      Table_put(valuetable,(void *) divstring,
		Doublelist_push(Table_get(valuetable,(void *) divstring),value));

      *label_totallength += labellength;
      Table_put(labeltable,(void *) divstring,
		List_push(Table_get(labeltable,(void *) divstring),label));

      lines = NULL;
      content_size = 0;
      if (restofheader != NULL) {
	lines = List_push(lines,(void *) restofheader);
	content_size += strlen(restofheader);
      }

    } else {
      tempstring = (char *) CALLOC(strlen(Buffer)+1,sizeof(char));
      strcpy(tempstring,Buffer);
      lines = List_push(lines,(void *) tempstring);
      content_size += strlen(Buffer);
    }
  }
  fprintf(stderr,"Finished reading FASTA file -- total entries: %d\n",nentries);

  /* Store as Univinterval_T now, but may need to change later */
  Table_put(intervaltable,(void *) divstring,
	    List_push(Table_get(intervaltable,(void *) divstring),
		      (void *) Univinterval_new(start,end,type)));

  lines = List_reverse(lines);
  if (restofheader == NULL && content_size > 0) {
    tempstring = (char *) CALLOC(2,sizeof(char));
    tempstring[0] = '\n';
    tempstring[1] = '\0';
    lines = List_push(lines,tempstring);
    content_size += 1;
  }
  *annot_totallength += content_size;
  Table_put(annottable,(void *) divstring,
	    List_push(Table_get(annottable,(void *) divstring),
		      (void *) concatenate_lines(lines,content_size)));
  List_free(&lines);

  if (seenp == true) {
    FREE(divstring);
  }
  
  fprintf(stderr,"Maximum coordinate: %llu\n",(unsigned long long) *max_coordinate);
  fprintf(stderr,"Total label length: %llu + %d separators\n",(unsigned long long) *label_totallength,nentries);
  fprintf(stderr,"Total annotation length: %llu + %d separators\n",(unsigned long long) *annot_totallength,nentries);
  *label_totallength += nentries;
  *annot_totallength += nentries;

  /* Reverse all lists */
  fprintf(stderr,"Saw %d distinct divisions/chromosomes\n",List_length(*divlist)-1);
  *divlist = List_reverse(*divlist);

  fprintf(stderr,"Saw %d distinct tags/types\n",List_length(*typelist));
  *typelist = List_reverse(*typelist);

  for (d = *divlist; d != NULL; d = List_next(d)) {
    divstring = (char *) List_head(d);
    Table_put(intervaltable,(void *) divstring,
	      List_reverse((List_T) Table_get(intervaltable,(void *) divstring)));
    Table_put(valuetable,(void *) divstring,
	      Doublelist_reverse((Doublelist_T) Table_get(valuetable,(void *) divstring)));
    Table_put(labeltable,(void *) divstring,
	      List_reverse((List_T) Table_get(labeltable,(void *) divstring)));
    Table_put(annottable,(void *) divstring,
	      List_reverse((List_T) Table_get(annottable,(void *) divstring)));
  }

  return;
}


static int
assign_columns (char **columns, char *Buffer, int maxfields) {
  char *token;
  int nfields = 0;
  
  columns[nfields++] = token = strtok(Buffer,"\t");
  while ((token = strtok(NULL,"\t")) != NULL && nfields < maxfields) {
    columns[nfields++] = token;
  }
  return nfields;
}


#define CHRCOLUMN 0
#define STARTCOLUMN 3
#define ENDCOLUMN 4
#define STRANDCOLUMN 6
#define FEATURECOLUMN 8
#define GFF3_COLUMNS 9

/* Modifies feature */
static char *
gff3_feature_id (char *feature, char *labelstr, int labelstrlen, int lineno) {
  char *token, *value, *p;

  token = strtok(feature,";");
  if (!strncmp(token,labelstr,labelstrlen)) {
    value = &(token[labelstrlen]);
    if (value[0] != '"') {
      return value;
    } else {
      value = &(value[1]);
      /* Quotation marks */
      if ((p = rindex(value,'"')) == NULL) {
	fprintf(stderr,"Error in line %d: Saw no matching quotation in %s\n",lineno,token);
	exit(9);
      } else {
	*p = '\0';
      }
      return value;
    }
  } else {
    while ((token = strtok(NULL,";")) != NULL) {
      if (!strncmp(token,labelstr,labelstrlen)) {
	value = &(token[labelstrlen]);
	if (value[0] != '"') {
	  return value;
	} else {
	  value = &(value[1]);
	  /* Quotation marks */
	  if ((p = rindex(value,'"')) == NULL) {
	    fprintf(stderr,"Error in line %d: Saw no matching quotation in %s\n",lineno,token);
	    exit(9);
	  } else {
	    *p = '\0';
	  }
	  return value;
	}
      }
    }
    return NULL;
  }
}

static bool
empty_line_p (char *line) {
  char *p = line;

  while (*p != '\0' && isspace(*p)) {
    p++;
  }
  if (*p == '\0') {
    return true;
  } else {
    return false;
  }
}

static void
parse_gff3 (List_T *divlist, Table_T intervaltable, Table_T labeltable, Table_T annottable,
	    FILE *fp, Tableint_T div_seenp) {
  char Buffer[LINELENGTH], Space[1000], *columns[GFF3_COLUMNS],
    *divstring, *label, *p, *chr, *line, *idptr;
  List_T d;
  Univcoord_T start, end;
  int nfields, lineno = 0, row = 0, labelstrlen;
  char strandchar;
  char *labelstr;

  labelstr = (char *) CALLOC(strlen(labelid) + strlen("=") + 1,sizeof(char));
  sprintf(labelstr,"%s=",labelid);
  labelstrlen = strlen(labelstr);

  while (fgets(Buffer,LINELENGTH,fp) != NULL) {
    lineno++;
    if (Buffer[0] == '#') {
      /* Skip comment */
    } else if (empty_line_p(Buffer) == true) {
      /* Skip empty line */
    } else {
      line = (char *) CALLOC(strlen(Buffer)+1,sizeof(char));
      strcpy(line,Buffer);

      if ((p = rindex(Buffer,'\n')) == NULL) {
	fprintf(stderr,"Line exceeds maximum length of %d\n",LINELENGTH);
	exit(9);
      } else {
	*p = '\0';
      }
      nfields = assign_columns(columns,Buffer,GFF3_COLUMNS); /* destroys Buffer */

      if (nfields < GFF3_COLUMNS-1) {
	/* Subract 1 to allow for an empty feature column */
	fprintf(stderr,"Skipping line %d with only %d fields: %s\n",lineno,nfields,line);
	FREE(line);
      } else {
	chr = columns[CHRCOLUMN];
	divstring = (char *) CALLOC(strlen(chr)+1,sizeof(char));
	sprintf(divstring,"%s",chr);
	
	if ((strandchar = columns[STRANDCOLUMN][0]) == '+') {
	  start = atof(columns[STARTCOLUMN]);
	  end = atof(columns[ENDCOLUMN]);
	} else if (strandchar == '-') {
	  start = atof(columns[ENDCOLUMN]);
	  end = atof(columns[STARTCOLUMN]);
	} else if (strandchar == '.' || strandchar == '?') {
	  start = atof(columns[STARTCOLUMN]);
	  end = atof(columns[ENDCOLUMN]);
	} else {
	  start = atof(columns[STARTCOLUMN]);
	  end = atof(columns[ENDCOLUMN]);
	}

	if (Tableint_get(div_seenp,(void *) divstring) == 0) {
	  Tableint_put(div_seenp,(void *) divstring,(int) true);
	  *divlist = List_push(*divlist,divstring);
	}

	/* Store Univinterval_T now, but may need to change later */
	Table_put(intervaltable,(void *) divstring,
		  List_push(Table_get(intervaltable,(void *) divstring),
			    (void *) Univinterval_new(start,end,/*type*/0)));

	if (nfields <= FEATURECOLUMN) {
	  sprintf(Space,"gff.%d",row);
	  label = (char *) CALLOC(strlen(Space)+1,sizeof(char));
	  strcpy(label,Space);
	} else if ((idptr = gff3_feature_id(columns[FEATURECOLUMN],labelstr,labelstrlen,lineno)) == NULL) {
	  sprintf(Space,"gff.%d",row);
	  label = (char *) CALLOC(strlen(Space)+1,sizeof(char));
	  strcpy(label,Space);
	} else {
	  label = (char *) CALLOC(strlen(idptr)+1,sizeof(char));
	  strcpy(label,idptr);
	}
	Table_put(labeltable,(void *) divstring,
		  List_push(Table_get(labeltable,(void *) divstring),label));
	Table_put(annottable,(void *) divstring,
		  List_push(Table_get(annottable,(void *) divstring),line));

	row++;
      }
    }
  }

  *divlist = List_reverse(*divlist);

  for (d = *divlist; d != NULL; d = List_next(d)) {
    divstring = (char *) List_head(d);
    Table_put(intervaltable,(void *) divstring,
	      List_reverse((List_T) Table_get(intervaltable,(void *) divstring)));
    Table_put(labeltable,(void *) divstring,
	      List_reverse((List_T) Table_get(labeltable,(void *) divstring)));
    Table_put(annottable,(void *) divstring,
	      List_reverse((List_T) Table_get(annottable,(void *) divstring)));
  }

  FREE(labelstr);

  return;
}


#ifdef __STRICT_ANSI__
int getopt (int argc, char *const argv[], const char *optstring);
#endif

int 
main (int argc, char *argv[]) {
  char *inputfile = NULL, *iitfile, *tempstring, *divstring, *typestring, *p;
  char firstchar;
  List_T d, l, templist = NULL, divlist = NULL, typelist = NULL, fieldlist = NULL;
  Doublelist_T valuelist;
  List_T newlist;
  FILE *fp;
  Univinterval_T univinterval;
  Interval_T interval;
  Tableint_T div_seenp, typetable;
  Table_T intervaltable, labeltable, valuetable, annottable;
  Chrom_T *chroms = NULL;
  int n_proper_divs = 0, i;
  bool coord_values_8p, label_pointers_8p, annot_pointers_8p, valuep = false;
  Univcoord_T order;
  Univcoord_T max_coordinate, label_totallength, annot_totallength;

  int opt;
  extern int optind;
  extern char *optarg;
  int long_option_index = 0;

  while ((opt = getopt_long(argc,argv,"o:1FGl:v:s:",
			    long_options,&long_option_index)) != -1) {
    switch (opt) {
    case 'o': outputfile = optarg; break;
    case '1': old_format_p = true; break;
    case 'F': fieldsp = true; break;
    case 'G': gff3_format_p = true; break;
    case 'l': labelid = optarg; break;
    case 'v': iit_version = atoi(optarg); break;
    case 's': 
      if (!strcmp(optarg,"none")) {
	divsort = NO_SORT;
      } else if (!strcmp(optarg,"alpha")) {
	divsort = ALPHA_SORT;
      } else if (!strcmp(optarg,"numeric-alpha")) {
	divsort = NUMERIC_ALPHA_SORT;
      } else if (!strcmp(optarg,"chrom")) {
	divsort = CHROM_SORT;
      } else {
	fprintf(stderr,"Don't recognize sort type %s.  Allowed values are none, alpha, or chrom.",optarg);
	exit(9);
      }
      break;
    case 'V': print_program_version(); exit(0);
    case '?': print_program_usage(); exit(0);
    default: exit(9);
    }
  }
  argc -= optind;
  argv += optind;

  if (outputfile == NULL) {
    fprintf(stderr,"Need to specify an output file with the -o flag\n");
    exit(9);
  } else if (iit_version > IIT_LATEST_VERSION_NOVALUES && iit_version > IIT_LATEST_VERSION_VALUES) {
    fprintf(stderr,"version %d requested, but this program can write only up to version %d or %d\n",
	    iit_version,IIT_LATEST_VERSION_NOVALUES,IIT_LATEST_VERSION_VALUES);
    exit(9);
  }

  if (argc < 1) {
    fp = stdin;
  } else {
    inputfile = argv[0];
    fp = FOPEN_READ_TEXT(inputfile);
    if (!fp) {
      fprintf(stderr,"Can't open file %s\n",inputfile);
      exit(9);
    }
  }

  div_seenp = Tableint_new(100,Table_string_compare,Table_string_hash);
  typetable = Tableint_new(100,Table_string_compare,Table_string_hash);
  intervaltable = Table_new(100,Table_string_compare,Table_string_hash);
  valuetable = Table_new(100,Table_string_compare,Table_string_hash);
  labeltable = Table_new(100,Table_string_compare,Table_string_hash);
  annottable = Table_new(100,Table_string_compare,Table_string_hash);
  
  /* The zeroth div is empty */
  divstring = (char *) CALLOC(1,sizeof(char));
  divstring[0] = '\0';
  divlist = List_push(NULL,divstring);

  /* The zeroth type is empty */
  typestring = (char *) CALLOC(1,sizeof(char));
  typestring[0] = '\0';
  typelist = List_push(NULL,typestring);

  if (gff3_format_p == true) {
    parse_gff3(&divlist,intervaltable,labeltable,annottable,fp,div_seenp);
  } else {
    fieldlist = parse_fieldlist(&firstchar,fp);
    parse_fasta(&valuep,&max_coordinate,&label_totallength,&annot_totallength,
		&divlist,&typelist,intervaltable,valuetable,labeltable,annottable,
		fp,div_seenp,typetable,firstchar);
  }

  if (inputfile != NULL) {
    fclose(fp);
  }


  if (iit_version == 0 && List_length(divlist) == 1) {
    /* No divs other than NULL */
    fprintf(stderr,"No divs/chromosomes provided, so storing as IIT version 1\n");
    iit_version = 1;
  }

  coord_values_8p = false;
  label_pointers_8p = false;
  annot_pointers_8p = false;
#ifdef HAVE_64_BIT
  if (max_coordinate > 4294967295U) {
    coord_values_8p = true;
  }
  if (iit_version == 0) {
    if (label_totallength > 4294967295U) {
      label_pointers_8p = true;
    }
    if (annot_totallength > 4294967295U) {
      annot_pointers_8p = true;
    }
    if (valuep == true) {
      iit_version = IIT_LATEST_VERSION_VALUES;
    } else {
      iit_version = IIT_LATEST_VERSION_NOVALUES;
    }

  } else if (iit_version == 4) {
    if (label_totallength > 4294967295U || annot_totallength > 4294967295U) {
      /* Both pointer types have to match */
      label_pointers_8p = true;
      annot_pointers_8p = true;
    }

  } else if (iit_version <= 3) {
    if (label_totallength > 4294967295U || annot_totallength > 4294967295U) {
      fprintf(stderr,"Need 8-byte pointers, which requires you to specify a version of 4 or greater\n");
      exit(9);
    }
  }
#else
  if (iit_version == 0) {
    if (valuep == true) {
      iit_version = IIT_LATEST_VERSION_VALUES;
    } else {
      iit_version = IIT_LATEST_VERSION_NOVALUES;
    }
  }
#endif

  if (iit_version == 1) {
    /* Will use Univinterval_T objects, which may print as UINT8 or UINT4 */
  } else if (coord_values_8p == true) {
    fprintf(stderr,"Cannot have large coordinates, except for chromosome IIT files\n");
    exit(9);
  } else {
    /* Convert all Univinterval_T objects to Interval_T objects */
    
    for (d = divlist; d != NULL; d = List_next(d)) {
      divstring = (char *) List_head(d);
      templist = (List_T) Table_get(intervaltable,(void *) divstring);
      newlist = (List_T) NULL;
      for (l = templist; l != NULL; l = List_next(l)) {
	univinterval = (Univinterval_T) List_head(l);
	if (Univinterval_sign(univinterval) < 0) {
	  newlist = List_push(newlist,
			      (void *) Interval_new(Univinterval_high(univinterval),
						    Univinterval_low(univinterval),
						    Univinterval_type(univinterval)));
	} else {
	  newlist = List_push(newlist,
			      (void *) Interval_new(Univinterval_low(univinterval),
						    Univinterval_high(univinterval),
						    Univinterval_type(univinterval)));
	}
	Univinterval_free(&univinterval);
      }
      Table_put(intervaltable,(void *) divstring,(void *) List_reverse(newlist));
      List_free(&templist);
    }
  }


  /* Figure out name of iit file */
  if (strlen(outputfile) < 4) {
    iitfile = (char *) CALLOC(strlen(outputfile)+strlen(".iit")+1,sizeof(char));
    sprintf(iitfile,"%s.iit",outputfile);
  } else {
    p = &(outputfile[strlen(outputfile)]);
    p -= 4;
    if (!strcmp(p,".iit")) {
      iitfile = (char *) CALLOC(strlen(outputfile)+1,sizeof(char));
      strcpy(iitfile,outputfile);
    } else {
      iitfile = (char *) CALLOC(strlen(outputfile)+strlen(".iit")+1,sizeof(char));
      sprintf(iitfile,"%s.iit",outputfile);
    }
  }

  order = 0;
  if ((n_proper_divs = List_length(divlist) - 1) > 0) {
    chroms = (Chrom_T *) CALLOC(n_proper_divs,sizeof(Chrom_T));
    for (l = divlist, i = 0; l != NULL; l = List_next(l)) {
      tempstring = List_head(l);
      if (tempstring[0] == '\0') {
	/* FREE(tempstring); -- Causes invalid read later in table_string_compare */
      } else {
	chroms[i++] = Chrom_from_string(tempstring,mitochondrial_string,order++,/*circularp*/false);
      }
    }
  }

  List_free(&divlist);

  switch (divsort) {
  case NO_SORT: qsort(chroms,n_proper_divs,sizeof(Chrom_T),Chrom_compare_order); break;
  case ALPHA_SORT: qsort(chroms,n_proper_divs,sizeof(Chrom_T),Chrom_compare_alpha); break;
  case NUMERIC_ALPHA_SORT: qsort(chroms,n_proper_divs,sizeof(Chrom_T),Chrom_compare_numeric_alpha); break;
  case CHROM_SORT: qsort(chroms,n_proper_divs,sizeof(Chrom_T),Chrom_compare_chrom); break;
  default: fprintf(stderr,"Don't recognize divsort type %d\n",divsort); abort();
  }
      
  /* The zeroth div is empty */
  divstring = (char *) CALLOC(1,sizeof(char));
  divstring[0] = '\0';
  divlist = List_push(NULL,divstring);

  for (i = 0; i < n_proper_divs; i++) {
    divlist = List_push(divlist,Chrom_string(chroms[i]));
  }
  divlist = List_reverse(divlist);

#if 0
  /* Causes invalid reads later on */
  for (i = 0; i < n_proper_divs; i++) {
    Chrom_free(&(chroms[i]));
  }
#endif

  FREE(chroms);


  if (iit_version == 1) {
    IIT_write_univ(iitfile,divlist,typelist,intervaltable,labeltable,annottable,
		   coord_values_8p,label_pointers_8p,annot_pointers_8p);
  } else if (valuep == false) {
    IIT_write(iitfile,divlist,typelist,fieldlist,intervaltable,/*valuetable*/NULL,labeltable,annottable,
	      divsort,iit_version,label_pointers_8p,annot_pointers_8p);
  } else {
    IIT_write(iitfile,divlist,typelist,fieldlist,intervaltable,valuetable,labeltable,annottable,
	      divsort,iit_version,label_pointers_8p,annot_pointers_8p);
  }
  FREE(iitfile);

  for (d = divlist; d != NULL; d = List_next(d)) {
    divstring = (char *) List_head(d);

    templist = (List_T) Table_get(annottable,(void *) divstring);
    for (l = templist; l != NULL; l = List_next(l)) {
      tempstring = (char *) List_head(l);
      FREE(tempstring);
    }
    List_free(&templist);
      
    valuelist = (Doublelist_T) Table_get(valuetable,(void *) divstring);
    Doublelist_free(&valuelist);

    templist = (List_T) Table_get(labeltable,(void *) divstring);
    for (l = templist; l != NULL; l = List_next(l)) {
      tempstring = (char *) List_head(l);
      FREE(tempstring);
    }
    List_free(&templist);

    templist = (List_T) Table_get(intervaltable,(void *) divstring);
    if (iit_version == 1) {
      for (l = templist; l != NULL; l = List_next(l)) {
	univinterval = (Univinterval_T) List_head(l);
	Univinterval_free(&univinterval);
      }
    } else {
      for (l = templist; l != NULL; l = List_next(l)) {
	interval = (Interval_T) List_head(l);
	Interval_free(&interval);
      }
    }
    List_free(&templist);
      
  }


  Table_free(&intervaltable);
  Table_free(&valuetable);
  Table_free(&labeltable);
  Table_free(&annottable);

  for (l = fieldlist; l != NULL; l = List_next(l)) {
    tempstring = (char *) List_head(l);
    FREE(tempstring);
  }
  List_free(&fieldlist);

  for (l = typelist; l != NULL; l = List_next(l)) {
    tempstring = (char *) List_head(l);
    FREE(tempstring);
  }
  List_free(&typelist);

  for (l = divlist; l != NULL; l = List_next(l)) {
    tempstring = (char *) List_head(l);
    FREE(tempstring);
  }
  List_free(&divlist);

  Tableint_free(&typetable);
  Tableint_free(&div_seenp);

  return 0;
}
