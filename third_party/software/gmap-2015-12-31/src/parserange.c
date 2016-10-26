static char rcsid[] = "$Id: parserange.c 145990 2014-08-25 21:47:32Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "access.h"
#include "parserange.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>		/* For isdigit */

#include "mem.h"
#include "univinterval.h"
#include "interval.h"
#include "iit-read-univ.h"
#include "iit-read.h"


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


/* Note that isnumber is a function in ctype.h on some systems */
bool
Parserange_iscoordp (Univcoord_T *result, char *string) {
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

bool
Parserange_islengthp (Chrpos_T *result, char *string) {
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

/* Returns coordinates as zero-based */
bool
Parserange_israngep (Univcoord_T *left, Chrpos_T *length, bool *revcomp, char *string) {
  bool result;
  char *copy, *startstring, *endstring;
  Univcoord_T start, end;

  copy = (char *) MALLOCA((strlen(string)+1) * sizeof(char));
  strcpy(copy,string);

  if (index(copy,'.')) {
    startstring = strtok(copy,"..");
    endstring = strtok(NULL,"..");
    if (!Parserange_iscoordp(&start,startstring) || !Parserange_iscoordp(&end,endstring)) {
      result = false;
    } else if (start <= end) {
      *length = end - start + 1;
      *left = start - 1;
      *revcomp = false;
      debug(printf("(..) "));
      result = true;
    } else {
      *length = start - end + 1;
      *left = end - 1;
      *revcomp = true;
      debug(printf("(..) "));
      result = true;
    }

  } else if (index(copy,'+')) {
    startstring = strtok(copy,"+");
    endstring = strtok(NULL,"+");
    if (!Parserange_iscoordp(&start,startstring)) {
      result = false;
    } else if (endstring[0] == '-' && Parserange_islengthp(&(*length),&(endstring[1]))) {
      *left = start - (*length);
      *revcomp = true;
      debug(printf("(-) "));
      result = true;
    } else if (!Parserange_islengthp(&(*length),endstring)) {
      result = false;
    } else {
      *left = start - 1;
      *revcomp = false;
      debug(printf("(+) "));
      result = true;
    }

  } else if (index(copy,'-')) {
    /* Old notation */
    startstring = strtok(copy,"--");
    endstring = strtok(NULL,"--");
    if (!Parserange_iscoordp(&start,startstring) || !Parserange_iscoordp(&end,endstring)) {
      result = false;
    } else if (start <= end) {
      *length = end - start + 1;
      *left = start - 1;
      *revcomp = false;
      debug(printf("(--) "));
      result = true;
    } else {
      *length = start - end + 1;
      *left = end - 1;
      *revcomp = true;
      debug(printf("(--) "));
      result = true;
    }

    /* Don't allow this yet ...
  } else if (index(copy,'-')) {
    startstring = strtok(copy,"-");
    endstring = strtok(NULL,"-");
    if (!Parserange_iscoordp(&start,startstring) || !Parserange_iscoordp(&end,endstring)) {
      result = false;
    } else if (end > start - 1) {
      result = false;
    } else {
      *left = start - 1 - end;
      *length = end;
      *revcomp = true;
      result = true;
    }
    */
    
  } else {
    result = false;
  }

  FREEA(copy);

  return result;
}


/* Retrieval functions */

static int
translate_chromosomepos_universal (Univcoord_T *genomicstart, Chrpos_T *genomiclength, 
				   char *chromosome, Univcoord_T left, Chrpos_T length,
				   Univ_IIT_T chromosome_iit) {
  int rc = 1, index;
  Univinterval_T interval;
#ifdef DEBUG
  bool allocp;
#endif
  
  if ((index = Univ_IIT_find_linear(chromosome_iit,chromosome)) >= 0) {
    debug(printf("chromosome %s => index %d\n",chromosome,index));
    interval = Univ_IIT_interval(chromosome_iit,index);
    debug(printf("  => label %s with interval low %u\n",
		 Univ_IIT_label(chromosome_iit,index,&allocp),Interval_low(interval)));
    *genomicstart = Univinterval_low(interval)+left;
    if (*genomicstart < Univinterval_low(interval)) {
      fprintf(stderr,"%lu + %lu = %lu (exceeds a 32-bit unsigned int)\n",
	      Univinterval_low(interval),left,*genomicstart);
      exit(9);
    }
    if (length == 0) {
      *genomiclength = Univinterval_length(interval)	/* - left (Why?) */; 
    } else {
      *genomiclength = length;
    }
    rc = 0;
  }
  
  return rc;
}


static int
translate_chromosomepos_segment (Univcoord_T *segmentstart, Chrpos_T *segmentlength, 
				 char *chromosome, Univcoord_T left, Chrpos_T length,
				 Univ_IIT_T chromosome_iit) {
  int rc = 1, index;
  Univinterval_T interval;
  
  if ((index = Univ_IIT_find_linear(chromosome_iit,chromosome)) >= 0) {
    interval = Univ_IIT_interval(chromosome_iit,index);
    *segmentstart = left;
    if (length == 0) {
      *segmentlength = Univinterval_length(interval) /* - left (Why?) */;
    } else {
      *segmentlength = length;
    }
    rc = 0;
  }
  
  return rc;
}


static int
translate_contig (Univcoord_T *genomicstart, Chrpos_T *genomiclength,
		  char *contig, Univcoord_T left, Chrpos_T length, Univ_IIT_T contig_iit) {
  int rc = 1, index;
  Univinterval_T interval;
  
  if ((index = Univ_IIT_find_one(contig_iit,contig)) >= 0) {
    interval = Univ_IIT_interval(contig_iit,index);
    *genomicstart = Univinterval_low(interval)+left;
    if (length == 0) {
      *genomiclength = Univinterval_length(interval) /* - left (Why?) */;
    } else {
      *genomiclength = length;
    }
    rc = 0;
  }

  return rc;
}


/* Assumes position is 0-based */
static char *
convert_to_chrpos (Chrpos_T *chrpos, char *genomesubdir, char *fileroot, Univcoord_T position) {
  char *chromosome, *filename;
  Univ_IIT_T chromosome_iit;
  
  filename = (char *) MALLOCA((strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
			       strlen(".chromosome.iit")+1) * sizeof(char));
  sprintf(filename,"%s/%s.chromosome.iit",genomesubdir,fileroot);
  if ((chromosome_iit = Univ_IIT_read(filename,/*readonlyp*/true,/*add_iit_p*/false)) == NULL) {
    fprintf(stderr,"Can't read IIT file %s\n",filename);
    exit(9);
  }
  FREEA(filename);

  /* Subtract 1 to make 0-based */
  chromosome = Univ_IIT_string_from_position(&(*chrpos),position-1U,chromosome_iit);
  *chrpos += 1U;		/* Bring back to 1-based */
  Univ_IIT_free(&chromosome_iit);

  return chromosome;
}

/* Assumes position is 0-based */
static char *
convert_to_chrpos_iit (Chrpos_T *chrpos, Univ_IIT_T chromosome_iit, Univcoord_T position) {
  char *chromosome;
  
  /* Subtract 1 to make 0-based */
  chromosome = Univ_IIT_string_from_position(&(*chrpos),position-1U,chromosome_iit);
  *chrpos += 1U;		/* Bring back to 1-based */
  Univ_IIT_free(&chromosome_iit);

  return chromosome;
}


static char *
find_div (int *div_strlen, char *string, int sep) {
  *div_strlen = 0;

  while (*string != '\0') {
    if (*string == sep) {
      return &(string[1]);
    } else {
      (*div_strlen) += 1;
      string++;
    }
  }

  return (char *) NULL;
}


bool
Parserange_query (char **divstring, Univcoord_T *coordstart, Univcoord_T *coordend, bool *revcomp,
		  char *query, char *filename) {
  char *coords;
  Univcoord_T result, left;
  Chrpos_T length;
  int div_strlen;
  IIT_T iit;
  
  *divstring = NULL;
  *revcomp = false;

  if ((coords = find_div(&div_strlen,query,':')) != NULL) {
    /* Query may have a div */
    *divstring = (char *) MALLOC((div_strlen+1) * sizeof(char)); /* Return value */
    strncpy(*divstring,query,div_strlen);

    debug(printf("Parsed query %s into divstring %s and coords %s\n",
		 query,*divstring,coords));

    if (IIT_read_divint(filename,*divstring,/*add_iit_p*/true) < 0) {
      fprintf(stderr,"Chromosome %s not found in IIT file\n",*divstring);
      debug(printf("  but divstring not found, so treat as label\n"));
      FREE(*divstring);	/* free only when returning false */
      return false;
    } else if (coords == NULL || *coords == '\0') {
      debug(printf("  entire div\n"));
      if ((iit = IIT_read(filename,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ONE,*divstring,
			  /*add_iit_p*/true,/*labels_read_p*/false)) == NULL) {
	if (Access_file_exists_p(filename) == false) {
	  fprintf(stderr,"Cannot read file %s\n",filename);
	} else {
	  fprintf(stderr,"File %s appears to be an invalid IIT file\n",filename);
	}
	exit(9);
      } else {
	*coordstart= 0;
	*coordend = IIT_divlength(iit,*divstring);
	debug(printf("  divlength is %u\n",*coordend));
	IIT_free(&iit);
      }
      return true;
    } else if (Parserange_iscoordp(&result,coords)) {
      debug(printf("  and coords %s as a number\n",coords));
      *coordstart = result;
      *coordend = result;
      return true;
    } else if (Parserange_israngep(&left,&length,&(*revcomp),coords)) {
      debug(printf("  and coords %s as a range starting at %u with length %u and revcomp = %d\n",
		   coords,left,length,*revcomp));
      *coordstart = left + 1;	/* Because isrange is 0-based */
      *coordend = left + length;
      return true;
    } else {
      debug(printf("  but coords %s is neither a number nor a range.  Interpret as a label.\n",coords));
      FREE(*divstring);		/* free only when returning false */
      return false;
    }

  } else {
    /* No div.  Query must be a number, range, or label */

    debug(printf("Parsed query %s without a div ",query));
    if (Parserange_iscoordp(&result,query)) {
      debug(printf("number\n"));
      *coordstart = result;
      *coordend = result;
      return true;
    } else if (Parserange_israngep(&left,&length,&(*revcomp),query)) {
      debug(printf("range\n"));
      *coordstart = left + 1;	/* Because isrange is 0-based */
      *coordend = left + length;
      return true;
    } else {
      debug(printf("label\n"));
      return false;
    }
  }
}



bool
Parserange_universal (char **div, bool *revcomp,
		      Univcoord_T *genomicstart, Chrpos_T *genomiclength,
		      Chrpos_T *chrstart, Chrpos_T *chrend,
		      Univcoord_T *chroffset, Chrpos_T *chrlength,
		      char *query, char *genomesubdir, char *fileroot) {
  char *coords, *filename;
  Univcoord_T result, left;
  Chrpos_T length;
  Univ_IIT_T chromosome_iit, contig_iit;
  Univinterval_T interval;
  int theindex;
  int rc;
  
  *revcomp = false;
  if (index(query,':')) {
    /* Segment must be a genome, chromosome, or contig */
    debug(printf("Parsed query %s into ",query));
    *div = strtok(query,":");
    if ((*div)[0] == '+') {
      *revcomp = false;
      *div = &((*div)[1]);
    } else if ((*div)[0] == '_') {
      *revcomp = true;
      *div = &((*div)[1]);
    }
    coords = strtok(NULL,":");
    debug(printf("segment %s and coords %s\n",*div,coords));

    /* Try chromosome first */
    filename = (char *) MALLOCA((strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
				 strlen(".chromosome.iit")+1) * sizeof(char));
    sprintf(filename,"%s/%s.chromosome.iit",genomesubdir,fileroot);
    chromosome_iit = Univ_IIT_read(filename,/*readonlyp*/true,/*add_iit_p*/false);
    FREEA(filename);

    debug(printf("Interpreting segment %s as a chromosome\n",*div));
    if (coords == NULL) {
      debug(printf("  entire chromosome\n"));
      rc = translate_chromosomepos_universal(&(*genomicstart),&(*genomiclength),*div,left=0,length=0,chromosome_iit);
    } else if (Parserange_iscoordp(&result,coords)) {
      debug(printf("  and coords %s as a number\n",coords));
      rc = translate_chromosomepos_universal(&(*genomicstart),&(*genomiclength),*div,left=result-1,length=1,chromosome_iit);
    } else if (Parserange_israngep(&left,&length,&(*revcomp),coords)) {
      debug(printf("  and coords %s as a range starting at %u with length %u and revcomp = %d\n",
		   coords,left,length,*revcomp));
      rc = translate_chromosomepos_universal(&(*genomicstart),&(*genomiclength),*div,left,length,chromosome_iit);
    } else {
      debug(printf("  but coords %s is neither a number nor a range\n",coords));
      rc = -1;
    }

    /* Compute chromosomal coordinates */
    *chrstart = left;
    *chrend = *chrstart + *genomiclength;
    *chrstart += 1U;		/* Make 1-based */

    /* Get chromosomal information */
    if ((theindex = Univ_IIT_find_one(chromosome_iit,*div)) < 0) {
      fprintf(stderr,"Cannot find chromosome %s in chromosome IIT file\n",*div);
      /* exit(9); */
    } else {
      interval = Univ_IIT_interval(chromosome_iit,theindex);
      *chroffset = Univinterval_low(interval);
      *chrlength = Univinterval_length(interval);
    }

    Univ_IIT_free(&chromosome_iit);

#if 0
    /* Contig IIT's are of type 1, which require some work to compute
       on current div-based scheme.  Just abandoning for now. */
    if (rc != 0) {
      /* Try contig */
      filename = (char *) MALLOCA((strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
				   strlen(".contig.iit")+1) * sizeof(char));
      sprintf(filename,"%s/%s.contig.iit",genomesubdir,fileroot);
      contig_iit = IIT_read(filename,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
			    /*divstring*/NULL,/*add_iit_p*/false,/*labels_read_p*/true);
      FREEA(filename);

      debug(printf("Interpreting segment %s as a contig\n",*div));
      if (coords == NULL) {
	debug(printf("  entire contig\n"));
	rc = translate_contig_universal(&(*genomicstart),&(*genomiclength),*div,left=0,length=0,chromosome_iit);
      } else if (Parserange_iscoordp(&result,coords)) {
	debug(printf("  and coords %s as a number\n",coords));
	rc = translate_contig(&(*genomicstart),&(*genomiclength),*div,left=result-1,length=1,contig_iit);
      } else if (Parserange_israngep(&left,&length,&(*revcomp),coords)) {
	debug(printf("  and coords %s as a range starting at %u with length %u and revcomp = %d\n",
		     coords,left,length,*revcomp));
	rc = translate_contig(&(*genomicstart),&(*genomiclength),*div,left,length,contig_iit);
      } else {
	debug(printf("  but coords %s is neither a number nor a range\n",coords));
	rc = -1;
      }

      IIT_free(&contig_iit);
    }
#endif

    if (rc != 0) {
      fprintf(stderr,"Can't find coordinates %s:%s\n",*div,coords);
      return false;
    } else {
      return true;
    }

  } else {
    /* Query must be a genomic position, genomic range, or contig */
    *chrstart = *chroffset = *chrlength = 0;

    debug(printf("Parsed query %s as atomic ",query));
    if (Parserange_iscoordp(&result,query)) {
      debug(printf("number\n"));
      *genomicstart = result-1;
      *genomiclength = 1;

    } else if (Parserange_israngep(&left,&length,&(*revcomp),query)) {
      debug(printf("range\n"));
      *genomicstart = left;
      *genomiclength = length;

    } else {

      debug(printf("contig\n"));
      return false;
#if 0
      filename = (char *) MALLOCA((strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
				   strlen(".contig.iit")+1) * sizeof(char));
      sprintf(filename,"%s/%s.contig.iit",genomesubdir,fileroot);
      contig_iit = IIT_read(filename,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
			    /*divstring*/NULL,/*add_iit_p*/false,/*labels_read_p*/true);
      FREEA(filename);

      rc = translate_contig(&(*genomicstart),&(*genomiclength),query,left=0,length=0,contig_iit);
      IIT_free(&contig_iit);
#endif
    }

    *div = convert_to_chrpos(&(*chrstart),genomesubdir,fileroot,*genomicstart);
    *chrend = *chrstart + *genomiclength;
    *chrstart += 1U;		/* Make 1-based */

    /* Try chromosome first */
    filename = (char *) MALLOCA((strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
				 strlen(".chromosome.iit")+1) * sizeof(char));
    sprintf(filename,"%s/%s.chromosome.iit",genomesubdir,fileroot);
    chromosome_iit = Univ_IIT_read(filename,/*readonlyp*/true,/*add_iit_p*/false);
    FREEA(filename);

    if ((theindex = Univ_IIT_find_one(chromosome_iit,*div)) < 0) {
      fprintf(stderr,"Cannot find chromosome %s in chromosome IIT file\n",*div);
      Univ_IIT_free(&chromosome_iit);
      return false;
    } else {
      interval = Univ_IIT_interval(chromosome_iit,theindex);
      *chroffset = Univinterval_low(interval);
      *chrlength = Univinterval_length(interval);
      Univ_IIT_free(&chromosome_iit);
      return true;
    }
  }
}


bool
Parserange_universal_iit (char **div, bool *revcomp,
			  Univcoord_T *genomicstart, Chrpos_T *genomiclength,
			  Chrpos_T *chrstart, Chrpos_T *chrend,
			  Univcoord_T *chroffset, Chrpos_T *chrlength,
			  char *query, Univ_IIT_T chromosome_iit, Univ_IIT_T contig_iit) {
  char *coords;
  Univcoord_T result, left;
  Chrpos_T length;
  Univinterval_T interval;
  int theindex;
  int rc;
  
  *revcomp = false;
  if (index(query,':')) {
    /* Segment must be a genome, chromosome, or contig */
    debug(printf("Parsed query %s into ",query));
    *div = strtok(query,":");
    if ((*div)[0] == '+') {
      *revcomp = false;
      *div = &((*div)[1]);
    } else if ((*div)[0] == '_') {
      *revcomp = true;
      *div = &((*div)[1]);
    }
    coords = strtok(NULL,":");
    debug(printf("segment %s and coords %s\n",*div,coords));


    debug(printf("Interpreting segment %s as a chromosome\n",*div));
    if (coords == NULL) {
      debug(printf("  entire chromosome\n"));
      rc = translate_chromosomepos_universal(&(*genomicstart),&(*genomiclength),*div,left=0,length=0,chromosome_iit);
    } else if (Parserange_iscoordp(&result,coords)) {
      debug(printf("  and coords %s as a number\n",coords));
      rc = translate_chromosomepos_universal(&(*genomicstart),&(*genomiclength),*div,left=result-1,length=1,chromosome_iit);
    } else if (Parserange_israngep(&left,&length,&(*revcomp),coords)) {
      debug(printf("  and coords %s as a range starting at %u with length %u and revcomp = %d\n",
		   coords,left,length,*revcomp));
      rc = translate_chromosomepos_universal(&(*genomicstart),&(*genomiclength),*div,left,length,chromosome_iit);
    } else {
      debug(printf("  but coords %s is neither a number nor a range\n",coords));
      rc = -1;
    }

    /* Compute chromosomal coordinates */
    *chrstart = left;
    *chrend = *chrstart + *genomiclength;
    *chrstart += 1U;		/* Make 1-based */

    /* Get chromosomal information */
    if ((theindex = Univ_IIT_find_one(chromosome_iit,*div)) < 0) {
      fprintf(stderr,"Cannot find chromosome %s in chromosome IIT file\n",*div);
      /* exit(9); */
    } else {
      interval = Univ_IIT_interval(chromosome_iit,theindex);
      *chroffset = Univinterval_low(interval);
      *chrlength = Univinterval_length(interval);
    }

    if (rc != 0) {
      /* Try contig */
      debug(printf("Interpreting segment %s as a contig\n",*div));
      if (coords == NULL) {
	debug(printf("  entire contig\n"));
	rc = translate_contig(&(*genomicstart),&(*genomiclength),*div,left=0,length=0,contig_iit);
      } else if (Parserange_iscoordp(&result,coords)) {
	debug(printf("  and coords %s as a number\n",coords));
	rc = translate_contig(&(*genomicstart),&(*genomiclength),*div,left=result-1,length=1,contig_iit);
      } else if (Parserange_israngep(&left,&length,&(*revcomp),coords)) {
	debug(printf("  and coords %s as a range starting at %u with length %u and revcomp = %d\n",
		     coords,left,length,*revcomp));
	rc = translate_contig(&(*genomicstart),&(*genomiclength),*div,left,length,contig_iit);
      } else {
	debug(printf("  but coords %s is neither a number nor a range\n",coords));
	rc = -1;
      }
    }

    if (rc != 0) {
      fprintf(stderr,"Can't find coordinates %s:%s\n",*div,coords);
      return false;
    } else {
      return true;
    }

  } else {
    /* Query must be a genomic position, genomic range, or contig */
    *chrstart = *chroffset = *chrlength = 0;

    debug(printf("Parsed query %s as atomic ",query));
    if (Parserange_iscoordp(&result,query)) {
      debug(printf("number\n"));
      *genomicstart = result-1;
      *genomiclength = 1;

    } else if (Parserange_israngep(&left,&length,&(*revcomp),query)) {
      debug(printf("range\n"));
      *genomicstart = left;
      *genomiclength = length;

    } else {

      debug(printf("contig\n"));
      return false;
#if 0
      rc = translate_contig(&(*genomicstart),&(*genomiclength),query,left=0,length=0,contig_iit);
      IIT_free(&contig_iit);
#endif
    }

    *div = convert_to_chrpos_iit(&(*chrstart),chromosome_iit,*genomicstart);
    *chrend = *chrstart + *genomiclength;
    *chrstart += 1U;		/* Make 1-based */

    /* Try chromosome first */
    if ((theindex = Univ_IIT_find_one(chromosome_iit,*div)) < 0) {
      fprintf(stderr,"Cannot find chromosome %s in chromosome IIT file\n",*div);
      return false;
    } else {
      interval = Univ_IIT_interval(chromosome_iit,theindex);
      *chroffset = Univinterval_low(interval);
      *chrlength = Univinterval_length(interval);
      return true;
    }
  }
}



bool
Parserange_simple (char **div, bool *revcomp, Chrpos_T *chrstart, Chrpos_T *chrend,
		   char *query) {
  char *coords;
  Univcoord_T result, left;
  Chrpos_T length;
  
  *revcomp = false;
  if (index(query,':')) {
    /* Segment must be a genome, chromosome, or contig */
    debug(printf("Parsed query %s into ",query));
    *div = strtok(query,":");
    if ((*div)[0] == '+') {
      *revcomp = false;
      *div = &((*div)[1]);
    } else if ((*div)[0] == '_') {
      *revcomp = true;
      *div = &((*div)[1]);
    }
    coords = strtok(NULL,":");
    debug(printf("segment %s and coords %s\n",*div,coords));

    debug(printf("Interpreting segment %s as a chromosome\n",*div));
    if (coords == NULL) {
      fprintf(stderr,"Need region after ':'\n");
      return false;
    } else if (Parserange_iscoordp(&result,coords)) {
      debug(printf("  and coords %s as a number\n",coords));
      left = result - 1;	/* Make 0-based */
      length = 1;
    } else if (Parserange_israngep(&left,&length,&(*revcomp),coords)) {
      debug(printf("  and coords %s as a range starting at %u with length %u and revcomp = %d\n",
		   coords,left,length,*revcomp));
    } else {
      fprintf(stderr,"Coordinates after ':' is neither a number nor a range\n");
      debug(printf("  but coords %s is neither a number nor a range\n",coords));
      return false;
    }

    /* Compute chromosomal coordinates */
    *chrstart = left;
    *chrend = *chrstart + length;
    *chrstart += 1U;		/* Make 1-based */

    return true;

  } else {
    fprintf(stderr,"Region %s does not contain ':'\n",query);
    return false;
  }
}

