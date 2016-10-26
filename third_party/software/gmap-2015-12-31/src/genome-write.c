static char rcsid[] = "$Id: genome-write.c 168395 2015-06-26 17:13:13Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "genome-write.h"

#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
#else
#include "littleendian.h"
#endif

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>		/* For rindex */
#include <ctype.h>		/* For isspace */
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>		/* For off_t */
#endif

#include "mem.h"
#include "types.h"
#include "fopen.h"
#include "interval.h"
#include "compress-write.h"
#include "iit-write.h"
#include "complement.h"
#include "genome.h"		/* For Genome_uncompress_memory */

#define CONTROLM 13		/* From PC */

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


/* Compressed genome file format: 32 characters in genome => 3
   bits/character = 96 bits, or 12 bytes, or 3 unsigned ints.  The
   first unsigned int (32 bits) has the high bits for these 32
   characters, the second unsigned int has the low bits, and the third
   unsigned int has the flags */

static void
find_positions (bool *revcompp, Univcoord_T *leftposition, Univcoord_T *rightposition,
		Univcoord_T *startposition, Univcoord_T *endposition, Chrpos_T *truelength,
		int *contigtype, char *accession, Univ_IIT_T contig_iit) {
  int index;
  Univinterval_T interval;
  char firstchar;

  if ((index = Univ_IIT_find_one(contig_iit,accession)) == -1) {
    fprintf(stderr,"Can't find accession %s in contig IIT file\n",
	    accession);
    exit(9);
  } else {	
    interval = Univ_IIT_interval(contig_iit,index);

    *leftposition = Univinterval_low(interval);
    *rightposition = Univinterval_high(interval);

#if 1
    /* contig IIT now always written with version 1 (universal) */
    firstchar = Univ_IIT_annotation_firstchar(contig_iit,index);
    if (firstchar == '-') {
      *revcompp = true;
      *startposition = Univinterval_high(interval) + 1U;
      *endposition = Univinterval_low(interval) + 1U;
    } else {
      *revcompp = false;
      *startposition = Univinterval_low(interval);
      *endposition = Univinterval_high(interval);
    }
#else
    /* Code for version 2 IITs */
    if (Interval_sign(interval) < 0) {
      *revcompp = true;
      *startposition = Interval_high(interval) + 1U;
      *endposition = Interval_low(interval) + 1U;
    } else {
      *revcompp = false;
      *startposition = Interval_low(interval);
      *endposition = Interval_high(interval);
    }
#endif

    *truelength = Univinterval_length(interval);
    *contigtype = Univinterval_type(interval);
    debug(fprintf(stderr,"revcompp = %d, leftposition = %llu, rightposition = %llu, startposition = %llu, endposition = %llu\n",
		  *revcompp,(unsigned long long) *leftposition,(unsigned long long) *rightposition,
		  (unsigned long long) *startposition,(unsigned long long) *endposition));

    return;
  }
}

static void
move_absolute (FILE *fp, Univcoord_T offset) {

#ifdef HAVE_FSEEKO
  if (fseeko(fp,(off_t) offset,SEEK_SET) < 0) {
    perror("Error in gmapindex, fseeko");
    exit(9);
  }
#else
  if (fseek(fp,(long) offset,SEEK_SET) < 0) {
    perror("Error in gmapindex, fseek");
    exit(9);
  }
#endif

  return;
}



#define WRITEBLOCK 1024
static char Empty[WRITEBLOCK];

static void
fill_zero (FILE *fp, Univcoord_T startpos, Univcoord_T endpos, bool uncompressedp,
	   int index1part) {
  Univcoord_T total;

  if (uncompressedp == true) {
    move_absolute(fp,startpos);
    total = endpos - startpos;
    while (total > WRITEBLOCK) {
      fwrite(Empty,sizeof(char),WRITEBLOCK,fp);
      total -= WRITEBLOCK;
    }
    if (total > 0) {
      fwrite(Empty,sizeof(char),total,fp);
    }
  } else {
    Compress_update_file(/*nbadchars*/0,fp,NULL,startpos,endpos,index1part);
  }
  return;
}

static void
fill_x (FILE *fp, Univcoord_T startpos, Univcoord_T endpos, bool uncompressedp,
	int index1part) {
  Univcoord_T total;

  if (uncompressedp == true) {
    move_absolute(fp,startpos);
    total = endpos - startpos;
    while (total > WRITEBLOCK) {
      fwrite(Empty,sizeof(char),WRITEBLOCK,fp);
      total -= WRITEBLOCK;
    }
    if (total > 0) {
      fwrite(Empty,sizeof(char),total,fp);
    }
  } else {
    Compress_update_file(/*nbadchars*/0,fp,NULL,startpos,endpos,index1part);
  }
  return;
}

static void
fill_x_memory (Genomecomp_T *genomecomp, Univcoord_T startpos, Univcoord_T endpos) {

  Compress_update_memory(/*nbadchars*/0,genomecomp,NULL,startpos,endpos);
  return;
}


#define BUFFERSIZE 8192

static char *
parse_accession (char *Header) {
  char *accession, *p;
  int strlength;

  p = &(Header[1]);		/* First character is '>' */
  while (*p != '\0' && *p != '\n' && *p != CONTROLM && !isspace((int) *p)) {
    p++;
  }
  *p = '\0';
  strlength = (p - &(Header[1]))/sizeof(char);
  accession = (char *) CALLOC(strlength+1,sizeof(char));
  strcpy(accession,&(Header[1]));

  return accession;
}

static void
make_reverse_buffered (char *reverse, char *sequence, Chrpos_T length) {
  Chrpos_T i, j;

  /* complement = (char *) CALLOC(length+1,sizeof(char)); */
  for (i = length-1, j = 0; j < length; i--, j++) {
    reverse[j] = sequence[i];
  }
  reverse[length] = '\0';
  return;
}

static void
make_complement_buffered (char *complement, char *sequence, Chrpos_T length) {
  char complCode[128] = COMPLEMENT_LC;
  Chrpos_T i, j;

  /* complement = (char *) CALLOC(length+1,sizeof(char)); */
  for (i = length-1, j = 0; j < length; i--, j++) {
    complement[j] = complCode[(int) sequence[i]];
  }
  complement[length] = '\0';
  return;
}

/* Puts reference genome into refgenome_fp (potentially compressed!),
   and puts alternate strain sequences into altstrain_iit. */
static void
genome_write_file (FILE *refgenome_fp, FILE *input, 
		   Univ_IIT_T contig_iit, IIT_T altstrain_iit, char *fileroot, bool uncompressedp,
		   int index1part, int nmessages) {
  char Buffer[BUFFERSIZE], Complement[BUFFERSIZE], *segment;
  char *accession, *p;
  Univcoord_T leftposition, rightposition, startposition, endposition,
    maxposition = 0, currposition = 0;
  Chrpos_T truelength;
  int contigtype;
  int i;
  bool revcompp;
#ifdef ALTSTRAIN
  int altstrain_index, altstrain_offset;
#endif
  int nbadchars = 0;
  int ncontigs = 0;

  for (i = 0; i < WRITEBLOCK; i++) {
    Empty[i] = 'X';
  }

  while (fgets(Buffer,BUFFERSIZE,input) != NULL) {
    if (Buffer[0] == '>') {
      /* HEADER */
      accession = parse_accession(Buffer);
      find_positions(&revcompp,&leftposition,&rightposition,&startposition,&endposition,
		     &truelength,&contigtype,accession,contig_iit);
      if (++ncontigs < nmessages) {
	if (revcompp == true) {
	  fprintf(stderr,"Writing contig %s to universal coordinates %llu..%llu in genome %s\n",
		  accession,(unsigned long long) startposition,(unsigned long long) endposition,fileroot);
	} else {
	  fprintf(stderr,"Writing contig %s to universal coordinates %llu..%llu in genome %s\n",
		  accession,(unsigned long long) startposition+1U,(unsigned long long) endposition+1U,fileroot);
	}
      } else if (ncontigs == nmessages) {
	fprintf(stderr,"More than %d contigs.  Will stop printing messages\n",nmessages);
      }

      if (contigtype > 0) {
#ifdef ALTSTRAIN
	fprintf(stderr," (alternate strain %s)",IIT_typestring(altstrain_iit,contigtype));
#endif
      }

      FREE(accession);

      if (contigtype > 0) {
#ifdef ALTSTRAIN
	/* Initialize file pointer for alternate strain */
	altstrain_index = IIT_get_exact(altstrain_iit,/*divstring*/NULL,leftposition,rightposition,contigtype);
	if (revcompp == true) {
	  altstrain_offset = rightposition + 1U - leftposition;
	} else {
	  altstrain_offset = 0;
	}
	debug(fprintf(stderr,"Setting altstrain_offset to be %d\n",altstrain_offset));
#endif
      }

      /* Handles case where true length is greater than provided
         coordinates.  This needs to be after call to IIT_get_exact */
      if (leftposition + truelength - 1U > rightposition) {
	debug(fprintf(stderr,"Extending endposition for truelength of %u\n",truelength));
	rightposition = leftposition + truelength;
	if (revcompp == true) {
	  endposition = startposition - truelength;
	} else {
	  endposition = startposition + truelength;
	}
      }

      /* In both cases, set file pointer for reference strain,
         although we won't write sequence of alternate strain.  For an
         alternate strain, ensure that we fill the reference strain
         with sufficient X's. */
      if (startposition > maxposition) {
	/* Start beyond end of file */
	debug(fprintf(stderr,"Filling with X's from %llu to %llu-1\n",
		      (unsigned long long) maxposition,(unsigned long long) startposition));
	fill_x(refgenome_fp,maxposition,startposition,uncompressedp,index1part);
	  
	if (contigtype > 0) {
#ifdef ALTSTRAIN
	  fill_x(refgenome_fp,leftposition,rightposition + 1,uncompressedp,index1part);
	  maxposition = currposition = rightposition + 1;
#endif
	} else {
	  maxposition = rightposition;
	  currposition = startposition;
	}

      } else {
	/* Start within file */
	if (contigtype > 0) {
#ifdef ALTSTRAIN
	  if (rightposition + 1 > maxposition) {
	    debug(fprintf(stderr,"Filling with X's from %llu to %llu-1\n",
			  (unsigned long long) maxposition,(unsigned long long) rightposition+1));
	    fill_x(refgenome_fp,maxposition,rightposition + 1,uncompressedp,index1part);
	    maxposition = currposition = rightposition + 1;
	  }
#endif
	} else {
	  debug(fprintf(stderr,"Moving to %llu\n",(unsigned long long) startposition));
	  if (uncompressedp == true) {
	    move_absolute(refgenome_fp,startposition);
	  }
	  currposition = startposition;
	}
      }

    } else {
      /* SEQUENCE */
      if ((p = rindex(Buffer,'\n')) != NULL) {
	*p = '\0';
      }
      if ((p = rindex(Buffer,CONTROLM)) != NULL) {
	*p = '\0';
      }
      if (revcompp == true) {
	make_complement_buffered(Complement,Buffer,strlen(Buffer));
	segment = Complement;
      } else {
	segment = Buffer;
      }

      if (contigtype > 0) {
#ifdef ALTSTRAIN
	/* Write alternate strain */
	if (revcompp == true) {
	  altstrain_offset -= strlen(segment);
	  debug(fprintf(stderr,"Writing alternate strain at %d\n",altstrain_offset));
	  IIT_backfill_sequence(altstrain_iit,altstrain_index,altstrain_offset,segment);
	} else {
	  debug(fprintf(stderr,"Writing alternate strain at %d\n",altstrain_offset));
	  IIT_backfill_sequence(altstrain_iit,altstrain_index,altstrain_offset,segment);
	  altstrain_offset += strlen(segment);
	}
#endif
      } else {
	/* Write reference strain */
	if (revcompp == true) {
	  debug(fprintf(stderr,"Filling with sequence from %llu-1 to %llu\n",
			(unsigned long long) currposition,(unsigned long long) (currposition-strlen(segment))));
	  currposition -= (Univcoord_T) strlen(segment);
	  if (uncompressedp == true) {
	    fwrite(segment,sizeof(char),strlen(segment),refgenome_fp);
	  } else {
	    nbadchars = Compress_update_file(nbadchars,refgenome_fp,segment,currposition,currposition + (Univcoord_T) strlen(segment),
					     index1part);
	  }
	} else {
	  debug(fprintf(stderr,"Filling with sequence from %llu to %llu-1\n",
			(unsigned long long) currposition,(unsigned long long) (currposition+strlen(segment))));
	  if (uncompressedp == true) {
	    fwrite(segment,sizeof(char),strlen(segment),refgenome_fp);
	  } else {
	    nbadchars = Compress_update_file(nbadchars,refgenome_fp,segment,currposition,currposition + (Univcoord_T) strlen(segment),
					     index1part);
	  }
	  currposition += (Univcoord_T) strlen(segment);
	  if (currposition > maxposition) {
	    maxposition = currposition;
	  }
	}
      }
    }
  }

  fprintf(stderr,"A total of %d non-ACGTNX characters were seen in the genome.\n",nbadchars);

  return;
}

/* Permits arbitrary ASCII characters.  Useful for storing numeric data */
static void
genome_writeraw_file (FILE *refgenome_fp, FILE *input, 
		      Univ_IIT_T contig_iit, IIT_T altstrain_iit, char *fileroot, int index1part,
		      int nmessages) {
  char Buffer[BUFFERSIZE], Reverse[BUFFERSIZE], *segment;
  char *accession, c;
  Univcoord_T leftposition, rightposition, startposition, endposition,
    maxposition = 0, currposition = 0;
  Chrpos_T truelength;
  int contigtype;
  int strlength;
  int i;
  bool revcompp;
#ifdef ALTSTRAIN
  int altstrain_index, altstrain_offset;
#endif
  int ncontigs = 0;

  for (i = 0; i < WRITEBLOCK; i++) {
    Empty[i] = 0;
  }

  while (fgets(Buffer,BUFFERSIZE,input) != NULL) {
    if (Buffer[0] != '>') {
      fprintf(stderr,"Expecting to see a new FASTA entry\n");
      fprintf(stderr,"Instead, saw %d (%c)\n",(int) Buffer[0],Buffer[0]);
      exit(9);
    } else {
      /* HEADER */
      accession = parse_accession(Buffer);
      find_positions(&revcompp,&leftposition,&rightposition,&startposition,&endposition,
		     &truelength,&contigtype,accession,contig_iit);
      if (++ncontigs < nmessages) {
	if (revcompp == true) {
	  fprintf(stderr,"Writing contig %s to universal coordinates %llu..%llu in genome %s\n",
		  accession,(unsigned long long) startposition,(unsigned long long) endposition,fileroot);
	} else {
	  fprintf(stderr,"Writing contig %s to universal coordinates %llu..%llu in genome %s\n",
		  accession,(unsigned long long) startposition+1U,(unsigned long long) endposition+1U,fileroot);
	}
      } else if (ncontigs == nmessages) {
	fprintf(stderr,"More than %d contigs.  Will stop printing messages\n",nmessages);
      }

      if (contigtype > 0) {
#ifdef ALTSTRAIN
	fprintf(stderr," (alternate strain %s)",IIT_typestring(altstrain_iit,contigtype));
#endif
      }
      FREE(accession);

      if (contigtype > 0) {
#ifdef ALTSTRAIN
	/* Initialize file pointer for alternate strain */
	altstrain_index = IIT_get_exact(altstrain_iit,/*divstring*/NULL,leftposition,rightposition,contigtype);
	if (revcompp == true) {
	  altstrain_offset = rightposition + 1U - leftposition;
	} else {
	  altstrain_offset = 0;
	}
	debug(fprintf(stderr,"Setting altstrain_offset to be %d\n",altstrain_offset));
#endif
      }

      /* Handles case where true length is greater than provided
         coordinates.  This needs to be after call to IIT_get_exact */
      if (leftposition + truelength - 1U > rightposition) {
	debug(fprintf(stderr,"Extending endposition for truelength of %u\n",truelength));
	rightposition = leftposition + truelength;
	if (revcompp == true) {
	  endposition = startposition - truelength;
	} else {
	  endposition = startposition + truelength;
	}
      }

      /* In both cases, set file pointer for reference strain,
         although we won't write sequence of alternate strain.  For an
         alternate strain, ensure that we fill the reference strain
         with sufficient X's. */
      if (startposition > maxposition) {
	/* Start beyond end of file */
	debug(fprintf(stderr,"Filling with zeroes from %llu to %llu-1\n",
		      (unsigned long long) maxposition,(unsigned long long) startposition));
	fill_zero(refgenome_fp,maxposition,startposition,/*uncompressedp*/true,index1part);
	  
	if (contigtype > 0) {
#ifdef ALTSTRAIN
	  fill_zero(refgenome_fp,leftposition,rightposition + 1,/*uncompressedp*/true,index1part);
	  maxposition = currposition = rightposition + 1;
#endif
	} else {
	  maxposition = rightposition;
	  currposition = startposition;
	}

      } else {
	/* Start within file */
	if (contigtype > 0) {
#ifdef ALTSTRAIN
	  if (rightposition + 1 > maxposition) {
	    debug(fprintf(stderr,"Filling with zeroes from %u to %u-1\n",maxposition,rightposition+1));
	    fill_zero(refgenome_fp,maxposition,rightposition + 1,/*uncompressedp*/true,index1part);
	    maxposition = currposition = rightposition + 1;
	  }
#endif
	} else {
	  debug(fprintf(stderr,"Moving to %llu\n",(unsigned long long) startposition));
	  move_absolute(refgenome_fp,startposition);
	  currposition = startposition;
	}
      }

      /* SEQUENCE */
      fprintf(stderr,"Processing %u characters\n",truelength);
      while (truelength > 0) {
	if (truelength > BUFFERSIZE) {
	  if ((strlength = fread(Buffer,sizeof(char),BUFFERSIZE,input)) < BUFFERSIZE) {
	    fprintf(stderr,"Expecting %u more characters, but saw end of file\n",truelength);
	    exit(9);
	  }
	  truelength -= BUFFERSIZE;
	} else {
	  if ((strlength = fread(Buffer,sizeof(char),truelength,input)) < truelength) {
	    fprintf(stderr,"Expecting %u more characters, but saw end of file\n",truelength);
	    exit(9);
	  }
	  truelength = 0;
	}

	if (revcompp == true) {
	  make_reverse_buffered(Reverse,Buffer,strlength);
	  segment = Reverse;
	} else {
	  segment = Buffer;
	}

	if (contigtype > 0) {
#ifdef ALTSTRAING
	  /* Write alternate strain.  It is likely that IIT commands
	     will fail because they depend on \0 to terminate the segment.  */
	  if (revcompp == true) {
	    altstrain_offset -= strlength;
	    debug(fprintf(stderr,"Writing alternate strain at %d\n",altstrain_offset));
	    IIT_backfill_sequence(altstrain_iit,altstrain_index,altstrain_offset,segment);
	  } else {
	    debug(fprintf(stderr,"Writing alternate strain at %d\n",altstrain_offset));
	    IIT_backfill_sequence(altstrain_iit,altstrain_index,altstrain_offset,segment);
	    altstrain_offset += strlength;
	  }
#endif
	} else {
	  /* Write reference strain */
	  if (revcompp == true) {
	    debug(fprintf(stderr,"Filling with sequence from %llu-1 to %llu\n",
			  (unsigned long long) currposition,(unsigned long long) (currposition-strlength)));
	    currposition -= strlength;
	    fwrite(segment,sizeof(char),strlength,refgenome_fp);
	  } else {
	    debug(fprintf(stderr,"Filling with sequence from %llu to %llu-1\n",
			  (unsigned long long) currposition,(unsigned long long) (currposition+strlength)));
	    fwrite(segment,sizeof(char),strlength,refgenome_fp);
	    currposition += strlength;
	    if (currposition > maxposition) {
	      maxposition = currposition;
	    }
	  }
	}
      }

      if ((c = fgetc(input)) != EOF && c != '\n') {
	fprintf(stderr,"Expecting linefeed at end of sequence.  Saw %d instead\n",
		c);
	exit(9);
      }

    }
  }

  return;
}

static void
fill_circular_chromosomes (UINT4 *genomecomp, Univ_IIT_T chromosome_iit, int circular_typeint) {
  int indx;
  Univinterval_T interval;
  char *segment, *chr;
  Univcoord_T alias_startpos, alias_endpos, orig_startpos, orig_endpos;
  Chrpos_T seglength;
  bool allocp;

  for (indx = 1; indx <= Univ_IIT_total_nintervals(chromosome_iit); indx++) {
    interval = Univ_IIT_interval(chromosome_iit,indx);
    if (Univinterval_type(interval) == circular_typeint) {
      orig_startpos = Univinterval_low(interval);
      orig_endpos = Univinterval_high(interval);
      seglength = orig_endpos - orig_startpos + 1U;

      alias_startpos = orig_endpos + 1U;
      alias_endpos = alias_startpos + seglength - 1U;

      chr = Univ_IIT_label(chromosome_iit,indx,&allocp);
      /* Add 1U to report 1-based coordinates */
      fprintf(stderr,"Chromosome %s is circular.  Copying %llu..%llu to %llu..%llu\n",
	      chr,(unsigned long long) orig_startpos+1U,(unsigned long long) orig_endpos+1U,
	      (unsigned long long) alias_startpos+1U,(unsigned long long) alias_endpos+1U);
      if (allocp) {
	FREE(chr);
      }

      segment = (char *) CALLOC(seglength+1U,sizeof(char));
      /* Add 1U because procedures below are expecting exclusive coordinates */
      Genome_uncompress_memory(segment,genomecomp,orig_startpos,orig_endpos+1U); /* not Genome_uncompress_mmap, which does bigendian conversion */
      Compress_update_memory(/*nbadchars*/0,genomecomp,segment,alias_startpos,alias_endpos+1U);
      FREE(segment);
    }
  }

  return;
}



/* Puts reference genome into refgenome_fp (assume compressed),
   and puts alternate strain sequences into altstrain_iit. */
static void
genome_write_memory (FILE *refgenome_fp, FILE *input, 
		     Univ_IIT_T contig_iit, IIT_T altstrain_iit, 
		     Univ_IIT_T chromosome_iit, int circular_typeint, Genomecomp_T *genomecomp,
		     size_t nuint4, char *fileroot, int nmessages) {
  char Buffer[BUFFERSIZE], Complement[BUFFERSIZE], *segment;
  char *accession, *p;
  Univcoord_T leftposition, rightposition, startposition, endposition,
    maxposition = 0, currposition = 0;
  Chrpos_T truelength;
  int contigtype;
  bool revcompp = false;
#ifdef ALTSTRAIN
  int altstrain_index, altstrain_offset;
#endif
  int nbadchars = 0;
  int ncontigs = 0;
  size_t i;
  
  while (fgets(Buffer,BUFFERSIZE,input) != NULL) {
    if (Buffer[0] == '>') {
      /* HEADER */
      accession = parse_accession(Buffer);
      find_positions(&revcompp,&leftposition,&rightposition,&startposition,&endposition,
		     &truelength,&contigtype,accession,contig_iit);
      if (++ncontigs < nmessages) {
	if (revcompp == true) {
	  fprintf(stderr,"Writing contig %s to universal coordinates %llu..%llu\n",
		  accession,(unsigned long long) startposition,(unsigned long long) endposition);
	} else {
	  fprintf(stderr,"Writing contig %s to universal coordinates %llu..%llu\n",
		  accession,(unsigned long long) startposition+1U,(unsigned long long) endposition+1U);
	}
      } else if (ncontigs == nmessages) {
	fprintf(stderr,"More than %d contigs.  Will stop printing messages\n",nmessages);
      }

      if (contigtype > 0) {
#ifdef ALTSTRAIN
	fprintf(stderr," (alternate strain %s)",IIT_typestring(altstrain_iit,contigtype));
#endif
      }
      FREE(accession);

      if (contigtype > 0) {
#ifdef ALTSTRAIN
	/* Initialize file pointer for alternate strain */
	altstrain_index = IIT_get_exact(altstrain_iit,/*divstring*/NULL,leftposition,rightposition,contigtype);
	if (revcompp == true) {
	  altstrain_offset = rightposition + 1U - leftposition;
	} else {
	  altstrain_offset = 0;
	}
	debug(fprintf(stderr,"Setting altstrain_offset to be %d\n",altstrain_offset));
#endif
      }

      /* Handles case where true length is greater than provided
         coordinates.  This needs to be after call to IIT_get_exact */
      if (leftposition + truelength - 1U > rightposition) {
	debug(fprintf(stderr,"Extending endposition for truelength of %u\n",truelength));
	rightposition = leftposition + truelength;
	if (revcompp == true) {
	  endposition = startposition - truelength;
	} else {
	  endposition = startposition + truelength;
	}
      }

      /* In both cases, set file pointer for reference strain,
         although we won't write sequence of alternate strain.  For an
         alternate strain, ensure that we fill the reference strain
         with sufficient X's. */
      if (startposition > maxposition) {
	/* Start beyond end of file */
	debug(fprintf(stderr,"Filling with X's from %llu to %llu-1\n",
		      (unsigned long long) maxposition,(unsigned long long) startposition));
	fill_x_memory(genomecomp,maxposition,startposition);
	  
	if (contigtype > 0) {
#ifdef ALTSTRAIN
	  fill_x_memory(genomecomp,leftposition,rightposition + 1);
	  maxposition = currposition = rightposition + 1;
#endif
	} else {
	  maxposition = rightposition;
	  currposition = startposition;
	}

      } else {
	/* Start within file */
	if (contigtype > 0) {
#ifdef ALTSTRAIN
	  if (rightposition + 1 > maxposition) {
	    debug(fprintf(stderr,"Filling with X's from %u to %u-1\n",maxposition,rightposition+1));
	    fill_x_memory(genomecomp,maxposition,rightposition + 1);
	    maxposition = currposition = rightposition + 1;
	  }
#endif
	} else {
	  debug(fprintf(stderr,"Moving to %llu\n",(unsigned long long) startposition));
	  currposition = startposition;
	}
      }

    } else {
      /* SEQUENCE */
      if ((p = rindex(Buffer,'\n')) != NULL) {
	*p = '\0';
      }
      if ((p = rindex(Buffer,CONTROLM)) != NULL) {
	*p = '\0';
      }
      if (revcompp == true) {
	make_complement_buffered(Complement,Buffer,strlen(Buffer));
	segment = Complement;
      } else {
	segment = Buffer;
      }

      if (contigtype > 0) {
#ifdef ALTSTRAIN
	/* Write alternate strain */
	if (revcompp == true) {
	  altstrain_offset -= strlen(segment);
	  debug(fprintf(stderr,"Writing alternate strain at %d\n",altstrain_offset));
	  IIT_backfill_sequence(altstrain_iit,altstrain_index,altstrain_offset,segment);
	} else {
	  debug(fprintf(stderr,"Writing alternate strain at %d\n",altstrain_offset));
	  IIT_backfill_sequence(altstrain_iit,altstrain_index,altstrain_offset,segment);
	  altstrain_offset += strlen(segment);
	}
#endif
      } else {
	/* Write reference strain */
	if (revcompp == true) {
	  debug(fprintf(stderr,"Filling with sequence from %llu-1 to %llu\n",
			(unsigned long long) currposition,(unsigned long long) currposition-strlen(segment)));
	  currposition -= strlen(segment);
	  nbadchars = Compress_update_memory(nbadchars,genomecomp,segment,currposition,currposition+strlen(segment));
	} else {
	  debug(fprintf(stderr,"Filling with sequence from %llu to %llu-1\n",
			(unsigned long long) currposition,(unsigned long long) currposition+strlen(segment)));
	  nbadchars = Compress_update_memory(nbadchars,genomecomp,segment,currposition,currposition+strlen(segment));
	  currposition += strlen(segment);
	  if (currposition > maxposition) {
	    maxposition = currposition;
	  }
	}

      }
    }
  }

  fill_circular_chromosomes(genomecomp,chromosome_iit,circular_typeint);

  move_absolute(refgenome_fp,0U);
  FWRITE_UINTS(genomecomp,nuint4,refgenome_fp);

  fprintf(stderr,"A total of %d non-ACGTNX characters were seen in the genome.\n",nbadchars);

  return;
}

void
Genome_write_comp32 (char *genomesubdir, char *fileroot, FILE *input, 
		     Univ_IIT_T contig_iit, IIT_T altstrain_iit, Univ_IIT_T chromosome_iit,
		     bool uncompressedp, bool rawp, bool writefilep,
		     Univcoord_T genomelength, int index1part, int nmessages) {
  size_t nuint4;
  FILE *refgenome_fp;
  char *filename;
  Genomecomp_T *genomecomp;
  int circular_typeint;

  fprintf(stderr,"Genome length is %llu nt\n",(unsigned long long) genomelength);
  if (uncompressedp == true) {
    filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			       strlen(fileroot)+strlen(".genome")+1,sizeof(char));
    sprintf(filename,"%s/%s.genome",genomesubdir,fileroot);
    if ((refgenome_fp = FOPEN_WRITE_BINARY(filename)) == NULL) {
      fprintf(stderr,"Can't write to file %s\n",filename);
      exit(9);
    }
    if (rawp == true) {
      genome_writeraw_file(refgenome_fp,input,contig_iit,altstrain_iit,fileroot,index1part,nmessages);
    } else {
      genome_write_file(refgenome_fp,input,contig_iit,altstrain_iit,fileroot,
			/*uncompressedp*/true,index1part,nmessages);
    }
    fclose(refgenome_fp);
    FREE(filename);

  } else if (writefilep == true) {
    filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			       strlen(fileroot)+strlen(".genomecomp")+1,sizeof(char));
    sprintf(filename,"%s/%s.genomecomp",genomesubdir,fileroot);
    fprintf(stderr,"User requested build of genome in file\n");

    if ((refgenome_fp = FOPEN_RW_BINARY(filename)) == NULL) {
      fprintf(stderr,"Can't open file %s for read/write\n",filename);
      exit(9);
    }
    genome_write_file(refgenome_fp,input,contig_iit,altstrain_iit,fileroot,
		      /*uncompressedp*/false,index1part,nmessages);
    fclose(refgenome_fp);
    FREE(filename);

  } else {
    filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			       strlen(fileroot)+strlen(".genomecomp")+1,sizeof(char));
    sprintf(filename,"%s/%s.genomecomp",genomesubdir,fileroot);

    nuint4 = ((genomelength + 31)/32U)*3;
    fprintf(stderr,"Trying to allocate %llu*%d bytes of memory...",
	    (unsigned long long) nuint4,(int) sizeof(Genomecomp_T));
    genomecomp = (Genomecomp_T *) CALLOC_NO_EXCEPTION(nuint4,sizeof(Genomecomp_T));
    if (genomecomp == NULL) {
      fprintf(stderr,"failed.  Building genome in file.\n");
      if ((refgenome_fp = FOPEN_RW_BINARY(filename)) == NULL) {
	fprintf(stderr,"Can't open file %s for read/write\n",filename);
	exit(9);
      }
      genome_write_file(refgenome_fp,input,contig_iit,altstrain_iit,fileroot,
			/*uncompressedp*/false,index1part,nmessages);
      fclose(refgenome_fp);

    } else {
      fprintf(stderr,"succeeded.  Building genome in memory.\n");
      /* Creates X's at end */
      genomecomp[nuint4-3] = 0xFFFFFFFF;
      genomecomp[nuint4-2] = 0xFFFFFFFF;
      genomecomp[nuint4-1] = 0xFFFFFFFF;

      if ((refgenome_fp = FOPEN_WRITE_BINARY(filename)) == NULL) {
	fprintf(stderr,"Can't open file %s for write\n",filename);
	exit(9);
      }
      circular_typeint = Univ_IIT_typeint(chromosome_iit,"circular");
      genome_write_memory(refgenome_fp,input,contig_iit,altstrain_iit,chromosome_iit,circular_typeint,
			  genomecomp,nuint4,fileroot,nmessages);
      fclose(refgenome_fp);
      FREE(genomecomp);
    }

    FREE(filename);
  }

  return;
}


