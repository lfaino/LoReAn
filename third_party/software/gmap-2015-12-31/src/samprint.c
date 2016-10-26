static char rcsid[] = "$Id: samprint.c 183725 2016-02-04 00:40:15Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "samprint.h"
#include "samflags.h"
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "mem.h"
#include "complement.h"
#include "mapq.h"
#include "assert.h"


#define SANGER_ILLUMINA_DIFF 31
/* #define PRINT_AMBIG_COORDS 1 */

/* BAM appears to truncate the H information on the ends of a cigar */
/* Also, this provides the information needed for getting term information */


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


/* compute_cigar */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif


/* print_md_string */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif


/* overlap */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif

/* compute_chrpos */
#ifdef DEBUG4
#define debug4(x) x
#else
#define debug4(x)
#endif


static bool add_paired_nomappers_p;
static bool paired_flag_means_concordant_p;
static bool quiet_if_excessive_p;
static int maxpaths_report;
static char *failedinput_root;
static bool fastq_format_p;
static bool hide_soft_clips_p;

static bool clip_overlap_p;
static bool merge_overlap_p;

static bool sam_multiple_primaries_p;
static bool force_xs_direction_p;
static bool md_lowercase_variant_p;
static IIT_T snps_iit;

static Univ_IIT_T chromosome_iit;
static Genome_T genome;

void
SAM_setup (bool add_paired_nomappers_p_in, bool paired_flag_means_concordant_p_in,
	   bool quiet_if_excessive_p_in, int maxpaths_report_in,
	   char *failedinput_root_in, bool fastq_format_p_in, bool hide_soft_clips_p_in,
	   bool clip_overlap_p_in, bool merge_overlap_p_in, bool sam_multiple_primaries_p_in,
	   bool force_xs_direction_p_in, bool md_lowercase_variant_p_in, IIT_T snps_iit_in,
	   Univ_IIT_T chromosome_iit_in, Genome_T genome_in) {
  add_paired_nomappers_p = add_paired_nomappers_p_in;
  paired_flag_means_concordant_p = paired_flag_means_concordant_p_in;
  quiet_if_excessive_p = quiet_if_excessive_p_in;
  failedinput_root = failedinput_root_in;
  fastq_format_p = fastq_format_p_in;
  hide_soft_clips_p = hide_soft_clips_p_in;
  clip_overlap_p = clip_overlap_p_in;
  merge_overlap_p = merge_overlap_p_in;
  maxpaths_report = maxpaths_report_in;
  sam_multiple_primaries_p = sam_multiple_primaries_p_in;
  force_xs_direction_p = force_xs_direction_p_in;
  md_lowercase_variant_p = md_lowercase_variant_p_in;
  snps_iit = snps_iit_in;

  chromosome_iit = chromosome_iit_in;
  genome = genome_in;

  return;
}


unsigned int
SAM_compute_flag (bool plusp, Stage3end_T mate, Resulttype_T resulttype,
		  bool first_read_p, int pathnum, int npaths, bool artificial_mate_p, int npaths_mate,
		  int absmq_score, int first_absmq, bool invertp, bool invert_mate_p) {
  unsigned int flag = 0U;

  debug(printf("Resulttype: %s\n",Resulttype_string(resulttype)));

  if (npaths == 0) {
    debug(printf("npaths = 0, so QUERY_UNMAPPED %d\n",QUERY_UNMAPPED));
    flag |= QUERY_UNMAPPED;
  } else if (plusp == invertp) {
    debug(printf("plusp %d and invertp %d, so QUERY_MINUSP %d\n",
		 plusp,invertp,QUERY_MINUSP));
    flag |= QUERY_MINUSP;
  }

  if (resulttype == SINGLEEND_NOMAPPING || resulttype == SINGLEEND_UNIQ || resulttype == SINGLEEND_TRANSLOC || resulttype == SINGLEEND_MULT) {
    /* No first or second read or mate */
  } else {
    debug(printf("PAIRED_READ %d\n",PAIRED_READ));
    flag |= PAIRED_READ;
    if (first_read_p == true) {
      debug(printf("FIRST_READ %d\n",FIRST_READ_P));
      flag |= FIRST_READ_P;
    } else {
      debug(printf("SECOND_READ %d\n",SECOND_READ_P));
      flag |= SECOND_READ_P;
    }
    if (artificial_mate_p == true) {
      debug(printf("MATE_UNMAPPED %d\n",MATE_UNMAPPED));
      flag |= MATE_UNMAPPED;

    } else if (npaths_mate == 0) {
      debug(printf("MATE_UNMAPPED %d\n",MATE_UNMAPPED));
      flag |= MATE_UNMAPPED;

    } else if (quiet_if_excessive_p && npaths_mate > maxpaths_report) {
      debug(printf("MATE_UNMAPPED %d\n",MATE_UNMAPPED));
      flag |= MATE_UNMAPPED;

    } else if (mate == NULL) {
      /* Unpaired; no mate.  Not clear if should be MATE_UNMAPPED. */

    } else if (npaths == 0) {
      /* Need to check npaths == 0 in case clipping of overlaps results in a nomapping */
      if (Stage3end_plusp(mate) == invert_mate_p) {
	debug(printf("MATE_MINUSP %d\n",MATE_MINUSP));
	flag |= MATE_MINUSP;
      }

    } else if (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_TRANSLOC || resulttype == CONCORDANT_MULT) {
      /* Can distinguish concordant mappings by presence of insert length */
      debug(printf("PAIRED_MAPPING %d\n",PAIRED_MAPPING));
      flag |= PAIRED_MAPPING;
      if (0 && Stage3end_chrnum(mate) == 0) {
	/* Splice without a direction.  But want the effective plusp anyway. */

      } else if (Stage3end_plusp(mate) == invert_mate_p) {
	debug(printf("MATE_MINUSP %d\n",MATE_MINUSP));
	flag |= MATE_MINUSP;
      }

    } else if (resulttype == PAIRED_UNIQ || resulttype == PAIRED_MULT) {
      /* Note: We are counting PAIRED_UNIQ and PAIRED_MULT as "paired" mappings.
	 However, we are no longer counting UNPAIRED_UNIQ as a "paired" mapping. */
      if (paired_flag_means_concordant_p == true) {
	/* Don't turn on paired flag */
      } else {
	debug(printf("PAIRED_MAPPING %d\n",PAIRED_MAPPING));
	flag |= PAIRED_MAPPING;
      }
      if (0 && Stage3end_chrnum(mate) == 0) {
	/* Splice without a direction.  But want the effective plusp anyway. */

      } else if (Stage3end_plusp(mate) == invert_mate_p) {
	debug(printf("MATE_MINUSP %d\n",MATE_MINUSP));
	flag |= MATE_MINUSP;
      }

    } else {
      if (0 && Stage3end_chrnum(mate) == 0) {
	/* Splice without a direction.  But want the effective plusp anyway. */

      } else if (Stage3end_plusp(mate) == invert_mate_p) {
	debug(printf("MATE_MINUSP %d\n",MATE_MINUSP));
	flag |= MATE_MINUSP;
      }
    }
  }

  if (pathnum > 1) {
    if (sam_multiple_primaries_p == false) {
      debug(printf("NOT_PRIMARY %d\n",NOT_PRIMARY));
      flag |= NOT_PRIMARY;
    } else if (absmq_score != first_absmq) {
      debug(printf("NOT_PRIMARY %d\n",NOT_PRIMARY));
      flag |= NOT_PRIMARY;
    } else {
      /* Just as good as first alignment, so don't mark as altloc */
    }
  }

  return flag;
}


Chrpos_T
SAM_compute_chrpos (int hardclip_low, int hardclip_high, Stage3end_T this, int querylength,
		    bool first_read_p) {
  Substring_T substring;
  Hittype_T hittype;

  if (this == NULL) {
    return 0U;

  } else if ((hittype = Stage3end_hittype(this)) == GMAP) {
    return Pair_genomicpos_low(hardclip_low,hardclip_high,Stage3end_pairarray(this),Stage3end_npairs(this),
			       querylength,/*watsonp*/Stage3end_plusp(this),hide_soft_clips_p);

  } else if (hittype == SAMECHR_SPLICE || hittype == TRANSLOC_SPLICE) {
    /* Want concordant substring */
    if (Stage3end_plusp(this) == true) {
      if (first_read_p == true) {
	/* Eventually want substringN */
	substring = Stage3end_substring2(this);
      } else {
	substring = Stage3end_substring1(this);
      }
    } else {
      if (first_read_p == true) {
	substring = Stage3end_substring1(this);
      } else {
	/* Eventually want substringN */
	substring = Stage3end_substring2(this);
      }
    }
    return Substring_compute_chrpos(substring,hardclip_low,hide_soft_clips_p);

  } else {
    /* Want low substring */
    substring = Stage3end_substring_low(this,hardclip_low);
    return Substring_compute_chrpos(substring,hardclip_low,hide_soft_clips_p);
  }
}

static void
print_chromosomal_pos (Filestring_T fp, Chrnum_T chrnum, Chrpos_T chrpos, Chrpos_T chrlength,
		       Univ_IIT_T chromosome_iit) {
  bool allocp;
  char *chr;

#if 0
  if (chrpos == 0U) {
    /* No mapping */
    FPRINTF(fp,"\t*\t0");
    return;
  }
#endif

  if (chrnum == 0) {
    /* Interchromosomal splice */
    fprintf(stderr,"Trying to print interchrosomal splice in one line\n");
    abort();

  } else {
    chr = Univ_IIT_label(chromosome_iit,chrnum,&allocp);

    /* chrpos already in 1-based coordinates */
    if (chrpos > chrlength) {
      FPRINTF(fp,"\t%s\t%u",chr,chrpos - chrlength /*+1U*/);
    } else {
      FPRINTF(fp,"\t%s\t%u",chr,chrpos /*+1U*/);
    }

    if (allocp == true) {
      FREE(chr);
    }

    return;
  }
}

static void
print_mate_chromosomal_pos (Filestring_T fp, Chrnum_T mate_chrnum, Chrnum_T mate_effective_chrnum,
			    Chrpos_T mate_chrpos, Chrpos_T mate_chrlength, Chrnum_T anchor_chrnum, Chrpos_T anchor_chrpos,
			    Univ_IIT_T chromosome_iit) {
  bool allocp;
  char *chr;

  if (mate_chrpos == 0U) {
    FPRINTF(fp,"\t*\t0");
    return;

  } else {
    if (mate_chrnum == 0) {
      /* Interchromosomal splice.  Choose effective chrnum. */
      mate_chrnum = mate_effective_chrnum;
    }
      
    if (anchor_chrpos > 0U && anchor_chrnum > 0 && mate_chrnum == anchor_chrnum) {
      /* chrpos already in 1-based coordinates */
      if (mate_chrpos > mate_chrlength) {
	FPRINTF(fp,"\t=\t%u",mate_chrpos - mate_chrlength /*+1U*/);
      } else {
	FPRINTF(fp,"\t=\t%u",mate_chrpos /*+1U*/);
      }

    } else {
      chr = Univ_IIT_label(chromosome_iit,mate_chrnum,&allocp);

      /* chrpos already in 1-based coordinates */
      if (mate_chrpos > mate_chrlength) {
	FPRINTF(fp,"\t%s\t%u",chr,mate_chrpos - mate_chrlength /*+1U*/);
      } else {
	FPRINTF(fp,"\t%s\t%u",chr,mate_chrpos /*+1U*/);
      }

      if (allocp == true) {
	FREE(chr);
      }
    }
    
    /* chrpos already in 1-based coordinates */
    return;
  }
}





static char complCode[128] = COMPLEMENT_LC;

static void
make_complement_buffered (char *complement, char *sequence, unsigned int length) {
  int i, j;

  /* complement = (char *) CALLOC(length+1,sizeof(char)); */
  for (i = length-1, j = 0; i >= 0; i--, j++) {
    complement[j] = complCode[(int) sequence[i]];
  }
  complement[length] = '\0';
  return;
}



/* npaths could be non-zero, if user selected --quiet-if-excessive */
void
SAM_print_nomapping (Filestring_T fp, char *abbrev, Shortread_T queryseq, Stage3end_T mate, char *acc1, char *acc2,
		     Univ_IIT_T chromosome_iit, Resulttype_T resulttype, bool first_read_p,
		     int pathnum, int npaths, bool artificial_mate_p, int npaths_mate,
		     Chrpos_T mate_chrpos, int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p) {
  unsigned int flag;


  /* 1. QNAME */
  if (acc2 == NULL) {
    FPRINTF(fp,"%s",acc1);
  } else {
    FPRINTF(fp,"%s,%s",acc1,acc2);
  }
  
  /* 2. FLAG */
  /* 3. RNAME: chr */
  /* 4. POS: chrpos */
  /* 5. MAPQ: Mapping quality.  Picard says MAPQ should be 0 for an unmapped read */
  /* 6. CIGAR */
  flag = SAM_compute_flag(/*plusp (NA)*/true,mate,resulttype,first_read_p,
			  /*pathnum*/0,/*npaths*/0,artificial_mate_p,npaths_mate,
			  /*absmq_score*/0,/*first_absmq*/0,invertp,invert_mate_p);
  FPRINTF(fp,"\t%u\t*\t0\t0\t*",flag);

  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  print_mate_chromosomal_pos(fp,Stage3end_chrnum(mate),Stage3end_effective_chrnum(mate),
			     mate_chrpos,Stage3end_chrlength(mate),
			     /*anchor_chrnum*/0,/*anchor_chrpos*/0U,chromosome_iit);


  /* 9. ISIZE: Insert size */
  FPRINTF(fp,"\t0");

  /* 10. SEQ: queryseq and 11. QUAL: quality scores */
  /* Since there is no mapping, we print the original query sequence. */
  if (invertp == false) {
    Shortread_print_chopped_sam(fp,queryseq,/*hardclip_low*/0,/*hardclip_high*/0);
    FPRINTF(fp,"\t");
    Shortread_print_quality(fp,queryseq,/*hardclip_low*/0,/*hardclip_high*/0,
			    quality_shift,/*show_chopped_p*/false);
  } else {
    Shortread_print_chopped_revcomp_sam(fp,queryseq,/*hardclip_low*/0,/*hardclip_high*/0);
    FPRINTF(fp,"\t");
    Shortread_print_quality_revcomp(fp,queryseq,/*hardclip_low*/0,/*hardclip_high*/0,
				    quality_shift,/*show_chopped_p*/false);
  }

  /* 12. TAGS: RG */
  if (sam_read_group_id != NULL) {
    FPRINTF(fp,"\tRG:Z:%s",sam_read_group_id);
  }
  
  /* 12. TAGS: NH */
  if (npaths > 0) {
    FPRINTF(fp,"\tNH:i:%d",npaths);
    if (add_paired_nomappers_p == true) {
      FPRINTF(fp,"\tHI:i:%d",pathnum);
    }
  }

  /* 12. TAGS: XB */
  Shortread_print_barcode(fp,queryseq);

  /* 12. TAGS: XP */
  Shortread_print_chop(fp,queryseq,invertp);

  /* 12. TAGS: XO */
  FPRINTF(fp,"\tXO:Z:%s",abbrev);

  FPRINTF(fp,"\n");

  return;
}


/* Derived from print_tokens_gff3 */
static void
print_tokens_sam (Filestring_T fp, List_T tokens) {
  List_T p;
  char *token;
  
  for (p = tokens; p != NULL; p = List_next(p)) {
    token = (char *) List_head(p);
    FPRINTF(fp,"%s",token);
    FREE(token);
  }

  return;
}

static List_T
push_token (List_T tokens, char *token) {
  char *copy;

  copy = (char *) CALLOC(strlen(token)+1,sizeof(char));
  strcpy(copy,token);
  return List_push(tokens,(void *) copy);
}


#if 0
/* Currently used for insertions and deletions */
static List_T
compute_cigar_old (List_T tokens, char type, int stringlength, int querypos, int querylength,
		   int hardclip_low, int hardclip_high, bool plusp, bool firstp, bool lastp) {
  char token[10];
  
  debug1(printf("\nEntering compute_cigar with type %c, stringlength %d, querypos %d, querylength %d, hardclip_low %d, hardclip_high %d, plusp %d\n",
		type,stringlength,querypos,querylength,hardclip_low,hardclip_high,plusp));

  if (firstp == true) {
    debug1(printf("firstp is true\n"));
    if (plusp == true) {
      if (hardclip_low > 0) {
	sprintf(token,"%dH",hardclip_low);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
      if (querypos > hardclip_low) {
	sprintf(token,"%dS",querypos - hardclip_low);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
    } else {
      if (hardclip_high > 0) {
	sprintf(token,"%dH",hardclip_high);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
      if (querypos < querylength - hardclip_high) {
	sprintf(token,"%dS",querypos - hardclip_high);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
    }
  }

  if (type == 'D' || type == 'N') {
    if (querypos < hardclip_low || querypos >= querylength - hardclip_high) {
      stringlength = 0;
    }

  } else if (plusp == true) {
    debug1(printf("Comparing querypos %d..%d against %d..%d\n",
		  querypos,querypos + stringlength,hardclip_low,querylength - hardclip_high));
    if (/* querypos < hardclip_low && */querypos + stringlength < hardclip_low) {
      /* Print nothing */
      stringlength = 0;
      debug1(printf("Case 1: stringlength 0\n"));
    } else if (querypos < hardclip_low) {
      if (querypos + stringlength < querylength - hardclip_high) {
	/* Print part after hardclip_low */
	stringlength = (querypos + stringlength) - hardclip_low;
	debug1(printf("Case 2: stringlength %d\n",stringlength));
      } else {
	/* Print part between hardclip_low and hardclip_high */
	stringlength = (querylength - hardclip_high) - hardclip_low;
	debug1(printf("Case 3: stringlength %d\n",stringlength));
      }
    } else if (querypos < querylength - hardclip_high) {
      if (querypos + stringlength >= querylength - hardclip_high) {
	/* Print up to hardclip_high */
	stringlength = (querylength - hardclip_high) - querypos;
	debug1(printf("Case 4: stringlength %d\n",stringlength));
      } else {
	/* Print full stringlength */
	debug1(printf("Case 5: stringlength %d\n",stringlength));
      }
    } else {
      /* Print nothing */
      stringlength = 0;
      debug1(printf("Case 6: stringlength 0\n"));
    }

  } else {
    debug1(printf("Comparing querypos %d..%d against %d..%d\n",
		  querypos,querypos - stringlength,hardclip_low,querylength - hardclip_high));
    if (/* querypos >= querylength - hardclip_high && */ querypos - stringlength >= querylength - hardclip_high) {
      /* Print nothing */
      stringlength = 0;
      debug1(printf("Case 1: stringlength 0\n"));
    } else if (querypos >= querylength - hardclip_high) {
      if (querypos - stringlength >= hardclip_low) {
	/* Print part after hardclip_high */
	stringlength = (querylength - hardclip_high) - (querypos - stringlength);
	debug1(printf("Case 2: stringlength %d\n",stringlength));
      } else {
	/* Print part between hardclip_low and hardclip_high */
	stringlength = (querylength - hardclip_high) - hardclip_low;
	debug1(printf("Case 3: stringlength %d\n",stringlength));
      }
    } else if (querypos >= hardclip_low) {
      if (querypos - stringlength < hardclip_low) {
	/* Print up to hardclip_low */
	stringlength = querypos - hardclip_low;
	debug1(printf("Case 4: stringlength %d\n",stringlength));
      } else {
	/* Print full stringlength */
	debug1(printf("Case 5: stringlength %d\n",stringlength));
      }
    } else {
      /* Print nothing */
      stringlength = 0;
      debug1(printf("Case 5: stringlength 0\n"));
    }
  }

  if (stringlength > 0) {
    sprintf(token,"%d%c",stringlength,type);
    debug1(printf("Pushing token %s\n",token));
    tokens = push_token(tokens,token);
  }

  if (lastp == true) {
    debug1(printf("lastp is true\n"));
    if (plusp == true) {
      querypos += stringlength;
      if (querypos < querylength - 1 - hardclip_high) {
	sprintf(token,"%dS",querylength - 1 - hardclip_high - querypos);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
      if (hardclip_high > 0) {
	sprintf(token,"%dH",hardclip_high);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
    } else {
      querypos -= stringlength;
      if (querypos > hardclip_low) {
	sprintf(token,"%dS",hardclip_low - querypos);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
      if (hardclip_low > 0) {
	sprintf(token,"%dH",hardclip_low);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
    }
  }

  return tokens;
}
#endif


/* Currently used for insertions and deletions */
static List_T
compute_cigar (List_T tokens, char type, int stringlength, int querypos, int querylength,
	       int hardclip_low, int hardclip_high, bool plusp, int lastp) {
  int matchlength = 0;
  int startpos, endpos;
  int cliplength = 0;
  char token[10];
  
  if (plusp == true) {
    debug1(printf("\nEntering compute_cigar with type %c, stringlength %d, querypos %d, querylength %d, hardclip_low %d, hardclip_high %d, plus\n",
		  type,stringlength,querypos,querylength,hardclip_low,hardclip_high));
    if (hardclip_low > querypos) { /* > not >= */
      startpos = hardclip_low;
      cliplength = hardclip_low;
    } else {
      startpos = querypos;
    }

    if (querylength - hardclip_high < querypos + stringlength) {
      endpos = querylength - hardclip_high;
      debug1(printf("  endpos %d = querylength %d - hardclip_high %d\n",endpos,querylength,hardclip_high));
    } else {
      endpos = querypos + stringlength;
      debug1(printf("  endpos %d = querypos %d + stringlength %d\n",endpos,querypos,stringlength));
    }

    debug1(printf("  new startpos %d, endpos %d, cliplength %d\n",startpos,endpos,cliplength));

    if (endpos >= startpos) {
      if (cliplength > 0) {
	debug1(printf("  Pushing initial %dH\n",cliplength));
	sprintf(token,"%dH",cliplength);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
      matchlength = endpos - startpos;
      if (matchlength > 0) {
	debug1(printf("  Pushing %d%c\n",matchlength,type));
	sprintf(token,"%d%c",matchlength,type);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
    }


    if (lastp == true) {
      /* cliplength = querypos + stringlength - endpos; */
      cliplength = querylength - endpos;
      if (cliplength > 0) {
	debug1(printf("  Pushing final %dH\n",cliplength));
	sprintf(token,"%dH",cliplength);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
    }

  } else {
    debug1(printf("\nEntering compute_cigar with type %c, stringlength %d, querypos %d, querylength %d, hardclip_low %d, hardclip_high %d, minus\n",
		  type,stringlength,querypos,querylength,hardclip_low,hardclip_high));

    if (querylength - hardclip_low < querypos) {
      startpos = querylength - hardclip_low;
      cliplength = hardclip_low;
    } else {
      startpos = querypos;
    }

    if (hardclip_high >= querypos - stringlength) {
      endpos = hardclip_high;
      debug1(printf("  endpos %d = hardclip_high %d\n",endpos,hardclip_high));
    } else {
      endpos = querypos - stringlength;
      debug1(printf("  endpos %d = querypos %d - stringlength %d\n",endpos,querypos,stringlength));
    }

    debug1(printf("  new startpos %d, endpos %d, cliplength %d\n",startpos,endpos,cliplength));

    if (endpos <= startpos) {
      if (cliplength > 0) {
	debug1(printf("  Pushing initial %dH\n",cliplength));
	sprintf(token,"%dH",cliplength);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
      matchlength = startpos - endpos;
      if (matchlength > 0) {
	debug1(printf("  Pushing %d%c\n",matchlength,type));
	sprintf(token,"%d%c",matchlength,type);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
    }


    if (lastp == true) {
      cliplength = endpos;
      if (cliplength > 0) {
	debug1(printf("  Pushing final %dH\n",cliplength));
	sprintf(token,"%dH",cliplength);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
    }
  }

  return tokens;
}


/* Modified from compute_cigar */
static Intlist_T
compute_cigar_types_only (Intlist_T types, char type, int stringlength, int querypos, int querylength,
			  int hardclip_low, int hardclip_high, bool plusp, int lastp) {
  int matchlength = 0;
  int startpos, endpos;
  int cliplength = 0;
  
  if (plusp == true) {
    debug1(printf("\nEntering compute_cigar_types_only with type %c, stringlength %d, querypos %d, querylength %d, hardclip_low %d, hardclip_high %d, plus\n",
		  type,stringlength,querypos,querylength,hardclip_low,hardclip_high));
    if (hardclip_low > querypos) { /* > not >= */
      startpos = hardclip_low;
      cliplength = hardclip_low;
    } else {
      startpos = querypos;
    }

    if (querylength - hardclip_high < querypos + stringlength) {
      endpos = querylength - hardclip_high;
      debug1(printf("  endpos %d = querylength %d - hardclip_high %d\n",endpos,querylength,hardclip_high));
    } else {
      endpos = querypos + stringlength;
      debug1(printf("  endpos %d = querypos %d + stringlength %d\n",endpos,querypos,stringlength));
    }

    debug1(printf("  new startpos %d, endpos %d, cliplength %d\n",startpos,endpos,cliplength));

    if (endpos >= startpos) {
      if (cliplength > 0) {
	debug1(printf("  Pushing initial %dH\n",cliplength));
	types = Intlist_push(types,'H');
      }
      matchlength = endpos - startpos;
      if (matchlength > 0) {
	debug1(printf("  Pushing %d%c\n",matchlength,type));
	types = Intlist_push(types,type);
      }
    }


    if (lastp == true) {
      /* cliplength = querypos + stringlength - endpos; */
      cliplength = querylength - endpos;
      if (cliplength > 0) {
	debug1(printf("  Pushing final %dH\n",cliplength));
	types = Intlist_push(types,'H');
      }
    }

  } else {
    debug1(printf("\nEntering compute_cigar with type %c, stringlength %d, querypos %d, querylength %d, hardclip_low %d, hardclip_high %d, minus\n",
		  type,stringlength,querypos,querylength,hardclip_low,hardclip_high));

    if (querylength - hardclip_low < querypos) {
      startpos = querylength - hardclip_low;
      cliplength = hardclip_low;
    } else {
      startpos = querypos;
    }

    if (hardclip_high >= querypos - stringlength) {
      endpos = hardclip_high;
      debug1(printf("  endpos %d = hardclip_high %d\n",endpos,hardclip_high));
    } else {
      endpos = querypos - stringlength;
      debug1(printf("  endpos %d = querypos %d - stringlength %d\n",endpos,querypos,stringlength));
    }

    debug1(printf("  new startpos %d, endpos %d, cliplength %d\n",startpos,endpos,cliplength));

    if (endpos <= startpos) {
      if (cliplength > 0) {
	debug1(printf("  Pushing initial %dH\n",cliplength));
	types = Intlist_push(types,'H');
      }
      matchlength = startpos - endpos;
      if (matchlength > 0) {
	debug1(printf("  Pushing %d%c\n",matchlength,type));
	types = Intlist_push(types,type);
      }
    }


    if (lastp == true) {
      cliplength = endpos;
      if (cliplength > 0) {
	debug1(printf("  Pushing final %dH\n",cliplength));
	types = Intlist_push(types,'H');
      }
    }
  }

  return types;
}


static void
print_cigar (Filestring_T fp, char type, int stringlength, int querypos, int querylength,
	     int hardclip_low, int hardclip_high, bool plusp, bool lastp, int trimlength) {
  int matchlength = 0;
  int startpos, endpos;
  int cliplength = 0;
  
  if (plusp == true) {
    debug1(printf("\nEntering print_cigar with type %c, stringlength %d, querypos %d, querylength %d, hardclip_low %d, hardclip_high %d, plus\n",
		  type,stringlength,querypos,querylength,hardclip_low,hardclip_high));
    if (hardclip_low > querypos) { /* > not >= */
      startpos = hardclip_low;
      cliplength = hardclip_low;
    } else {
      startpos = querypos;
    }

    if (querylength - hardclip_high < querypos + stringlength) {
      endpos = querylength - hardclip_high;
      debug1(printf("  endpos %d = querylength %d - hardclip_high %d\n",endpos,querylength,hardclip_high));
    } else {
      endpos = querypos + stringlength;
      debug1(printf("  endpos %d = querypos %d + stringlength %d\n",endpos,querypos,stringlength));
    }

    debug1(printf("  new startpos %d, endpos %d, cliplength %d\n",startpos,endpos,cliplength));

    if (endpos >= startpos) {
      if (cliplength > 0) {
	debug1(printf("  Pushing initial %dH\n",cliplength));
	FPRINTF(fp,"%dH",cliplength);
      }
      matchlength = endpos - startpos;
      if (matchlength <= 0) {
	/* Skip */
      } else if (type != 'E') {
	debug1(printf("  Pushing %d%c\n",matchlength,type));
	FPRINTF(fp,"%d%c",matchlength,type);
      } else if (matchlength == trimlength) {
	debug1(printf("  Pushing %dS\n",matchlength));
	FPRINTF(fp,"%dS",matchlength);
      } else {
	debug1(printf("  Pushing %dH\n",matchlength));
	FPRINTF(fp,"%dH",matchlength);
      }
    }


    if (lastp == true) {
      /* cliplength = querypos + stringlength - endpos; */
      cliplength = querylength - endpos;
      if (cliplength > 0) {
	debug1(printf("  Pushing final %dH\n",cliplength));
	FPRINTF(fp,"%dH",cliplength);
      }
    }

  } else {
    debug1(printf("\nEntering print_cigar with type %c, stringlength %d, querypos %d, querylength %d, hardclip_low %d, hardclip_high %d, minus\n",
		  type,stringlength,querypos,querylength,hardclip_low,hardclip_high));

    if (querylength - hardclip_low < querypos) {
      startpos = querylength - hardclip_low;
      cliplength = hardclip_low;
    } else {
      startpos = querypos;
    }

    if (hardclip_high >= querypos - stringlength) {
      endpos = hardclip_high;
      debug1(printf("  endpos %d = hardclip_high %d\n",endpos,hardclip_high));
    } else {
      endpos = querypos - stringlength;
      debug1(printf("  endpos %d = querypos %d - stringlength %d\n",endpos,querypos,stringlength));
    }

    debug1(printf("  new startpos %d, endpos %d, cliplength %d\n",startpos,endpos,cliplength));

    if (endpos <= startpos) {
      if (cliplength > 0) {
	debug1(printf("  Pushing initial %dH\n",cliplength));
	FPRINTF(fp,"%dH",cliplength);
      }
      matchlength = startpos - endpos;
      if (matchlength <= 0) {
	/* Skip */
      } else if (type != 'E') {
	debug1(printf("  Pushing %d%c\n",matchlength,type));
	FPRINTF(fp,"%d%c",matchlength,type);
      } else if (matchlength == trimlength) {
	debug1(printf("  Pushing %dS\n",matchlength));
	FPRINTF(fp,"%dS",matchlength);
      } else {
	debug1(printf("  Pushing %dH\n",matchlength));
	FPRINTF(fp,"%dH",matchlength);
      }
    }


    if (lastp == true) {
      cliplength = endpos;
      if (cliplength > 0) {
	debug1(printf("  Pushing final %dH\n",cliplength));
	FPRINTF(fp,"%dH",cliplength);
      }
    }
  }

  return;
}


static int
print_md_string (bool *printp, int *nmismatches_refdiff, int *nmismatches_bothdiff,
		 Filestring_T fp, int matchlength, char *genomicfwd_refdiff, char *genomicfwd_bothdiff,
		 int stringlength, int querypos, int querylength,
		 int hardclip_low, int hardclip_high, bool plusp, bool lastp) {
  int starti, endi, i;
  int local_nmismatches = 0;
  bool hardclip_end_p = false;

  if (plusp == true) {
    debug2(printf("\nEntering md_string with matchlength %d, querypos %d, querylength %d, hardclip_low %d, hardclip_high %d, plus: %s ref, %s both\n",
		  matchlength,querypos,querylength,hardclip_low,hardclip_high,genomicfwd_refdiff,genomicfwd_bothdiff));
    if (hardclip_low == 0) {
      starti = 0;
      hardclip_end_p = true;
    } else if (hardclip_low > querypos) {
      /* startpos = hardclip_low; */
      starti = hardclip_low - querypos;
      hardclip_end_p = true;
      debug2(printf("  Setting starti %d = hardclip_low %d - querypos %d\n",
		    starti,hardclip_low,querypos));
    } else {
      /* startpos = querypos; */
      starti = 0;
    }

    if (querylength - hardclip_high < querypos + stringlength) {
      /* endpos = querylength - hardclip_high; */
      endi = (querylength - hardclip_high) - querypos;
      debug2(printf("  Setting endi %d = (querylength %d - hardclip_high %d) - querypos %d\n",
		    endi,querylength,hardclip_high,querypos));
    } else {
      /* endpos = querypos + stringlength; */
      endi = stringlength;
    }

    debug2(printf("  Counting matches from %d to %d\n",starti,endi));

    if (genomicfwd_refdiff == NULL) {
      if (endi > starti) {
	matchlength += (endi - starti);
      }

    } else if (md_lowercase_variant_p == false) {
      for (i = starti; i < endi; i++) {
	if (isupper(genomicfwd_refdiff[i])) {
	  matchlength++;

	} else {
	  /* A true mismatch against both variants */
	  if (matchlength > 0 || hardclip_end_p == true) {
	    FPRINTF(fp,"%d",matchlength);
	    *printp = true;
	    hardclip_end_p = false;
	  }
	  FPRINTF(fp,"%c",toupper(genomicfwd_refdiff[i]));
	  *printp = true;
	  local_nmismatches += 1;
	  matchlength = 0;
	}
      }
      *nmismatches_refdiff += local_nmismatches;

    } else {
      for (i = starti; i < endi; i++) {
	if (isupper(genomicfwd_refdiff[i])) {
	  matchlength++;

	} else if (isupper(genomicfwd_bothdiff[i])) {
	  /* A mismatch against the reference only => alternate variant */
	  if (matchlength > 0 || hardclip_end_p == true) {
	    FPRINTF(fp,"%d",matchlength);
	    *printp = true;
	    hardclip_end_p = false;
	  }
	  FPRINTF(fp,"%c",genomicfwd_refdiff[i]); /* Leave as lower case */
	  *printp = true;
	  local_nmismatches += 1;
	  matchlength = 0;

	} else {
	  /* A true mismatch against both variants */
	  if (matchlength > 0 || hardclip_end_p == true) {
	    FPRINTF(fp,"%d",matchlength);
	    *printp = true;
	    hardclip_end_p = false;
	  }
	  FPRINTF(fp,"%c",toupper(genomicfwd_refdiff[i]));
	  *printp = true;
	  local_nmismatches += 1;
	  matchlength = 0;
	}
      }
      *nmismatches_refdiff += local_nmismatches;
    }

  } else {
    debug2(printf("\nEntering md_string with matchlength %d, querypos %d, querylength %d, hardclip_low %d, hardclip_high %d, minus: %s ref, %s both\n",
		  matchlength,querypos,querylength,hardclip_low,hardclip_high,genomicfwd_refdiff,genomicfwd_bothdiff));
    querypos = querylength - querypos - stringlength;
    debug2(printf("  Revising querypos to be %d\n",querypos));

    if (hardclip_low == 0) {
      starti = 0;
      hardclip_end_p = true;
    } else if (hardclip_low > querypos) {
      /* startpos = hardclip_low; */
      starti = hardclip_low - querypos;
      hardclip_end_p = true;
      debug2(printf("  Setting starti %d = hardclip_low %d - querypos %d\n",
		    starti,hardclip_low,querypos));
    } else {
      /* startpos = querypos; */
      starti = 0;
    }

    if (querylength - hardclip_high < querypos + stringlength) {
      /* endpos = querylength - hardclip_high; */
      endi = (querylength - hardclip_high) - querypos;
      debug2(printf("  Setting endi %d = (querylength %d - hardclip_high %d) - querypos %d\n",
		    endi,querylength,hardclip_high,querypos));
    } else {
      /* endpos = querypos + stringlength; */
      endi = stringlength;
    }

    debug2(printf("  Counting matches from %d to %d\n",starti,endi));

    if (genomicfwd_refdiff == NULL) {
      if (endi > starti) {
	matchlength += (endi - starti);
      }

    } else if (md_lowercase_variant_p == false) {
      for (i = starti; i < endi; i++) {
	if (isupper(genomicfwd_refdiff[i])) {
	  matchlength++;

	} else {
	  if (matchlength > 0 || hardclip_end_p == true) {
	    FPRINTF(fp,"%d",matchlength);
	    *printp = true;
	    hardclip_end_p = false;
	  }
	  FPRINTF(fp,"%c",toupper(genomicfwd_refdiff[i]));
	  *printp = true;
	  local_nmismatches += 1;
	  matchlength = 0;
	}
      }
      *nmismatches_refdiff += local_nmismatches;

    } else {
      for (i = starti; i < endi; i++) {
	if (isupper(genomicfwd_refdiff[i])) {
	  matchlength++;

	} else if (isupper(genomicfwd_bothdiff[i])) {
	  /* A mismatch against the reference only => alternate variant */
	  if (matchlength > 0 || hardclip_end_p == true) {
	    FPRINTF(fp,"%d",matchlength);
	    *printp = true;
	    hardclip_end_p = false;
	  }
	  FPRINTF(fp,"%c",genomicfwd_refdiff[i]); /* Leave as lower case */
	  *printp = true;
	  local_nmismatches += 1;
	  matchlength = 0;

	} else {
	  /* A true mismatch against both variants */
	  if (matchlength > 0 || hardclip_end_p == true) {
	    FPRINTF(fp,"%d",matchlength);
	    *printp = true;
	    hardclip_end_p = false;
	  }
	  FPRINTF(fp,"%c",toupper(genomicfwd_refdiff[i]));
	  *printp = true;
	  local_nmismatches += 1;
	  matchlength = 0;
	}
      }
      *nmismatches_refdiff += local_nmismatches;
    }
  }

  /* Update nmismatches_bothdiff */
  if (genomicfwd_bothdiff == NULL) {
    /* No change to nmismatches_bothdiff */
  } else if (genomicfwd_bothdiff == genomicfwd_refdiff) {
    *nmismatches_bothdiff += local_nmismatches;
  } else {
    for (i = starti; i < endi; i++) {
      if (!isupper(genomicfwd_bothdiff[i])) {
	*nmismatches_bothdiff += 1;
      }
    }
  }

  debug2(printf("  Ending with matchlength %d\n",matchlength));

  if (lastp == false) {
    return matchlength;
  } else if (matchlength > 0) {
    FPRINTF(fp,"%d",matchlength);
    *printp = true;
    return 0;
  } else {
    return 0;
  }
}


/* Copy also in pair.c for GMAP */
static bool
check_cigar_types (Intlist_T cigar_types) {
  Intlist_T p;
  int type;
  bool M_present_p = false;

  for (p = cigar_types; p != NULL; p = Intlist_next(p)) {
    type = Intlist_head(p);
    if (type == 'M') {
      M_present_p = true;
#if 0
    } else if (type == 'H' && last_type == 'S') {
      debug1(printf("check_cigar_types detects adjacent S and H, so returning false\n"));
      return false;
    } else if (type == 'S' && last_type == 'H') {
      debug1(printf("check_cigar_types detects adjacent S and H, so returning false\n"));
      return false;
#endif
    }
  }

  return M_present_p;
}



static void
print_substrings (Filestring_T fp, char *abbrev, Stage3end_T stage3end, Stage3end_T mate,
		  char *acc1, char *acc2, int pathnum, int npaths,
		  int absmq_score, int first_absmq, int second_absmq, int mapq_score,
		  Shortread_T queryseq, int pairedlength,
		  Chrpos_T chrpos, Chrpos_T mate_chrpos, int hardclip_low, int hardclip_high,
		  Resulttype_T resulttype, bool first_read_p, bool artificial_mate_p, int npaths_mate,
		  int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p,
		  bool circularp) {
  unsigned int flag = 0U;
  Substring_T substring, substringL, substringH, substringM;
  Junction_T post_junction;
  int type;
  int nindels = 0;

  List_T substrings_LtoH, junctions_LtoH;
  List_T startp, endp, startq, prevp, finalp, nextp, p, q;
  int substring_start, substring_length, matchlength;

  int nmismatches_refdiff = 0, nmismatches_bothdiff = 0, querylength;
  int sensedir;
  char *genomicfwd_refdiff, *genomicfwd_bothdiff, *genomicdir_refdiff, *genomicdir_bothdiff;
  char *deletion_string;
  bool plusp, lastp, printp;
  bool ambigL, ambigH;
  int n, i;
  Univcoord_T *ambcoords, splicecoord;
#ifdef PRINT_AMBIG_COORDS
  Univcoord_T chroffset;
#endif

  
  querylength = Shortread_fulllength(queryseq);
  plusp = Stage3end_plusp(stage3end);

  if ((sensedir = Stage3end_sensedir(stage3end)) == SENSE_NULL && mate != NULL) {
    sensedir = Stage3end_sensedir(mate);
  }
  /* sensep = (sensedir == SENSE_ANTI) ? false : true; */

  /* 1. QNAME */
  if (acc2 == NULL) {
    FPRINTF(fp,"%s",acc1);
  } else {
    FPRINTF(fp,"%s,%s",acc1,acc2);
  }

  /* 2. FLAG */
  flag = SAM_compute_flag(plusp,mate,resulttype,first_read_p,
			  pathnum,npaths,artificial_mate_p,npaths_mate,
			  absmq_score,first_absmq,invertp,invert_mate_p);
  FPRINTF(fp,"\t%u",flag);

  /* 3. RNAME: chr */
  /* 4. POS: chrpos */
  print_chromosomal_pos(fp,Stage3end_chrnum(stage3end),chrpos,Stage3end_chrlength(stage3end),chromosome_iit);


  /* 5. MAPQ: Mapping quality */
  FPRINTF(fp,"\t%d",mapq_score);

  /* 6. CIGAR */
  FPRINTF(fp,"\t");
  substrings_LtoH = Stage3end_substrings_LtoH(stage3end);
  junctions_LtoH = Stage3end_junctions_LtoH(stage3end);
  substringL = (Substring_T) List_head(substrings_LtoH);
  substringH = (Substring_T) List_last_value(substrings_LtoH);
  if (Substring_ambiguous_p(substringL) == true) {
    prevp = substrings_LtoH;
    startp = List_next(substrings_LtoH);
    startq = List_next(junctions_LtoH);
  } else {
    prevp = (List_T) NULL;
    startp = substrings_LtoH;
    startq = junctions_LtoH;
  }
  if (Substring_ambiguous_p(substringH) == true) {
    endp = List_last_item(substrings_LtoH);
  } else {
    endp = (List_T) NULL;
  }

  debug(printf("End has %d substrings\n",List_length(substrings_LtoH)));

  p = startp;
  q = startq;
  if (plusp == true) {
    /* Plus */
    while (p != endp && Substring_queryend((Substring_T) List_head(p)) < hardclip_low) {
      /* Skip, because substring entirely in hard-clipped region */
      debug(printf("Skipping %d..%d\n",Substring_querystart((Substring_T) List_head(p)),
		   Substring_queryend((Substring_T) List_head(p))));
      prevp = p;
      p = List_next(p);
      q = List_next(q);
    }

    substring = (Substring_T) List_head(p);
    if (List_next(p) == endp ||	Substring_queryend(substring) >= querylength - hardclip_high) {
      /* Single substring */
      debug(printf("Single substring %d..%d\n",Substring_querystart((Substring_T) List_head(p)),
		   Substring_queryend((Substring_T) List_head(p))));

      if (hide_soft_clips_p == true) {
	print_cigar(fp,/*type*/'M',
		    Substring_querystart(substring) + Substring_match_length(substring) +
		    (querylength - Substring_queryend(substring)),/*querypos*/0,querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true,/*trimlength*/0);
      } else {
	print_cigar(fp,/*type*/'S',Substring_querystart(substring),
		    /*querypos*/0,querylength,hardclip_low,hardclip_high,
		    /*plusp*/true,/*lastp*/false,/*trimlength*/0);
	print_cigar(fp,/*type*/'M',Substring_match_length(substring),
		    /*querypos*/Substring_querystart(substring),querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false,/*trimlength*/0);
	print_cigar(fp,/*type*/'S',querylength - Substring_queryend(substring),
		    /*querypos*/Substring_queryend(substring),querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true,/*trimlength*/0);
      }
      finalp = p;
      nextp = List_next(p);

    } else {
      /* First substring, plus */
      debug(printf("First substring, plus %d..%d\n",Substring_querystart((Substring_T) List_head(p)),
		   Substring_queryend((Substring_T) List_head(p))));

      post_junction = (Junction_T) List_head(q);

      if (hide_soft_clips_p == true) {
	print_cigar(fp,/*type*/'M',
		    Substring_querystart(substring) +
		    Substring_match_length(substring),
		    /*querypos*/0,querylength,hardclip_low,hardclip_high,
		    /*plusp*/true,/*lastp*/false,/*trimlength*/0);
      } else {
	print_cigar(fp,/*type*/'S',Substring_querystart(substring),
		    /*querypos*/0,querylength,hardclip_low,hardclip_high,
		    /*plusp*/true,/*lastp*/false,/*trimlength*/0);
	print_cigar(fp,/*type*/'M',Substring_match_length(substring),
		    /*querypos*/Substring_querystart(substring),querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false,/*trimlength*/0);
      }
      p = List_next(p);
      
      while (p != endp && Substring_queryend((Substring_T) List_head(p)) < querylength - hardclip_high) {
	if ((type = Junction_type(post_junction)) == DEL_JUNCTION) {
	  FPRINTF(fp,"%dD",Junction_nindels(post_junction));
	  nindels += Junction_nindels(post_junction);
	} else if (type == INS_JUNCTION) {
	  FPRINTF(fp,"%dI",Junction_nindels(post_junction));
	  nindels += Junction_nindels(post_junction);
	} else if (type == SPLICE_JUNCTION) {
	  FPRINTF(fp,"%uN",Junction_splice_distance(post_junction));
	}
	q = List_next(q);
	if (q == NULL) {
	} else {
	  post_junction = (Junction_T) List_head(q);
	}

	substring = (Substring_T) List_head(p);
	if (List_next(p) == endp) {
	  /* Last substring, plus, not hard-clipped */
	  debug(printf("Last substring, plus, not hard-clipped %d..%d\n",Substring_querystart((Substring_T) List_head(p)),
		       Substring_queryend((Substring_T) List_head(p))));
	  
	  if (hide_soft_clips_p == true) {
	    print_cigar(fp,/*type*/'M',
			Substring_match_length(substring) +
			(querylength - Substring_queryend(substring)),
			/*querypos*/Substring_querystart(substring),querylength,
			hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true,/*trimlength*/0);
	  } else {
	    print_cigar(fp,/*type*/'M',Substring_match_length(substring),
			/*querypos*/Substring_querystart(substring),querylength,
			hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false,/*trimlength*/0);
	    print_cigar(fp,/*type*/'S',querylength - Substring_queryend(substring),
			/*querypos*/Substring_queryend(substring),querylength,
			hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true,/*trimlength*/0);
	  }
	  finalp = p;
	  nextp = List_next(p);

	} else {
	  /* Middle substring, plus */
	  debug(printf("Middle substring, plus %d..%d\n",Substring_querystart((Substring_T) List_head(p)), 
		       Substring_queryend((Substring_T) List_head(p))));

	  print_cigar(fp,/*type*/'M',Substring_match_length(substring),
		      /*querypos*/Substring_querystart(substring),querylength,
		      hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false,/*trimlength*/0);
	}
	p = List_next(p);
      }
      
      if (p != endp) {
	if ((type = Junction_type(post_junction)) == DEL_JUNCTION) {
	  FPRINTF(fp,"%dD",Junction_nindels(post_junction));
	  nindels += Junction_nindels(post_junction);
	} else if (type == INS_JUNCTION) {
	  FPRINTF(fp,"%dI",Junction_nindels(post_junction));
	  nindels += Junction_nindels(post_junction);
	} else if (type == SPLICE_JUNCTION) {
	  FPRINTF(fp,"%uN",Junction_splice_distance(post_junction));
	}

	/* Last substring, plus, hard-clipped */
	substring = (Substring_T) List_head(p);
	debug(printf("Last substring, plus, hard-clipped %d..%d\n",Substring_querystart((Substring_T) List_head(p)),
		     Substring_queryend((Substring_T) List_head(p))));
	if (hide_soft_clips_p == true) {
	  print_cigar(fp,/*type*/'M',
		      Substring_match_length(substring) +
		      (querylength - Substring_queryend(substring)),
		      /*querypos*/Substring_querystart(substring),querylength,
		      hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true,/*trimlength*/0);
	} else {
	  print_cigar(fp,/*type*/'M',Substring_match_length(substring),
		      /*querypos*/Substring_querystart(substring),querylength,
		      hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false,/*trimlength*/0);
	  print_cigar(fp,/*type*/'S',querylength - Substring_queryend(substring),
		      /*querypos*/Substring_queryend(substring),querylength,
		      hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true,/*trimlength*/0);
	}
	finalp = p;
	nextp = List_next(p);

      }
    }

  } else {
    /* Minus */
    while (p != endp && Substring_querystart((Substring_T) List_head(p)) >= querylength - hardclip_low) {
      /* Skip, because substring entirely in hard-clipped region */
      debug(printf("Skipping %d..%d\n",Substring_querystart((Substring_T) List_head(p)),
		   Substring_queryend((Substring_T) List_head(p))));
      prevp = p;
      p = List_next(p);
      q = List_next(q);
    }

    substring = (Substring_T) List_head(p);
    if (List_next(p) == endp || Substring_querystart(substring) < hardclip_high) {
      /* Single substring */
      debug(printf("Single substring %d..%d\n",Substring_querystart((Substring_T) List_head(p)),
		   Substring_queryend((Substring_T) List_head(p))));

      if (hide_soft_clips_p == true) {
	print_cigar(fp,/*type*/'M',
		    (querylength - Substring_queryend(substring)) + 
		    Substring_match_length(substring) + Substring_querystart(substring),
		    /*querypos*/querylength,querylength,hardclip_low,hardclip_high,
		    /*plusp*/false,/*lastp*/true,/*trimlength*/0);
      } else {
	print_cigar(fp,/*type*/'S',querylength - Substring_queryend(substring),
		    /*querypos*/querylength,querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false,/*trimlength*/0);
	print_cigar(fp,/*type*/'M',Substring_match_length(substring),
		    /*querypos*/Substring_queryend(substring),querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false,/*trimlength*/0);
	print_cigar(fp,/*type*/'S',Substring_querystart(substring),
		    /*querypos*/Substring_querystart(substring),querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true,/*trimlength*/0);
      }
      finalp = p;
      nextp = List_next(p);

    } else {
      /* First substring, minus */
      debug(printf("First substring, minus %d..%d\n",Substring_querystart((Substring_T) List_head(p)),
		   Substring_queryend((Substring_T) List_head(p))));
    
      post_junction = (Junction_T) List_head(q);

      if (hide_soft_clips_p == true) {
	print_cigar(fp,/*type*/'M',
		    (querylength - Substring_queryend(substring)) +
		    Substring_match_length(substring),
		    /*querypos*/querylength,querylength,hardclip_low,hardclip_high,
		    /*plusp*/false,/*lastp*/false,/*trimlength*/0);
      } else {
	print_cigar(fp,/*type*/'S',querylength - Substring_queryend(substring),
		    /*querypos*/querylength,querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false,/*trimlength*/0);
	print_cigar(fp,/*type*/'M',Substring_match_length(substring),
		    /*querypos*/Substring_queryend(substring),querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false,/*trimlength*/0);
      }
      p = List_next(p);

      while (p != endp && Substring_querystart((Substring_T) List_head(p)) >= hardclip_high) {
	if ((type = Junction_type(post_junction)) == DEL_JUNCTION) {
	  FPRINTF(fp,"%dD",Junction_nindels(post_junction));
	  nindels += Junction_nindels(post_junction);
	} else if (type == INS_JUNCTION) {
	  FPRINTF(fp,"%dI",Junction_nindels(post_junction));
	  nindels += Junction_nindels(post_junction);
	} else if (type == SPLICE_JUNCTION) {
	  FPRINTF(fp,"%uN",Junction_splice_distance(post_junction));
	}
	q = List_next(q);
	if (q == NULL) {
	} else {
	  post_junction = (Junction_T) List_head(q);
	}

	substring = (Substring_T) List_head(p);
	if (List_next(p) == endp) {
	  /* Last substring, minus, not hard-clipped */
	  debug(printf("Last substring, minus, not hard-clipped %d..%d\n",Substring_querystart((Substring_T) List_head(p)),
		       Substring_queryend((Substring_T) List_head(p))));

	  if (hide_soft_clips_p == true) {
	    print_cigar(fp,/*type*/'M',
			Substring_match_length(substring) +
			Substring_querystart(substring),
			/*querypos*/Substring_queryend(substring),querylength,
			hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true,/*trimlength*/0);
	  } else {
	    print_cigar(fp,/*type*/'M',Substring_match_length(substring),
			/*querypos*/Substring_queryend(substring),querylength,
			hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false,/*trimlength*/0);
	    print_cigar(fp,/*type*/'S',Substring_querystart(substring),
			/*querypos*/Substring_querystart(substring),querylength,hardclip_low,hardclip_high,
			/*plusp*/false,/*lastp*/true,/*trimlength*/0);
	  }
	  finalp = p;
	  nextp = List_next(p);

	} else {
	  /* Middle substring, minus */
	  debug(printf("Middle substring, minus %d..%d\n",Substring_querystart((Substring_T) List_head(p)),
		       Substring_queryend((Substring_T) List_head(p))));

	  print_cigar(fp,/*type*/'M',Substring_match_length(substring),
		      /*querypos*/Substring_queryend(substring),querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false,/*trimlength*/0);
	}
	p = List_next(p);
      }

      if (p != endp) {
	if ((type = Junction_type(post_junction)) == DEL_JUNCTION) {
	  FPRINTF(fp,"%dD",Junction_nindels(post_junction));
	  nindels += Junction_nindels(post_junction);
	} else if (type == INS_JUNCTION) {
	  FPRINTF(fp,"%dI",Junction_nindels(post_junction));
	  nindels += Junction_nindels(post_junction);
	} else if (type == SPLICE_JUNCTION) {
	  FPRINTF(fp,"%uN",Junction_splice_distance(post_junction));
	}

	/* Last substring, minus, hard-clipped */
	substring = (Substring_T) List_head(p);
	debug(printf("Last substring, minus, hard-clipped %d..%d\n",Substring_querystart((Substring_T) List_head(p)),
		     Substring_queryend((Substring_T) List_head(p))));

	if (hide_soft_clips_p == true) {
	  print_cigar(fp,/*type*/'M',
		      Substring_match_length(substring) +
		      Substring_querystart(substring),
		      /*querypos*/Substring_queryend(substring),querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true,/*trimlength*/0);
	} else {
	  print_cigar(fp,/*type*/'M',Substring_match_length(substring),
		      /*querypos*/Substring_queryend(substring),querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false,/*trimlength*/0);
	  print_cigar(fp,/*type*/'S',Substring_querystart(substring),
		      /*querypos*/Substring_querystart(substring),querylength,hardclip_low,hardclip_high,
		      /*plusp*/false,/*lastp*/true,/*trimlength*/0);
	}
	finalp = p;
	nextp = List_next(p);

      }
    }
  }


  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  print_mate_chromosomal_pos(fp,Stage3end_chrnum(mate),Stage3end_effective_chrnum(mate),
			     mate_chrpos,Stage3end_chrlength(mate),
			     Stage3end_chrnum(stage3end),chrpos,chromosome_iit);


  /* 9. ISIZE: Insert size */
  if (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_TRANSLOC || resulttype == CONCORDANT_MULT) {
    if (plusp == invertp) {
      FPRINTF(fp,"\t%d",-pairedlength);
    } else {
      FPRINTF(fp,"\t%d",pairedlength);
    }
  } else if (mate_chrpos == 0) {
    FPRINTF(fp,"\t%d",pairedlength);
  } else if (chrpos < mate_chrpos) {
    FPRINTF(fp,"\t%d",pairedlength);
  } else if (chrpos > mate_chrpos) {
    FPRINTF(fp,"\t%d",-pairedlength);
  } else if (first_read_p == true) {
    FPRINTF(fp,"\t%d",pairedlength);
  } else {
    FPRINTF(fp,"\t%d",-pairedlength);
  }


  /* 10. SEQ: queryseq and 11. QUAL: quality scores */
  /* Queryseq has already been inverted, so just measure plusp relative to its current state */
  if (plusp == true) {
    Shortread_print_chopped_sam(fp,queryseq,hardclip_low,hardclip_high);
    FPRINTF(fp,"\t");
    Shortread_print_quality(fp,queryseq,hardclip_low,hardclip_high,
			   quality_shift,/*show_chopped_p*/false);
  } else {
    Shortread_print_chopped_revcomp_sam(fp,queryseq,hardclip_low,hardclip_high);
    FPRINTF(fp,"\t");
    Shortread_print_quality_revcomp(fp,queryseq,hardclip_low,hardclip_high,
				   quality_shift,/*show_chopped_p*/false);
  } 

  /* 12. TAGS: RG */
  if (sam_read_group_id != NULL) {
    FPRINTF(fp,"\tRG:Z:%s",sam_read_group_id);
  }

  /* 12. TAGS: XH and XI */
  if (hardclip_low > 0 || hardclip_high > 0) {
    FPRINTF(fp,"\tXH:Z:");
    if (plusp == true) {
      Shortread_print_chopped_end(fp,queryseq,hardclip_low,hardclip_high);
    } else {
      Shortread_print_chopped_end_revcomp(fp,queryseq,hardclip_low,hardclip_high);
    }

    if (Shortread_quality_string(queryseq) != NULL) {
      FPRINTF(fp,"\tXI:Z:");
      if (plusp == true) {
	Shortread_print_chopped_end_quality(fp,queryseq,hardclip_low,hardclip_high);
      } else {
	Shortread_print_chopped_end_quality_reverse(fp,queryseq,hardclip_low,hardclip_high);
      }
    }
  }

  /* 12. TAGS: XB */
  Shortread_print_barcode(fp,queryseq);

  /* 12. TAGS: XP.  Logically should be last in reconstructing a read. */
  Shortread_print_chop(fp,queryseq,invertp);

  /* 12. TAGS: MD */
  FPRINTF(fp,"\tMD:Z:");
  p = startp;
  q = startq;
  printp = false;

  if (plusp == true) {
    /* Plus */
    while (p != endp && Substring_queryend((Substring_T) List_head(p)) < hardclip_low) {
      /* Skip, because substring entirely in hard-clipped region */
      p = List_next(p);
      q = List_next(q);
    }

    substring = (Substring_T) List_head(p);
    if (List_next(p) == endp || Substring_queryend(substring) >= querylength - hardclip_high) {
      /* Single substring */
      if (hide_soft_clips_p == true) {
	substring_start = Substring_querystart_orig(substring);
	substring_length = Substring_match_length_orig(substring);
      } else {
	substring_start = Substring_querystart(substring);
	substring_length = Substring_match_length(substring);
      }

      if ((genomicfwd_bothdiff = Substring_genomic_bothdiff(substring)) == NULL) {
	/* matchlength = */ print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,/*matchlength*/0,
					    /*genomicfwd_refdiff*/NULL,/*genomicfwd_bothdiff*/NULL,
					    substring_length,/*querypos*/substring_start,querylength,
					    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);		    
      } else {
	genomicfwd_refdiff = Substring_genomic_refdiff(substring);
	/* matchlength = */ print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,/*matchlength*/0,
					    &(genomicfwd_refdiff[substring_start]),&(genomicfwd_bothdiff[substring_start]),
					    substring_length,/*querypos*/substring_start,querylength,
					    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
      }
	
    } else {
      /* First substring, plus */
      if (hide_soft_clips_p == true) {
	substring_start = Substring_querystart_orig(substring);
	substring_length = Substring_match_length_orig(substring);
      } else {
	substring_start = Substring_querystart(substring);
	substring_length = Substring_match_length(substring);
      }

      post_junction = (Junction_T) List_head(q);
      if ((type = Junction_type(post_junction)) == DEL_JUNCTION) {
	lastp = true;
      } else {
	lastp = false;
      }

      if ((genomicfwd_bothdiff = Substring_genomic_bothdiff(substring)) == NULL) {
	matchlength = print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,/*matchlength*/0,
				      /*genomicfwd_refdiff*/NULL,/*genomicfwd_bothdiff*/NULL,
				      substring_length,/*querypos*/substring_start,querylength,
				      hardclip_low,hardclip_high,/*plusp*/true,lastp);
      } else {
	genomicfwd_refdiff = Substring_genomic_refdiff(substring);
	matchlength = print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,/*matchlength*/0,
				      &(genomicfwd_refdiff[substring_start]),&(genomicfwd_bothdiff[substring_start]),
				      substring_length,/*querypos*/substring_start,querylength,
				      hardclip_low,hardclip_high,/*plusp*/true,lastp);
      }
      p = List_next(p);
      
      while (p != endp && Substring_queryend((Substring_T) List_head(p)) < querylength - hardclip_high) {
	if (type == DEL_JUNCTION) {
	  deletion_string = Junction_deletion_string(post_junction,genome,/*plusp*/true); 
	  FPRINTF(fp,"^%s",deletion_string);
	  FREE(deletion_string);
	}
	q = List_next(q);
	if (q == NULL) {
	  lastp = true;
	} else {
	  post_junction = (Junction_T) List_head(q);
	  if ((type = Junction_type(post_junction)) == DEL_JUNCTION) {
	    lastp = true;
	  } else {
	    lastp = false;
	  }
	}

	substring = (Substring_T) List_head(p);
	if (List_next(p) == endp) {
	  /* Last substring, plus, not hard-clipped */
	  if (hide_soft_clips_p == true) {
	    substring_start = Substring_querystart_orig(substring);
	    substring_length = Substring_match_length_orig(substring);
	  } else {
	    substring_start = Substring_querystart(substring);
	    substring_length = Substring_match_length(substring);
	  }
	  
	  if ((genomicfwd_bothdiff = Substring_genomic_bothdiff(substring)) == NULL) {
	    /* matchlength = */ print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,matchlength,
					       /*genomicfwd_refdiff*/NULL,/*genomicfwd_bothdiff*/NULL,
					       substring_length,/*querypos*/substring_start,querylength,
					       hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
	  } else {
	    genomicfwd_refdiff = Substring_genomic_refdiff(substring);
	    /* matchlength = */ print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,matchlength,
						&(genomicfwd_refdiff[substring_start]),&(genomicfwd_bothdiff[substring_start]),
						substring_length,/*querypos*/substring_start,querylength,
						hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
	  }

	} else {
	  /* Middle substring, plus */
	  substring_start = Substring_querystart(substring);
	  substring_length = Substring_match_length(substring);

	  if ((genomicfwd_bothdiff = Substring_genomic_bothdiff(substring)) == NULL) {
	    matchlength = print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,matchlength,
					  /*genomicfwd_refdiff*/NULL,/*genomicfwd_bothdiff*/NULL,
					  substring_length,/*querypos*/substring_start,querylength,
					  hardclip_low,hardclip_high,/*plusp*/true,lastp);
	  } else {
	    genomicfwd_refdiff = Substring_genomic_refdiff(substring);
	    matchlength = print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,matchlength,
					  &(genomicfwd_refdiff[substring_start]),&(genomicfwd_bothdiff[substring_start]),
					  substring_length,/*querypos*/substring_start,querylength,
					  hardclip_low,hardclip_high,/*plusp*/true,lastp);
	  }
	}
	p = List_next(p);
      }
      
      if (p != endp) {
	if (type == DEL_JUNCTION) {
	  deletion_string = Junction_deletion_string(post_junction,genome,/*plusp*/true); 
	  FPRINTF(fp,"^%s",deletion_string);
	  FREE(deletion_string);
	}

	/* Last substring, plus, hard-clipped */
	substring = (Substring_T) List_head(p);
	if (hide_soft_clips_p == true) {
	  substring_start = Substring_querystart_orig(substring);
	  substring_length = Substring_match_length_orig(substring);
	} else {
	  substring_start = Substring_querystart(substring);
	  substring_length = Substring_match_length(substring);
	}

	if ((genomicfwd_bothdiff = Substring_genomic_bothdiff(substring)) == NULL) {
	  /* matchlength = */ print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,matchlength,
					      /*genomicfwd_refdiff*/NULL,/*genomicfwd_bothdiff*/NULL,
					      substring_length,/*querypos*/substring_start,querylength,
					      hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
	} else {
	  genomicfwd_refdiff = Substring_genomic_refdiff(substring);
	  /* matchlength = */ print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,matchlength,
					      &(genomicfwd_refdiff[substring_start]),&(genomicfwd_bothdiff[substring_start]),
					      substring_length,/*querypos*/substring_start,querylength,
					      hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
	}
      }
    }

  } else {
    /* Minus */
    while (p != endp && Substring_querystart((Substring_T) List_head(p)) >= querylength - hardclip_low) {
      /* Skip, because substring entirely in hard-clipped region */
      p = List_next(p);
      q = List_next(q);
    }

    substring = (Substring_T) List_head(p);
    if (List_next(p) == endp ||	querylength - Substring_queryend(substring) >= querylength - hardclip_high) {
      /* Single substring */
      if (hide_soft_clips_p == true) {
	substring_start = Substring_querystart_orig(substring);
	substring_length = Substring_match_length_orig(substring);
      } else {
	substring_start = Substring_querystart(substring);
	substring_length = Substring_match_length(substring);
      }

      if ((genomicdir_bothdiff = Substring_genomic_bothdiff(substring)) == NULL) {
	/* matchlength = */ print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
					    fp,/*matchlength*/0,/*genomicfwd_refdiff*/NULL,/*genomicfwd_bothdiff*/NULL,
					    substring_length,/*querypos*/substring_start,querylength,
					    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      } else if ((genomicdir_refdiff = Substring_genomic_refdiff(substring)) == genomicdir_bothdiff) {
	genomicfwd_refdiff = (char *) MALLOCA((substring_length+1) * sizeof(char));
	make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
	/* matchlength = */ print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
					    fp,/*matchlength*/0,genomicfwd_refdiff,/*genomicfwd_bothdiff*/genomicfwd_refdiff,
					    substring_length,/*querypos*/substring_start,querylength,
					    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
	FREEA(genomicfwd_refdiff);
      } else {
	genomicfwd_refdiff = (char *) MALLOCA((substring_length+1) * sizeof(char));
	genomicfwd_bothdiff = (char *) MALLOCA((substring_length+1) * sizeof(char));
	make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
	make_complement_buffered(genomicfwd_bothdiff,&(genomicdir_bothdiff[substring_start]),substring_length);
	/* matchlength = */ print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
					    fp,/*matchlength*/0,genomicfwd_refdiff,genomicfwd_bothdiff,
					    substring_length,/*querypos*/substring_start,querylength,
					    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
	FREEA(genomicfwd_bothdiff);
	FREEA(genomicfwd_refdiff);
      }

    } else {
      /* First substring, minus */
      if (hide_soft_clips_p == true) {
	substring_start = Substring_querystart_orig(substring);
	substring_length = Substring_match_length_orig(substring);
      } else {
	substring_start = Substring_querystart(substring);
	substring_length = Substring_match_length(substring);
      }

      post_junction = (Junction_T) List_head(q);
      if ((type = Junction_type(post_junction)) == DEL_JUNCTION) {
	lastp = true;
      } else {
	lastp = false;
      }

      if ((genomicdir_bothdiff = Substring_genomic_bothdiff(substring)) == NULL) {
	matchlength = print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
				      fp,/*matchlength*/0,/*genomicfwd_refdiff*/NULL,/*genomicfwd_bothdiff*/NULL,
				      substring_length,/*querypos*/substring_start,querylength,
				      hardclip_low,hardclip_high,/*plusp*/false,lastp);
      } else if ((genomicdir_refdiff = Substring_genomic_refdiff(substring)) == genomicdir_bothdiff) {
	genomicfwd_refdiff = (char *) MALLOCA((substring_length+1) * sizeof(char));
	make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
	matchlength = print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
				      fp,/*matchlength*/0,genomicfwd_refdiff,/*genomicfwd_bothdiff*/genomicfwd_refdiff,
				      substring_length,/*querypos*/substring_start,querylength,
				      hardclip_low,hardclip_high,/*plusp*/false,lastp);
	FREEA(genomicfwd_refdiff);
      } else {
	genomicfwd_refdiff = (char *) MALLOCA((substring_length+1) * sizeof(char));
	genomicfwd_bothdiff = (char *) MALLOCA((substring_length+1) * sizeof(char));
	make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
	make_complement_buffered(genomicfwd_bothdiff,&(genomicdir_bothdiff[substring_start]),substring_length);
	matchlength = print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
				      fp,/*matchlength*/0,genomicfwd_refdiff,genomicfwd_bothdiff,
				      substring_length,/*querypos*/substring_start,querylength,
				      hardclip_low,hardclip_high,/*plusp*/false,lastp);
	FREEA(genomicfwd_bothdiff);
	FREEA(genomicfwd_refdiff);
      }
      p = List_next(p);

      while (p != endp && querylength - Substring_queryend((Substring_T) List_head(p)) < querylength - hardclip_high) {
	if (type == DEL_JUNCTION) {
	  deletion_string = Junction_deletion_string(post_junction,genome,/*plusp:true*/true); 
	  FPRINTF(fp,"^%s",deletion_string);
	  FREE(deletion_string);
	}
	q = List_next(q);
	if (q == NULL) {
	  lastp = true;
	} else {
	  post_junction = (Junction_T) List_head(q);
	  if ((type = Junction_type(post_junction)) == DEL_JUNCTION) {
	    lastp = true;
	  } else {
	    lastp = false;
	  }
	}

	substring = (Substring_T) List_head(p);
	if (List_next(p) == endp) {
	  /* Last substring, minus, not hard-clipped */
	  if (hide_soft_clips_p == true) {
	    substring_start = Substring_querystart_orig(substring);
	    substring_length = Substring_match_length_orig(substring);
	  } else {
	    substring_start = Substring_querystart(substring);
	    substring_length = Substring_match_length(substring);
	  }

	  if ((genomicdir_bothdiff = Substring_genomic_bothdiff(substring)) == NULL) {
	    /* matchlength = */ print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
						fp,matchlength,/*genomicfwd_refdiff*/NULL,/*genomicfwd_bothdiff*/NULL,
						substring_length,/*querypos*/substring_start,querylength,
						hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
	  } else if ((genomicdir_refdiff = Substring_genomic_refdiff(substring)) == genomicdir_bothdiff) {
	    genomicfwd_refdiff = (char *) MALLOCA((substring_length+1) * sizeof(char));
	    make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
	    /* matchlength = */ print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
						fp,matchlength,genomicfwd_refdiff,/*genomicfwd_bothdiff*/genomicfwd_refdiff,
						substring_length,/*querypos*/substring_start,querylength,
						hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
	    FREEA(genomicfwd_refdiff);
	  } else {
	    genomicfwd_refdiff = (char *) MALLOCA((substring_length+1) * sizeof(char));
	    genomicfwd_bothdiff = (char *) MALLOCA((substring_length+1) * sizeof(char));
	    make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
	    make_complement_buffered(genomicfwd_bothdiff,&(genomicdir_bothdiff[substring_start]),substring_length);
	    /* matchlength = */ print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
						fp,matchlength,genomicfwd_refdiff,genomicfwd_bothdiff,
						substring_length,/*querypos*/substring_start,querylength,
						hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
	    FREEA(genomicfwd_bothdiff);
	    FREEA(genomicfwd_refdiff);
	  }

	} else {
	  /* Middle substring, minus */
	  substring_start = Substring_querystart(substring);
	  substring_length = Substring_match_length(substring);

	  if ((genomicdir_bothdiff = Substring_genomic_bothdiff(substring)) == NULL) {
	    matchlength = print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
					  fp,matchlength,/*genomicfwd_refdiff*/NULL,/*genomicfwd_bothdiff*/NULL,
					  substring_length,/*querypos*/substring_start,querylength,
					  hardclip_low,hardclip_high,/*plusp*/false,lastp);
	  } else if ((genomicdir_refdiff = Substring_genomic_refdiff(substring)) == genomicdir_bothdiff) {
	    genomicfwd_refdiff = (char *) MALLOCA((substring_length+1) * sizeof(char));
	    make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
	    matchlength = print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
					  fp,matchlength,genomicfwd_refdiff,/*genomicfwd_bothdiff*/genomicfwd_refdiff,
					  substring_length,/*querypos*/substring_start,querylength,
					  hardclip_low,hardclip_high,/*plusp*/false,lastp);
	    FREEA(genomicfwd_refdiff);
	  } else {
	    genomicfwd_refdiff = (char *) MALLOCA((substring_length+1) * sizeof(char));
	    genomicfwd_bothdiff = (char *) MALLOCA((substring_length+1) * sizeof(char));
	    make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
	    make_complement_buffered(genomicfwd_bothdiff,&(genomicdir_bothdiff[substring_start]),substring_length);
	    matchlength = print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
					  fp,matchlength,genomicfwd_refdiff,genomicfwd_bothdiff,
					  substring_length,/*querypos*/substring_start,querylength,
					  hardclip_low,hardclip_high,/*plusp*/false,lastp);
	    FREEA(genomicfwd_bothdiff);
	    FREEA(genomicfwd_refdiff);
	  }
	}
	p = List_next(p);
      }

      if (p != endp) {
	if (type == DEL_JUNCTION) {
	  deletion_string = Junction_deletion_string(post_junction,genome,/*plusp:true*/true); 
	  FPRINTF(fp,"^%s",deletion_string);
	  FREE(deletion_string);
	}

	/* Last substring, minus, hard-clipped */
	substring = (Substring_T) List_head(p);
	if (hide_soft_clips_p == true) {
	  substring_start = Substring_querystart_orig(substring);
	  substring_length = Substring_match_length_orig(substring);
	} else {
	  substring_start = Substring_querystart(substring);
	  substring_length = Substring_match_length(substring);
	}

	if ((genomicdir_bothdiff = Substring_genomic_bothdiff(substring)) == NULL) {
	  /* matchlength = */ print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
					      fp,matchlength,/*genomicfwd_refdiff*/NULL,/*genomicfwd_bothdiff*/NULL,
					      substring_length,/*querypos*/substring_start,querylength,
					      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
	} else if ((genomicdir_refdiff = Substring_genomic_refdiff(substring)) == genomicdir_bothdiff) {
	  genomicfwd_refdiff = (char *) MALLOCA((substring_length+1) * sizeof(char));
	  make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
	  /* matchlength = */ print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
					      fp,matchlength,genomicfwd_refdiff,/*genomicfwd_bothdiff*/genomicfwd_refdiff,
					      substring_length,/*querypos*/substring_start,querylength,
					      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
	  FREEA(genomicfwd_refdiff);
	} else {
	  genomicfwd_refdiff = (char *) MALLOCA((substring_length+1) * sizeof(char));
	  genomicfwd_bothdiff = (char *) MALLOCA((substring_length+1) * sizeof(char));
	  make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
	  make_complement_buffered(genomicfwd_bothdiff,&(genomicdir_bothdiff[substring_start]),substring_length);
	  /* matchlength = */ print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
					      fp,matchlength,genomicfwd_refdiff,genomicfwd_bothdiff,
					      substring_length,/*querypos*/substring_start,querylength,
					      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
	  FREEA(genomicfwd_bothdiff);
	  FREEA(genomicfwd_refdiff);
	}
      }
    }
  }

  if (printp == false) {
    FPRINTF(fp,"0");
  }


  /* 12. TAGS: NH */
  /* 12. TAGS: HI */
  /* 12. TAGS: NM */
  FPRINTF(fp,"\tNH:i:%d\tHI:i:%d\tNM:i:%d",npaths,pathnum,nmismatches_refdiff + nindels);
  
  if (snps_iit) {
    /* 12. TAGS: XW and XV */
    FPRINTF(fp,"\tXW:i:%d",nmismatches_bothdiff);
    FPRINTF(fp,"\tXV:i:%d",nmismatches_refdiff - nmismatches_bothdiff);
  }

  /* 12. TAGS: SM */
  /* 12. TAGS: XQ */
  /* 12. TAGS: X2 */
  FPRINTF(fp,"\tSM:i:%d\tXQ:i:%d\tX2:i:%d",mapq_score,absmq_score,second_absmq);

  /* 12. TAGS: XO */
  FPRINTF(fp,"\tXO:Z:%s",abbrev);

  /* 12. TAGS: XS */
  if (sensedir == SENSE_FORWARD) {
    if (plusp == true) {
      FPRINTF(fp,"\tXS:A:+");
    } else {
      FPRINTF(fp,"\tXS:A:-");
    }
  } else if (sensedir == SENSE_ANTI) {
    if (plusp == true) {
      FPRINTF(fp,"\tXS:A:-");
    } else {
      FPRINTF(fp,"\tXS:A:+");
    }
#if 0
    /* Don't print XS field for SENSE_NULL */
  } else if (force_xs_direction_p == true) {
    FPRINTF(fp,"\tXS:A:+");
  } else {
    FPRINTF(fp,"\tXS:A:?");
#endif
  }

  /* 12. TAGS: XA */
  if (prevp == NULL) {
    /* substringL = (Substring_T) NULL; */
    ambigL = false;
  } else {
    substringL = (Substring_T) List_head(prevp);
    ambigL = Substring_ambiguous_p(substringL);
  }
  if (nextp == NULL) {
    ambigH = false;
  } else {
    substringH = (Substring_T) List_head(nextp);
    ambigH = Substring_ambiguous_p(substringH);
  }

  if (ambigL == true || ambigH == true) {
    FPRINTF(fp,"\tXA:Z:");

    if (ambigL == true) {
      ambcoords = Substring_ambcoords(substringL);
      n = Substring_nambcoords(substringL);
#ifdef PRINT_AMBIG_COORDS
      chroffset = Substring_chroffset(substringL);
      FPRINTF(fp,"%u",ambcoords[0] - chroffset + 1U);
      for (i = 1; i < n; i++) {
	FPRINTF(fp,",%u",ambcoords[i] - chroffset + 1U);
      }
#else
      substringM = (Substring_T) List_head(List_next(prevp));
      if (plusp == true) {
	splicecoord = Substring_alignstart(substringM);
      } else {
	splicecoord = Substring_alignend(substringM);
      }
      FPRINTF(fp,"%u",splicecoord - ambcoords[0]);
      for (i = 1; i < n; i++) {
	FPRINTF(fp,",%u",splicecoord - ambcoords[i]);
      }
#endif
    }
    FPRINTF(fp,"|");
    if (ambigH == true) {
      ambcoords = Substring_ambcoords(substringH);
      n = Substring_nambcoords(substringH);
#ifdef PRINT_AMBIG_COORDS
      chroffset = Substring_chroffset(substringH);
      FPRINTF(fp,"%u",ambcoords[0] - chroffset + 1U);
      for (i = 1; i < n; i++) {
	FPRINTF(fp,",%u",ambcoords[i] - chroffset + 1U);
      }
#else
      substringM = (Substring_T) List_head(finalp);
      if (plusp == true) {
	splicecoord = Substring_alignend(substringM);
      } else {
	splicecoord = Substring_alignstart(substringM);
      }
      FPRINTF(fp,"%u",ambcoords[0] - splicecoord);
      for (i = 1; i < n; i++) {
	FPRINTF(fp,",%u",ambcoords[i] - splicecoord);
      }
#endif
    }
  }

  /* 12. TAGS: XC */
  if (circularp == true) {
    FPRINTF(fp,"\tXC:A:+");
  }

  /* 12. TAGS: XG */
  if (Stage3end_sarrayp(stage3end) == true) {
    FPRINTF(fp,"\tXG:Z:A");
  }

  FPRINTF(fp,"\n");
  return;
}


static void
halfdonor_dinucleotide (char *donor1, char *donor2, Substring_T donor, int sensedir) {
  char *genomic;
  int substring_start, substring_end;

  genomic = Substring_genomic_refdiff(donor);
  if (sensedir == SENSE_FORWARD) {
    substring_end = Substring_queryend(donor);
    *donor1 = toupper(genomic[substring_end]);
    *donor2 = toupper(genomic[substring_end+1]);
  } else {
    substring_start = Substring_querystart(donor);
    *donor2 = toupper(complCode[(int) genomic[substring_start-2]]);
    *donor1 = toupper(complCode[(int) genomic[substring_start-1]]);
  }

  return;
}

static void
halfacceptor_dinucleotide (char *acceptor2, char *acceptor1, Substring_T acceptor, int sensedir) {
  char *genomic;
  int substring_start, substring_end;

  genomic = Substring_genomic_refdiff(acceptor);
  if (sensedir == SENSE_FORWARD) {
    substring_start = Substring_querystart(acceptor);
    *acceptor2 = toupper(genomic[substring_start-2]);
    *acceptor1 = toupper(genomic[substring_start-1]);
  } else {
    substring_end = Substring_queryend(acceptor);
    *acceptor1 = toupper(complCode[(int) genomic[substring_end]]);
    *acceptor2 = toupper(complCode[(int) genomic[substring_end+1]]);
  }

  return;
}



static void
print_halfdonor (Filestring_T fp, char *abbrev, Substring_T donor, Stage3end_T this, Stage3end_T mate,
		 char *acc1, char *acc2, int pathnum, int npaths, int absmq_score, int first_absmq, int second_absmq, int mapq_score,
		 Univ_IIT_T chromosome_iit, Shortread_T queryseq, int pairedlength,
		 Chrpos_T concordant_chrpos, Chrpos_T donor_chrpos, Chrpos_T acceptor_chrpos, Chrpos_T mate_chrpos,
		 int hardclip_low, int hardclip_high, Resulttype_T resulttype, bool first_read_p,
		 bool artificial_mate_p, int npaths_mate,
		 int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p,
		 bool use_hardclip_p, bool print_xt_p, int donor_sensedir, char donor_strand, char acceptor_strand,
		 char *donor_chr, char *acceptor_chr, char donor1, char donor2, char acceptor2, char acceptor1,
		 double donor_prob, double acceptor_prob, bool circularp) {
  unsigned int flag = 0U;
  int nmismatches_refdiff = 0, nmismatches_bothdiff = 0, querylength;
  bool sensep;
  char *genomicfwd_refdiff, *genomicfwd_bothdiff, *genomicdir_refdiff, *genomicdir_bothdiff;
  int substring_start, substring_length;
  int transloc_hardclip_low, transloc_hardclip_high;
  bool plusp, printp;
  bool start_ambig, end_ambig;
  int n, i;
  Univcoord_T *start_ambcoords, *end_ambcoords, splicecoord;
#ifdef PRINT_AMBIG_COORDS
  Univcoord_T chroffset;
#endif


  querylength = Shortread_fulllength(queryseq);
  plusp = Substring_plusp(donor);

  /* 1. QNAME */
  if (acc2 == NULL) {
    FPRINTF(fp,"%s",acc1);
  } else {
    FPRINTF(fp,"%s,%s",acc1,acc2);
  }

  /* 2. FLAG */
  flag = SAM_compute_flag(plusp,mate,resulttype,first_read_p,
			  pathnum,npaths,artificial_mate_p,npaths_mate,absmq_score,first_absmq,
			  invertp,invert_mate_p);
  FPRINTF(fp,"\t%u",flag);

  /* 3. RNAME: chr */
  /* 4. POS: chrpos */
  print_chromosomal_pos(fp,Substring_chrnum(donor),donor_chrpos,Substring_chrlength(donor),chromosome_iit);
  

  /* 5. MAPQ: Mapping quality */
  FPRINTF(fp,"\t%d",mapq_score);

  /* 6. CIGAR */
  FPRINTF(fp,"\t");
  if (Stage3end_sensedir(this) == SENSE_ANTI) {
    sensep = false;
  } else {
    sensep = true;
  }

  if (use_hardclip_p == true) {
    if (sensep == true) {
      if (plusp == true) {
	transloc_hardclip_low = 0;
	transloc_hardclip_high = querylength - Substring_queryend(donor);
      } else {
	transloc_hardclip_high = 0;
	transloc_hardclip_low = querylength - Substring_queryend(donor);
      }

    } else {
      if (plusp == true) {
	transloc_hardclip_high = 0;
	transloc_hardclip_low = Substring_querystart(donor);
      } else {
	transloc_hardclip_low = 0;
	transloc_hardclip_high = Substring_querystart(donor);
      }
    }

    if (transloc_hardclip_low > hardclip_low) {
      hardclip_low = transloc_hardclip_low;
    }
    if (transloc_hardclip_high > hardclip_high) {
      hardclip_high = transloc_hardclip_high;
    }
  }


  if (sensep == true) {
    assert(Substring_chimera_pos(donor) == Substring_queryend(donor));
    if (plusp == true) {
      /* sensep true, plusp true */
      /* FPRINTF(fp,"donor sensep true, plusp true\n"); */
      if (hide_soft_clips_p == true) {
	print_cigar(fp,/*type*/'M',
		    Substring_querystart(donor) + 
		    Substring_match_length(donor),
		    /*querypos*/0,querylength,hardclip_low,hardclip_high,
		    /*plusp*/true,/*lastp*/false,/*trimlength*/0);
	print_cigar(fp,/*type*/'E',querylength - Substring_queryend(donor),
		    /*querypos*/Substring_queryend(donor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true,
		    /*trimlength*/Substring_trim_right(donor));

      } else {
	print_cigar(fp,/*type*/'S',Substring_querystart(donor),
		    /*querypos*/0,querylength,hardclip_low,hardclip_high,
		    /*plusp*/true,/*lastp*/false,/*trimlength*/0);
	print_cigar(fp,/*type*/'M',Substring_match_length(donor),
		    /*querypos*/Substring_querystart(donor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false,
		    /*trimlength*/0);
	print_cigar(fp,/*type*/'E',querylength - Substring_queryend(donor),
		    /*querypos*/Substring_queryend(donor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true,
		    /*trimlength*/Substring_trim_right(donor));
      }

    } else {
      /* sensep true, plusp false */
      /* FPRINTF(fp,"donor sensep false, plusp false\n"); */
      if (hide_soft_clips_p == true) {
	print_cigar(fp,/*type*/'E',querylength - Substring_queryend(donor),
		    /*querypos*/querylength,querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false,
		    /*trimlength*/Substring_trim_right(donor));
	print_cigar(fp,/*type*/'M',
		    Substring_match_length(donor) +
		    Substring_querystart(donor),
		    /*querypos*/Substring_queryend(donor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true,
		    /*trimlength*/0);
      } else {
	print_cigar(fp,/*type*/'E',querylength - Substring_queryend(donor),
		    /*querypos*/querylength,querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false,
		    /*trimlength*/Substring_trim_right(donor));
	print_cigar(fp,/*type*/'M',Substring_match_length(donor),
		    /*querypos*/Substring_queryend(donor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false,
		    /*trimlength*/0);
	print_cigar(fp,/*type*/'S',Substring_querystart(donor),
		    /*querypos*/Substring_querystart(donor),querylength,hardclip_low,hardclip_high,
		    /*plusp*/false,/*lastp*/true,/*trimlength*/0);
      }
    }

  } else {
    assert(Substring_chimera_pos(donor) == Substring_querystart(donor));
    if (plusp == true) {
      /* sensep false, plusp true */
      /* FPRINTF(fp,"donor sensep false, plusp true\n"); */
      if (hide_soft_clips_p == true) {
	print_cigar(fp,/*type*/'E',Substring_querystart(donor),
		    /*querypos*/0,querylength,hardclip_low,hardclip_high,
		    /*plusp*/true,/*lastp*/false,/*trimlength*/Substring_trim_left(donor));
	print_cigar(fp,/*type*/'M',Substring_match_length(donor) + (querylength - Substring_queryend(donor)),
		    /*querypos*/Substring_querystart(donor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true,
		    /*trimlength*/0);
      } else {
	print_cigar(fp,/*type*/'E',Substring_querystart(donor),
		    /*querypos*/0,querylength,hardclip_low,hardclip_high,
		    /*plusp*/true,/*lastp*/false,/*trimlength*/Substring_trim_left(donor));
	print_cigar(fp,/*type*/'M',Substring_match_length(donor),
		    /*querypos*/Substring_querystart(donor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false,
		    /*trimlength*/0);
	print_cigar(fp,/*type*/'S',querylength - Substring_queryend(donor),
		    /*querypos*/Substring_queryend(donor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true,
		    /*trimlength*/0);
      }

    } else {
      /* sensep false, plusp false */
      /* FPRINTF(fp,"donor sensep true, plusp false\n"); */
      if (hide_soft_clips_p == true) {
	print_cigar(fp,/*type*/'M',(querylength - Substring_queryend(donor)) + Substring_match_length(donor),
		    /*querypos*/querylength,querylength,hardclip_low,hardclip_high,
		    /*plusp*/false,/*lastp*/false,/*trimlength*/0);
	print_cigar(fp,/*type*/'E',Substring_querystart(donor),
		    /*querypos*/Substring_querystart(donor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true,
		    /*trimlength*/Substring_trim_left(donor));

      } else {
	print_cigar(fp,/*type*/'S',querylength - Substring_queryend(donor),
		    /*querypos*/querylength,querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false,
		    /*trimlength*/0);
	print_cigar(fp,/*type*/'M',Substring_match_length(donor),
		    /*querypos*/Substring_queryend(donor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false,
		    /*trimlength*/0);
	print_cigar(fp,/*type*/'E',Substring_querystart(donor),
		    /*querypos*/Substring_querystart(donor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true,
		    /*trimlength*/Substring_trim_left(donor));
      }
    }
  }

  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  /* For anchor_chrnum, previously used Stage3end_chrnum(this), but this is 0 */
  print_mate_chromosomal_pos(fp,Stage3end_chrnum(mate),Stage3end_effective_chrnum(mate),
			     mate_chrpos,Stage3end_chrlength(mate),
			     /*anchor_chrnum*/Substring_chrnum(donor),donor_chrpos,chromosome_iit);


  /* 9. ISIZE: Insert size */
  if (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_TRANSLOC || resulttype == CONCORDANT_MULT) {
    if (plusp == invertp) {
      FPRINTF(fp,"\t%d",-pairedlength);
    } else {
      FPRINTF(fp,"\t%d",pairedlength);
    }
  } else if (mate_chrpos == 0) {
    FPRINTF(fp,"\t%d",pairedlength);
#if 0
  } else if (concordant_chrpos < mate_chrpos) {
    FPRINTF(fp,"\t%d",pairedlength);
  } else if (concordant_chrpos > mate_chrpos) {
    FPRINTF(fp,"\t%d",-pairedlength);
#endif
  } else if (first_read_p == true) {
    FPRINTF(fp,"\t%d",pairedlength);
  } else {
    FPRINTF(fp,"\t%d",-pairedlength);
  }


  /* 10. SEQ: queryseq and 11. QUAL: quality scores */
  /* Queryseq has already been inverted, so just measure plusp relative to its current state */
  if (plusp == true) {
    Shortread_print_chopped_sam(fp,queryseq,hardclip_low,hardclip_high);
    FPRINTF(fp,"\t");
    Shortread_print_quality(fp,queryseq,hardclip_low,hardclip_high,
			    quality_shift,/*show_chopped_p*/false);
  } else {
    Shortread_print_chopped_revcomp_sam(fp,queryseq,hardclip_low,hardclip_high);
    FPRINTF(fp,"\t");
    Shortread_print_quality_revcomp(fp,queryseq,hardclip_low,hardclip_high,
				    quality_shift,/*show_chopped_p*/false);
  }


  /* 12. TAGS: RG */
  if (sam_read_group_id != NULL) {
    FPRINTF(fp,"\tRG:Z:%s",sam_read_group_id);
  }

  /* 12. TAGS: XH and XI */
  if (hardclip_low > 0 || hardclip_high > 0) {
    FPRINTF(fp,"\tXH:Z:");
    if (plusp == true) {
      Shortread_print_chopped_end(fp,queryseq,hardclip_low,hardclip_high);
    } else {
      Shortread_print_chopped_end_revcomp(fp,queryseq,hardclip_low,hardclip_high);
    }

    if (Shortread_quality_string(queryseq) != NULL) {
      FPRINTF(fp,"\tXI:Z:");
      if (plusp == true) {
	Shortread_print_chopped_end_quality(fp,queryseq,hardclip_low,hardclip_high);
      } else {
	Shortread_print_chopped_end_quality_reverse(fp,queryseq,hardclip_low,hardclip_high);
      }
    }
  }

  /* 12. TAGS: XB */
  Shortread_print_barcode(fp,queryseq);

  /* 12. TAGS: XP.  Logically should be last in reconstructing a read. */
  Shortread_print_chop(fp,queryseq,invertp);

  /* 12. TAGS: MD */
  FPRINTF(fp,"\tMD:Z:");
  printp = false;

  if (hide_soft_clips_p == true) {
    substring_start = Substring_querystart_orig(donor);
    substring_length = Substring_match_length_orig(donor);
  } else {
    substring_start = Substring_querystart(donor);
    substring_length = Substring_match_length(donor);
  }

  if (use_hardclip_p == false) {
    genomicdir_refdiff = Substring_genomic_refdiff(donor);
    genomicdir_bothdiff = Substring_genomic_bothdiff(donor);
    if (plusp == true) {
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,/*matchlength*/0,
		      &(genomicdir_refdiff[substring_start]),&(genomicdir_bothdiff[substring_start]),
		      substring_length,/*querypos*/substring_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    } else if (genomicdir_bothdiff == genomicdir_refdiff) {
      genomicfwd_refdiff = (char *) MALLOCA((querylength+1) * sizeof(char));
      make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
		      fp,/*matchlength*/0,genomicfwd_refdiff,/*genomicfwd_bothdiff*/genomicfwd_refdiff,
		      substring_length,/*querypos*/substring_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      FREEA(genomicfwd_refdiff);
    } else {
      genomicfwd_refdiff = (char *) MALLOCA((querylength+1) * sizeof(char));
      genomicfwd_bothdiff = (char *) MALLOCA((querylength+1) * sizeof(char));
      make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
      make_complement_buffered(genomicfwd_bothdiff,&(genomicdir_bothdiff[substring_start]),substring_length);
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
		      fp,/*matchlength*/0,genomicfwd_refdiff,genomicfwd_bothdiff,
		      substring_length,/*querypos*/substring_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      FREEA(genomicfwd_bothdiff);
      FREEA(genomicfwd_refdiff);
    }

  } else if (sensep == true) {
    if (plusp == true) {
      genomicfwd_refdiff = Substring_genomic_refdiff(donor);
      genomicfwd_bothdiff = Substring_genomic_bothdiff(donor);
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,/*matchlength*/0,
		      &(genomicfwd_refdiff[substring_start]),&(genomicfwd_bothdiff[substring_start]),
		      substring_length,/*querypos*/substring_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    } else {
      genomicdir_refdiff = Substring_genomic_refdiff(donor);
      genomicdir_bothdiff = Substring_genomic_bothdiff(donor);
      if (genomicdir_bothdiff == genomicdir_refdiff) {
	genomicfwd_refdiff = (char *) MALLOCA((substring_length+1) * sizeof(char));
	make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
	print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
			fp,/*matchlength*/0,genomicfwd_refdiff,/*genomicfwd_bothdiff*/genomicfwd_refdiff,
			substring_length,/*querypos*/substring_start,querylength,
			hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
	FREEA(genomicfwd_refdiff);
      } else {
	genomicfwd_refdiff = (char *) MALLOCA((substring_length+1) * sizeof(char));
	genomicfwd_bothdiff = (char *) MALLOCA((substring_length+1) * sizeof(char));
	make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
	make_complement_buffered(genomicfwd_bothdiff,&(genomicdir_bothdiff[substring_start]),substring_length);
	print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
			fp,/*matchlength*/0,genomicfwd_refdiff,genomicfwd_bothdiff,
			substring_length,/*querypos*/substring_start,querylength,
			hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
	FREEA(genomicfwd_bothdiff);
	FREEA(genomicfwd_refdiff);
      }
    }

  } else {			/* sensep == false */
    if (plusp == true) {
      genomicfwd_refdiff = Substring_genomic_refdiff(donor);
      genomicfwd_bothdiff = Substring_genomic_bothdiff(donor);
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,/*matchlength*/0,
		      &(genomicfwd_refdiff[substring_start]),&(genomicfwd_bothdiff[substring_start]),
		      substring_length,/*querypos*/substring_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    } else {
      genomicdir_refdiff = Substring_genomic_refdiff(donor);
      genomicdir_bothdiff = Substring_genomic_refdiff(donor);
      if (genomicdir_bothdiff == genomicdir_refdiff) {
	genomicfwd_refdiff = (char *) MALLOCA((substring_length+1) * sizeof(char));
	make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
	print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
			fp,/*matchlength*/0,genomicfwd_refdiff,/*genomicfwd_bothdiff*/genomicfwd_refdiff,
			substring_length,/*querypos*/substring_start,querylength,
			hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
	FREEA(genomicfwd_refdiff);
      } else {
	genomicfwd_refdiff = (char *) MALLOCA((substring_length+1) * sizeof(char));
	genomicfwd_bothdiff = (char *) MALLOCA((substring_length+1) * sizeof(char));
	make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
	make_complement_buffered(genomicfwd_bothdiff,&(genomicdir_bothdiff[substring_start]),substring_length);
	print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
			fp,/*matchlength*/0,genomicfwd_refdiff,genomicfwd_bothdiff,
			substring_length,/*querypos*/substring_start,querylength,
			hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
	FREEA(genomicfwd_bothdiff);
	FREEA(genomicfwd_refdiff);
      }
    }
  }
  if (printp == false) {
    FPRINTF(fp,"0");
  }


  /* 12. TAGS: NH */
  /* 12. TAGS: HI */
  /* 12. TAGS: NM */
  FPRINTF(fp,"\tNH:i:%d\tHI:i:%d\tNM:i:%d",npaths,pathnum,nmismatches_refdiff);
  
  if (snps_iit) {
    /* 12. TAGS: XW and XV */
    FPRINTF(fp,"\tXW:i:%d",nmismatches_bothdiff);
    FPRINTF(fp,"\tXV:i:%d",nmismatches_refdiff - nmismatches_bothdiff);
  }

  /* 12. TAGS: SM */
  /* 12. TAGS: XQ */
  /* 12. TAGS: X2 */
  FPRINTF(fp,"\tSM:i:%d\tXQ:i:%d\tX2:i:%d",mapq_score,absmq_score,second_absmq);

  /* 12. TAGS: XO */
  FPRINTF(fp,"\tXO:Z:%s",abbrev);

  /* 12. TAGS: XS */
  assert(donor_sensedir != SENSE_NULL);
  FPRINTF(fp,"\tXS:A:%c",donor_strand);

  /* 12. TAGS: XA */
  if ((start_ambig = Stage3end_start_ambiguous_p(this)) == true ||
      (end_ambig = Stage3end_end_ambiguous_p(this)) == true) {
    FPRINTF(fp,"\tXA:Z:");

    if (plusp == true) {
      if ((n = Stage3end_start_nambcoords(this)) > 0) {
	assert(sensep == false);
	start_ambcoords = Stage3end_start_ambcoords(this);
	splicecoord = Substring_alignstart(donor);
#ifdef PRINT_AMBIG_COORDS
	chroffset = Substring_chroffset(donor);
	FPRINTF(fp,"%u",start_ambcoords[0] - chroffset + 1U);
	for (i = 1; i < n; i++) {
	  FPRINTF(fp,",%u",start_ambcoords[i] - chroffset + 1U);
	}
#else
	splicecoord = Substring_alignstart(donor);
	FPRINTF(fp,"%u",splicecoord - start_ambcoords[0]);
	for (i = 1; i < n; i++) {
	  FPRINTF(fp,",%u",splicecoord - start_ambcoords[i]);
	}
#endif
      }
      FPRINTF(fp,"|");
      if ((n = Stage3end_end_nambcoords(this)) > 0) {
	assert(sensep == true);
	end_ambcoords = Stage3end_end_ambcoords(this);
#ifdef PRINT_AMBIG_COORDS
	chroffset = Substring_chroffset(donor);
	FPRINTF(fp,"%u",end_ambcoords[0] - chroffset + 1U);
	for (i = 1; i < n; i++) {
	  FPRINTF(fp,",%u",end_ambcoords[i] - chroffset + 1U);
	}
#else
	splicecoord = Substring_alignend(donor);
	FPRINTF(fp,"%u",end_ambcoords[0] - splicecoord);
	for (i = 1; i < n; i++) {
	  FPRINTF(fp,",%u",end_ambcoords[i] - splicecoord);
	}
#endif
      }

    } else {
      if ((n = Stage3end_end_nambcoords(this)) > 0) {
	assert(sensep == true);
	end_ambcoords = Stage3end_end_ambcoords(this);
#ifdef PRINT_AMBIG_COORDS
	chroffset = Substring_chroffset(donor);
	FPRINTF(fp,"%u",end_ambcoords[0] - chroffset + 1U);
	for (i = 1; i < n; i++) {
	  FPRINTF(fp,",%u",end_ambcoords[i] - chroffset + 1U);
	}
#else
	splicecoord = Substring_alignend(donor);
	FPRINTF(fp,"%u",splicecoord - end_ambcoords[0]);
	for (i = 1; i < n; i++) {
	  FPRINTF(fp,",%u",splicecoord - end_ambcoords[i]);
	}
#endif
      }
      FPRINTF(fp,"|");
      if ((n = Stage3end_start_nambcoords(this)) > 0) {
	assert(sensep == false);
	start_ambcoords = Stage3end_start_ambcoords(this);
#ifdef PRINT_AMBIG_COORDS
	chroffset = Substring_chroffset(donor);
	FPRINTF(fp,"%u",start_ambcoords[0] - chroffset + 1U);
	for (i = 1; i < n; i++) {
	  FPRINTF(fp,",%u",start_ambcoords[i] - chroffset + 1U);
	}
#else
	splicecoord = Substring_alignstart(donor);
	FPRINTF(fp,"%u",start_ambcoords[0] - splicecoord);
	for (i = 1; i < n; i++) {
	  FPRINTF(fp,",%u",start_ambcoords[i] - splicecoord);
	}
#endif
      }
    }
  }

  /* 12. TAGS: XT */
  if (print_xt_p == true) {
    FPRINTF(fp,"\tXT:Z:%c%c-%c%c,%.2f,%.2f",donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob);
    FPRINTF(fp,",%c%s@%u..%c%s@%u",donor_strand,donor_chr,donor_chrpos,acceptor_strand,acceptor_chr,acceptor_chrpos);
  }

  /* 12. TAGS: XC */
  if (circularp == true) {
    FPRINTF(fp,"\tXC:A:+");
  }

  /* 12. TAGS: XG */
  if (Stage3end_sarrayp(this) == true) {
    FPRINTF(fp,"\tXG:Z:A");
  }

  FPRINTF(fp,"\n");
  return;
}



static void
print_halfacceptor (Filestring_T fp, char *abbrev, Substring_T acceptor, Stage3end_T this, Stage3end_T mate,
		    char *acc1, char *acc2, int pathnum, int npaths, int absmq_score, int first_absmq, int second_absmq, int mapq_score,
		    Univ_IIT_T chromosome_iit, Shortread_T queryseq, int pairedlength,
		    Chrpos_T concordant_chrpos, Chrpos_T donor_chrpos, Chrpos_T acceptor_chrpos, Chrpos_T mate_chrpos,
		    int hardclip_low, int hardclip_high, Resulttype_T resulttype, bool first_read_p, bool artificial_mate_p, int npaths_mate,
		    int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p,
		    bool use_hardclip_p, bool print_xt_p, int acceptor_sensedir, char donor_strand, char acceptor_strand,
		    char *donor_chr, char *acceptor_chr, char donor1, char donor2, char acceptor2, char acceptor1,
		    double donor_prob, double acceptor_prob, bool circularp) {
  unsigned int flag = 0U;
  int nmismatches_refdiff = 0, nmismatches_bothdiff = 0, querylength;
  bool sensep;
  char *genomicfwd_refdiff, *genomicfwd_bothdiff, *genomicdir_refdiff, *genomicdir_bothdiff;
  int substring_start, substring_length;
  int transloc_hardclip_low, transloc_hardclip_high;
  bool plusp, printp;
  bool start_ambig, end_ambig;
  int n, i;
  Univcoord_T *start_ambcoords, *end_ambcoords, splicecoord;
#ifdef PRINT_AMBIG_COORDS
  Univcoord_T chroffset;
#endif


  querylength = Shortread_fulllength(queryseq);
  plusp = Substring_plusp(acceptor);

  /* 1. QNAME */
  if (acc2 == NULL) {
    FPRINTF(fp,"%s",acc1);
  } else {
    FPRINTF(fp,"%s,%s",acc1,acc2);
  }

  /* 2. FLAG */
  flag = SAM_compute_flag(plusp,mate,resulttype,first_read_p,
			  pathnum,npaths,artificial_mate_p,npaths_mate,absmq_score,first_absmq,
			  invertp,invert_mate_p);
  FPRINTF(fp,"\t%u",flag);

  /* 3. RNAME: chr */
  /* 4. POS: chrpos */
  print_chromosomal_pos(fp,Substring_chrnum(acceptor),acceptor_chrpos,Substring_chrlength(acceptor),chromosome_iit);


  /* 5. MAPQ: Mapping quality */
  FPRINTF(fp,"\t%d",mapq_score);

  /* 6. CIGAR */
  FPRINTF(fp,"\t");
  if (Stage3end_sensedir(this) == SENSE_ANTI) {
    sensep = false;
  } else {
    sensep = true;
  }

  if (use_hardclip_p == true) {
    if (sensep == true) {
      if (plusp == true) {
	transloc_hardclip_high = 0;
	transloc_hardclip_low = Substring_querystart(acceptor);
      } else {
	transloc_hardclip_low = 0;
	transloc_hardclip_high = Substring_querystart(acceptor);
      }

    } else {
      if (plusp == true) {
	transloc_hardclip_low = 0;
	transloc_hardclip_high = querylength - Substring_queryend(acceptor);
      } else {
	transloc_hardclip_high = 0;
	transloc_hardclip_low = querylength - Substring_queryend(acceptor);
      }
    }

    if (transloc_hardclip_low > hardclip_low) {
      hardclip_low = transloc_hardclip_low;
    }
    if (transloc_hardclip_high > hardclip_high) {
      hardclip_high = transloc_hardclip_high;
    }
  }


  if (sensep == true) {
    assert(Substring_chimera_pos(acceptor) == Substring_querystart(acceptor));
    if (plusp == true) {
      /* sensep true, plusp true */
      /* FPRINTF(fp,"acceptor sensep true, plusp true\n"); */
      if (hide_soft_clips_p == true) {
	print_cigar(fp,/*type*/'M',Substring_querystart(acceptor) + Substring_match_length(acceptor),
		    /*querypos*/0,querylength,hardclip_low,hardclip_high,
		    /*plusp*/true,/*lastp*/false,/*trimlength*/0);
	print_cigar(fp,/*type*/'E',querylength - Substring_queryend(acceptor),
		    /*querypos*/Substring_queryend(acceptor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true,
		    /*trimlength*/Substring_trim_right(acceptor));
      } else {
	print_cigar(fp,/*type*/'S',Substring_querystart(acceptor),
		    /*querypos*/0,querylength,hardclip_low,hardclip_high,
		    /*plusp*/true,/*lastp*/false,/*trimlength*/0);
	print_cigar(fp,/*type*/'M',Substring_match_length(acceptor),
		    /*querypos*/Substring_querystart(acceptor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false,/*trimlength*/0);
	print_cigar(fp,/*type*/'E',querylength - Substring_queryend(acceptor),
		    /*querypos*/Substring_queryend(acceptor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true,
		    /*trimlength*/Substring_trim_right(acceptor));
      }

    } else {
      /* sensep true, plusp false */
      /* FPRINTF(fp,"acceptor sensep true, plusp false\n"); */
      if (hide_soft_clips_p == true) {
	print_cigar(fp,/*type*/'E',querylength - Substring_queryend(acceptor),
		    /*querypos*/querylength,querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false,
		    /*trimlength*/Substring_trim_right(acceptor));
	print_cigar(fp,/*type*/'M',Substring_match_length(acceptor) + Substring_querystart(acceptor),
		    /*querypos*/Substring_queryend(acceptor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true,
		    /*trimlength*/0);
      } else {
	print_cigar(fp,/*type*/'E',querylength - Substring_queryend(acceptor),
		    /*querypos*/querylength,querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false,
		    /*trimlength*/Substring_trim_right(acceptor));
	print_cigar(fp,/*type*/'M',Substring_match_length(acceptor),
		    /*querypos*/Substring_queryend(acceptor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false,
		    /*trimlength*/0);
	print_cigar(fp,/*type*/'S',Substring_querystart(acceptor),
		    /*querypos*/Substring_querystart(acceptor),querylength,hardclip_low,hardclip_high,
		    /*plusp*/false,/*lastp*/true,/*trimlength*/0);
      }
    }

  } else {
    /* sensep false, plusp true */
    assert(Substring_chimera_pos(acceptor) == Substring_queryend(acceptor));
    if (plusp == true) {
      /* FPRINTF(fp,"acceptor sensep false, plusp true\n"); */
      if (hide_soft_clips_p == true) {
	print_cigar(fp,/*type*/'E',Substring_querystart(acceptor),
		    /*querypos*/0,querylength,hardclip_low,hardclip_high,
		    /*plusp*/true,/*lastp*/false,/*trimlength*/Substring_trim_left(acceptor));
	print_cigar(fp,/*type*/'M',Substring_match_length(acceptor) + (querylength - Substring_queryend(acceptor)),
		    /*querypos*/Substring_querystart(acceptor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true,
		    /*trimlength*/0);
      } else {
	print_cigar(fp,/*type*/'E',Substring_querystart(acceptor),
		    /*querypos*/0,querylength,hardclip_low,hardclip_high,
		    /*plusp*/true,/*lastp*/false,/*trimlength*/Substring_trim_left(acceptor));
	print_cigar(fp,/*type*/'M',Substring_match_length(acceptor),
		    /*querypos*/Substring_querystart(acceptor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false,/*trimlength*/0);
	print_cigar(fp,/*type*/'S',querylength - Substring_queryend(acceptor),
		    /*querypos*/Substring_queryend(acceptor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true,
		    /*trimlength*/Substring_trim_right(acceptor));
      }

    } else {
      /* sensep false, plusp false */
      /* FPRINTF(fp,"acceptor sensep false, plusp false\n"); */
      if (hide_soft_clips_p == true) {
	print_cigar(fp,/*type*/'M',(querylength - Substring_queryend(acceptor)) + Substring_match_length(acceptor),
		    /*querypos*/querylength,querylength,hardclip_low,hardclip_high,
		    /*plusp*/false,/*lastp*/false,/*trimlength*/0);
	print_cigar(fp,/*type*/'E',Substring_querystart(acceptor),
		    /*querypos*/Substring_querystart(acceptor),querylength,hardclip_low,hardclip_high,
		    /*plusp*/false,/*lastp*/true,/*trimlength*/Substring_trim_left(acceptor));
      } else {
	print_cigar(fp,/*type*/'S',querylength - Substring_queryend(acceptor),
		    /*querypos*/querylength,querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false,
		    /*trimlength*/0);
	print_cigar(fp,/*type*/'M',Substring_match_length(acceptor),
		    /*querypos*/Substring_queryend(acceptor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false,
		    /*trimlength*/0);
	print_cigar(fp,/*type*/'E',Substring_querystart(acceptor),
		    /*querypos*/Substring_querystart(acceptor),querylength,hardclip_low,hardclip_high,
		    /*plusp*/false,/*lastp*/true,/*trimlength*/Substring_trim_left(acceptor));
      }
    }
  }


  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  /* For anchor_chrnum, previously used Stage3end_chrnum(this), but this is 0 */
  print_mate_chromosomal_pos(fp,Stage3end_chrnum(mate),Stage3end_effective_chrnum(mate),
			     mate_chrpos,Stage3end_chrlength(mate),
			     /*anchor_chrnum*/Substring_chrnum(acceptor),acceptor_chrpos,chromosome_iit);


  /* 9. ISIZE: Insert size */
  if (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_TRANSLOC || resulttype == CONCORDANT_MULT) {
    if (plusp == invertp) {
      FPRINTF(fp,"\t%d",-pairedlength);
    } else {
      FPRINTF(fp,"\t%d",pairedlength);
    }
  } else if (mate_chrpos == 0) {
    FPRINTF(fp,"\t%d",pairedlength);
#if 0
  } else if (concordant_chrpos < mate_chrpos) {
    FPRINTF(fp,"\t%d",pairedlength);
  } else if (concordant_chrpos > mate_chrpos) {
    FPRINTF(fp,"\t%d",-pairedlength);
#endif
  } else if (first_read_p == true) {
    FPRINTF(fp,"\t%d",pairedlength);
  } else {
    FPRINTF(fp,"\t%d",-pairedlength);
  }


  /* 10. SEQ: queryseq and 11. QUAL: quality scores */
  /* Queryseq has already been inverted, so just measure plusp relative to its current state */
  if (plusp == true) {
    Shortread_print_chopped_sam(fp,queryseq,hardclip_low,hardclip_high);
    FPRINTF(fp,"\t");
    Shortread_print_quality(fp,queryseq,hardclip_low,hardclip_high,
			    quality_shift,/*show_chopped_p*/false);
  } else {
    Shortread_print_chopped_revcomp_sam(fp,queryseq,hardclip_low,hardclip_high);
    FPRINTF(fp,"\t");
    Shortread_print_quality_revcomp(fp,queryseq,hardclip_low,hardclip_high,
				    quality_shift,/*show_chopped_p*/false);
  }


  /* 12. TAGS: RG */
  if (sam_read_group_id != NULL) {
    FPRINTF(fp,"\tRG:Z:%s",sam_read_group_id);
  }

  /* 12. TAGS: XH and XI */
  if (hardclip_low > 0 || hardclip_high > 0) {
    FPRINTF(fp,"\tXH:Z:");
    if (plusp == true) {
      Shortread_print_chopped_end(fp,queryseq,hardclip_low,hardclip_high);
    } else {
      Shortread_print_chopped_end_revcomp(fp,queryseq,hardclip_low,hardclip_high);
    }

    if (Shortread_quality_string(queryseq) != NULL) {
      FPRINTF(fp,"\tXI:Z:");
      if (plusp == true) {
	Shortread_print_chopped_end_quality(fp,queryseq,hardclip_low,hardclip_high);
      } else {
	Shortread_print_chopped_end_quality_reverse(fp,queryseq,hardclip_low,hardclip_high);
      }
    }
  }

  /* 12. TAGS: XB */
  Shortread_print_barcode(fp,queryseq);

  /* 12. TAGS: XP.  Logically should be last in reconstructing a read. */
  Shortread_print_chop(fp,queryseq,invertp);

  /* 12. TAGS: MD */
  FPRINTF(fp,"\tMD:Z:");
  printp = false;

  if (hide_soft_clips_p == true) {
    substring_start = Substring_querystart_orig(acceptor);
    substring_length = Substring_match_length_orig(acceptor);
  } else {
    substring_start = Substring_querystart(acceptor);
    substring_length = Substring_match_length(acceptor);
  }

  if (use_hardclip_p == false) {
    genomicdir_refdiff = Substring_genomic_refdiff(acceptor);
    genomicdir_bothdiff = Substring_genomic_bothdiff(acceptor);
    if (plusp == true) {
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,/*matchlength*/0,
		      &(genomicdir_refdiff[substring_start]),&(genomicdir_bothdiff[substring_start]),
		      substring_length,/*querypos*/substring_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    } else if (genomicdir_bothdiff == genomicdir_refdiff) {
      genomicfwd_refdiff = (char *) MALLOCA((querylength+1) * sizeof(char));
      make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
		      fp,/*matchlength*/0,genomicfwd_refdiff,/*genomicfwd_bothdiff*/genomicfwd_refdiff,
		      substring_length,/*querypos*/substring_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      FREEA(genomicfwd_refdiff);
    } else {
      genomicfwd_refdiff = (char *) MALLOCA((querylength+1) * sizeof(char));
      genomicfwd_bothdiff = (char *) MALLOCA((querylength+1) * sizeof(char));
      make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
      make_complement_buffered(genomicfwd_bothdiff,&(genomicdir_bothdiff[substring_start]),substring_length);
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
		      fp,/*matchlength*/0,genomicfwd_refdiff,genomicfwd_bothdiff,
		      substring_length,/*querypos*/substring_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      FREEA(genomicfwd_bothdiff);
      FREEA(genomicfwd_refdiff);
    }

  } else if (sensep == false) {
    if (plusp == true) {
      genomicfwd_refdiff = Substring_genomic_refdiff(acceptor);
      genomicfwd_bothdiff = Substring_genomic_bothdiff(acceptor);
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,/*matchlength*/0,
		      &(genomicfwd_refdiff[substring_start]),&(genomicfwd_bothdiff[substring_start]),
		      substring_length,/*querypos*/substring_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    } else {
      genomicdir_refdiff = Substring_genomic_refdiff(acceptor);
      genomicdir_bothdiff = Substring_genomic_bothdiff(acceptor);
      if (genomicdir_bothdiff == genomicdir_refdiff) {
	genomicfwd_refdiff = (char *) MALLOCA((substring_length+1) * sizeof(char));
	make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
	print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
			fp,/*matchlength*/0,genomicfwd_refdiff,/*genomicfwd_bothdiff*/genomicfwd_refdiff,
			substring_length,/*querypos*/substring_start,querylength,
			hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
	FREEA(genomicfwd_refdiff);
      } else {
	genomicfwd_refdiff = (char *) MALLOCA((substring_length+1) * sizeof(char));
	genomicfwd_bothdiff = (char *) MALLOCA((substring_length+1) * sizeof(char));
	make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
	make_complement_buffered(genomicfwd_bothdiff,&(genomicdir_bothdiff[substring_start]),substring_length);
	print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
			fp,/*matchlength*/0,genomicfwd_refdiff,genomicfwd_bothdiff,
			substring_length,/*querypos*/substring_start,querylength,
			hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
	FREEA(genomicfwd_bothdiff);
	FREEA(genomicfwd_refdiff);
      }

    }

  } else {			/* sensep true */
    if (plusp == true) {
      genomicfwd_refdiff = Substring_genomic_refdiff(acceptor);
      genomicfwd_bothdiff = Substring_genomic_bothdiff(acceptor);
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,/*matchlength*/0,
		      &(genomicfwd_refdiff[substring_start]),&(genomicfwd_bothdiff[substring_start]),
		      substring_length,/*querypos*/substring_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    } else {
      genomicdir_refdiff = Substring_genomic_refdiff(acceptor);
      genomicdir_bothdiff = Substring_genomic_bothdiff(acceptor);
      if (genomicdir_bothdiff == genomicdir_refdiff) {
	genomicfwd_refdiff = (char *) MALLOCA((substring_length+1) * sizeof(char));
	make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
	print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
			fp,/*matchlength*/0,genomicfwd_refdiff,/*genomicfwd_bothdiff*/genomicfwd_refdiff,
			substring_length,/*querypos*/substring_start,querylength,
			hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
	FREEA(genomicfwd_refdiff);
      } else {
	genomicfwd_refdiff = (char *) MALLOCA((substring_length+1) * sizeof(char));
	genomicfwd_bothdiff = (char *) MALLOCA((substring_length+1) * sizeof(char));
	make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
	make_complement_buffered(genomicfwd_bothdiff,&(genomicdir_bothdiff[substring_start]),substring_length);
	print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
			fp,/*matchlength*/0,genomicfwd_refdiff,genomicfwd_bothdiff,
			substring_length,/*querypos*/substring_start,querylength,
			hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
	FREEA(genomicfwd_bothdiff);
	FREEA(genomicfwd_refdiff);
      }
    }
  }
  if (printp == false) {
    FPRINTF(fp,"0");
  }


  /* 12. TAGS: NH */
  /* 12. TAGS: HI */
  /* 12. TAGS: NM */
  FPRINTF(fp,"\tNH:i:%d\tHI:i:%d\tNM:i:%d",npaths,pathnum,nmismatches_refdiff);
  
  if (snps_iit) {
    /* 12. TAGS: XW and XV */
    FPRINTF(fp,"\tXW:i:%d",nmismatches_bothdiff);
    FPRINTF(fp,"\tXV:i:%d",nmismatches_refdiff - nmismatches_bothdiff);
  }

  /* 12. TAGS: SM */
  /* 12. TAGS: XQ */
  /* 12. TAGS: X2 */
  FPRINTF(fp,"\tSM:i:%d\tXQ:i:%d\tX2:i:%d",mapq_score,absmq_score,second_absmq);

  /* 12. TAGS: XO */
  FPRINTF(fp,"\tXO:Z:%s",abbrev);

  /* 12. TAGS: XS */
  assert(acceptor_sensedir != SENSE_NULL);
  FPRINTF(fp,"\tXS:A:%c",acceptor_strand);

  /* 12. TAGS: XA */
  if ((start_ambig = Stage3end_start_ambiguous_p(this)) == true ||
      (end_ambig = Stage3end_end_ambiguous_p(this)) == true) {
    FPRINTF(fp,"\tXA:Z:");

    if (plusp == true) {
      if ((n = Stage3end_start_nambcoords(this)) > 0) {
	assert(sensep == true);
	start_ambcoords = Stage3end_start_ambcoords(this);
#ifdef PRINT_AMBIG_COORDS
	chroffset = Substring_chroffset(acceptor);
	FPRINTF(fp,"%u",start_ambcoords[0] - chroffset + 1U);
	for (i = 1; i < n; i++) {
	  FPRINTF(fp,",%u",start_ambcoords[i] - chroffset + 1U);
	}
#else
	splicecoord = Substring_alignstart(acceptor);
	FPRINTF(fp,"%u",splicecoord - start_ambcoords[0]);
	for (i = 1; i < n; i++) {
	  FPRINTF(fp,",%u",splicecoord - start_ambcoords[i]);
	}
#endif
      }
      FPRINTF(fp,"|");
      if ((n = Stage3end_end_nambcoords(this)) > 0) {
	assert(sensep == false);
	end_ambcoords = Stage3end_end_ambcoords(this);
#ifdef PRINT_AMBIG_COORDS
	chroffset = Substring_chroffset(acceptor);
	FPRINTF(fp,"%u",end_ambcoords[0] - chroffset + 1U);
	for (i = 1; i < n; i++) {
	  FPRINTF(fp,",%u",end_ambcoords[i] - chroffset + 1U);
	}
#else
	splicecoord = Substring_alignend(acceptor);
	FPRINTF(fp,"%u",end_ambcoords[0] - splicecoord);
	for (i = 1; i < n; i++) {
	  FPRINTF(fp,",%u",end_ambcoords[i] - splicecoord);
	}
#endif
      }

    } else {
      if ((n = Stage3end_end_nambcoords(this)) > 0) {
	assert(sensep == false);
	end_ambcoords = Stage3end_end_ambcoords(this);
#ifdef PRINT_AMBIG_COORDS
	chroffset = Substring_chroffset(acceptor);
	FPRINTF(fp,"%u",end_ambcoords[0] - chroffset + 1U);
	for (i = 1; i < n; i++) {
	  FPRINTF(fp,",%u",end_ambcoords[i] - chroffset + 1U);
	}
#else
	splicecoord = Substring_alignend(acceptor);
	FPRINTF(fp,"%u",splicecoord - end_ambcoords[0]);
	for (i = 1; i < n; i++) {
	  FPRINTF(fp,",%u",splicecoord - end_ambcoords[i]);
	}
#endif
      }
      FPRINTF(fp,"|");
      if ((n = Stage3end_start_nambcoords(this)) > 0) {
	assert(sensep == true);
	start_ambcoords = Stage3end_start_ambcoords(this);
#ifdef PRINT_AMBIG_COORDS
	chroffset = Substring_chroffset(acceptor);
	FPRINTF(fp,"%u",start_ambcoords[0] - chroffset + 1U);
	for (i = 1; i < n; i++) {
	  FPRINTF(fp,",%u",start_ambcoords[i] - chroffset + 1U);
	}
#else
	splicecoord = Substring_alignstart(acceptor);
	FPRINTF(fp,"%u",start_ambcoords[0] - splicecoord);
	for (i = 1; i < n; i++) {
	  FPRINTF(fp,",%u",start_ambcoords[i] - splicecoord);
	}
#endif
      }
    }
  }

  /* 12. TAGS: XT */
  if (print_xt_p == true) {
    FPRINTF(fp,"\tXT:Z:%c%c-%c%c,%.2f,%.2f",donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob);
    FPRINTF(fp,",%c%s@%u..%c%s@%u",donor_strand,donor_chr,donor_chrpos,acceptor_strand,acceptor_chr,acceptor_chrpos);
  }


  /* 12. TAGS: XC */
  if (circularp == true) {
    FPRINTF(fp,"\tXC:A:+");
  }

  /* 12. TAGS: XG */
  if (Stage3end_sarrayp(this) == true) {
    FPRINTF(fp,"\tXG:Z:A");
  }

  FPRINTF(fp,"\n");
  return;
}


/* Distant splicing, including scramble, inversion, translocation */
static void
print_exon_exon (Filestring_T fp, char *abbrev, Stage3end_T this, Stage3end_T mate,
		 char *acc1, char *acc2, int pathnum, int npaths,
		 int absmq_score, int first_absmq, int second_absmq, int mapq_score,
		 Univ_IIT_T chromosome_iit, Shortread_T queryseq, int pairedlength,
		 Chrpos_T mate_chrpos, int hardclip_low, int hardclip_high,
		 Resulttype_T resulttype, bool first_read_p, bool artificial_mate_p, int npaths_mate,
		 int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p) {
  Chrpos_T donor_chrpos, acceptor_chrpos, concordant_chrpos;
  Substring_T donor, acceptor;
  char *donor_chr, *acceptor_chr;
  char donor1, donor2, acceptor2, acceptor1;
  double donor_prob, acceptor_prob;
  int circularpos, querylength;
  char donor_strand, acceptor_strand;
  int sensedir, donor_sensedir, acceptor_sensedir;
  bool allocp;

  debug(printf("Entered print_exon_exon with hardclip_low %d, and hardclip_high %d\n",
	       hardclip_low,hardclip_high));

  sensedir = Stage3end_sensedir(this);
  donor = Stage3end_substring_donor(this);
  acceptor = Stage3end_substring_acceptor(this);

  querylength = Shortread_fulllength(queryseq);

  /* Shouldn't have any overlap on a distant splice */
  hardclip_low = hardclip_high = 0;

  donor_chrpos = Substring_compute_chrpos(donor,hardclip_low,hide_soft_clips_p);
  acceptor_chrpos = Substring_compute_chrpos(acceptor,hardclip_low,hide_soft_clips_p);

  halfdonor_dinucleotide(&donor1,&donor2,donor,sensedir);
  halfacceptor_dinucleotide(&acceptor2,&acceptor1,acceptor,sensedir);
  donor_chr = Univ_IIT_label(chromosome_iit,Substring_chrnum(donor),&allocp);
  acceptor_chr = Univ_IIT_label(chromosome_iit,Substring_chrnum(acceptor),&allocp);
  donor_prob = Substring_chimera_prob(donor);
  acceptor_prob = Substring_chimera_prob(acceptor);

  /* Code taken from that for XS tag for print_halfdonor and print_halfacceptor */
  /* For the donor and acceptor strands, use the substring sensedir and not the Stage3end_T sensedir */
  if ((donor_sensedir = Substring_chimera_sensedir(donor)) == SENSE_FORWARD) {
    if (Substring_plusp(donor) == true) {
      donor_strand = '+';
    } else {
      donor_strand = '-';
    }
  } else if (donor_sensedir == SENSE_ANTI) {
    if (Substring_plusp(donor) == true) {
      donor_strand = '-';
    } else {
      donor_strand = '+';
    }
  } else {
    abort();
  }

  if ((acceptor_sensedir = Substring_chimera_sensedir(acceptor)) == SENSE_FORWARD) {
    if (Substring_plusp(acceptor) == true) {
      acceptor_strand = '+';
    } else {
      acceptor_strand = '-';
    }
  } else if (acceptor_sensedir == SENSE_ANTI) {
    if (Substring_plusp(acceptor) == true) {
      acceptor_strand = '-';
    } else {
      acceptor_strand = '+';
    }
  } else {
    abort();
  }

  if (sensedir == SENSE_FORWARD) {

    /* NEEDS WORK: Need to decide whether to split halfdonor or halfacceptor */
    /* Not sure if circular chromosomes should participate in distant splicing anyway */
    if (0 && (circularpos = Stage3end_circularpos(this)) > 0) {
      print_halfdonor(fp,abbrev,donor,this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,
		      chromosome_iit,queryseq,pairedlength,
		      concordant_chrpos,donor_chrpos,acceptor_chrpos,mate_chrpos,
		      /*hardclip_low*/0,/*hardclip_high*/querylength-circularpos,
		      resulttype,first_read_p,artificial_mate_p,npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
		      donor_sensedir,donor_strand,acceptor_strand,donor_chr,acceptor_chr,
		      donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob,
		      /*circularp*/true);
      print_halfdonor(fp,abbrev,donor,this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,
		      chromosome_iit,queryseq,pairedlength,
		      /*concordant_chrpos*/1,/*donor_chrpos*/1,acceptor_chrpos,mate_chrpos,
		      /*hardclip_low*/circularpos,/*hardclip_high*/0,
		      resulttype,first_read_p,artificial_mate_p,npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
		      donor_sensedir,donor_strand,acceptor_strand,donor_chr,acceptor_chr,
		      donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob,
		      /*circularp*/true);
    } else {
      print_halfdonor(fp,abbrev,donor,this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,
		      chromosome_iit,queryseq,pairedlength,
		      concordant_chrpos,donor_chrpos,acceptor_chrpos,mate_chrpos,
		      hardclip_low,hardclip_high,resulttype,first_read_p,
		      artificial_mate_p,npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
		      donor_sensedir,donor_strand,acceptor_strand,donor_chr,acceptor_chr,
		      donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob,
		      /*circularp*/false);
    }

    if (0 && (circularpos = Stage3end_circularpos(this)) > 0) {
      print_halfacceptor(fp,abbrev,acceptor,this,mate,acc1,acc2,pathnum,npaths,
			 absmq_score,first_absmq,second_absmq,mapq_score,
			 chromosome_iit,queryseq,pairedlength,
			 concordant_chrpos,donor_chrpos,acceptor_chrpos,mate_chrpos,
			 /*hardclip_low*/0,/*hardclip_high*/querylength-circularpos,
			 resulttype,first_read_p,artificial_mate_p,npaths_mate,quality_shift,sam_read_group_id,
			 invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
			 acceptor_sensedir,donor_strand,acceptor_strand,donor_chr,acceptor_chr,
			 donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob,
			 /*circularp*/true);
      print_halfacceptor(fp,abbrev,acceptor,this,mate,acc1,acc2,pathnum,npaths,
			 absmq_score,first_absmq,second_absmq,mapq_score,
			 chromosome_iit,queryseq,pairedlength,
			 /*concordant_chrpos*/1,donor_chrpos,/*acceptor_chrpos*/1,mate_chrpos,
			 /*hardclip_low*/circularpos,/*hardclip_high*/0,
			 resulttype,first_read_p,artificial_mate_p,npaths_mate,quality_shift,sam_read_group_id,
			 invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
			 acceptor_sensedir,donor_strand,acceptor_strand,donor_chr,acceptor_chr,
			 donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob,
			 /*circularp*/true);
    } else {
      print_halfacceptor(fp,abbrev,acceptor,this,mate,acc1,acc2,pathnum,npaths,
			 absmq_score,first_absmq,second_absmq,mapq_score,
			 chromosome_iit,queryseq,pairedlength,
			 concordant_chrpos,donor_chrpos,acceptor_chrpos,mate_chrpos,
			 hardclip_low,hardclip_high,resulttype,first_read_p,
			 artificial_mate_p,npaths_mate,quality_shift,sam_read_group_id,
			 invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
			 acceptor_sensedir,donor_strand,acceptor_strand,donor_chr,acceptor_chr,
			 donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob,
			 /*circularp*/false);
    }

  } else if (Stage3end_sensedir(this) == SENSE_ANTI) {
    if (0 && (circularpos = Stage3end_circularpos(this)) > 0) {
      print_halfacceptor(fp,abbrev,acceptor,this,mate,acc1,acc2,pathnum,npaths,
			 absmq_score,first_absmq,second_absmq,mapq_score,
			 chromosome_iit,queryseq,pairedlength,
			 concordant_chrpos,donor_chrpos,acceptor_chrpos,mate_chrpos,
			 /*hardclip_low*/0,/*hardclip_high*/querylength-circularpos,
			 resulttype,first_read_p,artificial_mate_p,npaths_mate,quality_shift,sam_read_group_id,
			 invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
			 acceptor_sensedir,donor_strand,acceptor_strand,donor_chr,acceptor_chr,
			 donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob,
			 /*circularp*/true);
      print_halfacceptor(fp,abbrev,acceptor,this,mate,acc1,acc2,pathnum,npaths,
			 absmq_score,first_absmq,second_absmq,mapq_score,
			 chromosome_iit,queryseq,pairedlength,
			 /*concordant_chrpos*/1,donor_chrpos,/*acceptor_chrpos*/1,mate_chrpos,
			 /*hardclip_low*/circularpos,/*hardclip_high*/0,
			 resulttype,first_read_p,artificial_mate_p,npaths_mate,quality_shift,sam_read_group_id,
			 invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
			 acceptor_sensedir,donor_strand,acceptor_strand,donor_chr,acceptor_chr,
			 donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob,
			 /*circularp*/true);
    } else {
      print_halfacceptor(fp,abbrev,acceptor,this,mate,acc1,acc2,pathnum,npaths,
			 absmq_score,first_absmq,second_absmq,mapq_score,
			 chromosome_iit,queryseq,pairedlength,
			 concordant_chrpos,donor_chrpos,acceptor_chrpos,mate_chrpos,
			 hardclip_low,hardclip_high,resulttype,first_read_p,
			 artificial_mate_p,npaths_mate,quality_shift,sam_read_group_id,
			 invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
			 acceptor_sensedir,donor_strand,acceptor_strand,donor_chr,acceptor_chr,
			 donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob,
			 /*circularp*/false);
    }

    if (0 && (circularpos = Stage3end_circularpos(this)) > 0) {
      print_halfdonor(fp,abbrev,donor,this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,
		      chromosome_iit,queryseq,pairedlength,
		      concordant_chrpos,donor_chrpos,acceptor_chrpos,mate_chrpos,
		      /*hardclip_low*/0,/*hardclip_high*/querylength-circularpos,
		      resulttype,first_read_p,artificial_mate_p,npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
		      donor_sensedir,donor_strand,acceptor_strand,donor_chr,acceptor_chr,
		      donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob,
		      /*circularp*/true);
      print_halfdonor(fp,abbrev,donor,this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,
		      chromosome_iit,queryseq,pairedlength,
		      /*concordant_chrpos*/1,/*donor_chrpos*/1,acceptor_chrpos,mate_chrpos,
		      /*hardclip_low*/circularpos,/*hardclip_high*/0,
		      resulttype,first_read_p,artificial_mate_p,npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
		      donor_sensedir,donor_strand,acceptor_strand,donor_chr,acceptor_chr,
		      donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob,
		      /*circularp*/true);
    } else {
      print_halfdonor(fp,abbrev,donor,this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,
		      chromosome_iit,queryseq,pairedlength,
		      concordant_chrpos,donor_chrpos,acceptor_chrpos,mate_chrpos,
		      hardclip_low,hardclip_high,resulttype,first_read_p,
		      artificial_mate_p,npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
		      donor_sensedir,donor_strand,acceptor_strand,donor_chr,acceptor_chr,
		      donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob,
		      /*circularp*/false);
    }

  } else {
    /* SENSE_NULL (DNA distant chimera) */
    if (0 && (circularpos = Stage3end_circularpos(this)) > 0) {
      print_halfacceptor(fp,abbrev,acceptor,this,mate,acc1,acc2,pathnum,npaths,
			 absmq_score,first_absmq,second_absmq,mapq_score,
			 chromosome_iit,queryseq,pairedlength,
			 concordant_chrpos,donor_chrpos,acceptor_chrpos,mate_chrpos,
			 /*hardclip_low*/0,/*hardclip_high*/querylength-circularpos,
			 resulttype,first_read_p,artificial_mate_p,npaths_mate,quality_shift,sam_read_group_id,
			 invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
			 acceptor_sensedir,donor_strand,acceptor_strand,donor_chr,acceptor_chr,
			 donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob,
			 /*circularp*/true);
      print_halfacceptor(fp,abbrev,acceptor,this,mate,acc1,acc2,pathnum,npaths,
			 absmq_score,first_absmq,second_absmq,mapq_score,
			 chromosome_iit,queryseq,pairedlength,
			 /*concordant_chrpos*/1,donor_chrpos,/*acceptor_chrpos*/1,mate_chrpos,
			 /*hardclip_low*/circularpos,/*hardclip_high*/0,
			 resulttype,first_read_p,artificial_mate_p,npaths_mate,quality_shift,sam_read_group_id,
			 invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
			 acceptor_sensedir,donor_strand,acceptor_strand,donor_chr,acceptor_chr,
			 donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob,
			 /*circularp*/true);
    } else {
      print_halfacceptor(fp,abbrev,acceptor,this,mate,acc1,acc2,pathnum,npaths,
			 absmq_score,first_absmq,second_absmq,mapq_score,
			 chromosome_iit,queryseq,pairedlength,
			 concordant_chrpos,donor_chrpos,acceptor_chrpos,mate_chrpos,
			 hardclip_low,hardclip_high,resulttype,first_read_p,
			 artificial_mate_p,npaths_mate,quality_shift,sam_read_group_id,
			 invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
			 acceptor_sensedir,donor_strand,acceptor_strand,donor_chr,acceptor_chr,
			 donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob,
			 /*circularp*/false);
    }

    if (0 && (circularpos = Stage3end_circularpos(this)) > 0) {
      print_halfdonor(fp,abbrev,donor,this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,
		      chromosome_iit,queryseq,pairedlength,
		      concordant_chrpos,donor_chrpos,acceptor_chrpos,mate_chrpos,
		      /*hardclip_low*/0,/*hardclip_high*/querylength-circularpos,
		      resulttype,first_read_p,artificial_mate_p,npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
		      donor_sensedir,donor_strand,acceptor_strand,donor_chr,acceptor_chr,
		      donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob,
		      /*circularp*/true);
      print_halfdonor(fp,abbrev,donor,this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,
		      chromosome_iit,queryseq,pairedlength,
		      /*concordant_chrpos*/1,/*donor_chrpos*/1,acceptor_chrpos,mate_chrpos,
		      /*hardclip_low*/circularpos,/*hardclip_high*/0,
		      resulttype,first_read_p,artificial_mate_p,npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
		      donor_sensedir,donor_strand,acceptor_strand,donor_chr,acceptor_chr,
		      donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob,
		      /*circularp*/true);
    } else {
      print_halfdonor(fp,abbrev,donor,this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,
		      chromosome_iit,queryseq,pairedlength,
		      concordant_chrpos,donor_chrpos,acceptor_chrpos,mate_chrpos,
		      hardclip_low,hardclip_high,resulttype,first_read_p,
		      artificial_mate_p,npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
		      donor_sensedir,donor_strand,acceptor_strand,donor_chr,acceptor_chr,
		      donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob,
		      /*circularp*/false);
    }
  }

  if (allocp == true) {
    FREE(acceptor_chr);
    FREE(donor_chr);
  }

  return;
}



void
SAM_print (Filestring_T fp, Filestring_T fp_failedinput, char *abbrev,
	   Stage3end_T this, Stage3end_T mate, char *acc1, char *acc2, int pathnum, int npaths,
	   int absmq_score, int first_absmq, int second_absmq, int mapq_score, Univ_IIT_T chromosome_iit, Shortread_T queryseq,
	   Shortread_T queryseq_mate, int pairedlength, Chrpos_T chrpos, Chrpos_T mate_chrpos,
	   int clipdir, int hardclip5_low, int hardclip5_high, int hardclip3_low, int hardclip3_high,
	   Resulttype_T resulttype, bool first_read_p, bool artificial_mate_p, int npaths_mate,
	   int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p,
	   bool merge_samechr_p) {
  Hittype_T hittype;
  unsigned int flag;
  int circularpos, querylength;


  hittype = Stage3end_hittype(this);
  /* printf("hittype %s, chrpos %u\n",Stage3end_hittype_string(this),chrpos); */

  /* Test for nomapping was chrpos == 0, but we can actually align to chrpos 0 */
  /* Also, can use this test here because --quiet-if-excessive cases go directly to SAM_print_nomapping */
  if (npaths == 0) {
    SAM_print_nomapping(fp,abbrev,queryseq,mate,acc1,acc2,chromosome_iit,resulttype,first_read_p,
			/*pathnum*/0,/*npaths*/0,artificial_mate_p,npaths_mate,mate_chrpos,
			quality_shift,sam_read_group_id,invertp,invert_mate_p);

    if (fp_failedinput != NULL) {
      if (first_read_p == true) {
	Shortread_print_query_singleend(fp_failedinput,queryseq,/*headerseq*/queryseq);
      } else {
	Shortread_print_query_singleend(fp_failedinput,queryseq,/*headerseq*/queryseq_mate);
      }
    }

  } else if (hittype == GMAP) {
    /* Note: sam_paired_p must be true because we are calling GMAP only on halfmapping uniq */

    if (mate == NULL) {
      chrpos = SAM_compute_chrpos(/*hardclip_low*/0,/*hardclip_high*/hardclip5_high,
				  this,Shortread_fulllength(queryseq),/*first_read_p*/true);
      mate_chrpos = 0U;
      hardclip3_low = hardclip3_high = 0;

    } else if (first_read_p == true) {
      chrpos = SAM_compute_chrpos(/*hardclip_low*/hardclip5_low,/*hardclip_high*/hardclip5_high,
				  this,Shortread_fulllength(queryseq),/*first_read_p*/true);
      mate_chrpos = SAM_compute_chrpos(/*hardclip_low*/hardclip3_low,/*hardclip_high*/hardclip3_high,
				       mate,Shortread_fulllength(queryseq_mate),/*first_read_p*/false);
    } else {
      chrpos = SAM_compute_chrpos(/*hardclip_low*/hardclip3_low,/*hardclip_high*/hardclip3_high,
				  this,Shortread_fulllength(queryseq),/*first_read_p*/false);
      mate_chrpos = SAM_compute_chrpos(/*hardclip_low*/hardclip5_low,/*hardclip_high*/hardclip5_high,
				       mate,Shortread_fulllength(queryseq_mate),/*first_read_p*/true);
    }

    flag = SAM_compute_flag(Stage3end_plusp(this),mate,resulttype,first_read_p,
			    pathnum,npaths,artificial_mate_p,npaths_mate,
			    absmq_score,first_absmq,invertp,invert_mate_p);

    querylength = Shortread_fulllength(queryseq);
    if ((circularpos = Stage3end_circularpos(this)) > 0
#if 0
	&& Pair_check_cigar(Stage3end_pairarray(this),Stage3end_npairs(this),querylength,
			    /*clipdir*/+1,/*hardclip_low*/0,/*hardclip_high*/querylength-circularpos,
			    /*watsonp*/Stage3end_plusp(this),Stage3end_sensedir(this),
			    first_read_p,/*circularp*/true) == true &&
	Pair_check_cigar(Stage3end_pairarray(this),Stage3end_npairs(this),querylength,
			 /*clipdir*/+1,/*hardclip_low*/circularpos,/*hardclip_high*/0,
			 /*watsonp*/Stage3end_plusp(this),Stage3end_sensedir(this),
			 first_read_p,/*circularp*/true) == true
#endif
	) {
      Pair_print_sam(fp,abbrev,Stage3end_pairarray(this),Stage3end_npairs(this),
		     Stage3end_cigar_tokens(this),Stage3end_gmap_intronp(this),
		     acc1,acc2,Stage3end_chrnum(this),chromosome_iit,/*usersegment*/(Sequence_T) NULL,
		     Shortread_fullpointer(queryseq),Shortread_quality_string(queryseq),
		     /*clipdir*/+1,/*hardclip_low*/0,/*hardclip_high*/querylength-circularpos,Shortread_fulllength(queryseq),
		     /*watsonp*/Stage3end_plusp(this),Stage3end_sensedir(this),
		     /*chimera_part*/0,/*chimera*/NULL,quality_shift,first_read_p,
		     pathnum,npaths,absmq_score,first_absmq,second_absmq,chrpos,Stage3end_chrlength(this),
		     queryseq,resulttype,flag,/*pair_mapq_score*/mapq_score,/*end_mapq_score*/mapq_score,
		     Stage3end_chrnum(mate),Stage3end_effective_chrnum(mate),
		     mate_chrpos,Stage3end_chrlength(mate),/*mate_sensedir*/Stage3end_sensedir(mate),
		     pairedlength,sam_read_group_id,invertp,/*circularp*/true,/*merged_overlap_p*/false,
		     Stage3end_sarrayp(this));
      Pair_print_sam(fp,abbrev,Stage3end_pairarray(this),Stage3end_npairs(this),
		     Stage3end_cigar_tokens(this),Stage3end_gmap_intronp(this),
		     acc1,acc2,Stage3end_chrnum(this),chromosome_iit,/*usersegment*/(Sequence_T) NULL,
		     Shortread_fullpointer(queryseq),Shortread_quality_string(queryseq),
		     /*clipdir*/+1,/*hardclip_low*/circularpos,/*hardclip_high*/0,Shortread_fulllength(queryseq),
		     /*watsonp*/Stage3end_plusp(this),Stage3end_sensedir(this),
		     /*chimera_part*/0,/*chimera*/NULL,quality_shift,first_read_p,
		     pathnum,npaths,absmq_score,first_absmq,second_absmq,/*chrpos*/1,Stage3end_chrlength(this),
		     queryseq,resulttype,flag,/*pair_mapq_score*/mapq_score,/*end_mapq_score*/mapq_score,
		     Stage3end_chrnum(mate),Stage3end_effective_chrnum(mate),
		     mate_chrpos,Stage3end_chrlength(mate),/*mate_sensedir*/Stage3end_sensedir(mate),
		     pairedlength,sam_read_group_id,invertp,/*circularp*/true,/*merged_overlap_p*/false,
		     Stage3end_sarrayp(this));
    } else if (first_read_p == true) {
      Pair_print_sam(fp,abbrev,Stage3end_pairarray(this),Stage3end_npairs(this),
		     Stage3end_cigar_tokens(this),Stage3end_gmap_intronp(this),
		     acc1,acc2,Stage3end_chrnum(this),chromosome_iit,/*usersegment*/(Sequence_T) NULL,
		     Shortread_fullpointer(queryseq),Shortread_quality_string(queryseq),
		     clipdir,hardclip5_low,hardclip5_high,Shortread_fulllength(queryseq),
		     Stage3end_plusp(this),Stage3end_sensedir(this),
		     /*chimera_part*/0,/*chimera*/NULL,quality_shift,/*first_read_p*/true,
		     pathnum,npaths,absmq_score,first_absmq,second_absmq,chrpos,Stage3end_chrlength(this),
		     queryseq,resulttype,flag,/*pair_mapq_score*/mapq_score,/*end_mapq_score*/mapq_score,
		     Stage3end_chrnum(mate),Stage3end_effective_chrnum(mate),
		     mate_chrpos,Stage3end_chrlength(mate),/*mate_sensedir*/Stage3end_sensedir(mate),
		     pairedlength,sam_read_group_id,invertp,/*circularp*/false,/*merged_overlap_p*/false,
		     Stage3end_sarrayp(this));
    } else {
      Pair_print_sam(fp,abbrev,Stage3end_pairarray(this),Stage3end_npairs(this),
		     Stage3end_cigar_tokens(this),Stage3end_gmap_intronp(this),
		     acc1,acc2,Stage3end_chrnum(this),chromosome_iit,/*usersegment*/(Sequence_T) NULL,
		     Shortread_fullpointer(queryseq),Shortread_quality_string(queryseq),
		     clipdir,hardclip3_low,hardclip3_high,Shortread_fulllength(queryseq),
		     Stage3end_plusp(this),Stage3end_sensedir(this),
		     /*chimera_part*/0,/*chimera*/NULL,quality_shift,/*first_read_p*/false,
		     pathnum,npaths,absmq_score,first_absmq,second_absmq,chrpos,Stage3end_chrlength(this),
		     queryseq,resulttype,flag,/*pair_mapq_score*/mapq_score,/*end_mapq_score*/mapq_score,
		     Stage3end_chrnum(mate),Stage3end_effective_chrnum(mate),
		     mate_chrpos,Stage3end_chrlength(mate),/*mate_sensedir*/Stage3end_sensedir(mate),
		     pairedlength,sam_read_group_id,invertp,/*circularp*/false,/*merged_overlap_p*/false,
		     Stage3end_sarrayp(this));
    }

  } else if (hittype == TRANSLOC_SPLICE || (hittype == SAMECHR_SPLICE && merge_samechr_p == false)) {
    if (first_read_p == true) {
      print_exon_exon(fp,abbrev,this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		      mate_chrpos,hardclip5_low,hardclip5_high,resulttype,/*first_read_p*/true,
		      artificial_mate_p,npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p);
    } else {
      print_exon_exon(fp,abbrev,this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		      mate_chrpos,hardclip3_low,hardclip3_high,resulttype,/*first_read_p*/false,
		      artificial_mate_p,npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p);
    }

  } else {
    querylength = Shortread_fulllength(queryseq);
    if ((circularpos = Stage3end_circularpos(this)) > 0
#if 0
	&& check_cigar_single(hittype,this,querylength,/*hardclip_low*/0,/*hardclip_high*/querylength-circularpos) == true &&
	check_cigar_single(hittype,this,querylength,/*hardclip_low*/circularpos,/*hardclip_high*/0) == true
#endif
	) {
      print_substrings(fp,abbrev,this,mate,acc1,acc2,pathnum,npaths,
		       absmq_score,first_absmq,second_absmq,mapq_score,queryseq,pairedlength,
		       chrpos,mate_chrpos,/*hardclip_low*/0,/*hardclip_high*/querylength-circularpos,
		       resulttype,first_read_p,artificial_mate_p,npaths_mate,quality_shift,sam_read_group_id,
		       invertp,invert_mate_p,/*circularp*/true);
      print_substrings(fp,abbrev,this,mate,acc1,acc2,pathnum,npaths,
		       absmq_score,first_absmq,second_absmq,mapq_score,queryseq,pairedlength,
		       /*chrpos*/1,mate_chrpos,/*hardclip_low*/circularpos,/*hardclip_high*/0,
		       resulttype,first_read_p,artificial_mate_p,npaths_mate,quality_shift,sam_read_group_id,
		       invertp,invert_mate_p,/*circularp*/true);
    } else if (first_read_p == true) {
      print_substrings(fp,abbrev,this,mate,acc1,acc2,pathnum,npaths,
		       absmq_score,first_absmq,second_absmq,mapq_score,queryseq,pairedlength,
		       chrpos,mate_chrpos,hardclip5_low,hardclip5_high,resulttype,/*first_read_p*/true,
		       artificial_mate_p,npaths_mate,quality_shift,sam_read_group_id,
		       invertp,invert_mate_p,/*circularp*/false);
    } else {
      print_substrings(fp,abbrev,this,mate,acc1,acc2,pathnum,npaths,
		       absmq_score,first_absmq,second_absmq,mapq_score,queryseq,pairedlength,
		       chrpos,mate_chrpos,hardclip3_low,hardclip3_high,resulttype,/*first_read_p*/false,
		       artificial_mate_p,npaths_mate,quality_shift,sam_read_group_id,
		       invertp,invert_mate_p,/*circularp*/false);
    }
  }

  return;
}



void
SAM_print_paired (Filestring_T fp, Filestring_T fp_failedinput_1, Filestring_T fp_failedinput_2,
		  Result_T result, Resulttype_T resulttype, Univ_IIT_T chromosome_iit,
		  Shortread_T queryseq1, Shortread_T queryseq2, bool invert_first_p, bool invert_second_p,
		  bool nofailsp, bool failsonlyp, bool merge_samechr_p,
		  int quality_shift, char *sam_read_group_id) {
  Stage3pair_T *stage3pairarray, stage3pair;
  Stage3end_T *stage3array1, *stage3array2, stage3, mate, hit5, hit3;
  Chrpos_T chrpos, chrpos5, chrpos3;
  int npaths, npaths_max, npaths1, npaths2, pathnum;
  int first_absmq, second_absmq, first_absmq1, second_absmq1, first_absmq2, second_absmq2;
  int hardclip5_low = 0, hardclip5_high = 0, hardclip3_low = 0, hardclip3_high = 0, clipdir;
  char *acc1, *acc2;
  Pairtype_T pairtype;
  char *abbrev;

  struct Pair_T *pairarray;
  int npairs;
  char *queryseq_merged, *quality_merged;
  int querylength_merged;
  int flag;

  acc1 = Shortread_accession(queryseq1);
  acc2 = Shortread_accession(queryseq2); /* NULL, unless --allow-pe-name-mismatch is specified */

  if (resulttype == PAIREDEND_NOMAPPING) {
    if (nofailsp == true) {
      /* No output */
      return;
      
    } else 
      Filestring_set_split_output(fp,OUTPUT_NM);
      SAM_print_nomapping(fp,ABBREV_NOMAPPING_1,queryseq1,/*mate*/(Stage3end_T) NULL,
			  acc1,acc2,chromosome_iit,resulttype,
			  /*first_read_p*/true,/*pathnum*/0,/*npaths*/0,
			  /*artificial_mate_p*/false,/*npaths_mate*/0,
			  /*mate_chrpos*/0U,quality_shift,
			  sam_read_group_id,invert_first_p,invert_second_p);
      SAM_print_nomapping(fp,ABBREV_NOMAPPING_2,queryseq2,/*mate*/(Stage3end_T) NULL,
			  acc1,acc2,chromosome_iit,resulttype,
			  /*first_read_p*/false,/*pathnum*/0,/*npaths*/0,
			  /*artificial_mate_p*/false,/*npaths_mate*/0,
			  /*mate_chrpos*/0U,quality_shift,
			  sam_read_group_id,invert_second_p,invert_first_p);

      if (fp_failedinput_1 != NULL) {
	Shortread_print_query_pairedend(fp_failedinput_1,fp_failedinput_2,queryseq1,queryseq2);
      }

  } else {
    if (failsonlyp == true) {
      /* Unwanted success: skip */

    } else if (resulttype == CONCORDANT_UNIQ) {
      stage3pairarray = (Stage3pair_T *) Result_array(&npaths,&first_absmq,&second_absmq,result);
      /* Stage3pair_eval(stage3pairarray,npaths,maxpaths_report,queryseq1,queryseq2); */

      stage3pair = stage3pairarray[0];
      hit5 = Stage3pair_hit5(stage3pair);
      hit3 = Stage3pair_hit3(stage3pair);

      if (Stage3pair_circularp(stage3pair) == true) {
	/* Don't resolve overlaps on a circular alignment */
	clipdir = 0;
	hardclip5_low = hardclip5_high = hardclip3_low = hardclip3_high = 0;
	Filestring_set_split_output(fp,OUTPUT_CC);
	abbrev = ABBREV_CONCORDANT_CIRCULAR;

      } else if (clip_overlap_p == false && merge_overlap_p == false) {
	clipdir = 0;
	hardclip5_low = hardclip5_high = hardclip3_low = hardclip3_high = 0;
	Filestring_set_split_output(fp,OUTPUT_CU);
	abbrev = ABBREV_CONCORDANT_UNIQ;

      } else {
	clipdir = Stage3pair_overlap(&hardclip5_low,&hardclip5_high,&hardclip3_low,&hardclip3_high,stage3pair);
	debug3(printf("clipdir %d with hardclip5 = %d..%d, hardclip3 = %d..%d\n",
		      clipdir,hardclip5_low,hardclip5_high,hardclip3_low,hardclip3_high));
	Filestring_set_split_output(fp,OUTPUT_CU);
	abbrev = ABBREV_CONCORDANT_UNIQ;
      }

      chrpos5 = SAM_compute_chrpos(hardclip5_low,hardclip5_high,hit5,Shortread_fulllength(queryseq1),/*first_read_p*/true);
      chrpos3 = SAM_compute_chrpos(hardclip3_low,hardclip3_high,hit3,Shortread_fulllength(queryseq2),/*first_read_p*/false);

      if (merge_overlap_p == false || clipdir == 0) {
	/* print first end */
	SAM_print(fp,fp_failedinput_1,abbrev,hit5,/*mate*/hit3,
		  acc1,acc2,/*pathnum*/1,/*npaths*/1,
		  Stage3pair_absmq_score(stage3pair),first_absmq,/*second_absmq*/0,
		  Stage3pair_mapq_score(stage3pair),chromosome_iit,
		  /*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		  Stage3pair_pairlength(stage3pair),/*chrpos*/chrpos5,/*mate_chrpos*/chrpos3,
		  clipdir,hardclip5_low,hardclip5_high,hardclip3_low,hardclip3_high,
		  resulttype,/*first_read_p*/true,/*artificial_mate_p*/false,/*npaths_mate*/npaths,
		  quality_shift,sam_read_group_id,invert_first_p,invert_second_p,
		  merge_samechr_p);

	/* print second end */
	SAM_print(fp,fp_failedinput_2,abbrev,hit3,/*mate*/hit5,
		  acc1,acc2,/*pathnum*/1,/*npaths*/1,
		  Stage3pair_absmq_score(stage3pair),first_absmq,/*second_absmq*/0,
		  Stage3pair_mapq_score(stage3pair),chromosome_iit,
		  /*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		  Stage3pair_pairlength(stage3pair),/*chrpos*/chrpos3,/*mate_chrpos*/chrpos5,
		  clipdir,hardclip5_low,hardclip5_high,hardclip3_low,hardclip3_high,
		  resulttype,/*first_read_p*/false,/*artificial_mate_p*/false,/*npaths_mate*/npaths,
		  quality_shift,sam_read_group_id,invert_second_p,invert_first_p,
		  merge_samechr_p);

      } else {
	/* merge_overlap_p == true and overlap was found */
	pairarray = Stage3pair_merge(&npairs,&querylength_merged,&queryseq_merged,&quality_merged,
				     stage3pair,queryseq1,queryseq2,
				     /*querylength5*/Stage3end_querylength(hit5),
				     /*querylength3*/Stage3end_querylength(hit3),
				     clipdir,hardclip5_low,hardclip5_high,hardclip3_low,hardclip3_high);
	/* printf("queryseq_merged: %s\n",queryseq_merged); */
	if (clipdir >= 0) {
	  chrpos = chrpos5;
	} else {
	  chrpos = chrpos3;
	}
	/* merging changes resulttype from UNPAIRED_UNIQ to SINGLEEND_UNIQ */
	flag = SAM_compute_flag(Stage3end_plusp(hit5),/*mate*/NULL,/*resulttype*/SINGLEEND_UNIQ,/*first_read_p*/true,
				/*pathnum*/1,/*npaths*/1,/*artificial_mate_p*/false,/*npaths_mate*/0,
				Stage3pair_absmq_score(stage3pair),first_absmq,/*invertp*/false,
				/*invert_mate_p*/false);
	Filestring_set_split_output(fp,OUTPUT_UU);
	Pair_print_sam(fp,/*abbrev*/ABBREV_UNPAIRED_UNIQ,pairarray,npairs,/*cigar_tokens*/NULL,/*gmap_intronp*/false,
		       acc1,/*acc2*/NULL,Stage3end_chrnum(hit5),chromosome_iit,/*usersegment*/(Sequence_T) NULL,
		       /*queryseq_ptr*/queryseq_merged,/*quality_string*/quality_merged,
		       /*clipdir*/0,/*hardclip_low*/0,/*hardclip_high*/0,/*querylength*/querylength_merged,
		       Stage3end_plusp(hit5),Stage3end_sensedir(hit5),
		       /*chimera_part*/0,/*chimera*/NULL,quality_shift,/*first_read_p*/true,
		       /*pathnum*/1,/*npaths*/1,
#if 0
		       Stage3pair_absmq_score(stage3pair),first_absmq,/*second_absmq*/0,
#else
		       /*absmq_score*/MAX_QUALITY_SCORE,/*first_absmq*/MAX_QUALITY_SCORE,/*second_absmq*/0,
#endif
		       chrpos,Stage3end_chrlength(hit5),/*queryseq*/NULL,resulttype,flag,
		       /*pair_mapq_score*/MAX_QUALITY_SCORE,/*end_mapq_score*/MAX_QUALITY_SCORE,
		       /*mate_chrnum*/0,/*mate_effective_chrnum*/0,/*mate_chrpos*/0,/*mate_chrlength*/0,
		       /*mate_sensedir*/SENSE_NULL,/*pairedlength*/0,
		       sam_read_group_id,/*invertp*/false,/*circularp*/false,/*merged_overlap_p*/true,
		       Stage3end_sarrayp(hit5));
	if (quality_merged != NULL) {
	  FREE_OUT(quality_merged);
	}
	FREE_OUT(queryseq_merged);
	FREE_OUT(pairarray);
      }

    } else if (resulttype == CONCORDANT_TRANSLOC) {
      Filestring_set_split_output(fp,OUTPUT_CT);
      stage3pairarray = (Stage3pair_T *) Result_array(&npaths,&first_absmq,&second_absmq,result);

      if (quiet_if_excessive_p && npaths > maxpaths_report) {
	/* Print as nomapping, but send to fp_concordant_transloc */
	SAM_print_nomapping(fp,ABBREV_CONCORDANT_TRANSLOC,
			    queryseq1,/*mate*/(Stage3end_T) NULL,
			    acc1,acc2,chromosome_iit,resulttype,
			    /*first_read_p*/true,/*pathnum*/1,npaths,
			    /*artificial_mate_p*/false,/*npaths_mate*/npaths,/*mate_chrpos*/0U,
			    quality_shift,sam_read_group_id,invert_first_p,invert_second_p);
	SAM_print_nomapping(fp,ABBREV_CONCORDANT_TRANSLOC,
			    queryseq2,/*mate*/(Stage3end_T) NULL,
			    acc1,acc2,chromosome_iit,resulttype,
			    /*first_read_p*/false,/*pathnum*/1,npaths,
			    /*artificial_mate_p*/false,/*npaths_mate*/npaths,/*mate_chrpos*/0U,
			    quality_shift,sam_read_group_id,invert_second_p,invert_first_p);

      } else {
	/* Stage3pair_eval(stage3pairarray,npaths,maxpaths_report,queryseq1,queryseq2); */
	/* Note: We are ignoring merge_overlap for concordant_transloc */
	for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths_report; pathnum++) {

	  stage3pair = stage3pairarray[pathnum-1];
	  hit5 = Stage3pair_hit5(stage3pair);
	  hit3 = Stage3pair_hit3(stage3pair);

	  if (Stage3pair_circularp(stage3pair) == true) {
	    /* Don't resolve overlaps on a circular alignment */
	    clipdir = 0;
	    hardclip5_low = hardclip5_high = hardclip3_low = hardclip3_high = 0;

	  } else if (clip_overlap_p == false) {
	    clipdir = 0;
	    hardclip5_low = hardclip5_high = hardclip3_low = hardclip3_high = 0;

	  } else {
	    clipdir = Stage3pair_overlap(&hardclip5_low,&hardclip5_high,&hardclip3_low,&hardclip3_high,stage3pair);
	    debug3(printf("clipdir %d with hardclip5 = %d..%d, hardclip3 = %d..%d\n",
			  clipdir,hardclip5_low,hardclip5_high,hardclip3_low,hardclip3_high));
	  }

	  chrpos5 = SAM_compute_chrpos(hardclip5_low,hardclip5_high,hit5,Shortread_fulllength(queryseq1),/*first_read_p*/true);
	  chrpos3 = SAM_compute_chrpos(hardclip3_low,hardclip3_high,hit3,Shortread_fulllength(queryseq2),/*first_read_p*/false);

	  /* print first end */
	  SAM_print(fp,fp_failedinput_1,ABBREV_CONCORDANT_TRANSLOC,
		    hit5,/*mate*/hit3,acc1,acc2,pathnum,npaths,
		    Stage3pair_absmq_score(stage3pair),first_absmq,second_absmq,
		    Stage3pair_mapq_score(stage3pair),chromosome_iit,
		    /*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		    Stage3pair_pairlength(stage3pair),/*chrpos*/chrpos5,/*mate_chrpos*/chrpos3,
		    clipdir,hardclip5_low,hardclip5_high,hardclip3_low,hardclip3_high,
		    resulttype,/*first_read_p*/true,/*artificial_mate_p*/false,/*npaths_mate*/npaths,
		    quality_shift,sam_read_group_id,invert_first_p,invert_second_p,
		    merge_samechr_p);

	  /* print second end */
	  SAM_print(fp,fp_failedinput_2,ABBREV_CONCORDANT_TRANSLOC,
		    hit3,/*mate*/hit5,acc1,acc2,pathnum,npaths,
		    Stage3pair_absmq_score(stage3pair),first_absmq,second_absmq,
		    Stage3pair_mapq_score(stage3pair),chromosome_iit,
		    /*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		    Stage3pair_pairlength(stage3pair),/*chrpos*/chrpos3,/*mate_chrpos*/chrpos5,
		    clipdir,hardclip5_low,hardclip5_high,hardclip3_low,hardclip3_high,
		    resulttype,/*first_read_p*/false,/*artificial_mate_p*/false,/*npaths_mate*/npaths,
		    quality_shift,sam_read_group_id,invert_second_p,invert_first_p,
		    merge_samechr_p);
	}
      }
    
    } else if (resulttype == CONCORDANT_MULT) {
      stage3pairarray = (Stage3pair_T *) Result_array(&npaths,&first_absmq,&second_absmq,result);

      if (quiet_if_excessive_p && npaths > maxpaths_report) {
	/* Print as nomapping, but send to fp_concordant_mult_xs */
	Filestring_set_split_output(fp,OUTPUT_CX);
	SAM_print_nomapping(fp,ABBREV_CONCORDANT_MULT_XS,
			    queryseq1,/*mate*/(Stage3end_T) NULL,
			    acc1,acc2,chromosome_iit,resulttype,
			    /*first_read_p*/true,/*pathnum*/1,npaths,
			    /*artificial_mate_p*/false,/*npaths_mate*/npaths,
			    /*mate_chrpos*/0U,quality_shift,
			    sam_read_group_id,invert_first_p,invert_second_p);
	SAM_print_nomapping(fp,ABBREV_CONCORDANT_MULT_XS,
			    queryseq2,/*mate*/(Stage3end_T) NULL,
			    acc1,acc2,chromosome_iit,resulttype,
			    /*first_read_p*/false,/*pathnum*/1,npaths,
			    /*artificial_mate_p*/false,/*npaths_mate*/npaths,
			    /*mate_chrpos*/0U,quality_shift,
			    sam_read_group_id,invert_second_p,invert_first_p);

	if (fp_failedinput_1 != NULL) {
	  Shortread_print_query_pairedend(fp_failedinput_1,fp_failedinput_2,queryseq1,queryseq2);

	}

      } else {
	/* Stage3pair_eval(stage3pairarray,npaths,maxpaths_report,queryseq1,queryseq2); */
	Filestring_set_split_output(fp,OUTPUT_CM);
	for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths_report; pathnum++) {

	  stage3pair = stage3pairarray[pathnum-1];
	  hit5 = Stage3pair_hit5(stage3pair);
	  hit3 = Stage3pair_hit3(stage3pair);

	  if (Stage3pair_circularp(stage3pair) == true) {
	    /* Don't resolve overlaps on a circular alignment */
	    clipdir = 0;
	    hardclip5_low = hardclip5_high = hardclip3_low = hardclip3_high = 0;

	  } else if (clip_overlap_p == false && merge_overlap_p == false) {
	    clipdir = 0;
	    hardclip5_low = hardclip5_high = hardclip3_low = hardclip3_high = 0;

	  } else {
	    clipdir = Stage3pair_overlap(&hardclip5_low,&hardclip5_high,&hardclip3_low,&hardclip3_high,stage3pair);
	    debug3(printf("clipdir %d with hardclip5 = %d..%d, hardclip3 = %d..%d\n",
			  clipdir,hardclip5_low,hardclip5_high,hardclip3_low,hardclip3_high));
	  }

	  chrpos5 = SAM_compute_chrpos(hardclip5_low,hardclip5_high,hit5,Shortread_fulllength(queryseq1),/*first_read_p*/true);
	  chrpos3 = SAM_compute_chrpos(hardclip3_low,hardclip3_high,hit3,Shortread_fulllength(queryseq2),/*first_read_p*/false);

	  if (merge_overlap_p == false || clipdir == 0) {
	    /* print first end */
	    SAM_print(fp,fp_failedinput_1,ABBREV_CONCORDANT_MULT,
		      hit5,/*mate*/hit3,acc1,acc2,pathnum,npaths,
		      Stage3pair_absmq_score(stage3pair),first_absmq,second_absmq,
		      Stage3pair_mapq_score(stage3pair),chromosome_iit,
		      /*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		      Stage3pair_pairlength(stage3pair),/*chrpos*/chrpos5,/*mate_chrpos*/chrpos3,
		      clipdir,hardclip5_low,hardclip5_high,hardclip3_low,hardclip3_high,
		      resulttype,/*first_read_p*/true,/*artificial_mate_p*/false,/*npaths_mate*/npaths,
		      quality_shift,sam_read_group_id,invert_first_p,invert_second_p,
		      merge_samechr_p);

	    /* print second end */
	    SAM_print(fp,fp_failedinput_2,ABBREV_CONCORDANT_MULT,
		      hit3,/*mate*/hit5,acc1,acc2,pathnum,npaths,
		      Stage3pair_absmq_score(stage3pair),first_absmq,second_absmq,
		      Stage3pair_mapq_score(stage3pair),chromosome_iit,
		      /*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		      Stage3pair_pairlength(stage3pair),/*chrpos*/chrpos3,/*mate_chrpos*/chrpos5,
		      clipdir,hardclip5_low,hardclip5_high,hardclip3_low,hardclip3_high,
		      resulttype,/*first_read_p*/false,/*artificial_mate_p*/false,/*npaths_mate*/npaths,
		      quality_shift,sam_read_group_id,invert_second_p,invert_first_p,
		      merge_samechr_p);
	    
	  } else {
	    /* merge_overlap_p == true and overlap was found */
	    pairarray = Stage3pair_merge(&npairs,&querylength_merged,&queryseq_merged,&quality_merged,
					 stage3pair,queryseq1,queryseq2,
					 /*querylength5*/Stage3end_querylength(hit5),
					 /*querylength3*/Stage3end_querylength(hit3),
					 clipdir,hardclip5_low,hardclip5_high,hardclip3_low,hardclip3_high);
	    /* printf("queryseq_merged: %s\n",queryseq_merged); */
	    if (clipdir >= 0) {
	      chrpos = chrpos5;
	    } else {
	      chrpos = chrpos3;
	    }
	    /* merging changes resulttype from UNPAIRED_UNIQ to SINGLEEND_UNIQ */
	    flag = SAM_compute_flag(Stage3end_plusp(hit5),/*mate*/NULL,/*resulttype*/SINGLEEND_UNIQ,/*first_read_p*/true,
				    /*pathnum*/1,/*npaths*/1,/*artificial_mate_p*/false,/*npaths_mate*/0,
				    Stage3pair_absmq_score(stage3pair),first_absmq,/*invertp*/false,
				    /*invert_mate_p*/false);
	    Pair_print_sam(fp,ABBREV_CONCORDANT_MULT,pairarray,npairs,/*cigar_tokens*/NULL,/*gmap_intronp*/false,
			   acc1,/*acc2*/NULL,Stage3end_chrnum(hit5),chromosome_iit,
			   /*usersegment*/(Sequence_T) NULL,
			   /*queryseq_ptr*/queryseq_merged,/*quality_string*/quality_merged,
			   /*clipdir*/0,/*hardclip_low*/0,/*hardclip_high*/0,/*querylength*/querylength_merged,
			   Stage3end_plusp(hit5),Stage3end_sensedir(hit5),
			   /*chimera_part*/0,/*chimera*/NULL,quality_shift,/*first_read_p*/true,pathnum,npaths,
#if 0
			   Stage3pair_absmq_score(stage3pair),first_absmq,/*second_absmq*/0,
#else
			   /*absmq_score*/MAX_QUALITY_SCORE,/*first_absmq*/MAX_QUALITY_SCORE,/*second_absmq*/0,
#endif
			   chrpos,Stage3end_chrlength(hit5),/*queryseq*/NULL,resulttype,flag,
			   /*pair_mapq_score*/MAX_QUALITY_SCORE,/*end_mapq_score*/MAX_QUALITY_SCORE,
			   /*mate_chrnum*/0,/*mate_effective_chrnum*/0,/*mate_chrpos*/0,/*mate_chrlength*/0,
			   /*mate_sensedir*/SENSE_NULL,/*pairedlength*/0,
			   sam_read_group_id,/*invertp*/false,/*circularp*/false,/*merged_overlap_p*/true,
			   Stage3end_sarrayp(hit5));
	    if (quality_merged != NULL) {
	      FREE_OUT(quality_merged);
	    }
	    FREE_OUT(queryseq_merged);
	    FREE_OUT(pairarray);
	  }
	}
      }

    } else if (resulttype == PAIRED_UNIQ) {
      stage3pairarray = (Stage3pair_T *) Result_array(&npaths,&first_absmq,&second_absmq,result);
      /* Stage3pair_eval(stage3pairarray,npaths,maxpaths_report,queryseq1,queryseq2); */

      stage3pair = stage3pairarray[0];
      if (Stage3pair_circularp(stage3pair) == true) {
	Filestring_set_split_output(fp,OUTPUT_PC);
	abbrev = ABBREV_PAIRED_UNIQ_CIRCULAR;
      } else if ((pairtype = Stage3pair_pairtype(stage3pair)) == PAIRED_INVERSION) {
	Filestring_set_split_output(fp,OUTPUT_PI);
	abbrev = ABBREV_PAIRED_UNIQ_INV;
      } else if (pairtype == PAIRED_SCRAMBLE) {
	Filestring_set_split_output(fp,OUTPUT_PS);
	abbrev = ABBREV_PAIRED_UNIQ_SCR;
      } else if (pairtype == PAIRED_TOOLONG) {
	Filestring_set_split_output(fp,OUTPUT_PL);
	abbrev = ABBREV_PAIRED_UNIQ_LONG;
      } else {
	fprintf(stderr,"Unexpected pairtype %d\n",pairtype);
	abort();
      }

      hardclip5_low = hardclip5_high = hardclip3_low = hardclip3_high = 0;

      hit5 = Stage3pair_hit5(stage3pair);
      hit3 = Stage3pair_hit3(stage3pair);
      chrpos5 = SAM_compute_chrpos(hardclip5_low,hardclip5_high,hit5,Shortread_fulllength(queryseq1),/*first_read_p*/true);
      chrpos3 = SAM_compute_chrpos(hardclip3_low,hardclip3_high,hit3,Shortread_fulllength(queryseq2),/*first_read_p*/false);

      /* print first end */
      SAM_print(fp,fp_failedinput_1,abbrev,hit5,/*mate*/hit3,
		acc1,acc2,/*pathnum*/1,/*npaths*/1,
		Stage3pair_absmq_score(stage3pair),first_absmq,/*second_absmq*/0,
		Stage3pair_mapq_score(stage3pair),chromosome_iit,
		/*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		Stage3pair_pairlength(stage3pair),/*chrpos*/chrpos5,/*mate_chrpos*/chrpos3,
		/*clipdir*/0,hardclip5_low,hardclip5_high,hardclip3_low,hardclip3_high,
		resulttype,/*first_read_p*/true,/*artificial_mate_p*/false,/*npaths_mate*/npaths,
		quality_shift,sam_read_group_id,invert_first_p,invert_second_p,
		merge_samechr_p);

      /* print second end */
      SAM_print(fp,fp_failedinput_2,abbrev,hit3,/*mate*/hit5,
		acc1,acc2,/*pathnum*/1,/*npaths*/1,
		Stage3pair_absmq_score(stage3pair),first_absmq,/*second_absmq*/0,
		Stage3pair_mapq_score(stage3pair),chromosome_iit,
		/*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		Stage3pair_pairlength(stage3pair),/*chrpos*/chrpos3,/*mate_chrpos*/chrpos5,
		/*clipdir*/0,hardclip5_low,hardclip5_high,hardclip3_low,hardclip3_high,
		resulttype,/*first_read_p*/false,/*artificial_mate_p*/false,/*npaths_mate*/npaths,
		quality_shift,sam_read_group_id,invert_second_p,invert_first_p,
		merge_samechr_p);

    } else if (resulttype == PAIRED_MULT) {
      stage3pairarray = (Stage3pair_T *) Result_array(&npaths,&first_absmq,&second_absmq,result);

      if (quiet_if_excessive_p && npaths > maxpaths_report) {
	/* Print as nomapping, but send to fp_paired_mult */
	Filestring_set_split_output(fp,OUTPUT_PX);
	SAM_print_nomapping(fp,ABBREV_PAIRED_MULT_XS,
			    queryseq1,/*mate*/(Stage3end_T) NULL,
			    acc1,acc2,chromosome_iit,resulttype,
			    /*first_read_p*/true,/*pathnum*/1,npaths,
			    /*artificial_mate_p*/false,/*npaths_mate*/npaths,
			    /*mate_chrpos*/0U,quality_shift,
			    sam_read_group_id,invert_first_p,invert_second_p);
	SAM_print_nomapping(fp,ABBREV_PAIRED_MULT_XS,
			    queryseq2,/*mate*/(Stage3end_T) NULL,
			    acc1,acc2,chromosome_iit,resulttype,
			    /*first_read_p*/false,/*pathnum*/1,npaths,
			    /*artificial_mate_p*/false,/*npaths_mate*/npaths,
			    /*mate_chrpos*/0U,quality_shift,
			    sam_read_group_id,invert_second_p,invert_first_p);

	if (fp_failedinput_1 != NULL) {
	  Shortread_print_query_pairedend(fp_failedinput_1,fp_failedinput_2,queryseq1,queryseq2);
	}

      } else {
	/* Stage3pair_eval(stage3pairarray,npaths,maxpaths_report,queryseq1,queryseq2); */
	Filestring_set_split_output(fp,OUTPUT_PM);
	for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths_report; pathnum++) {

	  stage3pair = stage3pairarray[pathnum-1];
	  hardclip5_low = hardclip5_high = hardclip3_low = hardclip3_high = 0;

	  hit5 = Stage3pair_hit5(stage3pair);
	  hit3 = Stage3pair_hit3(stage3pair);
	  chrpos5 = SAM_compute_chrpos(hardclip5_low,hardclip5_high,hit5,Shortread_fulllength(queryseq1),/*first_read_p*/true);
	  chrpos3 = SAM_compute_chrpos(hardclip3_low,hardclip3_high,hit3,Shortread_fulllength(queryseq2),/*first_read_p*/false);

	  /* print first end */
	  SAM_print(fp,fp_failedinput_1,ABBREV_PAIRED_MULT,
		    hit5,/*mate*/hit3,acc1,acc2,pathnum,npaths,
		    Stage3pair_absmq_score(stage3pair),first_absmq,second_absmq,
		    Stage3pair_mapq_score(stage3pair),chromosome_iit,
		    /*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		    Stage3pair_pairlength(stage3pair),/*chrpos*/chrpos5,/*mate_chrpos*/chrpos3,
		    /*clipdir*/0,hardclip5_low,hardclip5_high,hardclip3_low,hardclip3_high,
		    resulttype,/*first_read_p*/true,/*artificial_mate_p*/false,/*npaths_mate*/npaths,
		    quality_shift,sam_read_group_id,invert_first_p,invert_second_p,
		    merge_samechr_p);

	  /* print second end */
	  SAM_print(fp,fp_failedinput_2,ABBREV_PAIRED_MULT,
		    hit3,/*mate*/hit5,acc1,acc2,pathnum,npaths,
		    Stage3pair_absmq_score(stage3pair),first_absmq,second_absmq,
		    Stage3pair_mapq_score(stage3pair),chromosome_iit,
		    /*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		    Stage3pair_pairlength(stage3pair),/*chrpos*/chrpos3,/*mate_chrpos*/chrpos5,
		    /*clipdir*/0,hardclip5_low,hardclip5_high,hardclip3_low,hardclip3_high,
		    resulttype,/*first_read_p*/false,/*artificial_mate_p*/false,/*npaths_mate*/npaths,
		    quality_shift,sam_read_group_id,invert_second_p,invert_first_p,
		    merge_samechr_p);
	}
      }

    } else if (resulttype == UNPAIRED_UNIQ) {
      /* Even though they are not related, we should print mate information in this situation */
      stage3array1 = (Stage3end_T *) Result_array(&npaths1,&first_absmq1,&second_absmq1,result);
      stage3array2 = (Stage3end_T *) Result_array2(&npaths2,&first_absmq2,&second_absmq2,result);

      hardclip5_low = hardclip5_high = hardclip3_low = hardclip3_high = 0;

      hit5 = stage3array1[0];
      hit3 = stage3array2[0];
      chrpos5 = SAM_compute_chrpos(hardclip5_low,hardclip5_high,hit5,Shortread_fulllength(queryseq1),/*first_read_p*/true);
      chrpos3 = SAM_compute_chrpos(hardclip3_low,hardclip3_high,hit3,Shortread_fulllength(queryseq2),/*first_read_p*/false);

      if (Stage3end_circularpos(hit5) > 0 || Stage3end_circularpos(hit3) > 0) {
	Filestring_set_split_output(fp,OUTPUT_UC);
	abbrev = ABBREV_UNPAIRED_CIRCULAR;
      } else {
	Filestring_set_split_output(fp,OUTPUT_UU);
	abbrev = ABBREV_UNPAIRED_UNIQ;
      }

      /* print first end */
      /* Stage3end_eval_and_sort(stage3array1,npaths1,maxpaths_report,queryseq1); */
      SAM_print(fp,fp_failedinput_1,abbrev,hit5,/*mate*/hit3,acc1,acc2,/*pathnum*/1,/*npaths*/1,
		Stage3end_absmq_score(stage3array1[0]),first_absmq1,/*second_absmq*/0,
		Stage3end_mapq_score(stage3array1[0]),chromosome_iit,
		/*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		/*pairedlength*/0U,/*chrpos*/chrpos5,/*mate_chrpos*/chrpos3,
		/*clipdir*/0,/*hardclip5_low*/0,/*hardclip5_high*/0,/*hardclip3_low*/0,/*hardclip3_high*/0,
		resulttype,/*first_read_p*/true,/*artificial_mate_p*/false,/*npaths_mate*/1,quality_shift,sam_read_group_id,
		invert_first_p,invert_second_p,merge_samechr_p);

      /* Note: Do not act on add_paired_nomappers_p, since the two reads are artificially paired up already */

      /* print second end */
      /* Stage3end_eval_and_sort(stage3array2,npaths2,maxpaths_report,queryseq2); */
      SAM_print(fp,fp_failedinput_2,abbrev,hit3,/*mate*/hit5,acc1,acc2,/*pathnum*/1,/*npaths*/1,
		Stage3end_absmq_score(stage3array2[0]),first_absmq2,/*second_absmq*/0,
		Stage3end_mapq_score(stage3array2[0]),chromosome_iit,
		/*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		/*pairedlength*/0U,/*chrpos*/chrpos3,/*mate_chrpos*/chrpos5,
		/*clipdir*/0,/*hardclip5_low*/0,/*hardclip5_high*/0,/*hardclip3_low*/0,/*hardclip3_high*/0,
		resulttype,/*first_read_p*/false,/*artificial_mate_p*/false,/*npaths_mate*/1,quality_shift,sam_read_group_id,
		invert_second_p,invert_first_p,merge_samechr_p);

    } else if (resulttype == UNPAIRED_MULT || resulttype == UNPAIRED_TRANSLOC) {
      if (resulttype == UNPAIRED_MULT) {
	if (quiet_if_excessive_p && npaths1 > maxpaths_report && npaths2 > maxpaths_report) {
	  Filestring_set_split_output(fp,OUTPUT_UX);
	} else {
	  Filestring_set_split_output(fp,OUTPUT_UM);
	}
	abbrev = ABBREV_UNPAIRED_MULT;
      } else {
	Filestring_set_split_output(fp,OUTPUT_UT);
	abbrev = ABBREV_UNPAIRED_TRANSLOC;
      }

      stage3array1 = (Stage3end_T *) Result_array(&npaths1,&first_absmq1,&second_absmq1,result);
      stage3array2 = (Stage3end_T *) Result_array2(&npaths2,&first_absmq2,&second_absmq2,result);

#if 0
      /* Do eval and sorting first */
      if (npaths1 == 1) {
	/* Stage3end_eval_and_sort(stage3array1,npaths1,maxpaths_report,queryseq1); */
      } else if (quiet_if_excessive_p && npaths1 > maxpaths_report) {
	/* Don't sort */
      } else {
	/* Stage3end_eval_and_sort(stage3array1,npaths1,maxpaths_report,queryseq1); */
      }
#endif

#if 0
      if (npaths2 == 1) {
	/* Stage3end_eval_and_sort(stage3array2,npaths2,maxpaths_report,queryseq2); */
      } else if (quiet_if_excessive_p && npaths2 > maxpaths_report) {
	/* Don't sort */
      } else {
	/* Stage3end_eval_and_sort(stage3array2,npaths2,maxpaths_report,queryseq2); */
      }
#endif

      if (add_paired_nomappers_p == true) {
	/* Artificially pair up results */
	npaths_max = (npaths1 > npaths2) ? npaths1 : npaths2;
	for (pathnum = 1; pathnum <= npaths1 && pathnum <= npaths2 && pathnum <= maxpaths_report; pathnum++) {
	  /* hardclip5_low = hardclip5_high = 0; */
	  chrpos5 = SAM_compute_chrpos(/*hardclip_low*/0,/*hardclip_high*/0,/*stage3*/stage3array1[pathnum-1],
				       Shortread_fulllength(queryseq1),/*first_read_p*/true);

	  /* hardclip3_low = hardclip3_high = 0; */
	  chrpos3 = SAM_compute_chrpos(/*hardclip_low*/0,/*hardclip_high*/0,/*stage3*/stage3array2[pathnum-1],
				       Shortread_fulllength(queryseq2),/*first_read_p*/false);

	  stage3 = stage3array1[pathnum-1];
	  SAM_print(fp,fp_failedinput_1,abbrev,stage3,/*mate*/stage3array2[pathnum-1],acc1,acc2,pathnum,npaths_max,
		    Stage3end_absmq_score(stage3),first_absmq1,second_absmq1,
		    Stage3end_mapq_score(stage3),chromosome_iit,
		    /*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		    /*pairedlength*/0U,/*chrpos*/chrpos5,/*mate_chrpos*/chrpos3,
		    /*clipdir*/0,/*hardclip5_low*/0,/*hardclip5_high*/0,/*hardclip3_low*/0,/*hardclip3_high*/0,
		    resulttype,/*first_read_p*/true,/*artificial_mate_p*/false,/*npaths_mate*/npaths_max,
		    quality_shift,sam_read_group_id,invert_first_p,invert_second_p,merge_samechr_p);

	  stage3 = stage3array2[pathnum-1];
	  SAM_print(fp,fp_failedinput_2,abbrev,stage3,/*mate*/stage3array1[pathnum-1],acc1,acc2,pathnum,npaths_max,
		    Stage3end_absmq_score(stage3),first_absmq2,second_absmq2,
		    Stage3end_mapq_score(stage3),chromosome_iit,
		    /*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		    /*pairedlength*/0U,/*chrpos*/chrpos3,/*mate_chrpos*/chrpos5,
		    /*clipdir*/0,/*hardclip5_low*/0,/*hardclip5_high*/0,/*hardclip3_low*/0,/*hardclip3_high*/0,
		    resulttype,/*first_read_p*/false,/*artificial_mate_p*/false,/*npaths_mate*/npaths_max,
		    quality_shift,sam_read_group_id,invert_second_p,invert_first_p,merge_samechr_p);
	}

	/* Print remaining results with non-mappers */
	if (npaths1 > npaths2) {
	  for ( ; pathnum <= npaths1 && pathnum <= maxpaths_report; pathnum++) {
	    stage3 = stage3array1[pathnum-1];
	    /* hardclip5_low = hardclip5_high = 0; */
	    chrpos5 = SAM_compute_chrpos(/*hardclip_low*/0,/*hardclip_high*/0,stage3,
					 Shortread_fulllength(queryseq1),/*first_read_p*/true);
	    chrpos3 = 0;

	    SAM_print(fp,fp_failedinput_1,abbrev,stage3,/*mate*/NULL,acc1,acc2,pathnum,npaths_max,
		      Stage3end_absmq_score(stage3),first_absmq1,second_absmq1,
		      Stage3end_mapq_score(stage3),chromosome_iit,
		      /*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		      /*pairedlength*/0U,/*chrpos*/chrpos5,/*mate_chrpos*/chrpos3,
		      /*clipdir*/0,/*hardclip5_low*/0,/*hardclip5_high*/0,/*hardclip3_low*/0,/*hardclip3_high*/0,
		      resulttype,/*first_read_p*/true,/*artificial_mate_p*/true,/*npaths_mate*/npaths_max,
		      quality_shift,sam_read_group_id,invert_first_p,invert_second_p,merge_samechr_p);

	    /* matching nomapper for second end */
	    SAM_print_nomapping(fp,abbrev,queryseq2,/*mate*/stage3,acc1,acc2,chromosome_iit,
				resulttype,/*first_read_p*/false,pathnum,npaths_max,
				/*artificial_mate_p*/false,/*npaths_mate*/npaths_max,/*mate_chrpos*/chrpos5,
				quality_shift,sam_read_group_id,invert_first_p,invert_second_p);
	  }

	} else if (npaths2 > npaths1) {
	  for ( ; pathnum <= npaths2 && pathnum <= maxpaths_report; pathnum++) {
	    stage3 = stage3array2[pathnum-1];
	    /* hardclip3_low = hardclip3_high = 0; */
	    chrpos5 = 0;
	    chrpos3 = SAM_compute_chrpos(/*hardclip_low*/0,/*hardclip_high*/0,stage3,
					 Shortread_fulllength(queryseq2),/*first_read_p*/false);

	    /* matching nomapper for first end */
	    SAM_print_nomapping(fp,abbrev,queryseq1,/*mate*/stage3,acc1,acc2,chromosome_iit,
				resulttype,/*first_read_p*/true,pathnum,npaths_max,
				/*artificial_mate_p*/false,/*npaths_mate*/npaths_max,/*mate_chrpos*/chrpos3,
				quality_shift,sam_read_group_id,invert_first_p,invert_second_p);

	    SAM_print(fp,fp_failedinput_2,abbrev,stage3,/*mate*/NULL,acc1,acc2,pathnum,npaths_max,
		      Stage3end_absmq_score(stage3),first_absmq2,second_absmq2,
		      Stage3end_mapq_score(stage3),chromosome_iit,
		      /*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		      /*pairedlength*/0U,/*chrpos*/chrpos3,/*mate_chrpos*/chrpos5,
		      /*clipdir*/0,/*hardclip5_low*/0,/*hardclip5_high*/0,/*hardclip3_low*/0,/*hardclip3_high*/0,
		      resulttype,/*first_read_p*/false,/*artificial_mate_p*/true,/*npaths_mate*/npaths_max,
		      quality_shift,sam_read_group_id,invert_second_p,invert_first_p,merge_samechr_p);
	  }
	}

      } else {
	/* print first end results */
	if (npaths2 == 0) {
	  mate = (Stage3end_T) NULL;
	  chrpos3 = 0U;
	} else if (quiet_if_excessive_p && npaths2 > maxpaths_report) {
	  mate = (Stage3end_T) NULL;
	  chrpos3 = 0U;
	} else {
	  mate = stage3array2[0];
	  hardclip3_low = hardclip3_high = 0;
	  chrpos3 = SAM_compute_chrpos(/*hardclip_low*/0,/*hardclip_high*/0,mate,Shortread_fulllength(queryseq2),/*first_read_p*/false);
	}

	if (npaths1 == 1) {
	  stage3 = stage3array1[0];
	  hardclip5_low = hardclip5_high = 0;
	  chrpos5 = SAM_compute_chrpos(/*hardclip_low*/0,/*hardclip_high*/0,stage3,Shortread_fulllength(queryseq1),/*first_read_p*/true);

	  SAM_print(fp,fp_failedinput_1,abbrev,stage3,mate,acc1,acc2,/*pathnum*/1,npaths1,
		    Stage3end_absmq_score(stage3),first_absmq1,second_absmq1,
		    Stage3end_mapq_score(stage3),chromosome_iit,
		    /*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		    /*pairedlength*/0U,/*chrpos*/chrpos5,/*mate_chrpos*/chrpos3,
		    /*clipdir*/0,/*hardclip5_low*/0,/*hardclip5_high*/0,/*hardclip3_low*/0,/*hardclip3_high*/0,
		    resulttype,/*first_read_p*/true,/*artificial_mate_p*/false,/*npaths_mate*/npaths2,
		    quality_shift,sam_read_group_id,invert_first_p,invert_second_p,merge_samechr_p);

	} else if (quiet_if_excessive_p && npaths1 > maxpaths_report) {
	  /* Just printing one end as nomapping */
	  SAM_print_nomapping(fp,abbrev,queryseq1,mate,acc1,acc2,chromosome_iit,
			      resulttype,/*first_read_p*/true,/*pathnum*/1,npaths1,
			      /*artificial_mate_p*/false,/*npaths_mate*/npaths2,/*mate_chrpos*/chrpos3,
			      quality_shift,sam_read_group_id,invert_first_p,invert_second_p);

	} else {
	  for (pathnum = 1; pathnum <= npaths1 && pathnum <= maxpaths_report; pathnum++) {
	    stage3 = stage3array1[pathnum-1];
	    hardclip5_low = hardclip5_high = 0;
	    chrpos5 = SAM_compute_chrpos(/*hardclip_low*/0,/*hardclip_high*/hardclip5_high,stage3,Shortread_fulllength(queryseq1),
					 /*first_read_p*/true);
	    
	    SAM_print(fp,fp_failedinput_1,abbrev,stage3,mate,acc1,acc2,pathnum,npaths1,
		      Stage3end_absmq_score(stage3),first_absmq1,second_absmq1,
		      Stage3end_mapq_score(stage3),chromosome_iit,
		      /*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		      /*pairedlength*/0U,/*chrpos*/chrpos5,/*mate_chrpos*/chrpos3,
		      /*clipdir*/0,/*hardclip5_low*/0,/*hardclip5_high*/0,/*hardclip3_low*/0,/*hardclip3_high*/0,
		      resulttype,/*first_read_p*/true,/*artificial_mate_p*/false,/*npaths_mate*/npaths2,
		      quality_shift,sam_read_group_id,invert_first_p,invert_second_p,merge_samechr_p);
	  }
	}
			  
	/* print second end results */
	if (npaths1 == 0) {
	  mate = (Stage3end_T) NULL;
	  chrpos5 = 0U;
	} else if (quiet_if_excessive_p && npaths1 > maxpaths_report) {
	  mate = (Stage3end_T) NULL;
	  chrpos5 = 0U;
	} else {
	  mate = stage3array1[0];
	  hardclip5_low = hardclip5_high = 0;
	  chrpos5 = SAM_compute_chrpos(/*hardclip_low*/0,/*hardclip_high*/0,mate,Shortread_fulllength(queryseq1),
				       /*first_read_p*/true);
	}

	if (npaths2 == 1) {
	  stage3 = stage3array2[0];
	  hardclip3_low = hardclip3_high = 0;
	  chrpos3 = SAM_compute_chrpos(/*hardclip_low*/0,/*hardclip_high*/0,stage3,Shortread_fulllength(queryseq2),
				       /*first_read_p*/false);
	  
	  SAM_print(fp,fp_failedinput_2,abbrev,stage3,mate,acc1,acc2,/*pathnum*/1,npaths2,
		    Stage3end_absmq_score(stage3),first_absmq2,second_absmq2,
		    Stage3end_mapq_score(stage3),chromosome_iit,
		    /*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		    /*pairedlength*/0U,/*chrpos*/chrpos3,/*mate_chrpos*/chrpos5,
		    /*clipdir*/0,/*hardclip5_low*/0,/*hardclip5_high*/0,/*hardclip3_low*/0,/*hardclip3_high*/0,
		    resulttype,/*first_read_p*/false,/*artificial_mate_p*/false,/*npaths_mate*/npaths1,
		    quality_shift,sam_read_group_id,invert_second_p,invert_first_p,merge_samechr_p);
	  
	} else if (quiet_if_excessive_p && npaths2 > maxpaths_report) {
	  /* Just printing one end as nomapping */
	  SAM_print_nomapping(fp,abbrev,queryseq2,mate,acc1,acc2,chromosome_iit,
			      resulttype,/*first_read_p*/false,/*pathnum*/1,npaths2,
			      /*artificial_mate_p*/false,/*npaths_mate*/npaths1,/*mate_chrpos*/chrpos5,
			      quality_shift,sam_read_group_id,invert_second_p,invert_first_p);
	  
	} else {
	  for (pathnum = 1; pathnum <= npaths2 && pathnum <= maxpaths_report; pathnum++) {
	    stage3 = stage3array2[pathnum-1];
	    hardclip3_low = hardclip3_high = 0;
	    chrpos3 = SAM_compute_chrpos(/*hardclip_low*/0,/*hardclip_high*/0,stage3,Shortread_fulllength(queryseq2),
					 /*first_read_p*/false);

	    SAM_print(fp,fp_failedinput_2,abbrev,stage3,mate,acc1,acc2,pathnum,npaths2,
		      Stage3end_absmq_score(stage3),first_absmq2,second_absmq2,
		      Stage3end_mapq_score(stage3),chromosome_iit,
		      /*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		      /*pairedlength*/0U,/*chrpos*/chrpos3,/*mate_chrpos*/chrpos5,
		      /*clipdir*/0,/*hardclip5_low*/0,/*hardclip5_high*/0,/*hardclip3_low*/0,/*hardclip3_high*/0,
		      resulttype,/*first_read_p*/false,/*artificial_mate_p*/false,/*npaths_mate*/npaths1,
		      quality_shift,sam_read_group_id,invert_second_p,invert_first_p,merge_samechr_p);
	  }
	}
      }

    } else {
      stage3array1 = (Stage3end_T *) Result_array(&npaths1,&first_absmq1,&second_absmq1,result);
      stage3array2 = (Stage3end_T *) Result_array2(&npaths2,&first_absmq2,&second_absmq2,result);

      if (resulttype == HALFMAPPING_UNIQ) {
	if (npaths1 == 1 && Stage3end_circularpos(stage3array1[0]) > 0) {
	  Filestring_set_split_output(fp,OUTPUT_HC);
	  abbrev = ABBREV_HALFMAPPING_CIRCULAR;
	} else if (npaths2 == 1 && Stage3end_circularpos(stage3array2[0]) > 0) {
	  Filestring_set_split_output(fp,OUTPUT_HC);
	  abbrev = ABBREV_HALFMAPPING_CIRCULAR;
	} else {
	  Filestring_set_split_output(fp,OUTPUT_HU);
	  abbrev = ABBREV_HALFMAPPING_UNIQ;
	}
      } else if (resulttype == HALFMAPPING_TRANSLOC) {
	Filestring_set_split_output(fp,OUTPUT_HT);
	abbrev = ABBREV_HALFMAPPING_TRANSLOC;
      } else if (resulttype == HALFMAPPING_MULT) {
	if (quiet_if_excessive_p == true && npaths1 > maxpaths_report && npaths2 > maxpaths_report) {
	  Filestring_set_split_output(fp,OUTPUT_HX);
	  abbrev = ABBREV_HALFMAPPING_MULT_XS;
	} else {
	  Filestring_set_split_output(fp,OUTPUT_HM);
	  abbrev = ABBREV_HALFMAPPING_MULT;
	}
      } else {
	abort();
      }

#if 0
      /* Do eval and sorting first */
      if (npaths1 == 0) {
	/* Nothing to sort */
      } else if (npaths1 == 1) {
	/* Stage3end_eval_and_sort(stage3array1,npaths1,maxpaths_report,queryseq1); */
      } else if (quiet_if_excessive_p && npaths1 > maxpaths_report) {
	/* Don't sort */
      } else {
	/* Stage3end_eval_and_sort(stage3array1,npaths1,maxpaths_report,queryseq1); */
      }
#endif

#if 0
      if (npaths2 == 0) {
	/* Nothing to sort */
      } else if (npaths2 == 1) {
	/* Stage3end_eval_and_sort(stage3array2,npaths2,maxpaths_report,queryseq2); */
      } else if (quiet_if_excessive_p && npaths2 > maxpaths_report) {
	/* Don't sort */
      } else {
	/* Stage3end_eval_and_sort(stage3array2,npaths2,maxpaths_report,queryseq2); */
      }
#endif


      /* print first end results */
      if (npaths2 == 0) {
	mate = (Stage3end_T) NULL;
	chrpos3 = 0U;
      } else if (quiet_if_excessive_p && npaths2 > maxpaths_report) {
	mate = (Stage3end_T) NULL;
	chrpos3 = 0U;
      } else {
	mate = stage3array2[0];
	hardclip3_low = hardclip3_high = 0;
	chrpos3 = SAM_compute_chrpos(/*hardclip_low*/0,/*hardclip_high*/0,mate,Shortread_fulllength(queryseq2),
				     /*first_read_p*/false);
      }

      if (npaths1 == 0) {
	/* just printing one end as nomapping */
	/* mate should be non-NULL here */
	if (add_paired_nomappers_p == true) {
	  /* Handle nomappers with each mapped mate */
	} else {
	  SAM_print_nomapping(fp,abbrev,queryseq1,mate,acc1,acc2,chromosome_iit,resulttype,
			      /*first_read_p*/true,/*pathnum*/0,npaths1,
			      /*artificial_mate_p*/false,/*npaths_mate*/npaths2,/*mate_chrpos*/chrpos3,
			      quality_shift,sam_read_group_id,invert_first_p,invert_second_p);
	}

      } else if (npaths1 == 1) {
	/* mate should be NULL here */

	stage3 = stage3array1[0];
	hardclip5_low = hardclip5_high = 0;
	chrpos5 = SAM_compute_chrpos(/*hardclip_low*/0,/*hardclip_high*/0,stage3,Shortread_fulllength(queryseq1),
				     /*first_read_p*/true);

	if (add_paired_nomappers_p == true) {
	  /* matching nomapper for second end */
	  npaths_max = npaths1;	/* since npaths2 == 0 */
	  SAM_print(fp,fp_failedinput_1,abbrev,stage3,mate,acc1,acc2,/*pathnum*/1,npaths_max,
		    Stage3end_absmq_score(stage3),first_absmq1,/*second_absmq1*/0,
		    Stage3end_mapq_score(stage3),chromosome_iit,
		    /*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		    /*pairedlength*/0U,/*chrpos*/chrpos5,/*mate_chrpos*/chrpos3,
		    /*clipdir*/0,/*hardclip5_low*/0,/*hardclip5_high*/0,/*hardclip3_low*/0,/*hardclip3_high*/0,
		    resulttype,/*first_read_p*/true,/*artificial_mate_p*/true,/*npaths_mate*/npaths_max,
		    quality_shift,sam_read_group_id,invert_first_p,invert_second_p,merge_samechr_p);
	  SAM_print_nomapping(fp,abbrev,queryseq2,/*mate*/stage3,acc1,acc2,chromosome_iit,
			      resulttype,/*first_read_p*/false,/*pathnum*/1,npaths1,
			      /*artificial_mate_p*/false,/*npaths_mate*/npaths_max,/*mate_chrpos*/chrpos5,
			      quality_shift,sam_read_group_id,invert_first_p,invert_second_p);
	} else {
	  SAM_print(fp,fp_failedinput_1,abbrev,stage3,mate,acc1,acc2,/*pathnum*/1,npaths1,
		    Stage3end_absmq_score(stage3),first_absmq1,/*second_absmq1*/0,
		    Stage3end_mapq_score(stage3),chromosome_iit,
		    /*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		    /*pairedlength*/0U,/*chrpos*/chrpos5,/*mate_chrpos*/chrpos3,
		    /*clipdir*/0,/*hardclip5_low*/0,/*hardclip5_high*/0,/*hardclip3_low*/0,/*hardclip3_high*/0,
		    resulttype,/*first_read_p*/true,/*artificial_mate_p*/false,/*npaths_mate*/npaths2,
		    quality_shift,sam_read_group_id,invert_first_p,invert_second_p,merge_samechr_p);
	}

      } else if (quiet_if_excessive_p && npaths1 > maxpaths_report) {
	/* Just printing one end as nomapping */
	/* mate should be NULL here */
	if (add_paired_nomappers_p == true) {
	  /* Handle nomappers with each mapped mate */
	} else {
	  SAM_print_nomapping(fp,abbrev,queryseq1,mate,acc1,acc2,chromosome_iit,resulttype,
			      /*first_read_p*/true,/*pathnum*/1,npaths1,
			      /*artificial_mate_p*/false,/*npaths_mate*/npaths2,/*mate_chrpos*/chrpos3,
			      quality_shift,sam_read_group_id,invert_first_p,invert_second_p);
	}

      } else {
	/* mate should be NULL here */
	for (pathnum = 1; pathnum <= npaths1 && pathnum <= maxpaths_report; pathnum++) {
	  stage3 = stage3array1[pathnum-1];
	  hardclip5_low = hardclip5_high = 0;
	  chrpos5 = SAM_compute_chrpos(/*hardclip_low*/0,/*hardclip_high*/0,stage3,Shortread_fulllength(queryseq1),
				       /*first_read_p*/true);

	  if (add_paired_nomappers_p == true) {
	    /* matching nomapper for second end */
	    npaths_max = npaths1; /* since npaths2 == 0 */
	    SAM_print_nomapping(fp,abbrev,queryseq2,/*mate*/stage3,acc1,acc2,chromosome_iit,
				resulttype,/*first_read_p*/false,pathnum,
				npaths1,/*artificial_mate_p*/false,/*npaths_mate*/npaths_max,/*mate_chrpos*/chrpos5,
				quality_shift,sam_read_group_id,invert_first_p,invert_second_p);
	    SAM_print(fp,fp_failedinput_1,abbrev,stage3,mate,acc1,acc2,pathnum,npaths_max,
		      Stage3end_absmq_score(stage3),first_absmq1,second_absmq1,
		      Stage3end_mapq_score(stage3),chromosome_iit,
		      /*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		      /*pairedlength*/0U,/*chrpos*/chrpos5,/*mate_chrpos*/chrpos3,
		      /*clipdir*/0,/*hardclip5_low*/0,/*hardclip5_high*/0,/*hardclip3_low*/0,/*hardclip3_high*/0,
		      resulttype,/*first_read_p*/true,/*artificial_mate_p*/true,/*npaths_mate*/npaths_max,quality_shift,sam_read_group_id,
		      invert_first_p,invert_second_p,merge_samechr_p);

	  } else {
	    SAM_print(fp,fp_failedinput_1,abbrev,stage3,mate,acc1,acc2,pathnum,npaths1,
		      Stage3end_absmq_score(stage3),first_absmq1,second_absmq1,
		      Stage3end_mapq_score(stage3),chromosome_iit,
		      /*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		      /*pairedlength*/0U,/*chrpos*/chrpos5,/*mate_chrpos*/chrpos3,
		      /*clipdir*/0,/*hardclip5_low*/0,/*hardclip5_high*/0,/*hardclip3_low*/0,/*hardclip3_high*/0,
		      resulttype,/*first_read_p*/true,/*artificial_mate_p*/false,/*npaths_mate*/npaths2,
		      quality_shift,sam_read_group_id,invert_first_p,invert_second_p,merge_samechr_p);
	  }
	}
      }
			  
      /* print second end results */
      if (npaths1 == 0) {
	mate = (Stage3end_T) NULL;
	chrpos5 = 0U;
      } else if (quiet_if_excessive_p && npaths1 > maxpaths_report) {
	mate = (Stage3end_T) NULL;
	chrpos5 = 0U;
      } else {
	mate = stage3array1[0];
	hardclip5_low = hardclip5_high = 0;
	chrpos5 = SAM_compute_chrpos(/*hardclip_low*/0,/*hardclip_high*/0,mate,Shortread_fulllength(queryseq1),
				     /*first_read_p*/true);
      }

      if (npaths2 == 0) {
	/* Just printing one end as nomapping */
	/* mate should be non-NULL here */
	if (add_paired_nomappers_p == true) {
	  /* Handle nomappers with each mapped mate */
	} else {
	  SAM_print_nomapping(fp,abbrev,queryseq2,mate,acc1,acc2,chromosome_iit,resulttype,
			      /*first_read_p*/false,/*pathnum*/0,npaths2,
			      /*artificial_mate_p*/false,/*npaths_mate*/npaths1,/*mate_chrpos*/chrpos5,
			      quality_shift,sam_read_group_id,invert_second_p,invert_first_p);
	}

      } else if (npaths2 == 1) {
	/* mate should be NULL here */

	stage3 = stage3array2[0];
	hardclip3_low = hardclip3_high = 0;
	chrpos3 = SAM_compute_chrpos(/*hardclip_low*/0,/*hardclip_high*/0,stage3,Shortread_fulllength(queryseq2),
				     /*first_read_p*/false);

	if (add_paired_nomappers_p == true) {
	  /* matching nomapper for first end */
	  npaths_max = npaths2;	/* since npaths1 == 0 */
	  SAM_print_nomapping(fp,abbrev,queryseq2,/*mate*/stage3,acc1,acc2,chromosome_iit,
			      resulttype,/*first_read_p*/true,/*pathnum*/1,
			      npaths2,/*artificial_mate_p*/false,/*npaths_mate*/npaths_max,/*mate_chrpos*/chrpos3,
			      quality_shift,sam_read_group_id,invert_first_p,invert_second_p);
	  SAM_print(fp,fp_failedinput_2,abbrev,stage3,mate,acc1,acc2,/*pathnum*/1,npaths_max,
		    Stage3end_absmq_score(stage3),first_absmq2,/*second_absmq2*/0,
		    Stage3end_mapq_score(stage3),chromosome_iit,
		    /*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		    /*pairedlength*/0U,/*chrpos*/chrpos3,/*mate_chrpos*/chrpos5,
		    /*clipdir*/0,/*hardclip5_low*/0,/*hardclip5_high*/0,/*hardclip3_low*/0,/*hardclip3_high*/0,
		    resulttype,/*first_read_p*/false,/*artificial_mate_p*/true,/*npaths_mate*/npaths_max,
		    quality_shift,sam_read_group_id,invert_second_p,invert_first_p,merge_samechr_p);

	} else {
	  SAM_print(fp,fp_failedinput_2,abbrev,stage3,mate,acc1,acc2,/*pathnum*/1,npaths2,
		    Stage3end_absmq_score(stage3),first_absmq2,/*second_absmq2*/0,
		    Stage3end_mapq_score(stage3),chromosome_iit,
		    /*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		    /*pairedlength*/0U,/*chrpos*/chrpos3,/*mate_chrpos*/chrpos5,
		    /*clipdir*/0,/*hardclip5_low*/0,/*hardclip5_high*/0,/*hardclip3_low*/0,/*hardclip3_high*/0,
		    resulttype,/*first_read_p*/false,/*artificial_mate_p*/false,/*npaths_mate*/npaths1,
		    quality_shift,sam_read_group_id,invert_second_p,invert_first_p,merge_samechr_p);
	}

      } else if (quiet_if_excessive_p && npaths2 > maxpaths_report) {
	/* Just printing one end as nomapping */
	/* mate should be NULL here */
	if (add_paired_nomappers_p == true) {
	  /* Handle nomappers with each mapped mate */
	} else {
	  SAM_print_nomapping(fp,abbrev,queryseq2,mate,acc1,acc2,chromosome_iit,resulttype,
			      /*first_read_p*/false,/*pathnum*/1,npaths2,
			      /*artificial_mate_p*/false,/*npaths_mate*/npaths1,/*mate_chrpos*/chrpos5,
			      quality_shift,sam_read_group_id,invert_second_p,invert_first_p);
	}

      } else {
	/* mate should be NULL here */
	for (pathnum = 1; pathnum <= npaths2 && pathnum <= maxpaths_report; pathnum++) {
	  stage3 = stage3array2[pathnum-1];
	  hardclip3_low = hardclip3_high = 0;
	  chrpos3 = SAM_compute_chrpos(/*hardclip_low*/0,/*hardclip_high*/0,stage3,Shortread_fulllength(queryseq2),
				       /*first_read_p*/false);

	  if (add_paired_nomappers_p == true) {
	    /* matching nomapper for first end */
	    npaths_max = npaths2;	/* since npaths1 == 0 */
	    SAM_print_nomapping(fp,abbrev,queryseq2,/*mate*/stage3,acc1,acc2,chromosome_iit,
				resulttype,/*first_read_p*/true,pathnum,
				npaths2,/*artificial_mate_p*/false,/*npaths_mate*/npaths_max,/*mate_chrpos*/chrpos3,
				quality_shift,sam_read_group_id,invert_first_p,invert_second_p);
	    SAM_print(fp,fp_failedinput_2,abbrev,stage3,mate,acc1,acc2,pathnum,npaths_max,
		      Stage3end_absmq_score(stage3),first_absmq2,second_absmq2,
		      Stage3end_mapq_score(stage3),chromosome_iit,
		      /*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		      /*pairedlength*/0U,/*chrpos*/chrpos3,/*mate_chrpos*/chrpos5,
		      /*clipdir*/0,/*hardclip5_low*/0,/*hardclip5_high*/0,/*hardclip3_low*/0,/*hardclip3_high*/0,
		      resulttype,/*first_read_p*/false,/*artificial_mate_p*/true,/*npaths_mate*/npaths_max,
		      quality_shift,sam_read_group_id,invert_second_p,invert_first_p,merge_samechr_p);

	  } else {
	    SAM_print(fp,fp_failedinput_2,abbrev,stage3,mate,acc1,acc2,pathnum,npaths2,
		      Stage3end_absmq_score(stage3),first_absmq2,second_absmq2,
		      Stage3end_mapq_score(stage3),chromosome_iit,
		      /*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		      /*pairedlength*/0U,/*chrpos*/chrpos3,/*mate_chrpos*/chrpos5,
		      /*clipdir*/0,/*hardclip5_low*/0,/*hardclip5_high*/0,/*hardclip3_low*/0,/*hardclip3_high*/0,
		      resulttype,/*first_read_p*/false,/*artificial_mate_p*/false,/*npaths_mate*/npaths1,
		      quality_shift,sam_read_group_id,invert_second_p,invert_first_p,merge_samechr_p);
	  }
	}
      }

    }
  }

  return;
}



