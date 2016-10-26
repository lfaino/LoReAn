static char rcsid[] = "$Id: output.c 183725 2016-02-04 00:40:15Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "output.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sequence.h"

#ifdef GSNAP
#include "shortread.h"
#include "samprint.h"
#include "stage3hr.h"
#endif

#include "samheader.h"
#include "samflags.h"		/* For output types */


/* For GSNAP, now handling --failsonlyp in sam_sort.  Still
   handling quiet-if-excessive in GMAP/GSNAP, because that changes the
   SAM line to a nomapping format.  */

static Univ_IIT_T chromosome_iit;
static bool nofailsp;
static bool failsonlyp;
static bool quiet_if_excessive_p;
static int maxpaths_report;
static int quality_shift;

#ifdef GSNAP
static bool output_sam_p;
static bool print_m8_p;
static bool invert_first_p;
static bool invert_second_p;

static bool merge_samechr_p;

#else
static Printtype_T printtype;
static int invertmode;
static int wraplength;
static int ngap;
static bool nointronlenp;
static bool sam_paired_p;

static int cds_startpos;
static bool fulllengthp;
static bool truncatep;
static bool strictp;
static bool checksump;

static Genome_T genome;
static Sequence_T usersegment;
static char *user_genomicseg;

static char *dbversion;
static char *chrsubset_name;
static Univ_IIT_T contig_iit;
static IIT_T altstrain_iit;
static bool chimeras_allowed_p;

static IIT_T map_iit;
static int *map_divint_crosstable;
static bool map_exons_p;
static bool map_bothstrands_p;
static int nflanking;
static bool print_comment_p;
#endif

static char *failedinput_root;
static char *sam_read_group_id;


void
Output_setup (Univ_IIT_T chromosome_iit_in,
	      bool nofailsp_in, bool failsonlyp_in, bool quiet_if_excessive_p_in, int maxpaths_report_in,
	      char *failedinput_root_in, int quality_shift_in,
#ifdef GSNAP
	      bool output_sam_p_in, bool print_m8_p_in,	bool invert_first_p_in, bool invert_second_p_in,
	      bool merge_samechr_p_in,
#else
	      Printtype_T printtype_in, int invertmode_in, int wraplength_in, int ngap_in,
	      bool nointronlenp_in, bool sam_paired_p_in, int cds_startpos_in,
	      bool fulllengthp_in, bool truncatep_in, bool strictp_in, bool checksump_in,

	      Genome_T genome_in, Sequence_T usersegment_in, char *user_genomicseg_in,
	      char *dbversion_in, char *chrsubset_name_in,
	      Univ_IIT_T contig_iit_in, IIT_T altstrain_iit_in, bool chimeras_allowed_p_in,
	      IIT_T map_iit_in, int *map_divint_crosstable_in, bool map_exons_p_in,
	      bool map_bothstrands_p_in, int nflanking_in, bool print_comment_p_in,
#endif
	      char *sam_read_group_id_in) {

  chromosome_iit = chromosome_iit_in;

  nofailsp = nofailsp_in;
  failsonlyp = failsonlyp_in;
  quiet_if_excessive_p = quiet_if_excessive_p_in;
  maxpaths_report = maxpaths_report_in;
  failedinput_root = failedinput_root_in;

  quality_shift = quality_shift_in;

#ifdef GSNAP
  output_sam_p = output_sam_p_in;
  print_m8_p = print_m8_p_in;
  invert_first_p = invert_first_p_in;
  invert_second_p = invert_second_p_in;

  merge_samechr_p = merge_samechr_p_in;

#else
  printtype = printtype_in;
  invertmode = invertmode_in;
  wraplength = wraplength_in;
  ngap = ngap_in;
  nointronlenp = nointronlenp_in;
  sam_paired_p = sam_paired_p_in;

  cds_startpos = cds_startpos_in;
  fulllengthp = fulllengthp_in;
  truncatep = truncatep_in;
  strictp = strictp_in;
  checksump = checksump_in;

  genome = genome_in;
  usersegment = usersegment_in;
  user_genomicseg = user_genomicseg_in;

  dbversion = dbversion_in;
  chrsubset_name = chrsubset_name_in;
  contig_iit = contig_iit_in;
  altstrain_iit = altstrain_iit_in;
  chimeras_allowed_p = chimeras_allowed_p_in;

  map_iit = map_iit_in;
  map_divint_crosstable = map_divint_crosstable_in;
  map_exons_p = map_exons_p_in;
  map_bothstrands_p = map_bothstrands_p_in;
  nflanking = nflanking_in;
  print_comment_p = print_comment_p_in;
#endif

  sam_read_group_id = sam_read_group_id_in;

  return;
}


#ifdef GSNAP
/************************************************************************
 *   Print routines and threads for GSNAP
 ************************************************************************/

/* Taken from print_result_sam from old outbuffer.c */
static Filestring_T
filestring_fromresult_sam (Filestring_T *fp_failedinput_1, Filestring_T *fp_failedinput_2,
			   Result_T result, Request_T request) {
  Filestring_T fp;
  Resulttype_T resulttype;
  Shortread_T queryseq1;
  Stage3end_T *stage3array, stage3;
  Chrpos_T chrpos;
  int npaths, pathnum, first_absmq, second_absmq;
  char *abbrev;

  fp = Filestring_new(Request_id(request));
  if (failedinput_root == NULL) {
    *fp_failedinput_1 = (Filestring_T) NULL;
  } else {
    *fp_failedinput_1 = Filestring_new(Request_id(request));
  }

  resulttype = Result_resulttype(result);
  if (resulttype == SINGLEEND_NOMAPPING) {
    *fp_failedinput_2 = (Filestring_T) NULL;
    queryseq1 = Request_queryseq1(request);
    if (nofailsp == true) {
      /* Skip */
    } else {
      Filestring_set_split_output(fp,OUTPUT_NM); /* Needs to go outside of nofailsp */
      SAM_print_nomapping(fp,ABBREV_NOMAPPING_1,
			  queryseq1,/*mate*/NULL,/*acc1*/Shortread_accession(queryseq1),
			  /*acc2*/NULL,chromosome_iit,resulttype,
			  /*first_read_p*/true,/*pathnum*/0,/*npaths*/0,
			  /*artificial_mate_p*/false,/*npaths_mate*/0,/*mate_chrpos*/0U,
			  quality_shift,sam_read_group_id,invert_first_p,invert_second_p);
      if (failedinput_root != NULL) {
	Shortread_print_query_singleend(*fp_failedinput_1,queryseq1,/*headerseq*/queryseq1);
      }
    }

  } else if (resulttype == SINGLEEND_UNIQ) {
    *fp_failedinput_2 = (Filestring_T) NULL;
    if (failsonlyp == true) {
      /* Skip */
    } else {
      queryseq1 = Request_queryseq1(request);

      stage3array = (Stage3end_T *) Result_array(&npaths,&first_absmq,&second_absmq,result);
      stage3 = stage3array[0];
      if (Stage3end_hittype(stage3) == SAMECHR_SPLICE || Stage3end_hittype(stage3) == TRANSLOC_SPLICE) {
	chrpos = 0;
      } else {
	chrpos = SAM_compute_chrpos(/*hardclip_low*/0,/*hardclip_high*/0,stage3,Shortread_fulllength(queryseq1),
				    /*first_read_p*/true);
      }
      if (Stage3end_circularpos(stage3) > 0) {
	Filestring_set_split_output(fp,OUTPUT_UC);
	abbrev = ABBREV_UNPAIRED_CIRCULAR;
      } else {
	Filestring_set_split_output(fp,OUTPUT_UU);
	abbrev = ABBREV_UNPAIRED_UNIQ;
      }
      SAM_print(fp,*fp_failedinput_1,abbrev,stage3,/*mate*/NULL,/*acc1*/Shortread_accession(queryseq1),/*acc2*/NULL,
		/*pathnum*/1,npaths,Stage3end_absmq_score(stage3array[0]),first_absmq,second_absmq,
		Stage3end_mapq_score(stage3array[0]),
		chromosome_iit,queryseq1,/*queryseq2*/NULL,
		/*pairedlength*/0,chrpos,/*mate_chrpos*/0U,
		/*clipdir*/0,/*hardclip5_low*/0,/*hardclip5_high*/0,/*hardclip3_low*/0,/*hardclip3_high*/0,
		resulttype,/*first_read_p*/true,/*artificial_mate_p*/false,/*npaths_mate*/0,quality_shift,
		sam_read_group_id,invert_first_p,invert_second_p,merge_samechr_p);
    }

  } else if (resulttype == SINGLEEND_TRANSLOC) {
    *fp_failedinput_2 = (Filestring_T) NULL;

    Filestring_set_split_output(fp,OUTPUT_UT);
    stage3array = (Stage3end_T *) Result_array(&npaths,&first_absmq,&second_absmq,result);
    if (failsonlyp == true) {
      /* Skip */
    } else if (quiet_if_excessive_p && npaths > maxpaths_report) {
      queryseq1 = Request_queryseq1(request);
      SAM_print_nomapping(fp,ABBREV_UNPAIRED_TRANSLOC,
			  queryseq1,/*mate*/NULL,/*acc1*/Shortread_accession(queryseq1),
			  /*acc2*/NULL,chromosome_iit,resulttype,
			  /*first_read_p*/true,/*pathnum*/1,npaths,
			  /*artificial_mate_p*/false,/*npaths_mate*/0,/*mate_chrpos*/0U,
			  quality_shift,sam_read_group_id,invert_first_p,invert_second_p);
      if (failedinput_root != NULL) {
	Shortread_print_query_singleend(*fp_failedinput_1,queryseq1,/*headerseq*/queryseq1);
      }

    } else {
      queryseq1 = Request_queryseq1(request);

      for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths_report; pathnum++) {
	stage3 = stage3array[pathnum-1];
	if (Stage3end_hittype(stage3) == SAMECHR_SPLICE || Stage3end_hittype(stage3) == TRANSLOC_SPLICE) {
	  chrpos = 0;
	} else {
	  chrpos = SAM_compute_chrpos(/*hardclip_low*/0,/*hardclip_high*/0,stage3,Shortread_fulllength(queryseq1),
				      /*first_read_p*/true);
	}
	SAM_print(fp,*fp_failedinput_1,ABBREV_UNPAIRED_TRANSLOC,
		  stage3,/*mate*/NULL,/*acc1*/Shortread_accession(queryseq1),
		  /*acc2*/NULL,pathnum,npaths,
		  Stage3end_absmq_score(stage3array[pathnum-1]),first_absmq,second_absmq,
		  Stage3end_mapq_score(stage3array[pathnum-1]),
		  chromosome_iit,queryseq1,/*queryseq2*/NULL,
		  /*pairedlength*/0,chrpos,/*mate_chrpos*/0U,
		  /*clipdir*/0,/*hardclip5_low*/0,/*hardclip5_high*/0,/*hardclip3_low*/0,/*hardclip3_high*/0,
		  resulttype,/*first_read_p*/true,/*artificial_mate_p*/false,/*npaths_mate*/0,quality_shift,
		  sam_read_group_id,invert_first_p,invert_second_p,merge_samechr_p);
      }
    }

  } else if (resulttype == SINGLEEND_MULT) {
    *fp_failedinput_2 = (Filestring_T) NULL;
    stage3array = (Stage3end_T *) Result_array(&npaths,&first_absmq,&second_absmq,result);

    if (failsonlyp == true) {
      /* Skip */
    } else if (quiet_if_excessive_p && npaths > maxpaths_report) {
      Filestring_set_split_output(fp,OUTPUT_UX);
      queryseq1 = Request_queryseq1(request);
      SAM_print_nomapping(fp,ABBREV_UNPAIRED_MULT_XS,
			  queryseq1,/*mate*/NULL,/*acc1*/Shortread_accession(queryseq1),
			  /*acc2*/NULL,chromosome_iit,resulttype,
			  /*first_read_p*/true,/*pathnum*/1,npaths,
			  /*artificial_mate_p*/false,/*npaths_mate*/0,/*mate_chrpos*/0U,
			  quality_shift,sam_read_group_id,invert_first_p,invert_second_p);
      if (failedinput_root != NULL) {
	Shortread_print_query_singleend(*fp_failedinput_1,queryseq1,/*headerseq*/queryseq1);
      }

    } else {
      Filestring_set_split_output(fp,OUTPUT_UM);
      queryseq1 = Request_queryseq1(request);
      for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths_report; pathnum++) {
	stage3 = stage3array[pathnum-1];
	if (Stage3end_hittype(stage3) == SAMECHR_SPLICE || Stage3end_hittype(stage3) == TRANSLOC_SPLICE) {
	  chrpos = 0;
	} else {
	  chrpos = SAM_compute_chrpos(/*hardclip_low*/0,/*hardclip_high*/0,stage3,Shortread_fulllength(queryseq1),
				      /*first_read_p*/true);
	}
	SAM_print(fp,*fp_failedinput_1,ABBREV_UNPAIRED_MULT,
		  stage3,/*mate*/NULL,/*acc1*/Shortread_accession(queryseq1),
		  /*acc2*/NULL,pathnum,npaths,
		  Stage3end_absmq_score(stage3array[pathnum-1]),first_absmq,second_absmq,
		  Stage3end_mapq_score(stage3array[pathnum-1]),
		  chromosome_iit,queryseq1,/*queryseq2*/NULL,
		  /*pairedlength*/0,chrpos,/*mate_chrpos*/0U,
		  /*clipdir*/0,/*hardclip5_low*/0,/*hardclip5_high*/0,/*hardclip3_low*/0,/*hardclip3_high*/0,
		  resulttype,/*first_read_p*/true,/*artificial_mate_p*/false,/*npaths_mate*/0,quality_shift,
		  sam_read_group_id,invert_first_p,invert_second_p,merge_samechr_p);
      }
    }

  } else {
    if (failedinput_root == NULL) {
      *fp_failedinput_2 = (Filestring_T) NULL;
    } else {
      *fp_failedinput_2 = Filestring_new(Request_id(request));
    }
    SAM_print_paired(fp,*fp_failedinput_1,*fp_failedinput_2,result,resulttype,chromosome_iit,
		     Request_queryseq1(request),Request_queryseq2(request),
		     invert_first_p,invert_second_p,nofailsp,failsonlyp,
		     merge_samechr_p,quality_shift,sam_read_group_id);
  }

  return fp;
}


static void
print_header_singleend (Filestring_T fp, Request_T request, bool translocationp, int npaths) {
  Shortread_T queryseq1;

  if (print_m8_p == false) {
    queryseq1 = Request_queryseq1(request);

    FPRINTF(fp,">");
    Shortread_print_oneline(fp,queryseq1);
    FPRINTF(fp,"\t%d",npaths);
    if (translocationp == true) {
      FPRINTF(fp," (transloc)");
    }

    /* No sequence inversion on single-end reads */
    if (Shortread_quality_string(queryseq1) != NULL) {
      FPRINTF(fp,"\t");
      Shortread_print_quality(fp,queryseq1,/*hardclip_low*/0,/*hardclip_high*/0,
			      quality_shift,/*show_chopped_p*/true);
    }

    FPRINTF(fp,"\t");
    Shortread_print_header(fp,queryseq1,/*queryseq2*/NULL);
    /* FPRINTF(fp,"\n"); -- included in header */
  }

  return;
}


/* Taken from print_result_gsnap from old outbuffer.c */
static Filestring_T
filestring_fromresult_gsnap (Filestring_T *fp_failedinput_1, Filestring_T *fp_failedinput_2,
			     Result_T result, Request_T request) {
  Filestring_T fp;
  Resulttype_T resulttype;
  Shortread_T queryseq1, queryseq2;
  Stage3end_T *stage3array, stage3;
  int npaths, pathnum, first_absmq, second_absmq;

  fp = Filestring_new(Request_id(request));
  if (failedinput_root == NULL) {
    *fp_failedinput_1 = (Filestring_T) NULL;
  } else {
    *fp_failedinput_1 = Filestring_new(Request_id(request));
  }

  resulttype = Result_resulttype(result);

  if (resulttype == SINGLEEND_NOMAPPING) {
    *fp_failedinput_2 = (Filestring_T) NULL;
    if (nofailsp == true) {
      /* Skip */
    } else if (print_m8_p) {
      /* Skip */
    } else {
      Filestring_set_split_output(fp,OUTPUT_NM);
      print_header_singleend(fp,request,/*translocationp*/false,/*npaths*/0);
      FPRINTF(fp,"\n");

      if (failedinput_root != NULL) {
	queryseq1 = Request_queryseq1(request);
	Shortread_print_query_singleend(*fp_failedinput_1,queryseq1,/*headerseq*/queryseq1);
      }
    }

  } else if (resulttype == SINGLEEND_UNIQ) {
    *fp_failedinput_2 = (Filestring_T) NULL;
    if (failsonlyp == true) {
      /* Skip */
    } else {
      Filestring_set_split_output(fp,OUTPUT_UU);
      queryseq1 = Request_queryseq1(request);

      stage3array = (Stage3end_T *) Result_array(&npaths,&first_absmq,&second_absmq,result);
      stage3 = stage3array[0];

      print_header_singleend(fp,request,/*translocationp*/false,/*npaths*/1);
      Stage3end_print(fp,stage3,Stage3end_score(stage3),
		      chromosome_iit,queryseq1,/*headerseq*/queryseq1,/*acc_suffix*/"",
		      invert_first_p,/*hit5*/(Stage3end_T) NULL,/*hit3*/(Stage3end_T) NULL,
		      /*pairlength*/0,/*pairscore*/0,/*pairtype*/UNPAIRED,
		      Stage3end_mapq_score(stage3));
      if (print_m8_p == false) {
	FPRINTF(fp,"\n");
      }
    }

  } else if (resulttype == SINGLEEND_TRANSLOC) {
    *fp_failedinput_2 = (Filestring_T) NULL;
    Filestring_set_split_output(fp,OUTPUT_UT);

    stage3array = (Stage3end_T *) Result_array(&npaths,&first_absmq,&second_absmq,result);

    if (failsonlyp == true) {
      /* Skip */

    } else if (quiet_if_excessive_p && npaths > maxpaths_report) {
      print_header_singleend(fp,request,/*translocationp*/true,npaths);
      FPRINTF(fp,"\n");

    } else {
      queryseq1 = Request_queryseq1(request);

      print_header_singleend(fp,request,/*translocationp*/true,npaths);
      for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths_report; pathnum++) {
	stage3 = stage3array[pathnum-1];
	Stage3end_print(fp,stage3,Stage3end_score(stage3),
			chromosome_iit,queryseq1,/*headerseq*/queryseq1,/*acc_suffix*/"",
			invert_first_p,/*hit5*/(Stage3end_T) NULL,/*hit3*/(Stage3end_T) NULL,
			/*pairlength*/0,/*pairscore*/0,/*pairtype*/UNPAIRED,
			Stage3end_mapq_score(stage3));
      }
      if (print_m8_p == false) {
	FPRINTF(fp,"\n");
      }
    }

  } else if (resulttype == SINGLEEND_MULT) {
    *fp_failedinput_2 = (Filestring_T) NULL;
    stage3array = (Stage3end_T *) Result_array(&npaths,&first_absmq,&second_absmq,result);

    if (failsonlyp == true) {
      /* Skip */

    } else if (quiet_if_excessive_p && npaths > maxpaths_report) {
      Filestring_set_split_output(fp,OUTPUT_UX);
      print_header_singleend(fp,request,/*translocationp*/false,npaths);
      if (print_m8_p == false) {
	FPRINTF(fp,"\n");
      }

    } else {
      queryseq1 = Request_queryseq1(request);

      Filestring_set_split_output(fp,OUTPUT_UM);
      print_header_singleend(fp,request,/*translocationp*/false,npaths);
      for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths_report; pathnum++) {
	stage3 = stage3array[pathnum-1];
	Stage3end_print(fp,stage3,Stage3end_score(stage3),
			chromosome_iit,queryseq1,/*headerseq*/queryseq1,/*acc_suffix*/"",
			invert_first_p,/*hit5*/(Stage3end_T) NULL,/*hit3*/(Stage3end_T) NULL,
			/*pairlength*/0,/*pairscore*/0,/*pairtype*/UNPAIRED,
			Stage3end_mapq_score(stage3));
      }
      if (print_m8_p == false) {
	FPRINTF(fp,"\n");
      }
    }

  } else if (resulttype == PAIREDEND_NOMAPPING) {
    if (failedinput_root == NULL) {
      *fp_failedinput_2 = (Filestring_T) NULL;
    } else {
      *fp_failedinput_2 = Filestring_new(Request_id(request));
    }

    if (nofailsp == true) {
      /* No output */

    } else {
      queryseq1 = Request_queryseq1(request);
      queryseq2 = Request_queryseq2(request);
      /* Stage3pair_print_end will call Filestring_set_split_output(), based on resulttype */

      /* First end */
      Stage3pair_print_end(fp,*fp_failedinput_1,result,resulttype,'>',/*firstp*/true,chromosome_iit,
			   /*queryseq*/queryseq1,/*headerseq1*/queryseq1,/*headerseq2*/queryseq2,
			   maxpaths_report,quiet_if_excessive_p,invert_first_p,quality_shift);

      /* Second end */
      Stage3pair_print_end(fp,*fp_failedinput_2,result,resulttype,'<',/*firstp*/false,chromosome_iit,
			   /*queryseq*/queryseq2,/*headerseq1*/queryseq1,/*headerseq2*/queryseq2,
			   maxpaths_report,quiet_if_excessive_p,invert_second_p,quality_shift);

      if (failedinput_root != NULL) {
	Shortread_print_query_pairedend(*fp_failedinput_1,*fp_failedinput_2,queryseq1,queryseq2);
      }
    }

  } else {
    if (failedinput_root == NULL) {
      *fp_failedinput_2 = (Filestring_T) NULL;
    } else {
      *fp_failedinput_2 = Filestring_new(Request_id(request));
    }

    if (failsonlyp == true) {
      /* Unwanted success: skip */
    
    } else {
      queryseq1 = Request_queryseq1(request);
      queryseq2 = Request_queryseq2(request);
      /* Stage3pair_print_end will call Filestring_set_split_output() based on resulttype */

      /* First end */
      Stage3pair_print_end(fp,*fp_failedinput_1,result,resulttype,'>',/*firstp*/true,chromosome_iit,
			   /*queryseq*/queryseq1,/*headerseq1*/queryseq1,/*headerseq2*/queryseq2,
			   maxpaths_report,quiet_if_excessive_p,invert_first_p,quality_shift);

      /* Second end */
      Stage3pair_print_end(fp,*fp_failedinput_2,result,resulttype,'<',/*firstp*/false,chromosome_iit,
			   /*queryseq*/queryseq2,/*headerseq1*/queryseq1,/*headerseq2*/queryseq2,
			   maxpaths_report,quiet_if_excessive_p,invert_second_p,quality_shift);
    }
  }

  return fp;
}

Filestring_T
Output_filestring_fromresult (Filestring_T *fp_failedinput_1, Filestring_T *fp_failedinput_2,
			      Result_T result, Request_T request) {
  if (output_sam_p == true) {
    return filestring_fromresult_sam(&(*fp_failedinput_1),&(*fp_failedinput_2),result,request);
  } else {
    return filestring_fromresult_gsnap(&(*fp_failedinput_1),&(*fp_failedinput_2),result,request);
  }
}

#else
/************************************************************************
 *   Print routines and threads for GMAP
 ************************************************************************/

static void
print_npaths (Filestring_T fp, int npaths, char *chrsubset_name, bool mergedp,
	      Chimera_T chimera, Failure_T failuretype) {

  if (npaths == 0) {
    FPRINTF(fp,"Paths (0):");
  } else if (mergedp == true) {
    FPRINTF(fp,"Paths (1):");
  } else {
    FPRINTF(fp,"Paths (%d):",npaths);
  }
  if (chrsubset_name != NULL) {
    FPRINTF(fp,"  [chrsubset: %s]",chrsubset_name);
  }
  if (failuretype == NO_FAILURE) {
    if (chimera != NULL) {
      Chimera_print(fp,chimera);
    }
  } else if (failuretype == EMPTY_SEQUENCE) {
    FPRINTF(fp," *** Empty sequence ***");
  } else if (failuretype == SHORT_SEQUENCE) {
    FPRINTF(fp," *** Short sequence < index oligo size ***");
  } else if (failuretype == POOR_SEQUENCE) {
    FPRINTF(fp," *** Poor sequence (use -p flag to change pruning behavior) ***");
  } else if (failuretype == REPETITIVE) {
    FPRINTF(fp," *** Repetitive sequence (use -p flag to change pruning behavior) ***");
  }
  FPRINTF(fp,"\n");
  if (npaths == 0) {
    FPRINTF(fp,"\n");
  }
  return;
}


/* Taken from Outbuffer_print_result */
Filestring_T
Output_filestring_fromresult (Filestring_T *fp_failedinput, Result_T result, Request_T request,
			      Sequence_T headerseq) {
  Filestring_T fp;
  char *abbrev;
  Sequence_T queryseq;
  Stage3_T *stage3array;
  int npaths, pathnum, effective_maxpaths, first_absmq, second_absmq;
  Chimera_T chimera = NULL;
  int chimerapos, chimeraequivpos, chimera_cdna_direction;
  int querylength;
  double donor_prob, acceptor_prob;
  bool mergedp = false;
  bool printp = true;

  fp = Filestring_new(Request_id(request));
  if (failedinput_root == NULL) {
    *fp_failedinput = (Filestring_T) NULL;
  } else {
    *fp_failedinput = Filestring_new(Request_id(request));
  }

  queryseq = Request_queryseq(request);
  querylength = Sequence_fulllength_given(queryseq);

  stage3array = Result_array(&npaths,&first_absmq,&second_absmq,result);

  chimerapos = chimeraequivpos = -1;
  chimera_cdna_direction = 0;
  donor_prob = acceptor_prob = 0.0;

  /* Translation */
  if (npaths == 0) {
    Filestring_set_split_output(fp,OUTPUT_NM);
    abbrev = ABBREV_NOMAPPING_1;
    effective_maxpaths = 0;
    if (nofailsp == true) {
      printp = false;
    }

    if (Result_failuretype(result) == POOR_SEQUENCE) {
      fprintf(stderr,"Accession %s skipped (poor sequence).  Use -p flag to change pruning behavior\n",Sequence_accession(headerseq));
    } else if (Result_failuretype(result) == REPETITIVE) {
      fprintf(stderr,"Accession %s skipped (repetitive sequence).  Use -p flag to change pruning behavior\n",Sequence_accession(headerseq));
    } else {
      fprintf(stderr,"No paths found for %s\n",Sequence_accession(headerseq));
    }

  } else if ((mergedp = Result_mergedp(result)) == true) {
    if (Stage3_circularpos(stage3array[0]) > 0) {
      Filestring_set_split_output(fp,OUTPUT_UC);
      abbrev = ABBREV_UNPAIRED_CIRCULAR;
    } else {
      Filestring_set_split_output(fp,OUTPUT_UU);
      abbrev = ABBREV_UNPAIRED_UNIQ;
    }
    effective_maxpaths = 1;
    if (failsonlyp == true) {
      printp = false;
    } else {
      Stage3_translate(stage3array[0],
#ifdef PMAP
		       queryseq,
#endif
		       querylength,fulllengthp,cds_startpos,truncatep,strictp);
    }

  } else if ((chimera = Result_chimera(result)) != NULL) {
    if (chimeras_allowed_p == true) {
      effective_maxpaths = 2;
    } else {
      effective_maxpaths = 0;
    }
    Filestring_set_split_output(fp,OUTPUT_UT);
    abbrev = ABBREV_UNPAIRED_TRANSLOC;

    if (failsonlyp == true) {
      printp = false;
    } else {
      chimerapos = Chimera_pos(chimera);
      chimeraequivpos = Chimera_equivpos(chimera);
      donor_prob = Chimera_donor_prob(chimera);
      acceptor_prob = Chimera_acceptor_prob(chimera);
      chimera_cdna_direction = Chimera_cdna_direction(chimera);

      Stage3_translate_chimera(stage3array[0],stage3array[1],
#ifdef PMAP
			       queryseq,
#endif
			       querylength,fulllengthp,cds_startpos,truncatep,strictp);
    }

  } else if (maxpaths_report == 0) {
    effective_maxpaths = 1;
    if (npaths > 1) {
      Filestring_set_split_output(fp,OUTPUT_UM);
      abbrev = ABBREV_UNPAIRED_MULT;
    } else if (Stage3_circularpos(stage3array[0]) > 0) {
      Filestring_set_split_output(fp,OUTPUT_UC);
      abbrev = ABBREV_UNPAIRED_CIRCULAR;
    } else {
      Filestring_set_split_output(fp,OUTPUT_UU);
      abbrev = ABBREV_UNPAIRED_UNIQ;
    }

    if (failsonlyp == true) {
      printp = false;
    } else if (quiet_if_excessive_p && npaths > maxpaths_report) {
      printp = false;
    } else {
      Stage3_translate(stage3array[0],
#ifdef PMAP
		       queryseq,
#endif
		       querylength,fulllengthp,cds_startpos,truncatep,strictp);
    }

  } else {
    if (npaths > 1) {
      Filestring_set_split_output(fp,OUTPUT_UM);
      abbrev = ABBREV_UNPAIRED_MULT;
    } else if (Stage3_circularpos(stage3array[0]) > 0) {
      Filestring_set_split_output(fp,OUTPUT_UC);
      abbrev = ABBREV_UNPAIRED_CIRCULAR;
    } else {
      Filestring_set_split_output(fp,OUTPUT_UU);
      abbrev = ABBREV_UNPAIRED_UNIQ;
    }

    if (npaths < maxpaths_report) {
      effective_maxpaths = npaths;
    } else {
      effective_maxpaths = maxpaths_report;
    }

    if (failsonlyp == true) {
      printp = false;
    } else if (quiet_if_excessive_p && npaths > maxpaths_report) {
      printp = false;
    } else {
      for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
	Stage3_translate(stage3array[pathnum-1],
#ifdef PMAP
			 queryseq,
#endif
			 querylength,fulllengthp,cds_startpos,truncatep,strictp);
      }
    }
  }

  /* Printing */
  if (printp == false) {
    /* No output, either because of --nofails or --quiet-if-excessive */

  } else {
    if (*fp_failedinput != NULL &&
	(npaths == 0 && quiet_if_excessive_p && npaths > maxpaths_report)) {
      PUTC('>',*fp_failedinput);
      Sequence_print_header(*fp_failedinput,headerseq,checksump);
      Sequence_print(*fp_failedinput,queryseq,/*uppercasep*/false,wraplength,/*trimmedp*/false);
    }

    if (printtype == SIMPLE || printtype == SUMMARY || printtype == ALIGNMENT) {
      /* Print header, even if no alignment is found */
      PUTC('>',fp);
      Sequence_print_header(fp,headerseq,checksump);

      if (npaths == 0) {
	print_npaths(fp,/*npaths*/0,chrsubset_name,
		     /*mergedp*/false,/*chimera*/NULL,Result_failuretype(result));


      } else {
	print_npaths(fp,npaths,chrsubset_name,mergedp,chimera,NO_FAILURE);
	for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
	  Stage3_print_pathsummary(fp,stage3array[pathnum-1],pathnum,
				   chromosome_iit,contig_iit,
				   altstrain_iit,queryseq,dbversion,/*maxmutations*/1000000);
	}
      }

      if (printtype != SIMPLE) {
	FPRINTF(fp,"Alignments:\n");
	for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
	  FPRINTF(fp,"  Alignment for path %d:\n\n",pathnum);
	  Stage3_print_alignment(fp,stage3array[pathnum-1],
				 genome,chromosome_iit,printtype,
				 /*continuousp*/false,/*continuous_by_exon_p*/false,
				 /*flipgenomep*/true,invertmode,nointronlenp,wraplength);
	}
      }

      if (map_iit != NULL) {
	FPRINTF(fp,"Maps:\n");
	for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
	  Stage3_print_map(fp,stage3array[pathnum-1],map_iit,map_divint_crosstable,
			   chromosome_iit,pathnum,map_exons_p,map_bothstrands_p,
			   nflanking,print_comment_p);
	}
      }

    } else if (printtype == COMPRESSED) {
      for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
	Stage3_print_compressed(fp,stage3array[pathnum-1],queryseq,chromosome_iit,
				dbversion,usersegment,pathnum,npaths,
				checksump,chimerapos,chimeraequivpos,
				donor_prob,acceptor_prob,chimera_cdna_direction);
      }

    } else if (printtype == CONTINUOUS) {
      PUTC('>',fp);
      Sequence_print_header(fp,headerseq,checksump);
      if (npaths == 0) {
	FPRINTF(fp,"\n\n\n");
      } else {
	Stage3_print_alignment(fp,stage3array[0],genome,chromosome_iit,printtype,
			       /*continuousp*/true,/*continuous_by_exon_p*/false,
			       /*flipgenomep*/true,invertmode,nointronlenp,wraplength);
      }

    } else if (printtype == CONTINUOUS_BY_EXON) {
      PUTC('>',fp);
      Sequence_print_header(fp,headerseq,checksump);
      print_npaths(fp,npaths,chrsubset_name,mergedp,chimera,NO_FAILURE);
      if (npaths == 0) {
	FPRINTF(fp,"\n\n\n");
      } else {
	Stage3_print_pathsummary(fp,stage3array[0],/*pathnum*/1,
				 chromosome_iit,contig_iit,
				 altstrain_iit,queryseq,
				 dbversion,/*maxmutations*/1000000);
	FPRINTF(fp,"Alignments:\n");
	FPRINTF(fp,"  Alignment for path %d:\n\n",/*pathnum*/1);
	Stage3_print_alignment(fp,stage3array[0],genome,chromosome_iit,printtype,
			       /*continuousp*/false,/*continuous_by_exon_p*/true,
			       /*flipgenomep*/true,invertmode,nointronlenp,wraplength);
      }

    } else if (printtype == EXONS_CDNA) {
      PUTC('>',fp);
      Sequence_print_header(fp,headerseq,checksump);
      for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
	FPRINTF(fp,"<path %d>\n",pathnum);
	Pair_print_exons(fp,Stage3_pairarray(stage3array[0]),Stage3_npairs(stage3array[0]),
			 wraplength,ngap,/*cdna*/true);
	FPRINTF(fp,"</path>\n");
      }

    } else if (printtype == EXONS_GENOMIC) {
      PUTC('>',fp);
      Sequence_print_header(fp,headerseq,checksump);
      for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
	FPRINTF(fp,"<path %d>\n",pathnum);
	Pair_print_exons(fp,Stage3_pairarray(stage3array[0]),Stage3_npairs(stage3array[0]),
			 wraplength,ngap,/*cdna*/false);
	FPRINTF(fp,"</path>\n");
      }

    } else if (printtype == CDNA) {
      for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
	PUTC('>',fp);
	Sequence_print_header(fp,headerseq,checksump);
	Stage3_print_cdna(fp,stage3array[pathnum-1],wraplength);
      }

    } else if (printtype == PROTEIN_GENOMIC) {
      for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
	PUTC('>',fp);
	Sequence_print_header(fp,headerseq,checksump);
	Stage3_print_protein_genomic(fp,stage3array[pathnum-1],wraplength);
      }

    } else if (printtype == PSL_NT) {
      for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
	Stage3_print_pslformat_nt(fp,stage3array[pathnum-1],
				  chromosome_iit,usersegment,queryseq);
      }

#ifdef PMAP
    } else if (printtype == PSL_PRO) {
      for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
	Stage3_print_pslformat_pro(fp,stage3array[pathnum-1],
				   chromosome_iit,usersegment,queryseq,strictp);
      }
#endif

    } else if (printtype == GFF3_GENE || printtype == GFF3_MATCH_CDNA ||
	       printtype == GFF3_MATCH_EST) {
      for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
	Stage3_print_gff3(fp,stage3array[pathnum-1],pathnum,
			  chromosome_iit,usersegment,queryseq,querylength,printtype,
			  /*sourcename*/usersegment ? user_genomicseg : dbversion);
      }

#ifndef PMAP
    } else if (printtype == SAM) {
      if (npaths == 0) {
	Pair_print_sam_nomapping(fp,abbrev,/*acc1*/Sequence_accession(headerseq),/*acc2*/NULL,
				 Sequence_fullpointer(queryseq),Sequence_quality_string(queryseq),
				 Sequence_fulllength(queryseq),quality_shift,
				 Sequence_firstp(queryseq),sam_paired_p,sam_read_group_id);

      } else if (quiet_if_excessive_p && npaths > maxpaths_report) {
	Pair_print_sam_nomapping(fp,abbrev,/*acc1*/Sequence_accession(headerseq),/*acc2*/NULL,
				 Sequence_fullpointer(queryseq),Sequence_quality_string(queryseq),
				 Sequence_fulllength(queryseq),quality_shift,
				 Sequence_firstp(queryseq),sam_paired_p,sam_read_group_id);

      } else if (mergedp == true) {
	Stage3_print_sam(fp,abbrev,stage3array[0],/*pathnum*/1,/*npaths*/1,
			 Stage3_absmq_score(stage3array[0]),first_absmq,second_absmq,
			 Stage3_mapq_score(stage3array[0]),
			 chromosome_iit,usersegment,queryseq,
			 /*chimera_part*/0,/*chimera*/NULL,quality_shift,sam_paired_p,
			 sam_read_group_id);

      } else if (chimera != NULL) {
	Stage3_print_sam(fp,abbrev,stage3array[0],/*pathnum*/1,npaths,
			 Stage3_absmq_score(stage3array[0]),first_absmq,second_absmq,
			 Stage3_mapq_score(stage3array[0]),
			 chromosome_iit,usersegment,queryseq,
			 /*chimera_part*/-1,chimera,quality_shift,sam_paired_p,
			 sam_read_group_id);
	Stage3_print_sam(fp,abbrev,stage3array[1],/*pathnum*/1,npaths,
			 Stage3_absmq_score(stage3array[0]),first_absmq,second_absmq,
			 Stage3_mapq_score(stage3array[0]),
			 chromosome_iit,usersegment,queryseq,
			 /*chimera_part*/+1,chimera,quality_shift,sam_paired_p,
			 sam_read_group_id);

      } else {
	for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
	  Stage3_print_sam(fp,abbrev,stage3array[pathnum-1],pathnum,npaths,
			   Stage3_absmq_score(stage3array[pathnum-1]),first_absmq,second_absmq,
			   Stage3_mapq_score(stage3array[pathnum-1]),
			   chromosome_iit,usersegment,queryseq,
			   /*chimera_part*/0,/*chimera*/NULL,quality_shift,sam_paired_p,
			   sam_read_group_id);
	}
      }
#endif

    } else if (printtype == COORDS) {
      for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
	FPRINTF(fp,">");
	Sequence_print_header(fp,headerseq,checksump);
	Stage3_print_coordinates(fp,stage3array[pathnum-1],chromosome_iit,invertmode);
      }

    } else if (printtype == SPLICESITES) {
      /* Print only best path */
      if (npaths > 0) {
	Stage3_print_splicesites(fp,stage3array[0],chromosome_iit,queryseq);
      }

    } else if (printtype == INTRONS) {
      /* Print only best path */
      if (npaths > 0) {
	Stage3_print_introns(fp,stage3array[0],chromosome_iit,queryseq);
      }

    } else if (printtype == MAP_RANGES) {
      for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
	Stage3_print_iit_map(fp,stage3array[pathnum-1],chromosome_iit,queryseq);
      }
      
    } else if (printtype == MAP_EXONS) {
      for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
	Stage3_print_iit_exon_map(fp,stage3array[pathnum-1],chromosome_iit,queryseq);
      }

    } else {
      fprintf(stderr,"Unexpected printtype %d\n",printtype);
      abort();

    }
  }

  return fp;
}

#endif


