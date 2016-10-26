static char rcsid[] = "$Id: uniqscan.c 173896 2015-09-12 00:11:40Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>		/* Needed to define pthread_t on Solaris */
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For strcpy */
#include <strings.h>		/* For rindex */
#include <ctype.h>
#include <math.h>		/* For rint */

#include "assert.h"
#include "except.h"
#include "mem.h"
#include "bool.h"
#include "types.h"
#include "fopen.h"

#include "mode.h"
#include "sequence.h"
#include "shortread.h"		/* For Shortread_setup */
#include "genome.h"
#include "genome128_hr.h"	/* For Genome_hr_setup */
#include "genome_sites.h"	/* For Genome_sites_setup */
#include "maxent_hr.h"		/* For Maxent_hr_setup */
#include "indexdb_hr.h"
#include "substring.h"
#include "stage3hr.h"
#include "spanningelt.h"
#include "splicestringpool.h"
#include "splicetrie_build.h"
#include "oligo.h"		/* For Oligo_setup */
#include "oligoindex_hr.h"	/* For Oligoindex_hr_setup */
#include "pairpool.h"
#include "diagpool.h"
#include "cellpool.h"
#include "stage2.h"		/* For Stage2_setup */
#include "indel.h"		/* For Indel_setup */
#include "dynprog.h"
#include "dynprog_single.h"
#include "dynprog_genome.h"
#include "dynprog_end.h"
#include "stage1hr.h"
#include "indexdb.h"
#include "intlist.h"
#include "list.h"
#include "listdef.h"
#include "iit-read-univ.h"
#include "iit-read.h"
#include "datadir.h"

#include "stage3.h"		/* To get EXTRAQUERYGAP */

#include "getopt.h"


#define MAX_QUERYLENGTH_FOR_ALLOC    100000
#define MAX_GENOMICLENGTH_FOR_ALLOC 1000000

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/************************************************************************
 *   GMAP parameters
 ************************************************************************/

static int gmap_mode = GMAP_PAIRSEARCH | GMAP_TERMINAL | GMAP_IMPROVEMENT;
static int gmap_min_nconsecutive = 20;
static int nullgap = 600;
static int maxpeelback = 11;
static int maxpeelback_distalmedial = 24;
static int extramaterial_end = 10;
static int extramaterial_paired = 8;
static int max_gmap_pairsearch = 3;
static int max_gmap_terminal = 3;
static int max_gmap_improvement = 3;

static double microexon_spliceprob = 0.90;
static int suboptimal_score_start = -1; /* Determined by simulations to have minimal effect */
static int suboptimal_score_end = 3;	/* Determined by simulations to have diminishing returns above 3 */

static int trigger_score_for_gmap = 5;
static int gmap_allowance = 3;
/* static int trigger_score_for_terminals = 5; -- obsolete */

static int min_intronlength = 9;
static int max_deletionlength = 50;

static int sufflookback = 60;
static int nsufflookback = 5;
static int maxintronlen = 200000; /* Was used previously in stage 1.  Now used only in stage 2 and Stage3_mergeable. */

static int extraband_single = 6; /* This is in addition to length2 -
				    length1.  If onesidegap is true in
				    dynprog.c, then this is equivalent
				    to extraband_single of 0.  Needs
				    to be > 0 to handle default
				    close_indels_mode. */
static int extraband_end = 6; /* Was 6.  Shouldn't differ from 0, since onesidegapp is true?
				 This is only on both sides of main diagonal */
static int extraband_paired = 14; /* This is in addition to length2 - length1 */
static int ngap = 3;


/************************************************************************
 *   Global parameters
 ************************************************************************/

static Univ_IIT_T chromosome_iit = NULL;
static int circular_typeint = -1;
static int nchromosomes = 0;
static bool *circularp = NULL;
static bool any_circular_p;
static Indexdb_T indexdb = NULL;
static Indexdb_T indexdb2 = NULL; /* For cmet or atoi */
static Genome_T genome = NULL;
static Genome_T genomealt = NULL;
static Genomecomp_T *genome_blocks = NULL;


/************************************************************************
 *   Program options
 ************************************************************************/

/* Input options */
static char *user_genomedir = NULL;
static char *dbroot = NULL;
static char *dbversion = NULL;

/* Compute options */
static bool from_right_p = false;

static bool query_unk_mismatch_p = false;
static bool genome_unk_mismatch_p = true;
static bool novelsplicingp = false;

static int trim_mismatch_score = -3;
static int trim_indel_score = -4;

static Access_mode_T offsetsstrm_access = USE_MMAP_ONLY;
static bool expand_offsets_p = false;

#ifdef HAVE_MMAP
static Access_mode_T positions_access = USE_MMAP_ONLY;
static Access_mode_T genome_access = USE_MMAP_ONLY;
#else
static Access_mode_T positions_access = USE_ALLOCATE;
static Access_mode_T genome_access = USE_ALLOCATE;
#endif

static int pairmax;
static int pairmax_dna = 1000;
static int pairmax_rna = 200000;
static int expected_pairlength = 200;
static int pairlength_deviation = 100;


/* static Masktype_T masktype = MASK_REPETITIVE; */
static int subopt_levels = 0;

/* If negative, then hasn't been specified by user.  If between 0 and
   1, then treated as a fraction of the querylength.  Else, treated as
   an integer */
static double user_maxlevel_float = 0.0;
static double user_mincoverage_float = 0.0;

/* Really have only one indel penalty */
static int indel_penalty_middle = 2;
static int indel_penalty_end = 2;

static bool allow_end_indels_p = true;
static int max_middle_insertions = 9;
static int max_middle_deletions = 30;
static int max_end_insertions = 3;
static int max_end_deletions = 6;
static int min_indel_end_matches = 4;
static Chrpos_T shortsplicedist = 200000;
static Chrpos_T shortsplicedist_known;
static Chrpos_T shortsplicedist_novelend = 50000;
static int localsplicing_penalty = 0;
static int distantsplicing_penalty = 100;
static int min_distantsplicing_end_matches = 20;
static double min_distantsplicing_identity = 0.95;
static int min_shortend = 2;
static int antistranded_penalty = 0; /* Most RNA-Seq is non-stranded */

static int index1part;
static int required_index1part = 0;
static int index1interval;
static int required_interval = 0;
static int spansize;
static int indexdb_size_threshold;


/* Genes IIT */
static char *genes_file = (char *) NULL;
static IIT_T genes_iit = NULL;
static int *genes_divint_crosstable = NULL;
static bool favor_multiexon_p = false;


/* Splicing IIT */
static bool knownsplicingp = false;
static bool distances_observed_p = false;
static char *user_splicingdir = (char *) NULL;
static char *splicing_file = (char *) NULL;
static IIT_T splicing_iit = NULL;
static bool amb_closest_p = false;

static int donor_typeint = -1;		/* for splicing_iit */
static int acceptor_typeint = -1;	/* for splicing_iit */

static int *splicing_divint_crosstable = NULL;
static Genomecomp_T *splicecomp = NULL;
static Univcoord_T *splicesites = NULL;
static Splicetype_T *splicetypes = NULL;
static Chrpos_T *splicedists = NULL; /* maximum observed splice distance for given splice site */
static List_T *splicestrings = NULL;
static Genomecomp_T *splicefrags_ref = NULL;
static Genomecomp_T *splicefrags_alt = NULL;
static int nsplicesites = 0;

/* Splicing via splicesites */
static int *nsplicepartners_skip = NULL;
static int *nsplicepartners_obs = NULL;
static int *nsplicepartners_max = NULL;

static bool splicetrie_precompute_p = true;
static Trieoffset_T *trieoffsets_obs = NULL;
static Triecontent_T *triecontents_obs = NULL;
static Trieoffset_T *trieoffsets_max = NULL;
static Triecontent_T *triecontents_max = NULL;


/* Dibase and CMET */
static char *user_cmetdir = NULL;
static char *user_atoidir = NULL;
static Mode_T mode = STANDARD;

/* SNPs IIT */
static char *user_snpsdir = NULL;
static char *snps_root = (char *) NULL;
static IIT_T snps_iit = NULL;
static int *snps_divint_crosstable = NULL;


/* Alignment options */
static bool uncompressedp = false;

/* getopt used alphabetically: AaDdEeGgiJjKklMmNnOoQqstVvwYyZz537 */

static struct option long_options[] = {
  /* Input options */
  {"dir", required_argument, 0, 'D'},	/* user_genomedir */
  {"db", required_argument, 0, 'd'}, /* dbroot */
  {"kmer", required_argument, 0, 'k'}, /* required_index1part, index1part */
  {"sampling", required_argument, 0, 'q'}, /* required_interval, index1interval */
  {"genomefull", no_argument, 0, 'G'}, /* uncompressedp */

  /* Compute options */
  {"from-5-end", no_argument, 0, '5'}, /* from_right_p */
  {"from-3-end", no_argument, 0, '3'}, /* from_right_p */

  {"pairmax-dna", required_argument, 0, 0}, /* pairmax_dna */
  {"pairmax-rna", required_argument, 0, 0}, /* pairmax_rna */

  {"query-unk-mismatch", required_argument, 0, 0}, /* query_unk_mismatch_p */
  {"genome-unk-mismatch", required_argument, 0, 0}, /* genome_unk_mismatch_p */

  {"trim-mismatch-score", required_argument, 0, 0}, /* trim_mismatch_score */
  {"trim-indel-score", required_argument, 0, 0}, /* trim_indel_score */
  {"novelsplicing", required_argument, 0, 'N'}, /* novelsplicingp */

  {"max-mismatches", required_argument, 0, 'm'}, /* user_maxlevel_float */

#if 0
  {"indel-penalty-middle", required_argument, 0, 'i'}, /* indel_penalty_middle */
  {"indel-penalty-end", required_argument, 0, 'I'}, /* indel_penalty_end */
#else
  {"indel-penalty", required_argument, 0, 'i'}, /* indel_penalty_middle, indel_penalty_end */
#endif

  {"indel-endlength", required_argument, 0, 0}, /* min_indel_end_matches, allow_end_indels_p */

  {"max-middle-insertions", required_argument, 0, 'y'}, /* max_middle_insertions */
  {"max-middle-deletions", required_argument, 0, 'z'}, /* max_middle_deletions */
  {"max-end-insertions", required_argument, 0, 'Y'}, /* max_end_insertions */
  {"max-end-deletions", required_argument, 0, 'Z'}, /* max_end_deletions */
  {"suboptimal-levels", required_argument, 0, 'M'}, /* subopt_levels */

  {"localsplicedist", required_argument, 0, 'w'}, /* shortsplicedist */
  {"splicingdir", required_argument, 0, 0},	  /* user_splicingdir */
  {"use-splicing", required_argument, 0, 's'}, /* splicing_iit, knownsplicingp */
  {"genes", required_argument, 0, 'g'}, /* genes_iit */
  {"favor-multiexon", no_argument, 0, 0}, /* favor_multiexon_p */

  {"local-splice-penalty", required_argument, 0, 'e'}, /* localsplicing_penalty */
  {"shortend-splice-endlength", required_argument, 0, 'l'}, /* min_shortend */
  {"antistranded-penalty", required_argument, 0, 0},	    /* antistranded_penalty */

  {"cmetdir", required_argument, 0, 0}, /* user_cmetdir */
  {"atoidir", required_argument, 0, 0}, /* user_atoidir */
  {"mode", required_argument, 0, 0}, /* mode */

  {"snpsdir", required_argument, 0, 'V'},   /* user_snpsdir */
  {"use-snps", required_argument, 0, 'v'}, /* snps_root */

  {"gmap-mode", required_argument, 0, 0}, /* gmap_mode */
  {"trigger-score-for-gmap", required_argument, 0, 0}, /* trigger_score_for_gmap */
  {"gmap-min-match-length", required_argument, 0, 0},      /* gmap_min_nconsecutive */
  {"gmap-allowance", required_argument, 0, 0}, /* gmap_allowance */
  {"max-gmap-pairsearch", required_argument, 0, 0}, /* max_gmap_pairsearch */
  {"max-gmap-terminal", required_argument, 0, 0}, /* max_gmap_terminal */
  {"max-gmap-improvement", required_argument, 0, 0}, /* max_gmap_improvement */
  {"microexon-spliceprob", required_argument, 0, 0}, /* microexon_spliceprob */
  {"stage2-start", required_argument, 0, 0},	     /* suboptimal_score_start */
  {"stage2-end", required_argument, 0, 0},	     /* suboptimal_score_end */

  /* Help options */
  {"version", no_argument, 0, 0}, /* print_program_version */
  {"help", no_argument, 0, 0}, /* print_program_usage */
  {0, 0, 0, 0}
};


static void
print_program_version () {
  fprintf(stdout,"\n");
  fprintf(stdout,"UNIQSCAN\n");
  fprintf(stdout,"Part of GMAP package, version %s\n",PACKAGE_VERSION);
  fprintf(stdout,"Build target: %s\n",TARGET);
  fprintf(stdout,"Features: ");
#ifdef HAVE_PTHREAD
  fprintf(stdout,"pthreads enabled, ");
#else
  fprintf(stdout,"no pthreads, ");
#endif
#ifdef HAVE_ZLIB
  fprintf(stdout,"zlib available, ");
#else
  fprintf(stdout,"no zlib, ");
#endif
#ifdef HAVE_MMAP
  fprintf(stdout,"mmap available, ");
#else
  fprintf(stdout,"no mmap, ");
#endif
#ifdef WORDS_BIGENDIAN
  fprintf(stdout,"bigendian, ");
#else
  fprintf(stdout,"littleendian, ");
#endif
#ifdef HAVE_SIGACTION
  fprintf(stdout,"sigaction available, ");
#else
  fprintf(stdout,"no sigaction, ");
#endif
#ifdef HAVE_64_BIT
  fprintf(stdout,"64 bits available");
#else
  fprintf(stdout,"64 bits not available");
#endif
  fprintf(stdout,"\n");

  fprintf(stdout,"Builtin functions:");
#ifdef HAVE_BUILTIN_CLZ
  fprintf(stdout," clz");
#endif
#ifdef HAVE_BUILTIN_CTZ
  fprintf(stdout," ctz");
#endif
#ifdef HAVE_BUILTIN_POPCOUNT
  fprintf(stdout," popcount");
#endif
  fprintf(stdout,"\n");

  fprintf(stdout,"Sizes: off_t (%d), size_t (%d), unsigned int (%d), long int (%d), long long int (%d)\n",
	  (int) sizeof(off_t),(int) sizeof(size_t),(int) sizeof(unsigned int),(int) sizeof(long int),(int) sizeof(long long int));
  fprintf(stdout,"Default gmap directory: %s\n",GMAPDB);
  fprintf(stdout,"Maximum read length: %d\n",MAX_READLENGTH);
  fprintf(stdout,"Thomas D. Wu, Genentech, Inc.\n");
  fprintf(stdout,"Contact: twu@gene.com\n");
  fprintf(stdout,"\n");
  return;
}

/* This flag is not well-supported, and therefore hidden, but
   kept for backwards compatibility */
/*  -R, --rel=STRING               Release\n\ */


static void
print_program_usage ();


/************************************************************************/



#define POOL_FREE_INTERVAL 200

static char digit[10];

static void
uniqueness_scan (bool from_right_p) {
  char *result, *subsequence, sequence[1024];
  int fulllength, sublength;
  Stage3end_T *stage3array;
  int npaths, first_absmq, second_absmq;

  Floors_T *floors_array;
  Shortread_T queryseq1;
  int i;

  /* For GMAP */
  Oligoindex_array_T oligoindices_major, oligoindices_minor;
  Dynprog_T dynprogL, dynprogM, dynprogR;
  Pairpool_T pairpool;
  Diagpool_T diagpool;
  Cellpool_T cellpool;

  oligoindices_major = Oligoindex_array_new_major(MAX_QUERYLENGTH_FOR_ALLOC,MAX_GENOMICLENGTH_FOR_ALLOC);
  oligoindices_minor = Oligoindex_array_new_minor(MAX_QUERYLENGTH_FOR_ALLOC,MAX_GENOMICLENGTH_FOR_ALLOC);
  dynprogL = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired,
			 /*doublep*/true);
  dynprogM = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired,
			 /*doublep*/false);
  dynprogR = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired,
			 /*doublep*/true);
  pairpool = Pairpool_new();
  diagpool = Diagpool_new();
  cellpool = Cellpool_new();

  floors_array = (Floors_T *) CALLOC(MAX_READLENGTH+1,sizeof(Floors_T));
  /* Except_stack_create(); -- requires pthreads */

  for (i = 0; i < 10; i++) {
    sprintf(&(digit[i]),"%d",i);
  }

  while (fgets(sequence,1024,stdin) != NULL) {
    fulllength = strlen(sequence) - 1; /* Ignore '\n' */
    result = (char *) CALLOC(fulllength+1,sizeof(char));
    for (i = 0; i < fulllength; i++) {
      result[i] = '.';
    }

    /* Handle full sequence */
    queryseq1 = Shortread_new(/*acc*/NULL,/*restofheader*/NULL,/*filterp*/false,sequence,/*long_sequence*/NULL,
			      /*sequence_length*/fulllength,/*quality*/NULL,/*long_quality*/NULL,/*quality_length*/0,
			      /*barcode_length*/0,/*invertp*/0,/*copy_acc_p*/false,/*skipp*/false);
    stage3array = Stage1_single_read(&npaths,&first_absmq,&second_absmq,
				     queryseq1,indexdb,indexdb2,indexdb_size_threshold,
				     floors_array,user_maxlevel_float,user_mincoverage_float,
				     indel_penalty_middle,indel_penalty_end,
				     allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
				     localsplicing_penalty,/*distantsplicing_penalty*/100,min_shortend,
				     oligoindices_major,oligoindices_minor,pairpool,diagpool,cellpool,
				     dynprogL,dynprogM,dynprogR,/*keep_floors_p*/true);

    /* printf("%d: %d\n",sublength,npaths); */
    if (from_right_p == true) {
      if (npaths < 10) {
	result[0] = digit[npaths];
      } else {
	result[0] = '*';
      }
    } else {
      if (npaths < 10) {
	result[fulllength-1] = digit[npaths];
      } else {
	result[fulllength-1] = '*';
      }
    }

    for (i = 0; i < npaths; i++) {
      Stage3end_free(&(stage3array[i]));
    }
    FREE_OUT(stage3array);
    Shortread_free(&queryseq1);
    

    if (npaths < 10) {
      subsequence = (char *) CALLOC(fulllength+1,sizeof(char));
      sublength = index1part+2;
      npaths = 2;
      while (sublength < fulllength && npaths > 1) {
	if (from_right_p == true) {
	  strncpy(subsequence,&(sequence[fulllength-sublength]),sublength);
	} else {
	  strncpy(subsequence,sequence,sublength);
	}
	queryseq1 = Shortread_new(/*acc*/NULL,/*restofheader*/NULL,/*filterp*/false,subsequence,/*long_sequence*/NULL,
				  /*sequence_length*/sublength,/*quality*/NULL,/*long_quality*/NULL,/*quality_length*/0,
				  /*barcode_length*/0,/*invertp*/0,/*copy_acc_p*/false,/*skipp*/false);
	stage3array = Stage1_single_read(&npaths,&first_absmq,&second_absmq,
					 queryseq1,indexdb,indexdb2,indexdb_size_threshold,
					 floors_array,user_maxlevel_float,user_mincoverage_float,
					 indel_penalty_middle,indel_penalty_end,
					 allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
					 localsplicing_penalty,/*distantsplicing_penalty*/100,min_shortend,
					 oligoindices_major,oligoindices_minor,pairpool,diagpool,cellpool,
					 dynprogL,dynprogM,dynprogR,/*keep_floors_p*/true);

	/* printf("%d: %d\n",sublength,npaths); */
	if (from_right_p == true) {
	  if (npaths < 10) {
	    result[fulllength-sublength] = digit[npaths];
	  } else {
	    result[fulllength-sublength] = '*';
	  }
	} else {
	  if (npaths < 10) {
	    result[sublength-1] = digit[npaths];
	  } else {
	    result[sublength-1] = '*';
	  }
	}

	for (i = 0; i < npaths; i++) {
	  Stage3end_free(&(stage3array[i]));
	}
	FREE_OUT(stage3array);
	Shortread_free(&queryseq1);
	sublength++;
      }

      FREE(subsequence);
    }

    printf("%s",sequence);
    printf("%s\n",result);
    FREE(result);

  }

  for (i = 0; i <= MAX_READLENGTH; i++) {
    if (floors_array[i] != NULL) {
      Floors_free_keep(&(floors_array[i]));
    }
  }
  FREE(floors_array);

  Cellpool_free(&cellpool);
  Diagpool_free(&diagpool);
  Pairpool_free(&pairpool);
  Dynprog_free(&dynprogR);
  Dynprog_free(&dynprogM);
  Dynprog_free(&dynprogL);
  Oligoindex_array_free(&oligoindices_minor);
  Oligoindex_array_free(&oligoindices_major);

  return;
}

static int
add_gmap_mode (char *string) {
  if (!strcmp(string,"none")) {
    gmap_mode = 0;
    return 0;
  } else {
    if (!strcmp(string,"improve")) {
      gmap_mode |= GMAP_IMPROVEMENT;
    } else if (!strcmp(string,"terminal")) {
      gmap_mode |= GMAP_TERMINAL;
    } else if (!strcmp(string,"pairsearch")) {
      gmap_mode |= GMAP_PAIRSEARCH;
    } else {
      fprintf(stderr,"Don't recognize gmap-mode type %s\n",string);
      fprintf(stderr,"Allowed values are: none, improve, terminal, pairsearch\n");
      exit(9);
    }
    return 1;
  }
}


static char *
check_valid_int (char *string) {
  char *p = string;

  if (*p == '+' || *p == '-') {
    p++;
  }

  if (!isdigit(*p)) {
    fprintf(stderr,"value %s is not a valid int\n",string);
    exit(9);
    return NULL;
  }
  while (*p != '\0' && isdigit(*p)) {
    p++;
  }

  if (*p == 'e') {
    p++;
    if (*p == '+') {
      p++;
    }
    if (!isdigit(*p)) {
      return false;
    }
    while (*p != '\0' && isdigit(*p)) {
      p++;
    }
  }

  if (*p == '\0') {
    return string;
  } else {
    fprintf(stderr,"value %s is not a valid int\n",string);
    exit(9);
    return NULL;
  }
}


static double
check_valid_float (char *string, const char *option) {
  double value;
  char *p = string;

  if (*p == '+' || *p == '-') {
    p++;
  }

  while (*p != '\0' && isdigit(*p)) {
    p++;
  }
  if (*p == '\0') {
    if ((value = atof(string)) > 1.0 || value < 0.0) {
      fprintf(stderr,"Value for option %s should be between 0.0 and 1.0\n",option);
      exit(9);
    } else {
      return value;
    }
  }

  if (*p == '.') {
    p++;
  }

  if (!isdigit(*p)) {
    fprintf(stderr,"Value %s for option %s is not a valid float\n",string,option);
    exit(9);
    return 0.0;
  }
  while (*p != '\0' && isdigit(*p)) {
    p++;
  }

  if (*p == 'e') {
    p++;
    if (*p == '+' || *p == '-') {
      p++;
    }
    if (!isdigit(*p)) {
      fprintf(stderr,"Value %s for option %s is not a valid float\n",string,option);
      exit(9);
      return 0.0;
    }
    while (*p != '\0' && isdigit(*p)) {
      p++;
    }
  }

  if (*p == '\0') {
    if ((value = atof(string)) > 1.0 || value < 0.0) {
      fprintf(stderr,"Value for option %s should be between 0.0 and 1.0\n",option);
      exit(9);
    } else {
      return value;
    }
  } else {
    fprintf(stderr,"Value %s for option %s is not a valid float\n",string,option);
    exit(9);
    return 0.0;
  }
}


static char *
check_valid_float_or_int (char *string) {
  char *p = string;

  if (*p == '+' || *p == '-') {
    p++;
  }

  while (*p != '\0' && isdigit(*p)) {
    p++;
  }
  if (*p == '\0') {
    return string;
  }

  if (*p == '.') {
    p++;
  }

  if (!isdigit(*p)) {
    fprintf(stderr,"value %s is not a valid float\n",string);
    exit(9);
    return NULL;
  }
  while (*p != '\0' && isdigit(*p)) {
    p++;
  }

  if (*p == 'e') {
    p++;
    if (*p == '+' || *p == '-') {
      p++;
    }
    if (!isdigit(*p)) {
      fprintf(stderr,"value %s is not a valid float\n",string);
      exit(9);
      return NULL;
    }
    while (*p != '\0' && isdigit(*p)) {
      p++;
    }
  }

  if (*p == '\0') {
    return string;
  } else {
    fprintf(stderr,"value %s is not a valid float\n",string);
    exit(9);
    return NULL;
  }
}



int
main (int argc, char *argv[]) {
  char *genomesubdir = NULL, *snpsdir = NULL, *modedir = NULL, *mapdir = NULL, *iitfile = NULL, *fileroot = NULL;

  int opt;
  extern int optind;
  extern char *optarg;
  int long_option_index = 0;
  const char *long_name;
  char *string;

  while ((opt = getopt_long(argc,argv,
			    "D:d:k:q:GN:M:m:i:y:Y:z:Z:w:e:l:g:S:s:V:v:53",
			    long_options, &long_option_index)) != -1) {
    switch (opt) {
    case 0:
      long_name = long_options[long_option_index].name;
      if (!strcmp(long_name,"version")) {
	print_program_version();
	exit(0);
      } else if (!strcmp(long_name,"help")) {
	print_program_usage();
	exit(0);

      } else if (!strcmp(long_name,"mode")) {
	if (!strcmp(optarg,"standard")) {
	  mode = STANDARD;
	} else if (!strcmp(optarg,"cmet-stranded")) {
	  mode = CMET_STRANDED;
	} else if (!strcmp(optarg,"cmet-nonstranded")) {
	  mode = CMET_NONSTRANDED;
	} else if (!strcmp(optarg,"atoi-stranded")) {
	  mode = ATOI_STRANDED;
	} else if (!strcmp(optarg,"atoi-nonstranded")) {
	  mode = ATOI_NONSTRANDED;
	} else if (!strcmp(optarg,"ttoc-stranded")) {
	  mode = TTOC_STRANDED;
	} else if (!strcmp(optarg,"ttoc-nonstranded")) {
	  mode = TTOC_NONSTRANDED;
	} else {
	  fprintf(stderr,"--mode must be standard, cmet-stranded, cmet-nonstranded, atoi-stranded, atoi-nonstranded, ttoc-stranded, or ttoc-nonstranded\n");
	  exit(9);
	}

      } else if (!strcmp(long_name,"cmetdir")) {
	user_cmetdir = optarg;
      } else if (!strcmp(long_name,"atoidir")) {
	user_atoidir = optarg;

      } else if (!strcmp(long_name,"splicingdir")) {
	user_splicingdir = optarg;

      } else if (!strcmp(long_name,"gmap-mode")) {
	gmap_mode = 0;
	string = strtok(optarg,",");
	if (add_gmap_mode(string) != 0) {
	  while ((string = strtok(NULL,",")) != NULL && add_gmap_mode(string) != 0) {
	  }
	}

      } else if (!strcmp(long_name,"trigger-score-for-gmap")) {
	trigger_score_for_gmap = atoi(check_valid_int(optarg));
      } else if (!strcmp(long_name,"gmap-min-match-length")) {
	gmap_min_nconsecutive = atoi(check_valid_int(optarg));
      } else if (!strcmp(long_name,"gmap-allowance")) {
	gmap_allowance = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"max-gmap-pairsearch")) {
	max_gmap_pairsearch = atoi(check_valid_int(optarg));
      } else if (!strcmp(long_name,"max-gmap-terminal")) {
	max_gmap_terminal = atoi(check_valid_int(optarg));
      } else if (!strcmp(long_name,"max-gmap-improvement")) {
	max_gmap_improvement = atoi(check_valid_int(optarg));
      } else if (!strcmp(long_name,"microexon-spliceprob")) {
	microexon_spliceprob = check_valid_float(optarg,long_name);
      } else if (!strcmp(long_name,"stage2-start")) {
	suboptimal_score_start = atoi(check_valid_int(optarg));
      } else if (!strcmp(long_name,"stage2-end")) {
	suboptimal_score_end = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"pairmax-dna")) {
	pairmax_dna = atoi(check_valid_int(optarg));
      } else if (!strcmp(long_name,"pairmax-rna")) {
	pairmax_rna = atoi(check_valid_int(optarg));
      } else if (!strcmp(long_name,"indel-endlength")) {
	min_indel_end_matches = atoi(check_valid_int(optarg));
	if (min_indel_end_matches > 14) {
	  allow_end_indels_p = false;
	}

      } else if (!strcmp(long_name,"antistranded-penalty")) {
	antistranded_penalty = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"favor-multiexon")) {
	favor_multiexon_p = true;

      } else if (!strcmp(long_name,"query-unk-mismatch")) {
	if (!strcmp(optarg,"1")) {
	  query_unk_mismatch_p = true;
	} else if (!strcmp(optarg,"0")) {
	  query_unk_mismatch_p = false;
	} else {
	  fprintf(stderr,"--query-unk-mismatch flag must be 0 or 1\n");
	  exit(9);
	}
      } else if (!strcmp(long_name,"genome-unk-mismatch")) {
	if (!strcmp(optarg,"1")) {
	  genome_unk_mismatch_p = true;
	} else if (!strcmp(optarg,"0")) {
	  genome_unk_mismatch_p = false;
	} else {
	  fprintf(stderr,"--genome-unk-mismatch flag must be 0 or 1\n");
	  exit(9);
	}

      } else if (!strcmp(long_name,"trim-mismatch-score")) {
	trim_mismatch_score = atoi(check_valid_int(optarg));
      } else if (!strcmp(long_name,"trim-indel-score")) {
	trim_indel_score = atoi(check_valid_int(optarg));

      } else {
	/* Shouldn't reach here */
	fprintf(stderr,"Don't recognize option %s.  For usage, run 'gsnap --help'",long_name);
	exit(9);
      }
      break;

    case 'D': user_genomedir = optarg; break;
    case 'd': 
      dbroot = (char *) CALLOC(strlen(optarg)+1,sizeof(char));
      strcpy(dbroot,optarg);
      break;
    case 'k':
      required_index1part = atoi(check_valid_int(optarg));
      if (required_index1part > 16) {
	fprintf(stderr,"The only values allowed for -k are 16 or less\n");
	exit(9);
      }
      break;
    case 'q': required_interval = atoi(check_valid_int(optarg)); break;
    case 'G': uncompressedp = true; break;

    case 'N':
      if (!strcmp(optarg,"1")) {
	novelsplicingp = true;
      } else if (!strcmp(optarg,"0")) {
	novelsplicingp = false;
      } else {
	fprintf(stderr,"Novel splicing (-N flag) must be 0 or 1\n");
	exit(9);
      }
      break;

#if 0
    case 'R': 
      if (!strcmp(optarg,"0")) {
	masktype = MASK_NONE;
      } else if (!strcmp(optarg,"1")) {
	masktype = MASK_FREQUENT;
      } else if (!strcmp(optarg,"2")) {
	masktype = MASK_REPETITIVE;
      } else if (!strcmp(optarg,"3")) {
	masktype = MASK_GREEDY_FREQUENT;
      } else if (!strcmp(optarg,"4")) {
	masktype = MASK_GREEDY_REPETITIVE;
      } else {
	fprintf(stderr,"Masking mode %s not recognized.\n",optarg);
	fprintf(stderr,"Mode 0 means no masking, mode 1 masks frequent oligomers;\n");
	fprintf(stderr,"  mode 2 masks frequent and repetitive oligomers;\n");
	fprintf(stderr,"  mode 3 does greedy masking of frequent oligomers,\n");
	fprintf(stderr,"    then no masking if necessary;\n");
	fprintf(stderr,"  mode 4 does greedy masking of frequent and repetitive oligomers,\n");
	fprintf(stderr,"    then no masking if necessary.\n");
	exit(9);
      }
      break;
#endif

    case 'M': subopt_levels = atoi(check_valid_int(optarg)); break;
    case 'm':
      user_maxlevel_float = atof(check_valid_float_or_int(optarg));
      if (user_maxlevel_float > 1.0 && user_maxlevel_float != rint(user_maxlevel_float)) {
	fprintf(stderr,"Cannot specify fractional value %f for --max-mismatches except between 0.0 and 1.0\n",user_maxlevel_float);
	exit(9);
      }
      break;

    case 'i': indel_penalty_middle = indel_penalty_end = atoi(check_valid_int(optarg)); break;

    case 'y': max_middle_insertions = atoi(check_valid_int(optarg)); break;
    case 'Y': max_end_insertions = atoi(check_valid_int(optarg)); break;
    case 'z': max_middle_deletions = atoi(check_valid_int(optarg)); break;
    case 'Z': max_end_deletions = atoi(check_valid_int(optarg)); break;

    case 'w': shortsplicedist = strtoul(optarg,NULL,10); break;

    case 'e': localsplicing_penalty = atoi(check_valid_int(optarg)); break;
    case 'l': min_shortend = atoi(check_valid_int(optarg)); break;

    case 'g': genes_file = optarg; break;
    case 's': splicing_file = optarg; knownsplicingp = true; break;

    case 'V': user_snpsdir = optarg; break;
    case 'v': snps_root = optarg; break;

    case '5': from_right_p = false; break;
    case '3': from_right_p = true; break;

    case '?': fprintf(stderr,"For usage, run 'gsnap --help'\n"); exit(9);
    default: exit(9);
    }
  }
  argc -= optind;
  argv += optind;


  Except_inactivate();

  if (dbroot == NULL) {
    fprintf(stderr,"Need to specify the -d flag.  For usage, run 'gsnap --help'\n");
    /* print_program_usage(); */
    exit(9);
  }

  if (novelsplicingp == true && knownsplicingp == true) {
    pairmax = pairmax_rna;
    shortsplicedist_known = shortsplicedist;

  } else if (knownsplicingp == true) {
    pairmax = pairmax_rna;
    shortsplicedist_known = shortsplicedist;

  } else if (novelsplicingp == true) {
    pairmax = pairmax_rna;
    shortsplicedist_known = 0;

  } else {
    /* Appears to be DNA-Seq */
    pairmax = pairmax_dna;
    shortsplicedist = shortsplicedist_known = 0U;
  }

  if (distantsplicing_penalty < localsplicing_penalty) {
    fprintf(stderr,"The distant splicing penalty %d cannot be less than local splicing penalty %d\n",
	    distantsplicing_penalty,localsplicing_penalty);
    exit(9);
  }

#if 0
  expand_offsets_p = false;
  offsetsstrm_access = USE_MMAP_ONLY;
  positions_access = USE_MMAP_ONLY;
  genome_access = USE_MMAP_ONLY;
#endif


  /* Prepare genomic data */

  genomesubdir = Datadir_find_genomesubdir(&fileroot,&dbversion,user_genomedir,dbroot);

  iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			    strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
  sprintf(iitfile,"%s/%s.chromosome.iit",genomesubdir,fileroot);
  if ((chromosome_iit = Univ_IIT_read(iitfile,/*readonlyp*/true,/*add_iit_p*/false)) == NULL) {
    fprintf(stderr,"IIT file %s is not valid\n",iitfile);
    exit(9);
  } else {
    nchromosomes = Univ_IIT_total_nintervals(chromosome_iit);
    circular_typeint = Univ_IIT_typeint(chromosome_iit,"circular");
    circularp = Univ_IIT_circularp(&any_circular_p,chromosome_iit);
  }
  FREE(iitfile);


  if (snps_root == NULL) {
    genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,/*genometype*/GENOME_OLIGOS,
			uncompressedp,genome_access,/*sharedp*/false);
    if (mode == CMET_STRANDED || mode == CMET_NONSTRANDED) {
      if (user_cmetdir == NULL) {
	modedir = genomesubdir;
      } else {
	modedir = user_cmetdir;
      }

      if ((indexdb = Indexdb_new_genome(&index1part,&index1interval,
					modedir,fileroot,/*idx_filesuffix*/"metct",/*snps_root*/NULL,
					required_index1part,required_interval,
					expand_offsets_p,offsetsstrm_access,positions_access,
					/*sharedp*/false)) == NULL) {
	fprintf(stderr,"Cannot find metct index file.  Need to run cmetindex first\n");
	exit(9);
      }

      if ((indexdb2 = Indexdb_new_genome(&index1part,&index1interval,
					 modedir,fileroot,/*idx_filesuffix*/"metga",/*snps_root*/NULL,
					 required_index1part,required_interval,
					 expand_offsets_p,offsetsstrm_access,positions_access,
					 /*sharedp*/false)) == NULL) {
	fprintf(stderr,"Cannot find metga index file.  Need to run cmetindex first\n");
	exit(9);
      }

    } else if (mode == ATOI_STRANDED || mode == ATOI_NONSTRANDED) {
      if (user_atoidir == NULL) {
	modedir = genomesubdir;
      } else {
	modedir = user_atoidir;
      }

      if ((indexdb = Indexdb_new_genome(&index1part,&index1interval,
					modedir,fileroot,/*idx_filesuffix*/"a2iag",/*snps_root*/NULL,
					required_index1part,required_interval,
					expand_offsets_p,offsetsstrm_access,positions_access,
					/*sharedp*/false)) == NULL) {
	fprintf(stderr,"Cannot find a2iag index file.  Need to run atoiindex first\n");
	exit(9);
      }

      if ((indexdb2 = Indexdb_new_genome(&index1part,&index1interval,
					 modedir,fileroot,/*idx_filesuffix*/"a2itc",/*snps_root*/NULL,
					 required_index1part,required_interval,
					 expand_offsets_p,offsetsstrm_access,positions_access,
					 /*sharedp*/false)) == NULL) {
	fprintf(stderr,"Cannot find a2itc index file.  Need to run atoiindex first\n");
	exit(9);
      }

    } else if (mode == TTOC_STRANDED || mode == TTOC_NONSTRANDED) {
      if (user_atoidir == NULL) {
	modedir = genomesubdir;
      } else {
	modedir = user_atoidir;
      }

      if ((indexdb = Indexdb_new_genome(&index1part,&index1interval,
					modedir,fileroot,/*idx_filesuffix*/"a2itc",/*snps_root*/NULL,
					required_index1part,required_interval,
					expand_offsets_p,offsetsstrm_access,positions_access,
					/*sharedp*/false)) == NULL) {
	fprintf(stderr,"Cannot find a2itc index file.  Need to run atoiindex first\n");
	exit(9);
      }

      if ((indexdb2 = Indexdb_new_genome(&index1part,&index1interval,
					 modedir,fileroot,/*idx_filesuffix*/"a2iag",/*snps_root*/NULL,
					 required_index1part,required_interval,
					 expand_offsets_p,offsetsstrm_access,positions_access,
					 /*sharedp*/false)) == NULL) {
	fprintf(stderr,"Cannot find a2iag index file.  Need to run atoiindex first\n");
	exit(9);
      }

    } else {
      /* Standard behavior */
      if ((indexdb = Indexdb_new_genome(&index1part,&index1interval,
					genomesubdir,fileroot,IDX_FILESUFFIX,/*snps_root*/NULL,
					required_index1part,required_interval,
					expand_offsets_p,offsetsstrm_access,positions_access,
					/*sharedp*/false)) == NULL) {
	fprintf(stderr,"Cannot find offsets file %s.%s*offsets, needed for GSNAP\n",fileroot,IDX_FILESUFFIX);
	exit(9);
      }
      indexdb2 = indexdb;
    }

  } else {
    if (user_snpsdir == NULL) {
      snpsdir = genomesubdir;
      mapdir = Datadir_find_mapdir(/*user_mapdir*/NULL,genomesubdir,fileroot);
    } else {
      snpsdir = user_snpsdir;
      mapdir = user_snpsdir;
    }

    /* SNPs */
    genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,/*genometype*/GENOME_OLIGOS,
			uncompressedp,genome_access,/*sharedp*/false);
    genomealt = Genome_new(snpsdir,fileroot,snps_root,/*genometype*/GENOME_OLIGOS,
			   uncompressedp,genome_access,/*sharedp*/false);

    if (mode == CMET_STRANDED || mode == CMET_NONSTRANDED) {
      if (user_cmetdir == NULL) {
	modedir = snpsdir;
      } else {
	modedir = user_cmetdir;
      }

      if ((indexdb = Indexdb_new_genome(&index1part,&index1interval,
					modedir,fileroot,/*idx_filesuffix*/"metct",snps_root,
					required_index1part,required_interval,
					expand_offsets_p,offsetsstrm_access,positions_access,
					/*sharedp*/false)) == NULL) {
	fprintf(stderr,"Cannot find metct index file.  Need to run cmetindex first\n");
	exit(9);
      }
      if ((indexdb2 = Indexdb_new_genome(&index1part,&index1interval,
					 modedir,fileroot,/*idx_filesuffix*/"metga",snps_root,
					 required_index1part,required_interval,
					 expand_offsets_p,offsetsstrm_access,positions_access,
					 /*sharedp*/false)) == NULL) {
	fprintf(stderr,"Cannot find metga index file.  Need to run cmetindex first\n");
	exit(9);
      }

    } else if (mode == ATOI_STRANDED || mode == ATOI_NONSTRANDED) {
      if (user_atoidir == NULL) {
	modedir = snpsdir;
      } else {
	modedir = user_atoidir;
      }

      if ((indexdb = Indexdb_new_genome(&index1part,&index1interval,
					modedir,fileroot,/*idx_filesuffix*/"a2iag",snps_root,
					required_index1part,required_interval,
					expand_offsets_p,offsetsstrm_access,positions_access,
					/*sharedp*/false)) == NULL) {
	fprintf(stderr,"Cannot find a2iag index file.  Need to run atoiindex first\n");
	exit(9);
      }
      if ((indexdb2 = Indexdb_new_genome(&index1part,&index1interval,
					 modedir,fileroot,/*idx_filesuffix*/"a2itc",snps_root,
					 required_index1part,required_interval,
					 expand_offsets_p,offsetsstrm_access,positions_access,
					 /*sharedp*/false)) == NULL) {
	fprintf(stderr,"Cannot find a2itc index file.  Need to run atoiindex first\n");
	exit(9);
      }

    } else if (mode == TTOC_STRANDED || mode == TTOC_NONSTRANDED) {
      if (user_atoidir == NULL) {
	modedir = snpsdir;
      } else {
	modedir = user_atoidir;
      }

      if ((indexdb = Indexdb_new_genome(&index1part,&index1interval,
					modedir,fileroot,/*idx_filesuffix*/"a2itc",snps_root,
					required_index1part,required_interval,
					expand_offsets_p,offsetsstrm_access,positions_access,
					/*sharedp*/false)) == NULL) {
	fprintf(stderr,"Cannot find a2itc index file.  Need to run atoiindex first\n");
	exit(9);
      }
      if ((indexdb2 = Indexdb_new_genome(&index1part,&index1interval,
					 modedir,fileroot,/*idx_filesuffix*/"a2iag",snps_root,
					 required_index1part,required_interval,
					 expand_offsets_p,offsetsstrm_access,positions_access,
					 /*sharedp*/false)) == NULL) {
	fprintf(stderr,"Cannot find a2iag index file.  Need to run atoiindex first\n");
	exit(9);
      }

    } else {
      indexdb = Indexdb_new_genome(&index1part,&index1interval,
				   snpsdir,fileroot,/*idx_filesuffix*/"ref",snps_root,
				   required_index1part,required_interval,
				   expand_offsets_p,offsetsstrm_access,positions_access,
				   /*sharedp*/false);
      if (indexdb == NULL) {
	fprintf(stderr,"Cannot find snps index file for %s in directory %s\n",snps_root,snpsdir);
	exit(9);
      }
      indexdb2 = indexdb;
    }

    iitfile = (char *) CALLOC(strlen(mapdir)+strlen("/")+
			      strlen(snps_root)+1,sizeof(char));
    sprintf(iitfile,"%s/%s",mapdir,snps_root);
    fprintf(stderr,"Reading SNPs file %s/%s...",mapdir,snps_root);
    if ((snps_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
			     /*divstring*/NULL,/*add_iit_p*/true,/*labels_read_p*/true)) == NULL) {
      fprintf(stderr,"SNPs file %s.iit not found in %s.\n",snps_root,mapdir);
      if (user_snpsdir == NULL) {
	fprintf(stderr,"Available files:\n");
	Datadir_list_directory(stderr,genomesubdir);
	fprintf(stderr,"Either install file %s or specify a full directory path\n",snps_root);
	fprintf(stderr,"using the -D flag to gsnap.\n");
	exit(9);
      }
    }

    snps_divint_crosstable = Univ_IIT_divint_crosstable(chromosome_iit,snps_iit);

    fprintf(stderr,"done\n");
    FREE(iitfile);
    if (user_snpsdir == NULL) {
      FREE(mapdir);
    }
  }
  genome_blocks = Genome_blocks(genome);


  Dynprog_init(/*mode*/STANDARD);
  Compoundpos_init_positions_free(Indexdb_positions_fileio_p(indexdb));
  Spanningelt_init_positions_free(Indexdb_positions_fileio_p(indexdb));
  Stage1_init_positions_free(Indexdb_positions_fileio_p(indexdb));
  indexdb_size_threshold = (int) (10*Indexdb_mean_size(indexdb,mode,index1part));
  debug(printf("Size threshold is %d\n",indexdb_size_threshold));


  FREE(genomesubdir);
  FREE(fileroot);
  FREE(dbroot);


  Genome_setup(genome,genomealt,/*mode*/STANDARD,circular_typeint);
  Genome_hr_setup(Genome_blocks(genome),/*snp_blocks*/genomealt ? Genome_blocks(genomealt) : NULL,
		  query_unk_mismatch_p,genome_unk_mismatch_p,mode);
  Genome_sites_setup(Genome_blocks(genome),/*snp_blocks*/genomealt ? Genome_blocks(genomealt) : NULL);
  Maxent_hr_setup(Genome_blocks(genome),/*snp_blocks*/genomealt ? Genome_blocks(genomealt) : NULL);
  Indexdb_setup(index1part);
  Indexdb_hr_setup(index1part);
  Oligo_setup(index1part);
  Splicetrie_setup(splicecomp,splicesites,splicefrags_ref,splicefrags_alt,
		   trieoffsets_obs,triecontents_obs,trieoffsets_max,triecontents_max,
		   /*snpp*/snps_iit ? true : false,amb_closest_p,/*amb_clip_p*/true,min_shortend);
  spansize = Spanningelt_setup(index1part,index1interval);
  Indel_setup(min_indel_end_matches,indel_penalty_middle);
  Stage1hr_setup(/*use_sarray_p*/false,/*use_only_sarray_p*/false,index1part,index1interval,
		 spansize,chromosome_iit,nchromosomes,
		 genome,genomealt,mode,/*maxpaths_search*/10,
		 splicesites,splicetypes,splicedists,nsplicesites,
		 novelsplicingp,knownsplicingp,/*find_dna_chimeras_p*/false,distances_observed_p,
		 subopt_levels,max_middle_insertions,max_middle_deletions,
		 shortsplicedist,shortsplicedist_known,shortsplicedist_novelend,min_intronlength,
		 min_distantsplicing_end_matches,min_distantsplicing_identity,
		 nullgap,maxpeelback,maxpeelback_distalmedial,
		 extramaterial_end,extramaterial_paired,gmap_mode,
		 trigger_score_for_gmap,gmap_allowance,max_gmap_pairsearch,
		 max_gmap_terminal,max_gmap_improvement,antistranded_penalty);
  Substring_setup(/*print_nsnpdiffs_p*/false,/*print_snplabels_p*/false,
		  /*show_refdiff_p*/false,snps_iit,snps_divint_crosstable,
		  genes_iit,genes_divint_crosstable,
		  splicing_iit,splicing_divint_crosstable,
		  donor_typeint,acceptor_typeint,trim_mismatch_score,
		  novelsplicingp,knownsplicingp,/*output_sam_p*/false,mode,
		  Univ_IIT_genomelength(chromosome_iit,/*with_circular_alias*/false));
  Dynprog_single_setup(/*homopolymerp*/false);
  Dynprog_genome_setup(novelsplicingp,splicing_iit,splicing_divint_crosstable,
		       donor_typeint,acceptor_typeint);
  Dynprog_end_setup(splicesites,splicetypes,splicedists,nsplicesites,
		    trieoffsets_obs,triecontents_obs,trieoffsets_max,triecontents_max);
  Oligoindex_hr_setup(Genome_blocks(genome),/*mode*/STANDARD);
  Stage2_setup(/*splicingp*/novelsplicingp == true || knownsplicingp == true,/*cross_species_p*/false,
	       suboptimal_score_start,suboptimal_score_end,sufflookback,nsufflookback,maxintronlen,
	       mode,/*snps_p*/snps_iit ? true : false);
  Pair_setup(trim_mismatch_score,trim_indel_score,/*gff3_separators_p*/false,/*sam_insert_0M_p*/false,
	     /*force_xs_direction_p*/false,/*md_lowercase_variant_p*/false,
	     /*snps_p*/snps_iit ? true : false,/*print_nsnpdiffs_p*/snps_iit ? true : false,
	     Univ_IIT_genomelength(chromosome_iit,/*with_circular_alias*/false));
  Stage3_setup(/*splicingp*/novelsplicingp == true || knownsplicingp == true,novelsplicingp,
	       /*require_splicedir_p*/false,splicing_iit,splicing_divint_crosstable,
	       donor_typeint,acceptor_typeint,
	       splicesites,min_intronlength,max_deletionlength,/*min_indel_end_matches*/6,
	       maxpeelback_distalmedial,nullgap,extramaterial_end,extramaterial_paired,
	       extraband_single,extraband_end,extraband_paired,
	       ngap,maxintronlen,/*output_sam_p*/false,/*homopolymerp*/false,
	       /*stage3debug*/NO_STAGE3DEBUG);
  Stage3hr_setup(/*invert_first_p*/false,/*invert_second_p*/false,genome,
		 chromosome_iit,nchromosomes,circular_typeint,genes_iit,genes_divint_crosstable,
		 /*tally_iit*/NULL,/*tally_divint_crosstable*/NULL,
		 /*runlength_iit*/NULL,/*runlength_divint_crosstable*/NULL,
		 distances_observed_p,pairmax,expected_pairlength,pairlength_deviation,
		 localsplicing_penalty,indel_penalty_middle,antistranded_penalty,
		 favor_multiexon_p,gmap_min_nconsecutive,index1part,index1interval,
		 novelsplicingp,/*merge_samechr_p*/false,circularp,/*failedinput_root*/NULL,
		 /*print_m8_p*/false,/*want_random_p*/true);

  uniqueness_scan(from_right_p);

  Dynprog_term();

  if (indexdb2 != indexdb) {
    Indexdb_free(&indexdb2);
  }
  if (indexdb != NULL) {
    Indexdb_free(&indexdb);
  }
  if (dbversion != NULL) {
    FREE(dbversion);
  }
  if (genomealt != NULL) {
    Genome_free(&genomealt);
    FREE(splicefrags_alt);	/* If genomealt == NULL, then splicefrags_alt == splicefrags_ref */
  }
  if (genome != NULL) {
    Genome_free(&genome);
  }

  if (nsplicesites > 0) {
    if (splicetrie_precompute_p == true) {
      FREE(triecontents_max);
      FREE(trieoffsets_max);
      FREE(triecontents_obs);
      FREE(trieoffsets_obs);
    } else {
      FREE(nsplicepartners_max);
      FREE(nsplicepartners_obs);
      FREE(nsplicepartners_skip);
      /* Splicestring_gc(splicestrings,nsplicesites); */
      FREE(splicestrings);
    }
    FREE(splicefrags_ref);
    FREE(splicedists);
    FREE(splicetypes);
    FREE(splicecomp);
    FREE(splicesites);
  }

  if (splicing_iit != NULL) {
    FREE(splicing_divint_crosstable);
    IIT_free(&splicing_iit);
  }

  if (genes_iit != NULL) {
    FREE(genes_divint_crosstable);
    IIT_free(&genes_iit);
  }

  if (snps_iit != NULL) {
    FREE(snps_divint_crosstable);
    IIT_free(&snps_iit);
  }

  if (circularp != NULL) {
    FREE(circularp);
  }

  if (chromosome_iit != NULL) {
    Univ_IIT_free(&chromosome_iit);
  }

  return 0;
}


static void
print_program_usage () {
  fprintf(stdout,"\
Usage: gsnap [OPTIONS...] <FASTA file>, or\n\
       cat <FASTA file> | gmap [OPTIONS...]\n\
\n\
");

  /* Input options */
  fprintf(stdout,"Input options (must include -d)\n");
  fprintf(stdout,"\
  -D, --dir=directory            Genome directory\n\
  -d, --db=STRING                Genome database\n\
  -k, --kmer=INT                 kmer size to use in genome database\n\
                                   (allowed values: 12-15).  If not specified, the program will\n\
                                   find the highest available kmer size in the genome database\n\
  -q, --part=INT/INT             Process only the i-th out of every n sequences\n\
                                   e.g., 0/100 or 99/100 (useful for distributing jobs\n\
                                   to a computer farm).\n\
  --input-buffer-size=INT        Size of input buffer (program reads this many sequences\n\
                                   at a time for efficiency) (default 1000)\n\
  --barcode-length=INT           Amount of barcode to remove from start of read\n\
                                   (default 0)\n\
  -o, --orientation=STRING       Orientation of paired-end reads\n\
                                   Allowed values: FR (fwd-rev, or typical Illumina; default),\n\
                                   RF (rev-fwd, for circularized inserts), or FF (fwd-fwd, same strand)\n\
  --fastq-id-start=INT           Starting position of identifier in FASTQ header, space-delimited (>= 1)\n\
  --fastq-id-end=INT             Ending position of identifier in FASTQ header, space-delimited (>= 1)\n\
                                 Examples:\n\
                                   @HWUSI-EAS100R:6:73:941:1973#0/1                           start=1, end=1 (default)\n\
                                   @SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36\n\
                                      start=1, end=1  => identifier is SRR001666.1\n\
                                      start=2, end=2  => identifier is 071112_SLXA-EAS1_s_7:5:1:817:345\n\
                                      start=1, end=2  => identifier is SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345\n\
");
#ifdef HAVE_ZLIB
  fprintf(stdout,"\
  --gunzip                       Uncompress gzipped input files\n\
");
#endif
  fprintf(stdout,"\n");

  /* Computation options */
  fprintf(stdout,"Computation options\n");
  fprintf(stdout,"\
  -5, --from-5-end               Compute successive substrings from 5' end (default)\n\
  -3, --from-3-end               Compute successive substrings from 3' end\n\
  -m, --max-mismatches=FLOAT     Maximum number of mismatches allowed (if not specified, then\n\
                                   defaults to the ultrafast level of ((readlength+2)/12 - 2))\n\
                                   If specified between 0.0 and 1.0, then treated as a fraction\n\
                                   of each read length.  Otherwise, treated as an integral number\n\
                                   of mismatches (including indel and splicing penalties)\n\
                                   For RNA-Seq, you may need to increase this value slightly\n\
                                   to align reads extending past the ends of an exon.\n\
  --query-unk-mismatch=INT       Whether to count unknown (N) characters in the query as a mismatch\n\
                                   (0=no (default), 1=yes)\n\
  --genome-unk-mismatch=INT      Whether to count unknown (N) characters in the genome as a mismatch\n\
                                   (0=no, 1=yes (default))\n\
");

#if 0
  fprintf(stdout,"\
  -i, --indel-penalty-middle=INT Penalty for an indel in middle of read (default 1).\n\
                                   Counts against mismatches allowed.  To find indels, make\n\
                                   indel-penalty less than or equal to max-mismatches\n\
  -I, --indel-penalty-end=INT    Penalty for an indel at end of read (default 2).\n\
                                   Counts against mismatches allowed.  To find indels, make\n\
                                   indel-penalty less than or equal to max-mismatches\n\
");
#else
  fprintf(stdout,"\
  --terminal-threshold=INT       Threshold for searching for a terminal alignment (from one end of the\n\
                                   read to the best possible position at the other end) (default 3).\n\
                                   To turn off terminal alignments, set this to a high value.\n\
  -i, --indel-penalty=INT        Penalty for an indel (default 2).\n\
                                   Counts against mismatches allowed.  To find indels, make\n\
                                   indel-penalty less than or equal to max-mismatches.\n\
                                   A value < 2 can lead to false positives at read ends\n\
");
#endif

#if 0
  /* No longer used */
  fprintf(stdout,"\
  -R, --masking=INT              Masking of frequent/repetitive oligomers to avoid spending time\n\
                                   on non-unique or repetitive reads\n\
                                   0 = no masking (will try to find non-unique or repetitive matches)\n\
                                   1 = mask frequent oligomers\n\
                                   2 = mask frequent and repetitive oligomers (fastest) (default)\n\
                                   3 = greedy frequent: mask frequent oligomers first, then\n\
                                       try no masking if alignments not found\n\
                                   4 = greedy repetitive: mask frequent and repetitive oligomers first, then\n\
                                       try no masking if alignments not found\n\
");
#endif

  fprintf(stdout,"\
  --indel-endlength=INT          Minimum length at end required for indel alignments (default 4)\n\
  -y, --max-middle-insertions=INT  Maximum number of middle insertions allowed (default 9)\n\
  -z, --max-middle-deletions=INT Maximum number of middle deletions allowed (default 30)\n\
  -Y, --max-end-insertions=INT   Maximum number of end insertions allowed (default 3)\n\
  -Z, --max-end-deletions=INT    Maximum number of end deletions allowed (default 6)\n\
  -M, --suboptimal-levels=INT    Report suboptimal hits beyond best hit (default 0)\n\
                                   All hits with best score plus suboptimal-levels are reported\n\
  --trim-mismatch-score=INT      Score to use for mismatches when trimming at ends (default is -3;\n\
                                   to turn off trimming, specify 0)\n\
  --trim-indel-score=INT         Score to use for mismatches when trimming at ends (default is -4;\n\
                                   to turn off trimming, specify 0)\n\
  -V, --snpsdir=STRING           Directory for SNPs index files (created using snpindex) (default is\n\
                                   location of genome index files specified using -D and -d)\n \
  -v, --use-snps=STRING          Use database containing known SNPs (in <STRING>.iit, built\n\
                                   previously using snpindex) for tolerance to SNPs\n\
  -cmetdir=STRING                Directory for methylcytosine index files (created using cmetindex)\n\
                                   (default is location of genome index files specified using -D, -V, and -d)\n\
  -atoidir=STRING                Directory for A-to-I RNA editing index files (created using atoiindex)\n\
                                   (default is location of genome index files specified using -D, -V, and -d)\n\
  --mode=STRING                  Alignment mode: standard (default), cmet-stranded, cmet-nonstranded,\n\
                                    atoi-stranded, or atoi-nonstranded.  Non-standard modes requires you\n\
                                    to have previously run the cmetindex or atoiindex programs on the genome\n\
  --tallydir=STRING              Directory for tally IIT file to resolve concordant multiple results (default is\n\
                                   location of genome index files specified using -D and -d).  Note: can\n\
                                   just give full path name to --use-tally instead.\n\
  --use-tally=STRING             Use this tally IIT file to resolve concordant multiple results\n\
  --runlengthdir=STRING          Directory for runlength IIT file to resolve concordant multiple results (default is\n\
                                   location of genome index files specified using -D and -d).  Note: can\n\
                                   just give full path name to --use-runlength instead.\n\
  --use-runlength=STRING         Use this runlength IIT file to resolve concordant multiple results\n\
");


#if 0
  fprintf(stdout,"\
  -2, --dibase                   Input is 2-base encoded (e.g., SOLiD), with database built\n\
                                   previously using dibaseindex)\n\
");
#endif

  fprintf(stdout,"\n");
  fprintf(stdout,"Options for GMAP alignment within GSNAP\n");
  fprintf(stdout,"\
  --gmap-mode=STRING             Cases to use GMAP for complex alignments containing multiple splices or indels\n\
                                 Allowed values: none, pairsearch, terminal, improve (or multiple,\n\
                                    separated by commas).  Default: pairsearch,terminal,improve\n\
  --trigger-score-for-gmap=INT   Try GMAP pairsearch on nearby genomic regions if best score (the total\n\
                                   of both ends if paired-end) exceeds this value (default 5)\n \
  --max-gmap-pairsearch=INT      Perform GMAP pairsearch on nearby genomic regions up to this many\n\
                                   many candidate ends (default 3).  Requires pairsearch in --gmap-mode\n\
  --max-gmap-terminal=INT        Perform GMAP terminal on nearby genomic regions up to this many\n\
                                   candidate ends (default 3).  Requires terminal in --gmap-mode\n\
  --max-gmap-improvement=INT     Perform GMAP improvement on nearby genomic regions up to this many\n\
                                   candidate ends (default 3).  Requires improve in --gmap-mode\n\
  --microexon-spliceprob=FLOAT   Allow microexons only if one of the splice site probabilities is\n\
                                   greater than this value (default 0.90)\n\
");
  fprintf(stdout,"\n");


  fprintf(stdout,"Genes options for RNA-Seq\n");
  fprintf(stdout,"\
  -g, --genes=STRING             Look for known genes in <STRING>.iit, to be used for resolving\n\
                                   multiple mapping reads.  See README instructions for the correct\n\
                                   formatting of a genes IIT file.\n\
  --favor-multiexon              In resolving multiple mapping reads, overlaps with known\n\
                                   multi-exon genes are favored over those with known single-exon\n\
                                   genes.  This favors spliced genes over psuedogenes.\n\
");
  fprintf(stdout,"\n");


  /* Splicing options */
  fprintf(stdout,"Splicing options for RNA-Seq\n");
  fprintf(stdout,"\
  -N, --novelsplicing=INT              Look for novel splicing (0=no (default), 1=yes)\n\
  -S, --splicesdir=STRING              Directory for splicing involving known sites or known introns,\n\
                                         as specified by the -s or --use-splices flag (default is\n\
                                         directory computed from -D and -d flags).  Note: can\n\
                                         just give full pathname to the -s flag instead.\n\
  -s, --use-splices=STRING             Look for splicing involving known sites or known introns\n\
                                         (in <STRING>.iit), at short or long distances\n\
                                         See README instructions for the distinction between known sites\n\
                                         and known introns\n\
  --novel-doublesplices                Allow GSNAP to look for two splices in a single-end involving novel\n\
                                         splice sites (default is not to allow this).  Caution: this option\n\
                                         can slow down the program considerably.  A better way to detect\n\
                                         double splices is with known splice sites, using the\n\
                                         --splicing option.\n\
  -w, --localsplicedist=INT            Definition of local novel splicing event (default 200000)\n\
  -e, --local-splice-penalty=INT       Penalty for a local splice (default 0).  Counts against mismatches allowed\n\
  -l, --shortend-splice-endlength=INT  Minimum length at end required for short-end spliced alignments (default 2,\n\
                                         but unless known splice sites are provided with the -s flag, GSNAP may still\n\
                                         need the end length to be the value of -k, or kmer size to find a given splice\n\
  --distant-splice-identity=FLOAT      Minimum identity at end required for distant spliced alignments (default 0.95)\n\
  --antistranded-penalty=INT           Penalty for antistranded splicing when using stranded RNA-Seq protocols.\n\
                                         A positive value, such as 1, expects antisense on the first read\n\
                                         and sense on the second read.  Default is 0, which treats sense and antisense\n\
                                         equally well\n\
");
  fprintf(stdout,"\n");


  /* Paired-end options */
  fprintf(stdout,"Options for paired-end reads\n");
  fprintf(stdout,"\
  --pairmax-dna=INT              Max total genomic length for DNA-Seq paired reads, or other reads\n\
                                   without splicing (default 1000).  Used if -N or -s is not specified.\n\
  --pairmax-rna=INT              Max total genomic length for RNA-Seq paired reads, or other reads\n\
                                   that could have a splice (default 200000).  Used if -N or -s is specified.\n\
                                   Should probably match the value for -w, --localsplicedist.\n\
");

  fprintf(stdout,"\n");

  /* Output options */
  fprintf(stdout,"Output options\n");
  fprintf(stdout,"\
  -Q, --quiet-if-excessive       If more than maximum number of paths are found,\n\
                                   then nothing is printed.\n\
  --print-snps                   Print detailed information about SNPs in reads (works only if -v also selected)\n\
                                   (not fully implemented yet)\n\
");

  fprintf(stdout,"\
  -A, --format=STRING            Another format type, other than default.\n\
                                   Currently implemented: sam\n\
                                   Also allowed, but not installed at compile-time: goby\n\
                                   (To install, need to re-compile with appropriate options)\n\
");

  fprintf(stdout,"\n");

  /* SAM options */
  fprintf(stdout,"Options for SAM output\n");
  fprintf(stdout,"\
  --no-sam-headers               Do not print headers beginning with '@'\n\
  --sam-headers-batch=INT        Print headers only for this batch, as specified by -q\n\
  --read-group-id=STRING         Value to put into read-group id (RG-ID) field\n\
  --read-group-name=STRING       Value to put into read-group name (RG-SM) field\n\
  --read-group-library=STRING    Value to put into read-group library (RG-LB) field\n\
  --read-group-platform=STRING   Value to put into read-group library (RG-PL) field\n\
");
  fprintf(stdout,"\n");


  /* Help options */
  fprintf(stdout,"Help options\n");
  fprintf(stdout,"\
  --version                      Show version\n\
  --help                         Show this help message\n\
");
  return;
}

