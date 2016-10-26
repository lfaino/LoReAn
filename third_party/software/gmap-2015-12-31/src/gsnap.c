static char rcsid[] = "$Id: gsnap.c 181923 2016-01-08 00:43:56Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef USE_MPI
#include <mpi.h>
#include "mpidebug.h"
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
#ifdef HAVE_SSE2
#include <emmintrin.h>
#endif
#ifdef HAVE_SSE4_1
#include <smmintrin.h>
#endif
#ifdef HAVE_POPCNT
#include <immintrin.h>
#endif
#if defined(HAVE_MM_POPCNT)
#include <nmmintrin.h>
#endif
#if defined(HAVE_LZCNT) || defined(HAVE_BMI1)
#include <immintrin.h>
#endif

#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif

#ifdef HAVE_ZLIB
#include <zlib.h>
#define GZBUFFER_SIZE 131072
#endif

#ifdef HAVE_BZLIB
#include "bzip2.h"
#endif

#include <signal.h>

#include "assert.h"
#include "except.h"
#include "mem.h"
#include "bool.h"
#include "types.h"
#include "fopen.h"

#include "mode.h"
#include "sequence.h"
#include "shortread.h"		/* For Shortread_setup */
#include "stopwatch.h"
#include "genome.h"
#include "genome128_hr.h"	/* For Genome_hr_setup */
#include "genome_sites.h"	/* For Genome_sites_setup */
#include "maxent_hr.h"		/* For Maxent_hr_setup */
#include "indexdb_hr.h"
#include "mapq.h"
#include "substring.h"
#include "stage3hr.h"
#include "spanningelt.h"
#include "splicestringpool.h"
#include "splicetrie_build.h"
#include "splice.h"		/* For Splice_setup */
#include "oligo.h"		/* For Oligo_setup */
#include "oligoindex_hr.h"	/* For Oligoindex_hr_setup */
#include "pairpool.h"
#include "diagpool.h"
#include "cellpool.h"
#include "stage2.h"		/* For Stage2_setup */
#ifndef LARGE_GENOMES
#include "sarray-read.h"
#endif
#include "indel.h"		/* For Indel_setup */
#include "dynprog.h"
#include "dynprog_single.h"
#include "dynprog_genome.h"
#include "dynprog_end.h"
#include "stage1hr.h"
#include "indexdb.h"
#include "resulthr.h"
#include "request.h"
#include "intlist.h"
#include "list.h"
#include "listdef.h"
#include "iit-read.h"
#include "datadir.h"
#include "samprint.h"		/* For SAM_setup */

#include "filestring.h"
#include "output.h"
#include "inbuffer.h"
#include "outbuffer.h"
#ifdef USE_MPI
#include "master.h"
#endif

#include "stage3.h"		/* To get EXTRAQUERYGAP */
#include "pair.h"
#include "getopt.h"


#define MIN_INDEXDB_SIZE_THRESHOLD 100

#define MAX_QUERYLENGTH_FOR_ALLOC    100000
#define MAX_GENOMICLENGTH_FOR_ALLOC 1000000


/* MPI Processing */
#ifdef DEBUGM
#define debugm(x) x
#else
#define debugm(x)
#endif

/* File open/close.  Want to turn on in shortread.c also. */
#ifdef DEBUGF
#define debugf(x) x
#else
#define debugf(x)
#endif


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/************************************************************************
 *   GMAP parameters
 ************************************************************************/

static int gmap_mode = GMAP_PAIRSEARCH | GMAP_INDEL_KNOWNSPLICE | GMAP_TERMINAL | GMAP_IMPROVEMENT;
static int gmap_min_nconsecutive = 20;
static int nullgap = 600;
static int maxpeelback = 20;	/* Now controlled by defect_rate */
static int maxpeelback_distalmedial = 24;
static int extramaterial_end = 10;
static int extramaterial_paired = 8;
static int sufflookback = 60;
static int nsufflookback = 5;
static int extraband_single = 3;
static int extraband_end = 3;  /* Shouldn't differ from 0, since onesidegapp is true? */
static int extraband_paired = 7;
static int ngap = 3;  /* 0? */

static int max_gmap_pairsearch = 50; /* Will perform GMAP on up to this many hits5 or hits3 */
static int max_gmap_terminal = 50;   /* Will perform GMAP on up to this many terminals5 or terminals3 */
static int max_gmap_improvement = 5;

static double microexon_spliceprob = 0.95;
static int suboptimal_score_start = -1; /* Determined by simulations to have minimal effect */
static int suboptimal_score_end = 3; /* Determined by simulations to have diminishing returns above 3 */

static int trigger_score_for_gmap = 5;
static int gmap_allowance = 3;
/* static int trigger_score_for_terminals = 5; -- obsolete */

static int min_intronlength = 9;
static int max_deletionlength = 50;


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
static Genome_T genomecomp = NULL;
static Genome_T genomecomp_alt = NULL;
static Genome_T genomebits = NULL;
static Genome_T genomebits_alt = NULL;


#ifdef LARGE_GENOMES
static bool use_sarray_p = false;
static bool use_only_sarray_p = false;
#else
static bool use_sarray_p = true; /* if present */
static bool use_only_sarray_p = false;
static Sarray_T sarray_fwd = NULL;
static Sarray_T sarray_rev = NULL;
#endif

#if 0
static char STANDARD_CHARTABLE[4] = {'A','C','G','T'};
static char CMET_FWD_CHARTABLE[4] = {'A','T','G','T'}; /* CT */
static char CMET_REV_CHARTABLE[4] = {'A','C','A','T'}; /* GA */
static char ATOI_FWD_CHARTABLE[4] = {'G','C','G','T'};     /* AG */
static char ATOI_REV_CHARTABLE[4] = {'A','C','G','C'};     /* TC */
#endif


static bool fastq_format_p = false;
static bool want_random_p = true; /* randomize among equivalent scores */
static Stopwatch_T stopwatch = NULL;

/************************************************************************
 *   Program options
 ************************************************************************/

/* Input options */
static char *user_genomedir = NULL;
static char *dbroot = NULL;
static char *dbversion = NULL;
static int part_modulus = 0;
static int part_interval = 1;
static int barcode_length = 0;
static bool invert_first_p = false;
static bool invert_second_p = true;
static int acc_fieldi_start = 0;
static int acc_fieldi_end = 0;
static bool force_single_end_p = false;
static bool filter_chastity_p = false;
static bool allow_paired_end_mismatch_p = false;
static bool filter_if_both_p = false;
static bool gunzip_p = false;
static bool bunzip2_p = false;

/* Compute options */
static bool chop_primers_p = false;
static bool query_unk_mismatch_p = false;
static bool genome_unk_mismatch_p = true;
static bool novelsplicingp = false;
static bool find_dna_chimeras_p = false;

static int trim_mismatch_score = -3;
static int trim_indel_score = -2; /* was -4 */


static bool sharedp = true;
static bool expand_offsets_p = false;

/* Note: sarray aux files (like lcpchilddc) are always allocated */
#ifdef HAVE_MMAP
/* Level 4 is now default */
static Access_mode_T offsetsstrm_access = USE_ALLOCATE;
static Access_mode_T positions_access = USE_ALLOCATE;
static Access_mode_T genome_access = USE_ALLOCATE;
static Access_mode_T sarray_access = USE_MMAP_PRELOAD;
static Access_mode_T lcp_access = USE_MMAP_PRELOAD;
static Access_mode_T guideexc_access = USE_ALLOCATE;
static Access_mode_T indexij_access = USE_ALLOCATE;
#else
static Access_mode_T offsetsstrm_access = USE_ALLOCATE;
static Access_mode_T positions_access = USE_ALLOCATE;
static Access_mode_T genome_access = USE_ALLOCATE;
static Access_mode_T sarray_access = USE_ALLOCATE;
static Access_mode_T lcp_access = USE_ALLOCATE;
static Access_mode_T guideexc_access = USE_ALLOCATE;
static Access_mode_T indexij_access = USE_ALLOCATE;
#endif

static int pairmax;
static int pairmax_dna = 1000;
static int pairmax_rna = 200000;
static int expected_pairlength = 200;
static int pairlength_deviation = 100;

#ifdef USE_MPI
static int nranks, n_slave_ranks, myid, provided;
static int exclude_ranks[1];
static MPI_Comm workers_comm;
static MPI_Group world_group, workers_group;
static int nthreads0;
static bool master_is_worker_p = false; /* default behavior */
#endif

#ifdef HAVE_PTHREAD
static pthread_t output_thread_id, *worker_thread_ids;
#ifdef USE_MPI
static pthread_t write_stdout_thread_id, parser_thread_id, mpi_interface_thread_id;
#endif
static pthread_key_t global_request_key;
static int nthreads = 1;	/* (int) sysconf(_SC_NPROCESSORS_ONLN) */
#else
static int nthreads = 0;	/* (int) sysconf(_SC_NPROCESSORS_ONLN) */
#endif

/* static Masktype_T masktype = MASK_REPETITIVE; */
static int subopt_levels = 0;

/* If negative, then hasn't been specified by user.  If between 0 and
   1, then treated as a fraction of the querylength.  Else, treated as
   an integer */
static double user_maxlevel_float = -1.0;
static double user_mincoverage_float = -1.0;

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
static int distantsplicing_penalty = 1;
static int min_distantsplicing_end_matches = 20;
static double min_distantsplicing_identity = 0.95;
static int min_shortend = 2;
/* static bool find_novel_doublesplices_p = true; */
static int antistranded_penalty = 0; /* Most RNA-Seq is non-stranded */

static Width_T index1part;
static Width_T required_index1part = 0;
static Width_T index1interval;
static Width_T required_index1interval = 0;
static Width_T spansize;
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
static bool amb_clip_p = true;

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


/* Cmet and AtoI */
static bool dibasep = false;
static char *user_cmetdir = NULL;
static char *user_atoidir = NULL;
static Mode_T mode = STANDARD;

/* SNPs IIT */
static char *user_snpsdir = NULL;
static char *snps_root = (char *) NULL;
static IIT_T snps_iit = NULL;
static int *snps_divint_crosstable = NULL;


/* Tally IIT */
static char *user_tallydir = NULL;
static char *tally_root = (char *) NULL;
static IIT_T tally_iit = NULL;
static int *tally_divint_crosstable = NULL;

static char *user_runlengthdir = NULL;
static char *runlength_root = (char *) NULL;
static IIT_T runlength_iit = NULL;
static int *runlength_divint_crosstable = NULL;


/* Output options */
static unsigned int output_buffer_size = 1000;
static bool output_sam_p = false;

/* For Illumina, subtract 64.  For Sanger, subtract 33. */
/* static int quality_score_adj = 64;  -- Stored in mapq.c */

static bool user_quality_score_adj = false;
static bool user_quality_shift = false;
static int quality_shift = 0;   /* For printing, may want -31 */

static bool exception_raise_p = true;
static bool add_paired_nomappers_p = false;
static bool paired_flag_means_concordant_p = false;
static bool quiet_if_excessive_p = false;
static int maxpaths_search = 1000;
static int maxpaths_report = 100;
static bool orderedp = false;
static bool failsonlyp = false;
static bool nofailsp = false;

static bool print_ncolordiffs_p = false;
static bool print_nsnpdiffs_p = false;
static bool print_snplabels_p = false;

static bool show_refdiff_p = false;
static bool clip_overlap_p = false;
static bool merge_overlap_p = false;
static bool merge_samechr_p = false;
static bool print_m8_p = false;

/* SAM */
static int sam_headers_batch = -1;
static bool sam_headers_p = true;
static bool sam_insert_0M_p = false;
static bool sam_multiple_primaries_p = false;
static char *sam_read_group_id = NULL;
static char *sam_read_group_name = NULL;
static char *sam_read_group_library = NULL;
static char *sam_read_group_platform = NULL;
static bool force_xs_direction_p = false;
static bool md_lowercase_variant_p = false;
static bool hide_soft_clips_p = false;


/* Input/output */
static char *split_output_root = NULL;
static char *output_file = NULL;
static char *failedinput_root = NULL;
static bool appendp = false;
static Outbuffer_T outbuffer;
static Inbuffer_T inbuffer;
static unsigned int inbuffer_nspaces = 1000;
static bool timingp = false;
static bool unloadp = false;


/* Alignment options */
static bool uncompressedp = false;

/* getopt used alphabetically: AaBDdEeGgiJjKklMmNnOoQqstVvwYyZz7 */

static struct option long_options[] = {
  /* Input options */
  {"dir", required_argument, 0, 'D'},	/* user_genomedir */
  {"db", required_argument, 0, 'd'}, /* dbroot */
  {"use-sarray", required_argument, 0, 0}, /* use_sarray_p, use_only_sarray_p */
  {"kmer", required_argument, 0, 'k'}, /* required_index1part, index1part */
  {"sampling", required_argument, 0, 0}, /* required_index1interval, index1interval */
  {"genomefull", no_argument, 0, 'G'}, /* uncompressedp */
  {"part", required_argument, 0, 'q'}, /* part_modulus, part_interval */
  {"orientation", required_argument, 0, 0}, /* invert_first_p, invert_second_p */
  {"input-buffer-size", required_argument, 0, 0}, /* inbuffer_nspaces */
  {"barcode-length", required_argument, 0, 0},	  /* barcode_length */
  {"fastq-id-start", required_argument, 0, 0},	  /* acc_fieldi_start */
  {"fastq-id-end", required_argument, 0, 0},	  /* acc_fieldi_end */
  {"force-single-end", no_argument, 0, 0},	  /* force_single_end_p */
  {"filter-chastity", required_argument, 0, 0},	/* filter_chastity_p, filter_if_both_p */
  {"allow-pe-name-mismatch", no_argument, 0, 0}, /* allow_paired_end_mismatch_p */

#ifdef HAVE_ZLIB
  {"gunzip", no_argument, 0, 0}, /* gunzip_p */
#endif

#ifdef HAVE_BZLIB
  {"bunzip2", no_argument, 0, 0}, /* bunzip2_p */
#endif

  /* Compute options */
  {"use-shared-memory", required_argument, 0, 0}, /* sharedp */
#ifdef HAVE_MMAP
  {"batch", required_argument, 0, 'B'}, /* offsetsstrm_access, positions_access, genome_access */
#endif
  {"expand-offsets", required_argument, 0, 0}, /* expand_offsets_p */
  {"pairmax-dna", required_argument, 0, 0}, /* pairmax_dna */
  {"pairmax-rna", required_argument, 0, 0}, /* pairmax_rna */
  {"pairexpect", required_argument, 0, 0},  /* expected_pairlength */
  {"pairdev", required_argument, 0, 0},  /* pairlength_deviation */

  {"nthreads", required_argument, 0, 't'}, /* nthreads */
  {"adapter-strip", required_argument, 0, 'a'},	/* chop_primers_p */

  {"query-unk-mismatch", required_argument, 0, 0}, /* query_unk_mismatch_p */
  {"genome-unk-mismatch", required_argument, 0, 0}, /* genome_unk_mismatch_p */

  {"trim-mismatch-score", required_argument, 0, 0}, /* trim_mismatch_score */
  {"trim-indel-score", required_argument, 0, 0}, /* trim_indel_score */
  {"novelsplicing", required_argument, 0, 'N'}, /* novelsplicingp */
  {"find-dna-chimeras", required_argument, 0, 0}, /* find_dna_chimeras */

  {"max-mismatches", required_argument, 0, 'm'}, /* user_maxlevel_float */
  {"min-coverage", required_argument, 0, 0}, /* user_mincoverage_float */

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
  {"novelend-splicedist", required_argument, 0, 0}, /* shortsplicedist_novelend */
  {"splicingdir", required_argument, 0, 0},	  /* user_splicingdir */
  {"use-splicing", required_argument, 0, 's'}, /* splicing_iit, knownsplicingp, find_dna_chimeras_p */
  {"ambig-splice-noclip", no_argument, 0, 0},  /* amb_clip_p */
  {"genes", required_argument, 0, 'g'}, /* genes_iit */
  {"favor-multiexon", no_argument, 0, 0}, /* favor_multiexon_p */

  {"local-splice-penalty", required_argument, 0, 'e'}, /* localsplicing_penalty */
  {"distant-splice-penalty", required_argument, 0, 'E'}, /* distantsplicing_penalty */
  {"distant-splice-endlength", required_argument, 0, 'K'}, /* min_distantsplicing_end_matches */
  {"shortend-splice-endlength", required_argument, 0, 'l'}, /* min_shortend */
  {"distant-splice-identity", required_argument, 0, 0}, /* min_distantsplicing_identity */
  {"antistranded-penalty", required_argument, 0, 0},	    /* antistranded_penalty */
  {"merge-distant-samechr", no_argument, 0, 0},		    /* merge_samechr_p */

  {"cmetdir", required_argument, 0, 0}, /* user_cmetdir */
  {"atoidir", required_argument, 0, 0}, /* user_atoidir */
  {"mode", required_argument, 0, 0}, /* mode */

  {"snpsdir", required_argument, 0, 'V'},   /* user_snpsdir */
  {"use-snps", required_argument, 0, 'v'}, /* snps_root */

  {"tallydir", required_argument, 0, 0},   /* user_tallydir */
  {"use-tally", required_argument, 0, 0}, /* tally_root */

  {"runlengthdir", required_argument, 0, 0},   /* user_runlengthdir */
  {"use-runlength", required_argument, 0, 0}, /* runlength_root */

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

  /* Output options */
  {"output-buffer-size", required_argument, 0, 0}, /* output_buffer_size */
  {"format", required_argument, 0, 'A'}, /* output_sam_p, print_m8_p */

  {"quality-protocol", required_argument, 0, 0}, /* quality_score_adj, quality_shift */
  {"quality-zero-score", required_argument, 0, 'J'}, /* quality_score_adj */
  {"quality-print-shift", required_argument, 0, 'j'}, /* quality_shift */
  {"sam-headers-batch", required_argument, 0, 0},	/* sam_headers_batch */
  {"no-sam-headers", no_argument, 0, 0},	/* sam_headers_p */
  {"sam-use-0M", no_argument, 0, 0},		/* sam_insert_0M_p */
  {"sam-multiple-primaries", no_argument, 0, 0}, /* sam_multiple_primaries_p */
  {"read-group-id", required_argument, 0, 0},	/* sam_read_group_id */
  {"read-group-name", required_argument, 0, 0},	/* sam_read_group_name */
  {"read-group-library", required_argument, 0, 0},	/* sam_read_group_library */
  {"read-group-platform", required_argument, 0, 0},	/* sam_read_group_platform */
  {"force-xs-dir", no_argument, 0, 0},			/* force_xs_direction_p */
  {"md-lowercase-snp", no_argument, 0, 0},		/* md_lowercase_variant_p */
  {"extend-soft-clips", no_argument, 0, 0},		/* hide_soft_clips_p */

  {"noexceptions", no_argument, 0, '0'}, /* exception_raise_p */
  {"maxsearch", required_argument, 0, 0}, /* maxpaths_search */
  {"npaths", required_argument, 0, 'n'}, /* maxpaths_report */
  {"add-paired-nomappers", no_argument, 0, 0}, /* add_paired_nomappers_p */
  {"paired-flag-means-concordant", required_argument, 0, 0}, /* paired_flag_means_concordant_p */
  {"quiet-if-excessive", no_argument, 0, 'Q'}, /* quiet_if_excessive_p */
  {"ordered", no_argument, 0, 'O'}, /* orderedp */
  {"clip-overlap", no_argument, 0, 0},	     /* clip_overlap_p */
  {"merge-overlap", no_argument, 0, 0},	     /* merge_overlap_p */
  {"show-refdiff", no_argument, 0, 0},	       /* show_refdiff_p */
  {"print-snps", no_argument, 0, 0}, /* print_snplabels_p */
  {"failsonly", no_argument, 0, 0}, /* failsonlyp */
  {"nofails", no_argument, 0, 0}, /* nofailsp */
  {"output-file", required_argument, 0, 'o'}, /* output_file */
  {"split-output", required_argument, 0, 0}, /* split_output_root */
  {"failed-input", required_argument, 0, 0}, /* failed_input_root */
  {"append-output", no_argument, 0, 0},	     /* appendp */

  {"order-among-best", required_argument, 0, 0}, /* want_random_p */

  /* Diagnostic options */
  {"time", no_argument, 0, 0},	/* timingp */
  {"unload", no_argument, 0, 0},	/* unloadp */

  /* Help options */
  {"check", no_argument, 0, 0}, /* check_compiler_assumptions */
  {"version", no_argument, 0, 0}, /* print_program_version */
  {"help", no_argument, 0, 0}, /* print_program_usage */
  {0, 0, 0, 0}
};


static void
print_program_version () {
  char *genomedir;

  fprintf(stdout,"\n");
  fprintf(stdout,"GSNAP: Genomic Short Nucleotide Alignment Program\n");
  fprintf(stdout,"Part of GMAP package, version %s\n",PACKAGE_VERSION);
  fprintf(stdout,"Build target: %s\n",TARGET);
  fprintf(stdout,"Features: ");
#ifdef HAVE_PTHREAD
  fprintf(stdout,"pthreads enabled, ");
#else
  fprintf(stdout,"no pthreads, ");
#endif
#ifdef HAVE_ALLOCA
  fprintf(stdout,"alloca available, ");
#else
  fprintf(stdout,"no alloca, ");
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

  fprintf(stdout,"Popcnt:");
#ifdef HAVE_POPCNT
  fprintf(stdout," popcnt/lzcnt/tzcnt");
#endif
#ifdef HAVE_MM_POPCNT
  fprintf(stdout," mm_popcnt");
#endif
#ifdef HAVE_BUILTIN_POPCOUNT
  fprintf(stdout," builtin_popcount");
#endif
  fprintf(stdout,"\n");

  fprintf(stdout,"Builtin functions:");
#ifdef HAVE_BUILTIN_CLZ
  fprintf(stdout," builtin_clz");
#endif
#ifdef HAVE_BUILTIN_CTZ
  fprintf(stdout," builtin_ctz");
#endif
#ifdef HAVE_BUILTIN_POPCOUNT
  fprintf(stdout," builtin_popcount");
#endif
  fprintf(stdout,"\n");

  fprintf(stdout,"SIMD functions:");
#ifdef HAVE_ALTIVEC
  fprintf(stdout," Altivec");
#endif
#ifdef HAVE_MMX
  fprintf(stdout," MMX");
#endif
#ifdef HAVE_SSE
  fprintf(stdout," SSE");
#endif
#ifdef HAVE_SSE2
  fprintf(stdout," SSE2");
#endif
#ifdef HAVE_SSE3
  fprintf(stdout," SSE3");
#endif
#ifdef HAVE_SSSE3
  fprintf(stdout," SSSE3");
#endif
#ifdef HAVE_SSE4_1
  fprintf(stdout," SSE4.1");
#endif
#ifdef HAVE_SSE4_2
  fprintf(stdout," SSE4.2");
#endif
#ifdef HAVE_AVX
  fprintf(stdout," AVX");
#endif
  fprintf(stdout,"\n");


  fprintf(stdout,"Sizes: off_t (%d), size_t (%d), unsigned int (%d), long int (%d), long long int (%d)\n",
	  (int) sizeof(off_t),(int) sizeof(size_t),(int) sizeof(unsigned int),(int) sizeof(long int),(int) sizeof(long long int));
  fprintf(stdout,"Default gmap directory (compiled): %s\n",GMAPDB);
  genomedir = Datadir_find_genomedir(/*user_genomedir*/NULL);
  fprintf(stdout,"Default gmap directory (environment): %s\n",genomedir);
  FREE(genomedir);
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


static void
check_compiler_assumptions () {
  unsigned int x = rand(), y = rand();
#ifdef HAVE_SSE2
  int z;
  __m128i a;
#ifdef HAVE_SSE4_1
  char negx, negy;
#endif
#endif


  fprintf(stderr,"Checking compiler assumptions for popcnt: ");
  fprintf(stderr,"%08X ",x);
#ifdef HAVE_LZCNT
  fprintf(stderr,"_lzcnt_u32=%d ",_lzcnt_u32(x));
#endif
#ifdef HAVE_BUILTIN_CLZ
  fprintf(stderr,"__builtin_clz=%d ",__builtin_clz(x));
#endif
#ifdef HAVE_BMI1
  fprintf(stderr,"_tzcnt_u32=%d ",_tzcnt_u32(x));
#endif
#ifdef HAVE_BUILTIN_CTZ
  fprintf(stderr,"__builtin_ctz=%d ",__builtin_ctz(x));
#endif

#ifdef HAVE_POPCNT
  fprintf(stderr,"_popcnt32=%d ",_popcnt32(x));
#endif
#if defined(HAVE_MM_POPCNT)
  fprintf(stderr,"_mm_popcnt_u32=%d ",_mm_popcnt_u32(x));
#endif
#if defined(HAVE_BUILTIN_POPCOUNT)
  fprintf(stderr,"__builtin_popcount=%d ",__builtin_popcount(x));
#endif

  fprintf(stderr,"\n");

#ifdef HAVE_SSE2
  /* With -mavx, compiler may use assembly instructions for _mm_set1_epi32 that don't work on non-AVX machines */
  fprintf(stderr,"Checking compiler assumptions for SSE2: ");
  fprintf(stderr,"%08X %08X",x,y);
  a = _mm_xor_si128(_mm_set1_epi32(x),_mm_set1_epi32(y));
  z = _mm_cvtsi128_si32(a);
  fprintf(stderr," xor=%08X\n",z);

#ifdef HAVE_SSE4_1
  if ((negx = (char) x) > 0) {
    negx = -negx;
  }
  if ((negy = (char) y) > 0) {
    negy = -negy;
  }

  fprintf(stderr,"Checking compiler assumptions for SSE4.1: ");
  fprintf(stderr,"%d %d",negx,negy);
  a = _mm_max_epi8(_mm_set1_epi8(negx),_mm_set1_epi8(negy));
  z = _mm_extract_epi8(a,0);
  fprintf(stderr," max=%d => ",z);
  if (negx > negy) {
    if (z == (int) negx) {
      fprintf(stderr,"compiler sign extends\n"); /* technically incorrect, but SIMD procedures behave properly */
    } else {
      fprintf(stderr,"compiler zero extends\n");
    }
  } else {
    if (z == (int) negy) {
      fprintf(stderr,"compiler sign extends\n"); /* technically incorrect, but SIMD procedures behave properly */
    } else {
      fprintf(stderr,"compiler zero extends\n");
    }
  }

#endif

#endif

  fprintf(stderr,"Finished checking compiler assumptions\n");

  return;
}


/************************************************************************/


static Filestring_T
process_request (Filestring_T *fp_failedinput_1, Filestring_T *fp_failedinput_2,
		 double *worker_runtime, Request_T request, Floors_T *floors_array,
		 Oligoindex_array_T oligoindices_major, Oligoindex_array_T oligoindices_minor,
		 Pairpool_T pairpool, Diagpool_T diagpool, Cellpool_T cellpool,
		 Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		 Stopwatch_T worker_stopwatch) {
  Filestring_T fp;
  Result_T result;
  int jobid;
  Shortread_T queryseq1, queryseq2;
  Stage3end_T *stage3array, *stage3array5, *stage3array3;
  Stage3pair_T *stage3pairarray;

  int npaths, npaths5, npaths3, i;
  int first_absmq, second_absmq, first_absmq5, second_absmq5, first_absmq3, second_absmq3;
  Pairtype_T final_pairtype;

  jobid = Request_id(request);
  queryseq1 = Request_queryseq1(request);
  queryseq2 = Request_queryseq2(request);

  Pairpool_reset(pairpool);
  Diagpool_reset(diagpool);
  Cellpool_reset(cellpool);

  /* printf("%s\n",Shortread_accession(queryseq1)); */

  if (worker_stopwatch != NULL) {
    Stopwatch_start(worker_stopwatch);
  }

  if (queryseq2 == NULL) {
    stage3array = Stage1_single_read(&npaths,&first_absmq,&second_absmq,
				     queryseq1,indexdb,indexdb2,indexdb_size_threshold,
				     floors_array,user_maxlevel_float,user_mincoverage_float,
				     indel_penalty_middle,indel_penalty_end,
				     allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
				     localsplicing_penalty,distantsplicing_penalty,min_shortend,
				     oligoindices_major,oligoindices_minor,
				     pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,
				     /*keep_floors_p*/true);

    result = Result_single_read_new(jobid,(void **) stage3array,npaths,first_absmq,second_absmq);
    fp = Output_filestring_fromresult(&(*fp_failedinput_1),&(*fp_failedinput_2),result,request);
    *worker_runtime = worker_stopwatch == NULL ? 0.00 : Stopwatch_stop(worker_stopwatch);
    Result_free(&result);
    return fp;

  } else if ((stage3pairarray = Stage1_paired_read(&npaths,&first_absmq,&second_absmq,&final_pairtype,
						   &stage3array5,&npaths5,&first_absmq5,&second_absmq5,
						   &stage3array3,&npaths3,&first_absmq3,&second_absmq3,
						   queryseq1,queryseq2,indexdb,indexdb2,indexdb_size_threshold,
						   floors_array,user_maxlevel_float,user_mincoverage_float,
						   indel_penalty_middle,indel_penalty_end,
						   allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
						   localsplicing_penalty,distantsplicing_penalty,min_shortend,
						   oligoindices_major,oligoindices_minor,
						   pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,
						   pairmax,/*keep_floors_p*/true)) != NULL) {
    /* Paired or concordant hits found */
    result = Result_paired_read_new(jobid,(void **) stage3pairarray,npaths,first_absmq,second_absmq,
				    final_pairtype);
    fp = Output_filestring_fromresult(&(*fp_failedinput_1),&(*fp_failedinput_2),result,request);
    *worker_runtime = worker_stopwatch == NULL ? 0.00 : Stopwatch_stop(worker_stopwatch);
    Result_free(&result);
    return fp;

  } else if (chop_primers_p == false || Shortread_chop_primers(queryseq1,queryseq2) == false) {
    /* No paired or concordant hits found, and no adapters found */
    /* Report ends as unpaired */
    result = Result_paired_as_singles_new(jobid,(void **) stage3array5,npaths5,first_absmq5,second_absmq5,
					  (void **) stage3array3,npaths3,first_absmq3,second_absmq3);
    fp = Output_filestring_fromresult(&(*fp_failedinput_1),&(*fp_failedinput_2),result,request);
    *worker_runtime = worker_stopwatch == NULL ? 0.00 : Stopwatch_stop(worker_stopwatch);
    Result_free(&result);
    return fp;

  } else {
    /* Try with potential primers chopped.  queryseq1 and queryseq2 altered by Shortread_chop_primers. */
    for (i = 0; i < npaths5; i++) {
      Stage3end_free(&(stage3array5[i]));
    }
    FREE_OUT(stage3array5);

    for (i = 0; i < npaths3; i++) {
      Stage3end_free(&(stage3array3[i]));
    }
    FREE_OUT(stage3array3);

    if ((stage3pairarray = Stage1_paired_read(&npaths,&first_absmq,&second_absmq,&final_pairtype,
					      &stage3array5,&npaths5,&first_absmq5,&second_absmq5,
					      &stage3array3,&npaths3,&first_absmq3,&second_absmq3,
					      queryseq1,queryseq2,indexdb,indexdb2,indexdb_size_threshold,
					      floors_array,user_maxlevel_float,user_mincoverage_float,
					      indel_penalty_middle,indel_penalty_end,
					      allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
					      localsplicing_penalty,distantsplicing_penalty,min_shortend,
					      oligoindices_major,oligoindices_minor,
					      pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,
					      pairmax,/*keep_floors_p*/false)) != NULL) {
      /* Paired or concordant hits found, after chopping adapters */
      result = Result_paired_read_new(jobid,(void **) stage3pairarray,npaths,first_absmq,second_absmq,
				      final_pairtype);

    } else {
      /* No paired or concordant hits found, after chopping adapters */
      result = Result_paired_as_singles_new(jobid,(void **) stage3array5,npaths5,first_absmq5,second_absmq5,
					    (void **) stage3array3,npaths3,first_absmq3,second_absmq3);
    }

    fp = Output_filestring_fromresult(&(*fp_failedinput_1),&(*fp_failedinput_2),result,request);
    *worker_runtime = worker_stopwatch == NULL ? 0.00 : Stopwatch_stop(worker_stopwatch);
    Result_free(&result);
    return fp;
  }
}



#ifdef HAVE_SIGACTION
static const Except_T sigfpe_error = {"SIGFPE--arithmetic exception"};
static const Except_T sigsegv_error = {"SIGSEGV--segmentation violation"};
static const Except_T sigtrap_error = {"SIGTRAP--hardware fault"};
static const Except_T misc_signal_error = {"Miscellaneous signal"};

static void
signal_handler (int sig) {
  Request_T request;
  Shortread_T queryseq1, queryseq2;

  switch (sig) {
  case SIGABRT: fprintf(stderr,"Signal received: SIGABRT\n"); break;
  case SIGFPE: fprintf(stderr,"Signal received: SIGFPE\n"); break;
  case SIGHUP: fprintf(stderr,"Signal received: SIGHUP\n"); break;
  case SIGILL:
    fprintf(stderr,"Signal received: SIGILL\n");
    fprintf(stderr,"An illegal instruction means that this program is being run on a computer\n");
    fprintf(stderr,"  with different features than the computer used to compile the program\n");
    fprintf(stderr,"You may need to re-compile the program with fewer features by doing something like\n");
    fprintf(stderr,"  ./configure --disable-simd\n");
    break;
  case SIGINT: fprintf(stderr,"Signal received: SIGINT\n"); break;
  case SIGPIPE: fprintf(stderr,"Signal received: SIGPIPE\n"); break;
  case SIGQUIT: fprintf(stderr,"Signal received: SIGQUIT\n"); break;
  case SIGSEGV: fprintf(stderr,"Signal received: SIGSEGV\n"); break;
  case SIGSYS: fprintf(stderr,"Signal received: SIGSYS\n"); break;
  case SIGTERM: fprintf(stderr,"Signal received: SIGTERM\n"); break;
  case SIGTRAP: fprintf(stderr,"Signal received: SIGTRAP\n"); break;
  case SIGXCPU: fprintf(stderr,"Signal received: SIGXCPU\n"); break;
  case SIGXFSZ: fprintf(stderr,"Signal received: SIGXFSZ\n"); break;
  }

  Access_emergency_cleanup();

#if 0
  /* Appears to hang */
#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif

#ifdef HAVE_PTHREAD
  request = (Request_T) pthread_getspecific(global_request_key);
  if (request == NULL) {
    /* fprintf(stderr,"Unable to retrieve request for thread\n"); */
  } else {
    queryseq1 = Request_queryseq1(request);
    queryseq2 = Request_queryseq2(request);
    if (queryseq1 == NULL) {
      fprintf(stderr,"Unable to retrieve queryseq for request\n");
    } else {
      fprintf(stderr,"Problem sequence: ");
      fprintf(stderr,"%s (%d bp)\n",Shortread_accession(queryseq1),Shortread_fulllength(queryseq1));
      if (queryseq2 == NULL) {
	Shortread_stderr_query_singleend_fasta(queryseq1,/*headerseq*/queryseq1);
      } else {
	Shortread_stderr_query_pairedend_fasta(queryseq1,queryseq2,invert_first_p,invert_second_p);
      }
    }
  }
#endif

  exit(9);

  return;
}
#endif


/* #define POOL_FREE_INTERVAL 200 */
#define POOL_FREE_INTERVAL 1


static void
single_thread () {
  Floors_T *floors_array;
  Request_T request;
  Filestring_T fp, fp_failedinput_1, fp_failedinput_2;
  Shortread_T queryseq1;
  int i;
  Stopwatch_T worker_stopwatch;

  /* For GMAP */
  Oligoindex_array_T oligoindices_major, oligoindices_minor;
  Dynprog_T dynprogL, dynprogM, dynprogR;
  Pairpool_T pairpool;
  Diagpool_T diagpool;
  Cellpool_T cellpool;
  int jobid = 0;
  double worker_runtime;

#ifdef MEMUSAGE
  long int memusage_constant = 0, memusage;
  char acc[100+1], comma0[20], comma1[20], comma2[20], comma3[20], comma4[20], comma5[20];
#endif

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
  worker_stopwatch = (timingp == true) ? Stopwatch_new() : (Stopwatch_T) NULL;

  floors_array = (Floors_T *) CALLOC(MAX_READLENGTH+1,sizeof(Floors_T));

  /* Except_stack_create(); -- requires pthreads */

#ifdef MEMUSAGE
  memusage_constant += Mem_usage_report_std_heap();
  Genomicpos_commafmt_fill(comma0,memusage_constant);
  Mem_usage_reset_heap_baseline(0);
#endif

  while ((request = Inbuffer_get_request(inbuffer)) != NULL) {
#ifdef USE_MPI
    debug(printf("rank %d, ",myid));
#endif
    debug(printf("single_thread got request %d\n",Request_id(request)));

#ifdef MEMUSAGE
    queryseq1 = Request_queryseq1(request);
    fprintf(stderr,"Single thread starting %s\n",Shortread_accession(queryseq1));
    Mem_usage_reset_stack_max();
    Mem_usage_reset_heap_max();
#endif

    TRY
      fp = process_request(&fp_failedinput_1,&fp_failedinput_2,&worker_runtime,
			   request,floors_array,oligoindices_major,oligoindices_minor,
			   pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,worker_stopwatch);
      if (timingp == true) {
        queryseq1 = Request_queryseq1(request);
        printf("%s\t%.6f\n",Shortread_accession(queryseq1),worker_runtime);
      }

    ELSE
      queryseq1 = Request_queryseq1(request);
      if (queryseq1 == NULL) {
	fprintf(stderr,"NULL");
      } else if (Shortread_accession(queryseq1) == NULL) {
	fprintf(stderr,"unnamed (%d bp)",Shortread_fulllength(queryseq1));
      } else {
	fprintf(stderr,"Problem sequence: ");
	fprintf(stderr,"%s (%d bp)",Shortread_accession(queryseq1),Shortread_fulllength(queryseq1));
      }
      fprintf(stderr,"\n");
      if (Request_queryseq2(request) == NULL) {
	Shortread_stderr_query_singleend_fasta(queryseq1,/*headerseq*/queryseq1);
      } else {
	Shortread_stderr_query_pairedend_fasta(queryseq1,Request_queryseq2(request),
					       invert_first_p,invert_second_p);
      }
      fprintf(stderr,"\n");
      fprintf(stderr,"To obtain a core dump, re-run program on problem sequence with the -0 [zero] flag\n");

      fprintf(stderr,"Exiting...\n");
      exit(9);
    RERAISE;
    END_TRY;

    Outbuffer_print_filestrings(fp,fp_failedinput_1,fp_failedinput_2);

    if (jobid % POOL_FREE_INTERVAL == 0) {
      Pairpool_free_memory(pairpool);
      Diagpool_free_memory(diagpool);
      Cellpool_free_memory(cellpool);
    }

#ifdef MEMUSAGE
    /* Copy acc before we free the request */
    queryseq1 = Request_queryseq1(request);
    strncpy(acc,Shortread_accession(queryseq1),100);
    acc[100] = '\0';
#endif

    Request_free(&request);

#ifdef MEMUSAGE
    Genomicpos_commafmt_fill(comma1,Mem_usage_report_std_heap_max());
    Genomicpos_commafmt_fill(comma2,Mem_usage_report_std_heap());
    Genomicpos_commafmt_fill(comma3,Mem_usage_report_keep());
    Genomicpos_commafmt_fill(comma4,Mem_usage_report_in());
    Genomicpos_commafmt_fill(comma5,Mem_usage_report_out());

    fprintf(stderr,"Acc %s: constant %s  max %s  std %s  keep %s  in %s  out %s\n",
	    acc,comma0,comma1,comma2,comma3,comma4,comma5);

    if ((memusage = Mem_usage_report_std_heap()) != 0) {
      fprintf(stderr,"Memory leak in single thread of %ld bytes\n",memusage);
      fflush(stdout);
      exit(9);
    }
#endif
  }

#ifdef MEMUSAGE
  Mem_usage_std_heap_add(memusage_constant);
#endif

  /* Except_stack_destroy(); -- requires pthreads */

  for (i = 0; i <= MAX_READLENGTH; i++) {
    if (floors_array[i] != NULL) {
      Floors_free_keep(&(floors_array[i]));
    }
  }
  FREE(floors_array);

  if (worker_stopwatch != NULL) {
    Stopwatch_free(&worker_stopwatch);
  }
  Cellpool_free(&cellpool);
  Diagpool_free(&diagpool);
  Pairpool_free(&pairpool);
  Dynprog_free(&dynprogR);
  Dynprog_free(&dynprogM);
  Dynprog_free(&dynprogL);
  Oligoindex_array_free(&oligoindices_minor);
  Oligoindex_array_free(&oligoindices_major);

#ifdef MEMUSAGE
  Mem_usage_set_threadname("main");
#endif

  return;
}


#ifdef HAVE_PTHREAD
static void *
worker_thread (void *data) {
  Floors_T *floors_array;
  Request_T request;
  Filestring_T fp, fp_failedinput_1, fp_failedinput_2;
  Shortread_T queryseq1;
  int i;
  Stopwatch_T worker_stopwatch;

  /* For GMAP */
  Oligoindex_array_T oligoindices_major, oligoindices_minor;
  Dynprog_T dynprogL, dynprogM, dynprogR;
  Pairpool_T pairpool;
  Diagpool_T diagpool;
  Cellpool_T cellpool;
  int worker_jobid = 0;
  double worker_runtime;
#if defined(DEBUG) || defined(MEMUSAGE)
  long int worker_id = (long int) data;
#endif

#ifdef MEMUSAGE
  long int memusage_constant = 0, memusage;
  char threadname[12];
  char acc[100+1], comma0[20], comma1[20], comma2[20], comma3[20], comma4[20], comma5[20];
  sprintf(threadname,"thread-%ld",worker_id);
  Mem_usage_set_threadname(threadname);
#endif

#ifdef USE_MPI
  debug(fprintf(stderr,"rank %d, ",myid));
#endif
  debug(fprintf(stderr,"worker_thread %ld starting\n",worker_id));

  /* Thread-specific data and storage */
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
  worker_stopwatch = (timingp == true) ? Stopwatch_new() : (Stopwatch_T) NULL;

  floors_array = (Floors_T *) CALLOC(MAX_READLENGTH+1,sizeof(Floors_T));

  Except_stack_create();

#ifdef MEMUSAGE
  memusage_constant += Mem_usage_report_std_heap();
  Genomicpos_commafmt_fill(comma0,memusage_constant);
  Mem_usage_reset_heap_baseline(0);
#endif

  while ((request = Inbuffer_get_request(inbuffer)) != NULL) {
#ifdef USE_MPI
    debug(fprintf(stderr,"rank %d, ",myid));
#endif
    debug(fprintf(stderr,"worker_thread %ld got request %d (%s)\n",
		  worker_id,Request_id(request),Shortread_accession(Request_queryseq1(request))));
    pthread_setspecific(global_request_key,(void *) request);

#ifdef MEMUSAGE
    queryseq1 = Request_queryseq1(request);
    fprintf(stderr,"Thread %d starting %s\n",worker_id,Shortread_accession(queryseq1));
    Mem_usage_reset_stack_max();
    Mem_usage_reset_heap_max();
#endif

    TRY
      fp = process_request(&fp_failedinput_1,&fp_failedinput_2,&worker_runtime,
			   request,floors_array,oligoindices_major,oligoindices_minor,
			   pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,worker_stopwatch);
      if (timingp == true) {
        queryseq1 = Request_queryseq1(request);
        printf("%s\t%.6f\n",Shortread_accession(queryseq1),worker_runtime);
      }

    ELSE
      queryseq1 = Request_queryseq1(request);
      if (queryseq1 == NULL) {
	fprintf(stderr,"NULL");
      } else if (Shortread_accession(queryseq1) == NULL) {
	fprintf(stderr,"unnamed (%d bp)",Shortread_fulllength(queryseq1));
      } else {
	fprintf(stderr,"Problem sequence: ");
	fprintf(stderr,"%s (%d bp)",Shortread_accession(queryseq1),Shortread_fulllength(queryseq1));
      }
      fprintf(stderr,"\n");
      if (Request_queryseq2(request) == NULL) {
	Shortread_stderr_query_singleend_fasta(queryseq1,/*headerseq*/queryseq1);
      } else {
	Shortread_stderr_query_pairedend_fasta(queryseq1,Request_queryseq2(request),
					       invert_first_p,invert_second_p);
      }
      fprintf(stderr,"\n");
      fprintf(stderr,"To obtain a core dump, re-run program on problem sequence with the -0 [zero] flag\n");

      fprintf(stderr,"Exiting...\n");
      exit(9);
    RERAISE;
    END_TRY;

    Outbuffer_put_filestrings(outbuffer,fp,fp_failedinput_1,fp_failedinput_2);

    if (worker_jobid % POOL_FREE_INTERVAL == 0) {
      Pairpool_free_memory(pairpool);
      Diagpool_free_memory(diagpool);
      Cellpool_free_memory(cellpool);
    }

#ifdef MEMUSAGE
    /* Copy acc before we free the request */
    queryseq1 = Request_queryseq1(request);
    strncpy(acc,Shortread_accession(queryseq1),100);
    acc[100] = '\0';
#endif

    Request_free(&request);

#ifdef MEMUSAGE
    Genomicpos_commafmt_fill(comma1,Mem_usage_report_std_heap_max());
    Genomicpos_commafmt_fill(comma2,Mem_usage_report_std_heap());
    Genomicpos_commafmt_fill(comma3,Mem_usage_report_keep());
    Genomicpos_commafmt_fill(comma4,Mem_usage_report_in());
    Genomicpos_commafmt_fill(comma5,Mem_usage_report_out());

    fprintf(stderr,"Acc %s, thread %d: constant %s  max %s  std %s  keep %s  in %s  out %s\n",
	    acc,worker_id,comma0,comma1,comma2,comma3,comma4,comma5);

    if ((memusage = Mem_usage_report_std_heap()) != 0) {
      fprintf(stderr,"Memory leak in worker thread %ld of %ld bytes\n",worker_id,memusage);
      fflush(stdout);
      exit(9);
    }
#endif
  }

#ifdef MEMUSAGE
  Mem_usage_std_heap_add(memusage_constant);
#endif

  Except_stack_destroy();

  for (i = 0; i <= MAX_READLENGTH; i++) {
    if (floors_array[i] != NULL) {
      Floors_free_keep(&(floors_array[i]));
    }
  }
  FREE(floors_array);

  if (worker_stopwatch != NULL) {
    Stopwatch_free(&worker_stopwatch);
  }
  Cellpool_free(&cellpool);
  Diagpool_free(&diagpool);
  Pairpool_free(&pairpool);
  Dynprog_free(&dynprogR);
  Dynprog_free(&dynprogM);
  Dynprog_free(&dynprogL);
  Oligoindex_array_free(&oligoindices_minor);
  Oligoindex_array_free(&oligoindices_major);

#ifdef MEMUSAGE
  Mem_usage_set_threadname("main");
#endif

#ifdef USE_MPI
  debug(fprintf(stderr,"rank %d, ",myid));
#endif
  debug(fprintf(stderr,"worker_thread %ld finished\n",worker_id));

  return (void *) NULL;
}
#endif


static void
parse_part (int *part_modulus, int *part_interval, char *string) {
  char *p = string;

  if (sscanf(p,"%d",&(*part_modulus)) < 1) {
    fprintf(stderr,"Cannot parse first integer from %s\n",string);
    exit(9);
  }

  while (*p != '\0' && isdigit(*p)) {
    p++;
  }
  while (*p != '\0' && !isdigit(*p)) {
    p++;
  }
  if (sscanf(p,"%d",&(*part_interval)) < 1) {
    fprintf(stderr,"Cannot parse first integer from %s\n",string);
    exit(9);
  }
  if ((*part_modulus) >= (*part_interval)) {
    fprintf(stderr,"In %s, batch number %d must be less than the number of batches %d\n",
	    string,*part_modulus,*part_interval);
    exit(9);
  }
  if (*part_interval == 0) {
    fprintf(stderr,"Bad batch specification %s.  Batch interval cannot be 0.\n",string);
    exit(9);
  }

  return;
}


static int
add_gmap_mode (char *string) {
  if (!strcmp(string,"none")) {
    gmap_mode = 0;
    return 0;
  } else if (!strcmp(string,"all")) {
    gmap_mode = (GMAP_IMPROVEMENT | GMAP_TERMINAL | GMAP_INDEL_KNOWNSPLICE | GMAP_PAIRSEARCH);
    return 1;
  } else {
    if (!strcmp(string,"improve")) {
      gmap_mode |= GMAP_IMPROVEMENT;
    } else if (!strcmp(string,"terminal")) {
      gmap_mode |= GMAP_TERMINAL;
    } else if (!strcmp(string,"indel_knownsplice")) {
      gmap_mode |= GMAP_INDEL_KNOWNSPLICE;
    } else if (!strcmp(string,"pairsearch")) {
      gmap_mode |= GMAP_PAIRSEARCH;
    } else {
      fprintf(stderr,"Don't recognize gmap-mode type %s\n",string);
      fprintf(stderr,"Allowed values are: none, all, improve, terminal, indel_knownsplice, pairsearch\n");
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


static int
parse_command_line (int argc, char *argv[], int optind) {
  int opt, c;
  extern char *optarg;
  int long_option_index = 0;
  const char *long_name;
  char **argstart;
  char *string;

  fprintf(stderr,"GSNAP version %s called with args:",PACKAGE_VERSION);
  argstart = &(argv[-optind]);
  for (c = 1; c < argc + optind; c++) {
      fprintf(stderr," %s",argstart[c]);
  }
  fprintf(stderr,"\n");

  while ((opt = getopt_long(argc,argv,
			    "D:d:k:Gq:o:a:N:M:m:i:y:Y:z:Z:w:E:e:J:K:l:g:s:V:v:B:t:A:j:0n:QO",
			    long_options, &long_option_index)) != -1) {
    switch (opt) {
    case 0:
      long_name = long_options[long_option_index].name;
      if (!strcmp(long_name,"version")) {
	print_program_version();
	return 1;
      } else if (!strcmp(long_name,"check")) {
	check_compiler_assumptions();
	return 1;
      } else if (!strcmp(long_name,"help")) {
	print_program_usage();
	return 1;

#ifdef LARGE_GENOMES
      } else if (!strcmp(long_name,"use-sarray")) {
	if (!strcmp(optarg,"0")) {
	  use_sarray_p = false;
	  use_only_sarray_p = false;
	} else {
	  fprintf(stderr,"--use-sarray flag for large genomes must be 0\n");
	  return 9;
	}

#else
      } else if (!strcmp(long_name,"use-sarray")) {
	if (!strcmp(optarg,"2")) {
	  use_sarray_p = true;
	  use_only_sarray_p = true;
	} else if (!strcmp(optarg,"1")) {
	  use_sarray_p = true;
	  use_only_sarray_p = false;
	} else if (!strcmp(optarg,"0")) {
	  use_sarray_p = false;
	  use_only_sarray_p = false;
	} else {
	  fprintf(stderr,"--use-sarray flag must be 0, 1, or 2\n");
	  return 9;
	}
#endif

      } else if (!strcmp(long_name,"use-shared-memory")) {
	if (!strcmp(optarg,"1")) {
	  sharedp = true;
	} else if (!strcmp(optarg,"0")) {
	  sharedp = false;
	} else {
	  fprintf(stderr,"--use-shared-memory flag must be 0 or 1\n");
	  return 9;
	}

      } else if (!strcmp(long_name,"expand-offsets")) {
	if (!strcmp(optarg,"1")) {
	  expand_offsets_p = true;
	} else if (!strcmp(optarg,"0")) {
	  expand_offsets_p = false;
	} else {
	  fprintf(stderr,"--expand-offsets flag must be 0 or 1\n");
	  return 9;
	}

      } else if (!strcmp(long_name,"sampling")) {
	required_index1interval = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"time")) {
	timingp = true;

      } else if (!strcmp(long_name,"unload")) {
	unloadp = true;

      } else if (!strcmp(long_name,"maxsearch")) {
	maxpaths_search = atoi(optarg);

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
	  return 9;
	}

      } else if (!strcmp(long_name,"cmetdir")) {
	user_cmetdir = optarg;

      } else if (!strcmp(long_name,"atoidir")) {
	user_atoidir = optarg;

      } else if (!strcmp(long_name,"novelend-splicedist")) {
	shortsplicedist_novelend = (Chrpos_T) strtoul(optarg,NULL,10);

      } else if (!strcmp(long_name,"splicingdir")) {
	user_splicingdir = optarg;

      } else if (!strcmp(long_name,"ambig-splice-noclip")) {
	amb_clip_p = false;

      } else if (!strcmp(long_name,"find-dna-chimeras")) {
	if (!strcmp(optarg,"1")) {
	  find_dna_chimeras_p = true;
	} else if (!strcmp(optarg,"0")) {
	  find_dna_chimeras_p = false;
	} else {
	  fprintf(stderr,"--find-dna-chimeras flag must be 0 or 1\n");
	  exit(9);
	}

      } else if (!strcmp(long_name,"tallydir")) {
	user_tallydir = optarg;

      } else if (!strcmp(long_name,"use-tally")) {
	tally_root = optarg;

      } else if (!strcmp(long_name,"runlengthdir")) {
	user_runlengthdir = optarg;

      } else if (!strcmp(long_name,"use-runlength")) {
	runlength_root = optarg;

      } else if (!strcmp(long_name,"gmap-mode")) {
	gmap_mode = 0;		/* Initialize */
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
	/* No longer used by stage 2 */
	suboptimal_score_start = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"stage2-end")) {
	suboptimal_score_end = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"input-buffer-size")) {
	inbuffer_nspaces = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"output-buffer-size")) {
	output_buffer_size = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"barcode-length")) {
	barcode_length = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"fastq-id-start")) {
	acc_fieldi_start = atoi(check_valid_int(optarg)) - 1;
	if (acc_fieldi_start < 0) {
	  fprintf(stderr,"Value for fastq-id-start must be 1 or greater\n");
	  return 9;
	}

      } else if (!strcmp(long_name,"fastq-id-end")) {
	acc_fieldi_end = atoi(check_valid_int(optarg)) - 1;
	if (acc_fieldi_end < 0) {
	  fprintf(stderr,"Value for fastq-id-end must be 1 or greater\n");
	  return 9;
	}

      } else if (!strcmp(long_name,"force-single-end")) {
	force_single_end_p = true;

      } else if (!strcmp(long_name,"filter-chastity")) {
	if (!strcmp(optarg,"off")) {
	  filter_chastity_p = false;
	  filter_if_both_p = false;
	} else if (!strcmp(optarg,"either")) {
	  filter_chastity_p = true;
	  filter_if_both_p = false;
	} else if (!strcmp(optarg,"both")) {
	  filter_chastity_p = true;
	  filter_if_both_p = true;
	} else {
	  fprintf(stderr,"--filter-chastity values allowed: off, either, both\n");
	  return 9;
	}

      } else if (!strcmp(long_name,"allow-pe-name-mismatch")) {
	allow_paired_end_mismatch_p = true;

#ifdef HAVE_ZLIB
      } else if (!strcmp(long_name,"gunzip")) {
	gunzip_p = true;
#endif

#ifdef HAVE_BZLIB
      } else if (!strcmp(long_name,"bunzip2")) {
	bunzip2_p = true;
#endif

      } else if (!strcmp(long_name,"orientation")) {
	if (!strcmp(optarg,"FR")) {
	  invert_first_p = false;
	  invert_second_p = true;
	} else if (!strcmp(optarg,"RF")) {
	  invert_first_p = true;
	  invert_second_p = false;
	} else if (!strcmp(optarg,"FF")) {
	  invert_first_p = invert_second_p = false;
	} else {
	  fprintf(stderr,"Currently allowed values for orientation: FR (fwd-rev), RF (rev-fwd) or FF (fwd-fwd)\n");
	  return 9;
	}

      } else if (!strcmp(long_name,"split-output")) {
	split_output_root = optarg;

      } else if (!strcmp(long_name,"failed-input")) {
	failedinput_root = optarg;

      } else if (!strcmp(long_name,"append-output")) {
	appendp = true;

      } else if (!strcmp(long_name,"order-among-best")) {
	if (!strcmp(optarg,"genomic")) {
	  want_random_p = false;
	} else if (!strcmp(optarg,"random")) {
	  want_random_p = true;
	} else {
	  fprintf(stderr,"--order-among-best values allowed: genomic, random (default)\n");
	  return 9;
	}

      } else if (!strcmp(long_name,"pairmax-dna")) {
	pairmax_dna = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"pairmax-rna")) {
	pairmax_rna = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"pairexpect")) {
	expected_pairlength = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"pairdev")) {
	pairlength_deviation = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"min-coverage")) {
	user_mincoverage_float = atof(check_valid_float_or_int(optarg));
	if (user_mincoverage_float > 1.0 && user_mincoverage_float != rint(user_mincoverage_float)) {
	  fprintf(stderr,"Cannot specify fractional value %f for --max-mismatches except between 0.0 and 1.0\n",user_mincoverage_float);
	  return 9;
	}

      } else if (!strcmp(long_name,"indel-endlength")) {
	min_indel_end_matches = atoi(check_valid_int(optarg));
	if (min_indel_end_matches > 14) {
	  allow_end_indels_p = false;
	}

      } else if (!strcmp(long_name,"antistranded-penalty")) {
	antistranded_penalty = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"favor-multiexon")) {
	favor_multiexon_p = true;

      } else if (!strcmp(long_name,"merge-distant-samechr")) {
	merge_samechr_p = true;

      } else if (!strcmp(long_name,"query-unk-mismatch")) {
	if (!strcmp(optarg,"1")) {
	  query_unk_mismatch_p = true;
	} else if (!strcmp(optarg,"0")) {
	  query_unk_mismatch_p = false;
	} else {
	  fprintf(stderr,"--query-unk-mismatch flag must be 0 or 1\n");
	  return 9;
	}

      } else if (!strcmp(long_name,"genome-unk-mismatch")) {
	if (!strcmp(optarg,"1")) {
	  genome_unk_mismatch_p = true;
	} else if (!strcmp(optarg,"0")) {
	  genome_unk_mismatch_p = false;
	} else {
	  fprintf(stderr,"--genome-unk-mismatch flag must be 0 or 1\n");
	  return 9;
	}

      } else if (!strcmp(long_name,"trim-mismatch-score")) {
	trim_mismatch_score = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"trim-indel-score")) {
	trim_indel_score = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"distant-splice-identity")) {
	min_distantsplicing_identity = check_valid_float(optarg,long_name);

      } else if (!strcmp(long_name,"force-xs-dir")) {
	force_xs_direction_p = true;

      } else if (!strcmp(long_name,"show-refdiff")) {
	show_refdiff_p = true;

      } else if (!strcmp(long_name,"clip-overlap")) {
	clip_overlap_p = true;

      } else if (!strcmp(long_name,"merge-overlap")) {
	merge_overlap_p = true;

      } else if (!strcmp(long_name,"no-sam-headers")) {
	sam_headers_p = false;

      } else if (!strcmp(long_name,"add-paired-nomappers")) {
	add_paired_nomappers_p = true;

      } else if (!strcmp(long_name,"paired-flag-means-concordant")) {
	if (!strcmp(optarg,"1")) {
	  paired_flag_means_concordant_p = true;
	} else if (!strcmp(optarg,"0")) {
	  paired_flag_means_concordant_p = false; /* Default */
	} else {
	  fprintf(stderr,"--paired-flag-means-concordant flag must be 0 or 1\n");
	  return 9;
	}

      } else if (!strcmp(long_name,"sam-headers-batch")) {
	sam_headers_batch = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"sam-use-0M")) {
	sam_insert_0M_p = true;

      } else if (!strcmp(long_name,"sam-multiple-primaries")) {
	sam_multiple_primaries_p = true;

      } else if (!strcmp(long_name,"quality-protocol")) {
	if (user_quality_score_adj == true) {
	  fprintf(stderr,"Cannot specify both -J (--quality-zero-score) and --quality-protocol\n");
	  return 9;
	} else if (user_quality_shift == true) {
	  fprintf(stderr,"Cannot specify both -j (--quality-print-shift) and --quality-protocol\n");
	  return 9;
	} else if (!strcmp(optarg,"illumina")) {
	  MAPQ_init(/*quality_score_adj*/64);
	  Pair_init(/*quality_score_adj*/64);
	  user_quality_score_adj = true;
	  quality_shift = -31;
	  user_quality_shift = true;
	} else if (!strcmp(optarg,"sanger")) {
	  MAPQ_init(/*quality_score_adj*/33);
	  Pair_init(/*quality_score_adj*/33);
	  user_quality_score_adj = true;
	  quality_shift = 0;
	  user_quality_shift = true;
	} else {
	  fprintf(stderr,"The only values allowed for --quality-protocol are illumina or sanger\n");
	  return 9;
	}

      } else if (!strcmp(long_name,"md-lowercase-snp")) {
	md_lowercase_variant_p = true;

      } else if (!strcmp(long_name,"extend-soft-clips")) {
	hide_soft_clips_p = true;

      } else if (!strcmp(long_name,"read-group-id")) {
	sam_read_group_id = optarg;

      } else if (!strcmp(long_name,"read-group-name")) {
	sam_read_group_name = optarg;

      } else if (!strcmp(long_name,"read-group-library")) {
	sam_read_group_library = optarg;

      } else if (!strcmp(long_name,"read-group-platform")) {
	sam_read_group_platform = optarg;

#ifdef USE_MPI
      } else if (!strcmp(long_name,"master-is-worker")) {
	if (!strcmp(optarg,"1")) {
	  master_is_worker_p = true;
	} else if (!strcmp(optarg,"0")) {
	  master_is_worker_p = false; /* Default */
	} else {
	  fprintf(stderr,"--master-is-worker flag must be 0 or 1\n");
	  return 9;
	}
#endif

      } else if (!strcmp(long_name,"print-snps")) {
	print_snplabels_p = true;

      } else if (!strcmp(long_name,"failsonly")) {
	if (nofailsp == true) {
	  fprintf(stderr,"Cannot specify both --nofails and --failsonly\n");
	  return 9;
	} else {
	  failsonlyp = true;
	}
      } else if (!strcmp(long_name,"nofails")) {
	if (failsonlyp == true) {
	  fprintf(stderr,"Cannot specify both --nofails and --failsonly\n");
	  return 9;
	} else {
	  nofailsp = true;
	}
      } else {
	/* Shouldn't reach here */
	fprintf(stderr,"Don't recognize option %s.  For usage, run 'gsnap --help'",long_name);
	return 9;
      }
      break;

    case 'D': user_genomedir = optarg; break;

    case 'd': dbroot = optarg; break;

    case 'k':
      required_index1part = atoi(check_valid_int(optarg));
      if (required_index1part > 16) {
	fprintf(stderr,"The value for k-mer size must be 16 or less\n");
	return 9;
      }
      break;

    case 'G': uncompressedp = true; break;

    case 'q': parse_part(&part_modulus,&part_interval,optarg); break;
    case 'o': output_file = optarg; break;

    case 'a': 
      if (!strcmp(optarg,"paired")) {
	chop_primers_p = true;
      } else if (!strcmp(optarg,"off")) {
	chop_primers_p = false;
      } else {
	fprintf(stderr,"Currently allowed values for adapter stripping (-a): off, paired\n");
	return 9;
      }
      break;

    case 'N':
      if (!strcmp(optarg,"1")) {
	novelsplicingp = true;
      } else if (!strcmp(optarg,"0")) {
	novelsplicingp = false;
      } else {
	fprintf(stderr,"Novel splicing (-N flag) must be 0 or 1\n");
	return 9;
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
	return 9;
      }
      break;
#endif

    case 'M': subopt_levels = atoi(check_valid_int(optarg)); break;
    case 'm':
      user_maxlevel_float = atof(check_valid_float_or_int(optarg));
      if (user_maxlevel_float > 1.0 && user_maxlevel_float != rint(user_maxlevel_float)) {
	fprintf(stderr,"Cannot specify fractional value %f for --max-mismatches except between 0.0 and 1.0\n",user_maxlevel_float);
	return 9;
      } else if (user_maxlevel_float > 0.10 && user_maxlevel_float < 1.0) {
	fprintf(stderr,"Your value %f for --max-mismatches implies more than 10%% mismatches, which does not make sense\n",
		user_maxlevel_float);
	return 9;
      }
      break;

    case 'i': indel_penalty_middle = indel_penalty_end = atoi(check_valid_int(optarg)); break;

    case 'y': max_middle_insertions = atoi(check_valid_int(optarg)); break;
    case 'Y': max_end_insertions = atoi(check_valid_int(optarg)); break;
    case 'z': max_middle_deletions = atoi(check_valid_int(optarg)); break;
    case 'Z': max_end_deletions = atoi(check_valid_int(optarg)); break;

    case 'w': shortsplicedist = strtoul(optarg,NULL,10); break;

    case 'E': distantsplicing_penalty = atoi(check_valid_int(optarg)); break;
    case 'e': localsplicing_penalty = atoi(check_valid_int(optarg)); break;
    case 'K': min_distantsplicing_end_matches = atoi(check_valid_int(optarg)); break;
    case 'l': min_shortend = atoi(check_valid_int(optarg)); break;

    case 'g': genes_file = optarg; break;

    case 's':
      splicing_file = optarg;
      knownsplicingp = true;
      break;

    case 'V': user_snpsdir = optarg; break;

    case 'v': snps_root = optarg; break;

    case 'B':
      if (!strcmp(optarg,"5")) {
#if 0
	/* Not true.  -B 5 allocates suffix array and suffix aux files */
	fprintf(stderr,"Note: Batch mode 5 is now the same as batch mode 4.\n");
	fprintf(stderr,"Expansion of offsets is now controlled separately by --expand-offsets (default=0).\n");
#endif
	offsetsstrm_access = USE_ALLOCATE; /* Doesn't matter */
	positions_access = USE_ALLOCATE;
	genome_access = USE_ALLOCATE;
	sarray_access = USE_ALLOCATE;
	lcp_access = USE_ALLOCATE;
	guideexc_access = USE_ALLOCATE;
	indexij_access = USE_ALLOCATE;

#ifdef HAVE_MMAP
      } else if (!strcmp(optarg,"4")) {
	offsetsstrm_access = USE_ALLOCATE;
	positions_access = USE_ALLOCATE;
	genome_access = USE_ALLOCATE;
	sarray_access = USE_MMAP_PRELOAD;
	lcp_access = USE_MMAP_PRELOAD;
	guideexc_access = USE_ALLOCATE;
	indexij_access = USE_ALLOCATE;

      } else if (!strcmp(optarg,"3")) {
	offsetsstrm_access = USE_ALLOCATE;
	positions_access = USE_ALLOCATE;
	genome_access = USE_MMAP_PRELOAD; /* was batch_genome_p = true */
	sarray_access = USE_MMAP_ONLY;
	lcp_access = USE_MMAP_PRELOAD;
	guideexc_access = USE_MMAP_PRELOAD;
	indexij_access = USE_ALLOCATE;

      } else if (!strcmp(optarg,"2")) {
	offsetsstrm_access = USE_ALLOCATE; /* was batch_offsets_p = true */
	positions_access = USE_MMAP_PRELOAD; /* was batch_positions_p = true */
	genome_access = USE_MMAP_PRELOAD; /* was batch_genome_p = true */
	sarray_access = USE_MMAP_ONLY;
	lcp_access = USE_MMAP_ONLY;
	guideexc_access = USE_MMAP_ONLY;
	indexij_access = USE_ALLOCATE;

      } else if (!strcmp(optarg,"1")) {
	offsetsstrm_access = USE_ALLOCATE; /* was batch_offsets_p = true */
	positions_access = USE_MMAP_PRELOAD; /* was batch_positions_p = true */
	genome_access = USE_MMAP_ONLY; /* was batch_genome_p = false */
	sarray_access = USE_MMAP_ONLY;
	guideexc_access = USE_MMAP_ONLY;
	indexij_access = USE_ALLOCATE;

      } else if (!strcmp(optarg,"0")) {
	offsetsstrm_access = USE_ALLOCATE; /* was batch_offsets_p = true */
	positions_access = USE_MMAP_ONLY; /* was batch_positions_p = false */
	genome_access = USE_MMAP_ONLY; /* was batch_genome_p = false */
	sarray_access = USE_MMAP_ONLY;
	lcp_access = USE_MMAP_ONLY;
	guideexc_access = USE_MMAP_ONLY;
	indexij_access = USE_ALLOCATE;

#endif
      } else {
#ifdef HAVE_MMAP
	fprintf(stderr,"Batch mode %s not recognized.  Only allow 0-5.  Run 'gsnap --help' for more information.\n",optarg);
#else
	fprintf(stderr,"Batch mode %s not recognized.  Only allow 4-5, since mmap is disabled.  Run 'gsnap --help' for more information.\n",optarg);
#endif
	return 9;
      }
      break;

#if defined(HAVE_PTHREAD)
    case 't': nthreads = atoi(check_valid_int(optarg)); break;
#else
    case 't': fprintf(stderr,"This version of GSNAP has pthreads disabled, so ignoring the value of %s for -t\n",optarg); break;
#endif

    case 'A':
      if (!strcmp(optarg,"sam")) {
	output_sam_p = true;
      } else if (!strcmp(optarg,"m8")) {
	print_m8_p = true;
      } else {
	fprintf(stderr,"Output format %s not recognized.  Allowed values: sam, m8\n",optarg);
	return 9;
      }
      break;

    case 'j':
      if (user_quality_shift == true) {
	fprintf(stderr,"Cannot specify both -j (--quality-print-shift) and --quality-protocol\n");
	return 9;
      } else {
	quality_shift = atoi(check_valid_int(optarg));
	user_quality_shift = true;
      }
      break;

    case 'J':
      if (user_quality_score_adj == true) {
	fprintf(stderr,"Cannot specify both -J (--quality-zero-score) and --quality-protocol\n");
	return 9;
      } else {
	MAPQ_init(/*quality_score_adj*/atoi(check_valid_int(optarg)));
	Pair_init(/*quality_score_adj*/atoi(check_valid_int(optarg)));
	user_quality_score_adj = true;
      }
      break;

    case '0': exception_raise_p = false; break; /* Allows signals to pass through */
    case 'n': maxpaths_report = atoi(check_valid_int(optarg)); break;
    case 'Q': quiet_if_excessive_p = true; break;

    case 'O': orderedp = true; break;

    case '?': fprintf(stderr,"For usage, run 'gsnap --help'\n"); return 9;
    default: return 9;
    }
  }

  /* Make inferences */
  if (dbroot == NULL) {
    fprintf(stderr,"Need to specify the -d flag.  For usage, run 'gsnap --help'\n");
    /* print_program_usage(); */
    return 9;
  }

  if (acc_fieldi_end < acc_fieldi_start) {
    fprintf(stderr,"--fastq-id-end must be equal to or greater than --fastq-id-start\n");
    return 9;
  }

  if (clip_overlap_p == true && merge_overlap_p == true) {
    fprintf(stderr,"Cannot specify both --clip-overlap and --merge-overlap.  Please choose one.\n");
    return 9;
  }

  if (novelsplicingp == true && knownsplicingp == true) {
    fprintf(stderr,"Novel splicing (-N) and known splicing (-s) both turned on => assume reads are RNA-Seq\n");
    find_dna_chimeras_p = false;
    pairmax = pairmax_rna;
    shortsplicedist_known = shortsplicedist;

  } else if (knownsplicingp == true) {
    fprintf(stderr,"Known splicing (-s) turned on => assume reads are RNA-Seq\n");
    find_dna_chimeras_p = false;
    pairmax = pairmax_rna;
    shortsplicedist_known = shortsplicedist;

  } else if (novelsplicingp == true) {
    fprintf(stderr,"Novel splicing (-N) turned on => assume reads are RNA-Seq\n");
    find_dna_chimeras_p = false;
    pairmax = pairmax_rna;
    shortsplicedist_known = 0;

  } else {
    /* Appears to be DNA-Seq */
    fprintf(stderr,"Neither novel splicing (-N) nor known splicing (-s) turned on => assume reads are DNA-Seq (genomic)\n");
    pairmax = pairmax_dna;
    shortsplicedist = shortsplicedist_known = 0U;
    shortsplicedist_novelend = 0U;
  }

  if (shortsplicedist_novelend > shortsplicedist) {
    fprintf(stderr,"The novelend-splicedist %d is greater than the localsplicedist %d.  Resetting novelend-splicedist to be %d\n",
	    shortsplicedist_novelend,shortsplicedist,shortsplicedist);
    shortsplicedist_novelend = shortsplicedist;
  }

  if (distantsplicing_penalty < localsplicing_penalty) {
    fprintf(stderr,"The distant splicing penalty %d cannot be less than local splicing penalty %d\n",
	    distantsplicing_penalty,localsplicing_penalty);
    return 9;
  }

  if (sam_headers_batch >= 0) {
    if (part_modulus == sam_headers_batch) {
      sam_headers_p = true;
    } else {
      sam_headers_p = false;
    }
  }

  if (sam_read_group_id == NULL && sam_read_group_name != NULL) {
    sam_read_group_id = sam_read_group_name;
  } else if (sam_read_group_id != NULL && sam_read_group_name == NULL) {
    sam_read_group_name = sam_read_group_id;
  }

  if (chop_primers_p == true) {
    if (invert_first_p == false && invert_second_p == true) {
      /* orientation FR */
    } else {
      fprintf(stderr,"Adapter stripping not currently implemented for given orientation\n");
      return 9;
    }
  }

#ifdef USE_MPI
  /* Code does allow for MPI output to stdout, but appears not to work
     yet, and may not work if rank 0 is also a worker */
  if (split_output_root == NULL && output_file == NULL) {
    fprintf(stderr,"For MPI version, need to specify either --split-output or --output-file\n");
    return 9;
  }
#endif

  return 0;
}


static bool
open_input_streams_parser (int *nextchar, int *nchars1, int *nchars2, char ***files, int *nfiles,
			   FILE **input, FILE **input2,
#ifdef HAVE_ZLIB
			   gzFile *gzipped, gzFile *gzipped2,
#endif
#ifdef HAVE_BZLIB
			   Bzip2_T *bzipped, Bzip2_T *bzipped2,
#endif
			   bool gunzip_p, bool bunzip2_p,
			   int argc, char **argv) {
  bool fastq_format_p = false;
  
  *input = *input2 = NULL;
#ifdef HAVE_ZLIB
  *gzipped = *gzipped2 = NULL;
#endif
#ifdef HAVE_BZLIB
  *bzipped = *bzipped2 = NULL;
#endif

  /* Open input stream and peek at first char */
  if (argc == 0) {
#ifdef USE_MPI
    fprintf(stderr,"For mpi_gsnap, cannot read from stdin\n");
    exit(9);
#else
    fprintf(stderr,"Reading from stdin\n");
    *input = stdin;
    *files = (char **) NULL;
    *nfiles = 0;
    *nextchar = Shortread_input_init(&(*nchars1),*input);
#endif
  } else {
    *files = argv;
    *nfiles = argc;

    if (gunzip_p == true) {
#ifdef HAVE_ZLIB
      if ((*gzipped = gzopen((*files)[0],"rb")) == NULL) {
	fprintf(stderr,"Cannot open gzipped file %s\n",(*files)[0]);
	exit(9);
      } else {
#ifdef HAVE_ZLIB_GZBUFFER
	gzbuffer(*gzipped,GZBUFFER_SIZE);
#endif
	*nextchar = Shortread_input_init_gzip(*gzipped);
      }
#endif

    } else if (bunzip2_p == true) {
#ifdef HAVE_BZLIB
      if ((*bzipped = Bzip2_new((*files)[0])) == NULL) {
	fprintf(stderr,"Cannot open bzipped file %s\n",(*files)[0]);
	exit(9);
      } else {
	*nextchar = Shortread_input_init_bzip2(*bzipped);
      }
#endif

    } else {
      if ((*input = FOPEN_READ_TEXT((*files)[0])) == NULL) {
	fprintf(stderr,"Cannot open file %s\n",(*files)[0]);
	exit(9);
      } else {
	debugf(fprintf(stderr,"Master opening file %s using fopen\n",(*files)[0]));
	*nextchar = Shortread_input_init(&(*nchars1),*input);
      }
    }

    (*files)++;
    (*nfiles)--;
  }

  /* Interpret first char to determine input type */
  if (*nextchar == EOF) {
    fprintf(stderr,"Input is empty\n");
    exit(9);

  } else if (*nextchar == '@') {
    /* Looks like a FASTQ file */
    if (*nfiles == 0 || force_single_end_p == true) {
#ifdef HAVE_ZLIB
      *gzipped2 = (gzFile) NULL;
#endif
#ifdef HAVE_BZLIB
      *bzipped2 = (Bzip2_T) NULL;
#endif
      *input2 = (FILE *) NULL;
    } else {
      if (gunzip_p == true) {
#ifdef HAVE_ZLIB
	if ((*gzipped2 = gzopen((*files)[0],"rb")) == NULL) {
	  fprintf(stderr,"Cannot open gzipped file %s\n",(*files)[0]);
	  exit(9);
	} else {
#ifdef HAVE_ZLIB_GZBUFFER
	  gzbuffer(*gzipped2,GZBUFFER_SIZE);
#endif
	  /* nextchar2 = */ Shortread_input_init_gzip(*gzipped2);
	}
#endif

      } else if (bunzip2_p == true) {
#ifdef HAVE_BZLIB
	if ((*bzipped2 = Bzip2_new((*files)[0])) == NULL) {
	  fprintf(stderr,"Cannot open bzip2 file %s\n",(*files)[0]);
	  exit(9);
	} else {
	  /* nextchar2 = */ Shortread_input_init_bzip2(*bzipped2);
	}
#endif

      } else {
	if ((*input2 = FOPEN_READ_TEXT((*files)[0])) == NULL) {
	  fprintf(stderr,"Cannot open file %s\n",(*files)[0]);
	  exit(9);
	} else {
	  debugf(fprintf(stderr,"Master opening file %s using fopen\n",(*files)[0]));
	  /* nextchar2 = */ Shortread_input_init(&(*nchars2),*input2);
	}
      }
      (*files)++;
      (*nfiles)--;
    }
    fastq_format_p = true;

  } else if (*nextchar == '>') {
    /* Looks like a FASTA file */

  } else {
    fprintf(stderr,"First char is %c.  Expecting either '>' for FASTA or '@' for FASTQ format.\n",*nextchar);
    exit(9);
  }

  return fastq_format_p;
}


#ifdef USE_MPI
static void
open_input_streams_worker (char ***files, int *nfiles,
#if defined(USE_MPI_FILE_INPUT)
			  MPI_File *input, MPI_File *input2, MPI_Comm workers_comm,
#else
			  FILE **input, FILE **input2,
#endif
#ifdef HAVE_ZLIB
			  gzFile *gzipped, gzFile *gzipped2,
#endif
#ifdef HAVE_BZLIB
			  Bzip2_T *bzipped, Bzip2_T *bzipped2,
#endif
			  bool gunzip_p, bool bunzip2_p, bool fastq_format_p,
			  int argc, char **argv) {

  *input = *input2 = NULL;
#ifdef HAVE_ZLIB
  *gzipped = *gzipped2 = NULL;
#endif
#ifdef HAVE_BZLIB
  *bzipped = *bzipped2 = NULL;
#endif

  /* Open input stream and peek at first char */
  if (argc == 0) {
    fprintf(stderr,"For mpi_gsnap, cannot read from stdin\n");
    exit(9);

  } else {
    *files = argv;
    *nfiles = argc;

    if (gunzip_p == true) {
#ifdef HAVE_ZLIB
      if ((*gzipped = gzopen((*files)[0],"rb")) == NULL) {
	fprintf(stderr,"Cannot open gzipped file %s\n",(*files)[0]);
	exit(9);
      } else {
#ifdef HAVE_ZLIB_GZBUFFER
	gzbuffer(*gzipped,GZBUFFER_SIZE);
#endif
      }
#endif

    } else if (bunzip2_p == true) {
#ifdef HAVE_BZLIB
      if ((*bzipped = Bzip2_new((*files)[0])) == NULL) {
	fprintf(stderr,"Cannot open bzipped file %s\n",(*files)[0]);
	exit(9);
      }
#endif

    } else {
#if defined(USE_MPI_FILE_INPUT)
      if ((*input = MPI_fopen((*files)[0],workers_comm)) == NULL) {
	fprintf(stderr,"Cannot open file %s\n",(*files)[0]);
	exit(9);
      }
      debugf(fprintf(stderr,"Slave opening file %s using MPI_File_open\n",(*files)[0]));
#else
      if ((*input = FOPEN_READ_TEXT((*files)[0])) == NULL) {
	fprintf(stderr,"Cannot open file %s\n",(*files)[0]);
	exit(9);
      }
      debugf(fprintf(stderr,"Slave opening file %s using fopen\n",(*files)[0]));
#endif
    }

    (*files)++;
    (*nfiles)--;
  }

  if (fastq_format_p == true) {
    /* Looks like a FASTQ file */
    if (*nfiles == 0 || force_single_end_p == true) {
#ifdef HAVE_ZLIB
      *gzipped2 = (gzFile) NULL;
#endif
#ifdef HAVE_BZLIB
      *bzipped2 = (Bzip2_T) NULL;
#endif
#if defined(USE_MPI_FILE_INPUT)
      *input2 = (MPI_File) NULL;
#else
      *input2 = (FILE *) NULL;
#endif
    } else {
      if (gunzip_p == true) {
#ifdef HAVE_ZLIB
	if ((*gzipped2 = gzopen((*files)[0],"rb")) == NULL) {
	  fprintf(stderr,"Cannot open gzipped file %s\n",(*files)[0]);
	  exit(9);
	} else {
#ifdef HAVE_ZLIB_GZBUFFER
	  gzbuffer(*gzipped2,GZBUFFER_SIZE);
#endif
	}
#endif

      } else if (bunzip2_p == true) {
#ifdef HAVE_BZLIB
	if ((*bzipped2 = Bzip2_new((*files)[0])) == NULL) {
	  fprintf(stderr,"Cannot open bzip2 file %s\n",(*files)[0]);
	  exit(9);
	}
#endif

      } else {
#if defined(USE_MPI_FILE_INPUT)
	if ((*input2 = MPI_fopen((*files)[0],workers_comm)) == NULL) {
	  fprintf(stderr,"Cannot open file %s\n",(*files)[0]);
	  exit(9);
	}
	debugf(fprintf(stderr,"Slave opening file %s using MPI_File_open\n",(*files)[0]));
#else
	if ((*input2 = FOPEN_READ_TEXT((*files)[0])) == NULL) {
	  fprintf(stderr,"Cannot open file %s\n",(*files)[0]);
	  exit(9);
	}
	debugf(fprintf(stderr,"Slave opening file %s using fopen\n",(*files)[0]));
#endif
      }
      (*files)++;
      (*nfiles)--;
    }
  }

  return;
}
#endif


static Univ_IIT_T
chromosome_iit_setup (int *nchromosomes, int *circular_typeint, bool *any_circular_p, bool **circularp,
		      char *genomesubdir, char *fileroot) {
  Univ_IIT_T chromosome_iit = NULL;
  char *iitfile = NULL;

  /* Prepare genomic data */

  iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			    strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
  sprintf(iitfile,"%s/%s.chromosome.iit",genomesubdir,fileroot);
  if ((chromosome_iit = Univ_IIT_read(iitfile,/*readonlyp*/true,/*add_iit_p*/false)) == NULL) {
    fprintf(stderr,"IIT file %s is not valid\n",iitfile);
    exit(9);
#ifdef LARGE_GENOMES
  } else if (Univ_IIT_coord_values_8p(chromosome_iit) == false) {
    fprintf(stderr,"This program gsnapl is designed for large genomes.\n");
    fprintf(stderr,"For small genomes of less than 2^32 (4 billion) bp, please run gsnap instead.\n");
    exit(9);
#endif
  } else {
    *nchromosomes = Univ_IIT_total_nintervals(chromosome_iit);
    *circular_typeint = Univ_IIT_typeint(chromosome_iit,"circular");
    *circularp = Univ_IIT_circularp(&(*any_circular_p),chromosome_iit);
  }
  FREE(iitfile);

  return chromosome_iit;
}


static void
worker_setup (char *genomesubdir, char *fileroot) {
  char *snpsdir = NULL, *modedir = NULL, *mapdir = NULL, *iitfile = NULL;
  Splicestringpool_T splicestringpool;

  if (snps_root == NULL) {
    genomecomp = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,/*genometype*/GENOME_OLIGOS,
			    uncompressedp,genome_access,sharedp);
    genomebits = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,/*genometype*/GENOME_BITS,
			    uncompressedp,genome_access,sharedp);

#ifndef LARGE_GENOMES
    if (use_sarray_p == true) {
      if (mode == STANDARD) {
	if ((sarray_fwd = Sarray_new(genomesubdir,fileroot,/*snps_root*/NULL,sarray_access,lcp_access,
				     guideexc_access,indexij_access,sharedp,mode,/*fwdp*/true)) == NULL) {
	  use_sarray_p = false;
	} else {
	  sarray_rev = sarray_fwd;
	}
      } else {
	if ((sarray_fwd = Sarray_new(genomesubdir,fileroot,/*snps_root*/NULL,sarray_access,lcp_access,
				     guideexc_access,indexij_access,sharedp,mode,/*fwdp*/true)) == NULL ||
	    (sarray_rev = Sarray_new(genomesubdir,fileroot,/*snps_root*/NULL,sarray_access,lcp_access,
				     guideexc_access,indexij_access,sharedp,mode,/*fwdp*/false)) == NULL) {
	  use_sarray_p = false;
	}
      }
    }
#endif

    if (use_only_sarray_p == true) {
      indexdb = indexdb2 = NULL;
      
    } else if (dibasep == true) {
      fprintf(stderr,"No longer supporting 2-base encoding\n");
      exit(9);
      if ((indexdb = Indexdb_new_genome(&index1part,&index1interval,
					genomesubdir,fileroot,/*idx_filesuffix*/"dibase",/*snps_root*/NULL,
					required_index1part,required_index1interval,
					expand_offsets_p,offsetsstrm_access,positions_access,sharedp)) == NULL) {
	fprintf(stderr,"Cannot find offsets file %s.%s*offsets, needed for GSNAP color mode\n",fileroot,"dibase");
	exit(9);
      }
      indexdb2 = indexdb;
      print_ncolordiffs_p = true;

    } else if (mode == CMET_STRANDED || mode == CMET_NONSTRANDED) {
      if (user_cmetdir == NULL) {
	modedir = genomesubdir;
      } else {
	modedir = user_cmetdir;
      }

      if ((indexdb = Indexdb_new_genome(&index1part,&index1interval,
					modedir,fileroot,/*idx_filesuffix*/"metct",/*snps_root*/NULL,
					required_index1part,required_index1interval,
					expand_offsets_p,offsetsstrm_access,positions_access,sharedp)) == NULL) {
	fprintf(stderr,"Cannot find metct index file.  Need to run cmetindex first\n");
	exit(9);
      }

      if ((indexdb2 = Indexdb_new_genome(&index1part,&index1interval,
					 modedir,fileroot,/*idx_filesuffix*/"metga",/*snps_root*/NULL,
					 required_index1part,required_index1interval,
					 expand_offsets_p,offsetsstrm_access,positions_access,sharedp)) == NULL) {
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
					required_index1part,required_index1interval,
					expand_offsets_p,offsetsstrm_access,positions_access,sharedp)) == NULL) {
	fprintf(stderr,"Cannot find a2iag index file.  Need to run atoiindex first\n");
	exit(9);
      }

      if ((indexdb2 = Indexdb_new_genome(&index1part,&index1interval,
					 modedir,fileroot,/*idx_filesuffix*/"a2itc",/*snps_root*/NULL,
					 required_index1part,required_index1interval,
					 expand_offsets_p,offsetsstrm_access,positions_access,sharedp)) == NULL) {
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
					required_index1part,required_index1interval,
					expand_offsets_p,offsetsstrm_access,positions_access,sharedp)) == NULL) {
	fprintf(stderr,"Cannot find a2itc index file.  Need to run atoiindex first\n");
	exit(9);
      }

      if ((indexdb2 = Indexdb_new_genome(&index1part,&index1interval,
					 modedir,fileroot,/*idx_filesuffix*/"a2iag",/*snps_root*/NULL,
					 required_index1part,required_index1interval,
					 expand_offsets_p,offsetsstrm_access,positions_access,sharedp)) == NULL) {
	fprintf(stderr,"Cannot find a2iag index file.  Need to run atoiindex first\n");
	exit(9);
      }


    } else {
      /* Standard behavior */
      if ((indexdb = Indexdb_new_genome(&index1part,&index1interval,
					genomesubdir,fileroot,IDX_FILESUFFIX,/*snps_root*/NULL,
					required_index1part,required_index1interval,
					expand_offsets_p,offsetsstrm_access,positions_access,sharedp)) == NULL) {
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
    genomecomp = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,/*genometype*/GENOME_OLIGOS,
			    uncompressedp,genome_access,sharedp);
    genomecomp_alt = Genome_new(snpsdir,fileroot,snps_root,/*genometype*/GENOME_OLIGOS,
				uncompressedp,genome_access,sharedp);
    genomebits = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,/*genometype*/GENOME_BITS,
			    uncompressedp,genome_access,sharedp);
    genomebits_alt = Genome_new(snpsdir,fileroot,snps_root,/*genometype*/GENOME_BITS,
				uncompressedp,genome_access,sharedp);

#ifndef LARGE_GENOMES
    if (use_sarray_p == true) {
      fprintf(stderr,"Note: Suffix arrays will bias against SNP-tolerant alignment.  For bias-free alignment, set --use-sarray=0\n");
      if (mode == STANDARD) {
	if ((sarray_fwd = Sarray_new(genomesubdir,fileroot,/*snps_root*/NULL,sarray_access,lcp_access,
				     guideexc_access,indexij_access,sharedp,mode,/*fwdp*/true)) == NULL) {
	  use_sarray_p = false;
	} else {
	  sarray_rev = sarray_fwd;
	}
      } else {
	if ((sarray_fwd = Sarray_new(genomesubdir,fileroot,/*snps_root*/NULL,sarray_access,lcp_access,
				     guideexc_access,indexij_access,sharedp,mode,/*fwdp*/true)) == NULL ||
	    (sarray_rev = Sarray_new(genomesubdir,fileroot,/*snps_root*/NULL,sarray_access,lcp_access,
				     guideexc_access,indexij_access,sharedp,mode,/*fwdp*/false)) == NULL) {
	  use_sarray_p = false;
	}
      }
    }
#endif

    if (dibasep == true) {
      fprintf(stderr,"Currently cannot combine SNPs with 2-base encoding\n");
      exit(9);
      print_ncolordiffs_p = true;

    } else if (mode == CMET_STRANDED || mode == CMET_NONSTRANDED) {
      if (user_cmetdir == NULL) {
	modedir = snpsdir;
      } else {
	modedir = user_cmetdir;
      }

      if ((indexdb = Indexdb_new_genome(&index1part,&index1interval,
					modedir,fileroot,/*idx_filesuffix*/"metct",snps_root,
					required_index1part,required_index1interval,
					expand_offsets_p,offsetsstrm_access,positions_access,sharedp)) == NULL) {
	fprintf(stderr,"Cannot find metct index file.  Need to run cmetindex first\n");
	exit(9);
      }
      if ((indexdb2 = Indexdb_new_genome(&index1part,&index1interval,
					 modedir,fileroot,/*idx_filesuffix*/"metga",snps_root,
					 required_index1part,required_index1interval,
					 expand_offsets_p,offsetsstrm_access,positions_access,sharedp)) == NULL) {
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
					required_index1part,required_index1interval,
					expand_offsets_p,offsetsstrm_access,positions_access,sharedp)) == NULL) {
	fprintf(stderr,"Cannot find a2iag index file.  Need to run atoiindex first\n");
	exit(9);
      }
      if ((indexdb2 = Indexdb_new_genome(&index1part,&index1interval,
					 modedir,fileroot,/*idx_filesuffix*/"a2itc",snps_root,
					 required_index1part,required_index1interval,
					 expand_offsets_p,offsetsstrm_access,positions_access,sharedp)) == NULL) {
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
					required_index1part,required_index1interval,
					expand_offsets_p,offsetsstrm_access,positions_access,sharedp)) == NULL) {
	fprintf(stderr,"Cannot find a2itc index file.  Need to run atoiindex first\n");
	exit(9);
      }
      if ((indexdb2 = Indexdb_new_genome(&index1part,&index1interval,
					 modedir,fileroot,/*idx_filesuffix*/"a2iag",snps_root,
					 required_index1part,required_index1interval,
					 expand_offsets_p,offsetsstrm_access,positions_access,sharedp)) == NULL) {
	fprintf(stderr,"Cannot find a2iag index file.  Need to run atoiindex first\n");
	exit(9);
      }

    } else {
      indexdb = Indexdb_new_genome(&index1part,&index1interval,
				   snpsdir,fileroot,/*idx_filesuffix*/"ref",snps_root,
				   required_index1part,required_index1interval,
				   expand_offsets_p,offsetsstrm_access,positions_access,sharedp);
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
	fprintf(stderr,"Either install file %s.iit or specify a directory for the IIT file\n",snps_root);
	fprintf(stderr,"using the -M flag.\n");
	exit(9);
      }
    }

    print_nsnpdiffs_p = true;
    snps_divint_crosstable = Univ_IIT_divint_crosstable(chromosome_iit,snps_iit);

    fprintf(stderr,"done\n");
    FREE(iitfile);
    if (user_snpsdir == NULL) {
      FREE(mapdir);
    }
  }

  if (min_distantsplicing_end_matches < index1part) {
    fprintf(stderr,"Minimum value for distant-splice-endlength is the value for -k (kmer size) %d\n",index1part);
    exit(9);
  }

  if (use_only_sarray_p == false) {
    Compoundpos_init_positions_free(Indexdb_positions_fileio_p(indexdb));
    Spanningelt_init_positions_free(Indexdb_positions_fileio_p(indexdb));
    Stage1_init_positions_free(Indexdb_positions_fileio_p(indexdb));

    indexdb_size_threshold = (int) (10*Indexdb_mean_size(indexdb,mode,index1part));
    debug(printf("Size threshold is %d\n",indexdb_size_threshold));
    if (indexdb_size_threshold < MIN_INDEXDB_SIZE_THRESHOLD) {
      indexdb_size_threshold = MIN_INDEXDB_SIZE_THRESHOLD;
    }
  }

  if (genes_file != NULL) {
    if ((genes_iit = IIT_read(genes_file,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
			      /*divstring*/NULL,/*add_iit_p*/true,/*labels_read_p*/true)) != NULL) {
      fprintf(stderr,"Reading genes file %s locally...",genes_file);
    } else {
      mapdir = Datadir_find_mapdir(/*user_mapdir*/NULL,genomesubdir,fileroot);
      iitfile = (char *) CALLOC(strlen(mapdir)+strlen("/")+
				strlen(genes_file)+1,sizeof(char));
      sprintf(iitfile,"%s/%s",mapdir,genes_file);
      if ((genes_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				/*divstring*/NULL,/*add_iit_p*/true,/*labels_read_p*/true)) != NULL) {
	fprintf(stderr,"Reading genes file %s...",iitfile);
	FREE(iitfile);
	FREE(mapdir);
      } else {
	fprintf(stderr,"Genes file %s.iit not found locally or in %s.  Available files:\n",genes_file,mapdir);
	Datadir_list_directory(stderr,mapdir);
	fprintf(stderr,"Either install file %s or specify a full directory path\n",genes_file);
	exit(9);
      }
    }
    genes_divint_crosstable = Univ_IIT_divint_crosstable(chromosome_iit,genes_iit);
  }


  if (splicing_file != NULL) {
    if (user_splicingdir == NULL) {
      if ((splicing_iit = IIT_read(splicing_file,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				   /*divstring*/NULL,/*add_iit_p*/true,/*labels_read_p*/true)) != NULL) {
	fprintf(stderr,"Reading splicing file %s locally...",splicing_file);
      }
    } else {
      iitfile = (char *) CALLOC(strlen(user_splicingdir)+strlen("/")+strlen(splicing_file)+1,sizeof(char));
      sprintf(iitfile,"%s/%s",user_splicingdir,splicing_file);
      if ((splicing_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				   /*divstring*/NULL,/*add_iit_p*/true,/*labels_read_p*/true)) != NULL) {
	fprintf(stderr,"Reading splicing file %s...",iitfile);
	FREE(iitfile);
      }
    }

    if (splicing_iit == NULL) {
      mapdir = Datadir_find_mapdir(/*user_mapdir*/NULL,genomesubdir,fileroot);
      iitfile = (char *) CALLOC(strlen(mapdir)+strlen("/")+
				strlen(splicing_file)+1,sizeof(char));
      sprintf(iitfile,"%s/%s",mapdir,splicing_file);
      if ((splicing_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				      /*divstring*/NULL,/*add_iit_p*/true,/*labels_read_p*/true)) != NULL) {
	fprintf(stderr,"Reading splicing file %s...",iitfile);
	FREE(iitfile);
	FREE(mapdir);
      } else {
	fprintf(stderr,"Splicing file %s.iit not found locally or in %s.  Available files:\n",splicing_file,mapdir);
	Datadir_list_directory(stderr,mapdir);
	fprintf(stderr,"Either install file %s or specify a full directory path\n",splicing_file);
	exit(9);
      }
    }

    splicing_divint_crosstable = Univ_IIT_divint_crosstable(chromosome_iit,splicing_iit);
    if ((donor_typeint = IIT_typeint(splicing_iit,"donor")) >= 0 && 
	(acceptor_typeint = IIT_typeint(splicing_iit,"acceptor")) >= 0) {
      fprintf(stderr,"found donor and acceptor tags, so treating as splicesites file\n");
      splicestringpool = Splicestringpool_new();
      splicesites = Splicetrie_retrieve_via_splicesites(&distances_observed_p,&splicecomp,&splicetypes,&splicedists,
							&splicestrings,&splicefrags_ref,&splicefrags_alt,
							&nsplicesites,splicing_iit,splicing_divint_crosstable,
							donor_typeint,acceptor_typeint,chromosome_iit,
							genomecomp,genomecomp_alt,shortsplicedist,splicestringpool);
      if (nsplicesites == 0) {
	fprintf(stderr,"\nWarning: No splicesites observed for genome %s.  Are you sure this splicesite file was built for this genome?  Please compare chromosomes below:\n",
		dbroot);
	fprintf(stderr,"Chromosomes in the genome: ");
	Univ_IIT_dump_labels(stderr,chromosome_iit);
	fprintf(stderr,"Chromosomes in the splicesites IIT file: ");
	IIT_dump_divstrings(stderr,splicing_iit);
	exit(9);

      } else {
	Splicetrie_npartners(&nsplicepartners_skip,&nsplicepartners_obs,&nsplicepartners_max,splicesites,splicetypes,splicedists,
			     splicestrings,nsplicesites,chromosome_iit,shortsplicedist,distances_observed_p);
	Splicetrie_build_via_splicesites(&triecontents_obs,&trieoffsets_obs,&triecontents_max,&trieoffsets_max,
					 nsplicepartners_skip,nsplicepartners_obs,nsplicepartners_max,splicetypes,
					 splicestrings,nsplicesites);
	FREE(nsplicepartners_max);
	FREE(nsplicepartners_obs);
	FREE(nsplicepartners_skip);
	/* Splicestring_gc(splicestrings,nsplicesites); */
	FREE(splicestrings);
      }
      Splicestringpool_free(&splicestringpool);

    } else {
      fprintf(stderr,"no donor or acceptor tags found, so treating as introns file\n");
      splicestringpool = Splicestringpool_new();
      splicesites = Splicetrie_retrieve_via_introns(&splicecomp,&splicetypes,&splicedists,
						    &splicestrings,&splicefrags_ref,&splicefrags_alt,
						    &nsplicesites,splicing_iit,splicing_divint_crosstable,
						    chromosome_iit,genomecomp,genomecomp_alt,splicestringpool);
      if (nsplicesites == 0) {
	fprintf(stderr,"\nWarning: No splicesites observed for genome %s.  Are you sure this splicesite file was built for this genome?  Please compare chromosomes below:\n",
		dbroot);
	fprintf(stderr,"Chromosomes in the genome: ");
	Univ_IIT_dump_labels(stderr,chromosome_iit);
	fprintf(stderr,"Chromosomes in the splicesites IIT file: ");
	IIT_dump_divstrings(stderr,splicing_iit);
	exit(9);
      } else {
	Splicetrie_build_via_introns(&triecontents_obs,&trieoffsets_obs,splicesites,splicetypes,
				     splicestrings,nsplicesites,chromosome_iit,splicing_iit,splicing_divint_crosstable);
	triecontents_max = (Triecontent_T *) NULL;
	trieoffsets_max =  (Trieoffset_T *) NULL;
	/* Splicestring_gc(splicestrings,nsplicesites); */
	FREE(splicestrings);
      }
      Splicestringpool_free(&splicestringpool);
      
    }

    /* For benchmarking purposes.  Can spend time/memory to load
       splicesites, but then not use them. */
    if (unloadp == true) {
      fprintf(stderr,"unloading...");

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
	FREE(splicesites);
	FREE(splicecomp);
	nsplicesites = 0;
      }

      FREE(splicing_divint_crosstable);
      IIT_free(&splicing_iit);
      splicing_iit = NULL;
      knownsplicingp = false;
      splicing_file = (char *) NULL;
    }

    fprintf(stderr,"done\n");
  }


  if (tally_root != NULL) {
    if (user_tallydir == NULL) {
      if ((tally_iit = IIT_read(tally_root,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				/*divstring*/NULL,/*add_iit_p*/true,/*labels_read_p*/true)) != NULL) {
	fprintf(stderr,"Reading tally file %s.iit locally...",tally_root);
      }
    } else {
      iitfile = (char *) CALLOC(strlen(user_tallydir)+strlen("/")+strlen(tally_root)+1,sizeof(char));
      sprintf(iitfile,"%s/%s",user_tallydir,tally_root);
      if ((tally_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				/*divstring*/NULL,/*add_iit_p*/true,/*labels_read_p*/true)) != NULL) {
	fprintf(stderr,"Reading tally file %s...",iitfile);
	FREE(iitfile);
      }
    }

    if (tally_iit == NULL) {
      iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(tally_root)+1,sizeof(char));
      sprintf(iitfile,"%s/%s",genomesubdir,tally_root);
      if ((tally_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				/*divstring*/NULL,/*add_iit_p*/true,/*labels_read_p*/true)) != NULL) {
	fprintf(stderr,"Reading tally file %s...",iitfile);
	FREE(iitfile);
      } else {
	fprintf(stderr,"Tally file %s.iit not found locally",tally_root);
	if (user_tallydir != NULL) {
	  fprintf(stderr," or in %s",user_tallydir);
	}
	fprintf(stderr," or in %s\n",genomesubdir);
	exit(9);
      }
    }

    tally_divint_crosstable = Univ_IIT_divint_crosstable(chromosome_iit,tally_iit);
    fprintf(stderr,"done\n");
  }


  if (runlength_root != NULL) {
    if (user_runlengthdir == NULL) {
      if ((runlength_iit = IIT_read(runlength_root,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				/*divstring*/NULL,/*add_iit_p*/true,/*labels_read_p*/true)) != NULL) {
	fprintf(stderr,"Reading runlength file %s.iit locally...",runlength_root);
      }
    } else {
      iitfile = (char *) CALLOC(strlen(user_runlengthdir)+strlen("/")+strlen(runlength_root)+1,sizeof(char));
      sprintf(iitfile,"%s/%s",user_runlengthdir,runlength_root);
      if ((runlength_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				/*divstring*/NULL,/*add_iit_p*/true,/*labels_read_p*/true)) != NULL) {
	fprintf(stderr,"Reading runlength file %s...",iitfile);
	FREE(iitfile);
      }
    }

    if (runlength_iit == NULL) {
      iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(runlength_root)+1,sizeof(char));
      sprintf(iitfile,"%s/%s",genomesubdir,runlength_root);
      if ((runlength_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				/*divstring*/NULL,/*add_iit_p*/true,/*labels_read_p*/true)) != NULL) {
	fprintf(stderr,"Reading runlength file %s...",iitfile);
	FREE(iitfile);
      } else {
	fprintf(stderr,"Runlength file %s.iit not found locally",runlength_root);
	if (user_runlengthdir != NULL) {
	  fprintf(stderr," or in %s",user_runlengthdir);
	}
	fprintf(stderr," or in %s\n",genomesubdir);
	exit(9);
      }
    }

    runlength_divint_crosstable = Univ_IIT_divint_crosstable(chromosome_iit,runlength_iit);
    fprintf(stderr,"done\n");
  }


  Genome_setup(genomecomp,genomecomp_alt,mode,circular_typeint);
#ifndef LARGE_GENOMES
  if (sarray_fwd != NULL && sarray_rev != NULL) {
    Sarray_setup(sarray_fwd,sarray_rev,genomecomp,mode,chromosome_iit,circular_typeint,circularp,
		 shortsplicedist,localsplicing_penalty,
		 max_deletionlength,max_end_deletions,max_middle_insertions,max_end_insertions,
		 splicesites,splicetypes,splicedists,nsplicesites);
  }
#endif

  if (genomebits == NULL) {
    fprintf(stderr,"This version of GSNAP requires the genomebits128 file\n");
    exit(9);
  } else {
    Genome_hr_setup(Genome_blocks(genomebits),/*snp_blocks*/genomebits_alt ? Genome_blocks(genomebits_alt) : NULL,
		    query_unk_mismatch_p,genome_unk_mismatch_p,mode);
  }
  Genome_sites_setup(Genome_blocks(genomecomp),/*snp_blocks*/genomecomp_alt ? Genome_blocks(genomecomp_alt) : NULL);
  Maxent_hr_setup(Genome_blocks(genomecomp),/*snp_blocks*/genomecomp_alt ? Genome_blocks(genomecomp_alt) : NULL);

  Dynprog_init(mode);			    /* Need Dynprog for sarray_gmap */
  Dynprog_single_setup(/*homopolymerp*/false);
  Dynprog_genome_setup(novelsplicingp,splicing_iit,splicing_divint_crosstable,
		       donor_typeint,acceptor_typeint);
  Dynprog_end_setup(splicesites,splicetypes,splicedists,nsplicesites,
		    trieoffsets_obs,triecontents_obs,trieoffsets_max,triecontents_max);
  Pair_setup(trim_mismatch_score,trim_indel_score,/*gff3_separators_p*/false,sam_insert_0M_p,
	     force_xs_direction_p,md_lowercase_variant_p,
	     /*snps_p*/snps_iit ? true : false,print_nsnpdiffs_p,
	     Univ_IIT_genomelength(chromosome_iit,/*with_circular_alias*/false));
  Stage3_setup(/*splicingp*/novelsplicingp == true || knownsplicingp == true,novelsplicingp,
	       /*require_splicedir_p*/true,splicing_iit,splicing_divint_crosstable,
	       donor_typeint,acceptor_typeint,
	       splicesites,min_intronlength,max_deletionlength,min_indel_end_matches,
	       maxpeelback_distalmedial,nullgap,extramaterial_end,extramaterial_paired,
	       extraband_single,extraband_end,extraband_paired,
	       ngap,/*maxintronlen*/shortsplicedist,
	       output_sam_p,/*homopolymerp*/false,/*stage3debug*/NO_STAGE3DEBUG);
  Oligoindex_hr_setup(Genome_blocks(genomecomp),mode);
  Stage2_setup(/*splicingp*/novelsplicingp == true || knownsplicingp == true,/*cross_species_p*/false,
	       suboptimal_score_start,suboptimal_score_end,sufflookback,nsufflookback,
	       /*maxintronlen*/shortsplicedist,mode,/*snps_p*/snps_iit ? true : false);

  if (use_only_sarray_p == true) {
    spansize = 1;
  } else {
    Indexdb_setup(index1part);
    Indexdb_hr_setup(index1part);
    Oligo_setup(index1part);
    spansize = Spanningelt_setup(index1part,index1interval);
  }

  Splicetrie_setup(splicecomp,splicesites,splicefrags_ref,splicefrags_alt,
		   trieoffsets_obs,triecontents_obs,trieoffsets_max,triecontents_max,
		   /*snpp*/snps_iit ? true : false,amb_closest_p,amb_clip_p,min_shortend);
  Splice_setup(min_shortend);
  Indel_setup(min_indel_end_matches,indel_penalty_middle);
  Stage1hr_setup(use_sarray_p,use_only_sarray_p,index1part,index1interval,spansize,chromosome_iit,nchromosomes,
		 genomecomp,genomecomp_alt,mode,maxpaths_search,
		 splicesites,splicetypes,splicedists,nsplicesites,
		 novelsplicingp,knownsplicingp,find_dna_chimeras_p,distances_observed_p,
		 subopt_levels,max_middle_insertions,max_middle_deletions,
		 shortsplicedist,shortsplicedist_known,shortsplicedist_novelend,min_intronlength,
		 min_distantsplicing_end_matches,min_distantsplicing_identity,
		 nullgap,maxpeelback,maxpeelback_distalmedial,
		 extramaterial_end,extramaterial_paired,gmap_mode,
		 trigger_score_for_gmap,gmap_allowance,max_gmap_pairsearch,
		 max_gmap_terminal,max_gmap_improvement,antistranded_penalty);
  Substring_setup(print_nsnpdiffs_p,print_snplabels_p,
		  show_refdiff_p,snps_iit,snps_divint_crosstable,
		  genes_iit,genes_divint_crosstable,
		  splicing_iit,splicing_divint_crosstable,
		  donor_typeint,acceptor_typeint,trim_mismatch_score,
		  novelsplicingp,knownsplicingp,output_sam_p,mode,
		  Univ_IIT_genomelength(chromosome_iit,/*with_circular_alias*/false));
  Stage3hr_setup(invert_first_p,invert_second_p,genomecomp,chromosome_iit,nchromosomes,circular_typeint,
		 genes_iit,genes_divint_crosstable,
		 tally_iit,tally_divint_crosstable,runlength_iit,runlength_divint_crosstable,
		 distances_observed_p,pairmax,
		 expected_pairlength,pairlength_deviation,
		 localsplicing_penalty,indel_penalty_middle,antistranded_penalty,
		 favor_multiexon_p,gmap_min_nconsecutive,index1part,index1interval,novelsplicingp,
		 merge_samechr_p,circularp,failedinput_root,print_m8_p,want_random_p);
  SAM_setup(add_paired_nomappers_p,paired_flag_means_concordant_p,
	    quiet_if_excessive_p,maxpaths_report,failedinput_root,fastq_format_p,hide_soft_clips_p,
	    clip_overlap_p,merge_overlap_p,sam_multiple_primaries_p,
	    force_xs_direction_p,md_lowercase_variant_p,snps_iit,chromosome_iit,genomecomp);
  Output_setup(chromosome_iit,nofailsp,failsonlyp,quiet_if_excessive_p,maxpaths_report,
	       failedinput_root,quality_shift,
	       output_sam_p,print_m8_p,invert_first_p,invert_second_p,
	       merge_samechr_p,sam_read_group_id);

  return;
}


static void
worker_cleanup () {

  if (use_only_sarray_p == false) {
    Stage1hr_cleanup();
  }

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
#ifndef LARGE_GENOMES
  if (sarray_fwd != NULL && sarray_rev != NULL) {
    if (mode == STANDARD) {
      Sarray_free(&sarray_fwd);
    } else {
      Sarray_free(&sarray_rev);
      Sarray_free(&sarray_fwd);
    }
  }
#endif
  if (genomecomp_alt != NULL) {
    Genome_free(&genomecomp_alt);
    Genome_free(&genomebits_alt);
    FREE(splicefrags_alt);	/* If genomealt == NULL, then splicefrags_alt == splicefrags_ref */
  }
  if (genomebits != NULL) {
    Genome_free(&genomebits);
  }
  if (genomecomp != NULL) {
    Genome_free(&genomecomp);
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

  if (runlength_iit != NULL) {
    FREE(runlength_divint_crosstable);
    IIT_free(&runlength_iit);
  }

  if (tally_iit != NULL) {
    FREE(tally_divint_crosstable);
    IIT_free(&tally_iit);
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

  Access_controlled_cleanup();

  return;
}


int
main (int argc, char *argv[]) {
  int nchars1 = 0, nchars2 = 0;
  bool multiple_sequences_p;
  int cmdline_status;

  char *genomesubdir, *fileroot, *dbversion;
  char **files;
  int nfiles;
#if defined(USE_MPI) && defined(USE_MPI_FILE_INPUT)
  MPI_File mpi_file_input, mpi_file_input_2;
#endif

#ifdef USE_MPI
  Master_T master;
  char **files_master;
  int nfiles_master;
  FILE *input_parser, *input2_parser;
#endif
  FILE *input, *input2;

#ifdef HAVE_ZLIB
#ifdef USE_MPI
  gzFile gzipped_master, gzipped2_master;
#endif
  gzFile gzipped, gzipped2;
#endif

#ifdef HAVE_BZLIB
#ifdef USE_MPI
  Bzip2_T bzipped_master, bzipped2_master;
#endif
  Bzip2_T bzipped, bzipped2;
#endif

  long int worker_id;

  int nread;
  int nextchar = '\0';
  double runtime;

#ifdef HAVE_PTHREAD
  int ret;
  pthread_attr_t thread_attr_join;
#ifdef USE_MPI
  pthread_attr_t thread_attr_detach;
#endif
#endif

#ifdef HAVE_SIGACTION
  struct sigaction signal_action;
#endif

  extern int optind;

#ifdef MEMUSAGE
  Mem_usage_init();
  Mem_usage_set_threadname("main");
#endif


  cmdline_status = parse_command_line(argc,argv,optind);
  argc -= optind;
  argv += optind;

  if (cmdline_status == 0) {
    /* okay to continue */
  } else if (cmdline_status == 1) {
    exit(0);
  } else {
    exit(cmdline_status);
  }

  check_compiler_assumptions();

  if (exception_raise_p == false) {
    fprintf(stderr,"Allowing signals and exceptions to pass through\n");
    Except_inactivate();
  } else {
#ifdef HAVE_SIGACTION
    signal_action.sa_handler = signal_handler;
    signal_action.sa_flags = 0;
    sigfillset(&signal_action.sa_mask); /* After first signal, block all other signals */

    /* Note: SIGKILL and SIGSTOP cannot be caught */

    sigaction(SIGABRT,&signal_action,NULL); /* abnormal termination (abort) */
    sigaction(SIGBUS,&signal_action,NULL);  /* bus error */
    sigaction(SIGFPE,&signal_action,NULL);  /* arithmetic exception */
    sigaction(SIGHUP,&signal_action,NULL);  /* hangup */
    sigaction(SIGILL,&signal_action,NULL);  /* illegal hardware instruction */
    sigaction(SIGINT,&signal_action,NULL);  /* terminal interruption (control-C) */
    sigaction(SIGPIPE,&signal_action,NULL);  /* write to pipe with no readers */
    sigaction(SIGQUIT,&signal_action,NULL);  /* terminal quit (control-backslash) */
    sigaction(SIGSEGV,&signal_action,NULL);  /* invalid memory reference */
    sigaction(SIGSYS,&signal_action,NULL);  /* invalid system call */
    sigaction(SIGTERM,&signal_action,NULL);  /* Unix kill command */
    sigaction(SIGTRAP,&signal_action,NULL);  /* hardware fault */
    sigaction(SIGXCPU,&signal_action,NULL);  /* CPU limit exceeded */
    sigaction(SIGXFSZ,&signal_action,NULL);  /* file size limit exceeded */
#endif
  }

#ifdef USE_MPI
  /* MPI_Init(&argc,&argv); */
  MPI_Init_thread(&argc,&argv,/*requested*/MPI_THREAD_MULTIPLE,&provided);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  MPI_Comm_size(MPI_COMM_WORLD,&nranks);
  MPI_Debug_setup(myid);

  nthreads0 = nthreads - 1;
  if (master_is_worker_p == false) {
    /* Default is to exclude master node from working */
    exclude_ranks[0] = 0;
    MPI_Comm_group(MPI_COMM_WORLD,&world_group);
    MPI_Group_excl(world_group,1,exclude_ranks,&workers_group);
    MPI_Comm_create(MPI_COMM_WORLD,workers_group,&workers_comm);
    MPI_Group_free(&workers_group);
    MPI_Group_free(&world_group);

  } else if (nthreads0 <= 0) {
    /* If insufficient threads, then also exclude master node from working */
    exclude_ranks[0] = 0;
    MPI_Comm_group(MPI_COMM_WORLD,&world_group);
    MPI_Group_excl(world_group,1,exclude_ranks,&workers_group);
    MPI_Comm_create(MPI_COMM_WORLD,workers_group,&workers_comm);
    MPI_Group_free(&workers_group);
    MPI_Group_free(&world_group);
    master_is_worker_p = false;

  } else {
    /* Include master rank 0 in workers group */
    MPI_Comm_group(MPI_COMM_WORLD,&world_group);
    MPI_Comm_create(MPI_COMM_WORLD,world_group,&workers_comm);
    MPI_Group_free(&world_group);
    /* master_is_worker_p = true; */
  }
  n_slave_ranks = nranks - 1;	/* Don't include master, even if it's a worker */

  if (myid == 0) {
    nthreads = nthreads0;
    fastq_format_p = open_input_streams_parser(&nextchar,&nchars1,&nchars2,
					       &files_master,&nfiles_master,&input_parser,&input2_parser,
#ifdef HAVE_ZLIB
					       &gzipped_master,&gzipped2_master,
#endif
#ifdef HAVE_BZLIB
					       &bzipped_master,&bzipped2_master,
#endif
					       gunzip_p,bunzip2_p,argc,argv);
    master = Master_new(n_slave_ranks,nextchar,nchars1,nchars2,
			input_parser,input2_parser,
#ifdef HAVE_ZLIB
			gzipped_master,gzipped2_master,
#endif
#ifdef HAVE_BZLIB
			bzipped_master,bzipped2_master,
#endif
			files_master,nfiles_master,inbuffer_nspaces,part_modulus,part_interval);
  }

  MPI_Bcast(&fastq_format_p,1,MPI_BOOL_T,/*root*/0,MPI_COMM_WORLD);
  MPI_Bcast(&nextchar,1,MPI_CHAR,/*root*/0,MPI_COMM_WORLD);

  /* If not using MPI_File, then master already has input and input2,
     and does not need mpi_file_input or mpi_file_input_2 (because of the workers_comm) */
  if (myid > 0 || master_is_worker_p == true) {
    open_input_streams_worker(&files,&nfiles,
#ifdef USE_MPI_FILE_INPUT
			      &mpi_file_input,&mpi_file_input_2,workers_comm,
#else
			      &input,&input2,
#endif
#ifdef HAVE_ZLIB
			      &gzipped,&gzipped2,
#endif
#ifdef HAVE_BZLIB
			      &bzipped,&bzipped2,
#endif
			      gunzip_p,bunzip2_p,fastq_format_p,argc,argv);

    /* Inbuffer_master_process skips to part_modulus, so workers need it set to 0 */
    Inbuffer_setup(filter_if_both_p,
#ifdef USE_MPI_FILE_INPUT
		   workers_comm,
#endif
		   /*part_modulus*/0,part_interval);

    inbuffer = Inbuffer_new(nextchar,myid,
#ifdef USE_MPI_FILE_INPUT
			    mpi_file_input,mpi_file_input_2,
#else
			    input,input2,
#endif
#ifdef HAVE_ZLIB
			    gzipped,gzipped2,
#endif
#ifdef HAVE_BZLIB
			    bzipped,bzipped2,
#endif
			    files,nfiles,inbuffer_nspaces);
  }

  Shortread_setup(acc_fieldi_start,acc_fieldi_end,force_single_end_p,filter_chastity_p,
		  allow_paired_end_mismatch_p,fastq_format_p,barcode_length,
		  invert_first_p,invert_second_p);
  multiple_sequences_p = true;


#else
  /* Non-MPI version */
  fastq_format_p = open_input_streams_parser(&nextchar,&nchars1,&nchars2,&files,&nfiles,&input,&input2,
#ifdef HAVE_ZLIB
					     &gzipped,&gzipped2,
#endif
#ifdef HAVE_BZLIB
					     &bzipped,&bzipped2,
#endif
					     gunzip_p,bunzip2_p,argc,argv);

  Inbuffer_setup(filter_if_both_p,part_modulus,part_interval);

  inbuffer = Inbuffer_new(nextchar,input,input2,
#ifdef HAVE_ZLIB
			  gzipped,gzipped2,
#endif
#ifdef HAVE_BZLIB
			  bzipped,bzipped2,
#endif
			  files,nfiles,inbuffer_nspaces);

  Shortread_setup(acc_fieldi_start,acc_fieldi_end,force_single_end_p,filter_chastity_p,
		  allow_paired_end_mismatch_p,fastq_format_p,barcode_length,
		  invert_first_p,invert_second_p);

  if ((nread = Inbuffer_fill_init(inbuffer)) > 1) {
    multiple_sequences_p = true;
  } else {
    multiple_sequences_p = false;
  }
#endif

  if (multiple_sequences_p == true) {
#if 0
    if (offsetsstrm_access != USE_ALLOCATE || genome_access != USE_ALLOCATE ||
	sarray_access != USE_ALLOCATE || lcp_access != USE_ALLOCATE) {
      fprintf(stderr,"Note: >1 sequence detected, so index files are being memory mapped.\n");
      fprintf(stderr,"  GSNAP can run slowly at first while the computer starts to accumulate\n");
      fprintf(stderr,"  pages from the hard disk into its cache.  To copy index files into RAM\n");
      fprintf(stderr,"  instead of memory mapping, use -B 3, -B 4, or -B 5, if you have enough RAM.\n");
#ifdef HAVE_PTHREAD
      fprintf(stderr,"  For more speed, also try multiple threads (-t <int>), if you have multiple processors or cores.");
#endif
      fprintf(stderr,"\n");
    }
#endif

  } else {
    /* fprintf(stderr,"Note: only 1 sequence detected.  Ignoring batch (-B) command\n"); */
    expand_offsets_p = false;
#ifdef HAVE_MMAP
    offsetsstrm_access = USE_MMAP_ONLY;
    positions_access = USE_MMAP_ONLY;
    genome_access = USE_MMAP_ONLY;
    sarray_access = USE_MMAP_ONLY;
    lcp_access = USE_MMAP_ONLY;
    guideexc_access = USE_MMAP_ONLY;
    indexij_access = USE_MMAP_ONLY;
#else
    /* No choice, since mmap is not available */
    offsetsstrm_access = USE_ALLOCATE;
    positions_access = USE_ALLOCATE;
    genome_access = USE_ALLOCATE;
    sarray_access = USE_ALLOCATE;
    lcp_access = USE_ALLOCATE;
    guideexc_access = USE_ALLOCATE;
    indexij_access = USE_ALLOCATE;
#endif
  }

  genomesubdir = Datadir_find_genomesubdir(&fileroot,&dbversion,user_genomedir,dbroot);
  FREE(dbversion);
  chromosome_iit = chromosome_iit_setup(&nchromosomes,&circular_typeint,&any_circular_p,&circularp,
					genomesubdir,fileroot);
  Outbuffer_setup(argc,argv,optind,chromosome_iit,any_circular_p,
		  nthreads,orderedp,quiet_if_excessive_p,
		  output_sam_p,sam_headers_p,sam_read_group_id,sam_read_group_name,
		  sam_read_group_library,sam_read_group_platform,
		  appendp,output_file,split_output_root,failedinput_root);

#if defined(USE_MPI) && defined(HAVE_PTHREAD)
  /* Needed for Master_parser and possibly Master_write_stdout, which never terminate */
  pthread_attr_init(&thread_attr_detach);
  if ((ret = pthread_attr_setdetachstate(&thread_attr_detach,PTHREAD_CREATE_DETACHED)) != 0) {
    fprintf(stderr,"ERROR: pthread_attr_setdetachstate returned %d\n",ret);
    exit(1);
  }
#endif

#ifdef USE_MPI
  if (myid == 0 && master_is_worker_p == false) {
    FREE(genomesubdir);
    FREE(fileroot);

    /* Master rank, which is not a worker */
    if (output_file != NULL) {
      fprintf(stderr,"Starting alignment.  Writing results to %s\n",output_file);
    } else if (split_output_root != NULL) {
      fprintf(stderr,"Starting alignment.  Writing results to %s.*\n",split_output_root);
    } else {
      fprintf(stderr,"Starting alignment\n");
    }

    stopwatch = Stopwatch_new();
    Stopwatch_start(stopwatch);

    if (split_output_root == NULL && output_file == NULL) {
      pthread_create(&write_stdout_thread_id,&thread_attr_detach,Master_write_stdout,(void *) NULL);
    }
    pthread_create(&parser_thread_id,&thread_attr_detach,Master_parser,(void *) master);
    Master_mpi_interface((void *) master); /* Can run as a normal procedure, not as a thread */

  } else {
    worker_setup(genomesubdir,fileroot);
    FREE(genomesubdir);
    FREE(fileroot);
    MPI_Barrier(workers_comm);

    outbuffer = Outbuffer_new(output_buffer_size,/*nread*/0);
    Inbuffer_set_outbuffer(inbuffer,outbuffer);
    /* MPI worker ranks continue on with creating output_thread and worker_threads below */

    if (myid == 0) {
      Inbuffer_set_master(inbuffer,master);

      if (output_file != NULL) {
	fprintf(stderr,"Starting alignment.  Writing results to %s\n",output_file);
      } else if (split_output_root != NULL) {
	fprintf(stderr,"Starting alignment.  Writing results to %s.*\n",split_output_root);
      } else {
	fprintf(stderr,"Starting alignment\n");
      }
      stopwatch = Stopwatch_new();
      Stopwatch_start(stopwatch);
    }

#else
  worker_setup(genomesubdir,fileroot);
  FREE(genomesubdir);
  FREE(fileroot);

  outbuffer = Outbuffer_new(output_buffer_size,nread);
  Inbuffer_set_outbuffer(inbuffer,outbuffer);

  if (output_file != NULL) {
    fprintf(stderr,"Starting alignment.  Writing results to %s\n",output_file);
  } else if (split_output_root != NULL) {
    fprintf(stderr,"Starting alignment.  Writing results to %s.*\n",split_output_root);
  } else {
    fprintf(stderr,"Starting alignment\n");
  }
  stopwatch = Stopwatch_new();
  Stopwatch_start(stopwatch);
#endif



#if !defined(HAVE_PTHREAD)
  /* Serial version */
  single_thread();

#else
  /* Pthreads version */
  if (nthreads == 0) {
    single_thread();

  } else if (multiple_sequences_p == false) {
    single_thread();

  } else {
    pthread_attr_init(&thread_attr_join);
    if ((ret = pthread_attr_setdetachstate(&thread_attr_join,PTHREAD_CREATE_JOINABLE)) != 0) {
      fprintf(stderr,"ERROR: pthread_attr_setdetachstate returned %d\n",ret);
      exit(1);
    }

#ifdef USE_MPI
    /* Master rank that is working or a Slave rank */
    if (myid == 0) {
      if (split_output_root == NULL && output_file == NULL) {
	pthread_create(&write_stdout_thread_id,&thread_attr_detach,Master_write_stdout,(void *) NULL);
      }
      pthread_create(&parser_thread_id,&thread_attr_detach,Master_parser,(void *) master);
      pthread_create(&mpi_interface_thread_id,&thread_attr_join,Master_mpi_interface,(void *) master);
    }
#endif

    worker_thread_ids = (pthread_t *) CALLOC(nthreads,sizeof(pthread_t));
    Except_init_pthread();
    pthread_key_create(&global_request_key,NULL);

    if (orderedp == true) {
      pthread_create(&output_thread_id,&thread_attr_join,Outbuffer_thread_ordered,
		     (void *) outbuffer);
    } else {
      pthread_create(&output_thread_id,&thread_attr_join,Outbuffer_thread_anyorder,
		     (void *) outbuffer);
    }

    for (worker_id = 0; worker_id < nthreads; worker_id++) {
      /* Need to have worker threads finish before we call Inbuffer_free() */
      pthread_create(&(worker_thread_ids[worker_id]),&thread_attr_join,worker_thread,(void *) worker_id);
    }
    
    pthread_join(output_thread_id,NULL);
    for (worker_id = 0; worker_id < nthreads; worker_id++) {
      pthread_join(worker_thread_ids[worker_id],NULL);
    }
#ifdef USE_MPI
    if (myid == 0) {
      pthread_join(mpi_interface_thread_id,NULL);
    }
#endif

    pthread_key_delete(global_request_key);
    /* Do not delete global_except_key, because worker threads might still need it */
    /* Except_term_pthread(); */

    FREE(worker_thread_ids);

  }
#endif /* HAVE_PTHREAD */

#ifdef USE_MPI
  /* MPI worker ranks finished with creating output_thread and worker_threads below */
  }
#endif


  /* Note: Shortread and Sequence procedures should close their own input files */
#ifdef USE_MPI
  if (myid == 0) {
    runtime = Stopwatch_stop(stopwatch);
    Stopwatch_free(&stopwatch);

    nread = Master_ntotal(master);
    fprintf(stderr,"Processed %u queries in %.2f seconds (%.2f queries/sec)\n",
	    nread,runtime,(double) nread/runtime);
    /* Master_free(&master); -- Master_parser thread still needs this */
  }

  if (myid > 0 || master_is_worker_p) {
    Outbuffer_free(&outbuffer);
    Inbuffer_free(&inbuffer);
  }

  Outbuffer_close_files();	/* All ranks have to close the files */

#else
  /* Single CPU or Pthreads version */
  runtime = Stopwatch_stop(stopwatch);
  Stopwatch_free(&stopwatch);

  nread = Outbuffer_nread(outbuffer);
  /* nbeyond = Outbuffer_nbeyond(outbuffer); */
  fprintf(stderr,"Processed %u queries in %.2f seconds (%.2f queries/sec)\n",
	  nread,runtime,(double) nread/runtime);
  
  Outbuffer_free(&outbuffer);
  Inbuffer_free(&inbuffer);

  Outbuffer_close_files();
#endif

  Outbuffer_cleanup();

#ifdef USE_MPI
  if (myid > 0 || master_is_worker_p == true) {
    worker_cleanup();
    MPI_Comm_free(&workers_comm);
  }
  MPI_Barrier(MPI_COMM_WORLD);	/* Make sure all processes have cleaned up */
  MPI_Finalize();
#else
  worker_cleanup();
#endif

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
  -D, --dir=directory            Genome directory.  Default (as specified by --with-gmapdb to the configure program) is\n\
                                   %s\n\
",GMAPDB);
  fprintf(stdout,"\
  -d, --db=STRING                Genome database\n\
  --use-sarray=INT               Whether to use a suffix array, which will give increased speed.\n\
                                   Allowed values: 0 (no), 1 (yes, plus GSNAP/GMAP algorithm, default),\n\
                                   or 2 (yes, and use only suffix array algorithm).\n\
                                   Note that suffix arrays will bias against SNP alleles in\n\
                                   SNP-tolerant alignment.\n\
  -k, --kmer=INT                 kmer size to use in genome database (allowed values: 16 or less)\n\
                                   If not specified, the program will find the highest available\n\
                                   kmer size in the genome database\n\
  --sampling=INT                 Sampling to use in genome database.  If not specified, the program\n\
                                   will find the smallest available sampling value in the genome database\n\
                                   within selected k-mer size\n\
  -q, --part=INT/INT             Process only the i-th out of every n sequences\n\
                                   e.g., 0/100 or 99/100 (useful for distributing jobs\n\
                                   to a computer farm).\n\
");
  fprintf(stdout,"\
  --input-buffer-size=INT        Size of input buffer (program reads this many sequences\n\
                                   at a time for efficiency) (default %d)\n\
",inbuffer_nspaces);
  fprintf(stdout,"\
  --barcode-length=INT           Amount of barcode to remove from start of read\n\
                                   (default %d)\n\
",barcode_length);
  fprintf(stdout,"\
  --orientation=STRING           Orientation of paired-end reads\n\
                                   Allowed values: FR (fwd-rev, or typical Illumina; default),\n\
                                   RF (rev-fwd, for circularized inserts), or FF (fwd-fwd, same strand)\n\
  --fastq-id-start=INT           Starting position of identifier in FASTQ header, space-delimited (>= 1)\n\
  --fastq-id-end=INT             Ending position of identifier in FASTQ header, space-delimited (>= 1)\n\
                                 Examples:\n\
                                   @HWUSI-EAS100R:6:73:941:1973#0/1\n\
                                      start=1, end=1 (default) => identifier is HWUSI-EAS100R:6:73:941:1973#0\n\
                                   @SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36\n\
                                      start=1, end=1  => identifier is SRR001666.1\n\
                                      start=2, end=2  => identifier is 071112_SLXA-EAS1_s_7:5:1:817:345\n\
                                      start=1, end=2  => identifier is SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345\n\
  --force-single-end             When multiple FASTQ files are provided on the command line, GSNAP assumes\n\
                                    they are matching paired-end files.  This flag treats each file as single-end.\n\
  --filter-chastity=STRING       Skips reads marked by the Illumina chastity program.  Expecting a string\n\
                                   after the accession having a 'Y' after the first colon, like this:\n\
                                         @accession 1:Y:0:CTTGTA\n\
                                   where the 'Y' signifies filtering by chastity.\n\
                                   Values: off (default), either, both.  For 'either', a 'Y' on either end\n\
                                   of a paired-end read will be filtered.  For 'both', a 'Y' is required\n\
                                   on both ends of a paired-end read (or on the only end of a single-end read).\n\
  --allow-pe-name-mismatch       Allows accession names of reads to mismatch in paired-end files\n\
");
#ifdef HAVE_ZLIB
  fprintf(stdout,"\
  --gunzip                       Uncompress gzipped input files\n\
");
#endif
#ifdef HAVE_BZLIB
  fprintf(stdout,"\
  --bunzip2                      Uncompress bzip2-compressed input files\n\
");
#endif
  fprintf(stdout,"\n");

  /* Computation options */
  fprintf(stdout,"Computation options\n");
  fprintf(stdout,"\
\n\
  Note: GSNAP has an ultrafast algorithm for calculating mismatches up to and including\n\
((readlength+2)/kmer - 2) (\"ultrafast mismatches\").  The program will run fastest if\n\
max-mismatches (plus suboptimal-levels) is within that value.\n\
Also, indels, especially end indels, take longer to compute, although the algorithm\n\
is still designed to be fast.\n\
\n\
");
#ifdef HAVE_MMAP
  fprintf(stdout,"\
  -B, --batch=INT                Batch mode (default = 2)\n\
                                 Mode     Offsets       Positions       Genome          Suffix array\n\
                                   0      see note      mmap            mmap            mmap\n\
                                   1      see note      mmap & preload  mmap            mmap\n\
                                   2      see note      mmap & preload  mmap & preload  mmap & preload\n\
                                   3      see note      allocate        mmap & preload  mmap & preload\n\
                      (default)    4      see note      allocate        allocate        mmap & preload\n\
                                   5      see note      allocate        allocate        allocate\n\
                           Note: For a single sequence, all data structures use mmap\n\
                           If mmap not available and allocate not chosen, then will use fileio (very slow)\n\
");
#else
  fprintf(stdout,"\
  -B, --batch=INT                Batch mode (default = 5, modes 0-4 disallowed because program configured without mmap)\n\
                                 Mode     Offsets       Positions       Genome         Suffix array\n\
                      (default)    5      see note      allocate        allocate       allocate\n\
");
#endif
  fprintf(stdout,"\
                       Note about offsets: Expansion of offsets can be controlled\n\
                       independently by the --expand-offsets flag.  However, offsets\n\
                       are accessed relatively fast in this version of GSNAP.\n\
\n\
  --use-shared-memory=INT        If 1 (default), then allocated memory is shared among all processes\n\
                                   on this node.  If 0, then each process has private allocated memory\n\
  --expand-offsets=INT           Whether to expand the genomic offsets index\n\
                                   Values: 0 (no, default), or 1 (yes).\n\
                                   Expansion gives faster alignment, but requires more memory\n\
");

  fprintf(stdout,"\
  -m, --max-mismatches=FLOAT     Maximum number of mismatches allowed (if not specified, then\n\
                                   defaults to the ultrafast level of ((readlength+index_interval-1)/kmer - 2))\n\
                                   (By default, the genome index interval is 3, but this can be changed\n \
                                   by providing a different value for -q to gmap_build when processing\n\
                                   the genome.)\n			\
                                   If specified between 0.0 and 1.0, then treated as a fraction\n\
                                   of each read length.  Otherwise, treated as an integral number\n\
                                   of mismatches (including indel and splicing penalties)\n\
                                   For RNA-Seq, you may need to increase this value slightly\n\
                                   to align reads extending past the ends of an exon.\n\
  --min-coverage=FLOAT           Minimum coverage required for an alignment.\n\
                                   If specified between 0.0 and 1.0, then treated as a fraction\n\
                                   of each read length.  Otherwise, treated as an integral number\n\
                                   of base pairs.  Default value is 0.0.\n\
  --query-unk-mismatch=INT       Whether to count unknown (N) characters in the query as a mismatch\n\
                                   (0=no (default), 1=yes)\n\
  --genome-unk-mismatch=INT      Whether to count unknown (N) characters in the genome as a mismatch\n\
                                   (0=no, 1=yes (default))\n\
");
  fprintf(stdout,"\
  --maxsearch=INT                Maximum number of alignments to find (default %d).\n\
                                   Must be larger than --npaths, which is the number to report.\n\
                                   Keeping this number large will allow for random selection among multiple alignments.\n\
                                   Reducing this number can speed up the program.\n\
",maxpaths_search);

#if 0
  fprintf(stdout,"\
  -i, --indel-penalty-middle=INT Penalty for an indel in middle of read (default %d).\n\
                                   Counts against mismatches allowed.  To find indels, make\n\
                                   indel-penalty less than or equal to max-mismatches\n\
",indel_penalty_middle);
  fprintf(stdout,"\
  -I, --indel-penalty-end=INT    Penalty for an indel at end of read (default %d).\n\
                                   Counts against mismatches allowed.  To find indels, make\n\
                                   indel-penalty less than or equal to max-mismatches\n\
",indel_penalty_end);
#endif

  fprintf(stdout,"\
  -i, --indel-penalty=INT        Penalty for an indel (default %d).\n\
                                   Counts against mismatches allowed.  To find indels, make\n\
                                   indel-penalty less than or equal to max-mismatches.\n\
                                   A value < 2 can lead to false positives at read ends\n\
",indel_penalty_middle);

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
  --indel-endlength=INT          Minimum length at end required for indel alignments (default %d)\n\
",min_indel_end_matches);
  fprintf(stdout,"\
  -y, --max-middle-insertions=INT  Maximum number of middle insertions allowed (default %d)\n\
",max_middle_insertions);
  fprintf(stdout,"\
  -z, --max-middle-deletions=INT Maximum number of middle deletions allowed (default %d)\n\
",max_middle_deletions);
  fprintf(stdout,"\
  -Y, --max-end-insertions=INT   Maximum number of end insertions allowed (default %d)\n\
",max_end_insertions);
  fprintf(stdout,"\
  -Z, --max-end-deletions=INT    Maximum number of end deletions allowed (default %d)\n\
",max_end_deletions);
  fprintf(stdout,"\
  -M, --suboptimal-levels=INT    Report suboptimal hits beyond best hit (default %d)\n\
                                   All hits with best score plus suboptimal-levels are reported\n\
",subopt_levels);
  fprintf(stdout,"\
  -a, --adapter-strip=STRING     Method for removing adapters from reads.  Currently allowed values: off, paired.\n\
                                   Default is \"off\".  To turn on, specify \"paired\", which removes adapters\n\
                                   from paired-end reads if they appear to be present.\n\
");
  fprintf(stdout,"\
  --trim-mismatch-score=INT      Score to use for mismatches when trimming at ends (default is %d;\n\
                                   to turn off trimming, specify 0).  Warning: turning trimming off\n\
                                   will give false positive mismatches at the ends of reads\n\
",trim_mismatch_score);
  fprintf(stdout,"\
  --trim-indel-score=INT         Score to use for indels when trimming at ends (default is %d;\n\
                                   to turn off trimming, specify 0).  Warning: turning trimming off\n\
                                   will give false positive indels at the ends of reads\n\
",trim_indel_score);
  fprintf(stdout,"\
  -V, --snpsdir=STRING           Directory for SNPs index files (created using snpindex) (default is\n\
                                   location of genome index files specified using -D and -d)\n \
  -v, --use-snps=STRING          Use database containing known SNPs (in <STRING>.iit, built\n\
                                   previously using snpindex) for tolerance to SNPs\n\
  --cmetdir=STRING               Directory for methylcytosine index files (created using cmetindex)\n\
                                   (default is location of genome index files specified using -D, -V, and -d)\n\
  --atoidir=STRING               Directory for A-to-I RNA editing index files (created using atoiindex)\n\
                                   (default is location of genome index files specified using -D, -V, and -d)\n\
  --mode=STRING                  Alignment mode: standard (default), cmet-stranded, cmet-nonstranded,\n\
                                    atoi-stranded, atoi-nonstranded, ttoc-stranded, or ttoc-nonstranded.\n\
                                    Non-standard modes requires you to have previously run the cmetindex\n\
                                    or atoiindex programs (which also cover the ttoc modes) on the genome\n\
");


#if 0
  fprintf(stdout,"\
  --tallydir=STRING              Directory for tally IIT file to resolve concordant multiple alignments (default is\n\
                                   location of genome index files specified using -D and -d).  Note: can\n\
                                   just give full path name to --use-tally instead.\n\
  --use-tally=STRING             Use this tally IIT file to resolve concordant multiple alignments\n\
  --runlengthdir=STRING          Directory for runlength IIT file to resolve concordant multiple alignments (default is\n\
                                   location of genome index files specified using -D and -d).  Note: can\n\
                                   just give full path name to --use-runlength instead.\n\
  --use-runlength=STRING         Use this runlength IIT file to resolve concordant multiple alignments\n\
");
#endif


#if 0
  fprintf(stdout,"\
  -2, --dibase                   Input is 2-base encoded (e.g., SOLiD), with database built\n\
                                   previously using dibaseindex)\n\
");
#endif

#ifdef HAVE_PTHREAD
  fprintf(stdout,"\
  -t, --nthreads=INT             Number of worker threads\n\
");
#else
  fprintf(stdout,"\
  -t, --nthreads=INT             Number of worker threads.  Flag is ignored in this version of GSNAP, which has pthreads disabled\n\
");
#endif

  fprintf(stdout,"\n");
  fprintf(stdout,"Options for GMAP alignment within GSNAP\n");
  fprintf(stdout,"\
  --gmap-mode=STRING             Cases to use GMAP for complex alignments containing multiple splices or indels\n\
                                 Allowed values: none, all, pairsearch, indel_knownsplice, terminal, improve\n\
                                   (or multiple values, separated by commas).\n\
                                   Default: all, i.e., pairsearch,indel_knownsplice,terminal,improve\n\
");
  fprintf(stdout,"\
  --trigger-score-for-gmap=INT   Try GMAP pairsearch on nearby genomic regions if best score (the total\n\
                                   of both ends if paired-end) exceeds this value (default %d)\n\
",trigger_score_for_gmap);
  fprintf(stdout,"\
  --gmap-min-match-length=INT    Keep GMAP hit only if it has this many consecutive matches (default %d)\n\
",gmap_min_nconsecutive);
  fprintf(stdout,"\
  --gmap-allowance=INT           Extra mismatch/indel score allowed for GMAP alignments (default %d)\n\
",gmap_allowance);
  fprintf(stdout,"\
  --max-gmap-pairsearch=INT      Perform GMAP pairsearch on nearby genomic regions up to this many\n\
                                   many candidate ends (default %d).  Requires pairsearch in --gmap-mode\n\
",max_gmap_pairsearch);
  fprintf(stdout,"\
  --max-gmap-terminal=INT        Perform GMAP terminal on nearby genomic regions up to this many\n\
                                   candidate ends (default %d).  Requires terminal in --gmap-mode\n\
",max_gmap_terminal);
  fprintf(stdout,"\
  --max-gmap-improvement=INT     Perform GMAP improvement on nearby genomic regions up to this many\n\
                                   candidate ends (default %d).  Requires improve in --gmap-mode\n\
",max_gmap_improvement);
  fprintf(stdout,"\
  --microexon-spliceprob=FLOAT   Allow microexons only if one of the splice site probabilities is\n\
                                   greater than this value (default %.2f)\n\
",microexon_spliceprob);
  fprintf(stdout,"\n");


#if 0
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
#endif


  /* Splicing options */
  fprintf(stdout,"Splicing options for DNA-Seq\n");
  fprintf(stdout,"\
  --find-dna-chimeras=INT              Look for distant splicing in DNA-Seq data (0=no (default), 1=yes)\n\
                                         Automatically inactivated for RNA-Seq data\n\
                                         if -N or -s are specified)\n\
");
  fprintf(stdout,"\n");

  /* Splicing options */
  fprintf(stdout,"Splicing options for RNA-Seq\n");
  fprintf(stdout,"\
  -N, --novelsplicing=INT              Look for novel splicing (0=no (default), 1=yes)\n\
  --splicingdir=STRING                 Directory for splicing involving known sites or known introns,\n\
                                         as specified by the -s or --use-splicing flag (default is\n\
                                         directory computed from -D and -d flags).  Note: can\n\
                                         just give full pathname to the -s flag instead.\n\
  -s, --use-splicing=STRING            Look for splicing involving known sites or known introns\n\
                                         (in <STRING>.iit), at short or long distances\n\
                                         See README instructions for the distinction between known sites\n\
                                         and known introns\n\
  --ambig-splice-noclip                For ambiguous known splicing at ends of the read, do not clip at the\n\
                                         splice site, but extend instead into the intron.  This flag makes\n\
                                         sense only if you provide the --use-splicing flag, and you are trying\n\
                                         to eliminate all soft clipping with --trim-mismatch-score=0\n\
");
  fprintf(stdout,"\
  -w, --localsplicedist=INT            Definition of local novel splicing event (default %d)\n\
",shortsplicedist);
  fprintf(stdout,"\
  --novelend-splicedist=INT            Distance to look for novel splices at the ends of reads (default %d)\n\
",shortsplicedist_novelend);
  fprintf(stdout,"\
  -e, --local-splice-penalty=INT       Penalty for a local splice (default %d).  Counts against mismatches allowed\n\
",localsplicing_penalty);
  fprintf(stdout,"\
  -E, --distant-splice-penalty=INT     Penalty for a distant splice (default %d).  A distant splice is one where\n\
                                         the intron length exceeds the value of -w, or --localsplicedist, or is an\n\
                                         inversion, scramble, or translocation between two different chromosomes\n\
                                         Counts against mismatches allowed\n\
",distantsplicing_penalty);
  fprintf(stdout,"\
  -K, --distant-splice-endlength=INT   Minimum length at end required for distant spliced alignments (default %d, min\n\
                                         allowed is the value of -k, or kmer size)\n\
",min_distantsplicing_end_matches);
  fprintf(stdout,"\
  -l, --shortend-splice-endlength=INT  Minimum length at end required for short-end spliced alignments (default %d,\n\
                                         but unless known splice sites are provided with the -s flag, GSNAP may still\n\
                                         need the end length to be the value of -k, or kmer size to find a given splice\n\
",min_shortend);
  fprintf(stdout,"\
  --distant-splice-identity=FLOAT      Minimum identity at end required for distant spliced alignments (default %.2f)\n\
",min_distantsplicing_identity);
  fprintf(stdout,"\
  --antistranded-penalty=INT           (Not currently implemented, since it leads to poor results)\n\
                                         Penalty for antistranded splicing when using stranded RNA-Seq protocols.\n\
                                         A positive value, such as 1, expects antisense on the first read\n\
                                         and sense on the second read.  Default is 0, which treats sense and antisense\n\
                                         equally well\n\
  --merge-distant-samechr              Report distant splices on the same chromosome as a single splice, if possible.\n\
                                         Will produce a single SAM line instead of two SAM lines, which is also done\n\
                                         for translocations, inversions, and scramble events\n\
");
  fprintf(stdout,"\n");


  /* Paired-end options */
  fprintf(stdout,"Options for paired-end reads\n");
  fprintf(stdout,"\
  --pairmax-dna=INT              Max total genomic length for DNA-Seq paired reads, or other reads\n\
                                   without splicing (default %d).  Used if -N or -s is not specified.\n\
",pairmax_dna);
  fprintf(stdout,"\
  --pairmax-rna=INT              Max total genomic length for RNA-Seq paired reads, or other reads\n\
                                   that could have a splice (default %d).  Used if -N or -s is specified.\n\
                                   Should probably match the value for -w, --localsplicedist.\n\
",pairmax_rna);
  fprintf(stdout,"\
  --pairexpect=INT               Expected paired-end length, used for calling splices in medial part\n\
                                   of paired-end reads (default %d).  Was turned off in previous versions, but reinstated.\n\
",expected_pairlength);
  fprintf(stdout,"\
  --pairdev=INT                  Allowable deviation from expected paired-end length, used for\n\
                                   calling splices in medial part of paired-end reads (default %d).\n\
                                   Was turned off in previous versions, but reinstated.\n\
",pairlength_deviation);
  fprintf(stdout,"\n");


  /* Quality score options */
  fprintf(stdout,"Options for quality scores\n");
  fprintf(stdout,"\
  --quality-protocol=STRING      Protocol for input quality scores.  Allowed values:\n\
                                   illumina (ASCII 64-126) (equivalent to -J 64 -j -31)\n\
                                   sanger   (ASCII 33-126) (equivalent to -J 33 -j 0)\n\
                                 Default is sanger (no quality print shift)\n\
                                 SAM output files should have quality scores in sanger protocol\n\
\n\
                                 Or you can customize this behavior with these flags:\n\
  -J, --quality-zero-score=INT   FASTQ quality scores are zero at this ASCII value\n\
                                   (default is 33 for sanger protocol; for Illumina, select 64)\n\
  -j, --quality-print-shift=INT  Shift FASTQ quality scores by this amount in output\n\
                                   (default is 0 for sanger protocol; to change Illumina input\n\
                                   to Sanger output, select -31)\n\
");

  /* Output options */
  fprintf(stdout,"Output options\n");
  fprintf(stdout,"\
  -n, --npaths=INT               Maximum number of paths to print (default %d).\n\
",maxpaths_report);
  fprintf(stdout,"\
  -Q, --quiet-if-excessive       If more than maximum number of paths are found,\n\
                                   then nothing is printed.\n\
  -O, --ordered                  Print output in same order as input (relevant\n\
                                   only if there is more than one worker thread)\n\
  --show-refdiff                 For GSNAP output in SNP-tolerant alignment, shows all differences\n\
                                   relative to the reference genome as lower case (otherwise, it shows\n\
                                   all differences relative to both the reference and alternate genome)\n\
  --clip-overlap                 For paired-end reads whose alignments overlap, clip the overlapping region.\n\
  --merge-overlap                For paired-end reads whose alignments overlap, merge the two ends into a single end (beta implementation)\n\
  --print-snps                   Print detailed information about SNPs in reads (works only if -v also selected)\n\
                                   (not fully implemented yet)\n\
  --failsonly                    Print only failed alignments, those with no results\n\
  --nofails                      Exclude printing of failed alignments\n\
");

  fprintf(stdout,"\
  -A, --format=STRING            Another format type, other than default.\n\
                                   Currently implemented: sam, m8 (BLAST tabular format)\n\
");

  fprintf(stdout,"\
  --split-output=STRING          Basename for multiple-file output, separately for nomapping,\n\
                                   halfmapping_uniq, halfmapping_mult, unpaired_uniq, unpaired_mult,\n\
                                   paired_uniq, paired_mult, concordant_uniq, and concordant_mult results\n\
  -o, --output-file=STRING       File name for a single stream of output results.\n\
  --failed-input=STRING          Print completely failed alignments as input FASTA or FASTQ format,\n\
                                    to the given file, appending .1 or .2, for paired-end data.\n\
                                    If the --split-output flag is also given, this file is generated\n\
                                    in addition to the output in the .nomapping file.\n\
  --append-output                When --split-output or --failed-input is given, this flag will append output\n\
                                    to the existing files.  Otherwise, the default is to create new files.\n\
  --order-among-best=STRING      Among alignments tied with the best score, order those alignments in this order.\n\
                                    Allowed values: genomic, random (default)\n\
");
  fprintf(stdout,"\
  --output-buffer-size=INT       Buffer size, in queries, for output thread (default %d).  When the number\n\
                                   of results to be printed exceeds this size, the worker threads are halted\n\
                                   until the backlog is cleared\n\
",output_buffer_size);
  fprintf(stdout,"\n");

  /* SAM options */
  fprintf(stdout,"Options for SAM output\n");
  fprintf(stdout,"\
  --no-sam-headers               Do not print headers beginning with '@'\n\
  --add-paired-nomappers         Add nomapper lines as needed to make all paired-end results alternate\n\
                                   between first end and second end\n\
  --paired-flag-means-concordant=INT  Whether the paired bit in the SAM flags means concordant only (1)\n\
                                 or paired plus concordant (0, default)\n\
  --sam-headers-batch=INT        Print headers only for this batch, as specified by -q\n\
  --sam-use-0M                   Insert 0M in CIGAR between adjacent insertions and deletions\n\
                                   Required by Picard, but can cause errors in other tools\n\
  --sam-multiple-primaries       Allows multiple alignments to be marked as primary if they\n\
                                   have equally good mapping scores\n\
  --force-xs-dir                 For RNA-Seq alignments, disallows XS:A:? when the sense direction\n\
                                   is unclear, and replaces this value arbitrarily with XS:A:+.\n\
                                   May be useful for some programs, such as Cufflinks, that cannot\n\
                                   handle XS:A:?.  However, if you use this flag, the reported value\n\
                                   of XS:A:+ in these cases will not be meaningful.\n\
  --md-lowercase-snp             In MD string, when known SNPs are given by the -v flag,\n\
                                   prints difference nucleotides as lower-case when they,\n\
                                   differ from reference but match a known alternate allele\n\
  --extend-soft-clips            Extends alignments through soft clipped regions\n\
  --action-if-cigar-error        Action to take if there is a disagreement between CIGAR length and sequence length\n\
                                   Allowed values: ignore, warning, noprint (default), abort\n\
  --read-group-id=STRING         Value to put into read-group id (RG-ID) field\n\
  --read-group-name=STRING       Value to put into read-group name (RG-SM) field\n\
  --read-group-library=STRING    Value to put into read-group library (RG-LB) field\n\
  --read-group-platform=STRING   Value to put into read-group library (RG-PL) field\n\
");
  fprintf(stdout,"\n");

#ifdef USE_MPI
  fprintf(stdout,"Options for MPI\n");
  fprintf(stdout,"\
  --master-is-worker=INT         Determines whether master node allocates threads for performing computation\n\
                                   in addition to coordinating input and output.  Number of worker threads\n\
                                   will be --nthreads minus 2\n\
                                   Values: 0 (no, default), 1 (yes if enough worker threads available)\n\
");
  fprintf(stdout,"\n");
#endif

  /* Help options */
  fprintf(stdout,"Help options\n");
  fprintf(stdout,"\
  --check                        Check compiler assumptions\n\
  --version                      Show version\n\
  --help                         Show this help message\n\
");
  return;
}

