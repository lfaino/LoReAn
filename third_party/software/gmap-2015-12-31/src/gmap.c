static char rcsid[] = "$Id: gmap.c 173190 2015-09-01 18:59:44Z twu $";
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

#include <signal.h>

#include "except.h"
#include "mem.h"
#include "bool.h"
#include "fopen.h"
#include "access.h"

#include "sequence.h"
#include "match.h"
#include "matchpool.h"
#include "pairpool.h"
#include "diagpool.h"
#include "cellpool.h"
#include "stopwatch.h"
#include "genome.h"
#include "genome-write.h"
#include "genome128_hr.h"	/* For Genome_hr_setup */
#include "genome_sites.h"	/* For Genome_sites_setup */
#include "compress-write.h"
#include "maxent_hr.h"		/* For Maxent_hr_setup */
#include "stage1.h"
#include "gregion.h"
#ifdef PMAP
#include "oligoindex_pmap.h"
#else
#include "oligoindex_hr.h"	/* For Oligoindex_hr_setup */
#endif
#include "stage2.h"
#include "splicestringpool.h"
#include "splicetrie_build.h"
#include "dynprog.h"
#include "dynprog_single.h"
#include "dynprog_genome.h"
#include "dynprog_end.h"
#include "pair.h"
#include "stage3.h"
#include "comp.h"
#include "chimera.h"
#ifdef PMAP
#include "oligop.h"		/* For Oligop_setup */
#include "backtranslation.h"
#else
#include "oligo.h"		/* For Oligo_setup */
#endif
#include "indexdb.h"
#include "result.h"
#include "request.h"
#include "intlist.h"
#include "list.h"
#include "listdef.h"
#include "iit-read-univ.h"
#include "iit-read.h"
#include "datadir.h"

#include "filestring.h"
#include "output.h"
#include "inbuffer.h"
#include "outbuffer.h"

#include "getopt.h"


#define MAX_QUERYLENGTH_FOR_ALLOC    100000
#define MAX_GENOMICLENGTH_FOR_ALLOC 1000000


#define POSSIBLE_OLIGOS 65536	/* 4^8 */
#define MAX_OLIGODEPTH 3.0
#define MAX_BADOLIGOS 0.30	/* Setting to 1.0 effectively turns this check off */
#define MAX_REPOLIGOS 0.40	/* Setting to 1.0 effectively turns this check off */

#define MAX_CHIMERA_ITER 3
#define CHIMERA_PENALTY 30	/* A small value for chimera_margin will reduce this  */
#define CHIMERA_IDENTITY 0.98
#define CHIMERA_PVALUE 0.01
#define CHIMERA_FVALUE 6.634897	/* qnorm(CHIMERA_PVALUE/2)^2 */
#define CHIMERA_SLOP 90	/* in nucleotides */

#define MIN_MATCHES 20


#define MAX_NALIGNMENTS 10


/* #define EXTRACT_GENOMICSEG 1 */


/* MPI Processing */
#ifdef DEBUGM
#define debugm(x) x
#else
#define debugm(x)
#endif


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Chimera detection */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Chimera detection, details */
#ifdef DEBUG2A
#define debug2a(x) x
#else
#define debug2a(x)
#endif

/* stage3list_remove_duplicates */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif



/************************************************************************
 *   Global variables
 ************************************************************************/

static Univ_IIT_T chromosome_iit = NULL;
static Univcoord_T genomelength;
static int circular_typeint = -1;
static int nchromosomes;
static bool *circularp = NULL;
static bool any_circular_p;
static Univ_IIT_T contig_iit = NULL;
static Genome_T genomecomp = NULL;
static Genome_T genomecomp_alt = NULL;
static Genome_T genomebits = NULL;
static Genome_T genomebits_alt = NULL;
static Genomecomp_T *genomecomp_blocks = NULL;
static Genomecomp_T *genomebits_blocks = NULL;

#ifdef PMAP
static Alphabet_T required_alphabet = AA0;
static Alphabet_T alphabet = AA20; /* Initialize in case we have a usersegment */
static int alphabet_size = 20;	   /* Initialize in case we have a usersegment */
static Width_T index1part_aa = 7;
#else
static Width_T index1part;
#endif

static Indexdb_T indexdb_fwd = NULL;
static Indexdb_T indexdb_rev = NULL;

static Width_T required_index1part = 0;
static Width_T index1interval;
static Width_T required_index1interval = 0;

static IIT_T altstrain_iit = NULL;

/* Cmet and AtoI */
static char *user_cmetdir = NULL;
static char *user_atoidir = NULL;
static Mode_T mode = STANDARD;


static char *user_snpsdir = NULL;
static char *snps_root = (char *) NULL;
static IIT_T map_iit = NULL;
static int *map_divint_crosstable = NULL;

#ifdef PMAP
#if 0
static Width_T minindexsize = 3;	/* In stage 2; in aa */
static Width_T maxindexsize = 6;	/* In stage 2; in aa */
#endif
/* Now controlled by defect_rate */
static int maxpeelback = 20;	/* Needs to be at least indexsize
				   because stage 2 jumps by indexsize.
				   Also should exceed length of
				   repeated nucleotides (e.g., a
				   string of consecutive T's) */
#else
/* Making minindexsize too small can lead to spurious exons in stage 2 */
/* FOOBAR */
#if 0
static Width_T minindexsize = 8;	/* In stage 2; in nt.  Used if sampling required in stage 1. */
static Width_T maxindexsize = 8;	/* In stage 2; in nt */
#endif
static int maxpeelback = 20;	/* Needs to be at least indexsize
				   because stage 2 jumps by indexsize.
				   Also should exceed length of
				   repeated nucleotides (e.g., a
				   string of consecutive T's) */
#endif
static int maxpeelback_distalmedial = 100; /* Needs to be longer to fix bad end exons */

/* static int stuttercycles = 2; */
static int stutterhits = 3;
static int sufflookback = 60;
static int nsufflookback = 5;

#if 0
static int maxoligohits = 400; /* Must be smaller than ALLOC in oligoindex.c */
#endif
static int nullgap = 600;
static int extramaterial_end = 10;
static int extramaterial_paired = 8; /* Should be at least indexsize in nt */
static int extraband_single = 6; /* This is in addition to length2 -
				    length1.  If onesidegap is true in
				    dynprog.c, then this is equivalent
				    to extraband_single of 0.  Needs
				    to be > 0 to handle default
				    close_indels_mode. */
static int extraband_end = 6; /* Was 6.  Shouldn't differ from 0, since onesidegapp is true?
				 This is only on both sides of main diagonal */
static int extraband_paired = 14; /* This is in addition to length2 - length1 */
static int minendexon = 9;

static Stopwatch_T stopwatch = NULL;


/************************************************************************
 *   Program options
 ************************************************************************/

/* Input options */
static char *user_genomedir = NULL;
static char *dbroot = NULL;
static char *dbversion = NULL;
static char *user_genomicseg = NULL;
static bool user_selfalign_p = false;
static bool user_pairalign_p = false;
static char *user_cmdline = NULL;
static Sequence_T global_usersegment = NULL;
static int part_modulus = 0;
static int part_interval = 1;

/* Compute options */
static int min_matches;

static bool sharedp = true;
static Access_mode_T offsetsstrm_access = USE_ALLOCATE;
static bool expand_offsets_p = false;

#ifdef HAVE_MMAP
static Access_mode_T positions_access = USE_MMAP_PRELOAD;
static Access_mode_T genome_access = USE_MMAP_PRELOAD;
#else
static Access_mode_T positions_access = USE_ALLOCATE;
static Access_mode_T genome_access = USE_ALLOCATE;
#endif

static int min_intronlength = 9;
static int max_deletionlength = 50;
static int maxtotallen_bound = 2400000;
static int maxintronlen = 200000; /* Was used previously in stage 1.  Now used only in stage 2 and Stage3_mergeable. */
static int maxextension = 1000000; /* Used in stage 1.  Not adjustable by user */
static int chimera_margin = 30;	/* Useful for finding readthroughs */
static int index1interval = 3; /* Stage 1 interval if user provides a genomic segment */
static char *referencefile = NULL;

#if 0
#ifndef PMAP
static bool literalrefp = false;
#endif
#endif

#ifdef USE_MPI
static int nprocs, n_worker_procs, proci, myid;
#endif


static bool altstrainp = false;
#ifdef HAVE_PTHREAD
static pthread_t output_thread_id, *worker_thread_ids;
static pthread_key_t global_request_key;
static int nworkers = 1;	/* (int) sysconf(_SC_NPROCESSORS_ONLN) */
#else
static int nworkers = 0;	/* (int) sysconf(_SC_NPROCESSORS_ONLN) */
#endif
#ifndef PMAP
static bool prune_poor_p = false;
static bool prune_repetitive_p = false;
#endif
static int canonical_mode = 1;
static bool cross_species_p = false;
static int homopolymerp = false;

static char *user_chrsubsetname = NULL;
static Univcoord_T chrsubset_start = 0;
static Univcoord_T chrsubset_end = -1;

static int close_indels_mode = +1;
static double microexon_spliceprob = 0.95;
static int suboptimal_score_start = -1; /* Determined by simulations to have minimal effect */
static int suboptimal_score_end = 3; /* Determined by simulations to have diminishing returns above 3 */

static int trim_mismatch_score = -3;
static int trim_indel_score = -2; /* was -4 */


/* Output options */
static unsigned int output_buffer_size = 1000;
static Printtype_T printtype = SIMPLE;
static bool exception_raise_p = true;
static bool debug_graphic_p = false;
static bool stage1debug = false;
static bool diag_debug = false;
static Stage3debug_T stage3debug = NO_STAGE3DEBUG;
static bool timingp = false;
static bool checkp = false;
static int maxpaths_report = 5;	/* 0 means 1 if nonchimeric, 2 if chimeric */
static bool quiet_if_excessive_p = false;
static int suboptimal_score = 1000000;
static bool require_splicedir_p = false;


/* GFF3 */
static bool gff3_separators_p = true;

/* SAM */
#ifndef PMAP
static bool sam_paired_p = false;
static bool user_quality_shift = false;
static int quality_shift = 0;
static bool sam_headers_p = true;
static char *sam_read_group_id = NULL;
static char *sam_read_group_name = NULL;
static char *sam_read_group_library = NULL;
static char *sam_read_group_platform = NULL;
#endif
static bool sam_insert_0M_p = false;

static bool orderedp = false;
static bool failsonlyp = false;
static bool nofailsp = false;
static bool checksump = false;
static int chimera_overlap = 0;
static bool force_xs_direction_p = false;
static bool md_lowercase_variant_p = false;

/* Map file options */
static char *user_mapdir = NULL;
static char *map_iitfile = NULL;
static bool map_exons_p = false;
static bool map_bothstrands_p = false;
static bool print_comment_p = false;
static int nflanking = 0;

/* Alignment options */
static bool fulllengthp = false;
static int cds_startpos = -1;
static bool truncatep = false;
static int sense_try = 0;		/* both */
static int sense_filter = 0;		/* both */
static double min_trimmed_coverage = 0.0;
static double min_identity = 0.0;
static bool strictp = true;
static int proteinmode = 1;
static bool uncompressedp = false;
static bool nointronlenp = false;
static int invertmode = 0;
static int ngap = 3;
static int wraplength = 50;


/* Splicing IIT */
static bool novelsplicingp = true; /* Can be disabled with --nosplicing flag */
static bool knownsplicingp = false;
static bool distances_observed_p = false;
static Chrpos_T shortsplicedist = 2000000;
static int min_extra_end;		      /* If knownsplicing, then equals shortsplicedist */
static char *user_splicingdir = (char *) NULL;
static char *splicing_file = (char *) NULL;
static IIT_T splicing_iit = NULL;
static bool amb_closest_p = false;

static int donor_typeint = -1;		/* for splicing_iit */
static int acceptor_typeint = -1;	/* for splicing_iit */

static int *splicing_divint_crosstable = NULL;
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


/* Input/output */
static char *split_output_root = NULL;
static char *failedinput_root = NULL;
static bool appendp = false;
static Inbuffer_T inbuffer = NULL;
static Outbuffer_T outbuffer = NULL;
static unsigned int inbuffer_nspaces = 1000;


#ifdef PMAP
/* Used alphabetically: 01235789ABbCcDdEefGgHIiKkLlMmNnOoPQRSstuVvwXxYZ */
#else
/* Used alphabetically: 01235789AaBbCcDdEeFfGgHIijKkLlMmNnOoPpQRSsTtuVvwXxYZ */
#endif

static struct option long_options[] = {
  /* Input options */
  {"dir", required_argument, 0, 'D'},	/* user_genomedir */
  {"db", required_argument, 0, 'd'}, /* dbroot */
#ifdef PMAP
  {"alphabet", required_argument, 0, 'a'}, /* required_alphabet */
#endif
  {"kmer", required_argument, 0, 'k'}, /* required_index1part, index1part */
  {"sampling", required_argument, 0, 0}, /* required_nterval, index1interval */
  {"genomefull", no_argument, 0, 'G'}, /* uncompressedp */
  {"gseg", required_argument, 0, 'g'}, /* user_genomicseg */
  {"selfalign", no_argument, 0, '1'}, /* user_selfalign_p */
  {"pairalign", no_argument, 0, '2'}, /* user_pairalign_p */
  {"cmdline", required_argument, 0, 0}, /* user_cmdline */
  {"part", required_argument, 0, 'q'}, /* part_modulus, part_interval */
  {"input-buffer-size", required_argument, 0, 0}, /* inbuffer_nspaces */

  /* Compute options */
#ifdef HAVE_MMAP
  {"batch", required_argument, 0, 'B'}, /* offsetsstrm_access, positions_access, genome_access */
#endif
  {"expand-offsets", required_argument, 0, 0}, /* expand_offsets_p */
  {"min-intronlength", required_argument, 0, 0}, /* min_intronlength */
  {"intronlength", required_argument, 0, 'K'}, /* maxintronlen */
  {"totallength", required_argument, 0, 'L'}, /* maxtotallen_bound */
  {"chimera-margin", required_argument, 0, 'x'}, /* chimera_margin */
  {"no-chimeras", no_argument, 0, 0},		 /* chimera_margin */
#if 0
  {"reference", required_argument, 0, 'w'}, /* referencefile */
#else
  {"localsplicedist", required_argument, 0, 'w'}, /* shortsplicedist */
#endif

  {"nthreads", required_argument, 0, 't'}, /* nworkers */
  {"splicingdir", required_argument, 0, 0}, /* user_splicingdir */
  {"nosplicing", no_argument, 0, 0},	    /* novelsplicingp */
  {"use-splicing", required_argument, 0, 's'}, /* splicing_iit, knownsplicingp (was previously altstrainp) */
  {"chrsubset", required_argument, 0, 'c'}, /* user_chrsubsetname */
  {"trimendexons", required_argument, 0, 'H'}, /* minendexon */
  {"canonical-mode", required_argument, 0, 0}, /* canonical_mode */
  {"cross-species", no_argument, 0, 0}, /* cross_species_p */
  {"homopolymer", no_argument, 0, 0},	/* homopolymerp */
#ifndef PMAP
  {"prunelevel", required_argument, 0, 'p'}, /* prune_poor_p, prune_repetitive_p */
#endif
  {"allow-close-indels", required_argument, 0, 0}, /* close_indels_mode, extraband_single */
  {"microexon-spliceprob", required_argument, 0, 0}, /* microexon_spliceprob */
  {"stage2-start", required_argument, 0, 0},	     /* suboptimal_score_start */
  {"stage2-end", required_argument, 0, 0},	     /* suboptimal_score_end */

  {"cmetdir", required_argument, 0, 0}, /* user_cmetdir */
  {"atoidir", required_argument, 0, 0}, /* user_atoidir */
  {"mode", required_argument, 0, 0}, /* mode */

  /* Output options */
  {"output-buffer-size", required_argument, 0, 0}, /* output_buffer_size */
  {"summary", no_argument, 0, 'S'}, /* printtype */
  {"align", no_argument, 0, 'A'}, /* printtype */
  {"continuous", no_argument, 0, '3'}, /* printtype */
  {"continuous-by-exon", no_argument, 0, '4'}, /* printtype */
  {"noexceptions", no_argument, 0, '0'}, /* exception_raise_p */
  {"graphic", no_argument, 0, '6'}, /* debug_graphic_p */
  {"stage3debug", required_argument, 0, '8'}, /* stage3debug */
  {"diagnostic", no_argument, 0, '9'}, /* checkp */
  {"npaths", required_argument, 0, 'n'}, /* maxpaths_report */
#if 0
  {"quiet-if-excessive", no_argument, 0, 0}, /* quiet_if_excessive_p */
#endif
  {"format", required_argument, 0, 'f'}, /* printtype */
  {"failsonly", no_argument, 0, 0}, /* failsonlyp */
  {"nofails", no_argument, 0, 0}, /* nofailsp */
  {"split-output", required_argument, 0, 0}, /* split_output_root */
  {"failed-input", required_argument, 0, 0}, /* failedinput_root */
  {"append-output", no_argument, 0, 0},	     /* appendp */
  {"suboptimal-score", required_argument, 0, 0}, /* suboptimal_score */
  {"require-splicedir", no_argument, 0, 0}, /* require_splicedir_p */

  {"gff3-add-separators", required_argument, 0, 0}, /* gff3_separators_p */

#ifndef PMAP
  {"quality-protocol", required_argument, 0, 0}, /* quality_shift */
  {"quality-print-shift", required_argument, 0, 'j'}, /* quality_shift */
  {"no-sam-headers", no_argument, 0, 0},	/* sam_headers_p */
  {"sam-use-0M", no_argument, 0, 0},		/* sam_insert_0M_p */
  {"read-group-id", required_argument, 0, 0},	/* sam_read_group_id */
  {"read-group-name", required_argument, 0, 0},	/* sam_read_group_name */
  {"read-group-library", required_argument, 0, 0}, /* sam_read_group_library */
  {"read-group-platform", required_argument, 0, 0}, /* sam_read_group_platform */
  {"force-xs-dir", no_argument, 0, 0},		    /* force_xs_direction_p */
  {"md-lowercase-snp", no_argument, 0, 0},	    /* md_lowercase_variant_p */
#endif

  {"compress", no_argument, 0, 'Z'}, /* printtype */
  {"ordered", no_argument, 0, 'O'}, /* orderedp */
  {"md5", no_argument, 0, '5'}, /* checksump */
  {"chimera-overlap", required_argument, 0, 'o'}, /* chimera_overlap */
  {"snpsdir", required_argument, 0, 'V'},   /* user_snpsdir */
  {"use-snps", required_argument, 0, 'v'}, /* snps_root */

  /* Map file options */
  {"mapdir", required_argument, 0, 'M'}, /* user_mapdir */
  {"map", required_argument, 0, 'm'},	/* map_iitfile */
  {"mapexons", no_argument, 0, 'e'}, /* map_exons_p */
  {"mapboth", no_argument, 0, 'b'}, /* map_bothstrands_p */
  {"nflanking", required_argument, 0, 'u'}, /* nflanking */
  {"print-comment", no_argument, 0, 0},	    /* print_comment_p */

  /* Alignment options */
  {"exons", required_argument, 0, 'E'}, /* printtype */
#ifdef PMAP
  {"protein_gen", no_argument, 0, 'P'}, /* printtype */
  {"nucleotide", no_argument, 0, 'Q'}, /* printtype */
#else
  {"protein_dna", no_argument, 0, 'P'}, /* printtype */
  {"protein_gen", no_argument, 0, 'Q'}, /* printtype */
  {"fulllength", no_argument, 0, 'F'}, /* fulllengthp */
  {"cdsstart", required_argument, 0, 'a'}, /* cds_startpos */
  {"truncate", no_argument, 0, 'T'}, /* truncatep */
  {"direction", required_argument, 0, 'z'}, /* sense_try, sense_filter */
#endif
  {"tolerant", no_argument, 0, 'Y'}, /* strictp */
  {"nolengths", no_argument, 0, 'N'},	/* nointronlenp */
  {"invertmode", required_argument, 0, 'I'}, /* invertmode */
  {"introngap", required_argument, 0, 'i'}, /* ngap */
  {"wraplength", required_argument, 0, 'l'}, /* wraplength */
  
  /* Filtering options */
  {"min-trimmed-coverage", required_argument, 0, 0}, /* min_trimmed_coverage */
  {"min-identity", required_argument, 0, 0},	/* min_identity */

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
#ifdef PMAP
  fprintf(stdout,"PMAP: Protein Mapping and Alignment Program\n");
#else
  fprintf(stdout,"GMAP: Genomic Mapping and Alignment Program\n");
#endif
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


#ifdef PMAP
  fprintf(stdout,"Stage 1 index size: %d aa\n",index1part_aa);
#endif
  fprintf(stdout,"Sizes: off_t (%d), size_t (%d), unsigned int (%d), long int (%d), long long int (%d)\n",
	  (int) sizeof(off_t),(int) sizeof(size_t),(int) sizeof(unsigned int),(int) sizeof(long int),(int) sizeof(long long int));
  fprintf(stdout,"Default gmap directory (compiled): %s\n",GMAPDB);
  genomedir = Datadir_find_genomedir(/*user_genomedir*/NULL);
  fprintf(stdout,"Default gmap directory (environment): %s\n",genomedir);
  FREE(genomedir);
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


/* Call before Stage1_compute */
static Diagnostic_T
evaluate_query (bool *poorp, bool *repetitivep, char *queryuc_ptr, int querylength,
		Oligoindex_T oligoindex) {
  Diagnostic_T diagnostic;

  diagnostic = Diagnostic_new();

#ifdef PMAP
  Oligoindex_set_inquery(&diagnostic->query_badoligos,&diagnostic->query_repoligos,
			 &diagnostic->query_trimoligos,&diagnostic->query_trim_start,
			 &diagnostic->query_trim_end,oligoindex,queryuc_ptr,
			 querylength);
  *poorp = false;
  *repetitivep = false;
#else
  diagnostic->query_oligodepth = 
    Oligoindex_set_inquery(&diagnostic->query_badoligos,&diagnostic->query_repoligos,
			   &diagnostic->query_trimoligos,&diagnostic->query_trim_start,
			   &diagnostic->query_trim_end,oligoindex,queryuc_ptr,
			   /*querystart*/0,/*queryend*/querylength,/*trimp*/true);

  debug2(printf("query_trimoligos %d, fraction badoligos %f = %d/%d, oligodepth %f, fraction repoligos %f = %d/%d\n",
		diagnostic->query_trimoligos,
		(double) diagnostic->query_badoligos/(double) diagnostic->query_trimoligos,
		diagnostic->query_badoligos,diagnostic->query_trimoligos,
		diagnostic->query_oligodepth,
		(double) diagnostic->query_repoligos/(double) diagnostic->query_trimoligos,
		diagnostic->query_repoligos,diagnostic->query_trimoligos));

  if (diagnostic->query_trimoligos == 0) {
    *poorp = true;
  } else if (((double) diagnostic->query_badoligos/(double) diagnostic->query_trimoligos > MAX_BADOLIGOS) ||
	     (diagnostic->query_trim_end - diagnostic->query_trim_start < 80 && diagnostic->query_badoligos > 0)) {
    *poorp = true;
  } else {
    *poorp = false;
  }

  if (diagnostic->query_trimoligos == 0) {
    *repetitivep = false;
  } else if (diagnostic->query_oligodepth > MAX_OLIGODEPTH || 
	     (double) diagnostic->query_repoligos/(double) diagnostic->query_trimoligos > MAX_REPOLIGOS) {
    *repetitivep = true;
  } else {
    *repetitivep = false;
  }
#endif

  return diagnostic;
}




static Stage3_T *
stage3array_from_list (int *npaths, int *first_absmq, int *second_absmq, List_T stage3list,
		       bool mergedp, bool chimerap, bool remove_overlaps_p) {
  Stage3_T *array1, *array0, x, y;
  bool *eliminate;
  int norig, i, j;
  int threshold_score;

  debug2(printf("Entering stage3array_from_list\n"));
  /* Stage3_recompute_goodness(stage3list); -- No longer necessary */
  Stage3_compute_mapq(stage3list);

  if ((norig = List_length(stage3list)) == 0) {
    *first_absmq = 0;
    *second_absmq = 0;
    return (Stage3_T *) NULL;

  } else if (mergedp == true) {
    debug2(printf("mergedp is true\n"));
    array0 = (Stage3_T *) List_to_array(stage3list,NULL);
    List_free(&stage3list);
    *first_absmq = 0;
    *second_absmq = 0;
    *npaths = norig;
    return array0;

  } else if (chimerap == true) {
    debug2(printf("chimerap is true\n"));
    array0 = (Stage3_T *) List_to_array(stage3list,NULL);
    List_free(&stage3list);
    *first_absmq = Stage3_absmq_score(array0[0]);
    if (norig <= 2) {
      *second_absmq = 0;
    } else {
      qsort(&(array0[2]),norig-2,sizeof(Stage3_T),Stage3_cmp);
      *second_absmq = Stage3_absmq_score(array0[2]);
    }
    *npaths = norig;
    return array0;

  } else if (remove_overlaps_p == false) {
    array0 = (Stage3_T *) List_to_array(stage3list,NULL);
    List_free(&stage3list);
    qsort(array0,norig,sizeof(Stage3_T),Stage3_cmp);

    threshold_score = Stage3_goodness(array0[0]) - suboptimal_score;
    i = 1;
    while (i < norig && Stage3_goodness(array0[i]) >= threshold_score) {
      i++;
    }
    *npaths = i;

    *first_absmq = Stage3_absmq_score(array0[0]);
    if (*npaths < 2) {
      *second_absmq = 0;
    } else {
      *second_absmq = Stage3_absmq_score(array0[1]);
    }

    return array0;

  } else {
    eliminate = (bool *) CALLOCA(norig,sizeof(bool));

    /* Initial sort to remove subsumed alignments */
    array0 = (Stage3_T *) MALLOCA(norig * sizeof(Stage3_T));
    List_fill_array_and_free((void **) array0,&stage3list);
    qsort(array0,norig,sizeof(Stage3_T),Stage3_cmp);
    for (i = 0; i < norig; i++) {
      x = array0[i];
      debug(printf("%d: chr %d:%u..%u, goodness %d, matches %d, npairs %d\n",
		   i,Stage3_chrnum(x),Stage3_chrstart(x),Stage3_chrend(x),Stage3_goodness(x),Stage3_matches(x),Stage3_npairs(x)));
      for (j = i+1; j < norig; j++) {
	y = array0[j];
	if (Stage3_overlap(x,y)) {
	  eliminate[j] = true;
	}
      }
    }

    *npaths = 0;
    for (i = 0; i < norig; i++) {
      if (eliminate[i] == false) {
	(*npaths)++;
      }
    }

    array1 = (Stage3_T *) MALLOC_OUT((*npaths) * sizeof(Stage3_T)); /* Return value */
    j = 0;
    for (i = 0; i < norig; i++) {
      x = array0[i];
      if (eliminate[i] == true) {
	Stage3_free(&x);
      } else {
	array1[j++] = x;
      }
    }
    FREEA(array0);
    FREEA(eliminate);

    threshold_score = Stage3_goodness(array1[0]) - suboptimal_score;
    i = 1;
    while (i < *npaths && Stage3_goodness(array1[i]) >= threshold_score) {
      i++;
    }
    *npaths = i;

    *first_absmq = Stage3_absmq_score(array1[0]);
    if (*npaths < 2) {
      *second_absmq = 0;
    } else {
      *second_absmq = Stage3_absmq_score(array1[1]);
    }
    return array1;
  }
}


static List_T
update_stage3list (List_T stage3list, Sequence_T queryseq,
#ifdef PMAP
		   Sequence_T queryntseq,
#endif
		   Sequence_T queryuc, Stage2_alloc_T stage2_alloc,
		   Oligoindex_array_T oligoindices_major, Oligoindex_array_T oligoindices_minor,
		   Pairpool_T pairpool, Diagpool_T diagpool, Cellpool_T cellpool, int straintype, char *strain,
		   Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
		   Chrpos_T chrstart, Chrpos_T chrend, bool watsonp, int genestrand,
		   Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		   Stopwatch_T worker_stopwatch) {
  int stage2_source, stage2_indexsize;

#ifdef PMAP
  Sequence_T genomicuc = NULL;
  char *genomicseg_ptr = NULL, *genomicuc_ptr = NULL;
#elif defined(EXTRACT_GENOMICSEG)
  Sequence_T genomicuc = NULL;
#endif
  List_T all_stage2results, p;
  Stage2_T stage2;
  Stage3_T stage3;

  struct Pair_T *pairarray;
  List_T pairs;
  int goodness;
  int npairs, cdna_direction, matches, unknowns, mismatches, qopens, qindels, topens, tindels,
    ncanonical, nsemicanonical, nnoncanonical;
  int sensedir;
  int nmatches_posttrim, max_match_length, ambig_end_length_5, ambig_end_length_3;
  Splicetype_T ambig_splicetype_5, ambig_splicetype_3;
  double ambig_prob_5, ambig_prob_3;
  double min_splice_prob;
  double stage3_runtime;
#ifdef PMAP
  int subseq_offset;
#endif


#ifdef PMAP_OLD
  /* Previously used for PMAP */
  if (user_genomicseg == NULL && uncompressedp == false && straintype == 0) {
    genomicuc = Sequence_alias(genomicseg);
  } else {
    genomicuc = Sequence_uppercase(genomicseg);
  }
  genomicseg_ptr = Sequence_fullpointer(genomicseg);
  genomicuc_ptr = Sequence_fullpointer(genomicuc);
#elif defined(EXTRACT_GENOMICSEG)
  if (user_genomicseg == NULL && uncompressedp == false && straintype == 0) {
    genomicuc = Sequence_alias(genomicseg);
  } else {
    genomicuc = Sequence_uppercase(genomicseg);
  }
  genomicseg_ptr = Sequence_fullpointer(genomicseg);
  genomicuc_ptr = Sequence_fullpointer(genomicuc);
#endif

#if 0
  if (canonical_mode == 0) {
    do_final_p = false;
  } else if (canonical_mode == 1) {
    do_final_p = true;
  } else if (lowidentityp == false) {
    do_final_p = false;
  } else {
    do_final_p = true;
  }
#endif

  debug2(printf("Beginning Stage2_compute with chrstart %u and chrend %u and query_subseq_offset %d\n",
		chrstart,chrend,Sequence_subseq_offset(queryseq)));
  all_stage2results = Stage2_compute(&stage2_source,&stage2_indexsize,
				     Sequence_trimpointer(queryseq),Sequence_trimpointer(queryuc),
				     Sequence_trimlength(queryseq),/*query_offset*/0,
				     chrstart,chrend,chroffset,chrhigh,/*plusp*/watsonp,genestrand,
				     stage2_alloc,oligoindices_major,/*proceed_pctcoverage*/0.3,
				     pairpool,diagpool,cellpool,
				     /*localp*/true,/*skip_repetitive_p*/true,
				     /*favor_right_p*/false,/*max_nalignments*/MAX_NALIGNMENTS,debug_graphic_p,
				     worker_stopwatch,diag_debug);

  debug(printf("End of Stage2_compute\n"));

  for (p = all_stage2results; p != NULL; p = List_next(p)) {
    stage2 = (Stage2_T) List_head(p);

    Stopwatch_start(worker_stopwatch);
#ifdef PMAP
    subseq_offset = Sequence_subseq_offset(queryseq); /* in nucleotides */
#endif
    pairarray = Stage3_compute(&pairs,&npairs,&goodness,&cdna_direction,&sensedir,
			       &matches,&nmatches_posttrim,&max_match_length,
			       &ambig_end_length_5,&ambig_end_length_3,
			       &ambig_splicetype_5,&ambig_splicetype_3,
			       &ambig_prob_5,&ambig_prob_3,&unknowns,&mismatches,&qopens,&qindels,&topens,&tindels,
			       &ncanonical,&nsemicanonical,&nnoncanonical,&min_splice_prob,
			       Stage2_middle(stage2),Stage2_all_starts(stage2),Stage2_all_ends(stage2),
#ifdef PMAP
			       /*queryaaseq_ptr*/Sequence_fullpointer(queryseq),
			       /*queryseq_ptr*/Sequence_subseq_pointer(queryntseq,subseq_offset),
			       /*queryuc_ptr*/Sequence_subseq_pointer(queryntseq,subseq_offset),
			       /*querylength*/Sequence_subseq_length(queryntseq,subseq_offset),
			       /*skiplength*/Sequence_skiplength(queryntseq),
			       /*query_subseq_offset*/subseq_offset,
#else
			       /*queryseq_ptr*/Sequence_fullpointer(queryseq),
			       /*queryuc_ptr*/Sequence_fullpointer(queryuc),
			       /*querylength*/Sequence_fulllength(queryseq),
			       /*skiplength*/Sequence_skiplength(queryseq),
			       /*query_subseq_offset*/Sequence_subseq_offset(queryseq),
#endif
			       chrnum,chroffset,chrhigh,
			       /*knownsplice_limit_low*/0U,/*knownsplice_limit_high*/-1U,
			       watsonp,genestrand,/*jump_late_p*/watsonp ? false : true,

			       maxpeelback,pairpool,dynprogL,dynprogM,dynprogR,
			       sense_try,sense_filter,oligoindices_minor,diagpool,cellpool);
    stage3_runtime = Stopwatch_stop(worker_stopwatch);
    if (pairarray == NULL) {
      /* Skip */
    } else if (matches < min_matches) {
      FREE_OUT(pairarray);
    } else if ((stage3 = Stage3_new(pairarray,pairs,npairs,goodness,cdna_direction,sensedir,
				    stage2_source,stage2_indexsize,matches,unknowns,mismatches,
				    qopens,qindels,topens,tindels,ncanonical,nsemicanonical,nnoncanonical,
				    chrnum,chroffset,chrhigh,chrlength,watsonp,
				    /*querylength*/Sequence_fulllength(queryseq),
				    /*skiplength*/Sequence_skiplength(queryseq),
				    /*trimlength*/Sequence_trimlength(queryseq),
				    stage3_runtime,straintype,strain,altstrain_iit)) != NULL) {
      stage3list = List_push(stage3list,(void *) stage3);
    }

    Stage2_free(&stage2);
  }

  List_free(&all_stage2results);

#ifdef PMAP_OLD
  Sequence_free(&genomicuc);
#elif defined(EXTRACT_GENOMICSEG)
  Sequence_free(&genomicuc);
#endif

  return stage3list;
}


#if 0
/* This code is duplicated in get-genome.c */
static int
index_compare (const void *a, const void *b) {
  int index1 = * (int *) a;
  int index2 = * (int *) b;
  int type1, type2;
  Chrpos_T pos1, pos2;

  type1 = Interval_type(IIT_interval(altstrain_iit,index1));
  type2 = Interval_type(IIT_interval(altstrain_iit,index2));
  
  if (type1 < type2) {
    return -1;
  } else if (type1 > type2) {
    return +1;
  } else {
    /* Store in descending genomic position, so right shifting works
       in Genome_patch_strain */
    pos1 = Interval_low(IIT_interval(altstrain_iit,index1));
    pos2 = Interval_low(IIT_interval(altstrain_iit,index2));

    if (pos1 > pos2) {
      return -1;
    } else if (pos1 < pos2) {
      return +1;
    } else {
      return 0;
    }
  }
}
#endif


static Stage3_T *
stage3_from_usersegment (int *npaths, int *first_absmq, int *second_absmq,
			 Sequence_T queryseq, Sequence_T queryuc,
#ifdef PMAP
			 Sequence_T queryntseq,
#endif
			 Sequence_T usersegment, Stage2_alloc_T stage2_alloc,
			 Oligoindex_array_T oligoindices_major, Oligoindex_array_T oligoindices_minor,
			 Pairpool_T pairpool, Diagpool_T diagpool, Cellpool_T cellpool,
			 Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
			 Stopwatch_T worker_stopwatch) {
  List_T stage3list;
  Univcoord_T chroffset, chrhigh;
  Chrpos_T chrlength, chrpos;
  Chrnum_T chrnum = 0;

#ifdef PMAP
  Sequence_T revcomp;
#endif
		    
  chroffset = chrpos = 0U;
  chrhigh = chrlength = Sequence_fulllength(usersegment);

  stage3list = update_stage3list(/*stage3list*/NULL,queryseq,
#ifdef PMAP
				 queryntseq,
#endif
				 queryuc,stage2_alloc,oligoindices_major,oligoindices_minor,
				 pairpool,diagpool,cellpool,/*straintype*/0,/*strain*/NULL,
				 chrnum,chroffset,chrhigh,chrlength,
				 /*chrstart*/0,/*chrend*/chrhigh,/*watsonp*/true,/*genestrand*/0,
				 dynprogL,dynprogM,dynprogR,worker_stopwatch);

#ifdef PMAP
  revcomp = Sequence_revcomp(usersegment);
#endif

  stage3list = update_stage3list(stage3list,queryseq,
#ifdef PMAP
				 queryntseq,
#endif
				 queryuc,stage2_alloc,oligoindices_major,oligoindices_minor,
				 pairpool,diagpool,cellpool,/*straintype*/0,/*strain*/NULL,
				 chrnum,chroffset,chrhigh,chrlength,
				 /*chrstart*/0,/*chrend*/chrhigh,/*watsonp*/false,/*genestrand*/0,
				 dynprogL,dynprogM,dynprogR,worker_stopwatch);

#ifdef PMAP
  Sequence_free(&revcomp);
#endif

  if (stage3list == NULL) {
    *npaths = 0;
    return NULL;
  } else {
    return stage3array_from_list(&(*npaths),&(*first_absmq),&(*second_absmq),stage3list,
				 /*mergedp*/false,/*chimerap*/false,/*remove_overlaps_p*/true);
  }
}


static List_T
stage3list_remove_duplicates (List_T stage3list) {
  List_T unique = NULL;
  Stage3_T *array;
  int best_score;
  Chrpos_T shortest_genomiclength;
  int n, besti, i, j, k;
  
  if ((n = List_length(stage3list)) == 0) {
    return (List_T) NULL;
  } else if (n == 1) {
    return stage3list;
  } else {
    array = (Stage3_T *) List_to_array(stage3list,NULL);
    List_free(&stage3list);
    qsort(array,n,sizeof(Stage3_T),Stage3_position_cmp);

    i = 0;
    while (i < n) {
      best_score = Stage3_goodness(array[i]);
      shortest_genomiclength = Stage3_genomiclength(array[i]);
      besti = i;
      debug3(printf("i = %d, score %d, genomiclength %u\n",
		    i,best_score,shortest_genomiclength));

      j = i + 1;
      while (j < n && Stage3_position_cmp(&(array[i]),&(array[j])) == 0) {
	debug3(printf("  j = %d, score %d, genomiclength %u\n",
		      j,Stage3_goodness(array[j]),Stage3_genomiclength(array[j])));

	if (Stage3_goodness(array[j]) < best_score) {
	  best_score = Stage3_goodness(array[j]);
	  shortest_genomiclength = Stage3_genomiclength(array[j]);
	  besti = j;

	} else if (Stage3_goodness(array[j]) == best_score &&
		   Stage3_genomiclength(array[j]) < shortest_genomiclength) {
	  best_score = Stage3_goodness(array[j]);
	  shortest_genomiclength = Stage3_genomiclength(array[j]);
	  besti = j;
	}

	j++;
      }
      debug3(printf("  => besti = %d, score %d, genomiclength %u\n",
		    besti,best_score,shortest_genomiclength));

      for (k = i; k < j; k++) {
	if (k == besti) {
	  unique = List_push(unique,(void *) array[besti]);
	} else {
	  Stage3_free(&(array[k]));
	}
      }

      i = j;
    }
    
    FREE(array);

    return unique;
  }
}


#if 0
static List_T
stage3list_remove_empties (List_T stage3list) {
  List_T nonempty = NULL, p;
  Stage3_T stage3;
  
  for (p = stage3list; p != NULL; p = List_next(p)) {
    stage3 = (Stage3_T) List_head(p);
    if (Stage3_pairs == NULL) {
      debug2(printf("Removing empty stage3 %p\n",stage3));
      Stage3_free(&stage3);
    } else {
      nonempty = List_push(nonempty,(void *) stage3);
    }
  }

  return nonempty;
}
#endif


static List_T
stage3list_sort (List_T stage3list) {
  List_T sorted = NULL;
  Stage3_T *array;
  int n, i;

  if ((n = List_length(stage3list)) == 0) {
    return (List_T) NULL;
  } else if (n == 1) {
    return stage3list;
  } else {
    array = (Stage3_T *) List_to_array(stage3list,NULL);
    List_free(&stage3list);
    qsort(array,n,sizeof(Stage3_T),Stage3_cmp);
    for (i = n-1; i >= 0; i--) {
      sorted = List_push(sorted,(void *) array[i]);
    }
    FREE(array);

    return sorted;
  }
}


static List_T
stage3list_filter_and_sort (Chimera_T *chimera, List_T stage3list) {
  List_T sorted = NULL;
  Stage3_T *array, stage3;
  int n, i;

  if ((n = List_length(stage3list)) == 0) {
    return (List_T) NULL;

  } else if (n == 1) {
    stage3 = (Stage3_T) List_head(stage3list);
    if (Stage3_passes_filter(stage3,min_trimmed_coverage,min_identity) == false) {
      Stage3_free(&stage3);
      List_free(&stage3list);
      return (List_T) NULL;
    } else {
      return stage3list;
    }

  } else if (*chimera == NULL) {
    array = (Stage3_T *) List_to_array(stage3list,NULL);
    List_free(&stage3list);
    qsort(array,n,sizeof(Stage3_T),Stage3_cmp);
    for (i = n-1; i >= 0; i--) {
      if (Stage3_passes_filter(array[i],min_trimmed_coverage,min_identity) == false) {
	Stage3_free(&(array[i]));
      } else {
	sorted = List_push(sorted,(void *) array[i]);
      }
    }
    FREE(array);
    return sorted;

  } else if (Stage3_passes_filter_chimera(*chimera,min_trimmed_coverage,min_identity) == true) {
    array = (Stage3_T *) List_to_array(stage3list,NULL);
    List_free(&stage3list);
    qsort(array,n,sizeof(Stage3_T),Stage3_cmp);
    for (i = n-1; i >= 0; i--) {
      if (Stage3_chimera_left_p(array[i]) == true) {
	sorted = List_push(sorted,(void *) array[i]);
      } else if (Stage3_chimera_right_p(array[i]) == true) {
	sorted = List_push(sorted,(void *) array[i]);
      } else if (Stage3_passes_filter(array[i],min_trimmed_coverage,min_identity) == false) {
	Stage3_free(&(array[i]));
      } else {
	sorted = List_push(sorted,(void *) array[i]);
      }
    }
    FREE(array);
    return sorted;

  } else {
    array = (Stage3_T *) List_to_array(stage3list,NULL);
    List_free(&stage3list);
    qsort(array,n,sizeof(Stage3_T),Stage3_cmp);
    for (i = n-1; i >= 0; i--) {
      if (Stage3_chimera_left_p(array[i]) == true) {
	Stage3_free(&(array[i]));
      } else if (Stage3_chimera_right_p(array[i]) == true) {
	Stage3_free(&(array[i]));
      } else if (Stage3_passes_filter(array[i],min_trimmed_coverage,min_identity) == false) {
	Stage3_free(&(array[i]));
      } else {
	sorted = List_push(sorted,(void *) array[i]);
      }
    }
    FREE(array);

    Chimera_free(&(*chimera));
    *chimera = (Chimera_T) NULL;

    return sorted;
  }
}


static List_T
stage3_from_gregions (List_T stage3list, List_T gregions,
		      Sequence_T queryseq, Sequence_T queryuc,
#ifdef PMAP
		      Sequence_T queryntseq,
#endif
		      Sequence_T usersegment, Stage2_alloc_T stage2_alloc,
		      Oligoindex_array_T oligoindices_major, Oligoindex_array_T oligoindices_minor,
		      Pairpool_T pairpool, Diagpool_T diagpool, Cellpool_T cellpool,
		      Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		      Stopwatch_T worker_stopwatch) {
  Gregion_T gregion, *array;
  int ngregions, ncovered, max_ncovered, stage2_source;
  int i;
#if 0
  int *indexarray, nindices, straintype, j;
#endif
  void *item;

#ifdef EXTRACT_GENOMICSEG
  genomicuc_ptr = Sequence_fullpointer(genomicuc);
  Sequence_T genomicseg = NULL, genomicuc = NULL;
#endif
		    
  if (usersegment == NULL && (ngregions = List_length(gregions)) > 0) {
    array = (Gregion_T *) List_to_array(gregions,NULL);
    List_free(&gregions);

    for (i = 0; i < ngregions; i++) {
      gregion = array[i];

#if defined(EXTRACT_GENOMICSEG)
      genomicseg = Genome_get_segment(genome,Gregion_genomicstart(gregion),Gregion_genomiclength(gregion),
				      /*chromosome_iit*/NULL,Gregion_revcompp(gregion));
      genomicuc = Sequence_uppercase(genomicseg);
      genomicuc_ptr = Sequence_fullpointer(genomicuc);
#endif
      ncovered = Stage2_scan(&stage2_source,Sequence_trimpointer(queryuc),Sequence_trimlength(queryseq),
			     Gregion_chrstart(gregion),Gregion_chrend(gregion),
			     Gregion_chroffset(gregion),Gregion_chrhigh(gregion),
			     /*plusp*/Gregion_revcompp(gregion) ? false : true,Gregion_genestrand(gregion),
			     stage2_alloc,oligoindices_major,diagpool,debug_graphic_p);
      Gregion_set_ncovered(gregion,ncovered,stage2_source);
#if defined(EXTRACT_GENOMICSEG)
      Sequence_free(&genomicuc);
      Sequence_free(&genomicseg);
#endif
    }
    qsort(array,ngregions,sizeof(Gregion_T),Gregion_cmp);
    max_ncovered = Gregion_ncovered(array[0]);
    debug(printf("max_ncovered of array[0] = %d\n",max_ncovered));
    if (max_ncovered < 0.10*Sequence_fulllength(queryseq)) {
      debug(printf("coverage is too short, so skipping\n"));
      for (i = 0; i < ngregions; i++) {
	Gregion_free(&(array[i]));
      }
      FREE(array);

    } else {
      gregions = (List_T) NULL;
      i = 0;
      while (i < ngregions && Gregion_ncovered(array[i]) > 0.25*max_ncovered) {
	debug(printf("Keeping %d ncovered relative to %d\n",Gregion_ncovered(array[i]),max_ncovered));
	gregions = List_push(gregions,(void *) array[i]);
	i++;
      }
      while (i < ngregions) {
	debug(printf("Discarding array %d with ncovered = %d\n",i,Gregion_ncovered(array[i])));
	Gregion_free(&(array[i]));
	i++;
      }
      FREE(array);
    }

    while (gregions != NULL) {
      gregions = List_pop(gregions,&item);
      gregion = (Gregion_T) item;

      /* if (Match_usep(match) == true) { */
      if (1) {
	if (usersegment != NULL) {
	  /* chrlength = Sequence_fulllength(usersegment); */
	  /* strain = NULL; */
	  stage3list = update_stage3list(stage3list,queryseq,
#ifdef PMAP
					 queryntseq,
#endif
					 queryuc,stage2_alloc,oligoindices_major,oligoindices_minor,
					 pairpool,diagpool,cellpool,
					 /*straintype*/0,/*strain*/NULL,Gregion_chrnum(gregion),
					 Gregion_chroffset(gregion),Gregion_chrhigh(gregion),Gregion_chrlength(gregion),
					 Gregion_chrstart(gregion),Gregion_chrend(gregion),
					 Gregion_plusp(gregion),Gregion_genestrand(gregion),
					 dynprogL,dynprogM,dynprogR,worker_stopwatch);
	} else {
	  stage3list = update_stage3list(stage3list,queryseq,
#ifdef PMAP
					 queryntseq,
#endif
					 queryuc,stage2_alloc,oligoindices_major,oligoindices_minor,
					 pairpool,diagpool,cellpool,
					 /*straintype*/0,/*strain*/NULL,Gregion_chrnum(gregion),
					 Gregion_chroffset(gregion),Gregion_chrhigh(gregion),Gregion_chrlength(gregion),
					 Gregion_chrstart(gregion),Gregion_chrend(gregion),
					 Gregion_plusp(gregion),Gregion_genestrand(gregion),
					 dynprogL,dynprogM,dynprogR,worker_stopwatch);

#if 0
	  /* We rely upon the fact that gbuffer1 still holds the genomic segment.  This code is duplicated in get-genome.c */
	  if (altstrain_iit != NULL) {
	    indexarray = IIT_get(&nindices,altstrain_iit,/*divstring*/NULL,Gregion_genomicstart(gregion)+1U,
				 Gregion_genomicstart(gregion)+Gregion_genomiclength(gregion)-1,/*sortp*/false);
	    if (nindices > 0) {
	      /* Sort according to type and genome position */
	      qsort(indexarray,nindices,sizeof(int),index_compare);
	      j = 0;
	      while (j < nindices) {
		i = j++;
		straintype = Interval_type(IIT_interval(altstrain_iit,indexarray[i]));
		/* strain = IIT_typestring(altstrain_iit,straintype); */
		while (j < nindices && Interval_type(IIT_interval(altstrain_iit,indexarray[j])) == straintype) {
		  j++;
		}
		/* Patch from i to j */
		genomicseg = Genome_patch_strain(&(indexarray[i]),j-i,altstrain_iit,
						 Gregion_genomicstart(gregion),Gregion_genomiclength(gregion),
						 Gregion_revcompp(gregion),
						 Gbuffer_chars1(gbuffer),Gbuffer_chars2(gbuffer),Gbuffer_chars3(gbuffer),
						 Gbuffer_gbufferlen(gbuffer));
		stage3list = update_stage3list(stage3list,queryseq,
#ifdef PMAP
					       queryntseq,
#endif					     
					       queryuc,stage2_alloc,oligoindices_major,oligoindices_minor,
					       pairpool,diagpool,straintype,strain,
					       Gregion_chrnum(gregion),Gregion_chroffset(gregion),
					       Gregion_chrhigh(gregion),Gregion_chrlength(gregion),
					       Gregion_chrstart(gregion),Gregion_chrend(gregion),
					       Gregion_plusp(gregion),Gregion_genestrand(gregion),
					       dynprogL,dynprogM,dynprogR,worker_stopwatch);
		Sequence_free(&genomicseg);
	      }
	      FREE(indexarray);
	    }
	  }
#endif

	}
      }
      Gregion_free(&gregion);
    }
  }
	
  return stage3list;
}


static bool
middle_piece_local_p (int *querystart, int *queryend,
		      Chrpos_T *chrstart, Chrpos_T *chrend,
		      Chrnum_T *chrnum, Univcoord_T *chroffset, Univcoord_T *chrhigh,
		      Chrpos_T *chrlength, bool *plusp, Stage3_T from, Stage3_T to) {

  debug2(printf("? middle_piece_local_p from [%p] %d..%d (%u..%u) -> to [%p] %d..%d (%u..%u) => ",
		from,Stage3_querystart(from),Stage3_queryend(from),
		Stage3_chrstart(from),Stage3_chrend(from),
		to,Stage3_querystart(to),Stage3_queryend(to),
		Stage3_chrstart(to),Stage3_chrend(to)));

  if (Stage3_chimera_right_p(from) == true) {
    debug2(printf("false, because from is already part of a chimera on its right\n"));
    return false;
    
  } else if (Stage3_chimera_left_p(to) == true) {
    debug2(printf("false, because to is already part of a chimera on its left\n"));
    return false;

  } else if ((*chrnum = Stage3_chrnum(from)) != Stage3_chrnum(to)) {
    /* Different chromosomes */
    debug2(printf("different chromosomes\n"));
    return false;

  } else if (Stage3_watsonp(from) != Stage3_watsonp(to)) {
    /* Different strands */
    debug2(printf("different strands\n"));
    return false;

  } else if (Stage3_querystart(to) <= Stage3_queryend(from) + CHIMERA_SLOP) {
    /* Already joinable */
    debug2(printf("wrong query order or already joinable\n"));
    return false;

  } else if ((*plusp = Stage3_watsonp(from)) == true) {
    if (Stage3_chrend(from) < Stage3_chrstart(to) &&
	Stage3_chrend(from) + 1000000 > Stage3_chrstart(to)) {
      debug2(printf("true, because %u < %u and %u + %u > %u\n",
		    Stage3_chrend(from),Stage3_chrstart(to),
		    Stage3_chrend(from),1000000,Stage3_chrstart(to)));
      Univ_IIT_interval_bounds(&(*chroffset),&(*chrhigh),&(*chrlength),chromosome_iit,
			       *chrnum,circular_typeint);
      *querystart = Stage3_queryend(from);
      *queryend = Stage3_querystart(to);
      *chrstart = Stage3_chrend(from);
      *chrend = Stage3_chrstart(to);
      return true;
    } else {
      debug2(printf("false, watsonp true, from_end %u, to start %u\n",
		    Stage3_chrend(from),Stage3_chrstart(to)));
      return false;
    }

  } else {
    if (Stage3_chrstart(to) < Stage3_chrend(from) &&
	Stage3_chrstart(to) + 1000000 > Stage3_chrend(from)) {
      debug2(printf("true, because %u < %u and %u + %u > %u\n",
		    Stage3_chrstart(to),Stage3_chrend(from),
		    Stage3_chrstart(to),1000000,Stage3_chrend(from)));
      Univ_IIT_interval_bounds(&(*chroffset),&(*chrhigh),&(*chrlength),chromosome_iit,
			       *chrnum,circular_typeint);
      *querystart = Stage3_queryend(from);
      *queryend = Stage3_querystart(to);
      *chrstart = Stage3_chrstart(to);
      *chrend = Stage3_chrend(from);
      return true;
    } else {
      debug2(printf("false, watsonp false, from_end %u, to start %u\n",
		    Stage3_chrend(from),Stage3_chrstart(to)));
      return false;
    }
  }
}


static bool
middle_piece_chimera_p (int *querystart, int *queryend, Stage3_T from, Stage3_T to) {

  debug2(printf("? middle_piece_chimera_p from [%p] %d..%d (%u..%u) -> to [%p] %d..%d (%u..%u) => ",
		from,Stage3_querystart(from),Stage3_queryend(from),
		Stage3_chrstart(from),Stage3_chrend(from),
		to,Stage3_querystart(to),Stage3_queryend(to),
		Stage3_chrstart(to),Stage3_chrend(to)));

  if (Stage3_chimera_right_p(from) == true) {
    debug2(printf("false, because from is already part of a chimera on its right\n"));
    return false;
    
  } else if (Stage3_chimera_left_p(to) == true) {
    debug2(printf("false, because to is already part of a chimera on its left\n"));
    return false;

  } else if (Stage3_querystart(to) <= Stage3_queryend(from) + CHIMERA_SLOP) {
    /* Already joinable */
    debug2(printf("wrong query order or already joinable\n"));
    return false;

  } else {
    *querystart = Stage3_queryend(from);
    *queryend = Stage3_querystart(to);
    return true;
  }
}


/* Returns nonjoinable */
static List_T
local_separate_paths (Stage3_T **stage3array_sub1, int *npaths_sub1, 
		      Stage3_T **stage3array_sub2, int *npaths_sub2,
		      List_T stage3list) {
  List_T nonjoinable = NULL, p;
  Stage3_T from, to, stage3;
  Stage3_T *by_queryend, *by_querystart;
  int npaths, i, j, k;
  int queryend;

  debug2(printf("local_separate_paths called with list length %d\n",List_length(stage3list)));

  if (stage3list == NULL) {
    *stage3array_sub1 = (Stage3_T *) NULL;
    *npaths_sub1 = 0;
    *stage3array_sub2 = (Stage3_T *) NULL;
    *npaths_sub2 = 0;
    return (List_T) NULL;
  } else {
    for (p = stage3list; p != NULL; p = List_next(p)) {
      stage3 = (Stage3_T) List_head(p);
      Stage3_clear_joinable(stage3);
    }
  }

  by_queryend = (Stage3_T *) List_to_array_n(&npaths,stage3list);
  qsort(by_queryend,npaths,sizeof(Stage3_T),Stage3_queryend_cmp);

  by_querystart = (Stage3_T *) List_to_array_n(&npaths,stage3list);
  qsort(by_querystart,npaths,sizeof(Stage3_T),Stage3_querystart_cmp);

#ifdef DEBUG2
  for (i = 0; i < npaths; i++) {
    stage3 = (Stage3_T) by_queryend[i];
    printf("from: %p query %d..%d, genomic %u..%u\t",
	   stage3,Stage3_querystart(stage3),Stage3_queryend(stage3),
	   Stage3_genomicstart(stage3),Stage3_genomicend(stage3));

    stage3 = (Stage3_T) by_querystart[i];
    printf("to: %p query %d..%d, genomic %u..%u\n",
	   stage3,Stage3_querystart(stage3),Stage3_queryend(stage3),
	   Stage3_genomicstart(stage3),Stage3_genomicend(stage3));
  }
#endif

  j = 0;
  for (i = 0; i < npaths; i++) {
    from = by_queryend[i];
    queryend = Stage3_queryend(from);

    while (j < npaths && Stage3_querystart(by_querystart[j]) < queryend + CHIMERA_SLOP) {
      j++;
    }
    j--;

    while (j >= 0 && Stage3_querystart(by_querystart[j]) > queryend - CHIMERA_SLOP) {
      j--;
    }
    j++;

    while (j < npaths && Stage3_querystart(by_querystart[j]) < queryend + CHIMERA_SLOP) {
      to = by_querystart[j];

      if (Chimera_local_join_p(from,to,CHIMERA_SLOP) == true) {
	debug2(printf("Found local join from %d to %d\n",i,j));
	Stage3_set_joinable_left(from);
	Stage3_set_joinable_right(to);
      }

      j++;
    }
  }

  FREE(by_querystart);
  FREE(by_queryend);


  *npaths_sub1 = *npaths_sub2 = 0;
  for (p = stage3list; p != NULL; p = List_next(p)) {
    stage3 = (Stage3_T) List_head(p);
    if (Stage3_joinable_left_p(stage3) == true) {
      debug2(printf("Putting stage3 %p into local sub1\n",stage3));
      (*npaths_sub1)++;
    }
    if (Stage3_joinable_right_p(stage3) == true) {
      debug2(printf("Putting stage3 %p into local sub2\n",stage3));
      (*npaths_sub2)++;
    }
  }

  if (*npaths_sub1 == 0 || *npaths_sub2 == 0) {
    *stage3array_sub1 = (Stage3_T *) NULL;
    *npaths_sub1 = 0;
    *stage3array_sub2 = (Stage3_T *) NULL;
    *npaths_sub2 = 0;
  } else {
    *stage3array_sub1 = (Stage3_T *) MALLOC((*npaths_sub1) * sizeof(Stage3_T)); /* Return value */
    *stage3array_sub2 = (Stage3_T *) MALLOC((*npaths_sub2) * sizeof(Stage3_T)); /* Return value */
    j = k = 0;
    for (p = stage3list; p != NULL; p = List_next(p)) {
      stage3 = (Stage3_T) List_head(p);
      if (Stage3_joinable_left_p(stage3) == false && Stage3_joinable_right_p(stage3) == false) {
	nonjoinable = List_push(nonjoinable,stage3);
      } else {
	/* Note: it is possible that the same stage3 object gets put into both lists */
	if (Stage3_joinable_left_p(stage3) == true) {
	  (*stage3array_sub1)[j++] = stage3;
	}
	if (Stage3_joinable_right_p(stage3) == true) {
	  (*stage3array_sub2)[k++] = stage3;
	}
      }
    }
  }

  return nonjoinable;
}


/* Returns nonjoinable */
static List_T
distant_separate_paths (Stage3_T **stage3array_sub1, int *npaths_sub1, 
			Stage3_T **stage3array_sub2, int *npaths_sub2,
			List_T stage3list) {
  List_T nonjoinable = NULL, p;
  Stage3_T from, to, stage3;
  Stage3_T *by_queryend, *by_querystart;
  int npaths, i, j, k;
  int queryend;

  debug2(printf("distant_separate_paths called with list length %d\n",List_length(stage3list)));

  if (stage3list == NULL) {
    *stage3array_sub1 = (Stage3_T *) NULL;
    *npaths_sub1 = 0;
    *stage3array_sub2 = (Stage3_T *) NULL;
    *npaths_sub2 = 0;
    return (List_T) NULL;
  } else {
    for (p = stage3list; p != NULL; p = List_next(p)) {
      stage3 = (Stage3_T) List_head(p);
      Stage3_clear_joinable(stage3);
    }
  }

  by_queryend = (Stage3_T *) List_to_array_n(&npaths,stage3list);
  qsort(by_queryend,npaths,sizeof(Stage3_T),Stage3_queryend_cmp);

  by_querystart = (Stage3_T *) List_to_array_n(&npaths,stage3list);
  qsort(by_querystart,npaths,sizeof(Stage3_T),Stage3_querystart_cmp);

  j = 0;
  for (i = 0; i < npaths; i++) {
    from = by_queryend[i];
    queryend = Stage3_queryend(from);

    while (j < npaths && Stage3_querystart(by_querystart[j]) < queryend + CHIMERA_SLOP) {
      j++;
    }
    j--;

    while (j >= 0 && Stage3_querystart(by_querystart[j]) > queryend - CHIMERA_SLOP) {
      j--;
    }
    j++;

    while (j < npaths && Stage3_querystart(by_querystart[j]) < queryend + CHIMERA_SLOP) {
      to = by_querystart[j];

      if (Chimera_distant_join_p(from,to,CHIMERA_SLOP) == true) {
	debug2(printf("Found distant join from %d to %d\n",i,j));
	Stage3_set_joinable_left(from);
	Stage3_set_joinable_right(to);
      }

      j++;
    }
  }

  FREE(by_querystart);
  FREE(by_queryend);


  *npaths_sub1 = *npaths_sub2 = 0;
  for (p = stage3list; p != NULL; p = List_next(p)) {
    stage3 = (Stage3_T) List_head(p);
    if (Stage3_joinable_left_p(stage3) == true) {
      (*npaths_sub1)++;
    }
    if (Stage3_joinable_right_p(stage3) == true) {
      (*npaths_sub2)++;
    }
  }

  if (*npaths_sub1 == 0 || *npaths_sub2 == 0) {
    *stage3array_sub1 = (Stage3_T *) NULL;
    *npaths_sub1 = 0;
    *stage3array_sub2 = (Stage3_T *) NULL;
    *npaths_sub2 = 0;
  } else {
    *stage3array_sub1 = (Stage3_T *) MALLOC((*npaths_sub1) * sizeof(Stage3_T)); /* Return value */
    *stage3array_sub2 = (Stage3_T *) MALLOC((*npaths_sub2) * sizeof(Stage3_T)); /* Return value */
    j = k = 0;
    for (p = stage3list; p != NULL; p = List_next(p)) {
      stage3 = (Stage3_T) List_head(p);
      if (Stage3_joinable_left_p(stage3) == false && Stage3_joinable_right_p(stage3) == false) {
	nonjoinable = List_push(nonjoinable,stage3);
      } else {
	/* Note: it is possible that the same stage3 object gets put into both lists */
	if (Stage3_joinable_left_p(stage3) == true) {
	  debug2(printf("Putting stage3 %p into distant sub1\n",stage3));
	  (*stage3array_sub1)[j++] = stage3;
	}
	if (Stage3_joinable_right_p(stage3) == true) {
	  debug2(printf("Putting stage3 %p into distant sub2\n",stage3));
	  (*stage3array_sub2)[k++] = stage3;
	}
      }
    }
  }

  return nonjoinable;
}


/* Returns a list with only one Stage3_T object */
static List_T
merge_left_and_right_readthrough (bool *mergedp, Stage3_T *stage3array_sub1, int npaths_sub1, int bestfrom,
				  Stage3_T *stage3array_sub2, int npaths_sub2, int bestto,
				  List_T nonjoinable, int breakpoint, int queryntlength,
#ifdef PMAP
				  char *queryaaseq_ptr,
#endif
				  char *queryseq_ptr, char *queryuc_ptr,
				  Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
				  Oligoindex_array_T oligoindices_minor, Diagpool_T diagpool, Cellpool_T cellpool) {
  List_T newstage3list, p;
  Stage3_T best0, best1, *array, last, freed0 = NULL, freed1 = NULL;
  int i, k;

  best0 = stage3array_sub1[bestfrom];
  best1 = stage3array_sub2[bestto];

  debug2(printf("\nEntering merge_left_and_right_readthrough with bestfrom %d: %p, bestto %d: %p, and nonjoinable %d\n",
		bestfrom,best0,bestto,best1,List_length(nonjoinable)));

#if 0
  /* Checked better by Stage3_mergeable */
  if (Stage3_sensedir(best0) != Stage3_sensedir(best1) &&
      Stage3_sensedir(best0) != SENSE_NULL && Stage3_sensedir(best1) != SENSE_NULL) {
    debug2(printf("sensedirs are not compatible: %d and %d\n",
		  Stage3_sensedir(best0),Stage3_sensedir(best1)));
    if (Stage3_npairs(best0) > Stage3_npairs(best1)) {
      newstage3list = (List_T) NULL;
      newstage3list = List_push(newstage3list,(void *) best0);
      freed1 = best1;
      Stage3_free(&best1);
    } else {
      newstage3list = (List_T) NULL;
      newstage3list = List_push(newstage3list,(void *) best1);
      freed0 = best0;
      Stage3_free(&best0);
    }
    *mergedp = false;

  }
#endif

  debug2(printf("Running Stage3_merge_local\n"));
  if (Stage3_merge_local(best0,best1,/*minpos1*/0,/*maxpos1*/breakpoint,
			 /*minpos2*/breakpoint+1,/*maxpos2*/queryntlength,
			 /*genestrand*/0,
#ifdef PMAP
			 queryaaseq_ptr,
#endif
			 queryseq_ptr,queryuc_ptr,
			 pairpool,dynprogL,dynprogM,dynprogR,
			 maxpeelback,oligoindices_minor,diagpool,cellpool) == false) {

    newstage3list = (List_T) NULL;
    newstage3list = List_push(newstage3list,(void *) best0);
    newstage3list = List_push(newstage3list,(void *) best1);
    for (p = nonjoinable; p != NULL; p = List_next(p)) {
      debug2(printf("1.  Pushing readthrough nonjoinable stage3 %p.  %d..%d\n",
		    List_head(p),Stage3_querystart(List_head(p)),Stage3_queryend(List_head(p))));
      if (List_head(p) == NULL) {
	debug2(printf("Unexpected: Have a NULL stage3 in nonjoinable\n"));
      } else {
	newstage3list = List_push(newstage3list,(void *) List_head(p));
      }
    }
    *mergedp = false;
    return List_reverse(newstage3list);

  } else {
    debug2(printf("done with Stage3_merge_local"));

    debug2(printf("Rearranging paths\n"));
    debug2(printf("Changing genomicend of merged stage3 from %u to %u\n",Stage3_genomicend(best0),Stage3_genomicend(best1)));
    Stage3_set_genomicend(best0,Stage3_genomicend(best1));
    newstage3list = (List_T) NULL;
    newstage3list = List_push(newstage3list,(void *) best0);
    debug2(printf("Freeing best1 %p\n",best1));
    freed1 = best1;
    Stage3_free(&best1);
    debug2(printf("Pushing stage3 %p: ",best0));
    debug2(Stage3_print_ends(best0));
    *mergedp = true;

    if (npaths_sub1 + npaths_sub2 > 2) {
      /* Push rest of results, taking care not to have duplicates */

      array = (Stage3_T *) MALLOCA((npaths_sub1 + npaths_sub2 - 2) * sizeof(Stage3_T));
      k = 0;
      for (i = 0; i < npaths_sub1; i++) {
	if (i != bestfrom) {
	  debug2(printf("array %d is now sub1 %d: %p\n",k,i,stage3array_sub1[i]));
	  array[k++] = stage3array_sub1[i];
	}
      }
      for (i = 0; i < npaths_sub2; i++) {
	if (i != bestto) {
	  debug2(printf("array %d is now sub2 %d: %p\n",k,i,stage3array_sub2[i]));
	  array[k++] = stage3array_sub2[i];
	}
      }
      qsort(array,npaths_sub1+npaths_sub2-2,sizeof(Stage3_T),Stage3_identity_cmp);

      last = (Stage3_T) NULL;
      for (i = 0; i < npaths_sub1+npaths_sub2-2; i++) {
	if (array[i] == last) {
	  /* Skip */
	  debug2(printf("array %d: Skipping stage3 %p, because just pushed, so duplicate\n",i,array[i]));
	} else if (array[i] == best0 || array[i] == best1) {
	  /* Skip */
	  debug2(printf("array %d: Skipping stage3 %p, because in chimera\n",i,array[i]));
	} else if (array[i] == freed0 || array[i] == freed1) {
	  /* Skip */
	  debug2(printf("array %d: Skipping stage3 %p, because already freed\n",i,array[i]));
	} else {
	  debug2(printf("array %d: Pushing stage3 %p\n",i,array[i]));
	  newstage3list = List_push(newstage3list,(void *) array[i]);
	  last = array[i];
	}
      }

      FREEA(array);
    }

    for (p = nonjoinable; p != NULL; p = List_next(p)) {
      debug2(printf("2.  Pushing readthrough nonjoinable stage3 %p.  %d..%d\n",
		    List_head(p),Stage3_querystart(List_head(p)),Stage3_queryend(List_head(p))));
      newstage3list = List_push(newstage3list,(void *) List_head(p));
    }

    return List_reverse(newstage3list);
  }
}


/* Returns a list with only two Stage3_T objects */
static List_T
merge_left_and_right_transloc (Stage3_T *stage3array_sub1, int npaths_sub1, int bestfrom,
			       Stage3_T *stage3array_sub2, int npaths_sub2, int bestto,
			       List_T nonjoinable) {
  List_T newstage3list, p;
  Stage3_T best0, best1, *array, last;
  int i, k;

  best0 = stage3array_sub1[bestfrom];
  best1 = stage3array_sub2[bestto];

  debug2(printf("\nEntering merge_left_and_right_transloc with bestfrom %d: %p, bestto %d: %p, and nonjoinable %d\n",
		bestfrom,best0,bestto,best1,List_length(nonjoinable)));

  debug2(printf("Before Stage3_merge_chimera, best0 is %p, query %d..%d\n",
		best0,Stage3_querystart(best0),Stage3_queryend(best0)));
  debug2(Stage3_print_ends(best0));
  debug2(printf("Before Stage3_merge_chimera, best1 is %p, query %d..%d\n",
		best1,Stage3_querystart(best1),Stage3_queryend(best1)));
  debug2(Stage3_print_ends(best1));

  debug2(printf("Rearranging paths\n"));
  newstage3list = (List_T) NULL;
  newstage3list = List_push(newstage3list,(void *) best0);
  newstage3list = List_push(newstage3list,(void *) best1);
  debug2(printf("Pushing stage3 %p, ",best0));
  debug2(Stage3_print_ends(best0));
  debug2(printf("Pushing stage3 %p, ",best1));
  debug2(Stage3_print_ends(best1));

  if (npaths_sub1 + npaths_sub2 > 2) {
    /* Push rest of results, taking care not to have duplicates */

    array = (Stage3_T *) MALLOCA((npaths_sub1 + npaths_sub2 - 2) * sizeof(Stage3_T));
    k = 0;
    for (i = 0; i < npaths_sub1; i++) {
      if (i != bestfrom) {
	array[k++] = stage3array_sub1[i];
      }
    }
    for (i = 0; i < npaths_sub2; i++) {
      if (i != bestto) {
	array[k++] = stage3array_sub2[i];
      }
    }
    qsort(array,npaths_sub1+npaths_sub2-2,sizeof(Stage3_T),Stage3_identity_cmp);

    last = (Stage3_T) NULL;
    for (i = 0; i < npaths_sub1+npaths_sub2-2; i++) {
      if (array[i] == last) {
	/* Skip */
	debug2(printf("Skipping stage3 %p, because just pushed\n",array[i]));
      } else if (array[i] == best0 || array[i] == best1) {
	/* Skip */
	debug2(printf("Skipping stage3 %p, because in chimera\n",array[i]));
      } else {
	debug2(printf("Pushing stage3 %p.  ",array[i]));
	debug2(Stage3_print_ends(array[i]));
	newstage3list = List_push(newstage3list,(void *) array[i]);
	last = array[i];
      }
    }

    FREEA(array);
  }

  for (p = nonjoinable; p != NULL; p = List_next(p)) {
    debug2(printf("Pushing transloc nonjoinable stage3 %p\n",List_head(p)));
    newstage3list = List_push(newstage3list,(void *) List_head(p));
  }

  return List_reverse(newstage3list);
}


static int
find_breakpoint (int *cdna_direction, int *chimerapos, int *chimeraequivpos, int *exonexonpos,
		 char *donor1, char *donor2, char *acceptor2, char *acceptor1,
		 bool *donor_watsonp, bool *acceptor_watsonp, double *donor_prob, double *acceptor_prob,
		 Stage3_T from, Stage3_T to,
#ifdef PMAP
		 Sequence_T queryntseq,
#endif
		 Sequence_T queryseq, Sequence_T queryuc,
		 int queryntlength, Genome_T genome, Genome_T genomealt,
		 Univ_IIT_T chromosome_iit, Pairpool_T pairpool) {
  int breakpoint, leftpos, rightpos, midpos;
  int maxpeelback_from, maxpeelback_to;
  int found_cdna_direction, try_cdna_direction;
  char comp;			/* Not really used anywhere */

  int queryjump;
  int genomejump;
  bool max_extend_p;
  Chrpos_T left_chrlength, right_chrlength;
  Univcoord_T chroffset, chrhigh;

  if (Stage3_queryend(from) < Stage3_querystart(to)) {
    /* Gap exists between the two parts */
    if ((leftpos = Stage3_queryend(from) - 8) < 0) {
      leftpos = 0;
    }
    if ((rightpos = Stage3_querystart(to) + 8) >= queryntlength) {
      rightpos = queryntlength - 1;
    }
    maxpeelback_from = 8;
    maxpeelback_to = 8;
    debug2(printf("overlap: leftpos %d, rightpos %d, queryntlength %d, maxpeelback_from %d, maxpeelback_to %d\n",
		  leftpos,rightpos,queryntlength,maxpeelback_from,maxpeelback_to));

    if (Stage3_watsonp(from) == true && Stage3_watsonp(to) == true) {
      queryjump = Stage3_querystart(to) - Stage3_queryend(from) - 1;
      genomejump = Stage3_genomicstart(to) - Stage3_genomicend(from) - 1U;
      max_extend_p = ((int) genomejump == queryjump) ? false : true;
      debug2(printf("gap exists: genomejump = %u, queryjump = %d, max_extend_p = %d\n",genomejump,queryjump,max_extend_p));
    } else if (Stage3_watsonp(from) == false && Stage3_watsonp(to) == false) {
      queryjump = Stage3_querystart(to) - Stage3_queryend(from) - 1;
      genomejump = Stage3_genomicend(from) - Stage3_genomicstart(to) - 1U;
      max_extend_p = ((int) genomejump == queryjump) ? false : true;
      debug2(printf("gap exists: genomejump = %u, queryjump = %d, max_extend_p = %d\n",genomejump,queryjump,max_extend_p));
    } else {
      max_extend_p = false;
    }
    
  } else {
    /* Two parts overlap */
    if ((leftpos = Stage3_querystart(to) - 8) < 0) {
      leftpos = 0;
    }
    if ((rightpos = Stage3_queryend(from) + 8) >= queryntlength) {
      rightpos = queryntlength - 1;
    }
    midpos = (leftpos+rightpos)/2;
    /* maxpeelback_from = rightpos - Stage3_querystart(to); */
    /* maxpeelback_to = Stage3_queryend(from) - leftpos; */
    maxpeelback_from = rightpos - midpos;
    maxpeelback_to = midpos - leftpos;
    debug2(printf("overlap: leftpos %d, rightpos %d, midpos %d, queryntlength %d, maxpeelback_from %d, maxpeelback_to %d\n",
		  leftpos,rightpos,midpos,queryntlength,maxpeelback_from,maxpeelback_to));
#if 0
    if (Stage3_watsonp(from) == true && Stage3_watsonp(to) == true) {
      queryjump = Stage3_queryend(from) - Stage3_querystart(to) - 1;
      genomejump = Stage3_genomicend(from) - Stage3_genomicstart(to) - 1U;
      max_extend_p = (genomejump == queryjump) ? false : true;
    } else if (Stage3_watsonp(from) == false && Stage3_watsonp(to) == false) {
      queryjump = Stage3_queryend(from) - Stage3_querystart(to) - 1;
      genomejump = Stage3_genomicstart(to) - Stage3_genomicend(from) - 1U;
      max_extend_p = (genomejump == queryjump) ? false : true;
    } else {
      max_extend_p = false;
    }
#else
    debug2(printf("parts overlap: max_extend_p is false\n"));
    max_extend_p = false;
#endif
  }

  debug2(printf("Before Stage3_extend_right, bestfrom is %p, query %d..%d\n",
		from,Stage3_querystart(from),Stage3_queryend(from)));
  debug2(Stage3_print_ends(from));
  debug2(printf("Before Stage3_extend_left, bestto is %p, query %d..%d\n",
		to,Stage3_querystart(to),Stage3_queryend(to)));
  debug2(Stage3_print_ends(to));
  
  Stage3_extend_right(from,/*goal*/rightpos,
#ifdef PMAP
		      /*querylength*/Sequence_fulllength(queryntseq),
		      /*queryseq_ptr*/Sequence_fullpointer(queryntseq),
		      /*queryuc_ptr*/Sequence_fullpointer(queryntseq),
#else
		      /*querylength*/Sequence_fulllength(queryseq),
		      /*queryseq_ptr*/Sequence_fullpointer(queryseq),
		      /*queryuc_ptr*/Sequence_fullpointer(queryuc),
#endif
		      max_extend_p,pairpool,maxpeelback_from);

  Stage3_extend_left(to,/*goal*/leftpos,
#ifdef PMAP
		     /*queryseq_ptr*/Sequence_fullpointer(queryntseq),
		     /*queryuc_ptr*/Sequence_fullpointer(queryntseq),
#else
		     /*queryseq_ptr*/Sequence_fullpointer(queryseq),
		     /*queryuc_ptr*/Sequence_fullpointer(queryuc),
#endif
		     max_extend_p,pairpool,maxpeelback_to);

  debug2(printf("Before Chimera_find_breakpoint, bestfrom is %p, query %d..%d\n",
		from,Stage3_querystart(from),Stage3_queryend(from)));
  debug2(Stage3_print_ends(from));
  debug2(printf("Before Chimera_find_breakpoint, bestto is %p, query %d..%d\n",
		to,Stage3_querystart(to),Stage3_queryend(to)));
  debug2(Stage3_print_ends(to));

  debug2(printf("Before Chimera_find_exonexon, bestfrom is %p, query %d..%d\n",
		from,Stage3_querystart(from),Stage3_queryend(from)));
  debug2(printf("Before Chimera_find_exonexon, bestto is %p, query %d..%d\n",
		to,Stage3_querystart(to),Stage3_queryend(to)));

  if ((*exonexonpos = Chimera_find_exonexon(&found_cdna_direction,&try_cdna_direction,
					    &(*donor1),&(*donor2),&(*acceptor2),&(*acceptor1),
					    &comp,&(*donor_watsonp),&(*acceptor_watsonp),&(*donor_prob),&(*acceptor_prob),
					    /*left_part*/from,/*right_part*/to,genome,genomealt ? genomealt : genome,
					    chromosome_iit,/*breakpoint_start*/Stage3_querystart(to),
					    /*breakpoint_end*/Stage3_queryend(from))) > 0) {
    breakpoint = *chimerapos = *chimeraequivpos = *exonexonpos;
    *cdna_direction = found_cdna_direction;
    debug2(printf("Exon-exon boundary found at %d, which is breakpoint.  Comp = %c\n",
		  *exonexonpos,comp));
    return breakpoint;

  } else {
    Univ_IIT_interval_bounds(&chroffset,&chrhigh,&left_chrlength,chromosome_iit,Stage3_chrnum(from),circular_typeint);
    Univ_IIT_interval_bounds(&chroffset,&chrhigh,&right_chrlength,chromosome_iit,Stage3_chrnum(to),circular_typeint);

    if ((*chimerapos = Chimera_find_breakpoint(&(*chimeraequivpos),&(*donor1),&(*donor2),&(*acceptor2),&(*acceptor1),
					       from,to,queryntlength,genome,left_chrlength,right_chrlength)) < 0) {
      /* TODO: Allow finding a breakpoint for DNA-Seq, which needs no donor or acceptor nucleotides */
      debug2(printf("Chimera_find_breakpoint returns no value\n"));
      *donor_prob = *acceptor_prob = 0.0;
      *donor_watsonp = *acceptor_watsonp = true;
      *cdna_direction = 0;
      return -1;

    } else {
      *donor_prob = *acceptor_prob = 0.0;
      *donor_watsonp = *acceptor_watsonp = true;
      
      debug2(printf("Chimera_find_breakpoint returns boundary at %d..%d (switch can occur at %d..%d)\n",
		    *chimerapos,*chimeraequivpos,(*chimerapos)-1,*chimeraequivpos));
      
      breakpoint = ((*chimerapos) + (*chimeraequivpos))/2;
      *cdna_direction = try_cdna_direction;
      debug2(printf("Exon-exon boundary not found, but setting breakpoint to be %d\n",breakpoint));
      return breakpoint;
    }
  }
}


static List_T
check_for_local (bool *mergedp, List_T stage3list, int effective_start, int effective_end,
		 Sequence_T queryseq, Sequence_T queryuc,
#ifdef PMAP
		 Sequence_T queryntseq,
#endif
		 int queryntlength, Sequence_T usersegment, Stage2_alloc_T stage2_alloc,
		 Oligoindex_array_T oligoindices_major, Oligoindex_array_T oligoindices_minor,
		 Matchpool_T matchpool, Pairpool_T pairpool, Diagpool_T diagpool, Cellpool_T cellpool,
		 Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR) {
  List_T gregions = NULL, nonjoinable = NULL, p;
  Stage3_T *stage3array_sub1 = NULL, *stage3array_sub2 = NULL, from, to, stage3;
  Sequence_T querysubseq = NULL, querysubuc = NULL;
  Diagnostic_T diagnostic;
  int bestfrom, bestto;
  int five_margin, three_margin, five_score = 0, three_score = 0;
  int extension;
  int npaths_sub1 = 0, npaths_sub2 = 0;
  bool lowidentityp, poorp, repetitivep;

  int max_single_goodness;
  int breakpoint, chimerapos, chimeraequivpos, exonexonpos;
  int chimera_cdna_direction;
  char donor1, donor2, acceptor2, acceptor1;
  bool donor_watsonp, acceptor_watsonp;
  double donor_prob, acceptor_prob;
  

#ifdef PMAP
  five_margin = effective_start - 3*Sequence_trim_start(queryseq);
  three_margin = 3*Sequence_trim_end(queryseq) - effective_end;
  debug2(printf("Margins are %d = %d - %d on the 5' end and %d = %d - %d on the 3' end\n",
		five_margin,effective_start,3*Sequence_trim_start(queryseq),
		three_margin,3*Sequence_trim_end(queryseq),effective_end));
#else
  five_margin = effective_start - Sequence_trim_start(queryseq);
  three_margin = Sequence_trim_end(queryseq) - effective_end;
  debug2(printf("Margins are %d = %d - %d on the 5' end and %d = %d - %d on the 3' end\n",
		five_margin,effective_start,Sequence_trim_start(queryseq),
		three_margin,Sequence_trim_end(queryseq),effective_end));
#endif

#ifdef DEBUG2A
  for (p = stage3list; p != NULL; p = List_next(p)) {
    stage3 = (Stage3_T) List_head(p);
    Pair_dump_array(Stage3_pairarray(stage3),Stage3_npairs(stage3),/*zerobasedp*/true);
    printf("\n");
  }
#endif

  /* Stage3_recompute_goodness(stage3list); */
  max_single_goodness = 0;
  for (p = stage3list; p != NULL; p = List_next(p)) {
    stage3 = (Stage3_T) List_head(p);
    if (Stage3_goodness(stage3) > max_single_goodness) {
      max_single_goodness = Stage3_goodness(stage3);
    }
  }
  debug2(printf("max single goodness = %d\n",max_single_goodness));


  /* List_free(&nonjoinable); */
  debug2(printf("Running local_separate_paths\n"));
  nonjoinable = local_separate_paths(&stage3array_sub1,&npaths_sub1,&stage3array_sub2,&npaths_sub2,
				     stage3list);
  debug2(printf("local: npaths_sub1 %d, npaths_sub2 %d, nonjoinable %d\n",
		npaths_sub1,npaths_sub2,List_length(nonjoinable)));

  if (npaths_sub1 == 0 && npaths_sub2 == 0) {
    /* Need to compute on margin explicitly */
    if (five_margin < chimera_margin && three_margin < chimera_margin) {
      debug2(printf("Insufficient margins\n"));
    } else if (five_margin > three_margin) {
#if 0
      /* extension makes it harder to find the other alignment.  The merging process will help fill in any gap. */
      extension = CHIMERA_SLOP;
      debug2(printf("Comparing extension %d with %d = (effective_start %d)/2\n",
		    extension,effective_start/2,effective_start));
      if (extension > effective_start/2) {
	/* Extension occupies more than 1/3 of sequence */
	debug2(printf("Proposed extension of %d is too long relative to effective_start %d\n",extension,effective_start));
	extension = effective_start/3;
      }
#else
      extension = 0;
#endif
      if ((querysubseq = Sequence_subsequence(queryseq,0,effective_start+extension)) != NULL) {
	if ((querysubuc = Sequence_subsequence(queryuc,0,effective_start+extension)) != NULL) {
	  debug2(printf("5 margin > 3 margin.  "));
	  debug2(printf("Beginning Stage1_compute on 5' margin from effective_start %d (%d..%d)\n",
			effective_start,0,effective_start+extension));
	  debug2a(Sequence_print(stdout,querysubseq,/*uppercasep*/true,wraplength,/*trimmedp*/true));

	  diagnostic = evaluate_query(&poorp,&repetitivep,Sequence_fullpointer(querysubuc),Sequence_fulllength(querysubuc),
				      Oligoindex_array_elt(oligoindices_major,0));
	  if (poorp == true || repetitivep == true) {
	    debug2(printf("Subsequence is poor or repetitive\n"));
	  } else {
	    gregions = Stage1_compute(&lowidentityp,querysubuc,indexdb_fwd,indexdb_rev,
				      /*indexdb_size_threshold*/100,chromosome_iit,
				      chrsubset_start,chrsubset_end,matchpool,
				      stutterhits,diagnostic,/*worker_stopwatch*/NULL,/*nbest*/10);
	    debug2(printf("A.  Performing Stage 3 starting with list length %d\n",List_length(stage3list)));
	    stage3list = stage3_from_gregions(stage3list,gregions,querysubseq,querysubuc,
#ifdef PMAP
					      queryntseq,
#endif
					      usersegment,stage2_alloc,oligoindices_major,oligoindices_minor,
					      pairpool,diagpool,cellpool,
					      dynprogL,dynprogM,dynprogR,/*worker_stopwatch*/NULL);
#ifdef DEBUG2
	    for (p = stage3list; p != NULL; p = List_next(p)) {
	      stage3 = (Stage3_T) List_head(p);
	      printf("%d..%d, %u..%u\n",
		     Stage3_querystart(stage3),Stage3_queryend(stage3),
		     Stage3_genomicstart(stage3),Stage3_genomicend(stage3));
	    }
#endif
	  }
	  Diagnostic_free(&diagnostic);

	  /* Above function frees gregions */
	  Sequence_free(&querysubuc);
	}
	Sequence_free(&querysubseq);
      }

      /* And recompute on original part, just in case stage 1 was led astray by the ends */
      if ((querysubseq = Sequence_subsequence(queryseq,effective_start,queryntlength)) != NULL) {
	if ((querysubuc = Sequence_subsequence(queryuc,effective_start,queryntlength)) != NULL) {
	  debug2(printf("Recomputing on original part.  "));
	  debug2(printf("Beginning Stage1_compute on 5' margin from effective_start %d (%d..%d)\n",
			effective_start,effective_start,queryntlength));
	  debug2a(Sequence_print(stdout,querysubseq,/*uppercasep*/true,wraplength,/*trimmedp*/true));

	  diagnostic = evaluate_query(&poorp,&repetitivep,Sequence_fullpointer(querysubuc),Sequence_fulllength(querysubuc),
				      Oligoindex_array_elt(oligoindices_major,0));
	  if (poorp == true || repetitivep == true) {
	    debug2(printf("Subsequence is poor or repetitive\n"));
	  } else {
	    gregions = Stage1_compute(&lowidentityp,querysubuc,indexdb_fwd,indexdb_rev,
				      /*indexdb_size_threshold*/100,chromosome_iit,
				      chrsubset_start,chrsubset_end,matchpool,
				      stutterhits,diagnostic,/*worker_stopwatch*/NULL,/*nbest*/10);
	    debug2(printf("B.  Performing Stage 3 starting with list length %d\n",List_length(stage3list)));
	    stage3list = stage3_from_gregions(stage3list,gregions,querysubseq,querysubuc,
#ifdef PMAP
					      queryntseq,
#endif
					      usersegment,stage2_alloc,oligoindices_major,oligoindices_minor,
					      pairpool,diagpool,cellpool,
					      dynprogL,dynprogM,dynprogR,/*worker_stopwatch*/NULL);
#ifdef DEBUG2
	    for (p = stage3list; p != NULL; p = List_next(p)) {
	      stage3 = (Stage3_T) List_head(p);
	      printf("%d..%d, %u..%u\n",
		     Stage3_querystart(stage3),Stage3_queryend(stage3),
		     Stage3_genomicstart(stage3),Stage3_genomicend(stage3));
	    }
#endif
	  }
	  Diagnostic_free(&diagnostic);

	  /* Above function frees gregions */
	  Sequence_free(&querysubuc);
	}
	Sequence_free(&querysubseq);
      }

      List_free(&nonjoinable);
      debug2(printf("Running local_separate_paths\n"));
      nonjoinable = local_separate_paths(&stage3array_sub1,&npaths_sub1,&stage3array_sub2,&npaths_sub2,
					 stage3list);
      debug2(printf("local: npaths_sub1 %d, npaths_sub2 %d, nonjoinable %d\n",
		    npaths_sub1,npaths_sub2,List_length(nonjoinable)));

    } else {
#if 0
      /* extension makes it harder to find the other alignment.  The merging process will help fill in any gap. */
      extension = CHIMERA_SLOP;
      debug2(printf("Comparing extension %d with %d = (queryntlength %d - effective_end %d)/2\n",
		    extension,(queryntlength-effective_end)/2,queryntlength,effective_end));
      if (extension > (queryntlength - effective_end)/2) {
	/* Extension occupies more than 1/3 of sequence */
	debug2(printf("Proposed extension of %d is too long relative to queryntlength %d and effective_end %d\n",
		      extension,queryntlength,effective_end));
	extension = (queryntlength - effective_end)/3;
      }
#else
      extension = 0;
#endif
      if ((querysubseq = Sequence_subsequence(queryseq,effective_end-extension,queryntlength)) != NULL) {
	if ((querysubuc = Sequence_subsequence(queryuc,effective_end-extension,queryntlength)) != NULL) {
	  debug2(printf("5 margin <= 3 margin.  "));
	  debug2(printf("Beginning Stage1_compute on 3' margin from effective_end %d (%d..%d) (extension %d)\n",
			effective_end,effective_end-extension,queryntlength,extension));
	  debug2(Sequence_stdout(querysubseq,/*uppercasep*/true,wraplength,/*trimmedp*/true));

	  diagnostic = evaluate_query(&poorp,&repetitivep,Sequence_fullpointer(querysubuc),Sequence_fulllength(querysubuc),
				      Oligoindex_array_elt(oligoindices_major,0));
	  if (poorp == true || repetitivep == true) {
	    debug2(printf("Subsequence is poor or repetitive\n"));
	  } else {
	    gregions = Stage1_compute(&lowidentityp,querysubuc,indexdb_fwd,indexdb_rev,
				      /*indexdb_size_threshold*/100,chromosome_iit,
				      chrsubset_start,chrsubset_end,matchpool,
				      stutterhits,diagnostic,/*worker_stopwatch*/NULL,/*nbest*/10);
	    debug2(printf("C.  Performing Stage 3 with list length %d\n",List_length(stage3list)));
	    stage3list = stage3_from_gregions(stage3list,gregions,querysubseq,querysubuc,
#ifdef PMAP
					      queryntseq,
#endif
					      usersegment,stage2_alloc,oligoindices_major,oligoindices_minor,
					      pairpool,diagpool,cellpool,
					      dynprogL,dynprogM,dynprogR,/*worker_stopwatch*/NULL);
#ifdef DEBUG2
	    for (p = stage3list; p != NULL; p = List_next(p)) {
	      stage3 = (Stage3_T) List_head(p);
	      printf("%d..%d, %u..%u\n",
		     Stage3_querystart(stage3),Stage3_queryend(stage3),
		     Stage3_genomicstart(stage3),Stage3_genomicend(stage3));
	    }
#endif
	  }
	  Diagnostic_free(&diagnostic);

	  /* Above function frees gregions */
	  Sequence_free(&querysubuc);
	}
	Sequence_free(&querysubseq);
      }

      /* And recompute on original part, just in case stage 1 was led astray by the ends */
      if ((querysubseq = Sequence_subsequence(queryseq,0,effective_end)) != NULL) {
	if ((querysubuc = Sequence_subsequence(queryuc,0,effective_end)) != NULL) {
	  debug2(printf("Recomputing on original part.  "));
	  debug2(printf("Beginning Stage1_compute on 3' margin from effective_end %d (%d..%d), extension %d\n",
			effective_end,0,effective_end,extension));
	  debug2(Sequence_stdout(querysubseq,/*uppercasep*/true,wraplength,/*trimmedp*/true));

	  diagnostic = evaluate_query(&poorp,&repetitivep,Sequence_fullpointer(querysubuc),Sequence_fulllength(querysubuc),
				      Oligoindex_array_elt(oligoindices_major,0));
	  if (poorp == true || repetitivep == true) {
	    debug2(printf("Subsequence is poor or repetitive\n"));
	  } else {
	    gregions = Stage1_compute(&lowidentityp,querysubuc,indexdb_fwd,indexdb_rev,
				      /*indexdb_size_threshold*/100,chromosome_iit,
				      chrsubset_start,chrsubset_end,matchpool,
				      stutterhits,diagnostic,/*worker_stopwatch*/NULL,/*nbest*/10);
	    debug2(printf("D.  Performing Stage 3 with list length %d\n",List_length(stage3list)));
	    stage3list = stage3_from_gregions(stage3list,gregions,querysubseq,querysubuc,
#ifdef PMAP
					      queryntseq,
#endif
					      usersegment,stage2_alloc,oligoindices_major,oligoindices_minor,
					      pairpool,diagpool,cellpool,
					      dynprogL,dynprogM,dynprogR,/*worker_stopwatch*/NULL);
#ifdef DEBUG2
	    for (p = stage3list; p != NULL; p = List_next(p)) {
	      stage3 = (Stage3_T) List_head(p);
	      printf("%d..%d, %u..%u\n",
		     Stage3_querystart(stage3),Stage3_queryend(stage3),
		     Stage3_genomicstart(stage3),Stage3_genomicend(stage3));
	    }
#endif
	  }
	  Diagnostic_free(&diagnostic);

	  /* Above function frees gregions */
	  Sequence_free(&querysubuc);

	}
	Sequence_free(&querysubseq);
      }

      List_free(&nonjoinable);
      debug2(printf("Running local_separate_paths\n"));
      nonjoinable = local_separate_paths(&stage3array_sub1,&npaths_sub1,&stage3array_sub2,&npaths_sub2,
					 stage3list);
      debug2(printf("local: npaths_sub1 %d, npaths_sub2 %d, nonjoinable %d\n",
		    npaths_sub1,npaths_sub2,List_length(nonjoinable)));
    }
  }

  *mergedp = false;
  if (npaths_sub1 == 0 || npaths_sub2 == 0) {
    /* Skip */

  } else if (Chimera_bestpath(&five_score,&three_score,&chimerapos,&chimeraequivpos,&bestfrom,&bestto,
			      stage3array_sub1,npaths_sub1,stage3array_sub2,npaths_sub2,queryntlength,
			      CHIMERA_SLOP,/*localp*/true) == false) {
    /* Skip */
    debug2(printf("Chimera_bestpath returns false\n"));

    FREE(stage3array_sub2);
    FREE(stage3array_sub1);

  } else {
    from = stage3array_sub1[bestfrom];
    to = stage3array_sub2[bestto];
    debug2(printf("Chimera_bestpath returns bestfrom %d (%d..%d, %u..%u) to bestto %d (%d..%d, %u..%u)\n",
		  bestfrom,Stage3_querystart(from),Stage3_queryend(from),Stage3_genomicstart(from),Stage3_genomicend(from),
		  bestto,Stage3_querystart(to),Stage3_queryend(to),Stage3_genomicstart(to),Stage3_genomicend(to)));

    breakpoint = find_breakpoint(&chimera_cdna_direction,&chimerapos,&chimeraequivpos,&exonexonpos,
				 &donor1,&donor2,&acceptor2,&acceptor1,
				 &donor_watsonp,&acceptor_watsonp,&donor_prob,&acceptor_prob,from,to,
#ifdef PMAP
				 queryntseq,
#endif
				 queryseq,queryuc,queryntlength,
				 genomecomp,genomecomp_alt,chromosome_iit,pairpool);
    debug2(printf("find_breakpoint returns %d\n",breakpoint));

    /* Check to see if we can merge chimeric parts */
    debug2(printf("Before Stage3_mergeable, bestfrom is %p, query %d..%d\n",
		  from,Stage3_querystart(from),Stage3_queryend(from)));
    debug2(printf("Before Stage3_mergeable, bestto is %p, query %d..%d\n",
		  to,Stage3_querystart(to),Stage3_queryend(to)));

    if (Stage3_mergeable(from,to,breakpoint,queryntlength) == true) {
      debug2(printf("Mergeable! -- Merging left and right as a readthrough\n"));
      List_free(&stage3list);
      stage3list = merge_left_and_right_readthrough(&(*mergedp),stage3array_sub1,npaths_sub1,bestfrom,
						    stage3array_sub2,npaths_sub2,bestto,
						    nonjoinable,breakpoint,queryntlength,
#ifdef PMAP
						    /*queryaaseq_ptr*/Sequence_fullpointer(queryseq),
						    /*queryseq_ptr*/Sequence_fullpointer(queryntseq),
						    /*queryuc_ptr*/Sequence_fullpointer(queryntseq),
#else
						    /*queryseq_ptr*/Sequence_fullpointer(queryseq),
						    /*queryuc_ptr*/Sequence_fullpointer(queryuc),
#endif
						    pairpool,dynprogL,dynprogM,dynprogR,
						    oligoindices_minor,diagpool,cellpool);
    }

    FREE(stage3array_sub2);
    FREE(stage3array_sub1);
  }

  List_free(&nonjoinable);

  debug2(printf("check_for_local returning list of length %d\n",List_length(stage3list)));

  /* stage3list = stage3list_remove_empties(stage3list); */

#if 0
  /* Should be handled by apply_stage3 loop */
  /* Needed after calls to stage3_from_gregions */
  Stage3_recompute_goodness(stage3list);
  stage3list = stage3list_remove_duplicates(stage3list);
#endif

  return stage3list;
}


static List_T
check_for_chimera (bool *mergedp, Chimera_T *chimera, List_T stage3list, int effective_start, int effective_end,
		   Sequence_T queryseq, Sequence_T queryuc,
#ifdef PMAP
		   Sequence_T queryntseq,
#endif
		   int queryntlength, Sequence_T usersegment, Stage2_alloc_T stage2_alloc,
		   Oligoindex_array_T oligoindices_major, Oligoindex_array_T oligoindices_minor,
		   Matchpool_T matchpool, Pairpool_T pairpool, Diagpool_T diagpool, Cellpool_T cellpool,
		   Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR) {
  List_T gregions = NULL, nonjoinable = NULL, p;
  Stage3_T *stage3array_sub1 = NULL, *stage3array_sub2 = NULL, from, to, stage3;
  Sequence_T querysubseq = NULL, querysubuc = NULL;
  Diagnostic_T diagnostic;
  int bestfrom, bestto;
  int five_margin, three_margin, five_score = 0, three_score = 0;
  int extension;
  int npaths_sub1 = 0, npaths_sub2 = 0;
  bool lowidentityp, poorp, repetitivep;

  int max_single_goodness, chimeric_goodness, penalty, matches0, matches1;
  int breakpoint, chimerapos, chimeraequivpos, exonexonpos;
  int chimera_cdna_direction;
  char donor1, donor2, acceptor2, acceptor1;
  bool donor_watsonp, acceptor_watsonp;
  double donor_prob, acceptor_prob;
  

#ifdef PMAP
  five_margin = effective_start - 3*Sequence_trim_start(queryseq);
  three_margin = 3*Sequence_trim_end(queryseq) - effective_end;
  debug2(printf("Margins are %d = %d - %d on the 5' end and %d = %d - %d on the 3' end\n",
		five_margin,effective_start,3*Sequence_trim_start(queryseq),
		three_margin,3*Sequence_trim_end(queryseq),effective_end));
#else
  five_margin = effective_start - Sequence_trim_start(queryseq);
  three_margin = Sequence_trim_end(queryseq) - effective_end;
  debug2(printf("Margins are %d = %d - %d on the 5' end and %d = %d - %d on the 3' end\n",
		five_margin,effective_start,Sequence_trim_start(queryseq),
		three_margin,Sequence_trim_end(queryseq),effective_end));
#endif

#ifdef DEBUG2A
  for (p = stage3list; p != NULL; p = List_next(p)) {
    stage3 = (Stage3_T) List_head(p);
    Pair_dump_array(Stage3_pairarray(stage3),Stage3_npairs(stage3),/*zerobasedp*/true);
    printf("\n");
  }
#endif

  /* Stage3_recompute_goodness(stage3list); */
  max_single_goodness = 0;
  for (p = stage3list; p != NULL; p = List_next(p)) {
    stage3 = (Stage3_T) List_head(p);
    if (Stage3_goodness(stage3) > max_single_goodness) {
      max_single_goodness = Stage3_goodness(stage3);
    }
  }
  debug2(printf("max single goodness = %d\n",max_single_goodness));


  /* List_free(&nonjoinable); */
  debug2(printf("Running distant_separate_paths\n"));
  nonjoinable = distant_separate_paths(&stage3array_sub1,&npaths_sub1,&stage3array_sub2,&npaths_sub2,
				       stage3list);
  debug2(printf("chimera: npaths_sub1 %d, npaths_sub2 %d, nonjoinable %d\n",
		npaths_sub1,npaths_sub2,List_length(nonjoinable)));

  if (npaths_sub1 == 0 && npaths_sub2 == 0) {
    /* Need to compute on margin explicitly */
    if (five_margin < chimera_margin && three_margin < chimera_margin) {
      debug2(printf("Insufficient margins\n"));
    } else if (five_margin > three_margin) {
      extension = CHIMERA_SLOP;
      debug2(printf("Comparing extension %d with %d = (effective_start %d)/2\n",
		    extension,effective_start/2,effective_start));
      if (extension > effective_start/2) {
	/* Extension occupies more than 1/3 of sequence */
	debug2(printf("Proposed extension of %d is too long relative to effective_start %d\n",extension,effective_start));
	extension = effective_start/3;
      }
      if ((querysubseq = Sequence_subsequence(queryseq,0,effective_start+extension)) != NULL) {
	if ((querysubuc = Sequence_subsequence(queryuc,0,effective_start+extension)) != NULL) {
	  debug2(printf("5 margin > 3 margin.  "));
	  debug2(printf("Beginning Stage1_compute on 5' margin from effective_start %d (%d..%d)\n",
			effective_start,0,effective_start+extension));
	  debug2a(Sequence_stdout(querysubseq,/*uppercasep*/true,wraplength,/*trimmedp*/true));

	  diagnostic = evaluate_query(&poorp,&repetitivep,Sequence_fullpointer(querysubuc),Sequence_fulllength(querysubuc),
				      Oligoindex_array_elt(oligoindices_major,0));
	  if (poorp == true || repetitivep == true) {
	    debug2(printf("Subsequence is poor or repetitive\n"));
	  } else {
	    gregions = Stage1_compute(&lowidentityp,querysubuc,indexdb_fwd,indexdb_rev,
				      /*indexdb_size_threshold*/100,chromosome_iit,
				      chrsubset_start,chrsubset_end,matchpool,
				      stutterhits,diagnostic,/*worker_stopwatch*/NULL,/*nbest*/10);
	    debug2(printf("A.  Performing Stage 3 starting with list length %d\n",List_length(stage3list)));
	    stage3list = stage3_from_gregions(stage3list,gregions,querysubseq,querysubuc,
#ifdef PMAP
					      queryntseq,
#endif
					      usersegment,stage2_alloc,oligoindices_major,oligoindices_minor,
					      pairpool,diagpool,cellpool,
					      dynprogL,dynprogM,dynprogR,/*worker_stopwatch*/NULL);
#ifdef DEBUG2
	    for (p = stage3list; p != NULL; p = List_next(p)) {
	      stage3 = (Stage3_T) List_head(p);
	      printf("%d..%d, %u..%u\n",
		     Stage3_querystart(stage3),Stage3_queryend(stage3),
		     Stage3_genomicstart(stage3),Stage3_genomicend(stage3));
	    }
#endif
	  }
	  Diagnostic_free(&diagnostic);

	  /* Above function frees gregions */
	  Sequence_free(&querysubuc);
	}
	Sequence_free(&querysubseq);
      }

      /* And recompute on original part, just in case stage 1 was led astray by the ends */
      if ((querysubseq = Sequence_subsequence(queryseq,effective_start,queryntlength)) != NULL) {
	if ((querysubuc = Sequence_subsequence(queryuc,effective_start,queryntlength)) != NULL) {
	  debug2(printf("Recomputing on original part.  "));
	  debug2(printf("Beginning Stage1_compute on 5' margin from effective_start %d (%d..%d)\n",
			effective_start,effective_start,queryntlength));
	  debug2a(Sequence_stdout(querysubseq,/*uppercasep*/true,wraplength,/*trimmedp*/true));

	  diagnostic = evaluate_query(&poorp,&repetitivep,Sequence_fullpointer(querysubuc),Sequence_fulllength(querysubuc),
				      Oligoindex_array_elt(oligoindices_major,0));
	  if (poorp == true || repetitivep == true) {
	    debug2(printf("Subsequence is poor or repetitive\n"));
	  } else {
	    gregions = Stage1_compute(&lowidentityp,querysubuc,indexdb_fwd,indexdb_rev,
				      /*indexdb_size_threshold*/100,chromosome_iit,
				      chrsubset_start,chrsubset_end,matchpool,
				      stutterhits,diagnostic,/*worker_stopwatch*/NULL,/*nbest*/10);
	    debug2(printf("B.  Performing Stage 3 starting with list length %d\n",List_length(stage3list)));
	    stage3list = stage3_from_gregions(stage3list,gregions,querysubseq,querysubuc,
#ifdef PMAP
					      queryntseq,
#endif
					      usersegment,stage2_alloc,oligoindices_major,oligoindices_minor,
					      pairpool,diagpool,cellpool,
					      dynprogL,dynprogM,dynprogR,/*worker_stopwatch*/NULL);
#ifdef DEBUG2
	    for (p = stage3list; p != NULL; p = List_next(p)) {
	      stage3 = (Stage3_T) List_head(p);
	      printf("%d..%d, %u..%u\n",
		     Stage3_querystart(stage3),Stage3_queryend(stage3),
		     Stage3_genomicstart(stage3),Stage3_genomicend(stage3));
	    }
#endif
	  }
	  Diagnostic_free(&diagnostic);

	  /* Above function frees gregions */
	  Sequence_free(&querysubuc);
	}
	Sequence_free(&querysubseq);
      }

      List_free(&nonjoinable);
      debug2(printf("Running distant_separate_paths\n"));
      nonjoinable = distant_separate_paths(&stage3array_sub1,&npaths_sub1,&stage3array_sub2,&npaths_sub2,
					   stage3list);
      debug2(printf("chimera: npaths_sub1 %d, npaths_sub2 %d, nonjoinable %d\n",
		    npaths_sub1,npaths_sub2,List_length(nonjoinable)));

    } else {
      extension = CHIMERA_SLOP;
      debug2(printf("Comparing extension %d with %d = (queryntlength %d - effective_end %d)/2\n",
		    extension,(queryntlength-effective_end)/2,queryntlength,effective_end));
      if (extension > (queryntlength - effective_end)/2) {
	/* Extension occupies more than 1/3 of sequence */
	debug2(printf("Proposed extension of %d is too long relative to queryntlength %d and effective_end %d\n",
		      extension,queryntlength,effective_end));
	extension = (queryntlength - effective_end)/3;
      }
      if ((querysubseq = Sequence_subsequence(queryseq,effective_end-extension,queryntlength)) != NULL) {
	if ((querysubuc = Sequence_subsequence(queryuc,effective_end-extension,queryntlength)) != NULL) {
	  debug2(printf("5 margin <= 3 margin.  "));
	  debug2(printf("Beginning Stage1_compute on 3' margin from effective_end %d (%d..%d)\n",
			effective_end,effective_end-extension,queryntlength));
	  debug2(Sequence_stdout(querysubseq,/*uppercasep*/true,wraplength,/*trimmedp*/true));

	  diagnostic = evaluate_query(&poorp,&repetitivep,Sequence_fullpointer(querysubuc),Sequence_fulllength(querysubuc),
				      Oligoindex_array_elt(oligoindices_major,0));
	  if (poorp == true || repetitivep == true) {
	    debug2(printf("Subsequence is poor or repetitive\n"));
	  } else {
	    gregions = Stage1_compute(&lowidentityp,querysubuc,indexdb_fwd,indexdb_rev,
				      /*indexdb_size_threshold*/100,chromosome_iit,
				      chrsubset_start,chrsubset_end,matchpool,
				      stutterhits,diagnostic,/*worker_stopwatch*/NULL,/*nbest*/10);
	    debug2(printf("C.  Performing Stage 3 with list length %d\n",List_length(stage3list)));
	    stage3list = stage3_from_gregions(stage3list,gregions,querysubseq,querysubuc,
#ifdef PMAP
					      queryntseq,
#endif
					      usersegment,stage2_alloc,oligoindices_major,oligoindices_minor,
					      pairpool,diagpool,cellpool,
					      dynprogL,dynprogM,dynprogR,/*worker_stopwatch*/NULL);
#ifdef DEBUG2
	    for (p = stage3list; p != NULL; p = List_next(p)) {
	      stage3 = (Stage3_T) List_head(p);
	      printf("%d..%d, %u..%u\n",
		     Stage3_querystart(stage3),Stage3_queryend(stage3),
		     Stage3_genomicstart(stage3),Stage3_genomicend(stage3));
	    }
#endif
	  }
	  Diagnostic_free(&diagnostic);

	  /* Above function frees gregions */
	  Sequence_free(&querysubuc);
	}
	Sequence_free(&querysubseq);
      }

      /* And recompute on original part, just in case stage 1 was led astray by the ends */
      if ((querysubseq = Sequence_subsequence(queryseq,0,effective_end)) != NULL) {
	if ((querysubuc = Sequence_subsequence(queryuc,0,effective_end)) != NULL) {
	  debug2(printf("Recomputing on original part.  "));
	  debug2(printf("Beginning Stage1_compute on 3' margin from effective_end %d (%d..%d)\n",
			effective_end,0,effective_end));
	  debug2(Sequence_stdout(querysubseq,/*uppercasep*/true,wraplength,/*trimmedp*/true));

	  diagnostic = evaluate_query(&poorp,&repetitivep,Sequence_fullpointer(querysubuc),Sequence_fulllength(querysubuc),
				      Oligoindex_array_elt(oligoindices_major,0));
	  if (poorp == true || repetitivep == true) {
	    debug2(printf("Subsequence is poor or repetitive\n"));
	  } else {
	    gregions = Stage1_compute(&lowidentityp,querysubuc,indexdb_fwd,indexdb_rev,
				      /*indexdb_size_threshold*/100,chromosome_iit,
				      chrsubset_start,chrsubset_end,matchpool,
				      stutterhits,diagnostic,/*worker_stopwatch*/NULL,/*nbest*/10);
	    debug2(printf("D.  Performing Stage 3 with list length %d\n",List_length(stage3list)));
	    stage3list = stage3_from_gregions(stage3list,gregions,querysubseq,querysubuc,
#ifdef PMAP
					      queryntseq,
#endif
					      usersegment,stage2_alloc,oligoindices_major,oligoindices_minor,
					      pairpool,diagpool,cellpool,
					      dynprogL,dynprogM,dynprogR,/*worker_stopwatch*/NULL);
#ifdef DEBUG2
	    for (p = stage3list; p != NULL; p = List_next(p)) {
	      stage3 = (Stage3_T) List_head(p);
	      printf("%d..%d, %u..%u\n",
		     Stage3_querystart(stage3),Stage3_queryend(stage3),
		     Stage3_genomicstart(stage3),Stage3_genomicend(stage3));
	    }
#endif
	  }
	  Diagnostic_free(&diagnostic);

	  /* Above function frees gregions */
	  Sequence_free(&querysubuc);

	}
	Sequence_free(&querysubseq);
      }

      List_free(&nonjoinable);
      debug2(printf("Running distant_separate_paths\n"));
      nonjoinable = distant_separate_paths(&stage3array_sub1,&npaths_sub1,&stage3array_sub2,&npaths_sub2,
					   stage3list);
      debug2(printf("chimera: npaths_sub1 %d, npaths_sub2 %d, nonjoinable %d\n",
		    npaths_sub1,npaths_sub2,List_length(nonjoinable)));
    }
  }

  *mergedp = false;
  *chimera = (Chimera_T) NULL;
  if (npaths_sub1 == 0 || npaths_sub2 == 0) {
    /* Skip */

  } else if (Chimera_bestpath(&five_score,&three_score,&chimerapos,&chimeraequivpos,&bestfrom,&bestto,
			      stage3array_sub1,npaths_sub1,stage3array_sub2,npaths_sub2,queryntlength,
			      CHIMERA_SLOP,/*localp*/false) == false) {
    /* Skip */
    debug2(printf("Chimera_bestpath returns false, so skipping\n"));
    FREE(stage3array_sub2);
    FREE(stage3array_sub1);

  } else {
    from = stage3array_sub1[bestfrom];
    to = stage3array_sub2[bestto];
    debug2(printf("Chimera_bestpath returns bestfrom %d (%d..%d, %u..%u) to bestto %d (%d..%d, %u..%u)\n",
		  bestfrom,Stage3_querystart(from),Stage3_queryend(from),Stage3_genomicstart(from),Stage3_genomicend(from),
		  bestto,Stage3_querystart(to),Stage3_queryend(to),Stage3_genomicstart(to),Stage3_genomicend(to)));

    chimeric_goodness = Stage3_chimeric_goodness(&matches0,&matches1,from,to,chimerapos);
    debug2(printf("chimeric goodness = %d\n",chimeric_goodness));
    
    penalty = CHIMERA_PENALTY;
    if (chimera_margin < penalty) {
      /* User is looking for higher sensitivity */
      penalty = chimera_margin;
    }

    if (chimeric_goodness < max_single_goodness + penalty) {
      debug2(printf("chimeric goodness not good enough relative to max_single_goodness %d and penalty %d\n",
		    max_single_goodness,penalty));

    } else if ((breakpoint = find_breakpoint(&chimera_cdna_direction,&chimerapos,&chimeraequivpos,&exonexonpos,
					     &donor1,&donor2,&acceptor2,&acceptor1,
					     &donor_watsonp,&acceptor_watsonp,&donor_prob,&acceptor_prob,from,to,
#ifdef PMAP
					     queryntseq,
#endif
					     queryseq,queryuc,queryntlength,
					     genomecomp,genomecomp_alt,chromosome_iit,pairpool)) < 0) {
      debug2(printf("find_breakpoint returns no value\n"));

    } else {
      debug2(printf("find_breakpoint returns %d\n",breakpoint));

      /* Check to see if we can merge chimeric parts */
      debug2(printf("Before Stage3_mergeable, bestfrom is %p, query %d..%d\n",
		    from,Stage3_querystart(from),Stage3_queryend(from)));
      debug2(printf("Before Stage3_mergeable, bestto is %p, query %d..%d\n",
		    to,Stage3_querystart(to),Stage3_queryend(to)));

      if (Stage3_mergeable(from,to,breakpoint,queryntlength) == false &&
	  Stage3_test_bounds(from,0,chimeraequivpos+chimera_overlap) == true &&
	  Stage3_test_bounds(to,chimerapos+1-chimera_overlap,queryntlength) == true &&
	  Stage3_merge_chimera(/*best0*/from,/*best1*/to,
			       /*minpos1*/0,/*maxpos1*/breakpoint,
			       /*minpos2*/breakpoint+1,/*maxpos2*/queryntlength,
#ifdef PMAP
			       Sequence_fullpointer(queryntseq),Sequence_fullpointer(queryntseq),
#else
			       Sequence_fullpointer(queryseq),Sequence_fullpointer(queryuc),
#endif
			       pairpool,dynprogL,dynprogR,maxpeelback) == true) {
	
	/* if maxpaths_report == 1, then don't want distant chimeras */
	if (maxpaths_report != 1) {
	  debug2(printf("Not mergeable -- Merging left and right as a transloc\n"));
	  *chimera = Chimera_new(from,to,chimerapos,chimeraequivpos,exonexonpos,chimera_cdna_direction,
				 donor1,donor2,acceptor2,acceptor1,donor_watsonp,acceptor_watsonp,
				 donor_prob,acceptor_prob);
	  List_free(&stage3list);

	  debug2(printf("Before merge_left_and_right_transloc, bestfrom is %p, query %d..%d\n",
			from,Stage3_querystart(from),Stage3_queryend(from)));
	  debug2(printf("Before merge_left_and_right_transloc, bestto is %p, query %d..%d\n",
			to,Stage3_querystart(to),Stage3_queryend(to)));
	  
	  stage3list = merge_left_and_right_transloc(stage3array_sub1,npaths_sub1,bestfrom,
						     stage3array_sub2,npaths_sub2,bestto,
						     nonjoinable);
	}
      }
    }

    FREE(stage3array_sub2);
    FREE(stage3array_sub1);
  }

  List_free(&nonjoinable);

  debug2(printf("check_for_chimera returning list of length %d\n",List_length(stage3list)));

#if 0
  /* Should be handled by apply_stage3 loop */
  /* Needed after calls to stage3_from_gregions */
  Stage3_recompute_goodness(stage3list);
  stage3list = stage3list_remove_duplicates(stage3list);
#endif

  return stage3list;
}


static List_T
merge_middlepieces (List_T stage3list, Stage3_T from, Stage3_T to,
		    List_T middlepieces, Stage3_T middle, bool mergeableAp, bool mergeableBp,
		    int breakpointA, int breakpointB, Sequence_T queryseq,
#ifdef PMAP
		    Sequence_T queryntseq,
#endif
		    Sequence_T queryuc, int queryntlength,
		    Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		    Oligoindex_array_T oligoindices_minor, Diagpool_T diagpool, Cellpool_T cellpool) {
  List_T newstage3list = NULL, merged;
  List_T nonjoinable, r;
  bool mergedAp, mergedBp;
  Stage3_T stage3;

  nonjoinable = (List_T) NULL;
  for (r = stage3list; r != NULL; r = List_next(r)) {
    stage3 = (Stage3_T) List_head(r);
    if (stage3 == from) {
      /* Skip */
    } else if (stage3 == to) {
      /* Skip */
    } else {
      nonjoinable = List_push(nonjoinable,(void *) stage3);
    }
  }
  

  if (mergeableAp == true && mergeableBp == true) {
    merged = merge_left_and_right_readthrough(&mergedAp,/*stage3array_sub1*/&from,/*npaths_sub1*/1,/*bestfrom*/0,
					      /*stage3array_sub2*/&middle,/*npaths_sub2*/1,/*bestto*/0,
					      /*nonjoinable*/NULL,breakpointA,queryntlength,
#ifdef PMAP
					      /*queryaaseq_ptr*/Sequence_fullpointer(queryseq),
					      /*queryseq_ptr*/Sequence_fullpointer(queryntseq),
					      /*queryuc_ptr*/Sequence_fullpointer(queryntseq),
#else
					      /*queryseq_ptr*/Sequence_fullpointer(queryseq),
					      /*queryuc_ptr*/Sequence_fullpointer(queryuc),
#endif
					      pairpool,dynprogL,dynprogM,dynprogR,
					      oligoindices_minor,diagpool,cellpool);
    List_free(&merged);

    newstage3list = merge_left_and_right_readthrough(&mergedBp,/*stage3array_sub1*/&from,/*npaths_sub1*/1,/*bestfrom*/0,
						     /*stage3array_sub2*/&to,/*npaths_sub2*/1,/*bestto*/0,
						     nonjoinable,breakpointB,queryntlength,
#ifdef PMAP
						     /*queryaaseq_ptr*/Sequence_fullpointer(queryseq),
						     /*queryseq_ptr*/Sequence_fullpointer(queryntseq),
						     /*queryuc_ptr*/Sequence_fullpointer(queryntseq),
#else
						     /*queryseq_ptr*/Sequence_fullpointer(queryseq),
						     /*queryuc_ptr*/Sequence_fullpointer(queryuc),
#endif
						     pairpool,dynprogL,dynprogM,dynprogR,
						     oligoindices_minor,diagpool,cellpool);

#ifndef PMAP
    Stage3_guess_cdna_direction(from);
#endif
    List_free(&stage3list);

  } else if (mergeableBp == true) {
    nonjoinable = List_push(nonjoinable,(void *) from);
    newstage3list = merge_left_and_right_readthrough(&mergedBp,/*stage3array_sub1*/&middle,/*npaths_sub1*/1,/*bestfrom*/0,
						     /*stage3array_sub2*/&to,/*npaths_sub2*/1,/*bestto*/0,
						     nonjoinable,breakpointB,queryntlength,
#ifdef PMAP
						     /*queryaaseq_ptr*/Sequence_fullpointer(queryseq),
						     /*queryseq_ptr*/Sequence_fullpointer(queryntseq),
						     /*queryuc_ptr*/Sequence_fullpointer(queryntseq),
#else
						     /*queryseq_ptr*/Sequence_fullpointer(queryseq),
						     /*queryuc_ptr*/Sequence_fullpointer(queryuc),
#endif
						     pairpool,dynprogL,dynprogM,dynprogR,
						     oligoindices_minor,diagpool,cellpool);
#ifndef PMAP
    Stage3_guess_cdna_direction(middle);
#endif
    List_free(&stage3list);

  } else if (mergeableAp == true) {
    nonjoinable = List_push(nonjoinable,(void *) to);
    newstage3list = merge_left_and_right_readthrough(&mergedAp,/*stage3array_sub1*/&from,/*npaths_sub1*/1,/*bestfrom*/0,
						     /*stage3array_sub2*/&middle,/*npaths_sub2*/1,/*bestto*/0,
						     nonjoinable,breakpointA,queryntlength,
#ifdef PMAP
						     /*queryaaseq_ptr*/Sequence_fullpointer(queryseq),
						     /*queryseq_ptr*/Sequence_fullpointer(queryntseq),
						     /*queryuc_ptr*/Sequence_fullpointer(queryntseq),
#else
						     /*queryseq_ptr*/Sequence_fullpointer(queryseq),
						     /*queryuc_ptr*/Sequence_fullpointer(queryuc),
#endif
						     pairpool,dynprogL,dynprogM,dynprogR,
						     oligoindices_minor,diagpool,cellpool);
#ifndef PMAP
    Stage3_guess_cdna_direction(from);
#endif
    List_free(&stage3list);
    
  } else {
    newstage3list = stage3list;	/* Contains all entries from nonjoinable */
    newstage3list = List_push(newstage3list,(void *) middle);
  }

  for (r = middlepieces; r != NULL; r = List_next(r)) {
    stage3 = (Stage3_T) List_head(r);
    if (stage3 == NULL) {
      /* Already freed */
    } else if (stage3 == middle) {
      /* Don't add again */
    } else {
      newstage3list = List_push(newstage3list,stage3);
    }
  }

  List_free(&nonjoinable);
  return newstage3list;
}



/* Returns stage3list */
static List_T
check_middle_piece_local (bool *foundp, List_T stage3list, Sequence_T queryseq, Sequence_T queryuc,
#ifdef PMAP
			  Sequence_T queryntseq,
#endif
			  int queryntlength, Stage2_alloc_T stage2_alloc,
			  Oligoindex_array_T oligoindices_major, Oligoindex_array_T oligoindices_minor,
			  Pairpool_T pairpool, Diagpool_T diagpool, Cellpool_T cellpool,
			  Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR) {
  Sequence_T querysubseq = NULL, querysubuc = NULL;
  int npaths, i, j;
  Stage3_T from = NULL, to = NULL, middle = NULL;
  Stage3_T *by_queryend, *by_querystart;
  List_T r;
  bool plusp;
  int querystart, queryend;
  Chrpos_T chrstart, chrend, chrlength;
  Univcoord_T chroffset, chrhigh;
  Chrnum_T chrnum;

  int breakpointA = 0, chimeraposA, chimeraequivposA, exonexonposA;
  char donorA1, donorA2, acceptorA2, acceptorA1;
  bool donor_watsonp_A, acceptor_watsonp_A;
  double donor_prob_A, acceptor_prob_A;

  int breakpointB = 0, chimeraposB, chimeraequivposB, exonexonposB;
  char donorB1, donorB2, acceptorB2, acceptorB1;
  bool donor_watsonp_B, acceptor_watsonp_B;
  double donor_prob_B, acceptor_prob_B;

  int chimera_cdna_direction_A, chimera_cdna_direction_B;
  bool mergeableAp, mergeableBp;

  List_T middlepieces;


#ifdef DEBUG2A
  for (p = stage3list; p != NULL; p = List_next(p)) {
    stage3 = (Stage3_T) List_head(p);
    Pair_dump_array(Stage3_pairarray(stage3),Stage3_npairs(stage3),/*zerobasedp*/true);
    printf("\n");
  }
#endif

  *foundp = false;

  by_queryend = (Stage3_T *) List_to_array_n(&npaths,stage3list);
  qsort(by_queryend,npaths,sizeof(Stage3_T),Stage3_queryend_cmp);

  by_querystart = (Stage3_T *) List_to_array_n(&npaths,stage3list);
  qsort(by_querystart,npaths,sizeof(Stage3_T),Stage3_querystart_cmp);

  j = 0;
  for (i = 0; i < npaths && *foundp == false; i++) {
    from = by_queryend[i];
    queryend = Stage3_queryend(from);

    while (j < npaths && Stage3_querystart(by_querystart[j]) < queryend) {
      j++;
    }
    j--;

    while (j >= 0 && Stage3_querystart(by_querystart[j]) > queryend) {
      j--;
    }
    j++;

    for ( ; j < npaths && *foundp == false; j++) {
      to = by_querystart[j];

      if (middle_piece_local_p(&querystart,&queryend,&chrstart,&chrend,
			       &chrnum,&chroffset,&chrhigh,&chrlength,&plusp,
			       from,to) == true) {
	debug2(printf("Found middle piece missing from %d to %d\n",i,j));

	if ((querysubseq = Sequence_subsequence(queryseq,querystart,queryend)) != NULL) {
	  if ((querysubuc = Sequence_subsequence(queryuc,querystart,queryend)) != NULL) {
	    debug2(printf("Performing Stage 3 on %d..%d against %u..%u\n",
			  querystart,queryend,chrstart,chrend));
	    if ((middlepieces = update_stage3list(/*stage3list*/NULL,querysubseq,
#ifdef PMAP
						  queryntseq,
#endif
						  querysubuc,stage2_alloc,oligoindices_major,oligoindices_minor,
						  pairpool,diagpool,cellpool,
						  /*straintype*/0,/*strain*/NULL,chrnum,
						  chroffset,chrhigh,chrlength,chrstart,chrend,plusp,/*genestrand*/0,
						  dynprogL,dynprogM,dynprogR,/*worker_stopwatch*/NULL)) != NULL) {
	      middlepieces = stage3list_sort(middlepieces);

	      /* 1.  Look first for middle piece that joins locally on both ends */
	      r = middlepieces;
	      mergeableAp = mergeableBp = false;
	      while (r != NULL && (mergeableAp == false || mergeableBp == false)) {
		middle = (Stage3_T) List_head(r);
		if (Chimera_local_join_p(from,middle,CHIMERA_SLOP) == true && Chimera_local_join_p(middle,to,CHIMERA_SLOP) == true) {
		  breakpointA = find_breakpoint(&chimera_cdna_direction_A,&chimeraposA,&chimeraequivposA,&exonexonposA,
						&donorA1,&donorA2,&acceptorA2,&acceptorA1,
						&donor_watsonp_A,&acceptor_watsonp_A,&donor_prob_A,&acceptor_prob_A,
						from,/*to*/middle,
#ifdef PMAP
						queryntseq,
#endif
						queryseq,queryuc,queryntlength,
						genomecomp,genomecomp_alt,chromosome_iit,pairpool);
		  breakpointB = find_breakpoint(&chimera_cdna_direction_B,&chimeraposB,&chimeraequivposB,&exonexonposB,
						&donorB1,&donorB2,&acceptorB2,&acceptorB1,
						&donor_watsonp_B,&acceptor_watsonp_B,&donor_prob_B,&acceptor_prob_B,
						/*from*/middle,to,
#ifdef PMAP
						queryntseq,
#endif
						queryseq,queryuc,queryntlength,
						genomecomp,genomecomp_alt,chromosome_iit,pairpool);

		  mergeableAp = Stage3_mergeable(from,/*to*/middle,breakpointA,queryntlength);
		  mergeableBp = Stage3_mergeable(/*from*/middle,to,breakpointB,queryntlength);
		}
		r = List_next(r);
	      }	/* End of while loop looking for dual merge */

	      if (mergeableAp == true && mergeableBp == true) {
		debug2(printf("Middle segment found and mergeable locally with both! -- Merging three as a readthrough.\n"));
		*foundp = true;
	      } else {
		/* 2.  Look for middle piece that joins locally on one end */
		r = middlepieces;
		mergeableAp = mergeableBp = false;
		while (r != NULL && mergeableAp == false && mergeableBp == false) {
		  middle = (Stage3_T) List_head(r);
		  if (Chimera_local_join_p(from,middle,CHIMERA_SLOP) == true && Chimera_local_join_p(middle,to,CHIMERA_SLOP) == true) {
		    breakpointA = find_breakpoint(&chimera_cdna_direction_A,&chimeraposA,&chimeraequivposA,&exonexonposA,
						  &donorA1,&donorA2,&acceptorA2,&acceptorA1,
						  &donor_watsonp_A,&acceptor_watsonp_A,&donor_prob_A,&acceptor_prob_A,
						  from,/*to*/middle,
#ifdef PMAP
						  queryntseq,
#endif
						  queryseq,queryuc,queryntlength,
						  genomecomp,genomecomp_alt,chromosome_iit,pairpool);
		    breakpointB = find_breakpoint(&chimera_cdna_direction_B,&chimeraposB,&chimeraequivposB,&exonexonposB,
						  &donorB1,&donorB2,&acceptorB2,&acceptorB1,
						  &donor_watsonp_B,&acceptor_watsonp_B,&donor_prob_B,&acceptor_prob_B,
						  /*from*/middle,to,
#ifdef PMAP
						  queryntseq,
#endif
						  queryseq,queryuc,queryntlength,
						  genomecomp,genomecomp_alt,chromosome_iit,pairpool);

		    mergeableAp = Stage3_mergeable(from,/*to*/middle,breakpointA,queryntlength);
		    mergeableBp = Stage3_mergeable(/*from*/middle,to,breakpointB,queryntlength);
		  }
		  r = List_next(r);
		} /* End of while loop looking for single merge */

		if (mergeableAp == true || mergeableBp == true) {
		  *foundp = true;
		}
	      }

	      stage3list = merge_middlepieces(stage3list,from,to,middlepieces,middle,mergeableAp,mergeableBp,
					      breakpointA,breakpointB,queryseq,
#ifdef PMAP
					      queryntseq,
#endif
					      queryuc,queryntlength,pairpool,dynprogL,dynprogM,dynprogR,
					      oligoindices_minor,diagpool,cellpool);
	      List_free(&middlepieces);
	    }

	    Sequence_free(&querysubuc);
	  }
	  Sequence_free(&querysubseq);
	}
      }
    }
  }

  FREE(by_querystart);
  FREE(by_queryend);

  return stage3list;
}


/* Returns stage3list */
static List_T
check_middle_piece_chimera (bool *foundp, List_T stage3list, Sequence_T queryseq, Sequence_T queryuc,
#ifdef PMAP
			    Sequence_T queryntseq,
#endif
			    int queryntlength, Sequence_T usersegment, Stage2_alloc_T stage2_alloc,
			    Oligoindex_array_T oligoindices_major, Oligoindex_array_T oligoindices_minor,
			    Matchpool_T matchpool, Pairpool_T pairpool, Diagpool_T diagpool, Cellpool_T cellpool,
			    Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR) {
  List_T newstage3list = NULL;
  Sequence_T querysubseq = NULL, querysubuc = NULL;
  int npaths, i, j;
  Stage3_T bestfrom, bestto, from, to, middle, stage3;
  Stage3_T *by_queryend, *by_querystart;
  List_T r;
  int querystart, queryend, maxdist, dist;

  int breakpointA, chimeraposA, chimeraequivposA, exonexonposA;
  char donorA1, donorA2, acceptorA2, acceptorA1;
  bool donor_watsonp_A, acceptor_watsonp_A;
  double donor_prob_A, acceptor_prob_A;

  int breakpointB, chimeraposB, chimeraequivposB, exonexonposB;
  char donorB1, donorB2, acceptorB2, acceptorB1;
  bool donor_watsonp_B, acceptor_watsonp_B;
  double donor_prob_B, acceptor_prob_B;

  int chimera_cdna_direction_A, chimera_cdna_direction_B;
  bool mergeableAp, mergeableBp, mergedAp, mergedBp;

  List_T nonjoinable = NULL, middlepieces = NULL;
  Diagnostic_T diagnostic;
  List_T gregions;
  bool lowidentityp, poorp, repetitivep;


#ifdef DEBUG2A
  for (p = stage3list; p != NULL; p = List_next(p)) {
    stage3 = (Stage3_T) List_head(p);
    Pair_dump_array(Stage3_pairarray(stage3),Stage3_npairs(stage3),/*zerobasedp*/true);
    printf("\n");
  }
#endif

  by_queryend = (Stage3_T *) List_to_array_n(&npaths,stage3list);
  qsort(by_queryend,npaths,sizeof(Stage3_T),Stage3_queryend_cmp);

  by_querystart = (Stage3_T *) List_to_array_n(&npaths,stage3list);
  qsort(by_querystart,npaths,sizeof(Stage3_T),Stage3_querystart_cmp);

  maxdist = 0;
  j = 0;
  for (i = 0; i < npaths; i++) {
    from = by_queryend[i];
    queryend = Stage3_queryend(from);

    while (j < npaths && Stage3_querystart(by_querystart[j]) < queryend) {
      j++;
    }
    j--;

    while (j >= 0 && Stage3_querystart(by_querystart[j]) > queryend) {
      j--;
    }
    j++;

    if (j < npaths) {
      /* Should have the first querystart just after queryend */
      to = by_querystart[j];

      if ((dist = Stage3_queryend(to) - Stage3_querystart(from)) > maxdist) {
	bestfrom = from;
	bestto = to;
	maxdist = dist;
      }
    }
  }

  FREE(by_querystart);
  FREE(by_queryend);


  *foundp = false;
  if (maxdist < CHIMERA_SLOP) {
    debug2(printf("maxdist %d < CHIMERA_SLOP %d\n",maxdist,CHIMERA_SLOP));
  } else {
    if (middle_piece_chimera_p(&querystart,&queryend,bestfrom,bestto) == true) {
      if ((querysubseq = Sequence_subsequence(queryseq,querystart,queryend)) != NULL) {
	if ((querysubuc = Sequence_subsequence(queryuc,querystart,queryend)) != NULL) {
	  debug2(printf("Performing Stage 3 on %d..%d\n",querystart,queryend));

	  diagnostic = evaluate_query(&poorp,&repetitivep,Sequence_fullpointer(querysubuc),
				      Sequence_fulllength(querysubuc),Oligoindex_array_elt(oligoindices_major,0));
	  if (poorp == true || repetitivep == true) {
	    debug2(printf("Subsequence is poor or repetitive\n"));
	  } else {
	    gregions = Stage1_compute(&lowidentityp,querysubuc,indexdb_fwd,indexdb_rev,
				      /*indexdb_size_threshold*/100,chromosome_iit,
				      chrsubset_start,chrsubset_end,matchpool,
				      stutterhits,diagnostic,/*worker_stopwatch*/NULL,/*nbest*/10);
	    debug2(printf("Performing Stage 3 starting with list length %d\n",List_length(stage3list)));
	    middlepieces = stage3_from_gregions(/*stage3list*/NULL,gregions,querysubseq,querysubuc,
#ifdef PMAP
						queryntseq,
#endif
						usersegment,stage2_alloc,oligoindices_major,oligoindices_minor,
						pairpool,diagpool,cellpool,
						dynprogL,dynprogM,dynprogR,/*worker_stopwatch*/NULL);
	  }
	  Diagnostic_free(&diagnostic);

	  /* Above function frees gregions */
	  Sequence_free(&querysubuc);
	}
	Sequence_free(&querysubseq);
      }
    }

    if (middlepieces != NULL) {
      middlepieces = stage3list_sort(middlepieces);

      r = middlepieces;
      mergeableAp = mergeableBp = false;
      while (r != NULL && mergeableAp == false && mergeableBp == false) {
	middle = (Stage3_T) List_head(r);
	if (middle != bestfrom && middle != bestto) {
	  if (Chimera_local_join_p(bestfrom,middle,CHIMERA_SLOP) == true) {
	    breakpointA = find_breakpoint(&chimera_cdna_direction_A,&chimeraposA,&chimeraequivposA,&exonexonposA,
					  &donorA1,&donorA2,&acceptorA2,&acceptorA1,
					  &donor_watsonp_A,&acceptor_watsonp_A,&donor_prob_A,&acceptor_prob_A,
					  bestfrom,/*to*/middle,
#ifdef PMAP
					  queryntseq,
#endif
					  queryseq,queryuc,queryntlength,
					  genomecomp,genomecomp_alt,chromosome_iit,pairpool);
	    mergeableAp = Stage3_mergeable(bestfrom,/*to*/middle,breakpointA,queryntlength);
	  }
	  if (Chimera_local_join_p(middle,bestto,CHIMERA_SLOP) == true) {
	    breakpointB = find_breakpoint(&chimera_cdna_direction_B,&chimeraposB,&chimeraequivposB,&exonexonposB,
					  &donorB1,&donorB2,&acceptorB2,&acceptorB1,
					  &donor_watsonp_B,&acceptor_watsonp_B,&donor_prob_B,&acceptor_prob_B,
					  /*from*/middle,to,
#ifdef PMAP
					  queryntseq,
#endif
					  queryseq,queryuc,queryntlength,
					  genomecomp,genomecomp_alt,chromosome_iit,pairpool);
	    mergeableBp = Stage3_mergeable(/*from*/middle,bestto,breakpointB,queryntlength);
	  }
	}
	r = List_next(r);
      }

      if (mergeableAp == true) {
	debug2(printf("Middle segment found and mergeable locally with from! -- Merging as a readthrough.  cdna_direction = %d\n",
		      chimera_cdna_direction_A));

	List_free(&nonjoinable);
	nonjoinable = (List_T) NULL;
	for (r = middlepieces; r != NULL; r = List_next(r)) {
	  stage3 = (Stage3_T) List_head(r);
	  if (stage3 == middle) {
	    /* Skip */
	  } else {
	    nonjoinable = List_push(nonjoinable,(void *) stage3);
	  }
	}
	List_free(&middlepieces);

	for (r = stage3list; r != NULL; r = List_next(r)) {
	  stage3 = (Stage3_T) List_head(r);
	  if (stage3 == bestfrom) {
	    /* Skip */
	  } else {
	    nonjoinable = List_push(nonjoinable,(void *) stage3);
	  }
	}

	newstage3list =
	  merge_left_and_right_readthrough(&mergedAp,/*stage3array_sub1*/&bestfrom,/*npaths_sub1*/1,/*bestfrom*/0,
					   /*stage3array_sub2*/&middle,/*npaths_sub2*/1,/*bestto*/0,
					   nonjoinable,breakpointA,queryntlength,
#ifdef PMAP
					   /*queryaaseq_ptr*/Sequence_fullpointer(queryseq),
					   /*queryseq_ptr*/Sequence_fullpointer(queryntseq),
					   /*queryuc_ptr*/Sequence_fullpointer(queryntseq),
#else
					   /*queryseq_ptr*/Sequence_fullpointer(queryseq),
					   /*queryuc_ptr*/Sequence_fullpointer(queryuc),
#endif
					   pairpool,dynprogL,dynprogM,dynprogR,
					   oligoindices_minor,diagpool,cellpool);
#ifndef PMAP
	Stage3_guess_cdna_direction(from);
#endif

	List_free(&nonjoinable);
	if (mergedAp == true) {
	  *foundp = true;
	}

      } else if (mergeableBp == true) {
	debug2(printf("Middle segment found and mergeable locally with to! -- Merging as a readthrough.  cdna_direction = %d\n",
		      chimera_cdna_direction_B));

	List_free(&nonjoinable);
	nonjoinable = (List_T) NULL;
	for (r = middlepieces; r != NULL; r = List_next(r)) {
	  stage3 = (Stage3_T) List_head(r);
	  if (stage3 == middle) {
	    /* Skip */
	  } else {
	    nonjoinable = List_push(nonjoinable,(void *) stage3);
	  }
	}
	List_free(&middlepieces);

	for (r = stage3list; r != NULL; r = List_next(r)) {
	  stage3 = (Stage3_T) List_head(r);
	  if (stage3 == bestto) {
	    /* Skip */
	  } else {
	    nonjoinable = List_push(nonjoinable,(void *) stage3);
	  }
	}

	newstage3list =
	  merge_left_and_right_readthrough(&mergedBp,/*stage3array_sub1*/&middle,/*npaths_sub1*/1,/*bestfrom*/0,
					   /*stage3array_sub2*/&bestto,/*npaths_sub2*/1,/*bestto*/0,
					   nonjoinable,breakpointB,queryntlength,
#ifdef PMAP
					   /*queryaaseq_ptr*/Sequence_fullpointer(queryseq),
					   /*queryseq_ptr*/Sequence_fullpointer(queryntseq),
					   /*queryuc_ptr*/Sequence_fullpointer(queryntseq),
#else
					   /*queryseq_ptr*/Sequence_fullpointer(queryseq),
					   /*queryuc_ptr*/Sequence_fullpointer(queryuc),
#endif
					   pairpool,dynprogL,dynprogM,dynprogR,
					   oligoindices_minor,diagpool,cellpool);

#ifndef PMAP
	Stage3_guess_cdna_direction(middle);
#endif

	List_free(&nonjoinable);
	if (mergedBp == true) {
	  *foundp = true;
	}

      } else {
	debug2(printf("Middle segment found but notmergeable\n"));
	for (r = middlepieces; r != NULL; r = List_next(r)) {
	  middle = (Stage3_T) List_head(r);
	  if (middle != NULL) {
	    Stage3_free(&middle);
	  }
	}
	List_free(&middlepieces);
      }

    }
  }

  if (newstage3list == NULL) {
    return stage3list;
  } else {
    List_free(&stage3list);
    return newstage3list;
  }
}



static List_T
apply_stage3 (bool *mergedp, Chimera_T *chimera, List_T gregions, Sequence_T queryseq, Sequence_T queryuc,
#ifdef PMAP
	      Sequence_T queryntseq,
#endif
	      Sequence_T usersegment, Stage2_alloc_T stage2_alloc,
	      Oligoindex_array_T oligoindices_major, Oligoindex_array_T oligoindices_minor,
	      Matchpool_T matchpool, Pairpool_T pairpool, Diagpool_T diagpool, Cellpool_T cellpool,
	      Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR, Stopwatch_T worker_stopwatch) {
#ifdef DEBUG2
  List_T p;
#endif
  List_T stage3list;
  Stage3_T nonchimericbest, chimera1, chimera2;
  bool testlocalp, testchimerap, foundp;
  int effective_start, effective_end;
  int queryntlength;
  int iter;

  
  *mergedp = false;
  *chimera = NULL;

  debug(printf("Calling stage3_from_gregions\n"));
  stage3list = stage3_from_gregions(/*stage3list*/(List_T) NULL,gregions,queryseq,queryuc,
#ifdef PMAP
				    queryntseq,
#endif
				    usersegment,stage2_alloc,oligoindices_major,oligoindices_minor,
				    pairpool,diagpool,cellpool,
				    dynprogL,dynprogM,dynprogR,worker_stopwatch);

  debug2(printf("Initial search gives stage3list of length %d\n",List_length(stage3list)));
#ifdef DEBUG2
  for (p = stage3list; p != NULL; p = List_next(p)) {
    Stage3_print_ends(List_head(p));
  }
#endif

  if (diag_debug == true) {
    return stage3list;		/* really diagonals */
  }

  queryntlength = Sequence_ntlength(queryseq);

  if (stage3list != NULL) {
    iter = 0;
    testlocalp = true;
    while (testlocalp == true && iter++ < MAX_CHIMERA_ITER) {
      debug2(printf("\n\n*** Testing for local on %d Stage3_T objects, iter %d ***\n",
		    List_length(stage3list),iter));

      /* Stage3_recompute_goodness(stage3list); */
      stage3list = stage3list_remove_duplicates(stage3list);
      stage3list = stage3list_sort(stage3list);

#ifdef DEBUG2
      for (p = stage3list; p != NULL; p = List_next(p)) {
	Stage3_print_ends(List_head(p));
      }
      printf("\n");
#endif
      nonchimericbest = (Stage3_T) List_head(stage3list);
      debug2(printf("nonchimericbest is %p\n",nonchimericbest));

#if 0
      if (List_length(stage3list) <= 1) {
	debug2(printf("Only 0 or 1 alignments, so won't look for local\n"));
	testlocalp = false;
      }
      else 
#endif

      if (Stage3_domain(nonchimericbest) < chimera_margin) {
	debug2(printf("Existing alignment is too short, so won't look for local\n"));
	testlocalp = false;

#if 0
      } else if (Stage3_fracidentity(nonchimericbest) < CHIMERA_IDENTITY &&
		 Chimera_alignment_break(&effective_start,&effective_end,nonchimericbest,Sequence_ntlength(queryseq),CHIMERA_FVALUE) >= chimera_margin
		 ) {
	debug2(printf("Break in alignment quality at %d..%d detected, so will look for local\n",
		      effective_start,effective_end));
	testlocalp = true;
#endif

      } else if (Stage3_largemargin(&effective_start,&effective_end,nonchimericbest,Sequence_ntlength(queryseq)) >= chimera_margin) {
	debug2(printf("Large margin at %d..%d detected (%d >= %d), so will look for local\n",
		      effective_start,effective_end,Stage3_largemargin(&effective_start,&effective_end,nonchimericbest,Sequence_ntlength(queryseq)),chimera_margin));
	testlocalp = true;
	
      } else {
	debug2(printf("Good alignment already with identity %f, so won't look for local\n",
		      Stage3_fracidentity(nonchimericbest)));
	testlocalp = false;
      }

      if (testlocalp == true) {
	testlocalp = false;
	debug2(printf("Checking for local, starting with list length %d, effective_start %d, effective_end %d\n",
		      List_length(stage3list),effective_start,effective_end));
	stage3list = check_for_local(&(*mergedp),stage3list,effective_start,effective_end,
				     queryseq,queryuc,
#ifdef PMAP
				     queryntseq,
#endif
				     queryntlength,usersegment,stage2_alloc,oligoindices_major,oligoindices_minor,
				     matchpool,pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR);
	
	if (*mergedp == true) {
	  testlocalp = true;	/* Local merge */
	} else {
	  debug2(printf("Checking for middle piece local, starting with list length %d\n",List_length(stage3list)));
	  stage3list = check_middle_piece_local(&foundp,stage3list,queryseq,queryuc,
#ifdef PMAP
						queryntseq,
#endif
						queryntlength,stage2_alloc,oligoindices_major,oligoindices_minor,
						pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR);
	  if (foundp == true) {
	    /* Iterate */
	    testlocalp = true;
	  }
	}
      }
    }
  }

  if (stage3list != NULL) {
    iter = 0;
    testchimerap = true;
    while (testchimerap == true && iter++ < MAX_CHIMERA_ITER) {
      debug2(printf("\n\n*** Testing for chimera on %d Stage3_T objects, iter %d ***\n",
		    List_length(stage3list),iter));

      /* Stage3_recompute_goodness(stage3list); */
      stage3list = stage3list_remove_duplicates(stage3list);
      stage3list = stage3list_sort(stage3list);

#ifdef DEBUG2
      for (p = stage3list; p != NULL; p = List_next(p)) {
	Stage3_print_ends(List_head(p));
      }
      printf("\n");
#endif
      nonchimericbest = (Stage3_T) List_head(stage3list);
      debug2(printf("nonchimericbest is %p\n",nonchimericbest));

      if (novelsplicingp == false) {
	testchimerap = false;

      } else if (chimera_margin <= 0) {
	debug2(printf("turned off\n"));
	testchimerap = false;

      } else if (maxpaths_report == 1) {
	debug2(printf("maxpaths set to 1\n"));
	testchimerap = false;

      } else if (Stage3_domain(nonchimericbest) < chimera_margin) {
	debug2(printf("Existing alignment is too short, so won't look for chimera\n"));
	testchimerap = false;

#if 0
      } else if (Stage3_fracidentity(nonchimericbest) < CHIMERA_IDENTITY &&
		 Chimera_alignment_break(&effective_start,&effective_end,nonchimericbest,Sequence_ntlength(queryseq),CHIMERA_FVALUE) >= chimera_margin
		 ) {
	debug2(printf("Break in alignment quality at %d..%d detected, so will look for chimera\n",
		      effective_start,effective_end));
	testchimerap = true;
#endif

      } else if (Stage3_largemargin(&effective_start,&effective_end,nonchimericbest,Sequence_ntlength(queryseq)) >= chimera_margin) {
	debug2(printf("Large margin at %d..%d detected (%d >= %d), so will look for chimera\n",
		      effective_start,effective_end,Stage3_largemargin(&effective_start,&effective_end,nonchimericbest,Sequence_ntlength(queryseq)),chimera_margin));
	testchimerap = true;
	
      } else {
	debug2(printf("Good alignment already with identity %f, so won't look for chimera\n",
		      Stage3_fracidentity(nonchimericbest)));
	testchimerap = false;
      }

      if (testchimerap == true) {
	testchimerap = false;
	debug2(printf("Checking for chimera, starting with list length %d, effective_start %d, effective_end %d\n",
		      List_length(stage3list),effective_start,effective_end));
	stage3list = check_for_chimera(&(*mergedp),&(*chimera),stage3list,effective_start,effective_end,
				       queryseq,queryuc,
#ifdef PMAP
				       queryntseq,
#endif
				       queryntlength,usersegment,stage2_alloc,oligoindices_major,oligoindices_minor,
				       matchpool,pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR);
	debug2(printf("chimera is %p\n",*chimera));
	if (*chimera != NULL) {
	  testchimerap = false;
	} else {
	  if (*mergedp == true) {
	    testchimerap = true;	/* Local merge */
	  } else {
	    debug2(printf("Checking for middle piece chimera, starting with list length %d\n",List_length(stage3list)));
	    stage3list = check_middle_piece_chimera(&foundp,stage3list,queryseq,queryuc,
#ifdef PMAP
						    queryntseq,
#endif
						    queryntlength,usersegment,stage2_alloc,oligoindices_major,oligoindices_minor,
						    matchpool,pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR);
	    if (foundp == true) {
	      /* Iterate */
	      testchimerap = true;
	    } else {
	      testchimerap = false;
	    }
	  }
	}
	debug2(printf("testchimerap is %d\n",testchimerap));
      }
    }
  }

  debug2(printf("apply_stage3 returning list of length %d\n",List_length(stage3list)));

  /* Needed after call to stage3_from_gregions */
  /* Stage3_recompute_goodness(stage3list); */

  /* Final call, so do both filtering and sorting */
  Stage3_recompute_coverage(stage3list,queryseq);
  stage3list = stage3list_filter_and_sort(&(*chimera),stage3list);

  if (*chimera != NULL && List_length(stage3list) > 2) {
    /* Compare chimera against non-chimeric alignments */
    chimera1 = (Stage3_T) List_head(stage3list);
    chimera2 = (Stage3_T) List_head(List_next(stage3list));
    nonchimericbest = (Stage3_T) List_head(List_next(List_next(stage3list)));
    debug2(printf("chimera1 %d, chimera2 %d\n",Stage3_goodness(chimera1),Stage3_goodness(chimera2)));
    debug2(printf("%p non-chimeric %d %d..%d\n",
		  nonchimericbest,Stage3_goodness(nonchimericbest),Stage3_querystart(nonchimericbest),Stage3_queryend(nonchimericbest)));

    if (Stage3_queryend(nonchimericbest) > (Stage3_querystart(chimera2) + Stage3_queryend(chimera2))/2 &&
	Stage3_querystart(nonchimericbest) < (Stage3_querystart(chimera1) + Stage3_queryend(chimera1))/2) {
      stage3list = List_pop(stage3list,(void **) &chimera1);
      stage3list = List_pop(stage3list,(void **) &chimera2);
      Stage3_free(&chimera1);
      Stage3_free(&chimera2);
      Chimera_free(&(*chimera));
      *chimera = (Chimera_T) NULL;
    }
  }

  return stage3list;
}


static Filestring_T
process_request (Filestring_T *fp_failedinput, double *worker_runtime, Request_T request, Sequence_T usersegment,
		 Matchpool_T matchpool, Pairpool_T pairpool, Diagpool_T diagpool, Cellpool_T cellpool,
		 Stage2_alloc_T stage2_alloc, Oligoindex_array_T oligoindices_major, Oligoindex_array_T oligoindices_minor,
		 Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		 Stopwatch_T worker_stopwatch) {
  Filestring_T fp;
  Result_T result;
  int jobid;
  Diagnostic_T diagnostic;
  Sequence_T queryseq, queryuc;
  Chimera_T chimera = NULL;
  bool mergedp, lowidentityp;
  bool repetitivep = false, poorp = false;

  List_T gregions = NULL, stage3list;
  Stage3_T *stage3array;
  int npaths, first_absmq, second_absmq;
#ifdef PMAP
  Sequence_T queryntseq;
#endif

  jobid = Request_id(request);
  queryseq = Request_queryseq(request);
  Matchpool_reset(matchpool);
  Pairpool_reset(pairpool);
  Diagpool_reset(diagpool);
  Cellpool_reset(cellpool);


  if (worker_stopwatch != NULL) {
    Stopwatch_start(worker_stopwatch);
  }

  if (Sequence_fulllength_given(queryseq) <= 0) {
    result = Result_new(jobid,/*mergedp*/false,(Chimera_T) NULL,(Stage3_T *) NULL,
			/*npaths*/0,/*first_absmq*/0,/*second_absmq*/0,/*diagnostic*/NULL,EMPTY_SEQUENCE);
      
  } else if (Sequence_fulllength_given(queryseq) < 
#ifdef PMAP
	     index1part_aa
#else
	     index1part
#endif
	     ) {
    result = Result_new(jobid,/*mergedp*/false,(Chimera_T) NULL,(Stage3_T *) NULL,
			/*npaths*/0,/*first_absmq*/0,/*second_absmq*/0,/*diagnostic*/NULL,SHORT_SEQUENCE);

  } else {			/* Sequence_fulllength_given(queryseq) > 0 */
    queryuc = Sequence_uppercase(queryseq);
#ifdef PMAP
    queryntseq = Sequence_convert_to_nucleotides(queryseq);
#endif

    diagnostic = evaluate_query(&poorp,&repetitivep,Sequence_fullpointer(queryuc),
				Sequence_fulllength(queryuc),Oligoindex_array_elt(oligoindices_major,0));

#ifndef PMAP
    if (poorp == true && prune_poor_p == true) {
      result = Result_new(jobid,/*mergedp*/false,(Chimera_T) NULL,(Stage3_T *) NULL,
			  /*npaths*/0,/*first_absmq*/0,/*second_absmq*/0,diagnostic,POOR_SEQUENCE);
    } else if (repetitivep == true && prune_repetitive_p == true) {
      result = Result_new(jobid,/*mergedp*/false,(Chimera_T) NULL,(Stage3_T *) NULL,
			  /*npaths*/0,/*first_absmq*/0,/*second_absmq*/0,diagnostic,REPETITIVE);
    }
#endif

    if (usersegment != NULL) {
#ifndef PMAP
#if 0
      /* Don't do Sequence_trim, because it affects sequences like NM_018406 */
      Sequence_trim(queryseq,diagnostic->query_trim_start,diagnostic->query_trim_end);
      Sequence_trim(queryuc,diagnostic->query_trim_start,diagnostic->query_trim_end);
#endif
#endif
      stage3array = stage3_from_usersegment(&npaths,&first_absmq,&second_absmq,queryseq,queryuc,
#ifdef PMAP
					    queryntseq,
#endif
					    usersegment,stage2_alloc,oligoindices_major,oligoindices_minor,
					    pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,worker_stopwatch);
      result = Result_new(jobid,/*mergedp*/false,(Chimera_T) NULL,stage3array,npaths,first_absmq,second_absmq,diagnostic,NO_FAILURE);

    } else {		/* Not user segment and not maponly */
#ifndef PMAP
#if 0
      /* Don't do Sequence_trim, because it affects sequences like NM_018406 */
      Sequence_trim(queryseq,diagnostic->query_trim_start,diagnostic->query_trim_end);
      Sequence_trim(queryuc,diagnostic->query_trim_start,diagnostic->query_trim_end);
#endif
#endif

      debug(printf("Calling stage 1\n"));
      if (mode == CMET_NONSTRANDED) {
	gregions = Stage1_compute_nonstranded(&lowidentityp,queryuc,indexdb_fwd,indexdb_fwd,
					      /*indexdb_size_threshold*/400,chromosome_iit,
					      chrsubset_start,chrsubset_end,matchpool,
					      stutterhits,diagnostic,worker_stopwatch,/*nbest*/10);
      } else {
	gregions = Stage1_compute(&lowidentityp,queryuc,indexdb_fwd,indexdb_rev,
				  /*indexdb_size_threshold*/100,chromosome_iit,
				  chrsubset_start,chrsubset_end,matchpool,
				  stutterhits,diagnostic,worker_stopwatch,/*nbest*/10);
      }
      debug(printf("Got %d gregions\n",List_length(gregions)));

      if (stage1debug == true) {
	result = Result_new_stage1debug(jobid,gregions,diagnostic,NO_FAILURE);
      } else {
	debug(printf("Applying stage 3\n"));
	stage3list = apply_stage3(&mergedp,&chimera,gregions,queryseq,queryuc,
#ifdef PMAP
				  queryntseq,
#endif
				  usersegment,stage2_alloc,oligoindices_major,oligoindices_minor,
				  matchpool,pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,worker_stopwatch);
	if (diag_debug == true) {
	  result = Result_new_diag_debug(jobid,/*diagonals*/stage3list,diagnostic,NO_FAILURE);
	} else if (stage3list == NULL) {
	  result = Result_new(jobid,mergedp,chimera,/*stage3array*/NULL,/*npaths*/0,/*first_absmq*/0,/*second_absmq*/0,
			      diagnostic,NO_FAILURE);
	} else if (chimera == NULL) {
	  stage3array = stage3array_from_list(&npaths,&first_absmq,&second_absmq,stage3list,mergedp,
					      /*chimerap*/false,/*remove_overlaps_p*/true);
	  result = Result_new(jobid,mergedp,/*chimera*/NULL,stage3array,npaths,first_absmq,second_absmq,
			      diagnostic,NO_FAILURE);
	} else {
	  stage3array = stage3array_from_list(&npaths,&first_absmq,&second_absmq,stage3list,mergedp,
					      /*chimerap*/true,/*remove_overlaps_p*/false);
	  result = Result_new(jobid,mergedp,chimera,stage3array,npaths,first_absmq,second_absmq,
			      diagnostic,NO_FAILURE);
	}
      }

      Oligoindex_clear_inquery(Oligoindex_array_elt(oligoindices_major,0),/*queryuc_ptr*/Sequence_fullpointer(queryuc),
			       /*querystart*/0,/*queryend*/Sequence_fulllength(queryuc));

    } /* Matches not user segment and not maponly */

#ifdef PMAP
    Sequence_free(&queryntseq);
#endif
    Sequence_free(&queryuc);
  } /* Matches sequence length > 0 */

  fp = Output_filestring_fromresult(&(*fp_failedinput),result,request,
				    /*headerseq*/user_pairalign_p == true ? usersegment : queryseq);
  *worker_runtime = worker_stopwatch == NULL ? 0.00 : Stopwatch_stop(worker_stopwatch);
  Result_free(&result);
  return fp;
}


#ifdef HAVE_SIGACTION
static const Except_T sigfpe_error = {"SIGFPE--arithmetic exception"};
static const Except_T sigsegv_error = {"SIGSEGV--segmentation violation"};
static const Except_T sigtrap_error = {"SIGTRAP--hardware fault"};
static const Except_T misc_signal_error = {"Miscellaneous signal"};

static void
signal_handler (int sig) {
  Request_T request;
  Sequence_T queryseq;

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

#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

#ifdef HAVE_PTHREAD
  request = (Request_T) pthread_getspecific(global_request_key);
  if (request == NULL) {
    /* fprintf(stderr,"Unable to retrieve request for thread\n"); */
  } else {
    queryseq = Request_queryseq(request);
    if (queryseq == NULL) {
      fprintf(stderr,"Unable to retrieve queryseq for request\n");
    } else {
      fprintf(stderr,"Problem sequence: ");
      fprintf(stderr,"%s (%d bp)\n",Sequence_accession(queryseq),Sequence_fulllength(queryseq));
    }
  }
#endif

  exit(9);

  return;
}
#endif


#define POOL_FREE_INTERVAL 200

#ifdef USE_MPI
static void
worker_mpi_process (int worker_id, Inbuffer_T inbuffer) {
  bool donep = false;
  int nread = 0;
  MPI_Status status;

  Stage2_alloc_T stage2_alloc;
  Oligoindex_array_T oligoindices_major, oligoindices_minor;
  Dynprog_T dynprogL, dynprogM, dynprogR;
  Matchpool_T matchpool;
  Pairpool_T pairpool;
  Diagpool_T diagpool;
  Cellpool_T cellpool;
  Stopwatch_T worker_stopwatch;
  Request_T request;
  Filestring_T fp, fp_failedinput;
  Sequence_T queryseq, usersegment, pairalign_segment;
  int filestringid, requestid, i;
  int ret;
  int worker_jobid = 0;
  double worker_runtime;

#ifdef MEMUSAGE
  Sequence_T queryseq;
  long int memusage_constant = 0, memusage, max_memusage;
  char procname[12];
  char acc[100+1], comma0[20], comma1[20], comma2[20], comma3[20], comma4[20], comma5[20];
  sprintf(procname,"proc-%ld",worker_id);
  Mem_usage_set_threadname(procname);
#endif

  stage2_alloc = Stage2_alloc_new(MAX_QUERYLENGTH_FOR_ALLOC);
  oligoindices_major = Oligoindex_array_new_major(MAX_QUERYLENGTH_FOR_ALLOC,MAX_GENOMICLENGTH_FOR_ALLOC);
  oligoindices_minor = Oligoindex_array_new_minor(MAX_QUERYLENGTH_FOR_ALLOC,MAX_GENOMICLENGTH_FOR_ALLOC);
  dynprogL = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired,
			 /*doublep*/true);
  dynprogM = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired,
			 /*doublep*/false);
  dynprogR = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired,
			 /*doublep*/true);
  matchpool = Matchpool_new();
  pairpool = Pairpool_new();
  diagpool = Diagpool_new();
  cellpool = Cellpool_new();
  worker_stopwatch = (timingp == true) ? Stopwatch_new() : (Stopwatch_T) NULL;

  usersegment = global_usersegment;

  /* Except_stack_create(); -- no worker threads, so no need to store request in global_request_key */

#ifdef MEMUSAGE
  memusage_constant += Mem_usage_report_std_heap();
  Genomicpos_commafmt_fill(comma0,memusage_constant);
  Mem_usage_reset_heap_baseline(0);
#endif

  /* Initial message to say that we are ready for a request */
  filestringid = -1;

  /* Use a synchronized send here to make sure outbuffer is ready */
  if ((ret = MPI_SSEND(&filestringid,1,MPI_INT,/*dest*/0,/*tag*/MPI_TAG_FILESTRING_AVAIL,MPI_COMM_WORLD)) != 0) {
    fprintf(stderr,"MPI_SSEND returns error %d\n",ret);
    MPI_Finalize();
    exit(9);
  }

  while (donep == false) {
    MPI_RECV(&requestid,1,MPI_INT,/*source*/0,/*tag*/MPI_ANY_TAG,MPI_COMM_WORLD,&status);
    debugm(printf("worker_id %ld got request %d\n",worker_id,requestid));

    while (nread < requestid &&
	   (queryseq = Inbuffer_read(&pairalign_segment,inbuffer,/*skipp*/true)) != NULL) {
      /* No need to free queryseq */
      nread++;
    }

    if (nread < requestid) {
      debugm(printf("because nread %d < requestid %d, worker_id %ld is done\n",nread,requestid,worker_id));
      donep = true;
    } else if ((queryseq = Inbuffer_read(&pairalign_segment,inbuffer,/*skipp*/false)) == NULL) {
      debugm(printf("because final read is NULL, worker_id %ld is done\n",worker_id));
      donep = true;
    } else {
      debugm(printf("worker_id %ld starting to process request %d\n",worker_id,requestid));
      request = Request_new(requestid,queryseq);
      nread++;

      if (user_pairalign_p == true) {
	genomecomp_blocks = Compress_create_blocks_comp(Sequence_fullpointer(usersegment),Sequence_fulllength(usersegment));
	genomebits_blocks = Compress_create_blocks_bits(genomecomp_blocks,Sequence_fulllength(usersegment));
	Genome_user_setup(genomecomp_blocks);
	Genome_hr_user_setup(genomebits_blocks,/*query_unk_mismatch_p*/false,
			     /*genome_unk_mismatch_p*/true,/*mode*/STANDARD);
	Genome_sites_setup(genomecomp_blocks,/*snp_blocks*/NULL);
	Maxent_hr_setup(genomecomp_blocks,/*genomealt_blocks*/genomecomp_blocks);
#ifdef PMAP
	Oligoindex_pmap_setup(genomecomp);
#else
	Oligoindex_hr_setup(genomecomp_blocks,mode);
#endif
	usersegment = pairalign_segment;
      }

#ifdef MEMUSAGE
      queryseq = Request_queryseq(request);
      fprintf(stderr,"Proc %d starting %s\n",worker_id,Sequence_accession(queryseq));
      Mem_usage_reset_stack_max();
      Mem_usage_reset_heap_max();
#endif

      TRY
        fp = process_request(&fp_failedinput,&worker_runtime,request,usersegment,
			     matchpool,pairpool,diagpool,cellpool,
			     stage2_alloc,oligoindices_major,oligoindices_minor,
			     dynprogL,dynprogM,dynprogR,worker_stopwatch);

      ELSE
	queryseq = Request_queryseq(request);
        if (Sequence_accession(queryseq) == NULL) {
	  fprintf(stderr,"Problem with unnamed sequence (%d bp)\n",Sequence_fulllength_given(queryseq));
	} else {
	  fprintf(stderr,"Problem with sequence %s (%d bp)\n",
		  Sequence_accession(queryseq),Sequence_fulllength_given(queryseq));
	}
	fprintf(stderr,"To obtain a core dump, re-run program on problem sequence with the -0 [zero] flag\n");
	fprintf(stderr,"Exiting...\n");
	exit(9);
      RERAISE;
      END_TRY;

      if (user_pairalign_p == true) {
	FREE(genomebits_blocks);
	FREE(genomecomp_blocks);
      }

      filestringid = Filestring_id(fp);
      debugm(printf("worker proc %d sending filestring %d...",worker_id,filestringid));

      /* Use a synchronized send here to make sure outbuffer is ready */
      if ((ret = MPI_SSEND(&filestringid,1,MPI_INT,/*dest*/0,/*tag*/MPI_TAG_FILESTRING_AVAIL,MPI_COMM_WORLD)) != 0) {
	fprintf(stderr,"MPI_SSEND returns error %d\n",ret);
	MPI_Finalize();
	exit(9);
      }
      Filestring_Send(fp,/*dest*/0,/*tag*/MPI_TAG_DEFAULT,MPI_COMM_WORLD);
      if (failedinput_root != NULL) {
	Filestring_Send(fp_failedinput,/*dest*/0,/*tag*/MPI_TAG_DEFAULT,MPI_COMM_WORLD);
      }
      debugm(printf("done with filestring %d\n",filestringid));

      if (worker_jobid % POOL_FREE_INTERVAL == 0) {
	Pairpool_free_memory(pairpool);
	Diagpool_free_memory(diagpool);
	Cellpool_free_memory(cellpool);
	Matchpool_free_memory(matchpool);
      }

#ifdef MEMUSAGE
      /* Copy acc before we free the request */
      queryseq = Request_queryseq(request);
      strncpy(acc,Sequence_accession(queryseq),100);
      acc[100] = '\0';
#endif

      Request_free(&request);

#ifdef MEMUSAGE
      Genomicpos_commafmt_fill(comma1,Mem_usage_report_std_heap_max());
      Genomicpos_commafmt_fill(comma2,Mem_usage_report_std_heap());
      Genomicpos_commafmt_fill(comma3,Mem_usage_report_keep());
      Genomicpos_commafmt_fill(comma4,Mem_usage_report_in());
      Genomicpos_commafmt_fill(comma5,Mem_usage_report_out());

      fprintf(stderr,"Acc %s, proc %d: constant %s  max %s  std %s  keep %s  in %s  out %s\n",
	      acc,worker_id,comma0,comma1,comma2,comma3,comma4,comma5);

      if ((memusage = Mem_usage_report_std_heap()) != 0) {
	fprintf(stderr,"Memory leak in proc of %ld bytes: %ld\n",worker_id,memusage);
	fflush(stdout);
	MPI_Finalize();
	exit(9);
      }
#endif
    }
  }

  /* Final message to say that we are done with all requests */
  debugm(printf("worker_id %ld sending final message to say it is done\n",worker_id));
  filestringid = -1;
  if ((ret = MPI_SSEND(&filestringid,1,MPI_INT,/*dest*/0,/*tag*/MPI_TAG_FILESTRING_AVAIL,MPI_COMM_WORLD)) != 0) {
    fprintf(stderr,"MPI_SSEND returns error %d\n",ret);
    MPI_Finalize();
    exit(9);
  }

#ifdef MEMUSAGE
  Mem_usage_std_heap_add(memusage_constant);
#endif

  /* Except_stack_destroy(); */

  Stopwatch_free(&worker_stopwatch);
  Cellpool_free(&cellpool);
  Diagpool_free(&diagpool);
  Pairpool_free(&pairpool);
  Matchpool_free(&matchpool);
  Dynprog_free(&dynprogR);
  Dynprog_free(&dynprogM);
  Dynprog_free(&dynprogL);
  Oligoindex_array_free(&oligoindices_minor);
  Oligoindex_array_free(&oligoindices_major);
  Stage2_alloc_free(&stage2_alloc);

#ifdef MEMUSAGE
  Mem_usage_set_threadname("main");
#endif

  debugm(printf("worker_id %ld is now returning\n",worker_id));
  return;
}
#endif


static void
single_thread () {
  Stage2_alloc_T stage2_alloc;
  Oligoindex_array_T oligoindices_major, oligoindices_minor;
  Dynprog_T dynprogL, dynprogM, dynprogR;
  Matchpool_T matchpool;
  Pairpool_T pairpool;
  Diagpool_T diagpool;
  Cellpool_T cellpool;
  Stopwatch_T worker_stopwatch;
  Request_T request;
  Sequence_T usersegment, pairalign_segment;
  Filestring_T fp, fp_failedinput;
  Sequence_T queryseq;
  int noutput = 0;
  int jobid = 0;
  double worker_runtime;

#ifdef MEMUSAGE
  long int memusage, memusage_constant = 0;
  char acc[100+1], comma0[20], comma1[20], comma2[20], comma3[20], comma4[20], comma5[20];
#endif

  stage2_alloc = Stage2_alloc_new(MAX_QUERYLENGTH_FOR_ALLOC);
  oligoindices_major = Oligoindex_array_new_major(MAX_QUERYLENGTH_FOR_ALLOC,MAX_GENOMICLENGTH_FOR_ALLOC);
  oligoindices_minor = Oligoindex_array_new_minor(MAX_QUERYLENGTH_FOR_ALLOC,MAX_GENOMICLENGTH_FOR_ALLOC);
  dynprogL = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired,
			 /*doublep*/true);
  dynprogM = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired,
			 /*doublep*/false);
  dynprogR = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired,
			 /*doublep*/true);
  matchpool = Matchpool_new();
  pairpool = Pairpool_new();
  diagpool = Diagpool_new();
  cellpool = Cellpool_new();
  worker_stopwatch = (timingp == true) ? Stopwatch_new() : (Stopwatch_T) NULL;

  usersegment = global_usersegment;

  /* Except_stack_create(); -- requires pthreads */

#ifdef MEMUSAGE
  memusage_constant += Mem_usage_report_std_heap();
  Genomicpos_commafmt_fill(comma0,memusage_constant);
  Mem_usage_reset_heap_baseline(0);
#endif

  while ((request = Inbuffer_get_request(&pairalign_segment,inbuffer)) != NULL) {

    if (user_pairalign_p == true) {
      genomecomp_blocks = Compress_create_blocks_comp(Sequence_fullpointer(usersegment),Sequence_fulllength(usersegment));
      genomebits_blocks = Compress_create_blocks_bits(genomecomp_blocks,Sequence_fulllength(usersegment));
      Genome_user_setup(genomecomp_blocks);
      Genome_hr_user_setup(genomebits_blocks,/*query_unk_mismatch_p*/false,
			   /*genome_unk_mismatch_p*/true,/*mode*/STANDARD);
      Genome_sites_setup(genomecomp_blocks,/*snp_blocks*/NULL);
      Maxent_hr_setup(genomecomp_blocks,/*genomealt_blocks*/genomecomp_blocks);
#ifdef PMAP
      Oligoindex_pmap_setup(genomecomp);
#else
      Oligoindex_hr_setup(genomecomp_blocks,mode);
#endif
      usersegment = pairalign_segment;
    }

#ifdef MEMUSAGE
    queryseq = Request_queryseq(request);
    fprintf(stderr,"Single thread starting %s\n",Sequence_accession(queryseq));
    Mem_usage_reset_stack_max();
    Mem_usage_reset_heap_max();
#endif

    TRY
      fp = process_request(&fp_failedinput,&worker_runtime,request,usersegment,
			   matchpool,pairpool,diagpool,cellpool,
			   stage2_alloc,oligoindices_major,oligoindices_minor,
			   dynprogL,dynprogM,dynprogR,worker_stopwatch);
      if (timingp == true) {
        queryseq = Request_queryseq(request);
        printf("%s\t%.6f\n",Sequence_accession(queryseq),worker_runtime);
      }

    ELSE
      queryseq = Request_queryseq(request);
      if (Sequence_accession(queryseq) == NULL) {
        fprintf(stderr,"Problem with unnamed sequence (%d bp)\n",Sequence_fulllength_given(queryseq));
      } else {
        fprintf(stderr,"Problem with sequence %s (%d bp)\n",
  	      Sequence_accession(queryseq),Sequence_fulllength_given(queryseq));
      }
      fprintf(stderr,"To obtain a core dump, re-run program on problem sequence with the -0 [zero] flag\n");
      fprintf(stderr,"Exiting...\n");
      exit(9);
    RERAISE;
    END_TRY;

    if (user_pairalign_p == true) {
      FREE(genomebits_blocks);
      FREE(genomecomp_blocks);
    }

    Outbuffer_print_filestrings(fp,fp_failedinput);

    if (jobid % POOL_FREE_INTERVAL == 0) {
      Pairpool_free_memory(pairpool);
      Diagpool_free_memory(diagpool);
      Cellpool_free_memory(cellpool);
      Matchpool_free_memory(matchpool);
    }

#ifdef MEMUSAGE
    /* Copy acc before we free the request */
    queryseq = Request_queryseq(request);
    strncpy(acc,Sequence_accession(queryseq),100);
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

  if (worker_stopwatch != NULL) {
    Stopwatch_free(&worker_stopwatch);
  }
  Cellpool_free(&cellpool);
  Diagpool_free(&diagpool);
  Pairpool_free(&pairpool);
  Matchpool_free(&matchpool);
  Dynprog_free(&dynprogR);
  Dynprog_free(&dynprogM);
  Dynprog_free(&dynprogL);
  Oligoindex_array_free(&oligoindices_minor);
  Oligoindex_array_free(&oligoindices_major);
  Stage2_alloc_free(&stage2_alloc);

#ifdef MEMUSAGE
  Mem_usage_set_threadname("main");
#endif

  return;
}


#ifdef HAVE_PTHREAD
static void *
worker_thread (void *data) {
  Stage2_alloc_T stage2_alloc;
  Oligoindex_array_T oligoindices_major, oligoindices_minor;
  Dynprog_T dynprogL, dynprogM, dynprogR;
  Matchpool_T matchpool;
  Pairpool_T pairpool;
  Diagpool_T diagpool;
  Cellpool_T cellpool;
  Stopwatch_T worker_stopwatch;
  Request_T request;
  Filestring_T fp, fp_failedinput;
  Sequence_T queryseq, usersegment, pairalign_segment;
  int worker_jobid = 0;
  double worker_runtime;
#if defined(DEBUG) || defined(MEMUSAGE)
  long int worker_id = (long int) data;
#endif

#ifdef MEMUSAGE
  long int memusage_constant = 0, memusage, max_memusage;
  char threadname[12];
  char acc[100+1], comma0[20], comma1[20], comma2[20], comma3[20], comma4[20], comma5[20];
  sprintf(threadname,"thread-%ld",worker_id);
  Mem_usage_set_threadname(threadname);
#endif

  /* Thread-specific data and storage */
  stage2_alloc = Stage2_alloc_new(MAX_QUERYLENGTH_FOR_ALLOC);
  oligoindices_major = Oligoindex_array_new_major(MAX_QUERYLENGTH_FOR_ALLOC,MAX_GENOMICLENGTH_FOR_ALLOC);
  oligoindices_minor = Oligoindex_array_new_minor(MAX_QUERYLENGTH_FOR_ALLOC,MAX_GENOMICLENGTH_FOR_ALLOC);
  dynprogL = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired,
			 /*doublep*/true);
  dynprogM = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired,
			 /*doublep*/false);
  dynprogR = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired,
			 /*doublep*/true);
  matchpool = Matchpool_new();
  pairpool = Pairpool_new();
  diagpool = Diagpool_new();
  cellpool = Cellpool_new();
  worker_stopwatch = (timingp == true) ? Stopwatch_new() : (Stopwatch_T) NULL;

  usersegment = global_usersegment;

  Except_stack_create();

#ifdef MEMUSAGE
  memusage_constant += Mem_usage_report_std_heap();
  Genomicpos_commafmt_fill(comma0,memusage_constant);
  Mem_usage_reset_heap_baseline(0);
#endif

  while ((request = Inbuffer_get_request(&pairalign_segment,inbuffer)) != NULL) {
    debug(printf("worker_thread %ld got request %d\n",worker_id,Request_id(request)));
    pthread_setspecific(global_request_key,(void *) request);

    if (user_pairalign_p == true) {
      genomecomp_blocks = Compress_create_blocks_comp(Sequence_fullpointer(usersegment),Sequence_fulllength(usersegment));
      genomebits_blocks = Compress_create_blocks_bits(genomecomp_blocks,Sequence_fulllength(usersegment));
      Genome_user_setup(genomecomp_blocks);
      Genome_hr_user_setup(genomebits_blocks,/*query_unk_mismatch_p*/false,
			   /*genome_unk_mismatch_p*/true,/*mode*/STANDARD);
      Genome_sites_setup(genomecomp_blocks,/*snp_blocks*/NULL);
      Maxent_hr_setup(genomecomp_blocks,/*genomealt_blocks*/genomecomp_blocks);
#ifdef PMAP
      Oligoindex_pmap_setup(genomecomp);
#else
      Oligoindex_hr_setup(genomecomp_blocks,mode);
#endif
      usersegment = pairalign_segment;
    }

#ifdef MEMUSAGE
    queryseq = Request_queryseq(request);
    fprintf(stderr,"Thread %d starting %s\n",worker_id,Sequence_accession(queryseq));
    Mem_usage_reset_stack_max();
    Mem_usage_reset_heap_max();
#endif

    TRY
      fp = process_request(&fp_failedinput,&worker_runtime,request,usersegment,
			   matchpool,pairpool,diagpool,cellpool,
			   stage2_alloc,oligoindices_major,oligoindices_minor,
			   dynprogL,dynprogM,dynprogR,worker_stopwatch);
      if (timingp == true) {
        queryseq = Request_queryseq(request);
        printf("%s\t%.6f\n",Sequence_accession(queryseq),worker_runtime);
      }

    ELSE
      queryseq = Request_queryseq(request);
      if (queryseq == NULL) {
	fprintf(stderr,"NULL");
      } else if (Sequence_accession(queryseq) == NULL) {
	fprintf(stderr,"unnamed (%d bp)",Sequence_fulllength_given(queryseq));
      } else {
	fprintf(stderr,"%s (%d bp)",Sequence_accession(queryseq),Sequence_fulllength_given(queryseq));
      }
      fprintf(stderr,"\n");
      fprintf(stderr,"To obtain a core dump, re-run program on problem sequence with the -0 [zero] flag\n");

      fprintf(stderr,"Exiting...\n");
      exit(9);
    RERAISE;
    END_TRY;

    if (user_pairalign_p == true) {
      FREE(genomebits_blocks);
      FREE(genomecomp_blocks);
    }

    debug(printf("worker_thread %ld putting filestring %d\n",worker_id,Filestring_id(fp)));
    Outbuffer_put_filestrings(outbuffer,fp,fp_failedinput);

    if (worker_jobid % POOL_FREE_INTERVAL == 0) {
      Pairpool_free_memory(pairpool);
      Diagpool_free_memory(diagpool);
      Cellpool_free_memory(cellpool);
      Matchpool_free_memory(matchpool);
    }

#ifdef MEMUSAGE
    /* Copy acc before we free the request */
    queryseq = Request_queryseq(request);
    strncpy(acc,Sequence_accession(queryseq),100);
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

  if (worker_stopwatch != NULL) {
    Stopwatch_free(&worker_stopwatch);
  }
  Cellpool_free(&cellpool);
  Diagpool_free(&diagpool);
  Pairpool_free(&pairpool);
  Matchpool_free(&matchpool);
  Dynprog_free(&dynprogR);
  Dynprog_free(&dynprogM);
  Dynprog_free(&dynprogL);
  Oligoindex_array_free(&oligoindices_minor);
  Oligoindex_array_free(&oligoindices_major);
  Stage2_alloc_free(&stage2_alloc);

#ifdef MEMUSAGE
  Mem_usage_set_threadname("main");
#endif

  return (void *) NULL;
}
#endif


#if 0

static void
align_relative (FILE *input, char **files, int nfiles, int nextchar,
		Sequence_T queryseq, Sequence_T referenceseq) {
  Stage2_alloc_T stage2_alloc;
  Oligoindex_array_T oligoindices_major, oligoindices_minor;
  Diagnostic_T diagnostic;
  bool lowidentityp;
#ifndef PMAP
  bool poorp, repetitivep;
#endif
  Dynprog_T dynprogL, dynprogM, dynprogR;
  Matchpool_T matchpool;
  Pairpool_T pairpool;
  Diagpool_T diagpool;
  Cellpool_T cellpool;
  Stopwatch_T stopwatch;

  Chrpos_T genomicstart, genomiclength;
  Sequence_T genomicseg, queryuc, referenceuc;
  int jobid = 0;

  Chimera_T chimera = NULL;
  List_T gregions, stage3list;
  Stage3_T *stage3array, stage3, stage3ref;
  int npaths, i;

  oligoindices_major = Oligoindex_array_new_major(&noligoindices_major);
  oligoindices_minor = Oligoindex_array_new_minor(&noligoindices_minor);
  dynprogL = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired);
  dynprogM = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired);
  dynprogR = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired);
  matchpool = Matchpool_new();
  pairpool = Pairpool_new();
  diagpool = Diagpool_new();
  cellpool = Cellpool_new();
  stopwatch = (timingp == true) ? Stopwatch_new() : (Stopwatch_T) NULL;

  Matchpool_reset(matchpool);
  Pairpool_reset(pairpool);
  Diagpool_reset(diagpool);
  Cellpool_reset(cellpool);

  referenceuc = Sequence_uppercase(referenceseq);

  /* Do not trim the mutation refseq */
  diagnostic = Diagnostic_new();
  Oligoindex_set_inquery(&diagnostic->query_badoligos,&diagnostic->query_repoligos,&diagnostic->query_trimoligos,
			 &diagnostic->query_trim_start,&diagnostic->query_trim_end,Oligoindex_array_elt(oligoindices_major,0),
			 Sequence_fullpointer(referenceuc),Sequence_fulllength(referenceuc),/*trimp*/false);
#ifndef PMAP
#if 0
  /* Don't do Sequence_trim, because it affects sequences like NM_018406 */
  Sequence_trim(referenceseq,diagnostic->query_trim_start,diagnostic->query_trim_end);
#endif
#endif
  gregions = Stage1_compute(&lowidentityp,referenceuc,indexdb_fwd,indexdb_rev,
			    /*indexdb_size_threshold*/100,chromosome_iit,
			    chrsubset_start,chrsubet_end,matchpool,
			    stutterhits,diagnostic,/*stopwatch*/NULL);
  stage3list = apply_stage3(&chimera,gregions,referenceseq,referenceuc,/*usersegment*/NULL,
			    oligoindices_major,oligoindices_minor,
			    matchpool,pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,stopwatch);
  if (stage3list == NULL) {
    npaths = 0;
    stage3array = (Stage3_T *) NULL;
  } else {
    stage3array = stage3array_from_list(&npaths,stage3list,/*mergedp*/false,/*chimerap*/false,/*remove_overlaps_p*/true);
  }

  Diagnostic_free(&diagnostic);

  /* chimera should be NULL */
  for (i = 1; i < npaths; i++) {
    stage3 = stage3array[i];
    Stage3_free(&stage3);
  }
  if (npaths > 0) {
    stage3ref = stage3array[0];
#ifdef PMAP
    Stage3_translate_cdna(stage3ref,queryseq,strictp);
    Stage3_backtranslate_cdna(stage3ref,/*diagnosticp*/false);
#else
    Stage3_translate_genomic(stage3ref,/*fulllengthp*/true,/*cds_startpos*/-1,
			     Sequence_fulllength_given(queryseq),/*truncatep*/false,strictp);
#endif
    FREE(stage3array);

    Stage3_genomicbounds(&genomicstart,&genomiclength,stage3ref);
    if (genomealt != NULL) {
      genomicseg = Genome_get_segment(genomealt,genomicstart,genomiclength,chromosome_iit,/*revcomp*/false);
    } else {
      genomicseg = Genome_get_segment(genome,genomicstart,genomiclength,chromosome_iit,/*revcomp*/false);
    }

    while (jobid == 0 || (queryseq = Sequence_read_multifile(&nextchar,&input,&files,&nfiles)) != NULL) {
      Matchpool_reset(matchpool);
      Pairpool_reset(pairpool);
      Diagpool_reset(diagpool);
      Cellpool_reset(cellpool);

      fprintf(fp,">");
      Sequence_print_header(stdout,queryseq,checksump);
      diagnostic = Diagnostic_new();
      if (Sequence_fulllength_given(queryseq) <= 0) {
	print_npaths(fp,0,diagnostic,/*usersegment*/NULL,chrsubset,/*chimera*/NULL,EMPTY_SEQUENCE);

      } else if (Sequence_fulllength_given(queryseq) <
#ifdef PMAP
		 index1part_aa
#else
		 index1part
#endif
		 ) {
	print_npaths(fp,0,diagnostic,/*usersegment*/NULL,chrsubset,/*chimera*/NULL,SHORT_SEQUENCE);

      } else {

	queryuc = Sequence_uppercase(queryseq);
#ifdef PMAP
	Oligoindex_set_inquery(&diagnostic->query_badoligos,&diagnostic->query_repoligos,
			       &diagnostic->query_trimoligos,&diagnostic->query_trim_start,
			       &diagnostic->query_trim_end,Oligoindex_array_elt(oligoindices_major,0),
			       Sequence_fullpointer(queryuc),Sequence_fulllength(queryuc),/*trimp*/false);
#else
	diagnostic->query_oligodepth = 
	  Oligoindex_set_inquery(&diagnostic->query_badoligos,&diagnostic->query_repoligos,
				 &diagnostic->query_trimoligos,&diagnostic->query_trim_start,
				 &diagnostic->query_trim_end,Oligoindex_array_elt(oligoindices_major,0),
				 Sequence_fullpointer(queryuc),/*querystart*/0,/*queryend*/Sequence_fulllength(queryuc),
				 /*trimp*/true);

	if (diagnostic->query_trimoligos == 0) {
	  poorp = true;
	} else if (((double) diagnostic->query_badoligos/(double) diagnostic->query_trimoligos > MAX_BADOLIGOS) ||
		   (diagnostic->query_trim_end - diagnostic->query_trim_start < 80 && diagnostic->query_badoligos > 0)) {
	  poorp = true;
	} else {
	  poorp = false;
	}
#if 0
	if (diagnostic->query_trimoligos == 0) {
	  repetitivep = false;
	} else if (diagnostic->query_oligodepth > MAX_OLIGODEPTH ||
		   (double) diagnostic->query_repoligos/(double) diagnostic->query_trimoligos > MAX_REPOLIGOS) {
	  repetitivep = true;
	} else {
	  repetitivep = false;
	}
#endif
	repetitivep = false;

	if (poorp == true && prune_poor_p == true) {
	  print_npaths(fp,0,diagnostic,/*usersegment*/NULL,chrsubset,/*chimera*/NULL,POOR_SEQUENCE);
	} else if (repetitivep == true && prune_repetitive_p == true) {
	  print_npaths(fp,0,diagnostic,/*usersegment*/NULL,chrsubset,/*chimera*/NULL,REPETITIVE);
	} else {
#endif /* PMAP */
	  stage3array = stage3_from_usersegment(&npaths,queryseq,queryuc,genomicseg,
						oligoindices_major,oligoindices_minor,
						pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,stopwatch);
	  
	  if (npaths == 0) {
	    print_npaths(fp,0,diagnostic,/*usersegment*/NULL,chrsubset,/*chimera*/NULL,NO_FAILURE);
	  } else if (printtype == COORDS) {
	    Stage3_print_coordinates(fp,stage3array[0],chromosome_iit,invertmode);

	  } else {
	    /* Usual output */
	    print_npaths(fp,1,diagnostic,/*usersegment*/NULL,chrsubset,/*chimera*/NULL,NO_FAILURE);
#ifndef PMAP
	    Stage3_translate_cdna_via_reference(stage3array[0],stage3ref,literalrefp);
#endif
	    Stage3_fix_cdna_direction(stage3array[0],stage3ref);
	    Stage3_print_mutations(stage3array[0],stage3ref,chromosome_iit,queryseq,
				   dbversion,printtype,proteinmode,
				   invertmode,nointronlenp,wraplength,
				   /*snps_p*/snp_blocks ? true : false,
				   /*maxmutations*/1000000);
	    for (i = 0; i < npaths; i++) {
	      stage3 = stage3array[i];
	      Stage3_free(&stage3);
	    }
	    FREE(stage3array);

	  }

#ifndef PMAP
	}
#endif

	Oligoindex_clear_inquery(Oligoindex_array_elt(oligoindices_major,0));

	Sequence_free(&queryuc);
      }
      Sequence_free(&queryseq);
      jobid++;
    }
    Sequence_free(&genomicseg);
    Stage3_free(&stage3ref);
  }

  Stopwatch_free(&stopwatch);
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

#endif

void
check_map_iit (IIT_T map_iit, Univ_IIT_T chromosome_iit) {
  char *typestring, *lookup, *p;
  int type, destranded_len;
  bool errorp = false;

  for (type = 1; type < IIT_ntypes(map_iit); type++) {
    lookup = typestring = IIT_typestring(map_iit,type);
    if ((p = rindex(typestring,'+')) != NULL) {
      destranded_len = (p - typestring)/sizeof(char);
      lookup = (char *) CALLOC(destranded_len+1,sizeof(char));
      strncpy(lookup,typestring,destranded_len);
    } else if ((p = rindex(typestring,'-')) != NULL) {
      destranded_len = (p - typestring)/sizeof(char);
      lookup = (char *) CALLOC(destranded_len+1,sizeof(char));
      strncpy(lookup,typestring,destranded_len);
    }

    if (Univ_IIT_find_one(chromosome_iit,lookup) < 0) {
      if (p != NULL) {
	fprintf(stderr,"Warning: In %s, type %s (without the %s) does not correspond to a known chromosome in %s.\n",
		map_iitfile,typestring,p,dbversion);
      } else {
	fprintf(stderr,"Warning: In %s, type %s does not correspond to a known chromosome in %s.\n",
		map_iitfile,typestring,dbversion);
      }
      errorp = true;
    }

    if (p != NULL) {
      FREE(lookup);
    }
  }
  if (errorp == true) {
    fprintf(stderr,"Known chromosomes: ");
    Univ_IIT_dump_labels(stderr,chromosome_iit);
  }
  return;
}


void
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


static int
parse_command_line (int argc, char *argv[], int optind) {
  int opt, c;
  extern char *optarg;
  int long_option_index = 0;
  const char *long_name;
  char **argstart;

  int len;
  int user_ngap = -1;


  fprintf(stderr,"GMAP version %s called with args:",PACKAGE_VERSION);
  argstart = &(argv[-optind]);
  for (c = 1; c < argc + optind; c++) {
    fprintf(stderr," %s",argstart[c]);
  }
  fprintf(stderr,"\n");

  while ((opt = getopt_long(argc,argv,
#ifdef PMAP
			    "q:D:a:d:k:Gg:2B:K:w:L:x:1t:s:c:H:SA03468:9n:f:ZO5o:V:v:M:m:ebu:E:PQYNI:i:l:",
#else
			    "q:D:d:k:Gg:2B:K:w:L:x:1t:s:c:H:p:SA03468:9n:f:ZO5o:V:v:M:m:ebu:E:PQFa:Tz:j:YNI:i:l:",
#endif
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

      } else if (!strcmp(long_name,"cmdline")) {
	user_cmdline = optarg;

      } else if (!strcmp(long_name,"suboptimal-score")) {
	suboptimal_score = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"require-splicedir")) {
	require_splicedir_p = true;

      } else if (!strcmp(long_name,"splicingdir")) {
	user_splicingdir = optarg;

      } else if (!strcmp(long_name,"nosplicing")) {
	novelsplicingp = false;

      } else if (!strcmp(long_name,"no-chimeras")) {
	chimera_margin = 0;

      } else if (!strcmp(long_name,"min-intronlength")) {
	min_intronlength = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"allow-close-indels")) {
	if (!strcmp(optarg,"0")) {
	  /* Disallow */
	  close_indels_mode = -1;
	  extraband_single = 0;
	} else if (!strcmp(optarg,"1")) {
	  /* Always allow */
	  close_indels_mode = +1;
	  extraband_single = 3;
	} else if (!strcmp(optarg,"2")) {
	  /* Allow for high-quality alignments */
	  close_indels_mode = 0;
	  extraband_single = 3;
	} else {
	  fprintf(stderr,"allow-close-indels argument %s not recognized.  Only allow 0, 1, or 2.  Run 'gsnap --help' for more information.\n",optarg);
	  return 9;
	}
      } else if (!strcmp(long_name,"microexon-spliceprob")) {
	microexon_spliceprob = check_valid_float(optarg,long_name);
      } else if (!strcmp(long_name,"stage2-start")) {
	suboptimal_score_start = atoi(check_valid_int(optarg));
      } else if (!strcmp(long_name,"stage2-end")) {
	suboptimal_score_end = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"canonical-mode")) {
	if (!strcmp(optarg,"0")) {
	  canonical_mode = 0;
	} else if (!strcmp(optarg,"1")) {
	  canonical_mode = 1;
	} else if (!strcmp(optarg,"2")) {
	  canonical_mode = 2;
	} else {
	  fprintf(stderr,"Canonical level %s not recognized.\n",optarg);
	  fprintf(stderr,"0=low reward for canonical introns, 1=high reward for canonical introns (default)\n");
	  fprintf(stderr,"2=low reward for high-identity seqs, high reward otherwise\n");
	  return 9;
	}

      } else if (!strcmp(long_name,"cross-species")) {
	cross_species_p = true;

      } else if (!strcmp(long_name,"homopolymer")) {
	homopolymerp = true;

      } else if (!strcmp(long_name,"cmetdir")) {
	user_cmetdir = optarg;

      } else if (!strcmp(long_name,"atoidir")) {
	user_atoidir = optarg;

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

      } else if (!strcmp(long_name,"min-trimmed-coverage")) {
	min_trimmed_coverage = check_valid_float(optarg,long_name);
      } else if (!strcmp(long_name,"min-identity")) {
	min_identity = check_valid_float(optarg,long_name);

      } else if (!strcmp(long_name,"input-buffer-size")) {
	inbuffer_nspaces = atoi(check_valid_int(optarg));
      } else if (!strcmp(long_name,"output-buffer-size")) {
	output_buffer_size = atoi(check_valid_int(optarg));
      } else if (!strcmp(long_name,"print-comment")) {
	print_comment_p = true;
      } else if (!strcmp(long_name,"failsonly")) {
	if (nofailsp == true) {
	  fprintf(stderr,"Cannot specify both --nofails and --failsonly\n");
	  return 9;
	} else {
	  failsonlyp = true;
	}
      } else if (!strcmp(long_name,"failed-input")) {
	failedinput_root = optarg;
#if 0
      } else if (!strcmp(long_name,"quiet-if-excessive")) {
	quiet_if_excessive_p = true;
#endif
      } else if (!strcmp(long_name,"nofails")) {
	if (failsonlyp == true) {
	  fprintf(stderr,"Cannot specify both --nofails and --failsonly\n");
	  return 9;
	} else {
	  nofailsp = true;
	}
      } else if (!strcmp(long_name,"split-output")) {
	split_output_root = optarg;
      } else if (!strcmp(long_name,"append-output")) {
	appendp = true;
      } else if (!strcmp(long_name,"gff3-add-separators")) {
	if (!strcmp(optarg,"1")) {
	  gff3_separators_p = true;
	} else if (!strcmp(optarg,"0")) {
	  gff3_separators_p = false;
	} else {
	  fprintf(stderr,"--gff3-add-separators flag must be 0 or 1\n");
	  return 9;
	}

#ifndef PMAP
      } else if (!strcmp(long_name,"no-sam-headers")) {
	sam_headers_p = false;
      } else if (!strcmp(long_name,"sam-use-0M")) {
	sam_insert_0M_p = true;
      } else if (!strcmp(long_name,"quality-protocol")) {
	if (user_quality_shift == true) {
	  fprintf(stderr,"Cannot specify both -j (--quality-print-shift) and --quality-protocol\n");
	  return 9;
	} else if (!strcmp(optarg,"illumina")) {
	  quality_shift = -31;
	  user_quality_shift = true;
	} else if (!strcmp(optarg,"sanger")) {
	  quality_shift = 0;
	  user_quality_shift = true;
	} else {
	  fprintf(stderr,"The only values allowed for --quality-protocol are illumina or sanger\n");
	  return 9;
	}
      } else if (!strcmp(long_name,"force-xs-dir")) {
	force_xs_direction_p = true;
      } else if (!strcmp(long_name,"md-lowercase-snp")) {
	md_lowercase_variant_p = true;
      } else if (!strcmp(long_name,"read-group-id")) {
	sam_read_group_id = optarg;
      } else if (!strcmp(long_name,"read-group-name")) {
	sam_read_group_name = optarg;
      } else if (!strcmp(long_name,"read-group-library")) {
	sam_read_group_library = optarg;
      } else if (!strcmp(long_name,"read-group-platform")) {
	sam_read_group_platform = optarg;
#endif
      } else {
	/* Shouldn't reach here */
	fprintf(stderr,"Don't recognize option %s.  For usage, run 'gsnap --help'",long_name);
	return 9;
      }
      break;

    case 'q': parse_part(&part_modulus,&part_interval,optarg); break;
    case 'D': user_genomedir = optarg; break;
    case 'd': 
      dbroot = (char *) CALLOC(strlen(optarg)+1,sizeof(char));
      strcpy(dbroot,optarg);
      break;
#ifdef PMAP
    case 'a':
      if ((required_alphabet = Alphabet_find(optarg)) == AA0) {
	return 9;
      }
      break;
    case 'k': required_index1part = atoi(check_valid_int(optarg)); break;
#else
    case 'k':
      required_index1part = atoi(check_valid_int(optarg));
      if (required_index1part > 16) {
	fprintf(stderr,"The value for k-mer size must be 16 or less\n");
	return 9;
      }
      break;
#endif
    case 'G': uncompressedp = true; break;
    case 'g': user_genomicseg = optarg; break;
    case '1': user_selfalign_p = true; break;
    case '2': user_pairalign_p = true; break;

    case 'B': 
      if (!strcmp(optarg,"5")) {
	fprintf(stderr,"Note: Batch mode 5 is now the same as batch mode 4.\n");
	fprintf(stderr,"Expansion of offsets is now controlled separately by --expand-offsets (default=1).\n");
	offsetsstrm_access = USE_ALLOCATE; /* Doesn't matter */
	positions_access = USE_ALLOCATE;
	genome_access = USE_ALLOCATE;
      } else if (!strcmp(optarg,"4")) {
	offsetsstrm_access = USE_ALLOCATE;
	positions_access = USE_ALLOCATE;
	genome_access = USE_ALLOCATE;
#ifdef HAVE_MMAP
      } else if (!strcmp(optarg,"3")) {
	offsetsstrm_access = USE_ALLOCATE;
	positions_access = USE_ALLOCATE;
	genome_access = USE_MMAP_PRELOAD; /* was batch_genome_p = true */
      } else if (!strcmp(optarg,"2")) {
	offsetsstrm_access = USE_ALLOCATE; /* was batch_offsets_p = true */
	positions_access = USE_MMAP_PRELOAD; /* was batch_positions_p = true */
	genome_access = USE_MMAP_PRELOAD; /* was batch_genome_p = true */
      } else if (!strcmp(optarg,"1")) {
	offsetsstrm_access = USE_ALLOCATE; /* was batch_offsets_p = true */
	positions_access = USE_MMAP_PRELOAD; /* was batch_positions_p = true */
	genome_access = USE_MMAP_ONLY; /* was batch_genome_p = false */
      } else if (!strcmp(optarg,"0")) {
	offsetsstrm_access = USE_ALLOCATE; /* was batch_offsets_p = true */
	positions_access = USE_MMAP_ONLY; /* was batch_positions_p = false */
	genome_access = USE_MMAP_ONLY; /* was batch_genome_p = false */
#endif

      } else {
#ifdef HAVE_MMAP
	fprintf(stderr,"Batch mode %s not recognized.  Only allow 0-5.  Run 'gmap --help' for more information.\n",optarg);
#else
	fprintf(stderr,"Batch mode %s not recognized.  Only allow 4-5, since mmap is disabled.  Run 'gmap --help' for more information.\n",optarg);
#endif
	return 9;
      }
      break;

    case 'K': maxintronlen = atoi(check_valid_int(optarg)); break;
    case 'w': shortsplicedist = strtoul(check_valid_int(optarg),NULL,10); break;

    case 'L': maxtotallen_bound = atoi(check_valid_int(optarg)); break;
    case 'x':
#ifdef PMAP
      chimera_margin = atoi(check_valid_int(optarg))/3; 
#else
      chimera_margin = atoi(check_valid_int(optarg)); 
#endif
      if (chimera_margin <= 0) {
	/* Disable finding of chimeras */
#if 0
      } else if (chimera_margin < CHIMERA_SLOP) {
	/* Not sure why chimera_margin should be tied to CHIMERA_SLOP */
	chimera_margin = CHIMERA_SLOP;
#endif
      }
      break;
      /* case 'w': referencefile = optarg; break; */

#ifdef HAVE_PTHREAD
    case 't': nworkers = atoi(check_valid_int(optarg)); break;
#else
    case 't': fprintf(stderr,"This version of GMAP has pthreads disabled, so ignoring the value of %s for -t\n",optarg); break;
#endif

    case 's': splicing_file = optarg; knownsplicingp = true; break;
    case 'c': user_chrsubsetname = optarg; break;
    case 'H': minendexon = atoi(check_valid_int(optarg)); break;

#ifndef PMAP
    case 'p': switch (atoi(check_valid_int(optarg))) {
      case 0: prune_poor_p = false, prune_repetitive_p = false; break;
      case 1: prune_poor_p = true; prune_repetitive_p = false; break;
      case 2: prune_poor_p = false; prune_repetitive_p = true; break;
      case 3: prune_poor_p = true; prune_repetitive_p = true; break;
      default: fprintf(stderr,"Prune level %s not recognized.\n",optarg);
	fprintf(stderr,"0=no pruning, 1=poor seqs, 2=repetitive seqs, 3=both poor and repetitive seqs (default)\n");
	return 9;
      }
      break;
#endif

    case 'S': printtype = SUMMARY; break;
    case 'A': printtype = ALIGNMENT; break;
    case '0': exception_raise_p = false; break; /* Allows signals to pass through */
    case '3': printtype = CONTINUOUS; break;
    case '4': printtype = CONTINUOUS_BY_EXON; break;
    case '6': debug_graphic_p = true; break;
    case '8':
      if (!strcmp(optarg,"stage1")) {
	stage1debug = true;
      } else if (!strcmp(optarg,"diag")) {
	diag_debug = true;
      } else if (!strcmp(optarg,"stage2")) {
	stage3debug = POST_STAGE2;
      } else if (!strcmp(optarg,"singles")) {
	stage3debug = POST_SINGLES;
      } else if (!strcmp(optarg,"introns")) {
	stage3debug = POST_INTRONS;
      } else if (!strcmp(optarg,"hmm")) {
	stage3debug = POST_HMM;
      } else if (!strcmp(optarg,"smoothing")) {
	stage3debug = POST_SMOOTHING;
      } else if (!strcmp(optarg,"dualintrons")) {
	stage3debug = POST_DUAL_INTRONS;
      } else if (!strcmp(optarg,"cycles")) {
	stage3debug = POST_CYCLES;
      } else if (!strcmp(optarg,"dualbreaks")) {
	stage3debug = POST_DUAL_BREAKS;
      } else if (!strcmp(optarg,"middle")) {
	stage3debug = POST_MIDDLE;
      } else if (!strcmp(optarg,"ends")) {
	stage3debug = POST_ENDS;
      } else if (!strcmp(optarg,"canonical")) {
	stage3debug = POST_CANONICAL;
      } else if (!strcmp(optarg,"trim")) {
	stage3debug = POST_CANONICAL;
      } else if (!strcmp(optarg,"changepoint")) {
	stage3debug = POST_CHANGEPOINT;
      } else if (!strcmp(optarg,"distalmedial")) {
	stage3debug = POST_DISTAL_MEDIAL;
      } else {
	fprintf(stderr,"Allowed arguments for -8 flag are stage2, smoothing, singles, introns, hmm, dualbreaks, cycles, canonical, changepoint, distalmedial\n");
	return 9;
      }
      break;
    case '9': checkp = true; break;
    case 'n':
      maxpaths_report = atoi(check_valid_int(optarg));
      if (maxpaths_report == 1) {
	fprintf(stderr,"Note: -n 1 will not report chimeric alignments.  If you want a single alignment plus chimeras, use -n 0 instead.\n");
      }
      break;
    case 'f':
      if (!strcmp(optarg,"1") || !strcmp(optarg,"psl_nt")) {
	printtype = PSL_NT;
#ifdef PMAP
      } else if (!strcmp(optarg,"0") || !strcmp(optarg,"psl_pro")) {
	printtype = PSL_PRO;
#else
      } else if (!strcmp(optarg,"psl")) {
	printtype = PSL_NT;
      } else if (!strcmp(optarg,"6") || !strcmp(optarg,"splicesites")) {
	printtype = SPLICESITES;
      } else if (!strcmp(optarg,"introns")) {
	printtype = INTRONS;
      } else if (!strcmp(optarg,"samse")) {
	printtype = SAM;
	sam_paired_p = false;
      } else if (!strcmp(optarg,"sampe")) {
	printtype = SAM;
	sam_paired_p = true;
#endif
      } else if (!strcmp(optarg,"2") || !strcmp(optarg,"gff3_gene")) {
	printtype = GFF3_GENE;
      } else if (!strcmp(optarg,"3") || !strcmp(optarg,"gff3_match_cdna")) {
	printtype = GFF3_MATCH_CDNA;
      } else if (!strcmp(optarg,"4") || !strcmp(optarg,"gff3_match_est")) {
	printtype = GFF3_MATCH_EST;
      } else if (!strcmp(optarg,"7") || !strcmp(optarg,"map_exons")) {
	printtype = MAP_EXONS;
      } else if (!strcmp(optarg,"8") || !strcmp(optarg,"map_ranges")) {
	printtype = MAP_RANGES;
      } else if (!strcmp(optarg,"9") || !strcmp(optarg,"coords")) {
	printtype = COORDS;
      } else {
	fprintf(stderr,"Output format \"%s\" not recognized.  Allowed formats are:\n",optarg);
	fprintf(stderr,"  psl_nt (1)\n");
#ifdef PMAP
	fprintf(stderr,"  psl_pro (0)\n");
#else
	fprintf(stderr,"  psl\n");
	fprintf(stderr,"  splicesites (6)\n");
	fprintf(stderr,"  introns\n");
	fprintf(stderr,"  samse\n");
	fprintf(stderr,"  sampe\n");
#endif
	fprintf(stderr,"  gff3_gene (2)\n");
	fprintf(stderr,"  gff3_match_cdna (3)\n");
	fprintf(stderr,"  gff3_match_est (4)\n");
	fprintf(stderr,"  map_exons (7)\n");
	fprintf(stderr,"  map_ranges (8)\n");
	fprintf(stderr,"  coords (9)\n");
	return 9;
      }
      break;
    case 'Z': printtype = COMPRESSED; break;
    case 'O': orderedp = true; break;
    case '5': checksump = true; break;
    case 'o': chimera_overlap = atoi(check_valid_int(optarg)); break;

    case 'V': user_snpsdir = optarg; break;
    case 'v': snps_root = optarg; break;

    case 'M': user_mapdir = optarg; break;
    case 'm': 
      map_iitfile = (char *) CALLOC(strlen(optarg)+1,sizeof(char));
      strcpy(map_iitfile,optarg);
      if ((len = strlen(map_iitfile)) > 4 && strcmp(&(map_iitfile[len-4]),".iit") == 0) {
	map_iitfile[len-4] = '\0';
      }
      break;

    case 'e': map_exons_p = true; break;
    case 'b': map_bothstrands_p = true; break;
    case 'u': nflanking = atoi(check_valid_int(optarg)); break;

    case 'E': 
      if (!strcmp(optarg,"cdna")) {
	printtype = EXONS_CDNA;
      } else if (!strncmp(optarg,"genomic",strlen("genomic"))) {
	printtype = EXONS_GENOMIC;
      } else {
	fprintf(stderr,"Argument to -E flag must be either \"cdna\" or \"genomic\"\n");
	return 9;
      }
      break;

#ifdef PMAP
    case 'P': printtype = PROTEIN_GENOMIC; break;
    case 'Q': printtype = CDNA; break; 
#else
    case 'P': printtype = CDNA; break;
    case 'Q': printtype = PROTEIN_GENOMIC; break;
    case 'F': fulllengthp = true; break;
    case 'a': cds_startpos = atoi(check_valid_int(optarg)); break;
    case 'T': truncatep = true; fulllengthp = true; break;
    case 'z':
      if (!strcmp(optarg,"sense_force")) {
	sense_try = +1;
	sense_filter = 0;
      } else if (!strcmp(optarg,"antisense_force")) {
	sense_try = -1;
	sense_filter = 0;
      } else if (!strcmp(optarg,"sense_filter")) {
	sense_try = 0;
	sense_filter = +1;
      } else if (!strcmp(optarg,"antisense_filter")) {
	sense_try = 0;
	sense_filter = -1;
      } else if (!strcmp(optarg,"auto")) {
	sense_try = 0;
	sense_filter = 0;
      } else {
	fprintf(stderr,"direction %s not recognized.  Must be sense_force, antisense_force, sense_filter, antisense_filter, or auto\n",optarg);
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

#endif
    case 'Y': strictp = false; break;
    case 'N': nointronlenp = true; break;
    case 'I': invertmode = atoi(check_valid_int(optarg)); break;
    case 'i': user_ngap = atoi(check_valid_int(optarg)); break;
    case 'l': wraplength = atoi(check_valid_int(optarg)); break;

    case '?': fprintf(stderr,"For usage, run 'gmap --help'\n"); return 9;
    default: return 9;
    }
  }

  if (printtype == SPLICESITES || printtype == INTRONS) {
    if (maxpaths_report > 1 || (sense_try != +1 && sense_filter != +1)) {
      fprintf(stderr,"For splicesites or introns output, you should probably add flags '-n 1' and either '-z sense_force' or '-z sense_filter'.\n");
    }
  }
  
  if (user_ngap >= 0) {
    ngap = user_ngap;
  } else if (printtype == EXONS_CDNA || printtype == EXONS_GENOMIC) {
    /* If user didn't specify, then set to zero */
    ngap = 0;
  };

  if (maxintronlen > maxtotallen_bound) {
    maxintronlen = maxtotallen_bound;
  }

#ifdef HAVE_PTHREAD
#ifdef USE_DIAGPOOL
  if (diag_debug == true && nworkers > 0) {
    fprintf(stderr,"For diag output, must specify 0 threads\n");
    exit(9);
  }
#endif
#endif

  if (user_cmdline != NULL) {
    part_modulus = 0;
    part_interval = 1;
    inbuffer_nspaces = 0;
    nchromosomes = 1;
    dbroot = (char *) NULL;
  } else if (user_selfalign_p == true) {
    nchromosomes = 1;
    dbroot = (char *) NULL;
  } else if (user_pairalign_p == true) {
    nchromosomes = 1;
    dbroot = (char *) NULL;
  } else if (user_genomicseg != NULL) {
    /* Ignore -D and -d flags */
    nchromosomes = 1;
    dbroot = (char *) NULL;
  } else if (dbroot == NULL) {
    fprintf(stderr,"Need to specify the -d, -g, -1, -2, or --cmdline flag\n");
    print_program_usage();
    return 9;
  } else if (!strcmp(dbroot,"?")) {
    Datadir_avail_gmap_databases(stdout,user_genomedir);
    return 1;
  }

#ifndef PMAP
  if (printtype == SAM) {
    if (sam_read_group_id == NULL && sam_read_group_name != NULL) {
      sam_read_group_id = sam_read_group_name;
    } else if (sam_read_group_id != NULL && sam_read_group_name == NULL) {
      sam_read_group_name = sam_read_group_id;
    }
  }
#endif

  return 0;
}


static Inbuffer_T
open_input_stream (int *nread, Sequence_T *usersegment, int argc, char **argv) {
  Inbuffer_T inbuffer;
  int nextchar = '\0';
  FILE *input = NULL;
  char **files;
  int nfiles;

  Request_T request;
  char *p;

  /* Read user segment before rest of sequences, because of shared usage of sequence.c */
  if (user_cmdline != NULL) {
    p = user_cmdline;
    while (*p != '\0' && *p != ',') {
      p++;
    }
    if (*p == '\0') {
      fprintf(stderr,"--cmdline requires two strings separated by a comma");
      exit(9);
    } else {
      *usersegment = global_usersegment = Sequence_genomic_new(user_cmdline,(int) (p - user_cmdline),/*copyp*/true);
      if ((min_matches = Sequence_fulllength(*usersegment)/2) > MIN_MATCHES) {
	min_matches = MIN_MATCHES;
      }
      p++;
    }

  } else if (user_selfalign_p == true) {
    /* usersegment will be assigned to query sequence below */

  } else if (user_pairalign_p == true) {
    /* Unfortunately, this procedure reads header of queryseq */
    *usersegment = Sequence_read_unlimited(&nextchar,stdin);
    if ((min_matches = Sequence_fulllength(*usersegment)/2) > MIN_MATCHES) {
      min_matches = MIN_MATCHES;
    }

  } else if (user_genomicseg != NULL) {
    if ((input = FOPEN_READ_TEXT(user_genomicseg)) == NULL) {
      fprintf(stderr,"Can't open file %s\n",user_genomicseg);
      exit(9);
    }
    if ((*usersegment = global_usersegment = Sequence_read_unlimited(&nextchar,input)) == NULL) {
      fprintf(stderr,"File %s is empty\n",user_genomicseg);
      exit(9);
    } else {
      genomelength = (Univcoord_T) Sequence_fulllength(*usersegment);
    }
    
    if ((min_matches = Sequence_fulllength(*usersegment)/2) > MIN_MATCHES) {
      min_matches = MIN_MATCHES;
    }
    fclose(input);

  } else {
    min_matches = MIN_MATCHES;
  }

  Inbuffer_setup(/*filter_if_both_p*/false,user_pairalign_p,global_usersegment,
		 part_modulus,part_interval);
  if (user_cmdline != NULL) {
    inbuffer = Inbuffer_cmdline(p,strlen(p));
    *nread = 1;

  } else if (user_selfalign_p == true) {
      input = stdin;
      files = (char **) NULL;
      nfiles = 0;

      /* Read in first batch of sequences */
      inbuffer = Inbuffer_new(nextchar,input,files,nfiles,inbuffer_nspaces);
      *nread = Inbuffer_fill_init(inbuffer);
      request = Inbuffer_first_request(inbuffer); /* Need usersegment, not the request itself */
      *usersegment = Request_queryseq(request);

  } else {
    /* Open input stream and peek at first char */
    if (user_pairalign_p == true) {
      input = stdin;
      files = (char **) NULL;
      nfiles = 0;
      inbuffer_nspaces = 1;
    } else if (argc == 0) {
      fprintf(stderr,"Reading from stdin\n");
      input = stdin;
      files = (char **) NULL;
      nfiles = 0;
    } else {
      input = NULL;
      files = argv;
      nfiles = argc;
    }

    /* Read in first batch of sequences */
    inbuffer = Inbuffer_new(nextchar,input,files,nfiles,inbuffer_nspaces);
#ifdef USE_MPI
    *nread = 0;
#else
    *nread = Inbuffer_fill_init(inbuffer);
#endif
  }

  return inbuffer;
}


int
main (int argc, char *argv[]) {
#ifdef USE_MPI
  int nbeyond;
#else
  bool multiple_sequences_p = false;
#endif
  int cmdline_status;

  char *genomesubdir = NULL, *snpsdir = NULL, *modedir = NULL, *mapdir = NULL, *iitfile = NULL, *fileroot = NULL;
  int divno;
  Univinterval_T interval;
  Sequence_T usersegment = NULL;

  bool showcontigp = true;
  int nread;
  double runtime;

  Splicestringpool_T splicestringpool;

#ifdef HAVE_PTHREAD
  int ret, i;
  pthread_attr_t thread_attr_join;
#ifdef WORKER_DETACH
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


#ifdef USE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

  if ((n_worker_procs = nprocs - 1) == 0) {
    if (myid == 0) {
      fprintf(stderr,"Need at least 2 processes for MPI version\n");
    }
    MPI_Finalize();
    exit(0);

  } else {
    MPI_Debug_setup(myid);
  }
#endif

  cmdline_status = parse_command_line(argc,argv,optind);
  argc -= optind;
  argv += optind;
      
  if (cmdline_status == 0) {
    /* okay to continue */
  } else if (cmdline_status == 1) {
    /* only information needed */
#ifdef USE_MPI
    MPI_Finalize();
#endif
    exit(0);
  } else {
#ifdef USE_MPI
    MPI_Finalize();
#endif
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

    sigaction(SIGFPE,&signal_action,NULL);
    sigaction(SIGSEGV,&signal_action,NULL);
    sigaction(SIGTRAP,&signal_action,NULL);
    sigaction(SIGUSR1,&signal_action,NULL);
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
  if (myid > 0) {
    inbuffer = open_input_stream(&nread,&usersegment,argc,argv);
  }

#else
  inbuffer = open_input_stream(&nread,&usersegment,argc,argv);

  if (nread > 1) {
    multiple_sequences_p = true;
#if 0
#ifdef HAVE_MMAP
    if (offsetsstrm_access != USE_ALLOCATE || genome_access != USE_ALLOCATE) {
      fprintf(stderr,"Note: >1 sequence detected, so index files are being memory mapped.\n");
      fprintf(stderr,"  GMAP can run slowly at first while the computer starts to accumulate\n");
      fprintf(stderr,"  pages from the hard disk into its cache.  To copy index files into RAM\n");
      fprintf(stderr,"  instead of memory mapping, use -B 3, -B 4, or -B 5, if you have enough RAM.\n");
#ifdef HAVE_PTHREAD
      fprintf(stderr,"  For more speed, also try multiple threads (-t <int>), if you have multiple processors or cores.");
#endif
      fprintf(stderr,"\n");
#endif
    }
#endif

  } else {
    /* multiple_sequences_p = false; */
    /* fprintf(stderr,"Note: only 1 sequence detected.  Ignoring batch (-B) command\n"); */
    expand_offsets_p = false;
#ifdef HAVE_MMAP
    offsetsstrm_access = USE_MMAP_ONLY;
    positions_access = USE_MMAP_ONLY;
    genome_access = USE_MMAP_ONLY;
#else
    offsetsstrm_access = USE_ALLOCATE;
    positions_access = USE_ALLOCATE;
    genome_access = USE_ALLOCATE;
#endif
  }

#endif


  if (dbroot != NULL) {
    /* Prepare genomic data */
    genomesubdir = Datadir_find_genomesubdir(&fileroot,&dbversion,user_genomedir,dbroot);

    iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			      strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
    sprintf(iitfile,"%s/%s.chromosome.iit",genomesubdir,fileroot);
    if ((chromosome_iit = Univ_IIT_read(iitfile,/*readonlyp*/true,/*add_iit_p*/false)) == NULL) {
      fprintf(stderr,"IIT file %s is not valid\n",iitfile);
      exit(9);
#ifdef LARGE_GENOMES
    } else if (Univ_IIT_coord_values_8p(chromosome_iit) == false) {
      fprintf(stderr,"This program gmapl is designed for large genomes.\n");
      fprintf(stderr,"For small genomes of less than 2^32 (4 billion) bp, please run gmap instead.\n");
      exit(9);
#endif
    } else {
      nchromosomes = Univ_IIT_total_nintervals(chromosome_iit);
      circular_typeint = Univ_IIT_typeint(chromosome_iit,"circular");
      circularp = Univ_IIT_circularp(&any_circular_p,chromosome_iit);
    }
    genomelength = Univ_IIT_genomelength(chromosome_iit,/*with_circular_alias*/false);

    FREE(iitfile);
  }

#ifdef USE_MPI
  /* Can prevent loading of files by rank 0 process */
#endif

  if (map_iitfile == NULL) {
    /* Skip */
  } else if (!strcmp(map_iitfile,"?")) {
    Datadir_avail_maps(stdout,user_mapdir,genomesubdir,fileroot);
    exit(0);
  } else {
    mapdir = Datadir_find_mapdir(user_mapdir,genomesubdir,fileroot);
    iitfile = (char *) CALLOC(strlen(mapdir)+strlen("/")+
			      strlen(map_iitfile)+strlen(".iit")+1,sizeof(char));
    sprintf(iitfile,"%s/%s.iit",mapdir,map_iitfile);
    if ((map_iit = IIT_read(iitfile,/*name*/map_iitfile,/*readonlyp*/true,/*divread*/READ_ALL,
			    /*divstring*/NULL,/*add_iit_p*/true,/*labels_read_p*/true)) == NULL) {
      fprintf(stderr,"Map file %s.iit not found in %s.  Available files:\n",map_iitfile,mapdir);
      Datadir_list_directory(stderr,mapdir);
      fprintf(stderr,"Either install file %s.iit or specify a directory for the IIT file\n",iitfile);
      fprintf(stderr,"using the -M flag.\n");
      exit(9);
    } else {
      map_divint_crosstable = Univ_IIT_divint_crosstable(chromosome_iit,map_iit);
    }

    check_map_iit(map_iit,chromosome_iit);

    FREE(iitfile);
    FREE(mapdir);
    FREE(map_iitfile);
  }

  if (splicing_file != NULL) {
    if (user_splicingdir == NULL) {
      if ((splicing_iit = IIT_read(splicing_file,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				   /*divstring*/NULL,/*add_iit_p*/false,/*labels_read_p*/true)) != NULL) {
	fprintf(stderr,"Reading splicing file %s locally...",splicing_file);
      } else {
	iitfile = (char *) CALLOC(strlen(user_splicingdir)+strlen("/")+strlen(splicing_file)+1,sizeof(char));
	sprintf(iitfile,"%s/%s",user_splicingdir,splicing_file);
	if ((splicing_iit = IIT_read(splicing_file,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				     /*divstring*/NULL,/*add_iit_p*/false,/*labels_read_p*/true)) != NULL) {
	  fprintf(stderr,"Reading splicing file %s locally...",splicing_file);
	  FREE(iitfile);
	}
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
  }

  /* Complement_init(); */
  Dynprog_init(mode);
#ifdef PMAP
  Backtranslation_init();
#endif

  if (user_pairalign_p == true) {
    showcontigp = false;
    /* maxpaths_report = 1; -- no; could have different paths against the user segment. */

    genomecomp = (Genome_T) NULL;
    genomebits = (Genome_T) NULL;
    genomecomp_alt = (Genome_T) NULL;
    genomebits_alt = (Genome_T) NULL;
    dbversion = (char *) NULL;
    /* Do for each usersegment */

  } else if (global_usersegment != NULL) {
    /* Map against user-provided genomic segment */
    showcontigp = false;
    /* maxpaths_report = 1; -- no; could have different paths against the user segment. */

    genomecomp = (Genome_T) NULL;
    genomebits = (Genome_T) NULL;
    genomecomp_alt = (Genome_T) NULL;
    genomebits_alt = (Genome_T) NULL;
    dbversion = (char *) NULL;
    genomecomp_blocks = Compress_create_blocks_comp(Sequence_fullpointer(global_usersegment),Sequence_fulllength(global_usersegment));
    genomebits_blocks = Compress_create_blocks_bits(genomecomp_blocks,Sequence_fulllength(global_usersegment));

    if (Sequence_fulllength(global_usersegment) > 1000000) {
      fprintf(stderr,"Genomic sequence is unusually long (%d bp).  GMAP handles genomes better when\n",
	      Sequence_fulllength(global_usersegment));
      fprintf(stderr,"  they are converted into gmap databases first using gmap_setup, and then accessed\n");
      fprintf(stderr,"  with the -d flag.\n");
    }

#ifdef PMAP
  } else {
    /* Map against genome */
    if (showcontigp == true) {
      iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
				strlen(fileroot)+strlen(".contig.iit")+1,sizeof(char));
      sprintf(iitfile,"%s/%s.contig.iit",genomesubdir,fileroot);
      if ((contig_iit = Univ_IIT_read(iitfile,/*readonlyp*/true,/*add_iit_p*/false)) == NULL) {
	fprintf(stderr,"IIT file %s is not valid\n",iitfile);
	exit(9);
      }
      FREE(iitfile);
    }
  
    genomecomp = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,/*genometype*/GENOME_OLIGOS,
			    uncompressedp,genome_access,sharedp);
    genomebits = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,/*genometype*/GENOME_BITS,
			    uncompressedp,genome_access,sharedp);
    if (snps_root == NULL) {
      genomecomp_alt = genomebits_alt = (Genome_T) NULL;
    } else {
      genomecomp_alt = Genome_new(genomesubdir,fileroot,snps_root,/*genometype*/GENOME_OLIGOS,
				  uncompressedp,genome_access,sharedp);
      genomebits_alt = Genome_new(genomesubdir,fileroot,snps_root,/*genometype*/GENOME_BITS,
				  uncompressedp,genome_access,sharedp);
    }

    indexdb_fwd = Indexdb_new_genome(&index1part_aa,&index1interval,
				     genomesubdir,fileroot,FWD_FILESUFFIX,/*snps_root*/NULL,
				     &alphabet,&alphabet_size,required_alphabet,
				     required_index1part,required_index1interval,
				     expand_offsets_p,offsetsstrm_access,positions_access,sharedp);
    indexdb_rev = Indexdb_new_genome(&index1part_aa,&index1interval,
				     genomesubdir,fileroot,REV_FILESUFFIX,/*snps_root*/NULL,
				     &alphabet,&alphabet_size,required_alphabet,
				     required_index1part,required_index1interval,
				     expand_offsets_p,offsetsstrm_access,positions_access,sharedp);

    if (indexdb_fwd == NULL || indexdb_rev == NULL) {
      fprintf(stderr,"Cannot find offsets file %s.%s*offsets or %s.%s*offsets.\n",
	      fileroot,FWD_FILESUFFIX,fileroot,REV_FILESUFFIX);
      fprintf(stderr,"You may need to run 'pmapindex -d %s' to build the indices needed for PMAP.\n",
	      fileroot);
      exit(9);
    } else if (user_chrsubsetname != NULL) {
      if ((divno = Univ_IIT_find_one(chromosome_iit,user_chrsubsetname)) < 0) {
	fprintf(stderr,"Cannot find chrsubset %s in chromosome IIT file.  Ignoring.\n",user_chrsubsetname);
      } else {
	interval = Univ_IIT_interval(chromosome_iit,divno);
	chrsubset_start = Univinterval_low(interval);
	chrsubset_end = Univinterval_high(interval);
      }
    }

    FREE(genomesubdir);
    FREE(fileroot);
    FREE(dbroot);

#else
  } else if (snps_root == NULL) {
    /* Map against genome without SNPs */
    if (showcontigp == true) {
      iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
				strlen(fileroot)+strlen(".contig.iit")+1,sizeof(char));
      sprintf(iitfile,"%s/%s.contig.iit",genomesubdir,fileroot);
      if ((contig_iit = Univ_IIT_read(iitfile,/*readonlyp*/true,/*add_iit_p*/false)) == NULL) {
	fprintf(stderr,"IIT file %s is not valid\n",iitfile);
	exit(9);
      }
      FREE(iitfile);
    }
  
    genomecomp = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,/*genometype*/GENOME_OLIGOS,
			    uncompressedp,genome_access,sharedp);
    genomecomp_blocks = Genome_blocks(genomecomp);
    if ((genomebits = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,/*genometype*/GENOME_BITS,
				 uncompressedp,genome_access,sharedp)) == NULL) {
      genomebits_blocks = (Genomecomp_T *) NULL;
    } else {
      genomebits_blocks = Genome_blocks(genomebits);
    }
    genomecomp_alt = (Genome_T) NULL;
    genomebits_alt = (Genome_T) NULL;

    if (mode == CMET_STRANDED || mode == CMET_NONSTRANDED) {
      if (user_cmetdir == NULL) {
	modedir = genomesubdir;
      } else {
	modedir = user_cmetdir;
      }

      if ((indexdb_fwd = Indexdb_new_genome(&index1part,&index1interval,
					    modedir,fileroot,/*idx_filesuffix*/"metct",/*snps_root*/NULL,
					    required_index1part,required_index1interval,
					    expand_offsets_p,offsetsstrm_access,positions_access,
					    sharedp)) == NULL) {
	fprintf(stderr,"Cannot find metct index file.  Need to run cmetindex first\n");
	exit(9);
      }

      if ((indexdb_rev = Indexdb_new_genome(&index1part,&index1interval,
					    modedir,fileroot,/*idx_filesuffix*/"metga",/*snps_root*/NULL,
					    required_index1part,required_index1interval,
					    expand_offsets_p,offsetsstrm_access,positions_access,
					    sharedp)) == NULL) {
	fprintf(stderr,"Cannot find metga index file.  Need to run cmetindex first\n");
	exit(9);
      }

    } else if (mode == ATOI_STRANDED || mode == ATOI_NONSTRANDED) {
      if (user_atoidir == NULL) {
	modedir = genomesubdir;
      } else {
	modedir = user_atoidir;
      }

      if ((indexdb_fwd = Indexdb_new_genome(&index1part,&index1interval,
					    modedir,fileroot,/*idx_filesuffix*/"a2iag",/*snps_root*/NULL,
					    required_index1part,required_index1interval,
					    expand_offsets_p,offsetsstrm_access,positions_access,
					    sharedp)) == NULL) {
	fprintf(stderr,"Cannot find a2iag index file.  Need to run atoiindex first\n");
	exit(9);
      }

      if ((indexdb_rev = Indexdb_new_genome(&index1part,&index1interval,
					    modedir,fileroot,/*idx_filesuffix*/"a2itc",/*snps_root*/NULL,
					    required_index1part,required_index1interval,
					    expand_offsets_p,offsetsstrm_access,positions_access,
					    sharedp)) == NULL) {
	fprintf(stderr,"Cannot find a2itc index file.  Need to run atoiindex first\n");
	exit(9);
      }

    } else if (mode == TTOC_STRANDED || mode == TTOC_NONSTRANDED) {
      if (user_atoidir == NULL) {
	modedir = genomesubdir;
      } else {
	modedir = user_atoidir;
      }

      if ((indexdb_fwd = Indexdb_new_genome(&index1part,&index1interval,
					    modedir,fileroot,/*idx_filesuffix*/"a2itc",/*snps_root*/NULL,
					    required_index1part,required_index1interval,
					    expand_offsets_p,offsetsstrm_access,positions_access,
					    sharedp)) == NULL) {
	fprintf(stderr,"Cannot find a2itc index file.  Need to run atoiindex first\n");
	exit(9);
      }

      if ((indexdb_rev = Indexdb_new_genome(&index1part,&index1interval,
					    modedir,fileroot,/*idx_filesuffix*/"a2iag",/*snps_root*/NULL,
					    required_index1part,required_index1interval,
					    expand_offsets_p,offsetsstrm_access,positions_access,
					    sharedp)) == NULL) {
	fprintf(stderr,"Cannot find a2iag index file.  Need to run atoiindex first\n");
	exit(9);
      }

    } else {
      /* Standard behavior */
      if ((indexdb_fwd = Indexdb_new_genome(&index1part,&index1interval,
					    genomesubdir,fileroot,IDX_FILESUFFIX,/*snps_root*/NULL,
					    required_index1part,required_index1interval,
					    expand_offsets_p,offsetsstrm_access,positions_access,
					    sharedp)) == NULL) {
	fprintf(stderr,"Cannot find offsets file %s.%s*offsets, needed for GSNAP\n",fileroot,IDX_FILESUFFIX);
	exit(9);
      }
      indexdb_rev = indexdb_fwd;
    }

    if (user_chrsubsetname != NULL) {
      if ((divno = Univ_IIT_find_one(chromosome_iit,user_chrsubsetname)) < 0) {
	fprintf(stderr,"Cannot find chrsubset %s in chromosome IIT file.  Ignoring.\n",user_chrsubsetname);
      } else {
	interval = Univ_IIT_interval(chromosome_iit,divno);
	chrsubset_start = Univinterval_low(interval);
	chrsubset_end = Univinterval_high(interval);
      }
    }

    FREE(genomesubdir);
    FREE(fileroot);
    FREE(dbroot);

  } else {
    /* Map against genome with SNPs */
    if (user_snpsdir == NULL) {
      snpsdir = genomesubdir;
    } else {
      snpsdir = user_snpsdir;
    }

    if (showcontigp == true) {
      iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
				strlen(fileroot)+strlen(".contig.iit")+1,sizeof(char));
      sprintf(iitfile,"%s/%s.contig.iit",genomesubdir,fileroot);
      if ((contig_iit = Univ_IIT_read(iitfile,/*readonlyp*/true,/*add_iit_p*/false)) == NULL) {
	fprintf(stderr,"IIT file %s is not valid\n",iitfile);
	exit(9);
      }
      FREE(iitfile);
    }

    genomecomp = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,/*genometype*/GENOME_OLIGOS,
			    uncompressedp,genome_access,sharedp);
    genomecomp_blocks = Genome_blocks(genomecomp);
    if ((genomebits = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,/*genometype*/GENOME_BITS,
				 uncompressedp,genome_access,sharedp)) == NULL) {
      genomebits_blocks = (Genomecomp_T *) NULL;
    } else {
      genomebits_blocks = Genome_blocks(genomebits);
    }
    genomecomp_alt = Genome_new(genomesubdir,fileroot,snps_root,/*genometype*/GENOME_OLIGOS,
				uncompressedp,genome_access,sharedp);
    genomebits_alt = Genome_new(genomesubdir,fileroot,snps_root,/*genometype*/GENOME_BITS,
				uncompressedp,genome_access,sharedp);

    if (mode == CMET_STRANDED || mode == CMET_NONSTRANDED) {
      if (user_cmetdir == NULL) {
	modedir = snpsdir;
      } else {
	modedir = user_cmetdir;
      }

      if ((indexdb_fwd = Indexdb_new_genome(&index1part,&index1interval,
					    modedir,fileroot,/*idx_filesuffix*/"metct",snps_root,
					    required_index1part,required_index1interval,
					    expand_offsets_p,offsetsstrm_access,positions_access,
					    sharedp)) == NULL) {
	fprintf(stderr,"Cannot find metct index file.  Need to run cmetindex first\n");
	exit(9);
      }
      if ((indexdb_rev = Indexdb_new_genome(&index1part,&index1interval,
					    modedir,fileroot,/*idx_filesuffix*/"metga",snps_root,
					    required_index1part,required_index1interval,
					    expand_offsets_p,offsetsstrm_access,positions_access,
					    sharedp)) == NULL) {
	fprintf(stderr,"Cannot find metga index file.  Need to run cmetindex first\n");
	exit(9);
      }

    } else if (mode == ATOI_STRANDED || mode == ATOI_NONSTRANDED) {
      if (user_atoidir == NULL) {
	modedir = snpsdir;
      } else {
	modedir = user_atoidir;
      }

      if ((indexdb_fwd = Indexdb_new_genome(&index1part,&index1interval,
					    modedir,fileroot,/*idx_filesuffix*/"a2iag",snps_root,
					    required_index1part,required_index1interval,
					    expand_offsets_p,offsetsstrm_access,positions_access,
					    sharedp)) == NULL) {
	fprintf(stderr,"Cannot find a2iag index file.  Need to run atoiindex first\n");
	exit(9);
      }
      if ((indexdb_rev = Indexdb_new_genome(&index1part,&index1interval,
					    modedir,fileroot,/*idx_filesuffix*/"a2itc",snps_root,
					    required_index1part,required_index1interval,
					    expand_offsets_p,offsetsstrm_access,positions_access,
					    sharedp)) == NULL) {
	fprintf(stderr,"Cannot find a2itc index file.  Need to run atoiindex first\n");
	exit(9);
      }

    } else if (mode == TTOC_STRANDED || mode == TTOC_NONSTRANDED) {
      if (user_atoidir == NULL) {
	modedir = snpsdir;
      } else {
	modedir = user_atoidir;
      }

      if ((indexdb_fwd = Indexdb_new_genome(&index1part,&index1interval,
					    modedir,fileroot,/*idx_filesuffix*/"a2itc",snps_root,
					    required_index1part,required_index1interval,
					    expand_offsets_p,offsetsstrm_access,positions_access,
					    sharedp)) == NULL) {
	fprintf(stderr,"Cannot find a2itc index file.  Need to run atoiindex first\n");
	exit(9);
      }
      if ((indexdb_rev = Indexdb_new_genome(&index1part,&index1interval,
					    modedir,fileroot,/*idx_filesuffix*/"a2iag",snps_root,
					    required_index1part,required_index1interval,
					    expand_offsets_p,offsetsstrm_access,positions_access,
					    sharedp)) == NULL) {
	fprintf(stderr,"Cannot find a2iag index file.  Need to run atoiindex first\n");
	exit(9);
      }

    } else {
      indexdb_fwd = Indexdb_new_genome(&index1part,&index1interval,
				       snpsdir,fileroot,/*idx_filesuffix*/"ref",snps_root,
				       required_index1part,required_index1interval,
				       expand_offsets_p,offsetsstrm_access,positions_access,
				       sharedp);
      if (indexdb_fwd == NULL) {
	fprintf(stderr,"Cannot find snps index file for %s in directory %s\n",snps_root,snpsdir);
	exit(9);
      }
      indexdb_rev = indexdb_fwd;
    }

    if (user_chrsubsetname != NULL) {
      if ((divno = Univ_IIT_find_one(chromosome_iit,user_chrsubsetname)) < 0) {
	fprintf(stderr,"Cannot find chrsubset %s in chromosome IIT file.  Ignoring.\n",user_chrsubsetname);
      } else {
	interval = Univ_IIT_interval(chromosome_iit,divno);
	chrsubset_start = Univinterval_low(interval);
	chrsubset_end = Univinterval_high(interval);
      }
    }

    FREE(genomesubdir);
    FREE(fileroot);
    FREE(dbroot);
#endif
  }

  if (splicing_file != NULL && genomecomp != NULL) {
    if (Genome_blocks(genomecomp) == NULL) {
      fprintf(stderr,"known splicing can be used only with compressed genome\n");
    } else {
      /* TODO: Handle case for observed distances */
      /* min_extra_end no longer used by gregion.c */
      min_extra_end = shortsplicedist;

      splicing_divint_crosstable = Univ_IIT_divint_crosstable(chromosome_iit,splicing_iit);
      if ((donor_typeint = IIT_typeint(splicing_iit,"donor")) >= 0 && 
	  (acceptor_typeint = IIT_typeint(splicing_iit,"acceptor")) >= 0) {
	fprintf(stderr,"found donor and acceptor tags, so treating as splicesites file\n");
	splicestringpool = Splicestringpool_new();
	splicesites = Splicetrie_retrieve_via_splicesites(&distances_observed_p,&splicetypes,&splicedists,
							  &splicestrings,&splicefrags_ref,&splicefrags_alt,
							  &nsplicesites,splicing_iit,splicing_divint_crosstable,
							  donor_typeint,acceptor_typeint,chromosome_iit,
							  genomecomp,genomecomp_alt/*can be NULL*/,shortsplicedist,
							  splicestringpool);
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
	splicesites = Splicetrie_retrieve_via_introns(&splicetypes,&splicedists,
						      &splicestrings,&splicefrags_ref,&splicefrags_alt,
						      &nsplicesites,splicing_iit,splicing_divint_crosstable,
						      chromosome_iit,genomecomp,genomecomp_alt/*can be NULL*/,
						      splicestringpool);
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
    }

    fprintf(stderr,"done\n");
  }

  if (user_pairalign_p == true) {
    /* Creation of genomebits/genomecomp and initialization done within single_thread() for each input sequence */

  } else if (usersegment != NULL) {
    Genome_user_setup(genomecomp_blocks);
    Genome_hr_user_setup(genomebits_blocks,/*query_unk_mismatch_p*/false,
			 /*genome_unk_mismatch_p*/true,/*mode*/STANDARD);
    Genome_sites_setup(genomecomp_blocks,/*snp_blocks*/NULL);
    Maxent_hr_setup(genomecomp_blocks,/*genomealt_blocks*/genomecomp_blocks);
#ifdef PMAP
    Oligoindex_pmap_setup(genomecomp);
#else
    Oligoindex_hr_setup(genomecomp_blocks,mode);
#endif

  } else if (genomecomp != NULL) {
    Genome_setup(genomecomp,genomecomp_alt/*can be NULL*/,mode,circular_typeint);
    if (genomebits_blocks == NULL) {
      fprintf(stderr,"This version of GMAP requires the genomebits128 file\n");
      exit(9);
    } else {
      Genome_hr_setup(genomebits_blocks,/*snp_blocks*/genomebits_alt ? Genome_blocks(genomebits_alt) : NULL,
		      /*query_unk_mismatch_p*/false,/*genome_unk_mismatch_p*/true,/*mode*/STANDARD);
    }
    Genome_sites_setup(Genome_blocks(genomecomp),/*snp_blocks*/genomecomp_alt ? Genome_blocks(genomecomp_alt) : NULL);
    Maxent_hr_setup(Genome_blocks(genomecomp),/*snp_blocks*/genomecomp_alt ? Genome_blocks(genomecomp_alt) : NULL);
#ifdef PMAP
    Alphabet_setup(alphabet,alphabet_size,index1part_aa);
    Oligoindex_pmap_setup(genomecomp);
    Oligop_setup(alphabet,alphabet_size,index1part_aa);
    Indexdb_setup(index1part_aa);
    Stage1_setup(index1part_aa,maxextension,maxtotallen_bound,min_extra_end,circular_typeint);
#else
    Oligoindex_hr_setup(Genome_blocks(genomecomp),mode);
    Oligo_setup(index1part);
    Indexdb_setup(index1part);
    Stage1_setup(index1part,maxextension,maxtotallen_bound,min_extra_end,circular_typeint);
#endif
  }

  Stage2_setup(/*splicingp*/novelsplicingp == true || knownsplicingp == true,cross_species_p,
	       suboptimal_score_start,suboptimal_score_end,sufflookback,nsufflookback,maxintronlen,mode,
	       /*snps_p*/genomecomp_alt ? true : false);
  Dynprog_single_setup(homopolymerp);
  Dynprog_genome_setup(novelsplicingp,splicing_iit,splicing_divint_crosstable,
		       donor_typeint,acceptor_typeint);
  Dynprog_end_setup(splicesites,splicetypes,splicedists,nsplicesites,
		    trieoffsets_obs,triecontents_obs,trieoffsets_max,triecontents_max);
  Pair_setup(trim_mismatch_score,trim_indel_score,gff3_separators_p,sam_insert_0M_p,
	     force_xs_direction_p,md_lowercase_variant_p,
	     /*snps_p*/genomecomp_alt ? true : false,
	     /*print_nsnpdiffs_p*/genomecomp_alt ? true : false,genomelength);
  Stage3_setup(/*splicingp*/novelsplicingp == true || knownsplicingp == true,novelsplicingp,
	       require_splicedir_p,splicing_iit,splicing_divint_crosstable,
	       donor_typeint,acceptor_typeint,
	       splicesites,min_intronlength,max_deletionlength,/*min_indel_end_matches*/6,
	       maxpeelback_distalmedial,nullgap,extramaterial_end,extramaterial_paired,
	       extraband_single,extraband_end,extraband_paired,
	       ngap,maxintronlen,/*output_sam_p*/printtype == SAM ? true : false,
	       homopolymerp,stage3debug);
  Splicetrie_setup(splicesites,splicefrags_ref,splicefrags_alt,
		   trieoffsets_obs,triecontents_obs,trieoffsets_max,triecontents_max,
		   /*snpp*/false,amb_closest_p,/*amb_clip_p*/true,/*min_shortend*/2);
  Output_setup(chromosome_iit,nofailsp,failsonlyp,quiet_if_excessive_p,maxpaths_report,
	       failedinput_root,quality_shift,
	       printtype,invertmode,wraplength,ngap,nointronlenp,sam_paired_p,cds_startpos,
	       fulllengthp,truncatep,strictp,checksump,genomecomp,usersegment,user_genomicseg,
	       dbversion,user_chrsubsetname,contig_iit,altstrain_iit,
	       /*chimeras_allowed_p*/chimera_margin > 0 ? true : false,
	       map_iit,map_divint_crosstable,map_exons_p,map_bothstrands_p,
	       nflanking,print_comment_p,sam_read_group_id);

#ifdef USE_MPI
  if (myid == 0) {
    Outbuffer_setup(argc,argv,optind,chromosome_iit,any_circular_p,
		    nworkers,orderedp,quiet_if_excessive_p,
		    printtype,usersegment,sam_headers_p,sam_read_group_id,sam_read_group_name,
		    sam_read_group_library,sam_read_group_platform,
		    appendp,/*output_file*/NULL,split_output_root,failedinput_root);
    outbuffer = Outbuffer_new(output_buffer_size,/*nread*/0);
    /* Inbuffer_set_outbuffer(inbuffer,outbuffer); */

    fprintf(stderr,"Starting alignment\n");
    stopwatch = Stopwatch_new();
    Stopwatch_start(stopwatch);
  }
#else
  Outbuffer_setup(argc,argv,optind,chromosome_iit,any_circular_p,
		  nworkers,orderedp,quiet_if_excessive_p,
		  printtype,usersegment,sam_headers_p,sam_read_group_id,sam_read_group_name,
		  sam_read_group_library,sam_read_group_platform,
		  appendp,/*output_file*/NULL,split_output_root,failedinput_root);
  outbuffer = Outbuffer_new(output_buffer_size,nread);
  Inbuffer_set_outbuffer(inbuffer,outbuffer);

  fprintf(stderr,"Starting alignment\n");
  stopwatch = Stopwatch_new();
  Stopwatch_start(stopwatch);
#endif


#ifdef USE_MPI
  /* MPI version */
  if (myid == 0) {
#ifdef WORKER_DETACH
    pthread_attr_init(&thread_attr_detach);
    if ((ret = pthread_attr_setdetachstate(&thread_attr_detach,PTHREAD_CREATE_DETACHED)) != 0) {
      fprintf(stderr,"ERROR: pthread_attr_setdetachstate %d\n",ret);
      exit(1);
    }
#endif
    pthread_attr_init(&thread_attr_join);
    if ((ret = pthread_attr_setdetachstate(&thread_attr_join,PTHREAD_CREATE_JOINABLE)) != 0) {
      fprintf(stderr,"ERROR: pthread_attr_setdetachstate %d\n",ret);
      exit(1);
    }

    Except_init_pthread();
    /* pthread_key_create(&global_request_key,NULL); */

    if (orderedp == true) {
      pthread_create(&output_thread_id,&thread_attr_join,Outbuffer_thread_ordered,
		     (void *) outbuffer);
    } else {
      pthread_create(&output_thread_id,&thread_attr_join,Outbuffer_thread_anyorder,
		     (void *) outbuffer);
    }

    Outbuffer_mpi_process(outbuffer,/*n_worker_procs*/nprocs - 1,part_modulus,part_interval);
    pthread_join(output_thread_id,NULL);

    /* pthread_key_delete(global_request_key); */
    /* Except_term_pthread(); */

  } else {
    worker_mpi_process(/*worker_id*/myid,inbuffer);
  }

#elif !defined(HAVE_PTHREAD)
  /* Serial version */
  single_thread();

#else
  /* Pthreads version */
  if (nworkers == 0) {
    single_thread();
    
  } else if (multiple_sequences_p == false) {
    single_thread();
    
  } else {
#ifdef WORKER_DETACH
    pthread_attr_init(&thread_attr_detach);
    if ((ret = pthread_attr_setdetachstate(&thread_attr_detach,PTHREAD_CREATE_DETACHED)) != 0) {
      fprintf(stderr,"ERROR: pthread_attr_setdetachstate %d\n",ret);
      exit(1);
    }
#endif
    pthread_attr_init(&thread_attr_join);
    if ((ret = pthread_attr_setdetachstate(&thread_attr_join,PTHREAD_CREATE_JOINABLE)) != 0) {
      fprintf(stderr,"ERROR: pthread_attr_setdetachstate %d\n",ret);
      exit(1);
    }
    
    worker_thread_ids = (pthread_t *) CALLOC(nworkers,sizeof(pthread_t));
    Except_init_pthread();
    pthread_key_create(&global_request_key,NULL);

    if (orderedp == true) {
      pthread_create(&output_thread_id,&thread_attr_join,Outbuffer_thread_ordered,
		     (void *) outbuffer);
    } else {
      pthread_create(&output_thread_id,&thread_attr_join,Outbuffer_thread_anyorder,
		     (void *) outbuffer);
    }

    for (i = 0; i < nworkers; i++) {
#ifdef WORKER_DETACH
      pthread_create(&(worker_thread_ids[i]),&thread_attr_detach,worker_thread,(void *) NULL);
#else
      /* Need to have worker threads finish before we call Inbuffer_free() */
      pthread_create(&(worker_thread_ids[i]),&thread_attr_join,worker_thread,(void *) NULL);
#endif
    }
    
    pthread_join(output_thread_id,NULL);
    for (i = 0; i < nworkers; i++) {
      pthread_join(worker_thread_ids[i],NULL);
    }

    pthread_key_delete(global_request_key);
    /* Do not delete global_except_key, because worker threads might still need it */
    /* Except_term_pthread(); */

    FREE(worker_thread_ids);

  }
#endif /* HAVE_PTHREAD */


#ifdef USE_MPI
  if (myid == 0) {
    runtime = Stopwatch_stop(stopwatch);
    Stopwatch_free(&stopwatch);

    nread = Outbuffer_nread(outbuffer);
    nbeyond = Outbuffer_nbeyond(outbuffer);
    fprintf(stderr,"Processed %u queries in %.2f seconds (%.2f queries/sec)\n",
	    nread-nbeyond,runtime,(double) nread/runtime);

    Outbuffer_free(&outbuffer);
    Inbuffer_free(&inbuffer);	/* Also closes inputs */
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
  Inbuffer_free(&inbuffer);	/* Also closes inputs */

  Outbuffer_close_files();
#endif

#ifdef PMAP
  Backtranslation_term();
#endif
  Dynprog_term();


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
  }

  if (splicing_iit != NULL) {
    FREE(splicing_divint_crosstable);
    IIT_free(&splicing_iit);
  }


#ifdef PMAP
 if (indexdb_rev != NULL) {
    Indexdb_free(&indexdb_rev);
  }
  if (indexdb_fwd != NULL) {
    Indexdb_free(&indexdb_fwd);
  }
#else
  if (indexdb_rev != indexdb_fwd) {
    Indexdb_free(&indexdb_rev);
  }
  if (indexdb_fwd != NULL) {
    Indexdb_free(&indexdb_fwd);
  }
#endif
  if (dbversion != NULL) {
    FREE(dbversion);
  }
  if (altstrain_iit != NULL) {
    IIT_free(&altstrain_iit);
  }
  if (genomecomp_alt != NULL) {
    Genome_free(&genomecomp_alt);
  }
  if (user_pairalign_p == true) {
    /* genomecomp_blocks freed within single_thread */
  } else if (usersegment != NULL) {
    FREE(genomecomp_blocks);
    FREE(genomebits_blocks);
  } else if (genomecomp != NULL) {
    Genome_free(&genomecomp);
  }
  if (genomebits != NULL) {
    Genome_free(&genomebits);
  }

  if (map_iit != NULL) {
    IIT_free(&map_iit);
  }
  if (contig_iit != NULL) {
    Univ_IIT_free(&contig_iit);
  }
  if (circularp != NULL) {
    FREE(circularp);
  }
  if (chromosome_iit != NULL) {
    Univ_IIT_free(&chromosome_iit);
  }

  if (user_selfalign_p == true) {
    /* Do not free usersegment */
  } else if (usersegment != NULL) {
    Sequence_free(&usersegment);
  }

  Outbuffer_cleanup();

  Access_controlled_cleanup();

#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);	/* Make sure all processes have cleaned up */
  MPI_Finalize();
#endif

  return 0;
}


static void
print_program_usage () {
#ifdef PMAP
    fprintf(stdout,"\
Usage: pmap [OPTIONS...] <FASTA files...>, or\n\
       cat <FASTA files...> | pmap [OPTIONS...]\n\
");
#else
    fprintf(stdout,"\
Usage: gmap [OPTIONS...] <FASTA files...>, or\n\
       cat <FASTA files...> | gmap [OPTIONS...]\n\
");
#endif
    fprintf(stdout,"\n");

    fprintf(stdout,"Input options (must include -d or -g)\n");
    fprintf(stdout,"\
  -D, --dir=directory            Genome directory.  Default (as specified by --with-gmapdb to the configure program) is\n \
                                   %s\n\
",GMAPDB);
    fprintf(stdout,"\
  -d, --db=STRING                Genome database.  If argument is '?' (with\n\
                                   the quotes), this command lists available databases.\n\
");
    fprintf(stdout,"\n");

#ifdef PMAP
    fprintf(stdout,"\
  -a, --alphabet=STRING          Alphabet to use in PMAP genome database\n\
                                   (allowed values in order of preference: 20, 15a, 12a).\n\
                                   If not specified, the program will find the first available\n\
                                   alphabet in the genome database in preference order\n\
");
#endif

    fprintf(stdout,"\
  -k, --kmer=INT                 kmer size to use in genome database (allowed values: 16 or less).\n\
                                   If not specified, the program will find the highest available\n\
                                   kmer size in the genome database\n\
  --sampling=INT                 Sampling to use in genome database.  If not specified, the program\n\
                                   will find the smallest available sampling value in the genome database\n\
                                   within selected k-mer size\n\
  -G, --genomefull               Use full genome (all ASCII chars allowed;\n\
                                   built explicitly during setup), not\n\
                                   compressed version\n\
  -g, --gseg=filename            User-supplied genomic segment\n\
  -1, --selfalign                Align one sequence against itself in FASTA format via stdin\n\
                                   (Useful for getting protein translation of a nucleotide sequence)\n\
  -2, --pairalign                Align two sequences in FASTA format via stdin, first one being\n\
                                   genomic and second one being cDNA\n\
  --cmdline=STRING,STRING        Align these two sequences provided on the command line,\n\
                                   first one being genomic and second one being cDNA\n\
  -q, --part=INT/INT             Process only the i-th out of every n sequences\n\
                                   e.g., 0/100 or 99/100 (useful for distributing jobs\n\
                                   to a computer farm).\n\
");
    fprintf(stdout,"\
  --input-buffer-size=INT        Size of input buffer (program reads this many sequences\n\
                                   at a time for efficiency) (default %d)\n\
",inbuffer_nspaces);
    fprintf(stdout,"\n");

    fprintf(stdout,"Computation options\n");
#ifdef HAVE_MMAP
    fprintf(stdout,"\
  -B, --batch=INT                Batch mode (default = 2)\n\
                                 Mode     Offsets       Positions       Genome\n\
                                   0      see note      mmap            mmap\n\
                                   1      see note      mmap & preload  mmap\n\
                      (default)    2      see note      mmap & preload  mmap & preload\n\
                                   3      see note      allocate        mmap & preload\n\
                                   4      see note      allocate        allocate\n\
                                   5      expand        allocate        allocate\n\
                           Note: For a single sequence, all data structures use mmap\n\
                           If mmap not available and allocate not chosen, then will use fileio (very slow)\n\
");
#else
    fprintf(stdout,"\
  -B, --batch=INT                Batch mode (default = 4, modes 0-3 disallowed because program configured without mmap)\n\
                                 Mode     Offsets       Positions       Genome\n\
                      (default)    4      see note      allocate        allocate\n\
                                   5      expand        allocate        allocate\n \
");
#endif
    fprintf(stdout,"\
                       Note about --batch and offsets: Expansion of offsets can be controlled\n\
                       independently by the --expand-offsets flag.  The --batch=5 option is equivalent\n\
                       to --batch=4 plus --expand-offsets=1\n\
\n\
  --expand-offsets=INT           Whether to expand the genomic offsets index\n\
                                   Values: 0 (no, default), or 1 (yes).\n\
                                   Expansion gives faster alignment, but requires more memory\n\
");

    fprintf(stdout,"\
  --nosplicing                   Turns off splicing (useful for aligning genomic sequences\n\
                                   onto a genome)\n\
");
    fprintf(stdout,"\
  --min-intronlength=INT         Min length for one internal intron (default %d).  Below this size,\n\
                                   a genomic gap will be considered a deletion rather than an intron.\n\
",min_intronlength);
    fprintf(stdout,"\
  -K, --intronlength=INT         Max length for one internal intron (default %d)\n\
",maxintronlen);
    fprintf(stdout,"\
  -w, --localsplicedist=INT      Max length for known splice sites at ends of sequence\n\
                                   (default %d)\n\
",shortsplicedist);
    fprintf(stdout,"\
  -L, --totallength=INT          Max total intron length (default %d)\n\
",maxtotallen_bound);
    fprintf(stdout,"\
  -x, --chimera-margin=INT       Amount of unaligned sequence that triggers\n\
                                   search for the remaining sequence (default %d).\n\
                                   Enables alignment of chimeric reads, and may help\n\
                                   with some non-chimeric reads.  To turn off, set to\n\
                                   zero.\n\
",chimera_margin);
    fprintf(stdout,"\
  --no-chimeras                  Turns off finding of chimeras.  Same effect as --chimera-margin=0\n\
");

#if 0
    fprintf(stdout,"\
  -w, --reference=filename       Reference cDNA sequence for relative alignment\n\
");
#endif

#ifdef HAVE_PTHREAD
    fprintf(stdout,"\
  -t, --nthreads=INT             Number of worker threads\n\
");
#else
  fprintf(stdout,"\
  -t, --nthreads=INT             Number of worker threads.  Flag is ignored in this version of GMAP, which has pthreads disabled\n\
");
#endif
    fprintf(stdout,"\
  -c, --chrsubset=string         Limit search to given chromosome\n\
  -z, --direction=STRING         cDNA direction (sense_force, antisense_force,\n\
                                   sense_filter, antisense_filter,or auto (default))\n\
");
    fprintf(stdout,"\
  -H, --trimendexons=INT         Trim end exons with fewer than given number of matches\n\
                                   (in nt, default %d)\n\
",minendexon);
    fprintf(stdout,"\
  --canonical-mode=INT           Reward for canonical and semi-canonical introns\n\
                                   0=low reward, 1=high reward (default), 2=low reward for\n\
                                   high-identity sequences and high reward otherwise\n\
  --cross-species                Use a more sensitive search for canonical splicing, which helps especially\n\
                                   for cross-species alignments and other difficult cases\n\
  --allow-close-indels=INT       Allow an insertion and deletion close to each other\n\
                                   (0=no, 1=yes (default), 2=only for high-quality alignments)\n\
");
    fprintf(stdout,"\
  --microexon-spliceprob=FLOAT   Allow microexons only if one of the splice site probabilities is\n\
                                   greater than this value (default %.2f)\n\
",microexon_spliceprob);

#if 0
    fprintf(stdout,"\
  --homopolymer                  In dynamic programming, favor indels in regions of homopolymers,\n\
                                   e.g., AAAAAA.  Useful for some platforms, such as Pacific Biosciences\n\
");
#endif

#ifndef PMAP
    fprintf(stdout,"\
  --cmetdir=STRING               Directory for methylcytosine index files (created using cmetindex)\n\
                                   (default is location of genome index files specified using -D, -V, and -d)\n\
  --atoidir=STRING               Directory for A-to-I RNA editing index files (created using atoiindex)\n\
                                   (default is location of genome index files specified using -D, -V, and -d)\n\
  --mode=STRING                  Alignment mode: standard (default), cmet-stranded, cmet-nonstranded,\n\
                                    atoi-stranded, atoi-nonstranded, ttoc-stranded, or ttoc-nonstranded.\n\
                                    Non-standard modes requires you to have previously run the cmetindex\n\
                                    or atoiindex programs (which also cover the ttoc modes) on the genome\n\
");
#endif

#if 0
    /* Causes seg faults, so do not advertise */
    fprintf(stdout,"\
  -s, --splicing=STRING          Look for splicing involving known sites\n\
                                   (in <STRING>.iit)\n\
");
#endif

#ifndef PMAP
    fprintf(stdout,"\
  -p, --prunelevel               Pruning level: 0=no pruning (default), 1=poor seqs,\n\
                                   2=repetitive seqs, 3=poor and repetitive\n\
");
#endif
    fprintf(stdout,"\n");

    fprintf(stdout,"\
Output types\n\
  -S, --summary                  Show summary of alignments only\n\
  -A, --align                    Show alignments\n\
  -3, --continuous               Show alignment in three continuous lines\n\
  -4, --continuous-by-exon       Show alignment in three lines per exon\n\
  -Z, --compress                 Print output in compressed format\n\
  -E, --exons=STRING             Print exons (\"cdna\" or \"genomic\")\n\
");

#ifdef PMAP    
    fprintf(stdout,"\
  -P, --protein_gen              Print protein sequence (genomic)\n\
  -Q, --nucleotide               Print inferred nucleotide sequence from protein\n\
");
#else
    fprintf(stdout,"\
  -P, --protein_dna              Print protein sequence (cDNA)\n\
  -Q, --protein_gen              Print protein sequence (genomic)\n\
");
#endif

#ifdef PMAP
    fprintf(stdout,"\
  -f, --format=INT               Other format for output (also note the -A and -S options\n\
                                   and other options listed under Output types):\n\
                                   psl_pro (or 0) = PSL format in protein coords,\n\
                                   psl_nt (or 1) = PSL format in nucleotide coords,\n\
                                   gff3_gene (or 2) = GFF3 gene format,\n\
                                   gff3_match_cdna (or 3) = GFF3 cDNA_match format,\n\
                                   gff3_match_est (or 4) = GFF3 EST_match format,\n\
                                   map_exons (or 7) = IIT FASTA exon map format,\n\
                                   map_ranges (or 8) = IIT FASTA range map format,\n\
                                   coords (or 9) = coords in table format\n\
");
#else
    fprintf(stdout,"\
  -f, --format=INT               Other format for output (also note the -A and -S options\n\
                                   and other options listed under Output types):\n\
                                   psl (or 1) = PSL (BLAT) format,\n\
                                   gff3_gene (or 2) = GFF3 gene format,\n\
                                   gff3_match_cdna (or 3) = GFF3 cDNA_match format,\n\
                                   gff3_match_est (or 4) = GFF3 EST_match format,\n\
                                   splicesites (or 6) = splicesites output (for GSNAP splicing file),\n\
                                   introns = introns output (for GSNAP splicing file),\n\
                                   map_exons (or 7) = IIT FASTA exon map format,\n\
                                   map_ranges (or 8) = IIT FASTA range map format,\n\
                                   coords (or 9) = coords in table format,\n\
                                   sampe = SAM format (setting paired_read bit in flag),\n\
                                   samse = SAM format (without setting paired_read bit)\n\
");
#endif
    fprintf(stdout,"\n");

    fprintf(stdout,"\
Output options\n\
  -n, --npaths=INT               Maximum number of paths to show (default %d).  If set to 1, GMAP\n\
                                   will not report chimeric alignments, since those imply\n\
                                   two paths.  If you want a single alignment plus chimeric\n\
                                   alignments, then set this to be 0.\n\
",maxpaths_report);
    fprintf(stdout,"\
  --suboptimal-score=INT         Report only paths whose score is within this value of the\n\
                                   best path.  By default, if this option is not provided,\n\
                                   the program prints all paths found.\n\
  -O, --ordered                  Print output in same order as input (relevant\n\
                                   only if there is more than one worker thread)\n\
  -5, --md5                      Print MD5 checksum for each query sequence\n\
  -o, --chimera-overlap          Overlap to show, if any, at chimera breakpoint\n\
  --failsonly                    Print only failed alignments, those with no results\n\
  --nofails                      Exclude printing of failed alignments\n\
\n\
  -V, --snpsdir=STRING           Directory for SNPs index files (created using snpindex) (default is\n\
                                   location of genome index files specified using -D and -d)\n \
  -v, --use-snps=STRING          Use database containing known SNPs (in <STRING>.iit, built\n\
                                   previously using snpindex) for tolerance to SNPs\n\
");

  fprintf(stdout,"\
  --split-output=STRING          Basename for multiple-file output, separately for nomapping,\n\
                                   uniq, mult, (and chimera, if --chimera-margin is selected)\n\
  --failed-input=STRING          Print completely failed alignments as input FASTA or FASTQ format\n\
                                   to the given file.  If the --split-output flag is also given, this file\n\
                                   is generated in addition to the output in the .nomapping file.\n\
  --append-output                When --split-output or --failedinput is given, this flag will append output\n\
                                   to the existing files.  Otherwise, the default is to create new files.\n\
");
  fprintf(stdout,"\
  --output-buffer-size=INT       Buffer size, in queries, for output thread (default %d).  When the number\n\
                                   of results to be printed exceeds this size, the worker threads are halted\n\
                                   until the backlog is cleared\n\
",output_buffer_size);


#ifdef PMAP    
    fprintf(stdout,"\
  -Y, --tolerant                 Translates genome with corrections for frameshifts\n\
");
#else
    fprintf(stdout,"\
  -F, --fulllength               Assume full-length protein, starting with Met\n\
  -a, --cdsstart=INT             Translate codons from given nucleotide (1-based)\n\
  -T, --truncate                 Truncate alignment around full-length protein, Met to Stop\n\
                                 Implies -F flag.\n\
  -Y, --tolerant                 Translates cDNA with corrections for frameshifts\n\
");
#endif

    fprintf(stdout,"\n");

#ifndef PMAP
  fprintf(stdout,"Options for GFF3 output\n");
  fprintf(stdout,"\
  --gff3-add-separators=INT      Whether to add a ### separator after each query sequence\n\
                                   Values: 0 (no), 1 (yes, default)\n\
");
  fprintf(stdout,"\n");

  fprintf(stdout,"Options for SAM output\n");
  fprintf(stdout,"\
  --no-sam-headers               Do not print headers beginning with '@'\n\
  --sam-use-0M                   Insert 0M in CIGAR between adjacent insertions and deletions\n\
                                   Required by Picard, but can cause errors in other tools\n\
  --force-xs-dir                 For RNA-Seq alignments, disallows XS:A:? when the sense direction\n\
                                   is unclear, and replaces this value arbitrarily with XS:A:+.\n\
                                   May be useful for some programs, such as Cufflinks, that cannot\n\
                                   handle XS:A:?.  However, if you use this flag, the reported value\n\
                                   of XS:A:+ in these cases will not be meaningful.\n\
  --md-lowercase-snp             In MD string, when known SNPs are given by the -v flag,\n\
                                   prints difference nucleotides as lower-case when they,\n\
                                   differ from reference but match a known alternate allele\n\
  --action-if-cigar-error        Action to take if there is a disagreement between CIGAR length and sequence length\n\
                                   Allowed values: ignore, warning (default), abort\n\
  --read-group-id=STRING         Value to put into read-group id (RG-ID) field\n\
  --read-group-name=STRING       Value to put into read-group name (RG-SM) field\n\
  --read-group-library=STRING    Value to put into read-group library (RG-LB) field\n\
  --read-group-platform=STRING   Value to put into read-group library (RG-PL) field\n\
");
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
                                 Or you can specify the print shift with this flag:\n\
  -j, --quality-print-shift=INT  Shift FASTQ quality scores by this amount in output\n\
                                   (default is 0 for sanger protocol; to change Illumina input\n\
                                   to Sanger output, select -31)\n\
");
#endif

    fprintf(stdout,"\
External map file options\n\
  -M, --mapdir=directory         Map directory\n\
  -m, --map=iitfile              Map file.  If argument is '?' (with the quotes),\n\
                                   this lists available map files.\n\
  -e, --mapexons                 Map each exon separately\n\
  -b, --mapboth                  Report hits from both strands of genome\n\
  -u, --flanking=INT             Show flanking hits (default 0)\n\
  --print-comment                Show comment line for each hit\n\
");
    fprintf(stdout,"\n");

    fprintf(stdout,"\
Alignment output options\n\
  -N, --nolengths                No intron lengths in alignment\n\
  -I, --invertmode=INT           Mode for alignments to genomic (-) strand:\n\
                                   0=Don't invert the cDNA (default)\n\
                                   1=Invert cDNA and print genomic (-) strand\n\
                                   2=Invert cDNA and print genomic (+) strand\n\
");
    fprintf(stdout,"\
  -i, --introngap=INT            Nucleotides to show on each end of intron (default %d)\n\
",ngap);
    fprintf(stdout,"\
  -l, --wraplength=INT           Wrap length for alignment (default %d)\n\
",wraplength);
    fprintf(stdout,"\n");

    fprintf(stdout,"\
Filtering output options\n\
  --min-trimmed-coverage=FLOAT   Do not print alignments with trimmed coverage less\n\
                                   this value (default=0.0, which means no filtering)\n\
                                   Note that chimeric alignments will be output regardless\n\
                                   of this filter\n\
  --min-identity=FLOAT           Do not print alignments with identity less\n\
                                   this value (default=0.0, which means no filtering)\n\
                                   Note that chimeric alignments will be output regardless\n\
                                   of this filter\n\
\n\
Help options\n\
  --check                        Check compiler assumptions\n\
  --version                      Show version\n\
  --help                         Show this help message\n\
");

    return;
}
