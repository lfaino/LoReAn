/*  A Program for Comparing a DNA sequence against a database of
    DNA or cDNA sequences (DDS, DNA-DNA Search)
    
    Proper attribution of the author as the source of the software would
    be appreciated:
    Huang, X.  (1996)
    Fast Comparison of a DNA Sequence with a Protein Sequence Database.
    Microbial & Comparative Genomics, 1(4): 281-291.
    
    Acknowledgments
    I thank the following people for discussions and suggestions:
    Mark Adams, Tony Kerlavage, Brendan Loftus, Steve Rounsley,
    Granger Sutton, and Jinghui Zhang.
    The integration of DDS with GAP2 was performaed at TIGR.
    
    The DDS program compares a DNA sequence to a DNA or cDNA database.
    The DDS enhances the existing methods by addressing the problems
    of gaps. DDS computes high-scoring chains of segment pairs, where
    there can be an intervening DNA sequence between adjacent segment
    pairs in a chain.
    
    We address the problem of introns by giving a constant penalty
    to a long gap.
    
    There are two gap length parameters qgaplen and lgaplen,
    qgaplen for the query and lgaplen for database sequence.
    Any query region of length <= qgaplen between two
    segment pairs is given a linear penalty, and any query
    region of length > qgaplen is penalized as a gap of length qgaplen.
    The same is used for penalizing a region of the database sequence
    between two segment pairs with lgaplen being the gap length parameter.
    If the query contains introns, then a small value should be used
    for qgaplen. If the database contains cDNA sequences, a large
    value should be used for lgaplen since cDNA sequences
    don't contain introns. This scheme does not heavily
    penalize an intron region between two segment pairs.
    The qgaplen parameter is specified using the -q option,
    and the lgaplen parameter using the -l option, the letter l.
    
    A chain of segment pairs is reported if the score of each segment
    pair is >= segment pair score cutoff, and the score of the chain
    is >= chain score cutoff. In addition, the percent identity
    of the chain must be >= percent identity cutoff, and the percent
    similarity of the chain must be >= percent similarity cutoff.
    Indels and substitutions involving N's count as similarities.
    
    The sensitivity of the program depends on the word size parameter W.
    Increased sensitivity comes at the expense of decreased speed.
    A hit is an exact occurrence of the DNA word in the query sequence.
    Each hit is extended in both directions.
    The extesion stops if the score drops by more than the D distance.
    
    For a segment pair s, let nstart(s) and nend(s) denote
    the starting and ending positions of the DNA segment in
    the DNA sequence, let astart(s) and aend(s) denote
    the starting and ending positions of the protein segment in
    the protein sequence, and let score(s) denote the score of s.
    The first antidiagonal of a segment pair s is defined to
    be antis(s) = nstart(s) + astart(s), and
    the last antidiagonal of s is defined to be
    antid(s) = nend(s) + aend(s).
    
    A chain of segment pairs is a list of segment pairs in increasing
    order of their last antidiagonal such that each segment pair
    is not far from its predecessor and adjacent segment pairs
    do not have a large overlap.
    Specifically, any two adjacent segment pairs s and s' in the list
    satisfy the requirement:
    
    antis(s')  -  antid(s)  <  A, astart(s') - aend(s) > -B,
    and nstart(s') - nend(s) > - B.
    
    for some nonnegative integers A and B.
    Here A is called the maximum number of antidiagonals between s and s',
    and B is the AAOVER parameter in the program.
    The value for A can be specified at the command line.
    A long intron requres use of a large value for A.
    
    The program shows at most C number of chains of segments.
    If there are extra chains to be reported, the program only
    prints out the chain headings without showing the segment alignments.
    
    The DDS program is written in C and runs under Unix systems.
    We think that the program is portable to many machines.
    
    To compare a DNA sequence with a cDNA database, use a command of form:
    
    Usage: dds DNA_Seq DNA_Database [options]
    
    DNA_Seq   file of one query DNA in FASTA format
    DNA_Database  file of DNA database in FASTA format
    Options (default values):
    -a  N  specify max number N of antidiagonals between segments (50000)
    -b|btab  emit btab output
    -bh      emit btab output with headers
    -c  N  specify max number N of chain alignments reported (300)
    -d  N  specify distance N for extension (30)
    -f  N  specify final chain score cutoff N (100)
    -i  N  specify initial segment score cutoff N (35)
    -l  N  specify gap length N for database sequence (1000)
    -m  N  specify match score N > 0 (2)
    -o  N  specify percent similarity cutoff N <= 100 (65)
    -p  N  specify percent identity cutoff N <= 100 (60)
    -q  N  specify gap length N for query sequence (20)
    -s  N  specify mismatch score N < 0 (-3)
    -w  N  specify nucleotide word size N >= 11 (11)
    
    A sample DNA file:
    >query DNA
    AACTGCTGCATTGCATCCTTGACCGAGAAATTCCCCTGAAAACTCAAACA
    ATGGAGAAATATCACTGAAAAAAGCATACCAAACTCATAAGAATCAATTC
    TATTCGGCCATTCAGTGCTTTTCAAATCCTCATTGTATCGTAGAGAATTG
    CCAAATATCTGCCTGATTCGTTCGACTAAATCCACATTTCTTGGCAATTC
    TTTACGAATCGGTTGAAACTTTCGACGCCAATCATTTTGATTTTGAAATC
    GATGTCTTTCGATTATTACTGAGGAAAATGGCCGGAAAACTCGACCGTGG
    CGAGTTACAAGGCGAGCCTGAAAGGAATTCAATTTTAATTTTTTTTTGTT
    GAATTTTCATAAACTTACAATTCGTAGAGCATTTTGTCGTGGAAATGGAG
    CCTTTTTGAGTACGACATGAACTCCTAAATTATAAGTGGCGGGAAATATT
    CTCGGCAGGAAGCGTTGGAAAATTGGTAGTCGTCCACCTTTTGCAACTAA
    CTAAAAAATTTCAATTTAACGATTATAAAAAAGAACGATTTTGAAGATAC
    CTCATTTGCGATTCGATATGCTGCTTTTCCGAATCGTTTCATAGACATTC
    TCG
    A sample cDNA database file:
    >gi|275717|gb|M89216|M89216 CEL18G12 Caenorhabditis elegans cDNA clone cm18g12 5'.
    CAAGTCTNTCAAGAAGAATGTCTATAAACGATTCGGAAAAGCAGCATATCGAATCGCAAATNAGTTAGTT
    GCAAAAGGTGGACGACTACCAATTTTCCAACGCTTCCTGCCGAGAATATTTCCCGCCACTTATAATTTAG
    GAGTTCATGTCGTACTCAAAAAGGCTCCATTTCCACGACAAAATGCTCTACGAATTGCTCGCCTTGTAAC
    TCGCCACGGTCGAGTTTTCCGGCCATTTTCCTCAGTAATAATCGAAAGACATCGATTTCAAAATCAAAAT
    GATTGGCGTCGAAAGTTTCAACCGATTCGTAAAGAATTGCCAAGAAATGTGGGATTTAGTCGAACGAATC
    AGGCAGATATTTGGCAATTCTCTACGANACAATAGGGATTTGAAANGCACCTAATGGCCNATAGANTNAT
    CCTTATGAGTTTGGGGGATTCCCGGNCAGGGTGCAT
    >gi|275718|gb|M89217|M89217 CEL18G2 Caenorhabditis elegans cDNA clone cm18g2 5'.
    AAATAACANGTTCGGTGGCCAATAAGGTCTGTCTCATCGTTATTGATGGATGGGGAGTTTCTGAAGATCC
    TTACGGTAACGCTATTCTCAACGCACAGACACCAGTTATGGACAAGCTGTGTTCGGGCAATTGGGCTCAA
    ATTGAGGCACATGGTCTTCATGTTGGTCTCCCAAANGGATTTATGGGAAATTCGGAAGTCGGACATTTGA
    ACATCGGAGCCGGACGNGTTATCTATCAAGACANTNTTCGTATTAATCTGGCAGTCAAGANCAACAAATT
    TGTGACTAATGAGAGCTTNGTGGATNCTTGCGATCGTGCTAAAAACGGAAATNGACGTCTTCATCTGGCC
    GGACTGGTTTCTACGGGGGTGTTCATTCNCATATTGANCACA
    >gi|275719|gb|M89218|M89218 CEL18G3 Caenorhabditis elegans cDNA clone cm18g3 5'.
    CGCACCACACAAGGTCTCTGTGACTGAATTAATGGAAGAAAACNGAGCAGTTGAAAGCTGAAGTTAAAGA
    TTTACAGCAAGAAATTGAAGAGATGCAGGACCAGTACCGTGAGAANGAAATCGAAGAATTCCGGGAGCTT
    CAGCANGAACTTGAGCTTAATGCAAAAAATTNTCGCGTTCTGCAGTTTAAGCTGAGAAAAACAGAAAGAA
    GTAGAGATCAAGCTGAAGCAGANAAAATGCACTCTGAGAAAAAACTAGATGAATACATGANCAGTTGCCC
    AGNGGCGGTANTGCCATCTATTAAATCAGACAGTGCAAAAGTGAAAGANCTTGAATATNAGATTCGAGTA
    GCAAAAGNGGGTTTCTGTT
*/

#include   <stdio.h>

#define  FCUTOFF     100	/* default value for fcutoff */       
#define  ICUTOFF      35	/* default value for icutoff */
#define  WORDSIZE     11   	/* default value for wordsize, must be >= 11 */
#define  DISTANCE     30        /* default value for distance */
#define  MINISEG  200000   	/* minimum value for segsize */
#define  ANTIDIS   50000   	/* default value for antidis */
#define  AAOVER       40   	/* max length of codon overlap between segments */
#define  MAXCHAINS   300   	/* default value for numchains */
#define  QGAPLEN      20   	/* default value for qgaplen */
#define  LGAPLEN    1000   	/* default value for lgaplen */
#define  MATCH         2   	/* default value for match */
#define  MISMAT       -3   	/* default value for mismh */
#define  SIMCUTOFF    65   	/* default value for simcutoff */
#define  IDNCUTOFF    60   	/* default value for idncutoff */

#define  DNASIZE       5	/* size of DNA alphabet: A, C, G, T, N */
#define  DNAS2  (DNASIZE * DNASIZE) /* DNASIZE square */
#define  DNAS3  (DNASIZE * DNAS2)   /* DNASIZE cube */
#define  DNAS4  (DNASIZE * DNAS3)   /* DNASIZE to the 4th power */
#define  MINF      -10000	/* negative infinity */
#define  ACCSIZE     300        // accession id buffer in fasta file
#define  FNAMESZ    1023        // longest output filename size       

// Version info
static float VersionS = 1.51f;
static char * BuildS = "$Revision: 1.17 $";

static int  wordsize;		/* size of word */
static int  diff8;		/* wordsize - 8 */
static int  distance;		/* stop extension if score drops by distance */
static int  maxchains; 		/* max number of chain alignments reported */
static int  numchains = -1; 		/* number of chain alignments */
static int  fcutoff; 		/* final cutoff on alignment score */       
static int  icutoff; 		/* initial cutoff on segment score */       
static int  antidis;		/* maximum no. of antidiagonals between segments */
static int  qgaplen;     /* gap length for constant gap penalty for query */
static int  lgaplen;     /* gap length for constant gap penalty for library */
static int  match;		/* match score */
static int  mismh;		/* mismatch score */
static int  simcutoff;		/* percent similarity cutoff */
static int  idncutoff;		/* percent identity cutoff */
static int  maxchno = 300000;   /* max number of chains reported */

static int  *frame1, *frame2, *frame3, *frame4;	/* 4 frames */
static int  *framept[4];	/* pointers to frames */
static int  *table;	/* lookup table for DNA */
static int  tsize;	/* table size */
static int  *list;   	/* lookup table for DNA */
static int  *cfm1;	/* one frame of cDNA sequence */
typedef struct MSP   {          
  int   dstart;		/* start position in codon sequence */
  int   dend;		/* end position in codon sequence */
  int   astart;		/* start position in cDNA sequence */
  int   aend;		/* end position in cDNA sequence */
  int   new;		/* 1, new; 0, old */
  int   score;		/* score of segment pair */
  int   total;		/* score of a list ending with the segment */
  int   pred;		/* predecessor of the segment in the list */
  int   first;		/* first segment in the list */
  int   antis;		/* start antidiagonal no. of the segment */
  int   antid;		/* end antidiagonal no. of the segment */
  int   diag;		/* diagonal no. of the segment */
  int   fnum;		/* frame number of DNA */
  int   *asp;		/* pointer to the start position in cfm */
  int   *dsp;		/* pointer to the srart position in frame */
  int   idn;		/* number of matches */
  int   sim;		/* number of similarities */
}  segtp, *segptr; /* type definition */
segptr segment; 		/* for holding segment pairs */
segptr *segstd;			/* a sorted list of segment pairs */
static int    segsize;		/* maximum size of segment */
static int    cseglen;		/* current size of segment */
static int    *diagonal;	/* diagonal flag */
typedef struct TUPT   {          
  int   first;		/* first segment in the chain */
  int   last;		/* last segment in the chain */
  int   total;		/* score of the chain */
  int   antid;		/* largest end antidiagonal seen */
}  tuptp, *tupptr; /* type definition */
tupptr tuplist;			/* a list of disjoint chains of segments */
static int    tupsize;		/* current size of tuplist */
static char   accn[300];		/* accession number */
static int    aaover;		/* aaover = -AAOVER, aaover can't be modified */
static int    fdst[AAOVER];	/* score of each 5' end of a segment */

static int    w[DNAS4+1][DNAS4+1]; /* score table for substitutions */
static int    da[DNAS4+1];	/* 1, no 'N'; 0, 'N' */
static int    di[128];		/* DNA letter to integer code */
static int    nncode;		/* DNASIZE - 1 */
static int    dnals2;		/* DNASIZE square */
static int    dnals3;		/* DNASIZE cube */
static int    dnals4;		/* DNASIZE raised to 4th power */
static int    dnalen;		/* length of DNA sequence */
static int    aalen;		/* maximum cDNA sequence length */
static int    inits;		/* initial segment score */

static char *name1, *name2;		/* names of sequence files    */
static int  q, r;       /* gap penalties */
static int  qpay;	/* constant-cost for long gap in query */
static int  lpay;	/* constant-cost for long gap in library sequence */
static int global_chain_num = 0; //bhaas added to track alignment number
static char match_header[ACCSIZE]; // mfs added to track db match header string
static char outfile[FNAMESZ];      // name of output file if not stdout
static int outfile_option = 0;     // whether output file option selected 
static int align = 1;              // whether align option selected 
static FILE * alignfp;             // alignment output file pointer

// btab option (mfs) 
static char btab_date[100];         // printable date  Mmm dd YYYY
static char btab_line[1024];        // line on which repeated btab info printed
static char * program = "";         // program name 
static char * type = "dds";         // blast type field
static int btab = 0;                // whether btab option selected 
static int btab_header = 0;         // whether btab option selected 
static FILE * btabfp;               // btab output file pointer
static char query_name[1024];       // buffer the Fasta query name
static char * headers[] = 
  {
    "Query name                      ",
    "Date         ",
    "Query length (bp)",
    "Program",
    "Database",
    "Accession",
    "dstart",
    "dstop",
    "astart",
    "astop",
    "Identity",
    "Similarity",
    "Bit score",
    "Chain number",
    "Exon number",
    "DB match",
    "?",
    "Orientation",
    "Match length (bp)",
    NULL
  };
static char * btab_formatS = 
"%s\t%s\t%d\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%.02lf\t%d\t%d\t%d\t%d\t%s\t%d\t%s\t%d";

// btab print macro (mfs)
#define BTAB_LINE()                                  \
{                                                    \
  fprintf(btabfp,                                    \
          btab_formatS,                              \
          query_name,                                \
          btab_date,                                 \
          dnalen,                                    \
          type,                                      \
          name2,                                     \
          accn,                                      \
          (ort == 1)? dstart+1 : dnalen - dstart,    \
          (ort == 1)? dend     : dnalen - dend + 1,  \
          astart+1,                                  \
          aend,                                      \
          ((double)100.0 * idn) / len,               \
          -1,                                        \
          score,                                     \
          global_chain_num,                          \
          exon_num,                                  \
          match_header,                              \
          -1,                                        \
          (ort == 1 ? "Plus" : "Minus"),             \
          psize                                      \
         );                                          \
  fprintf(btabfp,                                    \
          "\n"                                       \
         );                                          \
}

// debug option (mfs)
static int debug = 0;

void version()
{
  char f[32];
  int i, k=0;
  for (i = 0; i < strlen(BuildS) && i < 31; ++i)
    {
      if (isdigit((char)BuildS[i]) || BuildS[i] == '.')
	{
	  f[k++] = (char)BuildS[i];
	}
    }
  f[k] = '\0';
  fprintf(stderr,"%s  Version %1.02f (Build %s)\n", program, VersionS, f);
}

void usage()
{
  fprintf(stderr,"dds - Compare a DNA sequence to a DNA or cDNA database.\n\n");
  fprintf(stderr,"Usage: dds DNA_Seq DNA_Database [-btab] [-out <file>] [options]\n\n");
  fprintf(stderr,"  DNA_Seq       file of one query DNA in FASTA format\n");
  fprintf(stderr,"  DNA_Database  file of DNA database in FASTA format\n");
  fprintf(stderr,"  -b|btab       emit output in btab format\n");
  fprintf(stderr,"  -bh           emit output in btab format with column headers\n");
  fprintf(stderr,"  -out <file>   alignment report to <filename>.dds \n");
  fprintf(stderr,"                btab report to <filename>.dds.btab \n");
  fprintf(stderr,"  Options (default values):\n");
  fprintf(stderr,"  -a  N  specify max number N of antidiagonals between segments (50000)\n");
  fprintf(stderr,"  -c  N  specify max number N of chain alignments reported (300)\n");
  fprintf(stderr,"  -d  N  specify distance N for extension (30)\n");
  fprintf(stderr,"  -f  N  specify final chain score cutoff N (100)\n");
  fprintf(stderr,"  -i  N  specify initial segment score cutoff N (35)\n");
  fprintf(stderr,"  -l  N  specify gap length N for database sequence (1000)\n");
  fprintf(stderr,"  -m  N  specify match score N > 0 (2)\n");
  fprintf(stderr,"  -o  N  specify percent similarity cutoff N <= 100 (65)\n");
  fprintf(stderr,"  -p  N  specify percent identity cutoff N <= 100 (60)\n");
  fprintf(stderr,"  -q  N  specify gap length N for query sequence (20)\n");
  fprintf(stderr,"  -s  N  specify mismatch score N < 0 (-3)\n");
  fprintf(stderr,"  -w  N  specify nucleotide word size N >= 11 (11)\n");
  fprintf(stderr,"\n");
  fprintf(stderr,"See also:\n");
  fprintf(stderr,"  dps ext filter gap2 nap show\n");
  fprintf(stderr,"\n");
  fprintf(stderr,"Huang, X.  Fast Comparison of a DNA Sequence with a Protein Sequence Database\n");
  fprintf(stderr,"  Microbial & Comparative Genomics, 1(4): 281-291 (1996).\n");
  fprintf(stderr,"\n");
}

main(argc, argv) int argc; char *argv[];
{ char *D,  *P;				/* Storing two sequences */
 int  *DN, *PN;			/* two sequences of integer codes */
 int  symbol;				/* The next character	      */
 FILE *Bp, *Ap, *ckopen();
 char *ckalloc();			/* space-allocating function  */
 register int i, j, k;
 int  size;				/* size of alphabet */
 int  temp;
 char *dhead;				/* the name of the genomic DNA */
 long  cur_time;
 
 // static initializations
 wordsize = WORDSIZE;
 distance = DISTANCE;
 maxchains = MAXCHAINS;
 fcutoff = FCUTOFF;
 icutoff = ICUTOFF;
 antidis = ANTIDIS;
 aaover = - AAOVER;	/* don't change aaover */
 qgaplen = QGAPLEN;
 lgaplen = LGAPLEN;
 match = MATCH;
 mismh = MISMAT;
 simcutoff = SIMCUTOFF;
 idncutoff = IDNCUTOFF;
 
 // initialize operational parameters
 btab = 0;
 btab_header = 0;
 align = 1;
 
 // parse command line
 program = argv[0];
 
 if ( argc < 2 )
   { 
     usage();  // exit with usage error
     exit(1);
   }
 name1 = (argc > 1)? argv[1] : ""; 
 name2 = (argc > 2)? argv[2] : ""; 
 
 // parse options (mfs)
 k = 0;
 while (k++ < argc-1)
   {
     
     if (!strcmp(argv[k], "-help")  ||  !strcmp(argv[k], "-h"))
       {
	 usage();
	 exit(1);
       }
     if (!strcmp(argv[k], "-version")  ||  !strcmp(argv[k], "-V"))
       {
	 version();
	 exit(1);
       }
     if (!strcmp(argv[k], "-depend"))
       {
	 // No dependencies to print on stderr
	 exit(1);
       }
     if (!strcmp(argv[k], "-debug"))
       {
	 fprintf(stderr, "name1 = %s\n", name1); 
	 fprintf(stderr, "name2 = %s\n", name2); 
	 fprintf(stderr, "btab = %d\n", btab); 
	 fprintf(stderr, "antidis = %d\n", antidis); 
	 fprintf(stderr, "maxchains = %d\n", maxchains); 
	 fprintf(stderr, "distance = %d\n", distance); 
	 fprintf(stderr, "fcutoff = %d\n", fcutoff);
	 fprintf(stderr, "icutoff = %d\n", icutoff); 
	 fprintf(stderr, "lgaplen = %d\n", lgaplen); 
	 fprintf(stderr, "match = %d\n", match); 
	 fprintf(stderr, "simcutoff = %d\n", simcutoff); 
	 fprintf(stderr, "idncutoff = %d\n", idncutoff); 
	 fprintf(stderr, "qgaplen = %d\n", qgaplen); 
	 fprintf(stderr, "mismh = %d\n", mismh); 
	 fprintf(stderr, "wordsize = %d\n", wordsize); 
	 
	 debug = (k+1 < argc)? atoi(argv[k+1]) : 0;
	 exit(1);
       }
     
     if (!strcmp(argv[k], "-out"))
       {
	 char * name = (k+1 < argc)? argv[k+1] : "";
	 strncpy(outfile, name, FNAMESZ-5);
	 outfile_option = 1;
	 continue;
       }
     if (!strcmp(argv[k], "-bh"))
       {
	 btab = 1; 
	 btab_header = 1;
	 continue;
       }
     if (!strcmp(argv[k], "-b") || !strcmp(argv[k], "-btab"))
       {
	 btab = 1; 
	 continue;
       }
     if (!strcmp(argv[k], "-a"))
       {
	 antidis = (k+1 < argc)? atoi(argv[k+1]) : 0;
	 if (antidis <= 0)
	   {
	     fatal("antidis > 0");
	   } 
	 continue;
       }
     if (!strcmp(argv[k], "-c"))
       {
	 maxchains = (k+1 < argc)? atoi(argv[k+1]) : 0;
	 if (maxchains <= 0)
	   {
	     fatal("maxchains > 0");
	   } 
	 continue;
       }
     if (!strcmp(argv[k], "-d"))
       {
	 distance = (k+1 < argc)? atoi(argv[k+1]) : 0;
	 if (distance <= 0)
	   {
	     fatal("distance > 0");
	   } 
	 continue;
       }
     if (!strcmp(argv[k], "-f"))
       {
	 fcutoff = (k+1 < argc)? atoi(argv[k+1]) : 0;
	 if (fcutoff <= 0)
	   {
	     fatal("fcutoff > 0");
	   } 
	 continue;
       }
     if (!strcmp(argv[k], "-i"))
       {
	 icutoff = (k+1 < argc)? atoi(argv[k+1]) : 0;
	 if (icutoff <= 0)
	   {
	     fatal("icutoff > 0");
	   } 
	 continue;
       }
     if (!strcmp(argv[k], "-l"))
       {
	 lgaplen = (k+1 < argc)? atoi(argv[k+1]) : 0;
	 if (qgaplen <= 0)
	   {
	     fatal("qgaplen > 0");
	   } 
	 continue;
       }
     if (!strcmp(argv[k], "-m"))
       {
	 match = (k+1 < argc)? atoi(argv[k+1]) : 0;
	 if (match <= 0)
	   {
	     fatal("match > 0");
	   } 
	 continue;
       }
     if (!strcmp(argv[k], "-o"))
       {
	 simcutoff = (k+1 < argc)? atoi(argv[k+1]) : 0;
	 if (simcutoff < 0  ||  simcutoff > 100)
	   {
	     fatal("Value for percent similarity must be 0 <= simcutoff <= 100");
	   }
	 continue;
       }
     if (!strcmp(argv[k], "-p"))
       {
	 idncutoff = (k+1 < argc)? atoi(argv[k+1]) : 0;
	 if (idncutoff < 0  ||  idncutoff > 100)
	   {
	     fatal("Value for percent identify must be 0 <= idncutoff <= 100");
	   }
	 continue;
       }
     if (!strcmp(argv[k], "-q"))
       {
	 qgaplen = (k+1 < argc)? atoi(argv[k+1]) : 0;
	 if (qgaplen < 0)
	   {
	     fatal("qgaplen > 0");
	   } 
	 continue;
       }
     if (!strcmp(argv[k], "-s"))
       {
	 mismh = (k+1 < argc)? atoi(argv[k+1]) : 0;
	 if ( mismh >= 0 )
	   {
	     fatal("Number for mismatch score must be a negative integer");
	   }
	 continue;
       }
     if (!strcmp(argv[k], "-w"))
       {
	 wordsize = (k+1 < argc)? atoi(argv[k+1]) : 0;
	 if (wordsize < 11)
	   {
	     fatal("Word size must be at least 11");
	   }	                      
	 continue;
       }
     
     if (debug > 0) 
       {
	 
       }
   }
 
 // If -out option, send output to both .btab and .dds output files.
 // Otherwise, send output to stdout in btab format if -btab, otherwise
 // send output to stdout in align format (default if no options specified).
 // 
 if (outfile_option)
   {
     char * alignsuffix = ".dds";
     char * btabsuffix = ".dds.btab";
     char alignout[strlen(outfile) + strlen(alignsuffix) + 1];
     char btabout[strlen(outfile) + strlen(btabsuffix) + 1];
     strcpy(alignout, outfile);
     strcat(alignout, alignsuffix);
     strcpy(btabout, outfile);
     strcat(btabout, btabsuffix);
     
     align = 1;
     btab = 1;
     
     alignfp = ckopen(alignout, "w"); 
     btabfp = ckopen(btabout, "w");
     if (btab_header)
       {
	 int k = 0;
	 while (headers[k] != NULL)
	   {
	     fprintf(btabfp, "%s\t", headers[k]);
	     k++;
	   }
	 fprintf(btabfp, "\n");
       }
   }
 else
   {
     if (btab)
       { 
	 align = 0;        // btab on stdout only
	 btabfp = stdout; 
	 if (btab_header)
	   {
	     int k = 0;
	     btab = 1; 
	     while (headers[k] != NULL)
	       {
		 fprintf(btabfp, "%s\t", headers[k]);
		 k++;
	       }
	     fprintf(btabfp, "\n");
	   }
       }
     else
       {
	 alignfp = stdout; 
       }
   }
 
 time(&cur_time);
 //strftime(btab_date, sizeof(btab_date), "%b %e %Y", localtime(&cur_time));
 
 ComputeAiDi();
 ComputeDa();
 /* determine the sequence lengths */
 Ap = ckopen(name1, "r");
 if ( (symbol = getc(Ap) ) != '>' )
   fatal("The DNA sequence must be in the FASTA format, starting with '>'");
 for (k=0, size = 2; ( symbol = getc(Ap)) != EOF ; size++ )
   {
     if ( isspace(symbol) )
       {
	 break;
       }
     else
       {
	 query_name[k++] = (int)symbol;
       }
     query_name[k] = '\0';
   }
 for (dnalen = 0; ( symbol = getc(Ap)) != EOF ; )
   if ( symbol != '\n' )
     ++dnalen;
 (void) fclose(Ap);
 
 if ( maxchno < maxchains ) maxchno = maxchains;
 /* allocate space for D */
 D = ( char * ) ckalloc( (dnalen + 2) * sizeof(char));
 DN = ( int * ) ckalloc( (dnalen + 2) * sizeof(int));
 dhead = ( char * ) ckalloc( size * sizeof(char));
 i = dnalen / 4 + 3;
 framept[1] = frame1 = ( int * ) ckalloc( i * sizeof(int));
 framept[2] = frame2 = ( int * ) ckalloc( i * sizeof(int));
 framept[3] = frame3 = ( int * ) ckalloc( i * sizeof(int));
 framept[0] = frame4 = ( int * ) ckalloc( i * sizeof(int));
 /* read the DNA sequence into D */
 Ap = ckopen(name1, "r");
 for (i = 0; ( symbol = getc(Ap)) != EOF && symbol != '\n' ; )
   dhead[i++] = symbol;
 dhead[i] = '\0';      
 for (dnalen = 0; ( symbol = getc(Ap)) != EOF ; )
   if ( symbol != '\n' )
     {	D[++dnalen] = symbol;
     DN[dnalen] = di[symbol];
     }
 if ( dnalen < 4 )
   fatal("The DNA sequence is too short ");
 D[0] = D[dnalen+1] = '\200';
 DN[0] = DN[dnalen+1] = DNASIZE + 1000;
 tsize = dnals4 * dnals4;
 table = (int  * ) ckalloc(tsize * sizeof(int ));
 list = ( int  * ) ckalloc( (dnalen + 2) * sizeof(int ));
 
 q = 10;
 r = 1;
 qpay = q + r * qgaplen;
 lpay = q + r * lgaplen;
 
 inits = 8 * match;
 diff8 = wordsize - 8;
 
 SetUpW();
 /* determine the length of the longest cDNA sequence */
 Bp = ckopen(name2, "r");
 for (aalen = i = 0; ( symbol = getc(Bp)) != EOF ; )
   if ( symbol == '>' )
     { if ( i > aalen )
       aalen = i;
     i = 0;
     }
   else
     if ( symbol != '\n' )
       ++i;
 if ( i > aalen )
   aalen = i;
 (void) fclose(Bp);
 P = ( char * ) ckalloc( (aalen + 2) * sizeof(char));
 PN = ( int * ) ckalloc( (aalen + 2) * sizeof(int));
 P[0] = '\377';
 PN[0] = DNASIZE + 30000;
 i = aalen / 4;
 cfm1 = ( int * ) ckalloc( (i + 2) * sizeof(int));
 diagonal = ( int * ) ckalloc( (j = dnalen + aalen) * sizeof(int));
 for ( i = 0; i < j; i++ )
   diagonal[i] = 0;
 segsize = aalen > MINISEG ? aalen : MINISEG;
 segment = ( segptr ) ckalloc( (segsize + 1) * sizeof(segtp));
 segstd = ( segptr * ) ckalloc( (segsize + 1) * sizeof(segptr));
 tuplist = ( tupptr ) ckalloc( (segsize + 1) * sizeof(tuptp));
 if (align)
   {
     fprintf(alignfp, "Query sequence (top one on the alignment): %s\n", dhead);
     fprintf(alignfp, "Query sequence length: %d\n", dnalen);
     fprintf(alignfp, "Maximum database sequence length: %d\n\n", aalen);
     fprintf(alignfp, "       QStart     QEnd  Score  SPs  LStart  LEnd Accession Description\n");
   }
 Search(P, D, PN, DN, 1);
 size = dnalen / 2;
 for ( i = 1, j = dnalen; i <= size; i++, j-- )
   { temp = D[i];
   D[i] = D[j];
   D[j] = temp;
   }
 for ( i = 1; i <= dnalen; i++ )
   { switch ( D[i] )
     { case 'A' :
     case 'a' : D[i] = 'T'; break;
     case 'T' : 
     case 't' : D[i] = 'A'; break;
     case 'C' : 
     case 'c' : D[i] = 'G'; break;
     case 'G' : 
     case 'g' : D[i] = 'C'; break;
     default  : break;
     }
   DN[i] = di[D[i]];
   }
 Search(P, D, PN, DN, 0);
 
 // close open files
 if (outfile_option)
   {
     fclose(alignfp);
     if (btab)
       {
	 fclose(btabfp);
       }
   }
 
 exit(0);
 
}

Search(P, D, PN, DN, ort)
     char *P, *D;
     int  *PN, *DN, ort;
{ int psize;		/* size of cDNA sequence */
 int  symbol;		/* The next character */
 FILE *Bp, *ckopen();
 register int i, j;
 int  code;		/* codon code */
 int  succ,  pred;	/* segment index */
 int  astart, dstart, aend, dend; /* start and end positions */
 int  numsp;		/* number of segments in a chain */
 int  idn, sim;	/* numbers of matches and similarities in a segment */
 int  tolidn, tolsim;	/* numbers of matches and similarities in a segment */
 int  len;		/* length of a segment */
 int  cleng;		/* sum of lengths of segments in a chain */
 int  cfsz1; 		/* size of cfm1 */
 int  x, y;		/* temp variables */
 int *ap, *bp, *ep;   /* pointer vars */
 int  fsize1, fsize2, fsize3, fsize4; /* lengths of sequences */
 int exon_num = 0;     // track exon numbers -bhaas
 int score = 0;        // track bit score -mfs
 
 /* put the codons of the DNA sequence in 4 frames */
 code = DN[1] * dnals2 + DN[2] * DNASIZE + DN[3];
 fsize1 = fsize2 = fsize3 = fsize4 = 0;
 for ( i = 4; i <= dnalen; i++ )
   { code = ( code * DNASIZE + DN[i] ) % dnals4;
   if ( (j = i % 4) == 0 )
     frame1[++fsize1] = code;
   else
     if ( j == 1 )
       frame2[++fsize2] = code;
     else
       if ( j == 2 )
	 frame3[++fsize3] = code;
       else
	 frame4[++fsize4] = code;
   }
 frame1[0] = frame2[0] = frame3[0] = frame4[0] = dnals4;
 frame1[fsize1 + 1] = frame2[fsize2 + 1] = frame3[fsize3 + 1]
   = frame4[fsize4 + 1] = dnals4;
 /* make a lookup table */
 for ( i = 0; i < tsize; i++ )
   table[i] = 0;
 for ( code = len = 0, i = 1; i <= dnalen; i++ )
   { if ( (j = DN[i]) == nncode ) 
     len = -1;
   code = ( code * DNASIZE + j ) % tsize;
   if ( ++len >= 8 )
     { list[j = i - 7] = table[code];
     table[code] = j;
     }
   }
 
 Bp = ckopen(name2, "r");
 symbol = getc(Bp);
 /* process each cDNA sequence */
 while ( symbol != EOF )
   {
     if ( symbol == '>' )
       {
	 i = 0;
	 j = 0;
	 accn[0] = '\0';
	 match_header[0] = '\0';
	 while(i < ACCSIZE-1)
	   {
	     symbol = getc(Bp);
	     if (symbol == EOF)
	       {
		 break;
		 fatal("The cDNA database is not in the correct format");
	       }
	     if (isspace(symbol)) { break; }
	     accn[i++] = symbol;
	   }
	 accn[i] = '\0';
	 if (symbol != '\n') 
	   { 
	     // Grab the match header if it exists 
	     while((symbol = getc(Bp)) != EOF && symbol != '\n' && i+j < ACCSIZE-1)
	       {
		 match_header[j++] = (isspace(symbol))? ' ' : symbol;
	       }
	     match_header[j] = '\0';
	   }
	 
	 if (symbol != '\n')
	   {
	     while((symbol = getc(Bp)) != EOF && symbol != '\n')
	       {}
	   }
	 
       }
     else
       {
	 fatal("The cDNA database is not in the correct format");
       }
     for (psize = 0; ( symbol = getc(Bp)) != EOF && symbol != '>'; )
       if ( symbol != '\n' )
	 {  P[++psize] = symbol;
	 PN[psize] = di[symbol];
	 }
     if ( psize >= 8 )
       { /* pack the cDNA sequence */
	 cfsz1 = 0;
	 for ( i = 1; i <= psize - 3; i += 4 )
	   cfm1[++cfsz1] = PN[i] * dnals3 + PN[i+1] * dnals2
	     + PN[i+2] * DNASIZE + PN[i+3];
	 cfm1[0] = dnals4;
	 cfm1[cfsz1 + 1] = dnals4;
	 P[psize+1] = '\377';
	 PN[psize+1] = DNASIZE + 30000;
	 cseglen = tupsize = 0;
	 Compare(PN, DN, cfm1,cfsz1);
	 ChainSeg();
	 for ( i = 1; i <= tupsize; i++ )
	   { succ = numsp = 0;
	   for ( pred = tuplist[i].last; pred ; numsp++ )
	     { pred = segment[j = pred].pred;
	     segment[j].pred = succ;
	     succ = j;
	     }
	   pred = succ;
	   tolidn = tolsim = cleng = 0;
	   for ( ; succ ; succ = segment[succ].pred )
	     { astart = segment[succ].astart;
	     dstart = segment[succ].dstart;
	     aend = segment[succ].aend;
	     len = aend - astart + 1;
	     ap = PN + astart;
	     bp = DN + dstart;
	     ep = ap + len;
	     idn = sim = 0;
	     for ( ; ap < ep; )
	       { if ( (x = *ap++) == (y = *bp++) && x != nncode )
		 { idn++;
		 sim++;
		 }
	       if ( x == nncode || y == nncode )
		 sim++;
	       }
	     cleng += len;
	     tolidn += idn;
	     tolsim += sim;
	     segment[succ].idn = idn;
	     segment[succ].sim = sim;
	     }
	   if ( 100 * tolidn < idncutoff * cleng || 100 * tolsim < simcutoff * cleng )
	     continue;
	   if ( ort )
	     { dstart = segment[tuplist[i].first].dstart;
	     dend = segment[tuplist[i].last].dend;
	     }
	   else
	     { dstart = dnalen - segment[tuplist[i].first].dstart + 1;
	     dend = dnalen - segment[tuplist[i].last].dend + 1;
	     }
	   astart = segment[tuplist[i].first].astart;
	   aend = segment[tuplist[i].last].aend;
	   if ( ++numchains < maxchno )
	     {
	       global_chain_num++;
	       if (align)
		 {
		   fprintf(alignfp,
			   "\n//\n[Alignment_chain %d]\nChain%8d %8d %6d %4d %7d %5d %s\n",
			   global_chain_num,
			   dstart, 
			   dend, 
			   tuplist[i].total, 
			   numsp, 
			   astart, 
			   aend, 
			   accn
                           );
		 }
	     }
	   succ = pred;
	   exon_num = 0;
	   for ( ; succ ; succ = segment[succ].pred )
	     { astart = segment[succ].astart - 1;
	     dstart = segment[succ].dstart - 1;
	     aend = segment[succ].aend;
	     dend = segment[succ].dend;
	     segment[succ].new = 0;
	     len = aend - astart;
	     sim = segment[succ].sim;
	     idn = segment[succ].idn;
	     if ( numchains < maxchains )
	       { 
		 exon_num++;
		 score = segment[succ].score;
		 if (btab)
		   {
		     if (debug) { fprintf(stderr, "btab output at line %d\n", __LINE__); }
		     BTAB_LINE();
		   }
		 
		 if (align)
		   {
		     fprintf(alignfp,"\n[Segment %d]\nScore: %d Identity: %d/%d (%d%%) Similarity: %d%% QStrand: %s\n",
			     exon_num,
			     score, 
			     idn, 
			     len, 
			     (int) ( (100.0 * idn) / len ),
			     (int) ( (100.0 * sim) / len ), 
			     (ort == 1 ? "Plus" : "Minus") 
			     );
		     display2(P+astart,D+dstart,aend-astart,dend-dstart, astart+1,dstart+1,ort);
		   }
	       }
	     }
	   }
	 /* pick up remaining high-scoring segments */
	 for ( i = 1; i <= cseglen; i++ )
	   if ( segment[i].new && segment[i].score >= fcutoff ) 
	     { idn = sim = 0;
	     astart = segment[i].astart;
	     dstart = segment[i].dstart;
	     aend = segment[i].aend;
	     len = aend - astart + 1;
	     ap = PN + astart;
	     bp = DN + dstart;
	     ep = ap + len;
	     for ( ; ap < ep; )
	       { if ( (x = *ap++) == (y = *bp++) && x != nncode )
		 { idn++;
		 sim++;
		 }
	       if ( x == nncode || y == nncode )
		 sim++;
	       }
	     if ( 100 * idn < idncutoff * len || 100 * sim < simcutoff * len )
	       continue;
	     if ( ort )
	       { dstart = segment[i].dstart;
	       dend = segment[i].dend;
	       }
	     else
	       { dstart = dnalen - segment[i].dstart + 1;
	       dend = dnalen - segment[i].dend + 1;
	       }
	     if ( ++numchains < maxchno )
	       {
		 global_chain_num++;
		 exon_num = 0;
		 score = segment[i].score;
		 if (align)
		   {
		     fprintf(alignfp,
			     "\n//\n[Alignment_chain %d]\nChain%8d %8d %6d %4d %7d %5d %s\n",
			     global_chain_num,
                             dstart, 
                             dend, 
                             score, 
                             1, 
                             astart, 
                             aend, 
                             accn
			     );
		   }
	       }
	     if ( numchains < maxchains )
	       {  
		 astart--;
		 dstart = segment[i].dstart - 1;
		 dend = segment[i].dend;
		 exon_num++;
		 score = segment[i].score;
		 if (btab)
		   {
		     if (debug) { fprintf(stderr, "btab output at line %d\n", __LINE__); }
		     BTAB_LINE();
		   }
		 if (align)
		   {
		     fprintf(alignfp,
			     "\n[Segment %d]\nScore: %d Identity: %d/%d (%d%%) Similarity: %d%% QStrand: %s\n",
			     1,
			     score, 
			     idn, 
			     len, 
			     (int) ( (100.0 * idn) / len ),
			     (int) ( (100.0 * sim) / len ), 
			     (ort == 1 ? "Plus" : "Minus") 
                             );
		     display2(P+astart,D+dstart,aend-astart,dend-dstart, astart+1,dstart+1,ort);
		   }
	       }
	     }
       }
   }
 (void) fclose(Bp);
}

/* gap penalty for query */
#define gapq(k)  ((hh = k) <= 0 ? q : (hh <= qgaplen ? q+r*hh : qpay))

/* gap penalty for library sequence */
#define gapl(k)  ((hh = k) <= 0 ? q : (hh <= lgaplen ? q+r*hh : lpay))

ComputeAiDi()
{ int i;
 
 nncode = DNASIZE - 1;
 for ( i = 0; i < 128; i++ )
   di[i] = nncode; 	/* code of 'N' */
 di['a'] = di['A'] = 0;
 di['c'] = di['C'] = 1;
 di['g'] = di['G'] = 2;
 di['t'] = di['T'] = 3;
 di['n'] = di['N'] = 4;
}

/* Compute the da table */
ComputeDa()
{ int  i, j, k, m;	/* index variable */
 int  a, aa, aaa;	/* integer codes */
 
 nncode = DNASIZE - 1;
 dnals3 = (dnals2 = DNAS2) * DNASIZE;
 dnals4 = dnals3 * DNASIZE;
 for ( i = 0; i < dnals4; i++ )
   da[i] = 0;
 for ( i = 0; i < nncode; i++ )
   { a = i * DNASIZE;
   for ( j = 0; j < nncode; j++ )
     { aa = (a + j) * DNASIZE;
     for ( k = 0; k < nncode; k++ )
       { aaa = (aa + k) * DNASIZE;
       for ( m = 0; m < nncode; m++ )
	 da[aaa + m] = 1;
       }
     }
   }
}

/* Set up w */
SetUpW()
{ int  i, j;	/* index varialbes */
 int  *wa;	/* row of w */
 int  a1;	/* letter 1 code */
 int  a2;	/* letter 2 code */
 int  a3;	/* letter 3 code */
 int  a4;	/* letter 4 code */
 int  x;	/* temp */
 int  s;	/* score */
 
 /* compute the substitution scores for ijkm */
 for ( i = 0; i < dnals4; i++ )
   { wa = w[i];
   a1 = i / dnals3;
   a2 = ( x = i % dnals3 ) / dnals2; 
   a3 = ( x = x % dnals2 ) / DNASIZE; 
   a4 = x % DNASIZE;
   for ( j = 0; j < dnals4; j++ )
     { s = 0;
     if ( a1 == j / dnals3 && a1 != nncode )
       s += match;
     else
       s += mismh;
     if ( a2 == ( x = j % dnals3) / dnals2 && a2 != nncode )
       s += match;
     else
       s += mismh;
     if ( a3 == ( x = x % dnals2) / DNASIZE && a3 != nncode )
       s += match;
     else
       s += mismh;
     if ( a4 == x % DNASIZE && a4 != nncode )
       s += match;
     else
       s += mismh;
     wa[j] = s;
     }
   }
 /* handle cases where the end of an input sequence is reached. */
 wa = w[dnals4];
 for ( i = 0; i <= dnals4; i++ )
   { wa[i] = MINF;
   w[i][dnals4] = MINF;
   }
}

/* computes quickly scores of alignments between the two sequences */
Compare(PN,DN,cfm,cfsz)
     int PN[], DN[];
     int  cfm[], cfsz;
{ register int  value;		/* value of a word */
 register int  i, j, k, h, f,g;/* index variables */
 register int  ii;		/* index varialbes */
 register int  *aap;		/* pointer to element of cfm */
 register int  *codp;		/* pointer to element of frame */
 register int  maxs;		/* maximum score seen */
 register int  score;		/* current score */
 register int  *asp, *dsp, *aep, *dep;	/* start and end pointers */
 int  num;		/* number of segment pairs in segment */
 int  temp1, temp2, x;	/* temporary variables */
 int  *frame;		/* frame */
 int  fnum;		/* frame number */
 int  *dt, *pt;	/* pointers to DN and PN */
 char  y;
 
 num = cseglen;
 value = cfm[1];
 for ( ii = 1, i = 2; i <= cfsz ; i++, ii += 4 )
   { value = ( value * dnals4 + cfm[i] ) % tsize;
   for ( g = table[value]; g ; g = list[g] )
     { for ( x = 0, pt = PN + ii, dt = DN + g; *(--pt) == *(--dt); )
       if ( ++x >= diff8 )
	 break;
     if ( x < diff8 )
       { for ( pt = PN + ii + 8, dt = DN + g + 8; *pt++ == *dt++; )
	 if ( ++x >= diff8 )
	   break;
       if ( x < diff8 )
	 continue;
       }
     if ( (h = diagonal[k = g - ii + aalen])
	  && segment[h].aend >= i * 4 )
       continue;
     frame = framept[ fnum = g % 4 ];
     if ( !fnum ) fnum = 4;
     j = ( g - 1 ) / 4 + 1;
     aep = aap = cfm + i + 1;
     dep = codp = frame + j + 2;
     for ( maxs = score = inits; maxs - score < distance; )
       { score += w[*aap++][*codp++];
       if ( maxs < score )
	 { maxs = score;
	 aep = aap; 
	 dep = codp;
	 }
       }
     asp = aap = cfm + i - 2;
     dsp = codp = frame + j - 1;
     for ( score = maxs; maxs - score < distance; )
       { score += w[*aap--][*codp--];
       if ( maxs < score )
	 { maxs = score;
	 asp = aap; 
	 dsp = codp;
	 }
       }
     if ( maxs >= icutoff )
       { if ( h && segment[h].astart > (asp - cfm) * 4 + 1 )
	 f = h;
       else
	 if ( (f = ++num) > segsize )
	   DoubleSegSpace();
       if ( f <= segsize )
	 { temp1 = (asp - cfm) * 4 + 1;
	 temp2 = ( dsp - frame ) * 4 + fnum;
	 for ( ; (y = PN[--temp1]) == DN[--temp2] && y != nncode; )
	   maxs += match; 
	 for ( ; PN[++temp1] != DN[++temp2]; )
	   maxs -= mismh; 
	 segment[f].astart = temp1;
	 segment[f].asp = asp + 1;
	 segment[f].dstart = temp2;
	 segment[f].dsp = dsp + 1;
	 segment[f].antis = temp1 + temp2;
	 temp1 = (--aep - cfm) * 4;
	 temp2 = (--dep - frame ) * 4 + fnum - 1;
	 for ( ; (y = PN[++temp1]) == DN[++temp2] && y != nncode; )
	   maxs += match; 
	 for ( ; PN[--temp1] != DN[--temp2]; )
	   maxs -= mismh; 
	 segment[f].aend = temp1;
	 segment[f].dend = temp2;
	 segment[f].antid = temp1 + temp2;
	 segment[f].score = segment[f].total = maxs;
	 segment[f].fnum = fnum;
	 segment[f].new = 1;
	 segment[f].pred = 0;
	 segment[f].first = f;
	 segment[f].diag = k;
	 diagonal[k] = f;
	 }
       }
     }
   }
 for ( i = cseglen + 1; i <= num; i++ )
   diagonal[segment[i].diag] = 0;
 cseglen = num;
}

/* It sorts the segments by antidiagonals and chains close segments. */
ChainSeg()
{ int  i, j, d;		/* index variables */
 int  sorted;		/* flag */
 int  found;		/* flag */
 int  recent;		/* index of most recently visited element */
 int  temp, ss, hh;	/* temporary variables */
 int  first;		/* start segment in a chain */
 int  total;		/* score of the chain */
 segptr  p;
 int  *asp, *dsp;	/* pointers to PN and frame */
 int  iscore;		/* segstd[i]->score */
 int  is3;		/* score of 3' portion of a segment */
 int  am, dm;		/* gap lengths */
 
 /* make an array of pointers and sort the array */
 p = segment;
 for ( i = 1; i <= cseglen; )
   segstd[i++] = ++p;
 for ( i = 1; i < cseglen; i++ )
   { sorted = 1;
   for ( j = cseglen; j > i; j--)
     if ( segstd[j-1]->antid > segstd[j]->antid )
       { p = segstd[j];
       segstd[j] = segstd[j-1];
       segstd[j-1] = p;
       sorted = 0;
       }
   if ( sorted )
     break;
   }
 recent = 0;
 /* chain segments */
 for ( i = 1; i <= cseglen; i++ )
   { if ( (j = i - 1) && segstd[i]->antis - segstd[j]->antid < antidis )
     { asp = segstd[i]->asp;
     dsp = segstd[i]->dsp;
     total = (segstd[i]->aend - segstd[i]->astart - 6) / 4;
     for ( d = ss = 0, temp = w[*asp++][*dsp++]; d < AAOVER; d++ ) 
       if ( d <= total )
	 { fdst[d] = ss = temp;
	 temp += w[*asp++][*dsp++];
	 }
       else
	 fdst[d] = ss;
     }
   iscore = segstd[i]->score;
   for ( ; j && segstd[i]->antis - segstd[j]->antid < antidis; j-- )
     { is3 = iscore; 
     if ( (temp = (am = segstd[i]->astart - segstd[j]->aend)/4) > aaover &&
	  (d = (dm = segstd[i]->dstart - segstd[j]->dend) / 4) > aaover &&
	  ( ( temp > 0 && d > 0 ) ||
	    ( is3 = iscore - fdst[d < temp ? (-d) : (-temp)]) >= icutoff ) &&
	  (ss = segstd[j]->total + is3 - gapq(dm-1) - gapl(am-1)) > segstd[i]->total)
       { segstd[i]->first = segstd[j]->first;
       segstd[i]->pred = segstd[j] - segment;
       segstd[i]->total = ss;
       }
     }
   if ( (total = segstd[i]->total) >= fcutoff )
     { first = segstd[i]->first;
     found = 0;
     if ( recent && tuplist[recent].first == first )
       found = 1;
     else
       for ( d = tupsize; d ; d-- )
	 { if ( tuplist[d].first == first )
	   { found = 1;
	   recent = d;
	   break;
	   }
	 if ( segment[first].antis - tuplist[d].antid > antidis )
	   break;
	 }
     if ( found )
       { if ( tuplist[recent].total < total )
	 { tuplist[recent].total = total;
	 tuplist[recent].last = segstd[i] - segment;
	 ss = tuplist[recent].antid = segstd[i]->antid;
	 for (d = recent + 1; d <= tupsize; d++ )
	   tuplist[d].antid = ss;
	 }
       }
     else
       { tuplist[recent = ++tupsize].total = total;
       tuplist[recent].first = first;
       tuplist[recent].last = segstd[i] - segment;
       tuplist[recent].antid = segstd[i]->antid;
       }
     }
   }
}

/* Allocate more space to segment, segstd and tuplist */
DoubleSegSpace()
{ int    i;
 char   *ckalloc();		/* space-allocating function */
 int    oldsize;		/* old size */
 segptr oldsegment; 		/* for holding segment pairs */
 segptr *oldsegstd;		/* a sorted list of segment pairs */
 tupptr oldtuplist;	/* a list of disjoint chains of segments */
 segptr one, two;		/* temp */
 
 oldsize = segsize;
 segsize *= 2;
 oldsegment = segment;
 oldsegstd = segstd;
 oldtuplist = tuplist;
 segment = ( segptr ) ckalloc( (segsize + 1) * sizeof(segtp));
 segstd = ( segptr * ) ckalloc( (segsize + 1) * sizeof(segptr));
 tuplist = ( tupptr ) ckalloc( (segsize + 1) * sizeof(tuptp));
 one = segment;
 two = oldsegment;
 for ( i = 1; i <= oldsize; i++ )
   { one++;
   two++;
   one->dstart = two->dstart;
   one->dend = two->dend;
   one->astart = two->astart;
   one->aend = two->aend;
   one->new = two->new;
   one->score = two->score;
   one->total = two->total;
   one->pred = two->pred;
   one->first = two->first;
   one->antis = two->antis;
   one->antid = two->antid;
   one->diag = two->diag;
   one->fnum = two->fnum;
   one->asp = two->asp;
   one->dsp = two->dsp;
   one->idn = two->idn;
   one->sim = two->sim;
   }
 (void) free(oldsegment);
 (void) free(oldsegstd);
 (void) free(oldtuplist);
}

/* Alignment display routine */

static char ALINE[65];	/* cDNA line */
static char DLINE[65];	/* DNA line */
static char CLINE[65];	/* pair type indicator */ 

display2(A,B,M,N,AP,DP,orient)
     char A[], B[]; int M, N; int AP, DP, orient;
{ register char *a;	/* cDNA line pointer */ 
 register char *d;	/* DNA line pointer */
 register char *c;	/* center line pointer */
 register int  i, j;	/* index variables */
 int  dp;	/* DNA sequence position */
 int  ap;	/* cDNA sequence position */
 int mdp;	/* the DNA position relative to the plus strand */
 int  dp2;	/* DNA sequence position for right column */
 int  ap2;	/* cDNA sequence position for right column */
 int  mdp2;	/* the DNA position relative to the plus strand */
 char x, y;	/* temp */
 int  temp;   /* temp */
 
 i = j = 0;
 ap = AP;
 dp = DP;
 a = ALINE;
 d = DLINE;
 c = CLINE;
 while (i < M || j < N)
   { *a++ = x = A[++i];
   *d++ = y = B[++j];
   *c++ = ((temp = di[x]) == di[y] && temp != nncode) ? '|' : ' ';
   if (a >= ALINE+60 || i >= M && j >= N)
     { if ( i >= M && j >= N )
       for ( ; a < ALINE+60 ; )
	 *a++ = *d++ = ' ';
     *a = *d = *c = '\0';
     ap2 = AP + i - 1;
     dp2 = DP + j - 1;
     if ( orient )
       { mdp = dp;
       mdp2 = dp2;
       }
     else
       { mdp = dnalen - dp + 1;
       mdp2 = dnalen - dp2 + 1;
       }
     if (align)
       {
	 fprintf(alignfp,
		 "\n%8d %s %d\n         %s\n%8d %s %d\n",
		 mdp,
		 DLINE,
		 mdp2,
		 CLINE,
		 ap,
		 ALINE,
		 ap2
		 );
       }
     ap = AP + i;
     dp = DP + j;
     a = ALINE;
     d = DLINE;
     c = CLINE;
     }
   }
}

/* lib.c - library of C procedures. */

/* fatal - print message and die */
fatal(msg)
     char *msg;
{
  (void) fprintf(stderr, "%s\n", msg);
  exit(1);
}

/* fatalf - format message, print it, and die */
fatalf(msg, val)
     char *msg, *val;
{
  (void) fprintf(stderr, msg, val);
  (void) putc('\n', stderr);
  exit(1);
}

/* ckopen - open file; check for success */
FILE *ckopen(name, mode)
     char *name, *mode;
{
  FILE *fopen(), *fp;
  
  if ((fp = fopen(name, mode)) == NULL)
    fatalf("Cannot open %s.", name);
  return(fp);
}

/* ckalloc - allocate space; check for success */
char *ckalloc(amount)
     int amount;
{
  char *malloc(), *p;
  
  if ((p = malloc( (unsigned) amount)) == NULL)
    fatal("Ran out of memory.");
  return(p);
}
