/*  A Program for Comparing a DNA sequence against a protein sequence database
    (DPS, DNA Protein Search)

    Proper attribution of the author as the source of the software would
    be appreciated:
    Huang, X.  (1996)
    Fast Comparison of a DNA Sequence with a Protein Sequence Database.
    Microbial & Comparative Genomics, 1(4): 281-291.

    Acknowledgments
       I thank the following people for discussions and suggestions:
       Mark Adams, Tony Kerlavage, Brendan Loftus, Steve Rounsley,
       Granger Sutton, and Jinghui Zhang.
       The integration of DPS with NAP was performaed at TIGR.

    The DPS program compares a DNA sequence to a protein database.
    The DPS enhances the existing methods by addressing the problems
    of frameshifts and introns. DPS computes high-scoring chains of
    segment pairs, where segment pairs in a chain can be from different
    reading frames and there can be an intervening DNA sequence
    between adjacent segment pairs in a chain.

    The sensitivity of the program depends on the word size parameter W.
    Increased sensitivity comes at the expense of decreased speed.
    A hit is an exact coding of the amino acid word in the
    DNA sequence. Each hit is extended in both directions.
    The extesion stops if the score drops by more than the D distance.

    For a segment pair s, let nstart(s) and nend(s) denote
    the starting and ending positions of the DNA segment in
    the DNA sequence, let astart(s) and aend(s) denote
    the starting and ending positions of the protein segment in
    the protein sequence, and let score(s) denote the score of s.
    The first antidiagonal of a segment pair s is defined to
    be antis(s) = nstart(s) + 3 * astart(s), and
    the last antidiagonal of s is defined to be
    antid(s) = nend(s) + 3 * aend(s), where the protein
    position is scaled up by a factor of 3 since
    an amino acid corresponds to three nucleotides.

    A chain of segment pairs is a list of segment pairs in increasing
    order of their last antidiagonal such that each segment pair
    is not far from its predecessor and adjacent segment pairs
    do not have a large overlap.
    Specifically, any two adjacent segment pairs s and s' in the list
    satisfy the requirement:

       antis(s')  -  antid(s)  <  A, astart(s') - aend(s) > -B,
       and nstart(s') - nend(s) > -3 * B.

    for some nonnegative integers A and B.
    Here A is called the maximum number of antidiagonals between s and s',
    and B is the AAOVER parameter in the program.
    The value for A can be specified at the command line.
    A long intron requres use of a large value for A.

    For a chain of segment pairs to be reported, each segment pair
    in the chain must have a score no less than the initial segment score cutoff,
    and the chain must have a score no less than the final chain score cutoff.

    The program shows at most C number of chains of segments.
    If there are extra chains to be reported, the program only
    prints out the chain headings without showing the segment alignments.

    The DPS program is written in C and runs under Unix systems on
    Sun workstations and under DOS systems on PCs.
    We think that the program is portable to many machines.

    To compare a DNA sequence in file DNA_Seq and a protein
    database in file Protein_Database, use a command of form

        dps  DNA_Seq  Protein_Database  BLOSUM62  [options]  > result

    where 
       DNA_Seq is a file of a query DNA sequence in FASTA format,
       Protein_Database is a file of protein database in FASTA format,
       BLOSUM is a file of specially formatted BLOSUM matrix, and
       the options are:
	  -a  N  specify max number N of antidiagonals between segments,
          -b|btab  emit btab output
          -bh      emit btab output with headers
	  -c  N  specify max number N of chain alignments reported,
	  -d  N  specify distance N for extension,
	  -f  N  specify final chain score cutoff N,
	  -i  N  specify initial segment score cutoff N,
	  -w  N  specify amino acid word size of N <= 5.

	BLOSUM62 matrix:
ARNDCQEGHILKMFPSTWYVBZX
 4
-1  5
-2  0  6
-2 -2  1  6
 0 -3 -3 -3  9
-1  1  0  0 -3  5
-1  0  0  2 -4  2  5
 0 -2  0 -1 -3 -2 -2  6
-2  0  1 -1 -3  0  0 -2  8
-1 -3 -3 -3 -1 -3 -3 -4 -3  4
-1 -2 -3 -4 -1 -2 -3 -4 -3  2  4
-1  2  0 -1 -3  1  1 -2 -1 -3 -2  5
-1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5
-2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6
-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7
 1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4
 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5
-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11
-2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7
 0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4
-2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4
-1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4
 0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1
    A sample DNA file:
>DNA sequence
CTCATTTTTTCTTGCCCGTGTTTTAAATGTTTTCATCCACAGCATTTGAT
GGGATGATTGGAAGTGAGACGTTCGAGAAAATCCATATTTTGAGTCAAGA
ATTCAGATAATATACTGAGATGATTAGGTATGGCTGGGTTCTACAAAAAC
ACAAATATCCGGCTAGCAATGATCACTGAGCAAATTAAAGCGTTAACTCA
CTCATTATTGTAGCTTATGCGTTTCTCCTCCTCTCTTTTTTTCCTCGAAC
CGGAGTGGAAGATCCAATAACGTAATATTACTGATGTTGTTATTAAAGCT
GGCAAAAATAACATGAGGCGTAAAACCGCACTGCGGTAAGATGAGGGTAT
AAGGTGGAGATCAGGCGAACAAGCTGTTCTAAATCATACATATGTACAAT
GAGAACGTGTAACGATCCAATGAGCGTTTCATGATGCCATTGTTTAATCA
GAGTGATGAAAAAGAAATATTTGCGACCTTTTTTCGTTACATTGATCGTG
AAATTTTAATCAAAGATAATATAAGGACGTGAGATATTTATCTTTTTACT
TGAAATTAACAATAGAATTGCGCTAAGCGGAATAAGAGCTTTCGTAAACC
TTTCTATTTGCACCATTGCGTCAACGTATAAAATGGTATGACCTTTACAC
AAACGCATGCTTATAATCTTATGTTTTTCATAGGGTGTAATTTGGTTGAT
GACGTAGTCTAAATTTGATGCTATCTGCAATTGAGGTACATATAAGAGGT
CAATTTCGGGACCAACCCTTTTAATCGAAAAAAACGTAATTCACTAGGGC
AAGGGAGAACTTAGCAGCTAATATCGTAAACCTTTCATACTAAAAAAATG
CACTTACCATCAACAAAAAACTCAGGACCAATTTCCAAGCTTTTCTAGGT
GATTGCCTATAACACAAAAAGATTCGCTCATACATGAGATTTTTACATGT
AATAGCAATTTGTTCCGATCAGTTGAAGGTCATCAACGCACGGCAGGTAC
ATCCACACCTATCACAAAGCCCTTCAATAATTCACCTACGTAAAGTTATA
CCGAAACATGCAAAATCCATGAAAAATTCTGTATGATAACGATCATATCC
TTTTGTATTGGTGGTACGATGCTCAAAGATAGTTATTGTTGCACCTGAGG
CAAAAGCGGAAATGAAAAATCCAGATGGGGCCAAAAGCAGAAGTATTGTG
TACAACAATTGCTTCAGCAGTTTACCAAACCGTTTCCCAGCAATCATCAA
AAGTTGCTTTAGCCACATTTCCGCAAGATATCTTTGTGGCTCAACGAAGA
GGGCTATTCCAAATGCAATACAATACTAGCCGCTAGTGATCCATGCTTAT
AGGCAATTTGATTGATAACTGGCCGTTCTATTAAGGAGTCAATGCTAACC
ACATATAATGCATATAGGATTTGGCCTCTGCTGACCGTAATACTAGACAA
GGAATATAAAACAACAACGTAACCCAGCATAAAAACGATATAAATAAAAA
AAGAAACCAGATCATAAAGTTTGAGGGCACATCCCTCATGGTTTCAAAAT
CTCGTACATTGACTCAAACCTCGAGATCTTGTTTAAACGAACATAAGAAC
AGCGGTACCAAGTGACACATTGCGTACTGTGTAGTATGCGCCGTATATAA
CTTTTTTTTTCTGAAGGTACTTTGAATTACAATCTATTTTTTACAGTTCC
TATGGCAGGGGTTGAAGATATTTGGGTCTGAACCATAGCAGGATTAGTGT
TATAGTAGGTATGTGAATAGAAGCTAACAAAATGAGATGAACCTCATACA
AAGTCGTAGAGAAAACTGCTAACAGAAGAGCTGCGCCTTGAAATCGTATC
TCTAAGCTTATAATAAATTGAAAGGAAAAAATACGTGGTAAATGCAAGCG
ACCAAAAGGCTACGGCCCAACGCTAACCCGCCGATAGGTGCATAATCTAA
TTTACCTCCACCAGCAGGAGCCCTTTTTCTAAGTAATAAGCAAACCAGAT
AACTTACATCTTGCTGTAGGAAACAAAAGCCGGAATAATGGTTCACTCAT
ATTCTTCGTGTGAAACACAGAAGAAATCCAATATTTGCTTCAGTATTTAT
CTCTAAAAATTGGTCCTACATTGGAAACCATAAACCAATTATAACCGGTG
TACGAATTGTAAGCTAGTTCTGGAAATGTCATGTTGCGCAGGTAAAAGTG
GAGCTGAATTGTATATCTGTTTTGATCATTATTATCCCTCTGGGTGAGTG
GAAATATCAATAAAATGCAATGGCACATTTAATATCCTTCTCTTAATTCC
GTGATTTATAACATCTTGATGCCAGAAACACCTTTCGGATCCGGCAATAA
AGCGGAGATTAGCACGCTTTTCGCCGGTCCTACGGATTTAGTGTTGGCTA
TTGTTGAGATTAGTAATACGCAGAGAATTTTTCTACCGGTGAAGCGACCA
TCTCAGATTATTAGGTCAAGCAATAA
	A sample protein datbase file:
>P41260;  HEMOGLOBIN I (HB I).  GLB1_LUCPE 
SLEAAQKSNVTSSWAKASAAWGTAGPEFFMALFDAHDDVFAKFSGLFSGAAKGTVKNTPE
MAAQAQSFKGLVSNWVDNLDNAGALEGQCKTFAANHKARGISAGQLEAAFKVLSGFMKSY
GGDEGAWTAVAGALMGEIEPDM
>P41261;  HEMOGLOBIN II (HB II).  GLB2_LUCPE 
TTLTNPQKAAIRSSWSKFMDNGVSNGQGFYMDLFKAHPETLTPFKSLFGGLTLAQLQDNP
KMKAQSLVFCNGMSSFVDHLDDNMLVVLIQKMAKLHNNRGIRASDLRTAYDILIHYMEDH
NHMVGGAKDAWEVFVGFICKTLGDYMKELS
>P41262;  HEMOGLOBIN III (HB III).  GLB3_LUCPE 
SSGLTGPQKAALKSSWSRFMDNAVTNGTNFYMDLFKAYPDTLTPFKSLFEDVSFNQMTDH
PTMKAQALVFCDGMSSFVDNLDDHEVLVVLLQKMAKLHFNRGIRIKELRDGYGVLLRYLE
DHCHVEGSTKNAWEDFIAYICRVQGDFMKERL
>P36032;  HYPOTHETICAL 52.3 KD PROTEIN IN FRE2 5'REGION.  YKW1_YEAST
MSEERHEDHHRDVENKLNLNGKDDINGNTSISIEVPDGGYGWFILLAFILYNFSTWGANS
GYAIYLAHYLENNTFAGGSKLDYASIGGLAFSCGLFFAPVITWLYHIFSIQFIIGLGILF
QGAALLLAAFSVTLWEIYLTQGVLIGFGLAFIFIPSVTLIPLWFRNKRSLASGIGTAGSG
LGGIVFNLGMQSILQKRGVKWALIAQCIICTSLSTIALMLTRTTHQGLRQHKRSYKFELL
DYDVLSNFAVWLLFGFVSFAMLGYVVLLYSLSDFTVSLGYTSKQGSYVSCMVSVGSLLGR
PIVGHIADKYGSLTVGMILHLVMAILCWAMWIPCKNLATAIRFGLLVGSIMGTIWPTIAS
IVTRIVGLQKLPGTFGSTWIFMAAFALVAPIIGLELRSTDTNGNDYYRTAIFVGFAYFGV
SLCQWLLRGFIIARDEIAVREAYSADQNELHLNVKLSHMSKCLFRYKQLPRRV
*/



/* 

Improvement to hsp chaining algorithm:  08-08-2005 xqhuang@cs.iastate.edu

A mega-chain is used to address this problem. A mega-chain is an ordered list 
of normal chains, where each chain correspsonds to an alignment of a gene copy 
with the protein sequence. DPS computes a mega-chain with the maximum score.
Then the mega-chain is broken up to normal chains for output.

A mega-chain with two chains has a larger score than a mega-chain
with one chain. In other words, the correct solution seems to
be a mega-chain with the maximum score.

The minimum distance between chains in a mega-chain is controlled by the
-t option. 

*/

#include   <stdio.h>

#define  FCUTOFF      100	/* final cutoff on alignment score */       
#define  ICUTOFF      35	/* initial cutoff on segment score */       
#define  WORDSIZE      4   	/* degault value for wordsize */
#define  DISTANCE     20        /* default value for distance */
#define  MINISEG  200000   	/* minimum value for segsize */
#define  ANTIDIS   50000   	/* maximum no. of antidiagonals between segments */
#define  AAOVER       50   	/* maximum length of overlap between segments */
#define  MAXCHAINS   300        /* default value for numchains */

#define  DNASIZE       5	/* size of DNA alphabet: A, C, G, T, N */
#define  AASIZE       23	/* size of protein alphabet */
#define  DNAS2  (DNASIZE * DNASIZE) /* DNASIZE square */
#define  DNAS3  (DNASIZE * DNAS2)   /* DNASIZE cube */
#define  MINF      -10000	/* negative infinity */
#define  ACCSIZE     300        // accession id buffer in fasta file
#define  FNAMESZ    1023        // longest output filename size

// Version info
static float VersionS = 1.51f;
static char * BuildS = "$Revision: 1.20 $";

static int  wordsize;		/* size of words */
static int  distance;		/* stop extension if score drops by distance */
static int  maxchains;          /* max number of chain alignments reported */
static int  numchains = -1;          /* number of chain alignments */
static int  fcutoff;            /* final cutoff on alignment score */
static int  icutoff;            /* initial cutoff on segment score */
static int  antidis;            /* maximum no. of antidiagonals between segments */
static int  maxchno = 300000;   /* max number of chains reported */
static int  intergenedis = 500; /* minimum distance between genes */

static int  *frame1, *frame2, *frame3;	/* sequences of codon codes */
static int  *framept[3];		/* pointers to frames */
static int  *table;	/* lookup table for frames */
static int  tsize;	/* table size */
static int  *list;   	/* lookup list for frames */
typedef struct MSP   {          
       int   dstart;		/* start position in codon sequence */
       int   dend;		/* end position in codon sequence */
       int   astart;		/* start position in protein sequence */
       int   aend;		/* end position in protein sequence */
       int   inuse;		/* 1, output; 0, not used in output */
       int   score;		/* score of segment pair */
       int   total;		/* score of a list ending with the segment */
       int   pred;		/* predecessor of the segment in the list */
       int   first;		/* first segment in the list */
       int   antis;		/* start antidiagonal no. of the segment */
       int   antid;		/* end antidiagonal no. of the segment */
       int   diag;		/* diagonal no. of the segment */
       int   fnum;		/* frame number */
       int   *asp;		/* pointer to the start position in PN */
       int   *dsp;		/* pointer to the srart position in frame */
       short kind;		/* 0, intra-chain link; 1, inter-chain link */
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
static char   accn[ACCSIZE];	/* accession number */
static int    aaover;		/* aaover = -AAOVER, aaover can't be modified */
static int    fdst[AAOVER];	/* score of each 5' end of a segment */

static char   aa[AASIZE][4];	/* integer code to three-letter code */
static char   aa2[AASIZE][4];	/* integer code to one-letter code */
static int    w[AASIZE+1][DNAS3+1]; /* score table for substitutions */
static int    da[DNAS3];	/* codon integer code to aa integer code */
static int    ai[128];		/* aa letter to integer code */
static int    di[128];		/* DNA letter to integer code */
static int    pam[AASIZE+1][AASIZE+1]; /* pam matrix */
static int    ochre;		/* stop codon code */
static int    amber;		/* stop codon code */
static int    uga;		/* stop codon code */
static int    dnals2;		/* DNASIZE square */
static int    dnals3;		/* DNASIZE cube */
static int    gtcode;		/* the integer code of GT */
static int    agcode;		/* the integer code of AG */
static char   mt1[AASIZE][DNASIZE]; /* patial match indicator for base 1 */
static char   mt2[AASIZE][DNAS2];   /* patial match indicator for base 2 */
static char   mt3[AASIZE][DNAS3];   /* patial match indicator for base 3 */
static int    dnalen;		/* length of DNA sequence */
static int    aalen;		/* maximum protein sequence length */
static int    aalen3;		/* maximum protein sequence length */

static int match, mismh;		/* max and min substitution weights */
static char *name1, *name2;		/* names of sequence files    */
static char *name3;                     /* names of matrix file */
static int  q, r;       /* gap penalties */
static int  qr;         /* qr = q + r */
static int  gaplen;     /* minimum length for constant-cost insertion */
static int  pay;	/* constant-cost for long insertion */
static int global_chain_num = 0; //bhaas added to track seq alignment numbers.
static char match_header[ACCSIZE]; // mfs added to track db match header string
static char outfile[FNAMESZ];      // name of output file if not stdout
static int outfile_option = 0;     // whether output file option selected
static int align = 1;              // whether align option selected
static FILE * alignfp;             // alignment output file pointer


// btab option (mfs)
static char btab_date[100];         // printable date  Mmm dd YYYY
static char btab_line[1024];        // line on which repeated btab info printed
static char * program = "";         // program name
static char * type = "dps";         // blast type field
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
  fprintf(stderr,"dps - Compare a DNA sequence against a protein sequence database\n\n");
  fprintf(stderr,"Usage: dps DNA_Seq Prot_Database BLOSUM [-btab] [-out <file>] [options]\n\n");
  fprintf(stderr,"  DNA_Seq        file of one query DNA in FASTA format\n");
  fprintf(stderr,"  Prot_Database  file of protein database in FASTA format\n");
  fprintf(stderr,"  BLOSUM         file of specially formatted BLOSUM matrix\n");
  fprintf(stderr,"  -b|btab        emit output in btab format\n");
  fprintf(stderr,"  -bh            emit output in btab format with column headers\n");
  fprintf(stderr,"  -out <file>    alignment report to <filename>.dps \n");
  fprintf(stderr,"                 btab report to <filename>.dps.btab \n");
  fprintf(stderr,"  Options (default values):\n");
  fprintf(stderr,"  -a  N  specify max number N of antidiagonals between segments (50000)\n");
  fprintf(stderr,"  -c  N  specify max number N of chain alignments reported (300)\n");
  fprintf(stderr,"  -d  N  specify distance N for extension (20)\n");
  fprintf(stderr,"  -f  N  specify final chain score cutoff N (100)\n");
  fprintf(stderr,"  -i  N  specify initial segment score cutoff N (35)\n");
  fprintf(stderr,"  -t  N  specify minimum distance between tandem genes N (500)\n");
  fprintf(stderr,"  -w  N  specify amino acid word size N <= 5 (4)\n");
  fprintf(stderr,"\n");
  fprintf(stderr,"See also:\n");
  fprintf(stderr,"  dds ext filter gap2 nap show\n");
  fprintf(stderr,"\n");
  fprintf(stderr,"Huang, X.  Fast Comparison of a DNA Sequence with a Protein Sequence Database\n");
  fprintf(stderr,"  Microbial & Comparative Genomics, 1(4): 281-291 (1996).\n");
  fprintf(stderr,"\n");
}

main(argc, argv) int argc; char *argv[];
{ char *D,  *P;				/* Storing two sequences */
  int  *DN, *PN;			/* two sequences of integer codes */
  int  symbol;				/* The next character	      */
  int  ms;				/* User-supplied weights      */
  FILE *Bp, *Ap, *Sp, *ckopen();
  char *ckalloc();			/* space-allocating function  */
  register int i, j, k;
  char alph[129], *s;			/* alphabet */
  int  size;				/* size of alphabet */
  int  temp;
  char *dhead;                          /* the name of the genomic DNA */
  long  cur_time;

  // static initializations
  wordsize = WORDSIZE;
  distance = DISTANCE;
  maxchains = MAXCHAINS;
  fcutoff = FCUTOFF;
  icutoff = ICUTOFF;
  antidis = ANTIDIS;
  aaover = - AAOVER;      /* don't change aaover */

  // initialize operational parameters
  btab = 0;
  btab_header = 0;
  align = 1;

  // parse command line (mfs)
  program = argv[0];
  if ( argc < 2 )
  {
    usage();  // exit with usage error
    exit(1);
  }
  name1 = (argc > 1)? argv[1] : "";
  name2 = (argc > 2)? argv[2] : "";
  name3 = (argc > 3)? argv[3] : "";     // BLOSUM file or mismatch value

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
      fprintf(stderr, "name3 = %s\n", name3);
      fprintf(stderr, "btab = %d\n", btab);
      fprintf(stderr, "antidis = %d\n", antidis);
      fprintf(stderr, "maxchains = %d\n", maxchains);
      fprintf(stderr, "distance = %d\n", distance);
      fprintf(stderr, "fcutoff = %d\n", fcutoff);
      fprintf(stderr, "icutoff = %d\n", icutoff);
      fprintf(stderr, "wordsize = %d\n", wordsize);
      fprintf(stderr, "intergenedis = %d\n", intergenedis);
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
    if (!strcmp(argv[k], "-btab"))
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
    if (!strcmp(argv[k], "-w"))
    {
      wordsize = (k+1 < argc)? atoi(argv[k+1]) : 0;
      if (wordsize < 11)
      {
        fatal("Word size must be at least 11");
      }    
      continue;
    }
    if (!strcmp(argv[k], "-t"))
    {
      intergenedis = (k+1 < argc)? atoi(argv[k+1]) : 0;
      if (intergenedis <= 0)
      {
        fatal("intergenedis > 0");
      }
      continue;
    }

    if (debug > 0)
    {
     
    }




  }


  // If -out option, send output to both .btab and .dps output files.
  // Otherwise, send output to stdout in btab format if -btab, otherwise
  // send output to stdout in align format (default if no options specified).
  //
  if (outfile_option)
  {
    char * alignsuffix = ".dps";
    char * btabsuffix = ".dps.btab";
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
  // strftime(btab_date, sizeof(btab_date), "%b %e %Y", localtime(&cur_time));



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
	i = dnalen / 3 + 3;
	framept[1] = frame1 = ( int * ) ckalloc( i * sizeof(int));
	framept[2] = frame2 = ( int * ) ckalloc( i * sizeof(int));
	framept[0] = frame3 = ( int * ) ckalloc( i * sizeof(int));
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
	if ( dnalen < 3 )
	  fatal("The DNA sequence is too short ");
	for ( tsize = i = 1; i <= wordsize; i++ )
	   tsize *= AASIZE;
	table = (int  * ) ckalloc(tsize * sizeof(int ));
        list = ( int  * ) ckalloc( (dnalen + 1) * sizeof(int ));

	gaplen = 15;
	q = 15;
	r = 1;
	pay = q + r * gaplen;
	qr = q + r;

	/* check if the argument represents a negative integer */
	s = name3;  
	if ( *s == '-' ) s++;
	for ( ; *s >= '0' && *s <= '9' ; s++ );
	if ( *s == '\0' )
	  { (void) sscanf(name3,  "%d", &ms);
	    if ( ms >= 0 )
	       fatal("The mismatch weight is a negative integer");
	    match = 10;
	    mismh = ms;
	    /* set match and mismatch weights */
	    for ( i = 0; i < AASIZE ; i++ )
	      for ( j = 0; j < AASIZE ; j++ )
	         if (i == j )
	            pam[i][j] = 10;
	         else
	            pam[i][j] = mismh;
	  }
	else
	  { /* read a file containing alphabet and substitution weights */
	    Sp = ckopen(name3,   "r");
	    (void) fscanf(Sp, "%s", alph);
	    size = strlen(alph);
	    match = mismh = 0;
	    /* Initialize pam[][] */
	    for ( i = 0; i < AASIZE ; i++ )
	      for ( j = 0; j < AASIZE ; j++ )
	            pam[i][j] = 0;
	    for ( i = 0; i < size ; i++ )
	      for ( j = 0; j <= i ; j++ )
		{ (void) fscanf(Sp, "%d", &ms);
		  pam[ai[alph[i]]][ai[alph[j]]] = pam[ai[alph[j]]][ai[alph[i]]] = ms;
		  if ( ms > match ) match = ms;
		  if ( ms < mismh ) mismh = ms;
		}
	  }
	pam[AASIZE][AASIZE] = 0;
	SetUpW();
	SetUpMt();
	ComputeAa();
	/* determine the length of the longest protein sequence */
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
	aalen3 = 3 * aalen;
	diagonal = ( int * ) ckalloc( (j = dnalen + aalen3) * sizeof(int));
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

	exit(0);
}

Search(P, D, PN, DN, ort)
char *P, *D;
int  *PN, *DN, ort;
{ int psize;		/* size of protein sequence */
  int  symbol;		/* The next character */
  FILE *Bp, *ckopen();
  register int i, j;
  int  code;		/* codon code */
  int  succ,  pred;	/* segment index */
  int  astart, dstart, aend, dend; /* start and end positions */
  int  numsp;		/* number of segments in a chain */
  int  *asp, *dsp, *aep;/* pointers to PN and frame */
  int  idn;		/* number of matches in a segment */
  int  len;		/* length of a segment */
  int  fsize1, fsize2, fsize3;	/* lengths of sequences */
  int exon_num; //numbering of chain segments... bhaas addition
  int score = 0;        // track bit score -mfs

  //numchains = 0;
	/* put the codons of the DNA sequence in 3 frames */
	code = DN[1] * DNASIZE + DN[2];
	fsize1 = fsize2 = fsize3 = 0;
	for ( i = 3; i <= dnalen; i++ )
	 { code = ( code * DNASIZE + DN[i] ) % dnals3;
	   if ( (j = i % 3) == 0 )
	     frame1[++fsize1] = code;
	   else
	    if ( j == 1 )
	      frame2[++fsize2] = code;
	    else
	      frame3[++fsize3] = code;
	 }
	frame1[0] = frame2[0] = frame3[0] = dnals3;
	frame1[fsize1 + 1] = frame2[fsize2 + 1] = frame3[fsize3 + 1] = dnals3;
	/* set up tables for frames */
	for ( i = 0; i < tsize; i++ )
	    table[i] = 0;
        SetTable(frame1, fsize1, 1);
        SetTable(frame2, fsize2, 2);
        SetTable(frame3, fsize3, 3);

	Bp = ckopen(name2, "r");
	symbol = getc(Bp);
	/* process each protein sequence */
	PN[0] = AASIZE;
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
                fatal("The protein database is not in the correct format");
              }
              if (isspace(symbol)) { break; }
              accn[i++] = symbol;
            }
            accn[i] = '\0';
            if (symbol != '\n')
            {
              // Grab the match header if it exists
              while((symbol = getc(Bp)) != EOF && symbol != '\n' 
                                                  && i+j < ACCSIZE-1)
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
	    fatal("The protein database is not in the correct format");
          }
	  for (psize = 0; ( symbol = getc(Bp)) != EOF && symbol != '>'; )
	       if ( symbol != '\n' )
	         {  P[++psize] = symbol;
                    PN[psize] = ai[symbol];
	         }
	   if ( psize >= wordsize )
	    { PN[psize+1] = AASIZE;
	      cseglen = tupsize = 0;
	      Compare(PN,psize);
	      ChainSeg();
	      for ( i = 1; i <= tupsize; i++ )
	       { if ( ort )
		  { dstart = segment[tuplist[i].first].dstart;
		    dend = segment[tuplist[i].last].dend;
		  }
		 else
		  { dstart = dnalen - segment[tuplist[i].first].dstart + 1;
		    dend = dnalen - segment[tuplist[i].last].dend + 1;
		  }
		 astart = segment[tuplist[i].first].astart;
		 aend = segment[tuplist[i].last].aend;
	         succ = numsp = 0;
		 for ( pred = tuplist[i].last; pred ; numsp++ )
		  { pred = segment[j = pred].pred;
		    segment[j].pred = succ;
		    succ = j;
		  }
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
		 exon_num = 0;
		 for ( ; succ ; succ = segment[succ].pred )
		  { astart = segment[succ].astart - 1;
		    dstart = segment[succ].dstart - 1;
		    aend = segment[succ].aend;
		    dend = segment[succ].dend;
		    segment[succ].inuse = 1;
		    idn = 0;
		    len = aend - astart;
		    asp = segment[succ].asp;
		    dsp = segment[succ].dsp;
		    aep = asp + len;
		    for ( ; asp < aep; )
		      if ( *asp++ == da[*dsp++] )
			 idn++;
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
		        fprintf(alignfp,
                               "\n[Segment %d]\nFrame: %d  Score: %d  Identity: %d/%d (%d%%)\n", 
                               exon_num,
                               segment[succ].fnum * ( ort ? 1 : (-1) ), 
                               score,
		               idn,  
                               len, 
                               (int) ( 0.5 + (100.0 * idn) / len ) 
                              );
	                display2(P+astart,D+dstart,aend-astart,dend-dstart,
		            PN+astart,DN+dstart,astart+1,dstart+1,ort,aa);
                      }
                    }
                  }
	       }
	      /* pick remaining high-scoring segments */
  	      for ( i = 1; i <= cseglen; i++ )
	       if ( ! segment[i].inuse && segment[i].score > fcutoff ) 
	        { if ( ort )
		  { dstart = segment[i].dstart;
		    dend = segment[i].dend;
		  }
		 else
		  { dstart = dnalen - segment[i].dstart + 1;
		    dend = dnalen - segment[i].dend + 1;
		  }
		  astart = segment[i].astart;
		  aend = segment[i].aend;
		  if ( ++numchains < maxchno ) 
                  {
                    global_chain_num++;
                    exon_num = 1;
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
		  astart--;
		  dstart = segment[i].dstart - 1;
		  dend = segment[i].dend;
		  idn = 0;
		  len = aend - astart;
		  asp = segment[i].asp;
		  dsp = segment[i].dsp;
		  aep = asp + len;
		  for ( ; asp < aep; )
		    if ( *asp++ == da[*dsp++] )
		       idn++;
                  if ( numchains < maxchains )
                  { 
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
                              "\n[Segment %d]\nFrame: %d  Score: %d  Identity: %d/%d (%d%%)\n",
                              1,
                              segment[i].fnum * ( ort ? 1 : (-1) ),
                              score,
	                      idn, 
                              len, 
                              (int) ( 0.5 + (100.0 * idn) / len ) 
                             );
	               display2(P+astart,D+dstart,aend-astart,dend-dstart,
		                PN+astart,DN+dstart,astart+1,dstart+1,ort,aa);
                     }
                  }
	        }
	    }
	 }
	(void) fclose(Bp);
}


#define gap(k)  ((k) <= 0 ? q : q+r*(k))	/* k-symbol indel score */

ComputeAiDi()
{ int i;

	for ( i = 0; i < 128; i++ )
	 { di[i] = DNASIZE - 1; 	/* code of 'N' */
	   ai[i] = AASIZE - 1;	 	/* code of 'X' */
	 }
	di['a'] = di['A'] = 0;
	di['c'] = di['C'] = 1;
	di['g'] = di['G'] = 2;
	di['t'] = di['T'] = 3;
	di['n'] = di['N'] = 4;

	ai['a'] = ai['A'] = 0;
	ai['r'] = ai['R'] = 1;
	ai['n'] = ai['N'] = 2;
	ai['d'] = ai['D'] = 3;
	ai['c'] = ai['C'] = 4;
	ai['q'] = ai['Q'] = 5;
	ai['e'] = ai['E'] = 6;
	ai['g'] = ai['G'] = 7;
	ai['h'] = ai['H'] = 8;
	ai['i'] = ai['I'] = 9;
	ai['l'] = ai['L'] = 10;
	ai['k'] = ai['K'] = 11;
	ai['m'] = ai['M'] = 12;
	ai['f'] = ai['F'] = 13;
	ai['p'] = ai['P'] = 14;
	ai['s'] = ai['S'] = 15;
	ai['t'] = ai['T'] = 16;
	ai['w'] = ai['W'] = 17;
	ai['y'] = ai['Y'] = 18;
	ai['v'] = ai['V'] = 19;
	ai['b'] = ai['B'] = 20;
	ai['z'] = ai['Z'] = 21;
	ai['x'] = ai['X'] = 22;
}

/* Compute the da table */
ComputeDa()
{ int  i;	/* index variable */
  int  b;	/* one-symbol integer code */
  int  bb;	/* two-symbol integer code */
	
	dnals3 = (dnals2 = DNAS2) * DNASIZE;
	for ( i = 0; i < dnals3; i++ )
	  da[i] = AASIZE - 1;
	b = di['A'] * dnals2;
	bb = b + di['A'] * DNASIZE;
	da[bb + di['A']] = da[bb + di['G']] = 11;
	da[bb + di['C']] = da[bb + di['T']] = 2;
	bb = b + di['C'] * DNASIZE;
	da[bb + di['A']] = da[bb + di['G']] = 16;
	da[bb + di['C']] = da[bb + di['T']] = 16;
	bb = b + di['G'] * DNASIZE;
	da[bb + di['A']] = da[bb + di['G']] = 1;
	da[bb + di['C']] = da[bb + di['T']] = 15;
	bb = b + di['T'] * DNASIZE;
	da[bb + di['A']] = da[bb + di['C']] = da[bb + di['T']] = 9;
	da[bb + di['G']] = 12;

	b = di['C'] * dnals2;
	bb = b + di['A'] * DNASIZE;
	da[bb + di['A']] = da[bb + di['G']] = 5;
	da[bb + di['C']] = da[bb + di['T']] = 8;
	bb = b + di['C'] * DNASIZE;
	da[bb + di['A']] = da[bb + di['G']] = 14;
	da[bb + di['C']] = da[bb + di['T']] = 14;
	bb = b + di['G'] * DNASIZE;
	da[bb + di['A']] = da[bb + di['G']] = 1;
	da[bb + di['C']] = da[bb + di['T']] = 1;
	bb = b + di['T'] * DNASIZE;
	da[bb + di['A']] = da[bb + di['G']] = 10;
	da[bb + di['C']] = da[bb + di['T']] = 10;

	b = di['G'] * dnals2;
	bb = b + di['A'] * DNASIZE;
	da[bb + di['A']] = da[bb + di['G']] = 6;
	da[bb + di['C']] = da[bb + di['T']] = 3;
	bb = b + di['C'] * DNASIZE;
	da[bb + di['A']] = da[bb + di['G']] = 0;
	da[bb + di['C']] = da[bb + di['T']] = 0;
	bb = b + di['G'] * DNASIZE;
	da[bb + di['A']] = da[bb + di['G']] = 7;
	da[bb + di['C']] = da[bb + di['T']] = 7;
	bb = b + di['T'] * DNASIZE;
	da[bb + di['A']] = da[bb + di['G']] = 19;
	da[bb + di['C']] = da[bb + di['T']] = 19;

	b = di['T'] * dnals2;
	bb = b + di['A'] * DNASIZE;
	da[ochre = bb + di['A']] = da[amber = bb + di['G']] = AASIZE;
	da[bb + di['C']] = da[bb + di['T']] = 18;
	bb = b + di['C'] * DNASIZE;
	da[bb + di['A']] = da[bb + di['G']] = 15;
	da[bb + di['C']] = da[bb + di['T']] = 15;
	bb = b + di['G'] * DNASIZE;
	da[uga = bb + di['A']] = AASIZE;
	da[bb + di['G']] = 17;
	da[bb + di['C']] = da[bb + di['T']] = 4;
	bb = b + di['T'] * DNASIZE;
	da[bb + di['A']] = da[bb + di['G']] = 10;
	da[bb + di['C']] = da[bb + di['T']] = 13;
	gtcode = di['G'] * DNASIZE + di['T'];
	agcode = di['A'] * DNASIZE + di['G'];
}

/* Set up a lookup table for words of length wordsize in 3 frame.
   The value of a word is used as an index to the table.  */
SetTable(frame, size, fnum)
int *frame;		/* sequence of codon codes */
int size;		/* length of codon codes */
int fnum;		/* frame number */
{ int   value;			/* value of a word */
  int   len;			/* number of codons after the last one for 'X' */
  int   p;			/* index variables */
  register int  i;		/* index variables */
  int   xcode;			/* integer code of amino acid X */

	/* make a lookup table */
	xcode = AASIZE - 1;
	for ( value = len = 0, i = 1; i <= size; i++ )
	 { if ( ( p = da[frame[i]] ) >= xcode )
	      len = -1;
	   value = ( value * AASIZE + p ) % tsize;
	   if ( ++len >= wordsize )
	     { list[p = (i - wordsize) * 3 + fnum] = table[value];
	       table[value] = p;
	     }
	 }
}

/* Set up w */
SetUpW()
{ int  i, j, k;	/* index varialbes */
  int  row;	/* row number of w */
  int  *wa;	/* row of w */
  int  *va;	/* row of v */
  int  b;	/* one-symbol integer code */
  int  bb;	/* two-symbol integer code */
  int  bbb;	/* three-symbol integer code */
  int  ci;	/* number of codons in the i-indexed loop */
  int  si;	/* sum of substitution scores in the i-indexed loop */
  int  cj;	/* number of codons in the j-indexed loop */
  int  sj;	/* sum of substitution scores in the j-indexed loop */
  int  ck;	/* number of codons in the k-indexed loop */
  int  sk;	/* sum of substitution scores in the k-indexed loop */
  int  s;	/* substitution score */
  int  ds1;	/* DNASIZE -1 */
  int  NXX, XNX, NNX, NXN, XNN;	/* constants */	
	

	ds1 = DNASIZE - 1;
	NXX = ds1 * dnals2;
	XNX = ds1 * DNASIZE;
	NNX = NXX + XNX;
	NXN = NXX + ds1; 
	XNN = XNX + ds1;
	/* compute the substitution scores for ijk, ij'N' and i'N''N' */
	for ( row = 0; row < AASIZE; row++ )
	 { wa = w[row];
	   va = pam[row];
	   for ( i = 0; i < ds1; i++ )
	    { b = i * DNASIZE;
	      for ( sj = cj = j = 0; j < ds1; j++ )
	       { bb = (b + j) * DNASIZE;
		 for ( sk = ck = k = 0; k < ds1; k++ )
		  { bbb = bb + k;
		    if ( da[bbb] != AASIZE )
		     { s = wa[bbb] = va[da[bbb]];
		       sk += s;
		       sj += s;
		       cj++;
		       ck++;
		     }
		  }
		 bbb = bb + ds1;
		 wa[bbb] = ck > 0 ? (sk / ck) : 0;
	       }
              bbb = (b + ds1) * DNASIZE + ds1;
              wa[bbb] = cj > 0 ? (sj / cj) : 0;
	    }
	 }
	for ( row = 0; row < AASIZE; row++ )
	 { wa = w[row];
	   /* compute the substitution scores for 'N'jk and 'N''N'k */
	   for ( k = 0; k < ds1; k++ )
	    { for ( sj = cj = j = 0; j < ds1; j++ )
	       { bb = j * DNASIZE + k;
		 for ( si = ci = i = 0; i < ds1; i++ )
		  { bbb = i * dnals2 + bb;
		    if ( da[bbb] != AASIZE )
		     { si += (s = wa[bbb]);
		       sj += s;
		       ci++;
		       cj++;
		     }
		  }
		 bbb = NXX + bb;
		 wa[bbb] = ci > 0 ? (si / ci) : 0;
	       }
              bbb = NNX + k;
              wa[bbb] = cj > 0 ? (sj / cj) : 0;
	    }
	   /* compute the substitution scores for 'N'j'N' */
	   for ( j = 0; j < ds1; j++ )
	    { b = j * DNASIZE;
	      for ( si = ci = i = 0; i < ds1; i++ )
	       { bb = i * dnals2 + b;
		 for ( k = 0; k < ds1; k++ )
		  { bbb = bb + k;
		    if ( da[bbb] != AASIZE )
		     { si += wa[bbb];
		       ci++;
		     }
		  }
	       }
              bbb = NXN + b;
              wa[bbb] = ci > 0 ? (si / ci) : 0;
	    }
	   /* compute the substitution scores for i'N'k */
	   for ( i = 0; i < ds1; i++ )
	    { b = i * dnals2;
	      for ( k = 0; k < ds1; k++ )
	       { bb = b + k;
		 for ( sj = cj = j = 0; j < ds1; j++ )
		  { bbb = bb + j * DNASIZE;
		    if ( da[bbb] != AASIZE )
		     { sj += wa[bbb];
		       cj++;
		     }
		  }
                 bbb = b + XNX + k;
                 wa[bbb] = cj > 0 ? (sj / cj) : 0;
	       }
	    }
           wa[NXX + XNN] = -1;
	   /* compute the substitution scores for stop codons */
           wa[ochre] = mismh;
           wa[amber] = mismh;
           wa[uga] = mismh;
           wa[dnals3] = MINF;
	 }
	wa = w[AASIZE];
	for ( i = 0; i <= dnals3; i++ )
	  wa[i] = MINF;
}

/* Set up mt1, mt2 and mt3. */
SetUpMt()
{ int  i, j, k;	/* index varialbes */
  int  row;	/* row number of w */
  char *mr1;	/* row of mt1 */
  char *mr2;	/* row of mt2 */
  char *mr3;	/* row of mt3 */
  int  b;	/* one-symbol integer code */
  int  bb;	/* two-symbol integer code */
  int  bbb;	/* three-symbol integer code */
	
	/* set all mt to a blank */
	for ( row = 0; row < AASIZE; row++ )
	 { mr1 = mt1[row];
	   mr2 = mt2[row];
	   mr3 = mt3[row];
	   for ( i = 0; i < DNASIZE; i++ )
	    { mr1[i] = ' ';
	      b = i * DNASIZE;
	      for ( j = 0; j < DNASIZE; j++ )
	       { mr2[bb = b+j] = ' ';
		 bbb = bb * DNASIZE;
		 for ( k = 0; k < DNASIZE; k++ )
	            mr3[bbb + k] = ' ';
	       }
	    }
	 }
	mt1[ai['A']][di['G']] = '.';
	mt1[ai['R']][di['A']] = '.';
	mt1[ai['R']][di['C']] = '.';
	mt1[ai['N']][di['A']] = '.';
	mt1[ai['D']][di['G']] = '.';
	mt1[ai['C']][di['T']] = '.';
	mt1[ai['Q']][di['C']] = '.';
	mt1[ai['E']][di['G']] = '.';
	mt1[ai['G']][di['G']] = '.';
	mt1[ai['H']][di['C']] = '.';
	mt1[ai['I']][di['A']] = '.';
	mt1[ai['L']][di['C']] = '.';
	mt1[ai['L']][di['T']] = '.';
	mt1[ai['K']][di['A']] = '.';
	mt1[ai['M']][di['A']] = '.';
	mt1[ai['F']][di['T']] = '.';
	mt1[ai['P']][di['C']] = '.';
	mt1[ai['S']][di['A']] = '.';
	mt1[ai['S']][di['T']] = '.';
	mt1[ai['T']][di['A']] = '.';
	mt1[ai['W']][di['T']] = '.';
	mt1[ai['Y']][di['T']] = '.';
	mt1[ai['V']][di['G']] = '.';
	mt1[ai['B']][di['A']] = '.';
	mt1[ai['B']][di['G']] = '.';
	mt1[ai['Z']][di['C']] = '.';
	mt1[ai['Z']][di['G']] = '.';
	for ( i = 0; i < DNASIZE; i++ )
	 { b = i * DNASIZE;
	   mt2[ai['A']][b + di['C']] = '.';
	   mt2[ai['R']][b + di['G']] = '.';
	   mt2[ai['N']][b + di['A']] = '.';
	   mt2[ai['D']][b + di['A']] = '.';
	   mt2[ai['C']][b + di['G']] = '.';
	   mt2[ai['Q']][b + di['A']] = '.';
	   mt2[ai['E']][b + di['A']] = '.';
	   mt2[ai['G']][b + di['G']] = '.';
	   mt2[ai['H']][b + di['A']] = '.';
	   mt2[ai['I']][b + di['T']] = '.';
	   mt2[ai['L']][b + di['T']] = '.';
	   mt2[ai['K']][b + di['A']] = '.';
	   mt2[ai['M']][b + di['T']] = '.';
	   mt2[ai['F']][b + di['T']] = '.';
	   mt2[ai['P']][b + di['C']] = '.';
	   if ( i != di['A'] )
	     mt2[ai['S']][b + di['C']] = '.';
	   if ( i != di['T'] )
	     mt2[ai['S']][b + di['G']] = '.';
	   mt2[ai['T']][b + di['C']] = '.';
	   mt2[ai['W']][b + di['G']] = '.';
	   mt2[ai['Y']][b + di['A']] = '.';
	   mt2[ai['V']][b + di['T']] = '.';
	   mt2[ai['B']][b + di['A']] = '.';
	   mt2[ai['Z']][b + di['A']] = '.';
	   for ( j = 0; j < DNASIZE; j++ )
	    { bb = (b + j) * DNASIZE;
	      mr3 = mt3[ai['A']];
	      mr3[bb + di['A']] = mr3[bb + di['C']] = mr3[bb + di['G']] = '.';
	      mr3[bb + di['T']] = mr3[bb + di['N']] = '.';
	      mr3 = mt3[ai['R']];
	      mr3[bb + di['A']] = mr3[bb + di['G']] = '.';
	      if ( i != di['A'] )
		mr3[bb + di['C']] = mr3[bb + di['T']] = mr3[bb + di['N']] = '.';
	      mr3 = mt3[ai['N']];
	      mr3[bb + di['C']] = mr3[bb + di['T']] = '.';
	      mr3 = mt3[ai['D']];
	      mr3[bb + di['C']] = mr3[bb + di['T']] = '.';
	      mr3 = mt3[ai['C']];
	      mr3[bb + di['C']] = mr3[bb + di['T']] = '.';
	      mr3 = mt3[ai['Q']];
	      mr3[bb + di['A']] = mr3[bb + di['G']] = '.';
	      mr3 = mt3[ai['E']];
	      mr3[bb + di['A']] = mr3[bb + di['G']] = '.';
	      mr3 = mt3[ai['G']];
	      mr3[bb + di['A']] = mr3[bb + di['C']] = mr3[bb + di['G']] = '.';
	      mr3[bb + di['T']] = mr3[bb + di['N']] = '.';
	      mr3 = mt3[ai['H']];
	      mr3[bb + di['C']] = mr3[bb + di['T']] = '.';
	      mr3 = mt3[ai['I']];
	      mr3[bb + di['A']] = mr3[bb + di['C']] = mr3[bb + di['T']] = '.';
	      mr3 = mt3[ai['L']];
	      mr3[bb + di['A']] = mr3[bb + di['G']] = '.';
	      if ( i != di['T'] )
		mr3[bb + di['C']] = mr3[bb + di['T']] = mr3[bb + di['N']] = '.';
	      mr3 = mt3[ai['K']];
	      mr3[bb + di['A']] = mr3[bb + di['G']] = '.';
	      mr3 = mt3[ai['M']];
	      mr3[bb + di['G']] = '.';
	      mr3 = mt3[ai['F']];
	      mr3[bb + di['C']] = mr3[bb + di['T']] = '.';
	      mr3 = mt3[ai['P']];
	      mr3[bb + di['A']] = mr3[bb + di['C']] = mr3[bb + di['G']] = '.';
	      mr3[bb + di['T']] = mr3[bb + di['N']] = '.';
	      mr3 = mt3[ai['S']];
	      mr3[bb + di['C']] = mr3[bb + di['T']] = '.';
	      if ( i != di['A'] && ( i == di['T'] || j != di['G'] ) )
	       mr3[bb + di['A']] = mr3[bb + di['G']] = mr3[bb + di['N']] = '.';
	      mr3 = mt3[ai['T']];
	      mr3[bb + di['A']] = mr3[bb + di['C']] = mr3[bb + di['G']] = '.';
	      mr3[bb + di['T']] = mr3[bb + di['N']] = '.';
	      mr3 = mt3[ai['W']];
	      mr3[bb + di['G']] = '.';
	      mr3 = mt3[ai['Y']];
	      mr3[bb + di['C']] = mr3[bb + di['T']] = '.';
	      mr3 = mt3[ai['V']];
	      mr3[bb + di['A']] = mr3[bb + di['C']] = mr3[bb + di['G']] = '.';
	      mr3[bb + di['T']] = mr3[bb + di['N']] = '.';
	      mr3 = mt3[ai['B']];
	      mr3[bb + di['C']] = mr3[bb + di['T']] = '.';
	      mr3 = mt3[ai['Z']];
	      mr3[bb + di['A']] = mr3[bb + di['G']] = '.';
	    }
	 }
}

ComputeAa()
{  
  	strcpy(aa[0], "Ala");
	strcpy(aa[1], "Arg");
	strcpy(aa[2], "Asn");
	strcpy(aa[3], "Asp");
	strcpy(aa[4], "Cys");
	strcpy(aa[5], "Gln");
	strcpy(aa[6], "Glu");
	strcpy(aa[7], "Gly");
	strcpy(aa[8], "His");
	strcpy(aa[9], "Ile");
	strcpy(aa[10], "Leu");
	strcpy(aa[11], "Lys");
	strcpy(aa[12], "Met");
	strcpy(aa[13], "Phe");
	strcpy(aa[14], "Pro");
	strcpy(aa[15], "Ser");
	strcpy(aa[16], "Thr");
	strcpy(aa[17], "Trp");
	strcpy(aa[18], "Tyr");
	strcpy(aa[19], "Val");
	strcpy(aa[20], "Asx");
	strcpy(aa[21], "Glx");
	strcpy(aa[22], "Xxx");

  	strcpy(aa2[0], " A ");
	strcpy(aa2[1], " R ");
	strcpy(aa2[2], " N ");
	strcpy(aa2[3], " D ");
	strcpy(aa2[4], " C ");
	strcpy(aa2[5], " Q ");
	strcpy(aa2[6], " E ");
	strcpy(aa2[7], " G ");
	strcpy(aa2[8], " H ");
	strcpy(aa2[9], " I ");
	strcpy(aa2[10], " L ");
	strcpy(aa2[11], " K ");
	strcpy(aa2[12], " M ");
	strcpy(aa2[13], " F ");
	strcpy(aa2[14], " P ");
	strcpy(aa2[15], " S ");
	strcpy(aa2[16], " T ");
	strcpy(aa2[17], " W ");
	strcpy(aa2[18], " Y ");
	strcpy(aa2[19], " V ");
	strcpy(aa2[20], " B ");
	strcpy(aa2[21], " Z ");
	strcpy(aa2[22], " X ");
}

/* computes quickly scores of alignments between the two sequences */
Compare(PN,psize)
int  PN[], psize;
{ register int  value;		/* value of a word */
  register int  i, j, k, h, f, g; /* index variables */
  register int  *aap;		/* pointer to element of PN */
  register int  *codp;		/* pointer to element of frame */
  register int  maxs;		/* maximum score seen */
  register int  score;		/* current score */
  register int  s;		/* initial score */
  register int  *asp, *dsp, *aep, *dep;	/* start and end pointers */
  int  num;		/* number of segment pairs in segment */
  int  temp1, temp2;	/* temporary variables */
  int  *frame;		/* current frame */
  int  fnum;		/* current frame number */

	num = cseglen;
	for ( s = value = 0, i = 1; i < wordsize; i++ )
	  { value = value * AASIZE + (j = PN[i] );
	    s += pam[j][j];
	  }
	for ( ; i <= psize ; i++ )
	 { value = ( value * AASIZE + (j = PN[i] ) ) % tsize;
	   k = PN[i-wordsize];
	   s += pam[j][j] - pam[k][k];
	   for ( g = table[value]; g ; g = list[g] )
	     { if ( (h = diagonal[k = g - 3 * i + aalen3]) && segment[h].aend >= i )
		    continue;
	       frame = framept[fnum = g % 3];
	       if ( !fnum) fnum = 3;
	       j = ( g - 1 ) / 3 + 1;
	       aep = aap = PN + i + 1;
	       dep = codp = frame + j + wordsize;
	       for ( maxs = score = s; maxs - score < distance; )
		{ score += w[*aap++][*codp++];
		  if ( maxs < score )
		   { maxs = score;
		     aep = aap; 
		     dep = codp;
		   }
		}
	       asp = aap = PN + i - wordsize;
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
		{ if ( h && segment[h].astart > (asp - PN) )
		     f = h;
		  else
		   { if ( ( f = ++num ) > segsize )
                       DoubleSegSpace();
		   }
		  if ( f <= segsize )
		    { segment[f].astart = temp1 = ++asp - PN;
		      segment[f].asp = asp;
		      segment[f].dstart = temp2 = ( dsp - frame ) * 3 + fnum;
		      segment[f].dsp = dsp + 1;
		      segment[f].antis = 3 * temp1 + temp2;
		      segment[f].aend = temp1 =  --aep - PN;
		      segment[f].dend = temp2 = (--dep - frame ) * 3 + fnum - 1;
		      segment[f].antid = 3 * temp1 + temp2;
		      segment[f].score = segment[f].total = maxs;
		      segment[f].fnum = fnum;
		      segment[f].inuse = segment[f].pred = 0;
		      segment[f].first = f;
		      segment[f].diag = k;
		      segment[f].kind = 0;
		      diagonal[k] = f;
		    }
		}
	     }
	 }
	if ( num > segsize )
	   num = segsize;
	for ( i = cseglen + 1; i <= num; i++ )
	  diagonal[segment[i].diag] = 0;
	cseglen = num;
}

/* It sorts the segments by antidiagonals and chains close segments. */
ChainSeg()
{ int  i, j, d;		/* index variables */
  int  found;		/* flag */
  int  recent;		/* index of most recently visited element */
  int  temp, ss;	/* temporary variables */
  int  first;		/* start segment in a chain */
  int  total;		/* score of the chain */
  segptr  p;
  int  *asp, *dsp;	/* pointers to PN and frame */
  int  iscore;		/* segstd[i]->score */
  int  is3;		/* score of 3' portion of a segment */
  int  pred;		/* link to the previous segment pair */
  int  last;		/* last segment pair in a chain */

	/* make an array of pointers and sort the array */
	p = segment;
	for ( i = 1; i <= cseglen; )
	  segstd[i++] = ++p;
	QuickSort(1, cseglen);
      	recent = 0;
	/* chain segments */
	for ( i = 1; i <= cseglen; i++ )
	 { if ( (j = i - 1) && segstd[i]->antis - segstd[j]->antid < antidis )
	    { asp = segstd[i]->asp;
	      dsp = segstd[i]->dsp;
	      total = segstd[i]->aend - segstd[i]->astart;
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
	     if ( (temp = segstd[i]->astart - segstd[j]->aend) > aaover &&
		  (d = (segstd[i]->dstart - segstd[j]->dend) / 3) > aaover &&
		  ( ( temp > 0 && d > 0 ) ||
		   ( is3 = iscore - fdst[d < temp ? (-d) : (-temp)]) > icutoff ) &&
		 (ss = segstd[j]->total + is3 - gap(temp-1)) > segstd[i]->total )
              { segstd[i]->first = segstd[j]->first;
                segstd[i]->pred = segstd[j] - segment;
	        segstd[i]->total = ss;
	        segstd[i]->kind = 0;
	      }
	     else
	      if ( segstd[i]->dstart - segstd[j]->dend > intergenedis &&
	           segstd[j]->total > fcutoff &&
		   (ss = segstd[j]->total + iscore -
		    gap(segstd[i]->astart + segstd[j]->aend - 2)) > segstd[i]->total )
               { segstd[i]->first = segstd[j]->first;
                 segstd[i]->pred = segstd[j] - segment;
	         segstd[i]->total = ss;
	         segstd[i]->kind = 1;
	       }         /* build a mega-chain */
	    }
	   if ( (total = segstd[i]->total) > fcutoff )
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

/* break mega-chains into chains */
	recent = tupsize;
	for ( j = 1; j <= recent; j++ )
	 { for ( pred = tuplist[j].last; pred; )
	    { pred = segment[d = pred].pred;
	      if ( segment[d].kind )
	       { last = tuplist[j].last;
	         if ( ! pred ) fatal("ChainSeg: no previous link");
	         total = segment[last].total - segment[pred].total +
		        (is3 = gap(segment[pred].aend + segment[d].astart - 2) );
	         if ( total > fcutoff )
		  { tuplist[++tupsize].total = total;
		    tuplist[tupsize].first = d;
		    tuplist[tupsize].last = last;
		  }
	         for ( ss = last; ss; ss = segment[ss].pred )
		  { segment[ss].first = d;
		    segment[ss].total -= segment[pred].total - is3;
		    if ( ss == d ) break;
		  }
	         if ( ss != d ) fatal("ChainSeg: wrong link");
		 tuplist[j].last = pred;
		 tuplist[j].total = segment[pred].total;
		 segment[d].pred = 0;
		 segment[d].kind = 0;
	       }      /*  an inter-chain link a found */
	    }         /* for each link */
	 }            /* for each mega-chian */
}

/* It sorts the elements between positions p and r. */
QuickSort(p, r)
int  p, r; /* start and end positions */
{ int  i, j;		/* index variables */
  int  sorted;		/* flag */
  segptr  t;		/* pointer */
  int  pivot;	        /* pivot element */

     if ( p >= r ) return;
     if ( p + 20 > r )
      { for ( i = p; i < r; i++ )
	 { sorted = 1;
	   for ( j = r; j > i; j--)
	     if ( segstd[j-1]->antid > segstd[j]->antid )
	      { t = segstd[j];
		segstd[j] = segstd[j-1];
		segstd[j-1] = t;
		sorted = 0;
	      }
	   if ( sorted )
	      break;
	 }
	return;
      }
    pivot = segstd[p]->antid;
    i = p - 1;
    j = r + 1;
    while ( 1 )
     { for ( j--; segstd[j]->antid > pivot; j-- )
	  ;
       for ( i++; segstd[i]->antid < pivot; i++ )
	  ;
       if ( i < j )
	{ t = segstd[i];
	  segstd[i] = segstd[j];
	  segstd[j] = t;
        }
       else
	 break;
     }
    QuickSort(p, j);
    QuickSort(j+1, r);
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
           one->inuse = two->inuse;
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
	 }
 (void) free(oldsegment);
 (void) free(oldsegstd);
 (void) free(oldtuplist);
}

/* Alignment display routine */

static char ALINE[65];	/* protein line */
static char DLINE[65];	/* DNA line */
static char CLINE[65];	/* pair type indicator */ 

display2(A,B,M,N,AN,DN,AP,DP,orient,ttt)
char A[], B[], (*ttt)[4]; int M, N, AN[], DN[]; int AP, DP, orient;
{ register char *a;	/* protein line pointer */ 
  register char *d;	/* DNA line pointer */
  register char *c;	/* center line pointer */
  register char *e;	/* temp pointer */
  register int  i, j;	/* index variables */
           int  dp;	/* DNA sequence position */
           int  ap;	/* protein sequence position */
	   char *letter;/* pointer to three letter aa string */
	   int  bbb;	/* three-symbol integer code */
   int mdp;	/* the DNA position relative to the plus strand */
   int y, z;	/* current nucleotide integer code */
   int one;	/* current amino acid integer code */
   int  dp2;	/* DNA sequence position for right column */
   int  ap2;	/* protein sequence position for right column */
   int  mdp2;	/* the DNA position relative to the plus strand */

  i = j = 0;
  ap = AP;
  dp = DP;
  a = ALINE;
  d = DLINE;
  c = CLINE;
  while (i < M || j < N)
    { letter = ttt[one = AN[++i]];
      *d++ = B[++j];
      y = DN[j];
      *d++ = B[++j];
      z = DN[j];
      *d++ = B[++j];
      bbb = (z = y * DNASIZE + z) * DNASIZE + DN[j];
      *a++ = *letter++;
      *a++ = *letter++;
      *a++ = *letter++;
      if ( AN[i] == da[bbb] )
       { *c++ = ':';
         *c++ = ':';
         *c++ = ':';
       }
      else
       if ( da[bbb] == AASIZE )
	{ *c++ = '*';
          *c++ = '*';
          *c++ = '*';
	}
       else
        { *c++ = mt1[one][y];
          *c++ = mt2[one][z];
          *c++ = mt3[one][bbb];
        }
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
