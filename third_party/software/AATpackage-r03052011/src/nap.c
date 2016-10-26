/*  A global Nucleotide-Amino acid alignment Program (NAP):
    
This is a modified version of NAP for processing the output of DPS.
It also produces a Btab output. 

Proper attribution of the author as the source of the software would
be appreciated:
Huang, X. and Zhang, J.  (1996)
Methods for comparing a DNA sequence with a protein sequence,
CABIOS 12 (6), 497-506.

Acknowledgments
I thank the following people for discussions and suggestions:
Mark Adams, Tony Kerlavage, Brendan Loftus, Steve Rounsley,
Jinghui Zhang, and lixin Zhou.

The NAP program computes an optimal global alignment of a DNA sequence
and a protein sequence without penalizing terminal gaps.
NAP handles frameshifts and long introns in the DNA sequence.
It delivers the alignment in linear space, so long sequences can be aligned.
If DOUBLE = 1, then both strands of the DNA sequence are compared
with the protein sequence and one of the two alignments with the larger
score is reported. If DOUBLE = 0, only the plus strand is compared
with the protein sequence.

The NAP program is written in C and runs under Unix systems on
Sun workstations and under DOS systems on PCs.
We think that the program is portable to many machines.

Sequences to be analyzed are stored in separate files.
An input file contains all characters of a sequence, separated by
newline characters, in linear order. No other characters are allowed.
Sample DNA and protein sequence files are attached below.

To find the best alignment of a DNA sequence in file DNA and a protein
sequence in file Protein, use a command of form

nap  DNA  Protein  gs  BLOSUM62  q  r > result

where nap is the name of the object code, gs is the minimum length
of any deletion gap in the DNA sequence receiving a constant gap penalty,
BLOSUM62 is a specially formatted BLOSUM62 matrix, q and r are
non-negative integers specifying gap-open and gap-extend penalties,
respectively. The score of an i-symbol gap is -(q + r * i).
Any deletion gap of length greater than gs is given a score of
-(q + r * gs). The output alignment is saved in the file "result".
A set of possible values for the parameters are: gs = 10, q = 10, r = 2.

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
GCCAGTGCCAGAAGAGCCAAGGACAGGTACGGCTGTCATCACTTAGACCT
CACCCTGTGGAGCCACACCCTAGGGTTGGCCAATCTACTCCCAGGAGCAG
GGAGGGCAGGAGCCAGGGCTGGGCATAAAAGTCAGGGCAGAGCCATCTAT
TGCTTACATTTGCTTCTGACACAACTGTGTTCACTAGCAACCTCAAACAG
ACACCATGGTGCACCTGATCTGAGGAGAAGTCTGCCGTTACTGCCCTG
TGGGGCAAGGTGATACGTGGATGAAGTTGGTGGTGATATGGCCCTGGGCAGGTT
GGTATCAAGGTTACAAGACAGGTTTAAGGAGACCAATAGAAACTGGGCAT
GTGGAGACAGAGAAGACTCTTGGGTTTCTGATAGGCACTGACTCTCTCTG
CCTATTGGTCTATTTTCCCACCCTTAGGCTGCTGGTGGTCTACCCTTGGA
CCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTT
ATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTT
TAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACAC
TGAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGG
GTGAGTCTATGGGACCCTTGATGTTTTCTTTCCCCTTCTTTTCTATGGTT
AAGTTCATGTCATAGGAAGGGGAGAAGTAACAGGGTACAGTTTAGAATGG
GAAACAGACGAATGATTGCATCAGTGTGGAAGTCTCAGGATCGTTTTAGT
TTCTTTTATTTGCTGTTCATAACAATTGTTTTCTTTTGTTTAATTCTTGC
TTTCTTTTTTTTTCTTCTCCGCAATTTTTACTATTATACTTAATGCCTTA
ACATTGTGTATAACAAAAGGAAATATCTCTGAGATACATTAAGTAACTTA
AAAAAAAACTTTACACAGTCTGCCTAGTACATTACTATTTGGAATATATG
TGTGCTTATTTGCATATTCATAATCTCCCTACTTTATTTTCTTTTATTTT
TAATTGATACATAATCATTATACATATTTATGGGTTAAAGTGTAATGTTT
TAATATGTGTACACATATTGACCAAATCAGGGTAATTTTGCATTTGTAAT
TTTAAAAAATGCTTTCTTCTTTTAATATACTTTTTTGTTTATCTTATTTC
TAATACTTTCCCTAATCTCTTTCTTTCAGGGCAATAATGATACAATGTAT
CATGCCTCTTTGCACCATTCTAAAGAATAACAGTGATAATTTCTGGGTTA
AGGCAATAGCAATATTTCTGCATATAAATATTTCTGCATATAAATTGTAA
CTGATGTAAGAGGTTTCATATTGCTAATAGCAGCTACAATCCAGCTACCA
TTCTGCTTTTATTTTATGGTTGGGATAAGGCTGGATTATTCTGAGTCCAA
GCTAGGCCCTTTTGCTAATCATGTTCATACCTCTTATCTTCCTCCCACAG
CTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGA
ATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGG
CTAATGCCCTGGCCCACAAGTATCACTAAGCTCGCTTTCTTGCTGTCCAA
TTTCTATTAAAGGTTCCTTTGTTCCCTAAGTCCAACTACTAAACTGGGGG
ATATTATGAAGGGCCTTGAGCATCTGGATTCTGCCTAATAAAAAACATTT
ATTTTCATTGCAATGATGTATTTAAATTATTTCTGAATATTTTACTAAAA
AGGGAATGTGGGAGGTCAGTGCATTTAAAACATAAAGAAATGAAGAGCTA
GTTCAAACCTTGGGAAAATACACTATATCTTAAACTCCATGAAAGAAGGT
GAGGCTGCAAACAGCTAATGCACATTGGCAACAGCCCTGATGCCTATGCC
TTATTCATCCCTCAGAAAAGGATTCAAGTAG

A sample protein file:
VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKV
KAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENHFRLLGNVLVCVLAHHFGK
EFTPPVQAAYKVVAGVANALAHKYH
*/

#include   <stdio.h>

#define  DNASIZE       5	/* size of DNA alphabet: A, C, G, T, N */
#define  AASIZE       23	/* size of protein alphabet */
#define  DNAS2  (DNASIZE * DNASIZE) /* DNASIZE square */
#define  DNAS3  (DNASIZE * DNAS2)   /* DNASIZE cube */
#define  DOUBLE        1	/* 1, both strands; 0, given strand only */
#define  MINF      -10000	/* negative infinity */
#define  BONUS5    5	/* a bonus for a long gap starting with GT */
#define  BONUS3    5	/* a bonus for a long gap ending with AG */
#define  ACCLEN    400          /* maximum length of protein accession information */
#define  EXTRA     50           /* length of extra DNA region used */
#define  CSLEN     25           /* max length of 5' and 3' end of an exon examined */
#define  MINEXON    7           /* minimum exon length, must be < CSLEN */
#define  QGAPLEN   25	        /* default value for gaplen */
#define  GAPOPEN   15	        /* default value for q */
#define  GAPEXTN    1	        /* default value for r */
#define  LINELEN   60	        /* length of an output line */

// Version info
static float VersionS = 1.51f;
static char * BuildS = "$Revision: 1.9 $";

static char   aa[AASIZE][4];	/* integer code to three-letter code */
static char   aa2[AASIZE][4];	/* integer code to one-letter code with two blanks */
static char   aa3[AASIZE+1];	/* integer code to one-letter code */
static int    w[AASIZE][DNAS3]; /* score table for substitutions */
static int    da[DNAS3];	/* codon integer code to aa integer code */
static int    ai[128];		/* aa letter to integer code */
static int    di[128];		/* DNA letter to integer code */
static int    pam[AASIZE][AASIZE]; /* pam matrix */
static int    nnb[DNASIZE];	/* type 1: b to nnb */
static int    nbn[DNASIZE];	/* type 2: b to nbn */
static int    bnn[DNASIZE];	/* type 3: b to bnn */
static int    nbb[DNAS2];	/* type 4: bb to nbb */
static int    bnb[DNAS2];	/* type 5: bb to bnb */
static int    bbn[DNAS2];	/* type 6: bb to bbn */
/* type 7, bbb; 8, bnbb; 9, bbnb; 10 and 11, bbn...nb; 12 and 13, bn...nbb */
static int    ochre;		/* stop codon code */
static int    amber;		/* stop codon code */
static int    uga;		/* stop codon code */
static int    dnals2;		/* DNASIZE square */
static int    dnals3;		/* DNASIZE cube */
static int    gtcode;		/* the integer code of GT */
static int    agcode;		/* the integer code of AG */
static int    bonus5;		/* BONUS5 */
static int    bonus3;		/* BONUS3 */
static char   mt1[AASIZE][DNASIZE]; /* patial match indicator for base 1 */
static char   mt2[AASIZE][DNAS2];   /* patial match indicator for base 2 */
static char   mt3[AASIZE][DNAS3];   /* patial match indicator for base 3 */
static int  dnalen;		/* query sequence length */
static int  psize;		/* current protein length */
static int  minexon = MINEXON; /* min exon length */

/* Btab information */
static char *dhead;	        /* the name of the genomic DNA */
static char btab_date[100];	/* search date */
static char *program;           /* program name */
static char dblocus[ACCLEN+1];  /* locus name of database */
static int  dbstart;		/* start position of DB sequence */
static int  dbend;		/* end position of DB sequence */
static int  idnno;		/* number of identities */
static int  simno;		/* number of similarities */
static int  exlen;		/* length of an exon alignment */
static int  exscore;		/* length of an exon alignment */
static char dbdesp[ACCLEN+1];	/* DB sequence description */
static FILE *Btabp;		/* Btab output file pointer */
static char fname[ACCLEN+1];	/* Btab file name */

static int match, mismh;		/* max and min substitution weights */
static char *name1, *name2;		/* names of sequence files    */
static int  q, r;       /* gap penalties */
static int  qr;         /* qr = q + r */
static int  qr2;        /* q + 2 * r */
static int  qr3;        /* q + 3 * r */
static int  q2r2;       /* 2 * q + 2 * r */
static int  r2;         /* 2 * r */
static int  r3;         /* 3 * r */
static int  gaplen;     /* minimum length for constant-cost insertion */
static int  pay;	/* constant-cost for long insertion */

static int *CC, *PC, *DD;		/* saving matrix scores */
static int *RR, *SS;		 	/* for reverse pass */
static int *S;				/* saving operations for diff3 */
static int *ST;				/* saving substitution types */

/* The following definitions are for function diff3() */

int  diff3();
static int  zero = 0;				/* int type zero        */

#define gap(k)  ((k) <= 0 ? 0 : q+r*(k))	/* k-symbol indel score */

#define gap2(k)  ((k) <= 0 ? 0 : ((k) <= gaplen ? q+r*(k) : pay))
/* k-symbol insertion score */

static int *sapp;				/* Current script append ptr */
static int *stptr;				/* Substitution type ptr */
static int  last;				/* Last script op appended */

static int no_mat; 				/* number of matches */ 
static int no_mis; 				/* number of mismatches */ 
static int no_gap; 				/* number of indels */
static int term_gap5; 			/* no. of amono acids in 5' terminal gap */
static int term_gap3; 			/* no. of amono acids in 3' terminal gap */
/* Append "Delete k" op */
#define DEL(k)				\
{ no_gap += (3 * (k));			\
  if (last < 0)				\
    last = sapp[-1] -= (k);		\
  else					\
    last = *sapp++ = -(k);		\
  term_gap3 += (k);                     \
  if ( ! no_mat && ! no_mis )           \
    term_gap5 += (k);                   \
}
/* Append "Insert k" op */
#define INS(k)				\
{ no_gap += (k);			\
  if (last > 0)				\
    last = sapp[-1] += (k);		\
  else					\
    last = *sapp++ = (k);		\
}
/* Append "Replace" op */
#define REP(t) 				\
{ last = *sapp++ = 0; 			\
  *stptr++ = (t);                       \
  term_gap3 = 0;	                \
}

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

main(argc, argv) int argc; char *argv[];
{ int  M;				/* Sequence length */
 char *D, *D0,  *P;			/* Storing sequences */
 int  *PN, *DN, *DN0;			/* sequences of integer codes */
 int  symbol;				/* The next character	      */
 int  ms;				/* User-supplied weights      */
 FILE *Bp, *Ap, *Sp, *ckopen();
 char *ckalloc();			/* space-allocating function  */
 register int i, j;
 char  alph[129], *s;			/* alphabet */
 int  size;				/* size of alphabet */
 char  cline[ACCLEN+1];		/* a line of chain file */
 int  aalen;				/* maximu database sequence length */
 char accn[ACCLEN+1];			/* current protein accession */
 int  dstart, dend;			/* DNA positions */
 int  astart, aend;			/* protein positions */
 int  score;				/* alignment score */
 int  ort;				/* orientation */
 float  denst;				/* DNA length over protein length */
 int  off;				/* adjustment to length of DNA region */
 int  temp;				/* temp variable */
 typedef struct MSC   {
   int   dstart;            /* start position DNA sequence */
   int   dend;              /* end position DNA sequence */
   int   astart;            /* start position in protein sequence */
   int   aend;              /* end position in protein sequence */
   int   score;             /* score of chain */
   int   ort;               /* orientation */
   int   soff;              /* offset for dstart */
   int   eoff;              /* offset for dend */
   char  acc[ACCLEN];	/* protein accession */
   char  *aa;		/* pointer to protein seq */
   int   alen;		/* length of protein seq */
 }  chaintp, *chainptr;  /* type definition */
 chainptr chains;                /* for holding chains */
 chainptr  cpt;			/* temp pointer variable */
 int    chlen;			/* number of chains */
 char   *apt;			/* temp pointer var */
 char   *pp;			/* pointer to current protein seq */
 long  cur_time;
 
 program = argv[0];
 gaplen = QGAPLEN;
 q = GAPOPEN;
 r = GAPEXTN;
 
 
 if ( argc < 2 )
   { 
     fprintf(stderr,"Usage: %s DNA_Seq Prot_Database Chain_File BLOSUM [options]\n\n", argv[0]);
     fprintf(stderr,"  DNA_Seq        file of one query DNA in FASTA format\n");
     fprintf(stderr,"  Prot_Database  file of protein database in FASTA format\n");
     fprintf(stderr,"  Chain_File     file of coordinates produced by EXT or Filter\n");
     fprintf(stderr,"  BLOSUM         file of specially formatted BLOSUM matrix\n");
     fprintf(stderr,"Options (default values):\n");
     fprintf(stderr,"  -e  N  specify gap extension penalty N > 0 (1)\n");
     fprintf(stderr,"  -o  N  specify gap open penalty N >= 0 (15)\n");
     fprintf(stderr,"  -q  N  specify gap length N >= 5 for constant penalty (25)\n");
     fprintf(stderr,"  -x  N  specify min exon length in amino acids N > 0 and < %d (7)\n", CSLEN);
     exit(1);
   }
 
 if (!strcmp(argv[1], "-help")  ||  !strcmp(argv[1], "-h"))
   {
     fprintf(stderr,"Usage: %s DNA_Seq Prot_Database Chain_File BLOSUM [options]\n\n", argv[0]);
     fprintf(stderr,"  DNA_Seq        file of one query DNA in FASTA format\n");
     fprintf(stderr,"  Prot_Database  file of protein database in FASTA format\n");
     fprintf(stderr,"  Chain_File     file of coordinates produced by EXT or Filter\n");
     fprintf(stderr,"  BLOSUM         file of specially formatted BLOSUM matrix\n");
     fprintf(stderr,"Options (default values):\n");
     fprintf(stderr,"  -e  N  specify gap extension penalty N > 0 (1)\n");
     fprintf(stderr,"  -o  N  specify gap open penalty N >= 0 (15)\n");
     fprintf(stderr,"  -q  N  specify gap length N >= 5 for constant penalty (25)\n");
     fprintf(stderr,"  -x  N  specify min exon length in amino acids N > 0 and < %d (7)\n", CSLEN);
     exit(1);
   }
 if (!strcmp(argv[1], "-version")  ||  !strcmp(argv[1], "-V"))
   {
     version();
     exit(1);
   }
 if (!strcmp(argv[1], "-depend"))
   {
     // No dependencies to print on stderr
     exit(1);
   }
 if (!strcmp(argv[1], "-debug"))
   {
     fprintf(stderr, "default gap length = %d\n", gaplen); 
     fprintf(stderr, "default gap open penalty  = %d\n", q);
     fprintf(stderr, "default gap extension penalty = %d\n", r); 
     exit(1);
   }   
 
 if ( argc < 5 )
   { fprintf(stderr,"Usage: %s DNA_Seq Prot_Database Chain_File BLOSUM [options]\n\n", argv[0]);
   fprintf(stderr,"  DNA_Seq        file of one query DNA in FASTA format\n");
   fprintf(stderr,"  Prot_Database  file of protein database in FASTA format\n");
   fprintf(stderr,"  Chain_File     file of coordinates produced by EXT or Filter\n");
   fprintf(stderr,"  BLOSUM         file of specially formatted BLOSUM matrix\n");
   fprintf(stderr,"Options (default values):\n");
   fprintf(stderr,"  -e  N  specify gap extension penalty N > 0 (1)\n");
   fprintf(stderr,"  -o  N  specify gap open penalty N >= 0 (15)\n");
   fprintf(stderr,"  -q  N  specify gap length N >= 5 for constant penalty (25)\n");
   fprintf(stderr,"  -x  N  specify min exon length in amino acids N > 0 and < %d (7)\n", CSLEN);
   exit(1);
   }
 time(&cur_time);
 // strftime(btab_date, sizeof(btab_date), "%b %e %Y", localtime(&cur_time));
 
 ComputeAiDi();
 
 if ( argc > 5 )
   { for ( i = 5; i < argc - 1; i += 2 )
     { if ( argv[i][0] != '-' )
       fatal("Each option must begin with a dash (-)");
     (void) sscanf(argv[i+1],"%d", &temp);
     switch ( argv[i][1] )
       { case 'e' : if ( temp <= 0 )
	 fatal("Number for gap extension must be a positive integer");
	   r = temp;
	   break;
       case 'o' : if ( temp < 0 )
	 fatal("Number for gap open must be a non-negative integer");
	 q = temp;
	 break;
       case 'q' : if ( temp < 5 )
	 fatal("Number for gap length must be a positive integer >= 5");
	 gaplen = temp;
	 break;
       case 'x' : if ( temp >= CSLEN || temp < 1 )
	 fatalf("The value for min exon length must be > 0 and < %d", CSLEN);
	 minexon = temp;
	 break;
       default  : fatalf("Wrong option letter "); 
	 break;
       }
     }
   }
 
 /* read the chain file */
 Ap = ckopen(argv[3], "r");
 for ( chlen = 0; fgets(cline, ACCLEN, Ap) != NULL; )
   chlen++;
 if ( chlen < 2 ) {
   fprintf(stderr, "No chain was reported by DPS");
   exit(0);
 }
 chlen--;
 chains = ( chainptr ) ckalloc(chlen * sizeof(chaintp));
 (void) fclose(Ap);
 Ap = ckopen(argv[3], "r");
 fgets(cline, ACCLEN, Ap);
 sscanf(cline, "                     %8d %8d\n", &dnalen, &aalen);
 for ( cpt = chains; fgets(cline, ACCLEN, Ap) != NULL; cpt->aa = NULL, cpt++)
   sscanf(cline, "%8d %8d %6d %7d %5d %1d %5d %5d %s", &cpt->dstart,
	  &cpt->dend, &cpt->score, &cpt->astart, &cpt->aend, &cpt->ort,
	  &cpt->soff, &cpt->eoff, cpt->acc);
 (void) fclose(Ap);
 
 
 /* determine the heading length */
 Ap = ckopen(argv[1], "r");
 if ( (symbol = getc(Ap) ) != '>' )
   fatal("The DNA sequence must be in the FASTA format, starting with '>'");
 for (size = 2; ( symbol = getc(Ap)) != EOF ; size++ )
   if ( symbol == '\n' )
     break;
 (void) fclose(Ap);
 name1 = argv[1];
 
 /* allocate space for A */
 D = ( char * ) ckalloc( (dnalen + 1) * sizeof(char));
 DN = ( int * ) ckalloc( (dnalen + 1) * sizeof(int));
 dhead = ( char * ) ckalloc( size * sizeof(char));
 
 /* read the DNA sequence into D */
 Ap = ckopen(argv[1], "r");
 for (i = 0; ( symbol = getc(Ap)) != EOF && symbol != '\n' ; )
   dhead[i++] = symbol != '\t' ? symbol : ' ';
 dhead[i] = '\0';      
 for (M = 0; ( symbol = getc(Ap)) != EOF ; )
   if ( symbol != '\n' )
     {	D[++M] = symbol;
     DN[M] = di[symbol];
     }
 if ( M != dnalen )
   fatal("The DNA sequence for NAP must be identical to that for DPS");
 /* allocate space for D0 */
 D0 = ( char * ) ckalloc( (dnalen + 1) * sizeof(char));
 DN0 = ( int * ) ckalloc( (dnalen + 1) * sizeof(int));
 /* Set D0 and DN0 to the reverse complement of D and DN */
 for ( i = 1, j = M; i <= M; i++, j-- )
   switch( DN[i] )
     { case 0 : DN0[j] = 3;
	 D0[j] = 'T';
	 break;
     case 1 : DN0[j] = 2;
       D0[j] = 'G';
       break;
     case 2 : DN0[j] = 1;
       D0[j] = 'C';
       break;
     case 3 : DN0[j] = 0;
       D0[j] = 'A';
       break;
     case 4 : DN0[j] = 4;
       D0[j] = 'N';
       break;
     default : break;
     }
 
 pay = q + r * gaplen;
 qr = q + r;
 r2 = 2 * r;
 r3 = 3 * r;
 qr2 = q + r2;
 qr3 = q + r3;
 q2r2 = 2 * q + r2;
 
 /* check if the argument represents a negative integer */
 s = argv[4];
 if ( *s == '-' ) s++;
 for ( ; *s >= '0' && *s <= '9' ; s++ );
 if ( *s == '\0' )
   { (void) sscanf(argv[4],"%d", &ms);
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
     Sp = ckopen(argv[4], "r");
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
 ComputeDa();
 SetUpW();
 SetUpMt();
 ComputeAa();
 /* allocate space for all vectors */
 j = (dnalen + 1) * sizeof(int);
 CC = ( int * ) ckalloc(j);
 PC = ( int * ) ckalloc(j);
 DD = ( int * ) ckalloc(j);
 RR = ( int * ) ckalloc(j);
 SS = ( int * ) ckalloc(j);
 i = (aalen + 1) * sizeof(int);
 S = ( int * ) ckalloc(i + j);
 ST = ( int * ) ckalloc(i + j);
 
 /* allocate space for P */
 P = ( char * ) ckalloc( (aalen + 1) * sizeof(char));
 PN = ( int * ) ckalloc( (aalen + 1) * sizeof(int));
 Bp = ckopen(name2 = argv[2], "r");
 symbol = getc(Bp);
 (void) printf("Max Match   Min Mismatch   Gap-Open Penalty   Gap-Extension Penalty\n");
 (void) printf("   %d          %d              %d                  %d\n\n", match, mismh, q, r);
 (void) printf("                 Query Sequence : %s\n", name1);
 (void) printf("                 Database: %s\n", name2);
 /* process each protein sequence */
 while ( symbol != EOF )
   { if ( symbol == '>' )
     { for ( i = 0; ( symbol = getc(Bp)) != EOF && symbol != '\n'; )
       if ( i < ACCLEN && ( i != 0 || (symbol != ' ' && symbol != '\t') ) )
	 accn[i++] = symbol;
     if ( i >= ACCLEN )
       i = ACCLEN - 1;
     accn[i] = '\0';
     if ( symbol == EOF )
       fatal("The protein database is not in the correct format");
     }
   else
     fatal("The protein database is not in the correct format");
   for (psize = 0; ( symbol = getc(Bp)) != EOF && symbol != '>'; )
     if ( symbol != '\n' )
       P[++psize] = symbol;
   pp = NULL;
   for ( i = 0; i < chlen; i++ )
     { apt = chains[i].acc;
     for ( j = 0; apt[j] != '\0' && apt[j] == accn[j]; j++ )
       ;
     if ( apt[j] != '\0' )
       continue;
     if ( pp == NULL )
       { pp = ( char * ) ckalloc( (psize + 1) * sizeof(char));
       for ( j = 1; j <= psize; j++ )
	 pp[j] = P[j];
       }
     chains[i].aa = pp;
     chains[i].alen = psize;
     strcpy(chains[i].acc, accn);
     }
   }
 (void) fclose(Bp);
 sprintf(fname, "%s.nap.btab", argv[1]);
 Btabp = ckopen(fname, "w");
 for ( i = 0; i < chlen; i++ )
   { if ( (pp = chains[i].aa) == NULL )
     fatal("No protein seq was found for one chain.");
   psize = chains[i].alen;
   strcpy(dblocus, chains[i].acc);
   for ( j = 0; dblocus[j] != '\0'; j++ )
     if ( dblocus[j] == '\t' )
       dblocus[j] = ' ';
   for ( j = 0; dblocus[j] != '\0' && dblocus[j] != ' '; j++ )
     ;
   if ( dblocus[j] == '\0' )
     strcpy(dbdesp, " \0");
   else
     { strcpy(dbdesp, dblocus + j + 1);
     dblocus[j] = '\0';
     }
   for ( j = 1; j <= psize; j++ )
     { P[j] = symbol = pp[j];
     PN[j] = ai[symbol];
     }
   denst = (float) (chains[i].dend - chains[i].dstart + 1) /
     (float) (chains[i].aend - chains[i].astart + 1);
   if ( chains[i].ort )
     { dstart = chains[i].dstart;
     dend = chains[i].dend;
     }
   else
     { dstart = dnalen - chains[i].dend + 1;
     dend = dnalen - chains[i].dstart + 1;
     }
   if ( (temp = chains[i].astart) > 1 )
     { if ( (off = denst * temp + EXTRA) < chains[i].soff )
       off = chains[i].soff;
     dstart = (temp = dstart - off) > 0 ? temp : 1;
     }
   else
     { if ( dstart > EXTRA )
       dstart -= EXTRA;
     }
   if ( (temp = psize - chains[i].aend) > 0 )
     { if ( (off = denst * temp + EXTRA) < chains[i].eoff )
       off = chains[i].eoff;
     dend = (temp = dend + off) <= dnalen? temp : dnalen;
     }
   else
     { if ( dend + EXTRA <= dnalen )
       dend += EXTRA;
     }
   sapp = S;
   stptr = ST;
   last = 0;
   no_gap = 0;
   no_mat = 0;
   no_mis = 0;
   term_gap5 = 0;
   term_gap3 = 0;
   temp = dend-dstart+1;
   if ( chains[i].ort )
     { score = diff3(PN,DN+dstart-1,psize,temp,q,q,0,0,0,0);
     printf("\nAccession: %s\n", chains[i].acc);
     if ( term_gap5 + term_gap3 < psize )
       size = psize - term_gap5 - term_gap3;
     else
       size = psize;
     printf("Score: %d  Identity: %d/%d (%d%%)  Strand: %s\n",
	    score, no_mat, size, (100*no_mat)/size, "plus");
     display(P,D+dstart-1,psize,temp,PN,DN+dstart-1,S,ST,1,dstart,1,aa2, i);//add chain number i
     }
   else
     { score = diff3(PN,DN0+dstart-1,psize,temp,q,q,0,0,0,0);
     printf("\nAccession: %s\n", chains[i].acc);
     if ( term_gap5 + term_gap3 < psize )
       size = psize - term_gap5 - term_gap3;
     else
       size = psize;
     printf("Score: %d  Identity: %d/%d (%d%%)  Strand: %s\n",
	    score, no_mat, size, (100*no_mat)/size, "minus");
     display(P,D0+dstart-1,psize,temp,PN,DN0+dstart-1,S,ST,1,dstart,0,aa2, i);// add chain number
     }
   }
 
 exit(0);
 
}

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
 bonus5 = BONUS5;
 bonus3 = BONUS3;
}

/* Set up w, bbn, bnb, nbb, bnn, nbn, nnb */
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
 int  x;	/* index for bbn, bnb, nbb */
 int  NXX, XNX, NNX, NXN, XNN;	/* constants */	
 
 
 ds1 = DNASIZE - 1;
 NXX = ds1 * dnals2;
 XNX = ds1 * DNASIZE;
 NNX = NXX + XNX;
 NXN = NXX + ds1; 
 XNN = XNX + ds1;
 /* Set up bbn, bnb, nbb, bnn, nbn, nnb */
 for ( i = 0; i < DNASIZE; i++ )
   { bnn[i] = i * dnals2 + XNN; 
   b = i * DNASIZE;
   for ( j = 0; j < DNASIZE; j++ )
     { bbb = (x = b + j) * DNASIZE + ds1;
     bbn[x] = bbb; 
     }
   }
 for ( k = 0; k < DNASIZE; k++ )
   { nnb[k] = NNX + k; 
   for ( j = 0; j < DNASIZE; j++ )
     { x = j * DNASIZE + k; 
     nbb[x] = NXX + x;
     }
   }
 for ( j = 0; j < DNASIZE; j++ )
   nbn[j] = NXN + j * DNASIZE; 
 for ( i = 0; i < DNASIZE; i++ )
   { b = i * dnals2;
   x = i * DNASIZE;
   for ( k = 0; k < DNASIZE; k++ )
     bnb[x + k] = b + XNX + k;
   }
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
   wa[NXX + XNN] = 0;
   /* compute the substitution scores for stop codons */
   wa[ochre] = mismh;
   wa[amber] = mismh;
   wa[uga] = mismh;
   }
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
{  //int strcpy();
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
  
  strcpy(aa2[0], "A  ");
  strcpy(aa2[1], "R  ");
  strcpy(aa2[2], "N  ");
  strcpy(aa2[3], "D  ");
  strcpy(aa2[4], "C  ");
  strcpy(aa2[5], "Q  ");
  strcpy(aa2[6], "E  ");
  strcpy(aa2[7], "G  ");
  strcpy(aa2[8], "H  ");
  strcpy(aa2[9], "I  ");
  strcpy(aa2[10], "L  ");
  strcpy(aa2[11], "K  ");
  strcpy(aa2[12], "M  ");
  strcpy(aa2[13], "F  ");
  strcpy(aa2[14], "P  ");
  strcpy(aa2[15], "S  ");
  strcpy(aa2[16], "T  ");
  strcpy(aa2[17], "W  ");
  strcpy(aa2[18], "Y  ");
  strcpy(aa2[19], "V  ");
  strcpy(aa2[20], "B  ");
  strcpy(aa2[21], "Z  ");
  strcpy(aa2[22], "X  ");
  
  aa3[0] = 'A';
  aa3[1] = 'R';
  aa3[2] = 'N';
  aa3[3] = 'D';
  aa3[4] = 'C';
  aa3[5] = 'Q';
  aa3[6] = 'E';
  aa3[7] = 'G';
  aa3[8] = 'H';
  aa3[9] = 'I';
  aa3[10] = 'L';
  aa3[11] = 'K';
  aa3[12] = 'M';
  aa3[13] = 'F';
  aa3[14] = 'P';
  aa3[15] = 'S';
  aa3[16] = 'T';
  aa3[17] = 'W';
  aa3[18] = 'Y';
  aa3[19] = 'V';
  aa3[20] = 'B';
  aa3[21] = 'Z';
  aa3[22] = 'X';
  aa3[23] = '*';
}

/* diff3(AN,DN,M,N,tb,te,sc,sr,ec,er) returns the score of an optimum conversion
   between AN[1..M] and DN[1..N] that begins(ends) with a delete if tb(te) is zero
   and appends such a conversion to the current script. If sc = 0, then
   the beginning deletion is not penalized; if sr = 0, the beginning insertion is
   not penalized; if ec = 0, the ending deletion is not charged; if er = 0;
   then the ending insertion is not charged. Any insertion of length at least
   gaplen is given a constant cost */

int diff3(AN,DN,M,N,tb,te,sc,sr,ec,er)
     int *AN, *DN; int M, N; int tb, te, sc, sr, ec, er;
     
{ int  midi, midj, type;	/* Midpoint, type, and cost */
 int  midc;			/* C(midi, midj) */
 int  ss,cc;			/* terminal gap penalty indicators */
 
 { register int  i, j;		/* index variables */
 register int  c;		/* C(i, j) */
 register int  e;		/* E(i, j) */
 register int  d;		/* D(i, j) */
 register int  g;		/* G(i, j) */
 register int  s1;		/* C(i-1, j-1) */
 register int  s2;		/* C(i-1, j-2) */
 register int  s3;		/* C(i-1, j-3) */
 register int  s4;		/* C(i-1, j-4) */
 register int  f1;		/* D(i-1, j-1) */
 register int  f2;		/* D(i-1, j-2) */
 int  t;		/* D(i, 0) */
 int  *wa;		/* one row of w */
 int  temp;		/* temp variable */
 char *ckalloc();
 int  tp;		/* substitution type */
 int  codl;		/* length of the DNA codon */
 int  pen;		/* gap penalty */
 int  b;		/* one-symbol integer code */
 int  bb;		/* two-symbol integer code */
 int  bbb;		/* three-symbol integer code */
 int  bj;		/* bbb at midj */
 int  bdbb;		/* integer code of bbb */
 int  bbdb;		/* integer code of bbb */
 int  x, y, z;	/* temp variables */
 int  tbsc;		/* tb * sc */
 int  teec;		/* te * ec */
 int  ds;		/* D(i, j) */
 int  dp;		/* starting position for ds */
 int  hs1;		/* H(i, j) */
 int  hp1;		/* starting position for hs1 */
 int  hs2;		/* H'(i, j) */
 int  hp2;		/* starting position for hs2 */
 
 /* Boundary cases: M <= 1 or N == 0 */
 
 if (N <= 0)
   { if (M > 0) DEL(M)
		  if ( !sc || !ec )
		    return 0;
		  else
		    return - gap(3 * M);
   }
 if (M <= 1)
   { if (M <= 0)
     { INS(N);
     if ( !sr || !er )
       return 0;
     else
       return - gap2(N);
     }
   midc = - (sc * (tb + r3) + er * (z = gap2(N) ) );
   midj = -1;
   if ( midc < ( c =  - (ec * (te + r3) + sr * z ) ) )
     { midc = c;
     midj = 0;
     }
   f1 = ( N > 2 && gtcode == DN[1] * DNASIZE + DN[2] ) ? 1 : 0;
   f2 = ( N > 2 && agcode == DN[N-1] * DNASIZE + DN[N] ) ? 1 : 0;
   tbsc = tb * sc;
   teec = te * ec;
   wa = w[AN[1]];
   ds = hs1 = hs2 = MINF;
   dp = hp1 = hp2 = -1;
   for (j = 1; j <= N; j++)
     { b = DN[j];
     y = N - j;
     temp = er * gap2(y);
     if ( er && y > gaplen )
       { if ( gtcode == DN[j+1] * DNASIZE + DN[j+2] )
	 temp -= bonus5;
       if ( f2 )
	 temp -= bonus3;
       }
     y = j - 1;
     pen = sr * gap2(y) + temp + r2;
     if ( sr && y > gaplen )
       { if ( f1 )
	 pen -= bonus5;
       if ( agcode == DN[y-1] * DNASIZE + DN[y] )
	 pen -= bonus3;
       }
     z = j < N ? q : teec;
     x = y ? q : tbsc;
     if ( (c = wa[nnb[b]] - pen - x) > midc )
       { midc = c;
       midj = j;
       codl = 1;
       tp = 1;
       }
     if ( (c = wa[nbn[b]] - pen - x - z) > midc )
       { midc = c;
       midj = j;
       codl = 1;
       tp = 2;
       }
     if ( (c = wa[bnn[b]] - pen - z) > midc )
       { midc = c;
       midj = j;
       codl = 1;
       tp = 3;
       }
     if ( j > 1 )
       { bb = DN[j-1] * DNASIZE + b;
       y = j - 2;
       pen = sr * gap2(y) + temp + r;
       if ( sr && y > gaplen )
	 { if ( f1 )
	   pen -= bonus5;
	 if ( agcode == DN[y-1] * DNASIZE + DN[y] )
	   pen -= bonus3;
	 }
       x = y ? q : tbsc;
       if ( (c = wa[nbb[bb]] - pen - x) > midc )
	 { midc = c;
	 midj = j;
	 codl = 2;
	 tp = 4;
	 }
       if ( (c = wa[bnb[bb]] - pen - q) > midc )
	 { midc = c;
	 midj = j;
	 codl = 2;
	 tp = 5;
	 }
       if ( (c = wa[bbn[bb]] - pen - z) > midc )
	 { midc = c;
	 midj = j;
	 codl = 2;
	 tp = 6;
	 }
       if ( j > 2 )
	 { bbb = DN[j-2] * dnals2 + bb;
	 y = j - 3;
	 pen = sr * gap2(y) + temp;
	 if ( sr && y > gaplen )
	   { if ( f1 )
	     pen -= bonus5;
	   if ( agcode == DN[y-1] * DNASIZE + DN[y] )
	     pen -= bonus3;
	   }
	 if ( (c = wa[bbb] - pen) > midc )
	   { midc = c;
	   midj = j;
	   codl = 3;
	   tp = 7;
	   bj = bbb;
	   }
	 if ( j > 3 )
	   { bdbb = (x = DN[j-3] * dnals2) + bb;
	   y = j - 4;
	   s4 = sr * gap2(y);
	   if ( sr && y > gaplen )
	     { if ( f1 )
	       s4 -= bonus5;
	     if ( agcode == DN[y-1] * DNASIZE + DN[y] )
	       s4 -= bonus3;
	     }
	   pen = s4 + temp + qr;
	   if ( (c = wa[bdbb] - pen) > midc )
	     { midc = c;
	     midj = j;
	     codl = 4;
	     tp = 8;
	     }
	   bbdb = x + DN[j-2] * DNASIZE + b;
	   if ( (c = wa[bbdb] - pen) > midc )
	     { midc = c;
	     midj = j;
	     codl = 4;
	     tp = 9;
	     }
	   if ( (x = - (s4 + qr) ) > (ds = ds - r) )
	     { ds = x;
	     dp = j - 3;
	     }
	   bbdb = (y = DN[dp] * dnals2) + DN[dp+1] * DNASIZE + b;
	   if ( (c = ds + wa[bbdb] - temp) > midc )
	     { midc = c;
	     midj = j;
	     codl = j - dp + 1;
	     tp = 10;
	     }
	   bdbb = y + bb;
	   if ( (c = ds + wa[bdbb] - temp) > midc )
	     { midc = c;
	     midj = j;
	     codl = j - dp + 1;
	     tp = 12;
	     }
	   if ( (z = j - gaplen - 4) >= 0 )
	     { s1 = sr * gap2(z) + pay;
	     if ( sr && z > gaplen )
	       { if ( f1 )
		 s1 -= bonus5;
	       if ( agcode == DN[z-1] * DNASIZE + DN[z] )
		 s1 -= bonus3;
	       }
	     x = gtcode == (DN[z+3] * DNASIZE + DN[z+4]) ? bonus5 : 0;
	     if ( (y = x - s1 ) > hs2 )
	       { hs2 = y;
	       hp2 = z + 1;
	       }
	     bbdb = DN[hp2] * dnals2 + DN[hp2+1] * DNASIZE + b;
	     x = agcode == (DN[j-2] * DNASIZE + DN[j-1]) ? bonus3 : 0;
	     if ( (c = hs2 + x + wa[bbdb] - temp) > midc )
	       { midc = c;
	       midj = j;
	       codl = j - hp2 + 1;
	       tp = 10;
	       }
	     x = gtcode == (DN[z+2] * DNASIZE + DN[z+3]) ? bonus5 : 0;
	     if ( (y = x - s1 ) > hs1 )
	       { hs1 = y;
	       hp1 = z + 1;
	       }
	     bdbb = DN[hp1] * dnals2 + bb;
	     x = agcode == (DN[j-3] * DNASIZE + DN[j-2]) ? bonus3 : 0;
	     if ( (c = hs1 + x + wa[bdbb] - temp) > midc )
	       { midc = c;
	       midj = j;
	       codl = j - hp1 + 1;
	       tp = 12;
	       }
	     }
	   }
	 }
       }
     }
   if (midj == -1)
     { DEL(1) INS(N) }
   else
     if (midj == 0)
       { INS(N) DEL(1) }
     else
       { if (midj > codl) INS(midj-codl)
			    REP(tp)
			    if ( tp > 9 )
			      { INS(codl-3)
				  REP(tp+1)
				  }
			    else
			      no_gap += ( z = codl-3 ) >= 0 ? z : (-z);
       if ( tp == 7 && AN[1] == da[bj] )
	 no_mat += 1;
       else
	 no_mis += 1;
       if (midj < N) INS(N-midj)
		       }
   return midc;
   }
 
 /* Divide: Find optimum midpoint (midi,midj) of cost midc */
 
 midi = M/2;			/* Forward phase:                          */
 CC[0] = 0;			/*   Compute C(M/2,k) & D(M/2,k) for all k */
 t = - q * sr;
 if ( N <= gaplen )
   for (j = 1; j <= N; j++)
     { CC[j] = t = (t-r) * sr;
     DD[j] = t-q;
     }
 else
   { for (j = 1; j <= gaplen; j++)
     { CC[j] = t = (t-r) * sr;
     DD[j] = t-q;
     }
   x = gtcode == (DN[1] * DNASIZE + DN[2]) ? bonus5 : 0;
   x = (x - pay) * sr;
   for (j = gaplen+1; j <= N; j++)
     { y = (sr && agcode == (DN[j-1] * DNASIZE + DN[j])) ? bonus3 : 0;
     CC[j] = t = x + y;
     DD[j] = t - q;
     }
   }
 if ( !ec ) DD[N] += q;
 t = -tb * sc;
 s2 = s3 = s4 = f2 = 0;
 for (i = 1; i <= midi; i++)
   { s1 = PC[0] = CC[0];
   f1 = t;
   CC[0] = c = t = (t-r3) * sc;
   e = t-q;
   g = t - pay;
   wa = w[AN[i]];
   ds = hs1 = hs2 = t+MINF;
   dp = hp1 = hp2 = -1;
   for (j = 1; j <= N; j++)
     { b = DN[j];
     if ( j > 1 ) bb = DN[j-1] * DNASIZE + b;
     if ((c = c - qr) > (e = e - r)) e = c;
     if ( j == N && !ec )
       { if ((c = CC[j] ) > (d = DD[j] )) d = c;}
     else
       { if ((c = CC[j] - qr3) > (d = DD[j] - r3)) d = c;
       if ((c = f1 + (y = wa[nbn[b]]) - qr2) > d ) d = c;
       if ((x = wa[bnn[b]] - qr2) > (y = y - q2r2 ) ) y = x;
       if ((c = s1 + y) > d ) d = c;
       if ( j > 1 && (c = s2 + wa[bbn[bb]] - qr) > d ) d = c;
       }
     x = wa[nnb[b]] - r2;
     if ( (y = f1 + x) > (c = s1 + x - q) ) c = y;
     if ( j > 1 )
       { if ((x = wa[nbb[bb]]) > (y = wa[bnb[bb]]) ) y = x;
       if ( (y = s2 + y - qr) > c ) c = y;
       if ( (y = f2 + x - r) > c ) c = y;
       if ( j > 2 )
	 { bbb = DN[j-2] * dnals2 + bb;
	 if ( (x = s3 + wa[bbb]) > c ) c = x;
	 if ( j > 3 )
	   { bdbb = (y = DN[j-3] * dnals2) + bb;
	   bbdb = y + DN[j-2] * DNASIZE + b;
	   if ( (x = wa[bdbb]) > (y = wa[bbdb]) ) y = x;
	   if ( (x = s4 + y - qr) > c ) c = x;
	   if ( ( x = s4 - qr ) > ( ds = ds - r ) )
	     { ds = x;
	     dp = j - 3;
	     }
	   bdbb = (y = DN[dp] * dnals2) + bb;
	   bbdb = y + DN[dp+1] * DNASIZE + b;
	   if ( (x = wa[bdbb]) > (y = wa[bbdb]) ) y = x;
	   if ( (x = ds + y) > c ) c = x;
	   }
	 }
       }
     if (c < d) c = d;
     if (c < e) c = e;
     if ( (z = j - gaplen - 1) >= 0 )
       { x = gtcode == (DN[z+1] * DNASIZE + DN[z+2]) ? bonus5 : 0;
       if ( g < ( temp = CC[z] + x - pay ) )
	 g = temp;
       x = agcode == bb ? bonus3 : 0;
       if ( c < (y = g + x) ) c = y;
       if ( (z = z - 3) >= 0 )
	 { x = gtcode == (DN[z+2] * DNASIZE + DN[z+3]) ? bonus5 : 0;
	 if ( (y = x + (temp = PC[z] - pay) ) > hs1 )
	   { hs1 = y;
	   hp1 = z + 1;
	   }
	 bdbb = DN[hp1] * dnals2 + bb;
	 x = agcode == (DN[j-3] * DNASIZE + DN[j-2]) ? bonus3 : 0;
	 if ( (y = hs1 + x + wa[bdbb]) > c ) c = y;
	 x = gtcode == (DN[z+3] * DNASIZE + DN[z+4]) ? bonus5 : 0;
	 if ( (y = x + temp) > hs2 )
	   { hs2 = y;
	   hp2 = z + 1;
	   }
	 bbdb = DN[hp2] * dnals2 + DN[hp2+1] * DNASIZE + b;
	 x = agcode == (DN[j-2] * DNASIZE + DN[j-1]) ? bonus3 : 0;
	 if ( (y = hs2 + x + wa[bbdb]) > c ) c = y;
	 }
       }
     s4 = s3;
     s3 = s2;
     s2 = s1;
     s1 = PC[j] = CC[j];
     CC[j] = c;
     f2 = f1;
     f1 = DD[j];
     DD[j] = d;
     }
   }
 DD[0] = CC[0];
 
 RR[N] = 0;			/* Reverse phase:                          */
 t = -q * er;			/*   Compute R(M/2,k) & S(M/2,k) for all k */
 if ( N <= gaplen )
   for (j = N-1; j >= 0; j--)
     { RR[j] = t = (t-r) * er;
     SS[j] = t-q;
     }
 else
   { temp = N - gaplen;
   for (j = N-1; j >= temp; j--)
     { RR[j] = t = (t-r) * er;
     SS[j] = t-q;
     }
   x = agcode == (DN[N-1] * DNASIZE + DN[N]) ? bonus3 : 0;
   x = (x - pay) * er;
   for (j = temp-1; j >= 0; j--)
     { y = (er && gtcode == (DN[j+1] * DNASIZE + DN[j+2]) ) ? bonus5 : 0;
     RR[j] = t = x + y;
     SS[j] = t - q;
     }
   }
 if ( !sc ) SS[0] += q;
 t = -te * ec;
 for (i = M-1; i >= midi; i--)
   { s1 = PC[N] = RR[N];
   f1 = t;
   RR[N] = c = t = (t-r3) * ec;
   g = t - pay;
   e = t-q;
   wa = w[AN[i+1]];
   ds = hs1 = hs2 = t+MINF;
   dp = hp1 = hp2 = -1;
   for (j = N-1; j >= 0; j--)
     { b = DN[j+1];
     if ( j < N-1 ) bb = b * DNASIZE + DN[j+2];
     if ((c = c - qr) > (e = e - r)) e = c;
     if ( !j && !sc )
       { if ((c = RR[j] ) > (d = SS[j] )) d = c;}
     else
       { if ((c = RR[j] - qr3) > (d = SS[j] - r3)) d = c;
       if ((c = f1 + (y = wa[nbn[b]]) - qr2) > d ) d = c;
       if ((x = wa[nnb[b]] - qr2) > (y = y - q2r2 ) ) y = x;
       if ((c = s1 + y) > d ) d = c;
       if ( j < (N-1) && (c = s2 + wa[nbb[bb]] - qr) > d ) d = c;
       }
     x = wa[bnn[b]] - r2;
     if ( (y = f1 + x) > (c = s1 + x - q) ) c = y;
     if ( j < N - 1 )
       { if ((x = wa[bbn[bb]]) > (y = wa[bnb[bb]]) ) y = x;
       if ( (y = s2 + y - qr) > c ) c = y;
       if ( (y = f2 + x - r) > c ) c = y;
       if ( j < N - 2 )
	 { bbb = bb * DNASIZE + DN[j+3];
	 if ( (x = s3 + wa[bbb]) > c ) c = x;
	 if ( j < N - 3 )
	   { bdbb = (z = b * dnals2) + DN[j+3] * DNASIZE + (bj = DN[j+4]);
	   bbdb = (temp = bb * DNASIZE) + bj;
	   if ( (x = wa[bdbb]) > (y = wa[bbdb]) ) y = x;
	   if ( (x = s4 + y - qr) > c ) c = x;
	   if ( ( x = s4 - qr ) > ( ds = ds - r ) )
	     { ds = x;
	     dp = j + 4;
	     }
	   bdbb = z + DN[dp-1] * DNASIZE + (x = DN[dp]);
	   bbdb = temp + x;
	   if ( (x = wa[bdbb]) > (y = wa[bbdb]) ) y = x;
	   if ( (x = ds + y) > c ) c = x;
	   }
	 }
       }
     if (c < d) c = d;
     if (c < e) c = e;
     if ( (z = j + gaplen + 1) <= N )
       { x = agcode == (DN[z-1] * DNASIZE + DN[z]) ? bonus3 : 0;
       if ( g < ( temp = RR[z] + x - pay ) )
	 g = temp;
       x = gtcode == bb ? bonus5 : 0;
       if ( c < (y = g + x) ) c = y;
       if ( (z = z + 3) <= N )
	 { x = agcode == (DN[z-2] * DNASIZE + DN[z-1]) ? bonus3 : 0;
	 if ( (y = x + (temp = PC[z] - pay) ) > hs1 )
	   { hs1 = y;
	   hp1 = z;
	   }
	 bbdb =  bb * DNASIZE + DN[hp1];
	 x = gtcode == (DN[j+3] * DNASIZE + DN[j+4]) ? bonus5 : 0;
	 if ( (y = hs1 + x + wa[bbdb]) > c ) c = y;
	 x = agcode == (DN[z-3] * DNASIZE + DN[z-2]) ? bonus3 : 0;
	 if ( (y = x + temp) > hs2 )
	   { hs2 = y;
	   hp2 = z;
	   }
	 bdbb = b * dnals2 + DN[hp2-1] * DNASIZE + DN[hp2];
	 x = gtcode == (DN[j+2] * DNASIZE + DN[j+3]) ? bonus5 : 0;
	 if ( (y = hs2 + x + wa[bdbb]) > c ) c = y;
	 }
       }
     s4 = s3;
     s3 = s2;
     s2 = s1;
     s1 = PC[j] = RR[j];
     RR[j] = c;
     f2 = f1;
     f1 = SS[j];
     SS[j] = d;
     }
   }
 SS[N] = RR[N];
 
 midc = CC[0]+RR[0];		/* Find optimal midpoint */
 midj = 0;
 type = 1;
 for (j = 0; j <= N; j++)
   if ((c = CC[j] + RR[j]) >= midc)
     if (c > midc || CC[j] != DD[j] && RR[j] == SS[j])
       { midc = c;
       midj = j;
       }
 for (j = N; j >= 0; j--)
   { if ( j == N )
     d = q * ec;
   else
     if ( j == 0 )
       d = q * sc;
     else
       d = q;
   if ((c = DD[j] + SS[j] + d) > midc)
     { midc = c;
     midj = j;
     type = 2;
     }
   }
 }
 
 /* Conquer: recursively around midpoint */
 
 cc = midj == N ? ec : 1;
 ss = midj == 0 ? sc : 1;
 if (type == 1)
   { (void) diff3(AN,DN,midi,midj,tb,q,sc,sr,cc,1);
   (void) diff3(AN+midi,DN+midj,M-midi,N-midj,q,te,ss,1,ec,er);
   }
 else
   { (void) diff3(AN,DN,midi,midj,tb,zero,sc,sr,cc,1);
   (void) diff3(AN+midi,DN+midj,M-midi,N-midj,zero,te,ss,1,ec,er);
   }
 return midc;
}
/* Alignment display routine */

static char ALINE[LINELEN+5];	/* protein line */
static char DLINE[LINELEN+5];	/* DNA line */
static char CLINE[LINELEN+5];	/* pair type indicator */ 
static char PLINE[LINELEN+5];	/* predicated amino acids */ 

display(A,B,M,N,AN,DN,S,ST,AP,DP,orient,ttt, chain_num)
     char A[], B[], (*ttt)[4]; int M, N, chain_num, AN[], DN[]; int S[], ST[], AP, DP, orient;
{ register char *a;	/* protein line pointer */ 
 register char *d;	/* DNA line pointer */
 register char *c;	/* center line pointer */
 register char *e;	/* temp pointer */
 register int  i, j;	/* index variables */
 register int  op;	/* operation */
 int  lines;	/* line number */
 int  dp;	/* DNA sequence position */
 int  ap;	/* protein sequence position */
 char *letter;/* pointer to three letter aa string */
 int  bbb;	/* three-symbol integer code */
 char atemp[5]; /* for suffix of a line of length > LINELEN */
 char dtemp[5]; /* .. */
 char ctemp[5]; /* .. */
 char ptemp[5]; /* .. */
 int  ls;	/* length of the suffix */
 int  pst;    /* previous ST value */
 //int strcpy();
 int mdp;	/* the DNA position relative to the plus strand */
 int y, z;	/* current nucleotide integer code */
 int one;	/* current amino acid integer code */
 int nx, nxx, nnx;	/* 'N' part in the code */
 int isexon;		/* 1 if the current region is exon */
 int xstart;		/* start position of an exon */
 int xend;		/* end position of an exon */
 int exonid;		/* exon number */
 int f;		/* trinucleotide code for insertion */
 int g;		/* frame number */
 register char *pp;	/* translated amino acid pointer */
 int flen, tlen;	/* lengths of 5' and 3' ends of an exon */
 int fmc, tmc;	/* number of matches in the ends */
 int tring[CSLEN];	/* match-mismatch vector for 3' end */
 int tpos;		/* current position of tring */
 int gsize;		/* current gap length */
 int gcount;		/* length of a gap removed */
 int flag;		/* 1, print */
 int ncode, xnn;	/* 'N' part in the code */
 int prein;		/* 1, previous pair ends with an insertion */
 
 i = j = op = lines = isexon = exonid = gsize = gcount = 0;
 prein = 1;
 ap = AP;
 dp = DP;
 a = ALINE;
 d = DLINE;
 c = CLINE;
 pp = PLINE;
 pst = 1;
 nx = (ncode = DNASIZE - 1) * DNASIZE;
 nxx = nx * DNASIZE;
 nnx = nxx + nx;
 xnn = nx + ncode;
 chain_num++; // don't start at chain 0 in output.
 while (i < M || j < N)
   { if (op == 0 && *S == 0)
     { op = *S++;
     gsize = 0;
     if ( pst != 10 && pst != 12 )
       letter = ttt[one = AN[++i]];
     switch ( (pst = *ST++ ) )
       { case 1 : *d++ = ' ';
	   *d++ = ' ';
	   *d++ = B[++j];
	   *a++ = *letter++;
	   *a++ = *letter++;
	   *a++ = *letter++;
	   *c++ = '-';
	   *c++ = '-';
	   *c++ = mt3[one][z = nnx + DN[j]];
	   if ( ! isexon )
	     { isexon = 1;
	     xstart = DP + j - 1;
	     flen = tlen = fmc = tmc = tpos = 0;
	     dbstart = AP + i - 1;
	     exlen = idnno = simno = exscore = 0;
	     }
	   *pp++ = 'X';
	   *pp++ = ' ';
	   *pp++ = ' ';
	   if ( isexon )
	     { if ( flen < CSLEN )
	       flen++;
	     if ( tlen < CSLEN )
	       tlen++;
	     else
	       tmc -= tring[tpos];
	     tring[tpos] = 0;
	     tpos = (tpos+1) % CSLEN;
	     exlen++;
	     exscore += (y = w[one][z]) - q * prein - r2;
	     if ( y > 0 )
	       simno++;
	     }
	   prein = 1;
	   break;
       case 2 : *d++ = ' ';
	 *d++ = B[++j];
	 *d++ = ' ';
	 *a++ = *letter++;
	 *a++ = *letter++;
	 *a++ = *letter++;
	 *c++ = '-';
	 *c++ = mt2[one][z = nx + DN[j]];
	 *c++ = '-';
	 if ( ! isexon )
	   { isexon = 1;
	   xstart = DP + j - 1;
	   flen = tlen = fmc = tmc = tpos = 0;
	   dbstart = AP + i - 1;
	   exlen = idnno = simno = exscore = 0;
	   }
	 *pp++ = 'X';
	 *pp++ = ' ';
	 *pp++ = ' ';
	 if ( isexon )
	   { if ( flen < CSLEN )
	     flen++;
	   if ( tlen < CSLEN )
	     tlen++;
	   else
	     tmc -= tring[tpos];
	   tring[tpos] = 0;
	   tpos = (tpos+1) % CSLEN;
	   exlen++;
	   exscore += (y = w[one][z * DNASIZE + ncode]) - q * prein - qr2;
	   if ( y > 0 )
	     simno++;
	   }
	 prein = 0;
	 break;
       case 3 : *d++ = B[++j];
	 *d++ = ' ';
	 *d++ = ' ';
	 *a++ = *letter++;
	 *a++ = *letter++;
	 *a++ = *letter++;
	 *c++ = mt1[one][z = DN[j]];
	 *c++ = '-';
	 *c++ = '-';
	 if ( ! isexon )
	   { isexon = 1;
	   xstart = DP + j - 1;
	   flen = tlen = fmc = tmc = tpos = 0;
	   dbstart = AP + i - 1;
	   exlen = idnno = simno = exscore = 0;
	   }
	 *pp++ = 'X';
	 *pp++ = ' ';
	 *pp++ = ' ';
	 if ( isexon )
	   { if ( flen < CSLEN )
	     flen++;
	   if ( tlen < CSLEN )
	     tlen++;
	   else
	     tmc -= tring[tpos];
	   tring[tpos] = 0;
	   tpos = (tpos+1) % CSLEN;
	   exlen++;
	   exscore += (y = w[one][z * dnals2 + xnn]) - qr2;
	   if ( y > 0 )
	     simno++;
	   }
	 prein = 0;
	 break;
       case 4 : *d++ = ' ';
	 *d++ = B[++j];
	 y = DN[j];
	 *d++ = B[++j];
	 *a++ = *letter++;
	 *a++ = *letter++;
	 *a++ = *letter++;
	 *c++ = '-';
	 *c++ = mt2[one][z = nx + y];
	 *c++ = mt3[one][z = z * DNASIZE + DN[j]];
	 if ( ! isexon )
	   { isexon = 1;
	   xstart = DP + j - 2;
	   flen = tlen = fmc = tmc = tpos = 0;
	   dbstart = AP + i - 1;
	   exlen = idnno = simno = exscore = 0;
	   }
	 *pp++ = 'X';
	 *pp++ = ' ';
	 *pp++ = ' ';
	 if ( isexon )
	   { if ( flen < CSLEN )
	     flen++;
	   if ( tlen < CSLEN )
	     tlen++;
	   else
	     tmc -= tring[tpos];
	   tring[tpos] = 0;
	   tpos = (tpos+1) % CSLEN;
	   exlen++;
	   exscore += (y = w[one][z]) - q * prein - r;
	   if ( y > 0 )
	     simno++;
	   }
	 prein = 1;
	 break;
       case 5 : *d++ = B[++j];
	 y = DN[j];
	 *d++ = ' ';
	 *d++ = B[++j];
	 *a++ = *letter++;
	 *a++ = *letter++;
	 *a++ = *letter++;
	 *c++ = mt1[one][y];
	 *c++ = '-';
	 *c++ = mt3[one][z = y * dnals2 + nx + DN[j]];
	 if ( ! isexon )
	   { isexon = 1;
	   xstart = DP + j - 2;
	   flen = tlen = fmc = tmc = tpos = 0;
	   dbstart = AP + i - 1;
	   exlen = idnno = simno = exscore = 0;
	   }
	 *pp++ = 'X';
	 *pp++ = ' ';
	 *pp++ = ' ';
	 if ( isexon )
	   { if ( flen < CSLEN )
	     flen++;
	   if ( tlen < CSLEN )
	     tlen++;
	   else
	     tmc -= tring[tpos];
	   tring[tpos] = 0;
	   tpos = (tpos+1) % CSLEN;
	   exlen++;
	   exscore += (y = w[one][z]) - qr;
	   if ( y > 0 )
	     simno++;
	   }
	 prein = 1;
	 break;
       case 6 : *d++ = B[++j];
	 y = DN[j];
	 *d++ = B[++j];
	 *d++ = ' ';
	 *a++ = *letter++;
	 *a++ = *letter++;
	 *a++ = *letter++;
	 *c++ = mt1[one][y];
	 *c++ = mt2[one][z = y * DNASIZE + DN[j]];
	 *c++ = '-';
	 if ( ! isexon )
	   { isexon = 1;
	   xstart = DP + j - 2;
	   flen = tlen = fmc = tmc = tpos = 0;
	   dbstart = AP + i - 1;
	   exlen = idnno = simno = exscore = 0;
	   }
	 *pp++ = 'X';
	 *pp++ = ' ';
	 *pp++ = ' ';
	 if ( isexon )
	   { if ( flen < CSLEN )
	     flen++;
	   if ( tlen < CSLEN )
	     tlen++;
	   else
	     tmc -= tring[tpos];
	   tring[tpos] = 0;
	   tpos = (tpos+1) % CSLEN;
	   exlen++;
	   exscore += (y = w[one][z * DNASIZE + ncode]) - qr;
	   if ( y > 0 )
	     simno++;
	   }
	 prein = 0;
	 break;
       case 7 : *d++ = B[++j];
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
	   { *c++ = mt1[one][y];
	   *c++ = mt2[one][z];
	   *c++ = mt3[one][bbb];
	   }
	 if ( ! isexon )
	   { isexon = 1;
	   xstart = DP + j - 3;
	   flen = tlen = fmc = tmc = tpos = 0;
	   dbstart = AP + i - 1;
	   exlen = idnno = simno = exscore = 0;
	   }
	 *pp++ = aa3[da[bbb]];
	 *pp++ = ' ';
	 *pp++ = ' ';
	 if ( isexon )
	   { if ( flen < CSLEN )
	     { flen++;
	     if ( AN[i] == da[bbb] )
	       fmc += 1;
	     }
	   if ( tlen < CSLEN )
	     tlen++;
	   else
	     tmc -= tring[tpos];
	   if ( AN[i] == da[bbb] )
	     { tmc += 1;
	     tring[tpos] = 1;
	     }
	   else
	     tring[tpos] = 0;
	   tpos = (tpos+1) % CSLEN;
	   exlen++;
	   exscore += (y = w[one][bbb]);
	   if ( *(c-1) == ':' )
	     { idnno++;
	     simno++;
	     }
	   else
	     if ( y > 0 )
	       simno++;
	   }
	 prein = 1;
	 break;
       case 8 : *d++ = B[++j];
	 y = DN[j];
	 *d++ = B[++j];
	 *d++ = B[++j];
	 z = DN[j];
	 *d++ = B[++j];
	 *a++ = *letter++;
	 *a++ = ' ';
	 *a++ = *letter++;
	 *a++ = *letter++;
	 *c++ = mt1[one][y];
	 *c++ = '-';
	 *c++ = mt2[one][z = y * DNASIZE + z];
	 *c++ = mt3[one][z = z * DNASIZE + DN[j]];
	 if ( ! isexon )
	   { isexon = 1;
	   xstart = DP + j - 4;
	   flen = tlen = fmc = tmc = tpos = 0;
	   dbstart = AP + i - 1;
	   exlen = idnno = simno = exscore = 0;
	   }
	 *pp++ = aa3[da[z]];
	 *pp++ = ' ';
	 *pp++ = ' ';
	 *pp++ = ' ';
	 if ( isexon )
	   { if ( flen < CSLEN )
	     flen++;
	   if ( tlen < CSLEN )
	     tlen++;
	   else
	     tmc -= tring[tpos];
	   tring[tpos] = 0;
	   tpos = (tpos+1) % CSLEN;
	   exlen++;
	   exscore += (y = w[one][z]) - qr;
	   if ( y > 0 )
	     simno++;
	   }
	 prein = 1;
	 break;
       case 9 : *d++ = B[++j];
	 y = DN[j];
	 *d++ = B[++j];
	 z = DN[j];
	 *d++ = B[++j];
	 *d++ = B[++j];
	 *a++ = *letter++;
	 *a++ = *letter++;
	 *a++ = ' ';
	 *a++ = *letter++;
	 *c++ = mt1[one][y];
	 *c++ = mt2[one][z = y * DNASIZE + z];
	 *c++ = '-';
	 *c++ = mt3[one][z = z * DNASIZE + DN[j]];
	 if ( ! isexon )
	   { isexon = 1;
	   xstart = DP + j - 4;
	   flen = tlen = fmc = tmc = tpos = 0;
	   dbstart = AP + i - 1;
	   exlen = idnno = simno = exscore = 0;
	   }
	 *pp++ = aa3[da[z]];
	 *pp++ = ' ';
	 *pp++ = ' ';
	 *pp++ = ' ';
	 if ( isexon )
	   { if ( flen < CSLEN )
	     flen++;
	   if ( tlen < CSLEN )
	     tlen++;
	   else
	     tmc -= tring[tpos];
	   tring[tpos] = 0;
	   tpos = (tpos+1) % CSLEN;
	   exlen++;
	   exscore += (y = w[one][z]) - qr;
	   if ( y > 0 )
	     simno++;
	   }
	 prein = 1;
	 break;
       case 10 : *d++ = B[++j];
	 y = DN[j];
	 *d++ = B[++j];
	 *a++ = *letter++;
	 *a++ = *letter++;
	 *c++ = mt1[one][y];
	 *c++ = mt2[one][z = y * DNASIZE + DN[j]];
	 if ( isexon && *S >= gaplen )
	   { isexon = 0;
	   xend = DP + j - 1;
	   dbend = AP + i - 1;
	   if ( ! orient )
	     { xstart = dnalen - xstart + 1;
	     xend = dnalen - xend + 1;
	     }
	   if ( flen >= minexon )
	     { printf("\nExon %2d %8d %8d  Confidence: %3d %3d\n",
		      ++exonid, xstart, xend, (100*fmc)/flen, (100*tmc)/tlen);
	     fprintf(Btabp,"%s\t%s\t%d\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%f\t%f\t%d\t%d\t%d\t%s\t%d\t%s\t%d\n",
		     dhead+1, btab_date, dnalen, program,
		     name2, dblocus, xstart, xend, dbstart, dbend,
		     100.0 * idnno / exlen , 100.0 * simno / exlen,
		     exscore, chain_num, exonid, dbdesp, -1,
		     ( orient ? "Plus" : "Minus" ), psize);
	     }
	   }
	 z = z * DNASIZE + DN[*S+j+1];
	 *pp++ = aa3[da[z]];
	 *pp++ = ' ';
	 if ( isexon )
	   { exlen++;
	   exscore += (y = w[one][z]);
	   if ( y > 0 )
	     simno++;
	   }
	 prein = 1;
	 break;
       case 11 : *d++ = B[++j];
	 *a++ = *letter++;
	 *c++ = mt3[one][z];
	 if ( ! isexon )
	   { isexon = 1;
	   xstart = DP + j - 1;
	   flen = tlen = fmc = tmc = tpos = 0;
	   dbstart = AP + i - 1;
	   exlen = idnno = simno = exscore = 0;
	   }
	 *pp++ = ' ';
	 if ( isexon )
	   { if ( flen < CSLEN )
	     flen++;
	   if ( tlen < CSLEN )
	     tlen++;
	   else
	     tmc -= tring[tpos];
	   tring[tpos] = 0;
	   tpos = (tpos+1) % CSLEN;
	   }
	 prein = 1;
	 break;
       case 12 : *d++ = B[++j];
	 *a++ = *letter++;
	 *c++ = mt1[one][y = DN[j]];
	 if ( isexon && *S >= gaplen )
	   { isexon = 0;
	   xend = DP + j - 1;
	   dbend = AP + i - 1;
	   if ( ! orient )
	     { xstart = dnalen - xstart + 1;
	     xend = dnalen - xend + 1;
	     }
	   if ( flen >= minexon )
	     { printf("\nExon %2d %8d %8d  Confidence: %3d %3d\n",
		      ++exonid, xstart, xend, (100*fmc)/flen, (100*tmc)/tlen);
	     fprintf(Btabp,"%s\t%s\t%d\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%f\t%f\t%d\t%d\t%d\t%s\t%d\t%s\t%d\n",
		     dhead+1, btab_date, dnalen, program,
		     name2, dblocus, xstart, xend, dbstart, dbend,
		     100.0 * idnno / exlen , 100.0 * simno / exlen,
		     exscore, chain_num, exonid, dbdesp, -1,
		     ( orient ? "Plus" : "Minus" ), psize);
	     }
	   }
	 z = ( y = y * DNASIZE + DN[*S+j+1] ) * DNASIZE + DN[*S+j+2];
	 *pp++ = aa3[da[z]];
	 if ( isexon )
	   { exlen++;
	   exscore += (y = w[one][z]);
	   if ( y > 0 )
	     simno++;
	   }
	 prein = 1;
	 break;
       case 13 : *d++ = B[++j];
	 *d++ = B[++j];
	 *a++ = *letter++;
	 *a++ = *letter++;
	 *c++ = mt2[one][y];
	 *c++ = mt3[one][z];
	 if ( ! isexon )
	   { isexon = 1;
	   xstart = DP + j - 2;
	   flen = tlen = fmc = tmc = tpos = 0;
	   dbstart = AP + i - 1;
	   exlen = idnno = simno = exscore = 0;
	   }
	 *pp++ = ' ';
	 *pp++ = ' ';
	 if ( isexon )
	   { if ( flen < CSLEN )
	     flen++;
	   if ( tlen < CSLEN )
	     tlen++;
	   else
	     tmc -= tring[tpos];
	   tring[tpos] = 0;
	   tpos = (tpos+1) % CSLEN;
	   }
	 prein = 1;
	 break;
       default : break;
       }
     }
   else
     { if (op == 0)
       { op = *S++;
       gsize = 0;
       if ( isexon && (op >= gaplen || i == M || j == N) )
	 { isexon = 0;
	 xend = DP + j - 1;
	 dbend = AP + i - 1;
	 if ( ! orient )
	   { xstart = dnalen - xstart + 1;
	   xend = dnalen - xend + 1;
	   }
	 if ( flen >= minexon )
	   { printf("\nExon %2d %8d %8d  Confidence: %3d %3d\n",
		    ++exonid, xstart, xend, (100*fmc)/flen, (100*tmc)/tlen);
	   fprintf(Btabp,"%s\t%s\t%d\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%f\t%f\t%d\t%d\t%d\t%s\t%d\t%s\t%d\n",
		   dhead+1, btab_date, dnalen, program,
		   name2, dblocus, xstart, xend, dbstart, dbend,
		   100.0 * idnno / exlen , 100.0 * simno / exlen,
		   exscore, chain_num, exonid, dbdesp, -1,
		   ( orient ? "Plus" : "Minus" ), psize);
	   }
	 }
       g = 0;
       if ( isexon )
	 { if ( op > 0 )
	   { exscore -= q + r * op;
	   exlen += op / 3;
	   }
	 else
	   { exscore -= q * prein - r3 * op;
	   exlen -= op;
	   }
	 }
       }
     if (op > 0)
       { *d++ = B[++j];
       *a++ = ' ';
       *c++ = '-';
       op--;
       gsize++;
       if ( isexon && pst != 10 && pst != 12 )
	 { if ( ++g == 1 )
	   { if ( op > 1 )
	     { f = (DN[j] * DNASIZE + DN[j+1]) * DNASIZE + DN[j+2];
	     *pp++ = aa3[da[f]];
	     }
	   else
	     *pp++ = 'X';
	   }
	 else
	   *pp++ = ' ';
	 if ( g == 3 )
	   { g = 0;
	   if ( flen < CSLEN )
	     flen++;
	   if ( tlen < CSLEN )
	     tlen++;
	   else
	     tmc -= tring[tpos];
	   tring[tpos] = 0;
	   tpos = (tpos+1) % CSLEN;
	   }
	 }
       else
	 *pp++ = ' ';
       prein = 1;
       }
     else
       { letter = ttt[AN[++i]];
       *a++ = *letter++;
       *a++ = *letter++;
       *a++ = *letter++;
       *d++ = ' ';
       *d++ = ' ';
       *d++ = ' ';
       *c++ = '-';
       *c++ = '-';
       *c++ = '-';
       op++;
       *pp++ = ' ';
       *pp++ = ' ';
       *pp++ = ' ';
       gsize += 3;
       if ( isexon )
	 { if ( flen < CSLEN )
	   flen++;
	 if ( tlen < CSLEN )
	   tlen++;
	 else
	   tmc -= tring[tpos];
	 tring[tpos] = 0;
	 tpos = (tpos+1) % CSLEN;
	 }
       prein = 0;
       }
     }
   if (a >= ALINE+LINELEN || i >= M && j >= N)
     { *a = *d = *c = *pp = '\0';
     ls = a - ( ALINE+LINELEN );
     if ( ls > 0 )
       { strcpy(atemp, a = ALINE+LINELEN);
       strcpy(dtemp, DLINE+LINELEN);
       strcpy(ctemp, CLINE+LINELEN);
       strcpy(ptemp, PLINE+LINELEN);
       ALINE[LINELEN] = '\0';
       DLINE[LINELEN] = '\0';
       CLINE[LINELEN] = '\0';
       PLINE[LINELEN] = '\0';
       }
     flag = 0;
     if ( gsize < 2 * LINELEN )
       flag = 1;
     else
       if ( op >= 0 && op < LINELEN || op < 0 && -(3 * op) < LINELEN )
	 { if ( gcount > 0 )
	   { (void) printf("\n              ~~~ %d bp removed ~~~\n", gcount);
	   gcount = 0;
	   }
	 flag = 1;
	 }
       else
	 gcount += LINELEN;
     if ( flag )
       { (void) printf("\n%7d ",LINELEN*lines);
       for (e = ALINE+10; e <= a; e += 10)
	 (void) printf("    .    :");
       if (e <= a+5)
	 (void) printf("    .");
       mdp = orient ? dp : (dnalen - dp + 1);
       (void) printf("\nScript  %s\n%7d %s\n        %s\n%7d %s\n",
		     PLINE, mdp,DLINE,CLINE,ap,ALINE);
       }
     lines++;
     ap = AP + i;
     dp = DP + j;
     if ( ls > 0 )
       { for ( e = dtemp; *e ; )
	 if ( *e++ != ' ' )
	   dp--;
       for ( a = ALINE, e = atemp; *a++ = *e++; );
       for ( d = DLINE, e = dtemp; *d++ = *e++; );
       for ( c = CLINE, e = ctemp; *c++ = *e++; );
       for ( pp = PLINE, e = ptemp; *pp++ = *e++; );
       a--;
       d--;
       c--;
       pp--;
       if ( i >= M && j >= N)
	 { (void) printf("\n%7d ",LINELEN*lines++);
	 for (e = ALINE+10; e <= a; e += 10)
	   (void) printf("    .    :");
	 if (e <= a+5)
	   (void) printf("    .");
	 mdp = orient ? dp : (dnalen - dp + 1);
	 (void) printf("\nScript  %s\n%7d %s\n        %s\n%7d %s\n",
		       PLINE,mdp,DLINE,CLINE,ap,ALINE);
	 }
       }
     else
       { a = ALINE;
       d = DLINE;
       c = CLINE;
       pp = PLINE;
       }
     }
   }
 if ( isexon )
   { xend = DP + j - 1;
   dbend = AP + i - 1;
   if ( ! orient )
     { xstart = dnalen - xstart + 1;
     xend = dnalen - xend + 1;
     }
   if ( flen >= minexon )
     { printf("\nExon %2d %8d %8d  Confidence: %3d %3d\n",
	      ++exonid, xstart, xend, (100*fmc)/flen, (100*tmc)/tlen);
     fprintf(Btabp,"%s\t%s\t%d\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%f\t%f\t%d\t%d\t%d\t%s\t%d\t%s\t%d\n",
	     dhead+1, btab_date, dnalen, program,
	     name2, dblocus, xstart, xend, dbstart, dbend,
	     100.0 * idnno / exlen , 100.0 * simno / exlen,
	     exscore, chain_num, exonid, dbdesp, -1,
	     ( orient ? "Plus" : "Minus" ), psize);
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
