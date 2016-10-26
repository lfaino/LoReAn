/*  A GLOBAL ALIGNMENT PROGRAM (GAP2):
    
The GAP2 program is an improved version of GAP. GAP2 is designed
for comparing a genomic DNA sequence with a cDNA sequence.
GAP2 uses the GT and AG denucleotide patterns to recognize splice sites.
Deletion of any region of length > k of the DNA sequence is given a
constant penalty of q + k r, where q and r are gap-open and gap-extension
penalties. The GAP2 program computes a global alignment of two sequences
without penalizing terminal gaps. It delivers the alignment in
linear space, so long sequences can be aligned. 

Users supply scoring parameters. In the simplest form, users just
provide 3 integers: ms, q and r, where ms is the score of a mismatch
and the score of an i-symbol indel is -(q + r * i). Each match
automatically receives score 10. This simple scoring scheme may be
used for DNA sequences. NOTE: all scores are integers.

In general, users can define an alphabet of characters appearing
in the sequences and a matrix that gives the substitution score
for each pair of symbols in the alphabet. The 127 ASCII characters
are eligible. The alphabet and matrix are given in a file, where
the first line lists the characters in the alphabet and the lower
triangle of the matrix comes next. An example file looks as follows:

ARNDC	       
13
-15  19
-10 -22  11
-20 -10 -20  18
-10 -20 -10 -20  12

Here the -22 at position (3,2) is the score of replacing N by R.
This general scoring scheme is useful for protein sequences where the
set of protein characters and Dayhoff matrix are specified in the file.

The GAP2 program is written in C and runs under Unix systems on
Sun workstations and under DOS systems on PCs.
We think that the program is portable to many machines.

Sequences to be analyzed are stored in separate files.
An input file contains all characters of a sequence, separated by
newline characters, in linear order. No other characters are allowed.
Since upper case and lower case characters are different, use the same
case consistently. A sample sequence file of 4 lines is shown below.

GAATTCTAATCTCCCTCTCAACCCTACAGTCACCCATTTGGTATATTAAA
GATGTGTTGTCTACTGTCTAGTATCCCTCAAGTAGTGTCAGGAATTAGTC
ATTTAAATAGTCTGCAAGCCAGGAGTGGTGGCTCATGTCTGTAATTCCAG
CACTGGAGAGGTAGAAGTG

To find the best alignment of two sequences in files A and B,
use a command of form

gap2  A  B  gs  ms  q  r > result

where gap2 is the name of the object code, gs is the minimum length
of any deletion gap receiving a constant gap penalty,
ms is a negative integer specifying mismatch weight, q and r are
non-negative integers specifying gap-open and gap-extend penalties,
respectively. Output alignment is saved in the file "result".
If GAP2 is used to align DNA and cDNA sequences, then A must
be the genomic DNA sequence and B the cDNA sequence.
If GAP2 is used to align protein sequences of substantially
different lengths, then A must be the longer sequence.

For using a scoring matrix defined in file S, use a command of form

gap2  A  B  gs  S  q  r > result

Note that ms is replaced by the file S.

Acknowledgments
I thank the following people for discussions and suggestions:
Mark Adams, Tony Kerlavage, Brendan Loftus, Steve Rounsley,
Granger Sutton, and Lixin Zhou.
Phil Green and LaDeana Hillier suggested use of the splice
site information (GT and AG) in the GAP2 program.
*/

#include   <stdio.h>

#define  DNASIZE    5   /* size of DNA alphabet: A, C, G, T, N */
#define  BONUS5     3   /* a bonus for a long gap starting with GT */
#define  BONUS3     3   /* a bonus for a long gap ending with AG */
#define  LINELEN   60   /* lenght of an output line */
#define  ACCLEN   400   /* maximum length of protein accession information */
#define  EXTRA     50   /* length of extra DNA region used */
#define  CSLEN     75   /* max length of 5' and 3' end of an exon examined */
#define  MINEXON   21   /* minimum exon length, must be < CSLEN */
#define  QGAPLEN   20	/* default value for gaplen */
#define  GAPOPEN    6	/* default value for q */
#define  GAPEXTN    1	/* default value for r */
#define  MATCH      2   /* default value for match */
#define  MISMAT    -2   /* default value for mismh */

// Version info
static float VersionS = 1.51f;
static char * BuildS = "$Revision: 1.14 $";

static int    di[128];          /* DNA letter to integer code */
static int    gtcode;           /* the integer code of GT */
static int    agcode;           /* the integer code of AG */
static int    bonus5;           /* BONUS5 */
static int    bonus3;           /* BONUS3 */
static int    dnalen;           /* length of query seq */
static int    nncode;           /* code of 'N' */
static int   minexon = MINEXON; /* min exon length */

/* Btab information */
static char *dhead;	        /* the name of the genomic DNA */
static char btab_date[100];	/* search date */
static char *program;           /* program name */
static char dblocus[ACCLEN+1];  /* locus name of database */
static int  dbstart;		/* start position of DB sequence */
static int  dbend;		/* end position of DB sequence */
static int  dblen;		/* length of DB sequence */
static int  idnno;		/* number of identities */
static int  simno;		/* number of similarities */
static int  exlen;		/* length of an exon alignment */
static int  exscore;		/* length of an exon alignment */
static char dbdesp[ACCLEN+1];	/* DB sequence description */
static FILE *Btabp;		/* Btab output file pointer */
static char fname[ACCLEN+1];	/* Btab file name */

static int match, mismh;	/* max and min substitution weights */
static char *name1, *name2;	/* names of sequence files    */
static int v[128][128];	/* substitution scores */
static int  q, r;       /* gap penalties */
static int  qr;         /* qr = q + r */
static int  gaplen;     /* minimum length for constant-cost insertion */
static int  pay;	/* constant-cost for long insertion */

static int *CC, *DD;			/* saving matrix scores */
static int *RR, *SS;		 	/* saving start-points */
static int *S, *S0;			/* saving operations for diff */
void swap (int*, int*);  //swap coords


/* The following definitions are for function diff() */

int  diff();
static int  zero = 0;				/* int type zero        */

#define gap(k)  ((k) <= 0 ? 0 : q+r*(k))	/* k-symbol indel score */

#define gap2(k)  ((k) <= 0 ? 0 : ((k) <= gaplen ? q+r*(k) : pay))
/* k-symbol insertion score */

static int *sapp;				/* Current script append ptr */
static int  last;				/* Last script op appended */

static int no_mat; 				/* number of matches */ 
static int no_mis; 				/* number of mismatches */ 
static int al_len; 				/* length of alignment */
static int term_gap5; 			/* no. of amono acids in 5' terminal gap */
static int term_gap3; 			/* no. of amono acids in 3' terminal gap */
/* Append "Delete k" op */
#define DEL(k)				\
{ al_len += k;				\
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
{ al_len += k;				\
  if (last > 0)				\
    last = sapp[-1] += (k);		\
  else					\
    last = *sapp++ = (k);		\
}
/* Append "Replace" op */
#define REP 				\
{ last = *sapp++ = 0; 			\
  al_len += 1;				\
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
 char  *D, *D0, *P, *P0;		/* Storing two sequences */
 int   *DN, *DN0;			/* Integer code of the sequences */
 int  symbol;				/* The next character */
 FILE *Bp, *Ap, *ckopen();
 char *ckalloc();			/* space-allocating function  */
 register int i, j;
 int  size;				/* size of alphabet */
 char cline[ACCLEN+1];		        /* a line of chain file */
 int  aalen;				/* maximu database sequence length */
 int  psize;				/* current protein length */
 char accn[ACCLEN+1];			/* current protein accession */
 char accncpy[ACCLEN+1];	        /* a copy of the accession   */ 
 char *identifier;                     /* the identifier from the header of */
 /* of a cDNA sequence */  
 const char *delim;                    /* the delimiter for tokenizing the */
 /* cDNA sequence header */ 
 int  dstart, dend;			/* DNA positions */
 int  dstart0, dend0;			/* DNA positions */
 int  score;				/* alignment score */
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
 char   *aa;			/* pointer to current protein seq */
 int  mat0, ss0;		/* for other strand */
 int  term3, term5;		/* for other strand */
 long  cur_time;
 
 time(&cur_time);
 // strftime(btab_date, sizeof(btab_date), "%b %e %Y", localtime(&cur_time));
 program = argv[0];
 
 gaplen = QGAPLEN;
 q = GAPOPEN;
 r = GAPEXTN;
 match = MATCH;
 mismh = MISMAT;
 
 if ( argc < 2 )
   { fprintf(stderr,"Usage: %s DNA_Seq EST_Database Chain_File [Btab_file] [options]\n\n", argv[0]);
   fprintf(stderr,"  DNA_Seq        file of one query DNA in FASTA format\n");
   fprintf(stderr,"  EST_Database   file of EST database in FASTA format\n");
   fprintf(stderr,"  Chain_File     file of coordinates produced by EXT or Filter\n");
   fprintf(stderr,"  Btab_file      file with the btab output (default=DNA_Seq.gap2.btab)\n");
   fprintf(stderr,"Options (default values):\n");
   fprintf(stderr,"  -e  N  specify gap extension penalty N > 0 (1)\n");
   fprintf(stderr,"  -m  N  specify match score N > 0 (2)\n");
   fprintf(stderr,"  -o  N  specify gap open penalty N >= 0 (6)\n");
   fprintf(stderr,"  -q  N  specify gap length N >= 5 for constant penalty (20)\n");
   fprintf(stderr,"  -s  N  specify mismatch score of N < 0 (-2)\n");
   fprintf(stderr,"  -x  N  specify min exon length N > 0 and < %d (21)\n", CSLEN);
   exit(1);
   }
 
 if (!strcmp(argv[1], "-help")  ||  !strcmp(argv[1], "-h"))
   {
     fprintf(stderr,"Usage: %s DNA_Seq EST_Database Chain_File [Btab_file] [options]\n\n", argv[0]);
     fprintf(stderr,"  DNA_Seq        file of one query DNA in FASTA format\n");
     fprintf(stderr,"  EST_Database   file of EST database in FASTA format\n");
     fprintf(stderr,"  Chain_File     file of coordinates produced by EXT or Filter\n");
     fprintf(stderr,"  Btab_file      file with the btab output (default=DNA_Seq.gap2.btab)\n");
     fprintf(stderr,"Options (default values):\n");
     fprintf(stderr,"  -e  N  specify gap extension penalty N > 0 (1)\n");
     fprintf(stderr,"  -m  N  specify match score N > 0 (2)\n");
     fprintf(stderr,"  -o  N  specify gap open penalty N >= 0 (6)\n");
     fprintf(stderr,"  -q  N  specify gap length N >= 5 for constant penalty (20)\n");
     fprintf(stderr,"  -s  N  specify mismatch score of N < 0 (-2)\n");
     fprintf(stderr,"  -x  N  specify min exon length N > 0 and < %d (21)\n", CSLEN);
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
     fprintf(stderr, "default match score = %d\n", match); 
     fprintf(stderr, "default mismatch score = %d\n", mismh); 
     exit(1);
   }   
 
 if ( argc < 4 )
   { fprintf(stderr,"Usage: %s DNA_Seq EST_Database Chain_File [Btab_file] [options]\n\n", argv[0]);
   fprintf(stderr,"  DNA_Seq        file of one query DNA in FASTA format\n");
   fprintf(stderr,"  EST_Database   file of EST database in FASTA format\n");
   fprintf(stderr,"  Chain_File     file of coordinates produced by EXT or Filter\n");
   fprintf(stderr,"  Btab_file      file with the btab output (default=DNA_Seq.gap2.btab)\n");
   fprintf(stderr,"Options (default values):\n");
   fprintf(stderr,"  -e  N  specify gap extension penalty N > 0 (1)\n");
   fprintf(stderr,"  -m  N  specify match score N > 0 (2)\n");
   fprintf(stderr,"  -o  N  specify gap open penalty N >= 0 (6)\n");
   fprintf(stderr,"  -q  N  specify gap length N >= 5 for constant penalty (20)\n");
   fprintf(stderr,"  -s  N  specify mismatch score of N < 0 (-2)\n");
   fprintf(stderr,"  -x  N  specify min exon length N > 0 and < %d (21)\n", CSLEN);
   exit(1);
   }
 sprintf(fname, "%s.gap2.btab", argv[1]);
 if ( argc > 4 )
   { temp = 4;
   if ( argv[4][0] != '-' ) strcpy(fname, argv[temp++]);
   for ( i = temp; i < argc - 1; i += 2 )
     { if ( argv[i][0] != '-' )
       fatal("Each option must begin with a dash (-)");
     (void) sscanf(argv[i+1],"%d", &temp);
     switch ( argv[i][1] )
       { case 'e' : if ( temp <= 0 )
	 fatal("Number for gap extension must be a positive integer");
	   r = temp;
	   break;
       case 'm' : if ( temp <= 0 )
	 fatal("Number for match must be a positive integer");
	 match = temp;
	 break;
       case 'o' : if ( temp < 0 )
	 fatal("Number for gap open must be a non-negative integer");
	 q = temp;
	 break;
       case 'q' : if ( temp < 5 )
	 fatal("Number for gap length must be a positive integer >= 5");
	 gaplen = temp;
	 break;
       case 's' : if ( temp >= 0 )
	 fatal("Number for mismatch must be a negative integer");
	 mismh = temp;
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
   fprintf(stderr, "No chain was reported by DDS");
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
 
 nncode = DNASIZE - 1;         /* code of 'N' */
 for ( i = 0; i < 128; i++ )
   di[i] = nncode;
 di['a'] = di['A'] = 0;
 di['c'] = di['C'] = 1;
 di['g'] = di['G'] = 2;
 di['t'] = di['T'] = 3;
 gtcode = di['G'] * DNASIZE + di['T'];
 agcode = di['A'] * DNASIZE + di['G'];
 
 /* determine the heading length */
 Ap = ckopen(argv[1], "r");
 if ( (symbol = getc(Ap) ) != '>' )
   fatal("The DNA sequence must be in the FASTA format, starting with '>'");
 for (size = 2; ( symbol = getc(Ap)) != EOF ; size++ )
   if ( symbol == '\n' )
     break;
 (void) fclose(Ap);
 name1 = argv[1];
 
 /* allocate space for D */
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
   fatal("The DNA sequence for GAP2 must be identical to that for DDS");
 /* allocate space for D0 */
 D0 = ( char * ) ckalloc( (dnalen + 1) * sizeof(char));
 DN0 = ( int * ) ckalloc( (dnalen + 1) * sizeof(int));
 /* Set D0 and DN0 to the reverse complement of D and DN */
 for ( i = 1, j = M; i <= M; i++, j-- )
   switch( D[i] )
     { case 'A' : DN0[j] = 3; D0[j] = 'T'; break;
     case 'a' : DN0[j] = 3; D0[j] = 't'; break;
     case 'C' : DN0[j] = 2; D0[j] = 'G'; break;
     case 'c' : DN0[j] = 2; D0[j] = 'g'; break;
     case 'G' : DN0[j] = 1; D0[j] = 'C'; break;
     case 'g' : DN0[j] = 1; D0[j] = 'c'; break;
     case 'T' : DN0[j] = 0; D0[j] = 'A'; break;
     case 't' : DN0[j] = 0; D0[j] = 'a'; break;
     default  : DN0[j] = 4; D0[j] = 'N'; break;
     }
 
 pay = q + r * gaplen;
 qr = q + r;
 bonus5 = BONUS5;
 bonus3 = BONUS3;
 /* set match and mismatch weights */
 for ( i = 0; i < 128 ; i++ )
   for ( j = 0; j < 128 ; j++ )
     if (i == j )
       v[i][j] = match;
     else
       v[i][j] = mismh;
 v['N']['N'] = mismh;
 v['n']['n'] = mismh;
 v['A']['a'] = v['a']['A'] = match;
 v['C']['c'] = v['c']['C'] = match;
 v['G']['g'] = v['g']['G'] = match;
 v['T']['t'] = v['t']['T'] = match;
 
 /* allocate space for all vectors */
 j = (dnalen + 1) * sizeof(int);
 CC = ( int * ) ckalloc(j);
 DD = ( int * ) ckalloc(j);
 RR = ( int * ) ckalloc(j);
 SS = ( int * ) ckalloc(j);
 i = (aalen + 1) * sizeof(int);
 S = ( int * ) ckalloc(i + j);
 S0 = ( int * ) ckalloc(i + j);
 
 /* allocate space for P */
 P = ( char * ) ckalloc( (aalen + 1) * sizeof(char));
 P0 = ( char * ) ckalloc( (aalen + 1) * sizeof(char));
 Bp = ckopen(name2 = argv[2], "r");
 symbol = getc(Bp);
 (void) printf("Max Match   Min Mismatch   Gap-Open Penalty   Gap-Extension Penalty\n");
 (void) printf("   %d          %d              %d                  %d\n\n", match, mismh, q, r);
 (void) printf("                 Query Sequence : %s\n", name1);
 (void) printf("                 Database: %s\n", name2);
 /* process each database sequence */
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
   
   aa = NULL;
   for ( i = 0; i < chlen; i++ )
     { apt = chains[i].acc;
     //strcpy(accncpy, accn);
     // delim = " ";
     // identifier = (char *) strtok(accncpy, delim);
     /* Here is where the accessions in chain file are compared to the database sequences*/
     // if(strcmp(apt,identifier) != 0)
     // {
     //	continue;
     // }
     
     for ( j = 0; apt[j] != '\0' && apt[j] == accn[j]; j++ )
       ;
     if ( apt[j] != '\0')  //means there must have been character differences in comparison
       continue;
     // if at this point, apt-acc matched accn-acc for full length of apt-acc
     // still don't know if there's an extension on accn-acc that would make it a different
     // accession from apt-acc (ie. acc555 != acc5559).  Must make sure that accn-acc
     // trails with whitespace.
     
     if (! ((accn[j] == ' ') || (accn[j] == '\t') 
	    || (accn[j] == '\n') || (accn[j] == '\f')  
	    ||  accn[j] == '\0'))
       continue;
     
     if ( aa == NULL )
       { aa = ( char * ) ckalloc( (psize + 1) * sizeof(char));
       for ( j = 1; j <= psize; j++ )
	 aa[j] = P[j];
       }
     chains[i].aa = aa;
     chains[i].alen = psize;
     strcpy(chains[i].acc, accn);
     }
   }
 (void) fclose(Bp);
 Btabp = ckopen(fname, "w");
 for ( i = 0; i < chlen; i++ )
   { if ( (aa = chains[i].aa) == NULL )
     fatal("No protein seq was found for one chain.");
   dblen = psize = chains[i].alen;
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
   for ( j = 1, temp = psize; j <= psize; j++, temp-- )
     switch( P[j] = aa[j] )
       { case 'A' : P0[temp] = 'T'; break;
       case 'a' : P0[temp] = 't'; break;
       case 'C' : P0[temp] = 'G'; break;
       case 'c' : P0[temp] = 'g'; break;
       case 'G' : P0[temp] = 'C'; break;
       case 'g' : P0[temp] = 'c'; break;
       case 'T' : P0[temp] = 'A'; break;
       case 't' : P0[temp] = 'a'; break;
       default  : P0[temp] = 'N'; break;
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
   if ( chains[i].ort )
     { dstart0 = dnalen - dend + 1;
     dend0 = dnalen - dstart + 1;
     }
   else
     { dstart0 = dstart;
     dend0 = dend;
     dstart = dnalen - dend0 + 1;
     dend = dnalen - dstart0 + 1;
     }
   sapp = S;
   last = 0;
   al_len = 0;
   no_mat = 0;
   no_mis = 0;
   term_gap5 = 0;
   term_gap3 = 0;
   temp = dend-dstart+1;
   printf("\nAccession: %s\n", chains[i].acc);
   if ( chains[i].ort )
     { ss0 = diff2(P,D+dstart-1,psize,temp,DN+dstart-1,q,q,0,0,0,0);
     mat0 = no_mat;
     term5 = term_gap5;
     term3 = term_gap3;
     sapp = S0;
     last = 0;
     al_len = 0;
     no_mat = 0;
     no_mis = 0;
     term_gap5 = 0;
     term_gap3 = 0;
     score = diff2(P0,D0+dstart0-1,psize,temp,DN0+dstart0-1,q,q,0,0,0,0);
     if ( ss0 >= score )
       { if ( term5 + term3 < psize )
	 size = psize - term5 - term3;
       else
	 size = psize;
       printf("Score: %d  Identity: %d/%d (%d%%)  Qstrand: %s  Lstrand: %s\n",
	      ss0, mat0, size, (100*mat0)/size, "plus", "plus");
       display(P,D+dstart-1,psize,temp,S,1,dstart,1,1, i);
       }
     else
       { if ( term_gap5 + term_gap3 < psize )
	 size = psize - term_gap5 - term_gap3;
       else
	 size = psize;
       printf("Score: %d  Identity: %d/%d (%d%%)  Qstrand: %s  Lstrand: %s\n",
	      score, no_mat, size, (100*no_mat)/size, "minus", "minus");
       display(P0,D0+dstart0-1,psize,temp,S0,1,dstart0,0,0, i);
       }
     }
   else
     { ss0 = diff2(P,D0+dstart0-1,psize,temp,DN0+dstart0-1,q,q,0,0,0,0);
     mat0 = no_mat;
     term5 = term_gap5;
     term3 = term_gap3;
     sapp = S0;
     last = 0;
     al_len = 0;
     no_mat = 0;
     no_mis = 0;
     term_gap5 = 0;
     term_gap3 = 0;
     score = diff2(P0,D+dstart-1,psize,temp,DN+dstart-1,q,q,0,0,0,0);
     if ( ss0 >= score )
       { if ( term5 + term3 < psize )
	 size = psize - term5 - term3;
       else
	 size = psize;
       printf("Score: %d  Identity: %d/%d (%d%%)  Qstrand: %s  Lstrand: %s\n",
	      ss0, mat0, size, (100*mat0)/size, "minus", "plus");
       display(P,D0+dstart0-1,psize,temp,S,1,dstart0,1,0, i);
       }
     else
       { if ( term_gap5 + term_gap3 < psize )
	 size = psize - term_gap5 - term_gap3;
       else
	 size = psize;
       printf("Score: %d  Identity: %d/%d (%d%%)  Qstrand: %s  Lstrand: %s\n",
	      score, no_mat, size, (100*no_mat)/size, "plus", "minus");
       display(P0,D+dstart-1,psize,temp,S0,1,dstart,0,1, i);
       }
     }
   }
 
 exit(0);
 
}

/* diff2(A,B,M,N,DN,tb,te,sc,sr,ec,er) returns the score of an optimum conversion
   between A[1..M] and B[1..N] that begins(ends) with a delete if tb(te) is zero
   and appends such a conversion to the current script. If sc = 0, then
   the beginning deletion is not penalized; if sr = 0, the beginning insertion is
   not penalized; if ec = 0, the ending deletion is not charged; if er = 0;
   then the ending insertion is not charged. Any insertion of length at least
   gaplen is given a constant cost */

int diff2(A,B,M,N,DN,tb,te,sc,sr,ec,er)
     char *A, *B; int M, N, DN[]; int tb, te, sc, sr, ec, er;
{ int   midi, midj, type;	/* Midpoint, type, and cost */
 int midc;
 int  ss,cc;
 
 { register int   i, j;
 register int c, e, d, s;
 int t, *va;
 int  g, temp;
 int  x, y, z;
 
 /* Boundary cases: M <= 1 or N == 0 */
 
 if (N <= 0)
   { if (M > 0) DEL(M)
		  if ( !sc || !ec )
		    return 0;
		  else
		    return - gap(M);
   }
 if (M <= 1)
   { if (M <= 0)
     { INS(N);
     if ( !sr || !er )
       return 0;
     else
       return - gap2(N);
     }
   midc = - (sc * (tb + r) + er * (z = gap2(N) ) );
   midj = -1;
   if ( midc < ( c =  - (ec * (te + r) + sr * z ) ) )
     { midc = c;
     midj = 0;
     }
   if ( N > 2 && gtcode == DN[1] * DNASIZE + DN[2] )
     x = 1;
   else
     x = 0;
   if ( N > 2 && agcode == DN[N-1] * DNASIZE + DN[N] )
     z = 1;
   else
     z = 0;
   va = v[A[1]];
   for (j = 1; j <= N; j++)
     { y = N - j;
     temp = er * gap2(y);
     if ( er && y > gaplen )
       { if ( gtcode == DN[j+1] * DNASIZE + DN[j+2] )
	 temp -= bonus5;
       if ( z )
	 temp -= bonus3;
       }
     y = j - 1;
     temp += sr * gap2(y);
     if ( sr && y > gaplen )
       { if ( x )
	 temp -= bonus5;
       if ( agcode == DN[y-1] * DNASIZE + DN[y] )
	 temp -= bonus3;
       }
     c = va[B[j]] - temp;
     if (c > midc)
       { midc = c;
       midj = j;
       }
     }
   if (midj == -1)
     { DEL(1) INS(N) }
   else
     if (midj == 0)
       { INS(N) DEL(1) }
     else
       { if (midj > 1) INS(midj-1)
			 REP
			 if ( di[A[1]] == di[B[midj]] )
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
 for (i = 1; i <= midi; i++)
   { s = CC[0];
   CC[0] = c = t = (t-r) * sc;
   e = t-q;
   g = t - pay;
   va = v[A[i]];
   for (j = 1; j <= N; j++)
     { if ((c = c - qr) > (e = e - r)) e = c;
     if ( j == N && !ec )
       { if ((c = CC[j] ) > (d = DD[j] )) d = c;}
     else
       if ((c = CC[j] - qr) > (d = DD[j] - r)) d = c;
     c = s+va[B[j]];
     if (c < d) c = d;
     if (c < e) c = e;
     if ( (z = j - gaplen - 1) >= 0 )
       { x = gtcode == (DN[z+1] * DNASIZE + DN[z+2]) ? bonus5 : 0;
       if ( g < ( temp = CC[z] + x - pay ) )
	 g = temp;
       x = agcode == (DN[j-1] * DNASIZE + DN[j]) ? bonus3 : 0;
       if ( c < (y = g + x) ) c = y;
       }
     s = CC[j];
     CC[j] = c;
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
   { s = RR[N];
   RR[N] = c = t = (t-r) * ec;
   g = t - pay;
   e = t-q;
   va = v[A[i+1]];
   for (j = N-1; j >= 0; j--)
     { if ((c = c - qr) > (e = e - r)) e = c;
     if ( !j && !sc )
       { if ((c = RR[j] ) > (d = SS[j] )) d = c;}
     else
       if ((c = RR[j] - qr) > (d = SS[j] - r)) d = c;
     c =  s+va[B[j+1]];
     if (c < d) c = d;
     if (c < e) c = e;
     if ( (z = j + gaplen + 1) <= N )
       { x = agcode == (DN[z-1] * DNASIZE + DN[z]) ? bonus3 : 0;
       if ( g < ( temp = RR[z] + x - pay ) )
	 g = temp;
       x = gtcode == (DN[j+1] * DNASIZE + DN[j+2]) ? bonus5 : 0;
       if ( c < (y = g + x) ) c = y;
       }
     s = RR[j];
     RR[j] = c;
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
   { (void) diff2(A,B,midi,midj,DN,tb,q,sc,sr,cc,1);
   (void) diff2(A+midi,B+midj,M-midi,N-midj,DN+midj,q,te,ss,1,ec,er);
   }
 else
   { (void) diff2(A,B,midi,midj,DN,tb,zero,sc,sr,cc,1);
   (void) diff2(A+midi,B+midj,M-midi,N-midj,DN+midj,zero,te,ss,1,ec,er);
   }
 return midc;
}
/* Alignment display routine */

static char ALINE[LINELEN+1], BLINE[LINELEN+1], CLINE[LINELEN+1];

display(A,B,M,N,S,AP,BP,dbort,orient, chain_num)
     char A[], B[]; int M, N, chain_num; int S[], AP, BP, dbort, orient;
{ register char *a, *b, *c;
 register int   i,  j, op;
 int  lines, ap, bp;
 int  mbp;		/* position relative to the plus strand query */
 int  map;		/* position relative to the plus strand dbseq */
 int isexon;		/* 1 if the current region is exon */
 int xstart;		/* start position of an exon */
 int xend;		/* end position of an exon */
 int exonid;		/* exon number */
 int flen, tlen;	/* lengths of 5' and 3' ends of an exon */
 int fmc, tmc;		/* number of matches in the ends */
 int tring[CSLEN];	/* match-mismatch vector for 3' end */
 int tpos;		/* current position of tring */
 int temp, zz;		/* temp var */
 int gsize;		/* current gap length */
 int gcount;		/* length of a gap removed */
 int flag;		/* 1, print */
 
 i = j = op = lines = isexon = exonid = gsize = gcount = 0;
 ap = AP;
 bp = BP;
 a = ALINE;
 b = BLINE;
 c = CLINE;
 chain_num++; // start at 1 instead of zero.
 while (i < M || j < N)
   { if (op == 0 && *S == 0)
     { op = *S++;
     *a = A[++i];
     *b = B[++j];
     *c++ = ((temp = di[*a++]) == (zz = di[*b++]) && temp != nncode) ? '|' : ' ';
     gsize = 0;
     if ( ! isexon )
       { isexon = 1;
       xstart = BP + j - 1;
       dbstart = AP + i - 1;
       flen = tlen = fmc = tmc = tpos = exlen = idnno = simno = exscore = 0;
       }
     if ( isexon )
       { exlen++;
       if ( *(c-1) == '|' )
	 { idnno++;
	 simno++;
	 exscore += match;
	 }
       else
	 { exscore += mismh;
	 if ( temp == nncode ||  zz == nncode )
	   simno++;
	 }
       if ( flen < CSLEN )
	 { flen++;
	 if ( *(c-1) == '|' )
	   fmc += 1;
	 }
       if ( tlen < CSLEN )
	 tlen++;
       else
	 tmc -= tring[tpos];
       if ( *(c-1) == '|' )
	 { tmc += 1;
	 tring[tpos] = 1;
	 }
       else
	 tring[tpos] = 0;
       tpos = (tpos+1) % CSLEN;
       }
     }
   else
     { if (op == 0)
       { op = *S++;
       gsize = 0;
       if ( isexon && (op >= gaplen || i == M || j == N ) )
	 { isexon = 0;
	 xend = BP + j - 1;
	 dbend = AP + i - 1;
	 if ( ! orient )
	   { xstart = dnalen - xstart + 1;
	   xend = dnalen - xend + 1;
	   }
	 if ( ! dbort )
	   { dbstart = dblen - dbstart + 1;
	   dbend = dblen - dbend + 1;
	   }
	 if ( flen >= minexon )
	   { printf("\nExon %2d %8d %8d  Confidence: %3d %3d\n",
		    ++exonid, xstart, xend, (100*fmc)/flen, (100*tmc)/tlen);
	   //db coords should always be in the forward orientation.
	   if (dbstart > dbend) {
	     swap(&dbstart, &dbend);
	     swap(&xstart, &xend);
	   }
	   fprintf(Btabp,"%s\t%s\t%d\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%f\t%f\t%d\t%d\t%d\t%s\t%d\t%s\t%d\n",
		   dhead+1, btab_date, dnalen, program,
		   name2, dblocus, xstart, xend, dbstart, dbend,
		   100.0 * idnno / exlen , 100.0 * simno / exlen,
		   exscore, chain_num, exonid, dbdesp, -1,
		   ( orient ? "Plus" : "Minus" ), dblen);
	   }
	 }
       if ( isexon )
	 exscore -= q;
       }
     if (op > 0)
       { *a++ = ' ';
       *b++ = B[++j];
       op--;
       }
     else
       { *a++ = A[++i];
       *b++ = ' ';
       op++;
       }
     *c++ = '-';
     gsize++;
     if ( isexon )
       { exlen++;
       exscore -= r;
       simno++;
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
   if (a >= ALINE+LINELEN || i >= M && j >= N)
     { *a = *b = *c = '\0';
     flag = 0;
     if ( gsize < 2 * LINELEN )
       flag = 1;
     else
       if ( op >= 0 && op < LINELEN || op < 0 && -op < LINELEN )
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
       for (b = ALINE+10; b <= a; b += 10)
	 (void) printf("    .    :");
       if (b <= a+5)
	 (void) printf("    .");
       mbp = orient ? bp : (dnalen - bp + 1);
       map = dbort ? ap : (dblen - ap + 1);
       (void) printf("\n%7d %s\n        %s\n%7d %s\n",mbp,BLINE,CLINE,map,ALINE);
       }
     lines++;
     ap = AP + i;
     bp = BP + j;
     a = ALINE;
     b = BLINE;
     c = CLINE;
     }
   }
 if ( isexon )
   { xend = BP + j - 1;
   dbend = AP + i - 1;
   if ( ! orient )
     { xstart = dnalen - xstart + 1;
     xend = dnalen - xend + 1;
     }
   if ( ! dbort )
     { dbstart = dblen - dbstart + 1;
     dbend = dblen - dbend + 1;
     }
   if ( flen >= minexon )
     { printf("\nExon %2d %8d %8d  Confidence: %3d %3d\n",
	      ++exonid, xstart, xend, (100*fmc)/flen, (100*tmc)/tlen);
     //db coords should always be in the forward orientation.
     if (dbstart > dbend) {
       swap(&dbstart, &dbend);
       swap(&xstart, &xend);
     }
     fprintf(Btabp,"%s\t%s\t%d\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%f\t%f\t%d\t%d\t%d\t%s\t%d\t%s\t%d\n",
	     dhead+1, btab_date, dnalen, program,
	     name2, dblocus, xstart, xend, dbstart, dbend,
	     100.0 * idnno / exlen , 100.0 * simno / exlen,
	     exscore, chain_num, exonid, dbdesp, -1,
	     ( orient ? "Plus" : "Minus" ), dblen);
     }
   }
}


void swap (coord1, coord2)
     int *coord1, *coord2;
{
  int tmp_swap;
  tmp_swap = *coord1;
  *coord1 = *coord2;
  *coord2 = tmp_swap;
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
