/*  This program takes the pairwise alignments of GAP2 and NAP, and
    merge the pairwise alignments into a multiple alignment with
    respect to the query sequence.
    
    Acknowledgments
    I thank Mark Adams, Xiaoying Lin, and Tony Kerlavage for suggestions:
    
    Usage: show output_of_gap2/nap [ output_of_gap2/nap ]  > output       
    
    where
    output_of_gap2:  output file produced by GAP2
    output_of_nap:   output file produced by NAP
    
    Note that the output files of GAP2 and NAP must not be modified.
    If the query sequence is over 10 MB in length, then the value for
    FIELDS needs to be modified. The current value is for a query coordinate
    of at most 7 digits in length.
*/

#include   <stdio.h>

#define   ACCLEN  2000		/* max length of a line */
#define   FIELDS   7		/* max length of a line */
#define  LINELEN   60	/* length of an output line */
#define  NAMELEN   20	/* length of an accession */

// Version info
static float VersionS = 1.53f;
static char * BuildS = "$Revision: 1.13 $";
static char *program;           /* program name */

typedef struct EXN   {
  int   start;		/* start position of exon */
  int   end;		/* end position of exon */
  int   cons;		/* confidence for start position */
  int   cone;		/* confidence for end position */
  struct EXN  *next;	/* pointer to next exon node */
}   extp,  *exptr;	/* type definition for exon node */
typedef struct ALTYPE   {
  char  acc[ACCLEN];	/* db seq accession */
  int   score;             /* score of the alignment */
  int   idn;		/* number of identities */
  int   len;		/* db seq length */
  int   qort;		/* query strand */
  int   lort;		/* db seq strand */
  int   lend;		/* left end position of the alignment */
  int   rend;		/* right end position of the alignment */
  int   dstart;            /* start position of the query */
  int   dend;              /* end position of the query */
  int   astart;            /* start position of the db seq */
  int   aend;              /* end position of the db seq */
  int   new;		/* 1, not put into work yet */
  int   dna;		/* 1, db seq is a cDNA */
  int   pred;		/* index of previous db seq on the algnt */
  exptr first;		/* first exon node */
  exptr last;		/* last exon node */
}  alntp, *alnptr;	/* type definition */
alnptr *  alnmt;        /* record for 2 types of alignments */
int  lastin;		/* last index of db seq on the algnt */

int  topsize;		/* size of the top array */
char *top;		/* holding the top part of the alignment */
char *mid;		/* holding the middle part of the alignment */
char *bot;		/* holding the bot part of the alignment */
int  alnum; 		/* number of alignments */
int  *kind;		/* indicates if there is an alignment of that type */
char qname[ACCLEN];	/* query seq name */

typedef struct OUT
{ char  line[LINELEN+1];/* one print line for db seq */
  char  mark[LINELEN+1];/* one print line for marks */
  char  *a;		/* pointer to slot in line */
  char  *b;		/* pointer to slot in mark */
  int   id;		/* index of the alignment */
  int   loc;		/* position of next aln symbol in bot */
  int   rend;		/* right end position of aln in bot */
  int   start;		/* coodinate of first db seq symbol on line */
  int   cod;		/* coodinate of next db seq symbol */
  int   lort;		/* db seq orient */
  char  name[NAMELEN+2];/* name of segment */
  short done;		/* indicates if segment is done */
  exptr first;		/* first exon node */
  exptr last;		/* last exon node */
}  row, *rowptr;	/* node for segment */
rowptr *work;			/* a set of working segments */

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

int main(argc, argv) int argc; char *argv[];
{ FILE *Ap, *ckopen();
 int ckopenCheck(char* name, char*); 
 char *ckalloc();			/* space-allocating function  */
 register int i, j;
 int  score;		/* alignment score */
 int  idn;		/* number of identities */
 int  len;		/* db seq length */
 char qsd[ACCLEN];	/* query strand indicator string */  
 char dbsd[ACCLEN];	/* db seq strand indicator string */  
 int  pert;		/* percent identity */
 char tc, bc;		/* for omitted long gaps */
 char  *end;		/* file end flag */
 char line[ACCLEN+1];			/* a line */
 int  temp;				/* temp variable */
 int  flag;		/* file processing flag */
 int  llen;		/* length of current line */
 int  alen;		/* length of current alignment */
 int  fsize;		/* total length of all alignments in a file */
 int  xx;		/* 0, file is an EST; 1, file is a protein */
 int  num;		/* number of alignments in current file */
 char **dbname;	/* database names */
 alnptr cur;		/* pointer to current record */
 int  aind;		/* argument index */
 char accn[ACCLEN];	/* accession string */
 int  first;		/* indicates if the line is first */
 unsigned int strlen();
 exptr exon;		/* pointer to exon node */
 int* fileChecker;    /*determines status of file (1=ok, 0=skip)*/
 int startAlnmtNumber = 0;
 program = argv[0];
 
 if ( argc < 2 )
   { fprintf(stderr,"Usage: %s Output_of_GAP2/NAP [ Output_of_GAP2/NAP ]\n\n", argv[0]);
   fprintf(stderr,"It takes any number of output files by GAP2 or NAP\n\n", argv[0]);
   exit(1);
   } 
 
 
 if (!strcmp(argv[1], "-help")  ||  !strcmp(argv[1], "-h"))
   {
     fprintf(stderr,"Usage: %s Output_of_GAP2/NAP [ Output_of_GAP2/NAP ]\n\n", argv[0]);
     fprintf(stderr,"It takes any number of output files by GAP2 or NAP\n\n", argv[0]);
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
     exit(1);
   }   
 
 dbname = ( char ** ) ckalloc( argc * sizeof(char *) );
 kind = ( int * ) ckalloc( argc * sizeof(int) );
 for ( aind = 0; aind < argc; aind++ )
   { kind[aind] = 0;
   dbname[aind] = ( char * ) ckalloc( (ACCLEN + 2) * sizeof(char) );
   }
 topsize = alnum = 0;
 fileChecker = ( int * ) ckalloc( argc * sizeof(int) );
 
 for ( aind = 1; aind < argc; aind++ )
   { 
     fileChecker[aind] = 0; //initialize.
     if (! ckopenCheck(argv[aind], "r")) {
       fprintf (stderr, "ERROR: File: %s cannot be read or is empty. Ignoring this file.\n", argv[aind]);
       fileChecker[aind] = 0;
       continue;
     }
     //opening file
     Ap = ckopen(argv[aind], "r");
     //read in 7 lines for some reason.
     for ( i = 0; i < 7; i++ ) {
       if ( ( end = fgets(line, ACCLEN, Ap)) == NULL )
	 break;
     }
     
     if ( end == NULL )
       { (void) fclose(Ap);
       continue;
       }
     
     //finding first position of 'Accession: %s' in the file.
     if ( (flag = sscanf(line, "Accession: %s", accn) ) == 0 )
       fatalf("Alignment file No. %d is not in the correct format", aind);
     
     xx = 0;
     startAlnmtNumber = alnum;
     for ( ; flag > 0 && end != NULL; alnum++ )
       { fgets(line, ACCLEN, Ap);
       for ( alen = 0; (end = fgets(line, ACCLEN, Ap)) != NULL; )
	 {
	   // should be a blank line after the score info.
	   if ( line[0] != '\n' ) {
	     alnum = startAlnmtNumber;
	     fprintf (stderr, "ERROR: Input file: %s is corrupt. Ignoring rest of file.\n", argv[aind]);
	     goto END; //break to outermost loop. ...no continue-label in C, and no exception handling either.
	   }
	   fgets(line, ACCLEN, Ap);
	   if ( (flag = sscanf(line, "Accession: %s", accn) ) > 0 )
	     break;
	   if ( sscanf(line, "Exon %d", &temp) > 0 ||
		sscanf(line, "              ~~~ %d bp removed ~~~", &temp) > 0 )
	     continue;
	   sscanf(line, "%d", &alen);
	   fgets(line, ACCLEN, Ap);
	   if ( sscanf(line, "Script%c", &tc) > 0 )
	     { xx = 1;
	     fgets(line, ACCLEN, Ap);
	     }
	   fgets(line, ACCLEN, Ap);
	   fgets(line, ACCLEN, Ap);
	   llen = strlen(line);
	 }
       topsize += alen + llen;
       }
     kind[aind-1] = xx;
     fileChecker[aind] = 1; //made it this far, should be OK.
   END: //end of outermost loop.
     (void) fclose(Ap);
   }
 
 if ( alnum == 0 )
   { printf("No output produced\n");
   exit(0);
   }
 
 top = ( char * ) ckalloc( (topsize + 1) * sizeof(char));
 mid = ( char * ) ckalloc( (topsize + 1) * sizeof(char));
 bot = ( char * ) ckalloc( (topsize + 1) * sizeof(char));
 alnmt = ( alnptr * ) ckalloc( (alnum + 1) * sizeof(alnptr));
 
 num = fsize = 0;
 for ( aind = 1; aind < argc; aind++ )
   { 
     if (fileChecker[aind] == 0) { //avoid bad file.
       continue;
     }
     Ap = ckopen(argv[aind], "r");
     if ( ( end = fgets(line, ACCLEN, Ap)) == NULL )
       { (void) fclose(Ap);
       continue;
       }
     if ( ( end = fgets(line, ACCLEN, Ap)) == NULL )
       { (void) fclose(Ap);
       continue;
       }
     if ( ( end = fgets(line, ACCLEN, Ap)) == NULL )
       { (void) fclose(Ap);
       continue;
       }
     if ( ( end = fgets(line, ACCLEN, Ap)) == NULL )
       { (void) fclose(Ap);
       continue;
       }
     strcpy(qname, line);
     if ( ( end = fgets(line, ACCLEN, Ap)) == NULL )
       { dbname[aind-1][0] = '\0';
       (void) fclose(Ap);
       continue;
       }
     strcpy(dbname[aind-1], line);
     if ( ( end = fgets(line, ACCLEN, Ap)) == NULL )
       { (void) fclose(Ap);
       continue;
       }
     if ( ( end = fgets(line, ACCLEN, Ap)) == NULL )
       { (void) fclose(Ap);
       continue;
       }
     if ( (flag = sscanf(line, "Accession: %s", accn) ) == 0 )
       fatalf("Alignment file No. %d is not in the correct format", aind);
     for ( ; flag > 0 && end != NULL; num++ )
       { cur = alnmt[num] = ( alnptr ) ckalloc( (unsigned) sizeof(alntp) );
       strcpy(cur->acc, line + 11);
       fgets(line, ACCLEN, Ap);
       if ( kind[aind-1] == 0 )
	 { sscanf(line, "Score: %d  Identity: %d/%d (%d%%)  Qstrand: %s  Lstrand: %s",
		  &score, &idn, &len, &pert, qsd, dbsd);
	 cur->lort = (strcmp(dbsd, "plus") == 0) ? 1 : 0;
	 cur->dna = 1;
	 }
       else
	 { sscanf(line, "Score: %d  Identity: %d/%d (%d%%)  Strand: %s",
		  &score, &idn, &len, &pert, qsd);
	 cur->dna = 0;
	 cur->lort = 1;
	 }
       cur->score = score;
       cur->idn = idn;
       cur->len = len;
       cur->qort = ( strcmp(qsd, "plus") == 0 ) ? 1 : 0;
       cur->lend = fsize;
       cur->new = 1;
       cur->first = cur->last = NULL;
       first = -1;
       for ( alen = 0; (end = fgets(line, ACCLEN, Ap)) != NULL; )
	 { fgets(line, ACCLEN, Ap);
	 if ( (flag = sscanf(line, "Accession: %s", accn) ) > 0 )
	   break;
	 if ( sscanf(line, "Exon %d", &temp) > 0 )
	   {  exon = ( exptr ) ckalloc( (unsigned) sizeof(extp) );
	   sscanf(line, "Exon %d %d %d %s %d %d", &temp,
		  &exon->start, &exon->end, accn, &exon->cons, &exon->cone);
	   if ( cur->last == NULL )
	     { cur->first = cur->last = exon;
	     exon->next = NULL;
	     }
	   else
	     { cur->last->next = exon;
	     cur->last = exon;
	     exon->next = NULL;
	     }
	   continue;
	   }
	 if ( sscanf(line, "              ~~~ %d bp removed ~~~", &temp) > 0 )
	   { if ( top[fsize-1] != ' ' )
	     { tc = '*';
	     bc = ' ';
	     }
	   else
	     { tc = ' ';
	     bc = '*';
	     }
	   for ( i = fsize, j = 0; j < temp ; j++ )
	     { top[i] = tc;
	     mid[i] = '-';
	     bot[i++] = ( xx && bc == '*' && (j % 3) ) ? ' ' : bc;
	     }
	   fsize += temp;
	   continue;
	   }
	 fgets(line, ACCLEN, Ap);
	 if ( sscanf(line, "Script%c", &tc) > 0 )
	   fgets(line, ACCLEN, Ap);
	 sscanf(line, "%d", &temp);
	 if ( first < 0 )
	   cur->dstart = temp;
	 for ( i = fsize, j = FIELDS+1; line[j] != '\n'; )
	   top[i++] = line[j++];
	 fgets(line, ACCLEN, Ap);
	 for ( i = fsize, j = FIELDS+1; line[j] != '\n'; )
	   mid[i++] = line[j++];
	 fgets(line, ACCLEN, Ap);
	 sscanf(line, "%d", &temp);
	 if ( first < 0 )
	   { cur->astart = temp;
	   first = 1;
	   }
	 for ( i = fsize, j = FIELDS+1; line[j] != '\n'; )
	   bot[i++] = line[j++];
	 fsize = i;
	 }
       cur->rend = fsize - 1;
       }
     (void) fclose(Ap);
   }
 printf("%s", qname);
 for ( aind = 1; aind < argc; aind++ )
   if ( dbname[aind-1][0] != '\0' )
     printf("%s", dbname[aind-1]);
 printf("Plus (+) denotes forward strand, and minus (-) reverse strand.\n");
 printf("Asterisks (*) denote bases not shown on pair wise alignments.\n");
 Sort();
 Show();
 
 exit(0);
 
}

Sort()
{ int  sorted;
 int  i, j;
 alnptr cur;
 int  swap;
 int  x, y; 
 
 for ( i = 0; i < alnum - 1; i++ )
   { sorted = 1;
   for ( j = alnum - 1; j > i; j--)
     { swap = 0;
     if ( (x = alnmt[j-1]->qort) == (y = alnmt[j]->qort) )
       { if ( x == 1 )
	 { if ( alnmt[j-1]->dstart > alnmt[j]->dstart )
	   swap = 1;
	 }
       else
	 { if ( alnmt[j-1]->dstart < alnmt[j]->dstart )
	   swap = 1;
	 }
       }
     else
       if ( x == 0 && y == 1 )
	 swap = 1;
     if ( swap )
       { cur = alnmt[j];
       alnmt[j] = alnmt[j-1];
       alnmt[j-1] = cur;
       sorted = 0;
       }
     }
   if ( sorted )
     break;
   }
}

/* Display multiple alignments */
Show()
{ register  int   i, j, k; /* index variables */
 char *ckalloc();	/* space-allocating function */
 int   n;		/* number of working segments */
 int   limit;		/* number of slots in work */
 int   col;		/* number of output columns prepared */
 short done;		/* tells if current group is done */
 char  c;		/* temp variable */
 rowptr t;		/* temp pointer */
 int    x;		/* temp variables */
 int  qstart;	/* coodinate of first query seq symbol on line */
 int  qcod;	/* coodinate of next query seq symbol */
 char qline[LINELEN+1];	/* one print line for the query seq */
 char *qpt;	/* pointer to slot in qline */
 char qnm[NAMELEN+2];	/* query name */
 int  qort;	/* query orientation */
 char qlet;	/* next query symbol */
 int  adv;	/* indicates if qcod was incremented */
 int  subt;	/* 1, next columns of all aln involve a query symbol */
 int  indel;	/* index of an aln whose next column has no query symbol */
 int  numdone; /* number of entries done in work */
 int  group;	/* alignment number */
 int  dsign;	/* 1, the line is all '*'s */
 int  dsize;	/* number of '*'s */
 exptr exon;	/* pointer to exon node */
 
 
 
 work = ( rowptr * ) ckalloc( alnum * sizeof(rowptr));
 n = limit = col = group = 0;
 for ( i = 0; i < alnum; i++ )
   if ( alnmt[i]->new )
     { done = 0;
     lastin = -1;
     Enter(&limit, &n, i, col);
     adv = 1;
     indel = 0;
     qpt = qline;
     qstart = qcod = alnmt[i]->dstart;
     qort = alnmt[i]->qort;
     sprintf(qnm, "%s", "Query");
     qnm[j = 5] = qort ? '+' : '-';
     qnm[++j] = '\0';
     printf("\n                              Alignment %d\n", ++group);
     dsign = 1;
     dsize = 0;
     while ( ! done )
       { if ( adv )
	 { t = work[n-1];
	 for ( j = t->id + 1; j < alnum; j++ )
	   if ( alnmt[j]->new )
	     if ( alnmt[j]->qort == qort )
	       { if ( qort )
		 { if ( alnmt[j]->dstart <= qcod )
		   Enter(&limit, &n, j, col);
		 else
		   break;
		 }
	       else
		 { if ( alnmt[j]->dstart >= qcod )
		   Enter(&limit, &n, j, col);
		 else
		   break;
		 }
	       }
	     else
	       break;
	 adv = 0;
	 }
       subt = 1;
       for ( ; indel < n; indel++ )
	 { t = work[indel];
	 if ( t->done )
	   continue;
	 if ( top[t->loc] == ' ' )
	   { subt = 0;
	   break;
	   }
	 }
       numdone = 0;
       if ( subt )
	 { qlet = '*';
	 for ( j = 0; j < n; j++ )
	   { t = work[j];
	   if ( t->done )
	     { *t->a++ = ' ';
	     *t->b++ = ' ';
	     numdone++;
	     }
	   else
	     { if ( (*t->a++ = bot[t->loc]) != ' ' )
	       t->cod += (t->lort) ? 1 : (-1);
	     *t->b++ = mid[t->loc];
	     if ( (c = top[t->loc++]) != '*' )
	       qlet = c;
	     if ( t->loc > t->rend )
	       { t->done = 1;
	       numdone++;
	       }
	     
	     }
	   }
	 indel = 0;
	 qcod += qort ? 1 : (-1);
	 *qpt++ = qlet;
	 adv = 1;
	 if ( qlet != '*' )
	   dsign = 0;
	 }
       else
	 { for ( j = 0; j < n; j++ )
	   { t = work[j];
	   /*   if ( j == indel )    used before 12/13/98 */
	   if ( top[t->loc] == ' ' )
	     { if ( (*t->a++ = c = bot[t->loc]) != ' ' )
	       { t->cod += (t->lort) ? 1 : (-1);
	       if ( c != '*' )
		 dsign = 0;
	       }
	     *t->b++ = mid[t->loc++];
	     if ( t->loc > t->rend )
	       { t->done = 1;
	       numdone++;
	       }
	     }
	   else
	     { *t->a++ = ' ';
	     *t->b++ = ' ';
	     if ( t->done )
	       numdone++;
	     }
	   }
	 *qpt++ = ' ';
	 }
       if ( numdone == n )
	 done = 1;
       col++;
       if ( col == LINELEN || done )	/* output */
	 { col = 0;
	 if ( ! dsign ) 
	   { if ( dsize > 0 )
	     { (void) printf("\n                              ~~ %d bp removed ~~~\n", dsize);
	     dsize = 0;
	     }
	   
           (void) printf("\n                              ");
	   for ( j = 0; j < LINELEN; j += 10 )
	     (void)   printf("    .    :");
	   (void) printf("\n");
	   *qpt = '\0';
	   (void) printf("%-11s           %7d %s\n", qnm, qstart, qline);
	   }
	 else
	   { dsize += LINELEN;
	   }
	 qstart = qcod;
	 qpt = qline;
	 for ( j = 0; j < n; j++ )
	   { t = work[j];
	   *t->a = '\0';
	   *t->b = '\0';
	   if ( dsign == 0 )
	     {
	       (void) printf("                              %s\n", t->mark);
	       (void) printf("%-21s %7d %s\n", t->name, t->start, t->line);
	     }
	   t->a = t->line;
	   t->b = t->mark;
	   t->start = t->cod;
	   }
	 if ( dsign == 0 )
	   for ( j = 0; j < n; j++ )
	     { t = work[j];
	     if ( (exon = t->first) != NULL &&
		  ( qort && exon->end < qcod || !qort && exon->end > qcod))
	       { 
		 (void) printf("%-21s         Exon  %8d %8d  Confidence: %3d %3d\n",
			       t->name, exon->start, exon->end, exon->cons, exon->cone);
		 t->first = exon->next;
		 free(exon);
	       }
	     }
	 if ( ! done )
	   { for ( k = j = n - 1; j >= 0; j-- )
	     if ( work[j]->done )
	       { t = work[j];
	       for ( x = j; x < k; x++ )
		 work[x] = work[x+1];
	       work[k--] = t;
	       }
	   n = k + 1;
	   }
	 else
	   n = 0;
	 dsign = 1;
	 }
       }
     for ( k = -1, j = lastin; j != -1; )
       { j = alnmt[x = j]->pred;
       alnmt[x]->pred = k;
       k = x;
       }
     (void) printf("\n");
     for ( ; k != -1 ; k = alnmt[k]->pred )
       (void) printf("%s\n", alnmt[k]->acc);
     }
}

/* enter an alignment into working set */
Enter(b, d, id, pos) int  *b, *d, id, pos;
{ int  i;
 char *ckalloc();		/* space-allocating function */
 rowptr  t;	/* new record for work */
 alnptr  cur;	/* record for an alignment */
 char    x;	/* temp var */
 
 if ( *b <= *d )
   { work[*b] = ( rowptr ) ckalloc( (int ) sizeof(row));
   *b += 1;
   }
 t = work[*d];
 *d += 1;
 t->a = t->line;
 t->b = t->mark;
 for ( i = 0; i < pos; i++ )
   *t->a++ = *t->b++ = ' ';
 t->id = id;
 cur = alnmt[id];
 cur->new = 0;
 cur->pred = lastin;
 lastin = id;
 t->loc = cur->lend;
 t->rend = cur->rend;
 t->cod = t->start = cur->astart;
 t->lort = cur->lort;
 for ( i = 0; (x = cur->acc[i]) != '\n' && i < NAMELEN; i++ )
   t->name[i] = x;
 if ( cur->dna )
   t->name[i] = cur->lort ? '+' : '-';
 else
   t->name[i] = ' ';
 t->name[++i] = '\0';
 t->done = 0;
 t->first = cur->first;
 t->last = cur->last;
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
  long where;
  
  if ((fp = fopen(name, mode)) == NULL)
    fatalf("Cannot open %s.", name);
  
  if((fseek(fp, 0L, SEEK_SET) != 0) || (fseek(fp, 0L, SEEK_END) != 0)) {
    fclose(fp);
    fatalf("%s could not be seeked.", name);
  }
  else {
    where = ftell(fp); 
    if(where == 0L) {
      fclose(fp);
      fatalf("%s is empty.", name);
    }
    if(fseek(fp, 0L, SEEK_SET) != 0) {
      fclose(fp);
      fatalf("%s could not be seeked.", name);
    }
  }
  return(fp);
}

int ckopenCheck(name, mode)
     char *name, *mode;
{
  FILE *fopen(), *fp;
  long where;
  
  if ((fp = fopen(name, mode)) == NULL)
    return(0);
  
  if((fseek(fp, 0L, SEEK_SET) != 0) || (fseek(fp, 0L, SEEK_END) != 0)) {
    fclose(fp);
    return(0);
  }
  else {
    where = ftell(fp); 
    if(where == 0L) {
      fclose(fp);
      return(0);
    }
    if(fseek(fp, 0L, SEEK_SET) != 0) {
      fclose(fp);
      return(0);
    }
  }
  return(1); //file is just fine.
}




/* ckalloc - allocate space; check for success */
char *ckalloc(amount)
     unsigned int amount;
{
  char *malloc(), *p;
  
  if ((p = malloc( (unsigned) amount)) == NULL)
    fatal("Ran out of memory.");
  return(p);
}
