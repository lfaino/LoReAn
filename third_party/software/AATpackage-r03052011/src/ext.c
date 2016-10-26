/* This program (EXT) takes output of DPS or DDS and produces a summary of chain
   coordinates for use by GAP2 or NAP.
   It has an option to remove any match whose description contains
   a specified word in proper case such as "ALU".
*/

#define    ACCLEN  400	/* max length of a line */

#include   <stdio.h>

// Version info
static float VersionS = 1.51f;
static char * BuildS = "$Revision: 1.7 $";
static char * program = "";         // program name 

typedef struct BKTYPE           /* bucket type */
	{  char *line;          /* a line */
  	   struct  BKTYPE  *next; }  bbk, *bbkptr;
bbkptr   *bucket, item;	/* bucket */

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
  fprintf(stderr,"ext - Use dps or dds output to produce chain coordinates for use by gap2 or nap\n\n");
  fprintf(stderr,"Usage: ext {dds output | dps output} [options]\n\n");
  fprintf(stderr,"  dds output    file of dds output in alignment format\n");
  fprintf(stderr,"  dps output    file of dps output in alignment format\n");
  fprintf(stderr,"  Options (default values):\n");
  fprintf(stderr,"    -f  N     specify chain score cutoff (0)\n");
  fprintf(stderr,"    -r <word> specify a word in double quotes\n");
  fprintf(stderr,"              Any database entry whose description contains\n");
  fprintf(stderr,"              the specified word in proper case is removed.\n");
  fprintf(stderr,"              For example, use -r \"ALU\" to remove Alu hits.\n");
  fprintf(stderr,"\n");
  fprintf(stderr,"See also:\n");
  fprintf(stderr,"  dds dps filter gap2 nap show\n");
  fprintf(stderr,"\n");
  fprintf(stderr,"Huang, X.  Fast Comparison of a DNA Sequence with a Protein Sequence Database\n");
  fprintf(stderr,"  Microbial & Comparative Genomics, 1(4): 281-291 (1996).\n");
  fprintf(stderr,"\n");
}


main(argc, argv) int argc; char *argv[];
{ FILE *Ap, *ckopen();
 char  s[ACCLEN+1], *end;
 int   dstart, dend, astart, aend, score, sp;
 char  acc[ACCLEN];	/* accession */
 
 int   dnalen, aalen;	/* query length and max database seq length */
 int   i, j, k;	/* index var */
 char  ch[6];		/* "Chain" */
 int   cutoff;		/* score cutoff */
 char  *line;		/* for a line */
 int   sorted;		/* flag */
 char  *a, *b;		/* temp vars */
 char  *pt;		/* pointer */
 char *ckalloc();	/* space-allocating function  */
 char  *alu;		/* word */
 int   isword;		/* 1, if there is a word specified */
 int   first;		/* index for bucket */
 
 cutoff = 0;
 isword = 0;
 
 program = argv[0];
 if ( argc < 2 )
   { 
     usage();
     exit(1);
   }
 
 if (!strcmp(argv[1], "-help")  ||  !strcmp(argv[1], "-h"))
   {
     usage();
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
     fprintf(stderr,"default chain score cutoff = 0\n");
     exit(1);
   }
 
 if ( argc > 2 )
   { for ( i = 2; i < argc - 1; i += 2 )
     { if ( argv[i][0] != '-' )
       fatal("Each option must begin with a dash (-)");
     switch ( argv[i][1] )
       { case 'f' : (void) sscanf(argv[i+1],"%d", &cutoff);
	   if ( cutoff < 0 )
	     fatal("The cutoff must be a nonnegative integer");
	   break;
       case 'r' : for ( j = 0, a = argv[i+1]; a[j] != '\0'; j++ )
	 ;
	 alu = ( char * ) ckalloc( (j + 1) * sizeof(char));
	 strcpy(alu, argv[i+1]);
	 if ( j )
	   isword = 1;
	 break;
       default  : fatalf("Wrong option letter "); 
	 break;
       }
     }
   }
 strcpy(ch, "Chain");
 Ap = ckopen(argv[1], "r");
 fgets(s, ACCLEN, Ap);
 fgets(s, ACCLEN, Ap);
 sscanf(s, "Query sequence length: %d", &dnalen);
 fgets(s, ACCLEN, Ap);
 sscanf(s, "Maximum database sequence length: %d", &aalen);
 printf("                     %8d %8d\n", dnalen, aalen);
 bucket = ( bbkptr * ) ckalloc( (dnalen + 2) * sizeof(bbkptr) );
 for ( i = 0; i <= dnalen ; i++ )
   bucket[i] = NULL;
 for ( ; (end = fgets(s, ACCLEN, Ap) ) != NULL ; )
   { for ( i = 0; i < 5 && s[i] == ch[i]; i++ )
     ;
   if ( i < 5 )
     continue;
   sscanf(s, "Chain%d %d %d %d %d %d %s",
	  &dstart, &dend, &score, &sp, &astart, &aend, acc);
   for ( j = 5, k = 0; s[j] == ' ' | (s[j] >= '0' && s[j] <= '9'); j++ )
     if ( s[j] == ' ' && s[j-1] >= '0' && s[j-1] <= '9' )
       if ( ++k >= 6 )
	 break;
   for ( ; s[j] == ' ' ; j++ )
     ;
   strcpy(acc, s+j);
   if ( score < cutoff )
     continue;
   line = ( char * ) ckalloc( ACCLEN * sizeof(char));
   if ( dstart <= dend )
     { sprintf(line, "%8d %8d %6d %7d %5d %1d %5d %5d %s",
	       dstart, dend, score, astart, aend, 1, 0, 0, acc);
     first = dstart;
     }
   else
     { sprintf(line, "%8d %8d %6d %7d %5d %1d %5d %5d %s",
	       dend, dstart, score, astart, aend, 0, 0, 0, acc);
     first = dend;
     }
   if ( isword )
     { sorted = 1;
     for ( a = line; *a != '\0' && sorted; a++ )
       { for ( j = 0; a[j] != '\0' && alu[j] != '\0' && a[j] == alu[j]; j++ )
	 ;
       if ( alu[j] == '\0' )
	 sorted = 0; 
       }
     if ( sorted )
       { item = ( bbkptr ) ckalloc( (int) sizeof(bbk) );
       item->line = line;
       item->next = bucket[first];
       bucket[first] = item;
       }
     }
   else
     { item = ( bbkptr ) ckalloc( (int) sizeof(bbk) );
     item->line = line;
     item->next = bucket[first];
     bucket[first] = item;
     }
   }
 for ( i = 0; i <= dnalen; i++ )
   { for ( item = bucket[i]; item != NULL; item = item->next)
     printf("%s", item->line);
   }
 
 exit(0);
}


/* lib.c - library of C procedures. */

/* fatal - print message and die */
fatal(msg)
     char *msg;
{
  fprintf(stderr, "%s\n", msg);
  exit(1);
}

/* fatalf - format message, print it, and die */
fatalf(msg, val)
     char *msg, *val;
{
  fprintf(stderr, msg, val);
  putc('\n', stderr);
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

/* ckalloc - allocate space; check for success */
char *ckalloc(amount)
     long amount;
{
  char *malloc(), *p;
  
  if ((p = malloc( (unsigned) amount)) == NULL)
    fatal("Ran out of memory.");
  return(p);
}
