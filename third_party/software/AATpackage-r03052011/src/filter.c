
/* This filter program takes as input a file of coordinates produced by EXT.
   It lets at most N largest-scoring chains from each overlapping region
   go through.
   
   Acknowledgments
   I thank Xiaoying Lin for suggestions that improved this program.
   
   Usage: filter  Chain_file  [options]  >  out.filter
   
   Chain_file   file produced by EXT
   
   Options (default values):
   
   -c  N  specify max number N of chains from an overlapping region (3)
*/

#define    ACCLEN  400	/* max length of a line */
#define    DEPTH     5	/* default value for k */

#include   <stdio.h>

// Version info
static float VersionS = 1.51f;
static char * BuildS = "$Revision: 1.7 $";
static char * program = "";         // program name 

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
{ FILE *Ap, *ckopen();
 char  s[ACCLEN+1], *end;
 int   dstart, dend, astart, aend, score, ort;
 char  acc[ACCLEN];
 // int   strcpy();
 int   fcs, fce;
 int   rcs, rce;
 int   x, y;
 
 int   sorted;		/* flag */
 char  *pt;		/* pointer */
 char *ckalloc();			/* space-allocating function  */
 int   k;	/* max number of chains taken from an overlapping region */
 char  **outlist;	/* output list for this program */
 char  **flist;	/* chains taken from a plus strand region */
 int   *fscore;	/* their scores */
 int   *fend;		/* right end position of a chain on plus strand */
 char  **rlist;	/* chains taken from a minus strand region */
 int   *rscore;	/* their scores */
 int   *rend;		/* right end position of a chain on minus strand */
 int   insize;		/* size of input chains */
 int   outsize;	/* size of outlist */
 int   fsize;		/* size of flist */
 int   rsize;		/* size of rlist */
 char  *a, *b;		/* temp vars */
 char  *cc;		/* new chain */
 int   temp;		/* temp var */
 int   i, j;		/* index vars */
 int   fmin;		/* min score of chains saved in flist */
 int   fpos;		/* position of a chain with score fmin in flist */
 int   rmin;		/* min score of chains saved in rlist */
 int   rpos;		/* position of a chain with score rmin in rlist */
 
 k = DEPTH;
 program = argv[0];
 
 if ( argc < 2 )
   { fprintf(stderr,"Usage: %s  Chain_file  [options]\n\n", argv[0]);
   fprintf(stderr,"  Chain_file   file produced by EXT.\n");
   fprintf(stderr,"Options (default values):\n");
   fprintf(stderr,"  -c  N  specify max number N of chains from an overlapping region (5)\n");
   exit(1);
   }
 
 if (!strcmp(argv[1], "-help")  ||  !strcmp(argv[1], "-h"))
   {
     fprintf(stderr,"Usage: %s  Chain_file  [options]\n\n", argv[0]);
     fprintf(stderr,"  Chain_file   file produced by EXT.\n");
     fprintf(stderr,"Options (default values):\n");
     fprintf(stderr,"  -c  N  specify max number N of chains from an overlapping region (5)\n");
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
     fprintf(stderr,"default maximum number of chains from an overlapping region = 5\n");
     exit(1);
   }
 
 if ( argc > 2 )
   { for ( i = 2; i < argc - 1; i += 2 )
     { if ( argv[i][0] != '-' )
       fatal("Each option must begin with a dash (-)");
     (void) sscanf(argv[i+1],"%d", &temp);
     switch ( argv[i][1] )
       { case 'c' : if ( temp < 1 )
	 fatal("Value for max number of chains must be a positive integer");
	   k = temp;
	   break;
       default  : fatalf("Wrong option letter "); 
	 break;
       }
     }
   }
 Ap = ckopen(argv[1], "r");
 fgets(s, ACCLEN, Ap);
 for ( insize = 1 ; (end = fgets(s, ACCLEN, Ap) ) != NULL; )
   insize++;
 outlist = ( char ** ) ckalloc( insize * sizeof(char *));
 flist = ( char ** ) ckalloc( k * sizeof(char *));
 rlist = ( char ** ) ckalloc( k * sizeof(char *));
 fscore = ( int * ) ckalloc( k * sizeof(int));
 rscore = ( int * ) ckalloc( k * sizeof(int));
 fend = ( int * ) ckalloc( k * sizeof(int));
 rend = ( int * ) ckalloc( k * sizeof(int));
 fclose(Ap);
 
 fce = rce = fmin = rmin = outsize = fsize = rsize = 0;
 Ap = ckopen(argv[1], "r");
 fgets(s, ACCLEN, Ap);
 printf("%s", s);
 for ( ; (end = fgets(s, ACCLEN, Ap) ) != NULL; )
   { sscanf(s, "%d %d %d %d %d %d %d %d %s",
	    &dstart, &dend, &score, &astart, &aend, &ort, &x, &y, acc);
   if ( ort )
     { if ( dstart <= fce )
       { if ( fmin < score )
	 { if ( fsize < k )
	   { cc = ( char * ) ckalloc( ACCLEN * sizeof(char));
	   strcpy(cc, s);
	   flist[fsize] = cc; 
	   fend[fsize] = dend; 
	   fscore[fsize++] = score; 
	   }
	 else
	   { strcpy(flist[fpos], s);
	   fend[fpos] = dend; 
	   fscore[fpos] = score; 
	   }
	 if ( fce < dend )
	   fce = dend;
	 if ( fsize == k )
	   { fmin = fscore[fpos = 0];
	   fce = fend[0];
	   for ( i = 1; i < k; i++ )
	     { if ( fscore[i] < fmin )
	       fmin = fscore[fpos = i];
	     if ( fce < fend[i] )
	       fce = fend[i];
	     }
	   }
	 }
       }
     else
       { for ( i = 0; i < fsize; i++ )
	 outlist[outsize++] = flist[i];
       fcs = dstart;
       fce = dend;
       cc = ( char * ) ckalloc( ACCLEN * sizeof(char));
       strcpy(cc, s);
       flist[0] = cc; 
       fscore[0] = score; 
       fend[0] = dend; 
       fsize = 1;
       if ( k == 1 )
	 { fmin = score;
	 fpos = 0;
	 }
       else
	 fmin = 0;
       }
     }
   else
     { if ( dstart <= rce )
       { if ( rmin < score )
	 { if ( rsize < k )
	   { cc = ( char * ) ckalloc( ACCLEN * sizeof(char));
	   strcpy(cc, s);
	   rlist[rsize] = cc; 
	   rend[rsize] = dend; 
	   rscore[rsize++] = score; 
	   }
	 else
	   { strcpy(rlist[rpos], s);
	   rend[rpos] = dend; 
	   rscore[rpos] = score; 
	   }
	 if ( rce < dend )
	   rce = dend;
	 if ( rsize == k )
	   { rmin = rscore[rpos = 0];
	   rce = rend[0];
	   for ( i = 1; i < k; i++ )
	     { if ( rscore[i] < rmin )
	       rmin = rscore[rpos = i];
	     if ( rce < rend[i] )
	       rce = rend[i];
	     }
	   }
	 }
       }
     else
       { for ( i = 0; i < rsize; i++ )
	 outlist[outsize++] = rlist[i];
       rcs = dstart;
       rce = dend;
       cc = ( char * ) ckalloc( ACCLEN * sizeof(char));
       strcpy(cc, s);
       rlist[0] = cc; 
       rscore[0] = score; 
       rend[0] = dend; 
       rsize = 1;
       if ( k == 1 )
	 { rmin = score;
	 rpos = 0;
	 }
       else
	 rmin = 0;
       }
     }
   }
 
 if ( fsize )
   for ( i = 0; i < fsize; i++ )
     outlist[outsize++] = flist[i];
 if ( rsize )
   for ( i = 0; i < rsize; i++ )
     outlist[outsize++] = rlist[i];
 
 for ( i = 0; i < outsize - 1; i++ )
   { sorted = 1;
   for ( j = outsize - 1; j > i; j--)
     { a = outlist[j-1];
     b = outlist[j];
     for ( ; *a != '\n' && *b != '\n' && *a == *b; a++, b++ )
       ;
     if ( *b == '\n' || *a != '\n' && *b < *a )
       { pt = outlist[j];
       outlist[j] = outlist[j-1];
       outlist[j-1] = pt;
       sorted = 0;
       }
     }
   if ( sorted )
     break;
   }
 for ( i = 0; i < outsize; i++ )
   printf("%s", outlist[i]);
 
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
