static char rcsid[] = "$Id: filestring.c 162093 2015-03-26 18:54:22Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "filestring.h"
#include <stdlib.h>
#include <stdarg.h>
#include <ctype.h>		/* For isdigit() */
#include "assert.h"
#include "mem.h"
#include "list.h"


#define BLOCKSIZE 1024

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Simultaneous print to stdout */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif


#define T Filestring_T

struct T {
  int id;
  SAM_split_output_type split_output;

  List_T blocks;
  int nleft;
  char *ptr;

  char *string;
  int strlength;
};


int
Filestring_id (T this) {
  return this->id;
}

void
Filestring_set_split_output (T this, int split_output) {
  this->split_output = split_output;
  return;
}

SAM_split_output_type
Filestring_split_output (T this) {
  return this->split_output;
}


T
Filestring_new (int id) {
  T new = (T) MALLOC_OUT(sizeof(*new));

  new->id = id;
  new->split_output = OUTPUT_NONE;
  new->blocks = (List_T) NULL;
  new->nleft = 0;
  new->ptr = (char *) NULL;

  new->string = (char *) NULL;

  return new;
}

void
Filestring_free (T *old) {
  List_T p;
  char *block;

  if (*old) {
    if ((*old)->string != NULL) {
      FREE_OUT((*old)->string);
    }

    for (p = (*old)->blocks; p != NULL; p = List_next(p)) {
      block = (char *) List_head(p);
      FREE_OUT(block);
    }
    List_free_out(&(*old)->blocks);
    
    FREE_OUT(*old);
  }

  return;
}


void
Filestring_stringify (T this) {
  List_T p, next;
  char *ptr, *dest;
  int nblocks, i;

  if ((nblocks = List_length(this->blocks)) == 0) {
    this->string = (char *) NULL;
    this->strlength = -1;

  } else if (this->string != NULL) {
    /* Already stringified */

  } else {
    this->strlength = (nblocks - 1) * BLOCKSIZE + (BLOCKSIZE - this->nleft);
    dest = this->string = (char *) MALLOC_OUT((this->strlength + 1) * sizeof(char));

    p = this->blocks = List_reverse(this->blocks);

    next = List_next(p);
    while (next != NULL) {
      ptr = (char *) List_head(p);
      for (i = 0; i < BLOCKSIZE; i++) {
	*dest++ = *ptr++;
      }
      p = next;
      next = List_next(p);
    }

    ptr = (char *) List_head(p);
    for (i = 0; i < BLOCKSIZE - this->nleft; i++) {
      *dest++ = *ptr++;
    }

    *dest = '\0';
  }

  return;
}


/* Could assume that Filestring_stringify has been called */
void
Filestring_print (
#ifdef USE_MPI
		  MPI_File fp,
#else
		  FILE *fp,
#endif
		  T this) {
  List_T p, next;
  char *ptr;
  
  if (this == NULL) {
    return;

#ifdef USE_MPI
  } else if (fp == NULL) {
    /* This may not work if worker is from rank 0 */
    Filestring_send(this,/*dest*/0,/*tag*/MPI_TAG_WRITE_STDOUT,MPI_COMM_WORLD);
    
#endif

  } else if (this->string != NULL) {
    /* Already stringified */
#ifdef USE_MPI
    debug1(fwrite(this->string,sizeof(char),this->strlength,stdout));
    MPI_File_write_shared(fp,this->string,this->strlength,MPI_CHAR,MPI_STATUS_IGNORE);
#else
    fwrite(this->string,sizeof(char),this->strlength,fp);
#endif
	    
  } else if (this->blocks == NULL) {
    return;

  } else {
    p = this->blocks = List_reverse(this->blocks);

    next = List_next(p);
    while (next != NULL) {
      ptr = (char *) List_head(p);
#ifdef USE_MPI
      debug1(fwrite(ptr,sizeof(char),BLOCKSIZE,stdout));
      MPI_File_write_shared(fp,ptr,BLOCKSIZE,MPI_CHAR,MPI_STATUS_IGNORE);
#else
      fwrite(ptr,sizeof(char),BLOCKSIZE,fp);
#endif
      p = next;
      next = List_next(p);
    }

    ptr = (char *) List_head(p);
#ifdef USE_MPI
    debug1(fwrite(ptr,sizeof(char),BLOCKSIZE - this->nleft,stdout));
    MPI_File_write_shared(fp,ptr,BLOCKSIZE - this->nleft,MPI_CHAR,MPI_STATUS_IGNORE);
#else
    fwrite(ptr,sizeof(char),BLOCKSIZE - this->nleft,fp);
#endif
  }

  return;
}
      

static void
transfer_char (T this, char c) {
  char *block;

  if (this->nleft == 0) {
    block = (char *) MALLOC_OUT(BLOCKSIZE * sizeof(char));
    this->blocks = List_push_out(this->blocks,(void *) block);
    this->nleft = BLOCKSIZE;
    this->ptr = &(block[0]);
  }
  *this->ptr++ = c;
  this->nleft -= 1;

  return;
}

void
transfer_string (T this, char *string, int bufferlen) {
  char *block, *q;

  for (q = string; --bufferlen >= 0 && *q != '\0'; q++) {
    if (this->nleft == 0) {
      block = (char *) MALLOC_OUT(BLOCKSIZE * sizeof(char));
      this->blocks = List_push_out(this->blocks,(void *) block);
      this->nleft = BLOCKSIZE;
      this->ptr = &(block[0]);
    }
    *this->ptr++ = *q;
    this->nleft -= 1;
  }

  if (bufferlen < 0) {
    fprintf(stderr,"Overflowed buffer without seeing a terminating character\n");
    fprintf(stderr,"String was %s\n",q);
    abort();
  }

  return;
}



#define BUFFERLEN 1024

void
Filestring_put (T this, const char *format, ...) {
  va_list values;

  char BUFFER[BUFFERLEN];
  char *block;
  const char *p;
  char *q, c;
  int precision;

  va_start(values,format);

  p = format;
  debug(printf("format is %s\n",format));
  while (*p != '\0') {
    if ((c = *p) == '\\') {  /* escape */
      debug(printf("Saw an escape character\n"));
      switch (*++p) {
      case 't': transfer_char(this,'\t'); break; /* Actually \t shows up as an ASCII character */
      case '\\': transfer_char(this,'\\'); break;
      default: fprintf(stderr,"Cannot parse \\%c\n",*p);
      }

    } else if (c == '%') {  /* formatting */
      debug(printf("After formatting character saw %c\n",p[1]));
      switch (*++p) {
      case '%':			/* percent sign */
	transfer_char(this,'%');
	break;

      case 'c':			/* character */
	transfer_char(this,(char) va_arg(values, int));
	break;

      case 's': 		/* string */
	for (q = va_arg(values, char *); *q != '\0'; q++) {
	  transfer_char(this,*q);
	}
	break;

      case '.':			/* float or double */
	if (*++p == '*') {
	  precision = va_arg(values, int);
	  ++p;
	} else {
	  sscanf(p,"%d",&precision);
	  while (isdigit(*++p)) ;
	}
	switch (*p) {
	case 'f':
	  sprintf(BUFFER,"%.*f",precision,va_arg(values, double));
	  transfer_string(this,BUFFER,BUFFERLEN);
	  break;

	case 'e':
	  sprintf(BUFFER,"%.*e",precision,va_arg(values, double));
	  transfer_string(this,BUFFER,BUFFERLEN);
	  break;

	case 'g':
	  sprintf(BUFFER,"%.*g",precision,va_arg(values, double));
	  transfer_string(this,BUFFER,BUFFERLEN);
	  break;
	  
	case 's':
	  sprintf(BUFFER,"%.*s",precision,va_arg(values, char *));
	  transfer_string(this,BUFFER,BUFFERLEN);
	  break;

	default: fprintf(stderr,"Cannot parse %%.%d%c\n",precision,*p); abort();
	}
	break;

      case '*':			/* indirect int or string */
	precision = va_arg(values, int);
	debug(printf("format is %c\n",p[1]));
	switch (*++p) {
	case 'd':
	  sprintf(BUFFER,"%*d",precision,va_arg(values, int));
	  transfer_string(this,BUFFER,BUFFERLEN);
	  break;
	case 'u':
	  sprintf(BUFFER,"%*u",precision,va_arg(values, unsigned int));
	  transfer_string(this,BUFFER,BUFFERLEN);
	  break;
	case 's':
	  sprintf(BUFFER,"%*s",precision,va_arg(values, char *));
	  transfer_string(this,BUFFER,BUFFERLEN);
	  break;
	default: fprintf(stderr,"Cannot parse %%*%c\n",*p); abort();
	}
	break;

      case 'd':			/* int */
	sprintf(BUFFER,"%d",va_arg(values, int));
	transfer_string(this,BUFFER,BUFFERLEN);
	break;

      case 'f':			/* float */
	sprintf(BUFFER,"%f",va_arg(values, double));
	transfer_string(this,BUFFER,BUFFERLEN);
	break;

      case 'u':			/* unsigned int */
	sprintf(BUFFER,"%u",va_arg(values, unsigned int));
	transfer_string(this,BUFFER,BUFFERLEN);
	break;
      
      case 'l':
	switch (*++p) {
	case 'd':			/* long int */
	  sprintf(BUFFER,"%ld",va_arg(values, long int));
	  transfer_string(this,BUFFER,BUFFERLEN);
	  break;

	case 'u':			/* unsigned long */
	  sprintf(BUFFER,"%lu",va_arg(values, unsigned long));
	  transfer_string(this,BUFFER,BUFFERLEN);
	  break;

	case 'l':
	  switch (*++p) {
	  case 'd':			/* long long int */
	    sprintf(BUFFER,"%lld",va_arg(values, long long int));
	    transfer_string(this,BUFFER,BUFFERLEN);
	    break;

	  case 'u':			/* unsigned long long */
	    sprintf(BUFFER,"%llu",va_arg(values, unsigned long long));
	    break;

	  default: fprintf(stderr,"Cannot parse %%ll%c\n",*p); abort();
	  }
	  break;

	default: fprintf(stderr,"Cannot parse %%l%c\n",*p); abort();
	}
	break;

      default: fprintf(stderr,"Cannot parse %%%c\n",*p); abort();
      }

    } else {
      /* transfer_char(this,c); -- effectively inlined here */
      if (this->nleft == 0) {
	block = (char *) MALLOC_OUT(BLOCKSIZE * sizeof(char));
	this->blocks = List_push_out(this->blocks,(void *) block);
	this->nleft = BLOCKSIZE;
	this->ptr = &(block[0]);
      }
      *this->ptr++ = c;
      this->nleft -= 1;
    }

    p++;
  }

  va_end(values);

  return;
}

void
Filestring_putc (char c, T this) {
  char *block;

  if (this->nleft == 0) {
    block = (char *) MALLOC_OUT(BLOCKSIZE * sizeof(char));
    this->blocks = List_push_out(this->blocks,(void *) block);
    this->nleft = BLOCKSIZE;
    this->ptr = &(block[0]);
  }
  *this->ptr++ = c;
  this->nleft -= 1;
}


/* Modified from transfer_string */
void
Filestring_puts (T this, char *string, int strlength) {
  char *block, *q;

  for (q = string; --strlength >= 0; q++) {
    if (this->nleft == 0) {
      block = (char *) MALLOC_OUT(BLOCKSIZE * sizeof(char));
      this->blocks = List_push_out(this->blocks,(void *) block);
      this->nleft = BLOCKSIZE;
      this->ptr = &(block[0]);
    }
    *this->ptr++ = *q;
    this->nleft -= 1;
  }

  return;
}



#ifdef USE_MPI
char *
Filestring_extract (int *strlength, T this) {
  Filestring_stringify(this);
  if ((*strlength = this->strlength) == 0) {
    return (char *) NULL;
  } else {
    return this->string;
  }
}


void
Filestring_send (T this, int dest, int tag, MPI_Comm comm) {
  Filestring_stringify(this);
  MPI_SEND(&this->strlength,1,MPI_INT,dest,tag,comm);
  if (this->strlength > 0) {
    MPI_SEND(this->string,this->strlength+1,MPI_CHAR,dest,tag,comm);
  }
  return;
}


char *
Filestring_recv (int *strlength, int source, int tag, MPI_Comm comm) {
  char *string;
  MPI_Status status;

  MPI_RECV(&(*strlength),1,MPI_INT,source,tag,comm,&status);
  if (*strlength <= 0) {
    string = (char *) MALLOC(1 * sizeof(char));
    string[0] = '\0';
    *strlength = 0;
  } else {
    string = (char *) MALLOC(((*strlength) + 1) * sizeof(char));
    MPI_RECV(string,(*strlength) + 1,MPI_CHAR,source,tag,comm,&status);
  }

  return string;
}
#endif


  
