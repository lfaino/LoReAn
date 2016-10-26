static char rcsid[] = "$Id: reader.c 57095 2012-02-03 19:43:00Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "reader.h"
#include "mem.h"

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


#define T Reader_T
struct T {
  int querystart;
  int queryend;
  int querystart_save;
  int queryend_save;

  char *startinit;
  char *startptr;
  char *endptr;

  char *startbound;		/* Saved for reset */
  char *endbound;
};

int
Reader_querystart (T this) {
  return this->querystart;
}

int
Reader_queryend (T this) {
  return this->queryend;
}

int
Reader_startpos (T this) {
  return (this->startptr - this->startinit);
}

int
Reader_endpos (T this) {
  return (this->endptr - this->startinit);
}


void
Reader_reset_start (T this, int querypos) {
  char *sequence;

  sequence = this->startinit;
  this->startptr = &(sequence[querypos]);
  return;
}

void
Reader_reset_end (T this, int querypos) {
  char *sequence;

  sequence = this->startinit;
  this->endptr = &(sequence[querypos]);
  return;
}


void
Reader_reset_ends (T this) {
  /*
  char *sequence;
  sequence = this->startinit;
  this->startptr = &(sequence[this->querystart]);
  this->endptr = &(sequence[this->queryend-1]);
  */

  this->querystart = this->querystart_save;
  this->queryend = this->queryend_save;

  this->startptr = this->startbound;
  this->endptr = this->endbound;

  return;
}



T
Reader_new (char *sequence, int querystart, int queryend) {
  T new = (T) MALLOC(sizeof(*new));

  new->querystart_save = new->querystart = querystart;
  new->queryend_save = new->queryend = queryend;

  new->startinit = sequence;
  new->startbound = new->startptr = &(sequence[querystart]);
  new->endbound = new->endptr = &(sequence[queryend-1]);
  return new;
}

void
Reader_free (T *old) {
  if (*old) {
    FREE(*old);
  }
  return;
}

#if 0
char
Reader_getc_old (T this, cDNAEnd_T cdnaend) {
  debug(fprintf(stderr,"Read_getc has startptr %d and endptr %d\n",
		this->startptr-this->startinit,this->endptr-this->startinit));
  if (this->startptr - this->endptr >= this->blocksize) {
    return '\0';
  } else if (cdnaend == FIVE) {
    return *(this->startptr++);
  } else { 
    return *(this->endptr--);
  }
}
#endif

char
Reader_getc (T this, cDNAEnd_T cdnaend) {
  debug(printf("Read_getc has startptr %d and endptr %d\n",
	       this->startptr - this->startinit,this->endptr - this->startinit));
  if (cdnaend == FIVE) {
    if (this->startptr > this->endbound) {
      return '\0';
    } else {
      return *(this->startptr++);
    }
  } else {
    if (this->endptr < this->startbound) {
      return '\0';
    } else {
      return *(this->endptr--);
    }
  }
}


/* For testing */
/*
static void
process_input (FILE *input) {
  bool ssfile;
  int queryno = 0, seqlength, skiplength;
  char *query, initc;

  if ((initc = Reader_input_init(input)) == '>') {
    ssfile = false;
    query = Reader_input_header(input);
  } else {
    ssfile = true;
    query = (char *) CALLOC(strlen("NO_HEADER")+1,sizeof(char));
    strcpy(query,"NO_HEADER");
  }
  fprintf(stderr,"Read sequence %s\n",query);
  Reader_input_sequence(&seqlength,&skiplength,input,initc);
  print_contents(&(Sequence[0]),SEQUENCELEN);
  fprintf(stderr,"\nSeqlength = %d\n",seqlength);
  FREE(query);
  queryno++;

  while ((query = Reader_input_header(input)) != NULL) {
    fprintf(stderr,"Read sequence %s\n",query);
    Reader_input_sequence(&seqlength,&skiplength,input,initc);
    print_contents(&(Sequence[0]),SEQUENCELEN);
    fprintf(stderr,"\nSeqlength = %d\n",seqlength);
    FREE(query);
    queryno++;
  }

  return;
}

int
main (int argc, char *argv[]) {
  process_input(stdin);
  return 0;
}
*/
