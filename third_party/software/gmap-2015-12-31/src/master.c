static char rcsid[] = "$Id: master.c 162088 2015-03-26 18:29:04Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "master.h"
#include <stdio.h>
#include <stdlib.h>

#ifdef HAVE_PTHREAD
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>		/* Needed to define pthread_t on Solaris */
#endif
#include <pthread.h>
#endif

#include "mem.h"

#ifdef USE_MPI
#include "filestring.h"
#endif
#ifdef GSNAP
#include "shortread.h"
#endif


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


typedef struct RRlist_T *RRlist_T;
struct RRlist_T {
  int nextchar_start;
  int nextchar;
  int offset_start_1;
  int offset_start_2;
  int offset_end_1;
  int offset_end_2;
  Filestring_T filestring1;
  Filestring_T filestring2;
  bool donep;

  RRlist_T next;
};


#ifdef DEBUG
static void
RRlist_dump (RRlist_T head, RRlist_T tail) {
  RRlist_T this;

  fprintf(stdout,"head %p\n",head);
  for (this = head; this != NULL; this = this->next) {
    fprintf(stdout,"%p: offsets %d..%d and %d..%d, next %p\n",
	    this,this->offset_start_1,this->offset_end_1,this->offset_start_2,this->offset_end_2,this->next);
  }
  fprintf(stdout,"tail %p\n",tail);
  fprintf(stdout,"\n");
  return;
}
#endif


/* Returns new tail */
static RRlist_T
RRlist_push (RRlist_T *head, RRlist_T tail, int nextchar_start, int nextchar, int offset_start_1, int offset_start_2,
	     int offset_end_1, int offset_end_2, Filestring_T filestring1, Filestring_T filestring2, bool donep) {
  RRlist_T new;

  new = (RRlist_T) MALLOC_OUT(sizeof(*new));
  new->nextchar_start = nextchar_start;
  new->nextchar = nextchar;
  new->offset_start_1 = offset_start_1;
  new->offset_start_2 = offset_start_2;
  new->offset_end_1 = offset_end_1;
  new->offset_end_2 = offset_end_2;
  new->filestring1 = filestring1;
  new->filestring2 = filestring2;
  new->donep = donep;

  new->next = (RRlist_T) NULL;
  
  if (*head == NULL) {		/* Equivalent to tail == NULL, but using *head avoids having to set tail in RRlist_pop */
    *head = new;
  } else {
    tail->next = new;
  }

  return new;
}

/* Returns new head */
static RRlist_T
RRlist_pop (RRlist_T head, int *nextchar_start, int *nextchar, int *offset_start_1, int *offset_start_2,
	    int *offset_end_1, int *offset_end_2, Filestring_T *filestring1, Filestring_T *filestring2, bool *donep) {
  RRlist_T newhead;

  *nextchar_start = head->nextchar_start;
  *nextchar = head->nextchar;
  *offset_start_1 = head->offset_start_1;
  *offset_start_2 = head->offset_start_2;
  *offset_end_1 = head->offset_end_1;
  *offset_end_2 = head->offset_end_2;
  *filestring1 = head->filestring1;
  *filestring2 = head->filestring2;
  *donep = head->donep;

  newhead = head->next;

  FREE_OUT(head);
  return newhead;
}



#define T Master_T
struct T {
#ifdef HAVE_PTHREAD
  pthread_mutex_t lock;
  pthread_cond_t input_wanted_p;
  pthread_cond_t input_avail_p;
#endif

  /* Predicate for input_wanted_p */
  int nwanted;

  /* Predicate for input_avail_p */
  RRlist_T head;
  RRlist_T tail;
  
  int ntotal;

  int n_slave_ranks;
  int nextchar;
  int nchars1;
  int nchars2;

  FILE *input_parser;
  FILE *input2_parser;
#ifdef HAVE_ZLIB
  gzFile gzipped;
  gzFile gzipped2;
#endif
#ifdef HAVE_BZLIB
  Bzip2_T bzipped;
  Bzip2_T bzipped2;
#endif

  char **files;
  int nfiles;
  int nspaces;
  int part_modulus;
  int part_interval;
};


int
Master_ntotal (T this) {
  return this->ntotal;
}


/* Modeled after fill_buffer in inbuffer.c */
T
Master_new (int n_slave_ranks, int nextchar, int nchars1, int nchars2,
	    FILE *input_parser, FILE *input2_parser,
#ifdef HAVE_ZLIB
	    gzFile gzipped, gzFile gzipped2,
#endif
#ifdef HAVE_BZLIB
	    Bzip2_T bzipped, Bzip2_T bzipped2,
#endif
	    char **files, int nfiles, int nspaces, int part_modulus, int part_interval) {
  T new = (T) MALLOC(sizeof(*new));

#ifdef HAVE_PTHREAD
  pthread_mutex_init(&new->lock,NULL);
  pthread_cond_init(&new->input_wanted_p,NULL);
  pthread_cond_init(&new->input_avail_p,NULL);
#endif

  new->nwanted = 0;
  new->head = (RRlist_T) NULL;
  new->tail = (RRlist_T) NULL;

  new->ntotal = 0;

  new->n_slave_ranks = n_slave_ranks;
  new->nextchar = nextchar;
  new->nchars1 = nchars1;
  new->nchars2 = nchars2;

  new->input_parser = input_parser;
  new->input2_parser = input2_parser;

#ifdef HAVE_ZLIB
  new->gzipped = gzipped;
  new->gzipped2 = gzipped2;
#endif

#ifdef HAVE_BZLIB
  new->bzipped = bzipped;
  new->bzipped2 = bzipped2;
#endif

  new->files = files;
  new->nfiles = nfiles;
  new->nspaces = nspaces;
  new->part_modulus = part_modulus;
  new->part_interval = part_interval;

  return new;
}

void
Master_free (T *old) {
  pthread_cond_destroy(&(*old)->input_wanted_p);
  pthread_cond_destroy(&(*old)->input_avail_p);
  pthread_mutex_destroy(&(*old)->lock);

  FREE(*old);
  return;
}



/* Run only by rank 0, if output is designed to go to stdout */
void *
Master_write_stdout (void *data) {
  int strlength;
  char *string;

  while (true) {
    /* This may not work if worker is from rank 0 */
    string = Filestring_recv(&strlength,/*source*/MPI_ANY_SOURCE,/*tag*/MPI_TAG_WRITE_STDOUT,MPI_COMM_WORLD);
    fwrite(string,sizeof(char),strlength,stderr);
    fwrite(string,sizeof(char),strlength,stdout);
  }

  return (void *) NULL;
}


/* Run only by rank 0 */
/* Communicates below with Master_mpi_interface for ranks 1..n and with Master_self_interface for rank 0 */
void *
Master_parser (void *data) {
  T this = (T) data;

  int nextchar = this->nextchar;

  int nspaces = this->nspaces;
  int part_modulus = this->part_modulus;
  int part_interval = this->part_interval;

  bool donep = false;
  int nskip;
  Shortread_T queryseq1, queryseq2;
  Filestring_T filestring1, filestring2;
  int nextchar_start;
  int offset_start_1, offset_end_1, offset_start_2, offset_end_2;
  
#if 0
  /* For some reason, this doesn't work */
  offset_start_1 = nchars1;
  offset_start_2 = nchars2;
#else
  offset_start_1 = 0;
  offset_start_2 = 0;
#endif

  filestring1 = Filestring_new(/*id*/0);
  filestring2 = Filestring_new(/*id*/0);

  /* Skip to part_modulus */
  if (part_modulus > 0) {
    nskip = 0;
    while (nskip < part_modulus &&
	   (queryseq1 = Shortread_read(&nextchar,&offset_start_1,&offset_start_2,&queryseq2,
#ifdef USE_MPI
				       filestring1,filestring2,
#endif
				       &this->input_parser,&this->input2_parser,
#ifdef HAVE_ZLIB
				       &this->gzipped,&this->gzipped2,
#endif
#ifdef HAVE_BZLIB
				       &this->bzipped,&this->bzipped2,
#endif
				       &this->files,&this->nfiles,/*skipp*/true)) != NULL) {
      nskip++;
    }
    if (queryseq1 == NULL) {
      donep = true;
    }
  }

  Filestring_free(&filestring2);
  Filestring_free(&filestring1);



  nextchar_start = nextchar;
  offset_end_1 = offset_start_1;
  offset_end_2 = offset_start_2;

  while (true) {
    pthread_mutex_lock(&this->lock);
    debug(fprintf(stdout,"nwanted is %d\n",this->nwanted));
    while (this->nwanted == 0) {
      pthread_cond_wait(&this->input_wanted_p,&this->lock);
    }

    filestring1 = Filestring_new(/*id*/0);
    filestring2 = Filestring_new(/*id*/0);

    nskip = 0;
    while (nskip < nspaces * part_interval &&
	   (queryseq1 = Shortread_read(&nextchar,&offset_end_1,&offset_end_2,&queryseq2,
#ifdef USE_MPI
				       filestring1,filestring2,
#endif
				       &this->input_parser,&this->input2_parser,
#ifdef HAVE_ZLIB
				       &this->gzipped,&this->gzipped2,
#endif
#ifdef HAVE_BZLIB
				       &this->bzipped,&this->bzipped2,
#endif
				       &this->files,&this->nfiles,/*skipp*/true)) != NULL) {
      nskip++;
    }
    this->ntotal += (nskip + part_interval - 1)/part_interval;
    
    if (queryseq1 == NULL) {
      donep = true;
    }

    this->tail = RRlist_push(&this->head,this->tail,nextchar_start,nextchar,offset_start_1,offset_start_2,
			     offset_end_1,offset_end_2,filestring1,filestring2,donep);
    debug(RRlist_dump(this->head,this->tail));

    nextchar_start = nextchar;
    offset_start_1 = offset_end_1;
    offset_start_2 = offset_end_2;

    this->nwanted -= 1;
    pthread_cond_signal(&this->input_avail_p);
    pthread_mutex_unlock(&this->lock);
  }

  return (void *) NULL;
}


void
Master_self_interface (T this, int *nextchar_start, int *nextchar,
		       int *offset_start_1, int *offset_start_2,
		       int *offset_end_1, int *offset_end_2,
		       Filestring_T *filestring1, Filestring_T *filestring2,
		       bool *donep) {

  /* Get information from Master_parser */
  debug(printf("Master_self_interface called and locking master\n"));
  pthread_mutex_lock(&this->lock);
  this->nwanted += 1;
  pthread_cond_signal(&this->input_wanted_p);
  pthread_mutex_unlock(&this->lock);

  pthread_mutex_lock(&this->lock);
  while (this->head == NULL) {
    pthread_cond_wait(&this->input_avail_p,&this->lock);
  }
  this->head = RRlist_pop(this->head,&(*nextchar_start),&(*nextchar),&(*offset_start_1),&(*offset_start_2),&(*offset_end_1),&(*offset_end_2),
			  &(*filestring1),&(*filestring2),&(*donep));
  debug(RRlist_dump(this->head,this->tail));

  debug(printf("Master_self_interface unlocking master\n"));
  pthread_mutex_unlock(&this->lock);

  debug(fprintf(stdout,"Master_self_interface now returning\n"));
  return;
}



#ifdef USE_MPI
/* Run only by rank 0 */
/* Communicates below with fill_buffer procedure of ranks 1..n, and above with Master_parser */
/* Replaces Inbuffer_new for MPI master */
void *
Master_mpi_interface (void *data) {
  T this = (T) data;

  int n_slave_ranks = this->n_slave_ranks;
  int nfiles_slave;
  MPI_Status status;
  int ranki;
  int nextchar_start, nextchar = this->nextchar;
  int offset_start_1, offset_start_2, offset_end_1, offset_end_2;
  Filestring_T filestring1, filestring2;
  bool donep;

#ifdef HAVE_ZLIB
  gzFile gzipped = this->gzipped;
#endif
#ifdef HAVE_BZLIB
  Bzip2_T bzipped = this->bzipped;
#endif


  while (n_slave_ranks > 0) {
  /* Need to send nextchar_end (nextchar at end of block)
     because of the difference between filecontents end ('\0') and
     FILE * end (EOF) */

    MPI_RECV(&nfiles_slave,1,MPI_INT,/*source*/MPI_ANY_SOURCE,/*tag*/MPI_TAG_WANT_INPUT,MPI_COMM_WORLD,&status);
    ranki = status.MPI_SOURCE;

    /* Get information from Master_parser */
    debug(printf("Master_mpi_interface locking master\n"));
    pthread_mutex_lock(&this->lock);
    this->nwanted += 1;
    pthread_cond_signal(&this->input_wanted_p);
    pthread_mutex_unlock(&this->lock);

    pthread_mutex_lock(&this->lock);
    while (this->head == NULL) {
      pthread_cond_wait(&this->input_avail_p,&this->lock);
    }
    this->head = RRlist_pop(this->head,&nextchar_start,&nextchar,&offset_start_1,&offset_start_2,&offset_end_1,&offset_end_2,
			    &filestring1,&filestring2,&donep);
    debug(RRlist_dump(this->head,this->tail));
    debug(printf("Master_mpi_interface unlocking master\n"));
    pthread_mutex_unlock(&this->lock);

    debug(fprintf(stdout,"Master: received message from %d.  Sending nextchars %c..%c, donep %d",
		  ranki,nextchar_start,nextchar,donep));
    MPI_SEND(&nextchar_start,1,MPI_INT,/*dest*/ranki,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD);
    MPI_SEND(&nextchar,1,MPI_INT,/*dest*/ranki,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD);
    MPI_SEND(&donep,1,MPI_UNSIGNED_CHAR,/*dest*/ranki,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD);

#if defined(HAVE_ZLIB) && defined(HAVE_BZLIB)
    if (gzipped == NULL && bzipped == NULL) {
      debug(fprintf(stdout,"  Sending offsets %d..%d and %d..%d.\n",offset_start_1,offset_end_1,offset_start_2,offset_end_2));
      MPI_SEND(&offset_start_1,1,MPI_INT,/*dest*/ranki,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD);
      MPI_SEND(&offset_start_2,1,MPI_INT,/*dest*/ranki,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD);
      MPI_SEND(&offset_end_1,1,MPI_INT,/*dest*/ranki,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD);
      MPI_SEND(&offset_end_2,1,MPI_INT,/*dest*/ranki,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD);
    } else {
      debug(fprintf(stdout,"  Sending filestrings\n"));
      Filestring_send(filestring1,/*dest*/ranki,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD);
      Filestring_send(filestring2,/*dest*/ranki,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD);
    }

#elif defined(HAVE_ZLIB)
    if (gzipped == NULL) {
      debug(fprintf(stdout,"  Sending offsets %d..%d and %d..%d.\n",offset_start_1,offset_end_1,offset_start_2,offset_end_2));
      debug(fprintf(stdout,"  Sending end offsets %d and %d.\n",offset_end_1,offset_end_2));
      MPI_SEND(&offset_start_1,1,MPI_INT,/*dest*/ranki,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD);
      MPI_SEND(&offset_start_2,1,MPI_INT,/*dest*/ranki,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD);
      MPI_SEND(&offset_end_1,1,MPI_INT,/*dest*/ranki,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD);
      MPI_SEND(&offset_end_2,1,MPI_INT,/*dest*/ranki,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD);
    } else {
      debug(fprintf(stdout,"  Sending filestrings\n"));
      Filestring_send(filestring1,/*dest*/ranki,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD);
      Filestring_send(filestring2,/*dest*/ranki,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD);
    }

#elif defined(HAVE_BZLIB)
    if (bzipped == NULL) {
      debug(fprintf(stdout,"  Sending offsets %d..%d and %d..%d.\n",offset_start_1,offset_end_1,offset_start_2,offset_end_2));
      MPI_SEND(&offset_start_1,1,MPI_INT,/*dest*/ranki,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD);
      MPI_SEND(&offset_start_2,1,MPI_INT,/*dest*/ranki,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD);
      MPI_SEND(&offset_end_1,1,MPI_INT,/*dest*/ranki,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD);
      MPI_SEND(&offset_end_2,1,MPI_INT,/*dest*/ranki,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD);
    } else {
      debug(fprintf(stdout,"  Sending filestrings\n"));
      Filestring_send(filestring1,/*dest*/ranki,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD);
      Filestring_send(filestring2,/*dest*/ranki,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD);
    }

#else
    debug(fprintf(stdout,"  Sending offsets %d..%d and %d..%d.\n",offset_start_1,offset_end_1,offset_start_2,offset_end_2));
    MPI_SEND(&offset_start_1,1,MPI_INT,/*dest*/ranki,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD);
    MPI_SEND(&offset_start_2,1,MPI_INT,/*dest*/ranki,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD);
    MPI_SEND(&offset_end_1,1,MPI_INT,/*dest*/ranki,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD);
    MPI_SEND(&offset_end_2,1,MPI_INT,/*dest*/ranki,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD);
#endif

    Filestring_free(&filestring2);
    Filestring_free(&filestring1);

    debug(fprintf(stdout,"\n"));

    if (donep == true) {
      n_slave_ranks -= 1;
      debug(fprintf(stdout,"n_slave_ranks is now %d.  Rank %d knows we are done.\n",n_slave_ranks,ranki));
    }
  }

  debug(fprintf(stdout,"Master_mpi_interface now returning\n"));
  return (void *) NULL;
}

#endif

