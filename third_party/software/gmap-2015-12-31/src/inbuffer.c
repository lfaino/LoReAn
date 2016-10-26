static char rcsid[] = "$Id: inbuffer.c 175728 2015-09-30 15:08:16Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "inbuffer.h"
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


static bool filter_if_both_p;

#if defined(USE_MPI) && defined(USE_MPI_FILE_INPUT)
static MPI_Comm workers_comm;
#endif

#ifndef GSNAP
static bool user_pairalign_p;
static Sequence_T global_usersegment;
#endif

static int part_modulus;
static int part_interval;

void
Inbuffer_setup (bool filter_if_both_p_in, 
#if defined(USE_MPI) && defined(USE_MPI_FILE_INPUT)
		MPI_Comm workers_comm_in,
#endif
#ifndef GSNAP
		bool user_pairalign_p_in, Sequence_T global_usersegment_in,
#endif

		int part_modulus_in, int part_interval_in) {
  filter_if_both_p = filter_if_both_p_in;

#if defined(USE_MPI) && defined(USE_MPI_FILE_INPUT)
  workers_comm = workers_comm_in;
#endif

#ifndef GSNAP
  user_pairalign_p = user_pairalign_p_in;
  global_usersegment = global_usersegment_in;
#endif

  part_modulus = part_modulus_in;
  part_interval = part_interval_in;

  return;
}



#define T Inbuffer_T

struct T {
#ifdef USE_MPI
  Master_T master;
#endif
  Outbuffer_T outbuffer;

#if defined(USE_MPI) && defined(USE_MPI_FILE_INPUT)
  MPI_File input;
#ifdef GSNAP
  MPI_File input2;
#endif

#elif (defined(USE_MPI))
  FILE *input;
#ifdef GSNAP
  FILE *input2;
#endif

#else
  FILE *input;
#ifdef GSNAP
  FILE *input2;
#endif
#endif

#ifdef USE_MPI
  int myid;
  char *filecontents1_alloc;
  char *filecontents1;
  char *filecontents2_alloc;
  char *filecontents2;
#endif

#ifdef HAVE_ZLIB
  gzFile gzipped;
  gzFile gzipped2;
#else
  void *gzipped;
  void *gzipped2;
#endif

#ifdef HAVE_BZLIB
  Bzip2_T bzipped;
  Bzip2_T bzipped2;
#else
  void *bzipped;
  void *bzipped2;
#endif

  char **files;
  int nfiles;
  int nextchar;

#if defined(HAVE_PTHREAD)
  pthread_mutex_t lock;
#endif

#ifndef GSNAP
  Sequence_T pairalign_segment;
#endif
  Request_T *buffer;
  unsigned int nspaces;
  int ptr;
  int nleft;
  int inputid;
  int requestid;
};


#ifndef GSNAP
T
Inbuffer_cmdline (char *contents, int length) {
  T new = (T) MALLOC(sizeof(*new));

#if defined(USE_MPI) && defined(USE_MPI_FILE_INPUT)
  new->input = (MPI_File) NULL;
#else
  new->input = (FILE *) NULL;
#endif

#ifdef USE_MPI
  new->filecontents1_alloc = (char *) NULL;
  new->filecontents2_alloc = (char *) NULL;
#endif

  new->files = (char **) NULL;
  new->nfiles = 0;
  new->nextchar = '\0';

  new->pairalign_segment = (Sequence_T) NULL;
  new->buffer = (Request_T *) CALLOC(1,sizeof(Request_T));

  new->ptr = 0;
  new->nleft = 1;
  new->inputid = 0;
  new->requestid = 0;

  new->buffer[0] = Request_new(new->requestid++,Sequence_genomic_new(contents,length,/*copyp*/true));

#if defined(HAVE_PTHREAD)
  pthread_mutex_init(&new->lock,NULL);
#endif

  return new;
}
#endif


T
Inbuffer_new (int nextchar,
#ifdef USE_MPI
	      int myid,
#endif
#if defined(USE_MPI) && defined(USE_MPI_FILE_INPUT)
	      MPI_File input,
#else
	      FILE *input,
#endif
#ifdef GSNAP
#if defined(USE_MPI) && defined(USE_MPI_FILE_INPUT)
	      MPI_File input2,
#else
	      FILE *input2,
#endif
#ifdef HAVE_ZLIB
	      gzFile gzipped, gzFile gzipped2,
#endif
#ifdef HAVE_BZLIB
	      Bzip2_T bzipped, Bzip2_T bzipped2,
#endif
#endif
	      char **files, int nfiles, unsigned int nspaces) {

  T new = (T) MALLOC(sizeof(*new));

#ifdef USE_MPI
  new->myid = myid;
#endif

  new->input = input;
#ifdef GSNAP
  new->input2 = input2;
#ifdef HAVE_ZLIB
  new->gzipped = gzipped;
  new->gzipped2 = gzipped2;
#else
  new->gzipped = (void *) NULL;
  new->gzipped2 = (void *) NULL;
#endif
#ifdef HAVE_BZLIB
  new->bzipped = bzipped;
  new->bzipped2 = bzipped2;
#else
  new->bzipped = (void *) NULL;
  new->bzipped2 = (void *) NULL;
#endif
#endif

#ifdef USE_MPI
  new->filecontents1_alloc = (char *) NULL;
  new->filecontents2_alloc = (char *) NULL;
#endif

  new->files = files;
  new->nfiles = nfiles;
  new->nextchar = nextchar;

#if defined(HAVE_PTHREAD)
  pthread_mutex_init(&new->lock,NULL);
#endif

#ifndef GSNAP
  new->pairalign_segment = (Sequence_T) NULL;
#endif
  new->buffer = (Request_T *) CALLOC(nspaces,sizeof(Request_T));
  new->nspaces = nspaces;
  new->ptr = 0;
  new->nleft = 0;
  new->inputid = 0;
  new->requestid = 0;

  return new;
}

#ifdef USE_MPI
void
Inbuffer_set_master (T this, Master_T master) {
  this->master = master;
  return;
}
#endif

void
Inbuffer_set_outbuffer (T this, Outbuffer_T outbuffer) {
  this->outbuffer = outbuffer;
  return;
}

void
Inbuffer_free (T *old) {
  if (*old) {
    /* No need to close input, since done by Shortread and Sequence read procedures */

#ifdef USE_MPI
    FREE_IN((*old)->filecontents1_alloc);
    FREE_IN((*old)->filecontents2_alloc);
#endif

    FREE((*old)->buffer);
    
#if defined(HAVE_PTHREAD)
    pthread_mutex_destroy(&(*old)->lock);
#endif

    FREE(*old);
  }
  return;
}


#ifndef GSNAP
/* Can delete when we remove worker_mpi_process from gmap.c */
Sequence_T
Inbuffer_read (Sequence_T *pairalign_segment, T this, bool skipp) {
  Sequence_T queryseq;

  queryseq = Sequence_read_multifile(&this->nextchar,&this->input,&this->files,&this->nfiles);
  if (skipp == true) {
    Sequence_free(&queryseq);
  }

  if (user_pairalign_p == true) {
    /* assert(this->nspaces == 1) */
    if (this->pairalign_segment != NULL) {
      Sequence_free(&this->pairalign_segment);
    }
    this->pairalign_segment = Sequence_read_unlimited(&this->nextchar,stdin);
    debug(printf("  but first reading usersegment, got nextchar %c\n",this->nextchar));
  }

  this->inputid++;

  *pairalign_segment = this->pairalign_segment;
  return queryseq;
}

#endif


#ifdef USE_MPI
/* Used by rank 0 to communicate with Master_parser thread of rank 0 */
/* Returns number of requests read */
static unsigned int
fill_buffer_master (T this) {
  unsigned int nread = 0;
  Shortread_T queryseq1, queryseq2;
  Filestring_T filestring1, filestring2;
  bool skipp;
#if defined(USE_MPI_FILE_INPUT)
  MPI_Status status;
#endif

  int strlength1, strlength2;
  int offset_start_1, offset_end_1, offset_start_2, offset_end_2;
  int nextchar_end;
  bool donep;
#if 0
  int nchars1 = 0, nchars2 = 0;		/* Doesn't need to be saved as a field in Inbuffer_T. */
#endif

  /* Need to receive nextchar_end because of the difference between
     filecontents end ('\0') and FILE * end (EOF) */

  debug(fprintf(stdout,"Worker %d: accessing parser thread directly.  ",this->myid));
  Master_self_interface(this->master,&this->nextchar,&nextchar_end,
			&offset_start_1,&offset_start_2,&offset_end_1,&offset_end_2,
			&filestring1,&filestring2,&donep);

#if defined(HAVE_ZLIB) && defined(HAVE_BZLIB)
  if (this->gzipped == NULL && this->bzipped == NULL) {
    debug(fprintf(stdout,"Received offsets %d..%d and %d..%d, nextchars %c..%c, donep %d\n",
		  offset_start_1,offset_end_1,offset_start_2,offset_end_2,this->nextchar,nextchar_end,donep));

    FREE_IN(this->filecontents1_alloc);
    FREE_IN(this->filecontents2_alloc);
    this->filecontents1 = (char *) NULL;
    this->filecontents2 = (char *) NULL;

  } else {
    this->filecontents1 = this->filecontents1_alloc = Filestring_extract(&strlength1,filestring1);
    this->filecontents2 = this->filecontents2_alloc = Filestring_extract(&strlength2,filestring2);
    debug(fprintf(stdout,"Received filestrings of length %d and %d\n",strlength1,strlength2));
  }

#elif defined(HAVE_ZLIB)
  if (this->gzipped == NULL) {
    debug(fprintf(stdout,"Received offsets %d..%d and %d..%d, nextchars %c..%c, donep %d\n",
		  offset_start_1,offset_end_1,offset_start_2,offset_end_2,this->nextchar,nextchar_end,donep));

    FREE_IN(this->filecontents1_alloc);
    FREE_IN(this->filecontents2_alloc);
    this->filecontents1 = (char *) NULL;
    this->filecontents2 = (char *) NULL;

  } else {
    this->filecontents1 = this->filecontents1_alloc = Filestring_extract(&strlength1,filestring1);
    this->filecontents2 = this->filecontents2_alloc = Filestring_extract(&strlength2,filestring2);
    debug(fprintf(stdout,"Received filestrings of length %d and %d\n",strlength1,strlength2));
  }

#elif defined(HAVE_BZLIB)
  if (this->bzipped == NULL) {
    debug(fprintf(stdout,"Received offsets %d..%d and %d..%d, nextchars %c..%c, donep %d\n",
		  offset_start_1,offset_end_1,offset_start_2,offset_end_2,this->nextchar,nextchar_end,donep));

    FREE_IN(this->filecontents1_alloc);
    FREE_IN(this->filecontents2_alloc);
    this->filecontents1 = (char *) NULL;
    this->filecontents2 = (char *) NULL;

  } else {
    this->filecontents1 = this->filecontents1_alloc = Filestring_extract(&strlength1,filestring1);
    this->filecontents2 = this->filecontents2_alloc = Filestring_extract(&strlength2,filestring2);
    debug(fprintf(stdout,"Received filestrings of length %d and %d\n",strlength1,strlength2));
  }

#else
  debug(fprintf(stdout,"Received offsets %d..%d and %d..%d, nextchars %c..%c, donep %d\n",
		offset_start_1,offset_end_1,offset_start_2,offset_end_2,this->nextchar,nextchar_end,donep));

  FREE_IN(this->filecontents1_alloc);
  FREE_IN(this->filecontents2_alloc);
  this->filecontents1 = (char *) NULL;
  this->filecontents2 = (char *) NULL;
#endif

  Filestring_free(&filestring2);
  Filestring_free(&filestring1);


  if (this->filecontents1 == NULL) {
#if defined(USE_MPI_FILE_INPUT)
    MPI_File_seek(this->input,offset_start_1,MPI_SEEK_SET);
    this->filecontents1 = this->filecontents1_alloc = (char *) MALLOC_IN((offset_end_1 - offset_start_1 + 1) * sizeof(char));
    MPI_File_read(this->input,this->filecontents1,offset_end_1 - offset_start_1,MPI_CHAR,&status);
    this->filecontents1[offset_end_1 - offset_start_1] = '\0';

    if (this->input2 != NULL) {
      MPI_File_seek(this->input2,offset_start_2,MPI_SEEK_SET);
      this->filecontents2 = this->filecontents2_alloc = (char *) MALLOC_IN((offset_end_2 - offset_start_2 + 1) * sizeof(char));
      MPI_File_read(this->input2,this->filecontents2,offset_end_2 - offset_start_2,MPI_CHAR,&status);
      this->filecontents2[offset_end_2 - offset_start_2] = '\0';
    }
    
#else
#ifdef HAVE_FSEEKO
    fseeko(this->input,offset_start_1,SEEK_SET);
#else
    fseek(this->input,offset_start_1,SEEK_SET);
#endif
    this->filecontents1 = this->filecontents1_alloc = (char *) MALLOC_IN((offset_end_1 - offset_start_1 + 1) * sizeof(char));
    fread(this->filecontents1,offset_end_1 - offset_start_1,sizeof(char),this->input);
    this->filecontents1[offset_end_1 - offset_start_1] = '\0';
    if (this->input2 != NULL) {
#ifdef HAVE_FSEEKO
      fseeko(this->input2,offset_start_2,SEEK_SET);
#else
      fseek(this->input2,offset_start_2,SEEK_SET);
#endif
      this->filecontents2 = this->filecontents2_alloc = (char *) MALLOC_IN((offset_end_2 - offset_start_2 + 2) * sizeof(char));
      fread(this->filecontents2,offset_end_2 - offset_start_2,sizeof(char),this->input);
      this->filecontents2[offset_end_2 - offset_start_2] = '\0';
    }
#endif
  }

  /* Read from filecontents */
  while (nread < this->nspaces &&
	 (queryseq1 = Shortread_read_filecontents(&this->nextchar,&queryseq2,
						  &this->filecontents1,&this->filecontents2,&this->input,&this->input2,
#ifdef USE_MPI_FILE_INPUT
						  workers_comm,
#endif
						  &this->files,&this->nfiles,
						  skipp = (this->inputid % part_interval != part_modulus))) != NULL) {
    if (skipp) {
#if 0
      /* Shortread procedures won't allocate in this situation */
      Shortread_free(&queryseq1);
      if (queryseq2 != NULL) {
	Shortread_free(&queryseq2);
      }
#endif

    } else if (filter_if_both_p == true &&
	       Shortread_filterp(queryseq1) == true && (queryseq2 == NULL || Shortread_filterp(queryseq2) == true)) {
      Shortread_free(&queryseq1);
      if (queryseq2 != NULL) {
	Shortread_free(&queryseq2);
      }
      
    } else if (filter_if_both_p == false &&
	       (Shortread_filterp(queryseq1) == true || (queryseq2 != NULL && Shortread_filterp(queryseq2) == true))) {
      Shortread_free(&queryseq1);
      if (queryseq2 != NULL) {
	Shortread_free(&queryseq2);
      }
      
    } else {
      this->buffer[nread++] = Request_new(this->requestid++,queryseq1,queryseq2);
    }
    this->inputid++;
  }

  this->nleft = nread;
  this->ptr = 0;

  /* Need to set this to the FILE * end (EOF at end of file), and not the filecontents end (always '\0') */
  this->nextchar = nextchar_end;

#ifdef USE_MPI
  debug(printf("Worker %d: ",this->myid));
#endif
  debug(printf("this->nextchar (nextchar_end) is %c (%d)\n",this->nextchar,this->nextchar));

  return nread;
}



/* Used by ranks 1..n to communicate with Master_mpi_interface thread of rank 0 */
/* Returns number of requests read */
static unsigned int
fill_buffer_slave (T this) {
  unsigned int nread = 0;
  Shortread_T queryseq1, queryseq2;
  bool skipp, donep;

  int strlength1, strlength2;
  MPI_Status status;
  int offset_start_1, offset_end_1, offset_start_2, offset_end_2;
  int nextchar_end;
#if 0
  int nchars1 = 0, nchars2 = 0;		/* Doesn't need to be saved as a field in Inbuffer_T. */
#endif

  /* Need to receive nextchar_end because of the difference between
     filecontents end ('\0') and FILE * end (EOF) */

  debug(fprintf(stdout,"Worker %d: sending notification to master process.  ",this->myid));
  MPI_SEND(&this->nfiles,1,MPI_INT,/*dest*/0,/*tag*/MPI_TAG_WANT_INPUT,MPI_COMM_WORLD);
  MPI_RECV(&this->nextchar,1,MPI_INT,/*source*/0,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD,&status);
  MPI_RECV(&nextchar_end,1,MPI_INT,/*source*/0,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD,&status);
  MPI_RECV(&donep,1,MPI_UNSIGNED_CHAR,/*source*/0,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD,&status);

#if defined(HAVE_ZLIB) && defined(HAVE_BZLIB)
  if (this->gzipped == NULL && this->bzipped == NULL) {
    MPI_RECV(&offset_start_1,1,MPI_INT,/*source*/0,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD,&status);
    MPI_RECV(&offset_start_2,1,MPI_INT,/*source*/0,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD,&status);
    MPI_RECV(&offset_end_1,1,MPI_INT,/*source*/0,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD,&status);
    MPI_RECV(&offset_end_2,1,MPI_INT,/*source*/0,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD,&status);
    debug(fprintf(stdout,"Received offsets %d..%d and %d..%d, nextchars %c..%c, donep %d\n",
		  offset_start_1,offset_end_1,offset_start_2,offset_end_2,this->nextchar,nextchar_end,donep));

    FREE_IN(this->filecontents1_alloc);
    FREE_IN(this->filecontents2_alloc);
    this->filecontents1 = (char *) NULL;
    this->filecontents2 = (char *) NULL;

  } else {
    this->filecontents1 = this->filecontents1_alloc =
      Filestring_recv(&strlength1,/*source*/0,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD);
    this->filecontents2 = this->filecontents2_alloc =
      Filestring_recv(&strlength2,/*source*/0,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD);
    debug(fprintf(stdout,"Received filestrings of length %d and %d\n",strlength1,strlength2));
  }

#elif defined(HAVE_ZLIB)
  if (this->gzipped == NULL) {
    MPI_RECV(&offset_start_1,1,MPI_INT,/*source*/0,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD,&status);
    MPI_RECV(&offset_start_2,1,MPI_INT,/*source*/0,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD,&status);
    MPI_RECV(&offset_end_1,1,MPI_INT,/*source*/0,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD,&status);
    MPI_RECV(&offset_end_2,1,MPI_INT,/*source*/0,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD,&status);
    debug(fprintf(stdout,"Received offsets %d..%d and %d..%d, nextchars %c..%c, donep %d\n",
		  offset_start_1,offset_end_1,offset_start_2,offset_end_2,this->nextchar,nextchar_end,donep));
    debug(fprintf(stdout,"Received offsets %d..%d and %d..%d\n",offset_start_1,offset_end_1,offset_start_2,offset_end_2));

    FREE_IN(this->filecontents1_alloc);
    FREE_IN(this->filecontents2_alloc);
    this->filecontents1 = (char *) NULL;
    this->filecontents2 = (char *) NULL;

  } else {
    this->filecontents1 = this->filecontents1_alloc =
      Filestring_recv(&strlength1,/*source*/0,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD);
    this->filecontents2 = this->filecontents2_alloc =
      Filestring_recv(&strlength2,/*source*/0,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD);
    debug(fprintf(stdout,"Received filestrings of length %d and %d\n",strlength1,strlength2));
  }

#elif defined(HAVE_BZLIB)
  if (this->bzipped == NULL) {
    MPI_RECV(&offset_start_1,1,MPI_INT,/*source*/0,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD,&status);
    MPI_RECV(&offset_start_2,1,MPI_INT,/*source*/0,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD,&status);
    MPI_RECV(&offset_end_1,1,MPI_INT,/*source*/0,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD,&status);
    MPI_RECV(&offset_end_2,1,MPI_INT,/*source*/0,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD,&status);
    debug(fprintf(stdout,"Received offsets %d..%d and %d..%d, nextchars %c..%c, donep %d\n",
		  offset_start_1,offset_end_1,offset_start_2,offset_end_2,this->nextchar,nextchar_end,donep));
    debug(fprintf(stdout,"Received offsets %d..%d and %d..%d\n",offset_start_1,offset_end_1,offset_start_2,offset_end_2));

    FREE_IN(this->filecontents1_alloc);
    FREE_IN(this->filecontents2_alloc);
    this->filecontents1 = (char *) NULL;
    this->filecontents2 = (char *) NULL;

  } else {
    this->filecontents1 = this->filecontents1_alloc =
      Filestring_recv(&strlength1,/*source*/0,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD);
    this->filecontents2 = this->filecontents2_alloc =
      Filestring_recv(&strlength2,/*source*/0,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD);
    debug(fprintf(stdout,"Received filestrings of length %d and %d\n",strlength1,strlength2));
  }

#else
  MPI_RECV(&offset_start_1,1,MPI_INT,/*source*/0,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD,&status);
  MPI_RECV(&offset_start_2,1,MPI_INT,/*source*/0,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD,&status);
  MPI_RECV(&offset_end_1,1,MPI_INT,/*source*/0,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD,&status);
  MPI_RECV(&offset_end_2,1,MPI_INT,/*source*/0,/*tag*/MPI_TAG_GIVE_INPUT,MPI_COMM_WORLD,&status);
  debug(fprintf(stdout,"Received offsets %d..%d and %d..%d, nextchars %c..%c, donep %d\n",
		offset_start_1,offset_end_1,offset_start_2,offset_end_2,this->nextchar,nextchar_end,donep));

  FREE_IN(this->filecontents1_alloc);
  FREE_IN(this->filecontents2_alloc);
  this->filecontents1 = (char *) NULL;
  this->filecontents2 = (char *) NULL;
#endif


  if (this->filecontents1 == NULL) {
#if defined(USE_MPI_FILE_INPUT)
    MPI_File_seek(this->input,offset_start_1,MPI_SEEK_SET);
    this->filecontents1 = this->filecontents1_alloc = (char *) MALLOC_IN((offset_end_1 - offset_start_1 + 1) * sizeof(char));
    MPI_File_read(this->input,this->filecontents1,offset_end_1 - offset_start_1,MPI_CHAR,&status);
    this->filecontents1[offset_end_1 - offset_start_1] = '\0';

    if (this->input2 != NULL) {
      MPI_File_seek(this->input2,offset_start_2,MPI_SEEK_SET);
      this->filecontents2 = this->filecontents2_alloc = (char *) MALLOC_IN((offset_end_2 - offset_start_2 + 1) * sizeof(char));
      MPI_File_read(this->input2,this->filecontents2,offset_end_2 - offset_start_2,MPI_CHAR,&status);
      this->filecontents2[offset_end_2 - offset_start_2] = '\0';
    }
    
#else
#ifdef HAVE_FSEEKO
    fseeko(this->input,offset_start_1,SEEK_SET);
#else
    fseek(this->input,offset_start_1,SEEK_SET);
#endif
    this->filecontents1 = this->filecontents1_alloc = (char *) MALLOC_IN((offset_end_1 - offset_start_1 + 1) * sizeof(char));
    fread(this->filecontents1,offset_end_1 - offset_start_1,sizeof(char),this->input);
    this->filecontents1[offset_end_1 - offset_start_1] = '\0';
    if (this->input2 != NULL) {
#ifdef HAVE_FSEEKO
      fseeko(this->input2,offset_start_2,SEEK_SET);
#else
      fseek(this->input2,offset_start_2,SEEK_SET);
#endif
      this->filecontents2 = this->filecontents2_alloc = (char *) MALLOC_IN((offset_end_2 - offset_start_2 + 2) * sizeof(char));
      fread(this->filecontents2,offset_end_2 - offset_start_2,sizeof(char),this->input);
      this->filecontents2[offset_end_2 - offset_start_2] = '\0';
    }
#endif
  }

  /* Read from filecontents */
  while (nread < this->nspaces &&
	 (queryseq1 = Shortread_read_filecontents(&this->nextchar,&queryseq2,
						  &this->filecontents1,&this->filecontents2,&this->input,&this->input2,
#ifdef USE_MPI_FILE_INPUT
						  workers_comm,
#endif
						  &this->files,&this->nfiles,
						  skipp = (this->inputid % part_interval != part_modulus))) != NULL) {
    if (skipp) {
#if 0
      /* Shortread procedures won't allocate in this situation */
      Shortread_free(&queryseq1);
      if (queryseq2 != NULL) {
	Shortread_free(&queryseq2);
      }
#endif

    } else if (filter_if_both_p == true &&
	       Shortread_filterp(queryseq1) == true && (queryseq2 == NULL || Shortread_filterp(queryseq2) == true)) {
      Shortread_free(&queryseq1);
      if (queryseq2 != NULL) {
	Shortread_free(&queryseq2);
      }
      
    } else if (filter_if_both_p == false &&
	       (Shortread_filterp(queryseq1) == true || (queryseq2 != NULL && Shortread_filterp(queryseq2) == true))) {
      Shortread_free(&queryseq1);
      if (queryseq2 != NULL) {
	Shortread_free(&queryseq2);
      }
      
    } else {
      this->buffer[nread++] = Request_new(this->requestid++,queryseq1,queryseq2);
    }
    this->inputid++;
  }

  this->nleft = nread;
  this->ptr = 0;

  /* Need to set this to the FILE * end (EOF at end of file), and not the filecontents end (always '\0') */
  this->nextchar = nextchar_end;

#ifdef USE_MPI
  debug(printf("Worker %d: ",this->myid));
#endif
  debug(printf("this->nextchar (nextchar_end) is %c (%d)\n",this->nextchar,this->nextchar));

  return nread;
}

#elif defined(GSNAP)

/* Returns number of requests read */
static unsigned int
fill_buffer (T this) {
  unsigned int nread = 0;
  Shortread_T queryseq1, queryseq2;
  bool skipp;
  int nchars1 = 0, nchars2 = 0;		/* Returned only because MPI master needs it.  Doesn't need to be saved as a field in Inbuffer_T. */

  while (nread < this->nspaces &&
	 (queryseq1 = Shortread_read(&this->nextchar,&nchars1,&nchars2,&queryseq2,
				     &this->input,&this->input2,
#ifdef HAVE_ZLIB
				     &this->gzipped,&this->gzipped2,
#endif
#ifdef HAVE_BZLIB
				     &this->bzipped,&this->bzipped2,
#endif
				     &this->files,&this->nfiles,skipp = (this->inputid % part_interval != part_modulus))) != NULL) {
    if (skipp) {
#if 0
      /* Shortread procedures won't allocate in this situation */
      Shortread_free(&queryseq1);
      if (queryseq2 != NULL) {
	Shortread_free(&queryseq2);
      }
#endif
      
    } else if (filter_if_both_p == true &&
	       Shortread_filterp(queryseq1) == true && (queryseq2 == NULL || Shortread_filterp(queryseq2) == true)) {
      Shortread_free(&queryseq1);
      if (queryseq2 != NULL) {
	Shortread_free(&queryseq2);
      }
      
    } else if (filter_if_both_p == false &&
	       (Shortread_filterp(queryseq1) == true || (queryseq2 != NULL && Shortread_filterp(queryseq2) == true))) {
      Shortread_free(&queryseq1);
      if (queryseq2 != NULL) {
	Shortread_free(&queryseq2);
      }
      
    } else {
      this->buffer[nread++] = Request_new(this->requestid++,queryseq1,queryseq2);
    }
    this->inputid++;
  }

  this->nleft = nread;
  this->ptr = 0;

  return nread;
}

#else
	 
/* GMAP version */
/* Returns number of requests read */
static unsigned int
fill_buffer (T this) {
  unsigned int nread = 0;
#if 0
  unsigned int nchars = 0U;
#endif
  Sequence_T queryseq;

  while (nread < this->nspaces &&
	 (queryseq = Sequence_read_multifile(&this->nextchar,&this->input,
					     &this->files,&this->nfiles)) != NULL) {
    if (this->inputid % part_interval != part_modulus) {
      Sequence_free(&queryseq);
    } else {
      debug(printf("inbuffer creating request %d\n",this->requestid));
      this->buffer[nread++] = Request_new(this->requestid++,queryseq);
#if 0
      nchars += Sequence_fulllength(queryseq);
#endif
    }
    this->inputid++;
  }

  this->nleft = nread;
  this->ptr = 0;

  return nread;
}

#endif


#ifndef USE_MPI
/* No need to lock, since only main thread calls */
/* Returns nread to give to Outbuffer_new */
unsigned int
Inbuffer_fill_init (T this) {
  unsigned int nread;

  debug(printf("inbuffer filling initially\n"));
  nread = fill_buffer(this);
  debug(printf("inbuffer read %d sequences\n",nread));

  return nread;
}
#endif
  

Request_T
#ifdef GSNAP
Inbuffer_get_request (T this)
#else
Inbuffer_get_request (Sequence_T *pairalign_segment, T this) 
#endif
{
  Request_T request;
  unsigned int nread;

#if defined(HAVE_PTHREAD)
  pthread_mutex_lock(&this->lock);
#endif
  
  if (this->nleft > 0) {
    request = this->buffer[this->ptr++];
    this->nleft -= 1;

  } else if (this->nextchar == EOF) {
    /* ? Causes stall at end */
    /* Already know it is pointless to fill buffer */
    Outbuffer_add_nread(this->outbuffer,/*nread*/0);
    request = NULL;

  } else {
#ifdef USE_MPI
    debug(printf("Worker %d: ",this->myid));
#endif
    debug(printf("inbuffer filling with nextchar %c (%d)\n",this->nextchar,this->nextchar));

#ifndef GSNAP
    if (user_pairalign_p == true) {
      /* assert(this->nspaces == 1) */
      if (this->pairalign_segment != NULL) {
	Sequence_free(&this->pairalign_segment);
      }
      this->pairalign_segment = Sequence_read_unlimited(&this->nextchar,stdin);
      debug(printf("  but first reading usersegment, got nextchar %c\n",this->nextchar));
    }
#endif

#ifdef USE_MPI
    if (this->myid == 0) {
      nread = fill_buffer_master(this);
    } else {
      nread = fill_buffer_slave(this);
    }
#else
    nread = fill_buffer(this);
#endif

    Outbuffer_add_nread(this->outbuffer,nread);
#ifdef USE_MPI
    debug(printf("Worker %d: ",this->myid));
#endif
    debug(printf("inbuffer read %d sequences\n",nread));
    
    if (nread == 0) {
      /* Still empty */
      request = NULL;
    } else {
      request = this->buffer[this->ptr++];
      this->nleft -= 1;
    }
  }

#ifndef GSNAP
  *pairalign_segment = this->pairalign_segment;
#endif

#if defined(HAVE_PTHREAD)
  pthread_mutex_unlock(&this->lock);
#endif

  return request;
}


#ifndef GSNAP
/* Same as Inbuffer_get_request, but leaves sequence in buffer.  Used by GMAP for selfalign feature. */
Request_T
Inbuffer_first_request (T this) {
  Request_T request;
  unsigned int nread;

#if defined(HAVE_PTHREAD)
  pthread_mutex_lock(&this->lock);
#endif
  
  if (this->nleft > 0) {
    request = this->buffer[this->ptr/*++*/];
    /* this->nleft -= 1; */

  } else {
    debug(printf("inbuffer filling\n"));
    nread = fill_buffer(this);
    Outbuffer_add_nread(this->outbuffer,nread);
    debug(printf("inbuffer read %d sequences\n",nread));
    
    if (nread == 0) {
      /* Still empty */
      request = NULL;
    } else {
      request = this->buffer[this->ptr/*++*/];
      /* this->nleft -= 1; */
    }
  }

#if defined(HAVE_PTHREAD)
  pthread_mutex_unlock(&this->lock);
#endif

  return request;
}
#endif


