static char rcsid[] = "$Id: bzip2.c 83593 2013-01-16 22:59:40Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "bzip2.h"

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include "mem.h"


#ifdef HAVE_BZLIB
#include <bzlib.h>
#endif

#define BUFFER_SIZE 131072

#define T Bzip2_T
struct T {
#ifdef HAVE_BZLIB
  BZFILE *bzfile;
#endif

  bool eofp;

  char *buffer;
  int navail;
  int bufferi;
};


T
Bzip2_new (char *filename) {
  T new = (T) MALLOC(sizeof(*new));
  int error;
  FILE *fp;

#ifndef HAVE_BZLIB
  fprintf(stderr,"Bzip2_new called when BZLIB is unavailable\n");
  exit(9);

#else
  fp = fopen(filename,"r");
  if (fp == NULL) {
    fprintf(stderr,"Unable to open bzip2 file %s\n",filename);
    exit(9);
  }

  new->bzfile = BZ2_bzReadOpen(&error,fp,/*small*/0,/*verbosity*/0,/*unused*/NULL,/*nunused*/0);
  if (error != BZ_OK) {
    fprintf(stderr,"BZ2_bzReadOpen gave error %d\n",error);
    exit(9);
  }
#endif

  new->eofp = false;
  new->buffer = (char *) CALLOC(BUFFER_SIZE,sizeof(char));
  new->navail = 0;
  new->bufferi = 0;
  
  return new;
}

void
Bzip2_free (T *old) {
#ifdef HAVE_BZLIB
  int error;

  BZ2_bzReadClose(&error,(*old)->bzfile);
  if (error != BZ_OK) {
    fprintf(stderr,"BZ2_bzReadClose gave error %d\n",error);
    exit(9);
  }
#endif
  
  FREE((*old)->buffer);
  FREE(*old);
  return;
}


int
bzgetc (T this) {
#ifdef HAVE_BZLIB
  int error;
#endif

  if (this->navail > 0) {
    this->navail -= 1;
    return this->buffer[this->bufferi++];

  } else if (this->eofp == true) {
    return EOF;

  } else {
#ifdef HAVE_BZLIB
    this->navail = BZ2_bzRead(&error,this->bzfile,this->buffer,BUFFER_SIZE);
    if (error == BZ_STREAM_END) {
      this->eofp = true;
    } else if (error != BZ_OK) {
      fprintf(stderr,"BZ2_bzRead gave error %d\n",error);
      exit(9);
    }
#endif
    this->bufferi = 1;
    this->navail -= 1;
    return this->buffer[0];
  }
}


bool
bzeof (T this) {
  if (this->navail > 0) {
    return false;
  } else {
    return this->eofp;
  }
}

char *
bzgets (T this, char *buffer, int maxlength) {
  int n = 0;
  int c = '\0';

  if (this->navail == 0 && this->eofp == true) {
    return (char *) NULL;
  } else {
    while (c != '\n' && n < maxlength - 1 && (c = bzgetc(this)) != EOF) {
      buffer[n++] = c;
    }
    buffer[n] = '\0';
    return buffer;
  }
}


#if 0
int
main (int argc, char *argv[]) {
  char *filename;
  T bzip2;
  char string[100];

  filename = argv[1];
  bzip2 = Bzip2_new(filename);
  while (bzeof(bzip2) == false) {
    bzgets(bzip2,string,100);
    printf("%s",string);
  }
  Bzip2_free(&bzip2);
  return 0;
}
#endif

