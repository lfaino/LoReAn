static char rcsid[] = "$Id: except.c 47650 2011-09-17 04:45:25Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "except.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For strcat */
#include "assert.h"
#include "bool.h"

#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif

static bool raisep = true;
static bool threadedp = false;

void
Except_inactivate () {
  raisep = false;
}

/* Used only for non-threaded runs */
Except_Frame_T global_except_stack = NULL;

#ifdef HAVE_PTHREAD
static pthread_key_t global_except_key;

void
Except_init_pthread () {
  threadedp = true;
  pthread_key_create(&global_except_key,NULL);
  return;
}

void
Except_term_pthread () {
  /* Do not delete global_except_key, because worker threads might still need it */
  /* pthread_key_delete(global_except_key); */
  return;
}

void
Except_stack_create () {
  Except_Frame_T stackptr;

  stackptr = (Except_Frame_T) malloc(sizeof(*stackptr));
  pthread_setspecific(global_except_key,stackptr);
  return;
}

void
Except_stack_destroy () {
  Except_Frame_T stackptr;

  stackptr = (Except_Frame_T) pthread_getspecific(global_except_key);
  free(stackptr);
  return;
}
#endif

void
Except_link_stack (Except_Frame_T frameptr) {
#ifdef HAVE_PTHREAD
  Except_Frame_T stackptr;
#endif

  if (threadedp == false) {
    frameptr->prev = global_except_stack;
    global_except_stack = frameptr;
  } else {
#ifdef HAVE_PTHREAD
    stackptr = (Except_Frame_T) pthread_getspecific(global_except_key);
    frameptr->prev = stackptr;
    stackptr = frameptr;
#endif
  }

  return;
}

Except_Frame_T
Except_advance_stack () {
#ifdef HAVE_PTHREAD
  Except_Frame_T stackptr;

  if (threadedp == false) {
    global_except_stack = global_except_stack->prev;
    return global_except_stack;
  } else {
    stackptr = (Except_Frame_T) pthread_getspecific(global_except_key);
    stackptr = (stackptr)->prev;
    return stackptr;
  }

#else
  global_except_stack = global_except_stack->prev;
  return global_except_stack;
#endif
}

void
Except_raise (const Except_T *e, const char *file, int line) {
  Except_Frame_T frameptr;
  char message[512], piece[128];
#ifdef HAVE_PTHREAD
  Except_Frame_T stackptr;
#endif

  assert(e);
  message[0] = '\0';
  if (e->reason) {
    sprintf(piece," %s ", e->reason);
    strcat(message,piece);
  } else {
    sprintf(piece," at 0x%p",(void *) e);
    strcat(message,piece);
  }
  if (file && line > 0) {
    sprintf(piece," raised at %s:%d",file,line);
    strcat(message,piece);
  }
  fprintf(stderr,"Exception: %s\n",message);
  fflush(stderr);

  if (threadedp == false) {
    frameptr = global_except_stack;
  } else {
#ifdef HAVE_PTHREAD
    stackptr = (Except_Frame_T) pthread_getspecific(global_except_key);
    frameptr = stackptr;
#endif    
  }

  if (frameptr == NULL) {
    fprintf(stderr,"Uncaught exception: %s\n",message);
    fflush(stderr);
    abort();
  } else {
    frameptr->exception = e;
    frameptr->file = file;
    frameptr->line = line;
    if (raisep == true) {
      longjmp(frameptr->env,EXCEPT_RAISED);
    } else {
      fprintf(stderr,"Aborting...\n");
      abort();
    }
  }
  return;
}


