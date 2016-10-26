/* $Id: except.h 157225 2015-01-22 18:47:23Z twu $ */
#ifndef EXCEPT_INCLUDED
#define EXCEPT_INCLUDED
#ifdef HAVE_CONFIG_H
#include <config.h>		/* For HAVE_PTHREAD */
#endif

#include <setjmp.h>

#define T Except_T
typedef struct T {
  char *reason;
} T;

typedef struct Except_Frame_T *Except_Frame_T;
struct Except_Frame_T {
  Except_Frame_T prev;
  jmp_buf env;
  const char *file;
  int line;
  const T *exception;
};

enum {EXCEPT_ENTERED = 0, EXCEPT_RAISED, EXCEPT_HANDLED, EXCEPT_FINALIZED};

extern const Except_T Assert_Failed;

extern void
Except_inactivate ();

#ifdef HAVE_PTHREAD
extern void
Except_init_pthread ();
extern void
Except_term_pthread ();
extern void
Except_stack_create ();
extern void
Except_stack_destroy ();
#endif

extern void
Except_link_stack (Except_Frame_T frameptr);

extern Except_Frame_T
Except_advance_stack ();

extern void
Except_raise (const T *e, const char *file, int line);

#define RAISE(e) Except_raise(&(e), __FILE__, __LINE__)
#define RERAISE Except_raise(frame.exception, frame.file, frame.line)
#define RETURN switch (Except_advance_stack(),0) default: return

#define TRY \
do { \
    volatile int except_flag; \
    struct Except_Frame_T frame;      \
    Except_link_stack(&frame); \
    except_flag = setjmp(frame.env); \
    if (except_flag == EXCEPT_ENTERED) {

#define EXCEPT(e) \
      if (except_flag == EXCEPT_ENTERED) { \
        Except_advance_stack(); \
      } \
    } else if (frame.exception == &(e)) {	\
      except_flag = EXCEPT_HANDLED;

#define ELSE \
      if (except_flag == EXCEPT_ENTERED) { \
        Except_advance_stack(); \
      } \
    } else { \
      except_flag = EXCEPT_HANDLED;


#define FINALLY \
      if (except_flag == EXCEPT_ENTERED) { \
        Except_advance_stack(); \
      } \
    } { \
      if (except_flag == EXCEPT_ENTERED) { \
        except_flag = EXCEPT_FINALIZED; \
      }					  

#define END_TRY \
      if (except_flag == EXCEPT_ENTERED) { \
	Except_advance_stack(); \
     } \
   } \
   if (except_flag == EXCEPT_RAISED) { \
     RERAISE; \
   } \
} while (0)

#undef T
#endif

