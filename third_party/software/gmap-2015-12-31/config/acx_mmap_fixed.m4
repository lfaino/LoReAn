AC_DEFUN([ACX_MMAP_FIXED], [
AC_REQUIRE([AC_CANONICAL_HOST])
AC_LANG_SAVE
AC_LANG_C
acx_mmap_fixed_ok=no

AC_MSG_CHECKING([for working mmap with MAP_FIXED])

        AC_TRY_LINK([
/* Thanks to Mike Haertel and Jim Avera for this test.
   Here is a matrix of mmap possibilities:
	mmap private not fixed
	mmap private fixed at somewhere currently unmapped
	mmap private fixed at somewhere already mapped
	mmap shared not fixed
	mmap shared fixed at somewhere currently unmapped
	mmap shared fixed at somewhere already mapped
   For private mappings, we should verify that changes cannot be read()
   back from the file, nor mmap's back from the file at a different
   address.  (There have been systems where private was not correctly
   implemented like the infamous i386 svr4.0, and systems where the
   VM page cache was not coherent with the file system buffer cache
   like early versions of FreeBSD and possibly contemporary NetBSD.)
   For shared mappings, we should conversely verify that changes get
   propagated back to all the places they're supposed to be.

   Grep wants private fixed already mapped.
   The main things grep needs to know about mmap are:
   * does it exist and is it safe to write into the mmap'd area
   * how to use it (BSD variants)  */

#include <fcntl.h>
#include <sys/mman.h>
#include <stdlib.h>
#include <unistd.h>

#if !STDC_HEADERS && !HAVE_STDLIB_H
char *malloc ();
#endif

/* This mess was copied from the GNU getpagesize.h.  */
#if !HAVE_GETPAGESIZE
/* Assume that all systems that can run configure have sys/param.h.  */
# if !HAVE_SYS_PARAM_H
#  define HAVE_SYS_PARAM_H 1
# endif

# ifdef _SC_PAGESIZE
#  define getpagesize() sysconf(_SC_PAGESIZE)
# else /* no _SC_PAGESIZE */
#  if HAVE_SYS_PARAM_H
#   include <sys/param.h>
#   ifdef EXEC_PAGESIZE
#    define getpagesize() EXEC_PAGESIZE
#   else /* no EXEC_PAGESIZE */
#    ifdef NBPG
#     define getpagesize() NBPG * CLSIZE
#     ifndef CLSIZE
#      define CLSIZE 1
#     endif /* no CLSIZE */
#    else /* no NBPG */
#     ifdef NBPC
#      define getpagesize() NBPC
#     else /* no NBPC */
#      ifdef PAGESIZE
#       define getpagesize() PAGESIZE
#      endif /* PAGESIZE */
#     endif /* no NBPC */
#    endif /* no NBPG */
#   endif /* no EXEC_PAGESIZE */
#  else /* no HAVE_SYS_PARAM_H */
#   define getpagesize() 8192	/* punt totally */
#  endif /* no HAVE_SYS_PARAM_H */
# endif /* no _SC_PAGESIZE */

#endif /* no HAVE_GETPAGESIZE */
],
[
  char *data, *data2, *data3;
  size_t i, pagesize;
  int fd;

  pagesize = getpagesize ();

  /* First, make a file with some known garbage in it. */
  data = (char *) malloc (pagesize);
  if (!data)
    exit (1);
  for (i = 0; i < pagesize; ++i)
    *(data + i) = rand ();
  umask (0);
  fd = creat ("conftest.mmap", 0600);
  if (fd < 0)
    exit (1);
  if (write (fd, data, pagesize) != pagesize)
    exit (1);
  close (fd);

  /* Next, try to mmap the file at a fixed address which already has
     something else allocated at it.  If we can, also make sure that
     we see the same garbage.  */
  fd = open ("conftest.mmap", O_RDWR);
  if (fd < 0)
    exit (1);
  data2 = (char *) malloc (2 * pagesize);
  if (!data2)
    exit (1);
  data2 += (pagesize - ((size_t) data2 & (pagesize - 1))) & (pagesize - 1);
  if (data2 != mmap (data2, pagesize, PROT_READ | PROT_WRITE,
                     MAP_PRIVATE | MAP_FIXED, fd, 0L))
    exit (1);
  for (i = 0; i < pagesize; ++i)
    if (*(data + i) != *(data2 + i))
      exit (1);

  /* Finally, make sure that changes to the mapped area do not
     percolate back to the file as seen by read().  (This is a bug on
     some variants of i386 svr4.0.)  */
  for (i = 0; i < pagesize; ++i)
    *(data2 + i) = *(data2 + i) + 1;
  data3 = (char *) malloc (pagesize);
  if (!data3)
    exit (1);
  if (read (fd, data3, pagesize) != pagesize)
    exit (1);
  for (i = 0; i < pagesize; ++i)
    if (*(data + i) != *(data3 + i))
      exit (1);
  close (fd);
  exit (0);
],
	              [acx_mmap_fixed_ok=yes])

AC_MSG_RESULT($acx_mmap_fixed_ok)

AC_LANG_RESTORE
])dnl ACX_MMAP_FIXED
