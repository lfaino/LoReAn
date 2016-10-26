static char rcsid[] = "$Id: access.c 170327 2015-07-22 17:50:11Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "access.h"
#include "list.h"
#include "intlist.h"

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>		/* For strerror */
#include <errno.h>

/* <unistd.h> and <sys/types.h> included in access.h */
#include <sys/mman.h>		/* For mmap */

#define PROJECT_ID 42
#include <sys/ipc.h>
#include <sys/shm.h>		/* For shmat and shmdt */
#include <sys/sem.h>		/* For semaphores */

#ifdef USE_MPI
#include <mpi.h>
#endif

#ifdef HAVE_FCNTL_H
#include <fcntl.h>		/* For open */
#endif
#ifdef HAVE_SYS_STAT_H
#include <sys/stat.h>		/* For open and fstat */
#endif

#ifdef PAGESIZE_VIA_SYSCONF
#include <unistd.h>
#endif
#ifdef PAGESIZE_VIA_SYSCTL
#include <sys/sysctl.h>
#endif

#include "assert.h"
#include "mem.h"
#include "types.h"
#include "fopen.h"
#include "stopwatch.h"

#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
#else
#include "littleendian.h"
#endif


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


bool
Access_file_exists_p (char *filename) {
#ifdef HAVE_STRUCT_STAT64
  struct stat64 sb;
#else
  struct stat sb;
#endif

#ifdef HAVE_STAT64
  if (stat64(filename,&sb) == 0) {
    return true;
  } else {
    return false;
  }
#else
  if (stat(filename,&sb) == 0) {
    return true;
  } else {
    return false;
  }
#endif
}


off_t
Access_filesize (char *filename) {
#ifdef HAVE_STRUCT_STAT64
  struct stat64 sb;
#else
  struct stat sb;
#endif

#ifdef HAVE_STAT64
  stat64(filename,&sb);
#else
  stat(filename,&sb);
#endif
  debug(printf("filesize is %zu\n",sb.st_size));
  return sb.st_size;
}


size_t
Access_file_copy (char *dest_file, char *source_file) {
  size_t nbytes = 0;
  FILE *dest, *source;
  int c;

  if ((source = FOPEN_READ_BINARY(source_file)) == NULL) {
    fprintf(stderr,"Cannot open source file %s\n",source_file);
    return 0;

  } else if ((dest = FOPEN_WRITE_BINARY(dest_file)) == NULL) {
    fprintf(stderr,"Cannot open destination file %s\n",dest_file);
    fclose(source);
    return 0;

  } else {
    while ((c = fgetc(source)) != EOF) {
      fputc(c,dest);
      nbytes++;
    }
    fclose(dest);
    fclose(source);
    return nbytes;
  }
}


bool
Access_file_equal (char *file1, char *file2) {
  FILE *fp1, *fp2;
  int c1, c2;

  if ((fp1 = FOPEN_READ_BINARY(file1)) == NULL) {
    fprintf(stderr,"Cannot open file %s\n",file1);
    exit(9);

  } else if ((fp2 = FOPEN_READ_BINARY(file2)) == NULL) {
    fprintf(stderr,"Cannot open file %s\n",file2);
    fclose(fp1);
    exit(9);

  } else {
    c1 = fgetc(fp1);
    c2 = fgetc(fp2);
    while (c1 != EOF && c2 != EOF) {
      if (c1 != c2) {
	fclose(fp2);
	fclose(fp1);
	return false;
      }
      c1 = fgetc(fp1);
      c2 = fgetc(fp2);
    }
    fclose(fp2);
    fclose(fp1);

    if (c1 == EOF && c2 == EOF) {
      return true;
    } else {
      return false;
    }
  }
}


int
Access_fileio (char *filename) {
  int fd;

  if ((fd = open(filename,O_RDONLY,0764)) < 0) {
    fprintf(stderr,"Error: can't open file %s with open for reading\n",filename);
    exit(9);
  }
  return fd;
}

int
Access_fileio_rw (char *filename) {
  int fd;

  if ((fd = open(filename,O_RDWR | O_CREAT | O_TRUNC,0764)) < 0) {
    fprintf(stderr,"Error: can't open file %s with open for reading/writing\n",filename);
    exit(9);
  }
  return fd;
}


#ifndef WORDS_BIGENDIAN
/* Needed as a test on Macintosh machines */
static unsigned char
first_nonzero_char (size_t *i, char *filename) {
  size_t len;
  FILE *fp;
  unsigned char value = 0;
  void *p;

  len = (size_t) Access_filesize(filename);

  if ((fp = FOPEN_READ_BINARY(filename)) == NULL) {
    fprintf(stderr,"Error: can't open file %s with fopen\n",filename);
    exit(9);
  } else {
    *i = 0;
    p = (void *) &value;
    while ((size_t) *i < len && fread(p,sizeof(unsigned char),1,fp) > 0 &&
	   value == 0) {
      *i += 1;
    }
    if (value == 0) {
      *i = -1;
    }
    fclose(fp);
    return value;
  }
}

static UINT4
first_nonzero_uint (size_t *i, char *filename) {
  size_t len;
  FILE *fp;
  UINT4 value = 0;
  void *p;

  len = (size_t) Access_filesize(filename);

  if ((fp = FOPEN_READ_BINARY(filename)) == NULL) {
    fprintf(stderr,"Error: can't open file %s with fopen\n",filename);
    exit(9);
  } else {
    *i = 0;
    p = (void *) &value;
    while ((size_t) *i < len && fread(p,sizeof(UINT4),1,fp) > 0 &&
	   value == 0) {
      *i += 1;
    }
    if (value == 0) {
      *i = -1;
    }
    fclose(fp);
    return value;
  }
}

static UINT8
first_nonzero_uint8 (size_t *i, char *filename) {
  size_t len;
  FILE *fp;
  UINT8 value = 0;
  void *p;

  len = (size_t) Access_filesize(filename);

  if ((fp = FOPEN_READ_BINARY(filename)) == NULL) {
    fprintf(stderr,"Error: can't open file %s with fopen\n",filename);
    exit(9);
  } else {
    *i = 0;
    p = (void *) &value;
    while ((size_t) *i < len && fread(p,sizeof(UINT8),1,fp) > 0 &&
	   value == 0) {
      *i += 1;
    }
    if (value == 0) {
      *i = -1;
    }
    fclose(fp);
    return value;
  }
}
#endif



/************************************************************************
 *   Functions for shared memory and semaphores
 ************************************************************************/

static List_T shmem_memory = NULL;
static Intlist_T shmem_ids = NULL;
static Intlist_T semaphore_ids = NULL;


#define SEMAPHORE_NA 0	/* For commands, like removal, where semnum
			   argument is ignored */
#define SEMAPHORE_CREATION 0	/* -1 to lock, +1 to unlock */

/* See if item is already in shared memory */
static bool
shmem_exists_p (int *shmid, key_t key) {
  if ((*shmid = shmget(key,0,0)) == -1) {
    return false;
  } else {
    return true;
  }
}


static short
shmem_nattached (int shmid) {
  struct shmid_ds buf;

  if (shmctl(shmid,IPC_STAT,&buf) == -1) {
#if 0
    fprintf(stderr,"Error in shmem_nattached with shmctl.  Error %d: %s\n",
	    errno,strerror(errno));
#endif
    return -1;
  } else {
    return buf.shm_nattch;
  }
}


void
Access_shmem_remove (char *filename) {
  key_t key;
  int shmid, semid;
  struct shmid_ds *buf = NULL;

  key = ftok(filename,PROJECT_ID);
  if (shmem_exists_p(&shmid,key) == false) {
    /* Nothing to do */
  } else if (shmctl(shmid,IPC_RMID,buf) == -1) {
    fprintf(stderr,"Error with shmctl.  Error %d: %s\n",errno,strerror(errno));
  } else {
    fprintf(stderr,"Successfully removed existing memory\n");
  }

  if ((semid = semget(key,/*nsems*/0,0)) == -1) {
    /* Nothing to do */
  } else {
    fprintf(stderr,"Removing semaphore set %d\n",semid);
    semctl(semid,SEMAPHORE_NA,IPC_RMID,NULL);
  }

  return;
}

static void
semaphore_init (int semid, int sem_i, int value) {
  union semun {
    int val;
    struct semid_ds *buf;
    ushort *array;
  } argument;

  argument.val = value;
  semctl(semid,sem_i,SETVAL,argument);
  return;
}

#if 0
static int
semaphore_value (int semid, int sem_i) {
  return semctl(semid,sem_i,GETVAL,NULL);
}
#endif

#if 0
static int
semaphore_nwaiting (int semid, int sem_i) {
  printf("nwaiting = %d\n",semctl(semid,sem_i,GETNCNT,NULL));
  return semctl(semid,sem_i,GETNCNT,NULL);
}
#endif

#if 0
/* If already locked, then puts process to sleep */
static void
semaphore_lock (int semid) {
  struct sembuf op;

  /* printf("Process %d locking semaphore %d\n",getpid(),semid); */
  op.sem_num = SEMAPHORE_CREATION;
  op.sem_op = -1;
  op.sem_flg = SEM_UNDO;
  semop(semid,&op,1);

  return;
}
#endif

#if 0
static int
semaphore_unlock (int semid) {
  struct sembuf op;

  /* printf("Process %d unlocking semaphore %d\n",getpid(),semid); */
  op.sem_num = SEMAPHORE_CREATION;
  op.sem_op = +1;
  op.sem_flg = SEM_UNDO;
  semop(semid,&op,1);

  /* printf("%d processes still waiting\n",semctl(semid,SEMAPHORE_CREATION,GETNCNT,NULL)); */

  return semctl(semid,SEMAPHORE_CREATION,GETNCNT,NULL);
}
#endif


/* Same as semaphore_lock */
static void
semaphore_wait (int semid) {
  struct sembuf op;

  /* printf("Process %d locking semaphore %d\n",getpid(),semid); */
  op.sem_num = SEMAPHORE_CREATION;
  op.sem_op = -1;
  op.sem_flg = SEM_UNDO;
  semop(semid,&op,1);

  return;
}


void
Access_controlled_cleanup () {
  List_free(&shmem_memory);
  Intlist_free(&shmem_ids);
  Intlist_free(&semaphore_ids);
  return;
}


void
Access_emergency_cleanup () {
  List_T p;
  Intlist_T q;
  void *memory;
  int shmid, semid;
  int nattached;
  struct shmid_ds *buf = NULL;
  
  fprintf(stderr,"Calling Access_emergency_cleanup\n");
  for (p = shmem_memory, q = shmem_ids; p != NULL; p = List_next(p), q = Intlist_next(q)) {
    memory = List_head(p);
    shmid = Intlist_head(q);

    if (shmdt(memory) == -1) {
#if 0
      /* Somehow, shmdt forks and prints the error message and continues with the rest of the code */
      fprintf(stderr,"Error in Access_emergency_cleanup with shmdt on memory %p, shmid %d.  Error %d: %s\n",
	      memory,shmid,errno,strerror(errno));
#endif
    }

    if ((nattached = shmem_nattached(shmid)) > 0) {
      fprintf(stderr,"For shmid %d, %d other processes still attached\n",shmid,(int) nattached);

    } else if (shmctl(shmid,IPC_RMID,buf) == -1) {
#if 0
      fprintf(stderr,"Error in Access_emergency_cleanup with shmctl.  Error %d: %s\n",
	      errno,strerror(errno));
#endif
    } else {
      fprintf(stderr,"Removed existing memory for shmid %d\n",shmid);
    }
  }

  for (q = semaphore_ids; q != NULL; q = Intlist_next(q)) {
    /* Many of these removals will be for semaphores that don't exist */
    semid = Intlist_head(q);
    semctl(semid,/*semnum*/0,IPC_RMID,NULL);
  }
  Intlist_free(&semaphore_ids);

  if (shmem_memory != NULL) {
    fprintf(stderr,"\n");
    fprintf(stderr,"You may want to run 'ipcs -m' to see if any shared memory segments are still in use\n");
    fprintf(stderr,"You can remove a shared memory segment manually by doing 'ipcrm -m <shmid>'\n");
    fprintf(stderr,"\n");
    List_free(&shmem_memory);
    Intlist_free(&shmem_ids);
  }

  return;
}


void
Access_deallocate (void *memory, int shmid) {
  struct shmid_ds *buf = NULL;
  short nattached;

  if (shmdt(memory) == -1) {
#if 0      
      /* Somehow, shmdt forks and prints the error message and continues with the rest of the code */
      fprintf(stderr,"Error in Access_emergency_cleanup with shmdt on memory %p, shmid %d.  Error %d: %s\n",
	      memory,shmid,errno,strerror(errno));
#endif
  }

  if ((nattached = shmem_nattached(shmid)) > 0) {
    fprintf(stderr,"For shmid %d, %d processes still attached\n",shmid,(int) nattached);
  } else if (shmctl(shmid,IPC_RMID,buf) == -1) {
#if 0
    /* Somehow, shmctl forks and prints the error message and continues with the rest of the code */
    fprintf(stderr,"Error in Access_deallocate with shmctl.  Error %d: %s\n",
	    errno,strerror(errno));
#endif
  } else {
    fprintf(stderr,"Removed existing memory for shmid %d\n",shmid);
  }

  return;
}



#define FREAD_BATCH 100000000	/* 100 million elements at a time */

static void
copy_memory_from_file (void *memory, char *filename, size_t filesize, size_t eltsize) {
  FILE *fp;
  void *p;
  size_t i;

  if ((fp = FOPEN_READ_BINARY(filename)) == NULL) {
    fprintf(stderr,"Error: can't open file %s with fopen\n",filename);
    exit(9);
  }
  
  if (eltsize == 1) {
    for (i = 0; i + FREAD_BATCH < filesize/eltsize; i += FREAD_BATCH) {
      p = (void *) &(((unsigned char *) memory)[i]);
      fread(p,sizeof(unsigned char),FREAD_BATCH,fp);
    }
    
    if (i < filesize/eltsize) {
      p = (void *) &(((unsigned char *) memory)[i]);
      fread(p,sizeof(unsigned char),filesize/eltsize - i,fp);
    }

  } else if (eltsize == 4) {
    for (i = 0; i + FREAD_BATCH < filesize/eltsize; i += FREAD_BATCH) {
      p = (void *) &(((UINT4 *) memory)[i]);
      fread(p,sizeof(UINT4),FREAD_BATCH,fp);
    }

    if (i < filesize/eltsize) {
      p = (void *) &(((UINT4 *) memory)[i]);
      fread(p,sizeof(UINT4),filesize/eltsize - i,fp);
    }

  } else if (eltsize == 8) {
    for (i = 0; i + FREAD_BATCH < filesize/eltsize; i += FREAD_BATCH) {
      p = (void *) &(((UINT8 *) memory)[i]);
      fread(p,sizeof(UINT8),FREAD_BATCH,fp);
    }
    
    if (i < filesize/eltsize) {
      p = (void *) &(((UINT8 *) memory)[i]);
      fread(p,sizeof(UINT8),filesize/eltsize - i,fp);
    }

  } else {
    fprintf(stderr,"Access_allocated called with an element size of %d, which is not handled\n",(int) eltsize);
    exit(9);
  }
  fclose(fp);

  return;
}


static void *
shmem_attach (int *shmid, char *filename, off_t filesize, size_t eltsize) {
  void *memory = NULL;
  key_t key;
  int semid;

  key = ftok(filename,PROJECT_ID);
  if ((semid = semget(key,/*nsems*/1,IPC_CREAT | IPC_EXCL | 0666)) != -1) {
    /* Usually, we would set the value to be 1.  However, we can set
       the value to 0, because this process won't perform semaphore_wait */
    semaphore_init(semid,SEMAPHORE_CREATION,/*value*/0);

    /* Store semid in case we abort in the middle of this procedure */
    semaphore_ids = Intlist_push(semaphore_ids,semid);
  } else if ((semid = semget(key,0,0)) == -1) {
    fprintf(stderr,"Error in getting semaphore\n");
    abort();
  } else {
    semaphore_wait(semid);
  }

  /* The process tha created the semaphore will proceed, while the
     others wait.  They will be woken up when the semaphore is
     removed. */

  if ((*shmid = shmget(key,filesize,IPC_CREAT | IPC_EXCL |
#ifdef HAVE_SHM_NORESERVE
		       SHM_NORESERVE | 
#endif
		       0666)) != -1) {
    /* Created new shared memory */
    if ((memory = shmat(*shmid,NULL,0)) == (void *) -1) {
      fprintf(stderr,"Error with shmat.  Error %d: %s\n",errno,strerror(errno));
    } else {
      shmem_memory = List_push(shmem_memory,memory);
      shmem_ids = Intlist_push(shmem_ids,*shmid);
      copy_memory_from_file(memory,filename,filesize,eltsize);
      fprintf(stderr,"Attached new memory for %s...",filename);
    }

  } else if ((*shmid = shmget(key,0,0)) != -1) {
    /* Found existing shared memory */
    if ((memory = shmat(*shmid,NULL,0)) == (void *) -1) {
      fprintf(stderr,"Error with shmat.  Error %d: %s\n",errno,strerror(errno));
    } else {
      shmem_memory = List_push(shmem_memory,memory);
      shmem_ids = Intlist_push(shmem_ids,*shmid);
      fprintf(stderr,"Attached existing memory for %s...",filename);
    }
    
  } else {
    fprintf(stderr,"Using malloc instead of shmget for file %s\n",filename);
    memory = (void *) NULL;
  }

  /* The process that proceeded removes the semaphore here, allowing
     the other processes to continue after their waits.  The other
     processes will try to remove the semaphore too, yielding an
     error, which we simply ignore. */
  semctl(semid,SEMAPHORE_NA,IPC_RMID,NULL);

  return memory;
}


/* Bigendian conversion not needed after this */
void *
Access_allocate (int *shmid, size_t *len, double *seconds, char *filename, size_t eltsize, bool sharedp) {
  void *memory;
#ifdef CHECK
  void *memory2;
#endif
#if 0 && defined (USE_MPI)
  /* Does not work.  Gets ftruncate error */
  MPI_Comm comm;
  MPI_Win win;
#endif
  Stopwatch_T stopwatch;

  *len = (size_t) Access_filesize(filename);
  if (*len == 0) {
    *seconds = 0.0;
    return (void *) NULL;
  }

  Stopwatch_start(stopwatch = Stopwatch_new());

#ifdef CHECK
  memory2 = (void *) MALLOC(*len);
  if ((fp = FOPEN_READ_BINARY(filename)) == NULL) {
    fprintf(stderr,"Error: can't open file %s with fopen\n",filename);
    exit(9);
  }

  if (eltsize == 1) {
    FREAD_CHARS(memory2,(*len)/eltsize,fp);
  } else if (eltsize == 4) {
    FREAD_UINTS(memory2,(*len)/eltsize,fp);
  } else if (eltsize == 8) {
    FREAD_UINT8S(memory2,(*len)/eltsize,fp);
  } else {
    fprintf(stderr,"Access_allocated called with an element size of %d, which is not handled\n",(int) eltsize);
    exit(9);
  }
  fclose(fp);
#endif

  if (sharedp == true) {
#if 0 && defined(USE_MPI)
    /* Does not work.  Gives ftruncate error */
    MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,0,MPI_INFO_NULL,&comm);
    MPI_Win_allocate_shared(*len,/*disp_unit*/1,MPI_INFO_NULL,comm,&memory,&win);
    MPI_Win_free(&win);
#else
    if ((memory = shmem_attach(&(*shmid),filename,/*filesize*/*len,eltsize)) == NULL) {
      fprintf(stderr,"shm_attach not working on file %s, so using malloc instead on %lu bytes\n",
	      filename,*len);
      *shmid = 0;
      memory = (void *) MALLOC(*len);
      copy_memory_from_file(memory,filename,/*filesize*/*len,eltsize);
    }
#endif
  } else {
    *shmid = 0;
    memory = (void *) MALLOC(*len);
    copy_memory_from_file(memory,filename,/*filesize*/*len,eltsize);
  }

  /* Note: the following (old non-batch mode) requires conversion to bigendian later, as needed */
  /* fread(new->offsets,eltsize,sb.st_size/eltsize,fp); */

#ifdef CHECK
  for (i = 0; i < *len; i++) {
    if (((unsigned char *) memory)[i] != ((unsigned char *) memory2)[i]) {
      abort();
    }
  }
  FREE(memory2);
#endif

  *seconds = Stopwatch_stop(stopwatch);
  Stopwatch_free(&stopwatch);

  return memory;
}


#define PAGESIZE 1024*4

static int
get_pagesize () {

#ifdef PAGESIZE_VIA_SYSCTL
  int pagesize;
  size_t pagelen;
  int mib[2];
#endif

#ifdef __STRICT_ANSI__
  return PAGESIZE;
#elif defined(HAVE_GETPAGESIZE)
  return getpagesize();
#elif defined(PAGESIZE_VIA_SYSCONF)
  return (int) sysconf(_SC_PAGESIZE);
#elif defined(PAGESIZE_VIA_SYSCTL)
  pagelen = sizeof(pagesize);
  mib[0] = CTL_HW;
  mib[1] = HW_PAGESIZE;
  sysctl(mib,2,&pagesize,&pagelen,NULL,0);
  return pagesize;
#else
  return PAGESIZE;
#endif

}


#ifdef HAVE_MMAP
/* Returns NULL if mmap fails.  Bigendian conversion required */
#ifdef HAVE_CADDR_T
caddr_t
#else
void *
#endif
Access_mmap (int *fd, size_t *len, char *filename, size_t eltsize, bool randomp) {
  off_t length;
#ifdef HAVE_CADDR_T
  caddr_t memory;
#else
  void *memory;
#endif

  if ((*len = length = Access_filesize(filename)) == 0U) {
    fprintf(stderr,"Warning: file %s is empty\n",filename);
    *fd = open(filename,O_RDONLY,0764); /* Still need to initialize value */
    memory = (void *) NULL;

  } else if ((*fd = open(filename,O_RDONLY,0764)) < 0) {
    fprintf(stderr,"Error: can't open file %s with open for reading\n",filename);
    exit(9);

  } else if (sizeof(size_t) <= 4 && length > MAX32BIT) {
    debug(printf("Too big to mmap\n"));
    *len = 0;
    memory = (void *) NULL;

  } else {
    *len = (size_t) length;
    memory = mmap(NULL,length,PROT_READ,0
#ifdef HAVE_MMAP_MAP_SHARED
		  |MAP_SHARED
#endif
#ifdef HAVE_MMAP_MAP_FILE
		  |MAP_FILE
#endif
#ifdef HAVE_MMAP_MAP_VARIABLE
		  |MAP_VARIABLE
#endif
		  /*|MAP_NORESERVE*/
		  ,*fd,0);

    if (memory == MAP_FAILED) {
      fprintf(stderr,"Error in access.c (1): Got mmap failure on len %jd from length %jd.  Error %d: %s\n",
	      length,length,errno,strerror(errno));
      debug(printf("Got MAP_FAILED on len %jd from length %jd\n",length,length));
      memory = (void *) NULL;

    } else if (randomp == true) {
      debug(printf("Got mmap of %jd bytes at %p to %p\n",length,memory,memory+length-1));
#ifdef HAVE_MADVISE
#ifdef HAVE_MADVISE_MADV_RANDOM
      madvise(memory,*len,MADV_RANDOM);
#endif
#endif

    } else {
      debug(printf("Got mmap of %jd bytes at %p to %p\n",length,memory,memory+length-1));
#ifdef HAVE_MADVISE
#ifdef HAVE_MADVISE_MADV_DONTNEED
      madvise(memory,*len,MADV_DONTNEED);
#endif
#endif
    }
  }

  return memory;
}
#endif



#ifdef HAVE_MMAP
/* Returns NULL if mmap fails.  Bigendian conversion required */
#ifdef HAVE_CADDR_T
caddr_t
#else
void *
#endif
Access_mmap_offset (int *remainder, int fd, off_t offset, size_t length, size_t eltsize, bool randomp) {
#ifdef HAVE_CADDR_T
  caddr_t memory;
#else
  void *memory;
#endif

  if (length == 0) {
    abort();
  }

  *remainder = offset % get_pagesize();
  offset -= (off_t) *remainder;
  length += (size_t) *remainder;

  if (sizeof(size_t) <= 4 && length > MAX32BIT) {
    debug(printf("Too big to mmap\n"));
    memory = (void *) NULL;
  } else {
    memory = mmap(NULL,length,PROT_READ,0
#ifdef HAVE_MMAP_MAP_SHARED
		  |MAP_SHARED
#endif
#ifdef HAVE_MMAP_MAP_FILE
		  |MAP_FILE
#endif
#ifdef HAVE_MMAP_MAP_VARIABLE
		  |MAP_VARIABLE
#endif
		  /*|MAP_NORESERVE*/
		  ,fd,offset);

    if (memory == MAP_FAILED) {
      fprintf(stderr,"Error in access.c (2): Got mmap failure on fd %d, offset %jd, length %jd.  Error %d: %s\n",
	      fd,offset,length,errno,strerror(errno));
      debug(printf("Got MAP_FAILED on fd %d, offset %jd, length %zu\n",fd,offset,length));
      memory = (void *) NULL;

    } else if (randomp == true) {
      debug(printf("Got mmap of %jd bytes at %p to %p\n",length,memory,memory+length-1));
#ifdef HAVE_MADVISE
#ifdef HAVE_MADVISE_MADV_RANDOM
      madvise(memory,length,MADV_RANDOM);
#endif
#endif

    } else {
      debug(printf("Got mmap of %jd bytes at %p to %p\n",length,memory,memory+length-1));
#ifdef HAVE_MADVISE
#ifdef HAVE_MADVISE_MADV_DONTNEED
      madvise(memory,length,MADV_DONTNEED);
#endif
#endif
    }
  }

  return memory;
}
#endif



#ifdef HAVE_MMAP
/* Returns NULL if mmap fails.  Bigendian conversion required */
#ifdef HAVE_CADDR_T
caddr_t
#else
void *
#endif
Access_mmap_rw (int *fd, size_t *len, char *filename, size_t eltsize, bool randomp) {
  off_t length;
#ifdef HAVE_CADDR_T
  caddr_t memory;
#else
  void *memory;
#endif

  if ((*len = length = Access_filesize(filename)) == 0U) {
    fprintf(stderr,"Warning: file %s is empty\n",filename);
    *fd = open(filename,O_RDWR,0764); /* Still need to initialize value */
    memory = (void *) NULL;
  } else if ((*fd = open(filename,O_RDWR,0764)) < 0) {
    fprintf(stderr,"Error: can't open file %s with open for reading/writing\n",filename);
    exit(9);
  } else if (sizeof(size_t) <= 4 && length > MAX32BIT) {
    debug(printf("Too big to mmap\n"));
    *len = 0;
    memory = (void *) NULL;
  } else {
    *len = (size_t) length;
    memory = mmap(NULL,length,PROT_READ|PROT_WRITE,0
#ifdef HAVE_MMAP_MAP_SHARED
		  |MAP_SHARED
#endif
#ifdef HAVE_MMAP_MAP_FILE
		  |MAP_FILE
#endif
#ifdef HAVE_MMAP_MAP_VARIABLE
		  |MAP_VARIABLE
#endif
		  /*|MAP_NORESERVE*/
		  ,*fd,0);

    if (memory == MAP_FAILED) {
      fprintf(stderr,"Error in access.c (3): Got mmap failure on len %jd from length %jd.  Error %d: %s\n",
	      *len,length,errno,strerror(errno));
      debug(printf("Got MAP_FAILED on len %zu from length %jd\n",*len,length));
      memory = (void *) NULL;

    } else if (randomp == true) {
      debug(printf("Got mmap of %jd bytes at %p to %p\n",length,memory,memory+length-1));
#ifdef HAVE_MADVISE
#ifdef HAVE_MADVISE_MADV_RANDOM
      madvise(memory,*len,MADV_RANDOM);
#endif
#endif

    } else {
      debug(printf("Got mmap of %jd bytes at %p to %p\n",length,memory,memory+length-1));
#ifdef HAVE_MADVISE
#ifdef HAVE_MADVISE_MADV_DONTNEED
      madvise(memory,*len,MADV_DONTNEED);
#endif
#endif
    }
  }

  return memory;
}
#endif

#ifdef HAVE_MMAP
/* Returns NULL if mmap fails.  Bigendian conversion required */
#ifdef HAVE_CADDR_T
caddr_t
#else
void *
#endif
Access_mmap_offset_rw (int *remainder, int fd, off_t offset, size_t length, size_t eltsize, bool randomp) {
#ifdef HAVE_CADDR_T
  caddr_t memory;
#else
  void *memory;
#endif

  if (length == 0) {
    abort();
  }

  *remainder = offset % get_pagesize();
  offset -= (off_t) *remainder;
  length += (size_t) *remainder;

  if (sizeof(size_t) <= 4 && length > MAX32BIT) {
    debug(printf("Too big to mmap\n"));
    memory = (void *) NULL;
  } else {
    memory = mmap(NULL,length,PROT_READ|PROT_WRITE,0
#ifdef HAVE_MMAP_MAP_SHARED
		  |MAP_SHARED
#endif
#ifdef HAVE_MMAP_MAP_FILE
		  |MAP_FILE
#endif
#ifdef HAVE_MMAP_MAP_VARIABLE
		  |MAP_VARIABLE
#endif
		  /*|MAP_NORESERVE*/
		  ,fd,offset);

    if (memory == MAP_FAILED) {
      fprintf(stderr,"Error in access.c (4): Got mmap failure on offset %jd, length %jd.  Error %d: %s\n",
	      offset,length,errno,strerror(errno));
      debug(printf("Got MAP_FAILED on offset %jd, length %zu\n",offset,length));
      memory = (void *) NULL;

    } else if (randomp == true) {
      debug(printf("Got mmap of %zu bytes at %p to %p\n",length,memory,memory+length-1));
#ifdef HAVE_MADVISE
#ifdef HAVE_MADVISE_MADV_RANDOM
      madvise(memory,length,MADV_RANDOM);
#endif
#endif

    } else {
      debug(printf("Got mmap of %zu bytes at %p to %p\n",length,memory,memory+length-1));
#ifdef HAVE_MADVISE
#ifdef HAVE_MADVISE_MADV_DONTNEED
      madvise(memory,length,MADV_DONTNEED);
#endif
#endif
    }
  }

  return memory;
}
#endif




#ifdef HAVE_MMAP

#ifdef HAVE_CADDR_T
caddr_t
#else
void *
#endif
Access_mmap_and_preload (int *fd, size_t *len, int *npages, double *seconds, char *filename, size_t eltsize) {
  off_t length;
#ifdef HAVE_CADDR_T
  caddr_t memory;
#else
  void *memory;
#endif
  int pagesize, indicesperpage;
  size_t totalindices, i;	/* Needs to handle uncompressed genomes > 2 gigabytes */
  int nzero = 0, npos = 0;
  Stopwatch_T stopwatch;


  if ((*len = length = Access_filesize(filename)) == 0U) {
    fprintf(stderr,"Warning: file %s is empty\n",filename);
    *fd = open(filename,O_RDONLY,0764); /* Still need to initialize value */
    memory = (void *) NULL;

  } else if ((*fd = open(filename,O_RDONLY,0764)) < 0) {
    fprintf(stderr,"Error: can't open file %s with open for reading\n",filename);
    exit(9);

  } else if (sizeof(size_t) <= 4 && *len > MAX32BIT) {
    debug(printf("Too big to mmap\n"));
    *len = 0;
    *npages = 0;
    *seconds = 0.0;
    memory = (void *) NULL;

  } else {
    pagesize = get_pagesize();

    indicesperpage = pagesize/eltsize;
    
    Stopwatch_start(stopwatch = Stopwatch_new());

    memory = mmap(NULL,length,PROT_READ,0
#ifdef HAVE_MMAP_MAP_SHARED
		  |MAP_SHARED
#endif
#ifdef HAVE_MMAP_MAP_FILE
		  |MAP_FILE
#endif
#ifdef HAVE_MMAP_MAP_VARIABLE
		  |MAP_VARIABLE
#endif
		  /*|MAP_NORESERVE*/
		  ,*fd,0);

    if (memory == MAP_FAILED) {
      fprintf(stderr,"Error in access.c (5): Got mmap failure on len %jd from length %jd.  Error %d: %s\n",
	      *len,length,errno,strerror(errno));
      debug(printf("Got MAP_FAILED on len %jd from length %zu\n",*len,length));
      memory = (void *) NULL;
      Stopwatch_stop(stopwatch);
      Stopwatch_free(&stopwatch);

    } else {
      /* Touch all pages */
      debug(printf("Got mmap of %zu bytes at %p to %p\n",length,memory,memory+length-1));
#ifdef HAVE_MADVISE
#ifdef HAVE_MADVISE_MADV_WILLNEED
      madvise(memory,*len,MADV_WILLNEED);
#endif
#endif
      totalindices = (*len)/eltsize;
      for (i = 0; i < totalindices; i += indicesperpage) {
	if (((char *) memory)[i] == 0) {
	  nzero++;
	  if (i % 10000 == 0) {
	    fprintf(stderr,",");
	  }
	} else {
	  npos++;
	}
	if (i % 10000 == 0) {
	  fprintf(stderr,".");
	}
      }
      *npages = nzero + npos;
      *seconds = Stopwatch_stop(stopwatch);
      Stopwatch_free(&stopwatch);
    }
  }

  return memory;
}
#endif


