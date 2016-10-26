/* $Id: mpidebug.h 162091 2015-03-26 18:30:04Z twu $ */
#ifndef MPIDEBUG_INCLUDED
#define MPIDEBUG_INCLUDED
#include <mpi.h>


#define MPI_TAG_DEFAULT 0
#define MPI_TAG_WANT_INPUT 1
#define MPI_TAG_GIVE_INPUT 2
#define MPI_TAG_WRITE_STDOUT 3



/* #define MPIDEBUG 1 */
/* #define MPIDEBUG_I 1 */

extern void
MPI_Debug_setup (int myid_in);

extern MPI_File
MPI_fopen (char *filename, MPI_Comm comm);

extern int
MPI_Debug_Send (const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, const char *file, int line);
extern int
MPI_Debug_Recv (void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status, const char *file, int line);

#ifdef MPIDEBUG
#define MPI_SEND(buf,count,datatype,dest,tag,comm) MPI_Debug_Send(buf,count,datatype,dest,tag,comm,__FILE__,__LINE__)
#define MPI_RECV(buf,count,datatype,source,tag,comm,status) MPI_Debug_Recv(buf,count,datatype,source,tag,comm,status,__FILE__,__LINE__)

#else
#define MPI_SEND(buf,count,datatype,dest,tag,comm) MPI_Send(buf,count,datatype,dest,tag,comm)
#define MPI_RECV(buf,count,datatype,source,tag,comm,status) MPI_Recv(buf,count,datatype,source,tag,comm,status)
#endif


#ifdef MPIDEBUG_I
#define MPI_ISEND(buf,count,datatype,dest,tag,comm,req) MPI_Debug_Isend(buf,count,datatype,dest,tag,comm,req,__FILE__,__LINE__)
#define MPI_IRECV(buf,count,datatype,source,tag,comm,req) MPI_Debug_Irecv(buf,count,datatype,source,tag,comm,req,__FILE__,__LINE__)

#else
#define MPI_ISEND(buf,count,datatype,dest,tag,comm,req) MPI_Isend(buf,count,datatype,dest,tag,comm,req)
#define MPI_IRECV(buf,count,datatype,source,tag,comm,req) MPI_Irecv(buf,count,datatype,source,tag,comm,req)
#endif

#define MPI_SSEND(buf,count,datatype,dest,tag,comm) MPI_Ssend(buf,count,datatype,dest,tag,comm)

#endif


