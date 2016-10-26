static char rcsid[] = "$Id: mpidebug.c 162090 2015-03-26 18:29:53Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "mpidebug.h"

#include <stdio.h>
#include "bool.h"
#include "types.h"
#include "genomicpos.h"

static int myid;

void
MPI_Debug_setup (int myid_in) {
  myid = myid_in;
  return;
}


MPI_File
MPI_fopen (char *filename, MPI_Comm comm) {
  MPI_File mpi_input;
#if 0
  MPI_Datatype contig, filetype;
#endif

  MPI_File_open(comm,filename,MPI_MODE_RDONLY,MPI_INFO_NULL,&mpi_input);
#if 0
  MPI_Type_contiguous(endptr - startptr,MPI_BYTE,&contig);
  MPI_Type_create_resized(contig,/*lowerbound*/0,/*extent*/filesize - startptr,&filetype);
  MPI_Type_commit(&filetype);
  MPI_File_set_view(mpi_input,/*disp*/startptr,MPI_BYTE,filetype,"native",MPI_INFO_NULL);
#endif

  return mpi_input;
}



static char *
get_typename (MPI_Datatype datatype) {
  if (datatype == MPI_INT) {
    return "MPI_INT";
  } else if (datatype == MPI_BOOL_T) {
    return "MPI_BOOL_T";
  } else if (datatype == MPI_CHAR) {
    return "MPI_CHAR";
  } else if (datatype == MPI_DOUBLE) {
    return "MPI_DOUBLE";
  } else if (datatype == MPI_FLOAT) {
    return "MPI_FLOAT";
  } else {
    return "OTHER";
  }
}

static void
print_value (const void *buf, int count, MPI_Datatype datatype) {
  if (datatype == MPI_INT) {
    printf("%d",* (int *) buf);
  } else if (datatype == MPI_BOOL_T) {
    printf("%d",(int) (* (bool *) buf));
  } else if (datatype == MPI_CHAR) {
    printf("%s",(char *) buf);
  } else if (datatype == MPI_DOUBLE) {
    printf("%f",* (double *) buf);
  } else if (datatype == MPI_FLOAT) {
    printf("%f",* (float *) buf);
  } else {
    printf("??");
  }
}


int
MPI_Debug_Send (const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, const char *file, int line) {
  printf("MPI_Send (%d->%d) at %s:%d: proc %d sending count %d of datatype %s and tag %d to dest %d: ",
	 myid,dest,file,line,myid,count,get_typename(datatype),tag,dest);
  print_value(buf,count,datatype);
  printf("\n");
  return MPI_Send(buf,count,datatype,dest,tag,comm);
}

int
MPI_Debug_Recv (void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status,
		const char *file, int line) {
  int result;
  result = MPI_Recv(buf,count,datatype,source,tag,comm,&(*status));
  printf("MPI_Recv (%d<-%d) at %s:%d: proc %d receiving count %d of datatype %s and tag %d from source %d: ",
	 myid,(*status).MPI_SOURCE,file,line,myid,count,get_typename(datatype),tag,(*status).MPI_SOURCE);
  print_value(buf,count,datatype);
  printf("\n");
  return result;
}

int
MPI_Debug_Isend (const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
		 MPI_Request *req, const char *file, int line) {
  MPI_Status status;
  int result;

  printf("MPI_Isend (%d->%d) at %s:%d: proc %d sending count %d of datatype %s to dest %d: ",
	 myid,dest,file,line,myid,count,get_typename(datatype),dest);

  print_value(buf,count,datatype);
  printf("\n");

  result = MPI_Isend(buf,count,datatype,dest,tag,comm,&(*req));
  /* MPI_Wait(&(*req),&status); */
  return result;
}

int
MPI_Debug_Irecv (void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm,
		 MPI_Request *req, const char *file, int line) {
  MPI_Status status;
  int result;

  result = MPI_Irecv(buf,count,datatype,source,tag,comm,&(*req));
  MPI_Wait(&(*req),&status);
  printf("MPI_Irecv (%d<-%d) at %s:%d: proc %d receiving count %d of datatype %s from source %d: ",
	 myid,status.MPI_SOURCE,file,line,myid,count,get_typename(datatype),status.MPI_SOURCE);

  print_value(buf,count,datatype);
  printf("\n");
  return result;
}
