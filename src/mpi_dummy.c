/** \file mpi_dummy.c

Dummy MPI functions library for sequential built.
It is used instead of the true MPI library in the
case of a sequential built (see makefile).

*/

#include <stdio.h>
#include "mpi_dummy.h"

#pragma GCC diagnostic ignored "-Wmissing-prototypes"

void MPI_Comm_rank (int a, int *b) {*b = 0;} /* Only one process, with rank zero... */

void MPI_Comm_size (int a, int *b) {*b = 1;} /* Only one process in the world communicator... */

void MPI_Init (int *argc, char **argv[]) {
  fprintf (stderr, "\n       !!!! WARNING !!!!\n\n");
  fprintf (stderr, "This is a sequential built of the %s code\n", *argv[0]);
  fprintf (stderr, "If you planned to run the MPI-parallel version,\n");
  fprintf (stderr, "then you MUST rebuild the executable. In this case,\n");
  fprintf (stderr, " issue:\nmake PARALLEL=1 (or make para)\n");
  fprintf (stderr, "\nAny subsequent invocation of make will build an MPI version.\n");
}

void MPI_Allreduce (void *ptr, void *ptr2, int count, int type, int foo3, int foo4) {
  int i;
  for (i = 0; i < count; i++) {
    switch (type) {
    case MPI_FLOAT:
      *(((float *)ptr2)+i) = (float)(*(((float *)ptr)+i));
      break;
    case MPI_DOUBLE:
      *(((double *)ptr2)+i) = (double)(*(((double *)ptr)+i));
      break;
    case MPI_INT:
      *(((int *)ptr2)+i) = (int)(*(((int *)ptr)+i));
      break;
    }
  }
}

void MPI_Reduce (void *ptr, void *ptr2, int count, int type, int foo3, int foo4, int foo5) {
  int i;
  for (i = 0; i < count; i++) {
    switch (type) {
    case MPI_FLOAT:
      *(((float *)ptr2)+i) = (float)(*(((float *)ptr)+i));
      break;
    case MPI_DOUBLE:
      *(((double *)ptr2)+i) = (double)(*(((double *)ptr)+i));
      break;
    case MPI_INT:
      *(((int *)ptr2)+i) = (int)(*(((int *)ptr)+i));
      break;
    }
  }
}
void MPI_Allgather(void* ptr, int count, int type, void* ptr2, int f00, int f01, int f02) {
  int i;
  for (i = 0; i < count; i++) {
    switch (type) {
    case MPI_FLOAT:
      *(((float *)ptr2)+i) = (float)(*(((float *)ptr)+i));
      break;
    case MPI_DOUBLE:
      *(((double *)ptr2)+i) = (double)(*(((double *)ptr)+i));
      break;
    case MPI_INT:
      *(((int *)ptr2)+i) = (int)(*(((int *)ptr)+i));
      break;
    case MPI_CHAR:
      *(((char *)ptr2)+i) = (int)(*(((int *)ptr)+i));
      break;
    }
  }
}
void MPI_Allgatherv(void* ptr, int count, int type, void* ptr2, int f00, int f01, int f02, int f03) {
  int i;
  for (i = 0; i < count; i++) {
    switch (type) {
    case MPI_FLOAT:
      *(((float *)ptr2)+i) = (float)(*(((float *)ptr)+i));
      break;
    case MPI_DOUBLE:
      *(((double *)ptr2)+i) = (double)(*(((double *)ptr)+i));
      break;
    case MPI_INT:
      *(((int *)ptr2)+i) = (int)(*(((int *)ptr)+i));
      break;
    case MPI_CHAR:
      *(((char *)ptr2)+i) = (int)(*(((int *)ptr)+i));
      break;
    }
  }
}

//do nothing...
void MPI_Finalize(){}
void MPI_Bcast(){}
void MPI_Isend(){}
void MPI_Irecv(){}
void MPI_Send(){}
void MPI_Recv(){}
void MPI_Barrier(){}
void MPI_Wait(){}
void MPI_Scan(){} //In place scans require no special action
void MPI_Comm_split(){}
