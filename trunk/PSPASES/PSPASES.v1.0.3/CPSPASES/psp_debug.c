/* $Id: psp_debug.c,v 1.1 2005-01-05 16:51:31 paklein Exp $ */
#include <stdio.h>
#include "pspases_f2c.h"
#include "mpi.h"

/* BLAS routines */
void dsyrk_(char *UL, char *NT, integer *N, integer *K,
	doublereal *alpha, doublereal *A, integer *lda, 
	doublereal *beta, doublereal *C, integer *ldc, 
	ftnlen UL_len, ftnlen NT_len);

/* MPI routines */
int MPI_Isend(void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request);
int MPI_Get_count(MPI_Status *status, MPI_Datatype datatype, int *count);
 
void mydsyrk_(char *UL, char *NT, integer *N, integer *K,
	doublereal *alpha, doublereal *A, integer *lda, 
	doublereal *beta, doublereal *C, integer *ldc, 
	ftnlen UL_len, ftnlen NT_len)
{
	int i;

	/* report values */
	printf("dsyrk:");
	printf(" UT = \"");
	for (i = 0; i < UL_len; i++)
		printf("%c", UL[i]);
	printf("\"");
	printf(", N = %d, K = %d, lda = %d, ldc = %d", *N, *K, *lda, *ldc);
	printf(", NT = \"");
	for (i = 0; i < NT_len; i++)
		printf("%c", NT[i]);
	printf("\"\n");
	fflush(stdout);

	/* call */
	dsyrk_(UL, NT, N, K, alpha, A, lda, beta, C, ldc, UL_len, NT_len);
}

int myMPI_Isend(void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
{
	/* report values */
	printf("MPI_Isend: count = %d, datatype = %d, dest = %d, tag = %d\n", 
		count, datatype, dest, tag);
	fflush(stdout);
	
	/* call */
	return MPI_Isend(buf, count, datatype, dest, tag, comm, request);
}

int myMPI_Get_count(MPI_Status *status, MPI_Datatype datatype, int *count)
{
	/* report */
	printf("MPI_Get_count: status = %x, source = %d, tag = %d, datatype = %d\n",
		status, status->MPI_SOURCE, status->MPI_TAG, datatype);
	fflush(stdout);

	/* call */
	return MPI_Get_count(status, datatype, count);
}
