/* $Id: psp_debug.c,v 1.8 2005-01-15 16:31:56 paklein Exp $ */
#include <stdio.h>
#include "pspases_f2c.h"
#include "mpi.h"

/* BLAS routines */
void dsyrk_(char *UL, char *NT, integer *N, integer *K,
	doublereal *alpha, doublereal *A, integer *lda, 
	doublereal *beta, doublereal *C, integer *ldc, 
	ftnlen UL_len, ftnlen NT_len);

/* MPI routines */
int MPI_Irecv_d(void* buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request);
int MPI_Recv_d(void* buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status);

int MPI_Isend_d(void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request);
int MPI_Send_d(void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);

int MPI_Get_count_d(MPI_Status *status, MPI_Datatype datatype, int *count);
int MPI_Waitall_d(int count, MPI_Request* request, MPI_Status* status);
int MPI_Waitany_d(int count, MPI_Request* request, int* index, MPI_Status* status);
int MPI_Wait_d(MPI_Request* request, MPI_Status* status);

static const char* type_names[] = {"UNKNOWN", "MPI_BYTE", "MPI_INT", "MPIT_DOUBLE"};
const char* t2s(MPI_Datatype datatype) {

	if (datatype == MPI_BYTE)
		return type_names[1];
	else if (datatype == MPI_INT)
		return type_names[2];
	else if (datatype == MPI_DOUBLE)
		return type_names[3];
	else
		return type_names[0];
};
 
void dsyrk_d(char *UL, char *NT, integer *N, integer *K,
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

int MPI_Isend_d(void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
{
	int ret;

	/* report values */
	printf("MPI_Isend: count = %d, datatype = %s, dest = %d, tag = %d\n", 
		count, t2s(datatype), dest, tag);
	fflush(stdout);
	
	/* call */
	ret = MPI_Isend(buf, count, datatype, dest, tag, comm, request);

	printf("MPI_Isend: request = %x\n", *request);
	
	return ret;
}

int MPI_Send_d(void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
{
	/* report values */
	printf("MPI_Send: count = %d, datatype = %s, dest = %d, tag = %d\n", 
		count, t2s(datatype), dest, tag);
	fflush(stdout);
	
	/* call */
	return MPI_Send(buf, count, datatype, dest, tag, comm);
}

int MPI_Recv_d(void* buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status)
{
	/* report values */
	printf("MPI_Send: count = %d, datatype = %s, source = %d, tag = %d\n", 
		count, t2s(datatype), source, tag);
	fflush(stdout);
	
	/* call */
	return MPI_Recv(buf, count, datatype, source, tag, comm, status);
}

int MPI_Irecv_d(void* buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request)
{
	int ret;

	/* report values */
	printf("MPI_Irecv: count = %d, datatype = %s, dest = %d, tag = %d, request = %x\n", 
		count, t2s(datatype), source, tag, *request);
	fflush(stdout);
	
	/* call */
	ret = MPI_Irecv(buf, count, datatype, source, tag, comm, request);

	return ret;
}

int MPI_Get_count_d(MPI_Status *status, MPI_Datatype datatype, int *count)
{
	int ret;

	/* report */
	printf("MPI_Get_count: status = %x, source = %d, tag = %d, datatype = %s",
		status, status->MPI_SOURCE, status->MPI_TAG, t2s(datatype));
	fflush(stdout);

	/* call */
	ret = MPI_Get_count(status, datatype, count);

	/* result */
	printf(", count = %d\n", *count);
	fflush(stdout);

	return ret;
}

int MPI_Waitall_d(int count, MPI_Request* request, MPI_Status* status)
{
	int ret, i;

	printf("MPI_Waitall: count = %d\n", count);
	for (i = 0; i < count; i++)
		if (request[i] == MPI_REQUEST_NULL)
			printf("MPI_Waitall: request[%d] = NULL\n", i);
		else
			printf("MPI_Waitall: request[%d] = %x\n", i, request[i]);

	/* call */	
	ret = MPI_Waitall(count, request, status);

	for (i = 0; i < count; i++)
		printf("MPI_Waitall: status[%d]: status = %x, source = %d, tag = %d\n", i,
			status+i, status[i].MPI_SOURCE, status[i].MPI_TAG);
	
	return ret;
}

int MPI_Waitany_d(int count, MPI_Request *request, int *index, MPI_Status *status)
{
	int ret, i;

	printf("MPI_Waitany: count = %d\n", count);
	for (i = 0; i < count; i++)
		if (request[i] == MPI_REQUEST_NULL)
			printf("MPI_Waitany: request[%d] = NULL\n", i);
		else
			printf("MPI_Waitany: request[%d] = %x\n", i, request[i]);

	/* call */	
	ret = MPI_Waitany(count, request, index, status);

	printf("MPI_Waitany: status[%d]: status = %x, source = %d, tag = %d\n", *index,
		status + *index, status[*index].MPI_SOURCE, status[*index].MPI_TAG);
}
 
int MPI_Wait_d(MPI_Request* request, MPI_Status* status)
{
	int ret, i;

	printf("MPI_Wait: request = %x\n", *request);

	/* call */	
	ret = MPI_Wait(request, status);

	printf("MPI_Wait: source = %d, tag = %d\n", status->MPI_SOURCE, status->MPI_TAG);
	
	return ret;
}
