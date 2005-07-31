/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile: md_wrap_c.c,v $
 *
 * $Author: paklein $
 *
 * $Date: 2001-01-30 20:59:12 $
 *
 * $Revision: 1.1.1.1 $
 *
 * $Name: not supported by cvs2svn $
 *====================================================================*/
#ifndef lint
static char *cvs_wrapmpi_id =
  "$Id: md_wrap_c.c,v 1.1.1.1 2001-01-30 20:59:12 paklein Exp $";
#endif


/*******************************************************************************
 * Copyright 1995, Sandia Corporation.  The United States Government retains a *
 * nonexclusive license in this software as prescribed in AL 88-1 and AL 91-7. *
 * Export of this program may require a license from the United States         *
 * Government.                                                                 *
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Environment.h"

/* disable unused arguments warning */
#ifdef __MWERKS__
#pragma warn_unusedarg  off
#endif /* __MWERKS__ */           

/*******************************************************************************/
#ifndef __MPI__
/*******************************************************************************/

#define MPI_Request int

/* added to work around using MPI_ANY_SOURCE for Irecv */
int md_max(int my_value, int* max_value);
int md_allgather(char *send_buf, char *recv_buf, int bytes);
/* PAK (10/07/2000) */

/* function prototypes */
void get_parallel_info(int *proc, int *nprocs, int *dim);
int md_read(char *buf, int bytes, int *source, int *type, int *flag);
int md_write(char *buf, int bytes, int dest, int type, int *flag);
int md_wrap_iread(void *buf, int bytes, int *source, int *type,
                  MPI_Request *request);
int md_wrap_write(void *buf, int bytes, int dest, int type, int *flag);   
int md_wrap_wait(void *buf, int bytes, int *source, int *type, int *flag,
                 MPI_Request *request);

/* added to work around using MPI_ANY_SOURCE for Irecv */
int md_max(int my_value, int* max_value)
{
	*max_value = my_value;
	return 0;
}

int md_allgather(char *send_buf, char *recv_buf, int bytes)
{
	memcpy(recv_buf, send_buf, bytes);
	return 0;
}
/* PAK (10/07/2000) */
                 
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void get_parallel_info(int *proc, int *nprocs, int *dim)

{
  *proc   = 0;
  *nprocs = 1;
  *dim    = 0;

} /* get_parallel_info */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int md_read(char *buf, int bytes, int *source, int *type, int *flag)

{
  return bytes;

} /* md_read */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int md_write(char *buf, int bytes, int dest, int type, int *flag)

{
  return 0;

} /* md_write */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int md_wrap_iread(void *buf, int bytes, int *source, int *type,
                  MPI_Request *request)

/*******************************************************************************

  Machine dependent wrapped message-reading communication routine for the
  Intel.  This routine is a simple no-op but is used in order to provide
  compatibility with the MPI communication routine order.

  Author:          Scott A. Hutchinson, SNL, 9221
  =======

  Return code:     int
  ============

  Parameter list:
  ===============

  buf:             Beginning address of data to be sent.

  bytes:           Length of message in bytes.

  source:          Source processor number.

  type:            Message type

*******************************************************************************/

{

  return 0;

} /* md_wrap_iread */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int md_wrap_write(void *buf, int bytes, int dest, int type, int *flag)

/*******************************************************************************

  Machine dependent wrapped message-sending communication routine for the
  Intel.  This routine is exactly the same as md_write.

  Author:          Scott A. Hutchinson, SNL, 9221
  =======

  Return code:     int
  ============

  Parameter list:
  ===============

  buf:             Beginning address of data to be sent.

  bytes:           Length of message in bytes.

  dest:            Destination processor number.

  type:            Message type

  flag:

*******************************************************************************/

{

  return 0;

} /* md_wrap_write */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int md_wrap_wait(void *buf, int bytes, int *source, int *type, int *flag,
                 MPI_Request *request)

/*******************************************************************************

  Machine dependent wrapped message-wait communication routine for the Intel.
  This routine is identical to md_read but is put here in order to be compatible
  with the order required to do MPI communication.

  Author:          Scott A. Hutchinson, SNL, 9221
  =======

  Return code:     int
  ============

  Parameter list:
  ===============

  buf:             Beginning address of data to be sent.

  bytes:           Length of message in bytes.
  dest:            Destination processor number.

  type:            Message type

  flag:

*******************************************************************************/

{

  return bytes;

} /* md_wrap_wait */

/*******************************************************************************/
#else /* has MPI */
/*******************************************************************************/

#include "mpi.h"

/* added to work around using MPI_ANY_SOURCE for Irecv */
int md_max(int my_value, int* max_value);
int md_allgather(char *send_buf, char *recv_buf, int bytes);
/* PAK (10/07/2000) */

/* function prototypes */
void get_parallel_info(int *proc, int *nprocs, int *dim);
int md_read(char *buf, int bytes, int *source, int *type, int *flag);
int md_write(char *buf, int bytes, int dest, int type, int *flag);
int md_wrap_iread(void *buf, int bytes, int *source, int *type,
                  MPI_Request *request);
int md_wrap_write(void *buf, int bytes, int dest, int type, int *flag);   
int md_wrap_wait(void *buf, int bytes, int *source, int *type, int *flag,
                 MPI_Request *request);

int gl_rbuf = 3;
int gl_sbuf = 3;

/* added to work around using MPI_ANY_SOURCE for Irecv */
int md_max(int my_value, int* max_value)
{
	int err;
	err = MPI_Allreduce(&my_value, max_value, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	return err;
}

int md_allgather(char *send_buf, char *recv_buf, int bytes)
{
	int err;
	err = MPI_Allgather(send_buf, bytes, MPI_BYTE, recv_buf, bytes, MPI_BYTE, MPI_COMM_WORLD);
	return err;
}
/* PAK (10/07/2000) */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void get_parallel_info(int *proc, int *nprocs, int *dim)
{

  /* local variables */

  int i;

  MPI_Comm_size(MPI_COMM_WORLD, nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, proc);
  *dim = 0;

} /* get_parallel_info */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int md_read(char *buf, int bytes, int *source, int *type, int *flag)

{

  int        err, buffer = 1;
  MPI_Status status;

  if (*type   == -1) *type   = MPI_ANY_TAG;
  if (*source == -1) *source = MPI_ANY_SOURCE;

  if (bytes == 0) {
    err = MPI_Recv(&gl_rbuf, 1, MPI_BYTE, *source, *type, MPI_COMM_WORLD,
                   &status);
  }
  else {
    err = MPI_Recv(buf, bytes, MPI_BYTE, *source, *type, MPI_COMM_WORLD,
                   &status);
  }

  if (err != 0) (void) fprintf(stderr, "MPI_Recv error = %d\n", err);
  MPI_Get_count(&status,MPI_BYTE,&buffer);
  *source = status.MPI_SOURCE;
  *type   = status.MPI_TAG;
  if (bytes != 0) bytes = buffer;

  return bytes;

} /* md_read */


/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int md_write(char *buf, int bytes, int dest, int type, int *flag)

{

  int err;

  if (bytes == 0) {
    err = MPI_Send(&gl_sbuf, 1, MPI_BYTE, dest, type, MPI_COMM_WORLD);
  }
  else {
    err = MPI_Send(buf, bytes, MPI_BYTE, dest, type, MPI_COMM_WORLD);
  }

  if (err != 0) (void) fprintf(stderr, "MPI_Send error = %d\n", err);

  return 0;

} /* md_write */



/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int md_wrap_iread(void *buf, int bytes, int *source, int *type,
                  MPI_Request *request)


/*******************************************************************************

  Machine dependent wrapped message-reading communication routine for MPI.

  Author:          Scott A. Hutchinson, SNL, 9221
  =======

  Return code:     int
  ============

  Parameter list:
  ===============

  buf:             Beginning address of data to be sent.

  bytes:           Length of message in bytes.

  source:          Source processor number.

  type:            Message type

*******************************************************************************/

{

  int err = 0;

  if (*type   == -1) *type   = MPI_ANY_TAG;
  if (*source == -1) *source = MPI_ANY_SOURCE;

  if (bytes == 0) {
    err = MPI_Irecv(&gl_rbuf, 1, MPI_BYTE, *source, *type, MPI_COMM_WORLD,
                    request);
  }
  else {
    err = MPI_Irecv(buf, bytes, MPI_BYTE, *source, *type, MPI_COMM_WORLD,
                    request);
  }

  return err;

} /* md_wrap_iread */


/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int md_wrap_write(void *buf, int bytes, int dest, int type, int *flag)

/*******************************************************************************

  Machine dependent wrapped message-sending communication routine for MPI.

  Author:          Scott A. Hutchinson, SNL, 9221
  =======

  Return code:     int
  ============

  Parameter list:
  ===============

  buf:             Beginning address of data to be sent.

  bytes:           Length of message in bytes.

  dest:            Destination processor number.

  type:            Message type

  flag:

*******************************************************************************/

{

  int err = 0;

  if (bytes == 0) {
    err = MPI_Send(&gl_sbuf, 1, MPI_BYTE, dest, type, MPI_COMM_WORLD);
  }
  else {
    err = MPI_Send(buf, bytes, MPI_BYTE, dest, type, MPI_COMM_WORLD);
  }

  return err;

} /* md_wrap_write */



/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int md_wrap_wait(void *buf, int bytes, int *source, int *type, int *flag,
                 MPI_Request *request)

/*******************************************************************************

  Machine dependent wrapped message-wait communication routine for MPI.

  Author:          Scott A. Hutchinson, SNL, 9221
  =======

  Return code:     int
  ============

  Parameter list:
  ===============

  buf:             Beginning address of data to be sent.

  bytes:           Length of message in bytes.
  dest:            Destination processor number.

  type:            Message type

  flag:

*******************************************************************************/

{

  int        err, count;
  MPI_Status status;

  if ( MPI_Wait(request, &status) ) {
    (void) fprintf(stderr, "MPI_Wait error = %d\n", err);
    exit(-1);
  }

  MPI_Get_count(&status, MPI_BYTE, &count);
  *source = status.MPI_SOURCE;
  *type   = status.MPI_TAG;

  /* return the count, which is in bytes */

  return count;

} /* md_wrap_wait */

#endif /* has MPI */

/* reset unused arguments warning */
#ifdef __MWERKS__
#pragma warn_unusedarg  reset
#endif /* __MWERKS__ */           
