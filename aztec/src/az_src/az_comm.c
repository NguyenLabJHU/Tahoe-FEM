/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile: az_comm.c,v $
 *
 * $Author: paklein $
 *
 * $Date: 2001-01-30 20:59:13 $
 *
 * $Revision: 1.1.1.1 $
 *
 * $Name: not supported by cvs2svn $
 *====================================================================*/
#ifndef lint
static char rcsid[] = "$Id: az_comm.c,v 1.1.1.1 2001-01-30 20:59:13 paklein Exp $";
#endif


/*******************************************************************************
 * Copyright 1995, Sandia Corporation.  The United States Government retains a *
 * nonexclusive license in this software as prescribed in AL 88-1 and AL 91-7. *
 * Export of this program may require a license from the United States         *
 * Government.                                                                 *
 ******************************************************************************/

#include "az_aztec.h"
#include "f2c.h"
#undef abs

/* System Include files */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifndef TRUE
#define TRUE  1
#define FALSE 0
#endif

int AZ_sys_msg_type = AZ_MSG_TYPE;

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_exchange_bdry(double x[], int data_org[])

/*******************************************************************************

  Routine to locally exchange the components of the vector "x". This routine
  gathers the necessary components of the vector and then sends the required
  "num_neighbors" messages. The messages which are received are placed
  contiguously in the external_nodes part of x.

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  x:               Vector of unknowns defined on the current processor.
                   Indirect addressing will be used to gather those unknowns
                   that other processors need.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see file params.txt).

*******************************************************************************/

{

  /* local variables */

  char *message_send_add[AZ_MAX_NEIGHBORS];
  /* message_send_add[i] points to the beginning of the list of values to be
     sent to the ith neighbor (i.e. data_org[AZ_neighbors+i] or sometimes
     locally defined as proc_num_neighbor[i]). That is, *(message_send_add[i]+j)
     is the jth value to be sent to the ith neighbor. */
  int  message_send_length[AZ_MAX_NEIGHBORS];
  /* message_send_length[i] is the number of bytes to be sent to the ith
     neighbor (i.e. data_org[AZ_neighbors+i] or sometimes locally defined as
     proc_num_neighbor[i]). */
  char *message_recv_add[AZ_MAX_NEIGHBORS];
  /* message_recv_add[i] points to the beginning of the list of locations which
     are to receive values sent by the ith neighbor (i.e.
     data_org[AZ_neighbors+i] or sometimes locally defined as
     proc_num_neighbor[i]). That is, *(message_recv_add[i] + j) is the location
     where the jth value sent from the ith neighbor will be stored. */
  int  message_recv_length[AZ_MAX_NEIGHBORS];
  /* message_recv_length[i] is the number of bytes to be sent to the ith
     neighbor (i.e. data_org[AZ_neighbors+i] or sometimes locally defined as
     proc_num_neighbor[i]). */

  int              n;
  double          *ptr_send_list, *ptr_recv_list;
  register double *ptrd;
  register int    *ptr_int, i;
  int              size, num_send, num_recv;
  int              type;
  int              num_neighbors, *proc_num_neighbor, total_num_send_unknowns;
  int             *num_unknowns_send_neighbor, *list_send_unknowns;
  int              external_index, *num_unknowns_recv_neighbor;

  /**************************** execution begins ******************************/

  num_neighbors              = data_org[AZ_N_neigh];
  proc_num_neighbor          = &data_org[AZ_neighbors];
  total_num_send_unknowns    = data_org[AZ_total_send];
  num_unknowns_send_neighbor = &data_org[AZ_send_length];
  list_send_unknowns         = &data_org[AZ_send_list];
  external_index             = data_org[AZ_N_internal] + data_org[AZ_N_border];
  num_unknowns_recv_neighbor = &data_org[AZ_rec_length];

  type            = AZ_sys_msg_type;
  AZ_sys_msg_type = (AZ_sys_msg_type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;

  /* single processor case */

  if (num_neighbors == 0) return;

  /* Set up send messages: Gather send unknowns from "x" vector */

  ptrd = (double *) AZ_manage_memory(data_org[AZ_total_send]*sizeof(double),
                                     AZ_ALLOC, AZ_SYS, "ptrd", &n);
  ptr_send_list = ptrd;

  ptr_int = list_send_unknowns;
  for (i = total_num_send_unknowns; i > 0; i--) {
    *ptrd++ = x[*ptr_int++];
  }

  /* Define arrays for message passing */

  ptr_recv_list = &x[external_index];

  size = sizeof(double);
  for (n = 0; n < num_neighbors; n++) {
    num_send               = num_unknowns_send_neighbor[n];
    num_recv               = num_unknowns_recv_neighbor[n];
    message_send_add[n]    = (char *) ptr_send_list;
    message_recv_add[n]    = (char *) ptr_recv_list;
    message_send_length[n] = size * num_send;
    message_recv_length[n] = size * num_recv;
    ptr_send_list         += num_send;
    ptr_recv_list         += num_recv;
  }

  AZ_exchange_local_info(num_neighbors, proc_num_neighbor, message_send_add,
                         message_send_length, message_recv_add,
                         message_recv_length, type);

} /* AZ_exchange_bdry */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_exchange_local_info(int num_neighbors, int proc_num_neighbor[],
                            char *message_send_add[], int message_send_length[],
                            char *message_recv_add[], int message_recv_length[],
                            int type)

/*******************************************************************************

  Routine to communicate between a small number of processors by first writing
  all messages and then reading all messages. This generic routine can be called
  with the standard addresses, message lengths and message types for a wide
  variety of uses.

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  num_neighbors:             Total number of neighboring processors.

  proc_num_neighbor:         Array containing the processor id for each of the
                             current processor's neighbors.

  message_send_add:          message_send_add[i] points to the beginning of
                             the list of values to be sent to the ith neighbor
                             (i.e. data_org[AZ_neighbors+i] or sometimes locally
                             defined as proc_num_neighbor[i]). That is,
                             *(message_send_add[i] + j) is the jth value to be
                             sent to the ith neighbor.

  message_send_length:       message_send_length[i] is the number of bytes to
                             be sent to the ith neighbor (i.e.
                             data_org[AZ_neighbors+i] or sometimes locally
                             defined as proc_num_neighbor[i]).

  message_recv_add:          message_recv_add[i] points to the beginning of the
                             list of locations which are to receive values sent
                             by the ith neighbor (i.e.  data_org[AZ_neighbors+i]
                             or sometimes locally defined as
                             proc_num_neighbor[i]). That is,
                             *(message_recv_add[i] + j) is the location where
                             the jth value sent from the ith neighbor will be
                             stored.

  message_recv_length:       message_recv_length[i] is the number of bytes to
                             be sent to the ith neighbor (i.e.
                             data_org[AZ_neighbors+i] or sometimes locally
                             defined as proc_num_neighbor[i]).

  type:                      message type used when exchanging information.

*******************************************************************************/

{

  /* local declarations */

  register int n;
  int          rtype, st;
  int mesg_from;

  MPI_Request request[AZ_MAX_NEIGHBORS];  /* Message handle */

  /*********************** first executable statment *****************/

  /* post receives for all messages */

  for (n = 0; n < num_neighbors; n++) {
    rtype = type;
    (void) md_wrap_iread((void *) *(message_recv_add+n),
                         *(message_recv_length+n), proc_num_neighbor+n, &rtype,
                         request+n);
  }

  /* write out all messages */

  for (n = 0; n < num_neighbors; n++) {
    (void) md_wrap_write((void *) *(message_send_add+n),
                         *(message_send_length+n), *(proc_num_neighbor+n),
                         rtype, &st);
  }

  /* wait for all messages */

  for (n = 0; n < num_neighbors; n++) {
    rtype = type;
    (void) md_wrap_wait((void *) *(message_recv_add+n),
                        *(message_recv_length+n), proc_num_neighbor+n, &rtype,
                        &st, request+n);
  }

} /* AZ_exchange_local_info */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_gather_mesg_info(double x[] ,int data_org[], char *message_recv_add[],
                         char *message_send_add[], int message_recv_length[],
                         int message_send_length[])

/*******************************************************************************

  Routine to locally exchange the components of the vector "x". This routine
  gathers the necessary components of the vector and then sends the required
  "num_neighbors" messages. The messages which are received are placed
  contiguously in the external_nodes part of x.

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  x:                         Vector of unknowns defined on the current
                             processor. Indirect addressing will be used to
                             gather those unknowns that other processors need.

  data_org:                  Array containing information on the distribution of
                             the matrix to this processor as well as
                             communication parameters (see file params.txt).

  message_send_add:          message_send_add[i] points to the beginning of
                             the list of values to be sent to the ith neighbor
                             (i.e. data_org[AZ_neighbors+i] or sometimes locally
                             defined as proc_num_neighbor[i]). That is,
                             *(message_send_add[i] + j) is the jth value to be
                             sent to the ith neighbor.

  message_send_length:       message_send_length[i] is the number of bytes to
                             be sent to the ith neighbor (i.e.
                             data_org[AZ_neighbors+i] or sometimes locally
                             defined as proc_num_neighbor[i]).

  message_recv_add:          message_recv_add[i] points to the beginning of the
                             list of locations which are to receive values sent
                             by the ith neighbor (i.e.  data_org[AZ_neighbors+i]
                             or sometimes locally defined as
                             proc_num_neighbor[i]). That is,
                             *(message_recv_add[i] + j) is the location where
                             the jth value sent from the ith neighbor will be
                             stored.

  message_recv_length:       message_recv_length[i] is the number of bytes to
                             be sent to the ith neighbor (i.e.
                             data_org[AZ_neighbors+i] or sometimes locally
                             defined as proc_num_neighbor[i]).

  num_neighbors:             Total number of neighboring processors.


  Important external definitions:
  ===============================

  Num_Neighbors:             Total number of neighboring processors.
                             (type - int value)

  Proc_Num_Neighbor[]:       Array containing the processor id for each of the
                             current processor's neighbors.
                             (type - int vector with fixed length
                             AZ_MAX_NEIGHBORS)

  Total_Num_Send_Unknowns:   Total number of unknowns which are sent to
                             neighboring processors (size of mesg buff).
                             (type - int value)

  Num_Unknowns_Send_Neighbor[]:
                             Vector containing the number of unknowns to be sent
                             to each processor.
                             (type - int vector with fixed length
                             AZ_MAX_NEIGHBORS)

  List_Send_Unknowns[]:      Vector of local node numbers for the unknowns to be
                             sent to each neighbor.
                             (type - pointer to variable length int vector)

  Num_Unknowns_Recv_Neighbor[]:
                             Array of number of unknowns to be received from
                             each of the neighboring processors.
                             (type - int vector with fixed length
                             AZ_MAX_NEIGHBORS)

*******************************************************************************/

{

  /* local variables */

  double          *ptr_send_list, *ptr_recv_list;
  register double *ptrd;
  register int    *ptr_int, i;
  int              size, num_send, num_recv;
  int              n;

  int              Num_Neighbors;
  int             *Num_Unknowns_Recv_Neighbor, *Num_Unknowns_Send_Neighbor;

  /*********************** first executable statement *****************/

  Num_Neighbors              =  data_org[AZ_N_neigh];
  Num_Unknowns_Recv_Neighbor = &data_org[AZ_rec_length];
  Num_Unknowns_Send_Neighbor = &data_org[AZ_send_length];

  /* single processor case */

  if (Num_Neighbors == 0) return;

  /* Set up send messages:  Gather send unknowns from "x" vector */

  ptrd = (double *) AZ_manage_memory(data_org[AZ_total_send]*sizeof(double),
                                     AZ_ALLOC, AZ_SYS, "ptrd", &n);
  ptr_send_list = ptrd;

  ptr_int = &data_org[AZ_send_list];
  for (i = data_org[AZ_total_send]; i > 0; i--) {
    *ptrd++ = x[*ptr_int++];
  }

  /* Define arrays for message passing */

  ptr_recv_list = &x[ data_org[AZ_N_internal] + data_org[AZ_N_border] ];

  size = sizeof(double);
  for (n = 0; n < Num_Neighbors; n++) {
    num_send               = Num_Unknowns_Send_Neighbor[n];
    num_recv               = Num_Unknowns_Recv_Neighbor[n];
    message_send_add[n]    = (char *) ptr_send_list;
    message_recv_add[n]    = (char *) ptr_recv_list;
    message_send_length[n] = size * num_send;
    message_recv_length[n] = size * num_recv;
    ptr_send_list         += num_send;
    ptr_recv_list         += num_recv;
  }

} /* AZ_gather_mesg_info */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int AZ_gsum_int(int val, int proc_config[])

/*******************************************************************************

  Global integer sum.

  Author:
  =======

  Return code:     int, result of global sum.
  ============

  Parameter list:
  ===============

  val:             Individual processor value to be summed.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  /* local variables */

  int   type;             /* type of next message */
  int   partner;          /* processor I exchange with */
  int   mask;             /* bit pattern identifying partner */
  int   hbit;             /* largest nonzero bit in nprocs */
  int   nprocs_small;     /* largest power of 2 <= nprocs */
  int   val2;             /* arriving value to add */
  int   cflag;            /* dummy argument for compatability */
  int   node, nprocs;
  char *yo = "AZ_gsum_int: ";

  MPI_Request request;  /* Message handle */

  /*********************** first executable statment *****************/

  node            = proc_config[AZ_node];
  nprocs          = proc_config[AZ_N_procs];
  type            = AZ_sys_msg_type;
  AZ_sys_msg_type = (AZ_sys_msg_type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;

  /* Find next lower power of 2. */

  for (hbit = 0; (nprocs >> hbit) != 1; hbit++);

  nprocs_small = 1 << hbit;
  if (nprocs_small * 2 == nprocs) {
    nprocs_small *= 2;
    hbit++;
  }

  partner = node ^ nprocs_small;
  if (node+nprocs_small < nprocs) {

    /* post receives on the hypercube portion of the machine partition */

    if (md_wrap_iread((void *) &val2, sizeof(int), &partner, &type, &request)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_iread failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }
  else if (node & nprocs_small) {

    /*
     * Send messages from the portion of the machine partition "above" the
     * largest hypercube to the hypercube portion.
     */

    if (md_wrap_write((void *) &val, sizeof(int), partner, type, &cflag)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_write failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  if (node+nprocs_small < nprocs) {

    /* wait to receive the messages */

    if (md_wrap_wait((void *) &val2, sizeof(int), &partner, &type, &cflag,
                     &request) != sizeof(int)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_wait failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }

    /* sum values */

    val += val2;
  }

  /* Now do a binary exchange on nprocs_small nodes. */

  if (!(node & nprocs_small)) {
    for (mask = nprocs_small>>1; mask; mask >>= 1) {
      partner = node ^ mask;
      if (md_wrap_iread((void *) &val2, sizeof(int), &partner, &type,
                        &request)) {
        (void) fprintf(stderr, "%sERROR on node %d\nmd_iread failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (md_wrap_write((void *) &val, sizeof(int), partner, type, &cflag)) {
        (void) fprintf(stderr, "%sERROR on node %d\nmd_write failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (md_wrap_wait((void *) &val2, sizeof(int), &partner, &type, &cflag,
                       &request) != sizeof(int)) {
        (void) fprintf(stderr, "%sERROR on node %d\nmd_wait failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      val += val2;
    }
  }

  /* Finally, send message from lower half to upper half. */

  partner = node ^ nprocs_small;
  if (node & nprocs_small) {
    if (md_wrap_iread((void *) &val, sizeof(int), &partner, &type, &request)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_iread failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  else if (node+nprocs_small < nprocs ) {
    if (md_wrap_write((void *) &val, sizeof(int), partner, type, &cflag)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_write failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  if (node & nprocs_small) {
    if (md_wrap_wait((void *) &val, sizeof(int), &partner, &type, &cflag,
                     &request) != sizeof(int)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_wait failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  return val;

} /* AZ_gsum_int */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

#ifndef PUMA_GSUMD  /* Only compile the generic gsum_double function in */
                    /* if the "optimized" gsum_double for PUMA is not   */
                    /* desired.                                         */

double AZ_gsum_double(double val, int proc_config[])

/*******************************************************************************

  Global double sum.

  Author:
  =======

  Return code:     double, result of global sum.
  ============

  Parameter list:
  ===============

  val:             Individual processor value to be summed.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  /* local variables */

  int    type;          /* type of next message */
  int    partner;       /* processor I exchange with */
  int    mask;          /* bit pattern identifying partner */
  int    hbit;          /* largest nonzero bit in nprocs */
  int    nprocs_small;  /* largest power of 2 <= nprocs */
  double val2;          /* arriving value to add */
  int    cflag;         /* dummy argument for compatability */
  int    node, nprocs;
  char  *yo = "AZ_gsum_double: ";

  MPI_Request request;  /* Message handle */

  /**************************** execution begins ******************************/

  node   = proc_config[AZ_node];
  nprocs = proc_config[AZ_N_procs];

  type            = AZ_sys_msg_type;
  AZ_sys_msg_type = (AZ_sys_msg_type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;

  /* Find next lower power of 2. */

  for (hbit = 0; (nprocs >> hbit) != 1; hbit++);

  nprocs_small = 1 << hbit;
  if (nprocs_small*2 == nprocs) {
    nprocs_small *= 2;
    hbit++;
  }

  partner = node ^ nprocs_small;
  if (node+nprocs_small < nprocs) {

    /* post receives on the hypercube portion of the machine partition */

    if (md_wrap_iread((void *) &val2, sizeof(double), &partner, &type,
                      &request)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_iread failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }
  else if (node & nprocs_small) {

    /*
     * Send messages from the portion of the machine partition "above" the
     * largest hypercube to the hypercube portion.
     */

    if (md_wrap_write((void *) &val, sizeof(double), partner, type, &cflag)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_write failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  if (node+nprocs_small < nprocs) {

    /* wait to receive the messages */

    if (md_wrap_wait((void *) &val2, sizeof(double), &partner, &type, &cflag,
                     &request) != sizeof(double)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_wait failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }

    /* sum values */

    val += val2;
  }

  /* Now do a binary exchange on nprocs_small nodes. */

  if (!(node & nprocs_small)) {
    for (mask = nprocs_small>>1; mask; mask >>= 1) {
      partner = node ^ mask;

      if (md_wrap_iread((void *) &val2, sizeof(double), &partner, &type,
                        &request)) {
        (void) fprintf(stderr, "%sERROR on node %d\nmd_iread failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (md_wrap_write((void *) &val, sizeof(double), partner, type, &cflag)) {
        (void) fprintf(stderr, "%sERROR on node %d\nmd_write failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (md_wrap_wait((void *) &val2, sizeof(double), &partner, &type, &cflag,
                       &request) != sizeof(double)) {
        (void) fprintf(stderr, "%sERROR on node %d\nmd_wait failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      val += val2;
    }
  }

  /* Finally, send message from lower half to upper half. */

  partner = node ^ nprocs_small;
  if (node & nprocs_small) {
    if (md_wrap_iread((void *) &val, sizeof(double), &partner, &type,
                      &request)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_iread failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  else if (node+nprocs_small < nprocs ) {
    if (md_wrap_write((void *) &val, sizeof(double), partner, type, &cflag)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_write failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  if (node & nprocs_small) {
    if (md_wrap_wait((void *) &val, sizeof(double), &partner, &type, &cflag,
                     &request) != sizeof(double)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_wait failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  return val;

} /* AZ_gsum_double */

#endif    /* ifndef PUMA_GSUMD */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

double AZ_gmax_double(double val, int proc_config[])

/*******************************************************************************

  Global max of type double.

  Author:
  =======

  Return code:     double, maximum value across all processors.
  ============

  Parameter list:
  ===============

  val:             Individual processor value.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  /* local variables */

  int    type;          /* type of next message */
  int    partner;       /* processor I exchange with */
  int    mask;          /* bit pattern identifying partner */
  int    hbit;          /* largest nonzero bit in nprocs */
  int    nprocs_small;  /* largest power of 2 <= nprocs */
  double val2;          /* arriving value to add */
  int    cflag;         /* dummy argument for compatability */
  int    node, nprocs;
  char  *yo = "AZ_gmax_double: ";

  MPI_Request request;  /* Message handle */

  /**************************** execution begins ******************************/

  node   = proc_config[AZ_node];
  nprocs = proc_config[AZ_N_procs];

  type            = AZ_sys_msg_type;
  AZ_sys_msg_type = (AZ_sys_msg_type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;

  /* Find next lower power of 2. */

  for (hbit = 0; (nprocs >> hbit) != 1; hbit++);

  nprocs_small = 1 << hbit;
  if (nprocs_small*2 == nprocs) {
    nprocs_small *= 2;
    hbit++;
  }

  partner = node ^ nprocs_small;
  if (node+nprocs_small < nprocs) {

    /* post receives on the hypercube portion of the machine partition */

    if (md_wrap_iread((void *) &val2, sizeof(double), &partner, &type,
                      &request)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_iread failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }
  else if (node & nprocs_small) {

    /*
     * Send messages from the portion of the machine partition "above" the
     * largest hypercube to the hypercube portion.
     */

    if (md_wrap_write((void *) &val, sizeof(double), partner, type, &cflag)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_write failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  if (node+nprocs_small < nprocs) {

    /* wait to receive the messages */

    if (md_wrap_wait((void *) &val2, sizeof(double), &partner, &type, &cflag,
                     &request) != sizeof(double)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_wait failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }

    /* get max value */

    if (val2 > val) val = val2;
  }

  /* Now do a binary exchange on nprocs_small nodes. */

  if (!(node & nprocs_small)) {
    for (mask = nprocs_small>>1; mask; mask >>= 1) {
      partner = node ^ mask;

      if (md_wrap_iread((void *) &val2, sizeof(double), &partner, &type,
                        &request)) {
        (void) fprintf(stderr, "%sERROR on node %d\nmd_iread failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (md_wrap_write((void *) &val, sizeof(double), partner, type, &cflag)) {
        (void) fprintf(stderr, "%sERROR on node %d\nmd_write failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (md_wrap_wait((void *) &val2, sizeof(double), &partner, &type, &cflag,
                       &request) != sizeof(double)) {
        (void) fprintf(stderr, "%sERROR on node %d\nmd_wait failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (val2 > val)
        val = val2;
    }
  }

  /* Finally, send message from lower half to upper half. */

  partner = node ^ nprocs_small;
  if (node & nprocs_small) {
    if (md_wrap_iread((void *) &val, sizeof(double), &partner, &type,
                      &request)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_iread failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  else if (node+nprocs_small < nprocs ) {
    if (md_wrap_write((void *) &val, sizeof(double), partner, type, &cflag)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_write failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  if (node & nprocs_small) {
    if (md_wrap_wait((void *) &val, sizeof(double), &partner, &type, &cflag,
                     &request) != sizeof(double)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_wait failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  return val;

} /* AZ_gmax_double */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

double AZ_gmin_double(double val, int proc_config[])

/*******************************************************************************

  Global min of type double.

  Author:
  =======

  Return code:     double, minimum value across all processors.
  ============

  Parameter list:
  ===============

  val:             Individual processor value.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  /* local variables */

  int    type;               /* type of next message */
  int    partner;            /* processor I exchange with */
  int    mask;               /* bit pattern identifying partner */
  int    hbit;               /* largest nonzero bit in nprocs */
  int    nprocs_small;       /* largest power of 2 <= nprocs */
  double val2;               /* arriving value to add */
  int    cflag;              /* dummy argument for compatability */
  int    node, nprocs;
  char  *yo = "AZ_gmin_double: ";

  MPI_Request request;  /* Message handle */

  /**************************** execution begins ******************************/

  node   = proc_config[AZ_node];
  nprocs = proc_config[AZ_N_procs];

  type            = AZ_sys_msg_type;
  AZ_sys_msg_type = (AZ_sys_msg_type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;

  /* Find next lower power of 2. */

  for (hbit = 0; (nprocs >> hbit) != 1; hbit++);

  nprocs_small = 1 << hbit;
  if (nprocs_small*2 == nprocs) {
    nprocs_small *= 2;
    hbit++;
  }

  partner = node ^ nprocs_small;
  if (node+nprocs_small < nprocs) {

    /* post receives on the hypercube portion of the machine partition */

    if (md_wrap_iread((void *) &val2, sizeof(double), &partner, &type,
                      &request)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_iread failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }
  else if (node & nprocs_small) {

    /*
     * Send messages from the portion of the machine partition "above" the
     * largest hypercube to the hypercube portion.
     */

    if (md_wrap_write((void *) &val, sizeof(double), partner, type, &cflag)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_write failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  if (node+nprocs_small < nprocs) {

    /* wait to receive the messages */

    if (md_wrap_wait((void *) &val2, sizeof(double), &partner, &type, &cflag,
                     &request) != sizeof(double)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_wait failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }

    /* get max value */

    if (val2 < val) val = val2;
  }

  /* Now do a binary exchange on nprocs_small nodes. */

  if (!(node & nprocs_small)) {
    for (mask = nprocs_small>>1; mask; mask >>= 1) {
      partner = node ^ mask;

      if (md_wrap_iread((void *) &val2, sizeof(double), &partner, &type,
                        &request)) {
        (void) fprintf(stderr, "%sERROR on node %d\nmd_iread failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (md_wrap_write((void *) &val, sizeof(double), partner, type, &cflag)) {
        (void) fprintf(stderr, "%sERROR on node %d\nmd_write failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (md_wrap_wait((void *) &val2, sizeof(double), &partner, &type, &cflag,
                       &request) != sizeof(double)) {
        (void) fprintf(stderr, "%sERROR on node %d\nmd_wait failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (val2 < val)
        val = val2;
    }
  }

  /* Finally, send message from lower half to upper half. */

  partner = node ^ nprocs_small;
  if (node & nprocs_small) {
    if (md_wrap_iread((void *) &val, sizeof(double), &partner, &type,
                      &request)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_iread failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  else if (node+nprocs_small < nprocs ) {
    if (md_wrap_write((void *) &val, sizeof(double), partner, type, &cflag)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_write failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  if (node & nprocs_small) {
    if (md_wrap_wait((void *) &val, sizeof(double), &partner, &type, &cflag,
                     &request) != sizeof(double)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_wait failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  return val;

} /* AZ_gmin_double */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int AZ_gmax_int(int val, int proc_config[])

/*******************************************************************************

  Global max of type int.

  Author:
  =======

  Return code:     int, maximum value across all processors.
  ============

  Parameter list:
  ===============

  val:             Individual processor value.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  /* local variables */

  int   type;                     /* type of next message */
  int   partner;                  /* processor I exchange with */
  int   mask;                     /* bit pattern identifying partner */
  int   hbit;                     /* largest nonzero bit in nprocs */
  int   nprocs_small;             /* largest power of 2 <= nprocs */
  int   val2;                     /* arriving value to add */
  int   cflag;                    /* dummy argument for compatability */
  int   node, nprocs;
  char *yo = "AZ_gmax_int: ";

  MPI_Request request;  /* Message handle */

  /**************************** execution begins ******************************/

  node   = proc_config[AZ_node];
  nprocs = proc_config[AZ_N_procs];

  type            = AZ_sys_msg_type;
  AZ_sys_msg_type = (AZ_sys_msg_type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;

  /* Find next lower power of 2. */

  for (hbit = 0; (nprocs >> hbit) != 1; hbit++);

  nprocs_small = 1 << hbit;
  if (nprocs_small*2 == nprocs) {
    nprocs_small *= 2;
    hbit++;
  }

  partner = node ^ nprocs_small;
  if (node+nprocs_small < nprocs) {

    /* post receives on the hypercube portion of the machine partition */

    if (md_wrap_iread((void *) &val2, sizeof(int), &partner, &type, &request)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_iread failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }
  else if (node & nprocs_small) {

    /*
     * Send messages from the portion of the machine partition "above" the
     * largest hypercube to the hypercube portion.
     */

    if (md_wrap_write((void *) &val, sizeof(int), partner, type, &cflag)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_write failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  if (node+nprocs_small < nprocs) {

    /* wait to receive the messages */

    if (md_wrap_wait((void *) &val2, sizeof(int), &partner, &type, &cflag,
                     &request) != sizeof(int)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_wait failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }

    /* get max value */

    if (val2 > val) val += val2;
  }

  /* Now do a binary exchange on nprocs_small nodes. */

  if (!(node & nprocs_small)) {
    for (mask = nprocs_small>>1; mask; mask >>= 1) {
      partner = node ^ mask;
      if (md_wrap_iread((void *) &val2, sizeof(int), &partner, &type,
                        &request)) {
        (void) fprintf(stderr, "%sERROR on node %d\nmd_iread failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (md_wrap_write((void *) &val, sizeof(int), partner, type, &cflag)) {
        (void) fprintf(stderr, "%sERROR on node %d\nmd_write failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (md_wrap_wait((void *) &val2, sizeof(int), &partner, &type, &cflag,
                       &request) != sizeof(int)) {
        (void) fprintf(stderr, "%sERROR on node %d\nmd_wait failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (val2 > val) val = val2;
    }
  }

  /* Finally, send message from lower half to upper half. */

  partner = node ^ nprocs_small;
  if (node & nprocs_small) {
    if (md_wrap_iread((void *) &val, sizeof(int), &partner, &type, &request)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_iread failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  else if (node+nprocs_small < nprocs ) {
    if (md_wrap_write((void *) &val, sizeof(int), partner, type, &cflag)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_write failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  if (node & nprocs_small) {
    if (md_wrap_wait((void *) &val, sizeof(int), &partner, &type, &cflag,
                     &request) != sizeof(int)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_wait failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  return val;

} /* AZ_gmax_int */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

double AZ_gavg_double(double val, int proc_config[])

/*******************************************************************************

  Global average of type double.

  Author:
  =======

  Return code:     double, average value across all processors.
  ============

  Parameter list:
  ===============

  val:             Individual processor value.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  return (AZ_gsum_double(val, proc_config) / proc_config[AZ_N_procs]);

} /* AZ_gavg_double */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int AZ_gmin_int(int val, int proc_config[])

/*******************************************************************************

  Global min of type int.

  Author:
  =======

  Return code:     int, minimum value across all processors.
  ============

  Parameter list:
  ===============

  val:             Individual processor value.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  /* local variables */

  int   type;             /* type of next message */
  int   partner;          /* processor I exchange with */
  int   mask;             /* bit pattern identifying partner */
  int   hbit;             /* largest nonzero bit in nprocs */
  int   nprocs_small;     /* largest power of 2 <= nprocs */
  int   val2;             /* arriving value to add */
  int   cflag;            /* dummy argument for compatability */
  int   node, nprocs;
  char *yo = "AZ_gmin_int: ";

  MPI_Request request;  /* Message handle */

  /**************************** execution begins ******************************/

  node   = proc_config[AZ_node];
  nprocs = proc_config[AZ_N_procs];

  type            = AZ_sys_msg_type;
  AZ_sys_msg_type = (AZ_sys_msg_type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;

  /* Find next lower power of 2. */

  for (hbit = 0; (nprocs >> hbit) != 1; hbit++);

  nprocs_small = 1 << hbit;
  if (nprocs_small*2 == nprocs) {
    nprocs_small *= 2;
    hbit++;
  }

  partner = node ^ nprocs_small;
  if (node+nprocs_small < nprocs) {

    /* post receives on the hypercube portion of the machine partition */

    if (md_wrap_iread((void *) &val2, sizeof(int), &partner, &type, &request)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_iread failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }
  else if (node & nprocs_small) {

    /*
     * Send messages from the portion of the machine partition "above" the
     * largest hypercube to the hypercube portion.
     */

    if (md_wrap_write((void *) &val, sizeof(int), partner, type, &cflag)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_write failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  if (node+nprocs_small < nprocs) {

    /* wait to receive the messages */

    if (md_wrap_wait((void *) &val2, sizeof(int), &partner, &type, &cflag,
                     &request) != sizeof(int)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_wait failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }

    /* get min value */

    if (val2 < val) val = val2;
  }

  /* Now do a binary exchange on nprocs_small nodes. */

  if (!(node & nprocs_small)) {
    for (mask = nprocs_small>>1; mask; mask >>= 1) {
      partner = node ^ mask;
      if (md_wrap_iread((void *) &val2, sizeof(int), &partner, &type,
                        &request)) {
        (void) fprintf(stderr, "%sERROR on node %d\nmd_iread failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (md_wrap_write((void *) &val, sizeof(int), partner, type, &cflag)) {
        (void) fprintf(stderr, "%sERROR on node %d\nmd_write failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (md_wrap_wait((void *) &val2, sizeof(int), &partner, &type, &cflag,
                       &request) != sizeof(int)) {
        (void) fprintf(stderr, "%sERROR on node %d\nmd_wait failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (val2 < val) val = val2;
    }
  }

  /* Finally, send message from lower half to upper half. */

  partner = node ^ nprocs_small;
  if (node & nprocs_small) {
    if (md_wrap_iread((void *) &val, sizeof(int), &partner, &type, &request)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_iread failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  else if (node+nprocs_small < nprocs ) {
    if (md_wrap_write((void *) &val, sizeof(int), partner, type, &cflag)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_write failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  if (node & nprocs_small) {
    if (md_wrap_wait((void *) &val, sizeof(int), &partner, &type, &cflag,
                     &request) != sizeof(int)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_wait failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  return val;

} /* AZ_gmin_int */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_gdot_vec(int N, double dots[], double dots2[], int proc_config[])

/*******************************************************************************

  Author:
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  N:               Length of vectors 'dots' and 'dots2'.

  dots:            On output, AZ_gdot_vec() produces N global sums. Specifically
                   dots[i] (on all processors) = dots[i] (on proc 0) +
                                                 dots[i] (on proc 1) +
                                                      .
                                                      .
                                                      .
                                                 dots[i] (on proc proc_config[
                                                          AZ_N_procs])
                   for 0 <= i < N .

  dots2:           double precision work vector of length N.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  /* local variables */

  int   type;             /* type of next message */
  int   msg_size;         /* length of messages */
  int   partner;          /* processor I exchange with */
  int   mask;             /* bit pattern identifying partner */
  int   hbit;             /* largest nonzero bit in nprocs */
  int   nprocs_small;     /* largest power of 2 <= nprocs */
  int   cflag;            /* dummy argument for compatability */
  int   i;                /* loop counter */
  int   node, nprocs;
  char *yo = "AZ_gdot_vec: ";

  MPI_Request request;  /* Message handle */

  /**************************** execution begins ******************************/

  node   = proc_config[AZ_node];
  nprocs = proc_config[AZ_N_procs];

  type            = AZ_sys_msg_type;
  AZ_sys_msg_type = (AZ_sys_msg_type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;

  /* Find next lower power of 2. */

  for (hbit = 0; (nprocs >> hbit) != 1; hbit++);

  nprocs_small = 1 << hbit;
  if (nprocs_small * 2 == nprocs) {
    nprocs_small *= 2;
    hbit++;
  }

  msg_size = N * sizeof(double);
  partner  = node ^ nprocs_small;
  if (node+nprocs_small < nprocs) {

    /* post receives on the hypercube portion of the machine partition */

    if (md_wrap_iread((void *) dots2, msg_size, &partner, &type, &request)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_iread failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }
  else if (node & nprocs_small) {

    /*
     * Send messages from the portion of the machine partition "above" the
     * largest hypercube to the hypercube portion.
     */

    if (md_wrap_write((void *) dots, msg_size, partner, type, &cflag)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_write failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  if (node+nprocs_small < nprocs) {

    /* wait to receive the messages */

    if (md_wrap_wait((void *) dots2, msg_size, &partner, &type, &cflag,
                     &request) != msg_size) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_wait failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }

    /* sum values */

    for (i = 0; i < N; i++) dots[i] += dots2[i];
  }

  /* Now do a binary exchange on nprocs_small nodes. */

  if (!(node & nprocs_small)) {
    for (mask = nprocs_small>>1; mask; mask >>= 1) {
      partner = node ^ mask;

      if (md_wrap_iread((void *) dots2, msg_size, &partner, &type, &request)) {
        (void) fprintf(stderr, "%sERROR on node %d\nmd_iread failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (md_wrap_write((void *) dots, msg_size, partner, type, &cflag)) {
        (void) fprintf(stderr, "%sERROR on node %d\nmd_write failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (md_wrap_wait((void *) dots2, msg_size, &partner, &type, &cflag,
                       &request) != msg_size) {
        (void) fprintf(stderr, "%sERROR on node %d\nmd_wait failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      for (i = 0; i < N; i++) dots[i] += dots2[i];
    }
  }

  /* Finally, send message from lower half to upper half. */

  partner = node ^ nprocs_small;
  if (node & nprocs_small) {
    if (md_wrap_iread((void *) dots, msg_size, &partner, &type, &request)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_iread failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  else if (node+nprocs_small < nprocs ) {
    if (md_wrap_write((void *) dots, msg_size, partner, type, &cflag)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_write failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  if (node & nprocs_small) {
    if (md_wrap_wait((void *) dots, msg_size, &partner, &type, &cflag,
                     &request) != msg_size) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_wait failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

} /* AZ_gdot_vec */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_print_sync_start(int proc, int do_print_line)

/*******************************************************************************

  Routine to allow IO between print_sync_start and print_sync_end to be printed
  by each processor entirely before the next processor begins its IO.  The
  printing sequence is from proc = 0 to the last processor,
  number_of_procs = nprocs - 1.

  NOTE: THERE CAN BE NO COMMUNICATON BETWEEN THESE CALLS.

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  proc:            Current processor number.

  do_print_line:   Boolean variable.  If true, a line of # is printed to
                   indicate the start of a print_sync I/O block.

*******************************************************************************/

{

  /* local variables */

  int flag = 1, from, st, type;

  /**************************** execution begins ******************************/

  type = AZ_sys_msg_type;

  if (proc != 0) {
    from = proc - 1;

    if (md_read((char *) &flag, sizeof(int), &from, &type, &st) !=
        sizeof(int)) {
      (void) fprintf(stderr, "print_sync_start: ERROR on node %d\n", proc);
      (void) fprintf(stderr, "md_read failed, message type %d\n", type);
      exit (-1);
    }
  }

  else {
    if (do_print_line) {
      (void) printf("\n");
      for (flag = 0; flag < 37; flag++) (void) printf("#");
      (void) printf(" PRINT_SYNC_START ");
      for (flag = 0; flag < 25; flag++) (void) printf("#");
      (void) printf("\n");
    }
  }

} /* AZ_print_sync_start */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_print_sync_end(int proc, int nprocs, int do_print_line)

/*******************************************************************************

  Routine to allow IO between print_sync_start and print_sync_end to be printed
  by each processor entirely before the next processor begins its IO. The
  printing sequence is from proc = 0 to the last processor,
  number_of_procs = nprocs - 1.

  NOTE: THERE CAN BE NO COMMUNICATON BETWEEN THESE CALLS.

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  proc:            Current processor number.

  nprocs:          Number of processors in the current machine configuration.

  do_print_line:   Boolean variable.  If true, a line of # is printed to
                   indicate the start of a print_sync I/O block.

*******************************************************************************/

{

  /* local variables */

  int st, flag = 1, from, type, to;

  /**************************** execution begins ******************************/

  type            = AZ_sys_msg_type;
  AZ_sys_msg_type = (AZ_sys_msg_type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;

  if (proc < nprocs -1) to = proc + 1;
  else {
    to = 0;
    if (do_print_line) {
      (void) printf("\n");
      for (flag = 0; flag < 37; flag++) (void) printf("#");
      (void) printf(" PRINT_SYNC_END__ ");
      for (flag = 0; flag < 25; flag++) (void) printf("#");
      (void) printf("\n\n");
    }
  }

  if (md_write((char *) &flag, sizeof(int), to, type, &st) != 0) {
    (void) fprintf(stderr, "print_sync_end: ERROR on node %d\n", proc);
    (void) fprintf(stderr, "md_write failed, message type %d\n", type);
    exit (-1);
  }

  if (proc == 0) {
    from = nprocs -1;
    if (md_read((char *) &flag, sizeof(int), &from, &type, &st) !=
        sizeof(int)) {
      (void) fprintf(stderr, "print_sync_end: ERROR on node %d\n", proc);
      (void) fprintf(stderr, "md_read failed, message type %d/n", type);
      exit (-1);
    }
  }

  /*
   * Do a final sync amongst all the processors, so that all of the other
   * processors must wait for Proc 0 to receive the final message from
   * Proc (Num_Proc-1).
   */

  AZ_sync(proc, nprocs);

} /* AZ_print_sync_end */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_sync(int node, int nprocs)

/*******************************************************************************

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  node:            Current processor number.

  nprocs:          Number of processors in the current machine configuration.

*******************************************************************************/

{

  /* local variables */

  int   type;                     /* type of next message */
  int   partner;                  /* processor I exchange with */
  int   mask;                     /* bit pattern identifying partner */
  int   hbit;                     /* largest nonzero bit in nprocs */
  int   nprocs_small;             /* largest power of 2 <= nprocs */
  int   cflag;                    /* dummy argument for compatability */
  char *yo = "sync: ";

  MPI_Request request;  /* Message handle */

  /**************************** execution begins ******************************/

  type            = AZ_sys_msg_type;
  AZ_sys_msg_type = (AZ_sys_msg_type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;

  /*  Find next lower power of 2. */

  for (hbit = 0; (nprocs >> hbit) != 1; hbit++);

  nprocs_small = 1 << hbit;
  if (nprocs_small*2 == nprocs) {
    nprocs_small *= 2;
    hbit++;
  }

  partner = node ^ nprocs_small;
  if (node+nprocs_small < nprocs) {

    /* post receives on the hypercube portion of the machine partition */

    if (md_wrap_iread((void *) NULL, 0, &partner, &type, &request)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_iread failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }
  else if (node & nprocs_small) {

    /*
     * Send messages from the portion of the machine partition "above" the
     * largest hypercube to the hypercube portion.
     */

    if (md_wrap_write((void *) NULL, 0, partner, type, &cflag)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_write failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  if (node+nprocs_small < nprocs) {

    /*
     * Wait to receive the messages.  These messages will return length 1
     * because MPI will not necessarily send a zero-length message.
     */

    (void) md_wrap_wait((void *) NULL, 0, &partner, &type, &cflag, &request);
  }

  /*  Now do a binary exchange on nprocs_small nodes. */

  if (!(node & nprocs_small)) {
    for (mask = nprocs_small>>1; mask; mask >>= 1) {
      partner = node ^ mask;
      if (md_wrap_iread((void *) NULL, 0, &partner, &type, &request)) {
        (void) fprintf(stderr, "%sERROR on node %d\nmd_iread failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (md_wrap_write((void *) NULL, 0, partner, type, &cflag)) {
        (void) fprintf(stderr, "%sERROR on node %d\nmd_write failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      /*
       * Wait to receive the messages.  These messages will return length 1
       * because MPI will not necessarily send a zero-length message.
       */

      (void) md_wrap_wait((void *) NULL, 0, &partner, &type, &cflag, &request);
    }
  }

  /*  Finally, send message from lower half to upper half. */

  partner = node ^ nprocs_small;
  if (node & nprocs_small) {
    if (md_wrap_iread((void *) NULL, 0, &partner, &type, &request)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_iread failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  else if (node+nprocs_small < nprocs ) {
    if (md_wrap_write((void *) NULL, 0, partner, type, &cflag)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_write failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  /*
   * Wait to receive the messages.  These messages will return length 1
   * because MPI will not necessarily send a zero-length message.
   */

  if (node & nprocs_small) {
    (void) md_wrap_wait((void *) NULL, 0, &partner, &type, &cflag, &request);
  }

} /* AZ_sync */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_gsum_vec_int(int vals[], int vals2[], int length, int proc_config[])

/*******************************************************************************

  For each element in vals[], perform a global sum with the other processors.
  That is, on output vals[i] is equal to the sum of the input values in vals[i]
  on all the processors.

  Author:          Ray Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  vals:            On input, vals[i] on this processor is to be summed with
                   vals[i] on all the other processors.
                   On output, vals[i] is the sum of the input values in val[i]
                   defined on all processors.

  vals2:           Work space of size 'length'.

  node:            Current processor number.

  nprocs:          Number of processors in the current machine configuration.

  length:          Number of values in 'vals' (i.e. number of global sums).

*******************************************************************************/

{

  /* local variables */

  int   type;             /* type of next message */
  int   partner;          /* processor I exchange with */
  int   mask;             /* bit pattern identifying partner */
  int   hbit;             /* largest nonzero bit in nprocs */
  int   nprocs_small;     /* largest power of 2 <= nprocs */
  int   cflag;            /* dummy argument for compatability */
  int   k;
  int   node, nprocs;
  char *yo = "AZ_gsum_vec_int: ";

  MPI_Request request;  /* Message handle */

  /*********************** first executable statment *****************/

  node   = proc_config[AZ_node];
  nprocs = proc_config[AZ_N_procs];

  type            = AZ_sys_msg_type;
  AZ_sys_msg_type = (AZ_sys_msg_type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;

  /* Find next lower power of 2. */

  for (hbit = 0; (nprocs >> hbit) != 1; hbit++);

  nprocs_small = 1 << hbit;

  if (nprocs_small * 2 == nprocs) {
    nprocs_small *= 2;
    hbit++;
  }

  partner = node ^ nprocs_small;
  if (node+nprocs_small < nprocs) {

    /* post receives on the hypercube portion of the machine partition */

    if (md_wrap_iread((void *) vals2, length*sizeof(int), &partner, &type,
                      &request)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_iread failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }
  else if (node & nprocs_small) {

    /*
     * Send messages from the portion of the machine partition "above" the
     * largest hypercube to the hypercube portion.
     */

    if (md_wrap_write((void *) vals, length*sizeof(int), partner, type,
                      &cflag)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_write failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  if (node+nprocs_small < nprocs) {

    /* wait to receive the messages */

    if (md_wrap_wait((void *) vals2, length*sizeof(int), &partner, &type,
                     &cflag, &request) != length*sizeof(int)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_wait failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }

    /* sum values */

    for (k = 0; k < length; k++) vals[k] += vals2[k];
  }

  /* Now do a binary exchange on nprocs_small nodes. */

  if (!(node & nprocs_small)) {
    for (mask = nprocs_small >> 1; mask; mask >>= 1) {
      partner = node ^ mask;

      if (md_wrap_iread((void *) vals2, length*sizeof(int), &partner, &type,
                        &request)) {
        (void) fprintf(stderr, "%sERROR on node %d\nmd_iread failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (md_wrap_write((void *) vals, length*sizeof(int), partner, type,
                        &cflag)) {
        (void) fprintf(stderr, "%sERROR on node %d\nmd_write failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (md_wrap_wait((void *) vals2, length*sizeof(int), &partner, &type,
                       &cflag, &request) != length*sizeof(int)) {
        (void) fprintf(stderr, "%sERROR on node %d\nmd_wait failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      for (k = 0; k < length; k++) vals[k] += vals2[k];
    }
  }

  /* Finally, send message from lower half to upper half. */

  partner = node ^ nprocs_small;
  if (node & nprocs_small) {
    if (md_wrap_iread((void *) vals, length*sizeof(int), &partner, &type,
                      &request)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_iread failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  else if (node+nprocs_small < nprocs ) {
    if (md_wrap_write((void *) vals, length*sizeof(int), partner, type,
                      &cflag)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_write failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  if (node & nprocs_small) {
    if (md_wrap_wait((void *) vals, length*sizeof(int), &partner, &type, &cflag,
                     &request) != length*sizeof(int)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_wait failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

} /* AZ_gsum_vec_int */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_gappend(int vals[], int *cur_length, int total_length,
                int proc_config[])

/*******************************************************************************

  Take the list contained in vals[] that is defined on this processor and append
  it to the other lists (vals[]) defined on all the other processors.

  Author:          Ray Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  vals:            On input, vals[i] on this processor contains list to be
                   appended with vals[i] on all the other processors.
                   On output, vals[i] is the concatenation of the input lists
                   in val[i] defined on all processors.

  cur_length:      On input, initial length of vals[].
                   On output, length of new vals[] (after concatenation).

  total_length:    Maximum allowable length for vals[].

  node:            Current processor number.

  nprocs:          Number of processors in the current machine configuration.

*******************************************************************************/

{

  /* local variables */

  int   type;         /* type of next message */
  int   partner;      /* processor I exchange with */
  int   mask;         /* bit pattern identifying partner */
  int   hbit;         /* largest nonzero bit in nprocs */
  int   nprocs_small; /* largest power of 2 <= nprocs */
  int   cflag;        /* dummy argument for compatability */
  int   length;
  int   node, nprocs;
  char *yo = "AZ_gappend: ";

  MPI_Request request;  /* Message handle */

  /*********************** first executable statment *****************/

  node   = proc_config[AZ_node];
  nprocs = proc_config[AZ_N_procs];

  type            = AZ_sys_msg_type;
  AZ_sys_msg_type = (AZ_sys_msg_type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;

  /* Find next lower power of 2. */

  for (hbit = 0; (nprocs >> hbit) != 1; hbit++);

  nprocs_small = 1 << hbit;

  if (nprocs_small * 2 == nprocs) {
    nprocs_small *= 2;
    hbit++;
  }

  partner = node ^ nprocs_small;

  if (node+nprocs_small < nprocs) {

    /* post receives on the hypercube portion of the machine partition */

    if ( md_wrap_iread((void *) &(vals[*cur_length]),
                       (total_length - *cur_length) * sizeof(int), &partner,
                       &type, &request) ) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_iread failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }
  else if (node & nprocs_small) {

    /*
     * Send messages from the portion of the machine partition "above" the
     * largest hypercube to the hypercube portion.
     */

    if (md_wrap_write((void *) vals, (*cur_length)*sizeof(int), partner, type,
                      &cflag)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_write failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  if (node+nprocs_small < nprocs) {

    /* wait to receive the messages */

    length = md_wrap_wait((void *) &(vals[*cur_length]),
                          (total_length - *cur_length)*sizeof(int), &partner,
                          &type, &cflag, &request);
    (*cur_length) += (length / sizeof(int));
  }

  /* Now do a binary exchange on nprocs_small nodes. */

  if (!(node & nprocs_small)) {
    for (mask = nprocs_small >> 1; mask; mask >>= 1) {
      partner = node ^ mask;

      if (md_wrap_iread((void *) &(vals[*cur_length]),
                        (total_length - *cur_length)*sizeof(int), &partner,
                        &type, &request)) {
        (void) fprintf(stderr, "%sERROR on node %d\nmd_iread failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (md_wrap_write((void *) vals, *cur_length*sizeof(int), partner, type,
                        &cflag)) {
        (void) fprintf(stderr, "%sERROR on node %d\nmd_write failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      length = md_wrap_wait((void *) &(vals[*cur_length]),
                            (total_length - *cur_length)*sizeof(int), &partner,
                            &type, &cflag, &request);
      (*cur_length) += (length / sizeof(int));
    }
  }

  /* Finally, send message from lower half to upper half. */

  partner = node ^ nprocs_small;
  if (node & nprocs_small) {
    if (md_wrap_iread((void *) vals, total_length*sizeof(int), &partner, &type,
                      &request)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_iread failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  else if (node+nprocs_small < nprocs ) {
    if (md_wrap_write((void *) vals, *cur_length*sizeof(int), partner, type,
                      &cflag)) {
      (void) fprintf(stderr, "%sERROR on node %d\nmd_write failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  if (node & nprocs_small) {
    length = md_wrap_wait((void *) vals, total_length*sizeof(int), &partner,
                          &type, &cflag, &request);
    (*cur_length) = (length / sizeof(int));
  }

} /* AZ_gappend */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_broadcast(char *ptr, int length, int proc_config[], int action)

/*******************************************************************************

  Used to concatenate a buffer of information and to then broadcast this
  information from processor 0 to the other processors. The four possiblities
  are

    1) action == AZ_PACK  && proc_config[AZ_node] == 0
          Store ptr into internal buffer (brdcst_buffer).
    2) action == AZ_PACK  && proc_config[AZ_node] != 0
          Read from the internal buffer (brdcst_buffer) to ptr. If the
          internal buffer is empty, first receive broadcast information
    3) action == AZ_SEND  && proc_config[AZ_node] == 0
          Broadcast internal buffer (filled by previous AZ_broadcast
          calls) and then clear it.
    4) action == AZ_SEND  && proc_config[AZ_node] != 0
          Clear the internal buffer.


  Sample usage:

    if (proc_config[AZ_node] == 0) {
       a = 1; b = 2;
    }
    AZ_broadcast(&a, sizeof(int), proc_config, AZ_PACK);
    AZ_broadcast(&b, sizeof(int), proc_config, AZ_PACK);
    AZ_broadcast(NULL,         0, proc_config, AZ_SEND);

  Note: There can be no other communication calls between AZ_PACK
  and AZ_SEND calls to AZ_broadcast.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  ptr:             If action == AZ_PACK,
                   If proc==0, ptr is string to be broadcast to other nodes.
                   If proc~=0, on output ptr is string received from node 0.

  length:          Length of string to be broadcast/received.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

  action:          Determines what AZ_broadcast() does as described above.

*******************************************************************************/

{

  /* local variables */

  int i;
  static int   buffer_end = 0;
  static char *brdcst_buffer = 0;
  static int   buf_len = 0;
  static int   buffer_start = 0;
  char   *temp;
  int    *tt;

  /**************************** execution begins ******************************/

  if (action == AZ_PACK) {

    /* allocate the buffer */

    if (brdcst_buffer == 0) {
      buf_len = 1000;    /* Note: this is coordinated with the     */
      /* statement 'if (buf_len != 1000)' below */

      brdcst_buffer = (char *) calloc(buf_len,sizeof(int));
      if (brdcst_buffer == NULL) {

        (void) fprintf(stderr, "no space in AZ_broadcast: brdcst_buffer\n");
        exit(-1);
      }
    }

    /* If processor 0, pack buffer */

    if (proc_config[AZ_node] == 0) {

      if (buffer_end+length > buf_len)  {

        /* Buffer is not big enough. Allocate more space */

        buf_len += max(500, length);
        temp = (char *) calloc(buf_len, sizeof(int));
        if (temp == NULL) {
          (void) fprintf(stderr, "no space in AZ_broadcast: temp\n");
          exit(-1);
        }

        if (brdcst_buffer != 0) {
          for (i = 0; i < buffer_end; i++) temp[i] = brdcst_buffer[i];
          free(brdcst_buffer);
        }
        brdcst_buffer = temp;
      }

      if (brdcst_buffer == 0) {
        (void) fprintf(stderr,
                       "Error: Not enough space in AZ_broadcast_pack\n");
        exit(-1);
      }

      for (i = 0; i < length; i++) brdcst_buffer[i + buffer_end] = ptr[i];
      buffer_end += length;
    }

    /* For processors other than 0 ... */

    else {

      /*
       * If the buffer is empty, do a broadcast (i.e. post a read) to obtain
       * broadcast information.
       */

      if (buffer_end == 0) {

        buffer_end = AZ_broadcast_info(brdcst_buffer, proc_config, buf_len);

        /*
         * A single broadcasted integer, indicates that processor 0 increased
         * its buffer size during the packing stage. The size of this buffer is
         * in fact the value received.
         */

        if (buffer_end == sizeof(int)) {

          /* allocate a new buffer */

          tt = (int *) brdcst_buffer;
          buf_len = tt[0];
          free(brdcst_buffer);
          brdcst_buffer = (char *) calloc(buf_len, sizeof(char));
          if ( brdcst_buffer == NULL) {
            (void) fprintf(stderr,
                           "no space in AZ_broadcast: brdcst_buffer \n");
            exit(-1);

          }

          /* Do the broadcast again with the larger buffer */

          buffer_end = AZ_broadcast_info(brdcst_buffer, proc_config, buf_len);
        }
      }

      /* take the top elements off of 'buffer' and store them into 'ptr'. */

      for (i = 0; i < length; i++) ptr[i] = brdcst_buffer[buffer_start+i];
      buffer_start += length;
    }
  }

  else {
    if (proc_config[AZ_node] == 0) {

      /*
       * If additional buffer space was needed, processor 0 tells the others by
       * sending just 1 integer (the length of the new buffer)
       */

      if (buf_len != 1000)
        (void) AZ_broadcast_info((char *) &buffer_end,proc_config,sizeof(int));

      /*
       * If only 1 integer needs to be sent, we increase the number of bytes
       * sent to distinguish this situation from the case above where 1 integer
       * which is the length of the new buffer is sent
       */

      if (buffer_end == sizeof(int)) buffer_end++;

      /* broadcast the data */

      (void) AZ_broadcast_info(brdcst_buffer, proc_config, buffer_end);
    }

    /* clear the internal buffer */

    if (brdcst_buffer != (char *) NULL) free(brdcst_buffer);
    brdcst_buffer = (char *) NULL;
    buf_len = buffer_end = buffer_start = 0;
  }

} /* AZ_broadcast */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int AZ_broadcast_info(char buffer[], int proc_config[], int length)

/*******************************************************************************

  The information in 'buffer' on processor 0 is broadcast to all other
  processors.

  Author:          Ray Tuminaro, SNL, 1422
  =======

  Return code:     int, length of string broadcast.
  ============

  Parameter list:
  ===============

  buffer:          Buffer used to pack information to be broadcast.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

  length:          Length of string to be broadcast/received.

*******************************************************************************/

{

  /* local variables */

  int i;
  int type;         /* type of next message */
  int partner;      /* processor I exchange with */
  int hbit;         /* largest nonzero bit in nprocs */
  int my_lbit;      /* smallest nonzero bit in proc */
  int cflag;        /* dummy argument for compatability */
  int nprocs, proc;
  char *yo = "AZ_broadcast_info: ";

  MPI_Request request;  /* Message handle */

  /*********************** first executable statment *****************/

  nprocs = proc_config[AZ_N_procs];
  proc   = proc_config[AZ_node];

  type            = AZ_sys_msg_type;
  AZ_sys_msg_type = (AZ_sys_msg_type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;

  /* Find location of largest bit. */

  for (hbit = 0; ((nprocs - 1) >> hbit) != 0; hbit++);

  /* Find location of smallest bit corresponding to my processor name */

  if (proc != 0)
    for (my_lbit = 1; (proc | (1 << (my_lbit - 1))) != proc; my_lbit++);
  else my_lbit = hbit + 1;

  /* Zero out lowest bit in proc ... and receive from that processor */

  if (proc != 0) {
    partner = proc ^ (1 << (my_lbit - 1));
    (void) md_wrap_iread((void *) buffer, length, &partner, &type, &request);

    /* wait for messages */

    length  = md_wrap_wait((void *) buffer, length, &partner, &type, &cflag,
                           &request);
  }

  /*
   * Send to neighbors. The neighbors are defined by putting a 1 in one of the
   * locations in 'proc' to the right of the lowest nonzero bit.
   */

  for (i = my_lbit - 1; i > 0; i--) {
    partner = proc | (1 << (i - 1));
    if (partner < nprocs)
      (void) md_wrap_write((void *) buffer, length, partner, type, &cflag);
  }

  return length;

} /* AZ_broadcast_info */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

double AZ_sync_timer(int proc_config[])

/*******************************************************************************

  Use to synchronize the processors and take a timing.

  Author:          Ray Tuminaro, SNL, 1422
  =======

  Return code:     double, current time.
  ============

  Parameter list:
  ===============

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  int i = 0;

  i = AZ_gsum_int(i, proc_config);

  return AZ_second();

} /* AZ_sync_timer */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_splitup_big_msg(int num_neighbors, double *buffer, int *start_send_proc,
                        int *actual_send_length, int *actual_recv_length, int
                        *proc_num_neighbor, int type, int *total_num_recv,
                        int *proc_config)

/*******************************************************************************

   Author:        John Shadid, SNL, 1421 Date: 8/1/94
  =======

  Return code:    (void).
  ============

  Parameter list:
  ===============

  num_neighbors:      total number of neighbors to communicate with
  buffer:             on input  - buffer that holds the information to be sent
                      on output - buffer contains the recieve information
  start_send_proc:    contains index of start of information that is to be
                      sent to each processor
  actual_send_length: contains number of double precision entries to be sent to
                      each processor
  actual_recv_length: number of entries to be recieve from each processor
  proc_num_neighbor:  contains the processor number for each neighbor
  type:               the message type to be used in the mesage
  total_num_recv:     on output - total number of actual recvieved entried
  proc_config:        Machine configuration.  proc_config[AZ_node] is the node
                      number.  proc_config[AZ_N_procs] is the number
                      of processors.

*******************************************************************************/

{

  /*
   * This function handshakes big messages between all the neighbors.  The
   * length of the messages are calculated conservatively to not allow overflow
   * of the message buffers.
   */

  int     m, n, st, rtype, j, length, size, dummy_int;
  int     max_neighbors, messg_size_doubles, doubles_sent;
  int     total_doubles_to_send, dest, flag, messg_from, messg_type;
  int     total_doubles_to_recv, total_send_size;
  int     finished_send_messg[AZ_MAX_NEIGHBORS];
  int     finished_recv_messg[AZ_MAX_NEIGHBORS];
  int     number_of_messages, start_recv_proc[AZ_MAX_NEIGHBORS];
  int     allowed_buff_size, num_recv;
  int     max_buffer_size = 0, max_messg_size;
  double *send_buffer;
  char   *char_ptr;
  char   *yo = "AZ_splitup_big_msg ";
  int     split_up = FALSE;
  int     dummy_add;
  int     DEBUG = FALSE;
  
  
  
  MPI_Request request[AZ_MAX_NEIGHBORS];  /* Message handle */

  /**************************** execution begins ****************************/

  /* Compute the global maximum message buffer size needed */

  for (n = 0; n < num_neighbors; n++) {
    max_buffer_size += actual_recv_length[n];
  }
  max_buffer_size = AZ_gmax_int(max_buffer_size, proc_config);

  /* Determine if splitting of messages is necessary */

  if (max_buffer_size > (AZ_MAX_MSG_BUFF_SIZE / (2 * sizeof(double)))) {

     /* Too big for message buffers */

     split_up = TRUE;
  }

  if (split_up == TRUE) {

    /*
     * Compute maximum total message size in bytes that any processor will
     * recieve and the maximum number of neighbors that any proc must
     * communicate with. Also initalize some logical arrays.
     */

    max_messg_size = 0;
    for (n = 0; n < num_neighbors; n++) {
      max_messg_size = max(max_messg_size, actual_recv_length[n]);
      finished_send_messg[n] = finished_recv_messg[n] = AZ_FALSE;
    }
    max_messg_size = AZ_gmax_int(max_messg_size, proc_config);
    max_neighbors  = AZ_gmax_int(num_neighbors, proc_config);

    /*
     * Total received nonzeros and starting location for each processors
     * message that will be received
     */

    num_recv = 0;
    for (n = 0; n < num_neighbors; n++) {
      start_recv_proc[n] = num_recv;
      num_recv          += actual_recv_length[n];
    }
    *total_num_recv = num_recv;

    /*
     * Compute the global maximum allowed message size and the maximum number of
     * messages to send all the required information.
     */

    allowed_buff_size  = (int) floor((AZ_MAX_MSG_BUFF_SIZE /
                                      (3*sizeof(double))));

    messg_size_doubles = (int) floor(allowed_buff_size / max_neighbors);

    number_of_messages = (int) ceil((double) max_messg_size /
                                    (double) (messg_size_doubles));

    if (proc_config[AZ_node] == 0 && DEBUG == TRUE) {
      (void) printf("\n\t\tSplitting up messages in splitup_big_msg\n"
                    "\t\tmax_buffer_size required  (bytes): %d\n",
                    max_buffer_size*sizeof(double));
      (void) printf("\t\tmax_buffer_size allocated (bytes): %d\n",
                    allowed_buff_size*sizeof(double));
      (void) printf("\t\tindividual message size   (bytes): %d\n",
                    messg_size_doubles*sizeof(double));
      (void) printf("\t\ttotal number of split messages to be sent: %d\n\n",
                    number_of_messages);
    }

    /*
     * Allocate a temporary send buffer that can hold all out going messages.
     * Then copy all info to this buffer.
     */

    total_send_size = 0;
    for (n = 0; n < num_neighbors; n++) {
      total_send_size += actual_send_length[n];
    }

    send_buffer = (double *) calloc(total_send_size+1, sizeof(double));
    if (send_buffer == NULL) {
      (void) fprintf(stderr, "no space in AZ_splitup_big_msg: send_buffer \n");
      exit(-1);
    }
    for (n = 0; n < total_send_size ; n++) {
      send_buffer[n] = buffer[n];
    }

    /*
     * Send and receive messages in a series of communications. Each set of
     * exchanges is followed by a syncronization to not allow message buffers to
     * overflow.
     */

    doubles_sent = 0;

    for (m = 0; m < number_of_messages; m++) {

      /* post recieves for split messages */

      for (n = 0; n < num_neighbors; n++) {

        total_doubles_to_recv = actual_recv_length[n];
        messg_from            = proc_num_neighbor[n];
        dummy_int             = type;

        if (doubles_sent + messg_size_doubles < total_doubles_to_recv ) {


          /* read messg_size_doubles bytes */

          length = messg_size_doubles*sizeof(double);

          char_ptr = (char *) (&buffer[start_recv_proc[n]] + doubles_sent);



          (void) md_wrap_iread((void *) char_ptr, length, &messg_from, 
                               &dummy_int,  request+n);

        }

        else if (doubles_sent+messg_size_doubles >= total_doubles_to_recv &&
                 finished_recv_messg[n] == AZ_FALSE) {

          /* read actual_recv_length[n] - doubles_sent bytes */

          length = (total_doubles_to_recv - doubles_sent)*sizeof(double);

          char_ptr = (char *) (&buffer[start_recv_proc[n]] + doubles_sent);

          (void) md_wrap_iread((void *) char_ptr, length, &messg_from, 
                               &dummy_int,  request+n);

        }

        else if (finished_recv_messg[n] == AZ_TRUE) {

          /* read integer dummy message */

          length = sizeof(int);

         (void) md_wrap_iread((void *) &dummy_add, length, &messg_from, 
                               &dummy_int,  request+n);

        }
      }

      /* write split messages */

      for (n = 0; n < num_neighbors; n++) {

        total_doubles_to_send = actual_send_length[n];
        dest                  = proc_num_neighbor[n];

        if (doubles_sent + messg_size_doubles < total_doubles_to_send) {

          /* send out messg_size_doubles bytes */

          length = messg_size_doubles*sizeof(double);

          char_ptr = (char *) (&send_buffer[start_send_proc[n]] + doubles_sent);



          (void) md_wrap_write((void *) char_ptr, length, dest, type, &flag);

        }

        else if (doubles_sent + messg_size_doubles >= total_doubles_to_send &&
                 finished_send_messg[n] == AZ_FALSE) {

          /* send out actual_send_length[n] - doubles_sent bytes */

          length = (total_doubles_to_send - doubles_sent)*sizeof(double);

          char_ptr = (char *) (&send_buffer[start_send_proc[n]] + doubles_sent);

          (void) md_wrap_write((void *) char_ptr, length, dest, type, &flag);

          finished_send_messg[n] = AZ_TRUE;
        }

        else if (finished_send_messg[n] == AZ_TRUE) {

          /* send out integer dummy message */

          length = sizeof(int);

          (void) md_wrap_write((void *) &dummy_add, length, dest, type, &flag);

        }
      }

      /* read split messages */


      for (n = 0; n < num_neighbors; n++) {

        total_doubles_to_recv = actual_recv_length[n];
        messg_from            = proc_num_neighbor[n];
        messg_type            = type;

        if (doubles_sent + messg_size_doubles < total_doubles_to_recv ) {

          /* read messg_size_doubles bytes */

          length = messg_size_doubles*sizeof(double);

          char_ptr = (char *) (&buffer[start_recv_proc[n]] + doubles_sent);



          
          size =  md_wrap_wait((void *) char_ptr, length, &messg_from,
                             &messg_type, &flag, request+n); 

          if (length != size) {
           (void) fprintf(stderr,"%sERROR on node %d\nmd_wait failed, message "
                          "type = %d\n", yo, proc_config[AZ_node], messg_type);
           exit(-1);
          }

        }

        else if (doubles_sent+messg_size_doubles >= total_doubles_to_recv &&
                 finished_recv_messg[n] == AZ_FALSE) {

          /* read actual_recv_length[n] - doubles_sent bytes */

          length = (total_doubles_to_recv - doubles_sent)*sizeof(double);

          char_ptr = (char *) (&buffer[start_recv_proc[n]] + doubles_sent);

          size =  md_wrap_wait((void *) char_ptr, length, &messg_from,
                             &messg_type, &flag, request+n); 

          if (length != size) {
           (void) fprintf(stderr,"%sERROR on node %d\nmd_wait failed, message "
                          "type = %d\n", yo, proc_config[AZ_node], messg_type);
           exit(-1);
          }

          finished_recv_messg[n] = AZ_TRUE;
        }

        else if (finished_recv_messg[n] == AZ_TRUE) {

          /* read integer dummy message */

          length = sizeof(int);

          size =  md_wrap_wait((void *) &dummy_add, length, &messg_from,
                             &messg_type, &flag, request+n); 

          if (length != size) {
           (void) fprintf(stderr,"%sERROR on node %d\nmd_wait failed, message "
                          "type = %d\n", yo, proc_config[AZ_node], messg_type);
           exit(-1);
          }

        }
      }

      doubles_sent += messg_size_doubles;


      AZ_sync(proc_config[AZ_node], proc_config[AZ_N_procs] );
    }

    free(send_buffer);
    return;
  }

  else {
     
 
  

     /* Test to see if we can allocate enough space for send_buffer 
        to send the entire message at once                            */

     total_send_size = 0;
     for (n = 0; n < num_neighbors; n++) {
        total_send_size += actual_send_length[n];
     }
     send_buffer = (double *) calloc(total_send_size+1, sizeof(double));
     if (send_buffer == NULL) {
        (void) fprintf(stderr,"no space AZ_splitup_big_msg: send_buffer \n");
        exit(-1);
     }
    
   
     for (n = 0; n < total_send_size ; n++) {
        send_buffer[n] = buffer[n];
     }
     
     /* post receives for entire message from all neighbors */
     
     j = 0;
     for (n = 0; n < num_neighbors; n++) {

        messg_from = proc_num_neighbor[n];
        dummy_int = type;
        size      = actual_recv_length[n]*sizeof(double);

     
 
        
      
        (void) md_wrap_iread((void *) &buffer[j], size, 
                             &messg_from, &dummy_int, request+n);
        j += actual_recv_length[n];
        
     }


    /* send out entire message to each neighbor */

     for (n = 0; n < num_neighbors; n++) {

        size = actual_send_length[n]*sizeof(double);



        (void) md_wrap_write((void *) &send_buffer[start_send_proc[n]],
                             size, proc_num_neighbor[n], type, &st);

     }             


     /* wait for all messages */

    j = 0;
    for (n = 0; n < num_neighbors; n++) {
      messg_from = proc_num_neighbor[n];
      rtype     = type;
      size      = actual_recv_length[n]*sizeof(double);

      length =  md_wrap_wait((void *) &buffer[j], size, &messg_from,
                             &rtype, &st, request+n); 

      if ((length != size) && (size !=0) ) {
        (void) fprintf(stderr, "%sERROR on node %d\nmd_wait failed, message "
                       "type = %d\n", yo, proc_config[AZ_node] , rtype);
        exit(-1);
      }



      j += length / (sizeof(double));
    }

    *total_num_recv = j;

    free(send_buffer);
  }

} /* AZ_splitup_big_msg */
