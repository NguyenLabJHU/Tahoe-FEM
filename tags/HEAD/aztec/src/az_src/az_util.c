/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile: az_util.c,v $
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
static char rcsid[] = "$Id: az_util.c,v 1.1.1.1 2001-01-30 20:59:12 paklein Exp $";
#endif


/*******************************************************************************
 * Copyright 1995, Sandia Corporation.  The United States Government retains a *
 * nonexclusive license in this software as prescribed in AL 88-1 and AL 91-7. *
 * Export of this program may require a license from the United States         *
 * Government.                                                                 *
 ******************************************************************************/

#include "f2c.h"
#undef abs
#include "az_aztec.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifndef __MWERKS__
#include <malloc.h>
#endif /* __MWERKS__ */
#include <float.h>

#ifdef __MWERKS__
#include "utilsx.h" /* non-standard ANSI-like functions */
#endif /* __MWERKS__ */

/*
 * File containing utility functions for solvers.  Note: Some of the fem high
 * level solvers such as the nonlinear solver call these routines as well.
 */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_compute_residual(double val[], int indx[], int bindx[], int rpntr[],
                         int cpntr[], int bpntr[], double b[], double x[],
                         double r[], int data_org[])

/*******************************************************************************

  Compute the residual r = b - Ax.

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  val:             Array containing the nonzero entries of the matrix (see file
                   params.txt).

  indx,
  bindx,
  rpntr,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file params.txt).

  b:               Right hand side of linear system.

  x:               On input, contains the initial guess. On output contains the
                   solution to the linear system.

  r:               On output, residual vector.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see file params.txt).

*******************************************************************************/

{

  /* local variables */

  register int i;
  int N;

  /**************************** execution begins ******************************/

  N = data_org[AZ_N_internal] + data_org[AZ_N_border];

  AZ_matvec_mult(val, indx, bindx, rpntr, cpntr, bpntr, x, r, 1, data_org);

  for(i = 0; i < N; i++)
    r[i] = b[i] - r[i];

} /* AZ_compute_residual */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

double AZ_gmax_vec(int N, double vec[], int proc_config[])

/*******************************************************************************

  Routine to return the maximum element of a distributed vector "vec".

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     double, maximum value in vector 'vec'
  ============

  Parameter list:
  ===============

  N:               Length of vector 'vec'.

  vec:             Vector of length 'N'.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  /* local variables */

  register int i;
  double       rmax = 0.0;

  /**************************** execution begins ******************************/

  for (i = 0; i < N; i++) rmax = max(rmax, vec[i]);
  rmax = AZ_gmax_double(rmax, proc_config);

  return rmax;

} /* AZ_gmax_vec */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

double AZ_gdot(int N, double r[], double z[], int proc_config[])

/*******************************************************************************

  Routine to perform dot product of r and z with unit stride. This routine call
  the BLAS routine ddot to do the local vector dot product and then uses the
  global summation routine AZ_gsum_double to obtain the reguired global result.

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     double, dot product of vectors 'r' and 'z'
  ============

  Parameter list:
  ===============

  N:               Length of vector 'vec'.

  r, z:            Vectors of length 'N'.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_Procs] is the number of processors.

*******************************************************************************/

{

  static int one = 1;
  int        add_N;

  add_N = N;

  return AZ_gsum_double(ddot_(&add_N, r, &one, z, &one), proc_config);

} /* dot */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

#ifdef TIME_VB

void AZ_time_kernals(int gN, double gnzeros, double val[], int indx[],
                     int bindx[], int rpntr[], int cpntr[], int bpntr[],
                     double x[], double y[], int ntimes, int options[],
                     int data_org[], int proc_config[])

/*******************************************************************************

  Solve the system of equations given in the VBR format using an iterative
  method specified by 'options[AZ_solver]'. Store the result in 'x'.

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  gN:              Global order of the linear system of equations.

  gnzeros:         Global number of nonzeros in val.

  val:             Array containing the nonzero entries of the matrix (see file
                   params.txt).

  indx,
  bindx,
  rpntr,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file params.txt).

  x, y:            Vectors of length gN.

  ntimes:          Number of times to perform each operation.

  options:         Determines specific solution method and other parameters.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see file params.txt).

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  /* local variables */

  register int i, j;
  double       start_t, total_t;
  double       sprx_comm_t, sprx_comp_t, dot_comm_t, dot_comp_t;
  double       sprx_t, sprx_overlap_cop_t;
  double       daxpy_comp_t, dzpax_comp_t;
  double       Djunk;
  double      *ptr_vec1, *ptr_vec2, *ptr_vec3, *z1, *z2, *z3, alpha = 1.0;
  int          one = 1, bpntr_index;
  double       Mflop, Mflops_comp, Mflops_node_comp;
  double       Mflops_comm, Mflops_node_comm;
  double       sparax_overlap_border_comp_t, sparax_overlap_internal_comp_t;
  double       read_nonoverlap_t, read_nonoverlap_t_max, read_t_max;
  double       time_overlap_max;
  double       gather_t, write_t, read_t, time;
  double       max, min, avg;
  double       overall_t, start_overall_t;

  int          Num_Proc, Proc;
  int          Num_Internal_Blks;
  int          iout;

  char        *message_recv_add[AZ_MAX_NEIGHBORS];
  char        *message_send_add[AZ_MAX_NEIGHBORS];
  int          message_recv_length[AZ_MAX_NEIGHBORS];
  int          message_send_length[AZ_MAX_NEIGHBORS];

  /*
    message_send_add:
    message_send_add[i] points to the beginning of the list of
    values to be sent to the ith neighbor (i.e. data_org[AZ_neighbors+i] or
    sometimes locally defined as proc_num_neighbor[i]). That is,
    *(message_send_add[i] + j) is the jth value to be sent to the
    ith neighbor.

    message_send_length:
    message_send_length[i] is the number of bytes to be sent to
    the ith neighbor (i.e. data_org[AZ_neighbors+i] or sometimes
    locally defined as proc_num_neighbor[i]).

    message_recv_add:
    message_recv_add[i] points to the beginning of the list of
    locations which are to receive values sent by the ith neighbor (i.e.
    data_org[AZ_neighbors+i] or sometimes locally defined as
    proc_num_neighbor[i]). That is, *(message_recv_add[i] + j) is the
    location where the jth value sent from the ith neighbor will be stored.
    message_recv_length:
    message_recv_length[i] is the number of bytes to be sent to
    the ith neighbor (i.e. data_org[AZ_neighbors+i] or sometimes
    locally defined as proc_num_neighbor[i]).
    */

  int          temp1, temp2;

  /**************************** execution begins ******************************/

  Proc              = proc_config[AZ_node];
  Num_Proc          = proc_config[AZ_N_procs];
  Num_Internal_Blks = data_org[AZ_N_int_blk];
  iout              = options[AZ_output];

  if (data_org[AZ_matrix_type] != AZ_VBR_MATRIX) {
    (void) fprintf(stderr, "I'm not sure if the timing stuff works for \n"
                   "nonVBR matrices. For example I don't think that\n"
                   "we have overlapped communication/computations.\n"
                   "However, probably a lot of the stuff works?\n");
    exit(-1);
  }

  if (ntimes < 0)   ntimes = -ntimes;
  if (ntimes < 200) ntimes = 200;

  if (Proc == 0) {
    (void) printf("++++++++++ timing run ++++++++++\n");
    (void) printf("\nnprocs: %d ntimes: %d\n", Num_Proc, ntimes);
    (void) printf("N: %d\tNzeros: %e\t\n\n", gN, gnzeros);
  }

  /* sparax */

  if (iout > 0 && Proc == 0) {
    (void) printf("timing nonoverlapped matrix vector multiply\n");
  }

  start_t = AZ_second();
  for (i = 1; i <= ntimes; i++)
    AZ_matvec_mult(val, indx, bindx, rpntr, cpntr, bpntr, x, y, 1, data_org);

  sprx_t = AZ_second() - start_t;

  /* exchange boundary */

  if (iout > 0 && Proc ==0){
    (void) printf("timing unstructured communication\n");
  }

  start_t = AZ_second();
  for (i = 1; i <= ntimes; i++) {
    AZ_exchange_bdry(x, data_org);
    AZ_exchange_bdry(y, data_org);
  }

  sprx_comm_t = 0.5*(AZ_second() - start_t);
  sprx_comp_t = sprx_t - sprx_comm_t;

  /* overlapped communication */

  if (iout > 0 && Proc == 0) {
    (void) printf("timing overlapped sparse matrix vector multiply\n");
  }
  start_t = AZ_second();
  for (i = 1; i <= ntimes; i++)
    AZ_matvec_mult(val, indx, bindx, rpntr, cpntr, bpntr, x, y, 1, data_org);
  time = AZ_second() - start_t;

  gather_t                       = (double) 0.0;
  write_t                        = (double) 0.0;
  read_t                         = (double) 0.0;
  read_nonoverlap_t              = (double) 0.0;
  sparax_overlap_internal_comp_t = (double) 0.0;
  sparax_overlap_border_comp_t   = (double) 0.0;

  if (iout > 0 && Proc == 0) {
    (void) printf("time the individual routines for sparax_overlap\n");
  }

  start_overall_t = AZ_second();

  for (i = 1; i <= ntimes; i++) {

    /* time the individual routines for sparax_overlap */

    start_t = AZ_second();
    AZ_gather_mesg_info(x, data_org, message_recv_add, message_send_add,
                        message_recv_length, message_send_length);
    gather_t += AZ_second() - start_t;

    start_t = AZ_second();
    AZ_write_local_info(data_org, message_recv_add, message_send_add,
                        message_recv_length, message_send_length);
    write_t += AZ_second() - start_t;

    start_t = AZ_second();

    temp1                      = data_org[AZ_N_bord_blks];
    temp2                      = data_org[AZ_N_border];
    data_org[AZ_N_bord_blks] = 0;
    data_org[AZ_N_border] = 0;

    AZ_matvec_mult(val, indx, bindx, rpntr, cpntr, bpntr, x, y, 0, data_org);

    data_org[AZ_N_bord_blks] = temp1;
    data_org[AZ_N_border] = temp2;

    sparax_overlap_internal_comp_t += AZ_second() - start_t;

    start_t = AZ_second();
    AZ_read_local_info(data_org, message_recv_add, message_recv_length);
    read_t += AZ_second() - start_t;

    start_t = AZ_second();

    /* compute boundary portion of the sparse matrix - vector product */

    bpntr_index = bpntr[Num_Internal_Blks];

    temp1                      = data_org[AZ_N_int_blks];
    temp2                      = data_org[AZ_N_internal];
    data_org[AZ_N_int_blks]    = data_org[AZ_N_bord_blks];
    data_org[AZ_N_internal]    = data_org[AZ_N_border];
    data_org[AZ_N_bord_blks]   = 0;
    data_org[AZ_N_border]      = 0;

    AZ_matvec_mult(&val[indx[bpntr_index]], &indx[bpntr_index],
                   &bindx[bpntr_index], &rpntr[Num_Internal_Blks], cpntr,
                   &bpntr[Num_Internal_Blks], x, &y[rpntr[Num_Internal_Blks]],
                   0, data_org);

    data_org[AZ_N_bord_blks] = data_org[AZ_N_int_blks];
    data_org[AZ_N_border]    = data_org[AZ_N_internal];
    data_org[AZ_N_int_blks]  = temp1;
    data_org[AZ_N_internal]  = temp2;

    sparax_overlap_border_comp_t += AZ_second() - start_t;
  }

  overall_t = AZ_second() - start_overall_t;

  for (i = 1; i <= ntimes; i++) {

    /* time the reads in the nonoverlapped case */

    AZ_gather_mesg_info(x, data_org, message_recv_add, message_send_add,
                        message_recv_length, message_send_length);

    AZ_write_local_info(data_org, message_recv_add, message_send_add,
                        message_recv_length, message_send_length);

    start_t = AZ_second();
    AZ_read_local_info(data_org, message_recv_add, message_recv_length);
    read_nonoverlap_t += AZ_second() - start_t;
  }

  read_t_max            = AZ_gmax_double(read_t, proc_config);
  read_nonoverlap_t_max = AZ_gmax_double(read_nonoverlap_t, proc_config);
  time_overlap_max      = AZ_gmax_double(time, proc_config);

  /* dot */

  if (iout > 0 && Proc == 0) {
    (void) printf("time the individual routines for ddot\n");
  }

  start_t = AZ_second();
  for (i = 1; i <= ntimes; i++) Djunk = AZ_gdot(N, x, y, proc_config);
  total_t = AZ_second() - start_t;

  start_t = AZ_second();
  for (i = 1; i <= ntimes; i++) AZ_gsum_double(Djunk, proc_config);

  dot_comm_t = AZ_gmax_double(AZ_second() - start_t, proc_config);
  dot_comp_t = AZ_gmax_double(total_t - dot_comm_t, proc_config);

  /* daxpy */

  if (iout > 0 && Proc == 0) {
    (void) printf("time the individual routines for daxpy\n");
  }
  start_t = AZ_second();
  for (i = 1; i <= ntimes; i++) {
    Djunk = (double) i;
    daxpy_(&N, &Djunk, x, &one, y, &one);
  }
  daxpy_comp_t = AZ_gmax_double(AZ_second() - start_t, proc_config);

  if (Proc == 0) {              /* calculate and print results */
    (void) printf("\n********** sparax ***********\n\n");
    (void) printf("nonoverlapped\n");
    (void) printf("\t\tmax\t\tavg\t\tmin\n");
  }

  /* dzpax */

  if (iout > 0 && Proc == 0) {
    (void) printf("time the individual routines for dzpax\n");
  }

  z1 = (double *) malloc(N * sizeof(double));
  z2 = (double *) malloc(N * sizeof(double));
  z3 = (double *) malloc(N * sizeof(double));
  for (i = 0; i < N; i++) {
    z2[i] = (double) i;
    z3[i] = (double) 2*i;
  }

  start_t = AZ_second();
  for (i = 1; i <= ntimes; i++) {
    Djunk = 3.14 * (double) i;
    for (j = 0; j < N; j++)
      z1[j] = z2[j] + Djunk*z3[j];
    /*    dzpax_(&N, &Djunk, z3, &one, z2, &one, z1, &one);*/
  }
  dzpax_comp_t = AZ_gmax_double(AZ_second() - start_t, proc_config);

  free((void *) z1);
  free((void *) z2);
  free((void *) z3);

  if (Proc == 0) {              /* calculate and print results */
    (void) printf("\n********** sparax ***********\n\n");
    (void) printf("nonoverlapped\n");
    (void) printf("\t\tmax\t\tavg\t\tmin\n");
  }

  max = AZ_gmax_double(sprx_t, proc_config);
  avg = AZ_gavg_double(sprx_t, proc_config);
  min = AZ_gmin_double(sprx_t, proc_config);
  if (Proc == 0) (void) printf("total_time\t%e\t%e\t%e\n", max, avg, min);

  max = AZ_gmax_double(sprx_comp_t, proc_config);
  avg = AZ_gavg_double(sprx_comp_t, proc_config);
  min = AZ_gmin_double(sprx_comp_t, proc_config);
  if (Proc == 0) (void) printf("comp_time\t%e\t%e\t%e\n", max, avg, min);

  max = AZ_gmax_double(sprx_comm_t, proc_config);
  avg = AZ_gavg_double(sprx_comm_t, proc_config);
  min = AZ_gmin_double(sprx_comm_t, proc_config);
  if (Proc == 0) (void) printf("comm_time\t%e\t%e\t%e\n", max, avg, min);

  sprx_comp_t = AZ_gmax_double(sprx_comp_t, proc_config);
  sprx_t      = AZ_gmax_double(sprx_t, proc_config);

  if (Proc == 0) {
    Mflop             = (double) ntimes * (gnzeros + gnzeros) * 1.0e-6;
    Mflops_comp       = Mflop/sprx_comp_t;
    Mflops_node_comp  = Mflops_comp/(double) Num_Proc;
    Mflops_comm       = Mflop/sprx_t;
    Mflops_node_comm  = Mflops_comm/(double) Num_Proc;

    (void) printf("computation Mflops: %e \n", Mflops_comp);
    (void) printf("computation Mflops per node: %e \n", Mflops_node_comp);
    (void) printf("comp & comm Mflops: %e \n", Mflops_comm);
    (void) printf("comp & comm Mflops per node: %e \n\n", Mflops_node_comm);
  }

  /* statistics */

  if (Proc == 0) {
    (void) printf("overlapped\n\t\tmax\t\tavg\t\tmin\n");
  }

  max = AZ_gmax_double(time, proc_config);
  min = AZ_gmin_double(time, proc_config);
  avg = AZ_gavg_double(time, proc_config);
  if (Proc == 0) (void) printf("total_time\t%e\t%e\t%e\n", max, avg, min);

  max = AZ_gmax_double(gather_t, proc_config);
  min = AZ_gmin_double(gather_t, proc_config);
  avg = AZ_gavg_double(gather_t, proc_config);
  if (Proc == 0) (void) printf("gather_t\t%e\t%e\t%e\n", max, avg, min);

  max = AZ_gmax_double(write_t, proc_config);
  min = AZ_gmin_double(write_t, proc_config);
  avg = AZ_gavg_double(write_t, proc_config);
  if (Proc == 0) (void) printf("write_t \t%e\t%e\t%e\n", max, avg, min);

  max = AZ_gmax_double(sparax_overlap_internal_comp_t, proc_config);
  min = Az_gmin_double(sparax_overlap_internal_comp_t, proc_config);
  avg = AZ_gavg_double(sparax_overlap_internal_comp_t, proc_config);
  if (Proc == 0) (void) printf("internal_t\t%e\t%e\t%e\n", max, avg, min);

  max = AZ_gmax_double(read_t, proc_config);
  min = Az_gmin_double(read_t, proc_config);
  avg = AZ_gavg_double(read_t, proc_config);
  if (Proc == 0) (void) printf("read_t   \t%e\t%e\t%e\n", max, avg, min);

  max = AZ_gmax_double(sparax_overlap_border_comp_t, proc_config);
  min = AZ_gmin_double(sparax_overlap_border_comp_t, proc_config);
  avg = AZ_gavg_double(sparax_overlap_border_comp_t, proc_config);
  if (Proc == 0) (void) printf("border_t\t%e\t%e\t%e\n", max, avg, min);

  if (Proc == 0) {
    Mflops_comm      = Mflop/time_overlap_max;
    Mflops_node_comm = Mflops_comm/(double) Num_Proc;
    (void) printf("comp & comm Mflops: %e \n", Mflops_comm);
    (void) printf("comp & comm Mflops per node: %e \n\n", Mflops_node_comm);

    (void) printf("Ratio of overlapped/nonoverlapped sparax times: %e\n\n",
                  time_overlap_max/sprx_t);
    (void) printf("Ratio of overlapped/nonoverlapped read times: %e\n\n",
                  read_t_max/read_nonoverlap_t_max);
  }

  max = AZ_gmax_double(time, proc_config);
  if (time  ==  max) {
    (void) printf("max time proc\n");
    (void) printf("\nProc:%d\tNum_Neighbors: %d\n", Proc, data_org[AZ_N_neigh]);
    (void) printf("total_time: %e\n", time);
    (void) printf("gather_t  : %e\n", gather_t);
    (void) printf("write_t   : %e\n", write_t);
    (void) printf("internal_t: %e\n", sparax_overlap_internal_comp_t);
    (void) printf("read_t    : %e\n", read_t);
    (void) printf("border_t  : %e\n\n", sparax_overlap_border_comp_t);
  }

  max = AZ_gmax_double(overall_t, proc_config);
  if (overall_t  ==  max) {
    (void) printf("overall max time proc\n");
    (void) printf("\nProc:%d\tNum_Neighbors: %d\n", Proc, data_org[AZ_N_neigh]);
    (void) printf("overall total_time: %e\n", overall_t);
    (void) printf("total_time: %e\n", time);
    (void) printf("gather_t  : %e\n", gather_t);
    (void) printf("write_t   : %e\n", write_t);
    (void) printf("internal_t: %e\n", sparax_overlap_internal_comp_t);
    (void) printf("read_t    : %e\n", read_t);
    (void) printf("border_t  : %e\n\n", sparax_overlap_border_comp_t);
  }
  /*
    if (Proc == 0) {
    (void) printf("cop overlapped\n");
    (void) printf("\t\tmax\t\tavg\t\tmin\n");
    }

    max = AZ_gmax_double(sprx_overlap_cop_t, proc_config);
    min = AZ_gmin_double(sprx_overlap_cop_t, proc_config);
    avg = AZ_gavg_double(sprx_overlap_cop_t, proc_config);
    if (Proc == 0) (void) printf("total_time\t%e\t%e\t%e\n", max, avg, min);

    if (Proc == 0) {
    Mflops_comm = Mflop/max;
    Mflops_node_comm  = Mflops_comm/(double)Num_Proc;
    (void) printf("comp & comm Mflops: %e \n", Mflops_comm);
    (void) printf("comp & comm Mflops per node: %e \n\n", Mflops_node_comm);
    }
    */
  if (Proc == 0) {
    (void) printf("\n********* exchange **********\n\n");
    (void) printf("comm_time = %7.4e sec.\n", sprx_comm_t);

    (void) printf("\n************ dot ************\n\n");
    Mflop = (double) ntimes * 2.0 * ((double) gN * 1.0e-06);

    (void) printf("comp_time = %7.4e sec.\n", dot_comp_t);
    (void) printf("comm_time = %7.4e sec.\n", dot_comm_t);

    Mflops_comp      = Mflop/dot_comp_t;
    Mflops_node_comp = Mflops_comp/(double) Num_Proc;
    Mflops_comm      = Mflop/(dot_comm_t + dot_comp_t);
    Mflops_node_comm = Mflops_comm/(double) Num_Proc;

    (void) printf("computation Mflops: %e \n", Mflops_comp);
    (void) printf("computation Mflops per node: %e \n", Mflops_node_comp);
    (void) printf("comp & comm Mflops: %e \n", Mflops_comm);
    (void) printf("comp & comm Mflops per node: %e \n\n", Mflops_node_comm);

    (void) printf("\n*********** daxpy ***********\n\n");
    (void) printf("comp_time = %7.4e sec.\n", daxpy_comp_t);
    Mflops_comp      = Mflop/daxpy_comp_t;
    Mflops_node_comp = Mflops_comp/(double) Num_Proc;

    (void) printf("computation Mflops: %e \n", Mflops_comp);
    (void) printf("computation Mflops per node: %e \n", Mflops_node_comp);

    (void) printf("\n*********** dzpax ***********\n\n");
    (void) printf("comp_time = %7.4e sec.\n", dzpax_comp_t);
    Mflops_comp  = Mflop/dzpax_comp_t;
    Mflops_node_comp  = Mflops_comp/(double) Num_Proc;

    (void) printf("computation Mflops: %e \n", Mflops_comp);
    (void) printf("computation Mflops per node: %e \n", Mflops_node_comp);
  }

} /* AZ_time_kernals */

#endif

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int AZ_find_index(int key, int list[], int length)

/*******************************************************************************

  Find 'key' in 'list' and return the index number.

  Author:          Ray Tuminaro, SNL, 1422
  =======

  Return code:     int, -1 = key not found, i = list[i] = key
  ============

  Parameter list:
  ===============

  key:             Element to be search for in list.

  list:            List to be searched.

  length:          Length of list.

*******************************************************************************/

{

  /* local variables */

  int start, end;
  int mid;

  /**************************** execution begins ******************************/

  if (length == 0) return -1;

  start = 0;
  end   = length - 1;

  while (end - start > 1) {
    mid = (start + end) / 2;
    if (list[mid] < key) start = mid;
    else end = mid;
  }

  if (list[start] == key) return start;
  if (list[end] == key)   return end;
  return -1;

} /* AZ_find_index */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_free_memory(int label)

/*******************************************************************************

  Clear all internal memory that has been allocated by previous calls to
  AZ_manage_memory(*,*,label,*,*). This memory is typically includes
  preconditioner and scaling information.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  label:           An integer associated with memory that was allocated with
                   AZ_manage_memory(). On ouput, all memory allocated with this
                   integer will be freed.

*******************************************************************************/

{
  (void) AZ_manage_memory((int) NULL, AZ_CLEAR, label, (char *) NULL,
                          (int *) NULL);
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

double *AZ_manage_memory(int size, int action, int type, char *name,
                         int *status)

/*******************************************************************************

  AZ_manage_memory() either frees memory that was previously allocated or
  returns a pointer to memory for future use (this pointer can be a newly
  allocated block or a block that was previously allocated).


  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     double *,  pointer to allocated memory.
  ============

  Parameter list:
  ===============

  size:            integer variable. On input, size indicates the amount of
                   memory that is needed (NOTE: ignored if action == AZ_CLEAR).

  action:          On input, action indicates whether to allocate or free
                   memory.
                   = AZ_ALLOC: look for a chunk of memory (in the memory
                               management list pointed to by 'head') that has
                               already been allocated with a memory-type 'type',
                               a size in bytes of 'size', and labelled with the
                               string 'name'.

                               If this memory is found return a pointer to this
                               memory and set *status to AZ_OLD_ADDRESS.
                               If this memory is not found allocate the memory,
                               put it in the memory management list, return a
                               pointer to this piece of memory, and set *status
                               to AZ_NEW_ADDRESS.

                   = AZ_REALLOC:look for a chunk of memory (in the memory
                               management list pointed to by 'head') that has
                               already been allocated with a memory-type 'type',
                               and labelled with the string 'name'. Reallocate
                               this item with size 'size' and change the size
                               field.

                   = AZ_CLEAR: free all memory associated with any chunk of
                               memory that has the memory-type 'type'.


  type:            On input, a number associated with each piece of memory that
                   is allocated.  This number is used when checking to see if
                   a requested piece of memory is already allocated. This number
                   is also used when freeing up memory (see AZ_CLEAR).
                   NOTE: Aztec uses the type AZ_SYS for temporary internal
                         pieces of memory that can be deallocated after solving
                         a linear system.

  name:            On input, a character string associated with each piece of
                   memory that is allocated. This string is used when checking
                   to see if a requested piece of memory is already allocated.
                   Additionally, this string is printed when AZ_manage_memory()
                   fails to allocate memory.

  status:          On output,
                   *status = AZ_NEW_ADDRESS indicates that a new piece of memory
                                            has been allocated.
                   *status = AZ_OLD_ADDRESS indicates that a pointer to a
                                            previously allocated chunk of memory
                                            is returned.
                   NOTE: status is not used if action == AZ_CLEAR has been
                         allocated.


*******************************************************************************/

{

  /* local variables */

  struct mem_ptr {
    char   *name;
    double *address;
    int     type;
    int     size;
    struct mem_ptr *next;
  };

  static struct mem_ptr *head = NULL;
  struct mem_ptr        *current, *temp,*prev;
  int                   found = 0;

  /**************************** execution begins ******************************/

  current = head;

  if (action == -43) {

    /* print the list */

    while (current != NULL) {
      (void) printf("(%8d,%8d) ==> %s\n", current->type, current->size,
                    current->name);
      prev    = current;
      current = current->next;
    }
    return (double *) 0;
  }

  else if (action == AZ_ALLOC) {
    if (size == 0) return (double *) 0;

    /* first look for entry */

    while ( current != NULL) {
      if ( (current->size == size) && (current->type == type) &&
           (strcmp(current->name,name) == 0) ) {
        found = 1;
        break;
      }
      prev    = current;
      current = current->next;
    }

    if (found == 0) {

      /*
       * Put in a new entry if it has the type 0 we will put it at the front of
       * the list else put it at the end of the list This is done for efficiency
       * reasons as type = 0 is used a lot.
       */

      temp = (struct mem_ptr *) malloc(sizeof(struct mem_ptr));
      if (temp == NULL) {
        (void) fprintf(stderr, "Error: Not enough space to allocate\n");
        (void) fprintf(stderr, "       '%s' in memory_manage()\n", name);
        exit(-1);
      }

      temp->name    = (char *) az_strdup(name);
      temp->address = (double *) malloc(size);
      if (temp->address == NULL) {
        (void) fprintf(stderr, "Error:Not enough memory for '%s'\n", name);
        (void) fprintf(stderr, "      Asked for %d bytes. Perhaps\n", size);
        (void) fprintf(stderr, "      a smaller problem should be run\n");
        exit(-1);
      }

      temp->type = type;
      temp->size = size;
      if ( (type == AZ_SYS) || (head == NULL)) {
        temp->next = head;
        head       = temp;
      }
      else {
        prev->next = temp;
        temp->next = NULL;
      }

      *status = AZ_NEW_ADDRESS;
      return temp->address;
    }
    else {
      *status = AZ_OLD_ADDRESS;
      return current->address;
    }
  }

  else if (action == AZ_CLEAR) {
    prev = NULL;

    while (current != NULL) {
      if (current->type == type) {
        if (prev == NULL) head       = current->next;
        else              prev->next = current->next;

        temp = current->next;
        free(current->address);
        free(current->name);
        free(current);
        current = temp;
      }
      else {
        prev    = current;
        current = current->next;
      }
    }
    return (double *) 0;
  }
  else if (action == AZ_REALLOC) {
    if (size == 0)
      return (double *) 0;

    /* first look for entry */

    while ( current != NULL) {
      if ( (current->type == type) &&
           (strcmp(current->name,name) == 0) ) {
        found = 1;
        break;
      }
      prev    = current;
      current = current->next;
    }
    if (current == NULL) {
      (void) fprintf(stderr, "memory_management error: %s with type %d not",
                     name, type);
      (void) fprintf(stderr, " found during reallocation\n");
      exit(-1);
    }
    current->size = size;
    *status = AZ_OLD_ADDRESS;
    current->address = (double *) realloc((char *) current->address,size);
    if (current->address == NULL) {
      (void) fprintf(stderr, "Error:Not enough memory for '%s'\n", name);
      (void) fprintf(stderr, "      Asked for %d bytes. Perhaps\n", size);
      (void) fprintf(stderr, "      a smaller problem should be run\n");
      exit(-1);
    }
    return current->address;

  }

  else {
    (void) fprintf(stderr, "Error: Invalid action(%d) in AZ_manage_memory()\n",
                   action);
    exit(-1);
  }
  return((double *) NULL);

} /* AZ_manage_memory */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_mysleep(int i)

/*******************************************************************************

  A dummy routine that just wastes some time.  We essentially try and mimick the
  unix sleep() function which is not available on the nCUBE.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  i:               On input, indicates how much time to spend in this routine.
                   The larger the value of 'i', the more time is spent in this
                   routine. On the Ncube 2, this usually corresponds to about i
                   seconds.

*******************************************************************************/

{

  /* local variables */

  int    j, k, top, inner;
  double a;

  /**************************** execution begins ******************************/

  a = 10.01;
  top = 700 * i;
  inner = 10 * 10 * 10;

  for (k = 0; k < top; k++)
    for (j = 0; j < inner; j++)
      a = 2.0 / (a - 1);

} /* AZ_mysleep */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_p_error(char *str, int proc)

/*******************************************************************************

  Print the string if we are processor 0.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  str:             String to be printed.

  proc:            Current processor number.

*******************************************************************************/

{
  if (proc == 0) (void) fprintf(stdout, "%s", str);
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int AZ_get_new_eps(double *epsilon, double recursive, double actual,
                   int proc_config[])

/*******************************************************************************

  Routine which decides what to do when the computed residual has converged but
  not the true residual.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  epsilon:         Residual tolerance.

  recursive:       On input, the norm of the residual produced by an iterative
                   method via recursion (e.g. r = r_0 + delta_1 + delta_2 +
                   delta_3 ... )

  actual:          On input, the norm of the explicitly computed residual (i.e.
                   r = b - A x ). In the absence of rounding error, recursive
                   and actual should be the same value. However, due to rounding
                   errors these two might differ.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  /* local variables */

  double difference;

  /**************************** execution begins ******************************/

  difference = fabs(actual - recursive);
  if (difference > *epsilon)
    return AZ_QUIT;
  else {
    *epsilon = *epsilon - 1.5 * difference;
    while (*epsilon < 0.0) *epsilon += .1*difference;
   }

  if (proc_config[AZ_node] == 0)
    (void) printf("\n\t\tTrying to reduce actual residual "
                  "further\n\t\t     (recursive = %e, actual = %e)\n\n",
                  recursive, actual);

  return AZ_CONTINUE;

} /* AZ_get_new_eps */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_terminate_status_print(int situation, int iter, double status[],
                               double rec_residual, double params[],
                               double scaled_r_norm, double actual_residual,
                               int options[], int proc_config[])

/*******************************************************************************

  Routine to output conditions under which iterative solver was terminated if
  other than specified convergence.

  Author:          Scott A. Hutchinson, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  situation:       On input, indicates why iterative method terminated.
                   = AZ_normal   : iterative solver terminates normally (i.e.
                                   convergence criteria met).
                   = AZ_breakdown: iterative solver has broken down.
                   = AZ_loss     : iterative solver has terminated due to a lack
                                   of accuracy in the recursive residual
                                   (caused by rounding errors).
                   = AZ_maxits   : iterative solver has reached the maximum
                                   number of requested iterations without
                                   convergence.
                   = AZ_ill_cond : the upper triangular factor of the
                                   Hessenberg within GMRES is ill-conditioned.
                                   There is a possibility that the matrix
                                   is singular and we have the least-squares
                                   solution.

  iter:            Number of iterations completed.

  status:          !!!!!THIS MAY NOT BE ACCURATE FOR THIS ROUTINE!!!!!
                   On output, indicates termination status:
                    0:  terminated normally.
                   -1:  maximum number of iterations taken without achieving
                        convergence.
                   -2:  Breakdown. The algorithm can not proceed due to
                        numerical difficulties (usually a divide by zero).
                   -3:  Internal residual differs from the computed residual due
                        to a significant loss of precision.

  rec_residual:    On input, the norm of the residual produced by an iterative
                   method via recursion (e.g. r = r_0 + delta_1 + delta_2 +
                   delta_3 ... )

  params:          Drop tolerance and convergence tolerance info.

  scaled_r_norm:   Residual expression (requested by the user via
                   OPTIONS[AZ_conv] ... see user's guide and
                   AZ_compute_global_scalars()) which is used in printing and
                   compared with PARAMS[AZ_TOL] when determining convergence.

  actual_residual: On input, the norm of the explicitly computed residual (i.e.
                   r = b - A x ).  In the absence of rounding error, recursive
                   and actual should be the same value. However, due to rounding
                   errors these two might differ.

  options:         Determines specific solution method and other parameters.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  /* local variables */

  static int iterations = 0;
  char      *solver_name;
  int        solver_flag, conv_flag;
  double     eps;

  /**************************** execution begins ******************************/

  eps = params[AZ_tol];

  /* set status */

  if (scaled_r_norm < eps) situation = AZ_normal;

  status[AZ_its] = (double ) iter;
  status[AZ_why] = (double ) situation;
  status[AZ_r]   = actual_residual;

  if (actual_residual == -1.0) status[AZ_r] = rec_residual;

  status[AZ_rec_r]    = rec_residual;
  status[AZ_scaled_r] = scaled_r_norm;

  if (situation == AZ_normal) return; /* nothing else to do */

  if (options[AZ_output] ==  AZ_none) return;

  /* initialize */

  solver_flag = options[AZ_solver];
  conv_flag   = options[AZ_conv];
  if (!iterations) iterations = iter;

  switch (solver_flag) {
  case AZ_cg:
    solver_name = (char *) malloc(3*sizeof(char));
    (void) strcpy(solver_name, "cg");
    break;
  case AZ_gmres:
    solver_name = (char *) malloc(6*sizeof(char));
    (void) strcpy(solver_name, "gmres");
    break;
  case AZ_cgs:
    solver_name = (char *) malloc(4*sizeof(char));
    (void) strcpy(solver_name, "cgs");
    break;
  case AZ_tfqmr:
    solver_name = (char *) malloc(5*sizeof(char));
    (void) strcpy(solver_name, "tfqmr");
    break;
  case AZ_bicgstab:
    solver_name = (char *) malloc(9*sizeof(char));
    (void) strcpy(solver_name, "bicgstab");
    break;
  case AZ_symmlq:
    solver_name = (char *) malloc(7*sizeof(char));
    (void) strcpy(solver_name, "symmlq");
    break;
  case AZ_lu:
    solver_name = (char *) malloc(4*sizeof(char));
    (void) strcpy(solver_name, "lu");
    break;
  default:
    (void) fprintf(stderr,
                   "Error: invalid solver flag %d passed to terminate_status\n",
                   solver_flag);
  exit(-1);
  }

  if (proc_config[AZ_node] == 0) {
    (void) fprintf(stderr, "\n\n");
    (void) fprintf(stderr,"\t*************************************************"
                   "**************\n\n");

    switch (situation) {
    case AZ_ill_cond:
      (void) fprintf(stderr, "\tWarning: the GMRES Hessenberg matrix is "
                     "ill-conditioned.  This may \n\tindicate that the "
                     "application matrix is singular. In this case, GMRES\n"
                     "\tmay have a least-squares solution.\n");
    break;
    case AZ_breakdown:
      if (!solver_flag) {
        (void) fprintf(stderr, "\tWarning: direction vector is no longer "
                       "A-conjugate \n\tto previous vector but solution has "
                       "not converged.\n");
      }
      else {
        (void) fprintf(stderr, "\tWarning: a breakdown in this "
                       "method\n\thas occurred and solution has not "
                       "converged.\n");
      }
      break;

    case AZ_loss:
      (void) fprintf(stderr, "\tWarning: recursive residual indicates "
                     "convergence\n\tthough the true residual is too large.\n");
      (void) fprintf(stderr, "\n\tSometimes this occurs when storage is ");
      (void) fprintf(stderr, "overwritten (e.g. the\n\tsolution vector was not ");
      (void) fprintf(stderr, "dimensioned large enough to hold\n\texternal ");
      (void) fprintf(stderr, "variables). Other times, this is due to roundoff. ");
      (void) fprintf(stderr, "In\n\tthis case, the solution has either ");
      (void) fprintf(stderr, "converged to the accuracy\n\tof the machine or ");
      (void) fprintf(stderr, "intermediate roundoff errors ");
      (void) fprintf(stderr, "occurred\n\tpreventing full convergence. In the ");
      (void) fprintf(stderr, "latter case, try solving\n\tagain using the new ");
      (void) fprintf(stderr, "solution as an initial guess.\n");
      break;

    case AZ_maxits:
      (void) fprintf(stderr, "\tWarning: maximum number of iterations "
                     "exceeded\n\twithout convergence\n");
    break;

    default:
      (void) fprintf(stderr, "\tError: improper code passed from solver %s\n\n",
                     solver_name);
    (void) fprintf(stderr,"\t***********************************************",
                   "**********\n\n");
    exit(-1);
    }

    (void) fprintf(stdout,"\n\tSolver:\t\t\t%s\n", solver_name);
    (void) fprintf(stdout,"\tnumber of iterations:\t%d\n\n", iter);

    if (actual_residual != -1.0)
      (void) fprintf(stdout,"\tActual residual = %11.4e",actual_residual);

    (void) fprintf(stdout,"\tRecursive residual = %11.4e\n\n",rec_residual);
    (void) fprintf(stdout,"\tCalculated Norms\t\t\t\tRequested Norm\n");
    (void) fprintf(stdout,"\t--------------------------------------------");
    (void) fprintf(stdout,"\t--------------\n\n");

    switch (conv_flag) {
    case AZ_r0:
      (void) fprintf(stdout, "\t||r||_2 / ||r0||_2:\t\t%e\t%e\n", scaled_r_norm,
                     eps);
    break;
    case AZ_rhs:
      (void) fprintf(stdout, "\t||r||_2 / ||b||_2:\t\t%e\t%e\n", scaled_r_norm,
                     eps);
    break;
    case AZ_Anorm:
      (void) fprintf(stdout, "\t||r||_2 / ||A||_inf:\t\t%e\t%e\n",
                     scaled_r_norm, eps);
    break;
    case AZ_sol:
      (void) fprintf(stdout, "\t\t||r||_inf\n");
    (void) fprintf(stdout, "\t----------------------------- : %e\t%e\n",
                   scaled_r_norm, eps);
    (void) fprintf(stdout, "\t||A||_inf ||x||_1 + ||b||_inf\n");
    break;
    case AZ_weighted:
      (void) fprintf(stdout, "\t||r||_WRMS:\t\t%e\t%e\n", scaled_r_norm, eps);
    break;
    default:
      (void) fprintf(stderr, "terminate_status: ERROR: convergence test %d "
                     "not implemented\n", conv_flag);
    exit(-1);
    }

    (void) fprintf(stderr,"\n\t*********************************************"
                   "******************\n\n");
  }

  /* free memory */

  if (solver_name != NULL)
    free((void *) solver_name);

} /* AZ_terminate_status_print */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int AZ_breakdown_f(int N, double v[], double w[], double inner,
                   int proc_config[])

/*******************************************************************************

  Determine if the 2 vectors v and w are orthogonal. In particular,

       |<v,w>| <  100 ||v||_2 ||w||_2 DBL_EPSILON

  Author:          Ray Tuminaro, SNL, 1422
  =======

  Return code:     int
  ============

  Parameter list:
  ===============

  N:               Length of vectors v and w.

  v, w:            Two vectors of length N to be checked for othogonality.

  inner:           <v,w>

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  double v_norm, w_norm;

  v_norm = AZ_gvector_norm(N, 2, v, proc_config);
  w_norm = AZ_gvector_norm(N, 2, w, proc_config);

  return (fabs(inner) < 100.0 * v_norm * w_norm * DBL_EPSILON);

} /* breakdown */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_random_vector(double u[], int data_org[], int proc_config[])

/*******************************************************************************

  Set the the vector u to a random vector.

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  u:               Vector to be initialized.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see file params.txt).

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  /* local variables */

  static int seed = 493217272;
  int        i, N;

  /*********************** BEGIN EXECUTION *********************************/

  N    = (proc_config[AZ_node]+7) * (proc_config[AZ_node]+13) *
    (19 + proc_config[AZ_node]);

  seed = (int) (AZ_srandom1(&N)* (double) seed);
  N    = data_org[AZ_N_internal] + data_org[AZ_N_border];

  for (i = 0; i < N; i++) u[i] = AZ_srandom1(&seed);

} /* AZ_random_vector */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

double AZ_srandom1(int *seed)

/*******************************************************************************

  Random number generator.

  Author:
  =======

  Return code:     double, random number.
  ============

  Parameter list:
  ===============

  seed:            Random number seed.

*******************************************************************************/

{

  /* local variables */

  int    a = 16807;
  int    m = 2147483647;
  int    q = 127773;
  int    r = 2836;

  int    lo, hi, test;
  double rand_num;

  /**************************** execution begins ******************************/

  hi   = *seed / q;
  lo   = *seed % q;
  test = a * lo - r * hi;

  if (test > 0) *seed = test;
  *seed = test + m;

  rand_num = (double) *seed / (double) m;

  return rand_num;

} /* AZ_srandom1 */
