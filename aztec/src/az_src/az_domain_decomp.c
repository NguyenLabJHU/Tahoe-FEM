/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile: az_domain_decomp.c,v $
 *
 * $Author: paklein $
 *
 * $Date: 2001-01-30 20:59:14 $
 *
 * $Revision: 1.1.1.1 $
 *
 * $Name: not supported by cvs2svn $
 *====================================================================*/
#ifndef lint
static char rcsid[] = "$Id: az_domain_decomp.c,v 1.1.1.1 2001-01-30 20:59:14 paklein Exp $";
#endif


/*******************************************************************************
 * Copyright 1995, Sandia Corporation.  The United States Government retains a *
 * nonexclusive license in this software as prescribed in AL 88-1 and AL 91-7. *
 * Export of this program may require a license from the United States         *
 * Government.                                                                 *
 ******************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "az_aztec.h"

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_domain_decomp(double x[], double val[], int indx[], int rpntr[],
                   int cpntr[], int bindx[], int bpntr[], int options[],
                   int data_org[], int proc_config[], double params[])

/*******************************************************************************

  Driver for domain decomposition preconditioner.

  This routine preconditions x[] by backsolving 2 triangular systems. These
  triangular systems correspond to either an lu-decomposition, an
  ilu-decomposition, or a block ilu-decomposition. If
  options[AZ_pre_calc] < AZ_sys_reuse, it first computes the decomposition
  before performing the triangular solves. The specfic algorithm performed is
  determined by options[AZ_precond]. See file 'params.txt'.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  x:               On input, contains the residual(rhs) of the set of equations.
                   On output contains the preconditioned residual.

  val:             Array containing the nonzero entries of the matrix (see file
                   params.txt).

  indx,
  bindx,
  rpntr,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file params.txt).

  b:               Right hand side of linear system.

  options:         Determines specific solution method and other parameters.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see file params.txt).

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

  params:          Drop tolerance and convergence tolerance info.

*******************************************************************************/

{

  /* local variables */

  static int Nexp,      /* size of expanded system that will be factored */
             Mexp,      /* # of nonzeros in expanded system.             */
             new_blks;  /* # of blocks in the expanded system.           */

  double *buffer;       /* holds matrix rows received from neighbors
                           when overlapping domain decomposition is used */
  int    length;        /* length of buffer                              */
  int    N;             /* # of unknowns updated by this processor.      */
  int    i;

  double *little;
  int type,count,j, st, from, sym_flag, total;
  extern int AZ_sys_msg_type;
  MPI_Request request[AZ_MAX_NEIGHBORS];  /* Message handle */

  /**************************** execution begins ******************************/

  N = data_org[AZ_N_internal] + data_org[AZ_N_border];

  if (options[AZ_overlap] == AZ_sym_full) {
    options[AZ_overlap] = AZ_full;
    sym_flag = 1;
  }
  else sym_flag = 0;

  /* Set external points */

  if (options[AZ_overlap] == AZ_none)
    for (i = N; i < data_org[AZ_N_external] + N; i++) x[i] = 0.0;
  else AZ_exchange_bdry(x, data_org);

  /*
   * Compute the size and the number of nonzeros in the augmented system.
   * Note: the number of nonzeros is recomputed later if rows are exchanged
   */

  Nexp = N + data_org[AZ_N_external];

  if (N == 0)
    Mexp = 0;
  else if (data_org[AZ_matrix_type] == AZ_MSR_MATRIX)
    Mexp = data_org[AZ_N_external] + bindx[N];
  else
    Mexp = data_org[AZ_N_external] + indx[bpntr[data_org[AZ_N_int_blk] +
                                               data_org[AZ_N_bord_blk]]];

  /* Exchanging rows for the overlap points (if this option has been set) */


  if (options[AZ_pre_calc] <= AZ_recalc) {
    if (data_org[AZ_matrix_type] == AZ_VBR_MATRIX) {
      new_blks = 0;
      while (cpntr[new_blks] != Nexp) new_blks++;
    }

    if (options[AZ_overlap] == AZ_full) {
      AZ_exchange_rows(val,indx, &Mexp, Nexp, &buffer, &length, bindx, 
                       cpntr, bpntr, data_org, proc_config);
    }
  }
  if (options[AZ_precond] == AZ_ilu) {
    AZ_ilu_routine(Nexp, Mexp, x, length, buffer, val, indx, bindx, rpntr,
                   cpntr, bpntr, options, data_org);
  }
#ifdef eigen
  else if (options[AZ_precond] == AZ_slu) {
    AZ_lu_ng(Nexp, Mexp, x, val, bindx, options, data_org, params[AZ_drop]);
  }
#endif
  else if (options[AZ_precond] == AZ_lu) {
    AZ_lu_y12m(Nexp, Mexp, x, length, buffer, val, indx, bindx, rpntr, cpntr,
               bpntr, options, data_org, params[AZ_drop]);
  }
  else if (options[AZ_precond] == AZ_bilu) {
    AZ_block_ilu(val, indx, bindx, cpntr, bpntr, new_blks, Nexp, length,
                 buffer, x, options, data_org);
  }

  /* Add the values that are redundant. That is, add the external values 
   * to the border values that correspond to them. This will make the    
   * operator symmetric if the incomplete factorization used above was   
   * symmetric.                                                          */

  if (sym_flag == 1) {

     /* first send the external points to the neighbors */
 
     type            = AZ_sys_msg_type;
     AZ_sys_msg_type = (AZ_sys_msg_type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + 
                       AZ_MSG_TYPE;

     /* figure out what is the longest message to be */
     /* received and allocate space for it.          */
 
     total = 0;
     for ( i = 0 ; i < data_org[AZ_N_neigh] ; i++ )
        total += data_org[AZ_send_length+i];
     little = (double *) calloc(total, sizeof(double) ) ;

     /* post receives */

     count = 0;
     for ( i = 0 ; i < data_org[AZ_N_neigh] ; i++ ) {
        from = data_org[AZ_neighbors+i];
        (void) md_wrap_iread((void *) &(little[count]),
                  sizeof(double)*data_org[AZ_send_length+i],
                  &from, &type, request+i);
        count += data_org[AZ_send_length+i];
     }

     /* send messages */

     count = data_org[AZ_N_internal] + data_org[AZ_N_border];
     for ( i = 0 ; i < data_org[AZ_N_neigh] ; i++ ) {
        (void) md_wrap_write((void *) &(x[count]), data_org[AZ_rec_length+i]*
                     sizeof(double), data_org[AZ_neighbors+i], type, &st);
        count += data_org[AZ_rec_length+i];
     }
 
     /* receive messages and add recvd values to the send list */
 
     count = 0;
     for ( i = 0 ; i < data_org[AZ_N_neigh] ; i++ ) {
        from = data_org[AZ_neighbors+i];
        (void) md_wrap_wait((void *) &(little[count]),
                  sizeof(double)*data_org[AZ_send_length+i],
                  &from, &type, &st,request+i);
        count += data_org[AZ_send_length+i];
     }
     for ( j = 0 ; j < total; j++ ) 
        x[ data_org[AZ_send_list+j] ] += little[j];

     free(little);
     options[AZ_overlap] = AZ_sym_full;

  }

} /* AZ_domain_decomp */
