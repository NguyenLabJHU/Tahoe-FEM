/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile: az_converge.c,v $
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
static char rcsid[] = "$Id: az_converge.c,v 1.1.1.1 2001-01-30 20:59:14 paklein Exp $";
#endif


/*******************************************************************************
 * Copyright 1995, Sandia Corporation.  The United States Government retains a *
 * nonexclusive license in this software as prescribed in AL 88-1 and AL 91-7. *
 * Export of this program may require a license from the United States         *
 * Government.                                                                 *
 ******************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifndef __MWERKS__
#include <malloc.h>
#endif /* __MWERKS __ */
#include <float.h>
#include "az_aztec.h"

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_compute_global_scalars(double val[], int indx[], int bindx[],
                               int rpntr[], int cpntr[], int bpntr[],
                               double x[], double b[], double r[], double w[],
                               double *r_norm, double *scaled_r_norm,
                               int options[], int data_org[], int proc_config[],
                               int *r_avail, double v1[], double v2[],
                               double *value, int first_time)

/*******************************************************************************

  Routine to check against 'eps' for convergence. The type of norm use is
  determined by the variable 'conv_flag' as follows:

                0: ||r||2 / ||r0||2  < eps
                1: ||r||2 / ||b||2   < eps
                2: ||r||2 / ||A||inf < eps
                3: ||r||inf / (||A||inf * ||x||1 + ||b||inf)
                4: ||r/w||2

                where ||*||2 is the Euclidean norm, ||*||inf is the
                infinity norm and ||*|| is the sum norm.

  Author:          Scott A. Hutchinson, SNL, 1421
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

  x:               The current solution vector.

  b:               Right hand side of linear system.

  r:               The current residual vector.

  w:               Weighting vector for convergence norm #4.

  r_norm:          Norm of residual.

  scaled_r_norm:   Norm of residual scaled by norm of the rhs.

  options:         Determines specific solution method and other parameters.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see file params.txt).

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

  r_avail:         In general, this variable indicates whether or not the
                   residual is available or needs to be made available.
                   In particular,
                      first_time == TRUE     : real residual is available. The
                                               norm of this residual will be
                                               computed and stored in r_norm.
                      first_time == FALSE &&
                      r_avail    == TRUE     : real residual is available. The
                                               norm of this residual will be
                                               computed and stored in r_norm.
                      first_time == FALSE &&
                      r_avail    == FALSE    : real residual is not available.
                                               The norm of the residual is not
                                               computed. Instead, it is assumed
                                               that r_norm is an estimate of the
                                               residual norm.
                   All of this is done for gmres() and tfqmr() where we often
                   have estimates of the residual 2-norm without actually having
                   computed the residual.

                   IMPORTANT: if a convergence test requires a residual norm
                   other than the 2-norm, it is important that
                   AZ_compute_global_scalars() sets r_avail to TRUE. This tells
                   the iterative method (in particular gmres and tfqmr) that the
                   real residual must be computed (at additional cost)
                   and passed in to AZ_compute_global_scalars().

  v1,v2,value:     If v1 != NULL, *value = <v1,v2> where <.,.> is the standard
                   inner product, v1 and v2 are double precision vectors of
                   length data_org[AZ_num_int_unk] + data_org[AZ_num_bord_unk].

                   This is used so that 1 inner product can be grouped together
                   with the inner products required for convergence (to save
                   some messages).

  first_time:      Flag set AZ_TRUE if this is the first time this routine has
                   been called.  Set AZ_FALSE otherwise. See comments for
                   r_avail above for more information.

*******************************************************************************/

{

  /* local variables */

  register int  i;
  static double r0, amax, b_norm, *temp;
  static int    total_N;
  int           N, j;
  double        dots[5], tmp[5];
  int           count = 0, one = 1;

  /**************************** execution begins ******************************/

  N = data_org[AZ_N_internal] + data_org[AZ_N_border];

  if (v1 != NULL) dots[count++] = ddot_(&N, v1, &one, v2, &one);

  /* initialize */

  switch (options[AZ_conv]) {

  case AZ_r0:
    if ((*r_avail) || first_time) {
      dots[count] = ddot_(&N, r, &one, r, &one);

      AZ_gdot_vec(count + 1, dots, tmp, proc_config);
      *r_norm = sqrt(dots[count]);
    }
    else if (v1 != NULL) AZ_gdot_vec(1, dots, tmp, proc_config);

    if (v1 != NULL) *value = dots[0];

    if (first_time) {
      r0 = *r_norm;
      if (fabs(r0) < DBL_MIN) r0 = 1.0;
    }

    *scaled_r_norm = *r_norm / r0;
    break;

  case AZ_rhs:
    if (first_time) {
      dots[count++] = ddot_(&N, b, &one, b, &one);
      dots[count  ] = ddot_(&N, r, &one, r, &one);

      AZ_gdot_vec(count + 1, dots, tmp, proc_config);
      *r_norm = sqrt(dots[count]); count--;
      b_norm  = sqrt(dots[count]);
    }

    else if (*r_avail) {
      dots[count] = ddot_(&N, r, &one, r, &one);

      AZ_gdot_vec(count + 1, dots, tmp, proc_config);
      *r_norm = sqrt(dots[count]);
    }

    else if (v1 != NULL) AZ_gdot_vec(1, dots, tmp, proc_config);

    if (v1 != NULL) *value = dots[0];

    *scaled_r_norm = *r_norm / b_norm;
    break;

  case AZ_Anorm:
    if (first_time) {
      dots[count] = ddot_(&N, r, &one, r, &one);

      AZ_gdot_vec(count+1, dots, tmp, proc_config);
      *r_norm = sqrt(dots[count]);
      amax    = AZ_gmax_matrix_norm(val, indx, bindx, rpntr, cpntr, bpntr,
                                    proc_config, data_org);

      if (fabs(amax) < DBL_MIN) {
        AZ_p_error("Error: ||A||_infinity = 0\n\tSomething wrong with A\n",
                   proc_config[AZ_node]);
        exit(-1);
      }
    }

    else if (*r_avail) {
      dots[count] = ddot_(&N, r, &one, r, &one);

      AZ_gdot_vec(count + 1, dots, tmp, proc_config);
      *r_norm = sqrt(dots[count]);
    }

    else if (v1 != NULL) AZ_gdot_vec(1, dots, tmp, proc_config);

    if (v1 != NULL) *value = dots[0];

    *scaled_r_norm = *r_norm / amax;
    break;

  case AZ_sol:
    *r_avail = AZ_TRUE;

    if (v1 != NULL) {
      AZ_gdot_vec(1, dots, tmp, proc_config);
      *value = dots[0];
    }

    if (first_time) {
      amax = AZ_gmax_matrix_norm(val, indx, bindx, rpntr, cpntr, bpntr,
                                 proc_config, data_org);

      if (fabs(amax) < DBL_MIN) {
        AZ_p_error("Error: ||A||_infinity = 0\n\tSomething wrong with A\n",
                   proc_config[AZ_node]);
        exit(-1);
      }

      b_norm = AZ_gvector_norm(N, -1, b, proc_config);
    }

    *r_norm = AZ_gvector_norm(N, -1, r, proc_config);
    *scaled_r_norm = *r_norm / (amax *
                                AZ_gvector_norm(N, 1, x, proc_config) + b_norm);
    break;

  case AZ_weighted:
    *r_avail = AZ_TRUE;

    if (first_time) {
      temp = AZ_manage_memory((N + data_org[AZ_N_external]) * sizeof(double),
                              AZ_ALLOC, AZ_SYS,
                              "temp in AZ_compute_global_scalars", &j);
      total_N = AZ_gsum_int(N, proc_config);
    }

    for (i = 0; i < N; i++) temp[i] = r[i] / w[i];

    dots[count] = ddot_(&N, temp, &one, temp, &one);

    AZ_gdot_vec(count + 1, dots, tmp, proc_config);

    if (v1 != NULL) *value = dots[0];

    *scaled_r_norm = sqrt(dots[count] / (double) total_N);
    *r_norm        = *scaled_r_norm;
    break;

  default:
    if (proc_config[AZ_node] == 0) {
      (void) fprintf(stderr, "Error: Improper value, options[AZ_conv] = %d\n",
                     options[AZ_conv]);
    }
    exit(-1);
  }

} /* AZ_compute_global_scalars */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_scale_true_residual(double val[], int indx[], int bindx[], int rpntr[],
                            int cpntr[], int bpntr[], double x[], double b[],
                            double v[], double w[], double *actual_residual,
                            double *scaled_r_norm, int options[],
                            int data_org[], int proc_config[])

/*******************************************************************************

  Routine to check the true residual and decide whether or not actual
  convergence has been achieved based upon some initial convergence indicators.

  Author:          Scott A. Hutchinson, SNL, 1421
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

  x:               The current solution vector.

  b:               Right hand side of linear system.

  v:

  w:               Weighting vector for convergence norm #4.

  actual_residual: Norm of residual.

  scaled_r_norm:   Norm of residual scaled by norm of the rhs.

  options:         Determines specific solution method and other parameters.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see file params.txt).

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  /* local variables */

  int r_avail = AZ_TRUE;

  /**************************** execution begins ******************************/

  /* calculate true residual */

  AZ_compute_residual(val, indx, bindx, rpntr, cpntr, bpntr, b, x, v, data_org);

  /* compute scaled residual */

  AZ_compute_global_scalars(val, indx, bindx, rpntr, cpntr, bpntr, x, b, v, w,
                            actual_residual, scaled_r_norm, options, data_org,
                            proc_config,&r_avail, (double *) 0, (double *) 0,
                            (double *) 0, AZ_NOT_FIRST);

} /* AZ_scale_true_residual */
