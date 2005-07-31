/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile: az_cg.c,v $
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
static char rcsid[] = "$Id: az_cg.c,v 1.1.1.1 2001-01-30 20:59:14 paklein Exp $";
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
#include <float.h>
#include "az_aztec.h"

void AZ_pcg_f(double val[], int indx[], int bindx[], int rpntr[], int cpntr[],
              int bpntr[], double b[], double x[], double weight[],
              int options[], double params[], int data_org[], int proc_config[],
              double status[])

/*******************************************************************************

  Conjugate Gradient algorthm to solve the symmetric matrix problem Ax = b.

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

  weight:          Vector of weights for convergence norm #4.

  options:         Determines specific solution method and other parameters.

  params:          Drop tolerance and convergence tolerance info.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see file params.txt).

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

  status:          On output, indicates termination status:
                    0:  terminated normally.
                   -1:  maximum number of iterations taken without achieving
                        convergence.
                   -2:  Breakdown. The algorithm can not proceed due to
                        numerical difficulties (usually a divide by zero).
                   -3:  Internal residual differs from the computed residual due
                        to a significant loss of precision.

*******************************************************************************/

{

  /* local variables */

  register int i;
  int          N, NN, converged, one = 1, iter = 1, r_avail = AZ_TRUE, j;
  int          precond_flag, print_freq, proc, brkdown_will_occur = AZ_FALSE;
  double       alpha, beta = 0.0, nalpha, true_scaled_r;
  double      *r, *z, *p, *ap, actual_residual = -1.0;
  double       r_z_dot, r_z_dot_old, p_ap_dot, rec_residual;
  double       scaled_r_norm, epsilon, brkdown_tol = DBL_EPSILON;

  /**************************** execution begins ******************************/

  /* pull needed values out of parameter arrays */

  N            = data_org[AZ_N_internal] + data_org[AZ_N_border];
  precond_flag = options[AZ_precond];
  epsilon      = params[AZ_tol];
  proc         = proc_config[AZ_node];
  print_freq   = options[AZ_print_freq];

  /* allocate space for necessary vectors */

  NN = (N + data_org[AZ_N_external])*sizeof(double) + 1;
  /* +1: make sure everybody allocates something */

  z  = (double *) AZ_manage_memory(NN, AZ_ALLOC, AZ_SYS, "z in cg", &j);
  r  = (double *) AZ_manage_memory(NN, AZ_ALLOC, AZ_SYS, "r in cg", &j);
  p  = (double *) AZ_manage_memory(NN, AZ_ALLOC, AZ_SYS, "p in cg", &j);
  ap = (double *) AZ_manage_memory(NN, AZ_ALLOC, AZ_SYS, "ap in cg", &j);

  AZ_compute_residual(val, indx, bindx, rpntr, cpntr, bpntr, b, x, r, data_org);

  /*  z = M r */
  /*  p = 0   */

  dcopy_(&N, r, &one, z, &one);
  status[AZ_first_precond] = AZ_second();
  if (precond_flag)
    AZ_precondition(val, indx, bindx, rpntr, cpntr, bpntr, z, options, data_org,
                    proc_config, params);
  status[AZ_first_precond] = AZ_second() - status[AZ_first_precond];

  for (i = 0; i < N; i++ ) p[i] = 0.0;

  /* compute a few global scalars:                                 */
  /*     1) ||r||                corresponding to options[AZ_conv] */
  /*     2) scaled ||r||         corresponding to options[AZ_conv] */
  /*     3) r_z_dot = <z, r>                                       */

  AZ_compute_global_scalars(val, indx, bindx, rpntr, cpntr, bpntr, x, b, r,
                            weight, &rec_residual, &scaled_r_norm, options,
                            data_org, proc_config, &r_avail,r, z, &r_z_dot,
                            AZ_FIRST_TIME);
  true_scaled_r = scaled_r_norm;

  if ((options[AZ_output] != AZ_none) &&
      (options[AZ_output] != AZ_last) &&
      (options[AZ_output] != AZ_warnings) && (proc == 0))
    (void) fprintf(stdout, "\t\titer:    0\t\tresidual = %e\n", scaled_r_norm);

  converged = scaled_r_norm < epsilon;

  for (iter = 1; iter <= options[AZ_max_iter] && !converged; iter++ ) {

    /* p  = z + beta * p */
    /* ap = A p          */

    for (i = 0; i < N; i++) p[i] = z[i] + beta * p[i];
    AZ_matvec_mult(val, indx, bindx, rpntr, cpntr, bpntr, p, ap, 1, data_org);

    p_ap_dot = AZ_gdot(N, p, ap, proc_config);
    if (fabs(p_ap_dot) < brkdown_tol) {

      /* possible problem */

      if (AZ_breakdown_f(N, p, ap, p_ap_dot, proc_config)) {

        /* something wrong */

        AZ_scale_true_residual(val, indx, bindx, rpntr, cpntr, bpntr, x, b, ap,
                               weight, &actual_residual, &true_scaled_r,
                               options, data_org, proc_config);
        AZ_terminate_status_print(AZ_breakdown, iter, status, rec_residual,
                                  params, true_scaled_r, actual_residual,
                                  options, proc_config);
        return;
      }
      else brkdown_tol = 0.1 * fabs(p_ap_dot);
    }

    alpha  = r_z_dot / p_ap_dot;
    nalpha = -alpha;

    /* x = x + alpha*p  */
    /* r = r - alpha*Ap */
    /* z = M^-1 r       */

    daxpy_(&N, &alpha,  p,  &one, x, &one);
    daxpy_(&N, &nalpha, ap, &one, r, &one);
    dcopy_(&N, r, &one, z, &one);

    if (precond_flag) AZ_precondition(val, indx, bindx, rpntr, cpntr, bpntr, z,
                                      options, data_org, proc_config, params);

    r_z_dot_old = r_z_dot;

    /* compute a few global scalars:                                 */
    /*     1) ||r||                corresponding to options[AZ_conv] */
    /*     2) scaled ||r||         corresponding to options[AZ_conv] */
    /*     3) r_z_dot = <z, r>                                       */

    AZ_compute_global_scalars(val, indx, bindx, rpntr, cpntr, bpntr, x, b, r,
                              weight, &rec_residual, &scaled_r_norm, options,
                              data_org, proc_config, &r_avail, r, z, &r_z_dot,
                              AZ_NOT_FIRST);

    if (brkdown_will_occur) {
      AZ_scale_true_residual(val, indx, bindx, rpntr, cpntr, bpntr, x, b, ap,
                             weight, &actual_residual, &true_scaled_r, options,
                             data_org, proc_config);
      AZ_terminate_status_print(AZ_breakdown, iter, status, rec_residual,
                                params, true_scaled_r, actual_residual, options,
                                proc_config);
      return;
    }

    beta = r_z_dot / r_z_dot_old;

    if (fabs(r_z_dot) < brkdown_tol) {

      /* possible problem */

      if (AZ_breakdown_f(N, r, z, r_z_dot, proc_config))
        brkdown_will_occur = AZ_TRUE;
      else
        brkdown_tol = 0.1 * fabs(r_z_dot);
    }

    if ( (iter%print_freq == 0) && proc == 0 )
      (void) fprintf(stdout, "\t\titer: %4d\t\tresidual = %e\n", iter,
                     scaled_r_norm);

    /* convergence tests */

    if (scaled_r_norm < epsilon) {
      AZ_scale_true_residual(val, indx, bindx, rpntr, cpntr, bpntr, x, b, ap,
                             weight, &actual_residual, &true_scaled_r, options,
                             data_org, proc_config);

      converged = true_scaled_r < params[AZ_tol];

      /*
       * Note: epsilon and params[AZ_tol] may not be equal due to a previous
       * call to AZ_get_new_eps().
       */

      if (!converged &&
          (AZ_get_new_eps(&epsilon, scaled_r_norm, true_scaled_r,
                          proc_config) == AZ_QUIT)) {

        /*
         * Computed residual has converged, actual residual has not converged,
         * AZ_get_new_eps() has decided that it is time to quit.
         */

        AZ_terminate_status_print(AZ_loss, iter, status, rec_residual, params,
                                  true_scaled_r, actual_residual, options,
                                  proc_config);
        return;
      }
    }
  }

  iter--;
  if ( (iter%print_freq != 0) && (proc == 0) && (options[AZ_output] != AZ_none)
       && (options[AZ_output] != AZ_warnings) )
    (void) fprintf(stdout, "\t\titer: %4d\t\tresidual = %e\n", iter,
                   scaled_r_norm);

  /* check if we exceeded maximum number of iterations */

  if (converged) {
    i = AZ_normal; scaled_r_norm = true_scaled_r; }
  else
    i = AZ_maxits;

  AZ_terminate_status_print(i, iter, status, rec_residual, params,
                            scaled_r_norm, actual_residual, options,
                            proc_config);

} /* AZ_pcg */
