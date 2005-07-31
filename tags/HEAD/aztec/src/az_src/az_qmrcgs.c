/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile: az_qmrcgs.c,v $
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
static char rcsid[] = "$Id: az_qmrcgs.c,v 1.1.1.1 2001-01-30 20:59:12 paklein Exp $";
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

void AZ_pqmrs(double val[], int indx[], int bindx[], int rpntr[], int cpntr[],
              int bpntr[], double b[], double x[], double weight[],
              int options[], double params[], int data_org[], int proc_config[],
              double status[])

/*******************************************************************************

  Freund's transpose free QMR routine to solve the nonsymmetric matrix problem
  Ax = b. NOTE: this routine differs from Freund's paper in that we compute
  ubar (= M^-1 u ) and qbar (= M^-1 q) instead of u and q defined in Freund's
  paper.

  IMPORTANT NOTE: While an estimate of the 2-norm of the qmr residual is
  available (see comment below), the actual qmr residual is not normally
  computed as part of the qmr algorithm. Thus, if the user uses a convergence
  condition (see AZ_compute_global_scalars()) that is based on the 2-norm of the
  residual there is no need to compute the residual (i.e. r_avail = AZ_FALSE).
  However, if another norm of r is requested, AZ_compute_global_scalars() will
  set r_avail = AZ_TRUE and the algorithm will compute the residual.

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

  register int    i;
  int             N, NN, converged, one = 1, iter= 1,r_avail = AZ_FALSE, j;
  int             precond_flag, print_freq, proc;
  int             brkdown_will_occur = AZ_FALSE;
  double          alpha, beta = 0.0, true_scaled_r;
  double          *ubar, *v, *r_cgs, *rtilda, *Aubar, *qbar, *Aqbar, *d, *Ad;
  double          rhonm1, rhon, est_residual, actual_residual = -1.0;
  double          scaled_r_norm, sigma, epsilon, brkdown_tol = DBL_EPSILON;
  double          omega, c, norm_r_n_cgs, norm_r_nm1_cgs;
  double          tau_m, nu_m, eta_m, init_time;
  double          tau_mm1, nu_mm1 = 0.0, eta_mm1 = 0.0;
  register double dtemp;
  double          W_norm = 0.0;
  int             offset = 0;

  /**************************** execution begins ******************************/

  /* pull needed values out of parameter arrays */

  N            = data_org[AZ_N_internal] + data_org[AZ_N_border];
  precond_flag = options[AZ_precond];
  epsilon      = params[AZ_tol];
  proc         = proc_config[AZ_node];
  print_freq   = options[AZ_print_freq];

  /* allocate memory for required vectors */

  NN     = (N + data_org[AZ_N_external])*sizeof(double) + 1;
  /* +1 :make sure everyone allocates something */

  v      = (double *) AZ_manage_memory(NN, AZ_ALLOC, AZ_SYS, "v     in qmr",
                                       &j);
  ubar   = (double *) AZ_manage_memory(NN, AZ_ALLOC, AZ_SYS, "ubar  in qmr",
                                       &j);
  r_cgs  = (double *) AZ_manage_memory(NN, AZ_ALLOC, AZ_SYS, "r_cgs in qmr",
                                       &j);
  d      = (double *) AZ_manage_memory(NN, AZ_ALLOC, AZ_SYS, "d     in qmr",
                                       &j);
  qbar   = (double *) AZ_manage_memory(NN, AZ_ALLOC, AZ_SYS, "qbar  in qmr",
                                       &j);
  Aubar  = (double *) AZ_manage_memory(NN, AZ_ALLOC, AZ_SYS, "Aubar in qmr",
                                       &j);
  Aqbar  = (double *) AZ_manage_memory(NN, AZ_ALLOC, AZ_SYS, "Aqbar in qmr",
                                       &j);
  rtilda = (double *) AZ_manage_memory(NN, AZ_ALLOC, AZ_SYS, "rtilda in qmr",
                                       &j);

  AZ_compute_residual(val, indx, bindx, rpntr, cpntr, bpntr, b, x, r_cgs,
                      data_org);

  /* d, qbar, Aqbar, v = 0 */

  for (i = 0; i < N; i++)
    d[i] = qbar[i] = Aqbar[i] = v[i] = 0.0;

  /* set rtilda */

  if (options[AZ_aux_vec] == AZ_resid) dcopy_(&N, r_cgs, &one, rtilda, &one);
  else AZ_random_vector(rtilda, data_org, proc_config);

  /*
   * Compute a few global scalars:
   *     1) ||r_cgs||              corresponding to options[AZ_conv]
   *     2) scaled ||r_cgs||       corresponding to options[AZ_conv]
   *     3) rhon = <rtilda, r_cgs>
   * Note: step 1) is performed if r_avail = AZ_TRUE on entry or
   *       AZ_FIRST_TIME is passed in. Otherwise, ||r_cgs|| is taken as
   *       est_residual.
   */

  AZ_compute_global_scalars(val, indx, bindx, rpntr, cpntr, bpntr, x, b, r_cgs,
                            weight, &est_residual, &scaled_r_norm, options,
                            data_org, proc_config, &r_avail, r_cgs, rtilda,
                            &rhon, AZ_FIRST_TIME);
  true_scaled_r = scaled_r_norm;

  if ((options[AZ_output] != AZ_none) && (options[AZ_output] != AZ_last) &&
      (options[AZ_output] != AZ_warnings) && (proc == 0))
    (void) fprintf(stdout, "\t\titer:    0\t\tresidual = %e\n", scaled_r_norm);

  norm_r_nm1_cgs = est_residual;
  tau_mm1        = norm_r_nm1_cgs;
  rhonm1         = rhon;

  /* Set up aux-vector if we need to compute the qmr residual */

  if (r_avail) {
    Ad = (double *) AZ_manage_memory(NN, AZ_ALLOC, AZ_SYS, "Ad    in qmr", &j);
    for (i = 0; i < N; i++) Ad[i] = 0.0;
  }

  converged = scaled_r_norm < epsilon;

  for (iter = 1; iter <= options[AZ_max_iter] && !converged; iter++) {

    if (fabs(rhon) < brkdown_tol) { /* possible breakdown problem */

      if (AZ_breakdown_f(N, r_cgs, rtilda, rhon, proc_config))
        brkdown_will_occur = AZ_TRUE;
      else
        brkdown_tol = 0.1 * fabs(rhon);
    }

    /* ubar  = M^-1 r_cgs + beta*qbar               */
    /* Aubar = A ubar                               */
    /* v     = A ubar + beta ( A qbar + beta pnm1 ) */
    /*       = Aubar  + beta ( Aqbar +  beta v)     */

    dcopy_(&N, r_cgs, &one, ubar, &one);

    if (iter==1) init_time = AZ_second();

    if (precond_flag)
      AZ_precondition(val, indx, bindx, rpntr, cpntr, bpntr, ubar, options,
                      data_org, proc_config, params);

    if (iter==1) status[AZ_first_precond] = AZ_second() - init_time;

    for (i = 0; i < N; i++) ubar[i] = ubar[i] + beta * qbar[i];

    AZ_matvec_mult(val, indx, bindx, rpntr,cpntr, bpntr, ubar, Aubar, 1,
                   data_org);
    daxpy_(&N, &beta, v, &one, Aqbar, &one);
    for (i = 0; i < N; i++) v[i] = Aubar[i] + beta * Aqbar[i];

    sigma = AZ_gdot(N, rtilda, v, proc_config);

    if (fabs(sigma) < brkdown_tol) { /* possible problem */
      if (AZ_breakdown_f(N, rtilda, v, sigma, proc_config)) {

        /* break down */

        AZ_scale_true_residual(val, indx, bindx, rpntr, cpntr, bpntr, x, b, v,
                               weight, &actual_residual, &true_scaled_r,
                               options, data_org, proc_config);

        AZ_terminate_status_print(AZ_breakdown, iter, status, est_residual,
                                  params, true_scaled_r, actual_residual,
                                  options, proc_config);
        return;
      }
      else brkdown_tol = 0.1 * fabs(sigma);
    }

    alpha = rhon / sigma;

    /* qbar  = ubar - alpha* M^-1 v            */
    /* Aqbar = A qbar                          */
    /* r_cgs = r_cgs - alpha (A ubar + A qbar) */
    /*       = r_cgs - alpha (Aubar + Aqbar)   */

    dcopy_(&N, v, &one, qbar, &one);

    if (precond_flag)
      AZ_precondition(val, indx, bindx, rpntr, cpntr, bpntr, qbar, options,
                      data_org, proc_config, params);

    for (i = 0; i < N; i++) qbar[i] = ubar[i] - alpha * qbar[i];
    AZ_matvec_mult(val, indx, bindx, rpntr, cpntr, bpntr, qbar, Aqbar, 1,
                   data_org);
    for (i = 0; i < N; i++) r_cgs[i] = r_cgs[i] - alpha*(Aubar[i] + Aqbar[i]);

    /* QMRS scaling and iterates weights 5.11 */

    norm_r_n_cgs = sqrt(AZ_gdot(N, r_cgs, r_cgs, proc_config));

    /* m is odd in Freund's paper */

    omega = sqrt(norm_r_nm1_cgs * norm_r_n_cgs);
    nu_m  = omega / tau_mm1;
    c     = 1.0 / sqrt(1.0 + nu_m * nu_m);
    tau_m = tau_mm1 * nu_m * c;
    eta_m = c * c * alpha;

    if (brkdown_will_occur) {
      AZ_scale_true_residual(val, indx, bindx, rpntr, cpntr, bpntr, x, b, v,
                             weight, &actual_residual, &true_scaled_r, options,
                             data_org, proc_config);

      AZ_terminate_status_print(AZ_breakdown, iter, status, est_residual,
                                params, true_scaled_r, actual_residual, options,
                                proc_config);
      return;
    }

    dtemp = nu_mm1 *nu_mm1 * eta_mm1 / alpha;
    for (i = 0; i < N; i++) d[i] = ubar[i] + dtemp * d[i];
    daxpy_(&N, &eta_m, d, &one, x, &one); /* x = x - eta_m d  */

    if (r_avail) {
      for (i = 0; i < N; i++) Ad[i] = Aubar[i] + dtemp * Ad[i];
    }

    /* save some values */

    eta_mm1 = eta_m;  tau_mm1        = tau_m;
    nu_mm1  = nu_m;   norm_r_nm1_cgs = norm_r_n_cgs;

    /* m is even in Freund's paper */

    omega = norm_r_n_cgs;

    if (tau_mm1 == 0.0) nu_m = 0.0;
    else                nu_m  = omega / tau_mm1;

    c     = 1.0 / sqrt(1.0 + nu_m * nu_m);
    tau_m = tau_mm1 * nu_m * c;
    eta_m = c * c * alpha;

    dtemp = nu_mm1 * nu_mm1 * eta_mm1 / alpha;
    for (i = 0; i < N; i++) d[i] = qbar[i] + dtemp * d[i];
    daxpy_(&N, &eta_m, d, &one, x, &one); /* x = x - eta_m d  */

    if (r_avail) {
      for (i = 0; i < N; i++) Ad[i] = Aqbar[i] + dtemp * Ad[i];
    }

    /* save some values */

    eta_mm1 = eta_m;  tau_mm1        = tau_m;
    nu_mm1  = nu_m;   norm_r_nm1_cgs = norm_r_n_cgs;
    rhonm1  = rhon;

    if (r_avail) {
      for (i = 0; i < N; i++) Aubar[i] = r_cgs[i] - (eta_m - alpha) * Ad[i];

      /* Note: Aubar temporarily holds qmr residual */

    }
    else {

      /*
       * We want to estimate the 2-norm of the qmr residual. Freund gives the
       * bound ||r|| <= tau_m * sqrt(2*iter+1).  We use this bound until we get
       * close to the solution. At that point we compute the real residual norm
       * and use this to estimate the norm of ||W|| in Freund's paper.
       */

      dtemp = sqrt((double) (2 * iter + 1));

      if ((scaled_r_norm < epsilon * dtemp) && !offset) {
        AZ_scale_true_residual(val,indx, bindx, rpntr, cpntr, bpntr, x, b,
                               Aubar, weight, &actual_residual, &true_scaled_r,
                               options, data_org, proc_config);

        if (tau_m != 0.0) W_norm = actual_residual / tau_m;
        if (W_norm < 1.0) W_norm = 1.0;

        offset       = 2 * iter + 1;
        est_residual = actual_residual;
      }
      else est_residual = sqrt((double)(2 * iter + 1 - offset) +
                               W_norm * W_norm) * tau_m;
    }

    /*
     * Compute a few global scalars:
     *     1) ||r||                corresponding to options[AZ_conv]
     *     2) scaled ||r||         corresponding to options[AZ_conv]
     *     3) rhon = <rtilda, r_cgs>
     * Note: step 1) is performed if r_avail = AZ_TRUE or AZ_FIRST_TIME
     *       is passed in. Otherwise, ||r|| is taken as est_residual.
     */

    AZ_compute_global_scalars(val, indx, bindx, rpntr, cpntr, bpntr, x, b,
                              Aubar, weight, &est_residual, &scaled_r_norm,
                              options, data_org, proc_config, &r_avail, rtilda,
                              r_cgs, &rhon, AZ_NOT_FIRST);

    if ( (iter%print_freq == 0) && proc == 0 )
      (void) fprintf(stdout, "\t\titer: %4d\t\tresidual = %e\n", iter,
                     scaled_r_norm);

    /* convergence tests */

    if (scaled_r_norm < epsilon) {
      AZ_scale_true_residual(val, indx, bindx, rpntr, cpntr, bpntr, x, b, Aubar,
                             weight, &actual_residual, &true_scaled_r, options,
                             data_org, proc_config);

      converged = true_scaled_r < params[AZ_tol];

      /*
       * Note: epsilon and params[AZ_tol] may not be equal due to a previous
       * call to AZ_get_new_eps().
       */

      if (!converged  &&
          (AZ_get_new_eps(&epsilon, scaled_r_norm, true_scaled_r,
                          proc_config) == AZ_QUIT)) {

        /*
         * Computed residual has converged, actual residual has not converged,
         * AZ_get_new_eps() has decided that it is time to quit.
         */

        AZ_terminate_status_print(AZ_loss, iter, status, est_residual, params,
                                  true_scaled_r, actual_residual, options,
                                  proc_config);
        return;
      }
    }

    beta = rhon / rhonm1;
  }

  iter--;

  if ( (iter%print_freq != 0) && (proc == 0) && (options[AZ_output] != AZ_none)
       && (options[AZ_output] != AZ_warnings))
    (void) fprintf(stdout, "\t\titer: %4d\t\tresidual = %e\n", iter,
                   scaled_r_norm);

  /* check if we exceeded maximum number of iterations */

  if (converged) {
    i = AZ_normal;
    scaled_r_norm = true_scaled_r;
  }
  else
    i = AZ_maxits;

  AZ_terminate_status_print(i, iter, status, est_residual, params,
                            scaled_r_norm, actual_residual, options,
                            proc_config);

} /* pqmrs */
