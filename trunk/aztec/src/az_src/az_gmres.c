/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile: az_gmres.c,v $
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
static char rcsid[] = "$Id: az_gmres.c,v 1.1.1.1 2001-01-30 20:59:14 paklein Exp $";
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

void AZ_pgmres(double val[], int indx[], int bindx[], int rpntr[], int cpntr[],
            int bpntr[], double b[], double x[],double weight[], int options[],
            double params[], int data_org[], int proc_config[], double status[])

/*******************************************************************************

  This routine uses Saad's restarted Genralized Minimum Residual method to solve
  the nonsymmetric matrix problem Ax = b.

  IMPORTANT NOTE: While the 2-norm of the gmres residual is available, the
  actual residual is not normally computed as part of the gmres algorithm. Thus,
  if the user uses a convergence condition (see AZ_compute_global_scalars())
  that is based on the 2-norm of the residual there is no need to compute the
  residual (i.e. r_avail = AZ_FALSE). However, if another norm of r is
  requested, AZ_compute_global_scalars() sets r_avail = AZ_TRUE and the
  algorithm computes the residual.

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

  register int k;
  int          i, i1, k1;
  int          N, NN, converged, one = 1, iter = 1, r_avail = AZ_FALSE;
  int          precond_flag, print_freq, proc, kspace;
  double     **v, **hh, *c, *s, *rs, *dots, *tmp, *temp;
  double       *res, init_time;
  double       dble_tmp, r_2norm, first_2norm, epsilon;
  double       rec_residual, scaled_r_norm, true_scaled_r;
  double       actual_residual = -1.0;
  double       doubleone = 1.0, minusone = -1.0, *dummy = (double *) 0;
  char         str[30];

  /* condition number estimation variables */

  double *svbig, *svsml, *vectmp;
  double  big, cc, dble_tmp1, sestpr, small, ss;
  int     ijob;

  /**************************** execution begins ******************************/

  /* pull needed values out of parameter arrays */

  N            = data_org[AZ_N_internal] + data_org[AZ_N_border];
  precond_flag = options[AZ_precond];
  epsilon      = params[AZ_tol];
  proc         = proc_config[AZ_node];
  print_freq   = options[AZ_print_freq];
  kspace       = options[AZ_kspace];

  /* allocate memory for required vectors */

  NN    = (kspace + 1) * sizeof(double) + 1;
  /* +1: make sure everybody allocates something */

  dots  = AZ_manage_memory(NN, AZ_ALLOC, AZ_SYS, "dots in gmres", &i);
  tmp   = AZ_manage_memory(NN, AZ_ALLOC, AZ_SYS, "tmp  in gmres", &i);
  rs    = AZ_manage_memory(NN, AZ_ALLOC, AZ_SYS, "rs   in gmres", &i);
  v     = (double **) AZ_manage_memory(NN, AZ_ALLOC, AZ_SYS, "v   in gmres",
                                       &i);
  hh    = (double **) AZ_manage_memory(NN, AZ_ALLOC, AZ_SYS, "hh   in gmres",
                                       &i);

  NN    = kspace * sizeof(double) + 1;
  /* +1: make sure everybody allocates something */

  c     = (double *) AZ_manage_memory(NN, AZ_ALLOC, AZ_SYS, "c in gmres", &i);
  s     = (double *) AZ_manage_memory(NN, AZ_ALLOC, AZ_SYS, "s in gmres", &i);
  svbig = (double *) AZ_manage_memory(NN, AZ_ALLOC, AZ_SYS, "svbig in gmres",
                                      &i);
  svsml = (double *) AZ_manage_memory(NN, AZ_ALLOC, AZ_SYS, "svsml in gmres",
                                      &i);
  vectmp= (double *) AZ_manage_memory(NN, AZ_ALLOC, AZ_SYS, "vectmp in gmres",
                                      &i);

  for (k = 0; k < kspace + 1; k++) {
    sprintf(str, "hh[%d] in gmres", k);
    hh[k] = AZ_manage_memory(NN, AZ_ALLOC, AZ_SYS, str, &i);
  }

  NN = (N + data_org[AZ_N_external]) * sizeof(double) + 1;
  /* +1: make sure everybody allocates something */

  temp = AZ_manage_memory(NN, AZ_ALLOC, AZ_SYS, "temp in gmres", &i);

  for (k = 0; k < kspace+1; k++) {
    sprintf(str, "v[%d] in gmres", k);
    v[k] = AZ_manage_memory(NN, AZ_ALLOC, AZ_SYS, str, &i);
  }

  AZ_compute_residual(val, indx, bindx, rpntr, cpntr, bpntr, b, x, v[0],
                      data_org);

  /*
   * Compute a few global scalars:
   *     1) ||r||                corresponding to options[AZ_conv]
   *     2) scaled ||r||         corresponding to options[AZ_conv]
   *     3) r_2norm = <r,r>      corresponding to options[AZ_conv]
   */

  AZ_compute_global_scalars(val, indx, bindx, rpntr, cpntr, bpntr, x, b, v[0],
                            weight, &rec_residual, &scaled_r_norm, options,
                            data_org, proc_config, &r_avail, v[0], v[0],
                            &r_2norm, AZ_FIRST_TIME);
  true_scaled_r = scaled_r_norm;

  r_2norm   = sqrt(r_2norm);
  converged = scaled_r_norm < epsilon;

  if (r_avail)
    res = AZ_manage_memory(NN, AZ_ALLOC, AZ_SYS, "res in gmres", &i);
  else res = (double *) NULL;
  

  if ( (options[AZ_output] != AZ_none) && (options[AZ_output] != AZ_last) &&
       (options[AZ_output] != AZ_warnings) && (proc == 0) )
    (void) fprintf(stdout, "\t\titer:    0\t\tresidual = %e\n", scaled_r_norm);

  iter = 0;
  while (!converged && iter < options[AZ_max_iter]) {
    if (r_avail) dcopy_(&N, v[0], &one, res, &one);

    /* v1 = r0/beta */

    dble_tmp    = 1.0 / r_2norm;
    first_2norm = r_2norm;
    dscal_(&N, &dble_tmp, v[0], &one);

    rs[0] = r_2norm;  /* initialize 1st rhs term of H system */
    i     = 0;

    while (i < kspace && !converged && iter < options[AZ_max_iter]) {
      iter++;
      i1 = i + 1;

      /* v_i+1 = A M^-1 v_i */

      dcopy_(&N, v[i], &one, temp, &one);

      if (iter == 1) init_time = AZ_second();
 
      if (precond_flag) AZ_precondition(val, indx, bindx, rpntr, cpntr, bpntr,
                                        temp, options, data_org, proc_config,
                                        params);

      if (iter == 1) status[AZ_first_precond] = AZ_second() - init_time;

      AZ_matvec_mult(val, indx, bindx, rpntr, cpntr, bpntr, temp, v[i1], 1,
                     data_org);

      /* Gramm Schmidt orthogonalization */

      if (!options[AZ_orthog]) { /* classical */
        for (k = 0; k <= i; k++)
          dots[k] = ddot_(&N, v[k], &one, v[i1], &one);

        AZ_gdot_vec(i1, dots, tmp, proc_config);

        /* v_k+1 = inv(M)*A*v_k - sum_(1,k) h_(i,k)*v_i */

        for (k = 0; k <= i; k++) {
          hh[k][i] = dots[k];
          dble_tmp = -dots[k];
          daxpy_(&N, &dble_tmp, v[k], &one, v[i1], &one);
        }
      }
      else {                    /* modified */
        for (k = 0; k <= i; k++) {
          hh[k][i] = dble_tmp = AZ_gdot(N, v[k], v[i1], proc_config);
          dble_tmp = -dble_tmp;
          daxpy_(&N, &dble_tmp, v[k], &one, v[i1], &one);
        }
      }

      /* normalize vector */

      hh[i1][i] = dble_tmp = sqrt(AZ_gdot(N, v[i1], v[i1], proc_config));
      if (dble_tmp  > DBL_EPSILON*r_2norm)
        dble_tmp  = 1.0 / dble_tmp;
      else
        dble_tmp = 0.0;

      dscal_(&N, &dble_tmp, v[i1], &one);

      /* update factorization of hh by plane rotation */

      for (k = 1; k <= i; k++) {
        k1        = k - 1;
        dble_tmp  = hh[k1][i];
        hh[k1][i] =  c[k1]*dble_tmp + s[k1]*hh[k][i];
        hh[k][i]  = -s[k1]*dble_tmp + c[k1]*hh[k][i];
      }

      /* determine next plane rotation */

      dble_tmp = sqrt(hh[i][i] * hh[i][i] + hh[i1][i] * hh[i1][i]);

      /* Estimate condition number of the GMRES */
      /* least-squares problem using ICE.       */

      if (i == 0) {
        big = dble_tmp;
        small = big;
        svbig[0] = doubleone;
        svsml[0] = doubleone;
      }
      else {
        for (k = 0; k < i; k++) vectmp[k] = hh[k][i];
        vectmp[i] = dble_tmp;
        ijob = 1;
        dlaic1_(&ijob, &i, svbig, &big, vectmp, &vectmp[i], &sestpr, &ss, &cc);
        big = sestpr;
        dscal_(&i, &ss, svbig, &one);
        svbig[i] = cc;
        ijob = 2;
        dlaic1_(&ijob, &i, svsml, &small, vectmp, &vectmp[i], &sestpr, &ss,
                &cc);
        small = sestpr;
        dscal_(&i, &ss, svsml, &one);
        svsml[i] = cc;
      }

      if ((small == 0.0) || (dble_tmp  < DBL_EPSILON * r_2norm) ||
          (big/small > 1.0e+10) ) {

        /* take most recent solution and get out */

        for (k = 0; k < i1; k++) tmp[k] = rs[k];
        AZ_get_x_incr(val, indx, bindx, rpntr, cpntr, bpntr, options,
                      data_org, proc_config, params, i, hh, tmp, temp, v);
        daxpy_(&N, &doubleone, temp, &one, x, &one);

        AZ_scale_true_residual(val, indx, bindx, rpntr, cpntr, bpntr, x, b,
                               v[kspace], weight, &actual_residual,
                               &true_scaled_r, options, data_org, proc_config);

        if (dble_tmp  < DBL_EPSILON * r_2norm) i = AZ_breakdown;
        else                                   i = AZ_ill_cond;

        AZ_terminate_status_print(i, iter, status, rec_residual, params,
                                  true_scaled_r, actual_residual, options,
                                  proc_config);
        return;
      }

      dble_tmp = 1.0 / dble_tmp;
      c[i]     =  hh[i][i]  * dble_tmp;
      s[i]     =  hh[i1][i] * dble_tmp;
      rs[i1]   = -s[i] * rs[i];
      rs[i]   *=  c[i];

      if (r_avail) {
        dble_tmp  = c[i] * rs[i1];
        dble_tmp1 = s[i] * s[i];

        for (k = 0; k < N; k++)
          res[k] = dble_tmp*v[i1][k] + dble_tmp1*res[k];
      }

      /* determine residual norm & test convergence */

      hh[i][i]     = c[i] * hh[i][i] + s[i] * hh[i1][i];
      r_2norm      = fabs(rs[i1]);
      rec_residual = r_2norm;

      /*
       * Compute a few global scalars:
       *     1) ||r||                corresponding to options[AZ_conv]
       *     2) scaled ||r||         corresponding to options[AZ_conv]
       * NOTE: if r_avail = AZ_TRUE or AZ_FIRST is passed in, we perform
       * step 1), otherwise ||r|| is taken as rec_residual.
       */

      AZ_compute_global_scalars(val, indx, bindx, rpntr, cpntr, bpntr, x, b,
                                res, weight, &rec_residual, &scaled_r_norm,
                                options, data_org, proc_config, &r_avail, dummy,
                                dummy, dummy, AZ_NOT_FIRST);

      converged = scaled_r_norm < epsilon;

      if ( (iter%print_freq == 0) && proc == 0)
        (void) fprintf(stdout, "\t\titer: %4d\t\tresidual = %e\n", iter,
                       scaled_r_norm);

      i++;      /* subspace dim. counter dim(K) = i - 1 */

      if ( (i == kspace) || converged || iter == options[AZ_max_iter]) {

        /* update x and set temp to delta x */

        for (k = 0; k <= i1; k++) tmp[k] = rs[k];

        AZ_get_x_incr(val, indx, bindx, rpntr, cpntr, bpntr, options, data_org,
                      proc_config, params, i, hh, tmp, temp, v);
        daxpy_(&N, &doubleone, temp, &one, x, &one);
      }

      if (converged) {

        /* compute true residual using 'v[kspace]' as a temporary vector */

        AZ_scale_true_residual(val, indx, bindx, rpntr, cpntr, bpntr, x, b,
                               v[kspace], weight, &actual_residual,
                               &true_scaled_r, options, data_org, proc_config);

        converged = true_scaled_r < params[AZ_tol];

        if (!converged && (AZ_get_new_eps(&epsilon, scaled_r_norm,
                                          true_scaled_r,
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

        /* restore previous solution */

        if ( (!converged) && (i != kspace) )
          for (k = 0; k < N; k++) x[k] = x[k] - temp[k];
      }

      if ( (i == kspace) && !converged) {
        if (r_avail)
          for (k = 0; k < N; k++) v[0][k] = res[k];
        else {
          AZ_matvec_mult(val, indx, bindx, rpntr, cpntr, bpntr, temp, v[kspace],
                         1, data_org);
          dscal_(&N, &first_2norm, v[0], &one);
          daxpy_(&N, &minusone, v[kspace], &one, v[0], &one);
        }
      }
    }
  }

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

  AZ_terminate_status_print(i, iter, status, rec_residual, params,
                            scaled_r_norm, actual_residual, options,
                            proc_config);

} /* AZ_pgmres */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_get_x_incr(double val[], int indx[], int bindx[], int rpntr[],
                   int cpntr[], int bpntr[], int options[], int data_org[], int
                   proc_config[], double params[], int i, double **hh, double
                   *rs, double *temp, double **v)

/*******************************************************************************

  This routine is normally invoked from GMRES and is used to compute the
  increment to the solution that should be applied as a result of solving the
  upper hessenberg system produces by GMRES.

  Author:          Ray Tuminaro, SNL, 1422
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

  x:               On input, contains the initial guess. On output contains the
                   solution to the linear system.

  options:         Determines specific solution method and other parameters.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see file params.txt).

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

  params:          Drop tolerance and convergence tolerance info.

  i:               The size of the Krylov subspace that has be generated.

  hh:              The upper triangular matrix produced by gmres (after reducing
                   the upper hessenberg matrix) after sweeping out i steps.

  rs:              The projected right hand side produced by gmres.

  temp:            vector of length data_org[AZ_num_int_unk] +
                   data_org[AZ_num_bord_unk].  On output to AZ_get_x_incr(),
                   temp contains the increment that must be added the current
                   approximation to get the new gmres approximate solution.

  v:               Orthogonal vectors produced by Gramm-Schmidt process that
                   span the Krylov space swept out.

*******************************************************************************/

{

  /* local variables */

  int    ii, k, k1, j, N;
  double t;
  int    one = 1;

  /**************************** execution begins ******************************/

  N = data_org[AZ_N_internal] + data_org[AZ_N_border];

  /* solve upper triangular system, compute solution */

  i--;  /* set i = Krylov subspace dimension */

  rs[i] /= hh[i][i];
  for (ii = 1; ii <= i; ii++) {
    k  = i - ii;
    k1 = k + 1;
    t  = rs[k];
    for (j = k1; j <= i; j++) t -= hh[k][j] * rs[j];
    rs[k] = t / hh[k][k];
  }

  /*
   * done with back substitution; form linear combination to get solution
   */

  for (j = 0; j < N; j++) temp[j] = 0.0;

  for (j = 0; j <= i; j++) {
    t = rs[j];
    daxpy_(&N, &t, v[j], &one, temp, &one);
  }

  AZ_precondition(val, indx, bindx, rpntr, cpntr, bpntr, temp, options,
                  data_org, proc_config, params);

} /* AZ_get_x_incr */
