/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile: az_symmlq.c,v $
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
static char rcsid[] = "$Id: az_symmlq.c,v 1.1.1.1 2001-01-30 20:59:12 paklein Exp $";
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

void AZ_psymmlq(double val[], int indx[], int bindx[], int rpntr[], int cpntr[],
        int bpntr[], double b[], double x[], double weight[], int options[],
        double params[], int data_org[], int proc_config[], double status[])

/*******************************************************************************

  A symmlq-like algorithm to solve symmetric indefinite problems. This
  code was actually modified from R. Freund's complex  QMR where the 
  complex stuff was taken out.

  Author: 	Lydie Prevost Div 1422 SNL
  =======

  Return code:     void
  ============

  Paramter list:
  ==============

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

  
  double 	   beta_k, beta_kp1, inv_beta, scal;
  double           omega_km1, omega_k, omega_kp1;
  double           norm_xi_tilda, norm_xi;
  double           c_km2, c_km1, c_k;
  double 	   alpha, teta,minus_teta,  minus_nu, nu;
  double 	   xi, xi_tilda, tau, tau_tilda;
  double 	   s_km2, s_km1, s_k;
  double           skmag2;
  double 	   temp_tau;
  double 	   inv_xi;
  int              iter, Nsize;
  int              i, j, one = 1;
  double           d_one = 1.0;
  double           epsilon;
  double 	   *t_k, *t_kp1, *z;
  double 	   *v_k, *v_km1, *v_kp1, *a_v_k, *p_km2, *p_km1, *p_k;
  double 	   *residual,  *u;
int N;
int precond_flag, proc, print_freq, r_avail = AZ_TRUE;
double rec_residual,scaled_r_norm, r_z_dot;
int converged;
double *garbage, actual_residual, true_scaled_r;



  /**************************** execution begins ******************************/

  /* pull needed values out of parameter arrays */
 
  N            = data_org[AZ_N_internal] + data_org[AZ_N_border];
  precond_flag = options[AZ_precond];
  epsilon      = params[AZ_tol];
  proc         = proc_config[AZ_node];
  print_freq   = options[AZ_print_freq];

  if (proc == 0) {
     (void) fprintf(stderr,"Warning: Indefinte preconditioners can cause\n");
     (void) fprintf(stderr,"         this symmlq like algorithm to crash\n");
  }
  /* initialize */

  Nsize = (N + data_org[AZ_N_external])*sizeof(double) + 1; 

  z       = (double *) AZ_manage_memory(Nsize,AZ_ALLOC,AZ_SYS,"z in sqmr",&j);
  t_k     = (double *) AZ_manage_memory(Nsize,AZ_ALLOC,AZ_SYS,"t_k in sqmr",&j);
  t_kp1   = (double *) AZ_manage_memory(Nsize,AZ_ALLOC,AZ_SYS,"t_kp1 in sqmr",&j);
  v_k     = (double *) AZ_manage_memory(Nsize,AZ_ALLOC,AZ_SYS,"v_k in sqmr",&j);
  v_km1   = (double *) AZ_manage_memory(Nsize,AZ_ALLOC,AZ_SYS,"v_km1 in sqmr",&j);
  v_kp1   = (double *) AZ_manage_memory(Nsize,AZ_ALLOC,AZ_SYS,"v_kp1 in sqmr",&j);
  a_v_k   = (double *) AZ_manage_memory(Nsize,AZ_ALLOC,AZ_SYS,"a_v_k in sqmr",&j);
  p_km2   = (double *) AZ_manage_memory(Nsize,AZ_ALLOC,AZ_SYS,"p_km2 in sqmr",&j);
  p_km1   = (double *) AZ_manage_memory(Nsize,AZ_ALLOC,AZ_SYS,"p_km1 in sqmr",&j);
  p_k     = (double *) AZ_manage_memory(Nsize,AZ_ALLOC,AZ_SYS,"p_k in sqmr",&j);
  residual= (double *) AZ_manage_memory(Nsize,AZ_ALLOC,AZ_SYS,"residual in sqmr",&j);
  u       = (double *) AZ_manage_memory(Nsize,AZ_ALLOC,AZ_SYS,"u in sqmr",&j);
  garbage = (double *) AZ_manage_memory(Nsize,AZ_ALLOC,AZ_SYS,"garb in sqmr",&j);
  Nsize   = N + data_org[AZ_N_external];
  for (i = 0 ; i < Nsize ; i++ ) {
     t_k[i] = 0.0; t_kp1[i] = 0.0; v_k[i] = 0.0; v_km1[i] = 0.0; v_kp1[i] = 0.0;
     a_v_k[i] = 0.0; p_km2[i] = 0.0; p_km1[i] = 0.0; p_k[i] = 0.0; 
  }

  /* form residual: u  = b - A*x */

  AZ_compute_residual(val, indx, bindx, rpntr, cpntr, bpntr, b, x, u, data_org);

  dcopy_(&N, u, &one, residual, &one);

  /* compute z :  u = M*z  */

  status[AZ_first_precond] = AZ_second();
  dcopy_(&N, u, &one, z, &one);
  if (precond_flag) AZ_precondition(val, indx, bindx, rpntr, cpntr, bpntr, z, 
			            options, data_org, proc_config, params);
  status[AZ_first_precond] = AZ_second() - status[AZ_first_precond];

  /* compute a few global scalars:                                 */
  /*     1) ||r||                corresponding to options[AZ_conv] */
  /*     2) scaled ||r||         corresponding to options[AZ_conv] */
  /*     3) compute beta_k =     ( z .dot. u)                      */
 
  AZ_compute_global_scalars(val, indx, bindx, rpntr, cpntr, bpntr, x, b, residual, 
			   weight, &rec_residual, &scaled_r_norm, options, data_org,
                           proc_config, &r_avail, residual, z, &beta_k, AZ_FIRST_TIME);
  if ( beta_k < 0.0 ) {
     if (proc == 0) {
        (void)fprintf(stderr,"Error: Indefinite preconditioner has caused");
        (void)fprintf(stderr," problems.\n");
     }
     exit(-1);
  }
  beta_k = sqrt(beta_k);
 

  /* initialize vectors and scalers : vector v_km1 = 0
                                      vector p_km1 = 0
                                      c_km2 = 1
                                      c_km1 = 1
                                      s_km2 = 0
                                      s_km1= 0   */

  for (i = 0 ; i < N ; i++ ) v_km1[i]= 0.;
  dcopy_(&N, v_km1, &one, p_km2, &one);
  dcopy_(&N, v_km1, &one, p_km1, &one);

  omega_km1 = 0.;
  c_km2 = d_one;
  c_km1= d_one;
  s_km2= 0.;
  s_km1= 0.;


     
  /* compute v_k = u / beta_k */

  inv_beta = 1/beta_k;
  dcopy_(&N, u, &one, v_k, &one);
  dscal_(&N,&inv_beta, v_k, &one);

  /* compute t_k = z / beta_k */

  dcopy_(&N, z, &one, t_k, &one);
  dscal_(&N,&inv_beta, t_k, &one);

  /* compute omega_k = || v_k || and tau_tilda = omega_k * beta_k */

  omega_k = sqrt( AZ_gdot(N, v_k, v_k, proc_config) );
  tau_tilda = omega_k * beta_k;


  if ((options[AZ_output] != AZ_none) &&
      (options[AZ_output] != AZ_last) &&
      (options[AZ_output] != AZ_warnings) && (proc == 0))
    (void) fprintf(stdout, "\t\titer:    0\t\tresidual = %e\n", scaled_r_norm);
 
  converged = scaled_r_norm < epsilon;

  /* iteration loop */

  for(iter = 1; iter <= options[AZ_max_iter] && !converged; iter++ ) {

    /* compute  alpha_k = t_k . dot. A t_k  */

    AZ_matvec_mult(val, indx, bindx, rpntr, cpntr, bpntr, t_k, a_v_k, 1, data_org);

    alpha = AZ_gdot(N, t_k, a_v_k, proc_config);

    /* compute u = Atk - alpha*v_k - beta_k * v_km1  */
     
    dcopy_(&N, a_v_k, &one, u, &one);
    scal = - beta_k;
    daxpy_(&N, &scal, v_km1, &one, u, &one);
    scal = -alpha;
    daxpy_(&N, &scal, v_k, &one, u, &one);

    /* compute z : u = M*z */

    dcopy_(&N, u, &one, z, &one);
    if (precond_flag) AZ_precondition(val, indx, bindx, rpntr, cpntr, bpntr, z, 
			            options, data_org, proc_config, params);

    /* compute beta_kp1 = square root ( z .dot. u)  */

     beta_kp1 = sqrt( AZ_gdot(N, z, u, proc_config) );

     /* compute v_kp1 = u / beta_kp1 */
     /* compute t_kp1 = z / beta_kp1 */

     if ( beta_kp1 == 0. ) inv_beta = 1.;
     else inv_beta = 1./beta_kp1;

     dcopy_(&N, u, &one, v_kp1, &one);
     dscal_(&N,&inv_beta, v_kp1, &one);
     dcopy_(&N, z, &one, t_kp1, &one);
     dscal_(&N,&inv_beta, t_kp1, &one);
    
     /* compute omega_kp1 = || v_kp1 ||  */
       
     omega_kp1 = sqrt( AZ_gdot(N, v_kp1, v_kp1, proc_config) );

     /* compute teta = conjug(s_km2) * omega_km1 * beta_k
                nu   = c_km1 * c_km2 * omega_km1 * beta_k 
                     + conjug(s_km1) * omega_k * alpha     */      
     
     teta     = omega_km1 * s_km2 * beta_k;

     nu       = omega_k * s_km1*alpha + omega_km1 * c_km1 * c_km2 * beta_k;

     xi_tilda = c_km1 * omega_k * alpha - c_km2 * omega_km1 * s_km1 * beta_k;

     norm_xi_tilda = fabs(xi_tilda);
     norm_xi =  norm_xi_tilda*norm_xi_tilda 
                + ( omega_kp1*omega_kp1  * ( beta_kp1*beta_kp1) ) ;
     norm_xi = sqrt(norm_xi) ;
     if ( norm_xi_tilda != 0.0 )  { 
        xi  =  norm_xi * xi_tilda / norm_xi_tilda  ;
        c_k = norm_xi_tilda / norm_xi; }
     else {
        xi = norm_xi;
        c_k = 0;}
 
     /* compute s_k = (omega_kp1 * beta_k+1 ) / xi   */

     s_k =  beta_kp1/xi; 
     s_k = omega_kp1 * s_k;


     /* compute 
        vector p_k = (t_k -  nu * p_km1 - teta * p_km2 ) / xi    */

     dcopy_(&N,t_k, &one,p_k,&one);
     minus_nu = - nu;
     daxpy_(&N, &minus_nu, p_km1, &one, p_k, &one);
     minus_teta = - teta;
     daxpy_(&N, &minus_teta,p_km2 , &one, p_k, &one);
     inv_xi = 1.0/xi;
     dscal_(&N,&inv_xi,p_k, &one);

     /* compute tau = c_k * tau_tilda   */

     tau = c_k * tau_tilda;

     /* compute tau_tilda = -s_k * tau_tilda   */

     tau_tilda = s_k * tau_tilda;
     tau_tilda = - tau_tilda;


     /* compute x = x + tau * p_k   */

     daxpy_(&N, &tau, p_k, &one, x, &one);


    /* check convergence criteria */
    /* compute the current residual */
    /* residual = abs(s_k)**2 * residual + c_k*tau_tilda*v_kp1/omega_kp1  */


      skmag2 =  s_k * s_k;

      for (i = 0 ; i < N ; i++ ){
         temp_tau = tau_tilda * v_kp1[i];
         residual[i]  = skmag2 * residual[i] + c_k * temp_tau/ omega_kp1;
      }

    /* compute a few global scalars:                                 */
    /*     1) ||r||                corresponding to options[AZ_conv] */
    /*     2) scaled ||r||         corresponding to options[AZ_conv] */
    /*     3) r_z_dot = <z, r>                                       */
 
    AZ_compute_global_scalars(val, indx, bindx, rpntr, cpntr, bpntr, x, b, residual,
                           weight, &rec_residual, &scaled_r_norm, options,
                           data_org, proc_config, &r_avail, residual, z, &r_z_dot,
                           AZ_NOT_FIRST);


    /* save necessary old scalers */

    s_km2      = s_km1;
    s_km1     = s_k;
    beta_k    = beta_kp1;
    c_km2     = c_km1 ;
    c_km1     = c_k ;
    omega_km1 = omega_k ;
    omega_k   = omega_kp1 ;

    /* save necessary old vectors */

    dcopy_(&N, v_k, &one, v_km1, &one);
     dcopy_(&N, v_kp1, &one, v_k, &one);
     dcopy_(&N, t_kp1, &one, t_k, &one);
     dcopy_(&N, p_km1, &one, p_km2, &one);
     dcopy_(&N, p_k, &one, p_km1, &one);

    if ( (iter%print_freq == 0) && proc == 0 )
      (void) fprintf(stdout, "\t\titer: %4d\t\tresidual = %e\n", iter,
                     scaled_r_norm);
 
    /* convergence tests */

     

    if (scaled_r_norm < epsilon) {
      AZ_scale_true_residual(val, indx, bindx, rpntr, cpntr, bpntr, x, b, garbage,
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
  if ( (iter%print_freq != 0) && (proc == 0) && (options[AZ_output] != AZ_none))
    (void) fprintf(stdout, "\t\titer: %4d\t\tresidual = %e\n", iter,
                   scaled_r_norm);
 
  /* check if we exceeded maximum number of iterations */
 
  if (converged) { i = AZ_normal; scaled_r_norm = true_scaled_r; }
  else           i = AZ_maxits;
 
  AZ_terminate_status_print(i, iter, status, rec_residual, params, 
			scaled_r_norm , actual_residual, options, proc_config);

} /* pqmrs */

