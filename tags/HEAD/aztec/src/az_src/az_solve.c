/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile: az_solve.c,v $
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
static char rcsid[] = "$Id: az_solve.c,v 1.1.1.1 2001-01-30 20:59:12 paklein Exp $";
#endif


/*******************************************************************************
 * Copyright 1995, Sandia Corporation.  The United States Government retains a *
 * nonexclusive license in this software as prescribed in AL 88-1 and AL 91-7. *
 * Export of this program may require a license from the United States         *
 * Government.                                                                 *
 ******************************************************************************/


/*
 * All of these routines form the interface between the code drivers and the
 * user's choice of algorithm and PDE problem.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "az_aztec.h"


void AZ_solve(double x[], double b[], int options[], double params[],
              int indx[], int bindx[], int rpntr[], int cpntr[], int bpntr[],
              double val[], int data_org[], double status[], int proc_config[])

/*******************************************************************************

  Solve the system of equations given in the VBR format using an iterative
  method specified by 'options[AZ_solver]'. Store the result in 'x'.

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  x:               On input, contains the initial guess. On output contains the
                   solution to the linear system.

  b:               Right hand side of linear system.

  options:         Determines specific solution method and other parameters.

  params:          Drop tolerance and convergence tolerance info.

  indx,
  bindx,
  rpntr,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file params.txt).

  val:             Array containing the nonzero entries of the matrix (see file
                   params.txt).

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see file params.txt).

  status:          On output, indicates termination status:
                    0:  terminated normally.
                   -1:  maximum number of iterations taken without achieving
                        convergence.
                   -2:  Breakdown. The algorithm can not proceed due to
                        numerical difficulties (usually a divide by zero).
                   -3:  Internal residual differs from the computed residual due
                        to a significant loss of precision.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  /* local variables */

  int     i;
  int     gn;
  double  gnnz;
  double  start_t, total_time, MFlops;
  double *dummy = NULL;
  int     nonzeros;
  int     N, N_Blk, Proc, Num_Proc;
  int     old_pre_calc;

  /**************************** execution begins ******************************/

  N_Blk        = data_org[AZ_N_int_blk] + data_org[AZ_N_bord_blk];
  N            = data_org[AZ_N_internal] + data_org[AZ_N_border];
  Proc         = proc_config[AZ_node];
  Num_Proc     = proc_config[AZ_N_procs];
  old_pre_calc = options[AZ_pre_calc];

  /* output solver, scaling, and preconditioning options */

  AZ_print_call_iter_solve(options, Proc);

  /* check for inconsistent user input options */

  if (!AZ_check_options(options, Proc, data_org, Num_Proc, params)) {
    status[AZ_its]      = (double )  0.0;
    status[AZ_why]      = (double )  AZ_param;
    status[AZ_r]        = (double ) -1.0;
    status[AZ_rec_r]    = (double ) -1.0;
    status[AZ_scaled_r] = (double ) -1.0;

    (void) AZ_manage_memory(0, AZ_CLEAR, AZ_SYS, (char *) 0, (int *) 0);

    return;
  }

  /* If desired, print out the matrix and indexing arrays */

  if (options[AZ_output] == AZ_all)
     AZ_print_out((int *) NULL, (int *) NULL, (int *) NULL, (int *) NULL, 
        val, indx, bindx, rpntr, cpntr, bpntr, proc_config, AZ_input_form, 
        data_org[AZ_matrix_type], data_org[AZ_N_int_blk]+data_org[AZ_N_bord_blk], 
        data_org[AZ_N_ext_blk], 0);

  /* adjust print frequency according to user input */

  i = options[AZ_output];

  if  (i == AZ_warnings) options[AZ_print_freq] = options[AZ_max_iter] + 10;
  else if (i == AZ_none) options[AZ_print_freq] = options[AZ_max_iter] + 10;
  else if (i == AZ_all ) options[AZ_print_freq] = 1;
  else if (i == AZ_last) options[AZ_print_freq] = options[AZ_max_iter] + 1;
  else                   options[AZ_print_freq] = options[AZ_output];

  /* scale matrix, rhs and initial guess if required */

  AZ_scale_f(val, indx, bindx, rpntr, cpntr, bpntr, b, x, options, data_org,
             proc_config, AZ_SCALE_MAT);

  /* solve the system */

  AZ_sync(Proc, Num_Proc);
  start_t = AZ_second();

  fflush(stdout);
  switch (options[AZ_solver]) {

  case AZ_cg:

    /* conjugate gradient */

    AZ_pcg_f(val, indx, bindx, rpntr, cpntr, bpntr, b, x, &(params[AZ_weights]),
             options, params, data_org, proc_config, status);
    break;

  case AZ_gmres:

    /* GMRES */

    AZ_pgmres(val, indx, bindx, rpntr, cpntr, bpntr, b, x,
              &(params[AZ_weights]), options, params, data_org, proc_config,
              status);
    break;

  case AZ_cgs:

    /* conjugate gradient squared */

    AZ_pcgs(val, indx, bindx, rpntr, cpntr, bpntr, b, x, &(params[AZ_weights]),
            options, params, data_org, proc_config, status);
    break;

  case AZ_tfqmr:

    /* transpose-free quasi minimum residual */

    AZ_pqmrs(val, indx, bindx, rpntr, cpntr, bpntr, b, x, &(params[AZ_weights]),
             options, params, data_org, proc_config,status);
    break;

  case AZ_bicgstab:

    /* stabilized bi-conjugate gradient */

    AZ_pbicgstab(val, indx, bindx, rpntr, cpntr, bpntr, b, x,
                 &(params[AZ_weights]), options, params, data_org, proc_config,
                 status);
    break;

  case AZ_symmlq:

    /* conjugate gradient squared */

#ifdef eigen
    AZ_psymmlq(val,indx, bindx,rpntr,cpntr,bpntr, b, x, &(params[AZ_weights]),
            options, params, data_org, proc_config, status);
#else
    printf("symmlq not implemented in this version\n");
#endif
    break;

  case AZ_lu:

    /* direct sparse LU */

    for (i = 0; i < N + data_org[AZ_N_external]; i++) x[i] = b[i];

    if (data_org[AZ_matrix_type] == AZ_MSR_MATRIX) nonzeros = bindx[N];
    else nonzeros = indx[bpntr[N_Blk]];

    i                   = options[AZ_overlap];
    options[AZ_overlap] = AZ_none;
    AZ_lu_y12m(N, nonzeros, x, 0, dummy, val, indx, bindx, rpntr, cpntr, bpntr,
               options, data_org, 0.0);
    options[AZ_overlap] = i;
    status[AZ_its]      = (double ) 1.0;
    status[AZ_why]      = (double ) AZ_normal;
    status[AZ_r]        = (double ) 0.0;
    status[AZ_rec_r]    = (double ) 0.0;
    status[AZ_scaled_r] = (double ) 0.0;
    break;

  default:
    (void) fprintf(stderr,"ERROR: options[AZ_solver] has improper value(%d)\n",
                   options[AZ_solver]);
  exit(-1);
  }
  fflush(stdout);

  total_time = AZ_gmax_double(AZ_second() - start_t, proc_config);

  AZ_scale_f(val, indx, bindx, rpntr, cpntr, bpntr, b, x, options, data_org,
             proc_config, AZ_RESCALE_SOL);

  (void) AZ_manage_memory(0, AZ_CLEAR, AZ_SYS, (char *) 0, (int *) 0);

  options[AZ_pre_calc] = old_pre_calc;

  if ((options[AZ_output] != AZ_none) && (options[AZ_output] != AZ_warnings)) {
    if (proc_config[AZ_node] == 0) {
      (void) printf("\n\n\t\tSolution time: %f (sec.)\n", total_time);
      (void) printf("\t\ttotal iterations: %d\n", (int) status[AZ_its]);
    }
  }

  /*
   * calculate the solver MFlop/s
   */

#ifdef AZ_FLOP_CNTS
#endif
  if ( (options[AZ_output] != AZ_none) && (options[AZ_output] != AZ_warnings)) {
    gn = AZ_gsum_int(N, proc_config);
    if (data_org[AZ_matrix_type] == AZ_VBR_MATRIX)
      gnnz = AZ_gsum_double((double) (indx[bpntr[N_Blk]]-indx[0]), proc_config);
    else
      gnnz = AZ_gsum_double((double) (bindx[N]), proc_config);

    if (proc_config[AZ_node] == 0) {
      total_time += 1.0e-10;
      MFlops = AZ_calc_solve_flops(options, (int) status[AZ_its],total_time,gn,
                                   gnnz, data_org, proc_config);
      if (MFlops > 0.0) {
        (void) printf("\t\tSolver MFlop rate: %f MFlops/sec.\n\t\t", MFlops);
        (void) printf("Solver processor MFlop rate: %f MFlops/sec.\n\n",
                      MFlops/Num_Proc);
      }
    }
  }

#ifdef TIME_VB

  /* calculate the individual kernel times and Mflops */

  AZ_time_kernals(gn, gnnz, val, indx, bindx, rpntr, rpntr, bpntr, x, b,
                  (int) status[AZ_NUMBER_OF_ITS], options, data_org,
                  proc_config);

#endif

  fflush(stdout);

} /* AZ_solve */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_print_call_iter_solve(int options[], int az_proc)

/*******************************************************************************

  Print out the type of solver called, scaling and preconditioning information.

  Author:          SNL
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  options:         Determines specific solution method and other parameters.

  az_proc:         Current processor.

*******************************************************************************/

{
  int super_special = 0;

  /**************************** execution begins ******************************/

  if ( (options[AZ_output] == AZ_last) || (options[AZ_output] == AZ_none) ||
       (options[AZ_output] == AZ_warnings) ||  (az_proc != 0) ) return;

  (void) printf("\n\t\t****************************************"
                "***************\n\t\t***** ");

  /* first, print out chosen solver */

  switch (options[AZ_solver]) {

  case AZ_cg:
    (void) printf("Preconditioned CG"); break;
  case AZ_gmres:
    (void) printf("Preconditioned GMRES"); break;
  case AZ_cgs:
    (void) printf("Preconditioned CGS"); break;
  case AZ_tfqmr:
    (void) printf("Preconditioned TFQMR"); break;
  case AZ_bicgstab:
    (void) printf("Preconditioned BICGSTAB"); break;
  case AZ_symmlq:
    (void) printf("Preconditioned SYMMLQ-like"); break;
  case AZ_lu:
    (void) printf("LU");
  }

  (void) printf(" solution\n");

  /* next output preconditioning options */

  (void) printf("\t\t***** ");

  switch (options[AZ_precond]) {
  case AZ_none:
    (void) printf("No"); break;
  case AZ_Neumann:
    (void) printf("Order %d Neumann series polynomial", options[AZ_poly_ord]);
  break;
  case AZ_ls:
    (void) printf("Order %d least-squares polynomial", options[AZ_poly_ord]);
  break;
  case AZ_Jacobi:
    (void) printf("%d step block Jacobi", options[AZ_poly_ord]);
  break;
  case AZ_sym_GS:
    (void) printf("%d step symmetric Gauss-Seidel", options[AZ_poly_ord]);
  break;
  case AZ_bilu:
    (void) printf("BILU domain decomp. ");
  if      (options[AZ_overlap] == AZ_none) printf("without overlap");
  else if (options[AZ_overlap] == AZ_diag) printf("with diagonal overlap");
  else if (options[AZ_overlap] == AZ_full) printf("with full overlap");
  else if (options[AZ_overlap] == AZ_sym_full) printf("with sym. full overlap");
  super_special = 1;
  break;
  case AZ_ilu:
    (void) printf("ILU domain decomp. ");
  if      (options[AZ_overlap] == AZ_none) printf("without overlap");
  else if (options[AZ_overlap] == AZ_diag) printf("with diagonal overlap");
  else if (options[AZ_overlap] == AZ_full) printf("with full overlap");
  else if (options[AZ_overlap] == AZ_sym_full) printf("with sym. full overlap");
  super_special = 1;
  break;
#ifdef eigen
  case AZ_slu:
    (void) printf("Ng/LU domain decomp. ");
  if      (options[AZ_overlap] == AZ_none) printf("without overlap");
  else {
    (void) printf("not implemented with overlap\n"); exit(-1);
  }
  super_special = 1;
  break;
#endif
  case AZ_lu:
    (void) printf("LU domain decomp. ");
  if      (options[AZ_overlap] == AZ_none) printf("without overlap");
  else if (options[AZ_overlap] == AZ_diag) printf("with diagonal overlap");
  else if (options[AZ_overlap] == AZ_full) printf("with full overlap");
  else if (options[AZ_overlap] == AZ_sym_full) printf("with sym. full overlap");
  super_special = 1;
  break;
  case AZ_icc:
    (void) printf("incomplete Choleski decomposition");
  super_special = 1;
  }

  (void) printf(" preconditioning\n");
  (void) printf("\t\t***** ");

  /* lastly, print out the scaling information */

  switch (options[AZ_scaling]) {

  case AZ_none:
    (void) printf("No"); break;
  case AZ_Jacobi:
    (void) printf("block Jacobi"); break;
  case AZ_BJacobi:
    (void) printf("block Jacobi"); break;
  case AZ_row_sum:
    (void) printf("left row-sum"); break;
  case AZ_sym_diag:
    (void) printf("symmetric diagonal"); break;
  case AZ_sym_row_sum:
    (void) printf("symmetric row sum");
  }

  (void) printf(" scaling\n");

  if (super_special==1) {
    (void) printf("\t\t***** NOTE: convergence VARIES when the total number "
                  "of\n");
    (void) printf("\t\t*****       processors is changed.\n");
  }

  (void) printf("\t\t****************************************"
                "***************\n\n");

} /* AZ_print_call_iter_solver */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_output_matrix(double val[], int indx[], int bindx[], int rpntr[],
                      int cpntr[], int bpntr[], int proc_config[],
                      int data_org[])

/*******************************************************************************

  Routine to perform full matrix dump on each processor.

  Author:          Scott A. Hutchinson, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  a:               Array containing the nonzero entries of the matrix.

  indx,
  bindx,
  rpntr,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see file params.txt).

*******************************************************************************/

{

  /* local variables */

  int  iblk_row, i, j, ib1, ib2, n1, iblk, jblk, m1, ipoint, jpoint;
  int  ival = 0;
  int  k,num_nonzeros;
  int  num_total_nodes, N_external_nodes;
  int  Proc, Num_Proc;
  char str[5];
  char nstr[40];

  /********** execution begins **********/

  Proc               = proc_config[AZ_node];
  Num_Proc           = proc_config[AZ_N_procs];
  N_external_nodes = data_org[AZ_N_external];

  if (data_org[AZ_matrix_type] == AZ_VBR_MATRIX) {
    num_total_nodes    = data_org[AZ_N_int_blk]+data_org[AZ_N_bord_blk];
    N_external_nodes = data_org[AZ_N_ext_blk];

    /* Print out the VBR indexing information for the matrix */

    AZ_print_sync_start(Proc, AZ_TRUE);

    (void) printf("\n----- Proc: %d indx -----\n\n", Proc);
    for (iblk = 0; iblk < num_total_nodes; iblk++) {
      for (i = *(bpntr+iblk); i < *(bpntr+iblk+1); i++)
        (void) printf("%d ", *(indx+i));

      if (iblk == num_total_nodes - 1)
        (void) printf("%d\n",*(indx+i));
      else
        (void) printf("\n");
    }

    (void) printf("\n----- Proc: %d bindx -----\n\n", Proc);
    for (iblk = 0; iblk < num_total_nodes; iblk++) {
      for (i = *(bpntr+iblk); i < *(bpntr+iblk+1); i++)
        (void) printf("%d ", *(bindx+i));
      (void) printf("\n");
    }

    (void) printf("\n----- Proc: %d rpntr -----\n\n", Proc);
    for (i = 0; i < num_total_nodes + 1; i++)
      (void) printf("%d ", *(rpntr+i));
    (void) printf("\n");

    (void) printf("\n----- Proc: %d cpntr -----\n\n", Proc);
    for (i = 0; i < num_total_nodes + N_external_nodes + 1; i++)
      (void) printf("%d ", *(cpntr+i));
    (void) printf("\n");

    (void) printf("\n----- Proc: %d bpntr -----\n\n", Proc);
    for (i = 0; i < num_total_nodes + 1; i++)
      (void) printf("%d ", *(bpntr+i));
    (void) printf("\n");

    AZ_print_sync_end(Proc, Num_Proc, AZ_TRUE);

    AZ_print_sync_start(Proc, AZ_TRUE);
    (void) printf("AZ_solve debug output - full matrix dump: Processor %d\n",
                  Proc);

    /* loop over block rows */

    for (iblk_row = 0; iblk_row < num_total_nodes; iblk_row++) {

      /* number of rows in the current row block */

      m1 = rpntr[iblk_row+1] - rpntr[iblk_row];

      /* starting index of current row block */

      ival = indx[bpntr[iblk_row]];

      /* loop over all the blocks in the current block-row */

      for (j = bpntr[iblk_row]; j < bpntr[iblk_row+1]; j++) {
        jblk = bindx[j];

        /* the starting point column index of the current block */

        ib1 = cpntr[jblk];

        /* ending point column index of the current block */

        ib2 = cpntr[jblk+1];

        /* number of columns in the current block */

        n1 = ib2 - ib1;

        (void) printf("\nProc: %d Block Row: %d Block Column: %d "
                      "Row Pointer: %d Column Pointer: %d\n", Proc, iblk_row,
                      jblk, rpntr[iblk_row], rpntr[jblk]);
        (void) printf("----------------------------------------"
                      "----------------------------------------\n");

                      for (ipoint = 0; ipoint < m1; ipoint++) {
                        for (jpoint = 0; jpoint < n1; jpoint++)
                          (void) printf("a[%d]: %e ", ival+jpoint*m1+ipoint,
                                        val[ival+jpoint*m1+ipoint]);
                        (void) printf("\n");
                      }

                      ival += m1*n1;
      }
    }

    AZ_print_sync_end(Proc, Num_Proc, AZ_TRUE);
  }

  if (data_org[AZ_matrix_type] == AZ_MSR_MATRIX) {
    num_total_nodes = data_org[AZ_N_internal]+data_org[AZ_N_border];

    N_external_nodes = data_org[AZ_N_external];
    i                  = num_total_nodes + N_external_nodes;

    if      (i < 10    ) strcpy(str, "%1d");
    else if (i < 100   ) strcpy(str, "%2d");
    else if (i < 1000  ) strcpy(str, "%3d");
    else if (i < 10000 ) strcpy(str, "%4d");
    else if (i < 100000) strcpy(str," %5d");
    else strcpy(str, "%d");
    sprintf(nstr, "a(,%s)=%%8.1e ", str);

    AZ_print_sync_start(Proc, AZ_TRUE);

    (void) printf("\n----- Proc: %d -----\n\n", Proc);

    num_nonzeros = bindx[num_total_nodes];
    (void) printf("val:  ");
    for (i = 0; i < num_nonzeros; i++) {
      (void) printf("%9.1e", val[i]);
      if ((i%8) == 7) (void) printf("\n    ");
    }

    (void) printf("\nbindx:");
    for (i = 0; i < num_nonzeros; i++) {
      (void) printf("%9d", bindx[i]);
      if ((i%8) == 7) (void) printf("\n    ");
    }
    (void) printf("\n");

    for (i = 0; i < num_total_nodes; i++ ) {
      (void) printf("\nrow", Proc);
      (void) printf(str, i);
      (void) printf(":");
      (void) printf(nstr, i, val[i]);
      k = 0;
      for (j = bindx[i]; j < bindx[i+1]; j++ ) {
        (void) printf(nstr, bindx[j], val[j]);
        k++;
        if (((k%4) == 3) && (j != bindx[i+1]-1))
          (void) printf("\n      ");
      }
    }

    (void) printf("\n");
    AZ_print_sync_end(Proc, Num_Proc, AZ_TRUE);
  }

} /* AZ_output_matrix */
