/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile: az_precond.c,v $
 *
 * $Author: paklein $
 *
 * $Date: 2001-01-30 20:59:11 $
 *
 * $Revision: 1.1.1.1 $
 *
 * $Name: not supported by cvs2svn $
 *====================================================================*/
#ifndef lint
static char rcsid[] = "$Id: az_precond.c,v 1.1.1.1 2001-01-30 20:59:11 paklein Exp $";
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
#include <string.h>
#include "az_aztec.h"

/* missing prototypes - PAK (07/24/98) */
void jacobi(double val[], double x[], int data_org[]);

#ifdef eigen
extern void AZ_do_Jacobi(double val[], int indx[], int bindx[], int rpntr[],
                     int cpntr[], int bpntr[], double x[], double b[],
                     double temp[], int options[], int data_org[],
                     int proc_config[], double params[], int flag);
#endif
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_precondition(double val[], int indx[], int bindx[], int rpntr[],
                     int cpntr[], int bpntr[], double x[], int options[],
                     int data_org[], int proc_config[], double params[])

/*******************************************************************************

  This routine calls appropriate sparse matrix preconditioner.

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

  x:               On input, contains the current solution. On output contains
                   the preconditioned solution to the linear system.

  options:         Determines specific solution method and other parameters.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see file params.txt).

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

  params:          Drop tolerance and convergence tolerance info.

 * --------------------------------------------------------------------

 Related routines(VBR sparse format):

   scaling routines:
        AZ_block_diagonal_scaling -- block-diagonally scales sparse matrix
                                     problem.
        AZ_row_sum_scaling        -- row sum scales sparse matrix problem.
        sym_diagonal_scaling      -- diagonaly scales symm. sparse problem.
        sym_row_sum_scaling       -- row sum scales symmetric sparse problem.

   preconditioners:
        jacobi                 -- point Jacobi method.
        AZ_polynomial_expansion-- Polynomial expansion; Neumann series and
                                  least squares.
        domain decomposition   -- Block solvers (LU or ILU) used on each
                                  processor. The blocks are either
                                  non-overlapping or overlapping.
        icc                    -- incomplete sparse Choleski (symmetric
                                  version).

*******************************************************************************/

{

  /* local variables */

  int            ione = 1;
  double        *temp;
  int            m, N;
  int            i, step;
  static int    *d2_indx,*d2_bindx,*d2_rpntr,*d2_bpntr;
  static double *d2_inv;
  int            tsize;
  static int     previous_factors = -1;
  double        *v, *y;
  char          *yo = "precond: ";
#ifdef eigen
  double         *tb, *tr;
#endif

  /* prototype */

  void AZ_calc_blk_diag_inv(double *val, int *indx, int *bindx, int *rpntr,
                            int *cpntr, int *bpntr, double *d_inv, int *d_indx,
                            int *d_bindx, int *d_rpntr, int *d_bpntr,
                            int data_org[]);

  /**************************** execution begins ******************************/

  m    = data_org[AZ_N_int_blk] + data_org[AZ_N_bord_blk];
  N    = data_org[AZ_N_internal] + data_org[AZ_N_border];


  switch (options[AZ_precond]) {
  case AZ_none:

    /* no preconditioning, just return */

    break;

  case AZ_Jacobi:

    if (data_org[AZ_matrix_type] == AZ_MSR_MATRIX) {
      for (i = 0; i < N; i++)
        x[i] /= val[i];

      if (options[AZ_poly_ord] > 1) {
        v = AZ_manage_memory((N+data_org[AZ_N_external])*sizeof(double),
                             AZ_ALLOC, AZ_SYS, "v in precond", &i);
        y = AZ_manage_memory(N*sizeof(double), AZ_ALLOC, AZ_SYS,
                             "y in precond",&i);

        for (i = 0; i < N; i++)
          v[i] = x[i];

        for (step = 1; step < options[AZ_poly_ord]; step++) {
          AZ_matvec_mult(val, indx, bindx, rpntr, cpntr, bpntr, v, y, 1,
                         data_org);

          for(i = 0; i < N; i++)
            v[i] += x[i] - y[i] / val[i];
        }

        for (i = 0; i < N; i++)
          x[i] = v[i];
      }
    }

    else {

      /* block Jacobi preconditioning */

      if (options[AZ_pre_calc] < AZ_sys_reuse) {

        /*
         * First, compute the block-diagonal inverse (only if it hasn't already
         * been done).
         */

        tsize = 0;
        for (i = 0; i < m; i++)
          tsize += (rpntr[i+1] - rpntr[i]) * (cpntr[i+1] - cpntr[i]);

        d2_indx  = (int *) AZ_manage_memory((m+1)*sizeof(int), AZ_ALLOC,
                                            data_org[AZ_name],
                                            "d2_indx", &i);
        d2_bindx = (int *) AZ_manage_memory(m*sizeof(int), AZ_ALLOC,
                                            data_org[AZ_name],
                                            "d2_bindx", &i);
        d2_rpntr = (int *) AZ_manage_memory((m+1)*sizeof(int), AZ_ALLOC,
                                            data_org[AZ_name],
                                            "d2_rpntr", &i);
        d2_bpntr = (int *) AZ_manage_memory((m+1)*sizeof(int), AZ_ALLOC,
                                            data_org[AZ_name],
                                            "d2_bpntr", &i);
        d2_inv   = (double *) AZ_manage_memory(tsize*sizeof(double), AZ_ALLOC,
                                               data_org[AZ_name],
                                               "d2_inv", &i);

        if (options[AZ_pre_calc] != AZ_reuse) {
          AZ_calc_blk_diag_inv(val, indx, bindx, rpntr, cpntr, bpntr, d2_inv,
                               d2_indx, d2_bindx, d2_rpntr, d2_bpntr, data_org);
        }
        else if (i == AZ_NEW_ADDRESS) {
          (void) fprintf(stderr, "Error: options[AZ_pre_calc]==AZ_reuse and"
                         "previous factors\n       not found. Check"
                         "data_org[AZ_name].\n");
          exit(-1);
        }
      }

      else if (previous_factors != data_org[AZ_name]) {
        (void) fprintf(stderr, "Warning: Using a previous factorization as a"
                       "preconditioner\neven though matrix"
                       "(data_org[AZ_name]) has changed\n");
      }

      previous_factors = data_org[AZ_name];

      /* scale rhs */

      v = AZ_manage_memory((N+data_org[AZ_N_external])*sizeof(double),
                           AZ_ALLOC, AZ_SYS, "v in precond", &i);

      AZ_matvec_mult(d2_inv, d2_indx, d2_bindx, d2_rpntr, d2_rpntr, d2_bpntr,
                     x, v, 1, data_org);

#if defined (hp)
      vec_$dcopy(v, x, &N);
#else
      dcopy_(&N, v, &ione, x, &ione);
#endif

      if (options[AZ_poly_ord] > 1) {
        y = AZ_manage_memory((N+data_org[AZ_N_external])*sizeof(double),
                             AZ_ALLOC, AZ_SYS, "y in precond", &i);

        temp = AZ_manage_memory(N*sizeof(double), AZ_ALLOC, AZ_SYS,
                                "temp in precond", &i);

        for (step = 1; step < options[AZ_poly_ord]; step++) {
          AZ_matvec_mult(val, indx, bindx, rpntr, cpntr, bpntr, v, y, 1,
                         data_org);
          AZ_matvec_mult(d2_inv, d2_indx, d2_bindx, d2_rpntr, d2_rpntr,
                         d2_bpntr, y, temp, 1, data_org);

          for (i = 0; i < N; i++)
            v[i] += x[i] - temp[i];
        }

        for (i = 0; i < N; i++)
          x[i] = v[i];
      }
    }
    break;

  case AZ_sym_GS:

    /* symmetric Gauss-Seidel preconditioner only available on 1 proc */

    if (data_org[AZ_matrix_type] == AZ_VBR_MATRIX)
      AZ_sym_gauss_seidel();
    else
      AZ_sym_gauss_seidel_sl(val, bindx, x, data_org, options);
    break;

  case AZ_Neumann:
  case AZ_ls:

    /* polynomial preconditioners */

    if (!options[AZ_poly_ord]) return;

    AZ_polynomial_expansion(val, indx, bindx, rpntr, cpntr, bpntr, x, options,
                            data_org, proc_config);
    break;

  case AZ_ilu:
  case AZ_lu:
  case AZ_bilu:
    AZ_domain_decomp(x, val, indx, rpntr,cpntr,bindx,bpntr, options, data_org,
                     proc_config, params);
    break;

#ifdef eigen
  case AZ_slu:
    if (options[AZ_poly_ord] != 0) {
       tb = AZ_manage_memory((N+data_org[AZ_N_external])*sizeof(double),
                             AZ_ALLOC, AZ_SYS, "tb in precond", &i);
       tr = AZ_manage_memory((N+data_org[AZ_N_external])*sizeof(double),
                             AZ_ALLOC, AZ_SYS, "tr in precond", &i);
       for (i = 0 ; i < N ; i++ ) tb[i] = x[i];
       for (i = 0 ; i < N ; i++ ) x[i] = 0.0;
       AZ_do_Jacobi(val, indx, bindx, rpntr, cpntr, bpntr, x, tb, tr, options,
                 data_org, proc_config, params, 1);
       AZ_compute_residual(val, indx, bindx, rpntr, cpntr, bpntr, tb, x, tr,
                      data_org);
       AZ_domain_decomp(tr, val, indx, rpntr,cpntr,bindx,bpntr, options, 
                      data_org, proc_config, params);
       for (i = 0 ; i < N ; i++ ) x[i] += tr[i];
 
       /* do symmetric version */

       if ((options[AZ_solver] == AZ_cg) || (options[AZ_solver] == AZ_symmlq)) 
          AZ_do_Jacobi(val, indx, bindx, rpntr, cpntr, bpntr, x, tb, tr,
                       options,data_org, proc_config, params, 0);
   }
   else 
    AZ_domain_decomp(x, val, indx, rpntr,cpntr,bindx,bpntr, options, data_org,
                     proc_config, params);
    break;

#endif

  case AZ_icc:

    /* incomplete Cholesky factorization */

    (void) printf("not currently hooked in\n");
  /*
    icc();
    */
  break;

  default:
    (void) fprintf(stderr, "%sERROR: invalid preconditioning flag.\n"
                   "       options[AZ_precond] improperly set.\n", yo);
  exit(-1);

  }
  options[AZ_pre_calc] = AZ_sys_reuse;

} /* precond */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_calc_blk_diag_inv(double *val, int *indx, int *bindx, int *rpntr,
                          int *cpntr, int *bpntr, double *d_inv, int *d_indx,
                          int *d_bindx, int *d_rpntr, int *d_bpntr,
                          int *data_org)

/*******************************************************************************

  Routine to calculate the inverse of the block-diagonal portion of the sparse
  matrix in 'val' and the associated integer pointer vectors. This is used for
  scaling and/or preconditioning.

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

  d_inv:           Vector containing the inverses of the diagonal blocks.

  d_indx:          The 'indx' array corresponding to the inverse-block
                   diagonals.

  d_bindx:         The 'bindx' array corresponding to the inverse-block
                   diagonals.

  d_rpntr:         The 'rpntr' array corresponding to the inverse-block
                   diagonals.

  d_bpntr:         The 'bpntr' array corresponding to the inverse-block
                   diagonals.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see file params.txt).

*******************************************************************************/

{

  /* local variables */

  register int i, j, iblk_row, jblk, icount = 0, iblk_count = 0, ival;
  int          m1, n1, itemp;
  int          m;
  int          bpoff, idoff;
  int         *ipiv, info;
  double      *work;
  char        *yo = "AZ_calc_blk_diag_inv: ";

  /**************************** execution begins ******************************/

  m = data_org[AZ_N_int_blk] + data_org[AZ_N_bord_blk];

  if (m == 0) return;

  /* allocate vectors for lapack routines */

  ipiv = (int *)    malloc(rpntr[m]*sizeof(int));
  work = (double *) malloc(rpntr[m]*sizeof(double));

  /* offset of the first block */

  bpoff = *bpntr;
  idoff = *indx;

  /* loop over block rows */

  for (iblk_row = 0; iblk_row < m; iblk_row++) {

    /* number of rows in the current row block */

    m1 = rpntr[iblk_row+1] - rpntr[iblk_row];

    /* starting index of current row block */

    ival = indx[bpntr[iblk_row] - bpoff] - idoff;

    /* loop over column block numbers, looking for the diagonal block */

    for (j = bpntr[iblk_row] - bpoff; j < bpntr[iblk_row+1] - bpoff; j++) {
      jblk = bindx[j];

      /* determine the number of columns in this block */

      n1 = cpntr[jblk+1] - cpntr[jblk];

      itemp = m1*n1;

      if (jblk == iblk_row) {   /* diagonal block */

        /* error check */

        if (n1 != m1) {
          (void) fprintf(stderr, "%sERROR: diagonal blocks are not square\n.",
                         yo);
          exit(-1);
        }
        else {

          /* fill the vectors */

          d_indx[iblk_count]  = icount;
          d_rpntr[iblk_count] = rpntr[iblk_row];
          d_bpntr[iblk_count] = iblk_row;
          d_bindx[iblk_count] = iblk_row;

          for (i = 0; i < itemp; i++) d_inv[icount++] = val[ival + i];

          /* invert the dense matrix */

          dgetrf_(&m1, &m1, &d_inv[d_indx[iblk_count]], &m1, ipiv, &info);

          if (info < 0) {
            (void) fprintf(stderr, "%sERROR: argument %d is illegal.\n", yo,
                           -info);
            exit(-1);
          }

          else if (info > 0) {
            (void) fprintf(stderr, "%sERROR: the factorization has produced a "
                           "singular U with U[%d][%d] being exactly zero.\n",
                           yo, info, info);
            exit(-1);
          }

          dgetri_(&m1, &d_inv[d_indx[iblk_count]], &m1, ipiv, work, &m1, &info);

          if (info < 0) {
            (void) fprintf(stderr, "%sERROR: argument %d is illegal.\n", yo,
                           -info);
            exit(-1);
          }

          else if (info > 0) {
            (void) fprintf(stderr, "%sERROR: U[%d][%d] is exactly zero;\n", yo,
                           info, info);
            (void) fprintf(stderr, "the matrix is singular and its inverse "
                           "could not be computed.\n");
            exit(-1);
          }
          iblk_count++;
        }
        break;
      }
      else
        ival += itemp;
    }
  }

  d_indx[iblk_count]  = icount;
  d_rpntr[iblk_count] = rpntr[iblk_row];
  d_bpntr[iblk_count] = iblk_row;

  /* free vectors */

  free((void *) ipiv);
  free((void *) work);

} /* AZ_calc_blk_diag_inv */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void jacobi(double val[], double x[], int data_org[])

/*******************************************************************************

  Simple Jacobi iteration (undamped).  Not yet implemented for DVBR formatted
  matrices.

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

  x:               On input, contains the current solution to the linear system.
                   On output contains the Jacobi preconditioned solution.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see file params.txt).

*******************************************************************************/

{

  /* local variables */

  register int i;
  int          N;

  /**************************** execution begins ******************************/

  N = data_org[AZ_N_internal] + data_org[AZ_N_border];

  for (i = 0; i < N; i++)
    *x++ /=  *val++;

} /* jacobi */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

extern void AZ_sym_gauss_seidel(void)

/*******************************************************************************

  Symmetric Gauss-Siedel preconditioner

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:  void
  ===============

*******************************************************************************/

{

  /* local variables */

  /**************************** execution begins ******************************/

  (void) fprintf(stderr, "WARNING: sym Gauss-Seidel preconditioning not\n"
                 "         implemented for VBR matrices\n");
  exit(-1);

} /* AZ_sym_gauss_seidel */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_sym_gauss_seidel_sl(double val[], int bindx[], double x[],
                            int data_org[], int options[])

/*******************************************************************************

  Symmetric Gauss-Siedel preconditioner.

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

  x:               On input, contains the current solution to the linear system.
                   On output contains the Jacobi preconditioned solution.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see file params.txt).

  options:         Determines specific solution method and other parameters.

*******************************************************************************/

{

  /* local variables */

  register int    *bindx_ptr;
  register double sum, *ptr_val;
  int             i, bindx_row, j_last, N, step, ione = 1, j;
  double          *b, *ptr_b;

  /**************************** execution begins ******************************/

  N = data_org[AZ_N_internal] + data_org[AZ_N_border];

  b = AZ_manage_memory(N*sizeof(double), AZ_ALLOC, AZ_SYS, "b in sym GS", &i);

  dcopy_(&N, x, &ione, b, &ione);
  ptr_val = val;

  for (i = 0; i < N; i++) {
    (*ptr_val) = 1.0 / (*ptr_val);
    x[i]     = 0.0;
    ptr_val++;
  }

  for (step = 0; step < options[AZ_poly_ord]; step++) {
    AZ_exchange_bdry(x, data_org);

    bindx_row = bindx[0];
    bindx_ptr = &bindx[bindx_row];
    ptr_val   = &val[bindx_row];
    ptr_b   = b;

    for (i = 0; i < N; i++) {
      sum    = *ptr_b++;
      j_last = bindx[i+1] - bindx[i];

      for (j = 0; j < j_last; j++) {
        sum -= *ptr_val++ * x[*bindx_ptr++];
      }
      x[i] = sum * val[i];
    }

    bindx_row = bindx[N];
    bindx_ptr = &bindx[bindx_row-1];
    ptr_val   = &val[bindx_row-1];

    for (i = N - 1; i >= 0; i--) {
      sum = b[i];
      j_last  = bindx[i+1] - bindx[i];

      for (j = 0; j < j_last; j++) {
        sum -= *ptr_val-- * x[*bindx_ptr--];
      }
      x[i] = sum * val[i];
    }
  }

  for (i = 0; i < N; i++)
    val[i] = 1.0 / val[i];

} /* AZ_sym_gauss_seidel_sl */

#ifdef eigen
void AZ_do_Jacobi(double val[], int indx[], int bindx[], int rpntr[],
                     int cpntr[], int bpntr[], double x[], double b[],
                     double temp[], int options[], int data_org[],
                     int proc_config[], double params[], int flag)
{
double *v;
int i,step;
int N;
 
  N    = data_org[AZ_N_internal] + data_org[AZ_N_border];
 
    if (data_org[AZ_matrix_type] == AZ_MSR_MATRIX) {
 
      if ( (options[AZ_poly_ord] != 0) && (flag == 1) )
         for (i = data_org[AZ_N_internal]; i < N; i++) x[i] = b[i]/val[i];
 
      if (options[AZ_poly_ord] > flag) {
        v = AZ_manage_memory((N+data_org[AZ_N_external])*sizeof(double),
                             AZ_ALLOC, AZ_SYS, "v in do_jacobi", &i);
 
        for (i = 0; i < N; i++) v[i] = x[i];
 
        for (step = flag; step < options[AZ_poly_ord]; step++) {
          AZ_matvec_mult(val, indx, bindx, rpntr, cpntr, bpntr, v, temp, 1,
                         data_org);
          for(i = 0; i < N; i++) v[i] += (b[i] - temp[i]) / val[i];
        }
        for (i = 0; i < N; i++) x[i] = v[i];
      }
    }
    else {
       (void) fprintf(stderr,"AZ_slu with option[AZ_poly_ord] > 0 only \n");
       (void) fprintf(stderr,"implemented for MSR matrices.\n");
       exit(-1);
    }

}
#endif
