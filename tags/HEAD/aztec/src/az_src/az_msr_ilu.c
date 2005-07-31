/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile: az_msr_ilu.c,v $
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
static char rcsid[] = "$Id: az_msr_ilu.c,v 1.1.1.1 2001-01-30 20:59:12 paklein Exp $";
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

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void AZ_ilu_routine(int N, int M, double x[], int length, double buffer[],
                    double val[], int indx[], int bindx[], int rpntr[],
                    int cpntr[], int bpntr[], int options[], int data_org[])

/*******************************************************************************

  Routine preconditions the vector 'x' using an incomplete factorization.  If
  the factorization has not already been computed, this routine calculates the
  incomplete sparse lower and upper factors.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  N:               Order of linear system to be factorized (including any
                   external rows).

  M:               Number of nonzeros in sparse matrix 'a'.

  x:               On output, x[] is preconditioned by performing two backsolves
                   corresponding to incomplete ilu.

  length:          Length of buffer for exchanged rows.

  buffer:          Storage buffer for exchanged rows (i.e. external rows).

  val:             Array containing the nonzero entries of the matrix (see file
                   params.txt).

  indx,
  bindx,
  rpntr,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file params.txt).

  options:         Determines specific solution method and other parameters.  In
                   this routine, we are concerned with the options[AZ_overlap]:
                   options[AZ_overlap] == AZ_none
                   - nonoverlapping domain decomposition (don't use external
                     values).
                                       == AZ_diag
                   - use external values in conjunction with the external row,
                     however only keep the diagonal of the external row and set
                     up a Dirichlet-type boundary condition.
                                        == AZ_full
                   - duplicate external rows on this processor (only keeping the
                     coupling between variables contained within this
                     processor).

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see file params.txt).

*******************************************************************************/

{

  /* local variables */

  static double *msr_val;      /* matrix factors are stored in these arrays   */
  static int    *msr_bindx;    /* in MSR format.                              */
  static int    *extra;        /* extra[i] points to the start of the upper   */

  /* triangular elements in row i */

  static int previous_factors = -1;

  int    i, j;
  int    first2, last2;
  double sum;
  int    blks, oldN;
  int    st;

  /**************************** execution begins ******************************/

  oldN = data_org[AZ_N_internal] + data_org[AZ_N_border];
  blks = data_org[AZ_N_int_blk]  + data_org[AZ_N_bord_blk];

  if (oldN == 0) return;

  if (options[AZ_pre_calc] < AZ_sys_reuse) {

    /* allocate memory for new factors */

    extra = (int *) AZ_manage_memory((N+2)*sizeof(int), AZ_ALLOC,
                                     data_org[AZ_name], "extra in ilu", &st);

    /*
     * Use the last element of extra to store the size of M for future calls.
     */

    if (st == AZ_NEW_ADDRESS) extra[N+1] = M;
    else                      M          = extra[N+1];

    msr_bindx = (int *) AZ_manage_memory((M+3)*sizeof(int), AZ_ALLOC,
                                         data_org[AZ_name], "msr_bindx in ilu",
                                         &st);
    msr_val   = (double *) AZ_manage_memory((M+3)*sizeof(double), AZ_ALLOC,
                                            data_org[AZ_name], "msr_val in ilu",
                                            &st);
    if (msr_val == NULL) {
      (void) fprintf(stderr, "Error: not enough memory in ilu\n"
                     "       Run a smaller problem or another\n"
                     "       preconditioner.\n");
      exit(-1);
    }

    if (options[AZ_pre_calc] == AZ_reuse) {
      if (st == AZ_NEW_ADDRESS) {
        (void) fprintf(stderr, "Error: options[AZ_pre_calc] == AZ_reuse and "
                       "previous factors\n       not found. Check "
                       "data_org[AZ_name].\n");
        exit(-1);
      }
    }
  }

  else if (previous_factors != data_org[AZ_name] ) {
    (void) fprintf(stderr, "Warning: Using a previous factorization as a \n"
                   "preconditioner \neven though matrix"
                   "(data_org[AZ_name]) has changed\n");
  }

  previous_factors = data_org[AZ_name];

  if (options[AZ_pre_calc] <= AZ_recalc) {

    /* calculate the factors */

    for (i = 0; i < N+1; i++) extra[i] = 0;

    /*
     * Compute the l u factors corresponding to an incomplete factorization
     * first set up the msr matrix for this domain (including external rows).
     */

    if (data_org[AZ_matrix_type] == AZ_VBR_MATRIX) {
      AZ_vb2msr(blks, val, indx, bindx, rpntr, cpntr, bpntr, msr_val,
                msr_bindx);
      val   = msr_val;
      bindx = msr_bindx;
    }

    AZ_msr2ilu_setup(val, bindx, msr_val, msr_bindx, N, oldN, options,
                     length, buffer, data_org);

    /* Reorganize the matrix A such that
     *  msr_val[i]   = A(i,i)                               i < N
     *  msr_val[i]   = lower triangular nonzeros in row j
     *                                          msr_bindx[j]   <= i < extra[j]
     *  msr_bindx[i] = column index for msr_val[i]
     *                                          msr_bindx[j]   <= i < extra[j]
     *  msr_val[i]   = upper triangular nonzeros in col j
     *                                          extra[j] <= i < msr_bindx[j+1]
     *  msr_bindx[i] = row index for msr_val[i]
     *                                          extra[j] <= i < msr_bindx[j+1]
     */

    AZ_reorganize_matrix(msr_bindx, msr_val, N, extra);

    /* do factorization */

    for (i = 0; i < N; i++) {

      /* first do l for this row */

      for (j = msr_bindx[i]; j < extra[i]; j++) {

        /*
         * Compute the sum for l_i,msr_bindx[j]. That means compare row i with
         * column msr_bindx[j].
         */

        first2 = extra[msr_bindx[j]];
        last2  = msr_bindx[msr_bindx[j]+1] - 1;

        sum  = AZ_sum_comp(msr_bindx[i], j-1, first2, last2, msr_val,
                           msr_bindx);
        msr_val[j] = (msr_val[j] - sum) / msr_val[msr_bindx[j]];
      }

      /* compute the offdiag of u for this column */

      for (j = extra[i]; j < msr_bindx[i+1]; j++) {

        /*
         * Compute the sum for u_msr_bindx[j], i.  That means compare row
         * msr_bindx[j] with column i.
         */

        first2 = msr_bindx[msr_bindx[j]];
        last2  = extra[msr_bindx[j]] - 1;

        sum  = AZ_sum_comp(extra[i], j-1, first2, last2, msr_val, msr_bindx);
        msr_val[j] = (msr_val[j] - sum);
      }

      /* compute the u_ii term */

      sum  = AZ_sum_comp(msr_bindx[i], extra[i]-1, extra[i], msr_bindx[i+1]-1,
                         msr_val, msr_bindx);
      msr_val[i] -= sum;
    }
  }

  /*
   * Using lower triangular factor L (stored in msr_val,msr_bindx,extra) compute
   * L u = x where the solution u is overwritten in x.
   */

  AZ_lower(x, msr_val, msr_bindx, extra, N);

  /*
   * Using upper triangular factor U (stored in msr_val,msr_bindx,extra) compute
   * U u = x where the solution u is overwritten in x.
   */

  AZ_upper(x, msr_val, msr_bindx, extra, N);

} /* ilu */

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

double AZ_sum_comp(int first1, int last1, int first2, int last2, double val[],
                   int bindx[])

/*******************************************************************************

  Compute the inner product between the two sparse vectors:

     val[first1 ... last1]   and   val[first2 ... last2]

  where val[first ... last] represents a sparse vector 'v' defined as all zeros
  except for the following elements:

           v[bindx[k]] = val[k]    k = first to last

  IMPORTANT NOTE: it is assumed that

       bindx[i] < bindx[k]

  for  first <= i < k <= last

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     double, result of inner product.
  ============

  Parameter list:
  ===============

  first1:          'val[first1]' is first element of first vector.

  last1:           'val[last1]' is last element of first vector.

  first2:          'val[first1]' is first element of second vector.

  last2:           'val[last1]' is last element of second vector.

  val:              Holds nonzeros of vectors.

  bindx:           'bindx[k]' gives index of nonzero element 'val[k]'.

*******************************************************************************/

{

  /* local variables */

  int    flag;
  double sum = 0.0;

  /**************************** execution begins ******************************/

  if ((first1 > last1) || (first2 > last2)) flag = 0;
  else                                      flag = 1;

  while (flag == 1) {
    if (bindx[first2] < bindx[first1]) {
      first2++;

      if (first2 > last2) flag = 0;
    }
    else if (bindx[first2] > bindx[first1]) {
      first1++;

      if (first1 > last1) flag = 0;
    }
    else {
      sum += (val[first2]*val[first1]);
      first1++; first2++;

      if ((first2 > last2) || (first1 > last1)) flag = 0;
    }
  }

  return sum;

} /* AZ_sum_comp */

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void AZ_lower(double rhs[], double val[], int bindx[], int extra[], int n)

/*******************************************************************************

  Compute the solution to L u = rhs where

      u is written back over rhs[].
      L is the lower triangular matrix stored in val[], bindx[], extra[] and

         L(i,i) = 1.0
         L(j,i) = val[i]  where  i < j  and  bindx[j] <= i < extra[j]

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  rhs:             On input, right-hand side of system to be solved.
                   On output, solution of system.

  val:             Nonzeros of lower triangular matrix (see above).

  bindx:           Lower triangular matrix (see above comments).

  extra:           Lower triangular matrix (see above comments).

  n:               Order of linear system to be factorized (including any
                   external rows).

*******************************************************************************/

{

  /* local variables */

  int    k, j;
  double sum;

  /**************************** execution begins ******************************/

  for (k = 0; k < n; k++) {

    /* compute the sum */

    sum = 0.0;

    for (j = bindx[k]; j < extra[k]; j++) {
      if (bindx[j] > k )
        break;
      sum += val[j] * rhs[bindx[j]];
    }
    rhs[k] = (rhs[k] - sum);
  }

} /* lower */

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void AZ_upper(double rhs[], double val[], int bindx[], int extra[], int n)

/*******************************************************************************

  Compute the solution to U u = rhs where

      u is written back over rhs[].
      U is the upper triangular matrix stored in val[], bindx[], extra[] and

         U(i,i) = 1.0
         U(j,i) = val[i]  where  i < j  and  bindx[j] <= i < extra[j]

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  rhs:             On input, right-hand side of system to be solved.
                   On output, solution of system.

  val:             Nonzeros of lower triangular matrix (see above).

  bindx:           Lower triangular matrix (see above comments).

  extra:           Lower triangular matrix (see above comments).

  n:               Order of linear system to be factorized (including any
                   external rows).

*******************************************************************************/

{

  /* local variables */

  int    k, j;

  /**************************** execution begins ******************************/

  for (k = n-1; k >= 0; k--) {

    rhs[k] = rhs[k] / val[k];

    for (j = extra[k]; j < bindx[k+1]; j++) {
      rhs[bindx[j]] -= val[j] * rhs[k];
    }
  }

} /* upper */

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void AZ_reorganize_matrix(int bindx[], double val[], int N, int extra[])

/*******************************************************************************

  Reorganize the MSR matrix A such that

   val[i] = A(i,i)                                 i < N
   val[i] = lower triangular nonzeros for row j    bindx[j]   <= i < extra[j]
   val[i] = upper triangular nonzeros for column j extra[j]   <= i < bindx[j+1]

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  val, bindx:      On input, MSR matrix.
                   On output, reorganized as described above.

  n:               Order of linear system to be factorized (including any
                   external rows).

  extra:           Lower triangular matrix (see above comments).

*******************************************************************************/

{

  /* local variables */

  int    *extra2;
  int     i, j, k;
  int     offset1, offset2, buf_len, buf_head = 0, *buf_l;
  int    *buf_i;
  double *buf_a;
  int     location, ival, start;
  double  aval;
  int     first, last;
  int     temp_loc;

  /**************************** execution begins ******************************/

  extra2 = (int *) calloc(N + 1, sizeof(int));
  if (extra2 == NULL) {
    (void) fprintf(stderr, "Error: not enough memory for ilu\n"
                   "       Run a smaller problem or another\n"
                   "       preconditioner.\n");
    exit(-1);
  }

  /* Allocate a temporary buffer (as large as we can get away with) */

  buf_a   = NULL; buf_i = NULL; buf_l = NULL;
  buf_len = 2 * bindx[N];

  while (buf_a == NULL) {
    if (buf_i != NULL) free(buf_i);
    if (buf_l != NULL) free(buf_l);

    buf_len /= 2;

    if (buf_len == 0) {
      (void) fprintf(stderr, "Error: not enough memory for ilu\n"
                     "       Run a smaller problem or another\n"
                     "       preconditioner.\n");
      exit(-1);
    }

    buf_i = (int    *) calloc(buf_len, sizeof(int));
    buf_l = (int    *) calloc(buf_len, sizeof(int));
    buf_a = (double *) calloc(buf_len, sizeof(double));
  }

  start = bindx[0];

  /* calculate the number of nonzeros in the l factor for each row */
  /* and the number of nonzeros in the u factor for each column.   */

  for(i = 0; i < N; i++) {
    first = bindx[i];
    last  = bindx[i+1] - 1;
    AZ_sort( &(bindx[first]), last - first + 1, NULL, &(val[first]));

    for (j = first; j <= last; j++) {
      if (bindx[j] > i) extra[bindx[j]]++;
      else              extra2[i]++;
    }
  }

  /* computing starting location for each row of l and each column of u */

  first = bindx[0];
  for (i = 0; i < N; i++) {
    offset1   = extra2[i];
    extra2[i] = first;
    offset2   = extra[i];
    extra[i]  = extra2[i] + offset1;
    first     = extra[i] + offset2;
  }

  /* Reorganize the matrix A such that
   *  bindx[i] = A(i,i)                               i < N
   *  bindx[i] = lower triangular nonzeros of row j  bindx[j]   <= i <extra[j]
   *  val[i]   = upper triangular nonzeros of column j extra[j] <= i <bindx[j+1]
   */

  i = 0;
  for (j = bindx[0]; j < bindx[N]; j++) {
    while (bindx[i+1] <= j) i++;

    if (bindx[j] < i) {
      location = extra2[i];
      ival     = bindx[j];
      extra2[i]++;
    }
    else {
      location = extra[bindx[j]];
      ival     = i;
      extra[bindx[j]]++;
    }

    aval = val[j];

    /* try and store the new entry */

    if (location <= j) {
      val[location]   = aval;
      bindx[location] = ival;
    }
    else {

      /* store in the buffer if there is room */

      if (buf_head < buf_len) {
        buf_a[buf_head] = aval;
        buf_i[buf_head] = ival;
        buf_l[buf_head] = location;
        buf_head++;
      }
      else {

        /* empty what we can from the buffer */

        for (k = 0; k < buf_head; k++) {
          temp_loc = buf_l[k];

          if (temp_loc <= j ) {
            val[temp_loc]   = buf_a[k];
            bindx[temp_loc] = buf_i[k];
            buf_l[k]      = -1;
          }
        }

        /* compress the buffer */

        buf_head = 0;
        for (k = 0; k < buf_len; k++) {
          if (buf_l[k] != -1 ) {
            buf_l[buf_head] = buf_l[k];
            buf_a[buf_head] = buf_a[k];
            buf_i[buf_head] = buf_i[k];
            buf_head++;
          }
        }

        /* now try to store it in the buffer */

        if (buf_head < buf_len) {
          buf_a[buf_head] = aval;
          buf_i[buf_head] = ival;
          buf_l[buf_head] = location;
          buf_head++;
        }
        else {
          (void) fprintf(stderr,
                         "Error: buffer is full can not continue with ilu\n");
          (void) fprintf(stderr, "Error: not enough memory for ilu\n"
                         "       Run a smaller problem or\n"
                         "       another preconditioner.\n");
          exit(-1);
        }
      }
    }
  }

  /* empty the buffer */

  for (k = 0; k < buf_head; k++) {
    location      = buf_l[k];
    val[location]   = buf_a[k];
    bindx[location] = buf_i[k];
  }

  free(buf_a); free(buf_i); free(buf_l);

  for (i = N-1; i >= 0; i--) {
    bindx[i+1] = extra[i];
    extra[i] = extra2[i];
  }

  bindx[0] = start;

  free(extra2);

} /* AZ_reorganize_matrix */
