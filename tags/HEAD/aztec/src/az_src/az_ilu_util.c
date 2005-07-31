/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile: az_ilu_util.c,v $
 *
 * $Author: paklein $
 *
 * $Date: 2001-01-30 20:59:13 $
 *
 * $Revision: 1.1.1.1 $
 *
 * $Name: not supported by cvs2svn $
 *====================================================================*/
#ifndef lint
static char rcsid[] = "$Id: az_ilu_util.c,v 1.1.1.1 2001-01-30 20:59:13 paklein Exp $";
#endif


/*******************************************************************************
 * Copyright 1995, Sandia Corporation.  The United States Government retains a *
 * nonexclusive license in this software as prescribed in AL 88-1 and AL 91-7. *
 * Export of this program may require a license from the United States         *
 * Government.                                                                 *
 ******************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "az_aztec.h"

/* missing prototypes - PAK (07/24/98) */
void sort_blk_col_indx(int num_blks_row, int *bindx_start_row,
                       int *ordered_index);
void sort2(int n, int ra[], int rb[]);
void order_parallel(int M, double *val_old, double *val_new, int *bindx_old,
                    int *bindx_new, int *indx_old, int *indx_new,
                    int *bpntr_old, int *bpntr_new, int *diag_block);
void get_diag(int M, int *bindx, int *bpntr, int *diag_block);

/******************************************************************************/

void sort_blk_col_indx(int num_blks_row, int *bindx_start_row,
                       int *ordered_index)

/*******************************************************************************

  Routine to sort the block entires for a given block row into increasing order.

  Author:
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  num_blks_row:    Number of blocks in in this row.

  bindx_start_row: On Input:  starting address for the array of block column
                              indices for this row.
                   On Output: ordered block column indicies for this row.
                              (size >= num_blks_row).

  ordered_index:   On Input:  integer array.
                   On Output: an array of indices which describes the
                              reordering for the array "bindx_start_row"
                              (size >= num_blks_row).

*******************************************************************************/

{

  /* local variables */

  int i;

  /* externals */

  void sort2(int, int *, int*);

  /**************************** execution begins ******************************/

  /* Initialize the ordering index vector */

  for (i = 0; i < num_blks_row; i++) ordered_index[i] = i;

  /* Sort block column index array and produce the ordering index */

  sort2(num_blks_row, bindx_start_row-1, ordered_index-1);

} /* sort_blk_col_indx */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void sort2(int n, int ra[], int rb[])

/*******************************************************************************

  Numerical Recipes C source code modified to have first argument an integer
  array.

  Sorts the array ra[1,..,n] in ascending numerical order using heapsort
  algorithm, while making the corresponding rearrangement of the array
  rb[1,..,n].

  NOTE: The arrays start at 1 instead of 0, therefore you must pass call from C
  for a zero based array as:

                   sort(n, ra-1, rb-1);


  Author:          Modified by John N. Shadid, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  n:               Length of arrays ra and rb.

  ra:              Array to be sorted.

  rb               Second array order acording to the sorted ra array.

*******************************************************************************/

{

  /* local variables */

  int l, j, ir, i;
  int rra;
  int rrb;

  /**************************** execution begins ******************************/

  l  = (n >> 1) + 1;
  ir = n;

  for (;;) {
    if (l > 1) {
      rra = ra[--l];
      rrb = rb[l];
    }
    else {
      rra    = ra[ir];
      rrb    = rb[ir];
      ra[ir] = ra[1];
      rb[ir] = rb[1];

      if (--ir <= 1) {
        ra[1] = rra;
        rb[1] = rrb;
        return;
      }
    }

    i = l;
    j = l << 1;

    while (j <= ir) {
      if (j < ir && ra[j] < ra[j+1]) ++j;
      if (rra < ra[j]) {
        ra[i] = ra[j];
        rb[i] = rb[j];
        j    += (i = j);
      }
      else j = ir + 1;
    }

    ra[i] = rra;
    rb[i] = rrb;
  }

} /* sort2 */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_order(int M, double *val_old, double *val_new, int *bindx,
              int *indx_old, int *indx_new, int *bpntr, int *diag_block)

/*******************************************************************************

  For each row, reorders the blocks of the matrix (indices and values) in
  increasing order and constructs the array of pointers: diag_block to the
  diagonal blocks.

  Returns diag_block[i] = -1 if no diagonal block has been found at the i'th
  row.

  Author:          Lydie Prevost, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  M:               The number of (block) rows in the matrix.

  val_old:         Array containing the entries of the matrix before reordering.
                   The matrix is stored block-row-by-block-row. Each block entry
                   is dense and stored by columns (VBR).

  val_new:         On output, array containing the entries of the matrix after
                   reordering. The matrix is stored block-row-by-block-row. Each
                   block entry is dense and stored by columns (VBR).

  bindx:           On input, contains the block column indices of the non-zero
                   block entries of the matrix before reordering.
                   On output, contains the block column indices of the non-zero
                   block entries of the matrix after reordering.

  indx_old:        The ith element of indx_old points to the location in val_old
                   of the (0,0) entry of the ith block entry. The last element
                   is the number of nonzero entries of matrix plus one.

  indx_new:        The ith element of indx_new points to the location in val_new
                   of the (0,0) entry of the ith block entry. The last element
                   is the number of nonzero entries of matrix plus one.

  bpntr:           The ith element of bpntr points to the first block entry of
                   the ith row in bindx. The last element is the number of
                   nonzero blocks of matrix plus one.

  diag_block:      On output, array of size M points on each diagonal block.

*******************************************************************************/

{

  /* local variables */

  int     i, kk, j, ii;
  int     num_blks_row;
  int    *sort, old_blk_index, counter;
  int     new_blk;
  int    *temp_ind, size_temp_ind = 10, size_temp_val = 40;
  double *temp_val;
  int     total_vals;
  int     start, end;

  /**************************** execution begins ******************************/

  /* malloc array */

  temp_ind = (int *)    calloc(size_temp_ind, sizeof(int));
  temp_val = (double *) calloc(size_temp_val, sizeof(double));

  sort = (int *) malloc(sizeof(int) * (M));

  for (i = 0; i < M; i++) diag_block[i] = -1;

  for (i = 0; i < M; i++) {     /* loop over the rows */

    /* constructs the new array of block column indices in this row */

    num_blks_row = bpntr[i+1] - bpntr[i];

    if (num_blks_row+1 > size_temp_ind) {
      size_temp_ind = num_blks_row + 1;
      free(temp_ind);
      temp_ind = (int *) calloc(size_temp_ind, sizeof(int));
    }

    for (ii = bpntr[i]; ii <= bpntr[i+1]; ii++)
      temp_ind[ii - bpntr[i]] = indx_old[ii];

    total_vals = indx_old[bpntr[i+1]] - indx_old[bpntr[i]];

    sort_blk_col_indx(num_blks_row, bindx+bpntr[i], sort);

    /* for each block of this row computes the new indices and constructs the
       pointers on the diagonal block i */

    indx_new[0] = indx_old[0];
    for (kk = 0; kk < num_blks_row; kk++) {
      new_blk = kk + bpntr[i];

      /* index into old VBR matrix and get the block size */

      indx_new[new_blk+1] = indx_new[new_blk] +
        (temp_ind[sort[kk]+1] - temp_ind[sort[kk]]);

      if (bindx[new_blk] == i) diag_block[i] = new_blk;
    }

    /* constructs the new array containing the entries of the matrix */

    if (total_vals > size_temp_val) {
      size_temp_val = total_vals;
      free(temp_val);
      temp_val = (double *) calloc(size_temp_val, sizeof(double));
    }

    start   = indx_old[bpntr[i]];
    end     = indx_old[bpntr[i+1]];
    counter = 0;

    for (ii = start; ii < end; ii++)
      temp_val[counter++] = val_old[ii];

    for (kk = 0; kk < num_blks_row; kk++) {

      /*
       * Get old block index into the coefficient array and new block location.
       */

      old_blk_index = temp_ind[sort[kk]] - temp_ind[0];
      new_blk       = kk + bpntr[i];

      counter = 0;
      for (j = indx_new[new_blk]; j < indx_new[new_blk+1]; j++) {
        val_new[j] = temp_val[old_blk_index + counter++];
      }
    }
  }

  free((void *) sort);
  free((void *) temp_ind);
  free((void *) temp_val);

} /* order */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void order_parallel(int M, double *val_old, double *val_new, int *bindx_old,
                    int *bindx_new, int *indx_old, int *indx_new,
                    int *bpntr_old, int *bpntr_new, int *diag_block)

/*******************************************************************************

  For each row, reorders the blocks of the matrix (indices and values) in
  increasing order and constructs the array of pointers: diag_block to the
  diagonal blocks.

  Returns diag_block[i]=-1  if no diagonal block has been found at the i'th row.

  Author:          Lydie Prevost, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  M:               The number of (block) rows in the matrix.

  val_old:         Array containing the entries of the matrix before reordering.
                   The matrix is stored block-row-by-block-row. Each block entry
                   is dense and stored by columns (VBR).

  val_new:         On output, array containing the entries of the matrix after
                   reordering. The matrix is stored block-row-by-block-row. Each
                   block entry is dense and stored by columns (VBR).

  bindx_old:       Contains the block column indices of the non-zero block
                   entries of the matrix before reordering.

  bindx_new:       On output, contains the block column indices of the reordered
                   non-zero block entries of the matrix before reordering.

  indx_old:        The ith element of indx_old points to the location in val_old
                   of the (0,0) entry of the ith block entry. The last element
                   is the number of nonzero entries of matrix plus one.

  indx_new:        The ith element of indx_new points to the location in val_new
                   of the (0,0) entry of the ith block entry. The last element
                   is the number of nonzero entries of matrix plus one.

  bpntr_old:       The i'th element of bpntr points to the first block entry of
                   the i'th row in bindx_old. The last element is the number of
                   nonzero blocks of matrix plus one.

  bpntr_new:       On, output, the i'th element of bpntr points to the first
                   block entry of the i'th row in bindx_new. The last element
                   is the number of nonzero blocks of matrix plus one.

  diag_block:      On output, array of size M points on each diagonal block.

*******************************************************************************/

{

  /* local variables */

  int  i, kk, j;
  int  num_blks_row_old, num_blks_row_new;
  int *sort, ptr, compt, temp;

  /**************************** execution begins ******************************/

  /* Allocate work space */

  sort = (int *) malloc(sizeof(int) * (M));
  if (sort == NULL) {
    (void) fprintf(stderr, "Error: not enough memory inside order_parallel\n"
                   "       must run a smaller problem\n");
    exit(-1);
  }

  /* Initialize */

  for (i = 0; i < M; i++) diag_block[i] = -1;
  bpntr_new[0] = bindx_new[0] = 0;

  for (i = 0; i < M; i++) {     /* loop over the rows */

    /* constructs the new array of block column indices in this row */

    num_blks_row_old = bpntr_old[i+1] - bpntr_old[i];

    /* copy old block column index array and then sort it this defines the new
       block column index array */

    for (j = 0; j < num_blks_row_old; j++)
      bindx_new[bpntr_new[i] + j] = bindx_old[bpntr_old[i] + j];

    sort_blk_col_indx(num_blks_row_old, &bindx_new[bpntr_new[i]], sort);

    /* Count the blocks that multiply internal and border unknowns */

    num_blks_row_new = 0;
    for (j = 0; j < num_blks_row_old; j++) {
      if (bindx_new[bpntr_new[i] + j] >= M) break;

      num_blks_row_new++;
    }

    bpntr_new[i+1] = bpntr_new[i] + num_blks_row_new;

    /* for each block of this row compute the new index vector and construct the
       pointers to the diagonal block i */

    for (kk = bpntr_new[i]; kk < bpntr_new[i+1]; kk++) {

      /* Define new indx vector */

      if (kk - bpntr_new[i] == 0) {
        indx_new[0] = indx_old[0];
      }
      else {
        temp         = sort[kk-1-bpntr_old[i]] + bpntr_old[i];
        indx_new[kk] = indx_new[kk-1] + (indx_old[temp+1] - indx_old[temp]);
      }

      /* Get diagonal block pointers */

      if (bindx_new[kk] == i) diag_block[i] = kk;
    }

    /* constructs the new array containing the entries of the matrix */

    for (kk = bpntr_new[i]; kk < bpntr_new[i+1]; kk++) {
      ptr   = indx_old[sort[kk-bpntr_old[i]] + bpntr_old[i]];
      compt = -1;

      for (j = indx_new[kk]; j < indx_new[kk+1]; j++) {
        compt++;
        val_new[j] = val_old[ptr + compt];
      }
    }
  }

  free((void *) sort);

} /* order_parallel */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_lower_triang_vbr_solve(int N, int M, double *val, int *indx, int *bindx,
                            int *rpntr, int *cpntr, int *bpntr, int *diag_block,
                            double *y, double *b)

/*******************************************************************************

  Lower triangular solver for Ly = b with L stored in VBR matrix format.

  Note: In this version the diagonal blocks of L are assumed to be identity.

  Author:          Lydie Prevost, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  N:               Leading dimension of b and y.

  M:               Number of (block) rows in the matrix and L.

  val:             Array containing the nonzero entries of the matrix (see file
                   params.txt).

  indx,
  bindx,
  rpntr,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file params.txt).

  diag_block:      On output, array of size M points on each diagonal block.

  b:               Right hand side of linear system.

  y:               On output, contains the result vector y.

*******************************************************************************/

{

  /* local variables */

  register int  iblk_row, j, jblk;
  int           iy, m1, ival, ib1, ib2, n1;
  int           ione = 1;
  double        one = 1.0, minus_one=-1.0;
  char         *T = "N";
  double       *y_pntr;

  /**************************** execution begins ******************************/

  /* initialize the result vector */

  for (j = 0; j < N; j++) y[j] = b[j];

  /* loop over block rows */

  for (iblk_row = 1; iblk_row < M; iblk_row++) {

    /* starting point row of the current row block */

    iy = rpntr[iblk_row];

    /* set result pointer */

    y_pntr = y + iy;

    /* number of rows in the current row block */

    m1 = rpntr[iblk_row+1] - rpntr[iblk_row];

    /* starting index of current row block */

    ival = indx[bpntr[iblk_row]];

    /* loop over all the lower blocks in the current block-row */

    for (j = bpntr[iblk_row]; j < diag_block[iblk_row]; j++) {
      jblk = bindx[j];

      /* the starting point column index of the current block */

      ib1 = cpntr[jblk];

      /* ending point column index of the current block */

      ib2 = cpntr[jblk+1];

      /* number of columns in the current block */

      n1 = ib2 - ib1;

      /* dense matrix-vector multiplication */

      if (m1 == 1 && n1 == 1)
        *y_pntr -= (*(val+ival++) * *(y+ib1));
      else {
        if (m1 < 10)
          AZ_dgemv3(m1, n1, val+ival, y+ib1, y_pntr);
        else
          dgemv_(T, &m1, &n1, &minus_one, val+ival, &m1, y+ib1, &ione, &one,
                 y_pntr, &ione, strlen(T));
        ival += (m1*n1);
      }
    }
  }

} /* AZ_lower_triang_vbr_solve */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_upper_triang_vbr_solve(int N, int M, double *val, int *indx, int *bindx,
                               int *rpntr, int *cpntr, int *bpntr,
                               int *diag_block, double *y, double *b,
                               int *ipvt)

/*******************************************************************************

  Upper triangular solver for Uy = b with U stored in VBR matrix format.

  Author:          Lydie Prevost, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  N:               Leading dimension of b and y.

  M:               Number of (block) rows in the matrix and L.

  val:             Array containing the nonzero entries of the matrix (see file
                   params.txt).

  indx,
  bindx,
  rpntr,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file params.txt).

  diag_block:      On output, array of size M points on each diagonal block.

  b:               Right hand side of linear system.

  y:               On output, contains the result vector y.

  Routines called:
  ================

  AZ_dgemv3 or dgemv (blas2 dense matrix vector multiply)

*******************************************************************************/

{

  /* local variables */

  register int  iblk_row, j, jblk;
  int           iy, m1, ival, ib1, ib2, n1;
  int           ione = 1;
  int           info;
  double        one = 1.0, minus_one = -1.0;
  char         *T = "N";
  double       *y_pntr;

  /**************************** execution begins ******************************/

  /* initialize  the result vector */

  for (j = 0; j < N; j++) y[j] = b[j];

  /* loop over block rows */

  for (iblk_row = M - 1; iblk_row >= 0; iblk_row--) {

    /* starting point row of the current row block */

    iy = rpntr[iblk_row];

    /* set result pointer */

    y_pntr  = y + iy;

    /* number of rows in the current row block */

    m1 = rpntr[iblk_row+1] - rpntr[iblk_row];

    /* starting index of current row block */

    ival = indx[bpntr[iblk_row]];

    /* loop over all the upperblocks in the current block-row */

    for (j = diag_block[iblk_row] + 1; j < bpntr[iblk_row+1]; j++) {
      ival = indx[j];
      jblk = bindx[j];

      /* the starting point column index of the current block */

      ib1 = cpntr[jblk];

      /* ending point column index of the current block */

      ib2 = cpntr[jblk+1];

      /* number of columns in the current block */

      n1 = ib2 - ib1;

      /* dense matrix-vector multiplication and division */

      if (m1 == 1 && n1 == 1)
        *y_pntr -=  (*(y+ib1) * *(val+ival++) );
      else {
        if (m1 < 10)
          AZ_dgemv3(m1, n1, val+ival, y+ib1, y_pntr);
        else
          dgemv_(T, &m1, &n1, &minus_one, val+ival, &m1, y+ib1, &ione,
                 &one, y_pntr, &ione, strlen(T));
        ival += m1*n1;
      }
    }

    j    = diag_block[iblk_row];
    jblk = bindx[j];
    ib1  = cpntr[jblk];
    ib2  = cpntr[jblk+1];
    n1   = ib2 - ib1;
    ival = indx[j];
    dgetrs_(T, &n1, &ione, val+ival, &m1, &(ipvt[iy]), y_pntr, &m1, &info,
            strlen(T));
  }

} /* AZ_upper_triang_vbr_solve */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void get_diag(int M, int *bindx, int *bpntr, int *diag_block)

/*******************************************************************************

  Constructs the array of pointers: diag_block to the diagonal blocks.  Returns
  diag_block[i] = -1 if no diagonal block has been found at the ith row.

  Author:          Lydie Prevost, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  M:               Number of (block) rows in the matrix and L.

  bindx:           Contains the block column indices of the non-zero block
                   entries.

  bpntr:           The i'th element of bpntr points to the first block entry of
                   the i'th row in bindx. The last element is the number of
                   nonzero blocks of matrix plus one.

  diag_block:      On output, array of size M points on each diagonal block.

*******************************************************************************/

{
  int i, kk;

  for (i = 0; i < M; i++) diag_block[i] = -1;

  for (i = 0; i < M; i++) {
    for (kk = bpntr[i]; kk < bpntr[i+1]; kk++) {
      if (bindx[kk] == i)
        diag_block[i] = kk;
    }
  }

} /* get_diag */
