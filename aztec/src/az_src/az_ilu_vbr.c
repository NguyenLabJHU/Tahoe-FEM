/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile: az_ilu_vbr.c,v $
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
static char rcsid[] = "$Id: az_ilu_vbr.c,v 1.1.1.1 2001-01-30 20:59:13 paklein Exp $";
#endif


/*******************************************************************************
 * Copyright 1995, Sandia Corporation.  The United States Government retains a *
 * nonexclusive license in this software as prescribed in AL 88-1 and AL 91-7. *
 * Export of this program may require a license from the United States         *
 * Government.                                                                 *
 ******************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifndef __MWERKS__
#include <malloc.h>
#endif /* __MWERKS__ - PAK (07/24/98) */
#include <string.h>
#include "az_aztec.h"

/* missing prototypes - PAK (07/24/98) */
void AZ_funswill(int *trash);

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_block_ilu(double val2[], int indx2[], int bindx2[],
                  int cpntr2[], int bpntr2[], int new_blks, int new_N,
                  int length, double buffer[], double x[], int options[],
                  int data_org[])

/*******************************************************************************

  Routine preconditions the vector 'x' using an incomplete block factorization.
  If the factorization has not already been computed, this routine calculates
  the incomplete sparse lower and upper factors.

  Author:          Lydie Prevost, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  val2:            Array containing the entries of the matrix. The matrix is
                   stored block-row-by-block-row. Each block entry is dense and
                   stored by columns (VBR).

  indx2,
  bindx2,
  cpntr2,
  bpntr2:          Arrays used for DMSR and DVBR sparse matrix storage (see
                   file params.txt).

  new_blks:        Number of blocks rows in expanded system corresponding to
                   domain decomposition.

  new_N:           Number of rows in expanded system corresponding to domain
                   decomposition.

  length:          Length of buffer.

  buffer:          Buffer holding extra rows received from neighbors.

  x:               On input, contains the current solution to the linear system.

  options:         Determines specific solution method and other parameters.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see file params.txt).

*******************************************************************************/

{

  /* local variables */

  int     I, i, j, k, K_SUM, flag_block, ptr_j, diag_j;
  double *z;
  double *block1, *block2;
  int     indxsize, valsize;
  int     flag, point_col, block_row, block_col;
  int    *therows, row_count;
  int    *theptrs, space, *tbindx, *newbindx;
  int     kk;
  int     start, incr;
  int     value;
  int     N, blks, blk_size;
  int     info, m1, ival, iy;

  /* these arrays hold the Block ILU factors */

  static int *rpntr, *cpntr, *bpntr, *bindx, *indx, *diag_block, *ipvt;

  /*  diag_block == int work array of size n points to each diagonal block */

  static double *val;

  double        *x_tmp;
  int            st;
  static int     previous_factors = -1;

  /* external functions */

  extern void AZ_funswill(int *);

  /**************************** execution begins ******************************/

  N    = data_org[AZ_N_internal] + data_org[AZ_N_border];
  blks = data_org[AZ_N_int_blk] + data_org[AZ_N_bord_blk];

  if (N == 0) return;

  if (options[AZ_pre_calc] == AZ_reuse) {

    /*
     * Previous factorization is being used. Find pointers corresponding to this
     * previous factorization.
     */

    bpntr = (int *) AZ_manage_memory((new_blks+1)*sizeof(int), AZ_ALLOC,
                                     data_org[AZ_name], "bpntr", &st);
    if (st == AZ_NEW_ADDRESS) {
      (void) fprintf(stderr, "Error: options[AZ_pre_calc] == AZ_reuse and "
                     "previous factors\n       not found. "
                     "Check data_org[AZ_name].\n");
      exit(-1);
    }

    indxsize = bpntr[new_blks]+1;
    indx  = (int *) AZ_manage_memory((indxsize+1)*sizeof(int), AZ_ALLOC,
                                     data_org[AZ_name], "indx", &st);
    bindx = (int *) AZ_manage_memory((indxsize+1)*sizeof(int), AZ_ALLOC,
                                     data_org[AZ_name], "bindx", &st);
    rpntr = (int *) AZ_manage_memory((new_blks+1)*sizeof(int), AZ_ALLOC,
                                     data_org[AZ_name], "rpntr", &st);
    cpntr = (int *) AZ_manage_memory((new_blks+1)*sizeof(int), AZ_ALLOC,
                                     data_org[AZ_name], "cpntr", &st);

    valsize    = indx[bpntr[new_blks]]+1;
    val        = (double *) AZ_manage_memory((valsize+1)*sizeof(double),
                                             AZ_ALLOC, data_org[AZ_name], "val",
                                             &st);
    diag_block = (int *) AZ_manage_memory((new_blks+1)*sizeof(int), AZ_ALLOC,
                                          data_org[AZ_name], "diag_block", &st);
    ipvt       = (int *) AZ_manage_memory(rpntr[new_blks]*sizeof(int), AZ_ALLOC,
                                          data_org[AZ_name], "ipvt", &st);
  }

  else if (options[AZ_pre_calc] <= AZ_recalc) {

    /* compute a new factorization */

    indxsize = bpntr2[blks] + 1;
    valsize  = indx2[bpntr2[blks]] + 1;

    if (options[AZ_overlap] != AZ_full) {
      indxsize += new_blks - blks;
      for (i = blks; i < new_blks; i++) {
        valsize += ((cpntr2[i+1] - cpntr2[i]) * (cpntr2[i+1] - cpntr2[i]));
      }
    }
    else {

      /*
       * Figure out how much space the new entries will take.  That is, for each
       * entry: 1) figure out the block row and col, 2) add this block to the
       * new_blk list if it is a new block, 3) overwrite the point column
       * information with block information.
       */

      i    = 0;
      flag = 1;
      while (i < length) {
        if (buffer[i] == -1) flag = 1;

        while ( (i < length) && (buffer[i] == -1.0) ) i++;

        if (i < length) {
          point_col = (int) buffer[i];
          AZ_which_block_row(&block_col, cpntr2, point_col);

          if (flag == 1) {
            block_row = block_col;
            flag      = 0;
          }

          buffer[i+2] = (double) (buffer[i] - cpntr2[block_col]);
          buffer[i]   = (double) block_col;
          i += 3;
        }
      }

      /* make a rows array */

      therows = (int *) calloc(length/3+1, sizeof(int));
      i = row_count = 0;
      while (i < length) {
        while ( (i < length) && (buffer[i] == -1.0)) i++;

        if (i < length) {
          therows[row_count++] = (int) buffer[i];
          while ( (i < length) && (buffer[i] != -1.0) ) i += 3;
        }
      }

      AZ_sort(therows, row_count, NULL, NULL);

      j = 1;
      for (i = 0; i < row_count; i++) {
        if (therows[i] != therows[j-1]) {
          therows[j] = therows[i];
          j++;
        }
      }

      row_count = j;
      theptrs   = (int *) calloc(row_count + 1, sizeof(int));

      /* estimate the amount of space needed for each row */

      space  = (length - cpntr2[new_blks] + cpntr2[blks]) / 3 + 4 * row_count;
      tbindx = (int *) calloc(space+1, sizeof(int));

      for (i = 0; i <= space; i++) tbindx[i] = -1;

      j = theptrs[0] = 0;
      for (i = 1; i < row_count; i++)
        theptrs[i] = (space - theptrs[i-1])/(1 + row_count - i) + theptrs[i-1];
        theptrs[row_count] = space;

        i    = 0;
        flag = 1;
        while (i < length) {
          if (buffer[i] == -1) flag = 1;

          while ( (i < length) && (buffer[i] == -1.0) ) i++;

          if (i < length) {
            block_col = (int) buffer[i];
            if (flag == 1) {
              block_row = block_col;
              kk        = AZ_find_index(block_row, therows, row_count);
              flag      = 0;
            }

            start = theptrs[kk];
            while ( (tbindx[start] != -1) && (tbindx[start] != block_col) )
              start++;

            if (start >= theptrs[kk+1]) {
              incr     = 10;
              space   += incr;
              newbindx = (int *) malloc(space*sizeof(int));

              for (j = 0; j < space-incr; j++)
                newbindx[j] = tbindx[j];

              free(tbindx);
              tbindx = newbindx;

              /* move things around */

              for (j = space-incr-1; j >= theptrs[kk+1]; j--)
                tbindx[j+incr] = tbindx[j];

              for (j = theptrs[kk+1]; j < theptrs[kk+1]+incr; j++)
                tbindx[j] = -1;

              start = theptrs[kk+1];
              for (j = kk + 1; j <= row_count; j++)
                theptrs[j] += incr;
            }
            if (tbindx[start] == -1) tbindx[start] = block_col;
          }

          i += 3;
        }

        /* compress the storage */

        j = 0;
        for (i = 0; i < row_count; i++) {
          kk = theptrs[i];
          AZ_funswill(&kk);  /* This line is needed so that the puma compiler
                                (with optimization turned on) works properly. */

          value      = tbindx[kk];
          theptrs[i] = j;

          while ((value != -1) && (kk < theptrs[i+1])) {
            tbindx[j] = value;
            kk++; j++;
            value     = tbindx[kk];
          }
        }

        space              = j;
        theptrs[row_count] = space;
        indxsize          += space;
    }

    bpntr = (int *) AZ_manage_memory((new_blks+1)*sizeof(int), AZ_ALLOC,
                                     data_org[AZ_name], "bpntr", &st);
    bindx = (int *) AZ_manage_memory((indxsize+1)*sizeof(int), AZ_ALLOC,
                                     data_org[AZ_name], "bindx", &st);
    indx  = (int *) AZ_manage_memory((indxsize+1)*sizeof(int), AZ_ALLOC,
                                     data_org[AZ_name], "indx", &st);

    if (options[AZ_overlap] == AZ_full) {
      bpntr[blks] = bpntr2[blks];

      for (i = blks+1; i <= new_blks; i++) {
        bpntr[i] = bpntr[i-1] + theptrs[i-blks] - theptrs[i-blks-1];
      }

      for (i = bpntr[blks]; i < bpntr[new_blks]; i++)
        bindx[i] = tbindx[i-bpntr[blks]];

      kk       = bpntr[blks];
      indx[kk] = indx2[kk];

      for (i = 0; i < row_count; i++) {
        for (j = theptrs[i]; j < theptrs[i+1]; j++) {
          blk_size = (cpntr2[blks+i+1] - cpntr2[blks+i]) *
            (cpntr2[tbindx[j]+1] - cpntr2[tbindx[j]]);

          valsize   += blk_size;
          indx[kk+1] = indx[kk] + blk_size;
          kk++;
        }
      }

      free(tbindx); free(theptrs); free(therows);
    }

    rpntr = (int *) AZ_manage_memory((new_blks+1)*sizeof(int), AZ_ALLOC,
                                     data_org[AZ_name], "rpntr", &st);
    cpntr = (int *) AZ_manage_memory((new_blks+1)*sizeof(int), AZ_ALLOC,
                                     data_org[AZ_name], "cpntr", &st);
    val   = (double *) AZ_manage_memory((valsize+1)*sizeof(double),AZ_ALLOC,
                                        data_org[AZ_name], "val", &st);

    for (i = 0; i < valsize + 1; i++) val[i] = 0.0;

    diag_block = (int *) AZ_manage_memory((new_blks+1)*sizeof(int), AZ_ALLOC,
                                          data_org[AZ_name],
                                          "diag_block", &st);

    AZ_vb2bilu_setup(cpntr2, bpntr2, bindx2, indx2, val2, rpntr, cpntr,
                     bpntr, bindx, indx, val, new_N, blks, new_blks, length,
                     buffer, options, data_org);
    AZ_order(new_blks, val, val, bindx, indx, indx, bpntr, diag_block);

    N      = rpntr[new_blks];
    z      = (double *) malloc(sizeof(double) * (N));
    ipvt   = (int *)    malloc(sizeof(int)    * (N));
    block1 = (double *) malloc(sizeof(double) * (N));
    block2 = (double *) malloc(sizeof(double) * (N));

    for (i = 0; i < new_blks; i++) {  /* loop through block rows */

      /* compute the lower [L] sub-matrix */

      for (I = bpntr[i]; I < diag_block[i]; I++) {

        j = bindx[I];

        for (K_SUM = bpntr[i]; K_SUM < diag_block[i]; K_SUM++) {
          k = bindx[K_SUM];

          if (k < j) {
            flag_block = AZ_get_block(j, k, bindx, bpntr, &ptr_j);

            if (flag_block != -1) {
              AZ_update_block(I, K_SUM, ptr_j, val, indx, bindx, cpntr);
            }
          }
        }

        /* compute D[I] = D[I] / D(DIAG[j]) */

        diag_j = diag_block[j];
        AZ_divide_block(I, diag_j, val, indx, bindx, cpntr, z, block1, block2,
                        ipvt);
      }

      /* compute the upper [U] sub-matrix */

      for (I = diag_block[i]; I < bpntr[i+1]; I++) {
        j = bindx[I];

        for (K_SUM = bpntr[i]; K_SUM < diag_block[i]; K_SUM++) {
          k          = bindx[K_SUM];
          flag_block = AZ_get_block(j, k, bindx, bpntr, &ptr_j);

          if (flag_block != -1) {
            AZ_update_block(I, K_SUM, ptr_j, val, indx, bindx, cpntr);
          }
        }
      }
    }

    free((void *) z);
    free((void *) ipvt);
    free((void *) block1);
    free((void *) block2);

    ipvt = (int *) AZ_manage_memory(rpntr[new_blks]*sizeof(int), AZ_ALLOC,
                                    data_org[AZ_name], "ipvt", &st);

    for (i = new_blks - 1; i >= 0; i--) {
      m1   = rpntr[i+1] - rpntr[i];
      ival = indx[diag_block[i]];
      iy   = rpntr[i];
      dgetrf_( &m1, &m1, val+ival, &m1, &(ipvt[iy]), &info);
    }
  }

  else if (previous_factors != data_org[AZ_name] ) {
    (void) fprintf(stderr, "Warning: Using a previous factorization as a "
                   "preconditioner\neven though matrix "
                   "(data_org[AZ_name]) has changed\n");
  }

  previous_factors = data_org[AZ_name];

  x_tmp = (double *) calloc(new_N, sizeof(double));
  AZ_lower_triang_vbr_solve(new_N, new_blks, val, indx, bindx, rpntr, cpntr,
                            bpntr, diag_block, x_tmp, x);
  AZ_upper_triang_vbr_solve(new_N, new_blks, val, indx, bindx, rpntr, cpntr,
                            bpntr, diag_block, x, x_tmp, ipvt);

  free(x_tmp);

} /* ilu_vbr */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int AZ_get_block(int j, int k, int bindx[], int bpntr[], int *ptr_j)

/*******************************************************************************

  Return in *ptr_j the index in bindx[] corresponding to the (k,j)th block
  of the matrix. get_block() returns -1 if the block is not found.

  Author:          Lydie Prevost, SNL, 1422
  =======

  Return code:     int, 1 (block found) or -1 (block not found)
  ============

  Parameter list:
  ===============

  k, j:            On input, get_block() looks for the (k,j)th block in the
                   matrix.

  bindx,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file params.txt).

  ptr_j:           On output, index into bindx[] corresponding to the (k,j)th
                   block of the matrix.

*******************************************************************************/

{

  /* local variables */

  int kk;

  /**************************** execution begins ******************************/

  for (kk = bpntr[k]; kk < bpntr[k+1]; kk++){
    if (bindx[kk] == j) {
      *ptr_j = kk;
      return 1;
    }
  }

  return -1;

} /* AZ_get_block */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_update_block(int i, int k, int j, double val[], int indx[], int bindx[],
                     int cpntr[])

/*******************************************************************************

  Update block i for VBR sparse  format compute

        block[i] = block[i] - block[k]*block[j]

  Author:          Lydie Prevost, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  i, k, j:         Indices.

  val:             Array containing the nonzero entries of the matrix (see file
                   params.txt).

  indx,
  bindx,
  cpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file params.txt).

*******************************************************************************/

{

  /* local variables */

  int    nbelmt_i, nbelmt_k, nbelmt_j, colblock_i, colblock_k, colblock_j;
  int    ncol_i, nrow_i, ncol_k, nrow_k, ncol_j, nrow_j;
  int    m, n, kk, lda, ldb, ldc;
  double alpha = -1.0, beta = 1.0;
  char   T = 'N';
  char   None[2];

  /**************************** execution begins ******************************/

  strcpy(None, "N");

  nbelmt_i   = indx[i+1] - indx[i];
  nbelmt_k   = indx[k+1] - indx[k];
  nbelmt_j   = indx[j+1] - indx[j];
  colblock_i = bindx[i];
  colblock_k = bindx[k];
  colblock_j = bindx[j];
  ncol_i     = cpntr[colblock_i+1] - cpntr[colblock_i];
  ncol_k     = cpntr[colblock_k+1] - cpntr[colblock_k];
  ncol_j     = cpntr[colblock_j+1] - cpntr[colblock_j];
  nrow_i     = nbelmt_i / ncol_i;
  nrow_k     = nbelmt_k / ncol_k;
  nrow_j     = nbelmt_j / ncol_j;

  m   = nrow_k;
  n   = ncol_j;
  kk  = ncol_k;
  lda = nrow_k;
  ldb = nrow_j;
  ldc = nrow_i;

  dgemm_(&T, &T, &m, &n, &kk, &alpha, val+indx[k], &lda, val+indx[j], &ldb,
         &beta, val+indx[i], &ldc, strlen(None), strlen(None));

} /* AZ_update_block */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_divide_block(int i, int j, double val[], int indx[], int bindx[],
                     int cpntr[], double *z, double *blockj, double *blocki,
                     int *ipvt)

/*******************************************************************************

  Compute

     block[i] = block[i] / block[j]

  for  VBR sparse  format

  Author:          Lydie Prevost, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  i, j:            Indices.

  val:             Array containing the nonzero entries of the matrix (see file
                   params.txt).

  indx,
  bindx,
  cpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file params.txt).

  z:               Double precision work space array.

  blockj:          Points to block[j] in the above description.

  blocki:          Points to block[i] in the above description.

  ipvt:            Pivot vector??

*******************************************************************************/

{

  /* local variables */

  int    nbelmt_i, nbelmt_j, colblock_i, colblock_j;
  int    ncol_i, nrow_i, ncol_j, nrow_j;
  int    m, n, kk, lda, ldb, ldc, job;
  double det[2];
  double alpha = 1.0, beta = 0.0, rcond;
  char   T = 'N';
  char   None[2];

  /**************************** execution begins ******************************/

  strcpy(None, "N");

  nbelmt_i   = indx[i+1] - indx[i];
  nbelmt_j   = indx[j+1] - indx[j];
  colblock_i = bindx[i];
  colblock_j = bindx[j];
  ncol_i     = cpntr[colblock_i+1] - cpntr[colblock_i];
  ncol_j     = cpntr[colblock_j+1] - cpntr[colblock_j];
  nrow_i     = nbelmt_i / ncol_i;
  nrow_j     = nbelmt_j / ncol_j;

  /* compute the inverse of block j  */

  for (kk = 0; kk < nbelmt_j; kk++) {
    blockj[kk] = val[indx[j]+kk];
  }

  for (kk = 0; kk < nbelmt_i; kk++) {
    blocki[kk] = val[indx[i]+kk];
  }

  dgeco_(blockj, &nrow_j, &nrow_j, ipvt, &rcond, z);

  job = 01;
  dgedi_(blockj, &nrow_j, &nrow_j, ipvt, det, z, &job);

  m   = nrow_i;
  n   = ncol_j;
  kk  = ncol_i;
  lda = m;
  ldb = kk;
  ldc = m;

  dgemm_(&T, &T, &m, &n, &kk, &alpha, blocki, &lda, blockj, &ldb, &beta,
         val+indx[i], &ldc, strlen(None), strlen(None));

} /* AZ_divide_block */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_divide_block0(int i, int j, double val[], int indx[], int bindx[],
                      int cpntr[], int *ipvt)

/*******************************************************************************

  Compute

     block[i] = block[i] / block[j]

  for  VBR sparse  format

  Author:          Lydie Prevost, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  i, j:            Indices.

  val:             Array containing the nonzero entries of the matrix (see file
                   params.txt).

  indx,
  bindx,
  cpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file params.txt).

  z:

  ipvt:            Pivot vector??

*******************************************************************************/

{

  /* local variables */

  int    nbelmt_i, nbelmt_j, colblock_i, colblock_j;
  int    ncol_i, nrow_i, ncol_j, nrow_j;
  int    ival, info;
  int    ione = 1;
  char   *T = "N";

  /**************************** execution begins ******************************/

  nbelmt_i   = indx[i+1] - indx[i];
  nbelmt_j   = indx[j+1] - indx[j];
  colblock_i = bindx[i];
  colblock_j = bindx[j];
  ncol_i     = cpntr[colblock_i+1] - cpntr[colblock_i];
  ncol_j     = cpntr[colblock_j+1] - cpntr[colblock_j];
  nrow_i     = nbelmt_i / ncol_i;
  nrow_j     = nbelmt_j / ncol_j;
  ival       = indx[j];

  dgetrf_(&nrow_j, &colblock_j, val+ival, &nrow_j, ipvt, &info);

  for (i = 0; i < ncol_i; i++)
    dgetrs_(T, &nrow_i, &ione, val+ival, &nrow_j, ipvt, &(val[indx[i]]),
            &nrow_j, &info, strlen(T));

} /* divide_block0 */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_vb2bilu_setup(int cpntr2[], int bpntr2[], int bindx2[],
                      int indx2[], double val2[], int rpntr[], int cpntr[],
                      int bpntr[], int bindx[], int indx[], double val[],
                      int newN, int old_blks, int new_blks,int length, double
                      buffer[], int options[], int data_org[])

/*******************************************************************************

  Using the information in (cpntr2, bpntr2, bindx2, indx2, val2) and
  (buffer, length) create the extended matrix and store it in the arrays (rpntr,
  cpntr, bpntr, bindx, indx, val)

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  cpntr2,
  bpntr2,
  bindx2,
  indx2:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file params.txt).

  val2:            Array containing the nonzero entries of the matrix (see file
                   params.txt).

  rpntr,
  cpntr,
  bpntr,
  bindx,
  indx:            Arrays used for DMSR and DVBR sparse matrix storage (see
                   file params.txt).

  val:             Array containing the nonzero entries of the matrix (see file
                   params.txt).

  diag_block:      Work array of size n points to each diagonal block.

  oldN:            data_org[AZ_N_internal] + data_org[AZ_N_border].

  newN:            data_org[AZ_N_internal] + data_org[AZ_N_border] +
                   data_org[AZ_N_external]

  old_blks:        data_org[AZ_N_int_blks] + data_org[AZ_N_bord_blks]

  new_blks:        data_org[AZ_N_int_blks] + data_org[AZ_N_bord_blks] +
                   data_org[AZ_N_ext_blks]

  length:          Length of buffer.

  buffer:          Buffer holding extra rows received from neighbors.

  options:         Determines specific solution method and other parameters.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see file params.txt).

*******************************************************************************/

{

  /* local variables */

  int     i, j, ii, jj;
  double *diag;
  int     kk, block_col, blk_size;
  int     flag, col_offset, block_row, row_offset, num_blk_rows;

  /**************************** execution begins ******************************/

  /* first copy the internal submatrix */

  for (i = 0; i <= old_blks; i++) {
    bpntr[i] = bpntr2[i];
  }

  for (i = 0; i <= new_blks; i++) {
    cpntr[i] = cpntr2[i];
    rpntr[i] = cpntr2[i];
  }

  for (i = 0; i < bpntr2[old_blks]; i++) {
    indx[i]  = indx2[i];
    bindx[i] = bindx2[i];
  }

  i       = bpntr2[old_blks];
  indx[i] = indx2[i];

  for (i = 0; i < indx2[bpntr2[old_blks]]; i++) {
    val[i] = val2[i];
  }

  /* now add new entries */

  if (options[AZ_overlap] == AZ_full) {
    i = 0;

    while (i < length) {
      if (buffer[i] == -1) flag = 1;

      while ( (i < length) && (buffer[i] == -1.0) ) i++;

      if (i < length) {
        block_col  = (int) buffer[i];
        col_offset = (int) buffer[i+2];

        if (flag == 1) {
          block_row    = block_col;
          row_offset   = col_offset;
          flag         = 0;
          num_blk_rows = rpntr[block_row+1] - rpntr[block_row];
        }

        for (j = bpntr[block_row]; j < bpntr[block_row+1]; j++) {
          if (bindx[j] == block_col ) break;
        }

        val[indx[j] + col_offset*num_blk_rows + row_offset] = buffer[i+1];
      }
      i += 3;
    }
  }

  else if (options[AZ_overlap] == AZ_none) {

    /* add Dirichlet */

    for (i = old_blks; i < new_blks; i++) {
      bpntr[i+1] = bpntr[i] + 1;
      ii         = bpntr[i];
      indx[ii+1] = indx[ii] + (cpntr[i+1] - cpntr[i]) * (cpntr[i+1] - cpntr[i]);
      bindx[ii]  = i;
      ii         = indx[ii];

      for (j = cpntr[i]; j < cpntr[i+1]; j++) {
        for (jj = cpntr[i]; jj < cpntr[i+1]; jj++) {
          if (j == jj) val[ii++] = 1.0;
          else         val[ii++] = 0.0;
        }
      }
    }
  }

  else if (options[AZ_overlap] == AZ_diag) {
    diag = (double *) calloc(newN, sizeof(double));

    /* must find diagonal element for each real row */

    for (i = 0; i < old_blks; i++) {
      for (j = rpntr[i]; j < rpntr[i+1]; j++) {
        diag[j] = 0.0;

        for (kk = bpntr[i]; kk < bpntr[i+1]; kk++) {
          block_col = bindx[kk];
          blk_size  = cpntr[block_col+1] - cpntr[block_col];

          if (block_col == i)
            diag[j] = val[indx[kk]+(blk_size+1)*(j-rpntr[i])];
        }
      }
    }

    AZ_exchange_bdry(diag, data_org);

    for (i = old_blks; i < new_blks; i++) {
      bpntr[i+1] = bpntr[i] + 1;
      ii         = bpntr[i];
      indx[ii+1] = indx[ii] + (cpntr[i+1] - cpntr[i]) * (cpntr[i+1] - cpntr[i]);
      bindx[ii]  = i;
      ii         = indx[ii];

      for (j = cpntr[i]; j < cpntr[i+1]; j++) {
        for (jj = cpntr[i]; jj < cpntr[i+1]; jj++) {
          if (j == jj) val[ii++] = diag[j];
          else         val[ii++] = 0.0;
        }
      }
    }

    free(diag);
  }

} /* AZ_vb2bilu_setup */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_funswill(int *trash)

{
   /* Just some garbage which does nothing. This will fool the */
   /* lint compiler so we don't have lint warnings.            */

   *trash += 1;    *trash -= 1;
}
