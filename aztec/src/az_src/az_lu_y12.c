/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile: az_lu_y12.c,v $
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
static char rcsid[] = "$Id: az_lu_y12.c,v 1.1.1.1 2001-01-30 20:59:14 paklein Exp $";
#endif


/*******************************************************************************
 * Copyright 1995, Sandia Corporation.  The United States Government retains a *
 * nonexclusive license in this software as prescribed in AL 88-1 and AL 91-7. *
 * Export of this program may require a license from the United States         *
 * Government.                                                                 *
 ******************************************************************************/


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include "az_aztec.h"

void AZ_lu_y12m(int n, int nonzeros, double x[], int length, double buffer[],
                double val[], int indx[], int bindx[], int rpntr[], int cpntr[],
                int bpntr[], int options[], int data_org[],
                double drop_tolerance)

/*******************************************************************************

  C driver for sparse direct factorization solver. This routine takes a matrix
  (either msr or vbr) and translates it to a format suitable for the sparse
  solver. The matrix is factorized and the solved. On subsequent calls to lu(),
  the previously computed factors are used to backsolve the system.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  n:               Order of the new system (with external rows).

  nonzeros:        Number of nonzeros in original system.

  x:               On input, contains the rhs of subsystem to be solved. On
                   output contains the result.

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

  options:         Determines specific solution method and other parameters.
                   options[AZ_overlap] = AZ_none
                                         nonoverlapping domain decomposition
                                         (don't use external values).
                                       = AZ_diag
                                         Use external values in conjunction with
                                         the external row, however only keep the
                                         diagonal of the external row and set up
                                         a Dirichlet-type boundary condition.
                                       = AZ_full
                                         duplicate external rows on this
                                         processor (only keeping the coupling
                                         between variables contained within
                                         this processor).

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see file params.txt).

  drop_tolerance:  Drop tolerance used to compute an incomplete factorization.

*******************************************************************************/

{
  /* local variables */

  static double *newa, *aflag;
  static double *pivot;
  static int    *snr, *rnr, *iflag, z, nn1;
  static int     nn, iha;
  static int    *ha;
  static int     previous_factors = -1;

  int            i, j, ifail;
  int            block_col, kk, blk_size;
  int            tind;
  int            blks, oldN;
  int            st;
  int           *little;
  double        *diag;
  char          *yo = "AZ_lu_y12m: ";

  /**************************** execution begins ******************************/

  oldN = data_org[AZ_N_internal] + data_org[AZ_N_border];
  blks = data_org[AZ_N_int_blk]  + data_org[AZ_N_bord_blk];

  if (n == 0) return;

  if (options[AZ_pre_calc] <= AZ_reuse) {
    little = (int *) AZ_manage_memory(8*sizeof(int), AZ_ALLOC,
                                      data_org[AZ_name], "little", &st);

    iha   = n;
    pivot = (double *) AZ_manage_memory(n*sizeof(double), AZ_ALLOC,
                                        data_org[AZ_name], "pivot",
                                        &st);
    ha    = (int *) AZ_manage_memory(11*iha*sizeof(int), AZ_ALLOC,
                                     data_org[AZ_name], "ha", &st);

    if ((options[AZ_pre_calc] == AZ_reuse) && (st == AZ_NEW_ADDRESS)) {
      (void) fprintf(stderr, "%sERROR: options[AZ_pre_calc] == AZ_reuse and"
                     " previous factors\n       not found.  "
                     "Check data_org[AZ_name].\n", yo);
      exit(-1);
    }

    if (options[AZ_pre_calc] != AZ_reuse) {

      /* choose enough space to hopefully accomodate the factors */

      nn1 = nn = 8 * nonzeros;
      nn1 = nn = AZ_MAX_MEMORY_SIZE/(3*sizeof(int)+sizeof(double));


      if (2 * nonzeros < 0 ) {
        (void) fprintf(stderr, "%sERROR: not enough memory for lu\n"
                       "       Run a smaller problem, increase\n"
                       "       drop tolerance, or use another\n"
                       "       preconditioner.\n", yo);
        exit(-1);
      }

      snr  = (int    *) calloc(nn, sizeof(int ));
      rnr  = (int    *) calloc(nn1, sizeof(int ));
      newa = (double *) calloc(nn, sizeof(double));

      while (newa == NULL) {
        if (snr != NULL) free(snr);
        if (rnr != NULL) free(rnr);

        nn  = (int) (((double) nn) * .75);
        nn1 = nn;

        if (nn < 2 * nonzeros) {
          (void) fprintf(stderr, "%sERROR: not enough memory for lu\n"
                         "       Run a smaller problem, increase\n"
                         "       drop tolerance, or use another\n"
                         "       preconditioner.\n", yo);
          exit(-1);
        }

        snr  = (int    *) calloc(nn, sizeof(int ));
        rnr  = (int    *) calloc(nn1, sizeof(int ));
        newa = (double *) calloc(nn, sizeof(double));
      }

      nn  = (int) (((double) nn) * .96);
      nn1 = nn;

      little[0] = nn;
      little[1] = nn1;
      little[2] = iha;
      little[3] = z;

      free(snr); free(rnr); free(newa);
    }
    else {
      nn  = little[0];
      nn1 = little[1];
      iha = little[2];
      z   = little[3];
    }

    snr   = (int *)    AZ_manage_memory(nn*sizeof(int), AZ_ALLOC,
                                        data_org[AZ_name], "snr", &st);
    rnr   = (int *)    AZ_manage_memory(nn1*sizeof(int), AZ_ALLOC,
                                        data_org[AZ_name], "rnr", &st);
    newa  = (double *) AZ_manage_memory(nn*sizeof(double), AZ_ALLOC,
                                        data_org[AZ_name], "newa", &st);
    aflag = (double *) AZ_manage_memory(8*sizeof(double), AZ_ALLOC,
                                        data_org[AZ_name], "aflag",
                                        &st);
    iflag = (int *)    AZ_manage_memory(10*sizeof(int), AZ_ALLOC,
                                        data_org[AZ_name], "iflag",
                                        &st);
  }

  else if (previous_factors != data_org[AZ_name] ) {
    (void) fprintf(stderr, "%sWARNING: Using a previous factorization as a"
                   "preconditioner\neven though matrix"
                   "(data_org[AZ_name]) has changed\n", yo);
  }

  previous_factors = data_org[AZ_name];

  if (options[AZ_pre_calc] <= AZ_recalc) {

    /* calculate factorization */

    /* initialize constants here */

    z = 0;

    /* set up flags for the sparse factorization solver */

    iflag[0] = 1;    iflag[1] = 2;    iflag[2] = 1;     iflag[3] = 0;

    /* if matrix is pos def, iflag[2] = 2 */
    /* is cheaper */

    iflag[4] = 2;    aflag[0] = 16.0; aflag[2] = 1.0e8; aflag[3] = 1.0e-12;
    aflag[1] = drop_tolerance;

    /* translate the msr format to a new format suitable for the solver */

    if (data_org[AZ_matrix_type] == AZ_VBR_MATRIX)
      AZ_vb2lu(blks, val, indx, bindx, rpntr, cpntr, bpntr, newa, snr, rnr, &z);
    else AZ_msr2lu(oldN, val, newa, snr, rnr, &z, bindx);

    /* now add the new entries */

    if (options[AZ_overlap] == AZ_full) {
      i = 0;

      while (i < length) {
        while ((buffer[i] == -1.0) && (i < length)) i++;

        if (i < length) {
          tind     = (int) buffer[i];
          newa[z]  = buffer[i+1];
          rnr[z]   = tind + 1;
          snr[z++] = tind + 1;
          i       += 3;

          while (buffer[i] != -1.0) {
            newa[z]  = buffer[i+1];
            rnr[z]   = tind + 1;
            snr[z++] = ((int) buffer[i]) + 1;
            i       += 3;
          }
        }
      }

      if (buffer != 0) free((char *) buffer);
    }

    else if (options[AZ_overlap] == AZ_none) {

      /* add Dirichlet */

      for (i = oldN; i < n; i++) {
        newa[z]  = 1.0;
        rnr[z]   = i + 1;
        snr[z++] = i + 1;
      }
    }

    else if (options[AZ_overlap] == AZ_diag) {
      diag = (double *) calloc(n, sizeof(double));

      if (data_org[AZ_matrix_type] != AZ_VBR_MATRIX) {
        for (i = 0; i < oldN; i++)
          diag[i] = val[i];
      }

      else {

        /* must find diagonal element for each real row */

        for (i = 0; i < blks; i++) {
          for (j = rpntr[i]; j < rpntr[i+1]; j++) {
            diag[j] = 0.0;

            for (kk = bpntr[i]; kk < bpntr[i+1]; kk++) {
              block_col = bindx[kk];
              blk_size  = cpntr[block_col+1] - cpntr[block_col];

              if (block_col == i)
                diag[j] = val[indx[kk]+(blk_size+1) * (j-rpntr[i])];
            }
          }
        }
      }

      AZ_exchange_bdry(diag, data_org);

      for (i = oldN; i < n; i++) {
        newa[z] = diag[i];
        rnr[z]  = i + 1;
        snr[z]  = i + 1;
        z++;
      }
      free(diag);
    }

    if (z == 1) {
      newa[0] = 1./newa[0];
      ifail = 0;
      little[3] = z;
      x[0] *= newa[0];
    }
    else {
      (void) AZ_factandsolve(newa, aflag, pivot, x, snr, rnr, ha, iflag, &z,
                           &ifail, &nn, &n, &iha, &nn1);

      /* see if we can reduce the storage */

      if ( (3*iflag[7]/2 < nn) && (ifail ==0)) {
        nn   = iflag[7];
        nn1  = iflag[7];
        little[0] = nn;
        little[1] = nn1;

        snr  = (int    *) realloc(snr, nn*sizeof(int ));
        rnr  = (int    *) realloc(rnr, nn1*sizeof(int ));
        newa = (double *) realloc(newa, nn*sizeof(double));

        snr   = (int *)    AZ_manage_memory(nn*sizeof(int), AZ_REALLOC,
                                            data_org[AZ_name], "snr", &st);
        rnr   = (int *)    AZ_manage_memory(nn1*sizeof(int), AZ_REALLOC,
                                            data_org[AZ_name], "rnr", &st);
        newa  = (double *) AZ_manage_memory(nn*sizeof(double), AZ_REALLOC,
                                            data_org[AZ_name], "newa", &st);
      }
    }
  }

  else {
    if (z == 1) {
      x[0] *= newa[0];
      ifail = 0;
    }
    else
      (void) AZ_backsolve(newa, pivot,x, snr, ha, iflag, &ifail, &nn, &n, &iha);
  }

  if (ifail != 0) {
    (void) fprintf(stderr, "direct: ifail is not zero (%d)\n", ifail);
    switch(ifail) {
    case 4:
      (void) fprintf(stderr, "Large growth factor\n");
    break;
    case 3:
      (void) fprintf(stderr, "Matrix may be singular\n");
    break;
    case 5:
      (void) fprintf(stderr, "Allocated space not large enough\n");
    break;
    case 6:
      (void) fprintf(stderr, "Allocated space not large enough\n");
    break;
    case 7:
      (void) fprintf(stderr, "Either the matrix may be singular\n"
                     "or the drop tolerance may be too high\n");
    break;
    case 8:
      (void) fprintf(stderr, "Either the matrix may be singular\n"
                     "or the drop tolerance may be too high\n");
    break;
    case 11:
      (void) fprintf(stderr, "two elements in the same (i,j) position\n");
    break;
    case 12:
      (void) fprintf(stderr, "System has less than two rows\n");
    break;
    case 17:
      (void) fprintf(stderr, "A row without nonzero elements found\n");
    break;
    case 18:
      (void) fprintf(stderr, "A column without nonzero elements found\n");
    break;
    case 24:
      (void) fprintf(stderr, "A column index exceeds matrix dimension \n");
    break;
    case 25:
      (void) fprintf(stderr, "A row index exceeds matrix dimension \n");
    break;
    default:
      break;
    }
    (void) fprintf(stderr, "        Check y12m manual for more information.\n");
      exit(-1);
  }

} /* AZ_lu_y12m */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_factandsolve(double val[], double aflag[], double pivot[], double b[],
                     int snr[], int rnr[], int ha[], int iflag[], int *z,
                     int *ifail, int *nn, int *n, int *iha, int *nn1)

/*******************************************************************************

  Routines which call the correct sequence of routines for the sparse direct
  solver: y12m

  The subroutine 'AZ_factandsolve' first factorizes the matrix and then performs
  a backsolve. The second subroutine 'backsolve' uses precomputed factors to
  perform a back solve.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  val:             Array containing the nonzero entries of the matrix (see file
                   params.txt).

  aflag:

  pivot:

  b:

  snr, rnr:

  ha:

  iflag:

  z:

  ifail:

  nn:

  n:

  iha:

  nn1:

*******************************************************************************/

{

  *ifail = 0;

  y12mbf_(n, z, val, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail);

  if (*ifail != 0) return;

  y12mcf_(n, z, val, snr, nn, rnr, nn1, pivot, b, ha, iha, aflag, iflag, ifail);

  if (*ifail !=  0) return;

  y12mdf_(n, val, nn, b, pivot, snr, ha, iha, iflag, ifail);

} /* AZ_factandsolve */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_backsolve(double val[], double pivot[], double b[], int snr[], int ha[],
                  int iflag[], int *ifail, int *nn, int *n, int *iha)

/*******************************************************************************

  Routine which calls the correct sequence of routines for the sparse direct
  solver: y12m

  The subroutine 'backsolve' performs a backsolve using precomputed factors.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  val:             Array containing the nonzero entries of the matrix (see file
                   params.txt).

  pivot:

  b:

  snr:

  ha:

  iflag:

  ifail:

  nn:

  n:

  iha:

*******************************************************************************/

{

  iflag[4] = 3;
  y12mdf_(n, val, nn, b, pivot, snr, ha, iha, iflag, ifail);

} /* backsolve */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_msr2lu(int oldN, double val[], double newa[], int snr[], int rnr[],
               int *z, int bindx[])

/*******************************************************************************

  This routine converts an MSR matrix into a format suitable for the lu solver.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  oldN:

  val:             Array containing the nonzero entries of the matrix (see file
                   params.txt).

  newa:

  snr, rnr:

  z:

  bindx:           MSR index vector.

*******************************************************************************/

{

  /* local variables */

  int i, j, start, end;

  /**************************** execution begins ******************************/

  for (i = 0; i < oldN; i++) {
    newa[i] = val[i];
    snr[i]  = i + 1;
    rnr[i]  = i + 1;
  }

  *z = oldN;

  for (i = 0; i < oldN; i++) {
    start = bindx[i];
    end   = bindx[i+1];

    for (j = start; j < end; j++) {
      newa[*z] = val[j];
      rnr[*z]  = i + 1;
      snr[*z]  = bindx[j] + 1;
      (*z)++;
    }
  }

} /* AZ_msr2lu */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_vb2lu(int m, double val[], int indx[], int bindx[], int rpntr[],
              int cpntr[], int bpntr[], double newa[], int snr[], int rnr[],
              int *nz_ptr)

/*
 * This routine converts a VBR matrix into a format suitable for the lu solver.
 *
 * Author:         Ray Tuminaro, Div 1422 SNL
 * Date:           15 December 1994
 *
 *
 * Paramter list:
 *
 *      m       ==  Number of (block) rows in the matrix
 *      a       ==  matrix A in sparse format (VBR)
 *      indx    ==  The ith element of indx points to the location in
 *                   val of the (0,0) entry of the ith block entry. The
 *                   last element is the number of nonzero entries of
 *                   matrix A plus one.
 *      bindx   ==  Contains the block column indices of the non-zero
 *                   block entries.
 *      rpntr   ==  The ith element of rpntr indicates the first point
 *                   row in the ith block row. The last element is the
 *                   number of block rows plus one.
 *      cpntr   ==  The jth element of cpntr indicates the first point
 *                   column in the jth block column. The last element
 *                   is the number of block columns plus one.
 *      bpntr   ==   The ith element of bpntr points to the first block
 *                   entry of the ith row in bindx. The last element is
 *                   the number of nonzero blocks of matrix A plus one.
 *      rnr     ==   On output, rnr[i] is the row number of the ith element.
 *      snr     ==   On output, snr[i] is the column number of the ith element.
 *      newa    ==   On output, newa[i] is nonzero value of the ith element.
 */

{

  /* local variables */

  register int   indx_ptr, irow, iblk_row, jcol, jblk_col, icol_blk, iblk = 0;
  int            num_blk_rows, num_col_blks, num_blk_cols;
  int            realrow, realcol;

  /**************************** execution begins ******************************/

  /* loop over the block rows */

  for (iblk_row = 0; iblk_row < m; iblk_row++) {

    /* find out how many rows are in this block row */

    num_blk_rows = rpntr[iblk_row+1] - rpntr[iblk_row];

    /* find out how many block columns are in this block row */

    num_col_blks = bpntr[iblk_row+1] - bpntr[iblk_row];

    /* loop over all the rows in this block row */

    for (irow = 0; irow < num_blk_rows; irow++) {
      realrow = irow + rpntr[iblk_row];

      /* loop over all the column blocks in this block row */

      for (jblk_col = 0; jblk_col < num_col_blks; jblk_col++) {

        /* find out which column block this is */

        icol_blk = bindx[iblk];
        indx_ptr = indx[iblk++];

        /* find out how many columns are in this block */

        num_blk_cols = cpntr[icol_blk+1] - cpntr[icol_blk];

        /* loop over all the columns in this block */

        for (jcol = 0; jcol < num_blk_cols; jcol++) {
          realcol = jcol + cpntr[icol_blk];

          /* store a(realrow,realcol) into the MSR matrix */

          if (realcol == realrow) {
            newa[*nz_ptr] = val[indx_ptr + jcol*num_blk_rows + irow];
            snr[*nz_ptr]  = realrow+1;
            rnr[*nz_ptr]  = realrow+1;
            (*nz_ptr)++;
          }
          else {
            if (val[indx_ptr + jcol*num_blk_rows + irow] != 0.0) {
              newa[*nz_ptr] = val[indx_ptr + jcol*num_blk_rows + irow];
              snr[*nz_ptr]  = realcol + 1;
              rnr[*nz_ptr]  = realrow + 1;
              (*nz_ptr)++;
            }
          }
        }
      }

      iblk -= num_col_blks;
    }

    iblk += num_col_blks;
  }

} /* vb2msr */
