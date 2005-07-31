/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile: az_sparax.c,v $
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
static char rcsid[] = "$Id: az_sparax.c,v 1.1.1.1 2001-01-30 20:59:12 paklein Exp $";
#endif


/******************************************************************************
 * Copyright 1995, Sandia Corporation.  The United States Government retains a
 * nonexclusive license in this software as prescribed in AL 88-1 and AL 91-7.
 * Export of this program may require a license from the United States
 * Government.
 *****************************************************************************/


#include <stdio.h>
#include <float.h>
#include <string.h>
#include "az_aztec.h"

/* missing prototypes - PAK (07/24/98) */
void dvbr_sparax_basic(int m, double *val, int *indx, int *bindx, int *rpntr,
                       int *cpntr, int *bpntr, double *b, double *c,
                       int exchange_flag, int *data_org);
                       
void dvbr_sparax_basic(int m, double *val, int *indx, int *bindx, int *rpntr,
                       int *cpntr, int *bpntr, double *b, double *c,
                       int exchange_flag, int *data_org)

/******************************************************************************

  c = Ab:
  Sparse (square) matrix-vector multiply, using the variable block row (VBR)
  data structure (A = val).

  Author:          Scott A. Hutchinson, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  m:               Number of (block) rows in A.

  val:             Array containing the entries of the matrix. The matrix is
                   stored block-row-by-block-row. Each block entry is dense and
                   stored by columns (VBR).

  indx,
  bindx,
  rpntr,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file params.txt).

  b:               Right hand side of linear system.

  c:               On output contains the solution to the linear system.

  exchange_flag:   Flag which controls call to AZ_exchange_bdry() (ignored in
                   serial implementation).

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see file params.txt).

******************************************************************************/

{

  /* local variables */

  register double *x;
  register double *c_pntr;
  register int     iblk_row, j, jblk, iblk_size;
  int              m1, ib1, n1;
  int              bpoff, rpoff;
  int              ione = 1;
  int              irpntr, irpntr_next;
  int              ibpntr, ibpntr_next = 0;
  double           one = 1.0;
  double          *val_pntr;
  char            *N = "N";

  /**************************** execution begins *****************************/

  /* exchange boundary info */

  if (exchange_flag) AZ_exchange_bdry(b, data_org);

  /* offset of the first block */

  bpoff = *bpntr;
  rpoff = *rpntr;

  /* zero the result vector */

  for (j = 0; j < rpntr[m] - rpoff; c[j++] = 0.0);

  val_pntr    = val;
  irpntr_next = *rpntr++;
  bpntr++;
  c          -= rpoff;

  /* loop over block rows */

  for (iblk_row = 0; iblk_row < m; iblk_row++) {

    irpntr      = irpntr_next;
    irpntr_next = *rpntr++;

    ibpntr      = ibpntr_next;
    ibpntr_next = *bpntr++ - bpoff;

    /* set result pointer */

    c_pntr      = c + irpntr;

    /* number of rows in the current row block */

    m1          = irpntr_next - irpntr;

    /* loop over each block in the current row block */

    for (j = ibpntr; j < ibpntr_next; j++) {
      jblk = *(bindx+j);

      /* the starting point column index of the current block */

      ib1 = *(cpntr+jblk);

      /* number of columns in the current block */

      n1 = cpntr[jblk + 1] - ib1;
      iblk_size = m1*n1;

      /****************** Dense matrix-vector multiplication *****************/

      /*
       * Get base addresses
       */

      x = b + ib1;


      /*
       * Special case the m1 = n1 = 1 case
       */

      if (iblk_size == 1)
        *c_pntr += *val_pntr * *x;

      else if (m1 == n1) {

        /*
         * Inline small amounts of work
         */

        switch (m1) {

        case 2:
          c_pntr[0] += val_pntr[0]*x[0] + val_pntr[2]*x[1];
          c_pntr[1] += val_pntr[1]*x[0] + val_pntr[3]*x[1];
          break;

        case 3:
          c_pntr[0] += val_pntr[0]*x[0] + val_pntr[3]*x[1] + val_pntr[6]*x[2];
          c_pntr[1] += val_pntr[1]*x[0] + val_pntr[4]*x[1] + val_pntr[7]*x[2];
          c_pntr[2] += val_pntr[2]*x[0] + val_pntr[5]*x[1] + val_pntr[8]*x[2];
          break;

        case 4:
          c_pntr[0] += val_pntr[0]*x[0] + val_pntr[4]*x[1] + val_pntr[8] *x[2]
            + val_pntr[12]*x[3];
          c_pntr[1] += val_pntr[1]*x[0] + val_pntr[5]*x[1] + val_pntr[9] *x[2]
            + val_pntr[13]*x[3];
          c_pntr[2] += val_pntr[2]*x[0] + val_pntr[6]*x[1] + val_pntr[10]*x[2]
            + val_pntr[14]*x[3];

          c_pntr[3] += val_pntr[3]*x[0] + val_pntr[7]*x[1] + val_pntr[11]*x[2]
            + val_pntr[15]*x[3];
          break;

        case 5:
          c_pntr[0] += val_pntr[0]*x[0] + val_pntr[5]*x[1] + val_pntr[10]*x[2]
            + val_pntr[15]*x[3] + val_pntr[20]*x[4];
          c_pntr[1] += val_pntr[1]*x[0] + val_pntr[6]*x[1] + val_pntr[11]*x[2]
            + val_pntr[16]*x[3] + val_pntr[21]*x[4];
          c_pntr[2] += val_pntr[2]*x[0] + val_pntr[7]*x[1] + val_pntr[12]*x[2]
            + val_pntr[17]*x[3] + val_pntr[22]*x[4];
          c_pntr[3] += val_pntr[3]*x[0] + val_pntr[8]*x[1] + val_pntr[13]*x[2]
            + val_pntr[18]*x[3] + val_pntr[23]*x[4];
          c_pntr[4] += val_pntr[4]*x[0] + val_pntr[9]*x[1] + val_pntr[14]*x[2]
            + val_pntr[19]*x[3] + val_pntr[24]*x[4];
          break;

        case 6:
          c_pntr[0] += val_pntr[0]*x[0] + val_pntr[6] *x[1] + val_pntr[12]*x[2]
            + val_pntr[18]*x[3] + val_pntr[24]*x[4] + val_pntr[30]*x[5];
          c_pntr[1] += val_pntr[1]*x[0] + val_pntr[7] *x[1] + val_pntr[13]*x[2]
            + val_pntr[19]*x[3] + val_pntr[25]*x[4] + val_pntr[31]*x[5];
          c_pntr[2] += val_pntr[2]*x[0] + val_pntr[8] *x[1] + val_pntr[14]*x[2]
            + val_pntr[20]*x[3] + val_pntr[26]*x[4] + val_pntr[32]*x[5];
          c_pntr[3] += val_pntr[3]*x[0] + val_pntr[9] *x[1] + val_pntr[15]*x[2]
            + val_pntr[21]*x[3] + val_pntr[27]*x[4] + val_pntr[33]*x[5];
          c_pntr[4] += val_pntr[4]*x[0] + val_pntr[10]*x[1] + val_pntr[16]*x[2]
            + val_pntr[22]*x[3] + val_pntr[28]*x[4] + val_pntr[34]*x[5];
          c_pntr[5] += val_pntr[5]*x[0] + val_pntr[11]*x[1] + val_pntr[17]*x[2]
            + val_pntr[23]*x[3] + val_pntr[29]*x[4] + val_pntr[35]*x[5];
          break;

        default:

          /*
           * For most computers, a really well-optimized assembly-coded level 2
           * blas for small blocks sizes doesn't exist.  It's better to
           * optimize your own version, and take out all the overhead from the
           * regular dgemv call.  For large block sizes, it's also a win to
           * check for a column of zeroes; this is what dgemv_ does.  The
           * routine dgemvnsqr_() is a fortran routine that contains optimized
           * code for the hp, created from the optimizing preprocessor. Every
           * workstation will probably have an entry here eventually, since
           * this is a key optimization location.
           */

#ifdef hp
          dgemvnsqr_(&m1, val_pntr, x, c_pntr);
#else
          if (m1 < 10)
            AZ_dgemv2(m1, n1, val_pntr, x, c_pntr);
          else
            dgemv_(N, &m1, &n1, &one, val_pntr, &m1, x, &ione, &one, c_pntr,
                   &ione, strlen(N));
#endif

        }
      }

      /* nonsquare cases */

      else {
        if (m1 < 10)
          AZ_dgemv2(m1, n1, val_pntr, x, c_pntr);
        else
          dgemv_(N, &m1, &n1, &one, val_pntr, &m1, x, &ione, &one, c_pntr,
                 &ione, strlen(N));
      }

      val_pntr += iblk_size;
    }
  }

} /* dvbr_sparax_basic */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void AZ_matvec_mult(double *val, int *indx, int *bindx, int *rpntr, int *cpntr,
                    int *bpntr, double *b, register double *c,
                    int exchange_flag, int *data_org)

/******************************************************************************

  c = Ab:
  Sparse (square) overlapped matrix-vector multiply, using the distributed
  variable block row (DVBR) data structure (A = val).

  Author:          Scott A. Hutchinson, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  m:               Number of (block) rows in A.

  val:             Array containing the entries of the matrix. The matrix is
                   stored block-row-by-block-row. Each block entry is dense and
                   stored by columns (VBR).

  indx:            The ith element of indx points to the location in val of the
                   (0,0) entry of the ith block entry. The last element is the
                   number of nonzero entries of matrix A plus one.

  bindx:           Contains the block column indices of the non-zero block
                   entries.

  rpntr:           The ith element of rpntr indicates the first point row in
                   the i'th block row. The last element is the number of block
                   rows plus one.

  cpntr:           The jth element of cpntr indicates the first point column in
                   the jth block column. The last element is the number of
                   block columns plus one.

  bpntr:           The ith element of bpntr points to the first block entry of
                   the ith row in bindx. The last element is the number of
                   nonzero blocks of matrix A plus one.

  b:               Contains the vector b.

  c:               Contains the result vector b.

  exchange_flag:   Flag which controls call to exchange_bdry.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see file params.txt).

******************************************************************************/

{

  /* local variables */

  int          num_blks;
  register int j, k, irow, bindx_row;
  int          N, nzeros;

  /**************************** execution begins *****************************/

  if (data_org[AZ_matrix_type] == AZ_VBR_MATRIX) {
    num_blks = data_org[AZ_N_int_blk];

    /* multiple processors */

    if (exchange_flag )
      num_blks += data_org[AZ_N_bord_blk];

    /*
     * It is possible to overlap the communication and computation here for
     * performace gains.  The idea is to gather the messages each processor
     * needs to send (AZ_gather_mesg_info), write these out (send them),
     * perform a 'dvbr_sparax_basic' on the INTERNAL portion of the sparse
     * matrix-vector product, read the messages from the neighboring processors
     * and then perform a 'dvbr_sparax_basic' on the BOUNDARY portion of the
     * product.  We do not support this capability at this time.  SAH, 2/1996
     */

    /* perform the sparax - NOTE: the boundary exchange is done inside
       dvbr_sparax_basic */

    dvbr_sparax_basic(num_blks, val, indx, bindx, rpntr, cpntr, bpntr, b, c,
                      exchange_flag, data_org);
  }

  /* MSR version */

  else if (data_org[AZ_matrix_type] == AZ_MSR_MATRIX) {
    N = data_org[AZ_N_internal] + data_org[AZ_N_border];

    /* exchange boundary info */

    AZ_exchange_bdry(b, data_org);

    for (irow = 0; irow < N; irow++) {

      /* compute diagonal contribution */

      *c = val[irow] * b[irow];

      /* nonzero off diagonal contribution */

      bindx_row = bindx[irow];
      nzeros    = bindx[irow+1] - bindx_row;

      for (j = 0; j < nzeros; j++) {
        k   = bindx_row + j;
        *c += val[k] * b[bindx[k]];
      }
      c++;
    }
  }

} /* AZ_matvec_mult */
