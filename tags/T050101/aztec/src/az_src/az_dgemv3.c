/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile: az_dgemv3.c,v $
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
static char rcsid[] = "$Id: az_dgemv3.c,v 1.1.1.1 2001-01-30 20:59:12 paklein Exp $";
#endif


/*******************************************************************************
 * Copyright 1995, Sandia Corporation.  The United States Government retains a *
 * nonexclusive license in this software as prescribed in AL 88-1 and AL 91-7. *
 * Export of this program may require a license from the United States         *
 * Government.                                                                 *
 ******************************************************************************/

#include "az_aztec.h" /* missing include - PAK (07/24/98) */

void AZ_dgemv3(int m, int n, double *a, double *x, double *y)

/*******************************************************************************

  Perform the matrix-vector operation y = y - Ax where x and y are vectors and
  A is an m-by-n matrix.  This modeled on the blas 2 routine dgemv and used for
  small matrix sizes.  Also, there is a special case for m = n = 5, (3d fluids).

  Author:          Lydie Prevost, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  m:               Number of rows of the matrix a.

  n:               Number of columns of the matrix a.

  a:               The leading m-by-n part of the array a must contain the
                   matrix of coefficients.

  x:               The vector of length at least n to multiply by A.

  y:               The vector of length at least m which, on output, will
                   contain the result of the -Ax operation summed into it.
                   NOTE:  the result is the SUM of the entry y and -Ax.

*******************************************************************************/

{

  /* local variables */

  register int i = 0;

  /**************************** execution begins ******************************/

  switch (m) {
  case 5:
    while (i++ < n) {
      *y     -= *a++ * *x;
      *(y+1) -= *a++ * *x;
      *(y+2) -= *a++ * *x;
      *(y+3) -= *a++ * *x;
      *(y+4) -= *a++ * *x;

      x++;

    }

    break;

  default:

     while (i++ < n) {
      register int j = 0;
      while (j++ < m)
        *y++ -= *a++ * *x;
      y -= m; x++;
    }
  }

} /* AZ_dgemv3 */
