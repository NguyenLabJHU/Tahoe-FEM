/* reordx.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

/* debugging */
#undef __DO_DEBUG__
/* #define __DO_DEBUG__ 1 */

#include "mpi.h"
#include "pspases_f2c.h"
#include <stdio.h>

/* Table of constant values */

static integer c__13 = 13;

/* /+***************************************************************************+/ */
/* /+                                                                           +/ */
/* /+   (C) Copyright IBM Corporation, 1997                                     +/ */
/* /+   (C) Copyright Regents of the University of Minnesota, 1997              +/ */
/* /+                                                                           +/ */
/* /+   reordx.f                                                                +/ */
/* /+                                                                           +/ */
/* /+   Written by Mahesh Joshi, U of MN.                                       +/ */
/* /+                                                                           +/ */
/* /+***************************************************************************+/ */
/* /+                                                                           +/ */
/* /+ This code is meant to be used solely for educational, research, and       +/ */
/* /+ benchmarking purposes by non-profit institutions and US government        +/ */
/* /+ agencies only.  Use by any other organization requires prior written      +/ */
/* /+ permission from both IBM Corporation and the University of Minnesota.     +/ */
/* /+ The software may not be sold or redistributed.  One may make copies       +/ */
/* /+ of the software or modify it for their use provided that the copies,      +/ */
/* /+ modified or otherwise, are not sold or distributed, are used under the    +/ */
/* /+ same terms and conditions, and this notice and any part of the source     +/ */
/* /+ code that follows this notice are not separated.                          +/ */
/* /+                                                                           +/ */
/* /+ As unestablished research software, this code is provided on an           +/ */
/* /+ ``as is'' basis without warranty of any kind, either expressed or         +/ */
/* /+ implied, including but not limited to implied warranties of               +/ */
/* /+ merchantability and fitness for a particular purpose.  IBM does not       +/ */
/* /+ warrant that the functions contained in this software will meet the       +/ */
/* /+ user's requirements or that the operation of its routines will be         +/ */
/* /+ uninterrupted or error-free.  Acceptance and use of this program          +/ */
/* /+ constitutes the user's understanding that he/she will have no recourse    +/ */
/* /+ to IBM for any actual or consequential damages, including, but not        +/ */
/* /+ limited to, lost profits or savings, arising out of the use or inability  +/ */
/* /+ to use these libraries.  Even if the user informs IBM of the possibility  +/ */
/* /+ of such damages, IBM expects the user to accept the risk of any such      +/ */
/* /+ harm, or the user shall not attempt to use these libraries for any        +/ */
/* /+ purpose.                                                                  +/ */
/* /+                                                                           +/ */
/* /+ The downloading, compiling, or executing any part of this software        +/ */
/* /+ constitutes an implicit agreement to these terms.  These terms and        +/ */
/* /+ conditions are subject to change at any time without prior notice.        +/ */
/* /+                                                                           +/ */
/* /+***************************************************************************+/ */
/* /+ $Id: reordx.c,v 1.1 2004-12-28 18:11:12 paklein Exp $ +/ */
/* /+***************************************************************************+/ */
/*<    >*/
/* Subroutine */ int reordx_(integer *n, doublereal *x, integer *ldx, integer 
	*nrhs, doublereal *rx, integer *nown, integer *dsc, integer *dsd, 
	integer *drc, integer *drd, integer *iown, integer *rowdist, integer *
	mynnodes, integer *sloc, integer *rloc, doublereal *svals, doublereal 
	*rvals, integer *myid, MPI_Comm *comm)
{
    /* System generated locals */
    integer x_dim1, x_offset, rx_dim1, rx_offset, i__1, i__2;

    /* Local variables */
    integer i__, j;

/*<       implicit none >*/
/*<       include 'mpif.h' >*/
/*<       integer N,nown,dsc(0:*),dsd(0:*),ldx >*/
/* -*- fortran -*- */

/* double precision functions */

/*<       integer drc(0:*),drd(0:*),mynnodes >*/
/*<       integer sloc(0:*),iown(0:*),rloc(0:*),rowdist(0:*) >*/
/*<       integer i,j,nrhs,comm >*/
/*<       double precision x(0:ldx-1,0:*),rx(0:N-1,0:*) >*/
/*<       double precision svals(0:*),rvals(0:*) >*/
/*<       integer myid >*/
/*<       do j=0,nrhs-1 >*/
    /* Parameter adjustments */
    rx_dim1 = *n - 1 - 0 + 1;
    rx_offset = 0 + rx_dim1 * 0;
    rx -= rx_offset;
    x_dim1 = *ldx - 1 - 0 + 1;
    x_offset = 0 + x_dim1 * 0;
    x -= x_offset;

    /* Function Body */
    i__1 = *nrhs - 1;
    for (j = 0; j <= i__1; ++j) {
/*<         do i=0,nown-1 >*/
	i__2 = *nown - 1;
	for (i__ = 0; i__ <= i__2; ++i__) {
/*<           svals(sloc(i)*nrhs+j) = rx(iown(2*i),j) >*/
	    svals[sloc[i__] * *nrhs + j] = rx[iown[i__ * 2] + j * rx_dim1];
/*<         end do >*/
	}
/*<       end do >*/
    }
/*<    >*/
/*    call mpi_alltoallv(svals,dsc,dsd,MPI_DOUBLE_PRECISION,
     +                   rvals,drc,drd,MPI_DOUBLE_PRECISION,
     +                   comm,i) */
/*  mpi_alltoallv__(svals, dsc, dsd, &c__13, rvals, drc, drd, &c__13, comm, &i__); */
	MPI_Alltoallv(svals, dsc, dsd, MPI_DOUBLE, 
	              rvals, drc, drd, MPI_DOUBLE, *comm);

/*<       do j=0,nrhs-1 >*/
    i__1 = *nrhs - 1;
    for (j = 0; j <= i__1; ++j) {
/*<         do i=0,mynnodes-1 >*/
	i__2 = *mynnodes - 1;
	for (i__ = 0; i__ <= i__2; ++i__) {
/*<           x(rloc(i)-rowdist(myid),j) = rvals(i*nrhs+j) >*/
	    x[rloc[i__] - rowdist[*myid] + j * x_dim1] = rvals[i__ * *nrhs + 
		    j];
/*<         end do >*/
	}
/*<       end do >*/
    }
/*<       end >*/
    return 0;
} /* reordx_ */

