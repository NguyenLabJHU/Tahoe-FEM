/* $Id: preordxc.c,v 1.1 2004-12-28 18:50:41 paklein Exp $ */
/* preordxc.f -- translated by f2c (version 20030320).
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
static integer c__1 = 1;
static integer c__11 = 11;

/* /+***************************************************************************+/ */
/* /+                                                                           +/ */
/* /+   (C) Copyright IBM Corporation, 1997                                     +/ */
/* /+   (C) Copyright Regents of the University of Minnesota, 1997              +/ */
/* /+                                                                           +/ */
/* /+   preordxc.f                                                              +/ */
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
/* /+ $Id: preordxc.c,v 1.1 2004-12-28 18:50:41 paklein Exp $ +/ */
/* /+***************************************************************************+/ */
/*<    >*/
static integer lbit_shift(integer a, integer b) {
	return b >= 0 ? a << b : (integer)((uinteger)a >> -b);
};

/* Subroutine */ int preordxc_(integer *n, integer *iown, integer *rowdist, 
	integer *mynnodes, integer *isc, integer *isd, integer *irc, integer *
	ird, integer *nown, integer *slocx, integer *rlocx, integer *siprm, 
	integer *dd, integer *myid, integer *whichproc, MPI_Comm *comm)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer lbit_shift(integer, integer);

    /* Local variables */
    integer i__, j, k, pp, col;

/*<       implicit none >*/
/*<       include 'mpif.h' >*/
/*<       integer N,rowdist(0:*),nown,isc(0:*),isd(0:*) >*/
/* -*- fortran -*- */

/* double precision functions */
/*<       integer irc(0:*),ird(0:*),mynnodes,dd,comm >*/
/*<       integer whichproc(0:*),iown(0:*),slocx(0:*),rlocx(0:*) >*/
/*<       integer i,j,k,col,pp,siprm(0:*) >*/
/*<       integer myid >*/
/*<       pp = ishft(1,dd) >*/
    pp = lbit_shift(1, *dd);
/*<       do i=0,pp-1 >*/
    i__1 = pp - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
/*<         isc(i) = 0 >*/
	isc[i__] = 0;
/*<       end do >*/
    }
/*<       do i=0,nown-1 >*/
    i__1 = *nown - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
/*< 	col = iown(2*i+1) >*/
	col = iown[(i__ << 1) + 1];
/*< 	j = 1 >*/
	j = 1;
/*<         do while(col.ge.rowdist(j)) >*/
	while(col >= rowdist[j]) {
/*< 	  j = j+1 >*/
	    ++j;
/*<         end do >*/
	}
/*< 	k = j-1 >*/
	k = j - 1;
/*< 	whichproc(i) = j-1 >*/
	whichproc[i__] = j - 1;
/*< 	isc(k) = isc(k)+1 >*/
	++isc[k];
/*<       end do >*/
    }
/*<       call mpi_alltoall(isc,1,MPI_INTEGER,irc,1,MPI_INTEGER,comm,i) >*/
/*    mpi_alltoall__(isc, &c__1, &c__11, irc, &c__1, &c__11, comm, &i__); */
	MPI_Alltoall(isc,1,MPI_INT,irc,1,MPI_INT,*comm);

/*<       isd(0) = 0 >*/
    isd[0] = 0;
/*<       ird(0) = 0 >*/
    ird[0] = 0;
/*<       do i=1,pp-1 >*/
    i__1 = pp - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*< 	isd(i) = isd(i-1)+isc(i-1) >*/
	isd[i__] = isd[i__ - 1] + isc[i__ - 1];
/*< 	ird(i) = ird(i-1)+irc(i-1) >*/
	ird[i__] = ird[i__ - 1] + irc[i__ - 1];
/*<         isc(i-1) = 0 >*/
	isc[i__ - 1] = 0;
/*<       end do >*/
    }
/*<       isc(pp-1) = 0 >*/
    isc[pp - 1] = 0;
/*<       do i=0,nown-1 >*/
    i__1 = *nown - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
/*< 	k = whichproc(i) >*/
	k = whichproc[i__];
/*<         j = isd(k)+isc(k) >*/
	j = isd[k] + isc[k];
/*< 	slocx(i) = isd(k)+isc(k) >*/
	slocx[i__] = isd[k] + isc[k];
/*< 	siprm(j) = iown(2*i+1) >*/
	siprm[j] = iown[(i__ << 1) + 1];
/*< 	isc(k) = isc(k)+1 >*/
	++isc[k];
/*<       end do >*/
    }
/*<    >*/
/*   call mpi_alltoallv(siprm,isc,isd,MPI_INTEGER,
     +                   rlocx,irc,ird,MPI_INTEGER,
     +                   comm,i) */
/*   mpi_alltoallv__(siprm, isc, isd, &c__11, rlocx, irc, ird, &c__11, comm, &i__); */
	MPI_Alltoallv(siprm,isc,isd,MPI_INT,
	              rlocx,irc,ird,MPI_INT,*comm);

/*<       end >*/
    return 0;
} /* preordxc_ */
