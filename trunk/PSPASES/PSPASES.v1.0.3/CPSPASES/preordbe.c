/* $Id: preordbe.c,v 1.1 2004-12-29 05:57:17 paklein Exp $ */
/* preordbe.f -- translated by f2c (version 20030320).
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
/* /+   preordbe.f                                                              +/ */
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
/* /+ $Id: preordbe.c,v 1.1 2004-12-29 05:57:17 paklein Exp $ +/ */
/* /+***************************************************************************+/ */

static integer lbit_shift(integer a, integer b) {
	return b >= 0 ? a << b : (integer)((uinteger)a >> -b);
};
/*<    >*/
/* Subroutine */ int preordbe_(integer *n, integer *order, integer *ranmasks, 
	integer *nown, integer *isc, integer *isd, integer *irc, integer *ird,
	 integer *mynnodes, integer *dd, integer *myid, integer *lgblk, 
	integer *pgrid, integer *whichsnode, MPI_Comm *comm)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    integer i__, j, k, kc, kr, pp, col, ppc, ppg, ppr;
    integer cbits, rbits, pgrsize;

/*<       implicit none >*/
/*<       include 'mpif.h' >*/
/*<       integer N,order(0:*),ranmasks(5,0:*),nown,isc(0:*),isd(0:*) >*/
/* -*- fortran -*- */

/* Copyright (c) 2001-2002 The Trustees of Indiana University. */
/*                         All rights reserved. */
/* Copyright (c) 1998-2001 University of Notre Dame. */
/*                         All rights reserved. */
/* Copyright (c) 1994-1998 The Ohio State University. */
/*                         All rights reserved. */

/* double precision functions */
/*<       integer irc(0:*),ird(0:*),mynnodes,dd,lgblk,comm >*/
/*<       integer pgrid(0:*),whichsnode(0:*) >*/
/*<       integer rbits,cbits,pgrsize,ppr,ppc,ppg,kr,kc,i,j,k,col,pp >*/
/*<       integer myid >*/
/*<       pp = ishft(1,dd) >*/
    /* Parameter adjustments */
    --ranmasks;

    /* Function Body */
    pp = lbit_shift(1, *dd);
/*<       rbits = ishft(dd,-1) >*/
    rbits = lbit_shift(*dd, -1);
/*<       cbits = dd-rbits >*/
    cbits = *dd - rbits;
/*<       pgrsize = ishft(1,rbits) >*/
    pgrsize = lbit_shift(1, rbits);
/*<       ppr = 0 >*/
    ppr = 0;
/*<       ppc = pp >*/
    ppc = pp;
/*<       ppg = 2*pp >*/
    ppg = pp << 1;
/*<       do k=0,pp-1 >*/
    i__1 = pp - 1;
    for (k = 0; k <= i__1; ++k) {
/*<         kr = 0 >*/
	kr = 0;
/*<         i = ishft(k,-1) >*/
	i__ = lbit_shift(k, -1);
/*<         do j=0,rbits-1 >*/
	i__2 = rbits - 1;
	for (j = 0; j <= i__2; ++j) {
/*<           kr = ior(kr,ishft(iand(i,1),j)) >*/
	    kr |= lbit_shift(i__ & 1, j);
/*<           i = ishft(i,-2) >*/
	    i__ = lbit_shift(i__, -2);
/*<         end do >*/
	}
/*<         kc = 0 >*/
	kc = 0;
/*<         i = k >*/
	i__ = k;
/*<         do j=0,cbits-1 >*/
	i__2 = cbits - 1;
	for (j = 0; j <= i__2; ++j) {
/*<           kc = ior(kc,ishft(iand(i,1),j)) >*/
	    kc |= lbit_shift(i__ & 1, j);
/*<           i = ishft(i,-2) >*/
	    i__ = lbit_shift(i__, -2);
/*<         end do >*/
	}
/*<         pgrid(ppr+k) = kr >*/
	pgrid[ppr + k] = kr;
/*<         pgrid(ppc+k) = kc >*/
	pgrid[ppc + k] = kc;
/*<         pgrid(ppg+kc*pgrsize+kr) = k >*/
	pgrid[ppg + kc * pgrsize + kr] = k;
/*<       end do >*/
    }
/*<       do i=0,pp-1 >*/
    i__1 = pp - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
/*<         isc(i) = 0 >*/
	isc[i__] = 0;
/*<       end do >*/
    }
/*<       do i=0,mynnodes-1 >*/
    i__1 = *mynnodes - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
/*<         col = order(i) >*/
	col = order[i__];
/*<         j = 0 >*/
	j = 0;
/*<         do while (col.lt.ranmasks(1,j) .or. col.gt.ranmasks(2,j)) >*/
	while(col < ranmasks[j * 5 + 1] || col > ranmasks[j * 5 + 2]) {
/*<           j = j+1 >*/
	    ++j;
/*<         end do >*/
	}
/*<         whichsnode(i) = j >*/
	whichsnode[i__] = j;
/*<         k = ranmasks(3,j) >*/
	k = ranmasks[j * 5 + 3];
/*<         kr = pgrid(ppr+k)+iand(ishft(col,-lgblk),ranmasks(4,j)) >*/
	kr = pgrid[ppr + k] + (lbit_shift(col, -(*lgblk)) & ranmasks[j * 5 + 
		4]);
/*<         kc = pgrid(ppc+k)+iand(ishft(col,-lgblk),ranmasks(5,j)) >*/
	kc = pgrid[ppc + k] + (lbit_shift(col, -(*lgblk)) & ranmasks[j * 5 + 
		5]);
/*<         k = pgrid(ppg+kc*pgrsize+kr) >*/
	k = pgrid[ppg + kc * pgrsize + kr];
/*<         isc(k) = isc(k)+2 >*/
	isc[k] += 2;
/*<       end do >*/
    }
/*<       isd(0) = 0 >*/
    isd[0] = 0;
/*<       do i=1,pp-1 >*/
    i__1 = pp - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         isd(i) = isd(i-1)+isc(i-1) >*/
	isd[i__] = isd[i__ - 1] + isc[i__ - 1];
/*<       end do >*/
    }
/*<       call mpi_alltoall(isc,1,MPI_INTEGER,irc,1,MPI_INTEGER,comm,i) >*/
/*    mpi_alltoall__(isc, &c__1, &c__11, irc, &c__1, &c__11, comm, &i__); */
	MPI_Alltoall(isc,1,MPI_INT,irc,1,MPI_INT,*comm);

/*<       ird(0) = 0 >*/
    ird[0] = 0;
/*<       nown = 0 >*/
    *nown = 0;
/*<       do i=1,pp-1 >*/
    i__1 = pp - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         ird(i) = ird(i-1)+irc(i-1) >*/
	ird[i__] = ird[i__ - 1] + irc[i__ - 1];
/*<         nown = nown+irc(i-1) >*/
	*nown += irc[i__ - 1];
/*<       end do >*/
    }
/*<       nown = nown+irc(pp-1) >*/
    *nown += irc[pp - 1];
/*<       nown = ishft(nown,-1) >*/
    *nown = lbit_shift(*nown, -1);
/*<       end >*/
    return 0;
} /* preordbe_ */

