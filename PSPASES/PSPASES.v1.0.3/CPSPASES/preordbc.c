/* $Id: preordbc.c,v 1.1 2004-12-29 06:03:48 paklein Exp $ */
/* preordbc.f -- translated by f2c (version 20030320).
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
static integer c__11 = 11;

/* /+***************************************************************************+/ */
/* /+                                                                           +/ */
/* /+   (C) Copyright IBM Corporation, 1997                                     +/ */
/* /+   (C) Copyright Regents of the University of Minnesota, 1997              +/ */
/* /+                                                                           +/ */
/* /+   preordbc.f                                                              +/ */
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
/* /+ $Id: preordbc.c,v 1.1 2004-12-29 06:03:48 paklein Exp $ +/ */
/* /+***************************************************************************+/ */

static integer lbit_shift(integer a, integer b) {
	return b >= 0 ? a << b : (integer)((uinteger)a >> -b);
};

/*<    >*/
/* Subroutine */ int preordbc_(integer *n, integer *order, integer *ranmasks, 
	integer *nown, integer *isc, integer *isd, integer *irc, integer *ird,
	 integer *iown, integer *rowdist, integer *mynnodes, integer *sloc, 
	integer *siprm, integer *dd, integer *myid, integer *lgblk, integer *
	pgrid, integer *whichsnode, MPI_Comm *comm)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__, j, k, kc, kr, pp, col, ppc, ppg, ppr;
    integer pgrsize;

/*<       implicit none >*/
/*<       include 'mpif.h' >*/
/*<       integer N,order(0:*),ranmasks(5,0:*),nown,isc(0:*),isd(0:*) >*/
/* -*- fortran -*- */

/* double precision functions */
/*<       integer irc(0:*),ird(0:*),mynnodes,dd,lgblk,comm >*/
/*<       integer pgrid(0:*),whichsnode(0:*),rowdist(0:*) >*/
/*<       integer sloc(0:*),siprm(0:*),iown(0:*) >*/
/*<       integer pgrsize,ppr,ppc,ppg,kr,kc,i,j,k,col,pp >*/
/*<       integer myid >*/
/*<       pp = ishft(1,dd) >*/
    /* Parameter adjustments */
    --ranmasks;

    /* Function Body */
    pp = lbit_shift(1, *dd);
/*<       pgrsize = ishft(1,ishft(dd,-1)) >*/
    pgrsize = lbit_shift(1, lbit_shift(*dd, -1));
/*<       ppr = 0 >*/
    ppr = 0;
/*<       ppc = pp >*/
    ppc = pp;
/*<       ppg = 2*pp >*/
    ppg = pp << 1;
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
/*<         j = whichsnode(i) >*/
	j = whichsnode[i__];
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
/*<         sloc(i) = ishft((isd(k)+isc(k)),-1) >*/
	sloc[i__] = lbit_shift(isd[k] + isc[k], -1);
/*<         siprm(isd(k)+isc(k)) = col >*/
	siprm[isd[k] + isc[k]] = col;
/*<         siprm(isd(k)+isc(k)+1) = i+rowdist(myid) >*/
	siprm[isd[k] + isc[k] + 1] = i__ + rowdist[*myid];
/*<         isc(k) = isc(k)+2 >*/
	isc[k] += 2;
/*<       end do >*/
    }
/*<    >*/

/*   call mpi_alltoallv(siprm,isc,isd,MPI_INTEGER,
     +                   iown,irc,ird,MPI_INTEGER,
     +                   comm,i) */
/*  mpi_alltoallv__(siprm, isc, isd, &c__11, iown, irc, ird, &c__11, comm, &i__); */
	MPI_Alltoallv(siprm,isc,isd,MPI_INT,
	               iown,irc,ird,MPI_INT,*comm);

/*<       end >*/
    return 0;
} /* preordbc_ */
