/* $Id: alltoonev.c,v 1.2 2005-01-04 18:19:34 paklein Exp $ */
/* alltoonev.f -- translated by f2c (version 20030320).
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
static integer c__5 = 5;
static integer c__1 = 1;

/* /+***************************************************************************+/ */
/* /+                                                                           +/ */
/* /+   (C) Copyright IBM Corporation, 1997                                     +/ */
/* /+   (C) Copyright Regents of the University of Minnesota, 1997              +/ */
/* /+                                                                           +/ */
/* /+   alltoonev.f                                                             +/ */
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
/* /+ $Id: alltoonev.c,v 1.2 2005-01-04 18:19:34 paklein Exp $ +/ */
/* /+***************************************************************************+/ */

static integer lbit_shift(integer a, integer b) {
	return b >= 0 ? a << b : (integer)((uinteger)a >> -b);
};

/*<       subroutine alltoonev(tgt,tmp,size,psrc,lcsize,myidv,myid,comm) >*/
/* Subroutine */ int alltoonev_(doublereal *tgt, doublereal *tmp, integer *
	size, integer *psrc, integer *lcsize, integer *myidv, integer *myid, 
	MPI_Comm *comm)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    integer diffbits, i__, j, k, num, ierr;
    integer msglen, partner /*, mpistat[4] */;
	MPI_Status mpistat;

/*<       implicit none >*/
/*<       include 'mpif.h' >*/
/*<       integer lendp,itype >*/
/* -*- fortran -*- */

/* double precision functions */
/*<       parameter(lendp=8,itype=1) >*/
/*<       integer  size,psrc,lcsize,myid,myidv,partner,num,comm >*/
/*<       double precision tgt(size),tmp(size) >*/
/*<       integer diffbits,i,j,k,nbrecv,msglen,mids,midr >*/
/*<       integer mpistat(MPI_STATUS_SIZE),ierr >*/
/*<       diffbits = ieor(psrc,myidv) >*/
    /* Parameter adjustments */
    --tmp;
    --tgt;

    /* Function Body */
    diffbits = *psrc ^ *myidv;
/*<       num = lcsize >*/
    num = *lcsize;
/*<       if(diffbits.ne.0) then >*/
    if (diffbits != 0) {
/* the position of mostsignificant 1 in diffbits */
/* hints at how many stages of communication a */
/* processor will be involved in.... */
/*<         num = lcsize+1 >*/
	num = *lcsize + 1;
/*<         do while(diffbits.ne.0) >*/
	while(diffbits != 0) {
/*<           diffbits = ishft(diffbits,-1) >*/
	    diffbits = lbit_shift(diffbits, (ftnlen)-1);
/*<           num = num-1 >*/
	    --num;
/*<         end do >*/
	}
/*<       end if >*/
    }
/*<       msglen = size*lendp >*/
    msglen = *size << 3;
/*<       k = ishft(lcsize,1)-1 >*/
    k = (*lcsize << 1) - 1;
/*<       do i =1,num >*/
    i__1 = num;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         partner = ieor(myid,ishft(1,k)) >*/
	partner = *myid ^ lbit_shift((ftnlen)1, k);
/*<    >*/

/*        call mpi_sendrecv(tgt,msglen,MPI_BYTE,partner,itype,
     +                    tmp,msglen,MPI_BYTE,partner,itype,
     +                    comm,mpistat,ierr) */
/*	mpi_sendrecv__(&tgt[1], &msglen, &c__5, &partner, &c__1, &tmp[1], &
		msglen, &c__5, &partner, &c__1, comm, mpistat, &ierr); */

	MPI_Sendrecv(&tgt[1],msglen,MPI_BYTE,partner,1,
	             &tmp[1],msglen,MPI_BYTE,partner,1,*comm, &mpistat);

/*<         do j=1,size >*/
	i__2 = *size;
	for (j = 1; j <= i__2; ++j) {
/*<           tgt(j) = tgt(j) + tmp(j) >*/
	    tgt[j] += tmp[j];
/*<         end do >*/
	}
/*<         k = k-2 >*/
	k += -2;
/*<       end do >*/
    }
/*<       end >*/
    return 0;
} /* alltoonev_ */

