/* $Id: alltoallu_hc.c,v 1.1 2005-01-03 00:08:11 paklein Exp $ */
/* alltoallu_hc.f -- translated by f2c (version 20030320).
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
static integer c__1 = 1;

/* /+***************************************************************************+/ */
/* /+                                                                           +/ */
/* /+   (C) Copyright IBM Corporation, 1997                                     +/ */
/* /+   (C) Copyright Regents of the University of Minnesota, 1997              +/ */
/* /+                                                                           +/ */
/* /+   alltoallu_hc.f                                                          +/ */
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
/* /+ $Id: alltoallu_hc.c,v 1.1 2005-01-03 00:08:11 paklein Exp $ +/ */
/* /+***************************************************************************+/ */

static integer lbit_shift(integer a, integer b) {
	return b >= 0 ? a << b : (integer)((uinteger)a >> -b);
};

/*<    >*/
/* Subroutine */ int all_to_all_union_hc_(integer *src, integer *tgt, 
	integer *tmp, integer *mxsize, integer *srcsize, integer *tgtsize, 
	integer *myid, integer *dim, MPI_Comm *comm)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    integer ipartner, i__, j;
    extern /* Subroutine */ int mergelists_(integer *, integer *, integer *, 
	    integer *, integer *, integer *);
    integer ierr;
    integer /* mpistat[4], */ tmpsize;
	MPI_Status mpistat[4];

/*<       implicit none >*/
/*<       include 'mpif.h' >*/
/*<       integer itag >*/
/* -*- fortran -*- */

/* double precision functions */
/*<       parameter(itag=1) >*/
/*<       integer srcsize,tgtsize,mxsize,dim,myid,ipartner,comm >*/
/*<       integer src(*), tgt(*), tmp(*) >*/
/*<       integer i,j,tmpsize,nbrecv >*/
/*<       integer ierr,mpistat(MPI_STATUS_SIZE) >*/
/*<       do i=1,srcsize >*/
    /* Parameter adjustments */
    --tmp;
    --tgt;
    --src;

    /* Function Body */
    i__1 = *srcsize;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         tgt(i) = src(i) >*/
	tgt[i__] = src[i__];
/*<       end do >*/
    }
/*<       tgtsize = srcsize >*/
    *tgtsize = *srcsize;
/*<       do i=0,dim-1 >*/
    i__1 = *dim - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
/*<         ipartner = ieor(myid,ishft(1,i)) >*/
	ipartner = *myid ^ lbit_shift((ftnlen)1, i__);
/*<         if(ipartner.gt.myid) then >*/
	if (ipartner > *myid) {
/*<    >*/
	    MPI_Send(&tgt[1], *tgtsize, MPI_INT, ipartner, 1, *comm);
/*<    >*/
	    MPI_Recv(&tmp[1], *mxsize, MPI_INT, ipartner, 1, *comm, mpistat);
/*<         else >*/
	} else {
/*<    >*/
	    MPI_Recv(&tmp[1], *mxsize, MPI_INT, ipartner, 1, *comm, mpistat);
/*<    >*/
	    MPI_Send(&tgt[1], *tgtsize, MPI_INT, ipartner, 1, *comm);
/*<         end if >*/
	}
/*<         call mpi_get_count(mpistat,MPI_INTEGER,tmpsize,ierr) >*/
	MPI_Get_count__(mpistat, MPI_INT, &tmpsize);
/*<         call mergelists(tmp,tmpsize,tgt,tgtsize,src,srcsize) >*/
	mergelists_(&tmp[1], &tmpsize, &tgt[1], tgtsize, &src[1], srcsize);
/*<         do j=1,srcsize >*/
	i__2 = *srcsize;
	for (j = 1; j <= i__2; ++j) {
/*<           tgt(j) = src(j) >*/
	    tgt[j] = src[j];
/*<         end do >*/
	}
/*<         tgtsize = srcsize >*/
	*tgtsize = *srcsize;
/*<       end do >*/
    }
/*<       end >*/
    return 0;
} /* all_to_all_union_hc__ */

/*<       subroutine mergelists(as1,len1,as2,len2,at,lent) >*/
/* Subroutine */ int mergelists_(integer *as1, integer *len1, integer *as2, 
	integer *len2, integer *at, integer *lent)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__, j, k;

/*<       implicit none >*/
/*<       integer len1,len2,lent >*/
/*<       integer as1(*), as2(*), at(*) >*/
/*<       integer i,j,k >*/
/*<       i=1 >*/
    /* Parameter adjustments */
    --at;
    --as2;
    --as1;

    /* Function Body */
    i__ = 1;
/*<       j=1 >*/
    j = 1;
/*<       k=1 >*/
    k = 1;
/*<       do while( (i.le.len1) .and. (j.le.len2) ) >*/
    while(i__ <= *len1 && j <= *len2) {
/*<         if(as1(i) .lt. as2(j)) then >*/
	if (as1[i__] < as2[j]) {
/*<           at(k) = as1(i) >*/
	    at[k] = as1[i__];
/*<           i = i + 1 >*/
	    ++i__;
/*<         else >*/
	} else {
/*<           if(as1(i) .eq. as2(j)) then >*/
	    if (as1[i__] == as2[j]) {
/*<             at(k) = as1(i) >*/
		at[k] = as1[i__];
/*<             i = i + 1 >*/
		++i__;
/*<             j = j + 1 >*/
		++j;
/*<           else >*/
	    } else {
/*<             at(k) = as2(j) >*/
		at[k] = as2[j];
/*<             j = j + 1 >*/
		++j;
/*<           end if >*/
	    }
/*<         end if >*/
	}
/*<         k = k + 1 >*/
	++k;
/*<       end do >*/
    }
/*<       do i=i,len1 >*/
    i__1 = *len1;
    for (i__ = i__; i__ <= i__1; ++i__) {
/*<         at(k) = as1(i) >*/
	at[k] = as1[i__];
/*<         k = k + 1 >*/
	++k;
/*<       end do >*/
    }
/*<       do j=j,len2 >*/
    i__1 = *len2;
    for (j = j; j <= i__1; ++j) {
/*<         at(k) = as2(j) >*/
	at[k] = as2[j];
/*<         k = k + 1 >*/
	++k;
/*<       end do >*/
    }
/*<       lent = k-1 >*/
    *lent = k - 1;
/*<       end >*/
    return 0;
} /* mergelists_ */

