/* $Id: order.c,v 1.1 2004-12-11 01:44:14 paklein Exp $ */
/* order.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

/* debugging */
#undef __DO_DEBUG__
/* #define __DO_DEBUG__ 1 */

#include "mpi.h"
#include "pspases_f2c.h"

/* Table of constant values */
static integer c__3 = 3;
static integer c__1 = 1;
static integer c__9 = 9;
static integer c__0 = 0;

/* /+***************************************************************************+/ */
/* /+                                                                           +/ */
/* /+   (C) Copyright IBM Corporation, 1997                                     +/ */
/* /+   (C) Copyright Regents of the University of Minnesota, 1997              +/ */
/* /+                                                                           +/ */
/* /+   order.f                                                                 +/ */
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
/* /+ $Id: order.c,v 1.1 2004-12-11 01:44:14 paklein Exp $ +/ */
/* /+***************************************************************************+/ */
/*<    >*/
/* Subroutine */ int porder_(integer *rowdist, integer *aptrs, integer *ainds,
	 integer *order, integer *sizes, integer *myid, integer *pp, integer *
	serialorder, integer *dbgpp, MPI_Comm *comm)
/* , integer *xadj, integer *
	adjncy, integer *sorder, integer *counts, integer *mynnodes2, integer 
	*ioasize2)
*/
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer mynnodes, i__, j, k, l, m;
    extern /* Subroutine */ int  
    	parometisf_(integer *, integer *, integer *, integer *, integer *,
	     integer *, integer *, integer *);
    static integer ierr, opts[5], offdnz, ioasize;

/*<       implicit none >*/
/*<       include 'mpif.h' >*/
/*<       integer rowdist(0:*),aptrs(2,0:*),ainds(*),order(0:*),sizes(0:*) >*/
/* -*- fortran -*- */


/* double precision functions */

/*<       double precision MPI_WTIME, MPI_WTICK, PMPI_WTIME, PMPI_WTICK >*/
/*<       external MPI_WTIME, MPI_WTICK, PMPI_WTIME, PMPI_WTICK >*/
/*      integer, allocatable :: xadj(:),adjncy(:), sorder(:), counts(:) */
/*<       integer xadj, adjncy, sorder, counts >*/
/*<       integer mynnodes2, ioasize2 >*/
/*<       dimension xadj(0:mynnodes2) >*/
/*<       dimension adjncy(0:ioasize2 - mynnodes2 - 1) >*/
/*<       integer mynnodes,opts(5),i,j,k,l,m,ierr,ioasize,offdnz >*/
/*<       integer comm,myid,pp,serialorder >*/
/*<       integer dbgpp >*/
/*<       mynnodes = rowdist(myid+1)-rowdist(myid) >*/
    /* Parameter adjustments */
    --aptrs;
    --ainds;

    /* Function Body */
    mynnodes = rowdist[*myid + 1] - rowdist[*myid];
/*<       ioasize = aptrs(2,mynnodes-1)+aptrs(1,mynnodes-1)-1 >*/
    ioasize = aptrs[(mynnodes - 1 << 1) + 2] + aptrs[(mynnodes - 1 << 1) + 1] 
	    - 1;
/*      allocate(xadj(0:mynnodes),stat=i) */
	xadj = (integer*) malloc(mynnodes*sizeof(integer));

/*<       if(i.ne.0) then >*/
    if (!xadj) {
/*<         print *,myid,':memory allocation failure' >*/
		printf("%d: memory allocation failure", myid);
/*<         call mpi_abort(comm,0,ierr) >*/
		MPI_Abort(*comm, 0);
/*<       end if >*/
    }

/*      allocate(adjncy(0:ioasize-mynnodes-1),stat=i) */
	adjncy = (integer*) malloc((ioasize-mynnodes-1)*sizeof(integer));

/*<       if(i.ne.0) then >*/
    if (!adjncy) {
/*<         print *,myid,':memory allocation failure' >*/
		printf("%d: memory allocation failure", myid);
/*<         call mpi_abort(comm,0,ierr) >*/
		MPI_Abort(*comm, 0);
/*<       end if >*/
    }
/* form xadj and adjncy from aptrs, ainds */
/*<       m = rowdist(myid) >*/
    m = rowdist[*myid];
/*<       l = 0 >*/
    l = 0;
/*<       do i=0,mynnodes-1 >*/
    i__1 = mynnodes - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
/*<        xadj(i) = aptrs(1,i)-i-1 >*/
	xadj[i__] = aptrs[(i__ << 1) + 1] - i__ - 1;
/*<        do j=aptrs(1,i),aptrs(1,i)+aptrs(2,i)-1 >*/
	i__2 = aptrs[(i__ << 1) + 1] + aptrs[(i__ << 1) + 2] - 1;
	for (j = aptrs[(i__ << 1) + 1]; j <= i__2; ++j) {
/*<          k = ainds(j) >*/
	    k = ainds[j];
/*<          if(k.ne.m+i) then >*/
	    if (k != m + i__) {
/*<            adjncy(l) = k >*/
		adjncy[l] = k;
/*<            l = l+1 >*/
		++l;
/*<          end if >*/
	    }
/*<        end do >*/
	}
/*<       end do >*/
    }
/*<       xadj(mynnodes) = xadj(mynnodes-1)+aptrs(2,mynnodes-1)-1 >*/
    xadj[mynnodes] = xadj[mynnodes - 1] + aptrs[(mynnodes - 1 << 1) + 2] - 1;
/*      call mpi_allreduce(xadj(mynnodes),offdnz,1,MPI_INTEGER,MPI_SUM, */
/*     +                   comm,ierr); */
		MPI_Allreduce(xadj + mynnodes, &offdnz, 1, MPI_INTEGER, MPI_SUM, *comm);

/*<       if(offdnz.eq.0) then >*/
    if (offdnz == 0) {
/*<         do i=0,mynnodes-1 >*/
	i__1 = mynnodes - 1;
	for (i__ = 0; i__ <= i__1; ++i__) {
/*<           order(i) = i + m >*/
	    order[i__] = i__ + m;
/*<         end do >*/
	}
/*<         do i=0,pp-2 >*/
	i__1 = *pp - 2;
	for (i__ = 0; i__ <= i__1; ++i__) {
/*<           sizes(i) = rowdist(pp)/pp >*/
	    sizes[i__] = rowdist[*pp] / *pp;
/*<         end do >*/
	}
/*<         sizes(pp-1) = rowdist(pp) - (pp-1)*(rowdist(pp)/pp) >*/
	sizes[*pp - 1] = rowdist[*pp] - (*pp - 1) * (rowdist[*pp] / *pp);
/*<         do i=pp,2*pp-1 >*/
	i__1 = (*pp << 1) - 1;
	for (i__ = *pp; i__ <= i__1; ++i__) {
/*<           sizes(i) = 0 >*/
	    sizes[i__] = 0;
/*<         end do >*/
	}
/*<       else >*/
    } else {
/*<         opts(1) = 1 >*/
	opts[0] = 1;
/*<         opts(2) = 2 >*/
	opts[1] = 2;
/*<         opts(3) = dbgpp >*/
	opts[2] = *dbgpp;
/*<         opts(4) = 0 >*/
	opts[3] = 0;
/*<    >*/
	parometisf_(rowdist, xadj, adjncy, order, sizes, opts, serialorder, comm);
/*<       end if >*/
    }

/*      deallocate(xadj) */
/*      deallocate(adjncy) */
	free(xadj);
	free(adjncy);

/*<       end >*/
    return 0;
} /* porder_ */

