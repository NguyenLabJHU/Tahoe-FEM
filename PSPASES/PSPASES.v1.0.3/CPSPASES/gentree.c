/* $Id: gentree.c,v 1.1 2005-01-03 00:17:43 paklein Exp $ */
/* gentree.f -- translated by f2c (version 20030320).
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
/* /+   gentree.f                                                               +/ */
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
/* /+ $Id: gentree.c,v 1.1 2005-01-03 00:17:43 paklein Exp $ +/ */
/* /+***************************************************************************+/ */
/*<    >*/
/* Subroutine */ int gentree_(integer *n, integer *aptrs, integer *ainds, 
	doublereal *avals, integer *beginleafnode, integer *sizes, integer *
	parent, integer *tptrs, integer *tinds, integer *dd, integer *myid, 
	integer *pp, integer *wrkint, integer *asize, integer *recvsizs, 
	integer *nrsiz, MPI_Comm *comm)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__, j, k, l, m, pn, endnodepar, col;
    extern /* Subroutine */ int genleaftree_(integer *, integer *, integer *, 
	    integer *, integer *, integer *);
    integer ierr;
    integer endnode;

/*<       implicit none >*/
/*<       include 'mpif.h' >*/
/*<       integer N,beginleafnode,myid,dd,pp,asize,nrsiz,comm >*/
/* -*- fortran -*- */

/* double precision functions */
/*<       integer aptrs(2,0:*),ainds(*),sizes(0:*) >*/
/*<       integer parent(0:*),tptrs(3,0:*),tinds(*),recvsizs(0:*) >*/
/*<       integer wrkint(0:*),i,j,k,l,m,pn,ierr,endnode,endnodepar,col >*/
/*<       double precision avals(*) >*/
/* assuming tinds size is >= N, it is used as a temporary array here */
/*<       endnode = beginleafnode+sizes(myid)-1 >*/
    /* Parameter adjustments */
    --tinds;
    --tptrs;
    --avals;
    --ainds;
    --aptrs;

    /* Function Body */
    endnode = *beginleafnode + sizes[*myid] - 1;
/*<       endnodepar = parent(endnode) >*/
    endnodepar = parent[endnode];
/*<    >*/
    genleaftree_(beginleafnode, &sizes[*myid], parent, &aptrs[1], &ainds[1], &
	    tinds[1]);
/*<       do i=beginleafnode,endnode-1 >*/
    i__1 = endnode - 1;
    for (i__ = *beginleafnode; i__ <= i__1; ++i__) {
/*<         if(parent(i).eq.-1) then >*/
	if (parent[i__] == -1) {
/*<           parent(i) = endnode >*/
	    parent[i__] = endnode;
/*<         end if >*/
	}
/*<       end do >*/
    }
/*<       parent(endnode) = endnodepar >*/
    parent[endnode] = endnodepar;
/*<    >*/

/*      call mpi_allgather(beginleafnode,1,MPI_INTEGER,
     +                   wrkint,1,MPI_INTEGER,comm,ierr) */
/*    mpi_allgather__(beginleafnode, &c__1, &c__11, wrkint, &c__1, &c__11, comm, &ierr); */
    MPI_Allgather(beginleafnode, 1, MPI_INT, wrkint, 1, MPI_INT, *comm);
/*<    >*/

/*      call mpi_allgatherv(parent(beginleafnode),sizes(myid),
     +                    MPI_INTEGER,parent,sizes,wrkint,
     +                    MPI_INTEGER,comm,ierr) */
/*    mpi_allgatherv__(&parent[*beginleafnode], &sizes[*myid], &c__11, parent, 
	    sizes, wrkint, &c__11, comm, &ierr); */
    MPI_Allgatherv(&parent[*beginleafnode], sizes[*myid], MPI_INT, parent, sizes, wrkint, MPI_INT, *comm);

/*<       do i=0,N-1 >*/
    i__1 = *n - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
/*<         tptrs(2,i) = 0 >*/
	tptrs[i__ * 3 + 2] = 0;
/*<       end do >*/
    }
/* accumulate kid counts (tptrs(2,*)) */
/*<       do i=0,N-2 >*/
    i__1 = *n - 2;
    for (i__ = 0; i__ <= i__1; ++i__) {
/*<         j = parent(i) >*/
	j = parent[i__];
/*<         tptrs(2,j) = tptrs(2,j)+1 >*/
	++tptrs[j * 3 + 2];
/*<       end do >*/
    }
/* form displacements array (tptrs(1,*)) */
/*<       j = 1 >*/
    j = 1;
/*<       do i=0,N-1 >*/
    i__1 = *n - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
/*<         tptrs(1,i) = j >*/
	tptrs[i__ * 3 + 1] = j;
/*<         j = j+tptrs(2,i) >*/
	j += tptrs[i__ * 3 + 2];
/*<         tptrs(2,i) = 0 >*/
	tptrs[i__ * 3 + 2] = 0;
/*<       end do >*/
    }
/* fill tinds */
/*<       do i=0,N-2 >*/
    i__1 = *n - 2;
    for (i__ = 0; i__ <= i__1; ++i__) {
/*<         j = parent(i) >*/
	j = parent[i__];
/*<         tinds(tptrs(1,j)+tptrs(2,j)) = i >*/
	tinds[tptrs[j * 3 + 1] + tptrs[j * 3 + 2]] = i__;
/*<         tptrs(2,j) = tptrs(2,j)+1 >*/
	++tptrs[j * 3 + 2];
/*<       end do >*/
    }
/*<       pn = 1 >*/
    pn = 1;
/*<       i = 0 >*/
    i__ = 0;
/*<       do while (i.lt.nrsiz) >*/
    while(i__ < *nrsiz) {
/*<         col = recvsizs(i) >*/
	col = recvsizs[i__];
/*<         j = aptrs(2,col) >*/
	j = aptrs[(col << 1) + 2];
/*<         if(j.gt.0) then >*/
	if (j > 0) {
/*<           k = aptrs(1,col) >*/
	    k = aptrs[(col << 1) + 1];
/*<           l = 0 >*/
	    l = 0;
/*<           do while (ainds(k+l).lt.col) >*/
	    while(ainds[k + l] < col) {
/*<             l = l+1 >*/
		++l;
/*<           end do >*/
	    }
/*<           m = 0 >*/
	    m = 0;
/*<           do l=l,j-1 >*/
	    i__1 = j - 1;
	    for (l = l; l <= i__1; ++l) {
/*<             ainds(pn+m) = ainds(k+l) >*/
		ainds[pn + m] = ainds[k + l];
/*<             avals(pn+m) = avals(k+l) >*/
		avals[pn + m] = avals[k + l];
/*<             m = m+1 >*/
		++m;
/*<           end do >*/
	    }
/*<           aptrs(1,col) = pn >*/
	    aptrs[(col << 1) + 1] = pn;
/*<           aptrs(2,col) = m >*/
	    aptrs[(col << 1) + 2] = m;
/*<           pn = pn+m >*/
	    pn += m;
/*<         end if >*/
	}
/*<         i = i+2 >*/
	i__ += 2;
/*<       end do >*/
    }
/*<       asize = pn-1 >*/
    *asize = pn - 1;
/*<       end >*/
    return 0;
} /* gentree_ */

