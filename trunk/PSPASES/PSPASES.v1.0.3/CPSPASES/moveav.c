/* $Id: moveav.c,v 1.2 2004-12-12 23:00:41 paklein Exp $ */
/* moveav.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

/* debugging */
/* #undef __DO_DEBUG__ */
#define __DO_DEBUG__ 1

#include "mpi.h"
#include "pspases_f2c.h"

/* Table of constant values */
static integer c__9 = 9;
static integer c__1 = 1;
static integer c__11 = 11;
static integer c__13 = 13;

/* /+***************************************************************************+/ */
/* /+                                                                           +/ */
/* /+   (C) Copyright IBM Corporation, 1997                                     +/ */
/* /+   (C) Copyright Regents of the University of Minnesota, 1997              +/ */
/* /+                                                                           +/ */
/* /+   moveav.f                                                                +/ */
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
/* /+ $Id: moveav.c,v 1.2 2004-12-12 23:00:41 paklein Exp $ +/ */
/* /+***************************************************************************+/ */
/*<    >*/

integer lbit_shift(integer a, integer b) {
	return b >= 0 ? a << b : (integer)((uinteger)a >> -b);
};

/* Subroutine */ int moveav_(integer *n, integer *dd, integer *pp, integer *
	lgblk, integer *myid, integer *rowdista, integer *mynnodes, integer *
	order, integer *aptrs, integer *ainds, doublereal *avals, doublereal *
	pavals, integer *wrkint, integer *maxnzpercol, integer *ranmasks, 
	MPI_Comm *comm)

/* dummy declarations for f2c 
, integer *gorder, integer *whichsnode, integer *tainds, 
	doublereal *sendvals, integer *i2, integer *nsend2)
*/
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer beginrow, i__, j, k, l, m;
    extern /* Subroutine */ int ikeysortf_(integer *, integer *, integer *);
    static integer is1, col, ppc, ppg, ppr, row, ierr, proc, prcv, pscv, psdv,
	     prdv;
    static integer ptr_sendvals__, ptr_c__, nsend, ptr_r__;
    static integer bmaskc, bmaskr, fptr_r__;
    static integer itainds, pgrsize;

/*<       implicit none >*/
/*<       include 'mpif.h' >*/
/*<       integer*8 loc >*/
/* -*- fortran -*- */

/* double precision functions */

/*<       integer rowdista(0:*),order(0:*),aptrs(2,0:*),ainds(*) >*/
/*<       integer wrkint(0:*),ranmasks(5,0:*) >*/
/*<       integer N,dd,pp,lgblk,myid,mynnodes,maxnzpercol,comm >*/
/*<       double precision avals(*),pavals(*) >*/
/*      integer, allocatable :: gorder(:),whichsnode(:),tainds(:) */
	integer *gorder, *whichsnode, *tainds;
/*<       integer gorder, whichsnode, tainds >*/
/*<       integer i2 >*/
/*<       dimension gorder(0:N-1) >*/
/*<       dimension whichsnode(0:mynnodes-1) >*/
/*<       dimension tainds(2*i2) >*/
/*     double precision, allocatable :: sendvals(:) */
	double *sendvals;
/*<       double precision sendvals >*/
/*<       integer nsend2 >*/
/*<       dimension sendvals(0:nsend2-1) >*/
/*<       integer proc,pgrsize,ierr,bmaskr,bmaskc,row,col >*/
/*<       integer i,j,k,l,m,ptr_r,fptr_r,ptr_c,itainds >*/
/*<       integer is1,nsend,ptr_sendvals >*/
/*<       integer beginrow >*/
/*<       integer pscv,psdv,prcv,prdv,ppr,ppc,ppg >*/
/*<       pscv = 0 >*/
    /* Parameter adjustments */
    --aptrs;
    --ainds;
    --avals;
    --pavals;
    --ranmasks;
 /*   --tainds; */

    /* Function Body */
    pscv = 0;
/*<       psdv = pp >*/
    psdv = *pp;
/*<       prcv = 2*pp >*/
    prcv = *pp << 1;
/*<       prdv = 3*pp >*/
    prdv = *pp * 3;
/*<       ppr  = 8*pp >*/
    ppr = *pp << 3;
/*<       ppc  = 9*pp >*/
    ppc = *pp * 9;
/*<       ppg  = 10*pp >*/
    ppg = *pp * 10;
/*<       pgrsize = ishft(1,ishft(dd,-1)) >*/
    pgrsize = lbit_shift(1, lbit_shift(*dd, -1));

/*     allocate(gorder(0:N-1),stat=is1) */
	gorder = (integer*) malloc((*n)*sizeof(integer));
/*<       if(is1.ne.0) then >*/
    if (!gorder) {
/*<         print *,'Error in allocate' >*/
		printf("%d: Error in allocate", *myid);
/*<         call mpi_abort(comm,1,ierr) >*/
		MPI_Abort(*comm, 1);
/*<       end if >*/
    }
/*<       beginrow = rowdista(myid) >*/
    beginrow = rowdista[*myid];
/*<       do proc=0,pp-1 >*/
    i__1 = *pp - 1;
    for (proc = 0; proc <= i__1; ++proc) {
/*<         wrkint(pscv+proc) = rowdista(proc+1)-rowdista(proc) >*/
	wrkint[pscv + proc] = rowdista[proc + 1] - rowdista[proc];
/*<       end do >*/
    }
/*<    >*/

/*    call mpi_allgatherv(order,mynnodes,MPI_INTEGER,
     +                    gorder,wrkint(pscv),rowdista,MPI_INTEGER,
     +                    comm,ierr) */
/*  mpi_allgatherv__(order, mynnodes, &c__11, gorder, &wrkint[pscv], rowdista,
	     &c__11, comm, &ierr); */
	MPI_Allgatherv(order, *mynnodes, MPI_INT, gorder, &wrkint[pscv], rowdista, MPI_INT, *comm);

/*<       do proc=0,pp-1 >*/
    i__1 = *pp - 1;
    for (proc = 0; proc <= i__1; ++proc) {
/*<         wrkint(pscv+proc) = 0 >*/
	wrkint[pscv + proc] = 0;
/*<       end do >*/
    }

/*     allocate(whichsnode(0:mynnodes-1),stat=is1) */
	whichsnode = (integer*) malloc((*mynnodes)*sizeof(integer));
/*<       if(is1.ne.0) then >*/
    if (!whichsnode) {
/*<         print *,'Error in allocate' >*/
		printf("%d: Error in allocate", *myid);

/*<         call mpi_abort(comm,1,ierr) >*/
		MPI_Abort(*comm, 1);
/*<       end if >*/
    }
/*<       i = aptrs(1,mynnodes-1)+aptrs(2,mynnodes-1)-1 >*/
    i__ = aptrs[(*mynnodes - 1 << 1) + 1] + aptrs[(*mynnodes - 1 << 1) + 2] - 
	    1;

/*     allocate(tainds(2*i),stat=is1) */
	tainds = (integer*) malloc(2*i__*sizeof(integer));
	--tainds;
/*<       if(is1.ne.0) then >*/
    if (!tainds) {
/*<         print *,'Error in allocate' >*/
		printfd("%d: Error in allocate", *myid);
/*<         call mpi_abort(comm,1,ierr) >*/
		MPI_Abort(*comm, 1);
/*<       end if >*/
    }
/*<       itainds = 1 >*/
    itainds = 1;
/*<       do i=0,mynnodes-1 >*/
    i__1 = *mynnodes - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
/*<         col = gorder(beginrow+i) >*/
	col = gorder[beginrow + i__];
/*<         j = 0 >*/
	j = 0;
/*<         do while (col.lt.ranmasks(1,j) .or. col.gt.ranmasks(2,j)) >*/
	while(col < ranmasks[j * 5 + 1] || col > ranmasks[j * 5 + 2]) {
/*<           j = j+1 >*/
	    ++j;
/*<         end do >*/
	}
/*<         proc   = ranmasks(3,j) >*/
	proc = ranmasks[j * 5 + 3];
/*<         bmaskr = ranmasks(4,j) >*/
	bmaskr = ranmasks[j * 5 + 4];
/*<         bmaskc = ranmasks(5,j) >*/
	bmaskc = ranmasks[j * 5 + 5];
/*<         fptr_r = wrkint(ppr+proc) >*/
	fptr_r__ = wrkint[ppr + proc];
/*<         ptr_c  = wrkint(ppc+proc)+iand(ishft(col,-lgblk),bmaskc) >*/
	ptr_c__ = wrkint[ppc + proc] + (lbit_shift(col, -(*lgblk)) & bmaskc);
/*<         k = aptrs(1,i) >*/
	k = aptrs[(i__ << 1) + 1];
/*<         l = aptrs(2,i) >*/
	l = aptrs[(i__ << 1) + 2];
/*<         do j=0,l-1 >*/
	i__2 = l - 1;
	for (j = 0; j <= i__2; ++j) {
/*<           tainds(itainds+j) = gorder(ainds(k+j)) >*/
	    tainds[itainds + j] = gorder[ainds[k + j]];
/*<         end do >*/
	}
/*<         call ikeysortf(l,tainds(itainds),tainds(itainds+l)) >*/
	ikeysortf_(&l, &tainds[itainds], &tainds[itainds + l]);
/*<         m = 0 >*/
	m = 0;
/*<         do while(tainds(itainds+m).lt.col) >*/
	while(tainds[itainds + m] < col) {
/*<           m = m+1 >*/
	    ++m;
/*<         end do >*/
	}
/*<         whichsnode(i) = m >*/
	whichsnode[i__] = m;
/*<         do j=m,l-1 >*/
	i__2 = l - 1;
	for (j = m; j <= i__2; ++j) {
/*<           row = tainds(itainds+j) >*/
	    row = tainds[itainds + j];
/*<           ptr_r = fptr_r+iand(ishft(row,-lgblk),bmaskr) >*/
	    ptr_r__ = fptr_r__ + (lbit_shift(row, -(*lgblk)) & bmaskr);
/*<           proc = wrkint(ppg+ptr_c*pgrsize+ptr_r) >*/
	    proc = wrkint[ppg + ptr_c__ * pgrsize + ptr_r__];
/*<           wrkint(pscv+proc) = wrkint(pscv+proc)+1 >*/
	    ++wrkint[pscv + proc];
/*<           tainds(itainds+j) = proc >*/
	    tainds[itainds + j] = proc;
/*<         end do >*/
	}
/*<         itainds = itainds+2*l >*/
	itainds += l << 1;
/*<       end do  >*/
    }

/*     deallocate(gorder) */
	free(gorder);
/*<       wrkint(psdv) = 0 >*/
    wrkint[psdv] = 0;
/*<       do proc=1,pp-1 >*/
    i__1 = *pp - 1;
    for (proc = 1; proc <= i__1; ++proc) {
/*<         wrkint(psdv+proc) = wrkint(psdv+proc-1)+wrkint(pscv+proc-1) >*/
	wrkint[psdv + proc] = wrkint[psdv + proc - 1] + wrkint[pscv + proc - 
		1];
/*<         wrkint(pscv+proc-1) = 0 >*/
	wrkint[pscv + proc - 1] = 0;
/*<       end do >*/
    }
/*<       nsend = wrkint(psdv+pp-1)+wrkint(pscv+pp-1) >*/
    nsend = wrkint[psdv + *pp - 1] + wrkint[pscv + *pp - 1];
/*<       wrkint(pscv+pp-1) = 0 >*/
    wrkint[pscv + *pp - 1] = 0;

/*     allocate(sendvals(0:nsend-1),stat=is1) */
	sendvals = (double*) malloc(nsend*sizeof(double));
/*<       if(is1.ne.0) then >*/
    if (!sendvals) {
/*<         print *,'Error in allocate' >*/
		printf("%d: Error in allocate", *myid);
/*<         call mpi_abort(comm,1,ierr) >*/
		MPI_Abort(*comm, 1);
/*<       end if >*/
    }
/*<       itainds = 1 >*/
    itainds = 1;
/*<       do i=0,mynnodes-1 >*/
    i__1 = *mynnodes - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
/*<         k = aptrs(1,i) >*/
	k = aptrs[(i__ << 1) + 1];
/*<         l = aptrs(2,i) >*/
	l = aptrs[(i__ << 1) + 2];
/*<         m = whichsnode(i) >*/
	m = whichsnode[i__];
/*<         do j=m,l-1 >*/
	i__2 = l - 1;
	for (j = m; j <= i__2; ++j) {
/*<           proc = tainds(itainds+j) >*/
	    proc = tainds[itainds + j];
/*<           ptr_sendvals = wrkint(psdv+proc)+wrkint(pscv+proc) >*/
	    ptr_sendvals__ = wrkint[psdv + proc] + wrkint[pscv + proc];
/*<           sendvals(ptr_sendvals) = avals(k+tainds(itainds+l+j)) >*/
	    sendvals[ptr_sendvals__] = avals[k + tainds[itainds + l + j]];
/*<           wrkint(pscv+proc) = wrkint(pscv+proc)+1 >*/
	    ++wrkint[pscv + proc];
/*<         end do >*/
	}
/*<         itainds = itainds+2*l >*/
	itainds += l << 1;
/*<       end do  >*/
    }

/*     deallocate(whichsnode) */
/*     deallocate(tainds) */
	free(whichsnode);
	++tainds;
	free(tainds);
/*<    >*/
/*    call mpi_alltoall(wrkint(pscv),1,MPI_INTEGER,wrkint(prcv),
     +                  1,MPI_INTEGER,comm,ierr) */
/*   mpi_alltoall__(&wrkint[pscv], &c__1, &c__11, &wrkint[prcv], &c__1, &c__11,
	     comm, &ierr); */
	MPI_Alltall(&wrkint[pscv], 1, MPI_INT, &wrkint[prcv], 1, MPI_INT, *comm);

/*<       wrkint(prdv) = 0 >*/
    wrkint[prdv] = 0;
/*<       do proc=1,pp-1 >*/
    i__1 = *pp - 1;
    for (proc = 1; proc <= i__1; ++proc) {
/*<         wrkint(prdv+proc) = wrkint(prdv+proc-1)+wrkint(prcv+proc-1) >*/
	wrkint[prdv + proc] = wrkint[prdv + proc - 1] + wrkint[prcv + proc - 
		1];
/*<       end do >*/
    }
/*<    >*/
/*    call mpi_alltoallv(sendvals,wrkint(pscv),
     +                   wrkint(psdv),MPI_DOUBLE_PRECISION,pavals,
     +                   wrkint(prcv),wrkint(prdv),
     +                   MPI_DOUBLE_PRECISION,comm,ierr) */
/*    mpi_alltoallv__(sendvals, &wrkint[pscv], &wrkint[psdv], &c__13, &pavals[1]
	    , &wrkint[prcv], &wrkint[prdv], &c__13, comm, &ierr); */
	MPI_Alltoallv(sendvals, &wrkint[pscv], &wrkint[psdv], MPI_DOUBLE,
                &pavals[1], &wrkint[prcv], &wrkint[prdv], MPI_DOUBLE, *comm);

/*     deallocate(sendvals) */
	free(sendvals);
/*<       end >*/
    return 0;
} /* moveav_ */

