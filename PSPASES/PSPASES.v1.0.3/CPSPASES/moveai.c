/* $Id: moveai.c,v 1.2 2004-12-11 10:18:10 paklein Exp $ */
/* moveai.f -- translated by f2c (version 20030320).
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

/* /+***************************************************************************+/ */
/* /+                                                                           +/ */
/* /+   (C) Copyright IBM Corporation, 1997                                     +/ */
/* /+   (C) Copyright Regents of the University of Minnesota, 1997              +/ */
/* /+                                                                           +/ */
/* /+   moveai.f                                                                +/ */
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
/* /+ $Id: moveai.c,v 1.2 2004-12-11 10:18:10 paklein Exp $ +/ */
/* /+***************************************************************************+/ */
/*<    >*/

integer lbit_shift(integer a, integer b) {
	return b >= 0 ? a << b : (integer)((uinteger)a >> -b);
};

/* Subroutine */ int moveai_(integer *n, integer *dd, integer *pp, integer *
	lgblk, integer *myid, integer *mynnodes, integer *order, integer *
	paptrs, integer *recvsizs, integer *painds, integer *aptrs, integer *
	tainds, integer *wrkint, integer *ranmasks, integer *whichsnode, 
	integer *nsend, integer *nrecvsizs, MPI_Comm *comm)
	
/* dummy arguments needed for f2c
	, integer *sendinds, 
	integer *sendsizs, integer *iwillsend_sizs2__)
*/
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k;
    static integer is1, col, ppc, ppg, ppr, row, prci, psci, psdi, ierr, proc,
	     prdi, prcs, pscs, psds, prds;
    static integer ptr_sendinds__, ptr_sendsizs__, ptr_c__, ptr_r__;
    static integer bmaskc, bmaskr, fptr_r__, iwillsend_sizs__, pgrsize;

/*<       implicit none >*/
/*<       include 'mpif.h' >*/
/*<       integer*8 loc >*/
/* -*- fortran -*- */

/* double precision functions */

/*<       double precision MPI_WTIME, MPI_WTICK, PMPI_WTIME, PMPI_WTICK >*/
/*<       external MPI_WTIME, MPI_WTICK, PMPI_WTIME, PMPI_WTICK >*/
/*<       integer order(0:*),aptrs(2,0:*),paptrs(2,0:*),painds(*) >*/
/*<       integer tainds(*),wrkint(0:*),whichsnode(0:*),ranmasks(5,0:*) >*/
/*<       integer N,dd,pp,lgblk,myid,mynnodes,nsend,comm >*/
/*      integer, allocatable :: sendinds(:),sendsizs(:) */
	integer *sendinds, *sendsizs;
/*<       integer sendinds, sendsizs >*/
/*<       integer iwillsend_sizs2 >*/
/*<       dimension sendinds(0:nsend-1) >*/
/*<       dimension sendsizs(0:iwillsend_sizs2-1) >*/
/*<       integer recvsizs(0:*) >*/
/*<       integer proc,pgrsize,ierr,bmaskr,bmaskc,row,col >*/
/*<       integer i,j,k,l,m,ptr_r,fptr_r,ptr_c >*/
/*<       integer iwillsend_sizs,is1,ptr_sendsizs >*/
/*<       integer nrecvsizs,ptr_sendinds >*/
/*<       integer psci,psdi,prci,prdi,pscs,psds,prcs,prds,ppr,ppc,ppg >*/
/*<       psci = 0 >*/
    /* Parameter adjustments */
    --paptrs;
    --painds;
    --aptrs;
    --tainds;
    --ranmasks;

    /* Function Body */
    psci = 0;
/*<       psdi = pp >*/
    psdi = *pp;
/*<       prci = 2*pp >*/
    prci = *pp << 1;
/*<       prdi = 3*pp >*/
    prdi = *pp * 3;
/*<       pscs = 4*pp >*/
    pscs = *pp << 2;
/*<       psds = 5*pp >*/
    psds = *pp * 5;
/*<       prcs = 6*pp >*/
    prcs = *pp * 6;
/*<       prds = 7*pp >*/
    prds = *pp * 7;
/*<       ppr  = 8*pp >*/
    ppr = *pp << 3;
/*<       ppc  = 9*pp >*/
    ppc = *pp * 9;
/*<       ppg  = 10*pp >*/
    ppg = *pp * 10;
/*<       pgrsize = ishft(1,ishft(dd,-1)) >*/
    pgrsize = lbit_shift(1, lbit_shift(*dd, -1));
/*<       wrkint(psdi) = 0 >*/
    wrkint[psdi] = 0;
/*<       wrkint(psds) = 0 >*/
    wrkint[psds] = 0;
/*<       do proc=1,pp-1 >*/
    i__1 = *pp - 1;
    for (proc = 1; proc <= i__1; ++proc) {
/*<         wrkint(psdi+proc) = wrkint(psdi+proc-1)+wrkint(psci+proc-1) >*/
	wrkint[psdi + proc] = wrkint[psdi + proc - 1] + wrkint[psci + proc - 
		1];
/*<         wrkint(psds+proc) = wrkint(psds+proc-1)+wrkint(pscs+proc-1) >*/
	wrkint[psds + proc] = wrkint[psds + proc - 1] + wrkint[pscs + proc - 
		1];
/*<         wrkint(psci+proc-1) = 0 >*/
	wrkint[psci + proc - 1] = 0;
/*<         wrkint(pscs+proc-1) = 0 >*/
	wrkint[pscs + proc - 1] = 0;
/*<       end do >*/
    }
/*<       nsend = wrkint(psdi+pp-1)+wrkint(psci+pp-1) >*/
    *nsend = wrkint[psdi + *pp - 1] + wrkint[psci + *pp - 1];
/*<       iwillsend_sizs = wrkint(psds+pp-1)+wrkint(pscs+pp-1) >*/
    iwillsend_sizs__ = wrkint[psds + *pp - 1] + wrkint[pscs + *pp - 1];
/*<       wrkint(psci+pp-1) = 0 >*/
    wrkint[psci + *pp - 1] = 0;
/*<       wrkint(pscs+pp-1) = 0 >*/
    wrkint[pscs + *pp - 1] = 0;

/*      allocate(sendinds(0:nsend-1),stat=is1) */
	sendinds = (integer*) malloc((*nsend)*sizeof(integer));
/*<       if(is1.ne.0) then >*/
    if (!sendinds) {
/*<         print *,'Error in allocate' >*/
		printf("%d: Error in allocate", *myid);
/*<         call mpi_abort(comm,1,ierr) >*/
		MPI_Abort(*comm, 1);
/*<       end if >*/
    }

/*      allocate(sendsizs(0:iwillsend_sizs-1),stat=is1) */
	sendsizs = (integer*) malloc(iwillsend_sizs__*sizeof(integer));
/*<       if(is1.ne.0) then >*/
    if (!sendsizs) {
/*<         print *,'Error in allocate' >*/
		printf("%d: Error in allocate", *myid);
/*<         call mpi_abort(comm,1,ierr) >*/
		MPI_Abort(*comm, 1);
/*<       end if >*/
    }
/*<       do i=0,mynnodes-1 >*/
    i__1 = *mynnodes - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
/*<         col = order(i) >*/
	col = order[i__];
/*<         j = whichsnode(i) >*/
	j = whichsnode[i__];
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
/*<         do k=0,bmaskr >*/
	i__2 = bmaskr;
	for (k = 0; k <= i__2; ++k) {
/*<           proc = wrkint(ppg+ptr_c*pgrsize+fptr_r+k) >*/
	    proc = wrkint[ppg + ptr_c__ * pgrsize + fptr_r__ + k];
/*<           ptr_sendsizs = wrkint(psds+proc)+wrkint(pscs+proc) >*/
	    ptr_sendsizs__ = wrkint[psds + proc] + wrkint[pscs + proc];
/*<           sendsizs(ptr_sendsizs) = col >*/
	    sendsizs[ptr_sendsizs__] = col;
/*<           wrkint(prcs+proc) = ptr_sendsizs+1 >*/
	    wrkint[prcs + proc] = ptr_sendsizs__ + 1;
/*<           sendsizs(ptr_sendsizs+1) = 0 >*/
	    sendsizs[ptr_sendsizs__ + 1] = 0;
/*<           wrkint(pscs+proc) = wrkint(pscs+proc)+2 >*/
	    wrkint[pscs + proc] += 2;
/*<         end do >*/
	}
/*<         do j=aptrs(1,i),aptrs(1,i)+aptrs(2,i)-1 >*/
	i__2 = aptrs[(i__ << 1) + 1] + aptrs[(i__ << 1) + 2] - 1;
	for (j = aptrs[(i__ << 1) + 1]; j <= i__2; ++j) {
/*<           row = tainds(j) >*/
	    row = tainds[j];
/*<           ptr_r = fptr_r+iand(ishft(row,-lgblk),bmaskr) >*/
	    ptr_r__ = fptr_r__ + (lbit_shift(row, -(*lgblk)) & bmaskr);
/*<           proc = wrkint(ppg+ptr_c*pgrsize+ptr_r) >*/
	    proc = wrkint[ppg + ptr_c__ * pgrsize + ptr_r__];
/*<           ptr_sendinds = wrkint(psdi+proc)+wrkint(psci+proc) >*/
	    ptr_sendinds__ = wrkint[psdi + proc] + wrkint[psci + proc];
/*<           sendinds(ptr_sendinds) = row >*/
	    sendinds[ptr_sendinds__] = row;
/*<           sendsizs(wrkint(prcs+proc)) = sendsizs(wrkint(prcs+proc))+1 >*/
	    ++sendsizs[wrkint[prcs + proc]];
/*<           wrkint(psci+proc) = wrkint(psci+proc)+1 >*/
	    ++wrkint[psci + proc];
/*<         end do >*/
	}
/*<       end do  >*/
    }
/*<    >*/

/*< call mpi_alltoall(wrkint(pscs),1,MPI_INTEGER,wrkint(prcs),1,MPI_INTEGER,comm,ierr) >*/
/*  mpi_alltoall__(&wrkint[pscs], &c__1, &c__11, &wrkint[prcs], &c__1, &c__11, comm, &ierr);*/
	MPI_Alltoall(&wrkint[pscs], 1, MPI_INT, &wrkint[prcs], 1, MPI_INT, *comm);

/*<       wrkint(prds) = 0 >*/
    wrkint[prds] = 0;
/*<       do proc=1,pp-1 >*/
    i__1 = *pp - 1;
    for (proc = 1; proc <= i__1; ++proc) {
/*<         wrkint(prds+proc) = wrkint(prds+proc-1)+wrkint(prcs+proc-1) >*/
	wrkint[prds + proc] = wrkint[prds + proc - 1] + wrkint[prcs + proc - 
		1];
/*<       end do >*/
    }
/*<       nrecvsizs = wrkint(prds+pp-1)+wrkint(prcs+pp-1) >*/
    *nrecvsizs = wrkint[prds + *pp - 1] + wrkint[prcs + *pp - 1];
/*<    >*/
/*<      call mpi_alltoallv(sendinds,wrkint(psci),
     +                   wrkint(psdi),MPI_INTEGER,painds,
     +                   wrkint(prci),wrkint(prdi),
     +                   MPI_INTEGER,comm,ierr) >*/
/*    mpi_alltoallv__(sendinds, &wrkint[psci], &wrkint[psdi], &c__11, &painds[1]
	    , &wrkint[prci], &wrkint[prdi], &c__11, comm, &ierr); */
	MPI_Alltoallv(sendinds, &wrkint[psci], &wrkint[psdi], MPI_INT, 
	              &painds[1] , &wrkint[prci], &wrkint[prdi], MPI_INT, *comm);

/*<    >*/
/*<      call mpi_alltoallv(sendsizs,wrkint(pscs),
     +                   wrkint(psds),MPI_INTEGER,recvsizs,
     +                   wrkint(prcs),wrkint(prds),
     +                   MPI_INTEGER,comm,ierr) >*/
/*    mpi_alltoallv__(sendsizs, &wrkint[pscs], &wrkint[psds], &c__11, recvsizs, 
	    &wrkint[prcs], &wrkint[prds], &c__11, comm, &ierr); */
	MPI_Alltoallv(sendsizs, &wrkint[pscs], &wrkint[psds], MPI_INT, 
                  recvsizs, &wrkint[prcs], &wrkint[prds], MPI_INT, *comm);
	    
/*<       do i=0,N-1 >*/
    i__1 = *n - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
/*<         paptrs(1,i) = 1 >*/
	paptrs[(i__ << 1) + 1] = 1;
/*<         paptrs(2,i) = 0 >*/
	paptrs[(i__ << 1) + 2] = 0;
/*<       end do >*/
    }
/*<       i = 0 >*/
    i__ = 0;
/*<       j = 1 >*/
    j = 1;
/*<       do while (i.lt.nrecvsizs) >*/
    while(i__ < *nrecvsizs) {
/*<         col = recvsizs(i) >*/
	col = recvsizs[i__];
/*<         paptrs(2,col) = recvsizs(i+1) >*/
	paptrs[(col << 1) + 2] = recvsizs[i__ + 1];
/*<         i = i+2 >*/
	i__ += 2;
/*<         paptrs(1,col) = j >*/
	paptrs[(col << 1) + 1] = j;
/*<         j = j+paptrs(2,col) >*/
	j += paptrs[(col << 1) + 2];
/*<       end do >*/
    }
/*      deallocate(sendinds) */
/*      deallocate(sendsizs) */
	free(sendinds);
	free(sendsizs);

/*<       end >*/
    return 0;
} /* moveai_ */

