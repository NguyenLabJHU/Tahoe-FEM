/* $Id: premovea.c,v 1.1 2004-12-30 00:07:20 paklein Exp $ */
/* premovea.f -- translated by f2c (version 20030320).
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
static integer c__9 = 9;
static integer c__1 = 1;
static integer c__11 = 11;
static integer c__3 = 3;
static integer c__0 = 0;

/* /+***************************************************************************+/ */
/* /+                                                                           +/ */
/* /+   (C) Copyright IBM Corporation, 1997                                     +/ */
/* /+   (C) Copyright Regents of the University of Minnesota, 1997              +/ */
/* /+                                                                           +/ */
/* /+   premovea.f                                                              +/ */
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
/* /+ $Id: premovea.c,v 1.1 2004-12-30 00:07:20 paklein Exp $ +/ */
/* /+***************************************************************************+/ */
static int max(int a, int b) {
	return (a > b) ? a : b;
};

static integer lbit_shift(integer a, integer b) {
	return b >= 0 ? a << b : (integer)((uinteger)a >> -b);
};

/*<    >*/
/* Subroutine */ int premovea_(integer *n, integer *dd, integer *pp, integer *
	lgblk, integer *myid, integer *mynnodes, integer *rowdist, integer *
	order, integer *sizes, integer *mybeginleaf, integer *pasize, integer 
	*aptrs, integer *ainds, integer *tainds, integer *parent, integer *
	gorder, integer *temparr2, integer *wrkint, integer *ranmasks, 
	integer *whichsnode, integer *maxnzpercol, integer *checksymm, 
	integer *sortinds, MPI_Comm *comm)
	
/* dummy arguments needed for f2c
	, integer *sendinds, integer *
	sendsizs, integer *recvinds, integer *recvsizs, integer *
	iwillsend_inds2__, integer *iwillreceive_inds2__)
*/
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    integer sepbegin, sepindex, i__, j, k, l, m, nodebegin;
    extern /* Subroutine */ int ikeysortf_(integer *, integer *, integer *);
    integer iwillreceive_inds__, kc, ml, kr, is1, col, ppc, ppg, psi, ppr, 
	    pps, row;

    integer lend, kend, prci, psci, psdi, ierr, proc, prdi, psnb, tppc, prcs, 
	    pscs, psds, prds, tppr, psts;

    integer ptr_sendinds__, cbits, level, gbits, ptr_c__, gsize, rbits, 
	    ptr_r__, nseps;

    integer bmaskc, bmaskr;
    logical ifound;
    integer snodes, fptr_r__, nstree;

    integer iwillsend_inds__, nodeend, pgcsize, sepoffs, pgrsize;

/*<       implicit none >*/
/*<       include 'mpif.h' >*/
/*<       integer N,dd,pp,lgblk,myid,mynnodes,pasize,comm >*/
/* -*- fortran -*- */

/* double precision functions */
/*<       integer rowdist(0:*),order(0:*),sizes(0:*),ranmasks(5,0:*) >*/
/*<       integer aptrs(2,0:*),ainds(*),tainds(*),wrkint(0:*) >*/
/*<       integer gorder(0:N-1),temparr2(0:*),parent(0:*) >*/
/*<       integer whichsnode(0:*),maxnzpercol,checksymm,sortinds >*/
/*     integer, allocatable :: sendinds(:),sendsizs(:) */
	integer *sendinds, *sendsizs;

/*<       integer sendinds(*),sendsizs(*) >*/
/*<       dimension sendinds(0:iwillsend_inds2-1) >*/
/*<       dimension sendsizs(0:iwillsend_inds2-1)       >*/
/*     integer, allocatable :: recvinds(:),recvsizs(:) */
	integer *recvinds, *recvsizs;

/*<       integer recvinds(*),recvsizs(*) >*/
/*<       dimension recvinds(0:iwillreceive_inds2-1) >*/
/*<       dimension recvsizs(0:iwillreceive_inds2-1)       >*/
/*<       integer proc,rbits,cbits,pgrsize,pgcsize >*/
/*<       integer tppr,tppc >*/
/*<       integer ierr,level,nseps,sepindex,nodebegin,sepbegin,snodes >*/
/*<       integer sepoffs,gbits,gsize,bmaskr,bmaskc,row,col >*/
/*<       integer mybeginleaf,nodeend >*/
/*<       integer kr,kc,i,j,k,l,m,ptr_r,fptr_r,ptr_c >*/
/*<       integer iwillsend_inds,is1 >*/
/*<       integer iwillreceive_inds,ptr_sendinds >*/
/*<       integer lend,m1,ml,mk,kend,nstree >*/
/*<       logical ifound >*/
/*<       integer psci,psdi,prci,prdi,pscs,psds,prcs,prds >*/
/*<       integer ppr,ppc,ppg,psi,pps,psts,psnb >*/
/* non-overlapping wrkint storage */
/*<       psci = 0 >*/
    /* Parameter adjustments */
    --aptrs;
    --ainds;
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
/* overlapping wrkint storage */
/*<       psi  = psci >*/
    psi = psci;
/*<       pps  = psdi >*/
    pps = psdi;
/*<       psts = psdi >*/
    psts = psdi;
/*<       psnb = prdi >*/
    psnb = prdi;
/*<       maxnzpercol = 0 >*/
    *maxnzpercol = 0;
/*<       do i=0,mynnodes-1 >*/
    i__1 = *mynnodes - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
/*<         k = aptrs(1,i) >*/
	k = aptrs[(i__ << 1) + 1];
/*<         l = aptrs(2,i) >*/
	l = aptrs[(i__ << 1) + 2];
/*<         maxnzpercol = max(maxnzpercol,l) >*/
	*maxnzpercol = max(*maxnzpercol,l);
/*<         do j=0,l-1 >*/
	i__2 = l - 1;
	for (j = 0; j <= i__2; ++j) {
/*<           tainds(k+j) = ainds(k+j) >*/
	    tainds[k + j] = ainds[k + j];
/*<         end do >*/
	}
/*<         if(sortinds.eq.1 .and. checksymm.eq.1) then >*/
	if (*sortinds == 1 && *checksymm == 1) {
/*<           call ikeysortf(l,tainds(k),temparr2) >*/
	    ikeysortf_(&l, &tainds[k], temparr2);
/*<         end if >*/
	}
/*<       end do >*/
    }
/*<       if(checksymm.eq.1) then >*/
    if (*checksymm == 1) {
/*      check for symmetry of the input matrix. */
/*<       tppr =  0 >*/
	tppr = 0;
/*<       tppc =  mynnodes >*/
	tppc = *mynnodes;
/*<       do proc=0,pp-1 >*/
	i__1 = *pp - 1;
	for (proc = 0; proc <= i__1; ++proc) {
/*<         wrkint(psci+proc) = 0 >*/
	    wrkint[psci + proc] = 0;
/*<         do i=rowdist(proc),rowdist(proc+1)-1 >*/
	    i__2 = rowdist[proc + 1] - 1;
	    for (i__ = rowdist[proc]; i__ <= i__2; ++i__) {
/*<           temparr2(tppc+i) = proc >*/
		temparr2[tppc + i__] = proc;
/*<         end do >*/
	    }
/*<       end do >*/
	}
/*<       do i=0,mynnodes-1 >*/
	i__1 = *mynnodes - 1;
	for (i__ = 0; i__ <= i__1; ++i__) {
/*<         k = aptrs(1,i) >*/
	    k = aptrs[(i__ << 1) + 1];
/*<         l = aptrs(2,i) >*/
	    l = aptrs[(i__ << 1) + 2];
/*<         m = 0 >*/
	    m = 0;
/*<         do while (tainds(k+m).le.i+rowdist(myid)) >*/
	    while(tainds[k + m] <= i__ + rowdist[*myid]) {
/*<           m = m+1 >*/
		++m;
/*<         end do >*/
	    }
/*<         temparr2(tppr+i) = m >*/
	    temparr2[tppr + i__] = m;
/*<         do m=k+m,k+l-1 >*/
	    i__2 = k + l - 1;
	    for (m = k + m; m <= i__2; ++m) {
/*<           proc = temparr2(tppc+tainds(m)) >*/
		proc = temparr2[tppc + tainds[m]];
/*<           wrkint(psci+proc) = wrkint(psci+proc)+1 >*/
		++wrkint[psci + proc];
/*<         end do >*/
	    }
/*<       end do >*/
	}
/*<       wrkint(psdi) = 0 >*/
	wrkint[psdi] = 0;
/*<       iwillsend_inds = wrkint(psci) >*/
	iwillsend_inds__ = wrkint[psci];
/*<       do proc=1,pp-1 >*/
	i__1 = *pp - 1;
	for (proc = 1; proc <= i__1; ++proc) {
/*<         iwillsend_inds = iwillsend_inds+wrkint(psci+proc) >*/
	    iwillsend_inds__ += wrkint[psci + proc];
/*<         wrkint(psdi+proc) = wrkint(psdi+proc-1)+wrkint(psci+proc-1) >*/
	    wrkint[psdi + proc] = wrkint[psdi + proc - 1] + wrkint[psci + 
		    proc - 1];
/*<         wrkint(psci+proc-1) = 0 >*/
	    wrkint[psci + proc - 1] = 0;
/*<       end do >*/
	}
/*<       wrkint(psci+pp-1) = 0 >*/
	wrkint[psci + *pp - 1] = 0;
/*     allocate(sendinds(0:iwillsend_inds-1),stat=is1) */
	sendinds = (integer*) malloc(iwillsend_inds__*sizeof(integer));
	if (!sendinds) {
		printf("%d: Error in allocate 1", *myid);
	    MPI_Abort(*comm, 1);
	}

/*     allocate(sendsizs(0:iwillsend_inds-1),stat=is1) */
	sendsizs = (integer*) malloc(iwillsend_inds__*sizeof(integer));
	if (!sendsizs) {
		printf("%d: Error in allocate 2", *myid);
		MPI_Abort(*comm, 1);
	}

/*<       do i=0,mynnodes-1 >*/
	i__1 = *mynnodes - 1;
	for (i__ = 0; i__ <= i__1; ++i__) {
/*<         k = aptrs(1,i) >*/
	    k = aptrs[(i__ << 1) + 1];
/*<         l = aptrs(2,i) >*/
	    l = aptrs[(i__ << 1) + 2];
/*<         m = temparr2(tppr+i) >*/
	    m = temparr2[tppr + i__];
/*<         do m=k+m,k+l-1 >*/
	    i__2 = k + l - 1;
	    for (m = k + m; m <= i__2; ++m) {
/*<           proc = temparr2(tppc+tainds(m)) >*/
		proc = temparr2[tppc + tainds[m]];
/*<           ptr_sendinds = wrkint(psdi+proc)+wrkint(psci+proc) >*/
		ptr_sendinds__ = wrkint[psdi + proc] + wrkint[psci + proc];
/*<           sendinds(ptr_sendinds) = tainds(m) >*/
		sendinds[ptr_sendinds__] = tainds[m];
/*<           sendsizs(ptr_sendinds) = i+rowdist(myid) >*/
		sendsizs[ptr_sendinds__] = i__ + rowdist[*myid];
/*<           wrkint(psci+proc) = wrkint(psci+proc)+1 >*/
		++wrkint[psci + proc];
/*<         end do >*/
	    }
/*<       end do >*/
	}
/*<    >*/
/*    call mpi_alltoall(wrkint(psci),1,MPI_INTEGER,wrkint(prci),
     +                  1,MPI_INTEGER,comm,ierr) */
/* mpi_alltoall__(&wrkint[psci], &c__1, &c__11, &wrkint[prci], &c__1, &c__11, comm, &ierr); */
	MPI_Alltoall(&wrkint[psci],1,MPI_INT,&wrkint[prci],1,MPI_INT,*comm);

/*<       wrkint(prdi) = 0 >*/
	wrkint[prdi] = 0;
/*<       iwillreceive_inds = wrkint(prci) >*/
	iwillreceive_inds__ = wrkint[prci];
/*<       do proc=1,pp-1 >*/
	i__1 = *pp - 1;
	for (proc = 1; proc <= i__1; ++proc) {
/*<         iwillreceive_inds = iwillreceive_inds+wrkint(prci+proc) >*/
	    iwillreceive_inds__ += wrkint[prci + proc];
/*<         wrkint(prdi+proc) = wrkint(prdi+proc-1)+wrkint(prci+proc-1) >*/
	    wrkint[prdi + proc] = wrkint[prdi + proc - 1] + wrkint[prci + 
		    proc - 1];
/*<       end do >*/
	}
/*     allocate(recvinds(0:iwillreceive_inds-1),stat=is1) */
	recvinds = (integer*) malloc(iwillreceive_inds__*sizeof(integer));
	if (!recvinds) {
		printf("%d: Error in allocate 4", *myid);
		MPI_Abort(*comm, 1);
	}

/*     allocate(recvsizs(0:iwillreceive_inds-1),stat=is1) */
	recvsizs = (integer*) malloc(iwillreceive_inds__*sizeof(integer));
	if (!recvsizs) {
		printf("%d: Error in allocate 5", *myid);
		MPI_Abort(*comm, 1);
	}
/*<    >*/
/*       call mpi_alltoallv(sendinds,wrkint(psci),
     +                   wrkint(psdi),MPI_INTEGER,recvinds,
     +                   wrkint(prci),wrkint(prdi),
     +                   MPI_INTEGER,comm,ierr) */
/*  mpi_alltoallv__(sendinds, &wrkint[psci], &wrkint[psdi], &c__11, 
		recvinds, &wrkint[prci], &wrkint[prdi], &c__11, comm, &ierr); */
	MPI_Alltoallv(sendinds, &wrkint[psci], &wrkint[psdi], MPI_INT, 
		          recvinds, &wrkint[prci], &wrkint[prdi], MPI_INT, *comm);

/*<    >*/
/*    call mpi_alltoallv(sendsizs,wrkint(psci),
     +                   wrkint(psdi),MPI_INTEGER,recvsizs,
     +                   wrkint(prci),wrkint(prdi),
     +                   MPI_INTEGER,comm,ierr) */
/* 	mpi_alltoallv__(sendsizs, &wrkint[psci], &wrkint[psdi], &c__11, 
		recvsizs, &wrkint[prci], &wrkint[prdi], &c__11, comm, &ierr); */
	MPI_Alltoallv(sendsizs, &wrkint[psci], &wrkint[psdi], MPI_INT,
	              recvsizs, &wrkint[prci], &wrkint[prdi], MPI_INT, *comm);

/*<       do i=0,mynnodes-1 >*/
	i__1 = *mynnodes - 1;
	for (i__ = 0; i__ <= i__1; ++i__) {
/*<         temparr2(i) = 0 >*/
	    temparr2[i__] = 0;
/*<       end do >*/
	}
/*<       do i=0,iwillreceive_inds-1 >*/
	i__1 = iwillreceive_inds__ - 1;
	for (i__ = 0; i__ <= i__1; ++i__) {
/*<         j = recvinds(i)-rowdist(myid) >*/
	    j = recvinds[i__] - rowdist[*myid];
/*<         k = recvsizs(i) >*/
	    k = recvsizs[i__];
/*<         l = aptrs(1,j) >*/
	    l = aptrs[(j << 1) + 1];
/*<         m = temparr2(j) >*/
	    m = temparr2[j];
/*<         if(tainds(l+m).ne.k) then >*/
	    if (tainds[l + m] != k) {
	    	printf("%d: Matrix found to be unsymmetric in indices\n", *myid);
			printf("%d: Col index %d is not found in row %d", *myid, k, recvinds[i__]);
			MPI_Abort(*comm, 0);
	    }
/*<         temparr2(j) = temparr2(j)+1 >*/
	    ++temparr2[j];
/*<       end do >*/
	}
/*     deallocate(sendinds) */
/*     deallocate(sendsizs) */
/*     deallocate(recvinds) */
/*     deallocate(recvsizs) */
	free(sendinds);
	free(sendsizs);
	free(recvinds);
	free(recvsizs);

/*<       end if  >*/
    }
/*      collect global order with an allgather operation. */
/*<       do proc=0,pp-1 >*/
    i__1 = *pp - 1;
    for (proc = 0; proc <= i__1; ++proc) {
/*<         wrkint(proc) = rowdist(proc+1)-rowdist(proc) >*/
	wrkint[proc] = rowdist[proc + 1] - rowdist[proc];
/*<       end do >*/
    }
/*<    >*/
    mpi_allgatherv__(order, mynnodes, &c__11, gorder, wrkint, rowdist, &c__11,
	     comm, &ierr);
/*<       do i=0,N-1 >*/
    i__1 = *n - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
/*<         temparr2(gorder(i)) = i >*/
	temparr2[gorder[i__]] = i__;
/*<       end do >*/
    }
/*      adjust sizes array to ensure non-zero size supernodes. */
/*<       do i=0,N-1 >*/
    i__1 = *n - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
/*<         temparr2(N+i) = i >*/
	temparr2[*n + i__] = i__;
/*<       end do >*/
    }
/*<       do level = 0,dd >*/
    i__1 = *dd;
    for (level = 0; level <= i__1; ++level) {
/*<         wrkint(psi+level) = 0 >*/
	wrkint[psi + level] = 0;
/*<       end do >*/
    }
/*<       nseps = ishft(pp,1)-1 >*/
    nseps = (*pp << 1) - 1;
/*<       sepindex = 0 >*/
    sepindex = 0;
/*<       level = dd >*/
    level = *dd;
/*<       nodebegin = 0 >*/
    nodebegin = 0;
/*<       nstree = 0 >*/
    nstree = 0;
/*<       do snodes=0,nseps-1 >*/
    i__1 = nseps - 1;
    for (snodes = 0; snodes <= i__1; ++snodes) {
/*<         sepbegin = nseps-(ishft(1,level+1)-1) >*/
	sepbegin = nseps - (lbit_shift(1, level + 1) - 1);
/*<         sepoffs  = wrkint(psi+level) >*/
	sepoffs = wrkint[psi + level];
/*<         sepindex = sepbegin+sepoffs >*/
	sepindex = sepbegin + sepoffs;
/*<         nodeend = nodebegin+sizes(sepindex) >*/
	nodeend = nodebegin + sizes[sepindex];
/*<         wrkint(psts+sepindex) = nstree+sizes(sepindex) >*/
	wrkint[psts + sepindex] = nstree + sizes[sepindex];
/*<         nodebegin = nodeend >*/
	nodebegin = nodeend;
/*<         wrkint(psi+level) = sepoffs+1 >*/
	wrkint[psi + level] = sepoffs + 1;
/*<         if(mod(sepoffs+1,2).eq.0) then >*/
	if ((sepoffs + 1) % 2 == 0) {
/*<           level = level-1 >*/
	    --level;
/*<           nstree = wrkint(psts+sepindex)+wrkint(psts+sepindex-1) >*/
	    nstree = wrkint[psts + sepindex] + wrkint[psts + sepindex - 1];
/*<         else >*/
	} else {
/*<           level = dd >*/
	    level = *dd;
/*<           nstree = 0 >*/
	    nstree = 0;
/*<         end if >*/
	}
/*<       end do >*/
    }
/*<       level = 0 >*/
    level = 0;
/*<       j = 2*pp-2 >*/
    j = (*pp << 1) - 2;
/*<       wrkint(psnb+j) = N-1-sizes(j) >*/
    wrkint[psnb + j] = *n - 1 - sizes[j];
/*<       do i=2*pp-2,pp,-1 >*/
    i__1 = *pp;
    for (i__ = (*pp << 1) - 2; i__ >= i__1; --i__) {
/*<         if(sizes(i).eq.0) then >*/
	if (sizes[i__] == 0) {
/*<           kend = wrkint(psnb+i) >*/
	    kend = wrkint[psnb + i__];
/*<           lend = kend >*/
	    lend = kend;
/*<           sizes(i) = 1 >*/
	    sizes[i__] = 1;
/*<           wrkint(psnb+i) = wrkint(psnb+i)-1 >*/
	    --wrkint[psnb + i__];
/*<           ifound = .false. >*/
	    ifound = FALSE_;
/*<           m = i >*/
	    m = i__;
/*<           j = 2*pp-1-i  >*/
	    j = (*pp << 1) - 1 - i__;
/*<           do while (.not.ifound .and. m.ge.pp)  >*/
	    while(! ifound && m >= *pp) {
/*<             k = 2*(pp-j)-1  >*/
		k = (*pp - j << 1) - 1;
/*<             l = 2*(pp-j)-2  >*/
		l = (*pp - j << 1) - 2;
/*<             if(wrkint(psts+l).gt.wrkint(psts+k)) then  >*/
		if (wrkint[psts + l] > wrkint[psts + k]) {
/*<               m = l >*/
		    m = l;
/*<               lend = lend - wrkint(psts+k) >*/
		    lend -= wrkint[psts + k];
/*<             else >*/
		} else {
/*<               m = k >*/
		    m = k;
/*<             end if >*/
		}
/*<             wrkint(psts+m) = wrkint(psts+m)-1 >*/
		--wrkint[psts + m];
/*<             if(sizes(m).gt.0) then >*/
		if (sizes[m] > 0) {
/*<               ifound = .true. >*/
		    ifound = TRUE_;
/*<             else >*/
		} else {
/*<               j = 2*pp-1-m >*/
		    j = (*pp << 1) - 1 - m;
/*<             end if >*/
		}
/*<           end do >*/
	    }
/*<           if(m.lt.pp .and. .not.ifound) then >*/
	    if (m < *pp && ! ifound) {
/*<             if(myid.eq.0) then >*/
		if (*myid == 0) {
			printf("%d: Cannot solve this problem on %d processors\n", *myid, *pp);
			printf("%d: Insufficient nodes available to pull into zero-size separators\n", *myid);
		}
/*<             call mpi_barrier(comm,ierr) >*/
/*		mpi_barrier__(comm, &ierr); */
/*<             call mpi_abort(comm,0,ierr) >*/
/*		mpi_abort__(comm, &c__0, &ierr); */
/*<           end if >*/
		MPI_Barrier(*comm);
		MPI_Abort(*comm, 0);
	    }
/*<           ml = temparr2(N+lend) >*/
	    ml = temparr2[*n + lend];
/*<           do col = lend+1,kend >*/
	    i__2 = kend;
	    for (col = lend + 1; col <= i__2; ++col) {
/*<             temparr2(N+col-1) = temparr2(N+col) >*/
		temparr2[*n + col - 1] = temparr2[*n + col];
/*<           end do >*/
	    }
/*<           temparr2(N+kend) = ml >*/
	    temparr2[*n + kend] = ml;
/*<           sizes(m) = sizes(m)-1 >*/
	    --sizes[m];
/*<         end if >*/
	}
/*<         j = 2*pp-1-i  >*/
	j = (*pp << 1) - 1 - i__;
/*<         k = 2*(pp-j)-1  >*/
	k = (*pp - j << 1) - 1;
/*<         l = 2*(pp-j)-2  >*/
	l = (*pp - j << 1) - 2;
/*<         wrkint(psnb+k) = wrkint(psnb+i)-sizes(k) >*/
	wrkint[psnb + k] = wrkint[psnb + i__] - sizes[k];
/*<         wrkint(psnb+l) = wrkint(psnb+i)-wrkint(psts+k)-sizes(l) >*/
	wrkint[psnb + l] = wrkint[psnb + i__] - wrkint[psts + k] - sizes[l];
/*<       end do >*/
    }
/*<       do i=0,pp-1 >*/
    i__1 = *pp - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
/*< 	if(sizes(i).le.0) then >*/
	if (sizes[i__] <= 0) {
	
	if (*myid == 0) {
		printf("%d: Cannot solve this problem on %d processors\n", *myid);
		printf("%d: Reason: Zero size partition encountered. Probably the\n", *myid);
		printf("%d: matrix is too small, or it has a block that is too dense\n", *myid);
		prtinf("%d: to form required number of partitions\n", *myid);
	    }
/*<           call mpi_barrier(comm,ierr) >*/
/*	    mpi_barrier__(comm, &ierr); */
/*<           call mpi_abort(comm,0,ierr) >*/
/*	    mpi_abort__(comm, &c__0, &ierr); */
/*< 	end if >*/
	MPI_Barrier(*comm);
	MPI_Abort(*comm, 0);
	}
/*<       end do >*/
    }
/*<       do i=0,N-1 >*/
    i__1 = *n - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
/*<         j = temparr2(N+i) >*/
	j = temparr2[*n + i__];
/*<         k = temparr2(j) >*/
	k = temparr2[j];
/*<         gorder(k) = i >*/
	gorder[k] = i__;
/*<       end do >*/
    }
/* adjust order and inds array to reflect new numbering. */
/*<       do i=0,mynnodes-1 >*/
    i__1 = *mynnodes - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
/*<         order(i) = gorder(i+rowdist(myid)) >*/
	order[i__] = gorder[i__ + rowdist[*myid]];
/*<       end do >*/
    }
/*<       do i=0,mynnodes-1 >*/
    i__1 = *mynnodes - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
/*<         k = aptrs(1,i) >*/
	k = aptrs[(i__ << 1) + 1];
/*<         l = aptrs(2,i) >*/
	l = aptrs[(i__ << 1) + 2];
/*<         do j=k,k+l-1 >*/
	i__2 = k + l - 1;
	for (j = k; j <= i__2; ++j) {
/*<           tainds(j) = gorder(tainds(j)) >*/
	    tainds[j] = gorder[tainds[j]];
/*<         end do >*/
	}
/*<         call ikeysortf(l,tainds(k),temparr2) >*/
	ikeysortf_(&l, &tainds[k], temparr2);
/*<       end do >*/
    }
/* preprocessing for filling up some bookkeeping arrays. */
/*<       rbits = ishft(dd,-1) >*/
    rbits = lbit_shift(*dd, -1);
/*<       cbits = dd-rbits >*/
    cbits = *dd - rbits;
/*<       pgrsize = ishft(1,rbits) >*/
    pgrsize = lbit_shift(1, rbits);
/*<       pgcsize = ishft(1,cbits) >*/
    pgcsize = lbit_shift(1, cbits);
/*<       do proc = 0,pp-1 >*/
    i__1 = *pp - 1;
    for (proc = 0; proc <= i__1; ++proc) {
/*<         kr = 0 >*/
	kr = 0;
/*<         i = ishft(proc,-1) >*/
	i__ = lbit_shift(proc, -1);
/*<         do j = 0,rbits-1 >*/
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
/*<         i = proc >*/
	i__ = proc;
/*<         do j = 0,cbits-1 >*/
	i__2 = cbits - 1;
	for (j = 0; j <= i__2; ++j) {
/*<           kc = ior(kc,ishft(iand(i,1),j)) >*/
	    kc |= lbit_shift(i__ & 1, j);
/*<           i = ishft(i,-2) >*/
	    i__ = lbit_shift(i__, -2);
/*<         end do >*/
	}
/*<         wrkint(ppr+proc) = kr >*/
	wrkint[ppr + proc] = kr;
/*<         wrkint(ppc+proc) = kc >*/
	wrkint[ppc + proc] = kc;
/*<         wrkint(ppg+kc*pgrsize+kr) = proc >*/
	wrkint[ppg + kc * pgrsize + kr] = proc;
/*<       end do >*/
    }
/*<       do level = 0,dd >*/
    i__1 = *dd;
    for (level = 0; level <= i__1; ++level) {
/*<         wrkint(psi+level) = 0 >*/
	wrkint[psi + level] = 0;
/*<       end do >*/
    }
/*<       nseps = ishft(pp,1)-1 >*/
    nseps = (*pp << 1) - 1;
/*<       sepindex = 0 >*/
    sepindex = 0;
/*<       level = dd >*/
    level = *dd;
/*<       nodebegin = 0 >*/
    nodebegin = 0;
/*<       do snodes=0,nseps-1 >*/
    i__1 = nseps - 1;
    for (snodes = 0; snodes <= i__1; ++snodes) {
/*<         sepbegin = nseps-(ishft(1,level+1)-1) >*/
	sepbegin = nseps - (lbit_shift((ftnlen)1, level + 1) - 1);
/*<         sepoffs  = wrkint(psi+level) >*/
	sepoffs = wrkint[psi + level];
/*<         sepindex = sepbegin+sepoffs >*/
	sepindex = sepbegin + sepoffs;
/*<         gbits = dd-level >*/
	gbits = *dd - level;
/*<         gsize = ishft(1,gbits) >*/
	gsize = lbit_shift((ftnlen)1, gbits);
/*<         rbits = ishft(gbits,-1) >*/
	rbits = lbit_shift(gbits, (ftnlen)-1);
/*<         proc  = sepoffs*gsize >*/
	proc = sepoffs * gsize;
/*<         cbits = gbits-rbits >*/
	cbits = gbits - rbits;
/*<         bmaskr= ishft(1,rbits)-1 >*/
	bmaskr = lbit_shift((ftnlen)1, rbits) - 1;
/*<         bmaskc= ishft(1,cbits)-1 >*/
	bmaskc = lbit_shift((ftnlen)1, cbits) - 1;
/*<         if(level.eq.dd .and. sepoffs.eq.myid) mybeginleaf = nodebegin >*/
	if (level == *dd && sepoffs == *myid) {
	    *mybeginleaf = nodebegin;
	}
/*<         nodeend = nodebegin+sizes(sepindex) >*/
	nodeend = nodebegin + sizes[sepindex];
/*<         ranmasks(1,snodes) = nodebegin >*/
	ranmasks[snodes * 5 + 1] = nodebegin;
/*<         ranmasks(2,snodes) = nodeend-1 >*/
	ranmasks[snodes * 5 + 2] = nodeend - 1;
/*<         ranmasks(3,snodes) = proc >*/
	ranmasks[snodes * 5 + 3] = proc;
/*<         ranmasks(4,snodes) = bmaskr >*/
	ranmasks[snodes * 5 + 4] = bmaskr;
/*<         ranmasks(5,snodes) = bmaskc >*/
	ranmasks[snodes * 5 + 5] = bmaskc;
/*<         do col = nodebegin,nodeend-1 >*/
	i__2 = nodeend - 1;
	for (col = nodebegin; col <= i__2; ++col) {
/*<           parent(col) = col+1 >*/
	    parent[col] = col + 1;
/*<         end do >*/
	}
/*<         nodebegin = nodeend >*/
	nodebegin = nodeend;
/*<         wrkint(psi+level) = sepoffs+1 >*/
	wrkint[psi + level] = sepoffs + 1;
/*<         if(mod(sepoffs+1,2).eq.0) then >*/
	if ((sepoffs + 1) % 2 == 0) {
/*<           parent(wrkint(pps+level)) = nodeend >*/
	    parent[wrkint[pps + level]] = nodeend;
/*<           level = level-1 >*/
	    --level;
/*<         else >*/
	} else {
/*<           wrkint(pps+level) = nodeend-1 >*/
	    wrkint[pps + level] = nodeend - 1;
/*<           level = dd >*/
	    level = *dd;
/*<         end if >*/
	}
/*<       end do >*/
    }
/*<       parent(N-1) = -1 >*/
    parent[*n - 1] = -1;
/* estimate the size of arrays to be communicated */
/*<       do proc=0,pp-1 >*/
    i__1 = *pp - 1;
    for (proc = 0; proc <= i__1; ++proc) {
/*<         wrkint(psci+proc) = 0 >*/
	wrkint[psci + proc] = 0;
/*<         wrkint(pscs+proc) = 0 >*/
	wrkint[pscs + proc] = 0;
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
/*<           wrkint(psci+proc) = wrkint(psci+proc)+1 >*/
	    ++wrkint[psci + proc];
/*<         end do >*/
	}
/*<       end do  >*/
    }
/*<    >*/
/*      call mpi_alltoall(wrkint(psci),1,MPI_INTEGER,wrkint(prci),
     +                  1,MPI_INTEGER,comm,ierr) */
/*    mpi_alltoall__(&wrkint[psci], &c__1, &c__11, &wrkint[prci], &c__1, &c__11,
	     comm, &ierr); */
	MPI_Alltoall(&wrkint[psci], 1, MPI_INT, 
	             &wrkint[prci], 1, MPI_INT, *comm);

/*<       iwillreceive_inds = wrkint(prci) >*/
    iwillreceive_inds__ = wrkint[prci];
/*<       wrkint(prdi) = 0 >*/
    wrkint[prdi] = 0;
/*<       do proc=1,pp-1 >*/
    i__1 = *pp - 1;
    for (proc = 1; proc <= i__1; ++proc) {
/*<         iwillreceive_inds = iwillreceive_inds+wrkint(prci+proc) >*/
	iwillreceive_inds__ += wrkint[prci + proc];
/*<         wrkint(prdi+proc) = wrkint(prdi+proc-1)+wrkint(prci+proc-1) >*/
	wrkint[prdi + proc] = wrkint[prdi + proc - 1] + wrkint[prci + proc - 
		1];
/*<       end do >*/
    }
/*<       pasize = iwillreceive_inds >*/
    *pasize = iwillreceive_inds__;
/*<       end >*/
    return 0;
} /* premovea_ */

