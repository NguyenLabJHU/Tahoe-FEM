/* $Id: trisolve.c,v 1.2 2004-12-13 09:12:51 paklein Exp $ */
/* trisolve.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

/* debugging */
#undef __DO_DEBUG__
/* #define __DO_DEBUG__ 1 */

#include "mpi.h"
#include "pspases_f2c.h"

/* Table of constant values */
static integer c__9 = 9;
static integer c__1 = 1;
static integer c__112 = 112;
static integer c__3 = 3;
static integer c__113 = 113;
static integer c__111 = 111;
static integer c__0 = 0;

/* /+***************************************************************************+/ */
/* /+                                                                           +/ */
/* /+   (C) Copyright IBM Corporation, 1997                                     +/ */
/* /+   (C) Copyright Regents of the University of Minnesota, 1997              +/ */
/* /+                                                                           +/ */
/* /+   trisolve.f                                                              +/ */
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
/* /+ $Id: trisolve.c,v 1.2 2004-12-13 09:12:51 paklein Exp $ +/ */
/* /+***************************************************************************+/ */
/*<    >*/

static integer lbit_shift(integer a, integer b) {
	return b >= 0 ? a << b : (integer)((uinteger)a >> -b);
};

static integer min(integer a, integer b) {
	return (a < b) ? a : b;
}

static integer max(integer a, integer b) {
	return (a > b) ? a : b;
}
/* Subroutine */ int trisolve_(integer *n, integer *rowdist, integer *order, 
	integer *lptrs, integer *linds, doublereal *lvals, integer *tptrs, 
	integer *tinds, integer *sup, integer *tsind, integer *lc, integer *
	iptrs, integer *ifopts, integer *nrhs, integer *options, doublereal *
	rhso, integer *ldo, doublereal *rhsc, integer *ldc, integer *ranmasks,
	MPI_Comm *comm)

/* dummy variables introduced for f2c
, integer *hvbtemp, integer *lrud, integer *wrkord0, 
	integer *wrkord1, integer *wrkord2, doublereal *ty, doublereal *
	dworkmj, doublereal *ordvals, integer *mynnodes2, integer *bnrhs2, 
	integer *ordvalsiz2, integer *wsolvesize2, integer *trhsize2, integer 
	*wrkord1siz2, integer *wrkord2siz2, integer *ns2, integer *dd2, 
	integer *maxvsize2, integer *hvbsize2)
*/
{
    /* System generated locals */
    integer rhso_dim1, rhso_offset, rhsc_dim1, rhsc_offset, ty_dim1, 
	    ty_offset, i__1, i__2;

    /* Local variables */
    extern /* Subroutine */ int pbsolve1_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, MPI_Comm *), pfsolve1_(integer *, integer 
	    *, integer *, integer *, integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, doublereal *, integer *, MPI_Comm *), preordbc_(
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, MPI_Comm *), preordbe_(integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, MPI_Comm *), 
	    getmyhvb_(integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *), preordxc_(
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, integer *, integer *, doublereal *, MPI_Comm *), pbsolvem_(
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, MPI_Comm 
	    *);
    integer maxhsize, nsupnode, mynnodes;
    extern /* Subroutine */ int pfsolvem_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     doublereal *, integer *, MPI_Comm *);
    integer maxvsize, i__;
    integer ordvalsiz, dd, pp, ns, sr, wrkord1siz, wrkord2siz, is1, 
	    supindsize, wsolvesize, rnr, psv, prv;
    extern /* Subroutine */ int getmysnodes_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    integer pdrc, pdsc, pdsd, pdrd, pirc, pisc, ierr, myid, pisd, pird, ppgr, 
	    nown, psvx, prvx, lgblk, myidh, bnrhs, pdrcx, pdscx, piscx, pisdx,
	     myidv, pircx, pirdx, pdsdx, pdrdx, piown, psloc;
    extern /* Subroutine */ int reordb_(integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, MPI_Comm *), reordx_(integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, doublereal *, doublereal *, integer *, MPI_Comm *);
    integer pslocx, rhsptr, uvlptr, hvbsize, uindptr, recvptr, uvecptr, 
	    trhsize, supsize;

/*<       implicit none >*/
/*<       include 'mpif.h' >*/
/*<       double precision zero >*/
/* -*- fortran -*- */

/* double precision functions */

/*<       parameter(zero=0.d0) >*/
/*<       integer rowdist(0:*),order(0:*),lptrs(3,0:*),linds(*) >*/
/*<       integer tptrs(3,0:*),tinds(*),sup(*),tsind(*),lc(*),iptrs(2,0:*) >*/
/*<       integer ifopts(0:*),options(0:*),nrhs,ldo,ldc,ranmasks(5,0:*) >*/
/*<       double precision lvals(*),rhso(0:ldo-1,*),rhsc(0:ldc-1,*) >*/
/*<       integer N,dd,lgblk,wsolvesize,ns,supindsize,myid,myidh,myidv,pp >*/
/*<       integer comm,nsupnode,supsize,maxvsize,maxhsize,i,j,k,m,hvbsize >*/
/*<       integer uindptr,uvecptr,recvptr,rhsptr,trhsize,is1,nown,ierr >*/
/*<       integer mynnodes,uvlptr,bnrhs,sr,ir,rnr >*/
/*<       integer wrkord1siz,wrkord2siz,ordvalsiz >*/
/*<       integer pisc,pisd,pirc,pird,pdsc,pdsd,pdrc,pdrd >*/
/*<       integer piscx,pisdx,pircx,pirdx,pdscx,pdsdx,pdrcx,pdrdx,ppgr >*/
/*<       integer piown,psloc,pslocx >*/
/*<       integer psv,prv,psvx,prvx >*/
/*     integer, allocatable :: hvbtemp(:),lrud(:) */
	integer *hvbtemp, *lrud;
/*<       integer hvbtemp, lrud >*/
/*<       integer ns2, dd2, maxvsize2, hvbsize2 >*/
/*<       dimension hvbtemp(hvbsize2 + maxvsize2) >*/
/*<       dimension lrud(ns2+ 4*dd2) >*/
/*     integer, allocatable :: wrkord0(:), wrkord1(:), wrkord2(:) */
	integer *wrkord0, *wrkord1, *wrkord2;
/*<       integer wrkord0, wrkord1, wrkord2 >*/
/*<       integer mynnodes2, wrkord1siz2, wrkord2siz2 >*/
/*<       dimension wrkord0(0:mynnodes2-1) >*/
/*<       dimension wrkord1(0:wrkord1siz2-1) >*/
/*<       dimension wrkord2(0:wrkord2siz2-1) >*/
/*     double precision, allocatable :: ty(:,:),dworkmj(:) */
	doublereal *ty, *dworkmj;
/*<       double precision ty, dworkmj >*/
/*<       integer bnrhs2, wsolvesize2, trhsize2 >*/
/*<       dimension ty(0:N-1,bnrhs2) >*/
/*<       dimension dworkmj(wsolvesize2 + 4*trhsize2) >*/
/*     double precision, allocatable :: ordvals(:) */
	doublereal* ordvals;
/*<       double precision ordvals >*/
/*<       integer ordvalsiz2 >*/
/*<       dimension ordvals(0:ordvalsiz2-1) >*/
/*<       N = ifopts(0) >*/
    /* Parameter adjustments */
    --lptrs;
    --linds;
    --lvals;
    --tptrs;
    --tinds;
    --sup;
    --tsind;
    --lc;
    --iptrs;
    rhso_dim1 = *ldo - 1 - 0 + 1;
    rhso_offset = 0 + rhso_dim1;
    rhso -= rhso_offset;
    rhsc_dim1 = *ldc - 1 - 0 + 1;
    rhsc_offset = 0 + rhsc_dim1;
    rhsc -= rhsc_offset;
    --ranmasks;
/* pointers not set yet
   ty_dim1 = *n - 1 - 0 + 1;
    ty_offset = 0 + ty_dim1;
    ty -= ty_offset;
    --dworkmj;
    --lrud;
    --hvbtemp;
*/

    /* Function Body */
    *n = ifopts[0];
/*<       dd = ifopts(1) >*/
    dd = ifopts[1];
/*<       lgblk = ifopts(2) >*/
    lgblk = ifopts[2];
/*<       wsolvesize = ifopts(3)  >*/
    wsolvesize = ifopts[3];
/*<       ns = ifopts(4) >*/
    ns = ifopts[4];
/*<       supindsize = ifopts(5) >*/
    supindsize = ifopts[5];
/*<       myid = ifopts(6) >*/
    myid = ifopts[6];
/*<       myidh = ifopts(7) >*/
    myidh = ifopts[7];
/*<       myidv = ifopts(8) >*/
    myidv = ifopts[8];
/*<       supsize = ifopts(9) >*/
    supsize = ifopts[9];
/*<       bnrhs = min(options(0),nrhs) >*/
    bnrhs = min(options[0],*nrhs);
/*<       pp = ishft(1,dd) >*/
    pp = lbit_shift(1, dd);

/*     allocate(lrud(ns+4*dd),stat=is1) */
	lrud = (integer*) malloc((ns+4*dd)*sizeof(integer));
/*<       if(is1.ne.0) then >*/
    if (!lrud) {
/*<         print *,'Allocate error' >*/
		printf("%d: Allocate error", myid);
/*<         call mpi_abort(comm,112,ierr) >*/
		MPI_Abort(*comm, 112);
/*<       end if >*/
    }
	--lrud; /* pointer adjustment */
/*<    >*/
    i__1 = *n - 1;
    getmysnodes_(&i__1, &sup[1], &tinds[1], &tptrs[1], n, &supsize, &lrud[1], 
	    &nsupnode, &lrud[ns + 1], &dd, &maxhsize, &maxvsize, &ns, &myid);
/*<       i = sup(tptrs(3,lrud(nsupnode)))  >*/
    i__ = sup[tptrs[lrud[nsupnode] * 3 + 3]];
/*<       maxvsize = max(maxvsize,lptrs(2,i)) >*/
/* Computing MAX */
    i__1 = maxvsize, i__2 = lptrs[i__ * 3 + 2];
    maxvsize = max(i__1,i__2);
/*<       hvbsize=nsupnode*(10+maxvsize+3*((maxhsize+1)/2+(maxvsize+1)/2)) >*/
    hvbsize = nsupnode * (maxvsize + 10 + ((maxhsize + 1) / 2 + (maxvsize + 1)
	     / 2) * 3);

/*     allocate(hvbtemp(hvbsize+maxvsize),stat=is1) */
	hvbtemp = (integer*) malloc((hvbsize+maxvsize)*sizeof(integer));
/*<       if(is1.ne.0) then >*/
    if (!hvbtemp) {
/*<         print *,myid,': Cannot allocate memory for hvbtemp',hvbsize >*/
		printf("%d: Cannot allocate memory for hvbtemp", myid);
/*<         call mpi_abort(comm,113,ierr) >*/
		MPI_Abort(*comm, 113);
/*<       end if >*/
    }
	--hvbtemp; /* pointer offset */
/*<       uindptr = 1+hvbsize >*/
    uindptr = hvbsize + 1;
/*<    >*/
    getmyhvb_(&lrud[1], &nsupnode, &sup[1], &supsize, &tsind[1], &supindsize, 
	    &tptrs[1], &tinds[1], n, &dd, &lgblk, &hvbtemp[1], &hvbsize, &
	    myidh, &myidv, &myid);
/*<       wsolvesize = 4 *wsolvesize * nrhs >*/
    wsolvesize = (wsolvesize << 2) * *nrhs;
/*<       trhsize = maxvsize * nrhs >*/
    trhsize = maxvsize * *nrhs;

/*     allocate(dworkmj(wsolvesize+4*trhsize),stat=is1) */
	dworkmj = (doublereal*) malloc((wsolvesize+4*trhsize)*sizeof(doublereal));
/*<       if(is1.ne.0) then >*/
    if (!dworkmj) {
/*<         print *,myid,': Memory Allocation Failure' >*/
		printf("%d: Memory Allocation Failure", myid);
/*<         call mpi_abort(comm,111,ierr) >*/
		MPI_Abort(*comm, 111);
/*<       end if >*/
    }
	--dworkmj; /* pointer offset */
/*<       uvecptr = 1+wsolvesize >*/
    uvecptr = wsolvesize + 1;
/*<       uvlptr  = uvecptr+trhsize >*/
    uvlptr = uvecptr + trhsize;
/*<       recvptr = uvlptr+trhsize >*/
    recvptr = uvlptr + trhsize;
/*<       rhsptr  = recvptr+trhsize >*/
    rhsptr = recvptr + trhsize;

/*     allocate(ty(0:N-1,bnrhs),stat=is1) */
	ty = (doublereal*) malloc((*n)*bnrhs*sizeof(doublereal));
/*<       if(is1.ne.0) then >*/
    if (!ty) {
/*<         print *,'memory allocation error' >*/
		printf("%d: memory allocation error", myid);
/*<         call mpi_abort(comm,0,ierr) >*/
		MPI_Abort(*comm, 0);
/*<       end if >*/
    }
	ty_dim1 = *n - 1 - 0 + 1;
    ty_offset = 0 + ty_dim1;
    ty -= ty_offset;
/*<       mynnodes = rowdist(myid+1)-rowdist(myid) >*/
    mynnodes = rowdist[myid + 1] - rowdist[myid];

/*     allocate(wrkord0(0:mynnodes-1),stat=is1) */
	wrkord0 = (integer*) malloc(mynnodes*sizeof(integer));
/*<       if(is1.ne.0) then >*/
    if (!wrkord0) {
/*<         print *,'memory allocation error' >*/
		printf("%d: memory allocation error", myid);
/*<         call mpi_abort(comm,0,ierr) >*/
		MPI_Abort(*comm, 0);
/*<       end if >*/
    }
/*<       wrkord1siz = 19*pp >*/
    wrkord1siz = pp * 19;

/*     allocate(wrkord1(0:wrkord1siz-1),stat=is1) */
	wrkord1 = (integer*) malloc(wrkord1siz*sizeof(integer));
/*<       if(is1.ne.0) then >*/
    if (!wrkord1) {
/*<         print *,'memory allocation error' >*/
		printf("%d: memory allocation error", myid);
/*<         call mpi_abort(comm,0,ierr) >*/
		MPI_Abort(*comm, 0);
/*<       end if >*/
    }
/*<       pisc  = 0 >*/
    pisc = 0;
/*<       pisd  = pp >*/
    pisd = pp;
/*<       pirc  = 2*pp >*/
    pirc = pp << 1;
/*<       pird  = 2*pp+pp >*/
    pird = (pp << 1) + pp;
/*<       pdsc  = 4*pp >*/
    pdsc = pp << 2;
/*<       pdsd  = 4*pp+pp >*/
    pdsd = (pp << 2) + pp;
/*<       pdrc  = 4*pp+2*pp >*/
    pdrc = (pp << 2) + (pp << 1);
/*<       pdrd  = 4*pp+2*pp+pp >*/
    pdrd = (pp << 2) + (pp << 1) + pp;
/*<       piscx = 8*pp >*/
    piscx = pp << 3;
/*<       pisdx = 8*pp+pp >*/
    pisdx = (pp << 3) + pp;
/*<       pircx = 8*pp+2*pp >*/
    pircx = (pp << 3) + (pp << 1);
/*<       pirdx = 8*pp+2*pp+pp >*/
    pirdx = (pp << 3) + (pp << 1) + pp;
/*<       pdscx = 8*pp+4*pp >*/
    pdscx = (pp << 3) + (pp << 2);
/*<       pdsdx = 8*pp+4*pp+pp >*/
    pdsdx = (pp << 3) + (pp << 2) + pp;
/*<       pdrcx = 8*pp+4*pp+2*pp >*/
    pdrcx = (pp << 3) + (pp << 2) + (pp << 1);
/*<       pdrdx = 8*pp+4*pp+2*pp+pp >*/
    pdrdx = (pp << 3) + (pp << 2) + (pp << 1) + pp;
/*<       ppgr  = 16*pp >*/
    ppgr = pp << 4;
/*<    >*/
    preordbe_(n, order, &ranmasks[1], &nown, &wrkord1[pisc], &wrkord1[pisd], &
	    wrkord1[pirc], &wrkord1[pird], &mynnodes, &dd, &myid, &lgblk, &
	    wrkord1[ppgr], wrkord0, comm);
/*<       wrkord2siz = 3*nown+mynnodes >*/
    wrkord2siz = nown * 3 + mynnodes;

/*     allocate(wrkord2(0:wrkord2siz-1),stat=is1) */
	wrkord2 = (integer*) malloc(wrkord2siz*sizeof(integer));
/*<       if(is1.ne.0) then >*/
    if (!wrkord2) {
/*<         print *,'memory allocation error' >*/
		printf("%d: memory allocation error", myid);
/*<         call mpi_abort(comm,0,ierr) >*/
		MPI_Abort(*comm, 0);
/*<       end if >*/
    }
/*<       piown  = 0 >*/
    piown = 0;
/*<       psloc  = 2*nown >*/
    psloc = nown << 1;
/*<       pslocx = 2*nown+mynnodes >*/
    pslocx = (nown << 1) + mynnodes;
/*<       ordvalsiz = 2*max(mynnodes,nown)*bnrhs >*/
    ordvalsiz = (max(mynnodes,nown) << 1) * bnrhs;

/*     allocate(ordvals(0:ordvalsiz-1),stat=is1) */
	ordvals = (doublereal*) malloc(ordvalsiz*sizeof(doublereal));
/*<       if(is1.ne.0) then >*/
    if (!ordvals) {
/*<         print *,'memory allocation error' >*/
		printf("%d: memory allocation error", myid);
/*<         call mpi_abort(comm,0,ierr) >*/
		MPI_Abort(*comm, 0);
/*<       end if >*/
    }
/*<       psv  = 0 >*/
    psv = 0;
/*<       prv  = mynnodes*bnrhs >*/
    prv = mynnodes * bnrhs;
/*<       psvx = 0 >*/
    psvx = 0;
/*<       prvx = nown*bnrhs >*/
    prvx = nown * bnrhs;
/*<    >*/
    preordbc_(n, order, &ranmasks[1], &nown, &wrkord1[pisc], &wrkord1[pisd], &
	    wrkord1[pirc], &wrkord1[pird], &wrkord2[piown], rowdist, &
	    mynnodes, &wrkord2[psloc], &ordvals[psv], &dd, &myid, &lgblk, &
	    wrkord1[ppgr], wrkord0, comm);
/*<    >*/
    preordxc_(n, &wrkord2[piown], rowdist, &mynnodes, &wrkord1[piscx], &
	    wrkord1[pisdx], &wrkord1[pircx], &wrkord1[pirdx], &nown, &wrkord2[
	    pslocx], wrkord0, &ordvals[psvx], &dd, &myid, &ordvals[prvx], 
	    comm);
/*<       do i=0,pp-1 >*/
    i__1 = pp - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
/*<         wrkord1(pdsc+i) = ishft(wrkord1(pisc+i),-1)*bnrhs >*/
	wrkord1[pdsc + i__] = lbit_shift(wrkord1[pisc + i__], -1) * 
		bnrhs;
/*<         wrkord1(pdsd+i) = ishft(wrkord1(pisd+i),-1)*bnrhs >*/
	wrkord1[pdsd + i__] = lbit_shift(wrkord1[pisd + i__], -1) * 
		bnrhs;
/*<         wrkord1(pdrc+i) = ishft(wrkord1(pirc+i),-1)*bnrhs >*/
	wrkord1[pdrc + i__] = lbit_shift(wrkord1[pirc + i__], -1) * 
		bnrhs;
/*<         wrkord1(pdrd+i) = ishft(wrkord1(pird+i),-1)*bnrhs >*/
	wrkord1[pdrd + i__] = lbit_shift(wrkord1[pird + i__], -1) * 
		bnrhs;
/*<         wrkord1(pdscx+i) = wrkord1(piscx+i)*bnrhs >*/
	wrkord1[pdscx + i__] = wrkord1[piscx + i__] * bnrhs;
/*<         wrkord1(pdsdx+i) = wrkord1(pisdx+i)*bnrhs >*/
	wrkord1[pdsdx + i__] = wrkord1[pisdx + i__] * bnrhs;
/*<         wrkord1(pdrcx+i) = wrkord1(pircx+i)*bnrhs >*/
	wrkord1[pdrcx + i__] = wrkord1[pircx + i__] * bnrhs;
/*<         wrkord1(pdrdx+i) = wrkord1(pirdx+i)*bnrhs >*/
	wrkord1[pdrdx + i__] = wrkord1[pirdx + i__] * bnrhs;
/*<       end do >*/
    }
/*<       do sr=1,nrhs-bnrhs+1,bnrhs   >*/
    i__1 = *nrhs - bnrhs + 1;
    i__2 = bnrhs;
    for (sr = 1; i__2 < 0 ? sr >= i__1 : sr <= i__1; sr += i__2) {
/*<    >*/
	reordb_(n, &rhso[sr * rhso_dim1], ldo, &bnrhs, &ty[ty_offset], &nown, 
		&wrkord1[pdsc], &wrkord1[pdsd], &wrkord1[pdrc], &wrkord1[pdrd]
		, &wrkord2[piown], &mynnodes, &wrkord2[psloc], &ordvals[psv], 
		&ordvals[prv], comm);
/*<       if(bnrhs.ne.1) then >*/
	if (bnrhs != 1) {
/*<    >*/
	    pfsolvem_(&lrud[1], &nsupnode, &sup[1], &lptrs[1], &linds[1], &
		    lvals[1], &tptrs[1], &tinds[1], &myid, &myidh, &dd, &
		    lgblk, n, &ty[ty_offset], &bnrhs, &dworkmj[rhsptr], &
		    dworkmj[uvecptr], &dworkmj[uvlptr], &dworkmj[recvptr], &
		    hvbtemp[uindptr], &maxvsize, &lrud[ns + 1], &hvbtemp[1], &
		    hvbsize, &lc[1], &dworkmj[1], &iptrs[1], comm);
/*<       else >*/
	} else {
/*<    >*/
	    pfsolve1_(&lrud[1], &nsupnode, &sup[1], &lptrs[1], &linds[1], &
		    lvals[1], &tptrs[1], &tinds[1], &myid, &myidh, &dd, &
		    lgblk, n, &ty[ty_offset], &dworkmj[rhsptr], &dworkmj[
		    uvecptr], &dworkmj[uvlptr], &dworkmj[recvptr], &hvbtemp[
		    uindptr], &maxvsize, &lrud[ns + 1], &hvbtemp[1], &hvbsize,
		     &lc[1], &dworkmj[1], &iptrs[1], comm);
/*<       end if >*/
	}
/*<       if(bnrhs.ne.1) then >*/
	if (bnrhs != 1) {
/*<    >*/
	    pbsolvem_(&lrud[1], &nsupnode, &sup[1], &lptrs[1], &linds[1], &
		    lvals[1], &tptrs[1], &tinds[1], &myid, &myidh, &myidv, &
		    dd, &lgblk, n, &ty[ty_offset], &bnrhs, &dworkmj[uvecptr], 
		    &dworkmj[recvptr], &dworkmj[rhsptr], &hvbtemp[uindptr], &
		    maxvsize, &lrud[ns + 1], &hvbtemp[1], &hvbsize, &lc[1], &
		    dworkmj[1], &iptrs[1], comm);
/*<       else >*/
	} else {
/*<    >*/
	    pbsolve1_(&lrud[1], &nsupnode, &sup[1], &lptrs[1], &linds[1], &
		    lvals[1], &tptrs[1], &tinds[1], &myid, &myidh, &myidv, &
		    dd, &lgblk, n, &ty[ty_offset], &dworkmj[uvecptr], &
		    dworkmj[recvptr], &dworkmj[rhsptr], &hvbtemp[uindptr], &
		    maxvsize, &lrud[ns + 1], &hvbtemp[1], &hvbsize, &lc[1], &
		    dworkmj[1], &iptrs[1], comm);
/*<       end if >*/
	}
/*<    >*/
	reordx_(n, &rhsc[sr * rhsc_dim1], ldc, &bnrhs, &ty[ty_offset], &nown, 
		&wrkord1[pdscx], &wrkord1[pdsdx], &wrkord1[pdrcx], &wrkord1[
		pdrdx], &wrkord2[piown], rowdist, &mynnodes, &wrkord2[pslocx],
		 wrkord0, &ordvals[psvx], &ordvals[prvx], &myid, comm);
/*<       end do  >*/
    }
/*<       rnr = nrhs-sr+1 >*/
    rnr = *nrhs - sr + 1;
/*<       if(rnr.ne.0) then >*/
    if (rnr != 0) {
/*<       do i=0,pp-1 >*/
	i__2 = pp - 1;
	for (i__ = 0; i__ <= i__2; ++i__) {
/*<         wrkord1(pdsc+i) = ishft(wrkord1(pisc+i),-1)*rnr >*/
	    wrkord1[pdsc + i__] = lbit_shift(wrkord1[pisc + i__], -1) 
		    * rnr;
/*<         wrkord1(pdsd+i) = ishft(wrkord1(pisd+i),-1)*rnr >*/
	    wrkord1[pdsd + i__] = lbit_shift(wrkord1[pisd + i__], -1) 
		    * rnr;
/*<         wrkord1(pdrc+i) = ishft(wrkord1(pirc+i),-1)*rnr >*/
	    wrkord1[pdrc + i__] = lbit_shift(wrkord1[pirc + i__], -1) 
		    * rnr;
/*<         wrkord1(pdrd+i) = ishft(wrkord1(pird+i),-1)*rnr >*/
	    wrkord1[pdrd + i__] = lbit_shift(wrkord1[pird + i__], -1) 
		    * rnr;
/*<         wrkord1(pdscx+i) = wrkord1(piscx+i)*rnr >*/
	    wrkord1[pdscx + i__] = wrkord1[piscx + i__] * rnr;
/*<         wrkord1(pdsdx+i) = wrkord1(pisdx+i)*rnr >*/
	    wrkord1[pdsdx + i__] = wrkord1[pisdx + i__] * rnr;
/*<         wrkord1(pdrcx+i) = wrkord1(pircx+i)*rnr >*/
	    wrkord1[pdrcx + i__] = wrkord1[pircx + i__] * rnr;
/*<         wrkord1(pdrdx+i) = wrkord1(pirdx+i)*rnr >*/
	    wrkord1[pdrdx + i__] = wrkord1[pirdx + i__] * rnr;
/*<       end do >*/
	}
/*<    >*/
	reordb_(n, &rhso[sr * rhso_dim1], ldo, &rnr, &ty[ty_offset], &nown, &
		wrkord1[pdsc], &wrkord1[pdsd], &wrkord1[pdrc], &wrkord1[pdrd],
		 &wrkord2[piown], &mynnodes, &wrkord2[psloc], &ordvals[psv], &
		ordvals[prv], comm);
/*<       if(rnr.ne.1) then >*/
	if (rnr != 1) {
/*<    >*/
	    pfsolvem_(&lrud[1], &nsupnode, &sup[1], &lptrs[1], &linds[1], &
		    lvals[1], &tptrs[1], &tinds[1], &myid, &myidh, &dd, &
		    lgblk, n, &ty[ty_offset], &rnr, &dworkmj[rhsptr], &
		    dworkmj[uvecptr], &dworkmj[uvlptr], &dworkmj[recvptr], &
		    hvbtemp[uindptr], &maxvsize, &lrud[ns + 1], &hvbtemp[1], &
		    hvbsize, &lc[1], &dworkmj[1], &iptrs[1], comm);
/*<       else >*/
	} else {
/*<    >*/
	    pfsolve1_(&lrud[1], &nsupnode, &sup[1], &lptrs[1], &linds[1], &
		    lvals[1], &tptrs[1], &tinds[1], &myid, &myidh, &dd, &
		    lgblk, n, &ty[ty_offset], &dworkmj[rhsptr], &dworkmj[
		    uvecptr], &dworkmj[uvlptr], &dworkmj[recvptr], &hvbtemp[
		    uindptr], &maxvsize, &lrud[ns + 1], &hvbtemp[1], &hvbsize,
		     &lc[1], &dworkmj[1], &iptrs[1], comm);
/*<       end if >*/
	}
/*<       if(rnr.ne.1) then >*/
	if (rnr != 1) {
/*<    >*/
	    pbsolvem_(&lrud[1], &nsupnode, &sup[1], &lptrs[1], &linds[1], &
		    lvals[1], &tptrs[1], &tinds[1], &myid, &myidh, &myidv, &
		    dd, &lgblk, n, &ty[ty_offset], &rnr, &dworkmj[uvecptr], &
		    dworkmj[recvptr], &dworkmj[rhsptr], &hvbtemp[uindptr], &
		    maxvsize, &lrud[ns + 1], &hvbtemp[1], &hvbsize, &lc[1], &
		    dworkmj[1], &iptrs[1], comm);
/*<       else >*/
	} else {
/*<    >*/
	    pbsolve1_(&lrud[1], &nsupnode, &sup[1], &lptrs[1], &linds[1], &
		    lvals[1], &tptrs[1], &tinds[1], &myid, &myidh, &myidv, &
		    dd, &lgblk, n, &ty[ty_offset], &dworkmj[uvecptr], &
		    dworkmj[recvptr], &dworkmj[rhsptr], &hvbtemp[uindptr], &
		    maxvsize, &lrud[ns + 1], &hvbtemp[1], &hvbsize, &lc[1], &
		    dworkmj[1], &iptrs[1], comm);
/*<       end if >*/
	}
/*<    >*/
	reordx_(n, &rhsc[sr * rhsc_dim1], ldc, &rnr, &ty[ty_offset], &nown, &
		wrkord1[pdscx], &wrkord1[pdsdx], &wrkord1[pdrcx], &wrkord1[
		pdrdx], &wrkord2[piown], rowdist, &mynnodes, &wrkord2[pslocx],
		 wrkord0, &ordvals[psvx], &ordvals[prvx], &myid, comm);
/*<       end if  >*/
    }

	/* reset pointer offset */
    ty += ty_offset;
    ++dworkmj;
    ++lrud;
    ++hvbtemp;

/*      deallocate(lrud) */
/*      deallocate(hvbtemp) */
/*      deallocate(ty) */
/*      deallocate(dworkmj) */
/*      deallocate(wrkord0) */
/*      deallocate(wrkord1) */
/*      deallocate(wrkord2) */
/*      deallocate(ordvals) */

	free(lrud);
	free(hvbtemp);
	free(ty);
	free(dworkmj);
	free(wrkord0);
	free(wrkord1);
	free(wrkord2);
	free(ordvals);

/*<       end >*/
    return 0;
} /* trisolve_ */

