/* $Id: pbsolve1.c,v 1.6 2005-01-15 08:18:28 paklein Exp $ */
/* pbsolve1.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

/* debugging */
#undef __DO_DEBUG__
/* #define __DO_DEBUG__ 1 */

#include "mpi.h"
#include "pspases_int.h"
#include <stdio.h>

/* Table of constant values */
static doublereal c_b3 = 1.;
static integer c__1 = 1;
static integer c__5 = 5;
static integer c__2 = 2;
static integer c__11 = 11;

/* /+***************************************************************************+/ */
/* /+                                                                           +/ */
/* /+   (C) Copyright IBM Corporation, 1997                                     +/ */
/* /+   (C) Copyright Regents of the University of Minnesota, 1997              +/ */
/* /+                                                                           +/ */
/* /+   pbsolve1.f                                                              +/ */
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
/* /+ $Id: pbsolve1.c,v 1.6 2005-01-15 08:18:28 paklein Exp $ +/ */
/* /+***************************************************************************+/ */

static integer lbit_shift(integer a, integer b) {
	return b >= 0 ? a << b : (integer)((uinteger)a >> -b);
};

static integer max(integer a, integer b) {
	return (a > b) ? a : b;
}

/*<    >*/
/* Subroutine */ int pbsolve1_(integer *mysnodes, integer *nsupnode, integer *
	sup, integer *lptrs, integer *linds, doublereal *lvals, integer *
	tptrs, integer *tinds, integer *myid, integer *myidh, integer *myidv, 
	integer *dd, integer *lgblk, integer *n, doublereal *rhsc, doublereal 
	*uvec, doublereal *recvec, doublereal *rhs, integer *uinds, integer *
	maxvsize, integer *lrud, integer *hvbndry, integer *hvbsize, integer *
	lc, doublereal *w, integer *iptrs, MPI_Comm *comm)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    integer npending;

    integer klbotsiz, i__;

    integer npendings;
    integer ih, ik, ip, ir, is, iv, nhb, kid, bhp, mid, bip, nvb, bvp, /* req[2], */ bvs;
    MPI_Request req[2];
    
    integer vfac, itag, indk, kbip, rrec, hvbp, nbrp, ierr, kbvs, nvbr, nvbt, 
	    tstf, myup, ldalb, bnode, kbipi;
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    integer level, khvbp, ivlim, hsize, psrcv, vsize;
    extern /* Subroutine */ int dtrsv_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen);
    integer rhsst, level2;

    integer bmaskh, bmaskv, bhsize, szflag;
    extern /* Subroutine */ int bsolve_(integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, integer *, integer *, integer *, doublereal *);
    integer lhsize, myleft, nuinds, bvsize, nvinds, rowcol, bhstrt, lvalst, 
	    vsizek, supbot, mydown, lvsize, krhsst, vsizer, bvstrt, suptop, 
	    supsiz;
    extern /* Subroutine */ int putrhs1_(integer *, integer *, integer *, 
	    doublereal *, doublereal *);
    integer istflag, partner, knvinds, /* mpistat[4], */ myright, lbotsiz, mymaskv, 
	    ksupbot, kvsizer, ksupptr, ksupsiz;
	MPI_Status mpistat;
	
    extern /* Subroutine */ int bgetrhs1_(integer *, integer *, integer *, 
	    doublereal *, doublereal *);

/*<       implicit none >*/
/*<       include 'mpif.h' >*/
/*<       integer N,dd,lgblk,myid,myidh,myidv,maxvsize >*/
/* -*- fortran -*- */

/* double precision functions */
/*<       integer nsupnode,hvbsize,comm >*/
/*<       integer lptrs(3,0:*), linds(*), iptrs(2,0:*) >*/
/*<       integer tptrs(3,0:*), tinds(*) >*/
/*<       integer sup(*), mysnodes(*) >*/
/*<       integer lrud(*), hvbndry(*) >*/
/*<       integer uinds(*), lc(*) >*/
/*<       double precision lvals(*), w(*) >*/
/*<       double precision rhsc(0:*),uvec(*),rhs(*) >*/
/*<       double precision recvec(*) >*/
/*<       integer lendp,loglendp,itype >*/
/*<       parameter(lendp=8,loglendp=3,itype=1) >*/
/*<       double precision one,ngone,zero >*/
/*<       parameter(zero=0.d0,one=1.d0,ngone=-1.d0) >*/
/*<       integer level,level2,rowcol,bmaskh,bmaskv,bmaskhk,bmaskvk >*/
/*<       integer myleft,myright,myup,mydown >*/
/*<       integer supptr,supbot >*/
/*<       integer nhb,nvb,nbrp,hvbp,ih,iv,ivlim,bhp,bvp,bvs >*/
/*<       integer bhstrt,bhsize,bvstrt,bvsize >*/
/*<       integer uvecst,i,j,is,nbrecv,k1,k2 >*/
/*<       integer lhsize,lvsize,hsize,vsize,vfac >*/
/*<       integer psrcv,mymaskv,bnode >*/
/*<       integer vsizek,nuinds,partner,lvalst,ldalb,itag >*/
/*<       integer ir,ip,ik,indk,krhsst,kbipi,mid,npending,npendings >*/
/*<       integer ksupbot,klbotsiz,kbip,kvsizer,kbvs >*/
/*<       integer trhs,ksupptr,khvbp,kid,knvinds,nvbr,nvbt >*/
/*<       integer istflag,rhsst,lbotsiz,suptop,nvinds >*/
/*<       integer vsizer,bip,tstf,ksupsiz,supsiz,szflag,rrec >*/
/*<       integer mpistat(MPI_STATUS_SIZE),req(2),ierr >*/
/*<       req(1) = MPI_REQUEST_NULL >*/
    /* Parameter adjustments */
    --iptrs;
    --w;
    --lc;
    --hvbndry;
    --lrud;
    --uinds;
    --rhs;
    --recvec;
    --uvec;
    --tinds;
    --tptrs;
    --lvals;
    --linds;
    --lptrs;
    --sup;
    --mysnodes;

    /* Function Body */
    req[0] = MPI_REQUEST_NULL;
/*<       req(2) = MPI_REQUEST_NULL >*/
    req[1] = MPI_REQUEST_NULL;
/*<       level  = 0 >*/
    level = 0;
/*<       rowcol = mod(dd,2) >*/
    rowcol = *dd % 2;
/*<       bmaskh = (dd+1)/2 >*/
    bmaskh = (*dd + 1) / 2;
/*<       bmaskv = dd-bmaskh >*/
    bmaskv = *dd - bmaskh;
/*<       hsize  = ishft(1,bmaskh) >*/
    hsize = lbit_shift((ftnlen)1, bmaskh);
/*<       vsize  = ishft(1,bmaskv) >*/
    vsize = lbit_shift((ftnlen)1, bmaskv);
/*<       lhsize = bmaskh >*/
    lhsize = bmaskh;
/*<       lvsize = bmaskv >*/
    lvsize = bmaskv;
/*<       bmaskh = hsize-1 >*/
    bmaskh = hsize - 1;
/*<       bmaskv = vsize-1 >*/
    bmaskv = vsize - 1;
/*<       nbrp = 1+(dd-1)*4 >*/
    nbrp = (*dd - 1 << 2) + 1;
/*<       khvbp = hvbsize+1 >*/
    khvbp = *hvbsize + 1;
/*<       kid      = mysnodes(1) >*/
    kid = mysnodes[1];
/*<       ksupptr  = tptrs(3,kid) >*/
    ksupptr = tptrs[kid * 3 + 3];
/*<       ksupbot  = sup(ksupptr) >*/
    ksupbot = sup[ksupptr];
/*<       ksupsiz  = sup(ksupptr+1) >*/
    ksupsiz = sup[ksupptr + 1];
/*<       khvbp    = hvbndry(khvbp-1) >*/
    khvbp = hvbndry[khvbp - 1];
/*<       kbvs     = hvbndry(khvbp+1) >*/
    kbvs = hvbndry[khvbp + 1];
/*<       klbotsiz = hvbndry(khvbp+4) >*/
    klbotsiz = hvbndry[khvbp + 4];
/*<       kbip     = kbvs-3-klbotsiz >*/
    kbip = kbvs - 3 - klbotsiz;
/*<       kvsizer  = hvbndry(kbvs-3) >*/
    kvsizer = hvbndry[kbvs - 3];
/*<       knvinds  = hvbndry(kbip-1) >*/
    knvinds = hvbndry[kbip - 1];
/*<       kbipi    = kbip+knvinds-kvsizer >*/
    kbipi = kbip + knvinds - kvsizer;
/*<       do is=1,nsupnode-1  >*/
    i__1 = *nsupnode - 1;
    for (is = 1; is <= i__1; ++is) {
/*<         hvbp = khvbp >*/
	hvbp = khvbp;
/*<         myleft = lrud(nbrp)  >*/
	myleft = lrud[nbrp];
/*<         myright= lrud(nbrp+1)  >*/
	myright = lrud[nbrp + 1];
/*<         myup   = lrud(nbrp+2)  >*/
	myup = lrud[nbrp + 2];
/*<         mydown = lrud(nbrp+3)  >*/
	mydown = lrud[nbrp + 3];
/*<         suptop  = kid >*/
	suptop = kid;
/*<         supbot  = ksupbot >*/
	supbot = ksupbot;
/*<         lbotsiz = klbotsiz >*/
	lbotsiz = klbotsiz;
/*<         supsiz  = ksupsiz >*/
	supsiz = ksupsiz;
/*<         bvs     = kbvs >*/
	bvs = kbvs;
/*<         vsizer  = kvsizer >*/
	vsizer = kvsizer;
/*<         nvinds  = knvinds >*/
	nvinds = knvinds;
/*<         bip     = kbipi >*/
	bip = kbipi;
/*<         bnode   = ishft(suptop,-lgblk) >*/
	bnode = lbit_shift(suptop, -(*lgblk));
/*<         mymaskv = iand(myidv,bmaskv) >*/
	mymaskv = *myidv & bmaskv;
/*<         istflag = hvbndry(hvbp+3) >*/
	istflag = hvbndry[hvbp + 3];
/*<         nhb     = hvbndry(hvbp+5) >*/
	nhb = hvbndry[hvbp + 5];
/*<         nvb     = hvbndry(bvs-1) >*/
	nvb = hvbndry[bvs - 1];
/*<         nvbr    = hvbndry(bvs-2) >*/
	nvbr = hvbndry[bvs - 2];
/*<         bhp     = hvbp+6+(nhb-1)*3       >*/
	bhp = hvbp + 6 + (nhb - 1) * 3;
/*<         nvbt    = nvb-nvbr >*/
	nvbt = nvb - nvbr;
/*<         bvp     = bvs+3*(nvbt-1)  >*/
	bvp = bvs + (nvbt - 1) * 3;
/*<         psrcv   = iand(bnode,bmaskv) >*/
	psrcv = bnode & bmaskv;
/*<         vfac    = hsize/vsize >*/
	vfac = hsize / vsize;
/*<         iv      = nvbt >*/
	iv = nvbt;
/*<         szflag  = supsiz-lbotsiz >*/
	szflag = supsiz - lbotsiz;
/*<         if((nvb.ne.0).and.(nhb.ne.0)) then >*/
	if (nvb != 0 && nhb != 0) {
/*<         do ih = nhb,1,-1 >*/
	    for (ih = nhb; ih >= 1; --ih) {
/*<           bhstrt = hvbndry(bhp) >*/
		bhstrt = hvbndry[bhp];
/*<           bhsize = hvbndry(bhp+1) >*/
		bhsize = hvbndry[bhp + 1];
/*<           rhsst = nvinds-vsizer+1 >*/
		rhsst = nvinds - vsizer + 1;
/*<           tstf  = istflag >*/
		tstf = istflag;
/*<           ivlim = max(iv-vfac+1,1) >*/
/* Computing MAX */
		i__2 = iv - vfac + 1;
		ivlim = max(i__2,1);
/*<           if(bhsize.ne.0) then >*/
		if (bhsize != 0) {
/*<           ldalb  = lptrs(2,bhstrt) >*/
		    ldalb = lptrs[bhstrt * 3 + 2];
/*<           lvalst = lptrs(1,bhstrt)+ldalb-vsizer >*/
		    lvalst = lptrs[bhstrt * 3 + 1] + ldalb - vsizer;
/*<           do i=1,bhsize >*/
		    i__2 = bhsize;
		    for (i__ = 1; i__ <= i__2; ++i__) {
/*<             recvec(i) = zero >*/
			recvec[i__] = 0.;
/*<           end do >*/
		    }
/*<           if(vsizer.ne.0) then >*/
		    if (vsizer != 0) {
/*<    >*/
			dgemv_("t", &vsizer, &bhsize, &c_b3, &lvals[lvalst], &
				ldalb, &rhs[rhsst], &c__1, &c_b3, &recvec[1], 
				&c__1, (ftnlen)1);
/*<           end if >*/
		    }
/*<           if (szflag.ne.0) then >*/
		    if (szflag != 0) {
/*<    >*/
			alltoonev_(&recvec[1], &uvec[1], &bhsize, &psrcv, &
				lvsize, &mymaskv, myid, comm);
/*<           else >*/
		    } else {
/*<             szflag = 1 >*/
			szflag = 1;
/*<           end if >*/
		    }
/*<           do iv = iv,ivlim,-1 >*/
		    i__2 = ivlim;
		    for (iv = iv; iv >= i__2; --iv) {
/*<             bvstrt = hvbndry(bvp) >*/
			bvstrt = hvbndry[bvp];
/*<             bvsize = hvbndry(bvp+1) >*/
			bvsize = hvbndry[bvp + 1];
/*<             lvalst = lvalst-bvsize >*/
			lvalst -= bvsize;
/*<             rhsst  = rhsst-bvsize >*/
			rhsst -= bvsize;
/*<             bip    = bip-bvsize >*/
			bip -= bvsize;
/*<             itag   = MPI_ANY_TAG >*/
			itag = MPI_ANY_TAG;
/*<             if(bvsize.ne.0) then >*/
			if (bvsize != 0) {
/*<             if(tstf.eq.1 .or. tstf.eq.2) then >*/
			    if (tstf == 1 || tstf == 2) {
/*<               if(bhstrt.eq.bvstrt) then  >*/
				if (bhstrt == bvstrt) {
/*<               call bgetrhs1(N,hvbndry(bip),bvsize,rhs(rhsst),rhsc) >*/
				    bgetrhs1_(n, &hvbndry[bip], &bvsize, &rhs[
					    rhsst], rhsc);
/*<               do i=0,bvsize-1 >*/
				    i__3 = bvsize - 1;
				    for (i__ = 0; i__ <= i__3; ++i__) {
/*<                 rhs(rhsst+i) = rhs(rhsst+i) - recvec(i+1) >*/
					rhs[rhsst + i__] -= recvec[i__ + 1];
/*<               end do >*/
				    }
/*<    >*/
				    dtrsv_("l", "t", "n", &bvsize, &lvals[
					    lvalst], &ldalb, &rhs[rhsst], &
					    c__1, (ftnlen)1, (ftnlen)1, (
					    ftnlen)1);
/*<               itag = myid >*/
				    itag = *myid;
/*<    >*/
				    i__3 = bvsize << 3;

/*	call mpi_isend(rhs(rhsst),bvsize*lendp,MPI_BYTE,myleft,itag,comm,req(1),ierr) */
				    myMPI_Isend(&rhs[rhsst], i__3, MPI_BYTE, myleft, itag, *comm, req);

/*<               call putrhs1(N,hvbndry(bip),bvsize,rhs(rhsst),rhsc) >*/
				    putrhs1_(n, &hvbndry[bip], &bvsize, &rhs[rhsst], rhsc);
/*<               call mpi_wait(req(1),mpistat,ierr) >*/
				    MPI_Wait(req, &mpistat);
/*<               else >*/
				} else {
/*<                 npendings = 0 >*/
				    npendings = 0;
/*<    >*/
				    i__3 = bvsize << 3;
				    MPI_Irecv(&rhs[rhsst], i__3, MPI_BYTE, myright, itag, *comm, req);

/*<                 call mpi_wait(req(1),mpistat,ierr) >*/
				    MPI_Wait(req, &mpistat);
/*<                 itag = mpistat(MPI_TAG) >*/
				    itag = mpistat.MPI_TAG;
/*<                 if(itag.ne.myleft) then >*/
				    if (itag != myleft) {
/*<    >*/
					i__3 = bvsize << 3;
					myMPI_Isend(&rhs[rhsst], i__3, MPI_BYTE, myleft, itag, *comm, req);
/*<                   npendings = 1 >*/
					npendings = 1;
/*<                 end if >*/
				    }
/*<    >*/
				    dgemv_("t", &bvsize, &bhsize, &c_b3, &
					    lvals[lvalst], &ldalb, &rhs[rhsst]
					    , &c__1, &c_b3, &recvec[1], &c__1,
					     (ftnlen)1);
/*<                 if (myup.ne.myid) then >*/
				    if (myup != *myid) {
/*<    >*/
					i__3 = bhsize << 3;
					myMPI_Isend(&recvec[1], i__3, MPI_BYTE, myup, 1, *comm, &req[1]);
/*<                            npendings = npendings + 1 >*/
					++npendings;
/*<                 end if >*/
				    }
/*<                 do while (npendings.gt.0) >*/
				    while(npendings > 0) {
/*<                   call mpi_waitany(2,req,mid,mpistat,ierr) >*/
					MPI_Waitany(2, req, &mid, &mpistat);
					mid++; /* completed request indexed from 1 in FORTRAN */
					
/*<                   npendings = npendings-1 >*/
					--npendings;
/*<                 end do >*/
				    }
/*<               end if >*/
				}
/*<               tstf = 0 >*/
				tstf = 0;
/*<             else   >*/
			    } else {
/*<               if(bhstrt.eq.bvstrt) then  >*/
				if (bhstrt == bvstrt) {
/*<                 if(mydown.ne.myid) then >*/
				    if (mydown != *myid) {
/*<    >*/
					i__3 = bhsize << 3;
					MPI_Irecv(&recvec[1], i__3, MPI_BYTE, mydown, 1, *comm, &req[1]);
/*<                 end if >*/
				    }
/*<                 call bgetrhs1(N,hvbndry(bip),bvsize,rhs(rhsst),rhsc) >*/
				    bgetrhs1_(n, &hvbndry[bip], &bvsize, &rhs[
					    rhsst], rhsc);
/*<                 if(mydown.ne.myid) then >*/
				    if (mydown != *myid) {
/*<                   call mpi_wait(req(2),mpistat,ierr) >*/
					MPI_Wait(&req[1], &mpistat);
/*<                 end if >*/
				    }
/*<                 do i=0,bvsize-1 >*/
				    i__3 = bvsize - 1;
				    for (i__ = 0; i__ <= i__3; ++i__) {
/*<                   rhs(rhsst+i) = rhs(rhsst+i) - recvec(i+1) >*/
					rhs[rhsst + i__] -= recvec[i__ + 1];
/*<                 end do >*/
				    }
/*<    >*/
				    dtrsv_("l", "t", "n", &bvsize, &lvals[
					    lvalst], &ldalb, &rhs[rhsst], &
					    c__1, (ftnlen)1, (ftnlen)1, (
					    ftnlen)1);
/*<                 itag = myid >*/
				    itag = *myid;
/*<    >*/
				    i__3 = bvsize << 3;
				    myMPI_Isend(&rhs[rhsst], i__3, MPI_BYTE, myleft, itag, *comm, req);

/*<                 call putrhs1(N,hvbndry(bip),bvsize,rhs(rhsst),rhsc) >*/
				    putrhs1_(n, &hvbndry[bip], &bvsize, &rhs[
					    rhsst], rhsc);
/*<                 call mpi_wait(req(1),mpistat,ierr) >*/
				    MPI_Wait(req, &mpistat);
/*<               else    >*/
				} else {
/*<                 if(bvstrt.gt.bhstrt) then >*/
				    if (bvstrt > bhstrt) {
/*<                   npendings = 0 >*/
					npendings = 0;
/*<                   rrec = 0 >*/
					rrec = 0;
/*<    >*/
					i__3 = bvsize << 3;
					MPI_Irecv(&rhs[rhsst], i__3, MPI_BYTE, myright, itag, *comm, req);
/*<                   npending = 1 >*/
					npending = 1;
/*<                   if(mydown.ne.myid) then >*/
					if (mydown != *myid) {
/*<    >*/
					    i__3 = bhsize << 3;
					    MPI_Irecv(&recvec[1], i__3, MPI_BYTE, mydown, 1, *comm, &req[1]);
/*<                     npending = 2 >*/
					    npending = 2;
/*<                   end if >*/
					}
/*<                   do while(npending.gt.0) >*/
					while(npending > 0) {
/*<                     call mpi_waitany(2,req,mid,mpistat,ierr) >*/

					    MPI_Waitany(2, req, &mid, &mpistat);
					    mid++; /* completed requests numbered from 1 in FORTRAN */

/*<                     if(mid.eq.1 .and. rrec.eq.0) then >*/
					    if (mid == 1 && rrec == 0) {
/*<                       itag = mpistat(MPI_TAG) >*/
			  itag = mpistat.MPI_TAG;
/*<                       if(itag.ne.myleft) then >*/
			  if (itag != myleft) {
/*<    >*/
			      i__3 = bvsize << 3;
			      myMPI_Isend(&rhs[rhsst], i__3, MPI_BYTE, myleft, itag, *comm, req);
/*<                         npending = npending+1 >*/
			      ++npending;
/*<                       end if >*/
			  }
/*<                       rrec = 1 >*/
			  rrec = 1;
/*<                     end if >*/
					    }
/*<                     npending = npending-1 >*/
					    --npending;
/*<                   end do >*/
					}
/*<    >*/
					dgemv_("t", &bvsize, &bhsize, &c_b3, &
						lvals[lvalst], &ldalb, &rhs[
						rhsst], &c__1, &c_b3, &recvec[
						1], &c__1, (ftnlen)1);
/*<                   if(myup.ne.myid) then >*/
					if (myup != *myid) {
/*<    >*/
					    i__3 = bhsize << 3;
					    myMPI_Isend(&recvec[1], i__3, MPI_BYTE, myup, 1, *comm, &req[1]);
/*<                     call mpi_wait(req(2),mpistat,ierr) >*/
					    MPI_Wait(&req[1], &mpistat);
/*<                   end if >*/
					}
/*<                 else  >*/
				    } else {
/*<    >*/
					i__3 = bvsize << 3;
					MPI_Irecv(&rhs[rhsst], i__3, MPI_BYTE, myright, itag, *comm, req);
/*<                   call mpi_wait(req(1),mpistat,ierr) >*/
					MPI_Wait(req, &mpistat);
/*<                   itag = mpistat(MPI_TAG) >*/
					itag = mpistat.MPI_TAG;
/*<                   if(itag.ne.myleft) then >*/
					if (itag != myleft) {
/*<    >*/
					    i__3 = bvsize << 3;
					    myMPI_Isend(&rhs[rhsst], i__3, MPI_BYTE, myleft, itag, *comm, req);
/*<                     call mpi_wait(req(1),mpistat,ierr) >*/
					    MPI_Wait(req, &mpistat);
/*<                   end if >*/
					}
/*<                 end if  >*/
				    }
/*<               end if  >*/
				}
/*<             end if  >*/
			    }
/*<             else  >*/
			} else {
/*<               if(myup.ne.myid .and. bhstrt.lt.bvstrt) then >*/
			    if (myup != *myid && bhstrt < bvstrt) {
/*<                 if(tstf.eq.2) then >*/
				if (tstf == 2) {
/*<                   tstf = 0 >*/
				    tstf = 0;
/*<                 else >*/
				} else {
/*<    >*/
				    i__3 = bhsize << 3;
				    MPI_Irecv(&recvec[1], i__3, MPI_BYTE, mydown, 1, *comm, &req[1]);
/*<                   call mpi_wait(req(2),mpistat,ierr) >*/
				    MPI_Wait(&req[1], &mpistat);
/*<                 end if >*/
				}
/*<    >*/
				i__3 = bhsize << 3;
				myMPI_Isend(&recvec[1], i__3, MPI_BYTE, myup, 1, *comm, &req[1]);
/*<                 call mpi_wait(req(2),mpistat,ierr) >*/
				MPI_Wait(&req[1], &mpistat);
/*<               end if >*/
			    }
/*<             end if >*/
			}
/*<             vsizer = vsizer+bvsize >*/
			vsizer += bvsize;
/*<             bvp = bvp-3 >*/
			bvp += -3;
/*<           end do   >*/
		    }
/*<           else >*/
		} else {
/*<             do iv = iv,ivlim,-1 >*/
		    i__2 = ivlim;
		    for (iv = iv; iv >= i__2; --iv) {
/*<               bvsize = hvbndry(bvp+1) >*/
			bvsize = hvbndry[bvp + 1];
/*<               lvalst = lvalst-bvsize >*/
			lvalst -= bvsize;
/*<               rhsst  = rhsst-bvsize >*/
			rhsst -= bvsize;
/*<               bip    = bip-bvsize >*/
			bip -= bvsize;
/*<               itag   = MPI_ANY_TAG >*/
			itag = MPI_ANY_TAG;
/*<               if(bvsize.ne.0) then >*/
			if (bvsize != 0) {
/*<    >*/
			    i__3 = bvsize << 3;
			    MPI_Irecv(&rhs[rhsst], i__3, MPI_BYTE, myright, itag, *comm, req);
/*<                 call mpi_wait(req(1),mpistat,ierr) >*/
			    MPI_Wait(req, &mpistat);
/*<                 itag = mpistat(MPI_TAG) >*/
			    itag = mpistat.MPI_TAG;
/*<                 if(itag.ne.myleft) then >*/
			    if (itag != myleft) {
/*<    >*/
				i__3 = bvsize << 3;
				myMPI_Isend(&rhs[rhsst], i__3, MPI_BYTE, myleft, itag, *comm, req);
/*<                   call mpi_wait(req(1),mpistat,ierr) >*/
				MPI_Wait(req, &mpistat);
/*<                 end if >*/
			    }
/*<               end if >*/
			}
/*<               vsizer = vsizer+bvsize >*/
			vsizer += bvsize;
/*<               bvp = bvp-3 >*/
			bvp += -3;
/*<             end do >*/
		    }
/*<           end if >*/
		}
/*<           bhp = bhp-3 >*/
		bhp += -3;
/*<         end do  >*/
	    }
/*<         do iv = ivlim-1,1,-1 >*/
	    for (iv = ivlim - 1; iv >= 1; --iv) {
/*<           bvsize = hvbndry(bvp+1) >*/
		bvsize = hvbndry[bvp + 1];
/*<           if (bvsize.ne.0) then >*/
		if (bvsize != 0) {
/*<             rhsst  = rhsst-bvsize >*/
		    rhsst -= bvsize;
/*<             bip    = bip-bvsize >*/
		    bip -= bvsize;
/*<             itag   = MPI_ANY_TAG >*/
		    itag = MPI_ANY_TAG;
/*<    >*/
		    i__2 = bvsize << 3;
		    MPI_Irecv(&rhs[rhsst], i__2, MPI_BYTE, myright, itag, *comm, req);
/*<             call mpi_wait(req(1),mpistat,ierr) >*/
		    MPI_Wait(req, &mpistat);
/*<             itag = mpistat(MPI_TAG) >*/
		    itag = mpistat.MPI_TAG;
/*<             if(itag.ne.myleft) then >*/
		    if (itag != myleft) {
/*<    >*/
			i__2 = bvsize << 3;
			myMPI_Isend(&rhs[rhsst], i__2, MPI_BYTE, myleft, itag, *comm, req);
/*<               call mpi_wait(req(1),mpistat,ierr) >*/
			MPI_Wait(req, &mpistat);
/*<             end if >*/
		    }
/*<           end if >*/
		}
/*<           bvp = bvp-3 >*/
		bvp += -3;
/*<         end do >*/
	    }
/*<         else if((nvb.eq.0).and.(nhb.ne.0)) then >*/
	} else if (nvb == 0 && nhb != 0) {
/*<           if(szflag.ne.0) then  >*/
	    if (szflag != 0) {
/*<             bhsize = hvbndry(bhp+1) >*/
		bhsize = hvbndry[bhp + 1];
/*<             if(bhsize.ne.0) then >*/
		if (bhsize != 0) {
/*<             do i=1,bhsize >*/
		    i__2 = bhsize;
		    for (i__ = 1; i__ <= i__2; ++i__) {
/*<               uvec(i) = zero >*/
			uvec[i__] = 0.;
/*<             end do >*/
		    }
/*<    >*/
		    alltoonev_(&uvec[1], &recvec[1], &bhsize, &psrcv, &lvsize,
			     &mymaskv, myid, comm);
/*<             end if >*/
		}
/*<             bhp = bhp-3 >*/
		bhp += -3;
/*<           end if >*/
	    }
/*<           do ih=nhb-1,1,-1 >*/
	    for (ih = nhb - 1; ih >= 1; --ih) {
/*<             bhsize = hvbndry(bhp+1) >*/
		bhsize = hvbndry[bhp + 1];
/*<             if(bhsize.ne.0) then >*/
		if (bhsize != 0) {
/*<             do i=1,bhsize >*/
		    i__2 = bhsize;
		    for (i__ = 1; i__ <= i__2; ++i__) {
/*<               uvec(i) = zero >*/
			uvec[i__] = 0.;
/*<             end do >*/
		    }
/*<    >*/
		    alltoonev_(&uvec[1], &recvec[1], &bhsize, &psrcv, &lvsize,
			     &mymaskv, myid, comm);
/*<             end if >*/
		}
/*<             bhp = bhp-3 >*/
		bhp += -3;
/*<           end do >*/
	    }
/*<         else if((nhb.eq.0).and.(nvb.ne.nvbr)) then >*/
	} else if (nhb == 0 && nvb != nvbr) {
/*<           rhsst = nvinds-vsizer+1 >*/
	    rhsst = nvinds - vsizer + 1;
/*<           do iv = nvbt,1,-1 >*/
	    for (iv = nvbt; iv >= 1; --iv) {
/*<             bvsize = hvbndry(bvp+1) >*/
		bvsize = hvbndry[bvp + 1];
/*<             if(bvsize.ne.0) then >*/
		if (bvsize != 0) {
/*<               rhsst  = rhsst-bvsize >*/
		    rhsst -= bvsize;
/*<               bip    = bip-bvsize >*/
		    bip -= bvsize;
/*<               itag   = MPI_ANY_TAG >*/
		    itag = MPI_ANY_TAG;
/*<    >*/
		    i__2 = bvsize << 3;
		    MPI_Irecv(&rhs[rhsst], i__2, MPI_BYTE, myright, itag, *comm, req);
/*<               call mpi_wait(req(1),mpistat,ierr) >*/
		    MPI_Wait(req, &mpistat);
/*<               itag = mpistat(MPI_TAG) >*/
		    itag = mpistat.MPI_TAG;
/*<               if(itag.ne.myleft) then >*/
		    if (itag != myleft) {
/*<    >*/
			i__2 = bvsize << 3;
			myMPI_Isend(&rhs[rhsst], i__2, MPI_BYTE, myleft, itag, *comm, req);
/*<                 call mpi_wait(req(1),mpistat,ierr) >*/
			MPI_Wait(req, &mpistat);
/*<               end if >*/
		    }
/*<             end if >*/
		}
/*<             bvp = bvp-3 >*/
		bvp += -3;
/*<           end do >*/
	    }
/*<         end if >*/
	}
/*<         vsizek = vsize >*/
	vsizek = vsize;
/*<         if(tptrs(2,supbot).ne.1) then  >*/
	if (tptrs[supbot * 3 + 2] != 1) {
/*<           level2 = ishft(dd-level-1,-1) >*/
	    level2 = lbit_shift(*dd - level - 1, (ftnlen)-1);
/*<           if(rowcol.eq.1) then >*/
	    if (rowcol == 1) {
/*<             bmaskh = ieor(bmaskh,ishft(1,level2)) >*/
		bmaskh ^= lbit_shift((ftnlen)1, level2);
/*<             hsize  = ishft(hsize,-1) >*/
		hsize = lbit_shift(hsize, (ftnlen)-1);
/*<             lhsize = lhsize-1 >*/
		--lhsize;
/*<           else >*/
	    } else {
/*<             bmaskv = ieor(bmaskv,ishft(1,level2)) >*/
		bmaskv ^= lbit_shift((ftnlen)1, level2);
/*<             vsize  = ishft(vsize,-1) >*/
		vsize = lbit_shift(vsize, (ftnlen)-1);
/*<             lvsize = lvsize-1 >*/
		--lvsize;
/*<           end if >*/
	    }
/*<           level = level+1 >*/
	    ++level;
/*<           rowcol = 1-rowcol >*/
	    rowcol = 1 - rowcol;
/*<           nbrp = nbrp-4   >*/
	    nbrp += -4;
/*<         end if >*/
	}
/*<         bip = kbip >*/
	bip = kbip;
/*<         kid = mysnodes(is+1) >*/
	kid = mysnodes[is + 1];
/*<         if(is.ne.nsupnode-1) then >*/
	if (is != *nsupnode - 1) {
/*<           khvbp    = hvbndry(khvbp-1) >*/
	    khvbp = hvbndry[khvbp - 1];
/*<           ksupptr  = tptrs(3,kid) >*/
	    ksupptr = tptrs[kid * 3 + 3];
/*<           ksupbot  = sup(ksupptr) >*/
	    ksupbot = sup[ksupptr];
/*<           ksupsiz  = sup(ksupptr+1) >*/
	    ksupsiz = sup[ksupptr + 1];
/*<           kbvs     = hvbndry(khvbp+1) >*/
	    kbvs = hvbndry[khvbp + 1];
/*<           klbotsiz = hvbndry(khvbp+4) >*/
	    klbotsiz = hvbndry[khvbp + 4];
/*<           kbip     = kbvs-3-klbotsiz >*/
	    kbip = kbvs - 3 - klbotsiz;
/*<           kvsizer  = hvbndry(kbvs-3) >*/
	    kvsizer = hvbndry[kbvs - 3];
/*<           knvinds  = hvbndry(kbip-1) >*/
	    knvinds = hvbndry[kbip - 1];
/*<           kbipi    = kbip+knvinds-kvsizer >*/
	    kbipi = kbip + knvinds - kvsizer;
/*<           if(vsize.ne.vsizek) then   >*/
	    if (vsize != vsizek) {
/*<             partner = ieor(myid,ishft(1,dd-level)) >*/
		partner = *myid ^ lbit_shift((ftnlen)1, *dd - level);
/*<             if(partner.gt.myid) then >*/
		if (partner > *myid) {
/*<    >*/
		    MPI_Send(&hvbndry[bip], nvinds, MPI_INT, partner, 1, *comm);
/*<    >*/
		    MPI_Recv(&uinds[1], *maxvsize, MPI_INT, partner, 1, *comm, &mpistat);
/*<             else >*/
		} else {
/*<    >*/
		    MPI_Recv(&uinds[1], *maxvsize, MPI_INT, partner, 1, *comm, &mpistat);
/*<    >*/
		    MPI_Send(&hvbndry[bip], nvinds, MPI_INT, partner, 1, *comm);
/*<             end if >*/
		}
/*<             call mpi_get_count(mpistat,MPI_INTEGER,nuinds,ierr) >*/
		myMPI_Get_count(&mpistat, MPI_INT, &nuinds);
/*<             if(partner.gt.myid) then >*/
		if (partner > *myid) {
/*<    >*/
		    i__2 = nvinds << 3;
		    MPI_Send(&rhs[1], i__2, MPI_BYTE, partner, 1, *comm);
/*<    >*/
		    i__2 = *maxvsize << 3;
		    MPI_Recv(&recvec[1], i__2, MPI_BYTE, partner, 1, *comm, &mpistat);
/*<             else >*/
		} else {
/*<    >*/
		    i__2 = *maxvsize << 3;
		    MPI_Recv(&recvec[1], i__2, MPI_BYTE, partner, 1, *comm, &mpistat);
/*<    >*/
		    i__2 = nvinds << 3;
		    MPI_Send(&rhs[1], i__2, MPI_BYTE, partner, 1, *comm);
/*<             end if >*/
		}
/*<             ir = 1 >*/
		ir = 1;
/*<             ip = 1 >*/
		ip = 1;
/*<             do ik=1,kvsizer >*/
		i__2 = kvsizer;
		for (ik = 1; ik <= i__2; ++ik) {
/*<               indk = hvbndry(kbipi+ik-1) >*/
		    indk = hvbndry[kbipi + ik - 1];
/*<               do while(uinds(ir).lt.indk .and. ir.le.nuinds) >*/
		    while(uinds[ir] < indk && ir <= nuinds) {
/*<                 ir = ir+1 >*/
			++ir;
/*<               end do >*/
		    }
/*<               do while(hvbndry(bip+ip-1).lt.indk .and. ip.le.nvinds) >*/
		    while(hvbndry[bip + ip - 1] < indk && ip <= nvinds) {
/*<                 ip = ip+1 >*/
			++ip;
/*<               end do >*/
		    }
/*<               if(ip.le.nvinds .and. indk.eq.hvbndry(bip+ip-1)) then >*/
		    if (ip <= nvinds && indk == hvbndry[bip + ip - 1]) {
/*<                 uvec(ik) = rhs(ip) >*/
			uvec[ik] = rhs[ip];
/*<               end if >*/
		    }
/*<               if(ir.le.nuinds .and. indk.eq.uinds(ir)) then >*/
		    if (ir <= nuinds && indk == uinds[ir]) {
/*<                 uvec(ik) = recvec(ir) >*/
			uvec[ik] = recvec[ir];
/*<               end if >*/
		    }
/*<             end do >*/
		}
/*<           else  >*/
	    } else {
/*<             ip = 1 >*/
		ip = 1;
/*<             do ik=1,kvsizer >*/
		i__2 = kvsizer;
		for (ik = 1; ik <= i__2; ++ik) {
/*<               indk = hvbndry(kbipi+ik-1) >*/
		    indk = hvbndry[kbipi + ik - 1];
/*<               do while(hvbndry(bip+ip-1).lt.indk .and. ip.le.nvinds) >*/
		    while(hvbndry[bip + ip - 1] < indk && ip <= nvinds) {
/*<                 ip = ip+1 >*/
			++ip;
/*<               end do >*/
		    }
/*<               if(indk.eq.hvbndry(bip+ip-1)) then >*/
		    if (indk == hvbndry[bip + ip - 1]) {
/*<                 uvec(ik) = rhs(ip) >*/
			uvec[ik] = rhs[ip];
/*<               end if >*/
		    }
/*<             end do >*/
		}
/*<           end if >*/
	    }
/*<           krhsst = knvinds-kvsizer >*/
	    krhsst = knvinds - kvsizer;
/*<           do i=1,krhsst >*/
	    i__2 = krhsst;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/*<             rhs(i) = zero >*/
		rhs[i__] = 0.;
/*<           end do >*/
	    }
/*<           do i=1,kvsizer >*/
	    i__2 = kvsizer;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/*<             rhs(krhsst+i) = uvec(i) >*/
		rhs[krhsst + i__] = uvec[i__];
/*<           end do >*/
	    }
/*<         else  >*/
	} else {
/*<           i = tptrs(3,kid) >*/
	    i__ = tptrs[kid * 3 + 3];
/*<           knvinds = lptrs(2,sup(i)) >*/
	    knvinds = lptrs[sup[i__] * 3 + 2];
/*<           kvsizer = lptrs(2,kid)-1 >*/
	    kvsizer = lptrs[kid * 3 + 2] - 1;
/*<           kbipi   = lptrs(3,kid)+1 >*/
	    kbipi = lptrs[kid * 3 + 3] + 1;
/*<           krhsst = knvinds-kvsizer >*/
	    krhsst = knvinds - kvsizer;
/*<           ip = 1 >*/
	    ip = 1;
/*<           do ik=1,kvsizer >*/
	    i__2 = kvsizer;
	    for (ik = 1; ik <= i__2; ++ik) {
/*<             indk = linds(kbipi+ik-1) >*/
		indk = linds[kbipi + ik - 1];
/*<             do while(hvbndry(bip+ip-1).lt.indk .and. ip.le.nvinds) >*/
		while(hvbndry[bip + ip - 1] < indk && ip <= nvinds) {
/*<               ip = ip+1 >*/
		    ++ip;
/*<             end do >*/
		}
/*<             if(ip.le.nvinds .and. indk.eq.hvbndry(bip+ip-1)) then >*/
		if (ip <= nvinds && indk == hvbndry[bip + ip - 1]) {
/*<               w(krhsst+ik) = rhs(ip) >*/
		    w[krhsst + ik] = rhs[ip];
/*<             end if >*/
		}
/*<           end do >*/
	    }
/*<         end if >*/
	}
/*<       end do   >*/
    }
/*<    >*/
    bsolve_(n, &lvals[1], &linds[1], &lptrs[1], &tinds[1], &tptrs[1], &sup[1],
	     rhsc, &c__1, &mysnodes[*nsupnode], &lc[1], &iptrs[1], &w[1]);
/*<       end >*/
    return 0;
} /* pbsolve1_ */

