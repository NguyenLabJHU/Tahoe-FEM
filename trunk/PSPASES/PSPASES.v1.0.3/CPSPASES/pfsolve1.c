/* $Id: pfsolve1.c,v 1.2 2005-01-04 18:19:34 paklein Exp $ */
/* pfsolve1.f -- translated by f2c (version 20030320).
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
static integer c__5 = 5;
static integer c__11 = 11;
static integer c__2 = 2;
static integer c__3 = 3;
static doublereal c_b58 = -1.;
static doublereal c_b60 = 1.;

/* /+***************************************************************************+/ */
/* /+                                                                           +/ */
/* /+   (C) Copyright IBM Corporation, 1997                                     +/ */
/* /+   (C) Copyright Regents of the University of Minnesota, 1997              +/ */
/* /+                                                                           +/ */
/* /+   pfsolve1.f                                                              +/ */
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
/* /+ $Id: pfsolve1.c,v 1.2 2005-01-04 18:19:34 paklein Exp $ +/ */
/* /+***************************************************************************+/ */

static integer lbit_shift(integer a, integer b) {
	return b >= 0 ? a << b : (integer)((uinteger)a >> -b);
};

/*<    >*/
/* Subroutine */ int pfsolve1_(integer *mysnodes, integer *nsupnode, integer *
	sup, integer *lptrs, integer *linds, doublereal *lvals, integer *
	tptrs, integer *tinds, integer *myid, integer *myidh, integer *dd, 
	integer *lgblk, integer *n, doublereal *rhsc, doublereal *rhs, 
	doublereal *uvec, doublereal *uvl, doublereal *recvec, integer *uinds,
	 integer *maxvsize, integer *lrud, integer *hvbndry, integer *hvbsize,
	 integer *lc, doublereal *w, integer *iptrs, MPI_Comm *comm)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    integer diffbits, npending;
    integer mymaskhk, recvsize, i__, j;
    integer k1, k2;
    integer npendings, nb, ih;
    extern /* Subroutine */ int extend_op1__(doublereal *, integer *, 
	    doublereal *, integer *, integer *, integer *);
    integer npendings2, is, iv, nhb, kid, bhp, mid, bip, nvb, bvp, /* req[5], */
	    bvs, bip1;
	MPI_Request req[5];
    integer itag, kidl, iaml, pkid, indh, kidr, hvbp, nbrp, ierr, root, myup, 
	    ldalb, bhend, bnode, flagr, flags, bvend, ldauk, lpkid;
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    integer level, ldaup, pkido, lppar, hsize, usize;
    extern /* Subroutine */ int x10dad1_(doublereal *, integer *, doublereal *
	    , integer *, doublereal *, integer *, integer *, integer *, 
	    integer *);
    integer rhsst;
    extern /* Subroutine */ int dtrsv_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen);
    integer uvlst, level2;
    integer iamkid, bmaskh, iampar, nbrecv, cbsize, levelk, bmaskv, bhsize;
    extern /* Subroutine */ int fsolve_(integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, integer *, integer *, integer *, doublereal *);
    integer ntinds, myleft, nvinds, bvsize, rowcol, bhstrt, lvalst, uvecst, 
	    supbot, mydown, vsizer, bvstrt, suptop, supptr, supsiz;
    extern /* Subroutine */ int putrhs1_(integer *, integer *, integer *, 
	    doublereal *, doublereal *);
    integer bmaskhk, istflag, bmaskvk, /* statall[20], */ /* was [4][5] */ 
	    pparbot, cmpsize, partner, mymaskh, msizedp, /*mpistat[4],*/ myright, 
	    lbotsiz;
	MPI_Status mpistat, statall[5];
    extern /* Subroutine */ int bgetrhs1_(integer *, integer *, integer *, 
	    doublereal *, doublereal *);

/*<       implicit none >*/
/*<       include 'mpif.h' >*/
/*<       integer N,dd,lgblk,myid,myidh,maxvsize >*/
/* -*- fortran -*- */

/* double precision functions */
/*<       integer nsupnode,hvbsize,comm >*/
/*<       integer lptrs(3,0:*), linds(*), iptrs(2,0:*) >*/
/*<       integer tptrs(3,0:*), tinds(*) >*/
/*<       integer sup(*), mysnodes(*) >*/
/*<       integer lrud(*), hvbndry(*) >*/
/*<       integer uinds(*), lc(*) >*/
/*<       double precision lvals(*) >*/
/*<       double precision rhsc(0:N-1),w(*) >*/
/*<       double precision uvec(*),uvl(*) >*/
/*<       double precision rhs(*),recvec(*) >*/
/*<       integer lendp,loglendp,itype >*/
/*<       parameter(lendp=8,loglendp=3,itype=1) >*/
/*<       double precision one,ngone,zero >*/
/*<       parameter(zero=0.d0,one=1.d0,ngone=-1.d0) >*/
/*<       integer root >*/
/*<       integer level,level2,rowcol,bmaskh,bmaskv,bmaskhk,bmaskvk >*/
/*<       integer myleft,myright,myup,mydown >*/
/*<       integer supptr,suptop,supbot,supsiz,lbotsiz >*/
/*<       integer bvs,nb,nhb,nvb,nbrp,hvbp,ih,bhp >*/
/*<       integer bhstrt,bhsize,bhend,lvalst,ldalb >*/
/*<       integer rhsst,uvecst,uvlst,itag >*/
/*<       integer iv,bvp,bvstrt,bvsize,bvend >*/
/*<       integer recvsize,cmpsize >*/
/*<       integer npending,npendings,npendings2 >*/
/*<       integer mid,nbrecv >*/
/*<       integer lppar,kidtop,pparbot,kid >*/
/*<       integer diffbits,partner,ik,ldauk,lpkid,pkid >*/
/*<       integer i,j,k1,k2,is,usize >*/
/*<       integer mymaskh,mymaskhk,ldaup >*/
/*<       integer iamkid,iampar,kidr,kidl,pkido,bip,bip1,istflag >*/
/*<       integer vsizer,nbode,levelk,iaml,bnode,lvp,ivlim,nvinds,ntinds >*/
/*<       integer flagr,flags,indh,hsize >*/
/*<       integer cbsize,msizedp >*/
/*<       integer mpistat(MPI_STATUS_SIZE),req(5),ierr >*/
/*<       integer statall(MPI_STATUS_SIZE,5) >*/
/*<       do i=1,5 >*/
    /* Parameter adjustments */
    --mysnodes;
    --sup;
    --lptrs;
    --linds;
    --lvals;
    --tptrs;
    --tinds;
    --rhs;
    --uvec;
    --uvl;
    --recvec;
    --uinds;
    --lrud;
    --hvbndry;
    --lc;
    --w;
    --iptrs;

    /* Function Body */
    for (i__ = 1; i__ <= 5; ++i__) {
/*<         req(i) = MPI_REQUEST_NULL >*/
	req[i__ - 1] = MPI_REQUEST_NULL;
/*<       end do >*/
    }
/*<       root   = mysnodes(1) >*/
    root = mysnodes[1];
/*<       suptop = mysnodes(nsupnode) >*/
    suptop = mysnodes[*nsupnode];
/*<    >*/
    fsolve_(n, &lvals[1], &linds[1], &lptrs[1], &tinds[1], &tptrs[1], &sup[1],
	     rhsc, &c__1, &suptop, &lc[1], &iptrs[1], &w[1]);
/*<       supptr = tptrs(3,suptop) >*/
    supptr = tptrs[suptop * 3 + 3];
/*<       supbot = sup(supptr) >*/
    supbot = sup[supptr];
/*<       bhsize = sup(supptr+1) >*/
    bhsize = sup[supptr + 1];
/*<       ldalb  = lptrs(2,supbot) >*/
    ldalb = lptrs[supbot * 3 + 2];
/*<       lpkid  = lptrs(3,supbot) >*/
    lpkid = lptrs[supbot * 3 + 3];
/*<       ldauk  = ldalb-bhsize   >*/
    ldauk = ldalb - bhsize;
/*<       lpkid = lpkid+bhsize    >*/
    lpkid += bhsize;
/*<       do i=1,ldauk >*/
    i__1 = ldauk;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         uvl(i) = w(i+bhsize) >*/
	uvl[i__] = w[i__ + bhsize];
/*<       end do >*/
    }
/*<       level = dd   >*/
    level = *dd;
/*<       rowcol = 0   >*/
    rowcol = 0;
/*<       bmaskh = 0   >*/
    bmaskh = 0;
/*<       bmaskv = 0   >*/
    bmaskv = 0;
/*<       mymaskh = 0  >*/
    mymaskh = 0;
/*<       hsize  = 1   >*/
    hsize = 1;
/*<       hvbp = 1     >*/
    hvbp = 1;
/*<       nbrp = -3    >*/
    nbrp = -3;
/*<       do is = nsupnode-1,1,-1 >*/
    for (is = *nsupnode - 1; is >= 1; --is) {
/*<         bnode = ishft(suptop,-lgblk) >*/
	bnode = lbit_shift(suptop, -(*lgblk));
/*<         pkid = iand(bnode,bmaskh) >*/
	pkid = bnode & bmaskh;
/*<         suptop = mysnodes(is)         >*/
	suptop = mysnodes[is];
/*<         supptr = tptrs(3,suptop)       >*/
	supptr = tptrs[suptop * 3 + 3];
/*<         supbot = sup(supptr)             >*/
	supbot = sup[supptr];
/*<         supsiz = sup(supptr+1)       >*/
	supsiz = sup[supptr + 1];
/*<         bmaskhk = bmaskh             >*/
	bmaskhk = bmaskh;
/*<         bmaskvk = bmaskv             >*/
	bmaskvk = bmaskv;
/*<         mymaskhk = mymaskh             >*/
	mymaskhk = mymaskh;
/*<         levelk  = level                   >*/
	levelk = level;
/*<         if (tptrs(2,supbot).ne.1) then  >*/
	if (tptrs[supbot * 3 + 2] != 1) {
/*<           level = level-1 >*/
	    --level;
/*<           rowcol = 1-rowcol >*/
	    rowcol = 1 - rowcol;
/*<           level2 = ishft(dd-level-1,-1) >*/
	    level2 = lbit_shift(*dd - level - 1, (ftnlen)-1);
/*<           if(rowcol.eq.1) then   >*/
	    if (rowcol == 1) {
/*<             bmaskh = ior(bmaskh,ishft(1,level2)) >*/
		bmaskh |= lbit_shift((ftnlen)1, level2);
/*<             mymaskh = iand(myidh,bmaskh) >*/
		mymaskh = *myidh & bmaskh;
/*<             hsize = ishft(hsize,1) >*/
		hsize <<= 1;
/*<           else                   >*/
	    } else {
/*<             bmaskv = ior(bmaskv,ishft(1,level2)) >*/
		bmaskv |= lbit_shift((ftnlen)1, level2);
/*<           end if >*/
	    }
/*<           nbrp = nbrp + 4  >*/
	    nbrp += 4;
/*<         end if >*/
	}
/*<         bnode = ishft(supbot,-lgblk) >*/
	bnode = lbit_shift(supbot, -(*lgblk));
/*<         ldaup = lptrs(2,supbot)        >*/
	ldaup = lptrs[supbot * 3 + 2];
/*<         lppar = lptrs(3,supbot)       >*/
	lppar = lptrs[supbot * 3 + 3];
/*<         pparbot = iand(bnode,bmaskh) >*/
	pparbot = bnode & bmaskh;
/*<         iamkid = 0  >*/
	iamkid = 0;
/*<         iampar = 0  >*/
	iampar = 0;
/*<         if(mymaskhk.eq.pkid) iamkid = 1 >*/
	if (mymaskhk == pkid) {
	    iamkid = 1;
	}
/*<         if(mymaskh.eq.pparbot) iampar = 1 >*/
	if (mymaskh == pparbot) {
	    iampar = 1;
	}
/*<         kidl = tinds(tptrs(1,supbot))  >*/
	kidl = tinds[tptrs[supbot * 3 + 1]];
/*<         kidr = tinds(tptrs(1,supbot)+1)  >*/
	kidr = tinds[tptrs[supbot * 3 + 1] + 1];
/*<         iaml = 1  >*/
	iaml = 1;
/*<         if(iand(myid,ishft(1,dd-level-1)).ne.0) iaml = 0 >*/
	if ((*myid & lbit_shift((ftnlen)1, *dd - level - 1)) != 0) {
	    iaml = 0;
	}
/*<         cbsize  = ldauk*lendp  >*/
	cbsize = ldauk << 3;
/*<         msizedp = maxvsize*lendp >*/
	msizedp = *maxvsize << 3;
/*<         if(bmaskv.eq.bmaskvk) then >*/
	if (bmaskv == bmaskvk) {
/*<           partner  = myid >*/
	    partner = *myid;
/*<           if(iamkid.eq.1 .and. iampar.eq.1) then >*/
	    if (iamkid == 1 && iampar == 1) {
/*<             if(levelk.ne.level) then   >*/
		if (levelk != level) {
/*<               kid = kidr >*/
		    kid = kidr;
/*<               if(iaml.eq.0) kid = kidl >*/
		    if (iaml == 0) {
			kid = kidl;
		    }
/*<               bnode = ishft(kid,-lgblk) >*/
		    bnode = lbit_shift(kid, -(*lgblk));
/*<               pkido = iand(bnode,bmaskhk) >*/
		    pkido = bnode & bmaskhk;
/*<               diffbits = ieor(ieor(pkid,pkido),ishft(1,level2)) >*/
		    diffbits = pkid ^ pkido ^ lbit_shift((ftnlen)1, level2);
/*<               do j = 0,level2 >*/
		    i__1 = level2;
		    for (j = 0; j <= i__1; ++j) {
/*<                 partner = ieor(partner,ishft(iand(diffbits,1),j*2)) >*/
			partner ^= lbit_shift(diffbits & 1, j << 1);
/*<                 diffbits = ishft(diffbits,-1) >*/
			diffbits = lbit_shift(diffbits, (ftnlen)-1);
/*<               end do >*/
		    }
/*<    >*/
		    MPI_Irecv(&recvec[1], msizedp, MPI_BYTE, partner, 1, *comm, &req[1]);
/*		    mpi_irecv__(&recvec[1], &msizedp, &c__5, &partner, &c__1, comm, &req[1], &ierr); */
/*<    >*/
		    MPI_Irecv(&uinds[1], *maxvsize, MPI_INT, partner, 1, *comm, req);
/*		    mpi_irecv__(&uinds[1], maxvsize, &c__11, &partner, &c__1, comm, req, &ierr); */
/*<               npending = 2 >*/
		    npending = 2;
/*<               do while(npending.gt.0) >*/
		    while(npending > 0) {
/*<                 call mpi_waitany(2,req,mid,mpistat,ierr) >*/
			MPI_Waitany(2, req, &mid, &mpistat);
/*<                 if(mid.eq.1) then  >*/
			if (mid == 1) {
/*<                   call mpi_get_count(mpistat,MPI_INTEGER,usize,ierr) >*/
			    MPI_Get_count(&mpistat, MPI_INT, &usize);
/*<                 end if >*/
			}
/*<                 npending = npending-1 >*/
			--npending;
/*<               end do >*/
		    }
/*<    >*/
		    x10dad1_(&uvec[1], &ldaup, &uvl[1], &ldauk, &recvec[1], &
			    usize, &linds[lppar], &linds[lpkid], &uinds[1]);
/*<             else >*/
		} else {
/*<    >*/
		    extend_op1__(&uvec[1], &ldaup, &uvl[1], &ldauk, &linds[
			    lppar], &linds[lpkid]);
/*<             end if >*/
		}
/*<           end if >*/
	    }
/*<           if(iamkid.eq.0 .and. iampar.eq.1) then >*/
	    if (iamkid == 0 && iampar == 1) {
/*<             diffbits = ieor(pkid,mymaskhk) >*/
		diffbits = pkid ^ mymaskhk;
/*<             do j = 0,level2 >*/
		i__1 = level2;
		for (j = 0; j <= i__1; ++j) {
/*<               partner = ieor(partner,ishft(iand(diffbits,1),j*2)) >*/
		    partner ^= lbit_shift(diffbits & 1, j << 1);
/*<               diffbits = ishft(diffbits,-1) >*/
		    diffbits = lbit_shift(diffbits, (ftnlen)-1);
/*<             end do >*/
		}
/*<    >*/
		MPI_Irecv(&uvl[1], msizedp, MPI_BYTE, partner, 1, *comm, &req[2]);
/*		mpi_irecv__(&uvl[1], &msizedp, &c__5, &partner, &c__1, comm, &req[2], &ierr); */
/*<             npending = 1 >*/
		npending = 1;
/*<             if(levelk.ne.level) then >*/
		if (levelk != level) {
/*<               kid = kidr >*/
		    kid = kidr;
/*<               if(iaml.eq.0) kid = kidl >*/
		    if (iaml == 0) {
			kid = kidl;
		    }
/*<               bnode = ishft(kid,-lgblk) >*/
		    bnode = lbit_shift(kid, -(*lgblk));
/*<               pkido = iand(bnode,bmaskhk) >*/
		    pkido = bnode & bmaskhk;
/*<               diffbits = ieor(pparbot,pkido) >*/
		    diffbits = pparbot ^ pkido;
/*<               if(iaml.eq.1) diffbits = ieor(diffbits,ishft(1,level2)) >*/
		    if (iaml == 1) {
			diffbits ^= lbit_shift((ftnlen)1, level2);
		    }
/*<               partner = myid >*/
		    partner = *myid;
/*<               do j = 0,level2 >*/
		    i__1 = level2;
		    for (j = 0; j <= i__1; ++j) {
/*<                 partner = ieor(partner,ishft(iand(diffbits,1),j*2)) >*/
			partner ^= lbit_shift(diffbits & 1, j << 1);
/*<                 diffbits = ishft(diffbits,-1) >*/
			diffbits = lbit_shift(diffbits, (ftnlen)-1);
/*<               end do >*/
		    }
/*<    >*/
		    MPI_Irecv(&recvec[1], msizedp, MPI_BYTE, partner, 1, *comm, &req[1]);
/*<    >*/
		    MPI_Irecv(&uinds[1], *maxvsize, MPI_INT, partner, 1, *comm, req);
/*<               npending = npending+2 >*/
		    npending += 2;
/*<               do while(npending.gt.0) >*/
		    while(npending > 0) {
/*<                 call mpi_waitany(3,req,mid,mpistat,ierr) >*/
			MPI_Waitany(3, req, &mid, &mpistat);
/*<                 if(mid.eq.1) then >*/
			if (mid == 1) {
/*<                   call mpi_get_count(mpistat,MPI_INTEGER,usize,ierr) >*/
			    MPI_Get_count(&mpistat, MPI_INT, &usize);
/*<                 end if >*/
			}
/*<                 npending = npending-1 >*/
			--npending;
/*<               end do >*/
		    }
/*<    >*/
		    x10dad1_(&uvec[1], &ldaup, &uvl[1], &ldauk, &recvec[1], &
			    usize, &linds[lppar], &hvbndry[bip], &uinds[1]);
/*<             else >*/
		} else {
/*<               call mpi_wait(req(3),mpistat,ierr) >*/
		    MPI_Wait(&req[2], &mpistat);
/*<    >*/
		    extend_op1__(&uvec[1], &ldaup, &uvl[1], &ldauk, &linds[
			    lppar], &hvbndry[bip]);
/*<             end if >*/
		}
/*<           end if >*/
	    }
/*<           if(iamkid.eq.1 .and. iampar.eq.0) then >*/
	    if (iamkid == 1 && iampar == 0) {
/*<             diffbits = ieor(pparbot,mymaskh) >*/
		diffbits = pparbot ^ mymaskh;
/*<             i = diffbits >*/
		i__ = diffbits;
/*<             do j = 0,level2 >*/
		i__1 = level2;
		for (j = 0; j <= i__1; ++j) {
/*<               partner = ieor(partner,ishft(iand(i,1),j*2)) >*/
		    partner ^= lbit_shift(i__ & 1, j << 1);
/*<               i = ishft(i,-1) >*/
		    i__ = lbit_shift(i__, (ftnlen)-1);
/*<             end do >*/
		}
/*<    >*/
		MPI_Isend(&uvl[1], cbsize, MPI_BYTE, partner, 1, *comm, &req[3]);
/*<             npending = 1 >*/
		npending = 1;
/*<             if(iand(diffbits,ishft(1,level2)).ne.0) then  >*/
		if ((diffbits & lbit_shift((ftnlen)1, level2)) != 0) {
/*<    >*/
		    MPI_Isend(&linds[lpkid], ldauk, MPI_INT, partner, 1, *comm, &req[4]);
/*<               npending = 2 >*/
		    npending = 2;
/*<             end if >*/
		}
/*<             do while(npending.gt.0) >*/
		while(npending > 0) {
/*<               call mpi_waitany(2,req(4),mid,mpistat,ierr) >*/
		    MPI_Waitany(2, &req[3], &mid, &mpistat);
/*<               npending = npending-1 >*/
		    --npending;
/*<             end do >*/
		}
/*<           end if >*/
	    }
/*<         else >*/
	} else {
/*<           if(iamkid.eq.1 .and. iampar.eq.1) then >*/
	    if (iamkid == 1 && iampar == 1) {
/*<             partner = ieor(myid,ishft(1,dd-level-1)) >*/
		partner = *myid ^ lbit_shift((ftnlen)1, *dd - level - 1);
		MPI_Isend(&linds[lpkid], ldauk, MPI_INT, partner, 1, *comm, &req[3]);
		MPI_Isend(&uvl[1], cbsize, MPI_BYTE, partner, 1, *comm, &req[4]);
		MPI_Irecv(&uinds[1], *maxvsize, MPI_INT, partner, 1, *comm, req);
		MPI_Irecv(&recvec[1], msizedp, MPI_BYTE, partner, 1, *comm, &req[1]);
/*<             call mpi_wait(req(1),mpistat,ierr) >*/
		MPI_Wait(req, &mpistat);
/*<             call mpi_get_count(mpistat,MPI_INTEGER,usize,ierr) >*/
		MPI_Get_count(&mpistat, MPI_INT, &usize);
/*<             call mpi_wait(req(2),mpistat,ierr) >*/
		MPI_Wait(&req[1], &mpistat);
/*<    >*/
		x10dad1_(&uvec[1], &ldaup, &uvl[1], &ldauk, &recvec[1], &
			usize, &linds[lppar], &linds[lpkid], &uinds[1]);
/*<             call mpi_wait(req(4),mpistat,ierr) >*/
		MPI_Wait(&req[3], &mpistat);
/*<             call mpi_wait(req(5),mpistat,ierr) >*/
		MPI_Wait(&req[4], &mpistat);
/*<           end if >*/
	    }
/*<           if(iamkid.eq.0 .and. iampar.eq.1) then >*/
	    if (iamkid == 0 && iampar == 1) {
/*<             diffbits = ieor(pkid,mymaskhk) >*/
		diffbits = pkid ^ mymaskhk;
/*<             partner = myid >*/
		partner = *myid;
/*<             do j=0,level2 >*/
		i__1 = level2;
		for (j = 0; j <= i__1; ++j) {
/*<               partner = ieor(partner,ishft(iand(diffbits,1),j*2)) >*/
		    partner ^= lbit_shift(diffbits & 1, j << 1);
/*<               diffbits = ishft(diffbits,-1) >*/
		    diffbits = lbit_shift(diffbits, (ftnlen)-1);
/*<             end do >*/
		}
/*<    >*/
		MPI_Irecv(&uvl[1], msizedp, MPI_BYTE, partner, 1, *comm, &req[2]);
/*<             partner = ieor(myid,ishft(1,dd-level-1)) >*/
		partner = *myid ^ lbit_shift((ftnlen)1, *dd - level - 1);
/*<    >*/
		MPI_Irecv(&uinds[1], *maxvsize, MPI_INT, partner, 1, *comm, req);
/*<    >*/
		MPI_Irecv(&recvec[1], msizedp, MPI_BYTE, partner, 1, *comm, &req[1]);
/*<    >*/
		MPI_Isend(&hvbndry[bip], ldauk, MPI_INT, partner, 1, *comm, &req[4]);
/*<             call mpi_wait(req(3),mpistat,ierr) >*/
		MPI_Wait(&req[2], &mpistat);
/*<    >*/
		MPI_Isend(&uvl[1], cbsize, MPI_BYTE, partner, 1, *comm, &req[3]);
/*<             call mpi_wait(req(1),mpistat,ierr) >*/
		MPI_Wait(req, &mpistat);
/*<             call mpi_get_count(mpistat,MPI_INTEGER,usize,ierr) >*/
		MPI_Get_count(&mpistat, MPI_INT, &usize);
/*<             call mpi_wait(req(2),mpistat,ierr) >*/
		MPI_Wait(&req[1], &mpistat);
/*<    >*/
		x10dad1_(&uvec[1], &ldaup, &uvl[1], &ldauk, &recvec[1], &
			usize, &linds[lppar], &hvbndry[bip], &uinds[1]);
/*<             call mpi_wait(req(4),mpistat,ierr) >*/
		MPI_Wait(&req[3], &mpistat);
/*<             call mpi_wait(req(5),mpistat,ierr) >*/
		MPI_Wait(&req[4], &mpistat);
/*<           end if >*/
	    }
/*<           if(iamkid.eq.1 .and. iampar.eq.0) then >*/
	    if (iamkid == 1 && iampar == 0) {
/*<             diffbits = ieor(pparbot,mymaskh) >*/
		diffbits = pparbot ^ mymaskh;
/*<             partner = myid >*/
		partner = *myid;
/*<             do j=0,level2 >*/
		i__1 = level2;
		for (j = 0; j <= i__1; ++j) {
/*<               partner = ieor(partner,ishft(iand(diffbits,1),j*2)) >*/
		    partner ^= lbit_shift(diffbits & 1, j << 1);
/*<               diffbits = ishft(diffbits,-1) >*/
		    diffbits = lbit_shift(diffbits, (ftnlen)-1);
/*<             end do >*/
		}
/*<    >*/
		MPI_Isend(&uvl[1], cbsize, MPI_BYTE, partner, 1, *comm, &req[2]);
/*<             call mpi_wait(req(3),mpistat,ierr) >*/
		MPI_Wait(&req[2], &mpistat);
/*<           end if >*/
	    }
/*<         end if >*/
	}
/*<         bvs     = hvbndry(hvbp+1) >*/
	bvs = hvbndry[hvbp + 1];
/*<         nb      = hvbndry(hvbp+2) >*/
	nb = hvbndry[hvbp + 2];
/*<         istflag = hvbndry(hvbp+3) >*/
	istflag = hvbndry[hvbp + 3];
/*<         lbotsiz = hvbndry(hvbp+4) >*/
	lbotsiz = hvbndry[hvbp + 4];
/*<         nhb     = hvbndry(hvbp+5) >*/
	nhb = hvbndry[hvbp + 5];
/*<         vsizer  = hvbndry(bvs-3) >*/
	vsizer = hvbndry[bvs - 3];
/*<         bip     = hvbp+6+3*nhb+1 >*/
	bip = hvbp + 6 + nhb * 3 + 1;
/*<         bip1    = bip >*/
	bip1 = bip;
/*<         nvinds  = hvbndry(bip-1) >*/
	nvinds = hvbndry[bip - 1];
/*<         ntinds  = nvinds-vsizer >*/
	ntinds = nvinds - vsizer;
/*<         bip     = bip+ntinds >*/
	bip += ntinds;
/*<         myleft  = lrud(nbrp)             >*/
	myleft = lrud[nbrp];
/*<         myright = lrud(nbrp+1)       >*/
	myright = lrud[nbrp + 1];
/*<         myup    = lrud(nbrp+2)       >*/
	myup = lrud[nbrp + 2];
/*<         mydown  = lrud(nbrp+3)       >*/
	mydown = lrud[nbrp + 3];
/*<         uvlst   = 1                   >*/
	uvlst = 1;
/*<         ldauk   = vsizer             >*/
	ldauk = vsizer;
/*<         npending = 0  >*/
	npending = 0;
/*<         if(iampar.eq.0) then  >*/
	if (iampar == 0) {
/*<           do i=1,nvinds >*/
	    i__1 = nvinds;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/*<               uvec(i) = zero >*/
		uvec[i__] = 0.;
/*<             end do >*/
	    }
/*<         else >*/
	} else {
/*<           do i=1,nvinds >*/
	    i__1 = nvinds;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/*<             recvec(i) = zero >*/
		recvec[i__] = 0.;
/*<           end do >*/
	    }
/*<         end if >*/
	}
/*<         bhp = hvbp+6  >*/
	bhp = hvbp + 6;
/*<         if(nhb.ne.0) then >*/
	if (nhb != 0) {
/*<         if(mod(dd-level-1,2).eq.0) then >*/
	    if ((*dd - level - 1) % 2 == 0) {
/*<           i = bhp >*/
		i__ = bhp;
/*<           k1 = bip1 >*/
		k1 = bip1;
/*<           k2 = 1 >*/
		k2 = 1;
/*<           do ih=1,nhb >*/
		i__1 = nhb;
		for (ih = 1; ih <= i__1; ++ih) {
/*<             bhstrt = hvbndry(i) >*/
		    bhstrt = hvbndry[i__];
/*<             bhsize = hvbndry(i+1) >*/
		    bhsize = hvbndry[i__ + 1];
/*<             do while(hvbndry(k1).ne.bhstrt .and. k1-bip1.le.ntinds) >*/
		    while(hvbndry[k1] != bhstrt && k1 - bip1 <= ntinds) {
/*<               k1 = k1+1 >*/
			++k1;
/*<             end do >*/
		    }
/*<             if(bhstrt.eq.hvbndry(k1)) then >*/
		    if (bhstrt == hvbndry[k1]) {
/*<             do j=0,bhsize-1 >*/
			i__2 = bhsize - 1;
			for (j = 0; j <= i__2; ++j) {
/*<               uinds(k2+j) = hvbndry(k1+j) >*/
			    uinds[k2 + j] = hvbndry[k1 + j];
/*<             end do >*/
			}
/*<             k2 = k2+bhsize >*/
			k2 += bhsize;
/*<             end if >*/
		    }
/*<             i = i+3 >*/
		    i__ += 3;
/*<           end do >*/
		}
/*<           ntinds = k2-1 >*/
		ntinds = k2 - 1;
/*<           call bgetrhs1(N,uinds,ntinds,rhs,rhsc) >*/
		bgetrhs1_(n, &uinds[1], &ntinds, &rhs[1], rhsc);
/*<         else >*/
	    } else {
/*<           call bgetrhs1(N,hvbndry(bip1),ntinds,rhs,rhsc) >*/
		bgetrhs1_(n, &hvbndry[bip1], &ntinds, &rhs[1], rhsc);
/*<         end if >*/
	    }
/*<         rhsst = 1 >*/
	    rhsst = 1;
/*<         end if >*/
	}
/*<         do ih = 1,nhb >*/
	i__1 = nhb;
	for (ih = 1; ih <= i__1; ++ih) {
/*<           bhstrt = hvbndry(bhp)         >*/
	    bhstrt = hvbndry[bhp];
/*<           bhsize = hvbndry(bhp+1)       >*/
	    bhsize = hvbndry[bhp + 1];
/*<           bhend  = hvbndry(bhp+2)       >*/
	    bhend = hvbndry[bhp + 2];
/*<           lvalst = lptrs(1,bhstrt)   >*/
	    lvalst = lptrs[bhstrt * 3 + 1];
/*<           ldalb  = lptrs(2,bhstrt)   >*/
	    ldalb = lptrs[bhstrt * 3 + 2];
/*<           itag = MPI_ANY_TAG             >*/
	    itag = MPI_ANY_TAG;
/*<           npendings2 = 0             >*/
	    npendings2 = 0;
/*<           uvecst = 1                   >*/
	    uvecst = 1;
/*<           iv = 1 >*/
	    iv = 1;
/*<           bvp = hvbndry(hvbp+1) >*/
	    bvp = hvbndry[hvbp + 1];
/*<           nvb = hvbndry(bvp-1) >*/
	    nvb = hvbndry[bvp - 1];
/*<           do while (iv.le.nvb)  >*/
	    while(iv <= nvb) {
/*<             if(hvbndry(bvp).lt.bhstrt .or. hvbndry(bvp+1).eq.0) then >*/
		if (hvbndry[bvp] < bhstrt || hvbndry[bvp + 1] == 0) {
/*<               bvsize = hvbndry(bvp+1) >*/
		    bvsize = hvbndry[bvp + 1];
/*<               uvecst = uvecst+bvsize  >*/
		    uvecst += bvsize;
/*<               iv = iv+1 >*/
		    ++iv;
/*<               bvp = bvp+3 >*/
		    bvp += 3;
/*<             else >*/
		} else {
/*<               goto 10 >*/
		    goto L10;
/*<             end if >*/
		}
/*<           end do >*/
	    }
/*<  10          if(iv.le.nvb .and. bhsize.ne.0) then >*/
L10:
	    if (iv <= nvb && bhsize != 0) {
/*<           do while (iv.le.nvb) >*/
		while(iv <= nvb) {
/*<             bvstrt = hvbndry(bvp)       >*/
		    bvstrt = hvbndry[bvp];
/*<             bvsize = hvbndry(bvp+1)       >*/
		    bvsize = hvbndry[bvp + 1];
/*<             bvend  = hvbndry(bvp+2)       >*/
		    bvend = hvbndry[bvp + 2];
/*<             flags = 0 >*/
		    flags = 0;
/*<             if(ih+1 .gt. nhb) then >*/
		    if (ih + 1 > nhb) {
/*<               flags = 1 >*/
			flags = 1;
/*<             else >*/
		    } else {
/*<               if(hvbndry(bhp+3).gt.bvstrt) flags = 1 >*/
			if (hvbndry[bhp + 3] > bvstrt) {
			    flags = 1;
			}
/*<             end if >*/
		    }
/*<             flagr = 0 >*/
		    flagr = 0;
/*<             if(flags .eq. 1) then >*/
		    if (flags == 1) {
/*<               indh = bvend >*/
			indh = bvend;
/*<               if(indh.gt.suptop) indh = suptop >*/
			if (indh > suptop) {
			    indh = suptop;
			}
/*<               indh = iand(ishft(indh,-lgblk),bmaskh) >*/
			indh = lbit_shift(indh, -(*lgblk)) & bmaskh;
/*<    >*/
			if ((mymaskh + hsize - 1) % hsize != indh && bhstrt !=
				 supbot) {
			    flagr = 1;
			}
/*<               if(bhend.eq.suptop) flags = 0 >*/
			if (bhend == suptop) {
			    flags = 0;
			}
/*<             end if  >*/
		    }
/*<             if(bhstrt.eq.bvstrt) then  >*/
		    if (bhstrt == bvstrt) {
/*<               if(bhstrt.ne.supbot) then  >*/
			if (bhstrt != supbot) {
/*<    >*/
			    MPI_Irecv(&recvec[1], msizedp, MPI_BYTE, myleft, 1, *comm, req);
/*<                 call mpi_wait(req(1),mpistat,ierr) >*/
			    MPI_Wait(req, &mpistat);
/*<                 call mpi_get_count(mpistat,MPI_BYTE,nbrecv,ierr) >*/
			    MPI_Get_count(&mpistat, MPI_BYTE, &nbrecv);
/*<                 do i=0,bhsize-1 >*/
			    i__2 = bhsize - 1;
			    for (i__ = 0; i__ <= i__2; ++i__) {
/*<                   uvec(uvecst+i) = uvec(uvecst+i) + recvec(i+1) >*/
				uvec[uvecst + i__] += recvec[i__ + 1];
/*<                   rhs(rhsst+i) = rhs(rhsst+i) + uvec(uvecst+i) >*/
				rhs[rhsst + i__] += uvec[uvecst + i__];
/*<                 end do >*/
			    }
/*<               else >*/
			} else {
/*<                 do i=0,bhsize-1 >*/
			    i__2 = bhsize - 1;
			    for (i__ = 0; i__ <= i__2; ++i__) {
/*<                   rhs(rhsst+i) = rhs(rhsst+i) + uvec(uvecst+i) >*/
				rhs[rhsst + i__] += uvec[uvecst + i__];
/*<                 end do >*/
			    }
/*<               end if >*/
			}
/*<    >*/
			dtrsv_("l", "n", "n", &bhsize, &lvals[lvalst], &ldalb,
				 &rhs[rhsst], &c__1, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);
/*<               itag = myid >*/
			itag = *myid;
/*<    >*/
			if (mydown != *myid && ! (bvend == suptop && supsiz ==
				 lbotsiz)) {
/*<    >*/
			    i__2 = bhsize << 3;
			    MPI_Isend(&rhs[rhsst], i__2, MPI_BYTE, mydown, itag, *comm, &req[4]);
/*<                 npendings2 = 1 >*/
			    npendings2 = 1;
/*<               end if >*/
			}
/*<    >*/
			putrhs1_(n, &linds[lptrs[bhstrt * 3 + 3]], &bhsize, &
				rhs[rhsst], rhsc);
/*<               lvalst = lvalst+bhsize >*/
			lvalst += bhsize;
/*<               uvecst = uvecst+bhsize >*/
			uvecst += bhsize;
/*<             else  >*/
		    } else {
/*<               npendings = 0 >*/
			npendings = 0;
/*<               if(itag.eq.MPI_ANY_TAG) then >*/
			if (itag == MPI_ANY_TAG) {
/*<    >*/
			    MPI_Irecv(&rhs[rhsst], msizedp, MPI_BYTE, myup, itag, *comm, &req[3]);
/*<                 npendings = 1 >*/
			    npendings = 1;
/*<               end if >*/
			}
/*<               if(flagr.eq.1) then >*/
			if (flagr == 1) {
/*<                 if(npending.eq.1) call mpi_wait(req(1),mpistat,ierr) >*/
			    if (npending == 1) {
				MPI_Wait(req, &mpistat);
			    }
/*<    >*/
			    MPI_Irecv(&recvec[1], msizedp, MPI_BYTE, myleft, 1, *comm, req);
/*<               end if >*/
			}
/*<               if(bvstrt.gt.suptop) then  >*/
			if (bvstrt > suptop) {
/*<                 recvsize = 0  >*/
			    recvsize = 0;
/*<                 do while (iv.le.nvb) >*/
			    while(iv <= nvb) {
/*<                   recvsize = recvsize+hvbndry(bvp+1) >*/
				recvsize += hvbndry[bvp + 1];
/*<                   bvp = bvp+3 >*/
				bvp += 3;
/*<                   iv = iv+1 >*/
				++iv;
/*<                 end do >*/
			    }
/*<               else >*/
			} else {
/*<                 recvsize = bvsize >*/
			    recvsize = bvsize;
/*<               end if >*/
			}
/*<               if(npendings.eq.1) then >*/
			if (npendings == 1) {
/*<                 call mpi_wait(req(4),mpistat,ierr) >*/
			    MPI_Wait(&req[3], &mpistat);
/*<                 itag = mpistat(MPI_TAG) >*/
			    itag = mpistat.MPI_TAG;
/*<    >*/
			    if (itag != mydown && ! (bvend == suptop && 
				    lbotsiz == supsiz)) {
/*<    >*/
				i__2 = bhsize << 3;
				MPI_Isend(&rhs[rhsst], i__2, MPI_BYTE, mydown, itag, *comm, &req[4]);
/*<                   npendings2 = 1 >*/
				npendings2 = 1;
/*<                 end if >*/
			    }
/*<                 npendings = 0 >*/
			    npendings = 0;
/*<               end if >*/
			}
/*<    >*/
			dgemv_("n", &recvsize, &bhsize, &c_b58, &lvals[lvalst]
				, &ldalb, &rhs[rhsst], &c__1, &c_b60, &uvec[
				uvecst], &c__1, (ftnlen)1);
/*<               if(flagr.eq.1) then >*/
			if (flagr == 1) {
/*<                 call mpi_wait(req(1),mpistat,ierr) >*/
			    MPI_Wait(req, &mpistat);
/*<                 if(flags.eq.1) then >*/
			    if (flags == 1) {
/*<                   do i=0,recvsize-1 >*/
				i__2 = recvsize - 1;
				for (i__ = 0; i__ <= i__2; ++i__) {
/*<                     uvec(uvecst+i) = uvec(uvecst+i)+recvec(i+1) >*/
				    uvec[uvecst + i__] += recvec[i__ + 1];
/*<                   end do >*/
				}
/*<    >*/
				i__2 = recvsize << 3;
				MPI_Isend(&uvec[uvecst], i__2, MPI_BYTE, myright, 1, *comm, req);
/*<                   npending = 1 >*/
				npending = 1;
/*<                 else >*/
			    } else {
/*<                   do i = 0,recvsize-1 >*/
				i__2 = recvsize - 1;
				for (i__ = 0; i__ <= i__2; ++i__) {
/*<                     uvl(uvlst+i) = uvec(uvecst+i)+recvec(i+1) >*/
				    uvl[uvlst + i__] = uvec[uvecst + i__] + 
					    recvec[i__ + 1];
/*<                   end do >*/
				}
/*<                   uvlst = uvlst+recvsize >*/
				uvlst += recvsize;
/*<                 end if >*/
			    }
/*<               else >*/
			} else {
/*<                 if(flags.eq.1) then >*/
			    if (flags == 1) {
/*<                   if(npending.eq.1) call mpi_wait(req(1),mpistat,ierr) >*/
				if (npending == 1) {

				    MPI_Wait(req, &mpistat);
				}
/*<    >*/
				i__2 = recvsize << 3;
				MPI_Isend(&uvec[uvecst], i__2, MPI_BYTE, myright, 1, *comm, req);
/*<                   npending = 1 >*/
				npending = 1;
/*<                 else >*/
			    } else {
/*<                   if(bhend.eq.suptop) then >*/
				if (bhend == suptop) {
/*<                     do i = 0,recvsize-1 >*/
				    i__2 = recvsize - 1;
				    for (i__ = 0; i__ <= i__2; ++i__) {
/*<                       uvl(uvlst+i) = uvec(uvecst+i) >*/
					uvl[uvlst + i__] = uvec[uvecst + i__];
/*<                     end do >*/
				    }
/*<                     uvlst = uvlst+recvsize >*/
				    uvlst += recvsize;
/*<                   end if >*/
				}
/*<                 end if >*/
			    }
/*<               end if >*/
			}
/*<               uvecst = uvecst+recvsize >*/
			uvecst += recvsize;
/*<               lvalst = lvalst+recvsize >*/
			lvalst += recvsize;
/*<             end if >*/
		    }
/*<             bvp = bvp+3 >*/
		    bvp += 3;
/*<             iv = iv+1 >*/
		    ++iv;
/*<             do while (iv.le.nvb) >*/
		    while(iv <= nvb) {
/*<               if(hvbndry(bvp+1).eq.0) then >*/
			if (hvbndry[bvp + 1] == 0) {
/*<                 bvp = bvp+3 >*/
			    bvp += 3;
/*<                 iv = iv+1 >*/
			    ++iv;
/*<               else >*/
			} else {
/*<                 goto 30 >*/
			    goto L30;
/*<               end if >*/
			}
/*<             end do >*/
		    }
/*<  30       end do  >*/
L30:
		    ;
		}
/*<           else  >*/
	    } else {
/*<             if(iv.gt.nvb .and. bhsize.ne.0) then >*/
		if (iv > nvb && bhsize != 0) {
/*<               if(itag.eq.MPI_ANY_TAG .and. lbotsiz.ne.supsiz) then >*/
		    if (itag == MPI_ANY_TAG && lbotsiz != supsiz) {
/*<    >*/
			MPI_Irecv(&rhs[rhsst], msizedp, MPI_BYTE, myup, itag, *comm, &req[3]);
/*<                 call mpi_wait(req(4),mpistat,ierr) >*/
			MPI_Wait(&req[3], &mpistat);
/*<                 itag = mpistat(MPI_TAG) >*/
			itag = mpistat.MPI_TAG;
/*<                 if(itag.ne.mydown) then >*/
			if (itag != mydown) {
/*<    >*/
			    i__2 = bhsize << 3;
			    MPI_Isend(&rhs[rhsst], i__2, MPI_BYTE, mydown, itag, *comm, &req[4]);
/*<                   npendings2 = 1 >*/
			    npendings2 = 1;
/*<                 end if >*/
			}
/*<               end if >*/
		    }
/*<             else >*/
		} else {
/*<               do while (iv.le.nvb) >*/
		    while(iv <= nvb) {
/*<                 bvstrt = hvbndry(bvp)       >*/
			bvstrt = hvbndry[bvp];
/*<                 bvsize = hvbndry(bvp+1)       >*/
			bvsize = hvbndry[bvp + 1];
/*<                 bvend  = hvbndry(bvp+2)       >*/
			bvend = hvbndry[bvp + 2];
/*<                 flags = 0 >*/
			flags = 0;
/*<                 if(ih+1 .gt. nhb) then >*/
			if (ih + 1 > nhb) {
/*<                   flags = 1 >*/
			    flags = 1;
/*<                 else >*/
			} else {
/*<                   if(hvbndry(bhp+3).gt.bvstrt) flags = 1 >*/
			    if (hvbndry[bhp + 3] > bvstrt) {
				flags = 1;
			    }
/*<                 end if >*/
			}
/*<                 flagr = 0 >*/
			flagr = 0;
/*<                 if(flags .eq. 1) then >*/
			if (flags == 1) {
/*<                   indh = bvend >*/
			    indh = bvend;
/*<                   if(indh.gt.suptop) indh = suptop >*/
			    if (indh > suptop) {
				indh = suptop;
			    }
/*<                   indh = iand(ishft(indh,-lgblk),bmaskh) >*/
			    indh = lbit_shift(indh, -(*lgblk)) & bmaskh;
/*<                   if(mod(mymaskh+hsize-1,hsize).ne.indh) flagr=1 >*/
			    if ((mymaskh + hsize - 1) % hsize != indh) {
				flagr = 1;
			    }
/*<                 end if  >*/
			}
/*<                 if(flagr.eq.1) then >*/
			if (flagr == 1) {
/*<                   if(npending.eq.1) call mpi_wait(req(1),mpistat,ierr) >*/
			    if (npending == 1) {
				MPI_Wait(req, &mpistat);
			    }
/*<    >*/
			    MPI_Irecv(&recvec[1], msizedp, MPI_BYTE, myleft, 1, *comm, req);
/*<                   call mpi_wait(req(1),mpistat,ierr) >*/
			    MPI_Wait(req, &mpistat);
/*<                   call mpi_get_count(mpistat,MPI_BYTE,nbrecv,ierr) >*/
			    MPI_Get_count(&mpistat, MPI_BYTE, &nbrecv);
/*<                   recvsize = ishft(nbrecv,-loglendp) >*/
			    recvsize = lbit_shift(nbrecv, (ftnlen)-3);
/*<                   do i=0,recvsize-1 >*/
			    i__2 = recvsize - 1;
			    for (i__ = 0; i__ <= i__2; ++i__) {
/*<                     uvec(uvecst+i) = uvec(uvecst+i)+recvec(i+1) >*/
				uvec[uvecst + i__] += recvec[i__ + 1];
/*<                   end do >*/
			    }
/*<                   npending = 0 >*/
			    npending = 0;
/*<                   if(bvstrt.gt.suptop) then  >*/
			    if (bvstrt > suptop) {
/*<                     cmpsize = 0 >*/
				cmpsize = 0;
/*<                     do while(cmpsize.lt.recvsize) >*/
				while(cmpsize < recvsize) {
/*<                       cmpsize = cmpsize+hvbndry(bvp+1) >*/
				    cmpsize += hvbndry[bvp + 1];
/*<                       bvp = bvp+3 >*/
				    bvp += 3;
/*<                       iv = iv+1 >*/
				    ++iv;
/*<                     end do >*/
				}
/*<                   end if >*/
			    }
/*<                 else >*/
			} else {
/*<                   if(bvstrt.gt.suptop) then  >*/
			    if (bvstrt > suptop) {
/*<                     recvsize = 0  >*/
				recvsize = 0;
/*<                     do while (iv.le.nvb) >*/
				while(iv <= nvb) {
/*<                       recvsize = recvsize+hvbndry(bvp+1) >*/
				    recvsize += hvbndry[bvp + 1];
/*<                       bvp = bvp+3 >*/
				    bvp += 3;
/*<                       iv = iv+1 >*/
				    ++iv;
/*<                     end do >*/
				}
/*<                   else >*/
			    } else {
/*<                     recvsize = bvsize >*/
				recvsize = bvsize;
/*<                   end if >*/
			    }
/*<                 end if >*/
			}
/*<                 if(flags.eq.1) then >*/
			if (flags == 1) {
/*<    >*/
			    i__2 = recvsize << 3;
			    MPI_Isend(&uvec[uvecst], i__2, MPI_BYTE, myright, 1, *comm, req);
/*<                   npending = 1 >*/
			    npending = 1;
/*<                 end if >*/
			}
/*<                 uvecst = uvecst+recvsize >*/
			uvecst += recvsize;
/*<                 bvp = bvp+3 >*/
			bvp += 3;
/*<                 iv = iv+1 >*/
			++iv;
/*<                 do while (iv.le.nvb) >*/
			while(iv <= nvb) {
/*<                   if(hvbndry(bvp+1).eq.0) then >*/
			    if (hvbndry[bvp + 1] == 0) {
/*<                     bvp = bvp+3 >*/
				bvp += 3;
/*<                     iv = iv+1 >*/
				++iv;
/*<                   else >*/
			    } else {
/*<                     goto 40 >*/
				goto L40;
/*<                   end if >*/
			    }
/*<                 end do >*/
			}
/*<  40           end do >*/
L40:
			;
		    }
/*<             end if >*/
		}
/*<           end if  >*/
	    }
/*<           if(npendings2.eq.1) then >*/
	    if (npendings2 == 1) {
/*<             call mpi_wait(req(5),mpistat,ierr) >*/
		MPI_Wait(&req[4], &mpistat);
/*<             npendings2 = 0 >*/
		npendings2 = 0;
/*<           end if >*/
	    }
/*<           if(npending.eq.1) then >*/
	    if (npending == 1) {
/*<             call mpi_wait(req(1),mpistat,ierr) >*/
		MPI_Wait(req, &mpistat);
/*<             npending = 0 >*/
		npending = 0;
/*<           end if >*/
	    }
/*<           rhsst  = rhsst+bhsize >*/
	    rhsst += bhsize;
/*<           bhp = bhp+3 >*/
	    bhp += 3;
/*<         end do >*/
	}
/*<         call mpi_waitall(5,req,statall,ierr) >*/
	MPI_Waitall(5, req, statall);
/*<         lpkid = lptrs(3,suptop) >*/
	lpkid = lptrs[suptop * 3 + 3];
/*<         if(istflag.eq.1) lpkid = lpkid+1 >*/
	if (istflag == 1) {
	    ++lpkid;
	}
/*<         hvbp = hvbndry(hvbp) >*/
	hvbp = hvbndry[hvbp];
/*<       end do  >*/
    }
/*<       end >*/
    return 0;
} /* pfsolve1_ */

/*<    >*/
/* Subroutine */ int x10dad1_(doublereal *pvec, integer *psize, doublereal *
	kvec, integer *ksize, doublereal *rvec, integer *rsize, integer *
	indsp, integer *indsk, integer *indsr)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer ik, ip, ir;

/*<       integer psize,ksize,rsize,indsp(*),indsk(*),indsr(*) >*/
/*<       double precision pvec(*),kvec(*),rvec(*),zero >*/
/*<       parameter(zero=0.d0) >*/
/*<         ik = 1 >*/
    /* Parameter adjustments */
    --indsr;
    --indsk;
    --indsp;
    --rvec;
    --kvec;
    --pvec;

    /* Function Body */
    ik = 1;
/*<         ir = 1 >*/
    ir = 1;
/*<         do ip = 1,psize >*/
    i__1 = *psize;
    for (ip = 1; ip <= i__1; ++ip) {
/*<           do while(indsr(ir).lt.indsp(ip) .and. ir.le.rsize) >*/
	while(indsr[ir] < indsp[ip] && ir <= *rsize) {
/*<             ir = ir+1 >*/
	    ++ir;
/*<           end do >*/
	}
/*<           do while(indsk(ik).lt.indsp(ip) .and. ik.le.ksize) >*/
	while(indsk[ik] < indsp[ip] && ik <= *ksize) {
/*<             ik = ik+1 >*/
	    ++ik;
/*<           end do >*/
	}
/*<           if(ir.le.rsize .and. indsr(ir).eq.indsp(ip)) then >*/
	if (ir <= *rsize && indsr[ir] == indsp[ip]) {
/*<             pvec(ip) = rvec(ir) >*/
	    pvec[ip] = rvec[ir];
/*<           else >*/
	} else {
/*<             pvec(ip) = zero >*/
	    pvec[ip] = 0.;
/*<           end if >*/
	}
/*<           if(ik.le.ksize .and. indsk(ik).eq.indsp(ip)) then >*/
	if (ik <= *ksize && indsk[ik] == indsp[ip]) {
/*<             pvec(ip) = pvec(ip)+kvec(ik) >*/
	    pvec[ip] += kvec[ik];
/*<           end if >*/
	}
/*<         end do >*/
    }
/*<       end  >*/
    return 0;
} /* x10dad1_ */

/*<       subroutine extend_op1(pvec,psize,kvec,ksize,indsp,indsk) >*/
/* Subroutine */ int extend_op1__(doublereal *pvec, integer *psize, 
	doublereal *kvec, integer *ksize, integer *indsp, integer *indsk)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer ik, ip;

/*<       integer psize,ksize,indsp(*),indsk(*) >*/
/*<       double precision pvec(*),kvec(*) >*/
/*<         ik = 1 >*/
    /* Parameter adjustments */
    --indsk;
    --indsp;
    --kvec;
    --pvec;

    /* Function Body */
    ik = 1;
/*<         do ip = 1,psize >*/
    i__1 = *psize;
    for (ip = 1; ip <= i__1; ++ip) {
/*<           do while(indsk(ik).lt.indsp(ip) .and. ik.le.ksize) >*/
	while(indsk[ik] < indsp[ip] && ik <= *ksize) {
/*<             ik = ik+1 >*/
	    ++ik;
/*<           end do >*/
	}
/*<           if(ik.le.ksize .and. indsk(ik).eq.indsp(ip)) then >*/
	if (ik <= *ksize && indsk[ik] == indsp[ip]) {
/*<             pvec(ip) = kvec(ik) >*/
	    pvec[ip] = kvec[ik];
/*<           else >*/
	} else {
/*<             pvec(ip) = 0.d0 >*/
	    pvec[ip] = 0.;
/*<           end if >*/
	}
/*<         end do >*/
    }
/*<       end >*/
    return 0;
} /* extend_op1__ */

