/* $Id: pfsolvem.c,v 1.4 2005-01-05 16:51:31 paklein Exp $ */
/* pfsolvem.f -- translated by f2c (version 20030320).
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
static integer c__11 = 11;
static integer c__2 = 2;
static integer c__3 = 3;
static doublereal c_b50 = 1.;
static doublereal c_b58 = -1.;

/* /+***************************************************************************+/ */
/* /+                                                                           +/ */
/* /+   (C) Copyright IBM Corporation, 1997                                     +/ */
/* /+   (C) Copyright Regents of the University of Minnesota, 1997              +/ */
/* /+                                                                           +/ */
/* /+   pfsolvem.f                                                              +/ */
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
/* /+ $Id: pfsolvem.c,v 1.4 2005-01-05 16:51:31 paklein Exp $ +/ */
/* /+***************************************************************************+/ */

static integer lbit_shift(integer a, integer b) {
	return b >= 0 ? a << b : (integer)((uinteger)a >> -b);
};

/*<    >*/
/* Subroutine */ int pfsolvem_(integer *mysnodes, integer *nsupnode, integer *
	sup, integer *lptrs, integer *linds, doublereal *lvals, integer *
	tptrs, integer *tinds, integer *myid, integer *myidh, integer *dd, 
	integer *lgblk, integer *n, doublereal *rhsc, integer *nrhs, 
	doublereal *rhs, doublereal *uvec, doublereal *uvl, doublereal *
	recvec, integer *uinds, integer *maxvsize, integer *lrud, integer *
	hvbndry, integer *hvbsize, integer *lc, doublereal *w, integer *iptrs,
	 MPI_Comm *comm)
{
    /* System generated locals */
    integer rhsc_dim1, rhsc_offset, i__1, i__2, i__3;

    /* Local variables */
    integer diffbits, npending;
    integer mymaskhk, recvsize, i__, j;
    integer npendings;
    extern /* Subroutine */ int extend_op__(doublereal *, integer *, 
	    doublereal *, integer *, integer *, integer *, integer *);
    integer nb, ih, kd, npendings2, is, ks, iv, ks1, nhb, kid, bhp, mid, bip, 
	    nvb, bvp, /* req[5],*/ bvs;
	MPI_Request req[5];
    integer itag, kidl, iaml, pkid, indh, kidr, hvbp, nbrp, ierr, root, myup, 
	    ldalb, bhend;
    extern /* Subroutine */ int x10dad_(doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *);
    integer bnode;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    integer flagr, flags, bvend, ldauk, lpkid, level, ldaup, pkido, lppar, 
	    brhsr, hsize;
    extern /* Subroutine */ int dtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    integer usize, uvlst, level2;
    integer iamkid, bmaskh, iampar, nbrecv, cbsize, levelk, bmaskv, bhsize;
    extern /* Subroutine */ int fsolve_(integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, integer *, integer *, integer *, doublereal *);
    integer myleft, nvinds, bvsize, vcsize, rowcol, bhstrt, lvalst, uvecst, 
	    supbot, mydown, vsizer, bvstrt;
    extern /* Subroutine */ int putrhs_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *);
    integer suptop, supptr, supsiz, bmaskhk, ldauvec, istflag, bmaskvk;
    extern /* Subroutine */ int bgetrhs_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *);
    integer /* statall[20], */	/* was [4][5] */ pparbot, cmpsize, partner, mymaskh, 
	    msizedp, /* mpistat[4],*/ myright, lbotsiz;
	MPI_Status statall[5];
	MPI_Status mpistat;

/*<       implicit none >*/
/*<       include 'mpif.h' >*/
/*<       integer N,dd,lgblk,myid,myidh,maxvsize,nrhs >*/
/* -*- fortran -*- */

/* double precision functions */
/*<       integer nsupnode,hvbsize,comm >*/
/*<       integer lptrs(3,0:*), linds(*), iptrs(2,0:*) >*/
/*<       integer tptrs(3,0:*), tinds(*) >*/
/*<       integer sup(*), mysnodes(*) >*/
/*<       integer lrud(*), hvbndry(*) >*/
/*<       integer uinds(*), lc(*) >*/
/*<       double precision lvals(*) >*/
/*<       double precision rhsc(0:N-1,nrhs),w(*) >*/
/*<       double precision uvec(*),uvl(*) >*/
/*<       double precision rhs(*),recvec(*) >*/
/*<       integer lendp,loglendp,itype >*/
/*<       parameter(lendp=8,loglendp=3,itype=1) >*/
/*<       double precision one,ngone,zero >*/
/*<       parameter(one=1.d0,ngone=-1.d0,zero=0.d0) >*/
/*<       integer root >*/
/*<       integer level,level2,rowcol,bmaskh,bmaskv,bmaskhk,bmaskvk >*/
/*<       integer myleft,myright,myup,mydown >*/
/*<       integer supptr,suptop,supbot,supsiz,lbotsiz >*/
/*<       integer bvs,nb,nhb,nvb,nbrp,hvbp,ih,bhp >*/
/*<       integer bhstrt,bhsize,bhend,lvalst,ldalb >*/
/*<       integer uvecst,uvlst,itag >*/
/*<       integer iv,bvp,bvstrt,bvsize,bvend >*/
/*<       integer recvsize,cmpsize >*/
/*<       integer npending,npendings,npendings2 >*/
/*<       integer mid,nbrecv >*/
/*<       integer lppar,kidtop,pparbot,kid >*/
/*<       integer diffbits,partner,ik,ldauk,lpkid,pkid >*/
/*<       integer i,j,ks,kd,is,usize >*/
/*<       integer mymaskh,mymaskhk,ldaup >*/
/*<       integer iamkid,iampar,kidr,kidl,pkido,bip,istflag >*/
/*<       integer vsizer,nbode,levelk,iaml,bnode,lvp,ivlim,nvinds >*/
/*<       integer brhsr,msizedp,vcsize,cbsize >*/
/*<       integer flagr,flags,indh,hsize,ldauvec,ks1 >*/
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
    rhsc_dim1 = *n - 1 - 0 + 1;
    rhsc_offset = 0 + rhsc_dim1;
    rhsc -= rhsc_offset;
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
	     &rhsc[rhsc_offset], nrhs, &suptop, &lc[1], &iptrs[1], &w[1]);
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
/*<       ldauk  = ldalb-bhsize >*/
    ldauk = ldalb - bhsize;
/*<       lpkid = lpkid+bhsize >*/
    lpkid += bhsize;
/*<       do j=0,nrhs-1 >*/
    i__1 = *nrhs - 1;
    for (j = 0; j <= i__1; ++j) {
/*<         kd = j*ldauk >*/
	kd = j * ldauk;
/*<         ks = j*ldalb+bhsize >*/
	ks = j * ldalb + bhsize;
/*<         do i=1,ldauk >*/
	i__2 = ldauk;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<           uvl(kd+i) = w(ks+i) >*/
	    uvl[kd + i__] = w[ks + i__];
/*<         end do >*/
	}
/*<       end do >*/
    }
/*<       level = dd >*/
    level = *dd;
/*<       rowcol = 0 >*/
    rowcol = 0;
/*<       bmaskh = 0 >*/
    bmaskh = 0;
/*<       bmaskv = 0 >*/
    bmaskv = 0;
/*<       mymaskh = 0 >*/
    mymaskh = 0;
/*<       hsize  = 1 >*/
    hsize = 1;
/*<       brhsr = nrhs*lendp >*/
    brhsr = *nrhs << 3;
/*<       hvbp = 1 >*/
    hvbp = 1;
/*<       nbrp = -3  >*/
    nbrp = -3;
/*<       do is = nsupnode-1,1,-1 >*/
    for (is = *nsupnode - 1; is >= 1; --is) {
/*<         bnode = ishft(suptop,-lgblk) >*/
	bnode = lbit_shift(suptop, -(*lgblk));
/*<         pkid = iand(bnode,bmaskh) >*/
	pkid = bnode & bmaskh;
/*<         suptop = mysnodes(is) >*/
	suptop = mysnodes[is];
/*<         supptr = tptrs(3,suptop) >*/
	supptr = tptrs[suptop * 3 + 3];
/*<         supbot = sup(supptr) >*/
	supbot = sup[supptr];
/*<         supsiz = sup(supptr+1) >*/
	supsiz = sup[supptr + 1];
/*<         bmaskhk = bmaskh >*/
	bmaskhk = bmaskh;
/*<         bmaskvk = bmaskv >*/
	bmaskvk = bmaskv;
/*<         mymaskhk = mymaskh >*/
	mymaskhk = mymaskh;
/*<         levelk  = level >*/
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
/*<             hsize  = ishft(hsize,1) >*/
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
/*<         ldaup = lptrs(2,supbot) >*/
	ldaup = lptrs[supbot * 3 + 2];
/*<         lppar = lptrs(3,supbot) >*/
	lppar = lptrs[supbot * 3 + 3];
/*<         pparbot = iand(bnode,bmaskh) >*/
	pparbot = bnode & bmaskh;
/*<         iamkid = 0 >*/
	iamkid = 0;
/*<         iampar = 0 >*/
	iampar = 0;
/*<         if(mymaskhk.eq.pkid) iamkid = 1 >*/
	if (mymaskhk == pkid) {
	    iamkid = 1;
	}
/*<         if(mymaskh.eq.pparbot) iampar = 1 >*/
	if (mymaskh == pparbot) {
	    iampar = 1;
	}
/*<         kidl = tinds(tptrs(1,supbot)) >*/
	kidl = tinds[tptrs[supbot * 3 + 1]];
/*<         kidr = tinds(tptrs(1,supbot)+1) >*/
	kidr = tinds[tptrs[supbot * 3 + 1] + 1];
/*<         iaml = 1  >*/
	iaml = 1;
/*<         if(iand(myid,ishft(1,dd-level-1)).ne.0) iaml = 0 >*/
	if ((*myid & lbit_shift((ftnlen)1, *dd - level - 1)) != 0) {
	    iaml = 0;
	}
/*<         cbsize  = ldauk*brhsr >*/
	cbsize = ldauk * brhsr;
/*<         msizedp = maxvsize*brhsr >*/
	msizedp = *maxvsize * brhsr;
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
/*			mpi_waitany__(&c__2, req, &mid, mpistat, &ierr); */

/*<                 if(mid.eq.1) then >*/
			if (mid == 1) {
/*<                   call mpi_get_count(mpistat,MPI_INTEGER,usize,ierr) >*/

			    myMPI_Get_count(&mpistat, MPI_INT, &usize);
/*			    mpi_get_count__(mpistat, &c__11, &usize, &ierr); */

/*<                 end if >*/
			}
/*<                 npending = npending-1 >*/
			--npending;
/*<               end do >*/
		    }
/*<    >*/
		    x10dad_(&uvec[1], &ldaup, &uvl[1], &ldauk, &recvec[1], &
			    usize, nrhs, &linds[lppar], &linds[lpkid], &uinds[
			    1]);
/*<             else >*/
		} else {
/*<    >*/
		    extend_op__(&uvec[1], &ldaup, &uvl[1], &ldauk, nrhs, &
			    linds[lppar], &linds[lpkid]);
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
/*		    mpi_irecv__(&recvec[1], &msizedp, &c__5, &partner, &c__1, comm, &req[1], &ierr); */
/*<    >*/
		    MPI_Irecv(&uinds[1], *maxvsize, MPI_INT, partner, 1, *comm, req);
/*		    mpi_irecv__(&uinds[1], maxvsize, &c__11, &partner, &c__1, comm, req, &ierr); */
/*<               npending = npending+2 >*/
		    npending += 2;
/*<               do while(npending.gt.0) >*/
		    while(npending > 0) {
/*<                 call mpi_waitany(3,req,mid,mpistat,ierr) >*/

			MPI_Waitany(3, req, &mid, &mpistat);
/*			mpi_waitany__(&c__3, req, &mid, mpistat, &ierr); */

/*<                 if(mid.eq.1) then >*/
			if (mid == 1) {
/*<                   call mpi_get_count(mpistat,MPI_INTEGER,usize,ierr) >*/

			    myMPI_Get_count(&mpistat, MPI_INT, &usize);
/*			    mpi_get_count__(mpistat, &c__11, &usize, &ierr); */

/*<                 end if >*/
			}
/*<                 npending = npending-1 >*/
			--npending;
/*<               end do >*/
		    }
/*<    >*/
		    x10dad_(&uvec[1], &ldaup, &uvl[1], &ldauk, &recvec[1], &
			    usize, nrhs, &linds[lppar], &hvbndry[bip], &uinds[
			    1]);
/*<             else >*/
		} else {
/*<               call mpi_wait(req(3),mpistat,ierr) >*/

		    MPI_Wait(&req[2], &mpistat);
/*		    mpi_wait__(&req[2], mpistat, &ierr); */

/*<    >*/
		    extend_op__(&uvec[1], &ldaup, &uvl[1], &ldauk, nrhs, &
			    linds[lppar], &hvbndry[bip]);
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

		myMPI_Isend(&uvl[1], cbsize, MPI_BYTE, partner, 1, *comm, &req[3]);
/*		mpi_isend__(&uvl[1], &cbsize, &c__5, &partner, &c__1, comm, &req[3], &ierr); */

/*<             npending = 1 >*/
		npending = 1;
/*<             if(iand(diffbits,ishft(1,level2)).ne.0) then  >*/
		if ((diffbits & lbit_shift((ftnlen)1, level2)) != 0) {
/*<    >*/

		    myMPI_Isend(&linds[lpkid], ldauk, MPI_INT, partner, 1, *comm, &req[4]);
/*		    mpi_isend__(&linds[lpkid], &ldauk, &c__11, &partner, &c__1, comm, &req[4], &ierr); */

/*<               npending = 2 >*/
		    npending = 2;
/*<             end if >*/
		}
/*<             do while(npending.gt.0) >*/
		while(npending > 0) {
/*<               call mpi_waitany(2,req(4),mid,mpistat,ierr) >*/

		    MPI_Waitany(2, &req[3], &mid, &mpistat);
/*		    mpi_waitany__(&c__2, &req[3], &mid, mpistat, &ierr); */

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
/*<    >*/

		myMPI_Isend(&linds[lpkid], ldauk, MPI_INT, partner, 1, *comm, &req[3]);
/*		mpi_isend__(&linds[lpkid], &ldauk, &c__11, &partner, &c__1, comm, &req[3], &ierr); */

/*<    >*/
		myMPI_Isend(&uvl[1], cbsize, MPI_BYTE, partner, 1, *comm, &req[4]);
/*		mpi_isend__(&uvl[1], &cbsize, &c__5, &partner, &c__1, comm, &req[4], &ierr); */
/*<    >*/
		MPI_Irecv(&uinds[1], *maxvsize, MPI_INT, partner, 1, *comm, req);
/*		mpi_irecv__(&uinds[1], maxvsize, &c__11, &partner, &c__1, comm, req, &ierr); */
/*<    >*/
		MPI_Irecv(&recvec[1], msizedp, MPI_BYTE, partner, 1, *comm, &req[1]);
/*		mpi_irecv__(&recvec[1], &msizedp, &c__5, &partner, &c__1, comm, &req[1], &ierr); */
/*<             call mpi_wait(req(1),mpistat,ierr) >*/
		MPI_Wait(req, &mpistat);
/*		mpi_wait__(req, mpistat, &ierr); */
/*<             call mpi_get_count(mpistat,MPI_INTEGER,usize,ierr) >*/
		myMPI_Get_count(&mpistat, MPI_INT, &usize);
/*		mpi_get_count__(mpistat, &c__11, &usize, &ierr); */
/*<             call mpi_wait(req(2),mpistat,ierr) >*/
		MPI_Wait(&req[1], &mpistat);
/*		mpi_wait__(&req[1], mpistat, &ierr); */
/*<    >*/
		x10dad_(&uvec[1], &ldaup, &uvl[1], &ldauk, &recvec[1], &usize,
			 nrhs, &linds[lppar], &linds[lpkid], &uinds[1]);
/*<             call mpi_wait(req(4),mpistat,ierr) >*/
		MPI_Wait(&req[3], &mpistat);
/*		mpi_wait__(&req[3], mpistat, &ierr); */
/*<             call mpi_wait(req(5),mpistat,ierr) >*/
		MPI_Wait(&req[4], &mpistat);
/*		mpi_wait__(&req[4], mpistat, &ierr); */
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
/*		mpi_irecv__(&uvl[1], &msizedp, &c__5, &partner, &c__1, comm, &req[2], &ierr); */
/*<             partner = ieor(myid,ishft(1,dd-level-1)) >*/
		partner = *myid ^ lbit_shift((ftnlen)1, *dd - level - 1);
/*<    >*/
		MPI_Irecv(&uinds[1], *maxvsize, MPI_INT, partner, 1, *comm, req);
/*		mpi_irecv__(&uinds[1], maxvsize, &c__11, &partner, &c__1, comm, req, &ierr); */
/*<    >*/
		MPI_Irecv(&recvec[1], msizedp, MPI_BYTE, partner, 1, *comm, &req[1]);
/*		mpi_irecv__(&recvec[1], &msizedp, &c__5, &partner, &c__1, comm, &req[1], &ierr); */
/*<    >*/
		myMPI_Isend(&hvbndry[bip], ldauk, MPI_INT, partner, 1, *comm, &req[4]);
/*		mpi_isend__(&hvbndry[bip], &ldauk, &c__11, &partner, &c__1, comm, &req[4], &ierr); */
/*<             call mpi_wait(req(3),mpistat,ierr) >*/
		MPI_Wait(&req[2], &mpistat);
/*		mpi_wait__(&req[2], mpistat, &ierr); */
/*<    >*/
		myMPI_Isend(&uvl[1], cbsize, MPI_BYTE, partner, 1, *comm, &req[3]);
/*		mpi_isend__(&uvl[1], &cbsize, &c__5, &partner, &c__1, comm, &req[3], &ierr); */
/*<             call mpi_wait(req(1),mpistat,ierr) >*/
		MPI_Wait(req, &mpistat);
/*		mpi_wait__(req, mpistat, &ierr); */
/*<             call mpi_get_count(mpistat,MPI_INTEGER,usize,ierr) >*/
		myMPI_Get_count(&mpistat, MPI_INT, &usize);
/*		mpi_get_count__(mpistat, &c__11, &usize, &ierr); */
/*<             call mpi_wait(req(2),mpistat,ierr) >*/
		MPI_Wait(&req[1], &mpistat);
/*		mpi_wait__(&req[1], mpistat, &ierr); */
/*<    >*/
		x10dad_(&uvec[1], &ldaup, &uvl[1], &ldauk, &recvec[1], &usize,
			 nrhs, &linds[lppar], &hvbndry[bip], &uinds[1]);
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
		myMPI_Isend(&uvl[1], cbsize, MPI_BYTE, partner, 1, *comm, &req[2]);
/*		mpi_isend__(&uvl[1], &cbsize, &c__5, &partner, &c__1, comm, &req[2], &ierr); */
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
/*<         bip     = hvbp+6+3*nhb+1 >*/
	bip = hvbp + 6 + nhb * 3 + 1;
/*<         nvinds  = hvbndry(bip-1) >*/
	nvinds = hvbndry[bip - 1];
/*<         vsizer  = hvbndry(bvs-3) >*/
	vsizer = hvbndry[bvs - 3];
/*<         bip     = bip+nvinds-vsizer >*/
	bip = bip + nvinds - vsizer;
/*<         myleft  = lrud(nbrp) >*/
	myleft = lrud[nbrp];
/*<         myright = lrud(nbrp+1) >*/
	myright = lrud[nbrp + 1];
/*<         myup    = lrud(nbrp+2) >*/
	myup = lrud[nbrp + 2];
/*<         mydown  = lrud(nbrp+3) >*/
	mydown = lrud[nbrp + 3];
/*<         uvlst   = 1 >*/
	uvlst = 1;
/*<         ldauk   = vsizer >*/
	ldauk = vsizer;
/*<         ldauvec = nvinds >*/
	ldauvec = nvinds;
/*<         npending = 0  >*/
	npending = 0;
/*<         if(iampar.eq.0) then >*/
	if (iampar == 0) {
/*<           do i=1,nvinds*nrhs >*/
	    i__1 = nvinds * *nrhs;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/*<             uvec(i) = zero >*/
		uvec[i__] = 0.;
/*<           end do >*/
	    }
/*<         else >*/
	} else {
/*<           do i=1,nvinds*nrhs >*/
	    i__1 = nvinds * *nrhs;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/*<             recvec(i) = zero >*/
		recvec[i__] = 0.;
/*<           end do >*/
	    }
/*<         end if >*/
	}
/*<         bhp = hvbp+6  >*/
	bhp = hvbp + 6;
/*<         do ih = 1,nhb >*/
	i__1 = nhb;
	for (ih = 1; ih <= i__1; ++ih) {
/*<           bhstrt = hvbndry(bhp) >*/
	    bhstrt = hvbndry[bhp];
/*<           bhsize = hvbndry(bhp+1) >*/
	    bhsize = hvbndry[bhp + 1];
/*<           bhend  = hvbndry(bhp+2) >*/
	    bhend = hvbndry[bhp + 2];
/*<           lvalst = lptrs(1,bhstrt) >*/
	    lvalst = lptrs[bhstrt * 3 + 1];
/*<           ldalb  = lptrs(2,bhstrt) >*/
	    ldalb = lptrs[bhstrt * 3 + 2];
/*<           itag = MPI_ANY_TAG >*/
	    itag = MPI_ANY_TAG;
/*<           npendings2 = 0 >*/
	    npendings2 = 0;
/*<           vcsize = bhsize*brhsr >*/
	    vcsize = bhsize * brhsr;
/*<           uvecst = 1 >*/
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
/*<               uvecst = uvecst+hvbndry(bvp+1) >*/
		    uvecst += hvbndry[bvp + 1];
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
/*<             bvstrt = hvbndry(bvp) >*/
		    bvstrt = hvbndry[bvp];
/*<             bvsize = hvbndry(bvp+1) >*/
		    bvsize = hvbndry[bvp + 1];
/*<             bvend  = hvbndry(bvp+2) >*/
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
/*<             end if >*/
		    }
/*<             if(bhstrt.eq.bvstrt) then  >*/
		    if (bhstrt == bvstrt) {
/*<    >*/
			bgetrhs_(n, &linds[lptrs[bhstrt * 3 + 3]], &bhsize, &
				bhsize, nrhs, &rhs[1], &rhsc[rhsc_offset]);
/*<               if(bhstrt.ne.supbot) then  >*/
			if (bhstrt != supbot) {
/*<    >*/
			    MPI_Irecv(&recvec[1], msizedp, MPI_BYTE, myleft, 1, *comm, req);
/*			    mpi_irecv__(&recvec[1], &msizedp, &c__5, &myleft, &c__1, comm, req, &ierr); */
/*<                 call mpi_wait(req(1),mpistat,ierr) >*/
			    MPI_Wait(req, &mpistat);
/*<                 do j=0,nrhs-1 >*/
			    i__2 = *nrhs - 1;
			    for (j = 0; j <= i__2; ++j) {
/*<                   ks = j*bhsize+1 >*/
				ks = j * bhsize + 1;
/*<                   ks1 = uvecst+j*ldauvec >*/
				ks1 = uvecst + j * ldauvec;
/*<                   do i=0,bhsize-1 >*/
				i__3 = bhsize - 1;
				for (i__ = 0; i__ <= i__3; ++i__) {
/*<                     uvec(ks1+i) = uvec(ks1+i) + recvec(ks+i) >*/
				    uvec[ks1 + i__] += recvec[ks + i__];
/*<                     rhs(ks+i) = rhs(ks+i) + uvec(ks1+i) >*/
				    rhs[ks + i__] += uvec[ks1 + i__];
/*<                   end do >*/
				}
/*<                 end do >*/
			    }
/*<               else >*/
			} else {
/*<                 do j=0,nrhs-1 >*/
			    i__2 = *nrhs - 1;
			    for (j = 0; j <= i__2; ++j) {
/*<                   kd = j*bhsize+1 >*/
				kd = j * bhsize + 1;
/*<                   ks = uvecst+j*ldauvec >*/
				ks = uvecst + j * ldauvec;
/*<                   do i=0,bhsize-1 >*/
				i__3 = bhsize - 1;
				for (i__ = 0; i__ <= i__3; ++i__) {
/*<                     rhs(kd+i) = rhs(kd+i) + uvec(ks+i) >*/
				    rhs[kd + i__] += uvec[ks + i__];
/*<                   end do >*/
				}
/*<                 end do >*/
			    }
/*<               end if >*/
			}
/*<    >*/
			dtrsm_("l", "l", "n", "n", &bhsize, nrhs, &c_b50, &
				lvals[lvalst], &ldalb, &rhs[1], &bhsize, (
				ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
/*<               itag = myid >*/
			itag = *myid;
/*<    >*/
			if (mydown != *myid && ! (bvend == suptop && supsiz ==
				 lbotsiz)) {
/*<    >*/
			    myMPI_Isend(&rhs[1], vcsize, MPI_BYTE, mydown, itag, *comm, &req[4]);
/*			    mpi_isend__(&rhs[1], &vcsize, &c__5, &mydown, &itag, comm, &req[4], &ierr); */
/*<                 npendings2 = 1 >*/
			    npendings2 = 1;
/*<               end if >*/
			}
/*<    >*/
			putrhs_(n, &linds[lptrs[bhstrt * 3 + 3]], &bhsize, &
				bhsize, nrhs, &rhs[1], &rhsc[rhsc_offset]);
/*<               uvecst = uvecst+bhsize >*/
			uvecst += bhsize;
/*<               lvalst = lvalst+bhsize >*/
			lvalst += bhsize;
/*<             else  >*/
		    } else {
/*<               npendings = 0 >*/
			npendings = 0;
/*<               if(itag.eq.MPI_ANY_TAG) then >*/
			if (itag == MPI_ANY_TAG) {
/*<    >*/
			    MPI_Irecv(&rhs[1], msizedp, MPI_BYTE, myup, itag, *comm, &req[3]);
/*			    mpi_irecv__(&rhs[1], &msizedp, &c__5, &myup, &itag, comm, &req[3], &ierr); */
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
/*			    mpi_irecv__(&recvec[1], &msizedp, &c__5, &myleft, &c__1, comm, req, &ierr); */
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
				myMPI_Isend(&rhs[1], vcsize, MPI_BYTE, mydown, itag, *comm, &req[4]);
/*				mpi_isend__(&rhs[1], &vcsize, &c__5, &mydown, &itag, comm, &req[4], &ierr); */
/*<                   npendings2 = 1 >*/
				npendings2 = 1;
/*<                 end if >*/
			    }
/*<                 npendings = 0 >*/
			    npendings = 0;
/*<               end if >*/
			}
/*<    >*/
			dgemm_("n", "n", &recvsize, nrhs, &bhsize, &c_b58, &
				lvals[lvalst], &ldalb, &rhs[1], &bhsize, &
				c_b50, &uvec[uvecst], &ldauvec, (ftnlen)1, (
				ftnlen)1);
/*<               if(flagr.eq.1) then >*/
			if (flagr == 1) {
/*<                 call mpi_wait(req(1),mpistat,ierr) >*/
			    MPI_Wait(req, &mpistat);
/*<                 if(flags.eq.1) then >*/
			    if (flags == 1) {
/*<                   do j=0,nrhs-1 >*/
				i__2 = *nrhs - 1;
				for (j = 0; j <= i__2; ++j) {
/*<                     kd = uvecst+j*ldauvec >*/
				    kd = uvecst + j * ldauvec;
/*<                     ks = j*recvsize+1 >*/
				    ks = j * recvsize + 1;
/*<                     do i=0,recvsize-1 >*/
				    i__3 = recvsize - 1;
				    for (i__ = 0; i__ <= i__3; ++i__) {
/*<                       recvec(ks+i) = recvec(ks+i) + uvec(kd+i) >*/
					recvec[ks + i__] += uvec[kd + i__];
/*<                     end do >*/
				    }
/*<                   end do >*/
				}
/*<    >*/
				i__2 = recvsize * brhsr;
				myMPI_Isend(&recvec[1], i__2, MPI_BYTE, myright, 1, *comm, req);
/*				mpi_isend__(&recvec[1], &i__2, &c__5, &myright, &c__1, comm, req, &ierr); */
/*<                   npending = 1 >*/
				npending = 1;
/*<                 else >*/
			    } else {
/*<                   do j=0,nrhs-1 >*/
				i__2 = *nrhs - 1;
				for (j = 0; j <= i__2; ++j) {
/*<                     kd = j*ldauk+uvlst >*/
				    kd = j * ldauk + uvlst;
/*<                     ks = j*ldauvec+uvecst >*/
				    ks = j * ldauvec + uvecst;
/*<                     ks1 = j*recvsize+1 >*/
				    ks1 = j * recvsize + 1;
/*<                     do i = 0,recvsize-1 >*/
				    i__3 = recvsize - 1;
				    for (i__ = 0; i__ <= i__3; ++i__) {
/*<                       uvl(kd+i) = uvec(ks+i) + recvec(ks1+i) >*/
					uvl[kd + i__] = uvec[ks + i__] + 
						recvec[ks1 + i__];
/*<                     end do >*/
				    }
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
/*<                   do j=0,nrhs-1 >*/
				i__2 = *nrhs - 1;
				for (j = 0; j <= i__2; ++j) {
/*<                     kd = 1+recvsize*j >*/
				    kd = recvsize * j + 1;
/*<                     ks = uvecst+j*ldauvec >*/
				    ks = uvecst + j * ldauvec;
/*<                     do i=0,recvsize-1 >*/
				    i__3 = recvsize - 1;
				    for (i__ = 0; i__ <= i__3; ++i__) {
/*<                       recvec(kd+i) = uvec(ks+i) >*/
					recvec[kd + i__] = uvec[ks + i__];
/*<                     end do >*/
				    }
/*<                   end do >*/
				}
/*<    >*/
				i__2 = recvsize * brhsr;
				myMPI_Isend(&recvec[1], i__2, MPI_BYTE, myright, 1, *comm, req);
/*				mpi_isend__(&recvec[1], &i__2, &c__5, &myright, &c__1, comm, req, &ierr); */
/*<                   npending = 1 >*/
				npending = 1;
/*<                 else >*/
			    } else {
/*<                   if(bhend.eq.suptop) then >*/
				if (bhend == suptop) {
/*<                     do j=0,nrhs-1 >*/
				    i__2 = *nrhs - 1;
				    for (j = 0; j <= i__2; ++j) {
/*<                       kd = j*ldauk+uvlst >*/
					kd = j * ldauk + uvlst;
/*<                       ks = j*ldauvec+uvecst >*/
					ks = j * ldauvec + uvecst;
/*<                       do i = 0,recvsize-1 >*/
					i__3 = recvsize - 1;
					for (i__ = 0; i__ <= i__3; ++i__) {
/*<                         uvl(kd+i) = uvec(ks+i) >*/
					    uvl[kd + i__] = uvec[ks + i__];
/*<                       end do >*/
					}
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
			MPI_Irecv(&rhs[1], msizedp, MPI_BYTE, myup, itag, *comm, &req[3]);
/*			mpi_irecv__(&rhs[1], &msizedp, &c__5, &myup, &itag, comm, &req[3], &ierr); */
/*<                 call mpi_wait(req(4),mpistat,ierr) >*/
			MPI_Wait(&req[3], &mpistat);
/*<                 itag = mpistat(MPI_TAG) >*/
			itag = mpistat.MPI_TAG;

/*<                 call mpi_get_count(mpistat,MPI_BYTE,mpistat,ierr) >*/
/*<                 itag = mpistat(MPI_TAG) >*/
/*
 * NOTE: What's going on here? mpistat is being passed to receive the length of the
 * message?.
 *
 */
/*
			myMPI_Get_count(mpistat, MPI_BYTE, mpistat);
			itag = mpistat[MPI_TAG];
*/

/*<                 call mpi_get_count(mpistat,MPI_BYTE,nbrecv,ierr) >*/
			myMPI_Get_count(&mpistat, MPI_BYTE, &nbrecv);
/*<                 if(itag.ne.mydown) then >*/
			if (itag != mydown) {
/*<    >*/
			    myMPI_Isend(&rhs[1], nbrecv, MPI_BYTE, mydown, itag, *comm, &req[4]);
/*			    mpi_isend__(&rhs[1], &nbrecv, &c__5, &mydown, &itag, comm, &req[4], &ierr); */

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
/*<                 bvstrt = hvbndry(bvp)      >*/
			bvstrt = hvbndry[bvp];
/*<                 bvsize = hvbndry(bvp+1)    >*/
			bvsize = hvbndry[bvp + 1];
/*<                 bvend  = hvbndry(bvp+2)    >*/
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
/*<                 end if >*/
			}
/*<                 if(flagr.eq.1) then >*/
			if (flagr == 1) {
/*<                   if(npending.eq.1) call mpi_wait(req(1),mpistat,ierr) >*/
			    if (npending == 1) {
				MPI_Wait(req, &mpistat);
			    }
/*<    >*/
			    MPI_Irecv(&recvec[1], msizedp, MPI_BYTE, myleft, 1, *comm, req);
/*			    mpi_irecv__(&recvec[1], &msizedp, &c__5, &myleft, &c__1, comm, req, &ierr); */

/*<                   call mpi_wait(req(1),mpistat,ierr) >*/
			    MPI_Wait(req, &mpistat);
/*<                   call mpi_get_count(mpistat,MPI_BYTE,nbrecv,ierr) >*/
			    myMPI_Get_count(&mpistat, MPI_BYTE, &nbrecv);
/*<                   recvsize = ishft(nbrecv,-loglendp)/nrhs >*/
			    recvsize = lbit_shift(nbrecv, (ftnlen)-3) / *nrhs;
/*<                   do j=0,nrhs-1 >*/
			    i__2 = *nrhs - 1;
			    for (j = 0; j <= i__2; ++j) {
/*<                     kd = uvecst+j*ldauvec >*/
				kd = uvecst + j * ldauvec;
/*<                     ks = j*recvsize+1 >*/
				ks = j * recvsize + 1;
/*<                     do i=0,recvsize-1 >*/
				i__3 = recvsize - 1;
				for (i__ = 0; i__ <= i__3; ++i__) {
/*<                       uvec(kd+i) = uvec(kd+i) + recvec(ks+i) >*/
				    uvec[kd + i__] += recvec[ks + i__];
/*<                     end do >*/
				}
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
/*<                   if(npending.eq.1) call mpi_wait(req(1),mpistat,ierr) >*/
			    if (npending == 1) {
/*				mpi_wait__(req, mpistat, &ierr); */
				MPI_Wait(req, &mpistat);
			    }
/*<                   do j=0,nrhs-1 >*/
			    i__2 = *nrhs - 1;
			    for (j = 0; j <= i__2; ++j) {
/*<                     kd = 1+recvsize*j >*/
				kd = recvsize * j + 1;
/*<                     ks = uvecst+j*ldauvec >*/
				ks = uvecst + j * ldauvec;
/*<                     do i=0,recvsize-1 >*/
				i__3 = recvsize - 1;
				for (i__ = 0; i__ <= i__3; ++i__) {
/*<                       recvec(kd+i) = uvec(ks+i) >*/
				    recvec[kd + i__] = uvec[ks + i__];
/*<                     end do >*/
				}
/*<                   end do >*/
			    }
/*<    >*/
			    i__2 = recvsize * brhsr;
			    myMPI_Isend(&recvec[1], i__2, MPI_BYTE, myright, 1, *comm, req);
/*			    mpi_isend__(&recvec[1], &i__2, &c__5, &myright, &c__1, comm, req, &ierr); */

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
/*<             call mpi_wait(req(1),mpistat,ierr)  >*/
		MPI_Wait(req, &mpistat);
/*<             npending = 0 >*/
		npending = 0;
/*<           end if >*/
	    }
/*<           bhp = bhp+3 >*/
	    bhp += 3;
/*<         end do >*/
	}
/*<         call mpi_waitall(5,req,statall,ierr) >*/
	MPI_Waitall(5, req, statall);
/*	mpi_waitall__(&c__5, req, statall, &ierr); */
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
} /* pfsolvem_ */

/*<       subroutine extend_op(pvec,psize,kvec,ksize,nrhs,indsp,indsk) >*/
/* Subroutine */ int extend_op__(doublereal *pvec, integer *psize, doublereal 
	*kvec, integer *ksize, integer *nrhs, integer *indsp, integer *indsk)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    integer k, ik, ip;

/*<       integer psize,ksize,indsp(*),indsk(*),nrhs >*/
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
/*<             do k=0,nrhs-1 >*/
	    i__2 = *nrhs - 1;
	    for (k = 0; k <= i__2; ++k) {
/*<               pvec(psize*k+ip) = kvec(ksize*k+ik) >*/
		pvec[*psize * k + ip] = kvec[*ksize * k + ik];
/*<             end do >*/
	    }
/*<           else >*/
	} else {
/*<             do k=0,nrhs-1 >*/
	    i__2 = *nrhs - 1;
	    for (k = 0; k <= i__2; ++k) {
/*<               pvec(psize*k+ip) = 0.d0 >*/
		pvec[*psize * k + ip] = 0.;
/*<             end do >*/
	    }
/*<           end if >*/
	}
/*<         end do >*/
    }
/*<       end >*/
    return 0;
} /* extend_op__ */

/* @process opt(3) strict debug(inline) */
/*<    >*/
/* Subroutine */ int x10dad_(doublereal *pvec, integer *psize, doublereal *
	kvec, integer *ksize, doublereal *rvec, integer *rsize, integer *nrhs,
	 integer *indsp, integer *indsk, integer *indsr)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    integer k, ik, ip, ir, idst;

/*<       integer psize,ksize,rsize,nrhs,indsp(*),indsk(*),indsr(*) >*/
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
/*<             do k=0,nrhs-1 >*/
	    i__2 = *nrhs - 1;
	    for (k = 0; k <= i__2; ++k) {
/*<               pvec(psize*k+ip) = rvec(rsize*k+ir) >*/
		pvec[*psize * k + ip] = rvec[*rsize * k + ir];
/*<             end do >*/
	    }
/*<           else >*/
	} else {
/*<             do k=0,nrhs-1 >*/
	    i__2 = *nrhs - 1;
	    for (k = 0; k <= i__2; ++k) {
/*<               pvec(psize*k+ip) = zero >*/
		pvec[*psize * k + ip] = 0.;
/*<             end do >*/
	    }
/*<           end if >*/
	}
/*<           if(ik.le.ksize .and. indsk(ik).eq.indsp(ip)) then >*/
	if (ik <= *ksize && indsk[ik] == indsp[ip]) {
/*<             do k=0,nrhs-1 >*/
	    i__2 = *nrhs - 1;
	    for (k = 0; k <= i__2; ++k) {
/*<               idst = psize*k+ip >*/
		idst = *psize * k + ip;
/*<               pvec(idst) = pvec(idst)+kvec(ksize*k+ik) >*/
		pvec[idst] += kvec[*ksize * k + ik];
/*<             end do >*/
	    }
/*<           end if >*/
	}
/*<         end do >*/
    }
/*<       end >*/
    return 0;
} /* x10dad_ */

/* @process nosave */
/*     recursive */
/*<       subroutine comp_sty(root,sanity,tptrs,tinds,rhsc,N,nrhs) >*/
/* Subroutine */ int comp_sty__(integer *root, doublereal *sanity, integer *
	tptrs, integer *tinds, doublereal *rhsc, integer *n, integer *nrhs)
{
    /* System generated locals */
    integer rhsc_dim1, rhsc_offset, i__1, i__2;

    /* Local variables */
    integer i__, j, kid;
    integer size;

/*<       implicit none >*/
/*<       integer N,nrhs,i,j,kid,root,size >*/
/*<       integer tptrs(3,0:N-1),tinds(N) >*/
/*<       double precision rhsc(0:N-1,nrhs),sanity >*/
/*<       size = tptrs(2,root) >*/
    /* Parameter adjustments */
    --tinds;
    --tptrs;
    rhsc_dim1 = *n - 1 - 0 + 1;
    rhsc_offset = 0 + rhsc_dim1;
    rhsc -= rhsc_offset;

    /* Function Body */
    size = tptrs[*root * 3 + 2];
/*<       do i = 0,size-1 >*/
    i__1 = size - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
/*<         kid = tinds(tptrs(1,root)+i) >*/
	kid = tinds[tptrs[*root * 3 + 1] + i__];
/*<         call comp_sty_recursive(kid,sanity,tptrs,tinds,rhsc,N,nrhs) >*/
	comp_sty__(&kid, sanity, &tptrs[1], &tinds[1], &rhsc[
		rhsc_offset], n, nrhs);
/*<         do j=1,nrhs >*/
	i__2 = *nrhs;
	for (j = 1; j <= i__2; ++j) {
/*<          sanity = sanity + rhsc(kid,j) >*/
	    *sanity += rhsc[kid + j * rhsc_dim1];
/*<         end do >*/
	}
/*<       end do >*/
    }
/*<       end >*/
    return 0;
} /* comp_sty__ */

