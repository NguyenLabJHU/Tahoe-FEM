/* $Id: pparfact1i.c,v 1.3 2005-01-04 18:19:34 paklein Exp $ */
/* pparfact1i.f -- translated by f2c (version 20030320).
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
static integer c__3 = 3;
static integer c__1 = 1;
static integer c__9 = 9;
static integer c__0 = 0;
static integer c__11 = 11;
static integer c__5 = 5;
static integer c__2 = 2;

/* /+***************************************************************************+/ */
/* /+                                                                           +/ */
/* /+   (C) Copyright IBM Corporation, 1997                                     +/ */
/* /+   (C) Copyright Regents of the University of Minnesota, 1997              +/ */
/* /+                                                                           +/ */
/* /+   pparfact1i.f                                                            +/ */
/* /+                                                                           +/ */
/* /+   Written by Anshul Gupta, IBM Corp.                                      +/ */
/* /+   Modified by Mahesh Joshi, U of MN.                                      +/ */
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
/* /+ $Id: pparfact1i.c,v 1.3 2005-01-04 18:19:34 paklein Exp $ +/ */
/* /+***************************************************************************+/ */

static integer lbit_shift(integer a, integer b) {
	return b >= 0 ? a << b : (integer)((uinteger)a >> -b);
};

static integer max(integer a, integer b) {
	return (a > b) ? a : b;
}

static integer min(integer a, integer b) {
	return (a < b) ? a : b;
}

/*<    >*/
/* Subroutine */ int pparfact1_(integer *n, integer *aptrs, integer *ainds, 
	doublereal *avals, integer *lptrs, integer *linds, doublereal *lvals, 
	integer *tptrs, integer *tinds, integer *sup, integer *stak, integer *
	nstak, integer *root, integer *dd, integer *lgblk, integer *blk, 
	integer *myid, integer *cinfo, integer *supinds, integer *supindsize, 
	doublereal *dfopts, integer *ifopts, integer *dimstak, integer *
	wsize0, integer *wsize1, integer *ibuflen, integer *dbuflen, integer *
	iwspace, integer *node, integer *stakptr, integer *nptr, integer *lc, 
	integer *iptrs, integer *info, MPI_Comm *comm)
	
/* dummy arguments for f2c
	, doublereal *dbuf_s__, 
	doublereal *wmem, integer *ibuf_s__, integer *locinds, doublereal *
	wmem0, doublereal *wmem1, doublereal *dbuf_r__, integer *ibuf_r__)
*/
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    integer diagproc;
    integer csuptr_u__;
    extern /* Subroutine */ int parelimh1_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, integer *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, MPI_Comm *);
    integer rsuptr_u__, i__, j, k, l, m;
    extern /* Subroutine */ int parelimv1_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, integer *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, MPI_Comm *);
    integer ii, jj, kk, ll, halfbuflen, /* pw0, pw1, */ ldf;
    extern /* Subroutine */ int assimilate1_(doublereal *, integer *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *);
    extern integer loc_();
    integer ldu, hdim, locf, rank, udim, vdim, ierr;
    extern /* Subroutine */ int myfc_(integer *, integer *, integer *), mydc_(
	    integer *, doublereal *, doublereal *);
    integer fptr, uptr, myup, mask1, hcube[32], hcloc, inode, jnode, vcube[32]
	    , vcloc, ncols, nrows, /* pdbufr, */ hcsize, /* pibufr, */ myleft, vcsize, 
	    locptr, mydown, csuptr, rsuptr;
    extern /* Subroutine */ int factor6_(doublereal *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, integer *, integer *);
    integer ncols_u__, currdim, dbufptr, ibufptr, partner, /* mpistat[4], */
	    msgsize, myright, nrows_u__;
	MPI_Status mpistat[1];

/*<       integer N,root,dd,lgblk,blk,myid,nstak(*),supinds(*) >*/
/*<       integer aptrs(2,0:*),ainds(*),lptrs(3,0:*),linds(*) >*/
/*<       integer tptrs(3,0:*),tinds(*),sup(*),stak(3,*) >*/
/*<       integer supindsize,info,comm >*/
/*<       integer ifopts(*),dimstak(*) >*/
/*<       integer cinfo(0:*) >*/
/*<       double precision dfopts(*) >*/
/*<       double precision lvals(*),avals(*) >*/
/*<       integer AE_TYPE_I,AE_TYPE_D >*/
/*<       parameter (AE_TYPE_I=1,AE_TYPE_D=2) >*/
/*     double precision, allocatable:: dbuf_s(:),wmem(:) */
	doublereal* dbuf_s__ = NULL;
	doublereal* wmem = NULL;
	
/*<       double precision dbuf_s(*),wmem(*)       >*/
/*<       double precision wmem0(*),wmem1(*),dbuf_r(*) >*/
	doublereal* wmem0 = NULL;
	doublereal* wmem1 = NULL;
	doublereal* dbuf_r__ = NULL;

/*     integer, allocatable:: ibuf_s(:) */
	integer* ibuf_s__ = NULL;
/*<       integer ibuf_s(*) >*/
/*     integer, allocatable:: locinds(:) */
/*<       integer locinds(*) >*/
	integer* locinds = NULL;
/*<       dimension locinds(0:N-1) >*/
/*<       integer ibuf_r(*) >*/
	integer* ibuf_r__ = NULL;
/*     pointer (pdbufr,dbuf_r),(pibufr,ibuf_r),(pw0,wmem0),(pw1,wmem1) */
/* {pdbufr, pibufr, pw0, pw1} not needed because can work with pointers directly */
/*<       integer pdbufr, pibufr, pw0, pw1 >*/
/*<       integer rsuptr, csuptr, nptr, stakptr, rank, ibufptr, dbufptr >*/
/*<       integer hdim,vdim,dbuflen,halfbuflen,udim,ldu,wsize1,wsize0 >*/
/*<       integer fptr, uptr, rsuptr_u, csuptr_u, partner, ibuflen >*/
/*<       integer vcube(0:31),hcube(0:31) >*/
/*<       integer hcsize, vcsize, currdim, myleft, myright, myup, mydown >*/
/*<       integer diagproc, vcloc, hcloc >*/
/*<       integer lc(*),iptrs(*) !Cmj >*/
/*<       include 'mpif.h' >*/
/*<       integer mpistat(MPI_STATUS_SIZE),ierr,ip,level >*/
/* -*- fortran -*- */

/* double precision functions */
/*<       info = 0 !Cmj >*/
    /* Parameter adjustments */
    --aptrs;
    --ainds;
    --avals;
    --lptrs;
    --linds;
    --lvals;
    --tptrs;
    --tinds;
    --sup;
    stak -= 4;
    --nstak;
    --supinds;
    --dfopts;
    --ifopts;
    --dimstak;
    --lc;
    --iptrs;
/*    --dbuf_s__; */
/*    --wmem; */
/*    --ibuf_s__; */
/*    --wmem0; */
/*    --wmem1; */
/*    --dbuf_r__; */
/*    --ibuf_r__; */

    /* Function Body */
    *info = 0;
/*<       diagproc = 1 >*/
    diagproc = 1;
/*<       vcsize = 1 >*/
    vcsize = 1;
/*<       hcsize = 1 >*/
    hcsize = 1;
/*<       vcube(0) = myid >*/
    vcube[0] = *myid;
/*<       hcube(0) = myid >*/
    hcube[0] = *myid;
/*<       myleft = myid >*/
    myleft = *myid;
/*<       myright = myid >*/
    myright = *myid;
/*<       myup = myid >*/
    myup = *myid;
/*<       mydown = myid >*/
    mydown = *myid;
/*<       halfbuflen = ishft(ishft(dbuflen-2,-2),1) + 1 >*/
    halfbuflen = (lbit_shift(*dbuflen - 2, (ftnlen)-2) << 1) + 1;
/*     allocate(dbuf_s(dbuflen*2),stat=i) */
/*     allocate(ibuf_s(ibuflen*2),stat=j) */
	dbuf_s__ = (doublereal*) malloc(2*(*dbuflen)*sizeof(doublereal));
	ibuf_s__ = (integer*) malloc(2*(*ibuflen)*sizeof(integer));
	if (!dbuf_s__ || !ibuf_s__) {
		printf("%d: unable to allocate working storage (dbuf+ibuf)\n", *myid);
		MPI_Abort(*comm, 0);
	}
	
	/* adjust pointers */
	dbuf_s__--;
	ibuf_s__--;
    
/*<       pdbufr = loc(dbuf_s(dbuflen+1)) >*/
/*    pdbufr = loc_(&dbuf_s__[*dbuflen + 1]); */
	dbuf_r__ = dbuf_s__ + (*dbuflen);

/*<       pibufr = loc(ibuf_s(ibuflen+1)) >*/
/*    pibufr = loc_(&ibuf_s__[*ibuflen + 1]); */
    ibuf_r__ = ibuf_s__ + (*ibuflen);
    
/*     allocate(wmem(iwspace),stat=i) */
	wmem = (doublereal*) malloc((*iwspace)*sizeof(doublereal));
	if (!wmem) {
		printf("%d: unable to allocate working storage (wmem)\n", *myid);
		MPI_Abort(*comm, 0);
	}
	wmem--; /* adjust pointer */

/*<       pw0 = loc(wmem) >*/
/*    pw0 = loc_(&wmem[1]); */
	wmem0 = wmem;
/*<       ldu = 0 >*/
    ldu = 0;
/*<    >*/
    i__1 = *iwspace + 1 << 1;
    factor6_(&wmem0[1], &linds[1], &lptrs[1], &ainds[1], &aptrs[1], &avals[1],
	     &lvals[1], &tinds[1], &tptrs[1], &sup[1], &rank, &ldu, &udim, &
	    kk, node, &i__1, &nstak[*nptr], &lc[1], &iptrs[1], &dfopts[1], &
	    ifopts[1], info);
/*<       if (info.gt.0) return ! Cmj >*/
    if (*info > 0) {
	return 0;
    }
/*<       k = sup(tptrs(3,node)) >*/
    k = sup[tptrs[*node * 3 + 3]];
/*<       ncols_u = udim  >*/
    ncols_u__ = udim;
/*<       nrows_u = udim >*/
    nrows_u__ = udim;
/*<       csuptr_u = supindsize >*/
    csuptr_u__ = *supindsize;
/*<       rsuptr_u = supindsize + udim >*/
    rsuptr_u__ = *supindsize + udim;
/*<       call myfc(udim,linds(lptrs(3,k)),supinds(supindsize)) >*/
    myfc_(&udim, &linds[lptrs[k * 3 + 3]], &supinds[*supindsize]);
/*<       call myfc(udim,supinds(supindsize),supinds(rsuptr_u)) >*/
    myfc_(&udim, &supinds[*supindsize], &supinds[rsuptr_u__]);
/*<       uptr = kk + rank * (ldu + 1) >*/
    uptr = kk + rank * (ldu + 1);
/*<       if (wsize1 .lt. uptr) then >*/
    if (*wsize1 < uptr) {

/*<       pw1 = loc(wmem) >*/
/*	pw1 = loc_(&wmem[1]); */
	wmem1 = wmem;

/*<       pw0 = loc(wmem(wsize1+1)) >*/
/*	pw0 = loc_(&wmem[*wsize1 + 1]); */
	wmem0 = wmem + (*wsize1);

/*<         uptr = uptr - wsize1 >*/
	uptr -= *wsize1;
/*<       else  >*/
    } else {
/*<       j = kk + ldu*udim >*/
	j = kk + ldu * udim;
/*<       if (iwspace-j .ge. wsize1-1) then >*/
	if (*iwspace - j >= *wsize1 - 1) {

/*<         pw1 = loc(wmem(iwspace-wsize1+1)) >*/
/*	    pw1 = loc_(&wmem[*iwspace - *wsize1 + 1]); */
	wmem1 = wmem + (*iwspace) - (*wsize1);

/*<       else >*/
	} else {

/*<         pw1 = loc(wmem(wsize0+1)) >*/
/*	    pw1 = loc_(&wmem[*wsize0 + 1]); */
		wmem1 = wmem + (*wsize0);

/*<         fptr = 1 >*/
	    fptr = 1;
/*<         j = udim - rank  >*/
	    j = udim - rank;
/*<         k = j - 1 >*/
	    k = j - 1;
/*<         do i = 0, k >*/
	    i__1 = k;
	    for (i__ = 0; i__ <= i__1; ++i__) {
/*<           do l = i, k >*/
		i__2 = k;
		for (l = i__; l <= i__2; ++l) {
/*<             wmem0(fptr+l) = wmem0(uptr+l) >*/
		    wmem0[fptr + l] = wmem0[uptr + l];
/*<           end do >*/
		}
/*<           fptr = fptr + j >*/
		fptr += j;
/*<           uptr = uptr + ldu >*/
		uptr += ldu;
/*<         end do >*/
	    }
/*<         uptr = 1 >*/
	    uptr = 1;
/*<         ldu = j >*/
	    ldu = j;
/*<       end if >*/
	}
/*<       end if >*/
    }
/*<       nptr = nptr - 1 >*/
    --(*nptr);
/*<       hdim = 0 >*/
    hdim = 0;
/*<       vdim = 0 >*/
    vdim = 0;
/*<       currdim = 0 >*/
    currdim = 0;
/*<       mask1 = blk >*/
    mask1 = *blk;
/*<       ncols_u = ncols_u - rank >*/
    ncols_u__ -= rank;
/*<       csuptr_u = csuptr_u + rank >*/
    csuptr_u__ += rank;
/*<       if (ncols_u .gt. 0) then >*/
    if (ncols_u__ > 0) {
/*<         k = supinds(csuptr_u) >*/
	k = supinds[csuptr_u__];
/*<         do 110 i = rsuptr_u, rsuptr_u + nrows_u - 1 >*/
	i__1 = rsuptr_u__ + nrows_u__ - 1;
	for (i__ = rsuptr_u__; i__ <= i__1; ++i__) {
/*<           if (supinds(i) .ge. k) goto 120 >*/
	    if (supinds[i__] >= k) {
		goto L120;
	    }
/*< 110     continue >*/
/* L110: */
	}
/*< 120     nrows_u = nrows_u - i + rsuptr_u >*/
L120:
	nrows_u__ = nrows_u__ - i__ + rsuptr_u__;
/*<         rsuptr_u = i >*/
	rsuptr_u__ = i__;
/*<       else >*/
    } else {
/*<         nrows_u = 0 >*/
	nrows_u__ = 0;
/*<       end if >*/
    }

/*     allocate(locinds(0:N-1),stat=i) */
	locinds = (integer*) malloc((*n)*sizeof(integer));
	if (!locinds) {
		printf("%d: memory allocation error", *myid);
		MPI_Abort(*comm, 0);
	}

/*<       do while (currdim .ne. dd) >*/
    while(currdim != *dd) {
/*<         node = nstak(nptr) >*/
	*node = nstak[*nptr];
/*<         partner = ieor(myid,ishft(1,currdim)) >*/
	partner = *myid ^ lbit_shift((ftnlen)1, currdim);
/*<         stakptr = stakptr - 1 >*/
	--(*stakptr);
/*<         ncols = stak(2,stakptr) >*/
	ncols = stak[*stakptr * 3 + 2];
/*<         nrows = stak(3,stakptr) >*/
	nrows = stak[*stakptr * 3 + 3];
/*<         i = tptrs(3,node) + 2 >*/
	i__ = tptrs[*node * 3 + 3] + 2;
/*<         csuptr = sup(i) >*/
	csuptr = sup[i__];
/*<         rsuptr = csuptr + sup(i+1) >*/
	rsuptr = csuptr + sup[i__ + 1];
/*<         currdim = currdim + 1 >*/
	++currdim;
/*<         ldf = dimstak(currdim) >*/
	ldf = dimstak[currdim];
/*<         if (hdim .eq. vdim) then >*/
	if (hdim == vdim) {
/*<           j = ishft(1,currdim-1) >*/
	    j = lbit_shift((ftnlen)1, currdim - 1);
/*<           if (myid .lt. partner) then >*/
	    if (*myid < partner) {
/*<             do 310 i = 0, hcsize - 1 >*/
		i__1 = hcsize - 1;
		for (i__ = 0; i__ <= i__1; ++i__) {
/*<               hcube(hcsize + i) = ieor(hcube(i), j) >*/
		    hcube[hcsize + i__] = hcube[i__] ^ j;
/*< 310         continue >*/
/* L310: */
		}
/*<           else >*/
	    } else {
/*<             do 320 i = 0, hcsize - 1 >*/
		i__1 = hcsize - 1;
		for (i__ = 0; i__ <= i__1; ++i__) {
/*<               hcube(hcsize + i) = hcube(i) >*/
		    hcube[hcsize + i__] = hcube[i__];
/*<               hcube(i) = ieor(hcube(i), j) >*/
		    hcube[i__] ^= j;
/*< 320         continue >*/
/* L320: */
		}
/*<           end if >*/
	    }
/*<           hcsize = hcsize + hcsize >*/
	    hcsize += hcsize;
/*<           hcloc = 0 >*/
	    hcloc = 0;
/*<           do while (hcube(hcloc) .ne. myid) >*/
	    while(hcube[hcloc] != *myid) {
/*<             hcloc = hcloc + 1 >*/
		++hcloc;
/*<           end do >*/
	    }
/*<           if (hcloc .eq. 0) then >*/
	    if (hcloc == 0) {
/*<             myleft = hcube(hcsize - 1) >*/
		myleft = hcube[hcsize - 1];
/*<           else >*/
	    } else {
/*<             myleft = hcube(hcloc - 1) >*/
		myleft = hcube[hcloc - 1];
/*<           end if >*/
	    }
/*<           if (hcloc .eq. hcsize-1) then >*/
	    if (hcloc == hcsize - 1) {
/*<             myright = hcube(0) >*/
		myright = hcube[0];
/*<           else  >*/
	    } else {
/*<             myright = hcube(hcloc + 1) >*/
		myright = hcube[hcloc + 1];
/*<           end if >*/
	    }
/*<           ibuf_s(1) = nrows_u >*/
	    ibuf_s__[1] = nrows_u__;
/*<           call myfc(nrows_u,supinds(rsuptr_u),ibuf_s(2)) >*/
	    myfc_(&nrows_u__, &supinds[rsuptr_u__], &ibuf_s__[2]);
/*<           ibufptr = nrows_u + 2 >*/
	    ibufptr = nrows_u__ + 2;
/*<           hdim = hdim + 1 >*/
	    ++hdim;
/*<           fptr = min0(wsize1,wsize1-ldf*(ncols-1)-nrows+1) >*/
/* Computing MIN */
	    i__1 = *wsize1, i__2 = *wsize1 - ldf * (ncols - 1) - nrows + 1;
	    fptr = min(i__1,i__2);
/*<           i = rsuptr_u >*/
	    i__ = rsuptr_u__;
/*<           j = rsuptr_u + nrows_u >*/
	    j = rsuptr_u__ + nrows_u__;
/*<           k = rsuptr >*/
	    k = rsuptr;
/*<           l = rsuptr + nrows >*/
	    l = rsuptr + nrows;
/*<           locptr = -1 >*/
	    locptr = -1;
/*<           do while (i .lt. j .and. supinds(i) .lt. supinds(k)) >*/
	    while(i__ < j && supinds[i__] < supinds[k]) {
/*<             i = i + 1 >*/
		++i__;
/*<           end do >*/
	    }
/*<           do while (i .lt. j .and. k .lt. l) >*/
	    while(i__ < j && k < l) {
/*<             locptr = locptr + 2 >*/
		locptr += 2;
/*<             do while (i .lt. j .and. supinds(i) .eq. supinds(k)) >*/
		while(i__ < j && supinds[i__] == supinds[k]) {
/*<               i = i + 1 >*/
		    ++i__;
/*<               k = k + 1 >*/
		    ++k;
/*<             end do >*/
		}
/*<             locinds(locptr) = k - 1 >*/
		locinds[locptr] = k - 1;
/*<             if (i .lt. j) then >*/
		if (i__ < j) {
/*<               do while (supinds(i) .ne. supinds(k)) >*/
		    while(supinds[i__] != supinds[k]) {
/*<                 k = k + 1 >*/
			++k;
/*<               end do >*/
		    }
/*<               locinds(locptr+1) = k - 1 >*/
		    locinds[locptr + 1] = k - 1;
/*<             else >*/
		} else {
/*<               locinds(locptr+1) = l - 1 >*/
		    locinds[locptr + 1] = l - 1;
/*<             end if >*/
		}
/*<           end do >*/
	    }
/*<           if (locptr .eq. -1) then >*/
	    if (locptr == -1) {
/*<             locinds(1) = rsuptr - 1 >*/
		locinds[1] = rsuptr - 1;
/*<             locinds(2) = l - 1 >*/
		locinds[2] = l - 1;
/*<             locptr = 1 >*/
		locptr = 1;
/*<           end if >*/
	    }
/*<           dbufptr = 1 >*/
	    dbufptr = 1;
/*<           kk = csuptr >*/
	    kk = csuptr;
/*<           ll = csuptr + ncols >*/
	    ll = csuptr + ncols;
/*<           ncols_u = ncols_u + csuptr_u >*/
	    ncols_u__ += csuptr_u__;
/*<           nrows_u = nrows_u + rsuptr_u >*/
	    nrows_u__ += rsuptr_u__;
/*<           locf = fptr >*/
	    locf = fptr;
/*<           ii = rsuptr >*/
	    ii = rsuptr;
/*<           do while (kk .lt. ll) >*/
	    while(kk < ll) {
/*<    >*/
		while(csuptr_u__ < ncols_u__ && supinds[csuptr_u__] < supinds[
			kk]) {
/*<               do while (supinds(rsuptr_u) .lt. supinds(csuptr_u)) >*/
		    while(supinds[rsuptr_u__] < supinds[csuptr_u__]) {
/*<                 rsuptr_u = rsuptr_u + 1 >*/
			++rsuptr_u__;
/*<                 uptr = uptr + 1 >*/
			++uptr;
/*<               end do >*/
		    }
/*<               m = nrows_u - rsuptr_u >*/
		    m = nrows_u__ - rsuptr_u__;
/*<               ibuf_s(ibufptr) = supinds(csuptr_u) >*/
		    ibuf_s__[ibufptr] = supinds[csuptr_u__];
/*<               ibufptr = ibufptr + 1 >*/
		    ++ibufptr;
/*<               call mydc(m,wmem0(uptr),dbuf_s(dbufptr)) >*/
		    mydc_(&m, &wmem0[uptr], &dbuf_s__[dbufptr]);
/*<               uptr = uptr + ldu >*/
		    uptr += ldu;
/*<               dbufptr = dbufptr + m  >*/
		    dbufptr += m;
/*<               csuptr_u = csuptr_u + 1 >*/
		    ++csuptr_u__;
/*<             end do >*/
		}
/*<             if (csuptr_u .lt. ncols_u) then >*/
		if (csuptr_u__ < ncols_u__) {
/*<               if (supinds(csuptr_u) .eq. supinds(kk)) then >*/
		    if (supinds[csuptr_u__] == supinds[kk]) {
/*<                 do while (supinds(rsuptr_u) .lt. supinds(csuptr_u)) >*/
			while(supinds[rsuptr_u__] < supinds[csuptr_u__]) {
/*<                   rsuptr_u = rsuptr_u + 1 >*/
			    ++rsuptr_u__;
/*<                   uptr = uptr + 1 >*/
			    ++uptr;
/*<                 end do >*/
			}
/*<                 do while (supinds(ii) .lt. supinds(kk)) >*/
			while(supinds[ii] < supinds[kk]) {
/*<                   ii = ii + 1 >*/
			    ++ii;
/*<                   locf = locf + 1 >*/
			    ++locf;
/*<                 end do >*/
			}
/*<                 i = ii >*/
			i__ = ii;
/*<                 locf = locf - ii >*/
			locf -= ii;
/*<                 j = uptr >*/
			j = uptr;
/*<                 do 130 l = 1, locptr, 2 >*/
			i__1 = locptr;
			for (l = 1; l <= i__1; l += 2) {
/*<                   j = j - i >*/
			    j -= i__;
/*<                   do 140 k = i, locinds(l) >*/
			    i__2 = locinds[l];
			    for (k = i__; k <= i__2; ++k) {
/*<                     wmem1(locf+k) = wmem0(j+k) >*/
				wmem1[locf + k] = wmem0[j + k];
/*< 140               continue >*/
/* L140: */
			    }
/*<                   j = j + k >*/
			    j += k;
/*<                   do 150 i = k, locinds(l+1) >*/
			    i__2 = locinds[l + 1];
			    for (i__ = k; i__ <= i__2; ++i__) {
/*<                     wmem1(locf+i) = 0.d0 >*/
				wmem1[locf + i__] = 0.;
/*< 150               continue >*/
/* L150: */
			    }
/*< 130             continue >*/
/* L130: */
			}
/*<                 locf = locf + ii + ldf >*/
			locf = locf + ii + ldf;
/*<                 uptr = uptr + ldu >*/
			uptr += ldu;
/*<                 csuptr_u = csuptr_u + 1 >*/
			++csuptr_u__;
/*<               else >*/
		    } else {
/*<                 do while (supinds(ii) .lt. supinds(kk)) >*/
			while(supinds[ii] < supinds[kk]) {
/*<                   ii = ii + 1 >*/
			    ++ii;
/*<                   locf = locf + 1 >*/
			    ++locf;
/*<                 end do >*/
			}
/*<                 do 160 i = locf, locinds(locptr+1)-ii+locf >*/
			i__1 = locinds[locptr + 1] - ii + locf;
			for (i__ = locf; i__ <= i__1; ++i__) {
/*<                   wmem1(i) = 0.d0 >*/
			    wmem1[i__] = 0.;
/*< 160             continue >*/
/* L160: */
			}
/*<                 locf = locf + ldf >*/
			locf += ldf;
/*<               end if >*/
		    }
/*<               kk = kk + 1 >*/
		    ++kk;
/*<             else >*/
		} else {
/*<               do while (supinds(ii) .lt. supinds(kk)) >*/
		    while(supinds[ii] < supinds[kk]) {
/*<                 ii = ii + 1 >*/
			++ii;
/*<                 locf = locf + 1 >*/
			++locf;
/*<               end do >*/
		    }
/*<               do 170 i = locf, locinds(locptr+1)-ii+locf >*/
		    i__1 = locinds[locptr + 1] - ii + locf;
		    for (i__ = locf; i__ <= i__1; ++i__) {
/*<                 wmem1(i) = 0.d0 >*/
			wmem1[i__] = 0.;
/*< 170           continue >*/
/* L170: */
		    }
/*<               locf = locf + ldf >*/
		    locf += ldf;
/*<               kk = kk + 1 >*/
		    ++kk;
/*<             end if >*/
		}
/*<           end do >*/
	    }
/*<           do while (csuptr_u .lt. ncols_u) >*/
	    while(csuptr_u__ < ncols_u__) {
/*<             do while (supinds(rsuptr_u) .lt. supinds(csuptr_u)) >*/
		while(supinds[rsuptr_u__] < supinds[csuptr_u__]) {
/*<               rsuptr_u = rsuptr_u + 1 >*/
		    ++rsuptr_u__;
/*<               uptr = uptr + 1 >*/
		    ++uptr;
/*<             end do >*/
		}
/*<             m = nrows_u - rsuptr_u >*/
		m = nrows_u__ - rsuptr_u__;
/*<             ibuf_s(ibufptr) = supinds(csuptr_u) >*/
		ibuf_s__[ibufptr] = supinds[csuptr_u__];
/*<             ibufptr = ibufptr + 1 >*/
		++ibufptr;
/*<             call mydc(m,wmem0(uptr),dbuf_s(dbufptr)) >*/
		mydc_(&m, &wmem0[uptr], &dbuf_s__[dbufptr]);
/*<             uptr = uptr + ldu >*/
		uptr += ldu;
/*<             dbufptr = dbufptr + m >*/
		dbufptr += m;
/*<             csuptr_u = csuptr_u + 1 >*/
		++csuptr_u__;
/*<           end do >*/
	    }
/*<           ibuf_s(ibufptr) = -1 >*/
	    ibuf_s__[ibufptr] = -1;
/*<           if (partner .gt. myid) then >*/
	    if (partner > *myid) {
/*<             msgsize = ibufptr >*/
		msgsize = ibufptr;
/*<    >*/

		MPI_Send(&ibuf_s__[1], msgsize, MPI_INT, partner, 1, *comm);

/*<             msgsize = ishft(dbufptr-1,3) >*/
		msgsize = dbufptr - 1 << 3;
/*<    >*/

		MPI_Send(&dbuf_s__[1], msgsize, MPI_BYTE, partner, 2, *comm);

/*<             msgsize = ibuflen >*/
		msgsize = *ibuflen;
/*<    >*/

		MPI_Recv(&ibuf_r__[1], msgsize, MPI_INT, partner, 1, *comm, mpistat);

/*<    >*/

		i__1 = *dbuflen << 3;
		MPI_Recv(&dbuf_r__[1], i__1, MPI_BYTE, partner, 2, *comm, mpistat);

/*<           else >*/
	    } else {
/*<             msgsize = ibuflen >*/
		msgsize = *ibuflen;
/*<    >*/

		MPI_Recv(&ibuf_r__[1], msgsize, MPI_INT, partner, 1, *comm, mpistat);

/*<    >*/

		i__1 = *dbuflen << 3;
		MPI_Recv(&dbuf_r__[1], i__1, MPI_BYTE, partner, 2, *comm, mpistat);

/*<             msgsize = ibufptr >*/
		msgsize = ibufptr;
/*<    >*/

		MPI_Send(&ibuf_s__[1], msgsize, MPI_INT, partner, 1, *comm);

/*<             msgsize = ishft(dbufptr-1,3) >*/
		msgsize = dbufptr - 1 << 3;
/*<    >*/

		MPI_Send(&dbuf_s__[1], msgsize, MPI_BYTE, partner, 2, *comm);

/*<           end if >*/
	    }
/*<    >*/
	    i__1 = rsuptr + nrows;
	    assimilate1_(&wmem1[fptr], &ibuf_r__[1], &dbuf_r__[1], &supinds[1]
		    , locinds, &i__1, &csuptr, &rsuptr, &ldf, ibuflen, 
		    dbuflen);
/*<    >*/
	    i__1 = halfbuflen - 1;
	    parelimv1_(&wmem1[1], &dbuf_s__[1], &dbuf_s__[halfbuflen], &
		    dbuf_r__[1], &dbuf_r__[halfbuflen], &supinds[1], &tptrs[1]
		    , &tinds[1], cinfo, &stak[4], &nstak[1], &avals[1], &
		    aptrs[1], &ainds[1], &lvals[1], &lptrs[1], &linds[1], &
		    sup[1], stakptr, nptr, &mask1, &fptr, &ldf, &nrows, &
		    ncols, &rsuptr, &csuptr, myid, &myleft, &myright, &myup, &
		    mydown, &diagproc, &i__1, wsize1, locinds, n, &ibuf_s__[1]
		    , &dfopts[1], &ifopts[1], comm);
/*<         else >*/
	} else {
/*<           j = ishft(1,currdim-1) >*/
	    j = lbit_shift((ftnlen)1, currdim - 1);
/*<           if (myid .lt. partner) then >*/
	    if (*myid < partner) {
/*<             do 410 i = 0, vcsize - 1 >*/
		i__1 = vcsize - 1;
		for (i__ = 0; i__ <= i__1; ++i__) {
/*<               vcube(vcsize + i) = ieor(vcube(i), j) >*/
		    vcube[vcsize + i__] = vcube[i__] ^ j;
/*< 410         continue >*/
/* L410: */
		}
/*<           else >*/
	    } else {
/*<             do 420 i = 0, vcsize - 1 >*/
		i__1 = vcsize - 1;
		for (i__ = 0; i__ <= i__1; ++i__) {
/*<               vcube(vcsize + i) = vcube(i) >*/
		    vcube[vcsize + i__] = vcube[i__];
/*<               vcube(i) = ieor(vcube(i), j) >*/
		    vcube[i__] ^= j;
/*< 420         continue >*/
/* L420: */
		}
/*<           end if >*/
	    }
/*<           vcsize = vcsize + vcsize >*/
	    vcsize += vcsize;
/*<           vcloc = 0 >*/
	    vcloc = 0;
/*<           do while (vcube(vcloc) .ne. myid) >*/
	    while(vcube[vcloc] != *myid) {
/*<             vcloc = vcloc + 1 >*/
		++vcloc;
/*<           end do >*/
	    }
/*<           if (vcloc .eq. 0) then >*/
	    if (vcloc == 0) {
/*<             myup = vcube(vcsize - 1) >*/
		myup = vcube[vcsize - 1];
/*<           else >*/
	    } else {
/*<             myup = vcube(vcloc - 1) >*/
		myup = vcube[vcloc - 1];
/*<           end if >*/
	    }
/*<           if (vcloc .eq. vcsize-1) then >*/
	    if (vcloc == vcsize - 1) {
/*<             mydown = vcube(0) >*/
		mydown = vcube[0];
/*<           else  >*/
	    } else {
/*<             mydown = vcube(vcloc + 1) >*/
		mydown = vcube[vcloc + 1];
/*<           end if >*/
	    }
/*<           if (hcloc .eq. vcloc) then >*/
	    if (hcloc == vcloc) {
/*<             diagproc = 1 >*/
		diagproc = 1;
/*<           else >*/
	    } else {
/*<             diagproc = 0 >*/
		diagproc = 0;
/*<           end if >*/
	    }
/*<           vdim = vdim + 1 >*/
	    ++vdim;
/*<           fptr = min0(wsize0,wsize0-ldf*(ncols-1)-nrows+1) >*/
/* Computing MIN */
	    i__1 = *wsize0, i__2 = *wsize0 - ldf * (ncols - 1) - nrows + 1;
	    fptr = min(i__1,i__2);
/*<           nrows_u = rsuptr_u + nrows_u >*/
	    nrows_u__ = rsuptr_u__ + nrows_u__;
/*<           k = rsuptr >*/
	    k = rsuptr;
/*<           l = rsuptr + nrows - 1 >*/
	    l = rsuptr + nrows - 1;
/*<           ibufptr = 2 >*/
	    ibufptr = 2;
/*<           if (nrows .gt. 0) then >*/
	    if (nrows > 0) {
/*<             inode = supinds(k) >*/
		inode = supinds[k];
/*<           else >*/
	    } else {
/*<             inode = 1000000000 >*/
		inode = 1000000000;
/*<           end if >*/
	    }
/*<           do 210 i = rsuptr_u, nrows_u - 1 >*/
	    i__1 = nrows_u__ - 1;
	    for (i__ = rsuptr_u__; i__ <= i__1; ++i__) {
/*<             jnode = supinds(i) >*/
		jnode = supinds[i__];
/*<             do while (k .lt. l) >*/
		while(k < l) {
/*<               if (inode .lt. jnode) then >*/
		    if (inode < jnode) {
/*<                 k = k + 1 >*/
			++k;
/*<                 inode = supinds(k) >*/
			inode = supinds[k];
/*<               else >*/
		    } else {
/*<                 goto 220 >*/
			goto L220;
/*<               end if >*/
		    }
/*<             end do >*/
		}
/*< 220         if (inode .eq. jnode) then >*/
L220:
		if (inode == jnode) {
/*<               locinds(i-rsuptr_u) = 1 >*/
		    locinds[i__ - rsuptr_u__] = 1;
/*<             else >*/
		} else {
/*<               ibuf_s(ibufptr) = jnode >*/
		    ibuf_s__[ibufptr] = jnode;
/*<               locinds(i-rsuptr_u) = 0 >*/
		    locinds[i__ - rsuptr_u__] = 0;
/*<               ibufptr = ibufptr + 1 >*/
		    ++ibufptr;
/*<             end if >*/
		}
/*< 210       continue  >*/
/* L210: */
	    }
/*<           ibuf_s(1) = ibufptr - 2 >*/
	    ibuf_s__[1] = ibufptr - 2;
/*<           dbufptr = 1 >*/
	    dbufptr = 1;
/*<           locf = fptr >*/
	    locf = fptr;
/*<           k = csuptr  >*/
	    k = csuptr;
/*<           l = csuptr + ncols  >*/
	    l = csuptr + ncols;
/*<           kk = rsuptr >*/
	    kk = rsuptr;
/*<           m = rsuptr + nrows >*/
	    m = rsuptr + nrows;
/*<           ii = rsuptr_u  >*/
	    ii = rsuptr_u__;
/*<           if (ncols .gt. 0) then >*/
	    if (ncols > 0) {
/*<             inode = supinds(k) >*/
		inode = supinds[k];
/*<           else >*/
	    } else {
/*<             inode = 1000000000 >*/
		inode = 1000000000;
/*<           end if >*/
	    }
/*<           do 230 i = csuptr_u, csuptr_u + ncols_u - 1 >*/
	    i__1 = csuptr_u__ + ncols_u__ - 1;
	    for (i__ = csuptr_u__; i__ <= i__1; ++i__) {
/*<             jnode = supinds(i) >*/
		jnode = supinds[i__];
/*<             do while (supinds(ii) .lt. jnode) >*/
		while(supinds[ii] < jnode) {
/*<               ii = ii + 1 >*/
		    ++ii;
/*<               uptr = uptr + 1 >*/
		    ++uptr;
/*<             end do >*/
		}
/*<             do while (k .lt. l) >*/
		while(k < l) {
/*<               if (inode .lt. jnode) then >*/
		    if (inode < jnode) {
/*<                 locf = locf - kk >*/
			locf -= kk;
/*<                 do 260 j = kk, m - 1 >*/
			i__2 = m - 1;
			for (j = kk; j <= i__2; ++j) {
/*<                   wmem0(locf+j) = 0.d0 >*/
			    wmem0[locf + j] = 0.;
/*< 260             continue >*/
/* L260: */
			}
/*<                 k = k + 1 >*/
			++k;
/*<                 locf = locf + ldf + kk >*/
			locf = locf + ldf + kk;
/*<                 inode = supinds(k) >*/
			inode = supinds[k];
/*<                 do while (supinds(kk) .lt. inode) >*/
			while(supinds[kk] < inode) {
/*<                   locf = locf + 1 >*/
			    ++locf;
/*<                   kk = kk + 1 >*/
			    ++kk;
/*<                 end do >*/
			}
/*<               else >*/
		    } else {
/*<                 goto 240 >*/
			goto L240;
/*<               end if >*/
		    }
/*<             end do >*/
		}
/*<             kk = m >*/
		kk = m;
/*< 240         jj = 0 >*/
L240:
		jj = 0;
/*<             uptr = uptr - ii >*/
		uptr -= ii;
/*<             do 250 j = ii, nrows_u - 1 >*/
		i__2 = nrows_u__ - 1;
		for (j = ii; j <= i__2; ++j) {
/*<               do 275 jj = jj, m-kk-1 >*/
		    i__3 = m - kk - 1;
		    for (jj = jj; jj <= i__3; ++jj) {
/*<                 if (supinds(kk+jj) .lt. supinds(j)) then >*/
			if (supinds[kk + jj] < supinds[j]) {
/*<                   wmem0(locf+jj) = 0.d0 >*/
			    wmem0[locf + jj] = 0.;
/*<                 else >*/
			} else {
/*<                   goto 270 >*/
			    goto L270;
/*<                 end if >*/
			}
/*< 275           continue >*/
/* L275: */
		    }
/*< 270           if (locinds(j-rsuptr_u) .eq. 0) then >*/
L270:
		    if (locinds[j - rsuptr_u__] == 0) {
/*<                 dbuf_s(dbufptr) = wmem1(uptr+j) >*/
			dbuf_s__[dbufptr] = wmem1[uptr + j];
/*<                 dbufptr = dbufptr + 1 >*/
			++dbufptr;
/*<               else  >*/
		    } else {
/*<                 wmem0(locf+jj) = wmem1(uptr+j) >*/
			wmem0[locf + jj] = wmem1[uptr + j];
/*<                 jj = jj + 1 >*/
			++jj;
/*<               end if >*/
		    }
/*< 250         continue >*/
/* L250: */
		}
/*<             do 255 jj = jj, m-kk-1 >*/
		i__2 = m - kk - 1;
		for (jj = jj; jj <= i__2; ++jj) {
/*<               wmem0(locf+jj) = 0.d0 >*/
		    wmem0[locf + jj] = 0.;
/*< 255         continue >*/
/* L255: */
		}
/*<             if (inode .eq. jnode) then >*/
		if (inode == jnode) {
/*<               k = k + 1 >*/
		    ++k;
/*<               if (k .lt. l) then >*/
		    if (k < l) {
/*<                 locf = locf + ldf >*/
			locf += ldf;
/*<                 inode = supinds(k) >*/
			inode = supinds[k];
/*<                 do while (supinds(kk) .lt. inode) >*/
			while(supinds[kk] < inode) {
/*<                   locf = locf + 1 >*/
			    ++locf;
/*<                   kk = kk + 1 >*/
			    ++kk;
/*<                 end do >*/
			}
/*<               else >*/
		    } else {
/*<                 kk = m >*/
			kk = m;
/*<               end if >*/
		    }
/*<             end if >*/
		}
/*<             uptr = uptr + ldu + ii >*/
		uptr = uptr + ldu + ii;
/*<             ibuf_s(ibufptr) = jnode >*/
		ibuf_s__[ibufptr] = jnode;
/*<             ibufptr = ibufptr + 1 >*/
		++ibufptr;
/*< 230       continue >*/
/* L230: */
	    }
/*<           do 280 k = k, l-1, 1 >*/
	    i__1 = l - 1;
	    for (k = k; k <= i__1; ++k) {
/*<             do while (supinds(kk) .lt. supinds(k)) >*/
		while(supinds[kk] < supinds[k]) {
/*<               locf = locf + 1 >*/
		    ++locf;
/*<               kk = kk + 1 >*/
		    ++kk;
/*<             end do >*/
		}
/*<             locf = locf - kk >*/
		locf -= kk;
/*<             do 290 j = kk, m - 1 >*/
		i__2 = m - 1;
		for (j = kk; j <= i__2; ++j) {
/*<               wmem0(locf+j) = 0.d0 >*/
		    wmem0[locf + j] = 0.;
/*< 290         continue >*/
/* L290: */
		}
/*<             locf = locf + ldf + kk >*/
		locf = locf + ldf + kk;
/*< 280       continue >*/
/* L280: */
	    }
/*<           i = ibuf_s(1) >*/
	    i__ = ibuf_s__[1];
/*<           if (i .gt. 0) then >*/
	    if (i__ > 0) {
/*<             do while (ibuf_s(ibufptr-1) .gt. ibuf_s(i+1)) >*/
		while(ibuf_s__[ibufptr - 1] > ibuf_s__[i__ + 1]) {
/*<               ibufptr = ibufptr - 1 >*/
		    --ibufptr;
/*<             end do >*/
		}
/*<           end if >*/
	    }
/*<           ibuf_s(ibufptr) = -1  >*/
	    ibuf_s__[ibufptr] = -1;
/*<           if (partner .gt. myid) then >*/
	    if (partner > *myid) {
/*<             msgsize = ibufptr >*/
		msgsize = ibufptr;
/*<    >*/

		MPI_Send(&ibuf_s__[1], msgsize, MPI_INT, partner, 1, *comm);

/*<             msgsize = ishft(dbufptr-1,3) >*/
		msgsize = dbufptr - 1 << 3;
/*<    >*/

		MPI_Send(&dbuf_s__[1], msgsize, MPI_BYTE, partner, 2, *comm);
		
/*<             msgsize = ibuflen >*/
		msgsize = *ibuflen;
/*<    >*/

		MPI_Recv(&ibuf_r__[1], msgsize, MPI_INT, partner, 1, *comm, mpistat);

/*<    >*/

		i__1 = *dbuflen << 3;
		MPI_Recv(&dbuf_r__[1], i__1, MPI_BYTE, partner, 2, *comm, mpistat);

/*<           else >*/
	    } else {
/*<             msgsize = ibuflen  >*/
		msgsize = *ibuflen;
/*<    >*/

		MPI_Recv(&ibuf_r__[1], msgsize, MPI_INT, partner, 1, *comm, mpistat);

/*<    >*/

		i__1 = *dbuflen << 3;
		MPI_Recv(&dbuf_r__[1], i__1, MPI_BYTE, partner, 2, *comm, mpistat);

/*<             msgsize = ibufptr >*/
		msgsize = ibufptr;
/*<    >*/

		MPI_Send(&ibuf_s__[1], msgsize, MPI_INT, partner, 1, *comm);

/*<             msgsize = ishft(dbufptr-1,3) >*/
		msgsize = dbufptr - 1 << 3;
/*<    >*/

		MPI_Send(&dbuf_s__[1], msgsize, MPI_BYTE, partner, 2, *comm);

/*<           end if >*/
	    }
/*<    >*/
	    i__1 = rsuptr + nrows;
	    assimilate1_(&wmem0[fptr], &ibuf_r__[1], &dbuf_r__[1], &supinds[1]
		    , locinds, &i__1, &csuptr, &rsuptr, &ldf, ibuflen, 
		    dbuflen);
/*<    >*/
	    i__1 = halfbuflen - 1;
	    parelimh1_(&wmem0[1], &dbuf_s__[1], &dbuf_s__[halfbuflen], &
		    dbuf_r__[1], &dbuf_r__[halfbuflen], &supinds[1], &tptrs[1]
		    , &tinds[1], cinfo, &stak[4], &nstak[1], &avals[1], &
		    aptrs[1], &ainds[1], &lvals[1], &lptrs[1], &linds[1], &
		    sup[1], stakptr, nptr, &mask1, &fptr, &ldf, &nrows, &
		    ncols, &rsuptr, &csuptr, myid, &myleft, &myright, &myup, &
		    mydown, &diagproc, &i__1, wsize0, locinds, &ibuf_s__[1], &
		    dfopts[1], &ifopts[1], comm);
/*<           mask1 = ior(mask1,ishft(mask1,1)) >*/
	    mask1 |= mask1 << 1;
/*<         end if >*/
	}
/*<         uptr = fptr >*/
	uptr = fptr;
/*<         ldu = ldf >*/
	ldu = ldf;
/*<         ncols_u = ncols >*/
	ncols_u__ = ncols;
/*<         nrows_u = nrows >*/
	nrows_u__ = nrows;
/*<         csuptr_u = csuptr >*/
	csuptr_u__ = csuptr;
/*<         rsuptr_u = rsuptr >*/
	rsuptr_u__ = rsuptr;
/*<       end do >*/
    }
    
    /* restore pointers */
    dbuf_s__++;
    ibuf_s__++;
    wmem++;
    
/*     deallocate(dbuf_s) */
/*     deallocate(ibuf_s) */
/*     deallocate(wmem) */
/*     deallocate(locinds) */
	free(dbuf_s__);
    free(ibuf_s__);
    free(wmem);
	free(locinds);

/*<       return >*/
    return 0;

#if 0
/*< 111   print *,'Bad news in serial factor!' >*/
/* L111: */
    s_wsle(&io___57);
    do_lio(&c__9, &c__1, "Bad news in serial factor!", (ftnlen)26);
    e_wsle();
/*<       write(*,*) (ifopts(i),i=1,5) >*/
    s_wsle(&io___58);
    for (i__ = 1; i__ <= 5; ++i__) {
	do_lio(&c__3, &c__1, (char *)&ifopts[i__], (ftnlen)sizeof(integer));
    }
    e_wsle();
/*<       do i = 1, 7 >*/
    for (i__ = 1; i__ <= 7; ++i__) {
/*<         print *, dfopts(i) >*/
	s_wsle(&io___59);
	do_lio(&c__5, &c__1, (char *)&dfopts[i__], (ftnlen)sizeof(doublereal))
		;
	e_wsle();
/*<       end do >*/
    }
/*<       return >*/
    return 0;
/*<       end >*/
#endif
} /* pparfact1_ */
