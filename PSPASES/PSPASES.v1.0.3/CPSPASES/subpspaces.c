/* subpspaces.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

/* debugging */
#undef __DO_DEBUG__
/* #define __DO_DEBUG__ 1 */

#include "mpi.h"
#include "pspases_f2c.h"
#include <stdio.h>
#include <math.h>

/* Table of constant values */
static integer c__0 = 0;
static doublereal c_b57 = -1.;
static doublereal c_b58 = 1.;
static integer c__1 = 1;

/* /+***************************************************************************+/ */
/* /+                                                                           +/ */
/* /+   (C) Copyright IBM Corporation, 1997                                     +/ */
/* /+   (C) Copyright Regents of the University of Minnesota, 1997              +/ */
/* /+                                                                           +/ */
/* /+   subpspaces.f                                                            +/ */
/* /+                                                                           +/ */
/* /+   Written by Anshul Gupta, IBM Corp.                                      +/ */
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
/* /+ $Id: subpspaces.c,v 1.1 2004-12-28 17:46:01 paklein Exp $ +/ */
/* /+***************************************************************************+/ */
/* ------------------------------------------------------------------------------- */

/*<       subroutine DPACK(sbuf,dbuf,blklen,offset,blknum) >*/
/* Subroutine */ int dpack_(doublereal *sbuf, doublereal *dbuf, integer *
	blklen, integer *offset, integer *blknum)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__;
    extern /* Subroutine */ int mydc_(integer *, doublereal *, doublereal *);
    integer istart, jstart;

/*<       double precision sbuf(*), dbuf(*) >*/
/*<       integer blklen,offset,blknum >*/
/*<       jstart = 1 >*/
    /* Parameter adjustments */
    --dbuf;
    --sbuf;

    /* Function Body */
    jstart = 1;
/*<       istart = 1 >*/
    istart = 1;
/*<       do 10 i = 1, blknum >*/
    i__1 = *blknum;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         call mydc(blklen,sbuf(jstart),dbuf(istart)) >*/
	mydc_(blklen, &sbuf[jstart], &dbuf[istart]);
/*<         jstart = jstart + offset >*/
	jstart += *offset;
/*<         istart = istart + blklen >*/
	istart += *blklen;
/*< 10    continue >*/
/* L10: */
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* dpack_ */

/* ------------------------------------------------------------------------------- */
/*<    >*/
/* Subroutine */ int assimilate1_(doublereal *wmem, integer *ibuf, doublereal 
	*dbuf, integer *supinds, integer *locinds, integer *nrowlim, integer *
	csuptr, integer *rsuptr, integer *ldf, integer *ibuflen, integer *
	dbuflen)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    integer i__, k, l, jc, jr, locf, locptr, dbufptr, ibufptr, nextcol;

/*<       integer ibuf(*),supinds(*),locinds(*) >*/
/*<       integer nrowlim,csuptr,rsuptr,ldf,ibuflen,dbuflen >*/
/*<       double precision wmem(*), dbuf(*) >*/
/*<       integer ibufptr, dbufptr, locf, locptr, jc, jr >*/
/*<       i = ibuf(1) >*/
    /* Parameter adjustments */
    --locinds;
    --supinds;
    --dbuf;
    --ibuf;
    --wmem;

    /* Function Body */
    i__ = ibuf[1];
/*<       ibufptr = i + 2 >*/
    ibufptr = i__ + 2;
/*<       if (i .eq. 0 .or. ibuf(ibufptr) .eq. -1) then >*/
    if (i__ == 0 || ibuf[ibufptr] == -1) {
/*<       return >*/
	return 0;
/*<       end if >*/
    }
/*<       dbufptr = 1 >*/
    dbufptr = 1;
/*<       locf = 1 >*/
    locf = 1;
/*<       locptr = -1 >*/
    locptr = -1;
/*<       jc = csuptr >*/
    jc = *csuptr;
/*<       jr = rsuptr >*/
    jr = *rsuptr;
/*<       nextcol = ibuf(ibufptr) >*/
    nextcol = ibuf[ibufptr];
/*<       i = 2 >*/
    i__ = 2;
/*<       do while (ibuf(i) .lt. nextcol) >*/
    while(ibuf[i__] < nextcol) {
/*<         i = i + 1 >*/
	++i__;
/*<       end do >*/
    }
/*<       k = rsuptr >*/
    k = *rsuptr;
/*<       do while (i .lt. ibufptr .and. k .lt. nrowlim) >*/
    while(i__ < ibufptr && k < *nrowlim) {
/*<         locptr = locptr + 2 >*/
	locptr += 2;
/*<         do while (i .lt. ibufptr .and. ibuf(i) .eq. supinds(k)) >*/
	while(i__ < ibufptr && ibuf[i__] == supinds[k]) {
/*<           i = i + 1 >*/
	    ++i__;
/*<           k = k + 1 >*/
	    ++k;
/*<         end do >*/
	}
/*<         locinds(locptr) = k - 1 >*/
	locinds[locptr] = k - 1;
/*<         if (i .lt. ibufptr) then >*/
	if (i__ < ibufptr) {
/*<           do while (ibuf(i) .ne. supinds(k)) >*/
	    while(ibuf[i__] != supinds[k]) {
/*<             k = k + 1 >*/
		++k;
/*<           end do >*/
	    }
/*<           locinds(locptr+1) = k - 1 >*/
	    locinds[locptr + 1] = k - 1;
/*<         else >*/
	} else {
/*<           locinds(locptr+1) = nrowlim - 1 >*/
	    locinds[locptr + 1] = *nrowlim - 1;
/*<         end if >*/
	}
/*<       end do >*/
    }
/*<       do while (nextcol .ge. 0) >*/
    while(nextcol >= 0) {
/*<         do while (supinds(jc) .ne. nextcol) >*/
	while(supinds[jc] != nextcol) {
/*<           jc = jc + 1 >*/
	    ++jc;
/*<           locf = locf + ldf >*/
	    locf += *ldf;
/*<         end do >*/
	}
/*<         do while (supinds(jr) .lt. nextcol) >*/
	while(supinds[jr] < nextcol) {
/*<           jr = jr + 1 >*/
	    ++jr;
/*<           locf = locf + 1 >*/
	    ++locf;
/*<         end do >*/
	}
/*<         k = jr >*/
	k = jr;
/*<         locf = locf - jr >*/
	locf -= jr;
/*<         do 10 l = 1, locptr, 2 >*/
	i__1 = locptr;
	for (l = 1; l <= i__1; l += 2) {
/*<           dbufptr = dbufptr - k >*/
	    dbufptr -= k;
/*<           do 25 k = k, locinds(l) >*/
	    i__2 = locinds[l];
	    for (k = k; k <= i__2; ++k) {
/*<             wmem(locf+k) = wmem(locf+k) + dbuf(dbufptr+k) >*/
		wmem[locf + k] += dbuf[dbufptr + k];
/*< 25        continue >*/
/* L25: */
	    }
/*<           dbufptr = dbufptr + k >*/
	    dbufptr += k;
/*<           k = max0(k,locinds(l+1)+1) >*/
/* Computing MAX */
	    i__2 = k, i__3 = locinds[l + 1] + 1;
	    k = max(i__2,i__3);
/*< 10      continue >*/
/* L10: */
	}
/*<         locf = locf + jr >*/
	locf += jr;
/*<         ibufptr = ibufptr + 1 >*/
	++ibufptr;
/*<         nextcol = ibuf(ibufptr) >*/
	nextcol = ibuf[ibufptr];
/*<       end do >*/
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* assimilate1_ */

/* ------------------------------------------------------------------------------- */
/*<       SUBROUTINE MYFC (N,DX,DY) >*/
/* Subroutine */ int myfc_(integer *n, integer *dx, integer *dy)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__;

/*<       integer dx(*),dy(*) >*/
/*<       do i=1,n >*/
    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         dy(i) = dx(i) >*/
	dy[i__] = dx[i__];
/*<       end do >*/
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* myfc_ */

/* ------------------------------------------------------------------------------- */
/*<       SUBROUTINE MYDC (N,DX,DY) >*/
/* Subroutine */ int mydc_(integer *n, doublereal *dx, doublereal *dy)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__;

/*<       integer N >*/
/*<       double precision dx(*),dy(*) >*/
/*<       do i = 1, n >*/
    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       dy(i) = dx(i) >*/
	dy[i__] = dx[i__];
/*<       end do >*/
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* mydc_ */

/* ------------------------------------------------------------------------------- */
/*<       SUBROUTINE MYDDC (N,DX,DY,DZ) >*/
/* Subroutine */ int myddc_(integer *n, doublereal *dx, doublereal *dy, 
	doublereal *dz)
{
    extern /* Subroutine */ int mydc_(integer *, doublereal *, doublereal *);

/*<       integer N >*/
/*<       double precision DX, DY, DZ >*/
/*<       call mydc(n,dx,dy) >*/
    mydc_(n, dx, dy);
/*<       call mydc(n,dx,dz) >*/
    mydc_(n, dx, dz);
/*<       return >*/
    return 0;
/*<       end >*/
} /* myddc_ */

/* ------------------------------------------------------------------------------ */
/*<       subroutine a2i(i,a,LEN) >*/
/* Subroutine */ int a2i_(integer *i__, char *a, integer *len, ftnlen a_len)
{
    integer j, l;

/*<       integer i >*/
/*<       character a(*) >*/
/*      parameter(LEN = 64) */
/*<       i = 0 >*/
    /* Parameter adjustments */
    --a;

    /* Function Body */
    *i__ = 0;
/*<       j = 1 >*/
    j = 1;
/*<       l = 1 >*/
    l = 1;
/*<    >*/
    while(j < *len && *(unsigned char *)&a[j] != '0' && *(unsigned char *)&a[
	    j] != '1' && *(unsigned char *)&a[j] != '2' && *(unsigned char *)&
	    a[j] != '3' && *(unsigned char *)&a[j] != '4' && *(unsigned char *
	    )&a[j] != '5' && *(unsigned char *)&a[j] != '6' && *(unsigned 
	    char *)&a[j] != '7' && *(unsigned char *)&a[j] != '8' && *(
	    unsigned char *)&a[j] != '9') {
/*<         if (a(j) .eq. '-') l = -1 >*/
	if (*(unsigned char *)&a[j] == '-') {
	    l = -1;
	}
/*<         j = j + 1 >*/
	++j;
/*<       end do >*/
    }
/*<    >*/
    while(j < *len && (*(unsigned char *)&a[j] == '0' || *(unsigned char *)&a[
	    j] == '1' || *(unsigned char *)&a[j] == '2' || *(unsigned char *)&
	    a[j] == '3' || *(unsigned char *)&a[j] == '4' || *(unsigned char *
	    )&a[j] == '5' || *(unsigned char *)&a[j] == '6' || *(unsigned 
	    char *)&a[j] == '7' || *(unsigned char *)&a[j] == '8' || *(
	    unsigned char *)&a[j] == '9')) {
/*<         i = i * 10 >*/
	*i__ *= 10;
/*<         if (a(j) .eq. '1') i = i + 1 >*/
	if (*(unsigned char *)&a[j] == '1') {
	    ++(*i__);
	}
/*<         if (a(j) .eq. '2') i = i + 2 >*/
	if (*(unsigned char *)&a[j] == '2') {
	    *i__ += 2;
	}
/*<         if (a(j) .eq. '3') i = i + 3 >*/
	if (*(unsigned char *)&a[j] == '3') {
	    *i__ += 3;
	}
/*<         if (a(j) .eq. '4') i = i + 4 >*/
	if (*(unsigned char *)&a[j] == '4') {
	    *i__ += 4;
	}
/*<         if (a(j) .eq. '5') i = i + 5 >*/
	if (*(unsigned char *)&a[j] == '5') {
	    *i__ += 5;
	}
/*<         if (a(j) .eq. '6') i = i + 6 >*/
	if (*(unsigned char *)&a[j] == '6') {
	    *i__ += 6;
	}
/*<         if (a(j) .eq. '7') i = i + 7 >*/
	if (*(unsigned char *)&a[j] == '7') {
	    *i__ += 7;
	}
/*<         if (a(j) .eq. '8') i = i + 8 >*/
	if (*(unsigned char *)&a[j] == '8') {
	    *i__ += 8;
	}
/*<         if (a(j) .eq. '9') i = i + 9 >*/
	if (*(unsigned char *)&a[j] == '9') {
	    *i__ += 9;
	}
/*<         j = j + 1 >*/
	++j;
/*<       end do >*/
    }
/*<       i = i * l >*/
    *i__ *= l;
/*<       return >*/
    return 0;
/*<       end >*/
} /* a2i_ */

/* ------------------------------------------------------------------------------ */
/*<       subroutine compmysan(sanity,n,lptrs,lvals,linds,cinfo) >*/
/* Subroutine */ int compmysan_(doublereal *sanity, integer *n, integer *
	lptrs, doublereal *lvals, integer *linds, integer *cinfo)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    integer is1, is2;

/*      include 'mpif.h' */
/*<       integer n, lptrs(3,0:*) >*/
/*<       double precision sanity, lvals(*) >*/
/*<       integer cinfo(0:n-1) >*/
/*      call mpi_comm_rank(MPI_COMM_WORLD,myid,ierr) */
/*      write(10+myid,*)'cinfo:',(cinfo(i),i=0,n-1) */
/*<       sanity = 0.d0 >*/
    /* Parameter adjustments */
    --lptrs;
    --lvals;

    /* Function Body */
    *sanity = 0.;
/*<       do 1130 is1 = 0, N - 1 >*/
    i__1 = *n - 1;
    for (is1 = 0; is1 <= i__1; ++is1) {
/*<           if (cinfo(is1) .eq. 1) then >*/
	if (cinfo[is1] == 1) {
/*      write(10+myid,*)is1,':',(lvals(is2), */
/*     +			       linds(is2-lptrs(1,is1)+lptrs(3,is1)), */
/*     +			       is2 = lptrs(1,is1), */
/*     +				     lptrs(1,is1)+lptrs(2,is1)-1) */
/*<           do 1170 is2 = lptrs(1,is1), lptrs(1,is1)+lptrs(2,is1)-1 >*/
	    i__2 = lptrs[is1 * 3 + 1] + lptrs[is1 * 3 + 2] - 1;
	    for (is2 = lptrs[is1 * 3 + 1]; is2 <= i__2; ++is2) {
/*<             sanity = sanity + lvals(is2) >*/
		*sanity += lvals[is2];
/*<  1170     continue >*/
/* L1170: */
	    }
/*<         end if >*/
	}
/*< 1130  continue >*/
/* L1130: */
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* compmysan_ */

/* ------------------------------------------------------------------------------ */
/*    recursive */
/*    + */
/*<    >*/
/* Subroutine */ int symbolic6_(integer *aptrs, integer *ainds, integer *
	lptrs, integer *linds, integer *sup, integer *myscrach, integer *
	tptrs, integer *tinds, integer *nnz, integer *root, integer *scount, 
	integer *lspace, doublereal *opcount, integer *mystak, integer *lctr, 
	integer *nnzxtra, doublereal *opxtra, doublereal *ssthresh, integer *
	lsize, integer *info)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    integer isupsize, i__, j, k, ii, jj, kk, nnd, lds, jml;
    doublereal xop;
    integer node, istk, jlctr, itscr, nmiss;
    extern /* Subroutine */ int smerge_(integer *, integer *, integer *, 
	    integer *, integer *, integer *);
    integer imyscr, supbot, suptop, isupbot;

/*<       integer scrsize,mystak(*),lspace >*/
/*<       integer aptrs(2,0:*),ainds(*),lptrs(3,0:*),linds(*) >*/
/*<       integer sup(*),tptrs(3,0:*),tinds(*) >*/
/*<       integer jlctr,nnz,root,scount,lctr,nnzxtra >*/
/*<       double precision opcount,opxtra,xop,ssthresh >*/
/*<       integer myscrach(0:*) >*/
/*<       integer i,j,k,node,imyscr,itscr,supbot,suptop,nmiss >*/
/*<       integer lsize,info !Cmj >*/
/*<       if(info.eq.1) return !Cmj >*/
    /* Parameter adjustments */
    --mystak;
    --tinds;
    --tptrs;
    --sup;
    --linds;
    --lptrs;
    --ainds;
    --aptrs;

    /* Function Body */
    if (*info == 1) {
	return 0;
    }
/*<       istk = 1 >*/
    istk = 1;
/*<       node = root >*/
    node = *root;
/*<       k = tptrs(2,node) >*/
    k = tptrs[node * 3 + 2];
/*<       do while (k .eq. 1) >*/
    while(k == 1) {
/*<         mystak(istk) = node >*/
	mystak[istk] = node;
/*<         node = tinds(tptrs(1,node)) >*/
	node = tinds[tptrs[node * 3 + 1]];
/*<         k = tptrs(2,node) >*/
	k = tptrs[node * 3 + 2];
/*<         istk = istk + 1 >*/
	++istk;
/*<       end do >*/
    }
/*<       mystak(istk) = node >*/
    mystak[istk] = node;
/*<       if (k .eq. 0) then >*/
    if (k == 0) {
/*<       jj = aptrs(1,node) >*/
	jj = aptrs[(node << 1) + 1];
/*<         do 10 itscr = 0, aptrs(2,node) - 1 >*/
	i__1 = aptrs[(node << 1) + 2] - 1;
	for (itscr = 0; itscr <= i__1; ++itscr) {
/*<           linds(lctr+itscr) = ainds(jj+itscr) >*/
	    linds[*lctr + itscr] = ainds[jj + itscr];
/*< 10      continue >*/
/* L10: */
	}
/*<       else >*/
    } else {
/*<       do i = 0, k - 1 >*/
	i__1 = k - 1;
	for (i__ = 0; i__ <= i__1; ++i__) {
/*<    >*/
	    symbolic6_(&aptrs[1], &ainds[1], &lptrs[1], &linds[1], 
		    &sup[1], myscrach, &tptrs[1], &tinds[1], nnz, &tinds[
		    tptrs[node * 3 + 1] + i__], scount, lspace, opcount, &
		    mystak[istk + 1], lctr, nnzxtra, opxtra, ssthresh, lsize, 
		    info);
/*<           if(info.eq.1) return !Cmj >*/
	    if (*info == 1) {
		return 0;
	    }
/*<       end do >*/
	}
/*<       i = tinds(tptrs(1,node)) >*/
	i__ = tinds[tptrs[node * 3 + 1]];
/*<       j = tinds(tptrs(1,node)+1) >*/
	j = tinds[tptrs[node * 3 + 1] + 1];
/*<    >*/
	i__1 = lptrs[i__ * 3 + 2] - 1;
	i__2 = lptrs[j * 3 + 2] - 1;
	smerge_(&i__1, &linds[lptrs[i__ * 3 + 3] + 1], &i__2, &linds[lptrs[j *
		 3 + 3] + 1], &itscr, myscrach);
/*<         do i = 2, k - 2, 2 >*/
	i__1 = k - 2;
	for (i__ = 2; i__ <= i__1; i__ += 2) {
/*<         j = tinds(tptrs(1,node)+i) >*/
	    j = tinds[tptrs[node * 3 + 1] + i__];
/*<         ii = itscr >*/
	    ii = itscr;
/*<    >*/
	    i__2 = lptrs[j * 3 + 2] - 1;
	    smerge_(&ii, myscrach, &i__2, &linds[lptrs[j * 3 + 3] + 1], &
		    itscr, &linds[*lctr]);
/*<         j = tinds(tptrs(1,node)+i+1) >*/
	    j = tinds[tptrs[node * 3 + 1] + i__ + 1];
/*<         ii = itscr >*/
	    ii = itscr;
/*<    >*/
	    i__2 = lptrs[j * 3 + 2] - 1;
	    smerge_(&ii, &linds[*lctr], &i__2, &linds[lptrs[j * 3 + 3] + 1], &
		    itscr, myscrach);
/*<       end do >*/
	}
/*<         if (i .eq. k - 1) then >*/
	if (i__ == k - 1) {
/*<         j = tinds(tptrs(1,node)+i) >*/
	    j = tinds[tptrs[node * 3 + 1] + i__];
/*<         ii = itscr >*/
	    ii = itscr;
/*<    >*/
	    i__1 = lptrs[j * 3 + 2] - 1;
	    smerge_(&ii, myscrach, &i__1, &linds[lptrs[j * 3 + 3] + 1], &
		    itscr, &linds[*lctr]);
/*<         do j = 0, itscr - 1 >*/
	    i__1 = itscr - 1;
	    for (j = 0; j <= i__1; ++j) {
/*<           myscrach(j) = linds(lctr+j) >*/
		myscrach[j] = linds[*lctr + j];
/*<         end do >*/
	    }
/*<         end if >*/
	}
/*<         imyscr = itscr >*/
	imyscr = itscr;
/*<         itscr = 0 >*/
	itscr = 0;
/*<         ii = aptrs(1,node) >*/
	ii = aptrs[(node << 1) + 1];
/*<         kk = aptrs(2,node) + ii >*/
	kk = aptrs[(node << 1) + 2] + ii;
/*<         jj = 0 >*/
	jj = 0;
/*<         do while (ii .lt. kk .and. jj .lt. imyscr) >*/
	while(ii < kk && jj < imyscr) {
/*<           if (ainds(ii) .gt. myscrach(jj)) then >*/
	    if (ainds[ii] > myscrach[jj]) {
/*<             linds(lctr+itscr) = myscrach(jj) >*/
		linds[*lctr + itscr] = myscrach[jj];
/*<             jj = jj + 1 >*/
		++jj;
/*<           else >*/
	    } else {
/*<             if (ainds(ii) .eq. myscrach(jj)) then >*/
		if (ainds[ii] == myscrach[jj]) {
/*<               linds(lctr+itscr) = myscrach(jj) >*/
		    linds[*lctr + itscr] = myscrach[jj];
/*<               jj = jj + 1 >*/
		    ++jj;
/*<               ii = ii + 1 >*/
		    ++ii;
/*<             else >*/
		} else {
/*<               linds(lctr+itscr) = ainds(ii) >*/
		    linds[*lctr + itscr] = ainds[ii];
/*<               ii = ii + 1 >*/
		    ++ii;
/*<             end if >*/
		}
/*<           end if >*/
	    }
/*<           itscr = itscr + 1 >*/
	    ++itscr;
/*<         end do >*/
	}
/*<         do 444 ii = ii, kk - 1 >*/
	i__1 = kk - 1;
	for (ii = ii; ii <= i__1; ++ii) {
/*<           linds(lctr+itscr) = ainds(ii) >*/
	    linds[*lctr + itscr] = ainds[ii];
/*<           itscr = itscr + 1 >*/
	    ++itscr;
/*< 444     continue >*/
/* L444: */
	}
/*<         do 555 jj = jj, imyscr - 1 >*/
	i__1 = imyscr - 1;
	for (jj = jj; jj <= i__1; ++jj) {
/*<           linds(lctr+itscr) = myscrach(jj) >*/
	    linds[*lctr + itscr] = myscrach[jj];
/*<           itscr = itscr + 1 >*/
	    ++itscr;
/*< 555     continue >*/
/* L555: */
	}
/*<       end if >*/
    }
/*<       lptrs(1,node) = lspace >*/
    lptrs[node * 3 + 1] = *lspace;
/*<       lptrs(2,node) = itscr >*/
    lptrs[node * 3 + 2] = itscr;
/*<       lptrs(3,node) = lctr >*/
    lptrs[node * 3 + 3] = *lctr;
/*<       jlctr = lctr + 1 >*/
    jlctr = *lctr + 1;
/*<       lctr = lctr + itscr >*/
    *lctr += itscr;
/*<       if(lctr.gt.lsize) then >*/
    if (*lctr > *lsize) {
/*<       info = 1 >*/
	*info = 1;
/*<       return >*/
	return 0;
/*<       end if >*/
    }
/*<       nnz = nnz + itscr >*/
    *nnz += itscr;
/*<       opcount = opcount + dble(itscr * itscr) >*/
    *opcount += (doublereal) (itscr * itscr);
/*<       supbot = node >*/
    supbot = node;
/*<       isupbot = istk >*/
    isupbot = istk;
/*<       isupsize = 0 >*/
    isupsize = 0;
/*<       lds = itscr >*/
    lds = itscr;
/*<       lspace = lspace + lds >*/
    *lspace += lds;
/*<       do 30 j = istk-1, 1, -1 >*/
    for (j = istk - 1; j >= 1; --j) {
/*<         suptop = node >*/
	suptop = node;
/*<         isupsize = isupsize + 1 >*/
	++isupsize;
/*<         itscr = 0 >*/
	itscr = 0;
/*<         nmiss = 0 >*/
	nmiss = 0;
/*<       jml = -1 >*/
	jml = -1;
/*<         node = mystak(j) >*/
	node = mystak[j];
/*<         ii = jlctr >*/
	ii = jlctr;
/*<         jj = aptrs(1,node) >*/
	jj = aptrs[(node << 1) + 1];
/*<         k = aptrs(2,node) + jj >*/
	k = aptrs[(node << 1) + 2] + jj;
/*<         do while (jj .lt. k .and. ii .lt. lctr) >*/
	while(jj < k && ii < *lctr) {
/*<           if (ainds(jj) .gt. linds(ii)) then >*/
	    if (ainds[jj] > linds[ii]) {
/*<             linds(lctr+itscr) = linds(ii) >*/
		linds[*lctr + itscr] = linds[ii];
/*<             ii = ii + 1 >*/
		++ii;
/*<           else >*/
	    } else {
/*<             if (ainds(jj) .eq. linds(ii)) then >*/
		if (ainds[jj] == linds[ii]) {
/*<               linds(lctr+itscr) = linds(ii) >*/
		    linds[*lctr + itscr] = linds[ii];
/*<               ii = ii + 1 >*/
		    ++ii;
/*<               jj = jj + 1 >*/
		    ++jj;
/*<             else >*/
		} else {
/*<             if (nmiss .eq. 0) jml = itscr >*/
		    if (nmiss == 0) {
			jml = itscr;
		    }
/*<               linds(lctr+itscr) = ainds(jj) >*/
		    linds[*lctr + itscr] = ainds[jj];
/*<               jj = jj + 1 >*/
		    ++jj;
/*<               nmiss = nmiss + 1 >*/
		    ++nmiss;
/*<             end if >*/
		}
/*<           end if >*/
	    }
/*<           itscr = itscr + 1 >*/
	    ++itscr;
/*<         end do >*/
	}
/*<       if (jml .eq. -1 .and. k .gt. jj) jml = itscr >*/
	if (jml == -1 && k > jj) {
	    jml = itscr;
	}
/*<         nmiss = nmiss + k - jj >*/
	nmiss = nmiss + k - jj;
/*<         do 70 jj = jj, k - 1 >*/
	i__1 = k - 1;
	for (jj = jj; jj <= i__1; ++jj) {
/*<             linds(lctr+itscr) = ainds(jj) >*/
	    linds[*lctr + itscr] = ainds[jj];
/*<             itscr = itscr + 1 >*/
	    ++itscr;
/*< 70      continue >*/
/* L70: */
	}
/*<         do 80 ii = ii, lctr - 1 >*/
	i__1 = *lctr - 1;
	for (ii = ii; ii <= i__1; ++ii) {
/*<           linds(lctr+itscr) = linds(ii) >*/
	    linds[*lctr + itscr] = linds[ii];
/*<           itscr = itscr + 1 >*/
	    ++itscr;
/*< 80      continue >*/
/* L80: */
	}
/*<       if (jml .eq. -1) jml = itscr >*/
	if (jml == -1) {
	    jml = itscr;
	}
/*<         lptrs(2,node) = itscr >*/
	lptrs[node * 3 + 2] = itscr;
/*<         nnz = nnz + itscr >*/
	*nnz += itscr;
/*<         opcount = opcount + dble(itscr * itscr) >*/
	*opcount += (doublereal) (itscr * itscr);
/*<       kk = nmiss * isupsize >*/
	kk = nmiss * isupsize;
/*<       xop = dble(kk * kk) >*/
	xop = (doublereal) (kk * kk);
/*<       if (nmiss .gt. 0 .and. xop*(opxtra/opcount) .lt. ssthresh) then >*/
	if (nmiss > 0 && xop * (*opxtra / *opcount) < *ssthresh) {
/*<         lspace = lspace + kk >*/
	    *lspace += kk;
/*<         jj = 0 >*/
	    jj = 0;
/*<         do i = isupbot, j+1, -1 >*/
	    i__1 = j + 1;
	    for (i__ = isupbot; i__ >= i__1; --i__) {
/*<           nnd = mystak(i) >*/
		nnd = mystak[i__];
/*<           lptrs(1,nnd) = lptrs(1,nnd) + jj >*/
		lptrs[nnd * 3 + 1] += jj;
/*<           lptrs(2,nnd) = lptrs(2,nnd) + nmiss >*/
		lptrs[nnd * 3 + 2] += nmiss;
/*<           jj = jj + nmiss >*/
		jj += nmiss;
/*<         end do >*/
	    }
/*<         do i = jml, itscr - 1 >*/
	    i__1 = itscr - 1;
	    for (i__ = jml; i__ <= i__1; ++i__) {
/*<           linds(jlctr+i) = linds(lctr+i) >*/
		linds[jlctr + i__] = linds[*lctr + i__];
/*<         end do >*/
	    }
/*<         opxtra = opxtra + xop >*/
	    *opxtra += xop;
/*<         lds = lds + nmiss >*/
	    lds += nmiss;
/*<         lctr = jlctr + itscr >*/
	    *lctr = jlctr + itscr;
/*<         if(lctr.gt.lsize) then >*/
	    if (*lctr > *lsize) {
/*<           info = 1 >*/
		*info = 1;
/*<           return >*/
		return 0;
/*<         end if >*/
	    }
/*<         nnzxtra = nnzxtra + kk >*/
	    *nnzxtra += kk;
/*<         nmiss = 0 >*/
	    nmiss = 0;
/*<       end if >*/
	}
/*<         if (nmiss .ne. 0) then >*/
	if (nmiss != 0) {
/*<           lptrs(1,node) = lspace >*/
	    lptrs[node * 3 + 1] = *lspace;
/*<           lptrs(3,node) = lctr >*/
	    lptrs[node * 3 + 3] = *lctr;
/*<           lds = itscr >*/
	    lds = itscr;
/*<           jlctr = lctr + 1 >*/
	    jlctr = *lctr + 1;
/*<           lctr = lctr + itscr >*/
	    *lctr += itscr;
/*<         if(lctr.gt.lsize) then >*/
	    if (*lctr > *lsize) {
/*<           info = 1 >*/
		*info = 1;
/*<           return >*/
		return 0;
/*<         end if >*/
	    }
/*<           tptrs(3,supbot) = scount >*/
	    tptrs[supbot * 3 + 3] = *scount;
/*<           tptrs(3,suptop) = scount >*/
	    tptrs[suptop * 3 + 3] = *scount;
/*<           sup(scount) = supbot >*/
	    sup[*scount] = supbot;
/*<           sup(scount+1) = isupsize >*/
	    sup[*scount + 1] = isupsize;
/*<           sup(scount+2) = node >*/
	    sup[*scount + 2] = node;
/*<           scount = scount + 4 >*/
	    *scount += 4;
/*<           supbot = node >*/
	    supbot = node;
/*<         isupbot = j >*/
	    isupbot = j;
/*<           isupsize = 0 >*/
	    isupsize = 0;
/*<         else >*/
	} else {
/*<           lptrs(1,node) = lspace + isupsize >*/
	    lptrs[node * 3 + 1] = *lspace + isupsize;
/*<           lptrs(3,node) = jlctr >*/
	    lptrs[node * 3 + 3] = jlctr;
/*<           jlctr = jlctr + 1 >*/
	    ++jlctr;
/*<         end if >*/
	}
/*<         lspace = lspace + lds >*/
	*lspace += lds;
/*< 30    continue >*/
/* L30: */
    }
/*<       isupsize = isupsize + 1 >*/
    ++isupsize;
/*<       tptrs(3,supbot) = scount >*/
    tptrs[supbot * 3 + 3] = *scount;
/*<       tptrs(3,node) = scount >*/
    tptrs[node * 3 + 3] = *scount;
/*<       sup(scount) = supbot >*/
    sup[*scount] = supbot;
/*<       sup(scount+1) = isupsize >*/
    sup[*scount + 1] = isupsize;
/*<       sup(scount+2) = -1 >*/
    sup[*scount + 2] = -1;
/*<       scount = scount + 4 >*/
    *scount += 4;
/*<       return >*/
    return 0;
/*<       end  >*/
} /* symbolic6_ */

/* ------------------------------------------------------------------------------- */
/*<       subroutine smerge(i1,l1,i2,l2,i3,l3) >*/
/* Subroutine */ int smerge_(integer *i1, integer *l1, integer *i2, integer *
	l2, integer *i3, integer *l3)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__, ii, jj;

/*<       integer i1,i2,i3,l1(*),l2(*),l3(0:*) >*/
/*<       ii = 1 >*/
    /* Parameter adjustments */
    --l2;
    --l1;

    /* Function Body */
    ii = 1;
/*<       jj = 1 >*/
    jj = 1;
/*<       i3 = 0 >*/
    *i3 = 0;
/*<       do while (ii .le. i1 .and. jj .le. i2) >*/
    while(ii <= *i1 && jj <= *i2) {
/*<       if (l1(ii) .eq. l2(jj)) then >*/
	if (l1[ii] == l2[jj]) {
/*<         l3(i3) = l1(ii) >*/
	    l3[*i3] = l1[ii];
/*<         ii = ii + 1 >*/
	    ++ii;
/*<         jj = jj + 1 >*/
	    ++jj;
/*<         i3 = i3 + 1 >*/
	    ++(*i3);
/*<       else if (l1(ii) .gt. l2(jj)) then >*/
	} else if (l1[ii] > l2[jj]) {
/*<         l3(i3) = l2(jj) >*/
	    l3[*i3] = l2[jj];
/*<         jj = jj + 1 >*/
	    ++jj;
/*<         i3 = i3 + 1 >*/
	    ++(*i3);
/*<       else >*/
	} else {
/*<         l3(i3) = l1(ii) >*/
	    l3[*i3] = l1[ii];
/*<         ii = ii + 1 >*/
	    ++ii;
/*<         i3 = i3 + 1 >*/
	    ++(*i3);
/*<       end if >*/
	}
/*<       end do >*/
    }
/*<       do i = ii, i1 >*/
    i__1 = *i1;
    for (i__ = ii; i__ <= i__1; ++i__) {
/*<       l3(i3) = l1(i) >*/
	l3[*i3] = l1[i__];
/*<       i3 = i3 + 1 >*/
	++(*i3);
/*<       end do >*/
    }
/*<       do i = jj, i2 >*/
    i__1 = *i2;
    for (i__ = jj; i__ <= i__1; ++i__) {
/*<       l3(i3) = l2(i) >*/
	l3[*i3] = l2[i__];
/*<       i3 = i3 + 1 >*/
	++(*i3);
/*<       end do >*/
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* smerge_ */

/* ------------------------------------------------------------------------------- */
/*<    >*/
/* Subroutine */ int fsolvefsolve_(integer *n, doublereal *lvals, integer *
	linds, integer *lptrs, integer *tinds, integer *tptrs, integer *sup, 
	doublereal *rhs, integer *nrhs, integer *root, integer *lc, integer *
	iptrs, doublereal *w)
{
    /* System generated locals */
    integer rhs_dim1, rhs_offset;

    /* Local variables */
    integer i__, ldr;
    extern /* Subroutine */ int rfsolve_(integer *, integer *, doublereal *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *);

/*<       integer N,linds(*),lptrs(3,0:*),tinds(*),tptrs(3,0:*) >*/
/*<       integer sup(*),nrhs,root,lc(*),iptrs(2,0:*) >*/
/*<       double precision lvals(*),rhs(0:N-1,*),w(*) >*/
/*<    >*/
    /* Parameter adjustments */
    rhs_dim1 = *n - 1 - 0 + 1;
    rhs_offset = 0 + rhs_dim1;
    rhs -= rhs_offset;
    --lvals;
    --linds;
    --lptrs;
    --tinds;
    --tptrs;
    --sup;
    --lc;
    --iptrs;
    --w;

    /* Function Body */
    rfsolve_(n, root, &lvals[1], &linds[1], &lptrs[1], &tinds[1], &tptrs[1], &
	    sup[1], &rhs[rhs_offset], nrhs, &iptrs[1], &lc[1], &i__, &ldr, &w[
	    1]);
/*<       return >*/
    return 0;
/*<       end >*/
} /* fsolvefsolve_ */

/* ------------------------------------------------------------------------------ */
/*<    >*/
/* Subroutine */ int bsolve_(integer *n, doublereal *lvals, integer *linds, 
	integer *lptrs, integer *tinds, integer *tptrs, integer *sup, 
	doublereal *rhs, integer *nrhs, integer *root, integer *lc, integer *
	iptrs, doublereal *w)
{
    /* System generated locals */
    integer rhs_dim1, rhs_offset;

    /* Local variables */
    extern /* Subroutine */ int rbsolve_(integer *, integer *, doublereal *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *);

/*<       integer N,linds(*),lptrs(3,0:*),tinds(*),tptrs(3,0:*) >*/
/*<       integer sup(*),nrhs,root,lc(*),iptrs(2,0:*) >*/
/*<       double precision lvals(*),rhs(0:N-1,*),w(*) >*/
/*<    >*/
    /* Parameter adjustments */
    rhs_dim1 = *n - 1 - 0 + 1;
    rhs_offset = 0 + rhs_dim1;
    rhs -= rhs_offset;
    --lvals;
    --linds;
    --lptrs;
    --tinds;
    --tptrs;
    --sup;
    --lc;
    --iptrs;
    --w;

    /* Function Body */
    rbsolve_(n, root, &lvals[1], &linds[1], &lptrs[1], &tinds[1], &tptrs[1], &
	    sup[1], &rhs[rhs_offset], nrhs, &iptrs[1], &lc[1], &w[1], &c__0, &
	    w[1]);
/*<       return >*/
    return 0;
/*<       end >*/
} /* bsolve_ */

/* ------------------------------------------------------------------------------ */
/*     recursive */
/*    + */
/*<    >*/
/* Subroutine */ int rbsolve_(integer *n, integer *root, doublereal *lvals, 
	integer *linds, integer *lptrs, integer *tinds, integer *tptrs, 
	integer *sup, doublereal *rhs, integer *nrhs, integer *iptrs, integer 
	*lc, doublereal *w, integer *ldw, doublereal *wp)
{
    /* System generated locals */
    integer rhs_dim1, rhs_offset, i__1;

    /* Local variables */
    extern /* Subroutine */ int compress_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, integer *, integer *);
    integer i__;
    integer kid, ldj, ldr, jloc, rank, jbot;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen),
	     dgemv_(char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, ftnlen), dtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), dtrsv_(
	    char *, char *, char *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen);

/* why did f2c create this? */
/*    real riptrs; */

    extern /* Subroutine */ int putrhs_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *), bgetrhs_(
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, doublereal *);

/*<       integer N,root,linds(*),lptrs(3,0:*),tinds(*),tptrs(3,0:*) >*/
/*<       integer sup(*),iptrs(2,0:*),lc(*),ldr,ldw >*/
/*<       double precision lvals(*),w(0:*),rhs(0:N-1,*),wp(0:*) >*/
/*<       integer rank >*/
/*<       i = tptrs(3,root) >*/
    /* Parameter adjustments */
    rhs_dim1 = *n - 1 - 0 + 1;
    rhs_offset = 0 + rhs_dim1;
    rhs -= rhs_offset;
    --lvals;
    --linds;
    --lptrs;
    --tinds;
    --tptrs;
    --sup;
    --iptrs;
    --lc;

    /* Function Body */
    i__ = tptrs[*root * 3 + 3];
/*<       jbot = sup(i) >*/
    jbot = sup[i__];
/*<       rank = sup(i+1) >*/
    rank = sup[i__ + 1];
/*<       ldr = lptrs(2,jbot) >*/
    ldr = lptrs[jbot * 3 + 2];
/*<       ldj = ldr - rank >*/
    ldj = ldr - rank;
/*<       call compress(root,wp(rank),ldr,w,ldw,iptrs,lc,nrhs) >*/
    compress_(root, &wp[rank], &ldr, w, ldw, &iptrs[1], &lc[1], nrhs);
/*<       call bgetrhs(N,linds(lptrs(3,jbot)),rank,ldr,nrhs,wp,rhs) >*/
    bgetrhs_(n, &linds[lptrs[jbot * 3 + 3]], &rank, &ldr, nrhs, wp, &rhs[
	    rhs_offset]);
/*<       if (nrhs .gt. 1) then >*/
    if (*nrhs > 1) {
/*<    >*/
	dgemm_("T", "N", &rank, nrhs, &ldj, &c_b57, &lvals[lptrs[jbot * 3 + 1]
		 + rank], &ldr, &wp[rank], &ldr, &c_b58, wp, &ldr, (ftnlen)1, 
		(ftnlen)1);
/*<    >*/
	dtrsm_("L", "L", "T", "N", &rank, nrhs, &c_b58, &lvals[lptrs[jbot * 3 
		+ 1]], &ldr, wp, &ldr, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);
/*<       else >*/
    } else {
/*<    >*/
	i__1 = ldr - rank;
	dgemv_("T", &i__1, &rank, &c_b57, &lvals[lptrs[jbot * 3 + 1] + rank], 
		&ldr, &wp[rank], &c__1, &c_b58, wp, &c__1, (ftnlen)1);
/*<         call dtrsv('L','T','N',rank,lvals(lptrs(1,jbot)),ldr,wp,1) >*/
	dtrsv_("L", "T", "N", &rank, &lvals[lptrs[jbot * 3 + 1]], &ldr, wp, &
		c__1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
/*<       end if >*/
    }
/*<       call putrhs(N,linds(lptrs(3,jbot)),rank,ldr,nrhs,wp,rhs) >*/
    putrhs_(n, &linds[lptrs[jbot * 3 + 3]], &rank, &ldr, nrhs, wp, &rhs[
	    rhs_offset]);
/*<       jloc = ldr * nrhs >*/
    jloc = ldr * *nrhs;
/*<       do i = tptrs(1, jbot), tptrs(1, jbot) + tptrs(2,jbot) - 1 >*/
    i__1 = tptrs[jbot * 3 + 1] + tptrs[jbot * 3 + 2] - 1;
    for (i__ = tptrs[jbot * 3 + 1]; i__ <= i__1; ++i__) {
/*<         kid = tinds(i) >*/
	kid = tinds[i__];
/*<    >*/

/*
	rbsolve_(n, &kid, &lvals[1], &linds[1], &lptrs[1], &tinds[
		1], &tptrs[1], &sup[1], &riptrs, &lc[1], wp, &ldr, &wp[jloc]);
*/
	rbsolve_(n, &kid, &lvals[1], &linds[1], &lptrs[1], &tinds[
		1], &tptrs[1], &sup[1], &rhs[rhs_offset], nrhs, &iptrs[1], &lc[1], wp, &ldr, &wp[jloc]);


/*<       end do >*/
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* rbsolve_ */

/* ------------------------------------------------------------------------------ */
/*<       subroutine BGETRHS(N,linds,rank,ldr,nrhs,w,rhs) >*/
/* Subroutine */ int bgetrhs_(integer *n, integer *linds, integer *rank, 
	integer *ldr, integer *nrhs, doublereal *w, doublereal *rhs)
{
    /* System generated locals */
    integer w_dim1, w_offset, rhs_dim1, rhs_offset, i__1, i__2;

    /* Local variables */
    integer i__, j, k;

/*<       integer N,linds(*),rank,ldr,nrhs >*/
/*<       double precision w(ldr,*),rhs(0:N-1,*) >*/
/*<       do 10 i = 1, rank >*/
    /* Parameter adjustments */
    rhs_dim1 = *n - 1 - 0 + 1;
    rhs_offset = 0 + rhs_dim1;
    rhs -= rhs_offset;
    --linds;
    w_dim1 = *ldr;
    w_offset = 1 + w_dim1;
    w -= w_offset;

    /* Function Body */
    i__1 = *rank;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         j = linds(i) >*/
	j = linds[i__];
/*<         do 20 k = 1, nrhs >*/
	i__2 = *nrhs;
	for (k = 1; k <= i__2; ++k) {
/*<           w(i,k) = rhs(j,k) >*/
	    w[i__ + k * w_dim1] = rhs[j + k * rhs_dim1];
/*< 20      continue >*/
/* L20: */
	}
/*< 10    continue >*/
/* L10: */
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* bgetrhs_ */

/* ------------------------------------------------------------------------------ */
/*<       subroutine BGETRHS1(N,linds,rank,w,rhs) >*/
/* Subroutine */ int bgetrhs1_(integer *n, integer *linds, integer *rank, 
	doublereal *w, doublereal *rhs)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__;

/*<       integer N,linds(*),rank >*/
/*<       double precision w(*),rhs(0:*) >*/
/*<       do 10 i = 1, rank >*/
    /* Parameter adjustments */
    --w;
    --linds;

    /* Function Body */
    i__1 = *rank;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         w(i) = rhs(linds(i)) >*/
	w[i__] = rhs[linds[i__]];
/*< 10    continue >*/
/* L10: */
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* bgetrhs1_ */

/* ------------------------------------------------------------------------------ */
/*<       subroutine COMPRESS(node,dest,ldd,source,lds,iptrs,lc,nrhs) >*/
/* Subroutine */ int compress_(integer *node, doublereal *dest, integer *ldd, 
	doublereal *source, integer *lds, integer *iptrs, integer *lc, 
	integer *nrhs)
{
    /* System generated locals */
    integer dest_dim1, dest_offset, source_dim1, source_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    integer i__, j, k, m, nr;

/*<       integer node,ldd,lds,nrhs >*/
/*<       integer iptrs(2,0:*),lc(*) >*/
/*<       double precision dest(ldd,*),source(lds,*) >*/
/*<       m = 1 >*/
    /* Parameter adjustments */
    dest_dim1 = *ldd;
    dest_offset = 1 + dest_dim1;
    dest -= dest_offset;
    source_dim1 = *lds;
    source_offset = 1 + source_dim1;
    source -= source_offset;
    --iptrs;
    --lc;

    /* Function Body */
    m = 1;
/*<       j = 1 >*/
    j = 1;
/*<       do i = iptrs(1,node), iptrs(1,node) + iptrs(2,node) - 1, 2 >*/
    i__1 = iptrs[(*node << 1) + 1] + iptrs[(*node << 1) + 2] - 1;
    for (i__ = iptrs[(*node << 1) + 1]; i__ <= i__1; i__ += 2) {
/*<         do k = j, lc(i) + j - 1 >*/
	i__2 = lc[i__] + j - 1;
	for (k = j; k <= i__2; ++k) {
/*<           do nr = 1, nrhs >*/
	    i__3 = *nrhs;
	    for (nr = 1; nr <= i__3; ++nr) {
/*<             dest(m,nr) = source(k,nr) >*/
		dest[m + nr * dest_dim1] = source[k + nr * source_dim1];
/*<           end do >*/
	    }
/*<           m = m + 1 >*/
	    ++m;
/*<         end do >*/
	}
/*<         j = k + lc(i+1) >*/
	j = k + lc[i__ + 1];
/*<       end do >*/
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* compress_ */

/* ------------------------------------------------------------------------------ */
/*     recursive */
/*    + */
/*<    >*/
/* Subroutine */ int rfsolve_(integer *n, integer *root, doublereal *lvals, 
	integer *linds, integer *lptrs, integer *tinds, integer *tptrs, 
	integer *sup, doublereal *rhs, integer *nrhs, integer *iptrs, integer 
	*lc, integer *rank, integer *ldr, doublereal *w)
{
    /* System generated locals */
    integer rhs_dim1, rhs_offset, i__1;

    /* Local variables */
    extern /* Subroutine */ int mergerhs_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, integer *, integer *);
    integer i__;
    integer kid, ldj, jloc, jbot;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen),
	     dgemv_(char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, ftnlen);
    integer jrank;
    extern /* Subroutine */ int dtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), dtrsv_(
	    char *, char *, char *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen), getrhs_(integer 
	    *, integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *), putrhs_(integer *, integer *, integer *, integer *,
	     integer *, doublereal *, doublereal *);

/*<       integer N,root,linds(*),lptrs(3,0:*),tinds(N),tptrs(3,0:*) >*/
/*<       integer sup(*),iptrs(2,0:*),lc(*),rank,ldr >*/
/*<       double precision lvals(*),w(0:*),rhs(0:N-1,*) >*/
/*<       i = tptrs(3,root) >*/
    /* Parameter adjustments */
    rhs_dim1 = *n - 1 - 0 + 1;
    rhs_offset = 0 + rhs_dim1;
    rhs -= rhs_offset;
    --tinds;
    --lvals;
    --linds;
    --lptrs;
    --tptrs;
    --sup;
    --iptrs;
    --lc;

    /* Function Body */
    i__ = tptrs[*root * 3 + 3];
/*<       jbot = sup(i) >*/
    jbot = sup[i__];
/*<       rank = sup(i+1) >*/
    *rank = sup[i__ + 1];
/*<       ldr = lptrs(2,jbot) >*/
    *ldr = lptrs[jbot * 3 + 2];
/*<       jloc = ldr * nrhs >*/
    jloc = *ldr * *nrhs;
/*<       call getrhs(N,linds(lptrs(3,jbot)),rank,ldr,nrhs,w,rhs) >*/
    getrhs_(n, &linds[lptrs[jbot * 3 + 3]], rank, ldr, nrhs, w, &rhs[
	    rhs_offset]);
/*<       do i = tptrs(1,jbot), tptrs(1,jbot) + tptrs(2,jbot) - 1 >*/
    i__1 = tptrs[jbot * 3 + 1] + tptrs[jbot * 3 + 2] - 1;
    for (i__ = tptrs[jbot * 3 + 1]; i__ <= i__1; ++i__) {
/*<         kid = tinds(i) >*/
	kid = tinds[i__];
/*<    >*/
	rfsolve_(n, &kid, &lvals[1], &linds[1], &lptrs[1], &tinds[
		1], &tptrs[1], &sup[1], &rhs[rhs_offset], nrhs, &iptrs[1], &
		lc[1], &jrank, &ldj, &w[jloc]);
/*<         call mergerhs(ldr,w,ldj,w(jloc+jrank),nrhs,kid,iptrs,lc) >*/
	mergerhs_(ldr, w, &ldj, &w[jloc + jrank], nrhs, &kid, &iptrs[1], &lc[
		1]);
/*<       end do >*/
    }
/*<       if (nrhs .eq. 1) then >*/
    if (*nrhs == 1) {
/*<       call dtrsv('L','N','N',rank,lvals(lptrs(1,jbot)),ldr,w,1) >*/
	dtrsv_("L", "N", "N", rank, &lvals[lptrs[jbot * 3 + 1]], ldr, w, &
		c__1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
/*<    >*/
	i__1 = *ldr - *rank;
	dgemv_("N", &i__1, rank, &c_b57, &lvals[lptrs[jbot * 3 + 1] + *rank], 
		ldr, w, &c__1, &c_b58, &w[*rank], &c__1, (ftnlen)1);
/*<       else >*/
    } else {
/*<    >*/
	dtrsm_("L", "L", "N", "N", rank, nrhs, &c_b58, &lvals[lptrs[jbot * 3 
		+ 1]], ldr, w, ldr, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)
		1);
/*<    >*/
	i__1 = *ldr - *rank;
	dgemm_("N", "N", &i__1, nrhs, rank, &c_b57, &lvals[lptrs[jbot * 3 + 1]
		 + *rank], ldr, w, ldr, &c_b58, &w[*rank], ldr, (ftnlen)1, (
		ftnlen)1);
/*<       end if >*/
    }
/*<       call putrhs(N,linds(lptrs(3,jbot)),rank,ldr,nrhs,w,rhs) >*/
    putrhs_(n, &linds[lptrs[jbot * 3 + 3]], rank, ldr, nrhs, w, &rhs[
	    rhs_offset]);
/*<       return >*/
    return 0;
/*<       end >*/
} /* rfsolve_ */

/* ------------------------------------------------------------------------------ */
/*<       subroutine MERGERHS(ldr,w1,ldj,w2,nrhs,kid,iptrs,lc) >*/
/* Subroutine */ int mergerhs_(integer *ldr, doublereal *w1, integer *ldj, 
	doublereal *w2, integer *nrhs, integer *kid, integer *iptrs, integer *
	lc)
{
    /* System generated locals */
    integer w1_dim1, w1_offset, w2_dim1, w2_offset, i__1, i__2, i__3;

    /* Local variables */
    integer i__, j, k, l, m, nr;

/*<       integer ldr,ldj,nrhs,kid,iptrs(2,0:*),lc(*) >*/
/*<       double precision w1(ldr,*),w2(ldj,*) >*/
/*<       j = 1 >*/
    /* Parameter adjustments */
    w1_dim1 = *ldr;
    w1_offset = 1 + w1_dim1;
    w1 -= w1_offset;
    w2_dim1 = *ldj;
    w2_offset = 1 + w2_dim1;
    w2 -= w2_offset;
    --iptrs;
    --lc;

    /* Function Body */
    j = 1;
/*<       m = 1 >*/
    m = 1;
/*<       l = iptrs(1,kid) >*/
    l = iptrs[(*kid << 1) + 1];
/*<       do i = l, l+iptrs(2,kid)-1, 2 >*/
    i__1 = l + iptrs[(*kid << 1) + 2] - 1;
    for (i__ = l; i__ <= i__1; i__ += 2) {
/*<         do k = j, lc(i) + j - 1 >*/
	i__2 = lc[i__] + j - 1;
	for (k = j; k <= i__2; ++k) {
/*<           do nr = 1, nrhs >*/
	    i__3 = *nrhs;
	    for (nr = 1; nr <= i__3; ++nr) {
/*<             w1(k, nr) = w1(k,nr) + w2(m,nr) >*/
		w1[k + nr * w1_dim1] += w2[m + nr * w2_dim1];
/*<           end do >*/
	    }
/*<           m = m + 1 >*/
	    ++m;
/*<         end do >*/
	}
/*<         j = k + lc(i+1) >*/
	j = k + lc[i__ + 1];
/*<       end do >*/
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* mergerhs_ */

/* ------------------------------------------------------------------------------ */
/*<       subroutine GETRHS(N,linds,rank,ldr,nrhs,w,rhs) >*/
/* Subroutine */ int getrhs_(integer *n, integer *linds, integer *rank, 
	integer *ldr, integer *nrhs, doublereal *w, doublereal *rhs)
{
    /* System generated locals */
    integer w_dim1, w_offset, rhs_dim1, rhs_offset, i__1, i__2;

    /* Local variables */
    integer i__, j, k;

/*<       integer N,linds(*),rank,ldr,nrhs >*/
/*<       double precision  w(ldr,*),rhs(0:N-1,*) >*/
/*<       do i = 1, rank >*/
    /* Parameter adjustments */
    rhs_dim1 = *n - 1 - 0 + 1;
    rhs_offset = 0 + rhs_dim1;
    rhs -= rhs_offset;
    --linds;
    w_dim1 = *ldr;
    w_offset = 1 + w_dim1;
    w -= w_offset;

    /* Function Body */
    i__1 = *rank;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         j = linds(i) >*/
	j = linds[i__];
/*<         do k = 1, nrhs >*/
	i__2 = *nrhs;
	for (k = 1; k <= i__2; ++k) {
/*<           w(i,k) = rhs(j,k) >*/
	    w[i__ + k * w_dim1] = rhs[j + k * rhs_dim1];
/*<         end do >*/
	}
/*<       end do >*/
    }
/*<       do j = i, ldr >*/
    i__1 = *ldr;
    for (j = i__; j <= i__1; ++j) {
/*<         do k = 1, nrhs >*/
	i__2 = *nrhs;
	for (k = 1; k <= i__2; ++k) {
/*<           w(j,k) = 0.d0 >*/
	    w[j + k * w_dim1] = 0.;
/*<         end do >*/
	}
/*<       end do >*/
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* getrhs_ */

/* ------------------------------------------------------------------------------ */
/*<       subroutine PUTRHS(N,linds,rank,ldr,nrhs,w,rhs) >*/
/* Subroutine */ int putrhs_(integer *n, integer *linds, integer *rank, 
	integer *ldr, integer *nrhs, doublereal *w, doublereal *rhs)
{
    /* System generated locals */
    integer w_dim1, w_offset, rhs_dim1, rhs_offset, i__1, i__2;

    /* Local variables */
    integer i__, j, k;

/*<       integer N,linds(*),rank,ldr,nrhs >*/
/*<       double precision  w(ldr,*),rhs(0:N-1,*) >*/
/*<       do i = 1, rank >*/
    /* Parameter adjustments */
    rhs_dim1 = *n - 1 - 0 + 1;
    rhs_offset = 0 + rhs_dim1;
    rhs -= rhs_offset;
    --linds;
    w_dim1 = *ldr;
    w_offset = 1 + w_dim1;
    w -= w_offset;

    /* Function Body */
    i__1 = *rank;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         j = linds(i) >*/
	j = linds[i__];
/*<         do k = 1, nrhs >*/
	i__2 = *nrhs;
	for (k = 1; k <= i__2; ++k) {
/*<           rhs(j,k) = w(i,k) >*/
	    rhs[j + k * rhs_dim1] = w[i__ + k * w_dim1];
/*<         end do >*/
	}
/*<       end do >*/
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* putrhs_ */

/* ------------------------------------------------------------------------------ */
/*<    >*/
/* Subroutine */ int dchol2_(doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, integer *ldim, integer *lrank, integer *lup, integer *
	linds, doublereal *dfopts, integer *ifopts, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;

    /* Local variables */
    integer i__, j, k, k1;
    doublereal a11, a21, a22, a31, a41, a12, a32, a42, r21, t11, t21, t12, 
	    t22, x11, x21, x12, x22, x31, x41, x32, x42, rd1, rd2;
    integer lr1;
    doublereal ddmin, ddmax, dsmall;

/*<       implicit double precision (a-h,o-z) >*/
/*<       double precision a(0:lda-1,0:*),dfopts(7) >*/
/*<       double precision b(0:ldb-1,0:*) >*/
/*<       integer lda,ldim,lrank,lup,linds(0:*),ifopts(5),ldb >*/
/*<       parameter (one = 1.d0,zero = 0.d0) >*/
/*<       dsmall = dfopts(1) >*/
    /* Parameter adjustments */
    a_dim1 = *lda - 1 - 0 + 1;
    a_offset = 0 + a_dim1 * 0;
    a -= a_offset;
    b_dim1 = *ldb - 1 - 0 + 1;
    b_offset = 0 + b_dim1 * 0;
    b -= b_offset;
    --dfopts;
    --ifopts;

    /* Function Body */
    dsmall = dfopts[1];
/*<       ddmax = dfopts(5) >*/
    ddmax = dfopts[5];
/*<       ddmin = dfopts(6) >*/
    ddmin = dfopts[6];
/*<       lr1 = mod(lrank,2) >*/
    lr1 = *lrank % 2;
/*<       if (lr1 .eq. 1) then >*/
    if (lr1 == 1) {
/*<       a11 = a(0,0) >*/
	a11 = a[0];
/*<       if (a11 .le. 0.d0) go to 99 >*/
	if (a11 <= 0.) {
	    goto L99;
	}
/*<       a11 = dsqrt(a11) >*/
	a11 = sqrt(a11);
/*< 2111    if (a11 .gt. ddmax) ddmax = a11 >*/
/* L2111: */
	if (a11 > ddmax) {
	    ddmax = a11;
	}
/*<       if (a11 .lt. ddmin) ddmin = a11 >*/
	if (a11 < ddmin) {
	    ddmin = a11;
	}
/*< 3111    rd1 = one/a11 >*/
/* L3111: */
	rd1 = 1. / a11;
/*<         a(0,0) = a11 >*/
	a[0] = a11;
/*<       do i = 1, ldim - 4, 4 >*/
	i__1 = *ldim - 4;
	for (i__ = 1; i__ <= i__1; i__ += 4) {
/*<         a(i,0) = a(i,0) * rd1 >*/
	    a[i__] *= rd1;
/*<         a(i+1,0) = a(i+1,0) * rd1 >*/
	    a[i__ + 1] *= rd1;
/*<         a(i+2,0) = a(i+2,0) * rd1 >*/
	    a[i__ + 2] *= rd1;
/*<         a(i+3,0) = a(i+3,0) * rd1 >*/
	    a[i__ + 3] *= rd1;
/*<       end do >*/
	}
/*<       do j = i, ldim - 1 >*/
	i__1 = *ldim - 1;
	for (j = i__; j <= i__1; ++j) {
/*<         a(j,0) = a(j,0) * rd1 >*/
	    a[j] *= rd1;
/*<       end do >*/
	}
/*<       end if >*/
    }
/*<       do i = lr1, lrank - 1, 2 >*/
    i__1 = *lrank - 1;
    for (i__ = lr1; i__ <= i__1; i__ += 2) {
/*<       a11 = a(i,i) >*/
	a11 = a[i__ + i__ * a_dim1];
/*<       a21 = a(i+1,i) >*/
	a21 = a[i__ + 1 + i__ * a_dim1];
/*<       a22 = a(i+1,i+1) >*/
	a22 = a[i__ + 1 + (i__ + 1) * a_dim1];
/*<       do j = 0, i - 2, 2 >*/
	i__2 = i__ - 2;
	for (j = 0; j <= i__2; j += 2) {
/*<         t11 = a(i,j) >*/
	    t11 = a[i__ + j * a_dim1];
/*<         t21 = a(i+1,j) >*/
	    t21 = a[i__ + 1 + j * a_dim1];
/*<         t12 = a(i,j+1) >*/
	    t12 = a[i__ + (j + 1) * a_dim1];
/*<         t22 = a(i+1,j+1) >*/
	    t22 = a[i__ + 1 + (j + 1) * a_dim1];
/*<         a11 = a11 - t11 * t11 >*/
	    a11 -= t11 * t11;
/*<         a21 = a21 - t11 * t21 >*/
	    a21 -= t11 * t21;
/*<         a22 = a22 - t21 * t21 >*/
	    a22 -= t21 * t21;
/*<         a11 = a11 - t12 * t12 >*/
	    a11 -= t12 * t12;
/*<         a21 = a21 - t12 * t22 >*/
	    a21 -= t12 * t22;
/*<         a22 = a22 - t22 * t22 >*/
	    a22 -= t22 * t22;
/*<       end do >*/
	}
/*<       if (j .eq. i-1) then >*/
	if (j == i__ - 1) {
/*<         t11 = a(i,j) >*/
	    t11 = a[i__ + j * a_dim1];
/*<         t21 = a(i+1,j) >*/
	    t21 = a[i__ + 1 + j * a_dim1];
/*<         a11 = a11 - t11 * t11 >*/
	    a11 -= t11 * t11;
/*<         a21 = a21 - t11 * t21 >*/
	    a21 -= t11 * t21;
/*<         a22 = a22 - t21 * t21 >*/
	    a22 -= t21 * t21;
/*<       end if >*/
	}
/*<         if (a11 .le. 0.d0) go to 99 >*/
	if (a11 <= 0.) {
	    goto L99;
	}
/*<       a11 = dsqrt(a11) >*/
	a11 = sqrt(a11);
/*< 2112    if (a11 .gt. ddmax) ddmax = a11 >*/
/* L2112: */
	if (a11 > ddmax) {
	    ddmax = a11;
	}
/*<         if (a11 .lt. ddmin) ddmin = a11 >*/
	if (a11 < ddmin) {
	    ddmin = a11;
	}
/*< 3112    a(i,i) = a11 >*/
/* L3112: */
	a[i__ + i__ * a_dim1] = a11;
/*<       rd1 = one/a11 >*/
	rd1 = 1. / a11;
/*<       r21 = a21 * rd1 >*/
	r21 = a21 * rd1;
/*<       a(i+1,i) = r21 >*/
	a[i__ + 1 + i__ * a_dim1] = r21;
/*<       a22 = a22 - r21 * r21 >*/
	a22 -= r21 * r21;
/*<       if (a22 .le. 0.d0) go to 99 >*/
	if (a22 <= 0.) {
	    goto L99;
	}
/*<       a22 = dsqrt(a22) >*/
	a22 = sqrt(a22);
/*< 2113    if (a22 .gt. ddmax) ddmax = a22 >*/
/* L2113: */
	if (a22 > ddmax) {
	    ddmax = a22;
	}
/*<         if (a22 .lt. ddmin) ddmin = a22 >*/
	if (a22 < ddmin) {
	    ddmin = a22;
	}
/*< 3113    a(i+1,i+1) = a22 >*/
/* L3113: */
	a[i__ + 1 + (i__ + 1) * a_dim1] = a22;
/*<       rd2 = one/a22 >*/
	rd2 = 1. / a22;
/*<       do k = i+2, ldim - 4, 4 >*/
	i__2 = *ldim - 4;
	for (k = i__ + 2; k <= i__2; k += 4) {
/*<         a11 = a(k,i) >*/
	    a11 = a[k + i__ * a_dim1];
/*<         a21 = a(k+1,i) >*/
	    a21 = a[k + 1 + i__ * a_dim1];
/*<         a31 = a(k+2,i) >*/
	    a31 = a[k + 2 + i__ * a_dim1];
/*<         a41 = a(k+3,i) >*/
	    a41 = a[k + 3 + i__ * a_dim1];
/*<         a12 = a(k,i+1) >*/
	    a12 = a[k + (i__ + 1) * a_dim1];
/*<         a22 = a(k+1,i+1) >*/
	    a22 = a[k + 1 + (i__ + 1) * a_dim1];
/*<         a32 = a(k+2,i+1) >*/
	    a32 = a[k + 2 + (i__ + 1) * a_dim1];
/*<         a42 = a(k+3,i+1) >*/
	    a42 = a[k + 3 + (i__ + 1) * a_dim1];
/*<         do j = 0, i - 2, 2 >*/
	    i__3 = i__ - 2;
	    for (j = 0; j <= i__3; j += 2) {
/*<           t11 = a(i,j) >*/
		t11 = a[i__ + j * a_dim1];
/*<           t21 = a(i+1,j) >*/
		t21 = a[i__ + 1 + j * a_dim1];
/*<           t12 = a(i,j+1) >*/
		t12 = a[i__ + (j + 1) * a_dim1];
/*<           t22 = a(i+1,j+1) >*/
		t22 = a[i__ + 1 + (j + 1) * a_dim1];
/*<           x11 = a(k,j) >*/
		x11 = a[k + j * a_dim1];
/*<           x21 = a(k+1,j) >*/
		x21 = a[k + 1 + j * a_dim1];
/*<           x12 = a(k,j+1) >*/
		x12 = a[k + (j + 1) * a_dim1];
/*<           x22 = a(k+1,j+1) >*/
		x22 = a[k + 1 + (j + 1) * a_dim1];
/*<           a11 = a11 - x11 * t11 >*/
		a11 -= x11 * t11;
/*<           a21 = a21 - x21 * t11 >*/
		a21 -= x21 * t11;
/*<           a12 = a12 - x11 * t21 >*/
		a12 -= x11 * t21;
/*<           a22 = a22 - x21 * t21 >*/
		a22 -= x21 * t21;
/*<           a11 = a11 - x12 * t12 >*/
		a11 -= x12 * t12;
/*<           a21 = a21 - x22 * t12 >*/
		a21 -= x22 * t12;
/*<           a12 = a12 - x12 * t22 >*/
		a12 -= x12 * t22;
/*<           a22 = a22 - x22 * t22 >*/
		a22 -= x22 * t22;
/*<           x31 = a(k+2,j) >*/
		x31 = a[k + 2 + j * a_dim1];
/*<           x41 = a(k+3,j) >*/
		x41 = a[k + 3 + j * a_dim1];
/*<           x32 = a(k+2,j+1) >*/
		x32 = a[k + 2 + (j + 1) * a_dim1];
/*<           x42 = a(k+3,j+1) >*/
		x42 = a[k + 3 + (j + 1) * a_dim1];
/*<           a31 = a31 - x31 * t11 >*/
		a31 -= x31 * t11;
/*<           a41 = a41 - x41 * t11 >*/
		a41 -= x41 * t11;
/*<           a32 = a32 - x31 * t21 >*/
		a32 -= x31 * t21;
/*<           a42 = a42 - x41 * t21 >*/
		a42 -= x41 * t21;
/*<           a31 = a31 - x32 * t12 >*/
		a31 -= x32 * t12;
/*<           a41 = a41 - x42 * t12 >*/
		a41 -= x42 * t12;
/*<           a32 = a32 - x32 * t22 >*/
		a32 -= x32 * t22;
/*<           a42 = a42 - x42 * t22 >*/
		a42 -= x42 * t22;
/*<         end do >*/
	    }
/*<         if (j .eq. i-1) then >*/
	    if (j == i__ - 1) {
/*<           t11 = a(i,j) >*/
		t11 = a[i__ + j * a_dim1];
/*<           t21 = a(i+1,j) >*/
		t21 = a[i__ + 1 + j * a_dim1];
/*<           x11 = a(k,j) >*/
		x11 = a[k + j * a_dim1];
/*<           x21 = a(k+1,j)  >*/
		x21 = a[k + 1 + j * a_dim1];
/*<           x31 = a(k+2,j)  >*/
		x31 = a[k + 2 + j * a_dim1];
/*<           x41 = a(k+3,j)  >*/
		x41 = a[k + 3 + j * a_dim1];
/*<           a11 = a11 - x11 * t11 >*/
		a11 -= x11 * t11;
/*<           a21 = a21 - x21 * t11 >*/
		a21 -= x21 * t11;
/*<           a31 = a31 - x31 * t11 >*/
		a31 -= x31 * t11;
/*<           a41 = a41 - x41 * t11 >*/
		a41 -= x41 * t11;
/*<           a12 = a12 - x11 * t21 >*/
		a12 -= x11 * t21;
/*<           a22 = a22 - x21 * t21 >*/
		a22 -= x21 * t21;
/*<           a32 = a32 - x31 * t21 >*/
		a32 -= x31 * t21;
/*<           a42 = a42 - x41 * t21 >*/
		a42 -= x41 * t21;
/*<         end if >*/
	    }
/*<         a11 = a11 * rd1 >*/
	    a11 *= rd1;
/*<         a21 = a21 * rd1 >*/
	    a21 *= rd1;
/*<         a31 = a31 * rd1 >*/
	    a31 *= rd1;
/*<         a41 = a41 * rd1 >*/
	    a41 *= rd1;
/*<         a(k,i) = a11  >*/
	    a[k + i__ * a_dim1] = a11;
/*<         a(k+1,i) = a21  >*/
	    a[k + 1 + i__ * a_dim1] = a21;
/*<         a(k+2,i) = a31 >*/
	    a[k + 2 + i__ * a_dim1] = a31;
/*<         a(k+3,i) = a41 >*/
	    a[k + 3 + i__ * a_dim1] = a41;
/*<         a12 = a12 - a11 * r21 >*/
	    a12 -= a11 * r21;
/*<         a22 = a22 - a21 * r21 >*/
	    a22 -= a21 * r21;
/*<         a32 = a32 - a31 * r21 >*/
	    a32 -= a31 * r21;
/*<         a42 = a42 - a41 * r21 >*/
	    a42 -= a41 * r21;
/*<         a(k,i+1) = a12 * rd2 >*/
	    a[k + (i__ + 1) * a_dim1] = a12 * rd2;
/*<         a(k+1,i+1) = a22 * rd2 >*/
	    a[k + 1 + (i__ + 1) * a_dim1] = a22 * rd2;
/*<         a(k+2,i+1) = a32 * rd2 >*/
	    a[k + 2 + (i__ + 1) * a_dim1] = a32 * rd2;
/*<         a(k+3,i+1) = a42 * rd2 >*/
	    a[k + 3 + (i__ + 1) * a_dim1] = a42 * rd2;
/*<       end do >*/
	}
/*<       do k1 = k, ldim - 1 >*/
	i__2 = *ldim - 1;
	for (k1 = k; k1 <= i__2; ++k1) {
/*<         a11 = a(k1,i) >*/
	    a11 = a[k1 + i__ * a_dim1];
/*<         a12 = a(k1,i+1) >*/
	    a12 = a[k1 + (i__ + 1) * a_dim1];
/*<         do j = 0, i - 1 >*/
	    i__3 = i__ - 1;
	    for (j = 0; j <= i__3; ++j) {
/*<           t11 = a(i,j) >*/
		t11 = a[i__ + j * a_dim1];
/*<           t21 = a(i+1,j) >*/
		t21 = a[i__ + 1 + j * a_dim1];
/*<           x11 = a(k1,j) >*/
		x11 = a[k1 + j * a_dim1];
/*<           a11 = a11 - t11 * x11 >*/
		a11 -= t11 * x11;
/*<           a12 = a12 - t21 * x11 >*/
		a12 -= t21 * x11;
/*<         end do >*/
	    }
/*<         a11 = a11 * rd1 >*/
	    a11 *= rd1;
/*<         a(k1,i) = a11 >*/
	    a[k1 + i__ * a_dim1] = a11;
/*<         a12 = a12 - a11 * r21 >*/
	    a12 -= a11 * r21;
/*<         a(k1,i+1) = a12 * rd2 >*/
	    a[k1 + (i__ + 1) * a_dim1] = a12 * rd2;
/*<       end do       >*/
	}
/*<       end do >*/
    }
/*<       do i = lrank, lup - 2, 2 >*/
    i__1 = *lup - 2;
    for (i__ = *lrank; i__ <= i__1; i__ += 2) {
/*<       a11 = zero >*/
	a11 = 0.;
/*<       a21 = zero >*/
	a21 = 0.;
/*<       a22 = zero >*/
	a22 = 0.;
/*<       do j = 0, lrank - 2, 2 >*/
	i__2 = *lrank - 2;
	for (j = 0; j <= i__2; j += 2) {
/*<         t11 = a(i,j) >*/
	    t11 = a[i__ + j * a_dim1];
/*<         t21 = a(i+1,j) >*/
	    t21 = a[i__ + 1 + j * a_dim1];
/*<         t12 = a(i,j+1) >*/
	    t12 = a[i__ + (j + 1) * a_dim1];
/*<         t22 = a(i+1,j+1) >*/
	    t22 = a[i__ + 1 + (j + 1) * a_dim1];
/*<         a11 = a11 - t11 * t11 >*/
	    a11 -= t11 * t11;
/*<         a21 = a21 - t11 * t21 >*/
	    a21 -= t11 * t21;
/*<         a22 = a22 - t21 * t21 >*/
	    a22 -= t21 * t21;
/*<         a11 = a11 - t12 * t12 >*/
	    a11 -= t12 * t12;
/*<         a21 = a21 - t12 * t22 >*/
	    a21 -= t12 * t22;
/*<         a22 = a22 - t22 * t22 >*/
	    a22 -= t22 * t22;
/*<       end do >*/
	}
/*<       if (lr1 .eq. 1) then >*/
	if (lr1 == 1) {
/*<         t11 = a(i,j) >*/
	    t11 = a[i__ + j * a_dim1];
/*<         t21 = a(i+1,j) >*/
	    t21 = a[i__ + 1 + j * a_dim1];
/*<         a11 = a11 - t11 * t11 >*/
	    a11 -= t11 * t11;
/*<         a21 = a21 - t11 * t21 >*/
	    a21 -= t11 * t21;
/*<         a22 = a22 - t21 * t21 >*/
	    a22 -= t21 * t21;
/*<       end if >*/
	}
/*<       b(i,i) = a11 >*/
	b[i__ + i__ * b_dim1] = a11;
/*<       b(i+1,i) = a21 >*/
	b[i__ + 1 + i__ * b_dim1] = a21;
/*<       b(i+1,i+1) = a22 >*/
	b[i__ + 1 + (i__ + 1) * b_dim1] = a22;
/*<       do k = i+2, ldim - 4, 4 >*/
	i__2 = *ldim - 4;
	for (k = i__ + 2; k <= i__2; k += 4) {
/*<         a11 = zero >*/
	    a11 = 0.;
/*<         a21 = zero >*/
	    a21 = 0.;
/*<         a31 = zero >*/
	    a31 = 0.;
/*<         a41 = zero >*/
	    a41 = 0.;
/*<         a12 = zero >*/
	    a12 = 0.;
/*<         a22 = zero >*/
	    a22 = 0.;
/*<         a32 = zero >*/
	    a32 = 0.;
/*<         a42 = zero >*/
	    a42 = 0.;
/*<         do j = 0, lrank - 2, 2 >*/
	    i__3 = *lrank - 2;
	    for (j = 0; j <= i__3; j += 2) {
/*<           t11 = a(i,j) >*/
		t11 = a[i__ + j * a_dim1];
/*<           t21 = a(i+1,j) >*/
		t21 = a[i__ + 1 + j * a_dim1];
/*<           t12 = a(i,j+1) >*/
		t12 = a[i__ + (j + 1) * a_dim1];
/*<           t22 = a(i+1,j+1) >*/
		t22 = a[i__ + 1 + (j + 1) * a_dim1];
/*<           x11 = a(k,j) >*/
		x11 = a[k + j * a_dim1];
/*<           x21 = a(k+1,j) >*/
		x21 = a[k + 1 + j * a_dim1];
/*<           x12 = a(k,j+1) >*/
		x12 = a[k + (j + 1) * a_dim1];
/*<           x22 = a(k+1,j+1) >*/
		x22 = a[k + 1 + (j + 1) * a_dim1];
/*<           a11 = a11 - x11 * t11 >*/
		a11 -= x11 * t11;
/*<           a21 = a21 - x21 * t11 >*/
		a21 -= x21 * t11;
/*<           a12 = a12 - x11 * t21 >*/
		a12 -= x11 * t21;
/*<           a22 = a22 - x21 * t21 >*/
		a22 -= x21 * t21;
/*<           a11 = a11 - x12 * t12 >*/
		a11 -= x12 * t12;
/*<           a21 = a21 - x22 * t12 >*/
		a21 -= x22 * t12;
/*<           a12 = a12 - x12 * t22 >*/
		a12 -= x12 * t22;
/*<           a22 = a22 - x22 * t22 >*/
		a22 -= x22 * t22;
/*<           x31 = a(k+2,j) >*/
		x31 = a[k + 2 + j * a_dim1];
/*<           x41 = a(k+3,j) >*/
		x41 = a[k + 3 + j * a_dim1];
/*<           x32 = a(k+2,j+1) >*/
		x32 = a[k + 2 + (j + 1) * a_dim1];
/*<           x42 = a(k+3,j+1) >*/
		x42 = a[k + 3 + (j + 1) * a_dim1];
/*<           a31 = a31 - x31 * t11 >*/
		a31 -= x31 * t11;
/*<           a41 = a41 - x41 * t11 >*/
		a41 -= x41 * t11;
/*<           a32 = a32 - x31 * t21 >*/
		a32 -= x31 * t21;
/*<           a42 = a42 - x41 * t21 >*/
		a42 -= x41 * t21;
/*<           a31 = a31 - x32 * t12 >*/
		a31 -= x32 * t12;
/*<           a41 = a41 - x42 * t12 >*/
		a41 -= x42 * t12;
/*<           a32 = a32 - x32 * t22 >*/
		a32 -= x32 * t22;
/*<           a42 = a42 - x42 * t22 >*/
		a42 -= x42 * t22;
/*<         end do >*/
	    }
/*<         if (lr1 .eq. 1) then >*/
	    if (lr1 == 1) {
/*<           t11 = a(i,j) >*/
		t11 = a[i__ + j * a_dim1];
/*<           t21 = a(i+1,j) >*/
		t21 = a[i__ + 1 + j * a_dim1];
/*<           x11 = a(k,j) >*/
		x11 = a[k + j * a_dim1];
/*<           x21 = a(k+1,j)  >*/
		x21 = a[k + 1 + j * a_dim1];
/*<           x31 = a(k+2,j)  >*/
		x31 = a[k + 2 + j * a_dim1];
/*<           x41 = a(k+3,j)  >*/
		x41 = a[k + 3 + j * a_dim1];
/*<           a11 = a11 - x11 * t11 >*/
		a11 -= x11 * t11;
/*<           a21 = a21 - x21 * t11 >*/
		a21 -= x21 * t11;
/*<           a31 = a31 - x31 * t11 >*/
		a31 -= x31 * t11;
/*<           a41 = a41 - x41 * t11 >*/
		a41 -= x41 * t11;
/*<           a12 = a12 - x11 * t21 >*/
		a12 -= x11 * t21;
/*<           a22 = a22 - x21 * t21 >*/
		a22 -= x21 * t21;
/*<           a32 = a32 - x31 * t21 >*/
		a32 -= x31 * t21;
/*<           a42 = a42 - x41 * t21 >*/
		a42 -= x41 * t21;
/*<         end if >*/
	    }
/*<         b(k,i) = a11  >*/
	    b[k + i__ * b_dim1] = a11;
/*<         b(k+1,i) = a21  >*/
	    b[k + 1 + i__ * b_dim1] = a21;
/*<         b(k+2,i) = a31 >*/
	    b[k + 2 + i__ * b_dim1] = a31;
/*<         b(k+3,i) = a41 >*/
	    b[k + 3 + i__ * b_dim1] = a41;
/*<         b(k,i+1) = a12  >*/
	    b[k + (i__ + 1) * b_dim1] = a12;
/*<         b(k+1,i+1) = a22  >*/
	    b[k + 1 + (i__ + 1) * b_dim1] = a22;
/*<         b(k+2,i+1) = a32  >*/
	    b[k + 2 + (i__ + 1) * b_dim1] = a32;
/*<         b(k+3,i+1) = a42  >*/
	    b[k + 3 + (i__ + 1) * b_dim1] = a42;
/*<       end do >*/
	}
/*<       do k1 = k, ldim - 1 >*/
	i__2 = *ldim - 1;
	for (k1 = k; k1 <= i__2; ++k1) {
/*<         a11 = zero >*/
	    a11 = 0.;
/*<         a12 = zero >*/
	    a12 = 0.;
/*<         do j = 0, lrank - 1 >*/
	    i__3 = *lrank - 1;
	    for (j = 0; j <= i__3; ++j) {
/*<           t11 = a(i,j) >*/
		t11 = a[i__ + j * a_dim1];
/*<           t21 = a(i+1,j) >*/
		t21 = a[i__ + 1 + j * a_dim1];
/*<           x11 = a(k1,j) >*/
		x11 = a[k1 + j * a_dim1];
/*<           a11 = a11 - t11 * x11 >*/
		a11 -= t11 * x11;
/*<           a12 = a12 - t21 * x11 >*/
		a12 -= t21 * x11;
/*<         end do >*/
	    }
/*<         b(k1,i) = a11 >*/
	    b[k1 + i__ * b_dim1] = a11;
/*<         b(k1,i+1) = a12 >*/
	    b[k1 + (i__ + 1) * b_dim1] = a12;
/*<       end do       >*/
	}
/*<       end do >*/
    }
/*<       if (i .eq. lup-1) then >*/
    if (i__ == *lup - 1) {
/*<       do k = i, ldim - 1 >*/
	i__1 = *ldim - 1;
	for (k = i__; k <= i__1; ++k) {
/*<         a11 = zero >*/
	    a11 = 0.;
/*<         do j = 0, lrank - 1 >*/
	    i__2 = *lrank - 1;
	    for (j = 0; j <= i__2; ++j) {
/*<           a11 = a11 - a(k,j) * a(i,j) >*/
		a11 -= a[k + j * a_dim1] * a[i__ + j * a_dim1];
/*<         end do >*/
	    }
/*<         b(k,i) = a11 >*/
	    b[k + i__ * b_dim1] = a11;
/*<       end do >*/
	}
/*<       end if >*/
    }
/*<       dfopts(5) = ddmax       >*/
    dfopts[5] = ddmax;
/*<       dfopts(6) = ddmin       >*/
    dfopts[6] = ddmin;
/*<       info=0 !Cmj >*/
    *info = 0;
/*<       return >*/
    return 0;
/*< 99    dfopts(5) = ddmax >*/
L99:
    dfopts[5] = ddmax;
/*<       dfopts(6) = ddmin >*/
    dfopts[6] = ddmin;
/*<       info=1 !Cmj replaced "return 1" >*/
    *info = 1;
/*<       end >*/
    return 0;
} /* dchol2_ */

/* ----------------------------------------------------------------------------- */
/*<       subroutine PUTRHS1(N,linds,rank,w,rhs) >*/
/* Subroutine */ int putrhs1_(integer *n, integer *linds, integer *rank, 
	doublereal *w, doublereal *rhs)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__;

/*<       integer N,linds(*),rank,ldr >*/
/*<       double precision  w(*),rhs(0:*) >*/
/*<       do i = 1, rank >*/
    /* Parameter adjustments */
    --w;
    --linds;

    /* Function Body */
    i__1 = *rank;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         rhs(linds(i)) = w(i) >*/
	rhs[linds[i__]] = w[i__];
/*<       end do >*/
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* putrhs1_ */

