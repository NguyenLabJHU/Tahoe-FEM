/* $Id: eparfact1i.c,v 1.1 2005-01-03 06:39:34 paklein Exp $ */
/* eparfact1i.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

/* debugging */
#undef __DO_DEBUG__
/* #define __DO_DEBUG__ 1 */

#include "mpi.h"
#include "pspases_f2c.h"
#include <stdio.h>

/* /+***************************************************************************+/ */
/* /+                                                                           +/ */
/* /+   (C) Copyright IBM Corporation, 1997                                     +/ */
/* /+   (C) Copyright Regents of the University of Minnesota, 1997              +/ */
/* /+                                                                           +/ */
/* /+   eparfact1i.f                                                            +/ */
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
/* /+ $Id: eparfact1i.c,v 1.1 2005-01-03 06:39:34 paklein Exp $ +/ */
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
/* Subroutine */ int eparfact1_(integer *n, integer *aptrs, integer *ainds, 
	integer *lptrs, integer *linds, integer *tptrs, integer *tinds, 
	integer *sup, integer *stak, integer *nstak, integer *root, integer *
	dd, integer *lgblk, integer *myid, integer *supinds, integer *dimstak,
	 integer *wsize0, integer *wsize1, integer *ibuflen, integer *dbuflen,
	 integer *iwspace, integer *node, integer *stakptr, integer *nptr, 
	integer *lc, integer *iptrs, integer *lcsize, integer *wsolvesize, 
	integer *wa1, MPI_Comm *comm)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    integer j, i__, k, l, m, is3, is4, is5, idh, idv, hdim, vdim, itss, jtss, 
	    mask1, maskc, inode, maskr;
    extern /* Subroutine */ int gen_lc_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *);
    integer csuptr, iwstri, rsuptr, index_h__, index_v__, ncols_u__, 
	    nrows_u__;

/*<       integer KONSTANT >*/
/*<       parameter(KONSTANT=100000) >*/
/*<       integer N,root,dd,lgblk,blk,myid,nstak(*),supinds(*),comm >*/
/*<       integer aptrs(2,0:*),ainds(*),lptrs(3,0:*),linds(*) >*/
/*<       integer tptrs(3,0:*),tinds(*),sup(*),stak(3,*) >*/
/*<       integer dimstak(*),wa1(*) >*/
/*<       integer rsuptr, csuptr, nptr, stakptr, rank, ibufptr, dbufptr >*/
/*<       integer hdim,vdim,dbuflen,wsize1,wsize0 >*/
/*<       integer uptr, ibuflen >*/
/*<       integer lc(*),iptrs(*),lcsize,wsolvesize !Cmj >*/
/*<       include 'mpif.h' >*/
/*<       integer mpistat(MPI_STATUS_SIZE),ierr,ip,level >*/
/* -*- fortran -*- */

/* double precision functions */
/*<       vdim = ishft(dd,-1) >*/
    /* Parameter adjustments */
    --wa1;
    --iptrs;
    --lc;
    --dimstak;
    --supinds;
    --nstak;
    stak -= 4;
    --sup;
    --tinds;
    --tptrs;
    --linds;
    --lptrs;
    --ainds;
    --aptrs;

    /* Function Body */
    vdim = lbit_shift(*dd, (ftnlen)-1);
/*<       hdim = dd - vdim >*/
    hdim = *dd - vdim;
/*<       idv = 0 >*/
    idv = 0;
/*<       j = ishft(myid,-1) >*/
    j = lbit_shift(*myid, (ftnlen)-1);
/*<       do 50 i = 0, vdim - 1 >*/
    i__1 = vdim - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
/*<         idv = ior(idv,ishft(iand(j,1),i)) >*/
	idv |= lbit_shift(j & 1, i__);
/*<         j = ishft(j,-2) >*/
	j = lbit_shift(j, (ftnlen)-2);
/*< 50    continue >*/
/* L50: */
    }
/*<       idh = 0 >*/
    idh = 0;
/*<       j = myid >*/
    j = *myid;
/*<       do 60 i = 0, hdim - 1 >*/
    i__1 = hdim - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
/*<         idh = ior(idh,ishft(iand(j,1),i)) >*/
	idh |= lbit_shift(j & 1, i__);
/*<         j = ishft(j,-2) >*/
	j = lbit_shift(j, (ftnlen)-2);
/*< 60    continue >*/
/* L60: */
    }
/*<       wsize1 = 1 >*/
    *wsize1 = 1;
/*<       wsize0 = 1 >*/
    *wsize0 = 1;
/*<       mask1 = ishft(1,dd-1) >*/
    mask1 = lbit_shift((ftnlen)1, *dd - 1);
/*<       node = root >*/
    *node = *root;
/*<       nptr = 1 >*/
    *nptr = 1;
/*<       stakptr = 1 >*/
    *stakptr = 1;
/*<       do while ((hdim+vdim) .ne. 0) >*/
    while(hdim + vdim != 0) {
/*<         maskr = ishft(ishft(1,vdim)-1,lgblk) >*/
	maskr = lbit_shift(lbit_shift((ftnlen)1, vdim) - 1, *lgblk);
/*<         maskc = ishft(ishft(1,hdim)-1,lgblk) >*/
	maskc = lbit_shift(lbit_shift((ftnlen)1, hdim) - 1, *lgblk);
/*<         index_v = ishft(idv,lgblk) >*/
	index_v__ = lbit_shift(idv, *lgblk);
/*<         index_h = ishft(idh,lgblk) >*/
	index_h__ = lbit_shift(idh, *lgblk);
/*<         ncols_u = 0 >*/
	ncols_u__ = 0;
/*<         nrows_u = 0 >*/
	nrows_u__ = 0;
/*<         do while (tptrs(2,node) .eq. 1) >*/
	while(tptrs[*node * 3 + 2] == 1) {
/*<           nstak(nptr) = node >*/
	    nstak[*nptr] = *node;
/*<           nptr = nptr + 1 >*/
	    ++(*nptr);
/*<           k = tptrs(3,node) >*/
	    k = tptrs[*node * 3 + 3];
/*<           if (k .ne. 0) then >*/
	    if (k != 0) {
/*<             if (sup(k) .eq. node) then >*/
		if (sup[k] == *node) {
/*<               csuptr = sup(k+2) >*/
		    csuptr = sup[k + 2];
/*<               l = sup(k+3) >*/
		    l = sup[k + 3];
/*<               rsuptr = csuptr + l >*/
		    rsuptr = csuptr + l;
/*<               i = 0 >*/
		    i__ = 0;
/*<               j = 0 >*/
		    j = 0;
/*<               do 70 l = l-1, 0, -1 >*/
		    for (--l; l >= 0; --l) {
/*<                 if (iand(maskr,supinds(csuptr+l)).eq.index_v) goto 75 >*/
			if ((maskr & supinds[csuptr + l]) == index_v__) {
			    goto L75;
			}
/*< 70            continue >*/
/* L70: */
		    }
/*< 75            do 30 m = 0, l  >*/
L75:
		    i__1 = l;
		    for (m = 0; m <= i__1; ++m) {
/*<                 if (iand(maskc,supinds(csuptr+m)).eq.index_h) goto 35 >*/
			if ((maskc & supinds[csuptr + m]) == index_h__) {
			    goto L35;
			}
/*< 30            continue >*/
/* L30: */
		    }
/*< 35            do 10 m = m, l >*/
L35:
		    i__1 = l;
		    for (m = m; m <= i__1; ++m) {
/*<                 inode = supinds(csuptr+m) >*/
			inode = supinds[csuptr + m];
/*<                 if (iand(maskc,inode) .eq. index_h) then >*/
			if ((maskc & inode) == index_h__) {
/*<                   supinds(csuptr+i) = inode >*/
			    supinds[csuptr + i__] = inode;
/*<                   i = i + 1 >*/
			    ++i__;
/*<                 end if >*/
			}
/*<                 if (iand(maskr,inode) .eq. index_v) then >*/
			if ((maskr & inode) == index_v__) {
/*<                   supinds(rsuptr+j) = inode >*/
			    supinds[rsuptr + j] = inode;
/*<                   j = j + 1 >*/
			    ++j;
/*<                 end if >*/
			}
/*< 10            continue >*/
/* L10: */
		    }
/*<               stak(1,stakptr) = node >*/
		    stak[*stakptr * 3 + 1] = *node;
/*<               stak(2,stakptr) = i >*/
		    stak[*stakptr * 3 + 2] = i__;
/*<               stak(3,stakptr) = j >*/
		    stak[*stakptr * 3 + 3] = j;
/*<               stakptr = stakptr + 1 >*/
		    ++(*stakptr);
/*<               if (i .gt. ncols_u) ncols_u = i >*/
		    if (i__ > ncols_u__) {
			ncols_u__ = i__;
		    }
/*<               if (j .gt. nrows_u) nrows_u = j >*/
		    if (j > nrows_u__) {
			nrows_u__ = j;
		    }
/*<             end if >*/
		}
/*<           end if >*/
	    }
/*<           node = tinds(tptrs(1,node)) >*/
	    *node = tinds[tptrs[*node * 3 + 1]];
/*<         end do >*/
	}
/*<         nstak(nptr) = node >*/
	nstak[*nptr] = *node;
/*<         nptr = nptr + 1 >*/
	++(*nptr);
/*<         k = tptrs(3,node) >*/
	k = tptrs[*node * 3 + 3];
/*<         csuptr = sup(k+2) >*/
	csuptr = sup[k + 2];
/*<         l = sup(k+3) >*/
	l = sup[k + 3];
/*<         rsuptr = csuptr + l >*/
	rsuptr = csuptr + l;
/*<         i = 0 >*/
	i__ = 0;
/*<         j = 0 >*/
	j = 0;
/*<         do 80 l = l-1, 0, -1 >*/
	for (--l; l >= 0; --l) {
/*<           if (iand(maskr,supinds(csuptr+l)).eq.index_v) goto 85 >*/
	    if ((maskr & supinds[csuptr + l]) == index_v__) {
		goto L85;
	    }
/*< 80      continue >*/
/* L80: */
	}
/*< 85      do 40 m = 0, l  >*/
L85:
	i__1 = l;
	for (m = 0; m <= i__1; ++m) {
/*<           if (iand(maskc,supinds(csuptr+m)).eq.index_h) goto 45 >*/
	    if ((maskc & supinds[csuptr + m]) == index_h__) {
		goto L45;
	    }
/*< 40      continue >*/
/* L40: */
	}
/*< 45      do 20 m = m, l >*/
L45:
	i__1 = l;
	for (m = m; m <= i__1; ++m) {
/*<           inode = supinds(csuptr+m) >*/
	    inode = supinds[csuptr + m];
/*<           if (iand(maskc,inode) .eq. index_h) then >*/
	    if ((maskc & inode) == index_h__) {
/*<             supinds(csuptr+i) = inode >*/
		supinds[csuptr + i__] = inode;
/*<             i = i + 1 >*/
		++i__;
/*<           end if >*/
	    }
/*<           if (iand(maskr,inode) .eq. index_v) then >*/
	    if ((maskr & inode) == index_v__) {
/*<             supinds(rsuptr+j) = inode >*/
		supinds[rsuptr + j] = inode;
/*<             j = j + 1 >*/
		++j;
/*<           end if >*/
	    }
/*< 20      continue >*/
/* L20: */
	}
/*<         stak(1,stakptr) = node >*/
	stak[*stakptr * 3 + 1] = *node;
/*<         stak(2,stakptr) = i >*/
	stak[*stakptr * 3 + 2] = i__;
/*<         stak(3,stakptr) = j >*/
	stak[*stakptr * 3 + 3] = j;
/*<         stakptr = stakptr + 1 >*/
	++(*stakptr);
/*<         if (i .gt. ncols_u) ncols_u = i >*/
	if (i__ > ncols_u__) {
	    ncols_u__ = i__;
	}
/*<         if (j .gt. nrows_u) nrows_u = j >*/
	if (j > nrows_u__) {
	    nrows_u__ = j;
	}
/*<         dimstak(hdim+vdim) = nrows_u >*/
	dimstak[hdim + vdim] = nrows_u__;
/*<         if (iand(myid,mask1) .eq. 0) then >*/
	if ((*myid & mask1) == 0) {
/*<           node = tinds(tptrs(1,node)) >*/
	    *node = tinds[tptrs[*node * 3 + 1]];
/*<         else >*/
	} else {
/*<           node = tinds(1+tptrs(1,node)) >*/
	    *node = tinds[tptrs[*node * 3 + 1] + 1];
/*<         end if >*/
	}
/*<         mask1 = ishft(mask1,-1) >*/
	mask1 = lbit_shift(mask1, (ftnlen)-1);
/*<         if (hdim .ne. vdim) then >*/
	if (hdim != vdim) {
/*<           hdim = hdim - 1 >*/
	    --hdim;
/*<           idh = iand(idh,not(ishft(1,hdim))) >*/
	    idh &= ~ lbit_shift((ftnlen)1, hdim);
/*<           wsize1 = max0(wsize1,ncols_u*nrows_u) >*/
/* Computing MAX */
	    i__1 = *wsize1, i__2 = ncols_u__ * nrows_u__;
	    *wsize1 = max(i__1,i__2);
/*<         else >*/
	} else {
/*<           vdim = vdim - 1 >*/
	    --vdim;
/*<           idv = iand(idv,not(ishft(1,vdim))) >*/
	    idv &= ~ lbit_shift((ftnlen)1, vdim);
/*<           wsize0 = max0(wsize0,ncols_u*nrows_u) >*/
/* Computing MAX */
	    i__1 = *wsize0, i__2 = ncols_u__ * nrows_u__;
	    *wsize0 = max(i__1,i__2);
/*<         end if >*/
	}
/*<       end do >*/
    }
/*<       j = stak(1,stakptr-1) >*/
    j = stak[(*stakptr - 1) * 3 + 1];
/*<       jtss = lptrs(2,j) >*/
    jtss = lptrs[j * 3 + 2];
/*<       do while (jtss .eq. 0 .and. j .lt. root) >*/
    while(jtss == 0 && j < *root) {
/*<       j = j + 1 >*/
	++j;
/*<         jtss = lptrs(2,j) >*/
	jtss = lptrs[j * 3 + 2];
/*<       end do >*/
    }
/*<       itss = lptrs(2,sup(tptrs(3,node))) >*/
    itss = lptrs[sup[tptrs[*node * 3 + 3]] * 3 + 2];
/*<       m = itss >*/
    m = itss;
/*<       itss = itss * itss >*/
    itss *= itss;
/*<       if (jtss .gt. 0) then >*/
    if (jtss > 0) {
/*<         jtss = jtss - ishft(jtss, -2) >*/
	jtss -= lbit_shift(jtss, (ftnlen)-2);
/*<         jtss = max(jtss*jtss,ishft(itss,-1)) >*/
/* Computing MAX */
	i__1 = jtss * jtss, i__2 = lbit_shift(itss, (ftnlen)-1);
	jtss = max(i__1,i__2);
/*<       else >*/
    } else {
/*<       jtss = itss >*/
	jtss = itss;
/*<       end if >*/
    }
/*<       dbuflen = jtss + 36000 * (dd + 4) >*/
    *dbuflen = jtss + (*dd + 4) * 36000;
/*<       ibuflen = max(N+ishft(N,-1)+1, ishft(m,3)+640*(dd + 1)) >*/
/* Computing MAX */
    i__1 = *n + lbit_shift(*n, (ftnlen)-1) + 1, i__2 = (m << 3) + (*dd + 1) * 
	    640;
    *ibuflen = max(i__1,i__2);
/*<       is3 = 1 >*/
    is3 = 1;
/*<       iptrs(node+node+2) = 0  >*/
    iptrs[*node + *node + 2] = 0;
/*<    >*/
    gen_lc_(node, &lc[1], &linds[1], &lptrs[1], &tinds[1], &tptrs[1], &sup[1]
	    , &iptrs[1], lcsize, &wa1[1], &nstak[*nptr], &is3, &is5, &is4, &
	    iwstri);
/*<       lcsize = is3-1 >*/
    *lcsize = is3 - 1;
/*<       wsolvesize = iwstri >*/
    *wsolvesize = iwstri;
/*<       iwspace = is5 * is5 + is4 >*/
    *iwspace = is5 * is5 + is4;
/*<       wsize0 = max(wsize0,itss) >*/
    *wsize0 = max(*wsize0,itss);
/*<       iwspace = max(iwspace,wsize0+wsize1) >*/
/* Computing MAX */
    i__1 = *iwspace, i__2 = *wsize0 + *wsize1;
    *iwspace = max(i__1,i__2);
/*<       dbuflen = max(min((dbuflen*3+1)/2,iwspace),KONSTANT) >*/
/* Computing MAX */
/* Computing MIN */
    i__2 = (*dbuflen * 3 + 1) / 2;
    i__1 = min(i__2,*iwspace);
    *dbuflen = max(i__1,100000);
/*<       return >*/
    return 0;
/*<       end >*/
} /* eparfact1_ */

