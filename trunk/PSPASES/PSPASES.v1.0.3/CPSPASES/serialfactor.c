/* $Id: serialfactor.c,v 1.2 2004-12-28 18:01:25 paklein Exp $ */
/* serialfactor.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

/* debugging */
#undef __DO_DEBUG__
/* #define __DO_DEBUG__ 1 */

#include "pspases_f2c.h"

/* Table of constant values */

static doublereal c_b32 = 1.;
static doublereal c_b35 = -1.;

/* /+***************************************************************************+/ */
/* /+                                                                           +/ */
/* /+   (C) Copyright IBM Corporation, 1997                                     +/ */
/* /+   (C) Copyright Regents of the University of Minnesota, 1997              +/ */
/* /+                                                                           +/ */
/* /+   serialfactor.f                                                          +/ */
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
/* /+ $Id: serialfactor.c,v 1.2 2004-12-28 18:01:25 paklein Exp $ +/ */
/* /+***************************************************************************+/ */

static int max(int a, int b) {
	return (a > b) ? a : b;
};

/*     recursive */
/*    + */
/*<    >*/
/* Subroutine */ int factor6_(doublereal *wmem, integer *linds, integer *
	lptrs, integer *ainds, integer *aptrs, doublereal *avals, doublereal *
	lvals, integer *tinds, integer *tptrs, integer *sup, integer *rank, 
	integer *ldf, integer *fdim, integer *myfrontptr, integer *root, 
	integer *wm, integer *mystak, integer *lc, integer *iptrs, doublereal 
	*dfopts, integer *ifopts, integer *info)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;

    /* Local variables */
    integer frontptr, i__, j, k, i1, j1;
    integer jfrontptr, mystakptr, newstakptr, kid, ldl, kkk, ldu, myi, myj, 
	    locf, node;
    extern /* Subroutine */ int mydc_(integer *, doublereal *, doublereal *);
    integer udim, locu, kincr, mylim;
    extern /* Subroutine */ int dtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), dsyrk_(
	    char *, char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen),
	     dchol2_(doublereal *, integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, integer *), dpotrf_(char *, integer *, doublereal *, integer *,
	     integer *, ftnlen);
    integer locptr, myinds, locindi, locindj, kidincr;

/*<       integer lptrs(3,0:*),aptrs(2,0:*),ainds(*),tptrs(3,0:*) >*/
/*<       integer linds(*),tinds(*),sup(*),mystak(*),lc(*) >*/
/*<       integer myfrontptr,rank,ldf,fdim,root,wm,iptrs(2,0:*) >*/
/*<       integer ifopts(5) >*/
/*<       double precision wmem(*),avals(*),lvals(*),dfopts(7) >*/
/*<       integer frontptr,mystakptr,ldu,udim,locptr >*/
/*<       parameter (LIMDENSE = 50, LIMDRK = 9) >*/
/*<       mystakptr = 1 >*/
    /* Parameter adjustments */
    --ifopts;
    --dfopts;
    --iptrs;
    --lc;
    --mystak;
    --sup;
    --tptrs;
    --tinds;
    --lvals;
    --avals;
    --aptrs;
    --ainds;
    --lptrs;
    --linds;
    --wmem;

    /* Function Body */
    mystakptr = 1;
/*<       node = root >*/
    node = *root;
/*<       k = sup(tptrs(3,node)) >*/
    k = sup[tptrs[node * 3 + 3]];
/*<       ldf = max(ldf,sup(tptrs(3,root)+3)) >*/
/* Computing MAX */
    i__1 = *ldf, i__2 = sup[tptrs[*root * 3 + 3] + 3];
    *ldf = max(i__1,i__2);
/*<       do while (tptrs(2,k) .eq. 1) >*/
    while(tptrs[k * 3 + 2] == 1) {
/*<         mystak(mystakptr) = node >*/
	mystak[mystakptr] = node;
/*<         node = tinds(tptrs(1,k)) >*/
	node = tinds[tptrs[k * 3 + 1]];
/*<         mystakptr = mystakptr + 2 >*/
	mystakptr += 2;
/*<         k = sup(tptrs(3,node)) >*/
	k = sup[tptrs[node * 3 + 3]];
/*<       end do >*/
    }
/*<       jfrontptr = ldf * ldf + 1 >*/
    jfrontptr = *ldf * *ldf + 1;
/*<       fdim = lptrs(2,k) >*/
    *fdim = lptrs[k * 3 + 2];
/*<       do 11 i = mystakptr-2, 1, -2 >*/
    for (i__ = mystakptr - 2; i__ >= 1; i__ += -2) {
/*<         j = mystak(i) >*/
	j = mystak[i__];
/*<         j1 = lptrs(2,sup(tptrs(3,j))) >*/
	j1 = lptrs[sup[tptrs[j * 3 + 3]] * 3 + 2];
/*<         mystak(i+1) = jfrontptr - j1*ldf + ldf - j1 >*/
	mystak[i__ + 1] = jfrontptr - j1 * *ldf + *ldf - j1;
/*< 11    continue >*/
/* L11: */
    }
/*<       myfrontptr = jfrontptr - fdim * ldf + ldf - fdim >*/
    *myfrontptr = jfrontptr - *fdim * *ldf + *ldf - *fdim;
/*<       myinds = lptrs(3,k) >*/
    myinds = lptrs[k * 3 + 3];
/*<       mylim = myinds + fdim >*/
    mylim = myinds + *fdim;
/*<       frontptr = 1 >*/
    frontptr = 1;
/*<       ldu = ldf >*/
    ldu = *ldf;
/*<       do 20 i = 0, tptrs(2,k) - 1 >*/
    i__1 = tptrs[k * 3 + 2] - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
/*<       if (i .gt. 0) then >*/
	if (i__ > 0) {
/*<         frontptr = jfrontptr >*/
	    frontptr = jfrontptr;
/*<         ldu = 0 >*/
	    ldu = 0;
/*<         end if >*/
	}
/*<         kid = tinds(tptrs(1,k)+i) >*/
	kid = tinds[tptrs[k * 3 + 1] + i__];
/*<       info = 0  >*/
	*info = 0;
/*<    >*/
	i__2 = *wm - frontptr + 1;
	factor6_(&wmem[frontptr], &linds[1], &lptrs[1], &ainds[1], 
		&aptrs[1], &avals[1], &lvals[1], &tinds[1], &tptrs[1], &sup[1]
		, rank, &ldu, &udim, &kkk, &kid, &i__2, &mystak[mystakptr], &
		lc[1], &iptrs[1], &dfopts[1], &ifopts[1], info);
/*<       if (info.gt.0) goto 1111  >*/
	if (*info > 0) {
	    goto L1111;
	}
/*<         locptr = iptrs(2,kid) - 1 >*/
	locptr = iptrs[(kid << 1) + 2] - 1;
/*<         locu = frontptr + (ldu+1) * rank + kkk - 1 >*/
	locu = frontptr + (ldu + 1) * *rank + kkk - 1;
/*<         locf = myfrontptr - 1 >*/
	locf = *myfrontptr - 1;
/*<         kincr = 1 >*/
	kincr = 1;
/*<         kidincr = 1 >*/
	kidincr = 1;
/*<         myj = 1 >*/
	myj = 1;
/*<         locindj = 0 >*/
	locindj = 0;
/*<         if (i .eq. 0) then >*/
	if (i__ == 0) {
/*<           do 15 j1 = iptrs(1,kid), iptrs(1,kid) + locptr - 2, 2 >*/
	    i__2 = iptrs[(kid << 1) + 1] + locptr - 2;
	    for (j1 = iptrs[(kid << 1) + 1]; j1 <= i__2; j1 += 2) {
/*<           locindj = locindj + lc(j1) >*/
		locindj += lc[j1];
/*<             do 25 myj = myj, locindj >*/
		i__3 = locindj;
		for (myj = myj; myj <= i__3; ++myj) {
/*<               myi = kincr >*/
		    myi = kincr;
/*<             locindi = locindj - lc(j1) >*/
		    locindi = locindj - lc[j1];
/*<               do 35 i1 = j1, iptrs(1,kid) + locptr, 2 >*/
		    i__4 = iptrs[(kid << 1) + 1] + locptr;
		    for (i1 = j1; i1 <= i__4; i1 += 2) {
/*<                 locu = locu - myi >*/
			locu -= myi;
/*<             locindi = locindi + lc(i1) >*/
			locindi += lc[i1];
/*<                 do 46 myi = myi, locindi  >*/
			i__5 = locindi;
			for (myi = myi; myi <= i__5; ++myi) {
/*<                   wmem(locf+myi) = wmem(locu+myi) >*/
			    wmem[locf + myi] = wmem[locu + myi];
/*< 46              continue >*/
/* L46: */
			}
/*<                 locu = locu + myi >*/
			locu += myi;
/*<             locindi = locindi + lc(i1+1) >*/
			locindi += lc[i1 + 1];
/*<                 do 56 myi = myi, locindi >*/
			i__5 = locindi;
			for (myi = myi; myi <= i__5; ++myi) {
/*<                   wmem(locf+myi) = 0.d0 >*/
			    wmem[locf + myi] = 0.;
/*< 56              continue >*/
/* L56: */
			}
/*< 35            continue >*/
/* L35: */
		    }
/*<               locu = locu + rank + kidincr + ldu - udim >*/
		    locu = locu + *rank + kidincr + ldu - udim;
/*<               locf = locf + ldf >*/
		    locf += *ldf;
/*<               kincr = kincr + 1 >*/
		    ++kincr;
/*<               kidincr = kidincr + 1 >*/
		    ++kidincr;
/*< 25          continue >*/
/* L25: */
		}
/*<           locindj = locindj + lc(j1+1) >*/
		locindj += lc[j1 + 1];
/*<             do 65 myj = myj, locindj >*/
		i__3 = locindj;
		for (myj = myj; myj <= i__3; ++myj) {
/*<               do 76 i1 = kincr, fdim >*/
		    i__4 = *fdim;
		    for (i1 = kincr; i1 <= i__4; ++i1) {
/*<                 wmem(locf+i1) = 0.d0 >*/
			wmem[locf + i1] = 0.;
/*< 76            continue >*/
/* L76: */
		    }
/*<               locf = locf + ldf >*/
		    locf += *ldf;
/*<               kincr = kincr + 1 >*/
		    ++kincr;
/*< 65          continue >*/
/* L65: */
		}
/*< 15        continue >*/
/* L15: */
	    }
/*<         if (locf+kincr .ne. locu) then >*/
	    if (locf + kincr != locu) {
/*<           locindj = locindj + lc(j1) >*/
		locindj += lc[j1];
/*<             do 253 myj = myj, locindj >*/
		i__2 = locindj;
		for (myj = myj; myj <= i__2; ++myj) {
/*<               myi = kincr >*/
		    myi = kincr;
/*<             locindi = locindj - lc(j1) >*/
		    locindi = locindj - lc[j1];
/*<               do 353 i1 = j1, iptrs(1,kid) + locptr, 2 >*/
		    i__3 = iptrs[(kid << 1) + 1] + locptr;
		    for (i1 = j1; i1 <= i__3; i1 += 2) {
/*<                 locu = locu - myi >*/
			locu -= myi;
/*<             locindi = locindi + lc(i1) >*/
			locindi += lc[i1];
/*<                 do 463 myi = myi, locindi  >*/
			i__4 = locindi;
			for (myi = myi; myi <= i__4; ++myi) {
/*<                   wmem(locf+myi) = wmem(locu+myi) >*/
			    wmem[locf + myi] = wmem[locu + myi];
/*< 463             continue >*/
/* L463: */
			}
/*<                 locu = locu + myi >*/
			locu += myi;
/*<               locindi = locindi + lc(i1+1) >*/
			locindi += lc[i1 + 1];
/*<                 do 563 myi = myi, locindi >*/
			i__4 = locindi;
			for (myi = myi; myi <= i__4; ++myi) {
/*<                   wmem(locf+myi) = 0.d0 >*/
			    wmem[locf + myi] = 0.;
/*< 563             continue >*/
/* L563: */
			}
/*< 353           continue >*/
/* L353: */
		    }
/*<               locu = locu + rank + kidincr + ldu - udim >*/
		    locu = locu + *rank + kidincr + ldu - udim;
/*<               locf = locf + ldf >*/
		    locf += *ldf;
/*<               kincr = kincr + 1 >*/
		    ++kincr;
/*<               kidincr = kidincr + 1 >*/
		    ++kidincr;
/*< 253         continue >*/
/* L253: */
		}
/*<           locindj = locindj + lc(j1+1) >*/
		locindj += lc[j1 + 1];
/*<             do 653 myj = myj, locindj >*/
		i__2 = locindj;
		for (myj = myj; myj <= i__2; ++myj) {
/*<               do 763 i1 = kincr, fdim >*/
		    i__3 = *fdim;
		    for (i1 = kincr; i1 <= i__3; ++i1) {
/*<                 wmem(locf+i1) = 0.d0 >*/
			wmem[locf + i1] = 0.;
/*< 763           continue >*/
/* L763: */
		    }
/*<               locf = locf + ldf >*/
		    locf += *ldf;
/*<               kincr = kincr + 1 >*/
		    ++kincr;
/*< 653         continue >*/
/* L653: */
		}
/*<         end if >*/
	    }
/*<         else >*/
	} else {
/*<           do 115 j1 = iptrs(1,kid), iptrs(1,kid) + locptr, 2 >*/
	    i__2 = iptrs[(kid << 1) + 1] + locptr;
	    for (j1 = iptrs[(kid << 1) + 1]; j1 <= i__2; j1 += 2) {
/*<           locindj = locindj + lc(j1) >*/
		locindj += lc[j1];
/*<             do 125 myj = myj, locindj >*/
		i__3 = locindj;
		for (myj = myj; myj <= i__3; ++myj) {
/*<               myi = kincr >*/
		    myi = kincr;
/*<             locindi = locindj - lc(j1) >*/
		    locindi = locindj - lc[j1];
/*<               do 135 i1 = j1, iptrs(1,kid) + locptr, 2 >*/
		    i__4 = iptrs[(kid << 1) + 1] + locptr;
		    for (i1 = j1; i1 <= i__4; i1 += 2) {
/*<                 locu = locu - myi >*/
			locu -= myi;
/*<             locindi = locindi + lc(i1) >*/
			locindi += lc[i1];
/*<                 do 146 myi = myi, locindi >*/
			i__5 = locindi;
			for (myi = myi; myi <= i__5; ++myi) {
/*<                   wmem(locf+myi)=wmem(locf+myi)+wmem(locu+myi) >*/
			    wmem[locf + myi] += wmem[locu + myi];
/*< 146             continue >*/
/* L146: */
			}
/*<                 locu = locu + myi >*/
			locu += myi;
/*<             locindi = locindi + lc(i1+1) >*/
			locindi += lc[i1 + 1];
/*<                 myi = locindi + 1 >*/
			myi = locindi + 1;
/*< 135           continue >*/
/* L135: */
		    }
/*<               locu = locu + rank + kidincr + ldu - udim >*/
		    locu = locu + *rank + kidincr + ldu - udim;
/*<               locf = locf + ldf >*/
		    locf += *ldf;
/*<               kincr = kincr + 1 >*/
		    ++kincr;
/*<               kidincr = kidincr + 1 >*/
		    ++kidincr;
/*< 125         continue >*/
/* L125: */
		}
/*<           locindj = locindj + lc(j1+1) >*/
		locindj += lc[j1 + 1];
/*<             j = locindj + 1 - myj >*/
		j = locindj + 1 - myj;
/*<             locf = locf + ldf * j >*/
		locf += *ldf * j;
/*<             kincr = kincr + j >*/
		kincr += j;
/*<             myj = myj + j >*/
		myj += j;
/*< 115       continue >*/
/* L115: */
	    }
/*<         end if >*/
	}
/*< 20    continue >*/
/* L20: */
    }
/*<       rank = sup(tptrs(3,k)+1) >*/
    *rank = sup[tptrs[k * 3 + 3] + 1];
/*<       j = node >*/
    j = node;
/*<       if (tptrs(2,k) .eq. 0) then >*/
    if (tptrs[k * 3 + 2] == 0) {
/*<         locf = lptrs(1,node) >*/
	locf = lptrs[node * 3 + 1];
/*<         ldl = lptrs(2,k) >*/
	ldl = lptrs[k * 3 + 2];
/*<         do 80 j1 = 1, rank >*/
	i__1 = *rank;
	for (j1 = 1; j1 <= i__1; ++j1) {
/*<           myi = lptrs(3,j) >*/
	    myi = lptrs[j * 3 + 3];
/*<           locf = locf - myi >*/
	    locf -= myi;
/*<           do 90 i1 = aptrs(1,j), aptrs(1,j) + aptrs(2,j) - 1 >*/
	    i__2 = aptrs[(j << 1) + 1] + aptrs[(j << 1) + 2] - 1;
	    for (i1 = aptrs[(j << 1) + 1]; i1 <= i__2; ++i1) {
/*<             do while (linds(myi) .ne. ainds(i1)) >*/
		while(linds[myi] != ainds[i1]) {
/*<               lvals(locf + myi) = 0.d0 >*/
		    lvals[locf + myi] = 0.;
/*<               myi = myi + 1 >*/
		    ++myi;
/*<             end do >*/
		}
/*<             lvals(locf + myi) = avals(i1) >*/
		lvals[locf + myi] = avals[i1];
/*<             myi = myi + 1 >*/
		++myi;
/*< 90        continue >*/
/* L90: */
	    }
/*<           do 95 myi = myi, lptrs(3,j)+lptrs(2,j)-1 >*/
	    i__2 = lptrs[j * 3 + 3] + lptrs[j * 3 + 2] - 1;
	    for (myi = myi; myi <= i__2; ++myi) {
/*<             lvals(locf + myi) = 0.d0 >*/
		lvals[locf + myi] = 0.;
/*< 95        continue >*/
/* L95: */
	    }
/*<           locf = locf + lptrs(3,j) - ldl - 1 >*/
	    locf = locf + lptrs[j * 3 + 3] - ldl - 1;
/*<           j = tinds(tptrs(1,j)) >*/
	    j = tinds[tptrs[j * 3 + 1]];
/*< 80      continue >*/
/* L80: */
	}
/*<         locf = lptrs(1,k) >*/
	locf = lptrs[k * 3 + 1];
/*<       info = 0  >*/
	*info = 0;
/*<    >*/
	dchol2_(&lvals[locf], &ldl, &wmem[*myfrontptr], ldf, fdim, rank, fdim,
		 &linds[lptrs[k * 3 + 3]], &dfopts[1], &ifopts[1], info);
/*<       if (info.gt.0) goto 1111  >*/
	if (*info > 0) {
	    goto L1111;
	}
/*<       else >*/
    } else {
/*<         locf = myfrontptr + (ldf+1) * (rank-1) >*/
	locf = *myfrontptr + (*ldf + 1) * (*rank - 1);
/*<       newstakptr = mystakptr + rank >*/
	newstakptr = mystakptr + *rank;
/*<         do 130 j1 = rank, 1, -1 >*/
	for (j1 = *rank; j1 >= 1; --j1) {
/*<         mystak(newstakptr - j1) = j >*/
	    mystak[newstakptr - j1] = j;
/*<           myi = lptrs(3,j) >*/
	    myi = lptrs[j * 3 + 3];
/*<           locf = locf - lptrs(3,j) >*/
	    locf -= lptrs[j * 3 + 3];
/*<           do 140 i1 = aptrs(1,j), aptrs(1,j) + aptrs(2,j) - 1 >*/
	    i__1 = aptrs[(j << 1) + 1] + aptrs[(j << 1) + 2] - 1;
	    for (i1 = aptrs[(j << 1) + 1]; i1 <= i__1; ++i1) {
/*<             do while (linds(myi) .ne. ainds(i1)) >*/
		while(linds[myi] != ainds[i1]) {
/*<               myi = myi + 1 >*/
		    ++myi;
/*<             end do >*/
		}
/*<             wmem(locf + myi) = wmem(locf + myi) + avals(i1) >*/
		wmem[locf + myi] += avals[i1];
/*< 140       continue >*/
/* L140: */
	    }
/*<           locf = locf + lptrs(3,j) - ldf - 1 >*/
	    locf = locf + lptrs[j * 3 + 3] - *ldf - 1;
/*<           j = tinds(tptrs(1,j)) >*/
	    j = tinds[tptrs[j * 3 + 1]];
/*< 130     continue >*/
/* L130: */
	}
/*<         info = 0  >*/
	*info = 0;
/*<           call dpotrf('l',rank,wmem(myfrontptr),ldf,info) >*/
	dpotrf_("l", rank, &wmem[*myfrontptr], ldf, info, (ftnlen)1);
/*<         if (info.gt.0) goto 1111  >*/
	if (*info > 0) {
	    goto L1111;
	}
/*<    >*/
	i__1 = *fdim - *rank;
	dtrsm_("R", "L", "T", "N", &i__1, rank, &c_b32, &wmem[*myfrontptr], 
		ldf, &wmem[*myfrontptr + *rank], ldf, (ftnlen)1, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);
/*<    >*/
	i__1 = *fdim - *rank;
	dsyrk_("L", "N", &i__1, rank, &c_b35, &wmem[*myfrontptr + *rank], ldf,
		 &c_b32, &wmem[*myfrontptr + (*ldf + 1) * *rank], ldf, (
		ftnlen)1, (ftnlen)1);
/*<       locf = myfrontptr >*/
	locf = *myfrontptr;
/*<       do j1 = newstakptr - 1, mystakptr, -1 >*/
	i__1 = mystakptr;
	for (j1 = newstakptr - 1; j1 >= i__1; --j1) {
/*<         j = mystak(j1) >*/
	    j = mystak[j1];
/*<         call mydc(lptrs(2,j),wmem(locf),lvals(lptrs(1,j))) >*/
	    mydc_(&lptrs[j * 3 + 2], &wmem[locf], &lvals[lptrs[j * 3 + 1]]);
/*<         locf = locf + ldf + 1 >*/
	    locf = locf + *ldf + 1;
/*<       end do >*/
	}
/*<       end if >*/
    }
/*<       do 10 i = mystakptr - 2, 1, -2 >*/
    for (i__ = mystakptr - 2; i__ >= 1; i__ += -2) {
/*<         udim = fdim >*/
	udim = *fdim;
/*<         kid = node >*/
	kid = node;
/*<         frontptr = myfrontptr >*/
	frontptr = *myfrontptr;
/*<         node = mystak(i) >*/
	node = mystak[i__];
/*<         myfrontptr = mystak(i+1) >*/
	*myfrontptr = mystak[i__ + 1];
/*<         k = sup(tptrs(3,node)) >*/
	k = sup[tptrs[node * 3 + 3]];
/*<         myinds = lptrs(3,k) >*/
	myinds = lptrs[k * 3 + 3];
/*<         fdim = lptrs(2,k) >*/
	*fdim = lptrs[k * 3 + 2];
/*<         myj = myinds >*/
	myj = myinds;
/*<         mylim = myinds + fdim >*/
	mylim = myinds + *fdim;
/*<         locptr = iptrs(2,kid) - 1 >*/
	locptr = iptrs[(kid << 1) + 2] - 1;
/*<         locu = frontptr + (ldf+1) * rank >*/
	locu = frontptr + (*ldf + 1) * *rank;
/*<         locf = myfrontptr - 1 >*/
	locf = *myfrontptr - 1;
/*<         kincr = 1 >*/
	kincr = 1;
/*<         kidincr = 1 >*/
	kidincr = 1;
/*<         myj = 1 >*/
	myj = 1;
/*<         locindj = 0 >*/
	locindj = 0;
/*<         do 215 j1 = iptrs(1,kid), iptrs(1,kid) + locptr - 2, 2 >*/
	i__1 = iptrs[(kid << 1) + 1] + locptr - 2;
	for (j1 = iptrs[(kid << 1) + 1]; j1 <= i__1; j1 += 2) {
/*<         locindj = locindj + lc(j1) >*/
	    locindj += lc[j1];
/*<           do 225 myj = myj, locindj >*/
	    i__2 = locindj;
	    for (myj = myj; myj <= i__2; ++myj) {
/*<             myi = kincr >*/
		myi = kincr;
/*<           locindi = locindj - lc(j1) >*/
		locindi = locindj - lc[j1];
/*<             do 235 i1 = j1, iptrs(1,kid) + locptr, 2 >*/
		i__3 = iptrs[(kid << 1) + 1] + locptr;
		for (i1 = j1; i1 <= i__3; i1 += 2) {
/*<               locu = locu - myi >*/
		    locu -= myi;
/*<             locindi = locindi + lc(i1) >*/
		    locindi += lc[i1];
/*<               do 246 myi = myi, locindi >*/
		    i__4 = locindi;
		    for (myi = myi; myi <= i__4; ++myi) {
/*<                 wmem(locf+myi) = wmem(locu+myi) >*/
			wmem[locf + myi] = wmem[locu + myi];
/*< 246           continue >*/
/* L246: */
		    }
/*<               locu = locu + myi >*/
		    locu += myi;
/*<             locindi = locindi + lc(i1+1) >*/
		    locindi += lc[i1 + 1];
/*<               do 256, myi = myi, locindi >*/
		    i__4 = locindi;
		    for (myi = myi; myi <= i__4; ++myi) {
/*<                 wmem(locf+myi) = 0.d0 >*/
			wmem[locf + myi] = 0.;
/*< 256           continue >*/
/* L256: */
		    }
/*< 235         continue >*/
/* L235: */
		}
/*<             locu = locu + rank + kidincr + ldf - udim >*/
		locu = locu + *rank + kidincr + *ldf - udim;
/*<             locf = locf + ldf >*/
		locf += *ldf;
/*<             kincr = kincr + 1 >*/
		++kincr;
/*<             kidincr = kidincr + 1 >*/
		++kidincr;
/*< 225       continue >*/
/* L225: */
	    }
/*<         locindj = locindj + lc(j1+1) >*/
	    locindj += lc[j1 + 1];
/*<           do 265 myj = myj, locindj >*/
	    i__2 = locindj;
	    for (myj = myj; myj <= i__2; ++myj) {
/*<             do 276 i1 = kincr, fdim >*/
		i__3 = *fdim;
		for (i1 = kincr; i1 <= i__3; ++i1) {
/*<               wmem(locf+i1) = 0.d0 >*/
		    wmem[locf + i1] = 0.;
/*< 276         continue >*/
/* L276: */
		}
/*<             locf = locf + ldf >*/
		locf += *ldf;
/*<             kincr = kincr + 1 >*/
		++kincr;
/*< 265       continue >*/
/* L265: */
	    }
/*< 215     continue >*/
/* L215: */
	}
/*<         if ((locf+kincr) .ne. locu) then >*/
	if (locf + kincr != locu) {
/*<         locindj = locindj + lc(j1) >*/
	    locindj += lc[j1];
/*<           do 525 myj = myj, locindj >*/
	    i__1 = locindj;
	    for (myj = myj; myj <= i__1; ++myj) {
/*<             myi = kincr >*/
		myi = kincr;
/*<           locindi = locindj - lc(j1) >*/
		locindi = locindj - lc[j1];
/*<             do 535 i1 = j1, iptrs(1,kid) + locptr, 2 >*/
		i__2 = iptrs[(kid << 1) + 1] + locptr;
		for (i1 = j1; i1 <= i__2; i1 += 2) {
/*<               locu = locu - myi >*/
		    locu -= myi;
/*<             locindi = locindi + lc(i1) >*/
		    locindi += lc[i1];
/*<               do 546 myi = myi, locindi >*/
		    i__3 = locindi;
		    for (myi = myi; myi <= i__3; ++myi) {
/*<                 wmem(locf+myi) = wmem(locu+myi) >*/
			wmem[locf + myi] = wmem[locu + myi];
/*< 546           continue >*/
/* L546: */
		    }
/*<               locu = locu + myi >*/
		    locu += myi;
/*<             locindi = locindi + lc(i1+1) >*/
		    locindi += lc[i1 + 1];
/*<               do 555 myi = myi, locindi >*/
		    i__3 = locindi;
		    for (myi = myi; myi <= i__3; ++myi) {
/*<                 wmem(locf+myi) = 0.d0 >*/
			wmem[locf + myi] = 0.;
/*< 555           continue >*/
/* L555: */
		    }
/*< 535         continue >*/
/* L535: */
		}
/*<             locu = locu + rank + kidincr + ldf - udim >*/
		locu = locu + *rank + kidincr + *ldf - udim;
/*<             kincr = kincr + 1 >*/
		++kincr;
/*<             locf = locf + ldf >*/
		locf += *ldf;
/*<             kidincr = kidincr + 1 >*/
		++kidincr;
/*< 525       continue >*/
/* L525: */
	    }
/*<         locindj = locindj + lc(j1+1) >*/
	    locindj += lc[j1 + 1];
/*<           do 565 myj = myj, locindj >*/
	    i__1 = locindj;
	    for (myj = myj; myj <= i__1; ++myj) {
/*<             do 576 i1 = kincr, fdim >*/
		i__2 = *fdim;
		for (i1 = kincr; i1 <= i__2; ++i1) {
/*<               wmem(locf+i1) = 0.d0 >*/
		    wmem[locf + i1] = 0.;
/*< 576         continue >*/
/* L576: */
		}
/*<             locf = locf + ldf >*/
		locf += *ldf;
/* mj            kicr = kicr + 1 */
/*< 565       continue >*/
/* L565: */
	    }
/*<         end if >*/
	}
/*<         rank = sup(tptrs(3,k)+1) >*/
	*rank = sup[tptrs[k * 3 + 3] + 1];
/*<         j = node >*/
	j = node;
/*<         locf = myfrontptr + (ldf+1) * (rank-1) >*/
	locf = *myfrontptr + (*ldf + 1) * (*rank - 1);
/*<       newstakptr = mystakptr + rank >*/
	newstakptr = mystakptr + *rank;
/*<         do 30 j1 = rank, 1, -1 >*/
	for (j1 = *rank; j1 >= 1; --j1) {
/*<         mystak(newstakptr - j1) = j >*/
	    mystak[newstakptr - j1] = j;
/*<           myi = lptrs(3,j) >*/
	    myi = lptrs[j * 3 + 3];
/*<           locf = locf - myi >*/
	    locf -= myi;
/*<           do 40 i1 = aptrs(1,j), aptrs(1,j) + aptrs(2,j) - 1 >*/
	    i__1 = aptrs[(j << 1) + 1] + aptrs[(j << 1) + 2] - 1;
	    for (i1 = aptrs[(j << 1) + 1]; i1 <= i__1; ++i1) {
/*<             do while (linds(myi) .ne. ainds(i1)) >*/
		while(linds[myi] != ainds[i1]) {
/*<               myi = myi + 1 >*/
		    ++myi;
/*<             end do >*/
		}
/*<             wmem(locf + myi) = wmem(locf + myi) + avals(i1) >*/
		wmem[locf + myi] += avals[i1];
/*< 40        continue >*/
/* L40: */
	    }
/*<           locf = locf + lptrs(3,j) - ldf - 1 >*/
	    locf = locf + lptrs[j * 3 + 3] - *ldf - 1;
/*<           j = tinds(tptrs(1,j)) >*/
	    j = tinds[tptrs[j * 3 + 1]];
/*< 30      continue >*/
/* L30: */
	}
/*<         info = 0  >*/
	*info = 0;
/*<           call dpotrf('l',rank,wmem(myfrontptr),ldf,info) >*/
	dpotrf_("l", rank, &wmem[*myfrontptr], ldf, info, (ftnlen)1);
/*<           if (info.gt.0) goto 1111  >*/
	if (*info > 0) {
	    goto L1111;
	}
/*<    >*/
	i__1 = *fdim - *rank;
	dtrsm_("R", "L", "T", "N", &i__1, rank, &c_b32, &wmem[*myfrontptr], 
		ldf, &wmem[*myfrontptr + *rank], ldf, (ftnlen)1, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);
/*<    >*/
	i__1 = *fdim - *rank;
	dsyrk_("L", "N", &i__1, rank, &c_b35, &wmem[*myfrontptr + *rank], ldf,
		 &c_b32, &wmem[*myfrontptr + (*ldf + 1) * *rank], ldf, (
		ftnlen)1, (ftnlen)1);
/*<       locf = myfrontptr >*/
	locf = *myfrontptr;
/*<       do j1 = newstakptr - 1, mystakptr, -1 >*/
	i__1 = mystakptr;
	for (j1 = newstakptr - 1; j1 >= i__1; --j1) {
/*<         j = mystak(j1) >*/
	    j = mystak[j1];
/*<         call mydc(lptrs(2,j),wmem(locf),lvals(lptrs(1,j))) >*/
	    mydc_(&lptrs[j * 3 + 2], &wmem[locf], &lvals[lptrs[j * 3 + 1]]);
/*<         locf = locf + ldf + 1 >*/
	    locf = locf + *ldf + 1;
/*<       end do >*/
	}
/*< 10    continue >*/
/* L10: */
    }
/*<       info=0  >*/
    *info = 0;
/*<       return >*/
    return 0;
/*< 1111  continue >*/
L1111:
/*<       info=1  >*/
    *info = 1;
/*<       end >*/
    return 0;
} /* factor6_ */

