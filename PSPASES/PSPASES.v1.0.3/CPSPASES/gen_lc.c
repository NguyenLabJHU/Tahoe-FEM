/* $Id: gen_lc.c,v 1.1 2004-12-12 23:27:33 paklein Exp $ */
/* gen_lc.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

/* debugging */
#undef __DO_DEBUG__
/* #define __DO_DEBUG__ 1 */

#include "pspases_f2c.h"

/* /+***************************************************************************+/ */
/* /+                                                                           +/ */
/* /+   (C) Copyright IBM Corporation, 1997                                     +/ */
/* /+   (C) Copyright Regents of the University of Minnesota, 1997              +/ */
/* /+                                                                           +/ */
/* /+   gen_lc.f                                                                +/ */
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
/* /+ $Id: gen_lc.c,v 1.1 2004-12-12 23:27:33 paklein Exp $ +/ */
/* /+***************************************************************************+/ */
/*     recursive */
/*<    >*/
/* Subroutine */ int gen_lc_(integer *root, integer *lc, integer *linds, 
	integer *lptrs, integer *tinds, integer *tptrs, integer *sup, integer 
	*iptrs, integer *lcsize, integer *wa1, integer *wa2, integer *lctr, 
	integer *wsmy, integer *wstot, integer *wstri)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer maxwstri, i__, j, ln, it, kid, jbot, kbot, kptr;
    extern /* Subroutine */ int gen_l_(integer *, integer *, integer *, 
	    integer *, integer *, integer *);
    static integer itbot, itptr, kwstri;

/*<       integer root,lcsize,lctr,lc(*),linds(*),lptrs(3,0:*),sup(*) >*/
/*<       integer tinds(*),tptrs(3,0:*),iptrs(2,0:*), wsmy, wstot, wstri >*/
/*<       integer wa1(0:*),wa2(0:*) >*/
/*<       itbot = tptrs(3,root) >*/
    /* Parameter adjustments */
    --iptrs;
    --sup;
    --tptrs;
    --tinds;
    --lptrs;
    --linds;
    --lc;

    /* Function Body */
    itbot = tptrs[*root * 3 + 3];
/*<       jbot = sup(itbot) >*/
    jbot = sup[itbot];
/*<       wsmy = lptrs(2,jbot) >*/
    *wsmy = lptrs[jbot * 3 + 2];
/*<       wstri = wsmy >*/
    *wstri = *wsmy;
/*<       do while (tptrs(2,jbot) .eq. 1) >*/
    while(tptrs[jbot * 3 + 2] == 1) {
/*<         kid = tinds(tptrs(1,jbot)) >*/
	kid = tinds[tptrs[jbot * 3 + 1]];
/*<         kbot = sup(tptrs(3,kid)) >*/
	kbot = sup[tptrs[kid * 3 + 3]];
/*<         iptrs(1,kid) = lctr >*/
	iptrs[(kid << 1) + 1] = *lctr;
/*<         iptrs(1,kbot) = lctr >*/
	iptrs[(kbot << 1) + 1] = *lctr;
/*<         call gen_l(jbot,kid,ln,lc(lctr),linds,lptrs) >*/
	gen_l_(&jbot, &kid, &ln, &lc[*lctr], &linds[1], &lptrs[1]);
/*<         lctr = lctr + ln >*/
	*lctr += ln;
/*<         iptrs(2,kid) = ln >*/
	iptrs[(kid << 1) + 2] = ln;
/*<         iptrs(2,kbot) = ln >*/
	iptrs[(kbot << 1) + 2] = ln;
/*<         jbot = kbot >*/
	jbot = kbot;
/*<         j = lptrs(2,kbot) >*/
	j = lptrs[kbot * 3 + 2];
/*<         wsmy = max(wsmy,j) >*/
	*wsmy = max(*wsmy,j);
/*<         wstri = wstri + j >*/
	*wstri += j;
/*<       end do >*/
    }
/*<       wstot = 0 >*/
    *wstot = 0;
/*<       itptr = tptrs(2,jbot) >*/
    itptr = tptrs[jbot * 3 + 2];
/*<       maxwstri = 0 >*/
    maxwstri = 0;
/*<       do i = 0, itptr - 1 >*/
    i__1 = itptr - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
/*<         kid = tinds(tptrs(1,jbot)+i) >*/
	kid = tinds[tptrs[jbot * 3 + 1] + i__];
/*<         kbot = sup(tptrs(3,kid)) >*/
	kbot = sup[tptrs[kid * 3 + 3]];
/*<         iptrs(1,kid) = lctr >*/
	iptrs[(kid << 1) + 1] = *lctr;
/*<         iptrs(1,kbot) = lctr >*/
	iptrs[(kbot << 1) + 1] = *lctr;
/*<         call gen_l(jbot,kid,ln,lc(lctr),linds,lptrs) >*/
	gen_l_(&jbot, &kid, &ln, &lc[*lctr], &linds[1], &lptrs[1]);
/*<         lctr = lctr + ln  >*/
	*lctr += ln;
/*<         iptrs(2,kid) = ln >*/
	iptrs[(kid << 1) + 2] = ln;
/*<         iptrs(2,kbot) = ln >*/
	iptrs[(kbot << 1) + 2] = ln;
/*<    >*/
	gen_lc_(&kid, &lc[1], &linds[1], &lptrs[1], &tinds[1], &
		tptrs[1], &sup[1], &iptrs[1], lcsize, &wa1[itptr], &wa2[itptr]
		, lctr, &wa1[i__], &wa2[i__], &kwstri);
/*<         maxwstri = max(maxwstri,kwstri) >*/
	maxwstri = max(maxwstri,kwstri);
/*<       end do >*/
    }
/*<       wstri = wstri + maxwstri >*/
    *wstri += maxwstri;
/*<       if (itptr .gt. 0) then >*/
    if (itptr > 0) {
/*<        kptr = tptrs(1,jbot) >*/
	kptr = tptrs[jbot * 3 + 1];
/*<        do i = 1, itptr - 1 >*/
	i__1 = itptr - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<         if (wa1(i) .gt. wa1(0)) then >*/
	    if (wa1[i__] > wa1[0]) {
/*<           it = wa1(0) >*/
		it = wa1[0];
/*<           wa1(0) = wa1(i) >*/
		wa1[0] = wa1[i__];
/*<           wa1(i) = it >*/
		wa1[i__] = it;
/*<           it = wa2(0) >*/
		it = wa2[0];
/*<           wa2(0) = wa2(i) >*/
		wa2[0] = wa2[i__];
/*<           wa2(i) = it >*/
		wa2[i__] = it;
/*<           it = tinds(kptr) >*/
		it = tinds[kptr];
/*<           tinds(kptr) = tinds(kptr+i) >*/
		tinds[kptr] = tinds[kptr + i__];
/*<           tinds(kptr+i) = it >*/
		tinds[kptr + i__] = it;
/*<         end if >*/
	    }
/*<        end do >*/
	}
/*<        wsmy = max(wsmy, wa1(0)) >*/
	*wsmy = max(*wsmy,wa1[0]);
/*<        wstot = wa2(0) >*/
	*wstot = wa2[0];
/*<        do i = 1, itptr - 1 >*/
	i__1 = itptr - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<         wstot = max(wstot, wa1(i)*wa1(i)+wa2(i)) >*/
/* Computing MAX */
	    i__2 = *wstot, i__3 = wa1[i__] * wa1[i__] + wa2[i__];
	    *wstot = max(i__2,i__3);
/*<        end do >*/
	}
/*<       end if >*/
    }
/*<       sup(itbot+3) = wsmy        >*/
    sup[itbot + 3] = *wsmy;
/*<       return >*/
    return 0;
/*<       end >*/
} /* gen_lc__ */

/* ------------------------------------------------------------------------------ */
/*<       subroutine GEN_L(par,kid,ln,lc,linds,lptrs) >*/
/* Subroutine */ int gen_l_(integer *par, integer *kid, integer *ln, integer 
	*lc, integer *linds, integer *lptrs)
{
    static integer i__, j, k, kl, ilim, jlim;

/*<       integer par,kid,ln,lc(0:*),linds(*),lptrs(3,0:*) >*/
/*<       i = lptrs(3,par) >*/
    /* Parameter adjustments */
    --lptrs;
    --linds;

    /* Function Body */
    i__ = lptrs[*par * 3 + 3];
/*<       j = lptrs(3,kid)  >*/
    j = lptrs[*kid * 3 + 3];
/*<       ilim = i + lptrs(2,par) >*/
    ilim = i__ + lptrs[*par * 3 + 2];
/*<       jlim = j + lptrs(2,kid) >*/
    jlim = j + lptrs[*kid * 3 + 2];
/*<       j = j + 1 >*/
    ++j;
/*<       ln = 0 >*/
    *ln = 0;
/*<       if (j .eq. jlim) then >*/
    if (j == jlim) {
/*<         ln = 2 >*/
	*ln = 2;
/*<         lc(0) = 0 >*/
	lc[0] = 0;
/*<         lc(1) = ilim - i >*/
	lc[1] = ilim - i__;
/*<       else >*/
    } else {
/*<         do while (j .lt. jlim .and. i .lt. ilim) >*/
	while(j < jlim && i__ < ilim) {
/*<           k = i >*/
	    k = i__;
/*<           do while (j .lt. jlim .and. linds(j) .eq. linds(i)) >*/
	    while(j < jlim && linds[j] == linds[i__]) {
/*<             i = i + 1 >*/
		++i__;
/*<             j = j + 1 >*/
		++j;
/*<           end do  >*/
	    }
/*<           lc(ln) = i - k >*/
	    lc[*ln] = i__ - k;
/*<           if (j .lt. jlim) then >*/
	    if (j < jlim) {
/*<             k = i >*/
		k = i__;
/*<             kl = linds(j) >*/
		kl = linds[j];
/*<             do while (linds(i) .ne. kl) >*/
		while(linds[i__] != kl) {
/*<               i = i + 1 >*/
		    ++i__;
/*<             end do >*/
		}
/*<             lc(ln+1) = i - k >*/
		lc[*ln + 1] = i__ - k;
/*<           else >*/
	    } else {
/*<             lc(ln+1) = ilim - i >*/
		lc[*ln + 1] = ilim - i__;
/*<           end if >*/
	    }
/*<           ln = ln + 2 >*/
	    *ln += 2;
/*<         end do >*/
	}
/*<       end if >*/
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* gen_l__ */

