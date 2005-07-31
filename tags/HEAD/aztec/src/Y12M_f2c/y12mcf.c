/* PAK (07/24/98) */

/* y12mcf.f -- translated by f2c (version 19960611).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "y12m.h"

#ifdef __MWERKS__
#pragma warn_unusedarg off
#endif /* __MWERKS__ */

/* Subroutine */ int y12mcf(n, z__, a, snr, nn, rnr, nn1, pivot, b, ha, iha, 
	aflag, iflag, ifail)
integer *n, *z__;
doublereal *a;
integer *snr, *nn, *rnr, *nn1;
doublereal *pivot, *b;
integer *ha, *iha;
doublereal *aflag;
integer *iflag, *ifail;
{
    /* System generated locals */
    integer ha_dim1, ha_offset, i__1, i__2, i__3, i__4;
    doublereal d__1;

    /* Local variables */
    static integer slut, rrow, i__, j, l, k, r__;
    static doublereal u, v, t;
    static integer index, rcoll;
    static doublereal grmin;
    static integer c1, c2, i1, l1, l2, l3, l4, l5, l6, r2, n7, n8, r4, r5, r7,
	     r8, r9, r6, r3, r1, l7, r10, jj, kk;
    static doublereal td;
    static integer ll, nr, rr, zz, cr1, cr2, cr3, cr4;
    static doublereal td1;
    static integer rpivot, rr1, rr2, rr3, rr4, lfc, lfr;
    static doublereal tol1, tol2, tol3;


/*  systens of linear equations are solved by use of sparse matrix tech- 
*/
/*  nique and by gaussian elimination. */


/*  information which is necessary to begin the elimination is stored. */

    /* Parameter adjustments */
    --b;
    --pivot;
    --snr;
    --a;
    --rnr;
    ha_dim1 = *iha;
    ha_offset = ha_dim1 + 1;
    ha -= ha_offset;
    --aflag;
    --iflag;

    /* Function Body */
    *ifail = 0;
    if (iflag[1] != -1) {
	*ifail = 2;
    }
    if (aflag[1] < 1.) {
	aflag[1] = 1.0005;
    }
    if (aflag[3] < 1e5) {
	aflag[3] = 1e5;
    }
    if (aflag[4] < 0.) {
	aflag[4] = -aflag[4];
    }
    if (iflag[2] < 1) {
	*ifail = 19;
    }
    if (iflag[3] < 0 || iflag[3] > 2) {
	*ifail = 20;
    }
    if (iflag[5] < 1 || iflag[5] > 3) {
	*ifail = 21;
    }
    if (iflag[5] == 3) {
	*ifail = 22;
    }
    if (*ifail > 0) {
	goto L1110;
    }
    snr[*z__ + 1] = 0;
    rnr[*z__ + 1] = 0;
    n8 = *n + 1;
    n7 = *n - 1;
    u = aflag[1];
    grmin = aflag[4] * aflag[6];

/*  use the information about fill-ins if it is possible. */

    zz = *z__;
    nr = *n * *n;
    if (iflag[4] != 2) {
	goto L100;
    }
    if (iflag[10] > *nn) {
	goto L50;
    }
    l1 = iflag[10];
    l5 = l1 + 1;
    if (l5 <= *nn) {
	snr[l5] = 0;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l = n8 - i__;
	l2 = ha[l + ha_dim1 * 3] + 1;
	l3 = l2 - ha[l + ha_dim1];
	i__2 = l3;
	for (j = 1; j <= i__2; ++j) {
	    snr[l5 - j] = snr[l2 - j];
/* L10: */
	    a[l5 - j] = a[l2 - j];
	}
	ha[l + ha_dim1 * 3] = l1;
	ha[l + ha_dim1] = l5 - l3;
	l6 = l1 - l3;
	l5 -= ha[l + ha_dim1 * 9];
	if (l5 > l6) {
	    goto L30;
	}
	i__2 = l6;
	for (j = l5; j <= i__2; ++j) {
/* L20: */
	    snr[j] = 0;
	}
L30:
/* L40: */
	l1 = l5 - 1;
    }
L50:
    if (iflag[9] > *nn1) {
	goto L100;
    }
    l2 = iflag[9];
    l5 = l2 + 1;
    if (l5 <= *nn1) {
	rnr[l5] = 0;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l = n8 - i__;
	l1 = ha[l + ha_dim1 * 6] + 1;
	l4 = l1 - ha[l + (ha_dim1 << 2)];
	i__2 = l4;
	for (j = 1; j <= i__2; ++j) {
/* L60: */
	    rnr[l5 - j] = rnr[l1 - j];
	}
	ha[l + (ha_dim1 << 2)] = l5 - l4;
	ha[l + ha_dim1 * 6] = l2;
	l6 = l2 - l4;
	l5 -= ha[l + ha_dim1 * 10];
	if (l5 > l6) {
	    goto L80;
	}
	i__2 = l6;
	for (j = l5; j <= i__2; ++j) {
/* L70: */
	    rnr[j] = 0;
	}
L80:
/* L90: */
	l2 = l5 - 1;
    }
L100:
    r4 = ha[*n + ha_dim1 * 3];
    r5 = ha[*n + ha_dim1 * 6];
    aflag[7] = aflag[6];
    aflag[8] = aflag[6];
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	pivot[i__] = 0.;
	ha[i__ + (ha_dim1 << 1)] = ha[i__ + ha_dim1];
/* L110: */
	ha[i__ + ha_dim1 * 5] = ha[i__ + (ha_dim1 << 2)];
    }
    index = ha[*n + (ha_dim1 << 3)];

/*  start of gaussian elimination. */

    slut = ha[index + ha_dim1 * 3] - ha[index + (ha_dim1 << 1)] + 1;
    i__1 = n7;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rr3 = ha[i__ + (ha_dim1 << 1)];
	rr4 = ha[i__ + ha_dim1 * 3];
	c1 = ha[i__ + (ha_dim1 << 2)];
	cr4 = ha[i__ + ha_dim1 * 6];
	if (iflag[3] == 0) {
	    goto L350;
	}
	if (iflag[4] != 2) {
	    goto L120;
	}
	rrow = ha[i__ + ha_dim1 * 7];
	rcoll = ha[i__ + (ha_dim1 << 3)];
	goto L220;
L120:
	l4 = ha[i__ + (ha_dim1 << 3)];
	if (iflag[3] == 1) {
	    goto L130;
	}
	rrow = l4;
	rcoll = rrow;
	rpivot = i__;
	goto L170;
L130:
	r__ = nr;
	v = 0.;
	index = iflag[2];
	i__2 = index;
	for (kk = 1; kk <= i__2; ++kk) {
	    l1 = i__ - 1 + kk;
	    if (l1 > *n) {
		goto L170;
	    }
	    j = ha[l1 + (ha_dim1 << 3)];
	    r7 = ha[j + (ha_dim1 << 1)];
	    r8 = ha[j + ha_dim1 * 3];
	    r9 = r8 - r7;
	    t = 0.;
	    i__3 = r8;
	    for (k = r7; k <= i__3; ++k) {
		td = (d__1 = a[k], abs(d__1));
/* L140: */
		if (t < td) {
		    t = td;
		}
	    }
	    t /= u;
	    i__3 = r8;
	    for (k = r7; k <= i__3; ++k) {
		td = (d__1 = a[k], abs(d__1));
		if (td < t) {
		    goto L150;
		}
		r6 = snr[k];
		r3 = r9 * (ha[r6 + ha_dim1 * 6] - ha[r6 + ha_dim1 * 5]);
		if (r3 > r__) {
		    goto L150;
		}
		if (r3 < r__) {
		    goto L151;
		}
		if (v >= td) {
		    goto L150;
		}
L151:
		v = td;
		rrow = j;
		rcoll = r6;
		r__ = r3;
		rpivot = l1;
L150:
/* L160: */
		;
	    }
	}
L170:
	r3 = ha[rcoll + ha_dim1 * 10];
	ha[rcoll + ha_dim1 * 10] = ha[i__ + ha_dim1 * 10];
	ha[i__ + ha_dim1 * 10] = r3;
	r3 = ha[rrow + ha_dim1 * 9];
	ha[rrow + ha_dim1 * 9] = ha[i__ + ha_dim1 * 9];

/*  remove the pivot row of the list where the rows are ordered by */
/*  increasing numbers of non-zero elements. */

	ha[i__ + ha_dim1 * 9] = r3;
	l1 = 0;
	l = i__;
	l2 = ha[l4 + ha_dim1 * 3] - ha[l4 + (ha_dim1 << 1)] + 1;
L180:
	++l;
	if (l2 > l1) {
	    ha[l2 + ha_dim1 * 11] = l;
	}
	if (l > *n) {
	    goto L190;
	}
	l5 = ha[l + (ha_dim1 << 3)];
	l3 = ha[l5 + ha_dim1 * 3] - ha[l5 + (ha_dim1 << 1)] + 1;
	if (rpivot < l) {
	    goto L190;
	}
	ha[l4 + ha_dim1 * 7] = l;
	ha[l + (ha_dim1 << 3)] = l4;
	l4 = l5;
	l1 = l2;
	l2 = l3;
	l3 = n8;
	goto L180;
L190:
	if (l2 == l1) {
	    goto L200;
	}
	if (l3 == l2) {
	    goto L200;
	}
	ha[l2 + ha_dim1 * 11] = 0;
L200:
	l5 = ha[i__ + ha_dim1 * 7];
	if (rrow == i__) {
	    goto L210;
	}
	ha[l5 + (ha_dim1 << 3)] = rrow;
	ha[rrow + ha_dim1 * 7] = l5;
L210:
	ha[i__ + ha_dim1 * 7] = rrow;

/*  row interchanges. */

	ha[i__ + (ha_dim1 << 3)] = rcoll;
L220:
	if (rrow == i__) {
	    goto L290;
	}
	t = b[rrow];
	b[rrow] = b[i__];
	b[i__] = t;
	i__3 = rr4;
	for (j = rr3; j <= i__3; ++j) {
	    l1 = snr[j];
	    r__ = ha[l1 + ha_dim1 * 5] - 1;
	    r10 = ha[l1 + ha_dim1 * 6];
L240:
	    ++r__;
	    if (rnr[r__] != i__) {
		goto L240;
	    }
	    rnr[r__] = rnr[r10];
/* L250: */
	    rnr[r10] = rrow;
	}
	rr3 = ha[rrow + (ha_dim1 << 1)];
	rr4 = ha[rrow + ha_dim1 * 3];
	i__3 = rr4;
	for (j = rr3; j <= i__3; ++j) {
	    l1 = snr[j];
	    r__ = ha[l1 + ha_dim1 * 5] - 1;
L260:
	    ++r__;
	    if (rnr[r__] != rrow) {
		goto L260;
	    }
/* L270: */
	    rnr[r__] = i__;
	}
	for (j = 1; j <= 3; ++j) {
	    r3 = ha[rrow + j * ha_dim1];
	    ha[rrow + j * ha_dim1] = ha[i__ + j * ha_dim1];

/*  column interchanges. */

/* L280: */
	    ha[i__ + j * ha_dim1] = r3;
	}
L290:
	if (rcoll == i__) {
	    goto L350;
	}
	i__3 = cr4;
	for (j = c1; j <= i__3; ++j) {
	    l1 = rnr[j];
	    r__ = ha[l1 + (ha_dim1 << 1)] - 1;
	    r10 = ha[l1 + ha_dim1 * 3];
L300:
	    ++r__;
	    if (snr[r__] != i__) {
		goto L300;
	    }
	    t = a[r10];
	    a[r10] = a[r__];
	    a[r__] = t;
	    snr[r__] = snr[r10];
/* L310: */
	    snr[r10] = rcoll;
	}
	c1 = ha[rcoll + (ha_dim1 << 2)];
	cr4 = ha[rcoll + ha_dim1 * 6];
	i__3 = cr4;
	for (j = c1; j <= i__3; ++j) {
	    l1 = rnr[j];
	    r__ = ha[l1 + (ha_dim1 << 1)] - 1;
L320:
	    ++r__;
	    if (snr[r__] != rcoll) {
		goto L320;
	    }
/* L330: */
	    snr[r__] = i__;
	}
	for (j = 4; j <= 6; ++j) {
	    r3 = ha[rcoll + j * ha_dim1];
	    ha[rcoll + j * ha_dim1] = ha[i__ + j * ha_dim1];

/* end of the interchanges. */
/* the row ordered list and the column ordered list are prepared t
o */
/* begin step i of the elimination. */

/* L340: */
	    ha[i__ + j * ha_dim1] = r3;
	}
L350:
	r9 = rr4 - rr3;
	i__3 = rr4;
	for (rr = rr3; rr <= i__3; ++rr) {
	    if (snr[rr] == i__) {
		goto L370;
	    }
/* L360: */
	}
	*ifail = 9;
	goto L1110;
L370:
	v = a[rr];
	pivot[i__] = v;
	td = abs(v);
	if (td < aflag[8]) {
	    aflag[8] = td;
	}
	if (td >= grmin) {
	    goto L380;
	}
	*ifail = 3;
	goto L1110;
L380:
	r2 = ha[i__ + ha_dim1];
	a[rr] = a[rr3];
	snr[rr] = snr[rr3];
	a[rr3] = a[r2];
	snr[rr3] = snr[r2];
	snr[r2] = 0;
	--(*z__);
	++rr3;
	ha[i__ + (ha_dim1 << 1)] = rr3;
	ha[i__ + ha_dim1] = r2 + 1;
	cr3 = ha[i__ + ha_dim1 * 5];
	if (r9 <= 0) {
	    goto L431;
	}
	i__3 = rr4;
	for (j = rr3; j <= i__3; ++j) {
	    index = snr[j];
/* L430: */
	    pivot[index] = a[j];
	}
L431:
	r7 = cr4 - cr3 + 1;
	i__3 = r7;
	for (k = 1; k <= i__3; ++k) {
	    r1 = rnr[cr3 - 1 + k];
	    if (r1 == i__) {
		goto L870;
	    }
	    i1 = ha[r1 + ha_dim1];
	    rr1 = ha[r1 + (ha_dim1 << 1)];
	    rr2 = ha[r1 + ha_dim1 * 3];
	    l2 = rr2 - rr1 + 1;
	    l = rr1 - 1;
L390:
	    ++l;
	    if (snr[l] != i__) {
		goto L390;
	    }
	    t = a[l] / v;
	    if (iflag[5] == 2) {
		goto L400;
	    }
	    a[l] = a[i1];
	    snr[l] = snr[i1];
	    snr[i1] = 0;
	    ++i1;
	    ha[r1 + ha_dim1] = i1;
	    --(*z__);
	    goto L410;
L400:
	    a[l] = a[rr1];
	    a[rr1] = t;
	    r3 = snr[rr1];
	    snr[rr1] = snr[l];
	    snr[l] = r3;
L410:
	    ++rr1;
	    ha[r1 + (ha_dim1 << 1)] = rr1;
	    b[r1] -= b[i__] * t;
	    if (r9 <= 0) {
		goto L669;
	    }
	    r__ = rr1;
	    if (r__ > rr2) {
		goto L470;
	    }
	    i__2 = rr2;
	    for (l = r__; l <= i__2; ++l) {
		l1 = snr[l];
		td = pivot[l1];
		if (td == 0.) {
		    goto L450;
		}
		pivot[l1] = 0.;
		td = a[l] - td * t;
		a[l] = td;
		td1 = abs(td);
		if (td1 > aflag[7]) {
		    aflag[7] = td1;
		}

/*  too small element is created.remove it from the lists. */

		if (td1 > aflag[2]) {
		    goto L450;
		}
		--(*z__);
		a[l] = a[rr1];
		snr[l] = snr[rr1];
		a[rr1] = a[i1];
		snr[rr1] = snr[i1];
		snr[i1] = 0;
		++rr1;
		++i1;
		ha[r1 + (ha_dim1 << 1)] = rr1;
		ha[r1 + ha_dim1] = i1;
		r3 = ha[l1 + ha_dim1 * 5];
		r2 = r3 - 1;
		l4 = ha[l1 + (ha_dim1 << 2)];
		l5 = rnr[l4];
		l6 = rnr[r3];
L440:
		++r2;
		if (rnr[r2] != r1) {
		    goto L440;
		}
		rnr[r2] = l6;
		rnr[r3] = l5;
		rnr[l4] = 0;
		ha[l1 + ha_dim1 * 5] = r3 + 1;
		ha[l1 + (ha_dim1 << 2)] = l4 + 1;
L450:
/* L460: */
		;
	    }
L470:
	    i__2 = r9;
	    for (j = 1; j <= i__2; ++j) {
		r__ = rr3 - 1 + j;
		r2 = snr[r__];
		tol2 = pivot[r2];
		pivot[r2] = a[r__];
		if (tol2 == 0.) {
		    goto L740;
		}
		tol3 = -tol2 * t;
		tol1 = abs(tol3);
		if (tol1 < aflag[2]) {
		    goto L740;
		}
		c2 = ha[r2 + (ha_dim1 << 2)];
		cr2 = ha[r2 + ha_dim1 * 6];
		cr1 = ha[r2 + ha_dim1 * 5];
		lfr = rr2 - i1 + 2;
		lfc = cr2 - c2 + 2;
		if (iflag[4] != 1) {
		    goto L480;
		}
		if (lfr > ha[r1 + ha_dim1 * 9]) {
		    ha[r1 + ha_dim1 * 9] = lfr;
		}
		if (lfc > ha[r2 + ha_dim1 * 10]) {
		    ha[r2 + ha_dim1 * 10] = lfc;
		}
L480:
		if (i1 == 1) {
		    goto L490;
		}
		if (snr[i1 - 1] == 0) {
		    goto L600;
		}
L490:
		if (rr2 == *nn) {
		    goto L500;
		}
		if (snr[rr2 + 1] == 0) {
		    goto L580;
		}

/*  collection in row ordered list. */

L500:
		r10 = *nn - lfr;
		if (r10 >= r4) {
		    goto L560;
		}
		++iflag[6];
		i__4 = *n;
		for (jj = 1; jj <= i__4; ++jj) {
		    l1 = ha[jj + ha_dim1 * 3];
		    if (l1 < ha[jj + ha_dim1]) {
			goto L510;
		    }
		    ha[jj + ha_dim1 * 3] = snr[l1];
		    snr[l1] = -jj;
L510:
/* L520: */
		    ;
		}
		l3 = 0;
		l4 = 1;
		i__4 = r4;
		for (jj = 1; jj <= i__4; ++jj) {
		    if (snr[jj] == 0) {
			goto L540;
		    }
		    ++l3;
		    if (snr[jj] > 0) {
			goto L530;
		    }
		    l5 = -snr[jj];
		    snr[jj] = ha[l5 + ha_dim1 * 3];
		    ha[l5 + ha_dim1 * 3] = l3;
		    l6 = l4 + ha[l5 + (ha_dim1 << 1)] - ha[l5 + ha_dim1];
		    ha[l5 + (ha_dim1 << 1)] = l6;
		    ha[l5 + ha_dim1] = l4;
		    l4 = l3 + 1;
L530:
		    a[l3] = a[jj];
		    snr[l3] = snr[jj];
L540:
/* L550: */
		    ;
		}
		r4 = l3;
		snr[l3 + 1] = 0;
		rr3 = ha[i__ + (ha_dim1 << 1)];
		rr4 = ha[i__ + ha_dim1 * 3];
		i1 = ha[r1 + ha_dim1];
		rr1 = ha[r1 + (ha_dim1 << 1)];
		r__ = rr3 - 1 + j;
		if (r10 >= r4) {
		    goto L560;
		}
		*ifail = 5;

/* fill-in takes place in the row ordered list. */

		goto L1110;
L560:
		r8 = lfr - 1;
		rr2 = r4 + lfr;
		if (r8 <= 0) {
		    goto L579;
		}
		l3 = i1 - 1;
		i__4 = r8;
		for (ll = 1; ll <= i__4; ++ll) {
		    l4 = r4 + ll;
		    l5 = l3 + ll;
		    a[l4] = a[l5];
		    snr[l4] = snr[l5];
/* L570: */
		    snr[l5] = 0;
		}
L579:
		rr1 = r4 + rr1 - i1 + 1;
		ha[r1 + ha_dim1 * 3] = rr2;
		ha[r1 + (ha_dim1 << 1)] = rr1;
		i1 = r4 + 1;
		ha[r1 + ha_dim1] = i1;
		l1 = rr2;
		goto L590;
L580:
		++rr2;
		ha[r1 + ha_dim1 * 3] = rr2;
		l1 = rr2;
		if (rr2 <= r4) {
		    goto L610;
		}
L590:
		r4 = rr2;
		if (r4 < *nn) {
		    snr[r4 + 1] = 0;
		}
		goto L610;
L600:
		--rr1;
		--i1;
		ha[r1 + ha_dim1] = i1;
		ha[r1 + (ha_dim1 << 1)] = rr1;
		l1 = rr1;
		snr[i1] = snr[l1];
		a[i1] = a[l1];
L610:
		a[l1] = tol3;
		snr[l1] = snr[r__];
		td = (d__1 = a[l1], abs(d__1));
		if (td > aflag[7]) {
		    aflag[7] = td;
		}
		++(*z__);
		if (iflag[8] < *z__) {
		    iflag[8] = *z__;
		}
		if (c2 == 1) {
		    goto L620;
		}
		if (rnr[c2 - 1] == 0) {
		    goto L720;
		}
L620:
		if (cr2 == *nn1) {
		    goto L630;
		}
		if (rnr[cr2 + 1] == 0) {
		    goto L700;
		}

/*  collection in column ordered list. */

L630:
		r10 = *nn1 - lfc;
		if (r10 >= r5) {
		    goto L680;
		}
		++iflag[7];
		i__4 = *n;
		for (jj = i__; jj <= i__4; ++jj) {
		    l1 = ha[jj + ha_dim1 * 6];
		    ha[jj + ha_dim1 * 6] = rnr[l1];
/* L640: */
		    rnr[l1] = -jj;
		}
		l3 = 0;
		l4 = 1;
		i__4 = r5;
		for (jj = 1; jj <= i__4; ++jj) {
		    if (rnr[jj] == 0) {
			goto L660;
		    }
		    ++l3;
		    if (rnr[jj] > 0) {
			goto L650;
		    }
		    l5 = -rnr[jj];
		    rnr[jj] = ha[l5 + ha_dim1 * 6];
		    ha[l5 + ha_dim1 * 6] = l3;
		    l6 = l4 + ha[l5 + ha_dim1 * 5] - ha[l5 + (ha_dim1 << 2)];
		    ha[l5 + ha_dim1 * 5] = l6;
		    ha[l5 + (ha_dim1 << 2)] = l4;
		    l4 = l3 + 1;
L650:
		    rnr[l3] = rnr[jj];
L660:
/* L670: */
		    ;
		}
		r5 = l3;
		rnr[r5 + 1] = 0;
		c2 = ha[r2 + (ha_dim1 << 2)];
		cr3 = ha[i__ + ha_dim1 * 5];
		cr4 = ha[i__ + ha_dim1 * 6];
		cr1 = ha[r2 + ha_dim1 * 5];
		if (r10 >= r5) {
		    goto L680;
		}
		*ifail = 6;

/* fill-in takes place in the column ordered list. */

		goto L1110;
L680:
		r8 = lfc - 1;
		cr2 = r5 + lfc;
		if (r8 <= 0) {
		    goto L699;
		}
		l3 = c2 - 1;
		i__4 = r8;
		for (l = 1; l <= i__4; ++l) {
		    l4 = r5 + l;
		    l5 = l3 + l;
		    rnr[l4] = rnr[l5];
/* L690: */
		    rnr[l5] = 0;
		}
L699:
		cr1 = r5 + cr1 - c2 + 1;
		c2 = r5 + 1;
		ha[r2 + ha_dim1 * 6] = cr2;
		ha[r2 + (ha_dim1 << 2)] = c2;
		ha[r2 + ha_dim1 * 5] = cr1;
		r__ = cr2;
		goto L710;
L700:
		++cr2;
		ha[r2 + ha_dim1 * 6] = cr2;
		r__ = cr2;
		if (cr2 <= r5) {
		    goto L730;
		}
L710:
		r5 = cr2;
		if (r5 < *nn1) {
		    rnr[r5 + 1] = 0;
		}
		goto L730;
L720:
		--cr1;
		--c2;
		ha[r2 + (ha_dim1 << 2)] = c2;
		ha[r2 + ha_dim1 * 5] = cr1;
		r__ = cr1;
		rnr[c2] = rnr[r__];
L730:
		rnr[r__] = r1;
L740:
/* L750: */
		;
	    }
L669:
	    if (rr1 <= rr2) {
		goto L760;
	    }
	    *ifail = 7;

/*  update the information in the list where the rows are ordered 
by */
/*  increasing numbers of the non-zero elements. */

	    goto L1110;
L760:
	    if (iflag[4] == 2) {
		goto L870;
	    }
	    if (iflag[3] == 0) {
		goto L870;
	    }
	    l1 = rr2 - rr1 + 1;
	    if (l1 == l2) {
		goto L870;
	    }
	    l6 = ha[r1 + ha_dim1 * 7];
	    l4 = ha[l2 + ha_dim1 * 11];
	    if (l1 > l2) {
		goto L820;
	    }
	    if (l6 > l4) {
		goto L780;
	    }
	    if (l4 == *n) {
		goto L770;
	    }
	    l = ha[l4 + 1 + (ha_dim1 << 3)];
	    l5 = ha[l + ha_dim1 * 3] - ha[l + (ha_dim1 << 1)] + 1;
	    if (l5 == l2) {
		goto L790;
	    }
L770:
	    ha[l2 + ha_dim1 * 11] = 0;
	    goto L800;
L780:
	    l5 = ha[l4 + (ha_dim1 << 3)];
	    l3 = ha[l6 + (ha_dim1 << 3)];
	    ha[l4 + (ha_dim1 << 3)] = l3;
	    ha[l6 + (ha_dim1 << 3)] = l5;
	    ha[l5 + ha_dim1 * 7] = l6;
	    ha[l3 + ha_dim1 * 7] = l4;
	    l6 = l4;
L790:
	    ha[l2 + ha_dim1 * 11] = l4 + 1;
L800:
	    if (l4 == i__ + 1) {
		goto L810;
	    }
	    l = ha[l6 - 1 + (ha_dim1 << 3)];
	    l2 = ha[l + ha_dim1 * 3] - ha[l + (ha_dim1 << 1)] + 1;
	    l4 = ha[l2 + ha_dim1 * 11];
	    if (l1 < l2) {
		goto L780;
	    }
L810:
	    if (l1 != l2) {
		ha[l1 + ha_dim1 * 11] = l6;
	    }
	    goto L870;
L820:
	    if (l6 > l4) {
		goto L840;
	    }
	    if (l4 == *n) {
		goto L830;
	    }
	    l = ha[l4 + 1 + (ha_dim1 << 3)];
	    l5 = ha[l + ha_dim1 * 3] - ha[l + (ha_dim1 << 1)] + 1;
	    if (l5 == l2) {
		goto L840;
	    }
L830:
	    ha[l2 + ha_dim1 * 11] = 0;
L840:
	    ++l2;
	    if (l2 <= slut) {
		goto L850;
	    }
	    l3 = *n;
	    slut = l1;
	    l2 = l1;
	    goto L860;
L850:
	    l3 = ha[l2 + ha_dim1 * 11] - 1;
	    if (l3 == -1) {
		goto L840;
	    }
	    if (l2 > l1) {
		l2 = l1;
	    }
L860:
	    ha[l2 + ha_dim1 * 11] = l3;
	    l4 = ha[l3 + (ha_dim1 << 3)];
	    l7 = ha[l6 + (ha_dim1 << 3)];
	    ha[l3 + (ha_dim1 << 3)] = l7;
	    ha[l6 + (ha_dim1 << 3)] = l4;
	    ha[l7 + ha_dim1 * 7] = l3;
	    ha[l4 + ha_dim1 * 7] = l6;
	    l6 = l3;
	    if (l2 < l1) {
		goto L840;
	    }
L870:
/* L880: */
	    ;
	}
	if (r9 <= 0) {
	    goto L882;
	}
	i__3 = rr4;
	for (j = rr3; j <= i__3; ++j) {
	    index = snr[j];
/* L881: */
	    pivot[index] = 0.;
	}
L882:
	cr3 = ha[i__ + (ha_dim1 << 2)];
	i__3 = cr4;
	for (j = cr3; j <= i__3; ++j) {
/* L890: */
	    rnr[j] = 0;
	}
	if (r9 <= 0) {
	    goto L930;
	}
	l2 = ha[i__ + (ha_dim1 << 1)] - 1;
	i__3 = r9;
	for (ll = 1; ll <= i__3; ++ll) {
	    r__ = snr[l2 + ll];
	    r1 = ha[r__ + ha_dim1 * 5];
	    r2 = ha[r__ + ha_dim1 * 6];
	    if (r2 > r1) {
		goto L900;
	    }
	    *ifail = 8;
	    goto L1110;
L900:
	    ha[r__ + ha_dim1 * 5] = r1 + 1;
	    r3 = r1 - 1;
L910:
	    ++r3;
	    if (rnr[r3] != i__) {
		goto L910;
	    }
	    rnr[r3] = rnr[r1];
/* L920: */
	    rnr[r1] = i__;
	}
L930:
	aflag[5] = aflag[7] / aflag[6];
	if (aflag[5] < aflag[3]) {
	    goto L940;
	}
	*ifail = 4;
	goto L1110;
L940:

/*  preparation to begin the back substitution. */

/* L950: */
	;
    }
    index = ha[*n + (ha_dim1 << 1)];
    pivot[*n] = a[index];
    a[index] = 0.;
    td = (d__1 = pivot[*n], abs(d__1));
    if (td > aflag[7]) {
	aflag[7] = td;
    }
    if (td < aflag[8]) {
	aflag[8] = td;
    }
    if (td > grmin) {
	goto L960;
    }
    *ifail = 3;
    goto L1110;
L960:
    if (iflag[4] != 1) {
	goto L1060;
    }
    iflag[10] = ha[*n + ha_dim1 * 9];
    iflag[9] = ha[*n + ha_dim1 * 10];
    i__1 = n7;
    for (i__ = 1; i__ <= i__1; ++i__) {
	r1 = *n - i__;
	iflag[10] += ha[r1 + ha_dim1 * 9];
	iflag[9] += ha[r1 + ha_dim1 * 10];
	if (iflag[3] == 0) {
	    goto L980;
	}
	for (j = 9; j <= 10; ++j) {
	    r2 = ha[r1 + (j - 2) * ha_dim1];
	    r6 = ha[r2 + j * ha_dim1];
	    ha[r2 + j * ha_dim1] = ha[r1 + j * ha_dim1];
/* L970: */
	    ha[r1 + j * ha_dim1] = r6;
	}
L980:
/* L990: */
	;
    }
L1060:
    aflag[5] = aflag[7] / aflag[6];
    iflag[1] = -2;
L1110:
    *z__ = zz;
    return 0;
} /* y12mcf_ */

#ifdef __MWERKS__
#pragma warn_unusedarg reset
#endif /* __MWERKS__ */
