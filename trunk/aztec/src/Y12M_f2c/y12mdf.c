/* PAK (07/24/98 */

/* y12mdf.f -- translated by f2c (version 19960611).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "y12m.h"

#ifdef __MWERKS__
#pragma warn_unusedarg off
#endif /* __MWERKS__ */

/* Subroutine */ int y12mdf(n, a, nn, b, pivot, snr, ha, iha, iflag, ifail)
integer *n;
doublereal *a;
integer *nn;
doublereal *b, *pivot;
integer *snr, *ha, *iha, *iflag, *ifail;
{
    /* System generated locals */
    integer ha_dim1, ha_offset, i__1, i__2;

    /* Local variables */
    static integer mode, ipiv, i__, j;
    static doublereal t;
    static integer state, l1, r1, r2, n7, n8, rr1, rr2;

    /* Parameter adjustments */
    --pivot;
    --b;
    --snr;
    --a;
    ha_dim1 = *iha;
    ha_offset = ha_dim1 + 1;
    ha -= ha_offset;
    --iflag;

    /* Function Body */
    *ifail = 0;
    if (iflag[1] == -2) {
	goto L1000;
    }
    *ifail = 1;
    goto L1110;
L1000:
    mode = iflag[4];
    ipiv = iflag[3];
    n8 = *n + 1;
    n7 = *n - 1;
    state = iflag[5];

/*  solve the system with lower triangular matrix  l  (if the */
/*  lu-factorization is available). */

    if (state != 3) {
	goto L1051;
    }
    if (ipiv == 0) {
	goto L1020;
    }
    i__1 = n7;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l1 = ha[i__ + ha_dim1 * 7];
	t = b[l1];
	b[l1] = b[i__];
	b[i__] = t;
/* L1010: */
    }
L1020:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rr1 = ha[i__ + ha_dim1];
	rr2 = ha[i__ + (ha_dim1 << 1)] - 1;
	if (rr1 > rr2) {
	    goto L1040;
	}
	i__2 = rr2;
	for (j = rr1; j <= i__2; ++j) {
	    l1 = snr[j];
/* L1030: */
	    b[i__] -= a[j] * b[l1];
	}
L1040:
/* L1050: */
	;
    }

/*  solve the system with upper triagular matrix. */

L1051:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	r1 = n8 - i__;
	rr1 = ha[r1 + (ha_dim1 << 1)];
	rr2 = ha[r1 + ha_dim1 * 3];
	if (rr2 < rr1) {
	    goto L1080;
	}
	i__2 = rr2;
	for (j = rr1; j <= i__2; ++j) {
	    r2 = snr[j];
/* L1070: */
	    b[r1] -= a[j] * b[r2];
	}
L1080:
/* L1090: */
	b[r1] /= pivot[r1];
    }

/* if interchanges were used during the  elimination then a reordering in 
*/
/* lution vector is made. */

    if (ipiv == 0) {
	goto L1110;
    }
    i__1 = n7;
    for (i__ = 1; i__ <= i__1; ++i__) {
	r1 = *n - i__;
	r2 = ha[r1 + (ha_dim1 << 3)];
	t = b[r2];
	b[r2] = b[r1];
/* L1100: */
	b[r1] = t;
    }
L1110:
    return 0;
} /* y12mdf_ */

#ifdef __MWERKS__
#pragma warn_unusedarg reset
#endif /* __MWERKS__ */
