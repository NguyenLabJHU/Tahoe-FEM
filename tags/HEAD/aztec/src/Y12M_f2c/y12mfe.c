/* PAK (04/27/98) */

/* y12mfe.f -- translated by f2c (version 19960611).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "y12m.h"

#ifdef __MWERKS__
#pragma warn_unusedarg off
#endif /* __MWERKS__ */

/* Subroutine */ int y12mfe(n, a, snr, nn, rnr, nn1, a1, sn, nz, ha, iha, b, 
	b1, x, y, aflag, iflag, ifail)
integer *n;
real *a;
integer *snr, *nn, *rnr, *nn1;
real *a1;
integer *sn, *nz, *ha, *iha;
real *b, *b1, *x, *y, *aflag;
integer *iflag, *ifail;
{
    /* System generated locals */
    integer ha_dim1, ha_offset, i__1, i__2;
    real r__1;

    /* Local variables */
    static real dcor, dres;
    static integer nres;
    static real d__;
    static integer i__, j;
    static integer state, l1, l2, l3;
    static real dd;
    static doublereal er;
    static integer it;
    static real xm, xx;
    static doublereal er1, er2;
    static real gt1, gt2;
    static integer kit;


/*  store the non-zero elements,their column numbers,information about */
/*  row starts,information about row ends and the right-hand side. */

    /* Parameter adjustments */
    --y;
    --x;
    --b1;
    --b;
    --snr;
    --a;
    --rnr;
    --sn;
    --a1;
    ha_dim1 = *iha;
    ha_offset = ha_dim1 + 1;
    ha -= ha_offset;
    --aflag;
    --iflag;

    /* Function Body */
    *ifail = 0;
    nres = 0;
    dres = (float)0.;
    state = iflag[5];
    kit = 1;
    it = iflag[11];
    if (state == 1) {
	*ifail = 10;
    }
    if (it < 2) {
	*ifail = 23;
    }
    if (*ifail != 0) {
	goto L160;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	b1[i__] = b[i__];
    }
    if (state == 3) {
	goto L70;
    }
    y12mbe(n, nz, &a[1], &snr[1], nn, &rnr[1], nn1, &ha[ha_offset], iha, &
	    aflag[1], &iflag[1], ifail);
    if (*ifail != 0) {
	goto L160;
    }
    i__1 = *nz;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sn[i__] = snr[i__];
/* L20: */
	a1[i__] = a[i__];
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ha[i__ + ha_dim1 * 12] = ha[i__ + ha_dim1];
/* L30: */
	ha[i__ + ha_dim1 * 13] = ha[i__ + ha_dim1 * 3];
    }
    if (aflag[2] >= (float)0.) {
	goto L60;
    }
    gt1 = aflag[6];
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l1 = ha[i__ + ha_dim1];
	l2 = ha[i__ + ha_dim1 * 3];
	gt2 = (float)0.;
	i__2 = l2;
	for (j = l1; j <= i__2; ++j) {
	    d__ = (r__1 = a[j], dabs(r__1));
/* L40: */
	    if (gt2 < d__) {
		gt2 = d__;
	    }
	}
/* L50: */
	if (gt2 < gt1) {
	    gt1 = gt2;
	}
    }
    aflag[2] = -gt1 * aflag[2];

/*  find the first solution. */

L60:
    y12mce(n, nz, &a[1], &snr[1], nn, &rnr[1], nn1, &y[1], &b[1], &ha[
	    ha_offset], iha, &aflag[1], &iflag[1], ifail);
    if (*ifail != 0) {
	goto L160;
    }
L70:
    y12mde(n, &a[1], nn, &b[1], &y[1], &snr[1], &ha[ha_offset], iha, &iflag[
	    1], ifail);
    if (*ifail != 0) {
	goto L160;
    }

/*  prepare the data in order to begin the iterations. */

    dd = (float)0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = b[i__];
	xx = (r__1 = b[i__], dabs(r__1));
/* L80: */
	if (dd < xx) {
	    dd = xx;
	}
    }
    xm = dd;
    if (dd == (float)0.) {
	goto L160;
    }

/*  begin to iterate. */

L90:
    d__ = dd;
    dres = (float)0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	er = b1[i__];
	l1 = ha[i__ + ha_dim1 * 12];
	l2 = ha[i__ + ha_dim1 * 13];
	i__2 = l2;
	for (j = l1; j <= i__2; ++j) {
	    er1 = a1[j];
	    l3 = sn[j];
	    er2 = x[l3];
/* L100: */
	    er -= er1 * er2;
	}

/*  store residuals rounded to single precision. */

	b[i__] = er;
	xx = abs(er);
/* L110: */
	if (dres < xx) {
	    dres = xx;
	}
    }
    if (dres == (float)0.) {
	goto L160;
    }
    if (nres == 1) {
	goto L150;
    }
    if (dres > xm * (float)1e4) {
	goto L150;
    }
    ++kit;
    iflag[5] = 3;
    y12mde(n, &a[1], nn, &b[1], &y[1], &snr[1], &ha[ha_offset], iha, &iflag[
	    1], ifail);
    if (*ifail != 0) {
	goto L160;
    }

/*  compute the uniform norm of the current solution vector. */

    dd = (float)0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xx = (r__1 = b[i__], dabs(r__1));
/* L120: */
	if (dd < xx) {
	    dd = xx;
	}
    }
    if (dd == (float)0.) {
	goto L160;
    }

/*  check the convergence criterion. */

    if (dd > d__ && kit > 2) {
	goto L160;
    }

/*  calculate an improved solution. */

    dcor = dd;
    xm = (float)0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] += b[i__];
	xx = (r__1 = x[i__], dabs(r__1));
/* L130: */
	if (xx > xm) {
	    xm = xx;
	}
    }

/*  check the stopping criteria. */

    if (dd / xm + (float)10. == (float)10.) {
	goto L140;
    }
    if (kit < it) {
	goto L90;
    }

/*  end of the iterations. */

L140:
    nres = 1;
    goto L90;
L150:
    dd = dabs(dd);
L160:
    iflag[5] = state;
    iflag[12] = kit;
    aflag[9] = dd;
    aflag[10] = dres;
    aflag[11] = xm;
    return 0;
} /* y12mfe_ */

#ifdef __MWERKS__
#pragma warn_unusedarg reset
#endif /* __MWERKS__ */
