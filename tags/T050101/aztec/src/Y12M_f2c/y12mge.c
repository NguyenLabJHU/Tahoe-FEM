/* PAK (04/27/98) */

/* y12mge.f -- translated by f2c (version 19960611).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "y12m.h"

#ifdef __MWERKS__
#pragma warn_unusedarg off
#endif /* __MWERKS__ */

/* Subroutine */ int y12mge(n, nn, a, snr, w, pivot, anorm, rcond, iha, ha, 
	iflag, ifail)
integer *n, *nn;
real *a;
integer *snr;
real *w, *pivot, *anorm, *rcond;
integer *iha, *ha, *iflag, *ifail;
{
    /* System generated locals */
    integer ha_dim1, ha_offset, i__1, i__2;
    real r__1;

    /* Local variables */
    static integer i__, j, l;
    static real t;
    static integer l1, l2, l3;
    static real ynorm, znorm;
    static integer n7, n8;



/*   purpose. */
/*   -------- */

/*   this subroutine computes the number called    rcond    by */
/*   dongarra et al.(1979). this number is the reciprocal of the */
/*   condition number of matrix   a .   the subroutine can be */
/*   called immediately after the call of   y12mce.   the parameters */
/*   (except   rcond   and   anorm ) have the same meaning as the */
/*   corresponding parameters in the subroutines of package   y12m  (the 
*/
/*   subroutine can be  call only if the   lu   decomposition of matrix */
/*   a   computed by   y12mce   is not destroyed). subroutine  y12mhe */
/*   should be called before the call of subroutine   y12mce (this */
/*   subroutine calculates the   one-norm   of matrix   a   and */
/*   stores  it  in   anorm.   on successful exit   rcond   will */
/*   contain an approximation to the reciprocal of the condition */
/*   number of matrix   a.   more details, concerning the use */
/*   of similar subroutines for  dense matrices, can be found */
/*   in   j.dongarra, j.r.bunch, c.b.moler and g.w.stewart (1979): */
/*        "linpack - user's guide", siam, philadelphia. */




/*  declaration of the global variables and arrays. */




/*  declaration of the internal variables. */




/*   check whether the entry is correct or not. */


    /* Parameter adjustments */
    --pivot;
    --w;
    --snr;
    --a;
    ha_dim1 = *iha;
    ha_offset = ha_dim1 + 1;
    ha -= ha_offset;
    --iflag;

    /* Function Body */
    if (*ifail != 0) {
	goto L180;
    }
    if (iflag[5] != 1) {
	goto L10;
    }
    *ifail = 26;
    goto L180;


/*   no error detected.  the computations will be continued. */


L10:
    n8 = *n + 1;
    n7 = *n - 1;


/*   solve a system of the form    u1*w=e   where   u1   is the */
/*   transpose of matrix   u   in the   lu-factorization of matrix */
/*   a   and    e    is a vector whose components are equal to   +1 */
/*   or   -1. */


    w[1] = (float)1. / pivot[1];
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
/* L20: */
	w[i__] = (float)0.;
    }
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	l1 = ha[i__ + (ha_dim1 << 1)];
	l2 = ha[i__ + ha_dim1 * 3];
	if (l1 > l2) {
	    goto L40;
	}
	t = w[i__ - 1];
	i__2 = l2;
	for (j = l1; j <= i__2; ++j) {
	    l = snr[j];
/* L30: */
	    w[l] += t * a[j];
	}
L40:
	if (w[i__] > (float)0.) {
	    w[i__] += (float)1.;
	}
	if (w[i__] <= (float)0.) {
	    w[i__] += (float)-1.;
	}
/* L50: */
	w[i__] = -w[i__] / pivot[i__];
    }


/*   solve a system of the form   l1*y=w   where   l1   is the */
/*   transpose of matrix   l   in the   lu-factorization of */
/*   matrix   a .   the components of vector   y   are stored */
/*   array   w  (thus, the contents of array   w   are overwritten */
/*   by the components of vector   y ). */


    i__1 = n7;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l = *n - i__;
	l1 = ha[l + ha_dim1];
	l2 = ha[l + (ha_dim1 << 1)] - 1;
	if (l1 > l2) {
	    goto L70;
	}
	t = w[l + 1];
	i__2 = l2;
	for (j = l1; j <= i__2; ++j) {
	    l3 = snr[j];
/* L60: */
	    w[l3] -= t * a[j];
	}
L70:
/* L80: */
	;
    }


/*   calculate the one-norm of vector   y . */


    ynorm = (float)0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L90: */
	ynorm += (r__1 = w[i__], dabs(r__1));
    }


/*   compute the solution of    (lu)z=y .  this means that */
/*   two systems with triangular matrices are solved using the */
/*   same ideas as above. the components of the calculated solution */
/*   are stored in array   w . */


    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l1 = ha[i__ + ha_dim1];
	l2 = ha[i__ + (ha_dim1 << 1)] - 1;
	if (l1 > l2) {
	    goto L120;
	}
	i__2 = l2;
	for (j = l1; j <= i__2; ++j) {
	    l = snr[j];
/* L110: */
	    w[i__] -= a[j] * w[l];
	}
L120:
/* L130: */
	;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l3 = n8 - i__;
	l1 = ha[l3 + (ha_dim1 << 1)];
	l2 = ha[l3 + ha_dim1 * 3];
	if (l1 > l2) {
	    goto L150;
	}
	i__2 = l2;
	for (j = l1; j <= i__2; ++j) {
	    l = snr[j];
/* L140: */
	    w[l3] -= a[j] * w[l];
	}
L150:
/* L160: */
	w[l3] /= pivot[l3];
    }


/*   compute the one-norm of vector   z   (vector   z   is */
/*   the vector calculated above and stored in array   w . */


    znorm = (float)0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L170: */
	znorm += (r__1 = w[i__], dabs(r__1));
    }


/*   find the value of the required estimate for the reciprocal */
/*   of the condition number of matrix   a . */


    *rcond = ynorm / *anorm / znorm;


/*   end of the computations. */


L180:
    if (*ifail != 0) {
	*rcond = (float)-1.;
    }
    return 0;
} /* y12mge_ */

#ifdef __MWERKS__
#pragma warn_unusedarg reset
#endif /* __MWERKS__ */
