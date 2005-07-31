/* PAK (07/24/98) */

/* y12mhe.f -- translated by f2c (version 19960611).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "y12m.h"

#ifdef __MWERKS__
#pragma warn_unusedarg off
#endif /* __MWERKS__ */

/* Subroutine */ int y12mhe(n, nz, a, snr, work, anorm)
integer *n, *nz;
real *a;
integer *snr;
real *work, *anorm;
{
    /* System generated locals */
    integer i__1;
    real r__1;

    /* Local variables */
    static integer i__, l;



/*   purpose. */
/*   ------- */


/*   this subroutine computes the    one-norm   of a sparse */
/*   matrix   a.   all parameters  (except    anorm )   have  the */
/*   same meaning as in the other subroutines in  package   y12m. */
/*   on exit the   one-norm   of matrix   a   will be stored in */
/*   anorm. */




/*  declaration of the global variables and arrays. */




/*  declaration of the internal variables. */




/*  set all locations of array     work     equal to zero. */


    /* Parameter adjustments */
    --work;
    --snr;
    --a;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	work[i__] = (float)0.;
    }


/*  calculate the sums of the absolute values of the non-zero */
/*  elements in each row of matrix     a .     store these sums */
/*  in array     work . */


    i__1 = *nz;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l = snr[i__];
/* L20: */
	work[l] += (r__1 = a[i__], dabs(r__1));
    }


/*  calculate the one-norm of matrix     a . */


    *anorm = (float)0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L30: */
	if (work[i__] > *anorm) {
	    *anorm = work[i__];
	}
    }


/*  stop the computations. */


    return 0;
} /* y12mhe_ */

#ifdef __MWERKS__
#pragma warn_unusedarg reset
#endif /* __MWERKS__ */
