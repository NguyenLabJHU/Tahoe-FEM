/* PAK (07/24/98) */

/* y12mbf.f -- translated by f2c (version 19960611).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "y12m.h"

#ifdef __MWERKS__
#pragma warn_unusedarg off
#endif /* __MWERKS__ */

/* Subroutine */ int y12mbf(n, z__, a, snr, nn, rnr, nn1, ha, iha, aflag, 
	iflag, ifail)
integer *n, *z__;
doublereal *a;
integer *snr, *nn, *rnr, *nn1, *ha, *iha;
doublereal *aflag;
integer *iflag, *ifail;
{
    /* System generated locals */
    integer ha_dim1, ha_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer mode, i__, j, r__;
    static doublereal t;
    static integer index, l1, l2, l3, l4, l5;
    static doublereal gt1;



/*  the non-zero elements of a sparse matrix a are prepared  in order to 
*/
/*  solve the system ax=b by use of sparse matrix technique/ */


    /* Parameter adjustments */
    --snr;
    --a;
    --rnr;
    ha_dim1 = *iha;
    ha_offset = ha_dim1 + 1;
    ha -= ha_offset;
    --aflag;
    --iflag;

    /* Function Body */
    mode = iflag[4];
    *ifail = 0;
    if (*n < 2) {
	*ifail = 12;
    }
    if (*z__ <= 0) {
	*ifail = 13;
    }
    if (*nn < *z__ << 1) {
	*ifail = 5;
    }
    if (*nn1 < *z__) {
	*ifail = 6;
    }
    if (*ifail == 0 && *n > *z__) {
	*ifail = 14;
    }
    if (*iha < *n) {
	*ifail = 15;
    }
    if (mode < 0) {
	*ifail = 16;
    }
    if (mode > 2) {
	*ifail = 16;
    }
    if (*ifail != 0) {
	goto L22;
    }
    gt1 = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ha[i__ + (ha_dim1 << 1)] = 0;
	ha[i__ + ha_dim1 * 3] = 0;
/* L10: */
	ha[i__ + ha_dim1 * 6] = 0;
    }

/*  find the number of the non-zero elements in each row and column;move 
*/
/*  the non-zero elements in the end of the arrays a and snr;find the */
/*  largest non-zero element in a(in absolute value). */

    i__1 = *z__;
    for (i__ = 1; i__ <= i__1; ++i__) {
	t = (d__1 = a[i__], abs(d__1));
	l3 = rnr[i__];
	l4 = snr[i__];
	if (l4 > *n || l4 < 1) {
	    *ifail = 24;
	}
	if (l3 > *n || l3 < 1) {
	    *ifail = 25;
	}
	++ha[l3 + ha_dim1 * 3];
	++ha[l4 + ha_dim1 * 6];
	if (t > gt1) {
	    gt1 = t;
	}
	a[*z__ + i__] = a[i__];
/* L20: */
	snr[*z__ + i__] = snr[i__];
    }
    if (*ifail > 0) {
	goto L22;
    }

/*  store the information of the row starts(in ha(i,1))and of the column 
*/
/*  starts(in ha(i,4)). */

    l1 = 1;
    l2 = 1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l3 = ha[i__ + ha_dim1 * 3];
	l4 = ha[i__ + ha_dim1 * 6];
	if (l3 > 0) {
	    goto L21;
	}
	*ifail = 17;
	goto L22;
L21:
	if (l4 > 0) {
	    goto L23;
	}
	*ifail = 18;
	goto L22;
L23:
	if (mode == 2) {
	    goto L30;
	}
	ha[i__ + ha_dim1 * 9] = l3;
	ha[i__ + ha_dim1 * 10] = l4;
	ha[i__ + ha_dim1 * 11] = 0;
	++ha[l3 + (ha_dim1 << 1)];
	ha[i__ + ha_dim1 * 5] = l3;
L30:
	ha[i__ + ha_dim1] = l1;
	ha[i__ + (ha_dim1 << 2)] = l2;
	l1 += l3;
	l2 += l4;
	ha[i__ + ha_dim1 * 3] = 0;
/* L40: */
	ha[i__ + ha_dim1 * 6] = 0;
    }

/*  store the non-zero elements of matrix a(ordered in rows) in the */
/*  first z locations of the array a.do the same for their column numbers 
*/

    i__1 = *z__;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l1 = *z__ + i__;
	l3 = rnr[i__];
	l2 = ha[l3 + ha_dim1] + ha[l3 + ha_dim1 * 3];
	a[l2] = a[l1];
	snr[l2] = snr[l1];
/* L50: */
	++ha[l3 + ha_dim1 * 3];
    }

/*  store the row numbers of the non-zero elements ordered by columns in 
*/
/*  the first z locations of the array rnr. store information about row */
/*  ends(in ha(i,3)). */

    l4 = 1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (mode == 2) {
	    goto L60;
	}
	if (ha[i__ + (ha_dim1 << 1)] == 0) {
	    goto L60;
	}
	ha[i__ + ha_dim1 * 11] = l4;
	l4 += ha[i__ + (ha_dim1 << 1)];
	ha[i__ + (ha_dim1 << 1)] = ha[i__ + ha_dim1 * 11];
L60:
	ha[i__ + ha_dim1 * 3] = ha[i__ + ha_dim1] + ha[i__ + ha_dim1 * 3] - 1;
	l1 = ha[i__ + ha_dim1];
	l2 = ha[i__ + ha_dim1 * 3];
	i__2 = l2;
	for (j = l1; j <= i__2; ++j) {
	    l3 = snr[j];
	    r__ = ha[l3 + ha_dim1 * 6];
	    index = ha[l3 + (ha_dim1 << 2)] + r__;
	    rnr[index] = i__;
	    if (r__ == 0) {
		goto L70;
	    }
	    if (j == l1) {
		goto L70;
	    }
	    if (rnr[index - 1] != i__) {
		goto L70;
	    }
	    *ifail = 11;
	    goto L22;
L70:
	    ha[l3 + ha_dim1 * 6] = r__ + 1;
	}
    }
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	if (mode == 2) {
	    goto L80;
	}
	l3 = ha[i__ + ha_dim1 * 5];
	l5 = ha[l3 + (ha_dim1 << 1)];
	ha[l5 + (ha_dim1 << 3)] = i__;
	ha[i__ + ha_dim1 * 7] = l5;
	++ha[l3 + (ha_dim1 << 1)];
L80:
/* L90: */
	ha[i__ + ha_dim1 * 6] = ha[i__ + (ha_dim1 << 2)] + ha[i__ + ha_dim1 * 
		6] - 1;
    }
    aflag[6] = gt1;
    iflag[6] = 0;
    iflag[7] = 0;
    iflag[8] = *z__;
    iflag[1] = -1;
L22:
    return 0;
} /* y12mbf_ */

#ifdef __MWERKS__
#pragma warn_unusedarg reset
#endif /* __MWERKS__ */
