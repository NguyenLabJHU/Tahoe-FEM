/* PAK (04/27/98) */

/* y12maf.f -- translated by f2c (version 19960611).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "y12m.h"

#ifdef __MWERKS__
#pragma warn_unusedarg off
#endif /* __MWERKS__ */

/* Subroutine */ int y12maf(n, z__, a, snr, nn, rnr, nn1, pivot, ha, iha, 
	aflag, iflag, b, ifail)
integer *n, *z__;
doublereal *a;
integer *snr, *nn, *rnr, *nn1;
doublereal *pivot;
integer *ha, *iha;
doublereal *aflag;
integer *iflag;
doublereal *b;
integer *ifail;
{
    /* System generated locals */
    integer ha_dim1, ha_offset;

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
    aflag[1] = 16.;
    aflag[2] = 1e-12;
    aflag[3] = 1e16;
    aflag[4] = 1e-12;
    iflag[2] = 2;
    iflag[3] = 1;
    iflag[4] = 0;
    iflag[5] = 1;
    y12mbf(n, z__, &a[1], &snr[1], nn, &rnr[1], nn1, &ha[ha_offset], iha, &
	    aflag[1], &iflag[1], ifail);
    if (*ifail != 0) {
	goto L1;
    }
    y12mcf(n, z__, &a[1], &snr[1], nn, &rnr[1], nn1, &pivot[1], &b[1], &ha[
	    ha_offset], iha, &aflag[1], &iflag[1], ifail);
    if (*ifail != 0) {
	goto L1;
    }
    y12mdf(n, &a[1], nn, &b[1], &pivot[1], &snr[1], &ha[ha_offset], iha, &
	    iflag[1], ifail);
L1:
    return 0;
} /* y12maf_ */

#ifdef __MWERKS__
#pragma warn_unusedarg reset
#endif /* __MWERKS__ */
