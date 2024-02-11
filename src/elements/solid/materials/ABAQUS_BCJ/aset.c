/* $Id: aset.c,v 1.3 2004-08-01 20:42:35 paklein Exp $ */ 
/* aset.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

/* library support options */
#ifdef __F2C__

#include "f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Subroutine */ int aset_(doublereal *t, doublereal *cste, integer *n)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

/* ====================================================================== */
/* ---------------------------------------------------------------------- */
/* ---------  EXECUTION DE L'OPERATION TENSORIELLE: T = CSTE  ----------- */
/* ---------------------------------------------------------------------- */
/* ====================================================================== */
/* ----------- */
/*  ENTREES : */
/* ----------- */
/*  N 		: DIMENSION DU TENSEUR T */
/*  CSTE 	: CONSTANTE A AFFECTER AU TENSEUR T */
/* ----------- */
/*  SORTIES : */
/* ----------- */
/*  T		: TENSEUR CONSTANT */
/* ---------------------------------------------------------------------- */
/* ====================================================================== */
/* .1-----  Implicit, External */
/* .2-----  Parameter */
/* .3-----  Dimension */
/* .4-----  Real,Integer,Complex,Double precision,Logical,Character */
/* .5-----  Common */
/* .6-----  Equivalence */
/* .7-----  Save,Data */
/* .8-----  Definition de fonctions */
/* ======================================================================= */
    /* Parameter adjustments */
    --t;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	t[i__] = *cste;
    }
/* ======================================================================= */
    return 0;
} /* aset_ */

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* __F2C__ */
