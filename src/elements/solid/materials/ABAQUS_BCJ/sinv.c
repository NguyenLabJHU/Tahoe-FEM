/* $Id: sinv.c,v 1.3 2004-08-01 20:42:35 paklein Exp $ */ 
/* sinv.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

/* library support options */
#ifdef __F2C__

#include "f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Subroutine */ int sinv_(doublereal *stress, doublereal *sinv1, doublereal *
	sinv2, integer *ndi, integer *nshr)
{
    /* Initialized data */

    static doublereal zero = 0.;
    static doublereal deux = 2.;
    static doublereal trois = 3.;
    static doublereal tdemi = 1.5;

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal stresd[6];

/* ====================================================================== */
/* ---------------------------------------------------------------------- */
/* ------------  CALCUL DU PREMIER ET DEUXIEME INVARIANTS  -------------- */
/* ----------------  D'UN TENSEUR D'ORDRE 2 SYMETRIQUE  ----------------- */
/* ---------------------------------------------------------------------- */
/* ====================================================================== */
/* ------------- */
/*  EN ENTREE : */
/* ------------- */
/*  STRESS	: TENSEURS D'ORDRE DEUX */
/*  NDI		: NOMBRE DE COMPOSANTES DES CONTRAINTES DIRECTES */
/*  NSHR		: NOMBRE DE COMPOSANTES DES CONTRAINTES TANGENTIELLES */
/* ------------- */
/*  EN SORTIE : */
/* ------------- */
/*  SINV1	: PREMIER INVARIANT */
/*  SINV2	: DEUXIEME INVARIANT */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
/* .1-----  Implicit, External */
/* .2-----  Parameter */
/* .3-----  Dimension */
/* .4-----  Real,Integer,Complex,Double precision,Logical,Character */
/* .5-----  Common */
/* .6-----  Equivalence */
/* .7-----  Data */
    /* Parameter adjustments */
    --stress;

    /* Function Body */
/* .8-----  Definition de fonctions */
/* ======================================================================= */
    *sinv1 = zero;
    i__1 = *ndi;
    for (i__ = 1; i__ <= i__1; ++i__) {
	*sinv1 += stress[i__];
    }
    *sinv1 /= trois;
/* ------- */
    i__1 = *ndi;
    for (i__ = 1; i__ <= i__1; ++i__) {
	stresd[i__ - 1] = stress[i__] - *sinv1;
    }
    i__1 = *nshr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	stresd[i__ + *ndi - 1] = stress[i__ + *ndi];
    }
/* ------- */
    *sinv2 = zero;
    i__1 = *ndi;
    for (i__ = 1; i__ <= i__1; ++i__) {
	*sinv2 += stresd[i__ - 1] * stresd[i__ - 1];
    }
    i__1 = *nshr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	*sinv2 += deux * stresd[i__ + *ndi - 1] * stresd[i__ + *ndi - 1];
    }
    *sinv2 = sqrt(tdemi * *sinv2);
/* ======================================================================= */
    return 0;
} /* sinv_ */

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* __F2C__ */
