/* $Id: rotsig.c,v 1.3 2004-08-01 20:42:35 paklein Exp $ */ 
/* rotsig.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

/* library support options */
#ifdef __F2C__

#include "f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Subroutine */ int rotsig_(doublereal *s, doublereal *r__, doublereal *
	sprime, integer *lstr, integer *ndi, integer *nshr)
{
    /* Initialized data */

    static doublereal udemi = .5;
    static doublereal un = 1.;
    static doublereal deux = 2.;

    static doublereal sh, sh2, ssh3, ssh4, ssh5, ssh6;

/* ======================================================================= */
/* ----------------------------------------------------------------------- */
/* -------  ROTATION D'UN TENSEUR AVEC LA MATRICE ORIENTATION   ---------- */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
/* ----------- */
/*  ENTREES : */
/* ----------- */
/*  S		: TENSEUR NON TOURNE */
/*  R		: MATRICE INCREMENT DE ROTATION */
/*  LSTR		: SI LSTR=1 ---> S EST UN TENSEUR DES CONTRAINTES */
/* 		  SI LSTR=2 ---> S EST UN TENSEUR DES DEFORMATIONS */
/*  NDI		: NOMBRE DE COMPOSANTES DIRECTES */
/*  NSHR		: NOMBRE DE COMPOSANTES TANGEBTIELLES */
/* ----------- */
/*  SORTIES : */
/* ----------- */
/*  SPRIME	: TENSEUR TOURNE */
/* ----------------------------------------------------------------------* */
/* ======================================================================= */
/* ----------------------------------------------------------------------* */
/* .1-----  Implicit, External */
/* .2-----  Parameter */
/* .3-----  Dimension */
/* .4-----  Real,Integer,Complex,Double precision,Logical,Character */
/* .5-----  Common */
/* .6-----  Equivalence */
/* .7-----  Save,Data */
    /* Parameter adjustments */
    r__ -= 4;
    --sprime;
    --s;

    /* Function Body */
/* .8-----  Definition de fonctions */
/* ----------------------------------------------------------------------* */
/* ======================================================================= */
/* ----------------------------------------------------------------------* */
    if (*lstr == 1) {
	sh = un;
	sh2 = un;
    }
    if (*lstr == 2) {
	sh = udemi;
	sh2 = deux;
    }
/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------* */
    if (*ndi == 3 && *nshr == 1) {
/* ----------------------------------------------------------------------* */
/* -------  Etat Deformation Plane */

	ssh4 = sh * s[4];
	sprime[1] = r__[4] * (r__[4] * s[1] + r__[7] * ssh4) + r__[7] * (r__[
		4] * ssh4 + r__[7] * s[2]) + r__[10] * r__[10] * s[3];
	sprime[2] = r__[5] * (r__[5] * s[1] + r__[8] * ssh4) + r__[8] * (r__[
		5] * ssh4 + r__[8] * s[2]) + r__[11] * r__[11] * s[3];
	sprime[3] = r__[6] * (r__[6] * s[1] + r__[9] * ssh4) + r__[9] * (r__[
		6] * ssh4 + r__[9] * s[2]) + r__[12] * r__[12] * s[3];
	sprime[4] = r__[4] * (r__[5] * s[1] + r__[8] * ssh4) + r__[7] * (r__[
		5] * ssh4 + r__[8] * s[2]) + r__[10] * r__[11] * s[3];
	sprime[4] = sh2 * sprime[4];
/* ----------------------------------------------------------------------* */
    } else if (*ndi == 3 && *nshr == 3) {
/* ----------------------------------------------------------------------* */
/* -------  Etat 3D */

	ssh4 = sh * s[4];
	ssh5 = sh * s[5];
	ssh6 = sh * s[6];
	sprime[1] = r__[4] * (r__[4] * s[1] + r__[7] * ssh4 + r__[10] * ssh5) 
		+ r__[7] * (r__[4] * ssh4 + r__[7] * s[2] + r__[10] * ssh6) + 
		r__[10] * (r__[4] * ssh5 + r__[7] * ssh6 + r__[10] * s[3]);
	sprime[2] = r__[5] * (r__[5] * s[1] + r__[8] * ssh4 + r__[11] * ssh5) 
		+ r__[8] * (r__[5] * ssh4 + r__[8] * s[2] + r__[11] * ssh6) + 
		r__[11] * (r__[5] * ssh5 + r__[8] * ssh6 + r__[11] * s[3]);
	sprime[3] = r__[6] * (r__[6] * s[1] + r__[9] * ssh4 + r__[12] * ssh5) 
		+ r__[9] * (r__[6] * ssh4 + r__[9] * s[2] + r__[12] * ssh6) + 
		r__[12] * (r__[6] * ssh5 + r__[9] * ssh6 + r__[12] * s[3]);
	sprime[4] = r__[4] * (r__[5] * s[1] + r__[8] * ssh4 + r__[11] * ssh5) 
		+ r__[7] * (r__[5] * ssh4 + r__[8] * s[2] + r__[11] * ssh6) + 
		r__[10] * (r__[5] * ssh5 + r__[8] * ssh6 + r__[11] * s[3]);
	sprime[5] = r__[4] * (r__[6] * s[1] + r__[9] * ssh4 + r__[12] * ssh5) 
		+ r__[7] * (r__[6] * ssh4 + r__[9] * s[2] + r__[12] * ssh6) + 
		r__[10] * (r__[6] * ssh5 + r__[9] * ssh6 + r__[12] * s[3]);
	sprime[6] = r__[5] * (r__[6] * s[1] + r__[9] * ssh4 + r__[12] * ssh5) 
		+ r__[8] * (r__[6] * ssh4 + r__[9] * s[2] + r__[12] * ssh6) + 
		r__[11] * (r__[6] * ssh5 + r__[9] * ssh6 + r__[12] * s[3]);
	sprime[4] = sh2 * sprime[4];
	sprime[5] = sh2 * sprime[5];
	sprime[6] = sh2 * sprime[6];
/* ----------------------------------------------------------------------* */
    } else {
/* ----------------------------------------------------------------------* */
/* -------  Etat Contrainte Plane */

	ssh3 = sh * s[3];
	sprime[1] = r__[4] * (r__[4] * s[1] + r__[7] * ssh3) + r__[7] * (r__[
		4] * ssh3 + r__[7] * s[2]);
	sprime[2] = r__[5] * (r__[5] * s[1] + r__[8] * ssh3) + r__[8] * (r__[
		5] * ssh3 + r__[8] * s[2]);
	sprime[3] = r__[4] * (r__[5] * s[1] + r__[8] * ssh3) + r__[7] * (r__[
		5] * ssh3 + r__[8] * s[2]);
	sprime[3] = sh2 * sprime[3];
/* ----------------------------------------------------------------------* */
    }
/* ----------------------------------------------------------------------* */
/* ======================================================================= */
/* ----------------------------------------------------------------------* */
    return 0;
} /* rotsig_ */

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* __F2C__ */
