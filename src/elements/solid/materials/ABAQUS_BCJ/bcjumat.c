/* bcjumat.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

/* library support options */
#ifdef __F2C__

#include "f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Table of constant values */
static integer c__9 = 9;
static integer c__1 = 1;

/* ************************************************************************* */

/*              BAMMANN - CHIESA - JOHNSON  Plasticity Model */

/* ************************************************************************* */
/* ************************************************************************* */
/* **                                                                     ** */
/* **   WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING     ** */
/* **                                                                     ** */
/* **                    this is proprietary software                     ** */
/* **                                                                     ** */
/* **        recipients of this routine are not allowed to distribute     ** */
/* **        their source to other facilities.  rather, requests for      ** */
/* **        the source should be directed to the author at sandia.       ** */
/* **                                                                     ** */
/* **                        Michael L. Chiesa                            ** */
/* **                       Sandia National Labs                          ** */
/* **                         Mail Stop 9042                              ** */
/* **                            PO Box 969                               ** */
/* **                     Livermore, CA 94551-0969                        ** */
/* **                                                                     ** */
/* **                   email:    chiesa@sandia.gov                       ** */
/* **                                                                     ** */
/* **   WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING     ** */
/* **                                                                     ** */
/* ************************************************************************* */
/* ************************************************************************* */

/*                                Copyright */
/*                         Sandia Corporation 1995 */
/*                           All Rights Reserved */

/*     Notice: The Government is granted for itself and others acting on */
/*  its behalf a paid-up, nonexclusive, irrevocable worldwide license in */
/*  this data to reproduce, prepare derivative works, and perform publicly */
/*  and display publicly.  Beginning five (5) years after permission to */
/*  assert copyright was obtained, the Government is granted for itself */
/*  and others acting on its behalf a paid-up, nonexclusive, irrevocable */
/*  worldwide license in this data to reproduce, prepare derivative works, */
/*  distribute copies to the public, perform publicly and display publicly, */
/*  and to permit others to do so.  NEITHER THE UNITED STATES NOR THE */
/*  UNITED STATES DEPARTMENT OF ENERGY, NOR ANY OF THEIR EMPLOYEES, MAKES */
/*  ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR */
/*  RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY */
/*  INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS */
/*  THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

/* ************************************************************************* */

/* Subroutine */ int bcjumat_(doublereal *stress, doublereal *statev, doublereal 
	*ddsdde, doublereal *sse, doublereal *spd, doublereal *scd, 
	doublereal *rpl, doublereal *ddsddt, doublereal *drplde, doublereal *
	drpldt, doublereal *stran, doublereal *dstran, doublereal *time, 
	doublereal *dtime, doublereal *temp, doublereal *dtemp, doublereal *
	predef, doublereal *dpred, char *cmname, integer *ndi, integer *nshr, 
	integer *ntens, integer *nstatv, doublereal *props, integer *nprops, 
	doublereal *coords, doublereal *drot, doublereal *pnewdt, doublereal *
	celent, doublereal *dfgrd0, doublereal *dfgrd1, integer *noel, 
	integer *npt, integer *layer, integer *kspt, integer *kstep, integer *
	kinc, ftnlen cmname_len)
{
    /* Initialized data */

    static doublereal con1 = .816496580927726;
    static char name1[] = "bcj2.2  ";
    static char name2[] = "BCJ2.2  ";
    static integer iparam1 = 0;
    static integer iparam2 = 0;
    static integer iparam3 = 0;

    /* System generated locals */
    integer ddsdde_dim1, ddsdde_offset, i__1, i__2;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_wsle(cilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_wsle(void);
    /* Subroutine */ int s_stop(char *, ftnlen);
    double exp(doublereal), tanh(doublereal), sqrt(doublereal), log(
	    doublereal), sinh(doublereal), pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static doublereal g;
    static integer i__, j;
    static doublereal p, r__, h1, h2, z1, z2, z3, ak, c90, xi[6], cc1, cc2, 
	    cc3, cc4, cc5, cc6, cc7, cc8, cc9, ak2, rd1, rd2, rs1, rs2, bbb, 
	    cc10, cc11, cc12, cc13, cc14, cc15, cc16, cc17, cc18, cc19, cc20, 
	    cc21, ddd, blk, arg, phi, sig, dum, psi, tmp, sto, yyy, phi0, 
	    phi1, phi2, exp4, sto2;
    static integer iadd;
    static doublereal dgam, beta, davg, dsig, alpm, htcp, pold, expo, twog, 
	    dgam2, alpm2, temp0, term1, term2, term3, term4, term5, term6, 
	    term7, term8, term9, gamma, dalph, theta, ximag, tempi, slope, 
	    ximag2, ftheta, trialk, vtheta, ytheta, ztheta, epsdot, tratio;

    /* Fortran I/O blocks */
    static cilist io___7 = { 0, 6, 0, 0, 0 };
    static cilist io___8 = { 0, 6, 0, 0, 0 };
    static cilist io___9 = { 0, 6, 0, 0, 0 };
    static cilist io___10 = { 0, 6, 0, 0, 0 };
    static cilist io___16 = { 0, 6, 0, 0, 0 };
    static cilist io___77 = { 0, 6, 0, 0, 0 };
    static cilist io___104 = { 0, 6, 0, 0, 0 };





    /* Parameter adjustments */
    --time;
    --predef;
    --dpred;
    --dstran;
    --stran;
    --drplde;
    --ddsddt;
    ddsdde_dim1 = *ntens;
    ddsdde_offset = 1 + ddsdde_dim1;
    ddsdde -= ddsdde_offset;
    --stress;
    --statev;
    --props;
    --coords;
    drot -= 4;
    dfgrd0 -= 4;
    dfgrd1 -= 4;

    /* Function Body */

/* ************************************************************************* */

/*        Sandia strain rate and temperature dependent plasticity model */
/*          with damage/failure */
/*        Version 2.1 */
/*        Implemented by M.L. Chiesa - 07/21/97 for abaqus */
/*        This subroutine can be used for axisymmetric, plane strain and */
/*        3d solid elements. */

/* ************************************************************************* */

/*     *  ntens = number of non-zero stress components (4 for 2d, 6 for 3d) */
/*     *  nshr  = number of non-zero shear components (1 for 2d, 3 for 3d) */
/*     *  ndi   = number of non-zero normal stresses (always 3) */
/*     *  nstatv = number of state variables (11) */
/*     *  nprops = number of material parameters (30) */

/* ************************************************************************* */

/*     *  statev(1) = alpha-xx */
/*     *  statev(2) = alpha-yy */
/*     *  statev(3) = alpha-zz */
/*     *  statev(4) = alpha-xy */
/*     *  statev(5) = alpha-yz */
/*     *  statev(6) = alpha-zx */
/*     *  statev(7) = kappa */
/*     *  statev(8) = temperature */
/*     *  statev(9) = effective plastic strain */
/*     *  statev(10)= damage */
/*     *  statev(11)= damage growth rate */

/* ************************************************************************* */

/*     *  props(1) thru props(5) are constants for Johnson/Bammann */
/*     *  formulas for shear and bulk moduli */
/*     *  props(1) = mu zero         , props(2) = a */
/*     *  props(3) = K zero          , props(4) = b */
/*     *  props(5) = T melt          , props(6) = C1 */
/*     *  props(7) = C2              , props(8) = C3 */
/*     *  props(9) = C4              , props(10)= C5 */
/*     *  props(11)= C6              , props(12)= C7 */
/*     *  props(13)= C8              , props(14)= C9 */
/*     *  props(15)= C10             , props(16)= C11 */
/*     *  props(17)= C12             , props(18)= C13 */
/*     *  props(19)= C14             , props(20)= C15 */
/*     *  props(21)= C16             , props(22)= C17 */
/*     *  props(23)= C18             , props(24)= C19 */
/*     *  props(25)= C20             , props(26)= C21 */
/*     *  props(27)= initial temperature */
/*     *  props(28)= heat generation coefficient */
/*     *  props(29)= initial damage */
/*     *  props(30)= damage growth exponent */

/* ************************************************************************ */

/*   * iparam1 = 0 for linear return (uses total strain in recovery) */
/*   *         = 1 for quadratic return (uses plastic strain in recovery) */
/*   *           (only linear return implemented in this version) */

/*   * iparam2 = 0 for Simo tangent stiffness matrix */
/*   *         = 1 for Lathrop tangent stiffness matrix */
/*   *           (only Simo stiffness implemented in this version) */

/*   * iparam3 = 0 for trial kappa = kappa(n) */
/*   *         = 1 for trial kappa = kappa(n + 1/2) */

/* ************************************************************************ */

/* ---- check for correct umat version */

    if (0 && s_cmp(cmname, name1, (ftnlen)8, (ftnlen)8) != 0 && s_cmp(cmname, 
	    name2, (ftnlen)8, (ftnlen)8) != 0) {
	s_wsle(&io___7);
	do_lio(&c__9, &c__1, "   ********************************************"
		, (ftnlen)47);
	e_wsle();
	s_wsle(&io___8);
	do_lio(&c__9, &c__1, "   ** E R R O R - this is umat bcj2.2        **"
		, (ftnlen)47);
	e_wsle();
	s_wsle(&io___9);
	do_lio(&c__9, &c__1, "   ** your parameters may not be consistent  **"
		, (ftnlen)47);
	e_wsle();
	s_wsle(&io___10);
	do_lio(&c__9, &c__1, "   ********************************************"
		, (ftnlen)47);
	e_wsle();
	s_stop("", (ftnlen)0);
    }

/* ---- set model parameters */

    tempi = props[27];
    htcp = props[28];
    phi0 = props[29];
    expo = props[30];

/* ---- initialize state variables */
/*     temperature is set first by *initial condition command in abaqus */
/*     or by props(xx) */

    if (time[2] == 0.f) {
	statev[1] = 0.f;
	statev[2] = 0.f;
	statev[3] = 0.f;
	statev[4] = 0.f;
	statev[5] = 0.f;
	statev[6] = 0.f;
	statev[7] = 0.f;
	statev[9] = 0.f;
	statev[10] = phi0;
	statev[11] = 0.f;
	temp0 = tempi;
	if (*temp == 0.f) {
	    if (temp0 == 0.f) {
		s_wsle(&io___16);
		do_lio(&c__9, &c__1, " error - temperature is zero", (ftnlen)
			28);
		e_wsle();
		s_stop("", (ftnlen)0);
	    } else {
		statev[8] = temp0;
	    }
	} else {
	    statev[8] = *temp;
	}
    }

/* ****************************************************************** */

/* ---- set parameters */
/*     iadd is the number of parameters before the model constants */

    iadd = 5;
    cc1 = props[iadd + 1];
    cc2 = props[iadd + 2];
    cc3 = props[iadd + 3];
    cc4 = props[iadd + 4];
    cc5 = props[iadd + 5];
    cc6 = props[iadd + 6];
    cc7 = props[iadd + 7];
    cc8 = props[iadd + 8];
    cc9 = props[iadd + 9];
    cc10 = props[iadd + 10];
    cc11 = props[iadd + 11];
    cc12 = props[iadd + 12];
    cc13 = props[iadd + 13];
    cc14 = props[iadd + 14];
    cc15 = props[iadd + 15];
    cc16 = props[iadd + 16];
    cc17 = props[iadd + 17];
    cc18 = props[iadd + 18];
    cc19 = props[iadd + 19];
    cc20 = props[iadd + 20];
    cc21 = props[iadd + 21];

/* *********************************************************************** */

/* ---- determine elastic moduli */
/*     g = shear modulus   twog = 2*g   blk = bulk modulus */

    if (tempi == 0.f) {
	statev[8] = *temp;
    }
    theta = statev[8];
    if (props[5] == 0.f) {
	blk = props[3];
	g = props[1];
    } else {
	tratio = theta / props[5];
	tratio = min(tratio,.9999);
	blk = props[3] - props[4] * tratio;
	g = props[1] * (1.f - tratio * exp(props[2] * (1.f - 1.f / tratio)));
    }

    twog = g * 2.f;

/* ---- constants used in damage calculations */

    if (statev[10] == 0.f) {
	phi = 0.f;
	phi1 = 1.f;
	phi2 = 1.f;
    } else {
	phi = statev[10];
	phi1 = 1.f - phi;
/* Computing MIN */
	d__1 = 1.f, d__2 = *dtime * statev[11] / phi1;
	phi2 = 1.f - min(d__1,d__2);
    }

/* ---- calculate pressure */

    davg = (dstran[1] + dstran[2] + dstran[3]) * .33333333333333331;
    pold = (stress[1] + stress[2] + stress[3]) * .33333333333333331;
    p = pold * phi2 + blk * 3 * davg * phi1;

/* ---- check for melt */
/*     if material has melted set all deviatoric stresses and */
/*     state variables to zero. compressive hydrostatic pressure */
/*     is allowed. */

    if (props[5] != 0.f) {
	if (theta > props[5]) {
	    i__1 = *ntens;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		stress[i__] = 0.f;
/* L35: */
		statev[i__] = 0.f;
	    }
	    p = min(0.,p);
	    stress[1] = p;
	    stress[2] = p;
	    stress[3] = p;
	    statev[7] = 0.f;
	    statev[9] = 0.f;
	    statev[10] = phi0;
	    statev[11] = 0.f;
	    goto L315;
	}
    }

/* ************************************************************************* */

/* ---- compute function evaluations */
/*       theta = temperature */
/*       ytheta = static yield stress */
/*       vtheta,ftheta = functions to define rate dependence of yield */
/*       h1,h2 = plastic hardeing moduli */
/*       rs1,rs2 = static recovery functions */
/*       rd1,rd2 = dynamic recovery functions */

    ztheta = 1.f / theta;
    yyy = (tanh(cc19 * (cc20 - theta)) + 1.f) * .5f;
    if (cc19 == 0.f) {
	yyy = 1.;
    }
    exp4 = exp(cc4 * ztheta);
    ytheta = cc3 * exp4 * yyy / (cc21 * exp4 + 1.f);
    vtheta = cc1 * exp(-cc2 * ztheta);
    ftheta = cc5 * exp(-cc6 * ztheta);
    rd1 = cc7 * exp(-cc8 * ztheta);
/* Computing MAX */
    d__1 = cc9 - cc10 * theta, d__2 = cc9 * 1e-5f;
    h1 = max(d__1,d__2);
    rs1 = cc11 * exp(-cc12 * ztheta);
    rd2 = cc13 * exp(-cc14 * ztheta);
/* Computing MAX */
    d__1 = cc15 - cc16 * theta, d__2 = cc15 * 1e-5f;
    h2 = max(d__1,d__2);
    rs2 = cc17 * exp(-cc18 * ztheta);

/* *********************************************************************** */

/* ---- update alpha using abaqus rotation matrix */

    if (rs1 != 0.f || rd1 != 0.f || h1 != 0.f) {
	if (*ntens == 4) {
	    term1 = drot[4] * statev[1] + drot[7] * statev[4];
	    term2 = drot[4] * statev[4] + drot[7] * statev[2];
	    term3 = drot[5] * statev[1] + drot[8] * statev[4];
	    term4 = drot[5] * statev[4] + drot[8] * statev[2];
	    term5 = statev[1] + statev[2];
/* Computing 2nd power */
	    d__1 = drot[10];
	    statev[1] = drot[4] * term1 + drot[7] * term2 - term5 * (d__1 * 
		    d__1);
/* Computing 2nd power */
	    d__1 = drot[11];
	    statev[2] = drot[5] * term3 + drot[8] * term4 - term5 * (d__1 * 
		    d__1);
	    statev[3] = -(statev[1] + statev[2]);
	    statev[4] = drot[5] * term1 + drot[8] * term2 - term5 * drot[10] *
		     drot[11];
	} else {
	    term1 = drot[4] * statev[1] + drot[7] * statev[4] + drot[10] * 
		    statev[6];
	    term2 = drot[4] * statev[4] + drot[7] * statev[2] + drot[10] * 
		    statev[5];
	    term3 = drot[4] * statev[6] + drot[7] * statev[5] + drot[10] * 
		    statev[3];
	    term4 = drot[5] * statev[1] + drot[8] * statev[4] + drot[11] * 
		    statev[6];
	    term5 = drot[5] * statev[4] + drot[8] * statev[2] + drot[11] * 
		    statev[5];
	    term6 = drot[5] * statev[6] + drot[8] * statev[5] + drot[11] * 
		    statev[3];
	    term7 = drot[6] * statev[1] + drot[9] * statev[4] + drot[12] * 
		    statev[6];
	    term8 = drot[6] * statev[4] + drot[9] * statev[2] + drot[12] * 
		    statev[5];
	    term9 = drot[6] * statev[6] + drot[9] * statev[5] + drot[12] * 
		    statev[3];
	    statev[1] = term1 * drot[4] + term2 * drot[7] + term3 * drot[10];
	    statev[2] = term4 * drot[5] + term5 * drot[8] + term6 * drot[11];
	    statev[3] = term7 * drot[6] + term8 * drot[9] + term9 * drot[12];
	    statev[4] = term1 * drot[5] + term2 * drot[8] + term3 * drot[11];
	    statev[5] = term4 * drot[6] + term5 * drot[9] + term6 * drot[12];
	    statev[6] = term1 * drot[6] + term2 * drot[9] + term3 * drot[12];
	}
    }

/* ---- compute effective strain rate */

    if (*dtime != 0.f) {
/* Computing 2nd power */
	d__1 = dstran[4];
/* Computing 2nd power */
	d__2 = dstran[1] - dstran[2];
/* Computing 2nd power */
	d__3 = dstran[2] - dstran[3];
/* Computing 2nd power */
	d__4 = dstran[3] - dstran[1];
	dum = d__1 * d__1 * .33333333333333331 + (d__2 * d__2 + d__3 * d__3 + 
		d__4 * d__4) * .22222222222222221;
	if (*ntens == 6) {
/* Computing 2nd power */
	    d__1 = dstran[5];
/* Computing 2nd power */
	    d__2 = dstran[6];
	    dum += (d__1 * d__1 + d__2 * d__2) * .33333333333333331;
	}
	ddd = sqrt(dum) / *dtime;
    } else {
	ddd = 0.f;
    }

/* ---- calculate trial alpha, kappa and yield radius */

/* Computing 2nd power */
    d__1 = statev[1];
/* Computing 2nd power */
    d__2 = statev[2];
/* Computing 2nd power */
    d__3 = statev[3];
/* Computing 2nd power */
    d__4 = statev[4];
/* Computing 2nd power */
    d__5 = statev[5];
/* Computing 2nd power */
    d__6 = statev[6];
    alpm = con1 * sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + (d__4 * d__4 
	    + d__5 * d__5 + d__6 * d__6) * 2.f);
    sto = *dtime * rs1 * alpm;
    sto2 = *dtime * rs2;
    if (iparam1 == 0) {
	sto += *dtime * rd1 * ddd * alpm;
	sto2 += *dtime * rd2 * ddd;
    } else {
	s_wsle(&io___77);
	do_lio(&c__9, &c__1, " this option is not implemented", (ftnlen)31);
	e_wsle();
	s_stop("", (ftnlen)0);
    }
    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L50: */
	statev[i__] *= 1.f - sto;
    }
    if (iparam3 == 0) {
	trialk = statev[7];
    } else {
	trialk = (sqrt(sto2 * 2 * (statev[7] + h2 * .5f * ddd * *dtime) + 1.f)
		 - 1.f) / max(1e-30,sto2);
    }
    statev[7] -= sto2 * trialk * trialk;
/* Computing 2nd power */
    d__1 = ddd;
/* Computing 2nd power */
    d__2 = ftheta;
    ak = (vtheta * log((ddd + sqrt(d__1 * d__1 + d__2 * d__2)) / ftheta) + 
	    ytheta + statev[7]) * phi1;

/* ---- calculate trial elastic deviatoric stresses */

    for (i__ = 1; i__ <= 3; ++i__) {
/* L60: */
	stress[i__] = (stress[i__] - pold) * phi2 + phi1 * twog * (dstran[i__]
		 - davg);
    }
    i__1 = *ntens;
    for (i__ = 4; i__ <= i__1; ++i__) {
/* L70: */
	stress[i__] = stress[i__] * phi2 + phi1 * g * dstran[i__];
    }

/* ---- compute xi (deviatoric stress - 2/3 alpha) */

    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L80: */
	xi[i__ - 1] = stress[i__] - statev[i__] * .66666666666666663;
    }

/* ---- compute (magnitude of xi) squared */

/* Computing 2nd power */
    d__1 = xi[0];
/* Computing 2nd power */
    d__2 = xi[1];
/* Computing 2nd power */
    d__3 = xi[2];
/* Computing 2nd power */
    d__4 = xi[3];
    ximag2 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + d__4 * d__4 * 2.f;
    if (*ntens == 6) {
/* Computing 2nd power */
	d__1 = xi[4];
/* Computing 2nd power */
	d__2 = xi[5];
	ximag2 += (d__1 * d__1 + d__2 * d__2) * 2;
    }

/* ---- check for plasticity */

    ak2 = ximag2 - ak * .66666666666666663 * abs(ak);
    if (ak2 <= 0.f || ddd == 0.f) {
	goto L300;
    }

/* ---- plasticity process begins here */

    ximag = sqrt(ximag2);

/* ---- return trial stresses to yield surface, add pressure term */
/*      and update state variables */

    if (iparam1 == 0) {
	dgam = (ximag - con1 * ak) / (phi1 * twog + (h1 + h2 * phi1) * 
		.66666666666666663);
    }
    dgam2 = dgam / ximag;

    dsig = phi1 * twog * dgam2;
    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L90: */
	stress[i__] -= dsig * xi[i__ - 1];
    }
    stress[1] += p;
    stress[2] += p;
    stress[3] += p;

    statev[7] += dgam * con1 * h2;
    statev[7] = max(0.,statev[7]);

    dalph = h1 * dgam2;
    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L100: */
	statev[i__] += dalph * xi[i__ - 1];
    }

/* ---- update plastic strain (for output purposes only) */

    statev[9] += dgam * con1;

/* ---- update temperature for adiabatic problems */

    if (htcp > 0.f) {
	dum = 0.f;
	i__1 = *ntens;
	for (i__ = 4; i__ <= i__1; ++i__) {
/* L110: */
	    dum += stress[i__] * xi[i__ - 1];
	}
	statev[8] += htcp * dgam2 * (stress[1] * xi[0] + stress[2] * xi[1] + 
		stress[3] * xi[2] + dum * 2.f);
	theta = statev[8];
/*     if(props(5).gt.0.)then */
/*     tratio = theta/props(5) */
/*     tratio = min(tratio,0.9999) */
/*     g = props(1)*(1.-tratio*exp(props(2)*(1.-1./tratio))) */
/*     blk = props(3) - props(4)*tratio */
/*     twog = 2.0 * g */
/*     endif */
/*     h1 = max((1.e-5*cc9),(cc9-cc10*statev(8))) */
/*     h2 = max((1.e-5*cc15),(cc15-cc16*statev(8))) */
/*     rd1=cc7*exp(-cc8/statev(8)) */
/*     rd2=cc13*exp(-cc14/statev(8)) */
    }

/* ---- update damage */

    if (phi != 0.f && *dtime != 0.f) {
	epsdot = dgam * con1 / *dtime;
/* Computing 2nd power */
	d__1 = stress[1] - stress[2];
/* Computing 2nd power */
	d__2 = stress[2] - stress[3];
/* Computing 2nd power */
	d__3 = stress[3] - stress[1];
/* Computing 2nd power */
	d__4 = stress[4];
	dum = (d__1 * d__1 + d__2 * d__2 + d__3 * d__3) * .33333333333333331 
		+ d__4 * d__4 * 2;
	if (*ntens == 6) {
/* Computing 2nd power */
	    d__1 = stress[5];
/* Computing 2nd power */
	    d__2 = stress[6];
	    dum += (d__1 * d__1 + d__2 * d__2) * 2;
	}
/* Computing MAX */
	d__1 = 1e-15f, d__2 = sqrt(dum);
	sig = max(d__1,d__2);
/* Computing MIN */
	d__1 = 15.f, d__2 = p * 2 * (expo - .5f) / ((expo + .5f) * sig);
	arg = min(d__1,d__2);
	beta = sinh((max(0.,arg)));
	c90 = expo + 1.f;
/* Computing MIN */
	d__1 = 15.f, d__2 = beta * *dtime * epsdot * c90;
	psi = min(d__1,d__2);
/* Computing MAX */
	d__1 = 0.f, d__2 = (pow_dd(&phi1, &c90) - 1.f) * exp(psi) + 1.f;
	tmp = max(d__1,d__2);
/* Computing MIN */
	d__2 = 1.f / c90;
	d__1 = 1.f - pow_dd(&tmp, &d__2);
	phi = min(d__1,.999);
	d__1 = 1 - phi;
	statev[11] = beta * epsdot * (1.f / pow_dd(&d__1, &expo) - (1 - phi));
	statev[10] = phi;
	phi1 = 1.f - statev[10];
    }
/*     write(*,*)ak,phi,ximag,statev(7),statev(9) */

/* ---- form elastic-plastic constitutive matrix */
/*     if bbb=1 then conventional stiffness is calculated */
/*     otherwise simo stiffness is calculated */

    if (iparam2 == 0) {
	dum = 0.f;
	i__1 = *ntens;
	for (i__ = 4; i__ <= i__1; ++i__) {
/* L120: */
/* Computing 2nd power */
	    d__1 = statev[i__];
	    dum += d__1 * d__1;
	}
/* Computing 2nd power */
	d__1 = statev[1];
/* Computing 2nd power */
	d__2 = statev[2];
/* Computing 2nd power */
	d__3 = statev[3];
	alpm2 = con1 * sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + dum * 
		2.f);
/*     bbb = 1.0 */
	bbb = (con1 * (ak + h2 * con1 * dgam) + alpm2 - alpm) / ximag;
	bbb = min(1.,bbb);
/*     gamma = 1.0 / (1.0 + (h1 + h2)/3./g) - (1.0-bbb) */
	slope = h1 + h2 - rd1 * alpm2 * alpm2 - rd2 * statev[7] * statev[7];
/* Computing MAX */
	d__1 = slope, d__2 = (h1 + h2) * 1e-5f;
	slope = max(d__1,d__2);
	gamma = 1.f / (slope / 3.f / g + 1.f) - (1.f - bbb);
/* Computing 2nd power */
	d__1 = ximag;
	r__ = phi1 * (twog * gamma / (d__1 * d__1));
	z1 = phi1 * (blk + twog * .66666666666666663 * bbb);
	z2 = phi1 * (blk - twog * .33333333333333331 * bbb);
	z3 = phi1 * g * bbb;
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *ntens;
	    for (j = 1; j <= i__2; ++j) {
/* L140: */
		ddsdde[i__ + j * ddsdde_dim1] = -r__ * xi[i__ - 1] * xi[j - 1]
			;
	    }
	}
	ddsdde[ddsdde_dim1 + 1] = z1 + ddsdde[ddsdde_dim1 + 1];
	ddsdde[(ddsdde_dim1 << 1) + 2] = z1 + ddsdde[(ddsdde_dim1 << 1) + 2];
	ddsdde[ddsdde_dim1 * 3 + 3] = z1 + ddsdde[ddsdde_dim1 * 3 + 3];
	ddsdde[(ddsdde_dim1 << 1) + 1] = z2 + ddsdde[(ddsdde_dim1 << 1) + 1];
	ddsdde[ddsdde_dim1 * 3 + 1] = z2 + ddsdde[ddsdde_dim1 * 3 + 1];
	ddsdde[ddsdde_dim1 + 2] = z2 + ddsdde[ddsdde_dim1 + 2];
	ddsdde[ddsdde_dim1 + 3] = z2 + ddsdde[ddsdde_dim1 + 3];
	ddsdde[ddsdde_dim1 * 3 + 2] = z2 + ddsdde[ddsdde_dim1 * 3 + 2];
	ddsdde[(ddsdde_dim1 << 1) + 3] = z2 + ddsdde[(ddsdde_dim1 << 1) + 3];
	i__2 = *ntens;
	for (i__ = 4; i__ <= i__2; ++i__) {
/* L130: */
	    ddsdde[i__ + i__ * ddsdde_dim1] = z3 + ddsdde[i__ + i__ * 
		    ddsdde_dim1];
	}
    } else {
	s_wsle(&io___104);
	do_lio(&c__9, &c__1, " this option is not implemented", (ftnlen)31);
	e_wsle();
	s_stop("", (ftnlen)0);
    }
    goto L400;

/* ---- elastic process begins here */

L300:
    for (i__ = 1; i__ <= 3; ++i__) {
/* L310: */
	stress[i__] += p;
    }

/* ---- form elastic stiffness matrix */

L315:
    i__2 = *ntens;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__1 = *ntens;
	for (j = 1; j <= i__1; ++j) {
/* L320: */
	    ddsdde[i__ + j * ddsdde_dim1] = 0.f;
	}
    }
    z1 = phi1 * (blk + twog * .66666666666666663);
    z2 = phi1 * (blk - twog * .33333333333333331);
    ddsdde[ddsdde_dim1 + 1] = z1;
    ddsdde[(ddsdde_dim1 << 1) + 2] = z1;
    ddsdde[ddsdde_dim1 * 3 + 3] = z1;
    ddsdde[(ddsdde_dim1 << 1) + 1] = z2;
    ddsdde[ddsdde_dim1 * 3 + 1] = z2;
    ddsdde[ddsdde_dim1 + 2] = z2;
    ddsdde[ddsdde_dim1 + 3] = z2;
    ddsdde[ddsdde_dim1 * 3 + 2] = z2;
    ddsdde[(ddsdde_dim1 << 1) + 3] = z2;
    i__1 = *ntens;
    for (i__ = 4; i__ <= i__1; ++i__) {
/* L330: */
	ddsdde[i__ + i__ * ddsdde_dim1] = phi1 * g;
    }

/* ---- clean up */

L400:
    statev[7] = max(statev[7],0.);
    return 0;
} /* umat_ */

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* __F2C__ */
