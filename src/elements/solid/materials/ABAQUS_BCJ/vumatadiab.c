/* vumatadiab.f -- translated by f2c (version 19971204).
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
static doublereal c_b17 = .5;

/* small number */
static doublereal kSmall = 1.0e-12;

/* $Id: vumatadiab.c,v 1.2 2004-08-01 20:42:35 paklein Exp $ */

/* 23456789012345678901234567890123456789012345678901234567890123456789012 */
/* revision a  2000/12/6  correct problems */
/*              jlathrop */

/* Revision 4.0 1995/03/10 in progress kistler */
/* modified to interface with abaqus explicit v. 5.4 */

/* Revision 3.0 1995/03/09 jlathrop */
/* modified code to interface with jas3d or pronto3d */

/* Revision 1.1  2001/07/18 21:29:45  paklein */
/* initial check in of ABAQUS/Explicit VUMAT interface */

/* Revision 2.5  1993/05/07  17:06:57  fjmello */
/* corrected typo in variable types Scalar */

/* Revision 2.4  1992/10/15  18:39:46  mwheins */
/* fixed status error in mat19 */

/* Revision 2.3  1992/08/19  22:32:26  swattaw */
/* fixed problem with status */

/* Revision 2.2  1992/07/31  15:15:56  swattaw */
/* added status to sandia damage model */

/* Revision 2.1  1992/07/28  17:03:55  swattaw */
/* added comments to sandia damage model */

/* Revision 2.0  1992/05/01  19:18:59  swattaw */
/* added self-contact */

/* Revision 1.2  1992/02/19  15:28:39  swattaw */
/* new version from snll */

/* Revision 1.1  1991/12/23  19:56:36  swattaw */
/* Initial revision */

/* Subroutine */ int vumat_(integer *nblock, integer *ndir, integer *nshr, 
	integer *nstatev, integer *nfieldv, integer *nprops, integer *lanneal,
	 doublereal *steptime, doublereal *totaltime, doublereal *dt, char *
	cmname, doublereal *coordmp, doublereal *charlength, doublereal *
	props, doublereal *density, doublereal *straininc, doublereal *
	relspininc, doublereal *tempold, doublereal *stretchold, doublereal *
	defgradold, doublereal *fieldold, doublereal *sigold, doublereal *
	svold, doublereal *enerinternold, doublereal *enerinelasold, 
	doublereal *tempnew, doublereal *stretchnew, doublereal *defgradnew, 
	doublereal *fieldnew, doublereal *sig, doublereal *sv, doublereal *
	enerinternnew, doublereal *enerinelasnew, ftnlen cmname_len)
{
    /* System generated locals */
    integer sigold_dim1, sigold_offset, svold_dim1, svold_offset, sig_dim1, 
	    sig_offset, sv_dim1, sv_offset, straininc_dim1, straininc_offset, 
	    i__1, i__2;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);
    double exp(doublereal), tanh(doublereal), sqrt(doublereal), log(
	    doublereal), d_sign(doublereal *, doublereal *), sinh(doublereal),
	     pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static doublereal beta, davg, scle, dsig, sige;
    static integer iqcm;
    static doublereal bulk, pold, dgam2;
    static integer ndirnshr;
    static doublereal d__[6]	/* was [1][6] */, gamma, dalph, shear, theta, 
	    ximag, h1[256], h2[256], press, ximag2, ak, c90;
    static integer ie;
    static doublereal alpmag, pr, ym, ftheta[256], bratio, vtheta[256], 
	    ytheta[256], ak2, sratio, epsdot, rd1[256], rd2[256], props28, 
	    rs1[256], rs2[256], ddd, adj, arg, xi11, xi22, xi33, xi12, xi23, 
	    xi31, psi, tmp;
    static integer isv;
    static doublereal sto, skappan, phi1, phi2, constnu, sto2;

/* read only abaqus variables */
/* write only (modifiable) variables */

/* ----------------------------------------------------------------------- */
/*            january 1995 version */
/*     DESCRIPTION: SANDIA DAMAGE MODEL */
/*                  (Bammann - Chiesa  Model) */

/*     FORMAL PARAMETERS: */
/*       NBLOCK   INTEGER   Number of elements in this block */
/*       NSTATEV  INTEGER   Number of internal state variables */
/*       DT       REAL      Time step */
/*       PROPS    REAL      Material properties */
/*                          (1) = Youngs modulus */
/*                          (2) = Poissons ratio */
/*                          (3) = Initial temperature */
/* 			   (4) = Thermal expansion coeff */
/* 			   (5) = Heat coeff ( 1/ rho*Cv ) */
/* 			   (6) = c1 - V const */
/*                          (7) = c2 - V temp */
/* 			   (8) = c3 - Y const */
/* 			   (9) = c4 - Y temp */
/* 			   (10) = c5 - f const */
/* 			   (11) = c6 - f temp */
/* 			   (12) = c7 - rd const */
/* 			   (13) = c8 - rd temp */
/* 			   (14) = c9 - h constant */
/* 			   (15) = c10 - h temp */
/* 			   (16) = c11 - rs constant */
/* 			   (17) = c12 - rs temp */
/* 			   (18) = c13 - Rd constant */
/* 			   (19) = c14 - Rd temp */
/* 			   (20) = c15 - H constant */
/* 	 		   (21) = c16 - H temp */
/* 			   (22) = c17 - Rs constant */
/* 			   (23) = c18 - Rs temp */
/*                          (24) = c19 -   adjustment to Y */
/*                          (25) = c20 -   adjustment to Y */
/*                          (26) = c21 -   adjustment to Y */
/*                          (27) = casb - for adiabatic temp evolution */
/* 			   (28) = damage constant n */
/*                          (29) = initial damage value */
/*       SIGOLD REAL        OLD Stresses */
/*       SIG    REAL        NEW Stresses */
/*       D      REAL        Strain (increment tensor at ea point) --div. by DT */
/* 						to get strain rate D */
/*       SVOLD  REAL        Internal state variables:  (SV is NEW VALUES) */
/*                          (*,1) = backstress tensor, Axx */
/*                          (*,2) = backstress tensor, Ayy */
/*                          (*,3) = backstress tensor, Azz */
/*                          (*,4) = backstress tensor, Axy */
/*                          (*,5) = backstress tensor, Ayz */
/*                          (*,6) = backstress tensor, Axz */
/*                          (*,7) = hardening scalar */
/*                          (*,8) = temperature */
/*                          (*,9) = equivalent plastic strain rate */
/*                          (*,10) = damage term */
/*                          (*,11) = rate of change of damage term */
/*                          (*,12) = equivalent plastic strain */
/*                          (*,13) = Total damage (to retain consistency with the shell version) */
/*                        The following are used by JAS only */
/*                          (*,14) = Youngs Modulus at start of step (jas only ) */
/*                          (*,15) = Youngs Modulus at end of step        " */
/*                          (*,16) = Poisson's Ratio at start of step     " */
/*                          (*,17) = Poisson's Ratio at end of step       " */
/* 			 The following is used by ABAQUS EXPLICIT */
/* 			   (*,18) = failure flag (0=failed, 1=not) */


/* ----------------------------------------------------------------------- */

/* ******************************************************************* */


/*     Double precision version */
/* m_vec_Length = Maximum vector (group) block length. */
/* - for TAHOE only do one integration point at a time */
/* - so replace NBLOCK with NBLOCK1 */




    /* Parameter adjustments */
    --tempnew;
    sigold_dim1 = *nblock;
    sigold_offset = sigold_dim1 + 1;
    sigold -= sigold_offset;
    --tempold;
    straininc_dim1 = *nblock;
    straininc_offset = straininc_dim1 + 1;
    straininc -= straininc_offset;
    sig_dim1 = *nblock;
    sig_offset = sig_dim1 + 1;
    sig -= sig_offset;
    sv_dim1 = *nblock;
    sv_offset = sv_dim1 + 1;
    sv -= sv_offset;
    svold_dim1 = *nblock;
    svold_offset = svold_dim1 + 1;
    svold -= svold_offset;
    --props;
    --coordmp;
    --charlength;
    --density;
    --relspininc;
    --stretchold;
    --defgradold;
    --fieldold;
    --enerinternold;
    --enerinelasold;
    --stretchnew;
    --defgradnew;
    --fieldnew;
    --enerinternnew;
    --enerinelasnew;

    /* Function Body */
    ndirnshr = *ndir + *nshr;
    props28 = (props[28] * 2 - 1) * 2 / (props[28] * 2 + 1);

/* initialize sv */

/* 	write(*,987) */
/* 987	format('before start') */
    i__1 = *nblock;
    for (ie = 1; ie <= i__1; ++ie) {
	i__2 = *nstatev;
	for (isv = 1; isv <= i__2; ++isv) {
	    sv[ie + isv * sv_dim1] = svold[ie + isv * svold_dim1];
/* L2: */
	}
	if (*totaltime <= *dt && *steptime <= *dt) {
	    sv[ie + sv_dim1 * 18] = 1.f;
	    sv[ie + (sv_dim1 << 3)] = props[3];
	    sv[ie + sv_dim1 * 10] = props[29];
	    sv[ie + sv_dim1 * 14] = props[1];
	    sv[ie + sv_dim1 * 15] = props[1];
	    sv[ie + (sv_dim1 << 4)] = props[2];
	    sv[ie + sv_dim1 * 17] = props[2];
	    if (sv[ie + sv_dim1] != 0.f) {
		s_stop("vumatadiab: bad state variable", 0L);
	    }
	    if (sv[ie + (sv_dim1 << 1)] != 0.f) {
		s_stop("vumatadiab: bad state variable", 0L);
	    }
	    if (sv[ie + sv_dim1 * 3] != 0.f) {
		s_stop("vumatadiab: bad state variable", 0L);
	    }
	    if (sv[ie + (sv_dim1 << 2)] != 0.f) {
		s_stop("vumatadiab: bad state variable", 0L);
	    }
	    if (sv[ie + sv_dim1 * 5] != 0.f) {
		s_stop("vumatadiab: bad state variable", 0L);
	    }
	    if (sv[ie + sv_dim1 * 6] != 0.f) {
		s_stop("vumatadiab: bad state variable", 0L);
	    }
	    i__2 = ndirnshr;
	    for (isv = 1; isv <= i__2; ++isv) {
		if (sigold[ie + isv * sigold_dim1] != 0.f) {
		    s_stop("vumatadiab: bad state variable", 0L);
		}
/* L3: */
	    }
	}
	i__2 = ndirnshr;
	if (*dt > 1.0e-12)
		for (isv = 1; isv <= i__2; ++isv) {
			d__[ie + isv - 2] = straininc[ie + isv * straininc_dim1] / *dt;
/* L333: */
		}
	else /* allow dt -> 0 */
		for (isv = 1; isv <= i__2; ++isv) {
			d__[ie + isv - 2] = 0.0;
		}
/* L1: */
    }

/* ***pppp   the  following 7 lines for nontemp. dep. young's & poisson's */
/* 	write(*,991) */
/* 991	format('before ym definition') */
    iqcm = 0;
    ym = props[1];
    pr = props[2];
    shear = ym / (pr + 1.f);
    bulk = ym / ((1.f - pr * 2.f) * 3.f);
    sratio = 1.f;
    bratio = 1.f;
/* ***jjjj  the following 14 lines for jas3d modified for abaexp */
/* 	write(*,982) */
/* 982	format('before precalculated') */
    if (props[3] == 0.f) {
/*    PRECALCULATED MATERIAL PARAMETERS */
	i__1 = *nblock;
	for (ie = 1; ie <= i__1; ++ie) {
	    vtheta[ie - 1] = props[6];
	    ytheta[ie - 1] = props[8];
	    ftheta[ie - 1] = props[10];
	    rd1[ie - 1] = props[12];
	    h1[ie - 1] = props[14];
	    rs1[ie - 1] = props[16];
	    rd2[ie - 1] = props[18];
	    h2[ie - 1] = props[20];
	    rs2[ie - 1] = props[22];
/* L5: */
	}
    } else {
/* 	write(*,983) */
/* 983	format('afterprecalculated') */

/* 	write(*,994) */
/* 994	format('before do 10 definition') */
	i__1 = *nblock;
	for (ie = 1; ie <= i__1; ++ie) {
	    if (sv[ie + sv_dim1 * 18] < .01f) {
		goto L10;
	    }
/* ***  TEMPERATURE */
	    theta = sv[ie + (sv_dim1 << 3)];
	    if (theta == 0.f) {
		s_stop("vumatadiab: bad temperature", 0L);
	    }

/* ***  temperatured dependent parameters */
/* ***       yield */
	    vtheta[ie - 1] = props[6] * exp(-props[7] / theta);
/*          YTHETA(IE) = PROPS(8) * EXP( PROPS(9)/THETA) */
	    ytheta[ie - 1] = props[8] / (exp(-props[9] / theta) + props[26]);
	    ftheta[ie - 1] = props[10] * exp(-props[11] / theta);
/* ****       kinematic hardening */
	    rd1[ie - 1] = props[12] * exp(-props[13] / theta);
	    h1[ie - 1] = props[14] - theta * props[15];
	    rs1[ie - 1] = props[16] * exp(-props[17] / theta);
/* ****       isotropic hardening */
	    rd2[ie - 1] = props[18] * exp(-props[19] / theta);
	    h2[ie - 1] = props[20] - theta * props[21];
	    rs2[ie - 1] = props[22] * exp(-props[23] / theta);
/* ****       adjustment to y (if used) */
	    adj = 1.f;
	    if (props[24] > 0.f) {
/* Computing MAX */
		d__1 = 0.f, d__2 = props[25] - theta;
		adj = (tanh(props[24] * max(d__1,d__2)) + 1.f) * .5f;
	    }
	    ytheta[ie - 1] = adj * ytheta[ie - 1];
/* ****       adjustment to y (if used) */
/*         ADJ=1.0 */
/*         IF (PROPS(24) .GT. 0.0) ADJ =0.5*(1.+ TANH(PROPS(24)* */
/*    *                       MAX(0.,PROPS(25) - THETA))) */
/*         YTHETA(IE) = ADJ * YTHETA(IE) */
L10:
	    ;
	}
/* ***jjjj  this line for jas3d */
    }

    i__1 = *nblock;
    for (ie = 1; ie <= i__1; ++ie) {
	if (sv[ie + sv_dim1 * 18] < .01f) {
	    goto L20;
	}

/* ***jjjj   the following 8 lines for temperature dep. poisson & young's */
/*        YM = SV(IE,15) */
/*        PR = SV(IE,17) */
/*        SHEAR = YM / ( 1. + PR) */
/*        BULK = YM / (3. * (1. -2.*PR)) */
/*        SHEAR0 = SV(IE,14)/ (1. + SV(IE,16)) */
/*        SRATIO = SHEAR/SHEAR0 */
/*        BULK0 = SV(IE,14)/(3. *(1. - 2.*SV(IE,16))) */
/*        BRATIO = BULK/BULK0 */
/* volumetric strain */
	davg = (d__[ie - 1] + d__[ie] + d__[ie + 1]) * .33333333333333331;

/* strain rate mag */
	ddd = sv[ie + sv_dim1 * 9];
/* 	DDD = ROOT23 * SQRT((D(IE,1)-DAVG)**2 + */
/*     *                      (D(IE,2)-DAVG)**2 + */
/*     *                      (D(IE,3)-DAVG)**2 + */
/*     1         .5 * (D(IE,4)**2 + D(IE,5)**2 + D(IE,6)**2) ) */
/* c     2		/DT */

/* back stress mag */
/* Computing 2nd power */
	d__1 = sv[ie + sv_dim1];
/* Computing 2nd power */
	d__2 = sv[ie + (sv_dim1 << 1)];
/* Computing 2nd power */
	d__3 = sv[ie + sv_dim1 * 3];
/* Computing 2nd power */
	d__4 = sv[ie + (sv_dim1 << 2)];
/* Computing 2nd power */
	d__5 = sv[ie + sv_dim1 * 5];
/* Computing 2nd power */
	d__6 = sv[ie + sv_dim1 * 6];
	alpmag = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + (d__4 * d__4 
		+ d__5 * d__5 + d__6 * d__6) * 2.f) * .8164965809;

/*  recovery terms */
	sto = *dt * (rd1[ie - 1] * ddd + rs1[ie - 1]) * alpmag;
	sto = min(1.,sto);
/* 	  STO2 = DT * (RD2(IE) * DDD + RS2(IE)) * ABS(SV(IE,7)) */
	sto2 = *dt * (rd2[ie - 1] * ddd + rs2[ie - 1] * (d__1 = sv[ie + 
		sv_dim1 * 7], abs(d__1)));
	sto2 = min(1.,sto2);
/* backstress */
	sv[ie + sv_dim1] -= sto * sv[ie + sv_dim1];
	sv[ie + (sv_dim1 << 1)] -= sto * sv[ie + (sv_dim1 << 1)];
	sv[ie + sv_dim1 * 3] -= sto * sv[ie + sv_dim1 * 3];
	sv[ie + (sv_dim1 << 2)] -= sto * sv[ie + (sv_dim1 << 2)];
	sv[ie + sv_dim1 * 5] -= sto * sv[ie + sv_dim1 * 5];
	sv[ie + sv_dim1 * 6] -= sto * sv[ie + sv_dim1 * 6];

/* scalar hardening */
	skappan = sv[ie + sv_dim1 * 7];
	sv[ie + sv_dim1 * 7] -= sto2 * sv[ie + sv_dim1 * 7];
/* damage */
	phi1 = 1.f - sv[ie + sv_dim1 * 10];
/* Computing MIN */
	d__1 = 1.f, d__2 = *dt * sv[ie + sv_dim1 * 11] / phi1;
	phi2 = 1.f - min(d__1,d__2);

/* yield strength */
	if (ftheta[ie - 1] == 0.f) {
	    ftheta[ie - 1] = 1e-5;
	}
	d__1 = ddd / ftheta[ie - 1];
	ak = phi1 * (vtheta[ie - 1] * log(d__1 + sqrt(d__1 * d__1 + 1.f)) + 
		ytheta[ie - 1] + sv[ie + sv_dim1 * 7]);
/*   	  AK = PHI1 * ( VTHETA(IE) * */
/*    1            alog((ddd+sqrt(ddd**2.+ftheta**2.))/ftheta) + */
/* 	write(*,984)ddd,ytheta(ie),sv(ie,7) */
/* 984	format('VARIABLES...',4E15.5) */

/* old pressure */
	pold = (sigold[ie + sigold_dim1] + sigold[ie + (sigold_dim1 << 1)] + 
		sigold[ie + sigold_dim1 * 3]) * .33333333333333331;

/*  new pressure */
/* 	  PRESS = BRATIO*POLD * PHI2 + PHI1 * BULK * 3. * DAVG */
	press = bratio * pold * phi2 + *dt * phi1 * bulk * 3.f * davg;

/* predicted deviatoric part of the stress  (eq 11) */
	sig[ie + sig_dim1] = sratio * phi2 * (sigold[ie + sigold_dim1] - pold)
		 + phi1 * *dt * shear * (d__[ie - 1] - davg);
/*    *                PHI1*SHEAR*(D(IE,1)-DAVG) */
	sig[ie + (sig_dim1 << 1)] = sratio * phi2 * (sigold[ie + (sigold_dim1 
		<< 1)] - pold) + phi1 * *dt * shear * (d__[ie] - davg);
/*    *                PHI1*SHEAR*(D(IE,2)-DAVG) */
	sig[ie + sig_dim1 * 3] = sratio * phi2 * (sigold[ie + sigold_dim1 * 3]
		 - pold) + phi1 * *dt * shear * (d__[ie + 1] - davg);
/*    *                PHI1*SHEAR*(D(IE,3)-DAVG) */
	sig[ie + (sig_dim1 << 2)] = sratio * phi2 * sigold[ie + (sigold_dim1 
		<< 2)] + phi1 * *dt * shear * d__[ie + 2] / 2.f;
	sig[ie + sig_dim1 * 5] = sratio * phi2 * sigold[ie + sigold_dim1 * 5] 
		+ phi1 * *dt * shear * d__[ie + 3] / 2.f;
	sig[ie + sig_dim1 * 6] = sratio * phi2 * sigold[ie + sigold_dim1 * 6] 
		+ phi1 * *dt * shear * d__[ie + 4] / 2.f;
/* 	  SIG(IE,4) = SRATIO*PHI2*SIGOLD(IE,4) + PHI1*SHEAR* D(IE,4) */
/*   	  SIG(IE,5) = SRATIO*PHI2*SIGOLD(IE,5) + PHI1*SHEAR* D(IE,5) */
/* 	  SIG(IE,6) = SRATIO*PHI2*SIGOLD(IE,6) + PHI1*SHEAR* D(IE,6) */

/* equivalent stress (eq 12)  (sig-2/3 alpha) */
	xi11 = sig[ie + sig_dim1] - sv[ie + sv_dim1] * .66666666666666663;
	xi22 = sig[ie + (sig_dim1 << 1)] - sv[ie + (sv_dim1 << 1)] * 
		.66666666666666663;
	xi33 = sig[ie + sig_dim1 * 3] - sv[ie + sv_dim1 * 3] * 
		.66666666666666663;
	xi12 = sig[ie + (sig_dim1 << 2)] - sv[ie + (sv_dim1 << 2)] * 
		.66666666666666663;
	xi23 = sig[ie + sig_dim1 * 5] - sv[ie + sv_dim1 * 5] * 
		.66666666666666663;
	xi31 = sig[ie + sig_dim1 * 6] - sv[ie + sv_dim1 * 6] * 
		.66666666666666663;

/* equivalent stress magnitude */
/* Computing 2nd power */
	d__1 = xi11;
/* Computing 2nd power */
	d__2 = xi22;
/* Computing 2nd power */
	d__3 = xi33;
/* Computing 2nd power */
	d__4 = xi12;
/* Computing 2nd power */
	d__5 = xi23;
/* Computing 2nd power */
	d__6 = xi31;
	ximag2 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + (d__4 * d__4 + 
		d__5 * d__5 + d__6 * d__6) * 2.f;

/* check elastic assumption */
	ak2 = ximag2 - ak * .66666666666666663 * ak;

/*  scle = 0 for elastic,  1 for plastic */
	scle = d_sign(&c_b17, &ak2) + .5f;
/*  adjust so we dont divide by zero if there is no stress */
	ximag = sqrt(ximag2) + 1.f - scle;
	ddd = ddd + 1.f - scle;
/*  if plastic, calculate correction to the Integral of Plastic strain */

	gamma = scle * (ximag - ak * .8164965809) / (phi1 * shear + (h1[ie - 
		1] + h2[ie - 1] * phi1) * .66666666666666663);
/* 	write(*,993)PHI1,ROOT,SCLE,GAMMA */
/* 993	format('after PRESS def',4e15.5) */
	dgam2 = gamma / ximag;
	dsig = phi1 * shear * gamma / ximag;
/* correction to stress (eq 16) */
	sig[ie + sig_dim1] = sig[ie + sig_dim1] - dsig * xi11 + press;
	sig[ie + (sig_dim1 << 1)] = sig[ie + (sig_dim1 << 1)] - dsig * xi22 + 
		press;
	sig[ie + sig_dim1 * 3] = sig[ie + sig_dim1 * 3] - dsig * xi33 + press;
	sig[ie + (sig_dim1 << 2)] -= dsig * xi12;
	sig[ie + sig_dim1 * 5] -= dsig * xi23;
	sig[ie + sig_dim1 * 6] -= dsig * xi31;
/* update scalar hardening and backstress */
	sv[ie + sv_dim1 * 7] += gamma * .8164965809 * h2[ie - 1];

	dalph = h1[ie - 1] * dgam2;
	sv[ie + sv_dim1] += dalph * xi11;
	sv[ie + (sv_dim1 << 1)] += dalph * xi22;
	sv[ie + sv_dim1 * 3] += dalph * xi33;
	sv[ie + (sv_dim1 << 2)] += dalph * xi12;
	sv[ie + sv_dim1 * 5] += dalph * xi23;
	sv[ie + sv_dim1 * 6] += dalph * xi31;

/*   plastic strain rate */
	if (*dt > 1.0e-12)
		epsdot = gamma * .8164965809 / *dt;
	else /* allow dt -> 0 */
		epsdot = 0.0;
	sv[ie + sv_dim1 * 9] = epsdot;
/*   and save equivalent plastic strain */
	sv[ie + sv_dim1 * 12] += gamma * .8164965809;
/* vonmises stress */
/* Computing MAX */
/* Computing 2nd power */
	d__3 = sig[ie + sig_dim1] - sig[ie + (sig_dim1 << 1)];
/* Computing 2nd power */
	d__4 = sig[ie + (sig_dim1 << 1)] - sig[ie + sig_dim1 * 3];
/* Computing 2nd power */
	d__5 = sig[ie + sig_dim1 * 3] - sig[ie + sig_dim1];
/* Computing 2nd power */
	d__6 = sig[ie + (sig_dim1 << 2)];
/* Computing 2nd power */
	d__7 = sig[ie + sig_dim1 * 5];
/* Computing 2nd power */
	d__8 = sig[ie + sig_dim1 * 6];
	d__1 = 1e-15f, d__2 = sqrt((d__3 * d__3 + d__4 * d__4 + d__5 * d__5) *
		 .33333333333333331 + (d__6 * d__6 + d__7 * d__7 + d__8 * 
		d__8) * 2.f);
	sige = max(d__1,d__2);
/* 	write(*,992) */
/* 992	format('before damage definition') */
/* damage calculation */
/*  bounds of 0. and 15. below are unitless and merely */
/*  avoid potential numerical problems when there should not */
/*  be any damage */
/* 	  ARG = MIN( 15., PRESS*PROPS28/SIGE) */
/* JFL use yield strength if press/sige ratio too big */
	if (press / sige < 3.f) {
/* Computing MIN */
	    d__1 = 15.f, d__2 = press * props28 / sige;
	    arg = min(d__1,d__2);
	} else {
/* Computing MIN */
	    d__1 = 15.f, d__2 = press * props28 / ytheta[ie - 1];
	    arg = min(d__1,d__2);
	}
	beta = sinh((max(0.,arg)));

	c90 = props[28] + 1.f;
/* Computing MIN */
	d__1 = 15.f, d__2 = beta * *dt * epsdot * c90;
	psi = min(d__1,d__2);
/* Computing MAX */
	d__1 = 0.f, d__2 = (pow_dd(&phi1, &c90) - 1.f) * exp(psi) + 1.f;
	tmp = max(d__1,d__2);
/* damage (eq 18) bounded away from 1 - could use .9999 if desired */
/*                however this is adequate since damage greater than */
/*                .6 usually goes to 1 (or .99) in one step */
/* Computing MIN */
	d__2 = 1.f / c90;
	d__1 = 1.f - pow_dd(&tmp, &d__2);
	sv[ie + sv_dim1 * 10] = min(d__1,.99);
/*    and also store in total damage */
	sv[ie + sv_dim1 * 13] = sv[ie + sv_dim1 * 10];

/* abaqus - failure state variable- 1==not failed, 0==failed */
/*         where 0.01 below relates to 0.99 above to equal 1.0 */

	sv[ie + sv_dim1 * 18] = (real) (1 - (integer) (sv[ie + sv_dim1 * 10] 
		+ .0101f));
/* rate of damage (eq 7) */
	d__1 = 1.f - sv[ie + sv_dim1 * 10];
	sv[ie + sv_dim1 * 11] = beta * epsdot * (1.f / pow_dd(&d__1, &props[
		28]) - (1.f - sv[ie + sv_dim1 * 10]));

/* update the temp based on plastic work */
/* 	  SV(IE,8) = SV(IE,8) + */
/*     *    PROPS(5) * DGAM2 * */
/*     *   (SIG(IE,1)*XI11 + SIG(IE,2)*XI22 + SIG(IE,3)*XI33 + */
/*     * 2.*(SIG(IE,4)*XI12 + SIG(IE,5)*XI23 + SIG(IE,6)*XI31) ) */

/*         CONSTNU=(PONE+PR)/8.0D0 */
	if (props[27] < 1e-8) {
	    props[5] = 0.;
	    constnu = 1.;
	} else {
	    constnu = props[27];
	}
/* Computing 2nd power */
	d__1 = sv[ie + sv_dim1 * 7];
/* Computing 2nd power */
	d__2 = skappan;
	sv[ie + (sv_dim1 << 3)] += props[5] * (dgam2 * (sig[ie + sig_dim1] * 
		xi11 + sig[ie + (sig_dim1 << 1)] * xi22 + sig[ie + sig_dim1 * 
		3] * xi33 + (sig[ie + (sig_dim1 << 2)] * xi12 + sig[ie + 
		sig_dim1 * 5] * xi23 + sig[ie + sig_dim1 * 6] * xi31) * 2.) - 
		1. / (shear * constnu) * .5 * (d__1 * d__1 - d__2 * d__2));

L20:
	;
    }

    return 0;
} /* vumat_ */

#ifdef __cplusplus
} /* extern "C" */
#endif

#else
static void dummy(void)
{
	int a = 1;
}

#endif /* __F2C__ */
