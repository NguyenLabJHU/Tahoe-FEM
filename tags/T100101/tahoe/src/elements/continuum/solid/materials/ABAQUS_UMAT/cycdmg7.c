/* cycdmg7.f -- translated by f2c (version 19971204).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

/* library support options */
#ifdef __F2C__

#include "f2c.h"

/* Table of constant values */

static doublereal c_b2 = 2.;
static integer c__9 = 9;
static integer c__1 = 1;
static doublereal c_b19 = 3.;
static doublereal c_b20 = 1.5;
static doublereal c_b54 = .5;

/* ************************************************************************* */

/*              USCAR-USAMP CRADA property model developed at Sandia. */
/*              The BAMMANN - CHIESA - JOHNSON  Plasticity Model */
/*              was modified to include dendrite cell size effects and */
/*              damage nucleation, growth, and coalescence were */
/*              added by Mark F. Horstemeyer */

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
/* **                        Mark F. Horstemeyer                          ** */
/* **                       Sandia National Labs                          ** */
/* **                         Mail Stop 9405                              ** */
/* **                            PO Box 969                               ** */
/* **                     Livermore, CA 94551-0969                        ** */
/* **                                                                     ** */
/* **                   email:    mfhorst@sandia.gov                      ** */
/* **                                                                     ** */
/* **   WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING     ** */
/* **                                                                     ** */
/* ************************************************************************* */
/* ************************************************************************* */

/*                                Copyright */
/*                         Sandia Corporation 1998 */
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

/*        Sandia strain rate and temperature dependent plasticity model */
/*        was originally implemented by Mike Chiesa. */
/*        This subroutine can be used for axisymmetric, plane strain and */
/*        3d solid elements. */
/*        Porosity was added by Mark Horstemeyer 1/8/97. */
/*        Void Nucleation was added by Mark Horstemeyer 9/28/97. */
/*        Void coalescence was added by Mark Horstemeyer 3/1/98. */
/*        Dendrite cell size was added by Mark Horstemeyer 3/15/98. */

/* ************************************************************************* */

/*                          w a r n i n g */

/*                 this is proprietary software */

/*             recipients of this routine are not allowed to distribute */
/*             their source to other facilities.  rather, request for */
/*             the source should be directed to the author at snll. */
/*             phone: 510-294-2816 */


/* *********************************************************************** */

/*                          l e g a l   n o t i c e */

/*             this computer code material was prepared as an account */
/*             of work sponsored by the united states government. */
/*             neither the united states nor the united states */
/*             department of energy, nor any of their employees, */
/*             nor any of their contractors, subcontractors, or their */
/*             employees, makes any warranty, express or implied, or */
/*             assumes any legal liability or responsibility for the */
/*             accuracy, completeness or usefulness of any information, */
/*             apparatus, produce or process disclosed, or represents */
/*             that its use would not infringe privately-owned rights */

/* ************************************************************************* */
/* Subroutine */ int umat_(doublereal *stress, doublereal *statev, doublereal 
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

    static doublereal third = .333333333333333;
    static doublereal twothd = .666666666666667;
    static doublereal con1 = .81649658092773;
    static doublereal pi = 3.1415927;
    static integer iparam1 = 0;
    static integer iparam2 = 0;
    static integer iparam3 = 0;

    /* System generated locals */
    integer ddsdde_dim1, ddsdde_offset, i__1, i__2;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    /* Subroutine */ int s_stop(char *, ftnlen);
    double exp(doublereal), tanh(doublereal), sqrt(doublereal), log(
	    doublereal), sinh(doublereal);

    /* Local variables */
    static doublereal dgam, beta, davg, dsig, sige, alpm, vrad, htcp, pold, 
	    twog, dgam2, alpm2, zzzz, term1, term2, term3, term4, term5, 
	    term6, term7, term8, term9, g;
    static integer i__, j;
    static doublereal p, gamma, cacon, r__, dalph, theta, ximag, dterm, slope,
	     h1, h2, h3, r1, r2, r3, z1, z2, z3, ximag2, ca, cb, cf, ak, c90, 
	    damage, xi[6], ftheta, zz, trialk, cc1, cc2, cc3, cc4, cc5, cc6, 
	    cc7, cc8, cc9, cd1, cd2, tratio, dj2, dj3, vtheta, ytheta, rd1, 
	    rd2, rd3, ak2, epsdot, di1, rs1, rs2, rs3, zthird, abc, cc10, 
	    cc11, cc12, cc13, cc14, cc15, cc16, cc17, cc18, cc19, cc20, cc21, 
	    cc22, cc23, cc24, cc25, cc26, adj, blk, dcs, ds11, ds22, ds33, 
	    ds12, ds23, ds13, ddd, dum, arg, psi, bbb, tmp, sto, dam1, dam2, 
	    alphaxx, alphaxy, alphayy, alphayz, alphazz, dcs0, alphazx, 
	    zsecond, zzz, phi1, sto2;

    /* Fortran I/O blocks */
    static cilist io___8 = { 0, 6, 0, 0, 0 };


/* *************************************************** */

    /* Parameter adjustments */
    --predef;
    --dpred;
    --dstran;
    --stran;
    --drplde;
    --ddsddt;
    ddsdde_dim1 = *ntens;
    ddsdde_offset = ddsdde_dim1 + 1;
    ddsdde -= ddsdde_offset;
    --stress;
    --statev;
    --props;
    --coords;
    drot -= 4;
    dfgrd0 -= 4;
    dfgrd1 -= 4;

    /* Function Body */

/* *********************************************************************** */

/*     *  ntens = number of non-zero stress components (4 for 2d, 6 for 3d) */
/*     *  nshr  = number of non-zero shear components (1 for 2d, 3 for 3d) */
/*     *  ndi   = number of non-zero normal stresses (always 3) */
/*     *  nstatv = number of state variables (25) */
/*     *  nprops = number of material parameters (50) */

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
/*     *  statev(10) = McClintock void growth(second phase pores) */
/*     *  statev(11) = rate of change of M porosity */
/*     *  statev(12) = stress triaxiality */
/*     *  statev(13) = nucleation */
/*     *  statev(14) = damage */
/*     *  statev(15) = nucleation rate */
/*     *  statev(16) = damage rate */
/*     *  statev(17) = nucleation from previous time step */
/*     *  statev(18) = Cocks-Ashby void growth(large pores) */
/*     *  statev(19) = rate of change of CA porosity */
/*     *  statev(20) = alpha-xx long range */
/*     *  statev(21) = alpha-yy long range */
/*     *  statev(22) = alpha-zz long range */
/*     *  statev(23) = alpha-xy long range */
/*     *  statev(24) = alpha-yz long range */
/*     *  statev(25) = alpha-zx long range */

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
/*     *  props(25)= C20             , props(26)= CA */
/*     *  props(27)= CB */
/*     *  props(28)= initial temperature */
/*     *  props(29)= heat generation coefficient */
/*     *  props(30)= McClintock damage constant, n */
/*     *  props(31)= initial void radius */
/*     *  props(32)= torsional constant a in nucleation model */
/*     *  props(33)= tension/comp constant b in nucleation model */
/*     *  props(34)= triaxiality constant c in nucleation model */
/*     *  props(35)= coefficient constant in nucleation model */
/*     *  props(36)= fracture toughness, related to nucleation model */
/*     *  props(37)= ave size of particles, related to nucleation model */
/*     *  props(38)= particles vol fraction, related to nucleation model */
/*     *  props(39)= coalescence factor, D=nucleation*void volume*coal. */
/*     *  props(40)= coalescence factor, D=nucleation*void volume*coal. */
/*     *  props(41)= reference grain size or dendrite cell size,dcs0 */
/*     *  props(42)= grain size or dendrite cell size of material,dcs */
/*     *  props(43)= grain size or dendrite cell size exponent,zz */
/*     *  props(44)= initial void volume fraction for CA void growth */
/*     *  props(45)= C21             , props(46)= C22 */
/*     *  props(47)= C23             , props(48)= C24 */
/*     *  props(49)= C25             , props(50)= C26 */
/*     *  props(51)= nucleation temperature dependence */
/*     *  props(52)= coalescence temperature dependence */

/* ****************************************************************** */

/*   * iparam1 = 0 for linear return (uses total strain in recovery) */
/*   *         = 1 for quadratic return (uses plastic strain in recovery) */
/*   *           (only linear return implemented in this version) */

/*   * iparam2 = 0 for Simo tangent stiffness matrix */
/*   *         = 1 for Lathrop tangent stiffness matrix */
/*   *           (only Simo stiffness implemented in this version) */

/*   * iparam3 = 0 for trial kappa = kappa(n) */
/*   *         = 1 for trial kappa = kappa(n + 1/2) */

/* ************************************************************************ */
/* ---- initialize state variables */
/*     temperature is set first by *initial condition command in abaqus */
/*     or by props(28) */

    if (*time == 0.f) {
	statev[1] = 0.f;
	statev[2] = 0.f;
	statev[3] = 0.f;
	statev[4] = 0.f;
	statev[5] = 0.f;
	statev[6] = 0.f;
	statev[7] = 0.f;
	statev[8] = props[28];
	statev[9] = 0.f;
	statev[10] = pi * pow_dd(&props[31], &c_b2);
	statev[11] = 0.f;
	statev[12] = 0.f;
	statev[13] = props[35];
	statev[14] = props[35] * statev[10];
	statev[15] = 0.f;
	statev[16] = 0.f;
	statev[17] = 0.f;
	statev[18] = props[44];
	statev[19] = 0.f;
	statev[20] = 0.f;
	statev[21] = 0.f;
	statev[22] = 0.f;
	statev[23] = 0.f;
	statev[24] = 0.f;
	statev[25] = 0.f;
	if (*temp == 0.f) {
	    if (props[28] == 0.f) {
		s_wsle(&io___8);
		do_lio(&c__9, &c__1, " error - temperature is zero", 28L);
		e_wsle();
		s_stop("", 0L);
	    } else {
		statev[8] = props[28];
	    }
	} else {
	    statev[8] = *temp;
	}
    }

/* ****************************************************************** */

    cc1 = props[6];
    cc2 = props[7];
    cc3 = props[8];
    cc4 = props[9];
    cc5 = props[10];
    cc6 = props[11];
    cc7 = props[12];
    cc8 = props[13];
    cc9 = props[14];
    cc10 = props[15];
    cc11 = props[16];
    cc12 = props[17];
    cc13 = props[18];
    cc14 = props[19];
    cc15 = props[20];
    cc16 = props[21];
    cc17 = props[22];
    cc18 = props[23];
    cc19 = props[24];
    cc20 = props[25];
    cc21 = props[45];
    cc22 = props[46];
    cc23 = props[47];
    cc24 = props[48];
    cc25 = props[49];
    cc26 = props[50];
    ca = props[26];
    cb = props[27];
    cd1 = props[39];
    cd2 = props[40];
    dcs0 = props[41];
    dcs = props[42];
    zz = props[43];

    htcp = props[29];

/* ---- g = shear modulus   twog = 2*g   blk = bulk modulus */

    if (props[28] == 0.f) {
	statev[8] = *temp;
    }
    theta = statev[8];
    if (props[5] == 0.) {
	blk = props[3];
	g = props[1];
    } else {
	tratio = theta / props[5];
	tratio = min(tratio,.9999);
	g = props[1] * (1.f - tratio * exp(props[2] * (1.f - 1.f / tratio)));
    }
    twog = g * 2.f;
    blk = props[3] - props[4] * tratio;

/* ---- damage */
    dam1 = 1.f - statev[14];
/* Computing MIN */
    d__1 = 1.f, d__2 = *dtime * statev[16] / dam1;
    dam2 = 1.f - min(d__1,d__2);
    phi1 = 1.f - statev[18];

/* ---- calculate pressure */

    davg = third * (dstran[1] + dstran[2] + dstran[3]);
    pold = third * (stress[1] + stress[2] + stress[3]);
    p = pold * dam2 + dam1 * blk * davg * 3.f;

/* ---- check for melt */

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
	goto L315;
    }

/* ---- compute function evaluations */
/*       theta = temperature */
/*       ytheta = static yield stress */
/*       vtheta,ftheta = functions to define rate dependence of yield */
/*       h1,h2,h3 = plastic hardeing moduli */
/*       rs1,rs2,h3 = static recovery functions */
/*       rd1,rd2,h3 = dynamic recovery functions */

/* deviatoric stress */
    ds11 = stress[1] - p;
    ds22 = stress[2] - p;
    ds33 = stress[3] - p;
    ds12 = stress[4];
    ds23 = stress[5];
    ds13 = stress[6];
/* invariants of stress */
    dj2 = (pow_dd(&ds11, &c_b2) + pow_dd(&ds22, &c_b2) + pow_dd(&ds33, &c_b2) 
	    + (pow_dd(&ds12, &c_b2) + pow_dd(&ds23, &c_b2) + pow_dd(&ds13, &
	    c_b2)) * 2) * .5f;
    dj3 = ds11 * (ds22 * ds33 - pow_dd(&ds23, &c_b2)) - ds22 * (ds11 * ds33 - 
	    pow_dd(&ds13, &c_b2)) + ds33 * (ds22 * ds11 - pow_dd(&ds12, &c_b2)
	    );

/* stress state dependent material constants */
    adj = (tanh(cc19 * (cc20 - theta)) + 1.f) * .5f;
    if (cc19 == 0.f) {
	adj = 1.f;
    }
    vtheta = cc1 * exp(-cc2 / theta);
    ytheta = cc3 * exp(cc4 / theta) * adj;
    ftheta = cc5 * exp(-cc6 / theta);
    if (dj2 == 0.) {
	rd1 = cc7 * (1 - ca * .14814814814814814f) * exp(-cc8 / theta);
	h1 = (cc9 - cc10 * theta) * (ca * .14814814814814814f + 1);
	rs1 = cc11 * exp(-cc12 / theta);
	rd2 = cc13 * (1 - ca * .14814814814814814f) * exp(-cc14 / theta);
	h2 = (cc15 - cc16 * theta) * (ca * .14814814814814814f + 1);
	rs2 = cc17 * exp(-cc18 / theta);
	rd3 = cc21 * (1 - ca * .14814814814814814f) * exp(-cc22 / theta);
	h3 = (cc23 - cc24 * theta) * (ca * .14814814814814814f + 1);
	rs3 = cc25 * exp(-cc26 / theta);
    } else {
	rd1 = cc7 * (1 - ca * (.14814814814814814f - pow_dd(&dj3, &c_b2) / 
		pow_dd(&dj2, &c_b19)) - cb * dj3 / pow_dd(&dj2, &c_b20)) * 
		exp(-cc8 / theta);
	h1 = (cc9 - cc10 * theta) * (ca * (.14814814814814814f - pow_dd(&dj3, 
		&c_b2) / pow_dd(&dj2, &c_b19)) + 1 + cb * dj3 / pow_dd(&dj2, &
		c_b20));
	rs1 = cc11 * exp(-cc12 / theta);
	rd2 = cc13 * (1 - ca * (.14814814814814814f - pow_dd(&dj3, &c_b2) / 
		pow_dd(&dj2, &c_b19)) - cb * dj3 / pow_dd(&dj2, &c_b20)) * 
		exp(-cc14 / theta);
	h2 = (cc15 - cc16 * theta) * (ca * (.14814814814814814f - pow_dd(&dj3,
		 &c_b2) / pow_dd(&dj2, &c_b19)) + 1 + cb * dj3 / pow_dd(&dj2, 
		&c_b20));
	rs2 = cc17 * exp(-cc18 / theta);
	rd3 = cc21 * (1 - ca * (.14814814814814814f - pow_dd(&dj3, &c_b2) / 
		pow_dd(&dj2, &c_b19)) - cb * dj3 / pow_dd(&dj2, &c_b20)) * 
		exp(-cc22 / theta);
	h3 = (cc23 - cc24 * theta) * (ca * (.14814814814814814f - pow_dd(&dj3,
		 &c_b2) / pow_dd(&dj2, &c_b19)) + 1 + cb * dj3 / pow_dd(&dj2, 
		&c_b20));
	rs3 = cc25 * exp(-cc26 / theta);
    }

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
/* ---- update long range alpha using abaqus rotation matrix */

    if (rs3 != 0.f || rd3 != 0.f || h3 != 0.f) {
	if (*ntens == 4) {
	    term1 = drot[4] * statev[20] + drot[7] * statev[23];
	    term2 = drot[4] * statev[23] + drot[7] * statev[21];
	    term3 = drot[5] * statev[20] + drot[8] * statev[23];
	    term4 = drot[5] * statev[23] + drot[8] * statev[21];
	    term5 = statev[20] + statev[21];
/* Computing 2nd power */
	    d__1 = drot[10];
	    statev[20] = drot[4] * term1 + drot[7] * term2 - term5 * (d__1 * 
		    d__1);
/* Computing 2nd power */
	    d__1 = drot[11];
	    statev[21] = drot[5] * term3 + drot[8] * term4 - term5 * (d__1 * 
		    d__1);
	    statev[22] = -(statev[20] + statev[21]);
	    statev[23] = drot[5] * term1 + drot[8] * term2 - term5 * drot[10] 
		    * drot[11];
	} else {
	    term1 = drot[4] * statev[20] + drot[7] * statev[23] + drot[10] * 
		    statev[25];
	    term2 = drot[4] * statev[23] + drot[7] * statev[21] + drot[10] * 
		    statev[24];
	    term3 = drot[4] * statev[25] + drot[7] * statev[24] + drot[10] * 
		    statev[22];
	    term4 = drot[5] * statev[20] + drot[8] * statev[23] + drot[11] * 
		    statev[25];
	    term5 = drot[5] * statev[23] + drot[8] * statev[21] + drot[11] * 
		    statev[24];
	    term6 = drot[5] * statev[25] + drot[8] * statev[24] + drot[11] * 
		    statev[22];
	    term7 = drot[6] * statev[20] + drot[9] * statev[23] + drot[12] * 
		    statev[25];
	    term8 = drot[6] * statev[23] + drot[9] * statev[21] + drot[12] * 
		    statev[24];
	    term9 = drot[6] * statev[25] + drot[9] * statev[24] + drot[12] * 
		    statev[22];
	    statev[20] = term1 * drot[4] + term2 * drot[7] + term3 * drot[10];
	    statev[21] = term4 * drot[5] + term5 * drot[8] + term6 * drot[11];
	    statev[22] = term7 * drot[6] + term8 * drot[9] + term9 * drot[12];
	    statev[23] = term1 * drot[5] + term2 * drot[8] + term3 * drot[11];
	    statev[24] = term4 * drot[6] + term5 * drot[9] + term6 * drot[12];
	    statev[25] = term1 * drot[6] + term2 * drot[9] + term3 * drot[12];
	}
    }

/* ---- compute effective strain rate */

    if (*dtime != 0.f) {
	dum = 0.f;
	i__1 = *ntens;
	for (i__ = 4; i__ <= i__1; ++i__) {
/* L40: */
/* Computing 2nd power */
	    d__1 = dstran[i__];
	    dum += d__1 * d__1;
	}
/* Computing 2nd power */
	d__1 = dstran[1];
/* Computing 2nd power */
	d__2 = dstran[2];
/* Computing 2nd power */
	d__3 = dstran[3];
	ddd = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + dum * .5f) * 
		con1 / *dtime;
    } else {
	ddd = 0.f;
    }

/* ---- calculate trial alpha, kappa and yield radius */

    alphaxx = statev[1] + statev[20];
    alphayy = statev[2] + statev[21];
    alphazz = statev[3] + statev[22];
    alphaxy = statev[4] + statev[23];
    alphayz = statev[5] + statev[24];
    alphazx = statev[6] + statev[25];
/* Computing 2nd power */
    d__1 = alphaxx;
/* Computing 2nd power */
    d__2 = alphayy;
/* Computing 2nd power */
    d__3 = alphazz;
/* Computing 2nd power */
    d__4 = alphaxy;
/* Computing 2nd power */
    d__5 = alphayz;
/* Computing 2nd power */
    d__6 = alphazx;
    alpm = con1 * sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + (d__4 * d__4 
	    + d__5 * d__5 + d__6 * d__6) * 2.f);
    d__1 = dcs0 / dcs;
    sto = *dtime * rs1 * alpm * pow_dd(&d__1, &zz);
    d__1 = dcs0 / dcs;
    sto2 = *dtime * rs2 * pow_dd(&d__1, &zz);
    if (iparam1 == 0) {
	d__1 = dcs0 / dcs;
	sto = *dtime * (rs1 + rd1 * ddd + rs3 + rd3 * ddd) * alpm * pow_dd(&
		d__1, &zz);
	d__1 = dcs0 / dcs;
	sto2 = *dtime * (rs2 + rd2 * ddd) * pow_dd(&d__1, &zz);
    }
    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L50: */
	statev[i__] *= 1.f - sto;
    }
    i__1 = *ntens + 19;
    for (i__ = 19; i__ <= i__1; ++i__) {
/* L51: */
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
	    ytheta + statev[7]) * dam1;

/* ---- calculate trial elastic deviatoric stresses */

    for (i__ = 1; i__ <= 3; ++i__) {
/* L60: */
	stress[i__] = dam2 * (stress[i__] - pold) + dam1 * twog * (dstran[i__]
		 - davg);
    }
    i__1 = *ntens;
    for (i__ = 4; i__ <= i__1; ++i__) {
/* L70: */
	stress[i__] = dam2 * stress[i__] + dam1 * g * dstran[i__];
    }

/* ---- compute xi (deviatoric stress - 2/3 alpha) */

    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L80: */
	xi[i__ - 1] = stress[i__] - twothd * statev[i__];
    }

/* ---- compute (magnitude of xi) squared */

    dum = 0.f;
    i__1 = *ntens;
    for (i__ = 4; i__ <= i__1; ++i__) {
/* L85: */
/* Computing 2nd power */
	d__1 = xi[i__ - 1];
	dum += d__1 * d__1;
    }
/* Computing 2nd power */
    d__1 = xi[0];
/* Computing 2nd power */
    d__2 = xi[1];
/* Computing 2nd power */
    d__3 = xi[2];
    ximag2 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + dum * 2.f;

/* ---- check for plasticity */

    ak2 = ximag2 - twothd * ak * abs(ak);
    if (ak2 <= 0.f || ddd == 0.f) {
	goto L300;
    }

/* ---- plasticity process begins here */

    ximag = sqrt(ximag2);

/* ---- return trial stresses to yield surface, add pressure term */
/*      and update state variables */

    if (iparam1 == 0) {
	d__1 = dcs0 / dcs;
	dgam = (ximag - con1 * ak) / (dam1 * twog + twothd * pow_dd(&d__1, &
		zz) * (h1 + h3 + h2 * dam1));
    }
    dgam2 = dgam / ximag;

    dsig = dam1 * twog * dgam2;
    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L90: */
	stress[i__] -= dsig * xi[i__ - 1];
    }
    stress[1] += p;
    stress[2] += p;
    stress[3] += p;

    d__1 = dcs0 / dcs;
    statev[7] += dgam * con1 * h2 * pow_dd(&d__1, &zz);
    statev[7] = max(0.,statev[7]);

    dalph = (h1 + h3) * dgam2;
    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L100: */
	d__1 = dcs0 / dcs;
	statev[i__] += dalph * xi[i__ - 1] * pow_dd(&d__1, &zz);
    }

/* ---- update plastic strain (for output purposes only) */

    statev[9] += dgam * con1;

/* ---- update temperature for adiabatic problems */

    dum = 0.f;
    i__1 = *ntens;
    for (i__ = 4; i__ <= i__1; ++i__) {
/* L110: */
	dum += stress[i__] * xi[i__ - 1];
    }
    statev[8] += htcp * dgam2 * (stress[1] * xi[0] + stress[2] * xi[1] + 
	    stress[3] * xi[2] + dum * 2.f);

/* ---- update damage */

    epsdot = dgam * con1 / *dtime;
/* Computing MAX */
    d__3 = stress[1] - stress[2];
    d__4 = stress[2] - stress[3];
    d__5 = stress[3] - stress[1];
    d__1 = 1e-15f, d__2 = sqrt((pow_dd(&d__3, &c_b2) + pow_dd(&d__4, &c_b2) + 
	    pow_dd(&d__5, &c_b2) + (pow_dd(&stress[4], &c_b2) + pow_dd(&
	    stress[5], &c_b2) + pow_dd(&stress[6], &c_b2)) * 6) * .5f);
    sige = max(d__1,d__2);

/* Cocks-Ashby large pore growth term */
    cacon = (d__1 = vtheta / ytheta, abs(d__1));
    if (cacon < 8.f) {
	cacon = 8.f;
    }
    dterm = (cacon * 2 - 1) * 2 / (cacon * 2 + 1);
/* Computing MIN */
    d__1 = 15.f, d__2 = p * dterm / sige;
    arg = min(d__1,d__2);
    beta = sinh((max(0.,arg)));
    c90 = cacon + 1.f;
/* Computing MIN */
    d__1 = 15.f, d__2 = beta * *dtime * epsdot * c90;
    psi = min(d__1,d__2);
/* Computing MAX */
    d__1 = 0.f, d__2 = (pow_dd(&phi1, &c90) - 1.f) * exp(psi) + 1.f;
    tmp = max(d__1,d__2);
/* Computing MIN */
    d__2 = 1.f / c90;
    d__1 = 1.f - pow_dd(&tmp, &d__2);
    statev[18] = min(d__1,.99);
/* Cocks-Ashby void growth rate */
    d__1 = 1.f - statev[18];
    d__2 = vtheta / ytheta;
    statev[19] = beta * epsdot * (1.f / pow_dd(&d__1, &d__2) - (1.f - statev[
	    18]));
/* McClintock form of void growth */
    abc = pow_dd(&c_b19, &c_b54) / ((1.f - props[30]) * 2.f) * sinh(pow_dd(&
	    c_b19, &c_b54) * (1.f - props[30]) * p / sige);
    vrad = props[31] * exp(statev[9] * abc / con1);
    statev[10] = pi * pow_dd(&vrad, &c_b2);
    statev[11] = statev[10] * 3.f * abc * epsdot;
/* Nucleation of voids */
/* deviatoric stress */
    ds11 = stress[1] - p;
    ds22 = stress[2] - p;
    ds33 = stress[3] - p;
    ds12 = stress[4];
    ds23 = stress[5];
    ds13 = stress[6];
/* invariants of stress */
    di1 = p * 3.f;
    dj2 = (pow_dd(&ds11, &c_b2) + pow_dd(&ds22, &c_b2) + pow_dd(&ds33, &c_b2) 
	    + (pow_dd(&ds12, &c_b2) + pow_dd(&ds23, &c_b2) + pow_dd(&ds13, &
	    c_b2)) * 2) * .5f;
    dj3 = ds11 * (ds22 * ds33 - pow_dd(&ds23, &c_b2)) - ds22 * (ds11 * ds33 - 
	    pow_dd(&ds13, &c_b2)) + ds33 * (ds22 * ds11 - pow_dd(&ds12, &c_b2)
	    );
    if (dj2 <= 0.f) {
	r1 = 0.f;
	r2 = 0.f;
	r3 = 0.f;
    } else {
	r1 = .14814814814814814f - pow_dd(&dj3, &c_b2) / pow_dd(&dj2, &c_b19);
	r2 = dj3 / pow_dd(&dj2, &c_b20);
	r3 = di1 / pow_dd(&dj2, &c_b54);
    }
    r3 = abs(r3);
/* 	   if(r3.lt.0)then */
/*            r3=-r3 */
/*           endif */
    zzz = props[32] * r1 + props[33] * r2 + props[34] * r3;
    zzz = abs(zzz);
    zzzz = pow_dd(&props[37], &c_b54) / (props[36] * pow_dd(&props[38], &
	    third)) * zzz;
    statev[17] = statev[13];
    statev[13] = props[35] * exp(statev[9] * zzzz / con1) * exp(-props[51] / 
	    statev[8]);
/* added for nonmonotonic path sequences, statev(17) is old nucleation */
    if (statev[13] < statev[17]) {
	statev[13] = (d__1 = statev[17] + statev[13], abs(d__1));
    }
/* Coalescence factor */
    d__1 = dcs0 / dcs;
    cf = (cd1 + cd2 * statev[13] * statev[10]) * exp(props[52] * statev[8]) * 
	    pow_dd(&d__1, &zz);
/* Damage */
    damage = cf * (statev[13] * statev[10] + statev[18]);
    if (damage > .1f) {
	damage = .99f;
    }
    statev[14] = min(damage,.99);
/* Nucleation Rate */
    epsdot = abs(epsdot);
    statev[15] = zzzz * statev[13] * epsdot;
/* Damage Rate */
    zsecond = cf * (statev[15] * statev[10] + statev[13] * statev[11] + 
	    statev[19]);
    d__1 = dcs0 / dcs;
    zthird = (statev[13] * statev[10] + statev[18]) * cd2 * pow_dd(&d__1, &zz)
	     * exp(props[52] * statev[8]) * (statev[15] * statev[10] + statev[
	    13] * statev[11]);
    statev[16] = zsecond + zthird;

/* Triaxiality */
    statev[12] = p / sige;

/* ---- form elastic-plastic constitutive matrix */
/*     if bbb=1 then conventional stiffness is calculated */
/*     otherwise simo stiffness is calculated */

    if (iparam2 == 0) {
/* Computing 2nd power */
	d__1 = alphaxx;
/* Computing 2nd power */
	d__2 = alphayy;
/* Computing 2nd power */
	d__3 = alphazz;
	alpm2 = con1 * sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + (pow_dd(
		&alphaxy, &c_b2) + pow_dd(&alphayz, &c_b2) + pow_dd(&alphazx, 
		&c_b2)) * 2.f);
	bbb = (con1 * (ak + h2 * con1 * dgam) + alpm2 - alpm) / ximag;
	bbb = min(1.,bbb);
	slope = h1 + h3 + h2 - (rd1 + rd3) * alpm2 * alpm2 - rd2 * statev[7] *
		 statev[7];
/* Computing MAX */
	d__1 = (h1 + h2 + h3) * 1e-5f;
	slope = max(d__1,slope);
	gamma = 1.f / (slope / 3.f / g + 1.f) - (1.f - bbb);
/* Computing 2nd power */
	d__1 = ximag;
	r__ = dam1 * (twog * gamma / (d__1 * d__1));
	z1 = dam1 * (blk + twothd * twog * bbb);
	z2 = dam1 * (blk - third * twog * bbb);
	z3 = dam1 * (g * bbb);
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
    }
/* L2000: */
    goto L400;

/* ---- elastic process begins here */

L300:
    for (i__ = 1; i__ <= 3; ++i__) {
/* L310: */
	stress[i__] += p;
    }
    dgam = 0.f;

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
    z1 = dam1 * (blk + twothd * twog);
    z2 = dam1 * (blk - third * twog);
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
	ddsdde[i__ + i__ * ddsdde_dim1] = dam1 * g;
    }

/* ---- clean up */

L400:
    statev[7] = max(statev[7],0.);

    return 0;
} /* umat_ */

#endif /* __F2C__ */
