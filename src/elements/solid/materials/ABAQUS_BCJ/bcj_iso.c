/* $Id: bcj_iso.c,v 1.3 2004-08-01 20:42:35 paklein Exp $ */ 
/* bcj_iso.f -- translated by f2c (version 20030320).
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
static integer c__40 = 40;
static integer c__9 = 9;
static integer c__1 = 1;
static integer c__5 = 5;
static integer c__3 = 3;
static integer c__2 = 2;

/* Subroutine */ int bcj_iso_(doublereal *stress, doublereal *statev, doublereal 
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

    static doublereal zero = 0.;
    static doublereal un = 1.;
    static integer mpi = 6;

    /* Format strings */
    static char fmt_200[] = "(/1x,\002-------  SOLUTION ELASTIQUE  -------"
	    "-\002/)";
    static char fmt_201[] = "(/1x,\002-------  SOLUTION PLASTIQUE  -------"
	    "-\002/)";

    /* System generated locals */
    integer ddsdde_dim1, ddsdde_offset, i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    extern /* Subroutine */ int predelas_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, integer *, integer *, integer *);
    static doublereal e, g;
    static integer m;
    static doublereal v, alphaelas[6];
    extern /* Subroutine */ int elplastdo_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, logical *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *);
    static doublereal dc, hi, hk, pf, pr, xk, stresselas[6];
    static integer ncd;
    static doublereal rdi, phi, rdk;
    static integer nvi;
    static doublereal rsi, rsk, xip, xms, xnu, dlan, epsd;
    extern /* Subroutine */ int aset_(doublereal *, doublereal *, integer *);
    static doublereal teta, epsf, epsr;
    static integer irot;
    static doublereal epsx, alpha[6], eplas[6], yield;
    static integer nplan;
    extern /* Subroutine */ int recov_(doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *);
    static integer natur;
    extern /* Subroutine */ int store_(doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    static doublereal stvit[40], prelas, criter, vnewdt, dlambda, xlambda;
    extern /* Subroutine */ int calcrit_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *), solelas_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *)
	    ;
    static logical lonewdt;
    integer solve_done = 0;

    /* Fortran I/O blocks */
    static cilist io___42 = { 0, 0, 0, fmt_200, 0 };
    static cilist io___43 = { 0, 0, 0, fmt_201, 0 };


/* ====================================================================== */
/* ---------------------------------------------------------------------- */
/* -------------  MODELISATION D'UNE LOI DE COMPORTEMENT  --------------- */
/* -----  ELASTOPLASTIQUE AVEC ECROUISSAGE ISOTROPE ET CINEMATIQUE  ----- */
/* -------------------  COUPLEE A L'ENDOMMAGEMENT  ---------------------- */
/* ---------------------------------------------------------------------- */
/* ====================================================================== */
/* ------------- */
/*  EN ENTREE : */
/* ------------- */
/*  STRAN	: TABLEAU CONTENANT LES COMPOSANTES DE LA DEFORMATION */
/* 		  TOTALE EN DEBUT D'INCREMENT. */
/*  DSTRAN	: TABLEAU DES INCREMENTS DE DEFORMATION. */
/*  TIME(1)	: VALEUR DU PAS DE TEMPS EN DEBUT D'INCREMENT. */
/*  TIME(2)	: VALEUR DU TEMPS TOTAL EN DEBUT D'INCREMENT. */
/*  DTIME	: INCREMENT DE TEMPS. */
/*  TEMP		: TEMPERATURE EN DEBUT D'INCREMENT. */
/*  DTEMP	: INCREMENT DE TEMPERATURE. */
/*  PREDEF	: TABLEAU CONTENANT LES VALEURS INTERPOLEES DES */
/* 		  VARIABLES DE CHAMPS PREDEFINIES AU POINT D'INTEGRATION */
/* 		  EN DEBUT D'INCREMENT, BASEES SUR LES VALEURS PRISES */
/* 		  AUX NOEUDS DANS L'OPTION *FIELD. */
/*  DPRED	: TABLEAU CONTENANT LES INCREMENTS DES VARIABLES DE */
/* 		  CHAMPS PREDEFINIES. */
/*  CMNAME	: NOM DU MATERIAU DONNE PAR L'OPTION *MATERIAL. */
/*  NDI		: NOMBRE DE COMPOSANTES DES CONTRAINTES DIRECTES AU */
/* 		  POINT D'INTEGRATION. */
/*  NSHR		: NOMBRE DE COMPOSANTES DES CONTRAINTES TANGENTIELLES */
/* 		  POINT D'INTEGRATION. */
/*  NTENS	: TAILLE DU TABLEAU DES CONTRAINTES OU DEFORMATIONS */
/* 		  (NDI + NSHR). */
/*  NSTATV	: NOMBRE DE VARIABLES D'ETAT ASSOCIEES DEPENDANT DE LA */
/* 		  SOLUTION (TEL QU'IL EST DEFINI DANS L'OPTION *DEPVAR). */
/*  PROPS	: TABLEAU DES CONSTANTES CARACTERISTIQUES ENTREES DANS */
/* 		  L'OPTION *USER MATERIAL POUR CE MATERIAU. */
/*  NPROPS	: NOMBRE DES CONSTANTES CARACTERISTIQUES (LA VALEUR */
/* 		  EST DONNEE DANS LE PARAM£TRE 'CONSTANTS' DE L'OPTION */
/* 		  *USER MATERIAL). */
/*  COORDS(3)	: TABLEAU CONTENANT LES COORDONNEES COURANTES DU POINT. */
/*  DROT(3,3)	: MATRICE INCREMENTALE DE ROTATION (RIGIDE). */
/*  CELENT	: LONGUEUR D'ELEMENT CARACTERISQUE. */
/*  DFGRD0	: TABLEAU (3,3) CONTENANT LES GRADIENTS DE DEFORMATION */
/* 		  EN DEBUT D'INCREMENT. */
/*  DFGRD1	: TABLEAU (3,3) CONTENANT LES GRADIENTS DE DEFORMATION */
/* 		  EN FIN D'INCREMENT. */
/*  NOEL		: NUMERO D'ELEMENT. */
/*  NPT		: NUMERO DU POINT D'INTEGRATION. */
/*  LAYER	: */
/*  KSPT		: */
/*  KSTEP	: NUMERO DU PAS. */
/*  KINC		: NUMERO D'INCREMENT. */
/* ----------- */
/*  LOCALES : */
/* ----------- */
/*  EPLAS	: TENSEUR DES DEFORMATIONS PLASTIQUES */
/*  ALPHA(6)	: TENSEUR DES CONTRAINTES CINEMATIQUES */
/*  XLAMBDA	: DEFORMATION PLASTIQUE CUMULEE */
/*  PR		: VARIABLE D'ECROUISSAGE ISOTROPE */
/*  PHI		: DAMAGE */
/*  E		: MODULE D'YOUNG */
/*  XNU		: COEFFICIENT DE POISSON */
/*  G		: MODULE DE CISAILLEMENT */
/*  XK		: COEFFICIENT DE LAME */
/*  YIELD	: LIMITE D'ELASTICITE */
/*  V		: CONSTANTE VISCOPLASTIQUE */
/*  PF		: CONSTANTE VISCOPLASTIQUE */
/*  RD		: CONSTANTE D'ECROUISSAGE ISOTROPE */
/*  H		: CONSTANTE D'ECROUISSAGE ISOTROPE */
/*  RS		: CONSTANTE D'ECROUISSAGE ISOTROPE (RECOVERY) */
/*  PRD		: CONSTANTE D'ECROUISSAGE CINEMATIQUE */
/*  PH		: CONSTANTE D'ECROUISSAGE CINEMATIQUE */
/*  PRS		: CONSTANTE D'ECROUISSAGE CINEMATIQUE (RECOVERY) */
/*  XLD		: LENGTH SCALE PARAMETER */
/*  XFN		: VOLUME FRACTION */
/*  XKIC		: ROUGHNESS PARAMETER */
/*  XLF		: VOLUME FRACTION */
/*  XDN		: CONSTANTE XLD*XFN/(XKIC*XLF**(1/3)) */
/*  GS		: DUCTILE DAMAGE PARAMETERS */
/*  PS		: DUCTILE DAMAGE PARAMETERS */
/*  DC		: VALEUR D'ENDOMAGEMMENT CRITIQUE */
/*  CTE		: CONSTANTE DE TENSION */
/*  CCO		: CONSTANTE DE COMPRESSION */
/*  CTO		: CONSTANTE DE TORSION */
/*  VM		: VOID GROWTH INDICE M */
/*  DI		: VOID IMPINGEMENT CRITICAL DIAMETER RATIO */
/*  DS		: VOID SHEET CRITICAL DIAMETER RATIO */
/*  DC		: VALEUR DE L'ENDOMMAGEMENT CRITIQUE --> RUPTURE */
/*  TETA		: TETA-METHODE */
/*  EPSF		: PRECISION DE LA FONCTION DE CHARGE F */
/*  EPSG		: PRECISION DE LA FONCTION D'ENDOMMAGEMENT G */
/*  KCT		: KCT=0 --> MATRICE TANGENTE CONTINUE */
/*  		  KCT=1 --> MATRICE TANGENTE CONSISTENTE */
/*  NPLAN	: NPLAN=1 --> DEFORMATIONS PLANES */
/* 		  NPLAN=2 --> CONTRAINTES PLANES */
/*  NATUR	: NATUR=0 --> PETITES DEFORMATIONS */
/* 		  NATUR=1 --> GRANDES DEFORMATIONS */
/*  IROT		: IROT=0 --> PAS DE ROTATION */
/* 		  IROT=1 --> AVEC ROTATION */
/* 		: IROT=2 --> AVEC DERIVEE COROTATIONNELLE DE JAUMANN */
/* 		  IROT=3 --> AVEC ROTATION D'ABAQUS (ROTSIG) */
/*  NVI		: NOMBRE DE VARIABLE INTERNES STOCKEES (EE,R,X,J2(X),Y,...) */
/*  M		: COMMENTAIRES (PARAMETRE D'IMPRESSION) */
/*  MPI		: UNITE DE FICHIER DE SORTIE */
/* --------------------------------------- */
/*  VARIABLES DEVANT ETRE REACTUALISEES : */
/* --------------------------------------- */
/*  STRESS	: TENSEUR DES CONTRAINTES. */
/*  STATEV	: TABLEAU CONTENANT LES VARIABLES D'ETAT DEPENDANT DE */
/* 		  LA SOLUTION ET DEFINIES PAR L'UTILISATEUR. LA TAILLE */
/* 		  DE CE TABLEAU EST DEFINIE PAR L'OPTION *DEPVAR. */
/* 		  See Below for more information */
/*  SSE		: ENERGIE SPECIFIQUE ELASTIQUE DE DEFORMATION. */
/*  SPD		: DISSIPATION PLASTIQUE. */
/*  SCD		: "CREEP" DISSPATION */
/* ------------- */
/*  EN SORTIE : */
/* ------------- */
/*  DDSDDE	: MATRICE JACOBIENNE (NTENS,NTENS) DU MODELE THEORIQUE (DERIVEE */
/* 		  DE L'INCREMENT DE DEFORMATION PAR L'INCREMENT DES CONTRAINTES). */
/*  PNEWDT	: COEFFICIENT DU NOUVEAU TEMPS D'INCREMENT ALLANT ETRE UTILISE. */
/* 		  LES VARIABLES SUIVANTES SONT SEULEMENT DEFINIES SI ELLES SONT */
/* 		  UTILISEES AVEC L'ANALYSE TEMPERATURE-DEPLACEMENT TOTALEMENT COUPLEE. */
/*  RPL		: GENERATION VOLUMETRIQUE DE CHALEUR PAR UNITE DE TEMPS EN FIN */
/* 		  D'INCREMENT CAUSEE PAR LE TRAVAIL MECANIQUE DU MATERIAU. */
/*  DDSDDT	: VARIATION DES INCREMENTS DE CONTRAINTES PAR RAPPORT A LA */
/* 		  TEMPERATURE. */
/*  DRPLDE	: VARIATION DE RPL PAR RAPPORT AUX INCREMENTS DE DEFORMATION. */
/*  DRPLDT	: VARIATION DE RPL PAR RAPPORT A L'INCREMENT DE TEMPERATURE. */


/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* ------  Internal State Variable Array STATEV : */

/* --  Internal state variables required for stress integration: */

/* 	- Plastic Strain Tensor Ep */
/* 	- Kinematic Hardening (Deformation) Alpha (b in BCJ Model) */
/* 	- Isotropic Hardening (Deformation) r (R in BCJ Model) */
/* 	- Anisotropic Damage Tensor Endo (Phi in BCJ Model) */
/* 	- Plastic Equivalent Strain Rate Dlan at time t */
/* 	- Total Strain only for Green-Naghdi derivative (Irot = 2) */

/* --  Internal state variables for information: */

/* 	- Plastic Equivalent Strain Xlambda (PEEQ) */
/* 	- Stored Variables if Nvi > 0: */
/* 		- Isotropic Hardening R (Kappa in BCJ Model) */
/* 		- Kinematic Hardening Norm J2(X) (J2(Alpha) in BCJ Model) */
/* 		- Energy Release Rate Y */
/* 		- Triaxiality Yt */
/* 		- Reached interation number per increment */
/* 		- Nucleation Tensor Eta */
/* 		- Void Growth V */
/* 		- Coalescence Tensor c */

/*  *-------------------------------*-----------------*-----------------* */
/*  |     Component Definitions     |       2D        |       3D        | */
/*  *-------------------------------*-----------------*-----------------* */
/*  |      Plastic Strain Ep11      |    Statev(1)    |    Statev(1)    | */
/*  |      Plastic Strain Ep22      |    Statev(2)    |    Statev(2)    | */
/*  |      Plastic Strain Ep33      |    Statev(3)    |    Statev(3)    | */
/*  |      Plastic Strain Ep12      |    Statev(4)    |    Statev(4)    | */
/*  |      Plastic Strain Ep13      |                 |    Statev(5)    | */
/*  |      Plastic Strain Ep23      |                 |    Statev(6)    | */
/*  *-------------------------------*-----------------*-----------------* */
/*  |  Kinematic Hardening Alpha11  |    Statev(5)    |    Statev(7)    | */
/*  |  Kinematic Hardening Alpha22  |    Statev(6)    |    Statev(8)    | */
/*  |  Kinematic Hardening Alpha33  |    Statev(7)    |    Statev(9)    | */
/*  |  Kinematic Hardening Alpha12  |    Statev(8)    |    Statev(10)   | */
/*  |  Kinematic Hardening Alpha13  |                 |    Statev(11)   | */
/*  |  Kinematic Hardening Alpha23  |                 |    Statev(12)   | */
/*  *-------------------------------*-----------------*-----------------* */
/*  |   Plastic Equivalent Strain   |    Statev(9)    |    Statev(13)   | */
/*  *-------------------------------*-----------------*-----------------* */
/*  |     Isotropic Hardening r     |    Statev(10)   |    Statev(14)   | */
/*  *-------------------------------*-----------------*-----------------* */
/*  |      Istropic Damage Phi      |    Statev(11)   |    Statev(15)   | */
/*  *-------------------------------*-----------------*-----------------* */
/*  |           Precedent           |                 |                 | */
/*  |   Plastic Equivalent Strain   |    Statev(12)   |    Statev(16)   | */
/*  |          Rate  Dl/Dt          |                 |                 | */
/*  *-------------------------------*-----------------*-----------------* */
/*  | ------  Extra Internal State Variables for Post-Processing ------ | */
/*  *-------------------------------*-----------------*-----------------* */
/*  | Isotropic Hardening R (Kappa) |    Statev(13)   |    Statev(17)   | */
/*  *-------------------------------*-----------------*-----------------* */
/*  |   Kinematic Hardening Norm    |    Statev(14)   |    Statev(18)   | */
/*  |             J2(X)             |                 |                 | */
/*  *-------------------------------*-----------------*-----------------* */
/*  | Back Stress (Kin. Hard.) X11  |    Statev(15)   |    Statev(19)   | */
/*  | Back Stress (Kin. Hard.) X22  |    Statev(16)   |    Statev(20)   | */
/*  | Back Stress (Kin. Hard.) X33  |    Statev(17)   |    Statev(21)   | */
/*  | Back Stress (Kin. Hard.) X12  |    Statev(18)   |    Statev(22)   | */
/*  | Back Stress (Kin. Hard.) X13  |                 |    Statev(23)   | */
/*  | Back Stress (Kin. Hard.) X23  |                 |    Statev(24)   | */
/*  *-------------------------------*-----------------*-----------------* */
/*  |        Triaxiality Yt         |    Statev(19)   |    Statev(26)   | */
/*  *-------------------------------*-----------------*-----------------* */
/*  |   Reached Iteration Number    |    Statev(20)   |    Statev(27)   | */
/*  *-------------------------------*-----------------*-----------------* */
/*  |       Effective Strain        |    Statev(21)   |    Statev(28)   | */
/*  *-------------------------------*-----------------*-----------------* */

/* ------  END Internal State Variable Array STATEV */
/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* ======================================================================= */
/* ----------------------------------------------------------------------* */
/* .1-----  Include,Precision */
/*      INCLUDE 'ABA_PARAM.INC' */
/* .2-----  Parameter,Integer */
/* .3-----  Dimension */
/* .4-----  Real,Complex,Double precision,Logical,Character */
/* .5-----  Common */
/* .6-----  Equivalence */
/* .7-----  Save,Data */
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
/* .8-----  Definition de fonctions */
/* ----------------------------------------------------------------------* */
/* ======================================================================= */
/* ----------------------------------------------------------------------* */
/* -------  Recuperation de la solution a l'increment precedent, des */
/* 	  caracteristiques du materiau et des parametres de convergence. */

    recov_(&props[1], nprops, &statev[1], nstatv, &stress[1], &stran[1], &
	    dstran[1], eplas, alpha, &drot[4], &dfgrd0[4], &dfgrd1[4], &
	    xlambda, &pr, &phi, &dlan, &time[1], dtime, temp, dtemp, &e, &xnu,
	     &g, &xk, &yield, &v, &pf, &hi, &rdi, &rsi, &hk, &rdk, &rsk, &xip,
	     &xms, &dc, &teta, &epsf, &epsx, &epsr, &epsd, &vnewdt, &ncd, &
	    nplan, &natur, &irot, ntens, ndi, nshr, &nvi, &m, &mpi);

/* -------  FIN Recuperation de la solution a l'increment precedent, des */
/* 	  caracteristiques du materiau et des parametres de convergence. */
/* ----------------------------------------------------------------------* */
/* ------- Test : D = 1  --->  Pas d'increment et sortir de UMAT */

    if (phi == un) {
	i__1 = *ntens * *ntens;
	aset_(&ddsdde[ddsdde_offset], &zero, &i__1);
	return 0;
    }

/* ---------  FIN  Test : D = 1  --->  Pas d'increment et sortir de UMAT */
/* ----------------------------------------------------------------------* */
/* -------  Prediction elastique */

/* ---------  Calcul elastique */

    predelas_(&dstran[1], &stress[1], alpha, &pr, &phi, stresselas, alphaelas,
	     &prelas, &hi, &rsi, &hk, &rsk, &g, &xk, &nplan, ntens, ndi, nshr,
	     &m, &mpi);
/* -------  FIN Prediction elastique */
/* ----------------------------------------------------------------------* */
/* ------- Calcul du critere */
    calcrit_(&criter, stresselas, alphaelas, &prelas, &phi, &hk, &hi, &yield, 
	    ntens, ndi, nshr, &m, &mpi);
/* ------- FIN Calcul du critere */
/* ----------------------------------------------------------------------* */
/* ------- Test : F < 0  --->  SOLUTION PUREMENT ELASTIQUE */
/* 		F >=0  --->  SOLUTION AVEC CORRECTION PLASTIQUE */

/* ---------  Solution elastique or computing an elastic predictor */
    if (criter < zero || *kinc == -1) {
	if (m >= 3) {
	    io___42.ciunit = mpi;
	    s_wsfe(&io___42);
	    e_wsfe();
	}
	solelas_(&stress[1], &statev[1], &ddsdde[ddsdde_offset], stresselas, 
		alphaelas, &prelas, &phi, &xk, &g, &xip, &ncd, nstatv, ntens, 
		ndi, &nvi, &m, &mpi);
	return 0;
    }
/* ---------  FIN Solution elastique */
/* ----------------------------------------------------------------------* */

/* ----------------------------------------------------------------------* */
/* ---------  Solution plastique */
	solve_done = 0;
	while (!solve_done) 
	{
    	if (m >= 3) {
			io___43.ciunit = mpi;
			s_wsfe(&io___43);
			e_wsfe();
    	}

/* -----------  Calcul de la solution */

		elplastdo_(&ddsdde[ddsdde_offset], &stress[1], &stran[1], &dstran[1], 
			eplas, alpha, &xlambda, &dlambda, &pr, &phi, &dlan, stvit, &yield,
			&v, &pf, &hk, &rdk, &rsk, &hi, &rdi, &rsi, &dc, &xms, &g, &xk, &
			teta, &epsf, &epsx, &epsr, &epsd, dtime, &lonewdt, &ncd, &nplan, &
			natur, ntens, ndi, nshr, &c__40, &m, &mpi);

/* ---------  FIN Solution plastique */
/* ----------------------------------------------------------------------* */
/* -----------  Cas D > 1 ---> Recommencer l'increment avec DTIME*PNEWDT */

		/* signal to keep current time increment */
		*pnewdt = *dtime;
		
		/* found state (solved or unsolved) where damage exceeds max */
		if (ncd != 0 && phi > dc) {

			/* re-initialize the state */
		    recov_(&props[1], nprops, &statev[1], nstatv, &stress[1], &stran[1], &
				dstran[1], eplas, alpha, &drot[4], &dfgrd0[4], &dfgrd1[4], &
				xlambda, &pr, &phi, &dlan, &time[1], dtime, temp, dtemp, &e, &xnu,
				&g, &xk, &yield, &v, &pf, &hi, &rdi, &rsi, &hk, &rdk, &rsk, &xip,
				&xms, &dc, &teta, &epsf, &epsx, &epsr, &epsd, &vnewdt, &ncd, &
				nplan, &natur, &irot, ntens, ndi, nshr, &nvi, &m, &mpi);

			/* freeze damage at the maximum */
			ncd = 0;
			phi = dc;

			/* recompute elastic predictor */
			predelas_(&dstran[1], &stress[1], alpha, &pr, &phi, stresselas, alphaelas,
				&prelas, &hi, &rsi, &hk, &rsk, &g, &xk, &nplan, ntens, ndi, nshr,
				&m, &mpi);
		}
		/* solution failed even though damage has not exceeded max */
    	else if (lonewdt) {

			/* signal step cut */
			*pnewdt *= vnewdt;
			return 0;
    	}
    	else /* exit */
    		solve_done = 1;
    }
    
/* ---------  FIN  Cas D > 1 ---> Recommencer l'increment avec DTIME*PNEWDT */
/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* ---------  Enregistrement de la solution plastique */

    store_(&statev[1], nstatv, eplas, alpha, &xlambda, &pr, &phi, &dlambda, 
	    dtime, stvit, &c__40, ntens, ndi, &nvi, &ncd, &m, &mpi);

/* ---------  FIN Enregistrement de la solution plastique */
/* ----------------------------------------------------------------------* */
/* ======================================================================= */
/* ----------------------------------------------------------------------* */
/* ---------   FORMATS  -------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* ======================================================================= */
    return 0;
} /* umat_ */


/* $Id: bcj_iso.c,v 1.3 2004-08-01 20:42:35 paklein Exp $ */
/* Subroutine */ int addtens_(integer *ntens, doublereal *dt, doublereal *t)
{

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

/* ====================================================================== */
/* ---------------------------------------------------------------------- */
/* ---------  EXECUTION DE L'OPERATION TENSORIELLE: T = T + DT  --------- */
/* ---------------------------------------------------------------------- */
/* ====================================================================== */
/* ------------- */
/*  EN ENTREE : */
/* ------------- */
/*  NTENS 	: DIMENSION DES TENSEURS */
/*  T 		: TENSEUR INITIAL */
/*  DT 		: INCREMENT DE TENSEUR A RAJOUTER A T */
/* ------------- */
/*  EN SORTIE : */
/* ------------- */
/*  T 		: TENSEUR REACTUALISE */
/* ======================================================================= */
/* .1-----  Implicit, External */
/* .2-----  Parameter */
/* .3-----  Dimension */
/* .4-----  Real,Integer,Complex,Double precision,Logical,Character */
/* .5-----  Common */
/* .6-----  Equivalence */
/* .7-----  Save,Data */
    /* Parameter adjustments */
    --t;
    --dt;

    /* Function Body */
/* .8-----  Definition de fonctions */
/* ======================================================================= */
    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
	t[i__] += dt[i__];
    }
/* ======================================================================= */
    return 0;
} /* addtens_ */


/* Subroutine */ int affect_(integer *n, doublereal *v1, doublereal *v2)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

/* ====================================================================== */
/* ---------------------------------------------------------------------- */
/* ---------  EXECUTION DE L'OPERATION TENSORIELLE: T = T + DT  --------- */
/* ---------------------------------------------------------------------- */
/* ====================================================================== */
/* ------------- */
/*  EN ENTREE : */
/* ------------- */
/*  N	 	: DIMENSION DES TENSEURS */
/*  V1 		: TENSEUR INITIAL */
/* ------------- */
/*  EN SORTIE : */
/* ------------- */
/*  V2 		: TENSEUR REACTUALISE */
/* ---------------------------------------------------------------------- */
/* ====================================================================== */
/* ---------------------------------------------------------------------- */
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
    --v2;
    --v1;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	v2[i__] = v1[i__];
    }
/* ======================================================================= */
    return 0;
} /* affect_ */


/* Subroutine */ int calcoeff_(doublereal *zd, doublereal *tn, doublereal *
	zdn, doublereal *zdzd, doublereal *stres1d, doublereal *estar, 
	doublereal *p, doublereal *xkdh, doublereal *dsqhmp, doublereal *
	ryvshdl, doublereal *pold, doublereal *alpha1, doublereal *pr1, 
	doublereal *phi, doublereal *phiold, doublereal *dlambda, doublereal *
	h__, doublereal *unmd, doublereal *unmd1, doublereal *wr, doublereal *
	wx, doublereal *dmwx, doublereal *yield, doublereal *v, doublereal *
	pf, doublereal *g, doublereal *xk, doublereal *hk, doublereal *rdk, 
	doublereal *rsk, doublereal *xms, doublereal *hi, doublereal *rdi, 
	doublereal *rsi, doublereal *tdl, doublereal *untrwx, doublereal *
	phtrwx, doublereal *phtrwx2, doublereal *untrwr, doublereal *htrwr, 
	doublereal *htrwr2, doublereal *trswr, doublereal *rtdl, doublereal *
	teta, integer *ncd, integer *nitmax, integer *ntens, integer *ndi, 
	integer *nshr, integer *m, integer *mpi)
{
    /* Initialized data */

    static doublereal one = 1.;

    /* Format strings */
    static char fmt_101[] = "(1x,6(d20.14,2x))";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void), s_wsfe(cilist *), do_fio(integer *, char *, ftnlen),
	     e_wsfe(void);

    /* Local variables */
    static integer i__;
    static doublereal dphi, dlrd, dtlf, tntn;
    extern /* Subroutine */ int pdtsca_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *), calzdtn_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *);

    /* Fortran I/O blocks */
    static cilist io___55 = { 0, 6, 0, 0, 0 };
    static cilist io___56 = { 0, 6, 0, 0, 0 };
    static cilist io___57 = { 0, 6, 0, 0, 0 };
    static cilist io___58 = { 0, 6, 0, 0, 0 };
    static cilist io___59 = { 0, 6, 0, 0, 0 };
    static cilist io___60 = { 0, 6, 0, 0, 0 };
    static cilist io___61 = { 0, 6, 0, 0, 0 };
    static cilist io___63 = { 0, 6, 0, 0, 0 };
    static cilist io___64 = { 0, 0, 0, 0, 0 };
    static cilist io___66 = { 0, 0, 0, fmt_101, 0 };


/* ====================================================================== */
/* ---------------------------------------------------------------------- */
/* -------------------  CALCUL DES TERMES REPETITIFS  ------------------- */
/* ---------------------------------------------------------------------- */
/* ====================================================================== */
/* ------------- */
/*  EN ENTREE : */
/* ------------- */
/*  ESTAR	: TENSEUR DES DEFORMATIONS CONNUS A L'INSTANT T+DT */
/*  ALPHA1	: TENSEUR DES CONTRAINTES CINEMATIQUES CONNU A T */
/*  UNITDEV	: TENSEUR UNITE DEVIATORIQUE D'ORDRE 4 */
/*  PR1		: VARIABLE D'ECROUISSAGE ISOTROPE CONNU A T */
/*  DLAMBDA	: VITESSE DE DEFORMATION PLASTIQUE CUMULEE */
/*  PHI		: DAMAGE VARIABLE */
/*  H		: DAMAGE FUNCTIONS h1, h2 and h3 */
/*  YIELD	: LIMITE D'ELASTICITE */
/*  V		: CONSTANTE VISCOPLASTIQUE */
/*  PF		: CONSTANTE VISCOPLASTIQUE */
/*  RD		: CONSTANTE D'ECROUISSAGE ISOTROPE */
/*  H		: CONSTANTE D'ECROUISSAGE ISOTROPE */
/*  RS		: CONSTANTE D'ECROUISSAGE ISOTROPE (RECOVERY) */
/*  QS		: CONSTANTE D'ECROUISSAGE ISOTROPE (RECOVERY) */
/*  PRD		: CONSTANTE D'ECROUISSAGE CINEMATIQUE */
/*  PH		: CONSTANTE D'ECROUISSAGE CINEMATIQUE */
/*  PRS		: CONSTANTE D'ECROUISSAGE CINEMATIQUE (RECOVERY) */
/*  VM		: VOID GROWTH INDICE M */
/*  TETA		: TETA-METHODE */
/*  NDI		: NOMBRE DE COMPOSANTES DES CONTRAINTES DIRECTES */
/*  NSHR		: NOMBRE DE COMPOSANTES DES CONTRAINTES TANGENTIELLES */
/*  NTENS	: TAILLE DU TABLEAU DES CONTRAINTES OU DEFORMATIONS */
/* 		  (NDI + NSHR) */
/*  M		: COMMENTAIRES */
/*  MPI		: UNITE DE FICHIER DE SORTIE */
/* ------------- */
/*  EN SORTIE : */
/* ------------- */
/*  ZD		: TENSEUR DE DIRECTION NORMALE */
/*  TN		: TENSEUR NORMAL */
/*  ZDN		: NORME DU TENSEUR ZD (1ERE METHODE) */
/*  ZDZD		: NORME DU TENSEUR ZD (2NDE METHODE) */
/*  DM		: TENSEUR EFFET DU DOMMAGE */
/*  DMI		: TENSEUR EFFET DU DOMMAGE INVERSE */
/*  DMWX		: NORME DU TENSEUR Mt:(Alpha1+Tdl*N) */
/*  STRESS1D	: TENSEUR DEVIATORIQUE DES CONTRAINTES CONNU A T */
/*  UNMD		: EXEMPLE D'UN COEFFICIENT */
/*  NITMAX	: NOMBRE D'ITERATIONS DU SCHEMA DE NEWTON */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
/* .1-----  Implicit, External */
/* .2-----  Parameter,Integer */
/* .3-----  Dimension */
/* .----- */
/* .4-----  Real,Complex,Double precision,Logical,Character */
/* .5-----  Common */
/* .6-----  Equivalence */
/* .7-----  Save,Data */
    /* Parameter adjustments */
    --h__;
    --alpha1;
    --estar;
    --stres1d;
    --tn;
    --zd;

    /* Function Body */
/* .8-----  Definition de fonctions */
/* ---------------------------------------------------------------------* */
/* ====================================================================== */
/* ---------------------------------------------------------------------* */
/* -------  Number of iterations */

    ++(*nitmax);

/* -------  END Number of iterations */
/* ---------------------------------------------------------------------* */
/* -------  Calculation of coefficients */

    *unmd = one - *phi;
    *unmd1 = one / *unmd;
    dphi = *phi - *phiold;
/* ------- */
    *tdl = *teta * *dlambda;
    dlrd = *tdl * *unmd1;
/* ------- */
    dtlf = *dlambda / *pf;
/* ------- */
    *untrwx = one + *teta * (*dlambda * *rdk + *rsk) * *wx;
    *phtrwx = *hk / *untrwx;
    *phtrwx2 = *phtrwx / *untrwx;
/* ------- */
    *untrwr = one + *teta * (*dlambda * *rdi + *rsi) * *wr;
    *htrwr = *hi / *untrwr;
    *htrwr2 = *htrwr / *untrwr;
    *trswr = *teta * *rsi * *wr;
    *rtdl = *pr1 + *tdl;

/* -------  END  Calculation of coefficients */
/* ---------------------------------------------------------------------* */
/* -------  Tensor TN */

    calzdtn_(&zd[1], &tn[1], zdn, zdzd, dsqhmp, ryvshdl, &stres1d[1], &estar[
	    1], &alpha1[1], p, pr1, phi, &dphi, pold, xkdh, &h__[1], tdl, &
	    dtlf, unmd, dmwx, phtrwx, rtdl, htrwr, yield, v, pf, g, xk, xms, 
	    ncd, ntens, ndi, nshr, m, mpi);

/* -------  END Tensor TN */
/* ---------------------------------------------------------------------* */
/* -------  Comment Check Code */
    if (*m >= 4) {
	s_wsle(&io___55);
	do_lio(&c__9, &c__1, "CALCOEFF TDL     = ", (ftnlen)19);
	do_lio(&c__5, &c__1, (char *)&(*tdl), (ftnlen)sizeof(doublereal));
	e_wsle();
	s_wsle(&io___56);
	do_lio(&c__9, &c__1, "CALCOEFF TETA    = ", (ftnlen)19);
	do_lio(&c__5, &c__1, (char *)&(*teta), (ftnlen)sizeof(doublereal));
	e_wsle();
	s_wsle(&io___57);
	do_lio(&c__9, &c__1, "CALCOEFF DLAMBDA = ", (ftnlen)19);
	do_lio(&c__5, &c__1, (char *)&(*dlambda), (ftnlen)sizeof(doublereal));
	e_wsle();
	s_wsle(&io___58);
	do_lio(&c__9, &c__1, "CALCOEFF PHTRWX  = ", (ftnlen)19);
	do_lio(&c__5, &c__1, (char *)&(*phtrwx), (ftnlen)sizeof(doublereal));
	e_wsle();
	s_wsle(&io___59);
	do_lio(&c__9, &c__1, "CALCOEFF PHTRWX2 = ", (ftnlen)19);
	do_lio(&c__5, &c__1, (char *)&(*phtrwx2), (ftnlen)sizeof(doublereal));
	e_wsle();
	s_wsle(&io___60);
	do_lio(&c__9, &c__1, "CALCOEFF HTRWR   = ", (ftnlen)19);
	do_lio(&c__5, &c__1, (char *)&(*htrwr), (ftnlen)sizeof(doublereal));
	e_wsle();
	s_wsle(&io___61);
	do_lio(&c__9, &c__1, "CALCOEFF HTRWR2  = ", (ftnlen)19);
	do_lio(&c__5, &c__1, (char *)&(*htrwr2), (ftnlen)sizeof(doublereal));
	e_wsle();
	pdtsca_(ndi, nshr, &tn[1], &tn[1], &tntn);
	s_wsle(&io___63);
	do_lio(&c__9, &c__1, "CALCOEFF TNTN    = ", (ftnlen)19);
	do_lio(&c__5, &c__1, (char *)&tntn, (ftnlen)sizeof(doublereal));
	e_wsle();
	io___64.ciunit = *mpi;
	s_wsle(&io___64);
	do_lio(&c__9, &c__1, "CALCOEFF TN", (ftnlen)11);
	e_wsle();
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___66.ciunit = *mpi;
	    s_wsfe(&io___66);
	    do_fio(&c__1, (char *)&tn[i__], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }
/* -------  END Comment Check Code */
/* ----------------------------------------------------------------------- */
/* ---------   FORMATS  -------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* L102: */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
    return 0;
} /* calcoeff_ */


/* Subroutine */ int calcoeffi_(doublereal *zd, doublereal *tn, doublereal *
	zdn, doublereal *zdzd, doublereal *phi, doublereal *stres1d, 
	doublereal *estar, doublereal *p, doublereal *xkdh, doublereal *
	dsqhmp, doublereal *ryvshdl, doublereal *alpha1, doublereal *pold, 
	doublereal *pr1, doublereal *phiold, doublereal *unmd, doublereal *
	unmd1, doublereal *h__, doublereal *wr, doublereal *wx, doublereal *v,
	 doublereal *pf, doublereal *g, doublereal *xk, doublereal *hk, 
	doublereal *rdk, doublereal *rsk, doublereal *hi, doublereal *rdi, 
	doublereal *rsi, doublereal *xms, doublereal *dlambda, doublereal *
	tdl, doublereal *dmwx, doublereal *yield, doublereal *untrwx, 
	doublereal *phtrwx, doublereal *phtrwx2, doublereal *untrwr, 
	doublereal *htrwr, doublereal *htrwr2, doublereal *trswr, doublereal *
	rtdl, doublereal *teta, integer *nitmax, integer *ncd, integer *ntens,
	 integer *ndi, integer *nshr, integer *m, integer *mpi)
{
    /* Initialized data */

    static doublereal one = 1.;

    /* Format strings */
    static char fmt_101[] = "(1x,6(d20.14,2x))";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void), s_wsfe(cilist *), do_fio(integer *, char *, ftnlen),
	     e_wsfe(void);

    /* Local variables */
    static integer i__;
    static doublereal dphi, dtlf;
    extern /* Subroutine */ int calzdtn_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *);

    /* Fortran I/O blocks */
    static cilist io___71 = { 0, 0, 0, 0, 0 };
    static cilist io___73 = { 0, 0, 0, fmt_101, 0 };


/* ====================================================================== */
/* ---------------------------------------------------------------------- */
/* -------------------  CALCUL DES TERMES REPETITIFS  ------------------- */
/* ---------------------------------------------------------------------- */
/* ====================================================================== */
/* ------------- */
/*  EN ENTREE : */
/* ------------- */
/*  ESTAR	: TENSEUR DES DEFORMATIONS CONNUS A L'INSTANT T+DT */
/*  ALPHA1	: TENSEUR DES CONTRAINTES CINEMATIQUES CONNU A T */
/*  UNITDEV	: TENSEUR UNITE DEVIATORIQUE D'ORDRE 4 */
/*  PR1		: VARIABLE D'ECROUISSAGE ISOTROPE CONNU A T */
/*  DLAMBDA	: VITESSE DE DEFORMATION PLASTIQUE CUMULEE */
/*  PHI		: DAMAGE VARIABLE */
/*  H		: DAMAGE FUNCTIONS h1, h2 and h3 */
/*  YIELD	: LIMITE D'ELASTICITE */
/*  V		: CONSTANTE VISCOPLASTIQUE */
/*  PF		: CONSTANTE VISCOPLASTIQUE */
/*  RD		: CONSTANTE D'ECROUISSAGE ISOTROPE */
/*  H		: CONSTANTE D'ECROUISSAGE ISOTROPE */
/*  RS		: CONSTANTE D'ECROUISSAGE ISOTROPE (RECOVERY) */
/*  PRD		: CONSTANTE D'ECROUISSAGE CINEMATIQUE */
/*  PH		: CONSTANTE D'ECROUISSAGE CINEMATIQUE */
/*  PRS		: CONSTANTE D'ECROUISSAGE CINEMATIQUE (RECOVERY) */
/*  VM		: VOID GROWTH INDICE M */
/*  TETA		: TETA-METHODE */
/*  NDI		: NOMBRE DE COMPOSANTES DES CONTRAINTES DIRECTES */
/*  NSHR		: NOMBRE DE COMPOSANTES DES CONTRAINTES TANGENTIELLES */
/*  NTENS	: TAILLE DU TABLEAU DES CONTRAINTES OU DEFORMATIONS */
/* 			  (NDI + NSHR) */
/*  M		: COMMENTAIRES */
/*  MPI		: UNITE DE FICHIER DE SORTIE */
/* ------------- */
/*  EN SORTIE : */
/* ------------- */
/*  ZD		: TENSEUR DE DIRECTION NORMALE */
/*  TN		: TENSEUR NORMAL */
/*  ZDN		: NORME DU TENSEUR ZD (1ERE METHODE) */
/*  ZDZD		: NORME DU TENSEUR ZD (2NDE METHODE) */
/*  DM		: TENSEUR EFFET DU DOMMAGE */
/*  DMI		: TENSEUR EFFET DU DOMMAGE INVERSE */
/*  DMWX		: NORME DU TENSEUR Mt:(Alpha1+Tdl*N) */
/*  STRESS1D	: TENSEUR DEVIATORIQUE DES CONTRAINTES CONNU A T */
/*  DLAMBDA	: VITESSE DE DEFORMATION PLASTIQUE CUMULEE */
/*  NITMAX	: INITIALISATION DU NOMBRE D'ITERATIONS (A ZERO) */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
/* .1-----  Implicit, External */
/* .2-----  Parameter,Integer */
/* .3-----  Dimension */
/* .----- */
/* .4-----  Real,Complex,Double precision,Logical,Character */
/* .5-----  Common */
/* .6-----  Equivalence */
/* .7-----  Save,Data */
    /* Parameter adjustments */
    --h__;
    --alpha1;
    --estar;
    --stres1d;
    --tn;
    --zd;

    /* Function Body */
/* .8-----  Definition de fonctions */
/* ---------------------------------------------------------------------* */
/* ====================================================================== */
/* ---------------------------------------------------------------------* */
/* -------  Number of iterations */

    *nitmax = 0;

/* -------  END  Number of iterations */
/* ----------------------------------------------------------------------* */
/* -------  Calculation of coefficients */

    *tdl = *teta * *dlambda;
/* ------- */
    dtlf = *dlambda / *pf;
/* ------- */
    *untrwx = one + *teta * (*rdk * *dlambda + *rsk) * *wx;
    *phtrwx = *hk / *untrwx;
    *phtrwx2 = *phtrwx / *untrwx;
/* ------- */
    *untrwr = one + *teta * *rsi * *wr;
    *htrwr = *hi / *untrwr;
    *htrwr2 = *htrwr / *untrwr;
    *trswr = *teta * *rsi * *wr;
    *rtdl = *pr1 + *tdl;
/* ---------- */
    *unmd = one - *phi;
    *unmd1 = one / *unmd;
    dphi = *phi - *phiold;

/* -------  END  Calculation of coefficients */
/* ----------------------------------------------------------------------* */
/* -------  Tensors ZD and TN, and norms ZDN and ZDN */

    calzdtn_(&zd[1], &tn[1], zdn, zdzd, dsqhmp, ryvshdl, &stres1d[1], &estar[
	    1], &alpha1[1], p, pr1, phi, &dphi, pold, xkdh, &h__[1], tdl, &
	    dtlf, unmd, dmwx, phtrwx, rtdl, htrwr, yield, v, pf, g, xk, xms, 
	    ncd, ntens, ndi, nshr, m, mpi);

/* -------  END Tensors ZD and TN, norms ZDN and ZDN */
/* ----------------------------------------------------------------------* */
/* -------  Comment Check Code */
    if (*m >= 4) {
	io___71.ciunit = *mpi;
	s_wsle(&io___71);
	do_lio(&c__9, &c__1, "CALCOEFFI TN", (ftnlen)12);
	e_wsle();
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___73.ciunit = *mpi;
	    s_wsfe(&io___73);
	    do_fio(&c__1, (char *)&tn[i__], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }
/* -------  END Comment Check Code */
/* ----------------------------------------------------------------------* */
/* ---------   FORMATS  -------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* L102: */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
    return 0;
} /* calcoeffi_ */


/* Subroutine */ int calcrit_(doublereal *criter, doublereal *stress, 
	doublereal *alpha, doublereal *pr, doublereal *phi, doublereal *hk, 
	doublereal *hi, doublereal *yield, integer *ntens, integer *ndi, 
	integer *nshr, integer *m, integer *mpi)
{
    /* Initialized data */

    static doublereal un = 1.;

    /* Format strings */
    static char fmt_102[] = "(1x,a20,6(d20.14,2x))";
    static char fmt_101[] = "(1x,6(d20.14,2x))";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void),
	     s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    static doublereal sigmxeqv;
    static integer i__;
    static doublereal r__, x[6], unmd;
    extern /* Subroutine */ int sinv_(doublereal *, doublereal *, doublereal *
	    , integer *, integer *);
    static doublereal sigmx[6], sigmxh;

    /* Fortran I/O blocks */
    static cilist io___85 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___86 = { 0, 0, 0, 0, 0 };
    static cilist io___87 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___88 = { 0, 0, 0, fmt_102, 0 };


/* ====================================================================== */
/* ---------------------------------------------------------------------- */
/* ----------------  CALCUL DE LA FONCTION-SEUIL DE MISES  -------------- */
/* ---------------------------------------------------------------------- */
/* ====================================================================== */
/* ------------- */
/*  EN ENTREE : */
/* ------------- */
/*  STRESS	: CONTRAINTES A L'INSTANT N */
/*  ALPHA	: TENSEUR DES CONTRAINTES CINEMATIQUES A L'INSTANT N */
/*  PR		: VARIABLE D'ECROUISSAGE ISOTROPE */
/*  PHI		: DAMAGE */
/*  H		: PARAMETRE D'ECROUISSAGE ISOTROPE */
/*  PH		: PARAMETRE D'ECROUISSAGE CINEMATIQUE */
/*  YIELD	: MODULE D'YOUNG */
/*  NTENS	: LONGUEUR DES TABLEAU ALPHA, STRESS, ... */
/*  NDI		: NOMBRE DE CONTRAINTES DIRECTES (S11,S22,S33) */
/*  NSHR		: NOMBRE DE CONTRAINTES TANGENTIELLES (S12,S13,S23) */
/*  M		: COMMENTAIRES */
/*  MPI		: UNITE DE FICHIER DE SORTIE */
/* ----------- */
/*  LOCALES : */
/* ----------- */
/*  SIGMXT	: CONTRAINTES EFFECTIVES A L'INSTANT N */
/* ------------- */
/*  EN SORTIE : */
/* ------------- */
/*  CRITER	:  FONCTION CRITERE */
/* ----------------------------------------------------------------------* */
/* ======================================================================= */
/* ----------------------------------------------------------------------* */
/* .1-----  Implicit, External */
/* .2-----  Parameter,Integer */
/* .3-----  Dimension */
/* ------- */
/* .4-----  Real,Complex,Double precision,Logical,Character */
/* .5-----  Common */
/* .6-----  Equivalence */
/* .7-----  Save,Data */
    /* Parameter adjustments */
    --alpha;
    --stress;

    /* Function Body */
/* .8-----  Definition de fonctions */
/* ======================================================================= */
/* ----------------------------------------------------------------------* */
/* -------  Second Stress Invariant */

    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__ - 1] = *hk * alpha[i__];
	sigmx[i__ - 1] = stress[i__] - x[i__ - 1];
    }
    sinv_(sigmx, &sigmxh, &sigmxeqv, ndi, nshr);
/* ----------------------------------------------------------------------* */
/* -------  Isotropic Hardening */

    r__ = *hi * *pr;
/* ----------------------------------------------------------------------* */
/* -------  Yield Function */

    unmd = un - *phi;
    *criter = sigmxeqv / unmd - r__ - *yield;
/* ----------------------------------------------------------------------* */
/* -------  Comment Check Code */

    if (*m >= 2) {
	io___85.ciunit = *mpi;
	s_wsfe(&io___85);
	do_fio(&c__1, "CALCRIT CRITERE =", (ftnlen)17);
	do_fio(&c__1, (char *)&(*criter), (ftnlen)sizeof(doublereal));
	e_wsfe();
	if (*m >= 3) {
	    io___86.ciunit = *mpi;
	    s_wsle(&io___86);
	    do_lio(&c__9, &c__1, "CALCRIT SIGMX", (ftnlen)13);
	    e_wsle();
	    i__1 = *ntens;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		io___87.ciunit = *mpi;
		s_wsfe(&io___87);
		do_fio(&c__1, (char *)&sigmx[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	    io___88.ciunit = *mpi;
	    s_wsfe(&io___88);
	    do_fio(&c__1, "CALCRIT SIGMXEQV=", (ftnlen)17);
	    do_fio(&c__1, (char *)&sigmxeqv, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }
/* -------  END Comment Check Code */
/* ----------------------------------------------------------------------* */
/* ---------   FORMATS  -------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
    return 0;
} /* calcrit_ */


/* Subroutine */ int caldde_(doublereal *dlde, doublereal *ddde, doublereal *
	dnde, doublereal *unn, doublereal *tntn, doublereal *uzdnn, 
	doublereal *alde, doublereal *dpnde, doublereal *dzde, doublereal *
	dxide, doublereal *estar, doublereal *alpha1, doublereal *zd, 
	doublereal *tn, doublereal *tetn, doublereal *altn, doublereal *phi, 
	doublereal *zdzd, doublereal *dfdl, doublereal *dfdx, doublereal *
	dfdr, doublereal *dfdd, doublereal *dxdl, doublereal *dxdx, 
	doublereal *dxdd, doublereal *drdl, doublereal *drdr, doublereal *
	dhdl, doublereal *dhdx, doublereal *dhdr, doublereal *dhdd, 
	doublereal *dfds, doublereal *dxds, doublereal *drds, doublereal *
	dhds, doublereal *dedt, doublereal *p, doublereal *eyv, doublereal *
	dlambda, doublereal *tdl, doublereal *unmd, doublereal *unmd1, 
	doublereal *dmwx, doublereal *wx, doublereal *phtrwx, doublereal *
	phtrwx2, doublereal *h__, doublereal *g, doublereal *xk, doublereal *
	rdk, doublereal *rsk, doublereal *rsi, doublereal *xms, doublereal *
	teta, integer *ncd, integer *ntens, integer *ndi, integer *nshr, 
	integer *m, integer *mpi)
{
    /* Initialized data */

    static doublereal zero = 0.;
    static doublereal one = 1.;
    static doublereal two = 2.;
    static doublereal three = 3.;
    static doublereal xnine = 9.;
    static doublereal half = .5;
    static doublereal thalf = 1.5;

    /* Format strings */
    static char fmt_102[] = "(1x,a20,6(d20.14,2x))";
    static char fmt_101[] = "(1x,6(d20.13,2x))";

    /* System generated locals */
    integer dnde_dim1, dnde_offset, alde_dim1, alde_offset, dpnde_dim1, 
	    dpnde_offset, dzde_dim1, dzde_offset, dxide_dim1, dxide_offset, 
	    uzdnn_dim1, uzdnn_offset, unn_dim1, unn_offset, tntn_dim1, 
	    tntn_offset, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    double pow_dd(doublereal *, doublereal *), cosh(doublereal);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    static doublereal onethird, twothird, a[6], b[6];
    static integer i__, j;
    static doublereal tdlphdmwx;
    extern /* Subroutine */ int transpose_(integer *, doublereal *, 
	    doublereal *);
    static doublereal ad, bd, al, bl, pi, xi[6], d2g, adl, amm, dhs, dls, dxk,
	     d2gd, d2gt, xmm1, xms2, albd, eeal[6], d2gd3, dlca, drde[6], 
	    dxde[6], dsde[6], seff, aphi, cphi, dfrr, dfxx, dhxx, dh1dd, 
	    dh2dd, dh3dd, dhrr;
    extern /* Subroutine */ int aset_(doublereal *, doublereal *, integer *);
    static doublereal seff2, dpfde[6], dphde[6], fgxdd[6], cosha, dlcay, 
	    dqtdg, dpxde[6], unmdm, xkdgt, udrdr, xieqv, udxdx, uzdzd, trdwx, 
	    h1seff, trswx, xieqv2, adldkp, eealdv[6], dldxde[6], dlgdhx, 
	    dpndet[36];
    extern /* Subroutine */ int pdtsca_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal dxidet[36];
    extern /* Subroutine */ int prdmat_(integer *, doublereal *, doublereal *,
	     doublereal *), pdtmat_(integer *, doublereal *, doublereal *, 
	    doublereal *), pdtten_(doublereal *, doublereal *, doublereal *, 
	    integer *, integer *);
    static doublereal unmdxk, h2pksef, dseffde[6], dlcaxdk, xidxide[6];
    extern /* Subroutine */ int tensdev_(doublereal *, doublereal *, integer *
	    , integer *);
    static doublereal trddlrs;

    /* Fortran I/O blocks */
    static cilist io___121 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___122 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___168 = { 0, 0, 0, 0, 0 };
    static cilist io___169 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___170 = { 0, 0, 0, 0, 0 };
    static cilist io___171 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___172 = { 0, 0, 0, 0, 0 };
    static cilist io___173 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___174 = { 0, 0, 0, 0, 0 };
    static cilist io___175 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___176 = { 0, 0, 0, 0, 0 };
    static cilist io___177 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___178 = { 0, 0, 0, 0, 0 };
    static cilist io___179 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___180 = { 0, 0, 0, 0, 0 };
    static cilist io___181 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___182 = { 0, 0, 0, 0, 0 };
    static cilist io___183 = { 0, 0, 0, fmt_101, 0 };


/* ====================================================================== */
/* ---------------------------------------------------------------------- */
/* -----------------------  CALCUL DES DERIVEES  ------------------------ */
/* -------------  PAR RAPPORT AU TENSEUR DES DEFORMATIONS  -------------- */
/* ---------------------------------------------------------------------- */
/* ====================================================================== */
/* ------------- */
/*  EN ENTREE : */
/* ------------- */
/*  TETN		: TENSEUR = ESTAR - TDLRD*TN */
/*  ALTN		: TENSEUR = ALPHA1 + TDLRD*TN */
/*  PLTN		: PRODUCT P:|dEp/dl| */
/*  PLTN		: PRODUCT Q:|dEp/dl| */
/*  ESTAR	: TENSEUR e* */
/*  DEDT		: TENSEUR INCREMENT DE DEFORMATION : DSTRAN */
/*  ALPHA1	: TENSEUR D'ECROUISSAGE CINEMATIQUE CONNU */
/*  TN		: TENSEUR NORMAL */
/*  ZDZD		: NORME DU TENSEUR ZD = RACINE (2/3*ZD:ZD) */
/*  DFDL		: DERIVEE PARTIELLE DE F PAR RAPPORT A DLAMBDA */
/*  DFDX		: DERIVEE PARTIELLE DE F PAR RAPPORT A WX */
/*  DFDR		: DERIVEE PARTIELLE DE F PAR RAPPORT A WR */
/*  DFDD		: DERIVEE PARTIELLE DE F PAR RAPPORT A ENDO */
/*  DXDL		: DERIVEE PARTIELLE DE GX PAR RAPPORT A DLAMBDA */
/*  DXDX		: DERIVEE PARTIELLE DE GX PAR RAPPORT A WX */
/*  DXDD		: DERIVEE PARTIELLE DE GX PAR RAPPORT A ENDO */
/*  DRDL		: DERIVEE PARTIELLE DE GR PAR RAPPORT A DLAMBDA */
/*  DRDR		: DERIVEE PARTIELLE DE GR PAR RAPPORT A WR */
/*  DHDL		: DERIVEE PARTIELLE DE H PAR RAPPORT A DLAMBDA */
/*  DHDX		: DERIVEE PARTIELLE DE H PAR RAPPORT A WX */
/*  DHDR		: DERIVEE PARTIELLE DE H PAR RAPPORT A WR */
/*  DHDD		: DERIVEE PARTIELLE DE H PAR RAPPORT A ENDO */
/*  DFDS		: DERIVEE PARTIELLE DE F PAR RAPPORT A DLs */
/*  DXDS		: DERIVEE PARTIELLE DE FGX PAR RAPPORT A DLs */
/*  DRDS		: DERIVEE PARTIELLE DE FGR PAR RAPPORT A DLs */
/*  DHDS		: DERIVEE PARTIELLE DE H PAR RAPPORT A DLs */
/*  DMWX		: NORME DU TENSEUR Mt:(Alpha1+Tdl*N) */
/*  DLAMBDA	: INCREMENT DEFORMATION PLASTIQUE CUMULEE */
/*  YPS1		: VALEUR DE LA FONCTION YPS = (Y/S)**(s-1) */
/*  TDL		: TETA*DLAMBDA */
/*  DLRD		: TETA*DLAMBDA/UNMD */
/*  WR		: ISOTROPIC DYNAMIC RECOVERY TERM: WR = DSINH[QS*R] */
/*  WX		: KINEMATIC DYNAMIC RECOVERY TERM: WX = |DMI:X| */
/*  PHTRWX	: (2/3)*PH/(1+TETA*(PRD*DLAMBDA+PRS)*WX) */
/*  PHTRWX2	: (2/3)*PH/(1+TETA*(PRD*DLAMBDA+PRS)*WX)**2 */
/*  G		: MODULE DE CISAILLEMENT */
/*  XK		: COEFFICIENT DE LAME */
/*  RDK		: CONSTANTE D'ECROUISSAGE CINEMATIQUE */
/*  RSK		: CONSTANTE D'ECROUISSAGE CINEMATIQUE (RECOVERY) */
/*  RSI		: CONSTANTE D'ECROUISSAGE ISOTROPE (RECOVERY) */
/*  XMS		: DAMAGE RATE SENSITIVITY */
/*  TETA		: TETA-METHODE */
/*  NCD		: NCD=1 --> COUPLAGE ELASTOPLASTICITE-ENDOMMAGEMENT */
/*  		  NCD=0 --> ELASTOPLASTICITE SANS ENDOMMAGEMENT */
/*  NDI		: NOMBRE DE COMPOSANTES DES CONTRAINTES DIRECTES */
/*  NSHR		: NOMBRE DE COMPOSANTES DES CONTRAINTES TANGENTIELLES */
/*  NTENS	: NOMBRE TOTAL DE COMPOSANTES D'UN TENSEUR */
/*  M		: COMMENTAIRES */
/*  MPI		: UNITE DE FICHIER DE SORTIE */
/* ------------- */
/*  EN LOCAL  : */
/* ------------- */
/*  TNTN		: TENSEUR (D'ORDRE 4) = N*N */
/*  ALDE		: TENSEUR (D'ORDRE 4) = ALPHA1*DLDE */
/*  UDEVNN	: TENSEUR (D'ORDRE 4) = IDEV - 2/3 N*N */
/*  DPNDE	: TENSEUR (D'ORDRE 4) = 2G/ZDZD.(IDEV - 2/3 N*N) */
/*  DPXDE	: DERIVEE PARTIELLE DE LA FONCTION FGX PAR RAPPORT A DE */
/*  DPHDE	: DERIVEE PARTIELLE DE LA FONCTION H PAR RAPPORT A DE */
/* ------------- */
/*  EN SORTIE : */
/* ------------- */
/*  DLDE		: TENSEUR (D'ORDRE 2) DERIVEE DE DLAMBDA */
/*  DDDE		: TENSEUR (D'ORDRE 2) DERIVEE DE D */
/*  DNDE		: TENSEUR (D'ORDRE 4) DERIVEE DU TENSEUR NORMAL N */
/*  UDEVNN	: TENSEUR (D'ORDRE 4) = IDEV - 2/3 N*N */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
/* .1-----  Implicit, External */
/* .2-----  Parameter,Integer */
/* .3-----  Dimension */
/* ------- */
/* -----------------------------------------------------------* */
/* .4-----  Real,Complex,Double precision,Logical,Character */
/* .5-----  Common */
/* .6-----  Equivalence */
/* .7-----  Save,Data */
    /* Parameter adjustments */
    --h__;
    --dedt;
    --altn;
    --tetn;
    --tn;
    --zd;
    --alpha1;
    --estar;
    dxide_dim1 = *ntens;
    dxide_offset = 1 + dxide_dim1;
    dxide -= dxide_offset;
    dzde_dim1 = *ntens;
    dzde_offset = 1 + dzde_dim1;
    dzde -= dzde_offset;
    dpnde_dim1 = *ntens;
    dpnde_offset = 1 + dpnde_dim1;
    dpnde -= dpnde_offset;
    alde_dim1 = *ntens;
    alde_offset = 1 + alde_dim1;
    alde -= alde_offset;
    uzdnn_dim1 = *ntens;
    uzdnn_offset = 1 + uzdnn_dim1;
    uzdnn -= uzdnn_offset;
    tntn_dim1 = *ntens;
    tntn_offset = 1 + tntn_dim1;
    tntn -= tntn_offset;
    unn_dim1 = *ntens;
    unn_offset = 1 + unn_dim1;
    unn -= unn_offset;
    dnde_dim1 = *ntens;
    dnde_offset = 1 + dnde_dim1;
    dnde -= dnde_offset;
    --ddde;
    --dlde;

    /* Function Body */
/* .8-----  Definition de fonctions */
/* ---------------------------------------------------------------------* */
/* ====================================================================== */
/* ---------------------------------------------------------------------* */
/* -------  Initialization */

    onethird = one / three;
    twothird = two / three;
    d2g = two * *g;
    d2gd = d2g * *unmd;
    dqtdg = d2gd * twothird;
    uzdzd = one / *zdzd;
    d2gt = d2g / three;
    xkdgt = *xk - d2gt;

/* -------  END Initialization */
/* ---------------------------------------------------------------------* */
/* ---------------------------------------------------------------------* */
/* -------  Calcul des tenseurs UZDNN et DPNDE */

    pdtten_(&tn[1], &tn[1], &tntn[tntn_offset], ntens, ndi);
    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *ntens;
	for (j = 1; j <= i__2; ++j) {
	    unn[i__ + j * unn_dim1] = -twothird * tntn[i__ + j * tntn_dim1];
	}
	unn[i__ + i__ * unn_dim1] += one;
    }
    i__1 = *ndi;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *ndi;
	for (j = 1; j <= i__2; ++j) {
	    unn[i__ + j * unn_dim1] -= onethird;
	}
    }
/* ------- */
    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *ntens;
	for (j = 1; j <= i__2; ++j) {
	    uzdnn[i__ + j * uzdnn_dim1] = uzdzd * unn[i__ + j * unn_dim1];
	    dpnde[i__ + j * dpnde_dim1] = d2gd * uzdnn[i__ + j * uzdnn_dim1];
	}
    }

/* -------  END Calcul des tenseurs UZDNN et DPNDE */
/* ---------------------------------------------------------------------* */
/* ---------------------------------------------------------------------* */
/* -------  Tensor DPFDE */

    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dpfde[i__ - 1] = dqtdg * tn[i__];
    }

/* -------  END Tensor DPFDE */
/* ---------------------------------------------------------------------* */
/* ---------------------------------------------------------------------* */
/* -------  Tensor DPXDE */

    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
	altn[i__] = alpha1[i__] + *tdl * tn[i__];
    }
    transpose_(ntens, &dpnde[dpnde_offset], dpndet);
    pdtmat_(ntens, dpndet, &altn[1], fgxdd);
/* ------- */
    tdlphdmwx = *tdl * *phtrwx / *dmwx;
    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dpxde[i__ - 1] = -tdlphdmwx * fgxdd[i__ - 1];
    }

/* -------  END  Tensor DPXDE */
/* ---------------------------------------------------------------------* */
/* ---------------------------------------------------------------------* */
/* -------  Definition of the Tensor DSDE (Static Recovery) */

    if (*rsk == zero && *rsi == zero) {
	aset_(dsde, &zero, ntens);
    } else {
	pdtsca_(ndi, nshr, &tn[1], &dedt[1], &dls);
	if (dls == zero) {
	    dls = one;
	}
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dsde[i__ - 1] = tn[i__] / dls;
	}
    }

/* -------  END  Definition of the Tensor DSDE (Static Recovery) */
/* ---------------------------------------------------------------------* */
/* ---------------------------------------------------------------------* */
/* -------  Partial Derivative DPHDE */

    if (*ncd == 0) {
	aset_(dphde, &zero, ntens);
	goto L100;
    }

/* ------------------------------------------------------* */
/* ------------------------------------------------------* */
/* -------  Partial Derivative DPHDE */

    d2gd3 = onethird * d2gd;
    dlgdhx = *tdl * (d2gd + *phtrwx);
    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xi[i__ - 1] = zd[i__] - dlgdhx * tn[i__];
    }
/* -------------------------------------------* */
/* ---------  Norm |Xi| = Sqrt(3/2 * Xi:Xi) */
    xieqv = xi[0] * xi[0] + xi[1] * xi[1] + xi[2] * xi[2];
    i__1 = *ntens;
    for (i__ = *ndi + 1; i__ <= i__1; ++i__) {
	xieqv += two * xi[i__ - 1] * xi[i__ - 1];
    }
    xieqv = sqrt(thalf * xieqv);
    seff = sqrt(h__[1] * xieqv * xieqv + h__[2] * *p * *p);
/* -------------------------------------------* */
/* ---------  Partial Derivative DXIDE = dXi/dE */
    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *ntens;
	for (j = 1; j <= i__2; ++j) {
	    dxide[i__ + j * dxide_dim1] = -dlgdhx * dpnde[i__ + j * 
		    dpnde_dim1];
	}
	dxide[i__ + i__ * dxide_dim1] += d2gd;
    }
    i__1 = *ndi;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *ndi;
	for (j = 1; j <= i__2; ++j) {
	    dxide[i__ + j * dxide_dim1] -= d2gd3;
	}
    }
/* -------------------------------------------* */
/* ---------  Product Xi : dXi/dE */
    transpose_(ntens, &dxide[dxide_offset], dxidet);
    pdtmat_(ntens, dxidet, xi, xidxide);
/* -------------------------------------------* */
/* -------  Comment Check Code */
    if (*m >= 4) {
	io___121.ciunit = *mpi;
	s_wsfe(&io___121);
	do_fio(&c__1, "STRESS XI", (ftnlen)9);
	e_wsfe();
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___122.ciunit = *mpi;
	    s_wsfe(&io___122);
	    do_fio(&c__1, (char *)&xi[i__ - 1], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }
/* -----------------------------------------------------------* */
/* -----------------------------------------------------------* */
/* --------- Choice of Damage Evolution */
    if (*ncd == 1) {
/* -----------------------------------------------------------* */
/* ---------  Cocks [1989] */
/* -----------------------------------------------------------* */
	*eyv = onethird * *p / seff;
	xmm1 = *xms / (*xms + one);
	ad = two * (*xms + one) * (one + *phi) * pow_dd(unmd, &xmm1);
	aphi = xnine * *xms * *phi * *unmd * *eyv / ad;
/* -----------------------------------------------------------* */
	dh1dd = twothird;
	dh2dd = half * xmm1 / (*unmd1 * *unmd1);
	d__1 = xmm1 - one;
	dh3dd = xmm1 * pow_dd(unmd1, &d__1);
	dh1dd = zero;
	dh2dd = zero;
	dh3dd = zero;
/* -----------------------------------------------------------* */
/* ---------  Partial Derivative of DSEFFDE = d[1/Seff]/dE */
	dxk = onethird * *unmd * *xk;
	xieqv2 = xieqv * xieqv;
	seff2 = seff * seff;
	h1seff = thalf * h__[1] / seff2;
	h2pksef = dxk * h__[2] * *p / seff;
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dseffde[i__ - 1] = h1seff * xidxide[i__ - 1];
	}
	i__1 = *ndi;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dseffde[i__ - 1] += h2pksef;
	}
/* -------------------------------------------* */
/* ---------  Partial Derivative DPHDE = dFphi/dE */
	adl = aphi * *dlambda;
	pi = *p;
	if (*p == zero) {
	    pi = one;
	}
	adldkp = adl * dxk / pi;
/* ----------* */
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dphde[i__ - 1] = adl * dseffde[i__ - 1];
	}
	i__1 = *ndi;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dphde[i__ - 1] -= adldkp;
	}
/* -----------------------------------------------------------* */
    } else {
/* -----------------------------------------------------------* */
/* ---------  Cocks and Ashby [1980] */
/* -----------------------------------------------------------* */
	*eyv = onethird * *p / xieqv;
	xms2 = two * *xms;
	unmdm = pow_dd(unmd, xms);
	amm = two * (xms2 - one) / (xms2 + one);
	aphi = amm * *eyv;
	cphi = one / unmdm - *unmd;
/* -------------------------------------------* */
	xieqv2 = xieqv * xieqv;
	unmdxk = onethird * *unmd * *xk;
	cosha = cosh(aphi);
	dlca = *dlambda * cphi * cosha * amm;
	dlcay = thalf * dlca * *eyv / xieqv2;
	dlcaxdk = unmdxk * dlca / xieqv;
/* -------------------------------------------* */
/* ---------  Partial Derivative DHDE = dFphi/dE */
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dphde[i__ - 1] = dlcay * xidxide[i__ - 1];
	}
	i__1 = *ndi;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dphde[i__ - 1] -= dlcaxdk;
	}
/* -------------------------------------------* */
/* -----------------------------------------------------------* */
    }

/* -------  END Partial Derivative DPHDE = dH/dE */
/* ---------------------------------------------------------------------* */
/* ---------------------------------------------------------------------* */
/* -------  Tensor DDDE */

L100:
    udxdx = one / *dxdx;
    udxdx = one / *dxdx;
    udrdr = one / *drdr;
    dfxx = *dfdx * udxdx;
    dhxx = *dhdx * udxdx;
    dfrr = *dfdr * udrdr;
    dhrr = *dhdr * udrdr;
/* ------- */
    dls = *dfds - dfxx * *dxds - dfrr * *drds;
    dhs = *dhds - dhxx * *dxds - dhrr * *drds;
    al = *dfdl - dfxx * *dxdl - dfrr * *drdl;
    ad = *dfdd - dfxx * *dxdd;
    bl = *dhdl - dhxx * *dxdl - dhrr * *drdl;
    bd = *dhdd - dhxx * *dxdd;
/* ------------------------------------------* */
    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
	a[i__ - 1] = dpfde[i__ - 1] - dfxx * dpxde[i__ - 1] + dls * dsde[i__ 
		- 1];
	b[i__ - 1] = dphde[i__ - 1] - dhxx * dpxde[i__ - 1] + dhs * dsde[i__ 
		- 1];
    }
/* ------- */
    if (*ncd == 0) {
	ad = zero;
	aset_(&ddde[1], &zero, ntens);
    } else {
	albd = one / (al * bd - ad * bl);
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ddde[i__] = albd * (bl * a[i__ - 1] - al * b[i__ - 1]);
	}
    }

/* -------  END  Tensor DDDE */
/* ---------------------------------------------------------------------* */
/* ---------------------------------------------------------------------* */
/* -------  Tensor DLDE */

    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dlde[i__] = -(a[i__ - 1] + ad * ddde[i__]) / al;
    }

/* -------  END  Tensor DLDE */
/* ---------------------------------------------------------------------* */
/* ---------------------------------------------------------------------* */
/* -------  Tensors DXDE and DRDE */

    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dxde[i__ - 1] = *dxdl * dlde[i__] + *dxdd * ddde[i__] + *dxds * dsde[
		i__ - 1] + dpxde[i__ - 1];
	dxde[i__ - 1] = -udxdx * dxde[i__ - 1];
	drde[i__ - 1] = -udrdr * (*drdl * dlde[i__] + *drds * dsde[i__ - 1]);
    }

/* -------  END  Tensors DXDE and DRDE */
/* ---------------------------------------------------------------------* */
/* ---------------------------------------------------------------------* */
/* -------  Tensor DNDE */

    trdwx = *teta * *rdk * *wx;
    trswx = *teta * *rsk * *wx;
    trddlrs = *teta * *rdk * *dlambda + *rsk;
    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dldxde[i__ - 1] = trdwx * dlde[i__] + trddlrs * dxde[i__ - 1] + trswx 
		* dsde[i__ - 1];
    }
    pdtten_(&alpha1[1], dldxde, &alde[alde_offset], ntens, ndi);
/* ------- */
    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
	eeal[i__ - 1] = -d2g * estar[i__];
    }
    tensdev_(eeal, eealdv, ndi, nshr);
    pdtten_(eealdv, &ddde[1], &dzde[dzde_offset], ntens, ndi);
/* ------- */
    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *ntens;
	for (j = 1; j <= i__2; ++j) {
	    dzde[i__ + j * dzde_dim1] += *phtrwx2 * alde[i__ + j * alde_dim1];
	}
	dzde[i__ + i__ * dzde_dim1] += d2gd;
    }
    prdmat_(ntens, &uzdnn[uzdnn_offset], &dzde[dzde_offset], &dnde[
	    dnde_offset]);

/* -------  END  Tensor DNDE */
/* ---------------------------------------------------------------------* */
/* -------  Comment Check Code */

    if (*m >= 4) {
	io___168.ciunit = *mpi;
	s_wsle(&io___168);
	do_lio(&c__9, &c__1, "CALDDE DLDE", (ftnlen)11);
	e_wsle();
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___169.ciunit = *mpi;
	    s_wsfe(&io___169);
	    do_fio(&c__1, (char *)&dlde[i__], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	io___170.ciunit = *mpi;
	s_wsle(&io___170);
	do_lio(&c__9, &c__1, "CALDDE DXDE", (ftnlen)11);
	e_wsle();
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___171.ciunit = *mpi;
	    s_wsfe(&io___171);
	    do_fio(&c__1, (char *)&dxde[i__ - 1], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	io___172.ciunit = *mpi;
	s_wsle(&io___172);
	do_lio(&c__9, &c__1, "CALDDE DRDE", (ftnlen)11);
	e_wsle();
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___173.ciunit = *mpi;
	    s_wsfe(&io___173);
	    do_fio(&c__1, (char *)&drde[i__ - 1], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	io___174.ciunit = *mpi;
	s_wsle(&io___174);
	do_lio(&c__9, &c__1, "CALDDE DDDE", (ftnlen)11);
	e_wsle();
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___175.ciunit = *mpi;
	    s_wsfe(&io___175);
	    do_fio(&c__1, (char *)&ddde[i__], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	io___176.ciunit = *mpi;
	s_wsle(&io___176);
	do_lio(&c__9, &c__1, "CALDDE UZDNN", (ftnlen)12);
	e_wsle();
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___177.ciunit = *mpi;
	    s_wsfe(&io___177);
	    i__2 = *ntens;
	    for (j = 1; j <= i__2; ++j) {
		do_fio(&c__1, (char *)&uzdnn[i__ + j * uzdnn_dim1], (ftnlen)
			sizeof(doublereal));
	    }
	    e_wsfe();
	}
	io___178.ciunit = *mpi;
	s_wsle(&io___178);
	do_lio(&c__9, &c__1, "CALDDE DPNDE", (ftnlen)12);
	e_wsle();
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___179.ciunit = *mpi;
	    s_wsfe(&io___179);
	    i__2 = *ntens;
	    for (j = 1; j <= i__2; ++j) {
		do_fio(&c__1, (char *)&dpnde[i__ + j * dpnde_dim1], (ftnlen)
			sizeof(doublereal));
	    }
	    e_wsfe();
	}
	io___180.ciunit = *mpi;
	s_wsle(&io___180);
	do_lio(&c__9, &c__1, "CALDDE ALDE", (ftnlen)11);
	e_wsle();
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___181.ciunit = *mpi;
	    s_wsfe(&io___181);
	    i__2 = *ntens;
	    for (j = 1; j <= i__2; ++j) {
		do_fio(&c__1, (char *)&alde[i__ + j * alde_dim1], (ftnlen)
			sizeof(doublereal));
	    }
	    e_wsfe();
	}
	io___182.ciunit = *mpi;
	s_wsle(&io___182);
	do_lio(&c__9, &c__1, "CALDDE DNDE", (ftnlen)11);
	e_wsle();
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___183.ciunit = *mpi;
	    s_wsfe(&io___183);
	    i__2 = *ntens;
	    for (j = 1; j <= i__2; ++j) {
		do_fio(&c__1, (char *)&dnde[i__ + j * dnde_dim1], (ftnlen)
			sizeof(doublereal));
	    }
	    e_wsfe();
	}
    }
/* -------  END Comment Check Code */
/* ----------------------------------------------------------------------* */
/* ---------   FORMATS  -------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* ======================================================================= */
    return 0;
} /* caldde_ */


/* Subroutine */ int caldpf_(doublereal *dfdl, doublereal *dfdx, doublereal *
	dfdr, doublereal *dfdd, doublereal *dfds, doublereal *drds, 
	doublereal *dpdd, doublereal *dzdndl, doublereal *dzdndd, doublereal *
	alphzd, doublereal *dfrdl, doublereal *zd, doublereal *h__, 
	doublereal *p, doublereal *phi, doublereal *alpha1, doublereal *
	stres1d, doublereal *pr1, doublereal *zdn, doublereal *zdzd, 
	doublereal *dsqhmp, doublereal *phtrwx, doublereal *ryvshdl, 
	doublereal *dlambda, doublereal *tdl, doublereal *unmd, doublereal *
	unmd1, doublereal *xkdh, doublereal *wx, doublereal *wr, doublereal *
	phtrwx2, doublereal *htrwr, doublereal *htrwr2, doublereal *trswr, 
	doublereal *rtdl, doublereal *rdk, doublereal *rsk, doublereal *rdi, 
	doublereal *rsi, doublereal *yield, doublereal *g, doublereal *xk, 
	doublereal *v, doublereal *pf, doublereal *xms, doublereal *xmm1, 
	doublereal *teta, integer *ncd, integer *melpst, integer *ntens, 
	integer *ndi, integer *nshr, integer *m, integer *mpi)
{
    /* Initialized data */

    static doublereal zero = 0.;
    static doublereal one = 1.;
    static doublereal two = 2.;
    static doublereal three = 3.;
    static doublereal half = .5;

    /* Format strings */
    static char fmt_102[] = "(1x,a20,6(d20.14,2x))";

    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal), pow_dd(doublereal *, doublereal *);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static doublereal onethird, ondsqhmp, twothird, d2g, dtalphzdzd, h1h1, 
	    h2h1, h3h1, dh1dd, dh2dd, dh3dd, upsd1, trswx, dh2h1dd, dh3h1dd;
    extern /* Subroutine */ int pdtsca_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal dryvsh, streszd;

    /* Fortran I/O blocks */
    static cilist io___206 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___207 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___208 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___209 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___210 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___211 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___212 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___213 = { 0, 0, 0, fmt_102, 0 };


/* ====================================================================== */
/* ---------------------------------------------------------------------- */
/* ---------------  CALCUL DES DERIVEES PARTIELLES DE F  ---------------- */
/* --------------------  PAR RAPPORT DLAMBDA ET ENDO  ------------------- */
/* ---------------------------------------------------------------------- */
/* ====================================================================== */
/* ------------- */
/*  EN ENTREE : */
/* ------------- */
/*  ZD		: TENSEUR DE DIRECTION NORMALE */
/*  ENDO		: TENSEUR D'ENDOMMAGEMENT */
/*  ESTAR	: TENSEUR DES DEFORMATIONS CONNUS A T+DT */
/*  ALPHA1	: TENSEUR D'ECROUISSAGE CINEMATIQUE CONNU */
/*  PR1		: VARIABLE D'ECROUISSAGE ISOTROPE CONNUE */
/*  ZDN		: NORME DU TENSEUR DE DIRECTION NORMALE ZD */
/*  ZDZD		: NORME DU ZD = RACINE (2/3*ZD:ZD) */
/*  TDL		: TETA*DLAMBDA */
/*  WX		: KINEMATIC DYNAMIC RECOVERY TERM: WX = |DMI:X| */
/*  WR		: ISOTROPIC DYNAMIC RECOVERY TERM: WR = R */
/*  PHTRWX	: (2/3)*PH/(1+TETA*(RDK*DLAMBDA+RSK)*WX) */
/*  PHTRWX2	: (2/3)*PH/(1+TETA*(RDK*DLAMBDA+RSK)*WX)*RSK)*WX)**2 */
/*  HTRWR	: RT(2/3)*H/(1+TETA*(RDK*DLAMBDA+RS)) */
/*  HTRWR2	: RT(2/3)*H/(1+TETA*(RDK*DLAMBDA+RS))**2 */
/*  TRSWR	: TETA*RS*WR */
/*  RTDL		: PR1+TETA*DLAMBDA */
/*  RD		: CONSTANTE D'ECROUISSAGE ISOTROPE */
/*  RS		: CONSTANTE D'ECROUISSAGE ISOTROPE (RECOVERY) */
/*  RDK		: CONSTANTE D'ECROUISSAGE CINEMATIQUE */
/*  RSK		: CONSTANTE D'ECROUISSAGE CINEMATIQUE (RECOVERY) */
/*  YIELD	: LIMITE D'ELASTICITE */
/*  G		: MODULE DE CISAILLEMENT */
/*  V		: CONSTANTE VISCOPLASTIQUE */
/*  PF		: CONSTANTE VISCOPLASTIQUE */
/*  TETA		: TETA-METHODE */
/*  NCD		: NCD=1 --> COUPLAGE ELASTOPLASTICITE-ENDOMMAGEMENT */
/*  		  NCD=0 --> ELASTOPLASTICITE SANS ENDOMMAGEMENT */
/*  MELPST	: MELPST = 0 --> MATRICE CONSISTANTE AVEC V = 1 */
/*  NDI		: NOMBRE DE COMPOSANTES DES CONTRAINTES DIRECTES */
/*  NSHR		: NOMBRE DE COMPOSANTES DES CONTRAINTES TANGENTIELLES */
/*  NTENS	: TAILLE DU TABLEAU DES CONTRAINTES OU DEFORMATIONS */
/* 		  (NDI + NSHR). */
/*  M		: COMMENTAIRES */
/*  MPI		: UNITE DE FICHIER DE SORTIE */
/* ------------- */
/*  EN SORTIE : */
/* ------------- */
/*  DFDL		: DERIVEE PARTIELLE DE F PAR RAPPORT A DLAMBDA */
/*  DFDX		: DERIVEE PARTIELLE DE F PAR RAPPORT A WX */
/*  DFDR		: DERIVEE PARTIELLE DE F PAR RAPPORT A WR */
/*  DFDD		: DERIVEE PARTIELLE DE F PAR RAPPORT A ENDO */
/*  DZDNDL	: DERIVEE PARTIELLE DE ZDN PAR RAPPORT A DLAMBDA */
/*  DFRDL	: DERIVEE PARTIELLE DE R PAR RAPPORT A DLAMBDA */
/*  DMDDESTAR	: PRODUIT CONTRACTE DMDD : ESTAR */
/*  ALPHZD	: PRODUIT CONTRACTE ALPHA1 : ZD */
/*  DMDDALPH	: PRODUIT CONTRACTE DMDD : ALPHA1 */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
/* ----------------------------------------------------------------------- */
/* .1-----  Implicit, External */
/* .2-----  Parameter,Integer */
/* .3-----  Dimension */
/* ------- */
/* .4-----  Real,Integer,Complex,Double precision,Logical,Character */
/* .5-----  Common */
/* .6-----  Equivalence */
/* .7-----  Save,Data */
    /* Parameter adjustments */
    --h__;
    --stres1d;
    --alpha1;
    --zd;

    /* Function Body */
/* .8-----  Definition de fonctions */
/* ---------------------------------------------------------------------* */
/* ====================================================================== */
/* ---------------------------------------------------------------------* */
/* -------  Initialisation */
    onethird = one / three;
    twothird = two * onethird;
    d2g = two * *g;
    h2h1 = h__[2] / h__[1];
    h3h1 = h__[3] * h__[3] / h__[1];
/* -------  END Initialisation */
/* ---------------------------------------------------------------------* */
/* -------  Partial Derivative DZNDDL = dZDN/dlambda */
    if (*melpst == 0) {
	*v = zero;
    }
    ondsqhmp = one / *dsqhmp;
    *dfrdl = *teta * *htrwr2 * (one + *trswr - *pr1 * *rdi * *wr);
    dryvsh = *dfrdl + *v / sqrt(*pf * *pf + *dlambda * *dlambda);
    *dzdndl = twothird * *unmd * h3h1 * dryvsh * *ryvshdl * ondsqhmp + *teta *
	     (d2g * *unmd + *phtrwx - *tdl * *rdk * *wx * *phtrwx2);
/* -------  END Partial Derivative DFDL = df/dlambda */
/* ---------------------------------------------------------------------* */
/* -------  Partial Derivative DFDL = df/dlambda */
    pdtsca_(ndi, nshr, &alpha1[1], &zd[1], alphzd);
    dtalphzdzd = twothird * *alphzd / *zdzd;
    *dfdl = *teta * *rdk * *wx * *phtrwx2 * dtalphzdzd - *dzdndl;
/* -------  END Partial Derivative DFDL = df/dlambda */
/* ---------------------------------------------------------------------* */
/* -------  Partial Derivative DFDX = df/dWx */
    *dfdx = *teta * (*dlambda * *rdk + *rsk) * *phtrwx2 * (*tdl + dtalphzdzd);
/* -------  END Partial Derivative DFDX = df/dWx */
/* ---------------------------------------------------------------------* */
/* -------  Partial Derivative DFDR = df/dWr */
    *dfdr = *teta * (*dlambda * *rdi + *rsi) * *htrwr2 * *rtdl;
    *dfdr = twothird * *dfdr * *unmd * *ryvshdl * ondsqhmp;
/* -------  END Partial Derivative DFDR = df/dWr */
/* ---------------------------------------------------------------------* */
/* ---------  Partial Derivative DHxDD = dhx/dphi */
    if (*ncd == 1) {
	upsd1 = one / (one + *phi);
	dh1dd = twothird;
	dh2dd = half * *xmm1 * upsd1 * upsd1;
	d__1 = *xmm1 - one;
	dh3dd = *xmm1 * pow_dd(unmd1, &d__1);
    } else if (*ncd == 2) {
	dh1dd = zero;
	dh2dd = half / (*unmd1 * *unmd1);
	dh3dd = -one;
    }
    dh1dd = zero;
    dh2dd = zero;
    dh3dd = zero;
/* --------- */
    h1h1 = one / (h__[1] * h__[1]);
    h2h1 = h__[2] / h__[1];
    dh2h1dd = (h__[1] * dh2dd - h__[2] * dh1dd) / h1h1;
    dh3h1dd = (two * h__[1] * h__[3] * dh3dd - h__[3] * h__[3] * dh1dd) / 
	    h1h1;
/* ---------------------------------------------------------------------* */
/* -------  Partial Derivative DFDD = df/dD */
/*               Remark:  dZ/dD(t):Z = Z:dZ/dD */
    if (*ncd == 0) {
	*dfdd = zero;
    } else {
	pdtsca_(ndi, nshr, &stres1d[1], &zd[1], &streszd);
	*dpdd = -(*xkdh);
	*dzdndd = dh3h1dd * *ryvshdl * *ryvshdl - *p * *p * dh2h1dd - two * 
		h2h1 * *p * *dpdd;
	*dzdndd = half * *unmd * ondsqhmp * *dzdndd - *dsqhmp;
	*dzdndd = twothird * *dzdndd - d2g * *dlambda;
	*dfdd = -twothird * streszd / *zdzd - *dzdndd;
    }
/* -------  END Partial Derivative DFDD = df/dD */
/* ---------------------------------------------------------------------* */
/* -------  Partial Derivative DRDS = dR/dLs */
    *trswr = *teta * *rsi * *wr;
    *drds = *unmd * *trswr * *htrwr2 * *rtdl;
/* -------  END Partial Derivative DRDS = dR/dLs */
/* ---------------------------------------------------------------------* */
/* -------  Partial Derivative DFDS = df/dLs */
    trswx = *teta * *rsk * *wx;
    *dfds = twothird * *drds + trswx * *phtrwx2 * (*tdl + twothird * *alphzd /
	     *zdzd);
/* -------  END Partial Derivative DFDS = df/dLs */
/* ---------------------------------------------------------------------* */
/* ---------------------------------------------------------------------* */
/* -------  Comment Check Code */

    if (*m >= 4) {
	io___206.ciunit = *mpi;
	s_wsfe(&io___206);
	do_fio(&c__1, "CALDPF  DZDNDL = ", (ftnlen)17);
	do_fio(&c__1, (char *)&(*dzdndl), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___207.ciunit = *mpi;
	s_wsfe(&io___207);
	do_fio(&c__1, "CALDPF  DZDNDD = ", (ftnlen)17);
	do_fio(&c__1, (char *)&(*dzdndd), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___208.ciunit = *mpi;
	s_wsfe(&io___208);
	do_fio(&c__1, "CALDPF  DFDL   = ", (ftnlen)17);
	do_fio(&c__1, (char *)&(*dfdl), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___209.ciunit = *mpi;
	s_wsfe(&io___209);
	do_fio(&c__1, "CALDPF  DFDX   = ", (ftnlen)17);
	do_fio(&c__1, (char *)&(*dfdx), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___210.ciunit = *mpi;
	s_wsfe(&io___210);
	do_fio(&c__1, "CALDPF  DFDR   = ", (ftnlen)17);
	do_fio(&c__1, (char *)&(*dfdr), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___211.ciunit = *mpi;
	s_wsfe(&io___211);
	do_fio(&c__1, "CALDPF  DFDD   = ", (ftnlen)17);
	do_fio(&c__1, (char *)&(*dfdd), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___212.ciunit = *mpi;
	s_wsfe(&io___212);
	do_fio(&c__1, "CALDPF  DFDS   = ", (ftnlen)17);
	do_fio(&c__1, (char *)&(*dfds), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___213.ciunit = *mpi;
	s_wsfe(&io___213);
	do_fio(&c__1, "CALDPF  DRDS   = ", (ftnlen)17);
	do_fio(&c__1, (char *)&(*drds), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/* -------  END Comment Check Code */
/* ----------------------------------------------------------------------- */
/* ---------   FORMATS  -------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* L101: */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
    return 0;
} /* caldpf_ */


/* Subroutine */ int caldpg_(doublereal *dxdl, doublereal *dxdx, doublereal *
	dxdr, doublereal *dxdd, doublereal *dxds, doublereal *drdl, 
	doublereal *drdr, doublereal *dndl, doublereal *dndx, doublereal *
	dndr, doublereal *dndd, doublereal *dnds, doublereal *dfdr, 
	doublereal *dfrdl, doublereal *dmwx, doublereal *tn, doublereal *
	alpha1, doublereal *dlambda, doublereal *tdl, doublereal *wx, 
	doublereal *untrwx, doublereal *phtrwx, doublereal *phtrwx2, 
	doublereal *rdk, doublereal *rsk, doublereal *teta, integer *ncd, 
	integer *ntens, integer *ndi, integer *nshr, integer *m, integer *mpi)
{
    /* Initialized data */

    static doublereal zero = 0.;
    static doublereal one = 1.;
    static doublereal two = 2.;
    static doublereal three = 3.;
    static doublereal thalf = 1.5;

    /* Format strings */
    static char fmt_102[] = "(1x,a20,6(d20.14,2x))";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static doublereal dxdalpha, twothird;
    static integer i__;
    static doublereal altn[6], sqtt, dndln[6], udmwx, trswx;
    extern /* Subroutine */ int pdtsca_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal rddlrs;

    /* Fortran I/O blocks */
    static cilist io___228 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___229 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___230 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___231 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___232 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___233 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___234 = { 0, 0, 0, fmt_102, 0 };


/* ====================================================================== */
/* ---------------------------------------------------------------------- */
/* ------------  CALCUL DES DERIVEES PARTIELLES DE GX et GR  ------------ */
/* ----------------  PAR RAPPORT DLAMBDA, WX, WR ET ENDO  --------------- */
/* ---------------------------------------------------------------------- */
/* ====================================================================== */
/* ------------- */
/*  EN ENTREE : */
/* ------------- */
/*  DNDL		: DERIVEE PARTIELLE DE N PAR RAPPORT A DLAMBDA */
/*  DNDX		: DERIVEE PARTIELLE DE N PAR RAPPORT A WX */
/*  DNDX		: DERIVEE PARTIELLE DE N PAR RAPPORT A WR */
/*  DNDD		: DERIVEE PARTIELLE DE N PAR RAPPORT A ENDO */
/*  DNDX		: DERIVEE PARTIELLE DE N PAR RAPPORT A LS */
/*  DFDR		: DERIVEE PARTIELLE DE F PAR RAPPORT A WR */
/*  DFRDL	: DERIVEE PARTIELLE DE R PAR RAPPORT A DLAMBDA */
/*  DM		: TENSEUR EFFET D'ENDOMMAGEMENT (Mt) */
/*  DMDDALPH	: PRODUIT CONTRACTE DMDD : ALPHA1 */
/*  DMWX		: NORME DU TENSEUR Mt:(Alpha1+Tdl*N) */
/*  TN		: TENSEUR NORMAL */
/*  ALPHA1	: TENSEUR D'ECROUISSAGE CINEMATIQUE CONNU */
/*  DLAMBDA	: VITESSE DE DEFORMATION PLASTIQUE CUMULEE */
/*  PR1		: VARIABLE D'ECROUISSAGE ISOTROPE CONNUE */
/*  TDL		: TETA*DLAMBDA */
/*  WR		: ISOTROPIC DYNAMIC RECOVERY TERM: WR = DSINH[QS*R] */
/*  WX		: KINEMATIC DYNAMIC RECOVERY TERM: WX = |DMI:X| */
/*  UNTRWX	: 1+TETA*(PRD*DLAMBDA+PRS)*WX */
/*  PHTRWX	: (2/3)*PH/(1+TETA*(PRD*DLAMBDA+PRS)*WX) */
/*  PHTRWX2	: (2/3)*PH/(1+TETA*(PRD*DLAMBDA+PRS)*WX)**2 */
/*  PRD		: CONSTANTE D'ECROUISSAGE CINEMATIQUE */
/*  PRS		: CONSTANTE D'ECROUISSAGE CINEMATIQUE (RECOVERY) */
/*  TETA		: TETA-METHODE */
/*  NCD		: NCD=1 --> COUPLAGE ELASTOPLASTICITE-ENDOMMAGEMENT */
/*  		  NCD=0 --> ELASTOPLASTICITE SANS ENDOMMAGEMENT */
/*  NDI		: NOMBRE DE COMPOSANTES DES CONTRAINTES DIRECTES */
/*  NSHR		: NOMBRE DE COMPOSANTES DES CONTRAINTES TANGENTIELLES */
/*  NTENS	: TAILLE DU TABLEAU DES CONTRAINTES OU DEFORMATIONS */
/* 		  (NDI + NSHR). */
/*  M		: COMMENTAIRES */
/*  MPI		: UNITE DE FICHIER DE SORTIE */
/* ------------- */
/*  EN SORTIE : */
/* ------------- */
/*  DXDL		: DERIVEE PARTIELLE DE GX PAR RAPPORT A DLAMBDA */
/*  DXDX		: DERIVEE PARTIELLE DE GX PAR RAPPORT A WX */
/*  DXDD		: DERIVEE PARTIELLE DE GX PAR RAPPORT A ENDO */
/*  DRDL		: DERIVEE PARTIELLE DE GR PAR RAPPORT A DLAMBDA */
/*  DRDR		: DERIVEE PARTIELLE DE GR PAR RAPPORT A WR */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
/* ----------------------------------------------------------------------- */
/* .1-----  Implicit, External */
/* .2-----  Parameter,Integer */
/* .3-----  Dimension */
/* ------- */
/* .4-----  Real,Complex,Double precision,Logical,Character */
/* .5-----  Common */
/* .6-----  Equivalence */
/* .7-----  Save,Data */
    /* Parameter adjustments */
    --alpha1;
    --tn;
    --dnds;
    --dndd;
    --dndr;
    --dndx;
    --dndl;

    /* Function Body */
/* .8-----  Definition de fonctions */
/* ---------------------------------------------------------------------- */
/* ====================================================================== */
/* ---------------------------------------------------------------------* */
/* -------  Initialisation */

    twothird = two / three;
    sqtt = sqrt(twothird);
/* ------- */
    if (*dmwx == zero) {
	udmwx = one;
    } else {
	udmwx = one / *dmwx;
    }

/* -------  FIN Initialisation */
/* ----------------------------------------------------------------------* */
/* -------  Partial Derivative DXDL = dX/dlambda */

    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
	altn[i__ - 1] = alpha1[i__] + *tdl * tn[i__];
	dndln[i__ - 1] = *dlambda * dndl[i__] + tn[i__];
    }
    pdtsca_(ndi, nshr, altn, dndln, &dxdalpha);
    *dxdl = *teta * *phtrwx2 * (*rdk * *wx * *dmwx - *untrwx * dxdalpha * 
	    udmwx);

/* -------  FIN Partial Derivative DXDL = dGX/dlambda */
/* ---------------------------------------------------------------------* */
/* -------  Partial Derivative DXDX = dX/dWx */

    pdtsca_(ndi, nshr, altn, &dndx[1], &dxdalpha);
    rddlrs = *rdk * *dlambda + *rsk;
    *dxdx = one + *teta * rddlrs * *phtrwx2 * *dmwx - *tdl * *phtrwx * 
	    dxdalpha * udmwx;

/* -------  FIN Partial Derivative DXDX = dX/dWx */
/* ---------------------------------------------------------------------* */
/* -------  Partial Derivative DXDR = dX/dWr */

    pdtsca_(ndi, nshr, altn, &dndr[1], &dxdalpha);
    *dxdr = -(*tdl) * *phtrwx * dxdalpha * udmwx;
    *dxdr = zero;

/* -------  FIN Partial Derivative DXDR = dX/dWr */
/* ---------------------------------------------------------------------* */
/* -------  Partial Derivative DXDD = dX/dD */

    if (*ncd == 0) {
	*dxdd = zero;
    } else {
	pdtsca_(ndi, nshr, altn, &dndd[1], &dxdalpha);
	*dxdd = -(*tdl) * *phtrwx * dxdalpha * udmwx;
    }

/* -------  FIN Partial Derivative DXDD = dX/dD */
/* ---------------------------------------------------------------------* */
/* ---------------------------------------------------------------------* */
/* -------  Partial Derivative (Kinematic Hardening) DXDS = dX/dLs */
    pdtsca_(ndi, nshr, altn, &dnds[1], &dxdalpha);
    trswx = *teta * *rsk * *wx;
    *dxds = *phtrwx2 * (trswx * *dmwx - *untrwx * *tdl * dxdalpha * udmwx);
/* -------  END  Partial Derivative (Kinematic Hardening) DXDS = dX/dLs */
/* ---------------------------------------------------------------------* */
/* -------  Partial Derivative DRDL = dR/dlambda */
    *drdl = -(*dfrdl);
/* -------  FIN Partial Derivative DRDL = dR/dlambda */
/* ---------------------------------------------------------------------* */
/* -------  Partial Derivative DRDR = dR/dD */
    *drdr = one + thalf * *dfdr;
/* -------  FIN Partial Derivative DRDR = dR/dD */
/* ---------------------------------------------------------------------* */
/* ---------------------------------------------------------------------* */
/* -------  Comment Check Code */

    if (*m >= 4) {
	io___228.ciunit = *mpi;
	s_wsfe(&io___228);
	do_fio(&c__1, "CALDPG  DXDL =", (ftnlen)14);
	do_fio(&c__1, (char *)&(*dxdl), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___229.ciunit = *mpi;
	s_wsfe(&io___229);
	do_fio(&c__1, "CALDPG  DXDX =", (ftnlen)14);
	do_fio(&c__1, (char *)&(*dxdx), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___230.ciunit = *mpi;
	s_wsfe(&io___230);
	do_fio(&c__1, "CALDPG  DXDR =", (ftnlen)14);
	do_fio(&c__1, (char *)&(*dxdr), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___231.ciunit = *mpi;
	s_wsfe(&io___231);
	do_fio(&c__1, "CALDPG  DXDD =", (ftnlen)14);
	do_fio(&c__1, (char *)&(*dxdd), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___232.ciunit = *mpi;
	s_wsfe(&io___232);
	do_fio(&c__1, "CALDPG  DXDS =", (ftnlen)14);
	do_fio(&c__1, (char *)&(*dxds), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___233.ciunit = *mpi;
	s_wsfe(&io___233);
	do_fio(&c__1, "CALDPG  DRDL =", (ftnlen)14);
	do_fio(&c__1, (char *)&(*drdl), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___234.ciunit = *mpi;
	s_wsfe(&io___234);
	do_fio(&c__1, "CALDPG  DRDR =", (ftnlen)14);
	do_fio(&c__1, (char *)&(*drdr), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/* -------  FIN Comment Check Code */
/* ----------------------------------------------------------------------* */
/* ---------   FORMATS  -------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* L101: */
/* ----------------------------------------------------------------------* */
/* ======================================================================= */
    return 0;
} /* caldpg_ */


/* Subroutine */ int caldph_(doublereal *dhdl, doublereal *dhdx, doublereal *
	dhdr, doublereal *dhdd, doublereal *dhds, doublereal *dndl, 
	doublereal *dndx, doublereal *dndr, doublereal *dndd, doublereal *
	dnds, doublereal *dpdd, doublereal *zd, doublereal *tn, doublereal *p,
	 doublereal *phi, doublereal *stres1d, doublereal *alpha1, doublereal 
	*altn, doublereal *dlambda, doublereal *tdl, doublereal *h__, 
	doublereal *unmd, doublereal *unmd1, doublereal *wx, doublereal *
	phtrwx, doublereal *phtrwx2, doublereal *g, doublereal *rdk, 
	doublereal *rsk, doublereal *xms, doublereal *xmm1, doublereal *teta, 
	integer *ncd, integer *ntens, integer *ndi, integer *nshr, integer *m,
	 integer *mpi)
{
    /* Initialized data */

    static doublereal zero = 0.;
    static doublereal one = 1.;
    static doublereal two = 2.;
    static doublereal three = 3.;
    static doublereal xnine = 9.;
    static doublereal half = .5;
    static doublereal thalf = 1.5;

    /* Format strings */
    static char fmt_102[] = "(1x,a20,6(d20.14,2x))";
    static char fmt_101[] = "(1x,7(d20.14,2x))";

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal), pow_dd(doublereal *, doublereal *), sinh(
	    doublereal), cosh(doublereal);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static doublereal dxieqvdd, onethird, dxieqvdl, dxieqvdr, dxieqvds, 
	    dxieqvdx, twothird;
    static integer i__;
    static doublereal ad, xi[6], d2g, amm, eyv, xms2, dadd, dcdd, seff, aphi, 
	    cphi, deyv, dh1dd, dh2dd, dh3dd, d2gdl, seff2, dxidd[6], cosha, 
	    dlcay, dndln[6], sinha, dxidl[6], dxidr[6], dxids[6], tgdhx, 
	    dxidx[6], unmdm, xieqv, xieqv2, dlcaxi, dlgdhx, phrdwx, phrswx, 
	    dseffdd, dseffdl, dseffdr, dseffds, dseffdx, phrdsdl;

    /* Fortran I/O blocks */
    static cilist io___290 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___291 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___292 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___293 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___294 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___295 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___296 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___297 = { 0, 0, 0, fmt_102, 0 };


/* ====================================================================== */
/* ---------------------------------------------------------------------- */
/* ---------------  CALCUL DES DERIVEES PARTIELLES DE H  ---------------- */
/* ---------------  PAR RAPPORT DLAMBDA, WX, WR ET ENDO  ---------------- */
/* ---------------------------------------------------------------------- */
/* ====================================================================== */
/* ------------- */
/*  EN ENTREE : */
/* ------------- */
/*  DNDL		: DERIVEE PARTIELLE DE N PAR RAPPORT A DLAMBDA */
/*  DNDX		: DERIVEE PARTIELLE DE N PAR RAPPORT A WX */
/*  DNDD		: DERIVEE PARTIELLE DE N PAR RAPPORT A ENDO */
/*  TN		: TENSEUR NORMAL */
/*  ENDO		: TENSEUR D'ENDOMMAGEMENT */
/*  ALPHA1	: TENSEUR D'ECROUISSAGE CINEMATIQUE CONNU */
/*  ALTN		: TENSEUR = ALPHA1 + TDL*M(-t)*TN */
/*  PR1		: VARIABLE D'ECROUISSAGE ISOTROPE CONNUE */
/*  ZDN		: NORME DU TENSEUR DE DIRECTION NORMALE ZD */
/*  DFDD		: DERIVEE PARTIELLE DE F PAR RAPPORT A ENDO */
/*  DZDNDL	: DERIVEE PARTIELLE DE ZDN PAR RAPPORT A DLAMBDA */
/*  ALPHZD	: PRODUIT CONTRACTE ALPHA1 : ZD */
/*  DM		: TENSEUR D'EFFET D'ENDOMMAGEMENT  (Mt) */
/*  DMI		: TENSEUR D'EFFET D'ENDOMMAGEMENT INVERSE (M(-t)) */
/*  DN		: COMPONANTES DU VECTEUR ORIENTATION (NTENS,NTENS) */
/*  QD		: VOID COALESCENCE OPERATOR */
/*  DQD		: COALESCENCE DERIVATIVE HEAVISIDE */
/*  DLAMBDA	: INCREMENT DEFORMATION PLASTIQUE CUMULEE */
/*  YPS		: VALEUR DE LA FONCTION YPS = (Y/S)**s */
/*  YPS1		: VALEUR DE LA FONCTION YPS = (Y/S)**(s-1) */
/*  YV		: VALEUR DE LA FONCTION YV = Ye = 3/2*<Sh>/Seq */
/*  EYV		: VALEUR DE LA FONCTION Exp[YV] */
/*  UES		: VALEUR DE LA FONCTION Exp[YV]/Seq (if <Sh>.neq.0) */
/*  TDL		: TETA*DLAMBDA */
/*  UDHM		: VALUE 1/(1-Dh)**m */
/*  UDHM1	: VALUE 1/(1-Dh)**(m+1) */
/*  WR		: ISOTROPIC DYNAMIC RECOVERY TERM: WR = DSINH[QS*R] */
/*  WX		: KINEMATIC DYNAMIC RECOVERY TERM: WX = |DMI:X| */
/*  UNTRWX	: 1+TETA*(PRD*DLAMBDA+PRS)*WX */
/*  PHTRWX2	: (2/3)*PH/(1+TETA*(PRD*DLAMBDA+RT(2/3)*PRS)*WX)**2 */
/*  PHTRWX3	: (2/3)*PH/(1+TETA*(PRD*DLAMBDA+RT(2/3)*PRS)*WX)**3 */
/*  HTRWR3	: RT(2/3)*H/(1+TETA*(RD*DLAMBDA+RS)*WR)**3 */
/*  RTDL		: PR1+TETA*DLAMBDA */
/*  TRSWR	: TETA*RS*WR */
/*  G		: MODULE DE CISAILLEMENT */
/*  XK		: COEFFICIENT DE LAME */
/*  RD		: ISOTROPIC HARDENING DYNAMIC RECOVERY CONSTANT */
/*  H		: ISOTROPIC HARDENING MODULUS */
/*  RS		: ISOTROPIC HARDENING STATIC RECOVERY CONSTANT */
/*  PRD		: KINEMATIC HARDENING DYNAMIC RECOVERY CONSTANT */
/*  PH		: KINEMATIC HARDENING MODULUS */
/*  PRS		: KINEMATIC HARDENING STATIC RECOVERY CONSTANT */
/*  GS		: FONCTION D'ENDOMMAGEMENT DUCTILE */
/*  PS		: FONCTION D'ENDOMMAGEMENT DUCTILE */
/*  XDN		: CONSTANTE NUCLEATION XLD*XFN/(XKIC*XLF**(1/3)) */
/*  CTE		: CONSTANTE DE TENSION */
/*  CCO		: CONSTANTE DE COMPRESSION */
/*  CTO		: CONSTANTE DE TORSION */
/*  AV		: VOID GROWTH CONSTANT */
/*  VM		: VOID GROWTH CONSTANT */
/*  VP		: VOID GROWTH TEST (VC = 1 if Dh > Vc) */
/*  AC		: COALESCENCE GROWTH CONSTANT */
/*  TETA		: TETA-METHODE */
/*  MELPST	: MELPST = 1 --> TN = ZD/ZDN (RESOLUTION) */
/*                 SINON TN = ZD/ZDZD (MATRICE CONSISTANTE) */
/*  NDI		: NOMBRE DE COMPOSANTES DES CONTRAINTES DIRECTES */
/*  NSHR		: NOMBRE DE COMPOSANTES DES CONTRAINTES TANGENTIELLES */
/*  NTENS	: TAILLE DU TABLEAU DES CONTRAINTES OU DEFORMATIONS */
/* 		  (NDI + NSHR). */
/*  M		: COMMENTAIRES */
/*  MPI		: UNITE DE FICHIER DE SORTIE */
/* ------------- */
/*  EN SORTIE : */
/* ------------- */
/*  DHDL		: DERIVEE PARTIELLE DE H PAR RAPPORT A DLAMBDA */
/*  DHDX		: DERIVEE PARTIELLE DE H PAR RAPPORT A WX */
/*  DHDR		: DERIVEE PARTIELLE DE H PAR RAPPORT A WR */
/*  DHDD		: DERIVEE PARTIELLE DE H PAR RAPPORT A ENDO */
/*  UDEVNN	: TENSEUR (D'ORDRE 4) = Idev -(2/3) N*N */
/*  CSTRESD	: PRODUCT StresD : [ M : Idev : Mt ] */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
/* ----------------------------------------------------------------------- */
/* .1-----  Implicit, External */
/* .2-----  Parameter,Integer */
/* .3-----  Dimension */
/* ------- */
/* -----------------------------------------------------------* */
/* .4-----  Real,Complex,Double precision,Logical,Character */
/* .5-----  Common */
/* .6-----  Equivalence */
/* .7-----  Save,Data */
    /* Parameter adjustments */
    --h__;
    --altn;
    --alpha1;
    --stres1d;
    --tn;
    --zd;
    --dnds;
    --dndd;
    --dndr;
    --dndx;
    --dndl;

    /* Function Body */
/* .8-----  Definition de fonctions */
/* ---------------------------------------------------------------------* */
/* ====================================================================== */
/* ---------------------------------------------------------------------* */
/* ------- Initialization */
    d2g = two * *g;
    onethird = one / three;
    twothird = two / three;
/* ---------------------------------------------------------------------* */
    if (*ncd == 1 || *ncd == 2) {
/* -------------------------------------------* */
/* ---------  Tensors */
	tgdhx = *teta * (d2g * *unmd + *phtrwx);
	dlgdhx = tgdhx * *dlambda;
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dndln[i__ - 1] = *tdl * dndl[i__] + tn[i__];
	    altn[i__] = alpha1[i__] + *tdl * tn[i__];
	    xi[i__ - 1] = zd[i__] - dlgdhx * tn[i__];
	}
/* -------------------------------------------* */
/* ---------  Norm |Xi| = Sqrt(3/2 * Xi:Xi) */
	xieqv = xi[0] * xi[0] + xi[1] * xi[1] + xi[2] * xi[2];
	i__1 = *ntens;
	for (i__ = *ndi + 1; i__ <= i__1; ++i__) {
	    xieqv += two * xi[i__ - 1] * xi[i__ - 1];
	}
	xieqv = sqrt(thalf * xieqv);
	seff = sqrt(h__[1] * xieqv * xieqv + h__[2] * *p * *p);
/* -------------------------------------------* */
/* ---------  Partial Derivative DXIDx = dXi/dx */
	d2gdl = d2g * *tdl;
	phrdwx = *phtrwx2 * *rdk * *wx;
	phrswx = *phtrwx2 * *rsk * *wx;
	phrdsdl = *phtrwx2 * (*rdk * *dlambda + *rsk);
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dxidl[i__ - 1] = phrdwx * altn[i__] - tgdhx * dndln[i__ - 1];
	    dxidx[i__ - 1] = phrdsdl * altn[i__] - dlgdhx * dndx[i__];
	    dxidr[i__ - 1] = -dlgdhx * dndr[i__];
	    dxidd[i__ - 1] = -stres1d[i__] + d2gdl * tn[i__] - dlgdhx * dndd[
		    i__];
	    dxids[i__ - 1] = phrswx * altn[i__] - dlgdhx * dnds[i__];
	}
/* -------------------------------------------* */
/* ---------  Partial Derivative of DXIEQVDx = dXieq/dx (1st part) */
	dxieqvdl = xi[0] * dxidl[0] + xi[1] * dxidl[1] + xi[2] * dxidl[2];
	dxieqvdx = xi[0] * dxidx[0] + xi[1] * dxidx[1] + xi[2] * dxidx[2];
	dxieqvdr = xi[0] * dxidr[0] + xi[1] * dxidr[1] + xi[2] * dxidr[2];
	dxieqvdd = xi[0] * dxidd[0] + xi[1] * dxidd[1] + xi[2] * dxidd[2];
	dxieqvds = xi[0] * dxids[0] + xi[1] * dxids[1] + xi[2] * dxids[2];
	i__1 = *ntens;
	for (i__ = *ndi + 1; i__ <= i__1; ++i__) {
	    dxieqvdl += two * xi[i__ - 1] * dxidl[i__ - 1];
	    dxieqvdx += two * xi[i__ - 1] * dxidx[i__ - 1];
	    dxieqvdr += two * xi[i__ - 1] * dxidr[i__ - 1];
	    dxieqvdd += two * xi[i__ - 1] * dxidd[i__ - 1];
	    dxieqvds += two * xi[i__ - 1] * dxids[i__ - 1];
	}
    }
/* -----------------------------------------------------------* */
/* -----------------------------------------------------------* */
/* --------- Choice of Damage Evolution */
    if (*ncd == 1) {
/* -----------------------------------------------------------* */
/* ---------  Cocks [1989] */
/* -----------------------------------------------------------* */
	eyv = onethird * *p / seff;
	*xmm1 = *xms / (*xms + one);
	ad = two * (*xms + one) * (one + *phi) * pow_dd(unmd, xmm1);
	aphi = xnine * *xms * *phi * *unmd * eyv / ad;
/* -----------------------------------------------------------* */
	dh1dd = twothird;
	dh2dd = half * *xmm1 / (*unmd1 * *unmd1);
	d__1 = *xmm1 - one;
	dh3dd = *xmm1 * pow_dd(unmd1, &d__1);
	dh1dd = zero;
	dh2dd = zero;
	dh3dd = zero;
/* -----------------------------------------------------------* */
/* ---------  Partial Derivative of DSEFFDx = d[1/Seff]/dx */
	xieqv2 = xieqv * xieqv;
	seff2 = seff * seff;
	dseffdl = thalf * h__[1] * dxieqvdl / seff2;
	dseffdx = thalf * h__[1] * dxieqvdx / seff2;
	dseffdr = thalf * h__[1] * dxieqvdr / seff2;
	dseffdd = half * (dh1dd * xieqv2 + dh2dd * *p * *p) + thalf * h__[1] *
		 dxieqvdd + h__[2] * *p * *dpdd;
	dseffdd /= seff2;
	dseffds = thalf * h__[1] * dxieqvds / seff2;
/* -------------------------------------------* */
/* ---------  Partial Derivative DHDL = dFphi/dlambda */
	*dhdl = -aphi * (one - *dlambda * dseffdl);
/* -------------------------------------------* */
/* ---------  Partial Derivative DHDK = dFphi/dWx */
	*dhdx = aphi * *dlambda * dseffdx;
/* -------------------------------------------* */
/* ---------  Partial Derivative DHDB = dFphi/dWb */
	*dhdr = aphi * *dlambda * dseffdr;
/* -------------------------------------------* */
/* ---------  Partial Derivative DHDD = dFphi/dphi */
	dadd = (one + (*xmm1 - two) * *phi) / (*phi * *unmd) - one / (one + *
		phi);
	deyv = onethird * *dpdd / seff - eyv * dseffdd;
	*dhdd = aphi * dadd + xnine * *xms * *phi * *unmd * deyv / ad;
	*dhdd = one - *dlambda * *dhdd;
/* -------------------------------------------* */
/* ---------  Partial Derivative DHDS = dFphi/dls [Eq. 203] */
	*dhds = aphi * *dlambda * dseffds;
/* -----------------------------------------------------------* */
    } else if (*ncd == 2) {
/* -----------------------------------------------------------* */
/* ---------  Cocks and Ashby [1980] */
/* -----------------------------------------------------------* */
	eyv = onethird * *p / xieqv;
	eyv = (eyv > 0.0) ? eyv : 0.0; /* Macauley bracket */
	xms2 = two * *xms;
	unmdm = pow_dd(unmd, xms);
	amm = two * (xms2 - one) / (xms2 + one);
	aphi = amm * eyv;
	sinha = sinh(aphi);
	cphi = one / unmdm - *unmd;
/* -----------------------------------------------------------* */
/* ---------  Partial Derivative of DXIEQVDx = d[1/Xieq]/dx * Xieqv */
	xieqv2 = xieqv * xieqv;
	dxieqvdl = -thalf * dxieqvdl / xieqv2;
	dxieqvdx = -thalf * dxieqvdx / xieqv2;
	dxieqvdr = -thalf * dxieqvdr / xieqv2;
	dxieqvdd = -thalf * dxieqvdd / xieqv2;
	dxieqvds = -thalf * dxieqvds / xieqv2;
/* -------------------------------------------* */
	cosha = cosh(aphi);
	dlcay = *dlambda * cphi * cosha * amm * eyv;
	dlcaxi = *dlambda * cphi * cosha * amm / xieqv;
/* -------------------------------------------* */
/* ---------  Partial Derivative DHDL = dFphi/dlambda */
	*dhdl = -cphi * sinha - dlcay * dxieqvdl;
/* -------------------------------------------* */
/* ---------  Partial Derivative DHDX = dFphi/dWx */
	*dhdx = -dlcay * dxieqvdx;
/* -------------------------------------------* */
/* ---------  Partial Derivative DHDB = dFphi/dWr */
	*dhdr = -dlcay * dxieqvdr;
/* -------------------------------------------* */
/* ---------  Partial Derivative DHDD = dFphi/dphi */
	d__1 = *xms + one;
	dcdd = one + *xms * pow_dd(unmd1, &d__1);
	*dhdd = one - *dlambda * sinha * dcdd - dlcay * dxieqvdd - onethird * 
		dlcaxi * *dpdd;
/* -------------------------------------------* */
/* ---------  Partial Derivative DHDB = dFphi/dLs */
	*dhds = -dlcay * dxieqvds;
	if (*m >= 0)
	{
	io___290.ciunit = *mpi;
	s_wsfe(&io___290);
	do_fio(&c__1, "CALDPH DPDD = ", (ftnlen)14);
	do_fio(&c__1, (char *)&(*dpdd), (ftnlen)sizeof(doublereal));
	e_wsfe();
	}
/* -----------------------------------------------------------* */
    } else {
/* -----------------------------------------------------------* */
/* ---------  Plasticity without Damage */
/* -----------------------------------------------------------* */
	*dhdl = zero;
	*dhdx = zero;
	*dhdr = zero;
	*dhdd = zero;
	*dhds = zero;
/* -----------------------------------------------------------* */
    }
/* -----------------------------------------------------------* */
/* -------  END Partial Derivatives dFgamma/dphi  ----------------------* */
/* ---------------------------------------------------------------------* */
/* -------  Comment Check Code */
    if (*m >= 4) {
	if (*ncd == 1 || *ncd == 2) {
	    io___291.ciunit = *mpi;
	    s_wsfe(&io___291);
	    do_fio(&c__1, "STRESS XI", (ftnlen)9);
	    e_wsfe();
	    i__1 = *ntens;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		io___292.ciunit = *mpi;
		s_wsfe(&io___292);
		do_fio(&c__1, (char *)&xi[i__ - 1], (ftnlen)sizeof(doublereal)
			);
		e_wsfe();
	    }
	}
	io___293.ciunit = *mpi;
	s_wsfe(&io___293);
	do_fio(&c__1, "CALDPH DHDL = ", (ftnlen)14);
	do_fio(&c__1, (char *)&(*dhdl), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___294.ciunit = *mpi;
	s_wsfe(&io___294);
	do_fio(&c__1, "CALDPH DHDX = ", (ftnlen)14);
	do_fio(&c__1, (char *)&(*dhdx), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___295.ciunit = *mpi;
	s_wsfe(&io___295);
	do_fio(&c__1, "CALDPH DHDR = ", (ftnlen)14);
	do_fio(&c__1, (char *)&(*dhdr), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___296.ciunit = *mpi;
	s_wsfe(&io___296);
	do_fio(&c__1, "CALDPH DHDD = ", (ftnlen)14);
	do_fio(&c__1, (char *)&(*dhdd), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___297.ciunit = *mpi;
	s_wsfe(&io___297);
	do_fio(&c__1, "CALDPH DHDS = ", (ftnlen)14);
	do_fio(&c__1, (char *)&(*dhds), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
/* -------  END Comment Check Code */
/* ----------------------------------------------------------------------- */
/* ---------   FORMATS  -------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
    return 0;
} /* caldph_ */


/* Subroutine */ int caldpn_(doublereal *dndl, doublereal *dndx, doublereal *
	dndr, doublereal *dndd, doublereal *dnds, doublereal *zd, doublereal *
	alpha1, doublereal *stres1d, doublereal *zdzd, doublereal *zdn, 
	doublereal *dzdndl, doublereal *dzdndd, doublereal *dfdr, doublereal *
	dfdd, doublereal *alphzd, doublereal *unmd, doublereal *dlambda, 
	doublereal *tdl, doublereal *wx, doublereal *phtrwx, doublereal *
	phtrwx2, doublereal *rdk, doublereal *rsk, doublereal *teta, integer *
	ncd, integer *melpst, integer *ntens, integer *ndi, integer *nshr, 
	integer *m, integer *mpi)
{
    /* Initialized data */

    static doublereal zero = 0.;
    static doublereal one = 1.;
    static doublereal two = 2.;
    static doublereal three = 3.;

    /* Format strings */
    static char fmt_101[] = "(1x,6(d20.14,2x))";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void), s_wsfe(cilist *), do_fio(integer *, char *, ftnlen),
	     e_wsfe(void);

    /* Local variables */
    static doublereal thdlrpzd, twothird, trdhwxzd;
    static integer i__;
    static doublereal trshwxzd, trdhwazd3, thdlrpzd3, trshwazd3, uzd;
    extern /* Subroutine */ int aset_(doublereal *, doublereal *, integer *);
    static doublereal uzdn, uzdn2, uzdzd;
    extern /* Subroutine */ int pdtsca_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal thdlrp, trdhwx, dzdndd2, trshwx, dzdndl2, alphzd3, 
	    dzdzdd3, dzdndr2, dzdndx2, streszd;

    /* Fortran I/O blocks */
    static cilist io___326 = { 0, 0, 0, 0, 0 };
    static cilist io___327 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___328 = { 0, 0, 0, 0, 0 };
    static cilist io___329 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___330 = { 0, 0, 0, 0, 0 };
    static cilist io___331 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___332 = { 0, 0, 0, 0, 0 };
    static cilist io___333 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___334 = { 0, 0, 0, 0, 0 };
    static cilist io___335 = { 0, 0, 0, fmt_101, 0 };


/* ====================================================================== */
/* ---------------------------------------------------------------------- */
/* ---------------  CALCUL DES DERIVEES PARTIELLES DE N  ---------------- */
/* ---------------  PAR RAPPORT DLAMBDA, WX, XR et ENDO  ---------------- */
/* ---------------------------------------------------------------------- */
/* ====================================================================== */
/* ------------- */
/*  EN ENTREE : */
/* ------------- */
/*  ZD		: TENSEUR DE DIRECTION NORMALE */
/*  TN		: TENSEUR NORMAL */
/*  ESTAR	: TENSEUR DES DEFORMATIONS CONNUS A T+DT */
/*  ALPHA1	: TENSEUR D'ECROUISSAGE CINEMATIQUE CONNU */
/*  TETN		: TENSEUR = ESTAR - TDLRD*TN */
/*  ALTN		: TENSEUR = ALPHA1 + TDLRD*TN */
/*  PR1		: VARIABLE D'ECROUISSAGE ISOTROPE CONNUE */
/*  ZDN		: NORME DU TENSEUR DE DIRECTION NORMALE ZD */
/*  DFDD		: DERIVEE PARTIELLE DE F PAR RAPPORT A ENDO */
/*  DZDNDL	: DERIVEE PARTIELLE DE ZDN PAR RAPPORT A DLAMBDA */
/*  ALPHZD	: PRODUIT CONTRACTE ALPHA1 : ZD */
/*  DM		: TENSEUR D'EFFET D'ENDOMMAGEMENT (Mt) */
/*  DMI		: TENSEUR D'EFFET D'ENDOMMAGEMENT INVERSE (M(-t)) */
/*  DMDDESTAR	: PRODUIT CONTRACTE IDEV : DMDD : ESTAR */
/*  DMDDALPH	: PRODUIT CONTRACTE DMDD : ALPHA1 */
/*  DLAMBDA	: INCREMENT DEFORMATION PLASTIQUE CUMULEE */
/*  YPS		: VALEUR DE LA FONCTION YPS = (Y/S)**s */
/*  YPS1		: VALEUR DE LA FONCTION YPS = (Y/S)**(s-1) */
/*  TDL		: TETA*DLAMBDA */
/*  PHTRWX	: (2/3)*PH/(1+TETA*(PRD*DLAMBDA+RT(2/3)*PRS)*WX) */
/*  PHTRWX2	: (2/3)*PH/(1+TETA*(PRD*DLAMBDA+RT(2/3)*PRS)*WX)**2 */
/*  G		: MODULE DE CISAILLEMENT */
/*  XK		: COEFFICIENT DE LAME */
/*  PRD		: CONSTANTE D'ECROUISSAGE CINEMATIQUE */
/*  PRS		: CONSTANTE D'ECROUISSAGE CINEMATIQUE (RECOVERY) */
/*  NCD		: NCD=1 --> COUPLAGE ELASTOPLASTICITE-ENDOMMAGEMENT */
/*  		  NCD=0 --> ELASTOPLASTICITE SANS ENDOMMAGEMENT */
/*  TETA		: TETA-METHODE */
/*  MELPST	: MELPST = 1 --> TN = ZD/ZDN (RESOLUTION) */
/*                 SINON TN = ZD/ZDZD (MATRICE CONSISTANTE) */
/*  NDI		: NOMBRE DE COMPOSANTES DES CONTRAINTES DIRECTES */
/*  NSHR		: NOMBRE DE COMPOSANTES DES CONTRAINTES TANGENTIELLES */
/*  NTENS	: TAILLE DU TABLEAU DES CONTRAINTES OU DEFORMATIONS */
/* 		  (NDI + NSHR). */
/*  M		: COMMENTAIRES */
/*  MPI		: UNITE DE FICHIER DE SORTIE */
/* ------------- */
/*  EN SORTIE : */
/* ------------- */
/*  DNDL		: DERIVEE PARTIELLE DE N PAR RAPPORT A DLAMBDA */
/*  DNDX		: DERIVEE PARTIELLE DE N PAR RAPPORT A WX */
/*  DNDD		: DERIVEE PARTIELLE DE N PAR RAPPORT A ENDO */
/*  DNDS		: DERIVEE PARTIELLE DE N PAR RAPPORT A DLs */
/*  DZDDD	: DERIVEE PARTIELLE DE ZD PAR RAPPORT A ENDO */
/*  UDEVNN	: TENSEUR (D'ORDRE 4) = Idev -(2/3) N*N */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
/* .1-----  Implicit, External */
/* .2-----  Parameter,Integer */
/* .3-----  Dimension */
/* ------- */
/* .4-----  Real,Integer,Complex,Double precision,Logical,Character */
/* .5-----  Common */
/* .6-----  Equivalence */
/* .7-----  Save,Data */
    /* Parameter adjustments */
    --stres1d;
    --alpha1;
    --zd;
    --dnds;
    --dndd;
    --dndr;
    --dndx;
    --dndl;

    /* Function Body */
/* .8-----  Definition de fonctions */
/* ---------------------------------------------------------------------* */
/* ====================================================================== */
/* ---------------------------------------------------------------------* */
/* -------  Initialisation */
    twothird = two / three;

/* -------  FIN Initialisation */
/* ---------------------------------------------------------------------* */
/* -------  Calcul des derivees partielles dN/dlambda, dN/dWx et DNDD=dN/dD */

    trdhwx = *teta * *rdk * *wx * *phtrwx2;
    thdlrp = *teta * (*dlambda * *rdk + *rsk) * *phtrwx2;
/* ------- */
    if (*melpst == 1) {
/* ------- */
	uzdn = one / *zdn;
	uzdn2 = uzdn * uzdn;
	dzdndl2 = *dzdndl * uzdn2;
	trdhwxzd = uzdn * trdhwx;
	thdlrpzd = uzdn * thdlrp;
	dzdndx2 = uzdn * thdlrpzd * *tdl;
	dzdndr2 = *dfdr * uzdn2;
	dzdndd2 = *dzdndd * uzdn2;
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dndl[i__] = trdhwxzd * alpha1[i__] - dzdndl2 * zd[i__];
	    dndx[i__] = thdlrpzd * alpha1[i__] + dzdndx2 * zd[i__];
	    dndr[i__] = dzdndr2 * zd[i__];
	    dndd[i__] = -uzdn * stres1d[i__] - dzdndd2 * zd[i__];
	    dnds[i__] = zero;
	}
/* ------- */
    } else {
/* ------- */
	pdtsca_(ndi, nshr, &stres1d[1], &zd[1], &streszd);
	uzd = one / *zdzd;
	uzdzd = uzd * uzd;
	alphzd3 = twothird * *alphzd * uzdzd * uzdzd;
	trdhwxzd = uzdzd * trdhwx;
	trdhwazd3 = trdhwxzd * alphzd3;
	thdlrpzd = uzdzd * thdlrp;
	thdlrpzd3 = thdlrp * alphzd3;
	trshwx = *teta * *rsk * *wx * *phtrwx2;
	trshwxzd = uzdzd * trshwx;
	trshwazd3 = trshwxzd * alphzd3;
	dzdzdd3 = uzdzd * uzdzd * uzdzd * streszd;
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dndl[i__] = trdhwxzd * alpha1[i__] - trdhwazd3 * zd[i__];
	    dndx[i__] = thdlrpzd * alpha1[i__] - thdlrpzd3 * zd[i__];
	    dndr[i__] = zero;
	    dndd[i__] = -uzdzd * stres1d[i__] + dzdzdd3 * zd[i__];
	    dnds[i__] = trshwxzd * alpha1[i__] - trshwazd3 * zd[i__];
	}
/* ------- */
    }
/* ----------  Tensor DNDD if Ncd = 0 */
    if (*ncd == 0) {
	aset_(&dndd[1], &zero, ntens);
    }

/* -------  END Calcul des derivees partielles dN/dlambda, dN/dWx et dN/dD */
/* ---------------------------------------------------------------------* */
/* ---------------------------------------------------------------------* */
/* -------  Comment Check Code */

    if (*m >= 4) {
	io___326.ciunit = *mpi;
	s_wsle(&io___326);
	do_lio(&c__9, &c__1, "CALDPN DNDL", (ftnlen)11);
	e_wsle();
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___327.ciunit = *mpi;
	    s_wsfe(&io___327);
	    do_fio(&c__1, (char *)&dndl[i__], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	io___328.ciunit = *mpi;
	s_wsle(&io___328);
	do_lio(&c__9, &c__1, "CALDPN DNDX", (ftnlen)11);
	e_wsle();
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___329.ciunit = *mpi;
	    s_wsfe(&io___329);
	    do_fio(&c__1, (char *)&dndx[i__], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	if (*melpst == 1) {
	    io___330.ciunit = *mpi;
	    s_wsle(&io___330);
	    do_lio(&c__9, &c__1, "CALDPN DNDR", (ftnlen)11);
	    e_wsle();
	    i__1 = *ntens;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		io___331.ciunit = *mpi;
		s_wsfe(&io___331);
		do_fio(&c__1, (char *)&dndr[i__], (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	}
	io___332.ciunit = *mpi;
	s_wsle(&io___332);
	do_lio(&c__9, &c__1, "CALDPN DNDD", (ftnlen)11);
	e_wsle();
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___333.ciunit = *mpi;
	    s_wsfe(&io___333);
	    do_fio(&c__1, (char *)&dndd[i__], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	if (*melpst == 0) {
	    io___334.ciunit = *mpi;
	    s_wsle(&io___334);
	    do_lio(&c__9, &c__1, "CALDPN DNDS", (ftnlen)11);
	    e_wsle();
	    i__1 = *ntens;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		io___335.ciunit = *mpi;
		s_wsfe(&io___335);
		do_fio(&c__1, (char *)&dnds[i__], (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	}
    }

/* -------  END Comment Check Code */
/* ----------------------------------------------------------------------* */
/* ---------   FORMATS  -------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* L102: */
/* ----------------------------------------------------------------------* */
/* ======================================================================= */
    return 0;
} /* caldpn_ */


/* Subroutine */ int caldsig_(doublereal *ddsdde, doublereal *dnde, 
	doublereal *tnde, doublereal *snddde, doublereal *dlde, doublereal *
	ddde, doublereal *tn, doublereal *strant, doublereal *stres1d, 
	doublereal *dlambda, doublereal *unmd, doublereal *unmd1, doublereal *
	g, doublereal *xk, integer *ntens, integer *ndi, integer *m, integer *
	mpi)
{
    /* Initialized data */

    static doublereal udemi = .5;
    static doublereal deux = 2.;
    static doublereal trois = 3.;

    /* Format strings */
    static char fmt_101[] = "(1x,6(d12.5,2x))";

    /* System generated locals */
    integer ddsdde_dim1, ddsdde_offset, snddde_dim1, snddde_offset, dnde_dim1,
	     dnde_offset, tnde_dim1, tnde_offset, i__1, i__2;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void), s_wsfe(cilist *), do_fio(integer *, char *, ftnlen),
	     e_wsfe(void);

    /* Local variables */
    static doublereal xkstranh;
    static integer i__, j;
    static doublereal d2g, d2gdl, s1ddn[6];
    extern /* Subroutine */ int trace_(doublereal *, doublereal *, integer *, 
	    integer *);
    static doublereal unmd2g;
    extern /* Subroutine */ int pdtten_(doublereal *, doublereal *, 
	    doublereal *, integer *, integer *);
    static doublereal stranh, unmdxkg;

    /* Fortran I/O blocks */
    static cilist io___348 = { 0, 0, 0, 0, 0 };
    static cilist io___349 = { 0, 0, 0, fmt_101, 0 };


/* ====================================================================== */
/* ---------------------------------------------------------------------- */
/* ---------------------  CALCUL DU TENSEUR DDSDDE  --------------------- */
/* ---------------------------------------------------------------------- */
/* ====================================================================== */
/* ------------- */
/*  EN ENTREE : */
/* ------------- */
/*  DNDE		: TENSEUR (D'ORDRE 4) = DN/DE */
/*  TNDE		: TENSEUR (D'ORDRE 4) = TN*DLDE */
/*  SNDDDE 	: TENSEUR (D'ORDRE 4) = (XK*EH*I+S1D+2G*DLRD*TN)*DDDE */
/*  DDDE		: TENSEUR (D'ORDRE 2) = DD/DE */
/*  STRANT	: TENSEUR DES DEFORMATIONS TOTALES */
/*  TN		: TENSEUR NORMAL */
/*  STRES1D	: TENSEUR DES CONTRAINTES CONNU */
/*  DLAMBDA	: INCREMENT DEFORMATION PLASTIQUE CUMULEE */
/*  UNMD		: 1-D */
/*  RCUNMD	: RACINE(1-D) */
/*  DLRD		: DLAMBDA/RACINE(1-D) */
/*  G		: MODULE DE CISAILLEMENT */
/*  XK		: COEFFICIENT DE LAME */
/*  NDI		: NOMBRE DE COMPOSANTES DES CONTRAINTES DIRECTES */
/*  NTENS	: NOMBRE TOTAL DE COMPOSANTES D'UN TENSEUR */
/*  M		: COMMENTAIRES */
/*  MPI		: UNITE DE FICHIER DE SORTIE */
/* ------------- */
/*  EN LOCAL  : */
/* ------------- */
/*  ALDE		: TENSEUR (D'ORDRE 4) = TN*DLDE */
/* ------------- */
/*  EN SORTIE : */
/* ------------- */
/*  DDSDDE	: MATRICE TANGENTE CONSISTENTE */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
/* .1-----  Implicit, External */
/* .2-----  Parameter,Integer */
/* .3-----  Dimension */
/* .4-----  Real,Complex,Double precision,Logical,Character */
/* .5-----  Common */
/* .6-----  Equivalence */
/* .7-----  Save,Data */
    /* Parameter adjustments */
    --stres1d;
    --strant;
    --tn;
    --ddde;
    --dlde;
    snddde_dim1 = *ntens;
    snddde_offset = 1 + snddde_dim1;
    snddde -= snddde_offset;
    tnde_dim1 = *ntens;
    tnde_offset = 1 + tnde_dim1;
    tnde -= tnde_offset;
    dnde_dim1 = *ntens;
    dnde_offset = 1 + dnde_dim1;
    dnde -= dnde_offset;
    ddsdde_dim1 = *ntens;
    ddsdde_offset = 1 + ddsdde_dim1;
    ddsdde -= ddsdde_offset;

    /* Function Body */
/* .8-----  Definition de fonctions */
/* ---------------------------------------------------------------------- */
/* ====================================================================== */
/* ---------------------------------------------------------------------- */
/* -------  Initialization */

    d2g = deux * *g;
    unmd2g = *unmd * d2g;
    unmdxkg = *unmd * (*xk - d2g / trois);
    d2gdl = d2g * *dlambda;
    trace_(&strant[1], &stranh, ndi, ntens);
    xkstranh = *xk * stranh;

/* -------  END Initialisation */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* -------  Consistent Tangent DDSDDE */

    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s1ddn[i__ - 1] = stres1d[i__] - d2gdl * tn[i__];
    }
    i__1 = *ndi;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s1ddn[i__ - 1] += xkstranh;
    }
/* ------- */
    pdtten_(s1ddn, &ddde[1], &snddde[snddde_offset], ntens, ndi);
    pdtten_(&tn[1], &dlde[1], &tnde[tnde_offset], ntens, ndi);
/* ------- */
    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *ntens;
	for (j = 1; j <= i__2; ++j) {
	    ddsdde[i__ + j * ddsdde_dim1] = -snddde[i__ + j * snddde_dim1] - 
		    unmd2g * (tnde[i__ + j * tnde_dim1] + *dlambda * dnde[i__ 
		    + j * dnde_dim1]);
	}
	ddsdde[i__ + i__ * ddsdde_dim1] += unmd2g;
    }
/* ------- */
    i__1 = *ndi;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *ndi;
	for (j = 1; j <= i__2; ++j) {
	    ddsdde[i__ + j * ddsdde_dim1] += unmdxkg;
	}
    }
/* ------ */
    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *ntens;
	for (j = *ndi + 1; j <= i__2; ++j) {
	    ddsdde[i__ + j * ddsdde_dim1] = udemi * ddsdde[i__ + j * 
		    ddsdde_dim1];
	}
    }

/* -------  END Consistent Tangent DDSDDE */
/* ---------------------------------------------------------------------* */
/* ---------------------------------------------------------------------* */
/* -------  Check Comment Code */

    if (*m >= 4) {
	io___348.ciunit = *mpi;
	s_wsle(&io___348);
	do_lio(&c__9, &c__1, "CALDISG DDSDDE", (ftnlen)14);
	e_wsle();
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___349.ciunit = *mpi;
	    s_wsfe(&io___349);
	    i__2 = *ntens;
	    for (j = 1; j <= i__2; ++j) {
		do_fio(&c__1, (char *)&ddsdde[i__ + j * ddsdde_dim1], (ftnlen)
			sizeof(doublereal));
	    }
	    e_wsfe();
	}
    }

/* -------  END Check Comment Code */
/* ----------------------------------------------------------------------- */
/* ---------   FORMATS  -------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* L102: */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
    return 0;
} /* caldsig_ */


/* Subroutine */ int calfgh_(doublereal *ff, doublereal *fx, doublereal *fr, 
	doublereal *fd, doublereal *zdn, doublereal *zdzd, doublereal *tn, 
	doublereal *zd, doublereal *alpha1, doublereal *phi, doublereal *
	phiold, doublereal *unmd, doublereal *unmd1, doublereal *h__, 
	doublereal *p, doublereal *dlambda, doublereal *tdl, doublereal *eyv, 
	doublereal *wx, doublereal *dmwx, doublereal *phtrwx, doublereal *wr, 
	doublereal *rtdl, doublereal *htrwr, doublereal *g, doublereal *xms, 
	integer *ncd, integer *ntens, integer *ndi, integer *nshr, integer *m,
	 integer *mpi)
{
    /* Initialized data */

    static doublereal zero = 0.;
    static doublereal one = 1.;
    static doublereal two = 2.;
    static doublereal three = 3.;
    static doublereal xnine = 9.;
    static doublereal thalf = 1.5;

    /* Format strings */
    static char fmt_102[] = "(1x,a20,6(d20.14,2x))";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal), pow_dd(doublereal *, doublereal *), sinh(
	    doublereal);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static doublereal onethird;
    static integer i__;
    static doublereal xi[6], d2g, amm, xmm1, xms2, aphi, seff, cphi, altn[6], 
	    sinha, unmdm, xieqv, dlgdhx;

    /* Fortran I/O blocks */
    static cilist io___372 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___373 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___374 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___375 = { 0, 0, 0, fmt_102, 0 };


/* ====================================================================== */
/* ---------------------------------------------------------------------- */
/* ---------------  CALCUL DES FONCTIONS F, GX, GR ET H  ---------------- */
/* ---------------------------------------------------------------------- */
/* ====================================================================== */
/* ------------- */
/*  EN ENTREE : */
/* ------------- */
/*  ENDO 	: VARIABLE D'ENDOMMAGEMENT */
/*  ZDN		: NORME DU TENSEUR DE DIRECTION NORMALE ZD */
/*  ZDZD		: NORME DU TENSEUR DE DIRECTION NORMALE ZD */
/*  ENDO1	: VARIABLE D'ENDOMMAGEMENT CONNUE */
/*  TN		: PLASTIC UNIT NORMAL TENSOR */
/*  DMI		: INVERSE DAMAGE EFFECT TENSOR */
/*  DN		: COMPOSANTES DU VECTEUR ORIENTATION */
/*  TDL		: TETA*DLAMBDA */
/*  YPS		: VALEUR DE LA FONCTION YPS = (-Y/S)**s */
/*  EYV		: FUNCTION Ye = exp[3/2*<Sh>/Seq] */
/*  UDHM		: VALUE 1/(1-Dh)**m */
/*  UDHM1	: VALUE 1/(1-Dh)**(m+1) */
/*  WX		: KINEMATIC DYNAMIC RECOVERY UNKNOWN: WX = |DMI:X| */
/*  DMWX		: NORME DU TENSEUR Mt:(Alpha1+Tdl*N) */
/*  PHTRWX	: (2/3)*PH/(1+TETA*(PRD*DLAMBDA+PRS)*WX) */
/*  WR		: ISOTROPIC DYNAMIC RECOVERY UNKNOWN: WR = R */
/*  RTDL		: PR1+TETA*DLAMDA */
/*  HTRWR	: RT(2/3)*H/(1+TETA*(RD*DLAMBDA+RS)*WR) */
/*  XDN		: CONSTANTE NUCLEATION XLD*XFN/(XKIC*XLF**(1/3)) */
/*  CTE		: CONSTANTE DE TENSION */
/*  CCO		: CONSTANTE DE COMPRESSION */
/*  CTO		: CONSTANTE DE TORSION */
/*  AV		: VOID GROWTH CONSTANT */
/*  VM		: VOID GROWTH INDICE M */
/*  VC		: VOID GROWTH CRITICAL THRESHOLD */
/*  AC		: COALESCENCE CONSTANT */
/*  FKIC		: FRACTURE TOUGHNESS CONSTANT */
/*  DSCALE	: LENGTH SCALE PARAMETER */
/*  DI		: VOID IMPINGEMENT CRITICAL DIAMETER FRACTION */
/*  DS		: VOID SHEET CRITICAL DIAMETER FRACTION */
/*  NCD		: NCD=1 --> COUPLAGE ELASTOPLASTICITE-ENDOMMAGEMENT */
/*  		  NCD=0 --> ELASTOPLASTICITE SANS ENDOMMAGEMENT */
/*  NTENS	: TAILLE DU TABLEAU DES CONTRAINTES OU DEFORMATIONS */
/*  NDI		: NOMBRE DE CONTRAINTES DIRECTES (S11,S22,S33) */
/*  NSHR		: NOMBRE DE CONTRAINTES TANGENTIELLES (S12,S13,S23) */
/*  M		: COMMENTAIRES */
/*  MPI		: UNITE DE FICHIER DE SORTIE */
/* ------------- */
/*  EN SORTIE : */
/* ------------- */
/*  FF		: FONCTION EQUATION CRITERE */
/*  FGX		: FONCTION NORME DES CONTRAINTES EFFECTIVES Wx = |Xt| */
/*  FGR		: FONCTION SINUS HYPERBOLIQUE DE Wr = Sinh(Qs*R) */
/*  FH		: FONCTION EQUATION TENSORIELLE D'ENDOMMAGEMENT */
/*  DENS		: DENSITE DE DISTRIBUTION (TENSEUR) */
/*  ENDOEQ	: NORME DE LA DENSITE DE DISTRIBUTION */
/*  QD		: VOID COALESCENCE OPERATOR */
/*  DQD		: COALESCENCE DERIVATIVE HEAVISIDE */
/*  VP		: VOID GROWTH TEST (VC = 1 if Dh > Vc) */
/*  VT		: POSITIVE PART <Dh-Vc>  (Dh=1/3*Trace(Phi)) */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
/* ----------------------------------------------------------------------- */
/* .1-----  Implicit, External */
/* .2-----  Parameter,Integer */
/* .3-----  Dimension */
/* ------- */
/* .4-----  Real,Complex,Double precision,Logical,Character */
/* .5-----  Common */
/* .6-----  Equivalence */
/* .7-----  Save,Data */
    /* Parameter adjustments */
    --h__;
    --alpha1;
    --zd;
    --tn;

    /* Function Body */
/* .8-----  Definition of Functions */
/* ---------------------------------------------------------------------* */
/* ====================================================================== */
/* ---------------------------------------------------------------------* */
/* -------  Initialization */
    d2g = two * *g;
    onethird = one / three;
/* -------  END Initialization */
/* ---------------------------------------------------------------------* */
/* ---------------------------------------------------------------------* */
/* -------  Yield Function */
    *ff = *zdzd - *zdn;
/* -------  END Yield Function */
/* ---------------------------------------------------------------------* */
/* -------  Kinematic Hardening Function */
    *fx = *wx - *phtrwx * *dmwx;
/* -------  END Kinematic Hardening Function */
/* ---------------------------------------------------------------------* */
/* -------  Isotropic Hardening Function */
    *fr = *wr - *htrwr * *rtdl;
/* -------  END Isotropic Hardening Function */
/* ---------------------------------------------------------------------* */
    if (*ncd == 1 || *ncd == 2) {
/* -------  Damage Function */
/* -------------------------------------------* */
/* ---------  Stress Tensors */
	dlgdhx = *tdl * (d2g * *unmd + *phtrwx);
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    altn[i__ - 1] = alpha1[i__] + *tdl * tn[i__];
	    xi[i__ - 1] = zd[i__] - dlgdhx * tn[i__];
	}
/* -------------------------------------------* */
/* ---------  Norm |Xi| = Sqrt(3/2 * Xi:Xi) */
	xieqv = xi[0] * xi[0] + xi[1] * xi[1] + xi[2] * xi[2];
	i__1 = *ntens;
	for (i__ = *ndi + 1; i__ <= i__1; ++i__) {
	    xieqv += two * xi[i__ - 1] * xi[i__ - 1];
	}
	xieqv = sqrt(thalf * xieqv);
/* ---------  Norm Seff = Sqrt(h1*Xieqv^2+h2*p^2) */
	seff = sqrt(h__[1] * xieqv * xieqv + h__[2] * *p * *p);
    }
/* -------------------------------------------* */
    if (*ncd == 1) {
	*eyv = onethird * *p / seff;
	xmm1 = *xms / (*xms + one);
	aphi = two * (*xms + one) * (one + *phi) * pow_dd(unmd, &xmm1);
	aphi = xnine * *xms * *phi * *unmd * *eyv / aphi;
	*fd = *phi - *phiold - *tdl * aphi;
    } else if (*ncd == 2) {
	*eyv = onethird * *p / xieqv;
	*eyv = (*eyv > 0.0) ? *eyv : 0.0; /* Macauley bracket */
	xms2 = two * *xms;
	unmdm = pow_dd(unmd, xms);
	amm = two * (xms2 - one) / (xms2 + one);
	aphi = amm * *eyv;
	sinha = sinh(aphi);
	cphi = one / unmdm - *unmd;
	*fd = *phi - *phiold - *tdl * cphi * sinha;
    } else if (*ncd == 0) {
	*fd = zero;
    }
/* -------  END Damage Function */
/* ---------------------------------------------------------------------* */
/* ---------------------------------------------------------------------* */
/* -------  Comment Check Code */
    if (*m >= 2) {
	io___372.ciunit = *mpi;
	s_wsfe(&io___372);
	do_fio(&c__1, "CALFG  FF = ", (ftnlen)12);
	do_fio(&c__1, (char *)&(*ff), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___373.ciunit = *mpi;
	s_wsfe(&io___373);
	do_fio(&c__1, "CALFG  FX = ", (ftnlen)12);
	do_fio(&c__1, (char *)&(*fx), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___374.ciunit = *mpi;
	s_wsfe(&io___374);
	do_fio(&c__1, "CALFG  FR = ", (ftnlen)12);
	do_fio(&c__1, (char *)&(*fr), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___375.ciunit = *mpi;
	s_wsfe(&io___375);
	do_fio(&c__1, "CALFG  FD = ", (ftnlen)12);
	do_fio(&c__1, (char *)&(*fd), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
/* -------  END Comment Check Code */
/* ----------------------------------------------------------------------* */
/* ---------   FORMATS  -------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* L101: */
/* ----------------------------------------------------------------------* */
/* ======================================================================= */
    return 0;
} /* calfgh_ */


/* Subroutine */ int calfghi_(doublereal *ff, doublereal *fx, doublereal *fr, 
	doublereal *fd, doublereal *zdn, doublereal *zdzd, doublereal *phi, 
	doublereal *phiold, doublereal *unmd, doublereal *wx, doublereal *
	dmwx, doublereal *phtrwx, doublereal *wr, doublereal *rtdl, 
	doublereal *htrwr, doublereal *xms, doublereal *xmm1, integer *ncd, 
	integer *m, integer *mpi)
{
    /* Initialized data */

    static doublereal zero = 0.;

    /* Format strings */
    static char fmt_102[] = "(1x,a20,6(d20.14,2x))";

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Fortran I/O blocks */
    static cilist io___377 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___378 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___379 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___380 = { 0, 0, 0, fmt_102, 0 };


/* ====================================================================== */
/* ---------------------------------------------------------------------- */
/* -------------  CALCUL DES FONCTIONS Fi, GXi, GRi ET Hi  -------------- */
/* ---------------------------------------------------------------------- */
/* ====================================================================== */
/* ------------- */
/*  EN ENTREE : */
/* ------------- */
/*  ENDO 	: TENSEUR D'ENDOMMAGEMENT */
/*  ZDN		: NORME DU TENSEUR DE DIRECTION NORMALE ZD */
/*  ZDZD		: NORME DU TENSEUR DE DIRECTION NORMALE ZD */
/*  ENDO1	: TENSEUR D'ENDOMMAGEMENT CONNU */
/*  WX		: KINEMATIC DYNAMIC RECOVERY TERM: WX = |DMI:X| */
/*  DMWX		: NORME DU TENSEUR Mt:(ALPHA1+TDL*TN) */
/*  PHTRWX	: (2/3)*PH/(1+TETA*(PRD*DLAMBDA+PRS)*WX) */
/*  WR		: ISOTROPIC DYNAMIC RECOVERY TERM: WR = R */
/*  RTDL		: PR1+TETA*DLAMDA */
/*  HTRWR	: RT(2/3)*H/(1+TETA*(RD*DLAMBDA+RS)*WR) */
/*  VC		: VOID GROWTH CRITICAL THRESHOLD */
/*  FKIC		: FRACTURE TOUGHNESS CONSTANT */
/*  DSCALE	: LENGTH SCALE PARAMETER */
/*  DI		: VOID IMPINGEMENT CRITICAL DIAMETER FRACTION */
/*  DS		: VOID SHEET CRITICAL DIAMETER FRACTION */
/*  NCD		: NCD=1 --> COUPLAGE ELASTOPLASTICITE-ENDOMMAGEMENT */
/*  		  NCD=0 --> ELASTOPLASTICITE SANS ENDOMMAGEMENT */
/*  M		: COMMENTAIRES */
/*  MPI		: UNITE DE FICHIER DE SORTIE */
/* ------------- */
/*  EN SORTIE : */
/* ------------- */
/*  FF		: YIELD FUNCTION */
/*  FX		: KINMETIC */
/*  FR		: ISOTROPIC HARDENING STRESS Wr = R */
/*  FD		: DAMAGE FUNCTION FH */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
/* ----------------------------------------------------------------------- */
/* .1-----  Implicit, External */
/* .2-----  Parameter,Integer */
/* .3-----  Dimension */
/* .4-----  Real,Complex,Double precision,Logical,Character */
/* .5-----  Common */
/* .6-----  Equivalence */
/* .7-----  Save,Data */
/* .8-----  Definition de fonctions */
/* ---------------------------------------------------------------------* */
/* ====================================================================== */
/* ---------------------------------------------------------------------* */
/* -------  Yield Function Ff */

    *ff = *zdzd - *zdn;

/* -------  END Yield Function Ff */
/* ---------------------------------------------------------------------* */
/* -------  Kinematic Hardening Function Fx */

    *fx = *wx - *phtrwx * *dmwx;

/* -------  END Kinematic Hardening Function Fx */
/* ---------------------------------------------------------------------* */
/* -------  Isotropic Hardening Function Fr */

    *fr = *wr - *htrwr * *rtdl;

/* -------  END Isotropic Hardening Function Fr */
/* ---------------------------------------------------------------------* */
/* -------  Damage Functions Fh */

    if (*ncd == 0) {
	*fd = zero;
    } else {
	*fd = *phi - *phiold;
    }

/* -------  END Damage Functions Fh */
/* ---------------------------------------------------------------------* */
/* -------  Comment Check Code */
    if (*m >= 4) {
	io___377.ciunit = *mpi;
	s_wsfe(&io___377);
	do_fio(&c__1, "CALFGI  FF = ", (ftnlen)13);
	do_fio(&c__1, (char *)&(*ff), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___378.ciunit = *mpi;
	s_wsfe(&io___378);
	do_fio(&c__1, "CALFGI  FX = ", (ftnlen)13);
	do_fio(&c__1, (char *)&(*fx), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___379.ciunit = *mpi;
	s_wsfe(&io___379);
	do_fio(&c__1, "CALFGI  FR = ", (ftnlen)13);
	do_fio(&c__1, (char *)&(*fr), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___380.ciunit = *mpi;
	s_wsfe(&io___380);
	do_fio(&c__1, "CALFGI  FH = ", (ftnlen)13);
	do_fio(&c__1, (char *)&(*fd), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
/* -------  END Comment Check Code */
/* ----------------------------------------------------------------------* */
/* ---------   FORMATS  -------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* L101: */
/* ----------------------------------------------------------------------* */
/* ======================================================================= */
    return 0;
} /* calfghi_ */


/* Subroutine */ int calzdtn_(doublereal *zd, doublereal *tn, doublereal *zdn,
	 doublereal *zdzd, doublereal *dsqhmp, doublereal *ryvshdl, 
	doublereal *stres1d, doublereal *estar, doublereal *alpha1, 
	doublereal *p, doublereal *pr1, doublereal *phi, doublereal *dphi, 
	doublereal *pold, doublereal *xkdh, doublereal *h__, doublereal *tdl, 
	doublereal *dtlf, doublereal *unmd, doublereal *dmwx, doublereal *
	phtrwx, doublereal *rtdl, doublereal *htrwr, doublereal *yield, 
	doublereal *v, doublereal *pf, doublereal *g, doublereal *xk, 
	doublereal *xms, integer *ncd, integer *ntens, integer *ndi, integer *
	nshr, integer *m, integer *mpi)
{
    /* Initialized data */

    static doublereal zero = 0.;
    static doublereal half = .5;
    static doublereal one = 1.;
    static doublereal two = 2.;
    static doublereal three = 3.;

    /* Format strings */
    static char fmt_102[] = "(1x,a20,6(d20.14,2x))";
    static char fmt_101[] = "(1x,6(d20.14,2x))";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal), pow_dd(doublereal *, doublereal *), log(
	    doublereal);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void),
	     s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    static doublereal twothird;
    static integer i__;
    static doublereal d2g, dasinhdtlf, hmp, xmm1, xms1, altn[6], upsd, sqtt;
    extern /* Subroutine */ int pdtsca_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal gestar[6];
    extern /* Subroutine */ int tensdev_(doublereal *, doublereal *, integer *
	    , integer *);

    /* Fortran I/O blocks */
    static cilist io___397 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___398 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___399 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___400 = { 0, 0, 0, 0, 0 };
    static cilist io___401 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___402 = { 0, 0, 0, 0, 0 };
    static cilist io___403 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___404 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___405 = { 0, 0, 0, fmt_101, 0 };


/* ====================================================================== */
/* ---------------------------------------------------------------------- */
/* --------------------  CALCUL DU TENSEUR NORMAL  ---------------------- */
/* ---------------------------------------------------------------------- */
/* ====================================================================== */
/* ------------- */
/*  EN ENTREE : */
/* ------------- */
/*  ESTAR	: TENSEUR DES DEFORMATIONS CONNUS A L'INSTANT T+DT */
/*  ALPHA1	: TENSEUR DES CONTRAINTES CINEMATIQUES CONNU A T */
/*  PR1		: VARIABLE D'ECROUISSAGE ISOTROPE CONNU A T */
/*  XKDEH	: Xk*Tr(d) */
/*  YIELD	: LIMITE D'ELASTICITE */
/*  TDL		: TETA*DLAMBDA */
/*  DTLF		: RT(2/3)*DLAMBDA/PF */
/*  PHTRWX	: (2/3)*PH/(1+TETA*(PRD*DLAMBDA+RT(2/3)*PRS)*WX) */
/*  RTDL		: PR1+TETA*DLAMDA */
/*  HTRWR	: RT(2/3)*H/(1+(TETA*RD*DLAMBDA+RS)*WR) */
/*  YIELD	: LIMITE D'ELASTICITE */
/*  V		: CONSTANTE VISCOPLASTIQUE */
/*  PF		: CONSTANTE VISCOPLASTIQUE */
/*  G		: MODULE DE CISAILLEMENT */
/*  PH		: CONSTANTE D'ECROUISSAGE CINEMATIQUE */
/*  NDI		: NOMBRE DE COMPOSANTES DES CONTRAINTES DIRECTES */
/*  NSHR		: NOMBRE DE COMPOSANTES DES CONTRAINTES TANGENTIELLES */
/*  NTENS	: TAILLE DU TABLEAU DES CONTRAINTES OU DEFORMATIONS */
/* 		  (NDI + NSHR) */
/*  M		: COMMENTAIRES */
/*  MPI		: UNITE DE FICHIER DE SORTIE */
/* ------------- */
/*  EN SORTIE : */
/* ------------- */
/*  ZD		: TENSEUR DE DIRECTION NORMALE */
/*  TN		: TENSEUR NORMAL */
/*  ZDN		: NORME DU TENSEUR ZD */
/*  ZDZD		: NORME DU TENSEUR ZD EGALE A RACINE(2/3*ZD:ZD) */
/*  P		: HYDROSTATIC STRESS */
/*  POLD		: HYDROSTATIC STRESS AT TIME T */
/*  DMWX		: NORME DU TENSEUR Mt:(Alpha1+Tdl*N) */
/*  STRESS1D	: TENSEUR DEVIATORIQUE DES CONTRAINTES CONNU A T */
/*  H		: DAMAGE FUNCTIONS h1, h2 and h3 */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
/* .1-----  Implicit, External */
/* .2-----  Parameter,Integer */
/* .3-----  Dimension */
/* ------- */
/* .4-----  Real,Complex,Double precision,Logical,Character */
/* .5-----  Common */
/* .6-----  Equivalence */
/* .7-----  Save,Data */
    /* Parameter adjustments */
    --h__;
    --alpha1;
    --estar;
    --stres1d;
    --tn;
    --zd;

    /* Function Body */
/* .8-----  Definition de fonctions */
/* ---------------------------------------------------------------------* */
/* ====================================================================== */
/* ---------------------------------------------------------------------* */
/* -------  Initialisation */
    twothird = two / three;
    sqtt = sqrt(twothird);
/* -------  FIN  Initialisation */
/* ---------------------------------------------------------------------* */
/* -------  Hydrostatic Stress */
/*      P = POLD+UNMD*XKDEH+XK*DPHI */
    *p = *unmd * *xkdh;
/* -------  END  Hydrostatic Stress */
/* ---------------------------------------------------------------------* */
/* -------  Sigma1D */
    d2g = two * *g;
    tensdev_(&estar[1], gestar, ndi, nshr);
/* ------- */
    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
	stres1d[i__] = d2g * gestar[i__ - 1];
    }
/* -------  END  Sigma1D */
/* ---------------------------------------------------------------------* */
/* -------  Tensor ZD */
    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
	zd[i__] = *unmd * stres1d[i__] - *phtrwx * alpha1[i__];
    }
/* -------  END  Tensor ZD */
/* ---------------------------------------------------------------------* */
/* -------  Marin and McDowell Damage Functions H */
    xms1 = *xms + one;
    xmm1 = one / xms1;
    upsd = one + *phi;
    if (*ncd == 1) {
	h__[1] = one + twothird * *phi;
	h__[2] = half * xmm1 * *phi / upsd;
	h__[3] = pow_dd(unmd, &xmm1);
    } else if (*ncd == 2) {
	h__[1] = one;
	h__[2] = half * *phi / upsd;
	h__[3] = *unmd;
    } else {
	h__[1] = one;
	h__[2] = zero;
	h__[3] = one;
    }
    h__[1] = one;
    h__[2] = zero;
    h__[3] = one;
/* -------------------------------------------* */
/* -------  First Norm of Tensor ZD with the Sine Function */
    dasinhdtlf = log(*dtlf + sqrt(one + *dtlf * *dtlf));
    *ryvshdl = *htrwr * *rtdl + *yield + *v * dasinhdtlf;
    hmp = (h__[3] * h__[3] * *ryvshdl * *ryvshdl - h__[2] * *p * *p) / h__[1];
    *dsqhmp = sqrt(hmp);
    *zdn = twothird * *unmd * *dsqhmp + *tdl * (*unmd * d2g + *phtrwx);
/* -------  Second Norm of Tensor ZD with |Zd*Zd| */
    pdtsca_(ndi, nshr, &zd[1], &zd[1], zdzd);
    *zdzd = sqrt(twothird * *zdzd);
/* -------  END  Calcul de la norme du tenseur ZD suivant deux methodes */
/* ---------------------------------------------------------------------* */
/* -------  Tensor TN */
    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
	tn[i__] = zd[i__] / *zdn;
    }
/* -------  END  Tensor TN */
/* ---------------------------------------------------------------------* */
/* -------  Norm DMWX = |Mt:(Alpha1+Tdl*N)| */
    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
	altn[i__ - 1] = alpha1[i__] + *tdl * tn[i__];
    }
    pdtsca_(ndi, nshr, altn, altn, dmwx);
    *dmwx = sqrt(*dmwx);
/* -------  END  Norm DMWX = |Mt:(Alpha1+Tdl*N)| */
/* ---------------------------------------------------------------------* */
/* -------  Comment Check Code */

    if (*m >= 4) {
	io___397.ciunit = *mpi;
	s_wsfe(&io___397);
	do_fio(&c__1, "CALZDTN ZDZD = ", (ftnlen)15);
	do_fio(&c__1, (char *)&(*zdzd), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___398.ciunit = *mpi;
	s_wsfe(&io___398);
	do_fio(&c__1, "CALZDTN ZDN  = ", (ftnlen)15);
	do_fio(&c__1, (char *)&(*zdn), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___399.ciunit = *mpi;
	s_wsfe(&io___399);
	do_fio(&c__1, "CALZDTN DMWX = ", (ftnlen)15);
	do_fio(&c__1, (char *)&(*dmwx), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___400.ciunit = *mpi;
	s_wsle(&io___400);
	do_lio(&c__9, &c__1, "CALZDTN ZD", (ftnlen)10);
	e_wsle();
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___401.ciunit = *mpi;
	    s_wsfe(&io___401);
	    do_fio(&c__1, (char *)&zd[i__], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	io___402.ciunit = *mpi;
	s_wsle(&io___402);
	do_lio(&c__9, &c__1, "CALZDTN TN", (ftnlen)10);
	e_wsle();
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___403.ciunit = *mpi;
	    s_wsfe(&io___403);
	    do_fio(&c__1, (char *)&tn[i__], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	io___404.ciunit = *mpi;
	s_wsfe(&io___404);
	do_fio(&c__1, "CALZDTN H", (ftnlen)9);
	e_wsfe();
	for (i__ = 1; i__ <= 3; ++i__) {
	    io___405.ciunit = *mpi;
	    s_wsfe(&io___405);
	    do_fio(&c__1, (char *)&h__[i__], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }

/* -------  FIN Comment Check Code */
/* ----------------------------------------------------------------------* */
/* ---------   FORMATS  -------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* ======================================================================= */
    return 0;
} /* calzdtn_ */


/* Subroutine */ int dpolr_(doublereal *rot, doublereal *f, integer *m, 
	integer *mpi)
{
    /* Initialized data */

    static doublereal un = 1.;

    /* Format strings */
    static char fmt_101[] = "(1x,3(d18.12,2x))";

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void), s_wsfe(cilist *), do_fio(integer *, char *, ftnlen),
	     e_wsfe(void);

    /* Local variables */
    static doublereal d__[3];
    static integer i__, j;
    static doublereal v[9]	/* was [3][3] */;
    extern /* Subroutine */ int transpose_(integer *, doublereal *, 
	    doublereal *);
    static doublereal tf[9]	/* was [3][3] */, vn[9]	/* was [3][3] */, rd1,
	     rd2, rd3, dv1[9]	/* was [3][3] */, ftf[9]	/* was [3][3] 
	    */, rtr[9]	/* was [3][3] */, trot[9]	/* was [3][3] */;
    extern /* Subroutine */ int prdmat_(integer *, doublereal *, doublereal *,
	     doublereal *), valprop_(doublereal *, doublereal *, doublereal *,
	     integer *, integer *);

    /* Fortran I/O blocks */
    static cilist io___422 = { 0, 0, 0, 0, 0 };
    static cilist io___423 = { 0, 0, 0, fmt_101, 0 };


/* ====================================================================== */
/* ---------------------------------------------------------------------- */
/* ---------------------  POLAR DECOMPOSITION OF F  --------------------- */
/* ---------------------------------------------------------------------- */
/* ====================================================================== */
/* ----------- */
/*  ENTREES : */
/* ----------- */
/*  F		: DEFORMATION GRADIENT */
/*  M		: PARAMETRE D'IMPRESSION */
/*  MPI		: UNITE DE FICHIER D'IMPRESSION */
/* ----------- */
/*  SORTIES : */
/* ----------- */
/*  ROT		: MATRICE ROTATION PROPRE */
/* ----------------------------------------------------------------------* */
/* ======================================================================= */
/* ----------------------------------------------------------------------* */
/* .1-----  Precision */
/* .2-----  Parameter */
/* .3-----  Dimension */
/* .4-----  Real,Integer,Complex,Double precision,Logical,Character */
/* .5-----  Common */
/* .6-----  Equivalence */
/* .7-----  Save,Data */
    /* Parameter adjustments */
    f -= 4;
    rot -= 4;

    /* Function Body */
/* .8-----  Definition de fonctions */
/* ----------------------------------------------------------------------* */
/* ======================================================================= */
/* ----------------------------------------------------------------------* */
/* -------  Transposition de la matrice F et calcul de TF*F */

    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    tf[i__ + j * 3 - 4] = f[j + i__ * 3];
	}
    }
    prdmat_(&c__3, tf, &f[4], ftf);

/* -------  FIN  Transposition des matrices F0 et F1 */
/* ----------------------------------------------------------------------* */
/* -------  Calcul des valeurs et vecteurs propres de FTF */

    valprop_(ftf, d__, v, &c__3, &c__3);
/* -------  FIN  Calcul des valeurs et vecteurs propres de FTF */
/* ----------------------------------------------------------------------* */
/* -------  Calcul de l'inverse du tenseur de deformations pure gauche DV */


/* -------  Calcul de l'inverse des racines des valeurs propres */

    rd1 = un / sqrt(d__[0]);
    rd2 = un / sqrt(d__[1]);
    rd3 = un / sqrt(d__[2]);

/* -------  Calcul des vecteurs propres nI, nII et nIII */

    vn[0] = rd1 * (f[4] * v[0] + f[7] * v[1] + f[10] * v[2]);
    vn[1] = rd1 * (f[5] * v[0] + f[8] * v[1] + f[11] * v[2]);
    vn[2] = rd1 * (f[6] * v[0] + f[9] * v[1] + f[12] * v[2]);
/* ------- */
    vn[3] = rd2 * (f[4] * v[3] + f[7] * v[4] + f[10] * v[5]);
    vn[4] = rd2 * (f[5] * v[3] + f[8] * v[4] + f[11] * v[5]);
    vn[5] = rd2 * (f[6] * v[3] + f[9] * v[4] + f[12] * v[5]);
/* ------- */
    vn[6] = rd3 * (f[4] * v[6] + f[7] * v[7] + f[10] * v[8]);
    vn[7] = rd3 * (f[5] * v[6] + f[8] * v[7] + f[11] * v[8]);
    vn[8] = rd3 * (f[6] * v[6] + f[9] * v[7] + f[12] * v[8]);

/* -------  Calcul de DV1 inverse de DV */

    dv1[0] = rd1 * vn[0] * vn[0] + rd2 * vn[3] * vn[3] + rd3 * vn[6] * vn[6];
    dv1[3] = rd1 * vn[0] * vn[1] + rd2 * vn[3] * vn[4] + rd3 * vn[6] * vn[7];
    dv1[6] = rd1 * vn[0] * vn[2] + rd2 * vn[3] * vn[5] + rd3 * vn[6] * vn[8];
    dv1[4] = rd1 * vn[1] * vn[1] + rd2 * vn[4] * vn[4] + rd3 * vn[7] * vn[7];
    dv1[7] = rd1 * vn[1] * vn[2] + rd2 * vn[4] * vn[5] + rd3 * vn[7] * vn[8];
    dv1[8] = rd1 * vn[2] * vn[2] + rd2 * vn[5] * vn[5] + rd3 * vn[8] * vn[8];
    dv1[1] = dv1[3];
    dv1[2] = dv1[6];
    dv1[5] = dv1[7];

/* -------  FIN  Calcul de l'inverse du tenseur de deformations pure DV */
/* ----------------------------------------------------------------------* */
/* -------  Calcul de ROT */

    prdmat_(&c__3, dv1, &f[4], &rot[4]);
    transpose_(&c__3, &rot[4], trot);
    prdmat_(&c__3, &rot[4], trot, rtr);

/* -------  FIN  Calcul de DR */
/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* -------  Comment Check Code */

    if (*m >= 4) {
	io___422.ciunit = *mpi;
	s_wsle(&io___422);
	do_lio(&c__9, &c__1, "DPOLR ROT", (ftnlen)9);
	e_wsle();
	for (i__ = 1; i__ <= 3; ++i__) {
	    io___423.ciunit = *mpi;
	    s_wsfe(&io___423);
	    for (j = 1; j <= 3; ++j) {
		do_fio(&c__1, (char *)&rot[i__ + j * 3], (ftnlen)sizeof(
			doublereal));
	    }
	    e_wsfe();
	}
    }

/* -------  FIN Comment Check Code */
/* ----------------------------------------------------------------------* */
/* ======================================================================= */
/* ----------------------------------------------------------------------* */
/* ---------   FORMATS  -------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
    return 0;
} /* dpolr_ */


/* Subroutine */ int elplastdo_(doublereal *ddsdde, doublereal *stress, 
	doublereal *stran, doublereal *dstran, doublereal *eplas, doublereal *
	alpha, doublereal *xlambda, doublereal *dlambda, doublereal *pr, 
	doublereal *phi, doublereal *dlan, doublereal *stvit, doublereal *
	yield, doublereal *v, doublereal *pf, doublereal *hk, doublereal *rdk,
	 doublereal *rsk, doublereal *hi, doublereal *rdi, doublereal *rsi, 
	doublereal *dc, doublereal *xms, doublereal *g, doublereal *xk, 
	doublereal *teta, doublereal *epsf, doublereal *epsx, doublereal *
	epsr, doublereal *epsd, doublereal *dtime, logical *lonewdt, integer *
	ncd, integer *nplan, integer *natur, integer *ntens, integer *ndi, 
	integer *nshr, integer *nstvi, integer *m, integer *mpi)
{
    /* Initialized data */

    static doublereal zero = 0.;
    static doublereal un = 1.;
    static integer melpst = 1;

    /* Format strings */
    static char fmt_300[] = "(1x,\002SUB ELPLASTDO: DAMAGE COMPONENT GREATER"
	    " THAN 1\002)";
    static char fmt_400[] = "(1x,\002*--------------------------------------"
	    "--------*\002/1x,\002| *----*  REACHED ITERATION NUMBER: \002,i3,"
	    "\002  *----* |\002,/1x,\002*------------------------------------"
	    "----------*\002)";
    static char fmt_200[] = "(1x,\002SUB ELPLASTDO: MAXIMUM ITERATION\002,"
	    "1x,\002NUMBER HAS BEEN REACHED\002)";

    /* System generated locals */
    integer ddsdde_dim1, ddsdde_offset, i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), s_wsle(cilist *), do_lio(integer *
	    , integer *, char *, ftnlen), e_wsle(void), do_fio(integer *, 
	    char *, ftnlen);

    /* Local variables */
    extern /* Subroutine */ int calcoeff_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *), solution_(doublereal 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *);
    static doublereal h__[3], p;
    extern /* Subroutine */ int calcoeffi_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *);
    static doublereal fd, ff, fr, zd[6], fx, tn[6], wr, wx;
    extern /* Subroutine */ int verifgence_(logical *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *);
    static doublereal pr1, tdl, zdn, eyv;
    static logical convergence;
    static doublereal phi1, xmm1, dfdd, dhdd, dndd[6], dfdl, dhdl, dpdd, dfdr,
	     dedt[6], dndl[6], dfds, dxdd, dhdr, dfdx, drdl, dndr[6], dnds[6],
	     dhdx, dhds, dxdl, altn[6], dndx[6], pold, dxdr, unmd, dxds, rtdl,
	     drdr, dxdx, drds, tetn[6], zdzd, xkdh;
    extern /* Subroutine */ int aset_(doublereal *, doublereal *, integer *);
    static doublereal dmwx, unmd1, dfrdl, estar[6], alpha1[6], htrwr, trswr, 
	    eplas1[6];
    extern /* Subroutine */ int calfgh_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, integer *, integer *, integer *);
    static doublereal htrwr2;
    extern /* Subroutine */ int caldpf_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *), caldpg_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, integer *, integer *, integer *), caldph_(
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, integer *, integer *, integer *), caldpn_(
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *)
	    ;
    static doublereal dzdndd, dzdndl, alphzd, gestar[6];
    extern /* Subroutine */ int ktimen_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, integer *, integer *, integer *);
    static doublereal dsqhmp;
    static integer nitmax;
    static doublereal strant[6], phtrwx, untrwr, untrwx, stres1d[6];
    extern /* Subroutine */ int calfghi_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *), jddsdde_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *);
    static doublereal phtrwx2, ryvshdl;
    extern /* Subroutine */ int solsyst_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, integer *, integer *, integer *);

    /* Fortran I/O blocks */
    static cilist io___496 = { 0, 0, 0, fmt_300, 0 };
    static cilist io___497 = { 0, 0, 0, 0, 0 };
    static cilist io___498 = { 0, 0, 0, 0, 0 };
    static cilist io___500 = { 0, 0, 0, fmt_400, 0 };
    static cilist io___501 = { 0, 0, 0, fmt_200, 0 };


/* ====================================================================== */
/* ---------------------------------------------------------------------- */
/* ---------------  SOLUTION PLASTIQUE AVEC RETOUR RADIAL  -------------- */
/* ---------------------------------------------------------------------- */
/* ====================================================================== */
/* ----------- */
/*  ENTREES : */
/* ----------- */
/*  STRAN	: TENSEUR DES DEFORMATIONS TOTALES (ESTRAN) */
/*  DSTRAN	: INCREMENT DE DEFORMATION TOTALE */
/*  UNITDEV	: TENSEUR UNITE DEVIATORIQUE Idev=I-1/3(1*1) */
/*  E		: MODULE D'YOUNG */
/*  XNU		: COEFFICIENT DE POISSON */
/*  G		: MODULE DE CISAILLEMENT */
/*  XK		: COEFFICIENT DE LAME */
/*  YIELD	: LIMITE D'ELASTICITE */
/*  V		: CONSTANTE VISCOPLASTIQUE */
/*  PF		: CONSTANTE VISCOPLASTIQUE */
/*  RD		: CONSTANTE D'ECROUISSAGE ISOTROPE */
/*  H		: CONSTANTE D'ECROUISSAGE ISOTROPE */
/*  RS		: CONSTANTE D'ECROUISSAGE ISOTROPE (RECOVERY) */
/*  PRD		: CONSTANTE D'ECROUISSAGE CINEMATIQUE */
/*  PH		: CONSTANTE D'ECROUISSAGE CINEMATIQUE */
/*  PRS		: CONSTANTE D'ECROUISSAGE CINEMATIQUE (RECOVERY) */
/*  XIP		: INITIAL POROSITY */
/*  DC		: VALEUR DE L'ENDOMMAGEMENT CRITIQUE --> RUPTURE */
/*  XIP		: INITIAL POROSITY */
/*  XLF		: VOLUME FRACTION */
/*  XDN		: CONSTANTE NUCLEATION XLD*XFN/(XKIC*XLF**(1/3)) */
/*  GS		: FONCTION D'ENDOMMAGEMENT DUCTILE */
/*  PS		: FONCTION D'ENDOMMAGEMENT DUCTILE */
/*  DC		: VALEUR D'ENDOMAGEMMENT CRITIQUE */
/*  CTE		: CONSTANTE DE TENSION */
/*  CCO		: CONSTANTE DE COMPRESSION */
/*  CTO		: CONSTANTE DE TORSION */
/*  AV		: VOID GROWTH CONSTANT */
/*  VM		: VOID GROWTH INDICE M */
/*  VC		: VOID GROWTH CRITICAL THRESHOLD */
/*  AC		: COALESCENCE CONSTANT */
/*  DSCALE	: LENGTH SCALE PARAMETER */
/*  DI		: VOID IMPINGEMENT CRITICAL DIAMETER RATIO */
/*  DS		: VOID SHEET CRITICAL DIAMETER RATIO */
/*  TETA		: TETA-METHODE */
/*  EPSD		: PRECISION DE LA FONCTION DE CHARGE H */
/*  EPSG		: PRECISION DE LA FONCTION D'ENDOMMAGEMENT G */
/*  KCT		: KCT=0 --> MATRICE TANGENTE CONTINUE */
/*  		  KCT=1 --> MATRICE TANGENTE CONSISTENTE */
/*  NPLAN	: NPLAN=1 --> DEFORMATIONS PLANES */
/* 		  NPLAN=2 --> CONTRAINTES PLANES */
/*  NATUR	: NATUR=0 --> PETITES DEFORMATIONS */
/* 		  NATUR=1 --> GRANDES DEFORMATIONS */
/*  IROT		: IROT=0 --> PAS DE ROTATION */
/* 		  IROT=1 --> AVEC ROTATION */
/* 		: IROT=2 --> AVEC DERIVEE COROTATIONNELLE DE JAUMANN */
/* 		  IROT=3 --> AVEC ROTATION D'ABAQUS (ROTSIG) */
/*  NCD		: NCD=1 --> COUPLAGE ELASTOPLASTICITE-ENDOMMAGEMENT */
/*  		  NCD=0 --> ELASTOPLASTICITE SANS ENDOMMAGEMENT */
/*  NDI		: NOMBRE DE COMPOSANTES DES CONTRAINTES DIRECTES */
/*  NSHR		: NOMBRE DE COMPOSANTES DES CONTRAINTES TANGENTIELLES */
/*  NTENS	: TAILLE DU TABLEAU DES CONTRAINTES OU DEFORMATIONS */
/* 		  (NDI + NSHR). */
/*  NSTVI	: TAILLE DU TABLEAU STVIT DE VARIABLES INTERNES */
/*  M		: COMMENTAIRES (PARAMETRE D'IMPRESSION) */
/*  MPI		: UNITE DE FICHIER DE SORTIE */
/* ----------- */
/*  LOCALES : */
/* ----------- */
/*  ZD		: TENSEUR DE DIRECTION NORMALE */
/*  TN		: TENSEUR NORMAL */
/*  STRANT	: TENSEUR DES DEFORMATIONS TOTALES : STRAN+DSTRAN */
/*  DEDT		: TENSEUR INCREMENT DE DEFORMATION : DTIME*DSTRAN */
/*  STRES1D	: TENSEUR DEVIATORIQUE DES CONTRAINTES CONNU A T */
/*  EPLAS1	: TENSEUR DES DEFORMATIONS PLASTIQUES CONNU A T */
/*  ALPHA1	: TENSEUR DES CONTRAINTES CINEMATIQUES CONNU A T */
/*  PR1		: VARIABLE D'ECROUISSAGE ISOTROPE CONNU A T */
/* --------------------------------------- */
/*  VARIABLES DEVANT ETRE REACTUALISEES : */
/* --------------------------------------- */
/*  STRESS	: TENSEUR DES CONTRAINTES */
/*  EPLAS	: TENSEUR DES DEFORMATIONS PLASTIQUES */
/*  ALPHA	: TENSEUR DES CONTRAINTES CINEMATIQUES */
/*  XLAMBDA	: DEFORMATION PLASTIQUE CUMULEE */
/*  PR		: VARIABLE D'ECROUISSAGE ISOTROPE */
/*  ENDO		: TENSEUR D'ENDOMMAGEMENT */
/*  DM		: TENSEUR D'EFFET D'ENDOMMAGEMENT (Mt) */
/*  DMI		: TENSEUR D'EFFET D'ENDOMMAGEMENT INVERSE (M(-t)) */
/* ----------- */
/*  SORTIES : */
/* ----------- */
/*  DDSDDE	: MATRICE TANGENTE CONSISTENTE (DSTRAN/DSTRESS) */
/*  STVIT	: TABLEAU DES VARIABLE INTERNES (R,X,J2(X),Y,...) */
/* ----------------------------------------------------------------------* */
/* ======================================================================= */
/* ----------------------------------------------------------------------* */
/* .1-----  Precision */
/* .2-----  Parameter,Integer */
/* .3-----  Dimension */
/* ------- */
/* .4-----  Real,Complex,Double precision,Logical,Character */
/* .5-----  Common */
/* .6-----  Equivalence */
/* .7-----  Save,Data */
    /* Parameter adjustments */
    --alpha;
    --eplas;
    --dstran;
    --stran;
    --stress;
    ddsdde_dim1 = *ntens;
    ddsdde_offset = 1 + ddsdde_dim1;
    ddsdde -= ddsdde_offset;
    --stvit;

    /* Function Body */
/* .8-----  Definition de fonctions */
/* ----------------------------------------------------------------------* */
/* ======================================================================= */
/* ----------------------------------------------------------------------* */
/* --------------------------   INITIALISATION  -------------------------* */
/* ----------------------------------------------------------------------* */
/* -------  LoNewdt --> False si D < 1 and LoNewdt --> True si D > 1 */
    *lonewdt = FALSE_;
/* --------- END LoNewdt --> False si D < 1 and LoNewdt --> True si D > 1 */
/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* -------  Calcul des constantes E*, e*, alpha1, r1, d1. */

    ktimen_(&stress[1], &stran[1], &dstran[1], &eplas[1], &alpha[1], &p, pr, 
	    phi, &xkdh, dlambda, dlan, strant, dedt, estar, gestar, stres1d, 
	    eplas1, alpha1, &pr1, &phi1, &pold, &wr, &wx, xk, g, hk, rdk, rsk,
	     hi, rdi, rsi, xms, teta, dtime, ncd, ntens, ndi, nshr, m, mpi);
/* --------- END Calcul des  constantes E*, e*, alpha1, r1, d1. */
/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* ---------  Initial non constant repetitive coefficients */

    calcoeffi_(zd, tn, &zdn, &zdzd, phi, stres1d, estar, &p, &xkdh, &dsqhmp, &
	    ryvshdl, alpha1, &pold, &pr1, &phi1, &unmd, &unmd1, h__, &wr, &wx,
	     v, pf, g, xk, hk, rdk, rsk, hi, rdi, rsi, xms, dlambda, &tdl, &
	    dmwx, yield, &untrwx, &phtrwx, &phtrwx2, &untrwr, &htrwr, &htrwr2,
	     &trswr, &rtdl, teta, &nitmax, ncd, ntens, ndi, nshr, m, mpi);
/* ----------- END  Initial non constant repetitive coefficients */
/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* -------  Initial functions Fi, GXi, GRi and Hi */

    calfghi_(&ff, &fx, &fr, &fd, &zdn, &zdzd, phi, &phi1, &unmd, &wx, &dmwx, &
	    phtrwx, &wr, &rtdl, &htrwr, xms, &xmm1, ncd, m, mpi);

/* --------- END Initial functions Fi, GXi, GRi and Hi */
/* ----------------------------------------------------------------------* */
/* ------------------------   END  INITIALISATION  ----------------------* */
/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------* */

/* -------  Solution Z, DLAMBDA and ENDO */

    convergence = FALSE_;
    while(! convergence && nitmax < 100) {

/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* -------  Partial derivatives of df/dlambda, df/dWx, df/dWr and df/dD */

	caldpf_(&dfdl, &dfdx, &dfdr, &dfdd, &dfds, &drds, &dpdd, &dzdndl, &
		dzdndd, &alphzd, &dfrdl, zd, h__, &p, phi, alpha1, stres1d, &
		pr1, &zdn, &zdzd, &dsqhmp, &phtrwx, &ryvshdl, dlambda, &tdl, &
		unmd, &unmd1, &xkdh, &wx, &wr, &phtrwx2, &htrwr, &htrwr2, &
		trswr, &rtdl, rdk, rsk, rdi, rsi, yield, g, xk, v, pf, xms, &
		xmm1, teta, ncd, &melpst, ntens, ndi, nshr, m, mpi);

/* ----------- END  Partial derivatives of de f */
/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* -------  Partial derivatives of dN/dlambda, dN/dWx and dN/dD */

	caldpn_(dndl, dndx, dndr, dndd, dnds, zd, alpha1, stres1d, &zdzd, &
		zdn, &dzdndl, &dzdndd, &dfdr, &dfdd, &alphzd, &unmd, dlambda, 
		&tdl, &wx, &phtrwx, &phtrwx2, rdk, rsk, teta, ncd, &melpst, 
		ntens, ndi, nshr, m, mpi);

/* ----------- END    Partial derivatives of de N */
/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* -------  Partial derivatives of dgX/dlambda, dgX/dWx and dgX/dD */
/*                                        dgR/dlambda and dgR/dWr */

	caldpg_(&dxdl, &dxdx, &dxdr, &dxdd, &dxds, &drdl, &drdr, dndl, dndx, 
		dndr, dndd, dnds, &dfdr, &dfrdl, &dmwx, tn, alpha1, dlambda, &
		tdl, &wx, &untrwx, &phtrwx, &phtrwx2, rdk, rsk, teta, ncd, 
		ntens, ndi, nshr, m, mpi);

/* ----------- END  Partial derivatives of de gX and de gR */
/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* -------  Partial derivatives of dh/dlambda, dh/dWx, dh/dWr and dh/dD */

	caldph_(&dhdl, &dhdx, &dhdr, &dhdd, &dhds, dndl, dndx, dndr, dndd, 
		dnds, &dpdd, zd, tn, &p, phi, stres1d, alpha1, altn, dlambda, 
		&tdl, h__, &unmd, &unmd1, &wx, &phtrwx, &phtrwx2, g, rdk, rsk,
		 xms, &xmm1, teta, ncd, ntens, ndi, nshr, m, mpi);

/* ----------- END    Partial derivatives of de h */
/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* -------  System solution clambda, cwx, cwr and cd */

	solsyst_(dlambda, phi, &wx, &wr, &ff, &fx, &fr, &fd, &dfdl, &dfdx, &
		dfdr, &dfdd, &dxdl, &dxdx, &dxdr, &dxdd, &drdl, &drdr, &dhdl, 
		&dhdx, &dhdr, &dhdd, ncd, ntens, ndi, nshr, m, mpi);

/* ----------- END  System solution clambda, cwx, cwr and cd */
/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* -------  Test: if D > 1 then decrease of loading time increment */
	if (*phi > un) {
		if (*m >= 0) {
		    io___496.ciunit = *mpi;
		    s_wsfe(&io___496);
		    e_wsfe();
		    io___497.ciunit = *mpi;
		    s_wsle(&io___497);
		    do_lio(&c__9, &c__1, "SUB ELPLASTDO PHI = ", (ftnlen)20);
		    do_lio(&c__5, &c__1, (char *)&(*phi), (ftnlen)sizeof(doublereal));
		    e_wsle();
		    io___498.ciunit = *mpi;
		    s_wsle(&io___498);
		    do_lio(&c__9, &c__1, "ELPLASTDO NITMAX =", (ftnlen)18);
		    do_lio(&c__3, &c__1, (char *)&nitmax, (ftnlen)sizeof(integer));
		    e_wsle();
		}
	    *lonewdt = TRUE_;
		return 0;
	}
/* ----------- END  Test: if D > 1 then decrease of loading time increment */
/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* ---------  Update of non constant repetitive coefficients */

	calcoeff_(zd, tn, &zdn, &zdzd, stres1d, estar, &p, &xkdh, &dsqhmp, &
		ryvshdl, &pold, alpha1, &pr1, phi, &phi1, dlambda, h__, &unmd,
		 &unmd1, &wr, &wx, &dmwx, yield, v, pf, g, xk, hk, rdk, rsk, 
		xms, hi, rdi, rsi, &tdl, &untrwx, &phtrwx, &phtrwx2, &untrwr, 
		&htrwr, &htrwr2, &trswr, &rtdl, teta, ncd, &nitmax, ntens, 
		ndi, nshr, m, mpi);
/* ----------- END  Update of non constant repetitive coefficients */
/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* -------  Functions F, FX, FR and FH */

	calfgh_(&ff, &fx, &fr, &fd, &zdn, &zdzd, tn, zd, alpha1, phi, &phi1, &
		unmd, &unmd1, h__, &p, dlambda, &tdl, &eyv, &wx, &dmwx, &
		phtrwx, &wr, &rtdl, &htrwr, g, xms, ncd, ntens, ndi, nshr, m, 
		mpi);
/* ----------- END  Functions F, GX, GR and H */
/* ----------------------------------------------------------------------* */
/* -------  Verifivation of convergence criteria */

	verifgence_(&convergence, &ff, &fx, &fr, &fd, epsf, epsx, epsr, epsd, 
		ntens, ndi, nshr, m, mpi);
/* ----------- END  Verifivation of convergence criteria */
/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/*      PAUSE */
    }
    if (*m >= 1) {
	io___500.ciunit = *mpi;
	s_wsfe(&io___500);
	do_fio(&c__1, (char *)&nitmax, (ftnlen)sizeof(integer));
	e_wsfe();
    }
/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* -------  Test: if NITMAX >100 then return and restart */
	if (nitmax >= 100) {
		if (*m >= 0) {
			io___501.ciunit = *mpi;
			s_wsfe(&io___501);
			e_wsfe();
		}
		*lonewdt = TRUE_;
		return 0;
    }
/* -------  END  Test: if NITMAX >100 */

/* --------- END Solution Z, DLAMBDA and ENDO */
/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* -------  Tensors EPLAS(n+1), ALPHA(n+1), PR(n+1) and STRESS(n+1) */

    solution_(&stress[1], &eplas[1], &alpha[1], xlambda, pr, strant, phi, zd, 
	    tn, stres1d, eplas1, alpha1, &pr1, dlambda, dlan, &tdl, &unmd, &
	    untrwx, &untrwr, &rtdl, &eyv, g, xk, hk, hi, xms, teta, dtime, &
	    stvit[1], &nitmax, nstvi, ntens, ndi, nshr, ncd, m, mpi);

/* -------  END  Tensors EPLAS(n+1), ALPHA(n+1), PR(n+1) and STRESS(n+1) */
/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* -------  Test: if D > DC then D = 1, STRESS and DDSDDE equal to zero */

/* freeze damage, don't allow bulk failure */
#if 0
    if (*phi > *dc) {
	*phi = un;
	aset_(&stress[1], &zero, ntens);
	i__1 = *ntens * *ntens;
	aset_(&ddsdde[ddsdde_offset], &zero, &i__1);
	i__1 = *ntens + 2;
	aset_(&stvit[1], &zero, &i__1);
	return 0;
    }
#endif

/* ----------- END  Test: if D > DC then D = 1, STRESS and DDSDDE equal to zero */
/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* -------  Tangent consistent matrix J=DDSDDE */

    jddsdde_(&ddsdde[ddsdde_offset], tn, zd, strant, stres1d, estar, alpha1, 
	    tetn, altn, dedt, &zdn, &zdzd, &p, &dsqhmp, &ryvshdl, &eyv, phi, 
	    h__, &pr1, &xkdh, dlambda, &tdl, &unmd, &unmd1, &dmwx, &wx, &wr, &
	    untrwx, &phtrwx, &phtrwx2, &htrwr, &htrwr2, &rtdl, &trswr, yield, 
	    v, pf, g, xk, rdk, rsk, rdi, rsi, xms, &xmm1, teta, ncd, ntens, 
	    ndi, nshr, m, mpi);

/* -------  END  Tangent consistent matrix J=DDSDDE */
/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* -------  Comment Check Code */
/* -------  END Comment Check Code */
/* ----------------------------------------------------------------------* */
/* ======================================================================= */
/* ----------------------------------------------------------------------* */
/* ---------   FORMATS  -------------------------------------------------- */
/* ----------------------------------------------------------------------* */
/* L101: */
/* L102: */
/* ----------------------------------------------------------------------* */
/* ======================================================================= */
    return 0;
} /* elplastdo_ */


/* Subroutine */ int jddsdde_(doublereal *ddsdde, doublereal *tn, doublereal *
	zd, doublereal *strant, doublereal *stres1d, doublereal *estar, 
	doublereal *alpha1, doublereal *tetn, doublereal *altn, doublereal *
	dedt, doublereal *zdn, doublereal *zdzd, doublereal *p, doublereal *
	dsqhmp, doublereal *ryvshdl, doublereal *eyv, doublereal *phi, 
	doublereal *h__, doublereal *pr1, doublereal *xkdh, doublereal *
	dlambda, doublereal *tdl, doublereal *unmd, doublereal *unmd1, 
	doublereal *dmwx, doublereal *wx, doublereal *wr, doublereal *untrwx, 
	doublereal *phtrwx, doublereal *phtrwx2, doublereal *htrwr, 
	doublereal *htrwr2, doublereal *rtdl, doublereal *trswr, doublereal *
	yield, doublereal *v, doublereal *pf, doublereal *g, doublereal *xk, 
	doublereal *rdk, doublereal *rsk, doublereal *rdi, doublereal *rsi, 
	doublereal *xms, doublereal *xmm1, doublereal *teta, integer *ncd, 
	integer *ntens, integer *ndi, integer *nshr, integer *m, integer *mpi)
{
    /* Initialized data */

    static integer melpst = 0;

    /* Format strings */
    static char fmt_101[] = "(1x,6(d20.14,2x))";

    /* System generated locals */
    integer ddsdde_dim1, ddsdde_offset, i__1, i__2;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void), s_wsfe(cilist *), do_fio(integer *, char *, ftnlen),
	     e_wsfe(void);

    /* Local variables */
    static integer i__, j;
    static doublereal unn[36], ddde[6], dfdd, dhdd, alde[36], dlde[6], dndd[
	    36], dnde[36], dfdl, dhdl, dpdd, dfdr, dndl[6], dfds, dxdd, dhdr, 
	    dfdx, dzde[36], dndr[6], dnds[6], drdl, tnde[36], dxdl, drdr, 
	    dndx[6], dxdr, dxds, drds, dhdx, dxdx, dhds, tntn[36], dpnde[36], 
	    dfrdl, dxide[36], uzdnn[36];
    extern /* Subroutine */ int caldde_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, integer *, integer *, integer *), caldpf_(
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *), caldpg_(
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, integer *), caldph_(doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *);
    static doublereal snddde[36];
    extern /* Subroutine */ int caldpn_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *);
    static doublereal dzdndd, dzdndl, alphzd;
    extern /* Subroutine */ int caldsig_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, integer *);

    /* Fortran I/O blocks */
    static cilist io___549 = { 0, 0, 0, 0, 0 };
    static cilist io___551 = { 0, 0, 0, fmt_101, 0 };


/* ======================================================================= */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* -------------  CALCUL DE LA MATRICE TANGENTE CONSISTENTE  ------------- */
/* --------------------  J = DDSDDE = DSTRESS/DSTRAN  -------------------- */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
/* ------------- */
/*  EN ENTREE : */
/* ------------- */
/*  TN		: TENSEUR NORMAL */
/*  ZD		: TENSEUR DE DIRECTION NORMALE */
/*  STRESS	: TENSEUR DES CONTRAINTES */
/*  STRANT	: TENSEUR DES DEFORMATIONS TOTALES */
/*  DEDT		: TENSEUR INCREMENT DE DEFORMATION : DTIME*DSTRAN */
/*  ESTAR	: TENSEUR e* DES DEFORMATIONS CONNUES */
/*  ENDO		: TENSEUR D'ENDOMMAGEMENT */
/*  ALPHA1	: TENSEUR D'ECROUISSAGE CINEMATIQUE CONNUE */
/*  TETN		: TENSEUR = ESTAR - TDLRD*TN */
/*  ALTN		: TENSEUR = ALPHA1 + TDLRD*TN */
/*  QD		: VOID COALESCENCE OPERATOR */
/*  DQD		: COALESCENCE DERIVATIVE HEAVISIDE */
/*  PLTN		: PRODUIT P:|dEp/dl| */
/*  DM		: TENSEUR D'EFFET D'ENDOMMAGEMENT (Mt) */
/*  DMI		: TENSEUR D'EFFET D'ENDOMMAGEMENT INVERSE (M(-t)) */
/*  DN		: COMPONANTES DU VECTEUR ORIENTATION (NTENS,NTENS) */
/*  DENS		: DENSITE DE DISTRIBUTION DU DOMMAGE */
/*  UNITDEV	: TENSEUR UNITE DEVIATORIQUE */
/*  ZDN		: NORME DU TENSEUR ZD (1ERE METHODE) */
/*  ZDZD		: NORME DU TENSEUR ZD = RACINE (2/3*ZD:ZD) */
/*  PR1		: VARIABLE D'ECROUISSAGE ISOTROPE CONNUE */
/*  DLAMBDA	: VITESSE DE DEFORMATION PLASTIQUE CUMULEE */
/*  YPS		: VALEUR DE LA FONCTION YPS = (-Y/S)**s */
/*  YPS1		: VALEUR DE LA FONCTION YPS = (-Y/S)**(s-1) */
/*  YV		: VALEUR DE LA FONCTION YV = Ye = 3/2*<Sh>/Seq */
/*  EYV		: VALEUR DE LA FONCTION Exp[YV] */
/*  UES		: VALEUR DE LA FONCTION 1/Seq (if <Sh>.neq.0) */
/*  TDL		: TETA*DLAMBDA */
/*  UDHM		: VALUE 1/(1-Dh)**m */
/*  UDHM1	: VALUE 1/(1-Dh)**(m+1) */
/*  WR		: ISOTROPIC DYNAMIC RECOVERY TERM: WR = DSINH[QS*R] */
/*  WX		: KINEMATIC DYNAMIC RECOVERY TERM: WX = |DMI:X| */
/*  UNTRWX	: 1+TETA*(PRD*DLAMBDA+PRS)*WX */
/*  PHTRWX	: (2/3)*PH/(1+TETA*(PRD*DLAMBDA+PRS)*WX) */
/*  PHTRWX2	: (2/3)*PH/(1+TETA*(PRD*DLAMBDA+RT(2/3)*PRS)*WX)**2 */
/*  HTRWR3	: RT(2/3)*H/(1+TETA*(RD*DLAMBDA+RS)*WR)**3 */
/*  RTDL		: PR1+TETA*DLAMBDA */
/*  TRSWR	: TETA*RS*WR */
/*  G		: MODULE DE CISAILLEMENT */
/*  XK		: COEFFICIENT DE LAME */
/*  RD		: CONSTANTE D'ECROUISSAGE ISOTROPE */
/*  H		: CONSTANTE D'ECROUISSAGE ISOTROPE */
/*  RS		: CONSTANTE D'ECROUISSAGE ISOTROPE (RECOVERY) */
/*  PRD		: CONSTANTE D'ECROUISSAGE CINEMATIQUE */
/*  PH		: CONSTANTE D'ECROUISSAGE CINEMATIQUE */
/*  PRS		: CONSTANTE D'ECROUISSAGE CINEMATIQUE (RECOVERY) */
/*  GS		: FONCTION D'ENDOMMAGEMENT DUCTILE */
/*  PS		: FONCTION D'ENDOMMAGEMENT DUCTILE */
/*  XDN		: CONSTANTE NUCLEATION XLD*XFN/(XKIC*XLF**(1/3)) */
/*  AV		: VOID GROWTH CONSTANT */
/*  VM		: VOID GROWTH CONSTANT */
/*  VC		: VOID GROWTH TEST (VC = 1 if Dh > Vc) */
/*  AC		: COALESCENCE CONSTANT */
/*  TETA		: TETA-METHODE */
/*  NCD		: NCD=1 --> COUPLAGE ELASTOPLASTICITE-ENDOMMAGEMENT */
/*  		  KCT=1 --> MATRICE TANGENTE CONSISTENTE */
/*  NDI		: NOMBRE DE COMPOSANTES DES CONTRAINTES DIRECTES */
/*  NSHR		: NOMBRE DE COMPOSANTES DES CONTRAINTES TANGENTIELLES */
/*  NTENS	: TAILLE DU TABLEAU DES CONTRAINTES OU DEFORMATIONS */
/* 		  (NDI + NSHR). */
/*  M		: COMMENTAIRES */
/*  MPI		: UNITE DE FICHIER DE SORTIE */
/* ------------- */
/*  EN SORTIE : */
/* ------------- */
/*  DDSDDE	: MATRICE TANGENTE CONSISTANTE */
/* ----------------------------------------------------------------------* */
/* ======================================================================= */
/* ----------------------------------------------------------------------* */
/* .1-----  Implicit, External */
/* .2-----  Parameter,Integer */
/* .3-----  Dimension */
/* -------- */
/* .4-----  Real,Complex,Double precision,Logical,Character */
/* .5-----  Common */
/* .6-----  Equivalence */
/* .7-----  Save,Data */
    /* Parameter adjustments */
    --stres1d;
    --h__;
    --dedt;
    --altn;
    --tetn;
    --alpha1;
    --estar;
    --strant;
    --zd;
    --tn;
    ddsdde_dim1 = *ntens;
    ddsdde_offset = 1 + ddsdde_dim1;
    ddsdde -= ddsdde_offset;

    /* Function Body */
/* .8-----  Definition of fonctions */
/* ----------------------------------------------------------------------* */
/* ======================================================================= */
/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* -------  Partial derivatives df/dlambda, df/dWx, df/dWr, df/dD and df/dLs */

    caldpf_(&dfdl, &dfdx, &dfdr, &dfdd, &dfds, &drds, &dpdd, &dzdndl, &dzdndd,
	     &alphzd, &dfrdl, &zd[1], &h__[1], p, phi, &alpha1[1], &stres1d[1]
	    , pr1, zdn, zdzd, dsqhmp, phtrwx, ryvshdl, dlambda, tdl, unmd, 
	    unmd1, xkdh, wx, wr, phtrwx2, htrwr, htrwr2, trswr, rtdl, rdk, 
	    rsk, rdi, rsi, yield, g, xk, v, pf, xms, xmm1, teta, ncd, &melpst,
	     ntens, ndi, nshr, m, mpi);

/* ----------- END  Partial derivatives of f */
/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* -------  Partial derivatives dN/dlambda, dN/dWx, dN/dWr, dN/dD and dN/ds */

    caldpn_(dndl, dndx, dndr, dndd, dnds, &zd[1], &alpha1[1], &stres1d[1], 
	    zdzd, zdn, &dzdndl, &dzdndd, &dfdr, &dfdd, &alphzd, unmd, dlambda,
	     tdl, wx, phtrwx, phtrwx2, rdk, rsk, teta, ncd, &melpst, ntens, 
	    ndi, nshr, m, mpi);

/* ----------- END Partial derivatives of N */
/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* -------  Partial derivatives dgX/dlambda, dgX/dWx and dgX/dD */
/*                                        dgR/dlambda and dgR/dWr */

    caldpg_(&dxdl, &dxdx, &dxdr, &dxdd, &dxds, &drdl, &drdr, dndl, dndx, dndr,
	     dndd, dnds, &dfdr, &dfrdl, dmwx, &tn[1], &alpha1[1], dlambda, 
	    tdl, wx, untrwx, phtrwx, phtrwx2, rdk, rsk, teta, ncd, ntens, ndi,
	     nshr, m, mpi);

/* ----------- END  Partial derivatives of gX and of gR */
/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* -------  Partial derivatives dh/dlambda, dh/dWx, dh/dWr and dh/dD */

    caldph_(&dhdl, &dhdx, &dhdr, &dhdd, &dhds, dndl, dndx, dndr, dndd, dnds, &
	    dpdd, &zd[1], &tn[1], p, phi, &stres1d[1], &alpha1[1], &altn[1], 
	    dlambda, tdl, &h__[1], unmd, unmd1, wx, phtrwx, phtrwx2, g, rdk, 
	    rsk, xms, xmm1, teta, ncd, ntens, ndi, nshr, m, mpi);

/* ----------- END  Partial derivatives of h */
/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* -------  Partial derivatives in respect with  E = STRANT */

    caldde_(dlde, ddde, dnde, unn, tntn, uzdnn, alde, dpnde, dzde, dxide, &
	    estar[1], &alpha1[1], &zd[1], &tn[1], &tetn[1], &altn[1], phi, 
	    zdzd, &dfdl, &dfdx, &dfdr, &dfdd, &dxdl, &dxdx, &dxdd, &drdl, &
	    drdr, &dhdl, &dhdx, &dhdr, &dhdd, &dfds, &dxds, &drds, &dhds, &
	    dedt[1], p, eyv, dlambda, tdl, unmd, unmd1, dmwx, wx, phtrwx, 
	    phtrwx2, &h__[1], g, xk, rdk, rsk, rsi, xms, teta, ncd, ntens, 
	    ndi, nshr, m, mpi);

/* -------- END  Partial derivatives in respect with  E = STRANT */
/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* -------  Tangent consistent matrix DDSDDE */

    caldsig_(&ddsdde[ddsdde_offset], dnde, tnde, snddde, dlde, ddde, &tn[1], &
	    strant[1], &stres1d[1], dlambda, unmd, unmd1, g, xk, ntens, ndi, 
	    m, mpi);

/* -------  END  Tangent consistent matrix DDSDDE */
/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* -------  Comment Check Code */

    if (*m >= 4) {
	io___549.ciunit = *mpi;
	s_wsle(&io___549);
	do_lio(&c__9, &c__1, "JDDSDDE DDSDDE", (ftnlen)14);
	e_wsle();
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___551.ciunit = *mpi;
	    s_wsfe(&io___551);
	    i__2 = *ntens;
	    for (j = 1; j <= i__2; ++j) {
		do_fio(&c__1, (char *)&ddsdde[i__ + j * ddsdde_dim1], (ftnlen)
			sizeof(doublereal));
	    }
	    e_wsfe();
	}
    }

/* -------  END Comment Check Code */
/* -----------------------------------------------------------------------* */
/* ---------   FORMATS  --------------------------------------------------* */
/* -----------------------------------------------------------------------* */
/* L102: */
/* -----------------------------------------------------------------------* */
/* ======================================================================== */
    return 0;
} /* jddsdde_ */


/* Subroutine */ int ktimen_(doublereal *stress, doublereal *stran, 
	doublereal *dstran, doublereal *eplas, doublereal *alpha, doublereal *
	p, doublereal *pr, doublereal *phi, doublereal *xkdh, doublereal *
	dlambda, doublereal *dlan, doublereal *strant, doublereal *dedt, 
	doublereal *estar, doublereal *gestar, doublereal *stres1d, 
	doublereal *eplas1, doublereal *alpha1, doublereal *pr1, doublereal *
	phi1, doublereal *pold, doublereal *wr, doublereal *wx, doublereal *
	xk, doublereal *g, doublereal *hk, doublereal *rdk, doublereal *rsk, 
	doublereal *hi, doublereal *rdi, doublereal *rsi, doublereal *xms, 
	doublereal *teta, doublereal *dtime, integer *ncd, integer *ntens, 
	integer *ndi, integer *nshr, integer *m, integer *mpi)
{
    /* Initialized data */

    static doublereal zero = 0.;
    static doublereal one = 1.;
    static doublereal two = 2.;
    static doublereal three = 3.;
    static doublereal half = .5;
    static doublereal thalf = 1.5;

    /* Format strings */
    static char fmt_101[] = "(1x,6(d20.14,2x))";
    static char fmt_102[] = "(1x,a20,6(d20.14,2x))";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void), s_wsfe(cilist *), do_fio(integer *, char *, ftnlen),
	     e_wsfe(void);

    /* Local variables */
    static doublereal onethird, stressdv[6];
    static integer i__;
    static doublereal r__, z__[6], unmtetadt, tn[6], zz, d2g, rwr, 
	    unmtetadtdl, unmd, tdzz, xkdeh;
    extern /* Subroutine */ int pdtsca_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal unmdhk, rdlswx;
    extern /* Subroutine */ int tensdev_(doublereal *, doublereal *, integer *
	    , integer *);
    static doublereal unmteta;

    /* Fortran I/O blocks */
    static cilist io___576 = { 0, 0, 0, 0, 0 };
    static cilist io___577 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___578 = { 0, 0, 0, 0, 0 };
    static cilist io___579 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___580 = { 0, 0, 0, 0, 0 };
    static cilist io___581 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___582 = { 0, 0, 0, 0, 0 };
    static cilist io___583 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___584 = { 0, 0, 0, 0, 0 };
    static cilist io___585 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___586 = { 0, 0, 0, 0, 0 };
    static cilist io___587 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___588 = { 0, 0, 0, 0, 0 };
    static cilist io___589 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___590 = { 0, 0, 0, 0, 0 };
    static cilist io___591 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___592 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___593 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___594 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___595 = { 0, 0, 0, fmt_102, 0 };


/* ======================================================================= */
/* ----------------------------------------------------------------------- */
/* -------------  CALCUL DES TERMES E*, e*, ALPHA1, r1, D1  -------------- */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
/* ------------- */
/*  EN ENTREE : */
/* ------------- */
/*  STRESS	: TENSEUR DES CONTRAINTES */
/*  STRAN	: DEFORMATION TOTALE */
/*  EPLAS	: TENSEUR DES DEFORMATIONS PLASTIQUES */
/*  ALPHA	: TENSEUR DES CONTRAINTES CINEMATIQUES */
/*  PR		: VARIABLE D'ECROUISSAGE ISOTROPE */
/*  PHI		: VARIABLE D'ENDOMMAGEMENT */
/*  DLAN		: INCREMENT DEFORMATION PLASTIQUE CUMULEE PRECEDENT */
/*  NTENS	: LONGUEUR DES TABLEAU DSTRAN, STRESS, ... */
/*  G		: MODULE DE CISAILLEMENT */
/*  RD		: CONSTANTE D'ECROUISSAGE ISOTROPE */
/*  H		: CONSTANTE D'ECROUISSAGE ISOTROPE */
/*  RS		: CONSTANTE D'ECROUISSAGE ISOTROPE (RECOVERY) */
/*  PRD		: CONSTANTE D'ECROUISSAGE CINEMATIQUE */
/*  PH		: CONSTANTE D'ECROUISSAGE CINEMATIQUE */
/*  PRS		: CONSTANTE D'ECROUISSAGE CINEMATIQUE (RECOVERY) */
/*  XDN		: CONSTANTE XLD*XFN/(XKIC*XLF**(1/3)) */
/*  GS		: FONCTION D'ENDOMMAGEMENT DUCTILE */
/*  PS		: FONCTION D'ENDOMMAGEMENT DUCTILE */
/*  DC		: VALEUR D'ENDOMAGEMMENT CRITIQUE */
/*  CTE		: CONSTANTE DE TENSION */
/*  CCO		: CONSTANTE DE COMPRESSION */
/*  CTO		: CONSTANTE DE TORSION */
/*  TETA		: TETA-METHODE */
/*  DTIME	: INCREMENT DE TEMPS */
/*  NCD		: NCD=1 --> COUPLAGE ELASTOPLASTICITE-ENDOMMAGEMENT */
/*  		  NCD=0 --> ELASTOPLASTICITE SANS ENDOMMAGEMENT */
/*  NDI		: NOMBRE DE COMPOSANTES DES CONTRAINTES DIRECTES */
/*  NSHR		: NOMBRE DE COMPOSANTES DES CONTRAINTES TANGENTIELLES */
/*  NTENS	: TAILLE DU TABLEAU DES CONTRAINTES OU DEFORMATIONS */
/* 		  (NDI + NSHR). */
/*  M		: COMMENTAIRES */
/*  MPI		: UNITE DE FICHIER DE SORTIE */
/* ------------- */
/*  EN SORTIE : */
/* ------------- */
/*  STRANT	: TENSEUR DES DEFORMATIONS TOTALES : STRAN+DSTRAN */
/*  DEDT		: TENSEUR INCREMENT DE DEFORMATION : DSTRAN */
/*  ESTAR	: TENSEUR DES DEFORMATIONS CONNUS A T+DT */
/*  GESTAR	: TENSEUR DEVIATORIQUE DE ESTAR */
/*  STRES1D	: TENSEUR DEVIATORIQUE DES CONTRAINTES CONNU A T */
/*  EPLAS1	: TENSEUR DES DEFORMATIONS PLASTIQUES CONNUS A T+DT */
/*  ALPHA1	: TENSEUR DES CONTRAINTES CINEMATIQUES CONNUS A T+DT */
/*  PR1		: VARIABLE D'ECROUISSAGE ISOTROPE CONNUS A T+DT */
/*  PHI1		: TENSEUR D'ENDOMMAGEMENT CONNUS A T+DT */
/*  P		: HYDROSTATIC STRESS */
/*  POLD		: HYDROSTATIC STRESS AT TIME T */
/*  PLTN		: PRODUIT P:|dEp/dl| */
/*  DN		: COMPONANTES DU VECTEUR ORIENTATION */
/*  DLAMBDA	: INCREMENT DE DEFORMATION PLASTIQUE EQUIVALENTE */
/*  WR		: ISOTROPIC DYNAMIC RECOVERY TERM: WR = R */
/*  WX		: KINEMATIC DYNAMIC RECOVERY TERM: WX = |DMI:X| */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
/* .1-----  Implicit, External */
/* .2-----  Parameter,Integer */
/* .3-----  Dimension */
/* ------- */
/* .4-----  Real,Complex,Double precision,Logical,Character */
/* .5-----  Common */
/* .6-----  Equivalence */
/* .7-----  Save,Data */
    /* Parameter adjustments */
    --alpha1;
    --eplas1;
    --stres1d;
    --gestar;
    --estar;
    --dedt;
    --strant;
    --alpha;
    --eplas;
    --dstran;
    --stran;
    --stress;

    /* Function Body */
/* .8-----  Definition de fonctions */
/* ---------------------------------------------------------------------* */
/* ====================================================================== */
/* ---------------------------------------------------------------------* */
/* -------  Initialisation */
    onethird = one / three;
    unmteta = one - *teta;
    unmtetadt = unmteta * *dtime;
    unmtetadtdl = unmtetadt * *dlan;
/* -------  FIN Initialisation */
/* ---------------------------------------------------------------------* */
/* -------  Prediction Increment de deformation plastique cumulee */
    if (*teta == one) {
	*dlambda = *dlan * *dtime;
    } else {
	*dlambda = zero;
    }
/* -------  END  Prediction Increment de deformation plastique cumulee */
/* ---------------------------------------------------------------------* */
/* -------  Total Strain Tensor STRANT = STRAN + DSTRAN */

    i__1 = *ndi;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dedt[i__] = dstran[i__];
	strant[i__] = stran[i__] + dedt[i__];
    }
    i__1 = *ntens;
    for (i__ = *ndi + 1; i__ <= i__1; ++i__) {
	dedt[i__] = half * dstran[i__];
	strant[i__] = half * stran[i__] + dedt[i__];
    }
/* -------  END  Total Strain Tensor STRANT = STRAN + DSTRAN */
/* ---------------------------------------------------------------------* */
/* -------  Tensors Z and TN at time t */
    unmd = one - *phi;
    unmdhk = unmd * *hk;
    tensdev_(&stress[1], stressdv, ndi, nshr);
/* ------- */
    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
	z__[i__ - 1] = stressdv[i__ - 1] - *hk * alpha[i__];
    }
/* ------- */
    pdtsca_(ndi, nshr, z__, z__, &zz);
    if (zz == zero) {
	zz = one;
    }
    tdzz = thalf / sqrt(thalf * zz);
    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
	tn[i__ - 1] = tdzz * z__[i__ - 1];
    }
/* -------  FIN Tensors Z and TN at time t */
/* ---------------------------------------------------------------------* */
/* -------  Norm WX */
    pdtsca_(ndi, nshr, &alpha[1], &alpha[1], wx);
    *wx = unmdhk * sqrt(*wx);
/* -------  END  Norm WX */
/* ---------------------------------------------------------------------* */
/* -------  Tensor Ep1 */
    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
	eplas1[i__] = eplas[i__] + unmtetadtdl * tn[i__ - 1];
    }
/* -------  END  Tensor Ep1 */
/* ---------------------------------------------------------------------* */
/* -------  Tensor e* */
    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
	estar[i__] = strant[i__] - eplas1[i__];
    }
/* -------  END  Tensor e* */
/* ---------------------------------------------------------------------* */
/* -------  Tensor E* (Deviatoric tensor of e*) */
    tensdev_(&estar[1], &gestar[1], ndi, nshr);
/* -------  END  Tensor E* */
/* ---------------------------------------------------------------------* */
/* -------  Tensor Sigma1D */
    d2g = two * *g;
/* ------- */
    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
	stres1d[i__] = d2g * gestar[i__];
    }
/* -------  END  Tensor Sigma1D */
/* ---------------------------------------------------------------------* */
/* -------  Tensor ALPHA1 */
    rdlswx = (*rdk * *dlan + *rsk) * *wx;
/* ------- */
    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
	alpha1[i__] = alpha[i__] + unmtetadt * (*dlan * tn[i__ - 1] - rdlswx *
		 alpha[i__]);
    }
/* -------  FIN Tensor ALPHA1 */
/* ---------------------------------------------------------------------* */
/* -------  Hydrostratic Stress */
    *pold = onethird * (stress[1] + stress[2] + stress[3]);
    xkdeh = *xk * (dstran[1] + dstran[2] + dstran[3]);
    *xkdh = *xk * (strant[1] + strant[2] + strant[3]);
    *p = *pold + unmd * xkdeh;
/* -------  END  Tensor Sigma1D */
/* ---------------------------------------------------------------------* */
/* -------  Calcul de r1 */
    r__ = *hi * *pr;
    *wr = r__;
    rwr = *pr * *wr;
    *pr1 = *pr + unmtetadtdl * (one - *rdi * rwr) - unmtetadt * *rsi * rwr;
/* -------  END  Calcul de r1 */
/* ---------------------------------------------------------------------* */
/* -------  Damage Phi1 */

	/* keep damage even if it is not evolving */
	*phi1 = *phi;

#if 0
    if (*ncd == 0) {
	*phi1 = zero;
    } else {
	*phi1 = *phi;
    }
#endif

/* -------  END  Damage Phi1 */
/* ---------------------------------------------------------------------* */
/* -------  Comment Check Code */
    if (*m >= 4) {
	io___576.ciunit = *mpi;
	s_wsle(&io___576);
	do_lio(&c__9, &c__1, "KTIMEN TN", (ftnlen)9);
	e_wsle();
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___577.ciunit = *mpi;
	    s_wsfe(&io___577);
	    do_fio(&c__1, (char *)&tn[i__ - 1], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	io___578.ciunit = *mpi;
	s_wsle(&io___578);
	do_lio(&c__9, &c__1, "KTIMEN STRAN", (ftnlen)12);
	e_wsle();
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___579.ciunit = *mpi;
	    s_wsfe(&io___579);
	    do_fio(&c__1, (char *)&stran[i__], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	io___580.ciunit = *mpi;
	s_wsle(&io___580);
	do_lio(&c__9, &c__1, "KTIMEN STRANT", (ftnlen)13);
	e_wsle();
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___581.ciunit = *mpi;
	    s_wsfe(&io___581);
	    do_fio(&c__1, (char *)&strant[i__], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	io___582.ciunit = *mpi;
	s_wsle(&io___582);
	do_lio(&c__9, &c__1, "KTIMEN ESTAR", (ftnlen)12);
	e_wsle();
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___583.ciunit = *mpi;
	    s_wsfe(&io___583);
	    do_fio(&c__1, (char *)&estar[i__], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	io___584.ciunit = *mpi;
	s_wsle(&io___584);
	do_lio(&c__9, &c__1, "KTIMEN GESTAR", (ftnlen)13);
	e_wsle();
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___585.ciunit = *mpi;
	    s_wsfe(&io___585);
	    do_fio(&c__1, (char *)&gestar[i__], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	io___586.ciunit = *mpi;
	s_wsle(&io___586);
	do_lio(&c__9, &c__1, "KTIMEN EPLAS1", (ftnlen)13);
	e_wsle();
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___587.ciunit = *mpi;
	    s_wsfe(&io___587);
	    do_fio(&c__1, (char *)&eplas1[i__], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	io___588.ciunit = *mpi;
	s_wsle(&io___588);
	do_lio(&c__9, &c__1, "KTIMEN ALPHA1", (ftnlen)13);
	e_wsle();
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___589.ciunit = *mpi;
	    s_wsfe(&io___589);
	    do_fio(&c__1, (char *)&alpha1[i__], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	io___590.ciunit = *mpi;
	s_wsle(&io___590);
	do_lio(&c__9, &c__1, "KTIMEN STRES1D", (ftnlen)14);
	e_wsle();
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___591.ciunit = *mpi;
	    s_wsfe(&io___591);
	    do_fio(&c__1, (char *)&stres1d[i__], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	io___592.ciunit = *mpi;
	s_wsfe(&io___592);
	do_fio(&c__1, "KTIMEN PR1  =", (ftnlen)13);
	do_fio(&c__1, (char *)&(*pr1), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___593.ciunit = *mpi;
	s_wsfe(&io___593);
	do_fio(&c__1, "KTIMEN PHI1 =", (ftnlen)13);
	do_fio(&c__1, (char *)&(*phi1), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___594.ciunit = *mpi;
	s_wsfe(&io___594);
	do_fio(&c__1, "KTIMEN WR   =", (ftnlen)13);
	do_fio(&c__1, (char *)&(*wr), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___595.ciunit = *mpi;
	s_wsfe(&io___595);
	do_fio(&c__1, "KTIMEN WX   =", (ftnlen)13);
	do_fio(&c__1, (char *)&(*wx), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
/* -------  FIN Comment Check Code */
/* ----------------------------------------------------------------------- */
/* ---------   FORMATS  -------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
    return 0;
} /* ktimen_ */


/* Subroutine */ int pdtmat_(integer *n, doublereal *a, doublereal *b, 
	doublereal *ab)
{
    /* Initialized data */

    static doublereal zero = 0.;

    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j;

/* ====================================================================== */
/* ---------------------------------------------------------------------- */
/* -------  CALCUL DU PRODUIT CONTRACTE D'UN TENSEUR D'ORDRE 4  --------- */
/* --------------  AVEC UN TENSEUR D'ORDRE 2 SYMETRIQUE  ---------------- */
/* ---------------------------------------------------------------------- */
/* ====================================================================== */
/* ------------- */
/*  EN ENTREE : */
/* ------------- */
/*  A(N,N)	: TENSEURS D'ORDRE QUATRE */
/*  B(N)		: TENSEURS D'ORDRE DEUX */
/* ------------- */
/*  EN SORTIE : */
/* ------------- */
/*  AB(N)	: PRODUIT CONTRACTE A:B */
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
    --ab;
    --b;
    a_dim1 = *n;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
/* .8-----  Definition de fonctions */
/* ======================================================================= */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ab[i__] = zero;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    ab[i__] += a[i__ + j * a_dim1] * b[j];
	}
    }
/* ======================================================================= */
    return 0;
} /* pdtmat_ */


/* Subroutine */ int pdtsca_(integer *ndi, integer *nshr, doublereal *a, 
	doublereal *b, doublereal *ab)
{
    /* Initialized data */

    static doublereal zero = 0.;
    static doublereal deux = 2.;

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

/* ====================================================================== */
/* ---------------------------------------------------------------------- */
/* ----------  CALCUL DU PRODUIT CONTRACTE DE DEUX TENSEURS   ----------- */
/* --------------------   D'ORDRE 2 SYMETRIQUES  ------------------------ */
/* ---------------------------------------------------------------------- */
/* ====================================================================== */
/* ------------- */
/*  EN ENTREE : */
/* ------------- */
/*  NDI		:  NOMBRE DE COMPOSANTES DES CONTRAINTES DIRECTES. */
/*  NSHR		:  NOMBRE DE COMPOSANTES DES CONTRAINTES TANGENTIELLES. */
/*  A,B		:  TENSEURS A ET B */
/* ------------- */
/*  EN SORTIE : */
/* ------------- */
/*  AB		:  PRODUIT CONTRACTE A:B */

/* ======================================================================= */
/* .1-----  Implicit, External */
/* .2-----  Parameter */
/* .3-----  Dimension */
/* .4-----  Real,Integer,Complex,Double precision,Logical,Character */
/* .5-----  Common */
/* .6-----  Equivalence */
/* .7-----  Data */
    /* Parameter adjustments */
    --b;
    --a;

    /* Function Body */
/* .8-----  Definition de fonctions */
/* ======================================================================= */
    *ab = zero;
    i__1 = *ndi;
    for (i__ = 1; i__ <= i__1; ++i__) {
	*ab += a[i__] * b[i__];
    }
    i__1 = *nshr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	*ab += deux * a[*ndi + i__] * b[*ndi + i__];
    }
/* ======================================================================= */
    return 0;
} /* pdtsca_ */


/* Subroutine */ int pdtten_(doublereal *a, doublereal *b, doublereal *atb, 
	integer *ntens, integer *ndi)
{
    /* Initialized data */

    static doublereal deux = 2.;

    /* System generated locals */
    integer atb_dim1, atb_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j;

/* ====================================================================== */
/* ---------------------------------------------------------------------- */
/* ----------  CALCUL DU PRODUIT TENSORIEL DE DEUX TENSEURS   ----------- */
/* --------------------   D'ORDRE 2 SYMETRIQUES  ------------------------ */
/* ---------------------------------------------------------------------- */
/* ====================================================================== */
/* ------------- */
/*  EN ENTREE : */
/* ------------- */
/*  NTENS	: DIMENSION DES TENSEURS */
/*  A et B	: TENSEURS D'ORDRE DEUX */
/* ------------- */
/*  EN SORTIE : */
/* ------------- */
/*  ATB		: ATB = A * B --->  TENSEUR D'ORDRE QUATRE */
/* ======================================================================= */
/* .1-----  Implicit, External */
/* .2-----  Parameter */
/* .3-----  Dimension */
/* .4-----  Real,Integer,Complex,Double precision,Logical,Character */
/* .5-----  Common */
/* .6-----  Equivalence */
/* .7-----  Data */
    /* Parameter adjustments */
    atb_dim1 = *ntens;
    atb_offset = 1 + atb_dim1;
    atb -= atb_offset;
    --b;
    --a;

    /* Function Body */
/* .8-----  Definition de fonctions */
/* ======================================================================= */
    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *ntens;
	for (j = 1; j <= i__2; ++j) {
	    atb[i__ + j * atb_dim1] = a[i__] * b[j];
	}
    }
/* ------- */
    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *ntens;
	for (j = *ndi + 1; j <= i__2; ++j) {
	    atb[i__ + j * atb_dim1] = deux * atb[i__ + j * atb_dim1];
	}
    }
/* ======================================================================= */
    return 0;
} /* pdtten_ */


/* Subroutine */ int prdmat_(integer *ntens, doublereal *vm1, doublereal *vm2,
	 doublereal *vm3)
{
    /* Initialized data */

    static doublereal zero = 0.;

    /* System generated locals */
    integer vm1_dim1, vm1_offset, vm2_dim1, vm2_offset, vm3_dim1, vm3_offset, 
	    i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k;
    static doublereal r__;

/* ====================================================================== */
/* ---------------------------------------------------------------------- */
/* -----------------------  PRODUIT DE DEUX MATRICES  ------------------- */
/* ---------------------------------------------------------------------- */
/* ====================================================================== */
/* ------------- */
/*  EN ENTREE : */
/* ------------- */
/*  VM1		: MATRICE OU TENSEUR D'ORDRE 4 (NDIM1,NDIM2) */
/*  VM2		: MATRICE OU TENSEUR D'ORDRE 4 (NDIM2,NDIM3) */
/* ------------- */
/*  EN SORTIE : */
/* ------------- */
/*  VM3		: MATRICE OU TENSEUR D'ORDRE 4 (NDIM1,NDIM3) */
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
    vm3_dim1 = *ntens;
    vm3_offset = 1 + vm3_dim1;
    vm3 -= vm3_offset;
    vm2_dim1 = *ntens;
    vm2_offset = 1 + vm2_dim1;
    vm2 -= vm2_offset;
    vm1_dim1 = *ntens;
    vm1_offset = 1 + vm1_dim1;
    vm1 -= vm1_offset;

    /* Function Body */
/* .8-----  Definition de fonctions */
/* ======================================================================= */
    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *ntens;
	for (j = 1; j <= i__2; ++j) {
	    r__ = zero;
	    i__3 = *ntens;
	    for (k = 1; k <= i__3; ++k) {
		r__ += vm1[i__ + k * vm1_dim1] * vm2[k + j * vm2_dim1];
	    }
	    vm3[i__ + j * vm3_dim1] = r__;
	}
    }
/* ======================================================================= */
    return 0;
} /* prdmat_ */


/* Subroutine */ int predelas_(doublereal *dstran, doublereal *stress, 
	doublereal *alpha, doublereal *pr, doublereal *phi, doublereal *
	stresselas, doublereal *alphaelas, doublereal *prelas, doublereal *hi,
	 doublereal *rsi, doublereal *hk, doublereal *rsk, doublereal *g, 
	doublereal *xk, integer *nplan, integer *ntens, integer *ndi, integer 
	*nshr, integer *m, integer *mpi)
{
    /* Initialized data */

    static doublereal one = 1.;
    static doublereal two = 2.;

    /* Format strings */
    static char fmt_101[] = "(1x,6(d20.14,2x))";
    static char fmt_102[] = "(1x,a20,6(d20.14,2x))";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void), s_wsfe(cilist *), do_fio(integer *, char *, ftnlen),
	     e_wsfe(void);

    /* Local variables */
    static doublereal dstrandv[6];
    static integer i__;
    static doublereal d2g, d2gd, dehk, unmd, dmwx;
    extern /* Subroutine */ int sinv_(doublereal *, doublereal *, doublereal *
	    , integer *, integer *);
    static doublereal dmwxh, unrswx, dstran2[6], dstranh;
    extern /* Subroutine */ int tensdev_(doublereal *, doublereal *, integer *
	    , integer *);

    /* Fortran I/O blocks */
    static cilist io___624 = { 0, 0, 0, 0, 0 };
    static cilist io___625 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___626 = { 0, 0, 0, 0, 0 };
    static cilist io___627 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___628 = { 0, 0, 0, fmt_102, 0 };


/* ====================================================================== */
/* ---------------------------------------------------------------------- */
/* -------------  CALCUL D'UNE SOLUTION PUREMENT ELASTIQUE  ------------- */
/* ---------------------------------------------------------------------- */
/* ====================================================================== */
/* ------------- */
/*  EN ENTREE : */
/* ------------- */
/*  DSTRAN	: INCREMENT DE DEFORMATION */
/*  STRESS	: CONTRAINTES A L'INSTANT N */
/*  ALPHA	: TENSEUR DES CONTRAINTES CINEMATIQUES */
/*  PR		: VARIABLE D'ECROUISSAGE ISOTROPE */
/*  PHI		: DAMAGE */
/*  H		: CONSTANTE D'ECROUISSAGE ISOTROPE */
/*  RS		: CONSTANTE D'ECROUISSAGE ISOTROPE (RECOVERY) */
/*  PH		: CONSTANTE D'ECROUISSAGE CINEMATIQUE */
/*  PRS		: CONSTANTE D'ECROUISSAGE CINEMATIQUE (RECOVERY) */
/*  G		: MODULE DE CISAILLEMENT */
/*  XK		: COEFFICIENT DE LAME */
/*  XNU		: COEFFICIENT DE POISSON */
/*  NTENS	: LONGUEUR DES TABLEAU DSTRAN, STRESS, ... */
/*  NDI		: NOMBRE DE COMPOSANTES DIRECTES D'UN TENSEUR */
/*  NSHR		: NOMBRE DE COMPOSANTES TANGENTIELLES D'UN TENSEUR */
/* 		  POINT D'INTEGRATION. */
/*  M		: COMMENTAIRES */
/*  MPI		: UNITE DE FICHIER DE SORTIE */
/* ------------- */
/*  EN SORTIE : */
/* ------------- */
/*  STRESSELAS	:  CONTRAINTES A L'INSTANT N+1 AVEC PREDICTION ELASTIQUE */
/*  ALPHAELAS	:  ECROUISSAGE CINEMATIQUE A L'INSTANT N+1 (TRIAL) */
/*  PRELAS	:  ECROUISSAGE ISOTROPE A L'INSTANT N+1 (TRIAL) */
/*  UNIT		:  TENSEUR UNITE 1*1 */
/*  UNITDEV	:  TENSEUR UNITE DEVIATORIQUE Idev=I-1/3(1*1) */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
/* ----------------------------------------------------------------------- */
/* .1-----  Implicit, External */
/* .2-----  Parameter,Integer */
/* .3-----  Dimension */
/*     DOUBLE PRECISION UNITDEV(NTENS,NTENS) */
/* ------- */
/* .4-----  Real,Complex,Double precision,Logical,Character */
/* .5-----  Common */
/* .6-----  Equivalence */
/* .7-----  Save,Data */
    /* Parameter adjustments */
    --alphaelas;
    --stresselas;
    --alpha;
    --stress;
    --dstran;

    /* Function Body */
/* .8-----  Definition de fonctions */
/* ---------------------------------------------------------------------* */
/* ====================================================================== */
/* ---------------------------------------------------------------------* */
/* -------  Division par deux des deformations tangentielles */

    i__1 = *ndi;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dstran2[i__ - 1] = dstran[i__];
    }
    i__1 = *ntens;
    for (i__ = *ndi + 1; i__ <= i__1; ++i__) {
	dstran2[i__ - 1] = dstran[i__] / two;
    }

/* -------  FIN Division par deux des deformations tangentielles */
/* ---------------------------------------------------------------------* */
/* -------  Calcul des contraintes modifiees */

/*      IF(NPLAN.EQ.2) DSTRAN2(3)=-XNU/(ONE-XNU)*(DSTRAN2(1)+DSTRAN2(2)) */
/* ------- */
    unmd = one - *phi;
    d2g = two * *g;
    d2gd = unmd * d2g;
    dstranh = dstran[1] + dstran[2] + dstran[3];
    dehk = unmd * *xk * dstranh;
    tensdev_(dstran2, dstrandv, ndi, nshr);
    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
	stresselas[i__] = stress[i__] + d2gd * dstrandv[i__ - 1];
    }
    i__1 = *ndi;
    for (i__ = 1; i__ <= i__1; ++i__) {
	stresselas[i__] += dehk;
    }

/* -------  FIN Calcul des contraintes modifiees */
/* ---------------------------------------------------------------------* */
/* -------  Calcul des contraintes cinematiques modifiees */

    sinv_(&alpha[1], &dmwxh, &dmwx, ndi, nshr);
    dmwx = *hk * dmwx;
    unrswx = one - *rsk * dmwx;
    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
	alphaelas[i__] = unrswx * alpha[i__];
    }

/* -------  FIN Calcul des contraintes cinematiques modifiees */
/* ---------------------------------------------------------------------* */
/* -------  Calcul des contraintes isotropes modifiees */

    *prelas = (one - *rsi * *hi * *pr) * *pr;

/* -------  FIN Calcul des contraintes isotropes modifiees */
/* ---------------------------------------------------------------------* */
/* -------  Comment Check Code */
    if (*m >= 4) {
	io___624.ciunit = *mpi;
	s_wsle(&io___624);
	do_lio(&c__9, &c__1, "PREDELAS STRESSELAS", (ftnlen)19);
	e_wsle();
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___625.ciunit = *mpi;
	    s_wsfe(&io___625);
	    do_fio(&c__1, (char *)&stresselas[i__], (ftnlen)sizeof(doublereal)
		    );
	    e_wsfe();
	}
	io___626.ciunit = *mpi;
	s_wsle(&io___626);
	do_lio(&c__9, &c__1, "PREDELAS ALPHAELAS", (ftnlen)18);
	e_wsle();
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___627.ciunit = *mpi;
	    s_wsfe(&io___627);
	    do_fio(&c__1, (char *)&alphaelas[i__], (ftnlen)sizeof(doublereal))
		    ;
	    e_wsfe();
	}
	io___628.ciunit = *mpi;
	s_wsfe(&io___628);
	do_fio(&c__1, "PREDELAS  PRELAS  = ", (ftnlen)20);
	do_fio(&c__1, (char *)&(*prelas), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
/* -------  FIN Comment Check Code */
/* ----------------------------------------------------------------------* */
/* ---------   FORMATS  -------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* ======================================================================= */
    return 0;
} /* predelas_ */


/* Subroutine */ int recov_(doublereal *props, integer *nprops, doublereal *
	statev, integer *nstatv, doublereal *stress, doublereal *stran, 
	doublereal *dstran, doublereal *eplas, doublereal *alpha, doublereal *
	drot, doublereal *dfgrd0, doublereal *dfgrd1, doublereal *xlambda, 
	doublereal *pr, doublereal *phi, doublereal *dlan, doublereal *time, 
	doublereal *dtime, doublereal *temp, doublereal *dtemp, doublereal *e,
	 doublereal *xnu, doublereal *g, doublereal *xk, doublereal *yield, 
	doublereal *v, doublereal *pf, doublereal *hi, doublereal *rdi, 
	doublereal *rsi, doublereal *hk, doublereal *rdk, doublereal *rsk, 
	doublereal *xip, doublereal *xms, doublereal *dc, doublereal *teta, 
	doublereal *epsf, doublereal *epsx, doublereal *epsr, doublereal *
	epsd, doublereal *vnewdt, integer *ncd, integer *nplan, integer *
	natur, integer *irot, integer *ntens, integer *ndi, integer *nshr, 
	integer *nvi, integer *m, integer *mpi)
{
    /* Initialized data */

    static doublereal zero = 0.;
    static doublereal un = 1.;
    static doublereal deux = 2.;
    static doublereal trois = 3.;
    static doublereal udemi = .5;

    /* Format strings */
    static char fmt_101[] = "(1x,3(d20.14,2x))";
    static char fmt_102[] = "(1x,\002XLAMBDA=\002,d20.14,2x,\002PHI=\002,d20"
	    ".14,2x,\002PR=\002,d20.14,2x,\002DLAN=\002,d20.14,2x)";
    static char fmt_103[] = "(\002*-----------------------------------------"
	    "------------*\002//,\002*-----   Material Constants   -----*\002"
	    "/,\002  V(T)      = \002,e12.5/,\002  Y(T)      = \002,e12.5/"
	    ",\002  f(T)      = \002,e12.5/,\002  rd(T)     = \002,e12.5/,"
	    "\002  h(T)      = \002,e12.5/,\002  rs(T)     = \002,e12.5/,\002"
	    "  Rd(T)     = \002,e12.5/,\002  H(T)      = \002,e12.5/,\002  Rs"
	    "(T)     = \002,e12.5/)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal), tanh(doublereal);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void), s_wsfe(cilist *), do_fio(integer *, char *, ftnlen),
	     e_wsfe(void);

    /* Local variables */
    static doublereal tunmdxnu, c__;
    static integer i__, j;
    static doublereal t0;
    extern /* Subroutine */ int transpose_(integer *, doublereal *, 
	    doublereal *);
    static doublereal dr[9]	/* was [3][3] */, ya, yb, ut0, rho, cts[22], 
	    rot0[9]	/* was [3][3] */, rot1[9]	/* was [3][3] */, 
	    rttd, trot0[9]	/* was [3][3] */, dtier;
    extern /* Subroutine */ int dpolr_(doublereal *, doublereal *, integer *, 
	    integer *);
    static doublereal tdrot[9]	/* was [3][3] */, utier;
    extern /* Subroutine */ int prdmat_(integer *, doublereal *, doublereal *,
	     doublereal *);
    static doublereal drotgn[9]	/* was [3][3] */;
    extern /* Subroutine */ int rotsig_(doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *);
    static doublereal stressi[6];

    /* Fortran I/O blocks */
    static cilist io___654 = { 0, 0, 0, 0, 0 };
    static cilist io___655 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___656 = { 0, 0, 0, 0, 0 };
    static cilist io___657 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___658 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___659 = { 0, 0, 0, fmt_103, 0 };


/* ====================================================================== */
/* ---------------------------------------------------------------------- */
/* ------  EXTRACTION DES VARIABLES UTILISATEURS ET DES CONSTANTES  ----- */
/* ---------------------------------------------------------------------- */
/* ====================================================================== */
/* ----------- */
/*  ENTREES : */
/* ----------- */
/*  STAVEV(NSTATEV) : VARIABLES UTILISATEUR POUR LE MODELE UMAT */
/*  PROPS(NPROPS)   : CONSTANTES DEFINIES DANS *USER MATERIAL */
/* ----------- */
/*  SORTIES : */
/* ----------- */
/*  STRAN	: TENSEUR DES DEFORMATIONS TOTALES */
/*  DSTRAN	: TENSEUR INCREMENT DES DEFORMATIONS */
/*  EPLAS	: TENSEUR DES DEFORMATIONS PLASTIQUES */
/*  ALPHA	: TENSEUR DES CONTRAINTES CINEMATIQUES */
/*  DROT		: INCREMENT DE MATRICE ROTATION (WINGET) */
/*  DFGRD0	: GRADIENT DES DEFORMATIONS EN DEBUT D'INCREMENT */
/*  DFGRD1	: GRADIENT DES DEFORMATIONS EN FIN D'INCREMENT */
/*  XLAMBDA	: DEFORMATION PLASTIQUE CUMULEE */
/*  PR		: VARIABLE D'ECROUISSAGE ISOTROPE */
/*  ENDO		: TENSEUR D'ENDOMMAGEMENT */
/*  DLAN		: INCREMENT DE DEFORMATION PLASTIQUE CUMULEE PRECEDENT */
/*  TIME		: STEP TIME (#1) AND TOTAL TIME (#2) AT BEGINNING */
/* 		  OF THE INCREMENT */
/*  DTIME	: TIME INCREMENT */
/*  TEMP		: TEMPERATURE AT BEGINNING OF INCREMENT */
/*  DTEMP	: TEMPERATURE INCREMENT */
/*  T0		: TEMPERATURE INITIALE */
/*  RHO		: DENSITE VOLUMIQUE */
/*  C		: COEFFICIENT ADIABATIQUE */
/*  E		: MODULE D'YOUNG */
/*  G		: MODULE DE CISAILLEMENT */
/*  XK		: COEFFICIENT DE LAME */
/*  XNU		: COEFFICIENT DE POISSON */
/*  YIELD	: LIMITE D'ELASTICITE */
/*  V		: CONSTANTE VISCOPLASTIQUE */
/*  PF		: CONSTANTE VISCOPLASTIQUE */
/*  RD		: CONSTANTE D'ECROUISSAGE ISOTROPE */
/*  H		: CONSTANTE D'ECROUISSAGE ISOTROPE */
/*  RS		: CONSTANTE D'ECROUISSAGE ISOTROPE (RECOVERY) */
/*  PRD		: CONSTANTE D'ECROUISSAGE CINEMATIQUE */
/*  PH		: CONSTANTE D'ECROUISSAGE CINEMATIQUE */
/*  PRS		: CONSTANTE D'ECROUISSAGE CINEMATIQUE (RECOVERY) */
/*  XIP		: INITIAL POROSITY */
/*  DC		: VALEUR DE L'ENDOMMAGEMENT CRITIQUE --> RUPTURE */
/*  XMS		: DAMAGE RATE SENSIVITY */
/*  TETA		: TETA-METHODE */
/*  EPSF		: PRECISION DE LA FONCTION DE CHARGE F */
/*  EPSX		: PRECISION DE LA FONCTION D'ENDOMMAGEMENT OMEGA-X */
/*  EPSR		: PRECISION DE LA FONCTION D'ENDOMMAGEMENT OMEGA-R */
/*  EPSD		: PRECISION DE LA FONCTION TENSORIELLE D'ENDOMMAGEMENT PHI */
/*  NCD		: NCD=1 --> COUPLAGE ELASTOPLASTICITE-ENDOMMAGEMENT */
/*  		  NCD=0 --> ELASTOPLASTICITE SANS ENDOMMAGEMENT */
/*  NPLAN	: NPLAN=1 --> DEFORMATIONS PLANES */
/* 		  NPLAN=2 --> CONTRAINTES PLANES */
/*  NATUR	: NATUR=0 --> PETITES DEFORMATIONS */
/* 		  NATUR=1 --> GRANDES DEFORMATIONS */
/*  IROT		: IROT=0 --> PAS DE ROTATION */
/* 		  IROT=1 --> AVEC ROTATION */
/* 		: IROT=2 --> AVEC DERIVEE COROTATIONNELLE DE JAUMANN */
/* 		  IROT=3 --> AVEC ROTATION D'ABAQUS (ROTSIG) */
/*  NVI		: NOMBRE DE VARIABLE INTERNES STOCKEES (R,X,J2(X),Y,...) */
/*  M		: COMMENTAIRES (PARAMETRE D'IMPRESSION) */
/*  MPI		: UNITE DE FICHIER DE SORTIE */
/* ----------------------------------------------------------------------* */
/* ======================================================================* */
/* ----------------------------------------------------------------------* */
/* .1-----  Precision */
/* .2-----  Parameter,Integer */
/* ------- */
/* .3-----  Dimension */
/* ------- */
/* .4-----  Real,Complex,Double precision,Logical,Character */
/* .5-----  Common */
/* .6-----  Equivalence */
/* .7-----  Save,Data */
    /* Parameter adjustments */
    --props;
    --statev;
    drot -= 4;
    dfgrd0 -= 4;
    dfgrd1 -= 4;
    --time;
    --alpha;
    --eplas;
    --dstran;
    --stran;
    --stress;

    /* Function Body */
/* .8-----  Definition de fonctions */
/* ----------------------------------------------------------------------* */
/* ======================================================================* */
/* ----------------------------------------------------------------------* */
/* -------  Eclatement du tableau PROPS(NPROPS) */
    *m = (integer) props[42];
    t0 = props[1];
/*      T0     = TEMP+DTEMP */
    rho = props[2];
    c__ = props[3];
    *e = props[4];
    *xnu = props[5];
    *xip = props[28];
    *dc = props[29];
    *xms = props[30];
    *teta = props[31];
    *epsf = props[32];
    *epsx = props[33];
    *epsr = props[34];
    *epsd = props[35];
    *vnewdt = props[36];
    *ncd = (integer) props[37];
    *nplan = (integer) props[38];
    *natur = (integer) props[39];
    *irot = (integer) props[40];
    *nvi = (integer) props[41];
/* ----------- Constantes du modele BCJ (c1 --> c22) */
    for (i__ = 1; i__ <= 22; ++i__) {
	cts[i__ - 1] = props[i__ + 5];
    }
/* -------  FIN Eclatement du tableau PROPS(NPROPS) */
/* ----------------------------------------------------------------------* */
/* -------  Calcul du module de cisaillement G et de K */

    *g = udemi * *e / (un + *xnu);
    tunmdxnu = trois * (un - deux * *xnu);
    *xk = *e / tunmdxnu;

/* -------  FIN  Calcul du module de cisaillement G et de K */
/* ----------------------------------------------------------------------* */
/* -------  Calcul de la constante de nucleation */

    utier = un / trois;
    dtier = deux * utier;
    rttd = sqrt(trois * udemi);

/* -------  FIN  Calcul de la constante de nucleation */
/* ----------------------------------------------------------------------* */
/* -------  Calcul des constantes elastoplastiques */

    ut0 = un / t0;
    ya = cts[2] / (cts[20] + exp(-cts[3] * ut0));
    if (cts[18] == 0.0) { /* BCJ convention */
    	yb = 1.0;
    } else {
    	yb = udemi*(un + tanh(cts[18] * (cts[19] - t0)));
    }
/* ------- */
    *v = cts[0] * exp(-cts[1] * ut0);
    *yield = ya * yb;
/*      YIELD = CTS(3)*DEXP(CTS(4)*UT0) */
    *pf = *dtime * cts[4] * exp(-cts[5] * ut0);
    *hk = dtier * (cts[8] - cts[9] * t0);
    *rdk = rttd * cts[6] * exp(-cts[7] * ut0);
    *rsk = *dtime * rttd * cts[10] * exp(-cts[11] * ut0);
    *hi = cts[14] - cts[15] * t0;
    *rdi = cts[12] * exp(-cts[13] * ut0);
    *rsi = *dtime * cts[16] * exp(-cts[17] * ut0);
/* ------- */
/*      HK    = ZERO */
/*      RDK   = ZERO */
/*      RSK   = ZERO */
/*      HI    = ZERO */
/*      RDI   = ZERO */
/*      RSI   = ZERO */
/*      V     = ZERO */

/* -------  FIN  Calcul des constantes elastoplastiques */
/* ----------------------------------------------------------------------* */
/* -------  Eclatement du tableau STATEV(NSTATV) */
/* ----------  Deformations totales, plastiques et elastiques */
/*            Ecrouissage cinematique */
    if (*natur == 0) {
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    eplas[i__] = statev[i__];
	    alpha[i__] = statev[i__ + *ntens];
	}
    } else {
	if (*irot == 1) {
	    rotsig_(&statev[1], &drot[4], &eplas[1], &c__1, ndi, nshr);
	    rotsig_(&statev[*ntens + 1], &drot[4], &alpha[1], &c__1, ndi, 
		    nshr);
	} else {
	    i__1 = *ntens;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		stressi[i__ - 1] = stress[i__];
	    }
	    transpose_(&c__3, &drot[4], tdrot);
	    dpolr_(rot0, &dfgrd0[4], m, mpi);
	    dpolr_(rot1, &dfgrd1[4], m, mpi);
	    transpose_(&c__3, rot0, trot0);
	    prdmat_(&c__3, rot1, trot0, dr);
	    for (i__ = 1; i__ <= 3; ++i__) {
		for (j = 1; j <= 3; ++j) {
		    dr[i__ + j * 3 - 4] = udemi * dr[i__ + j * 3 - 4];
		}
		dr[i__ + i__ * 3 - 4] = udemi + dr[i__ + i__ * 3 - 4];
	    }
	    prdmat_(&c__3, dr, tdrot, drotgn);
	    rotsig_(stressi, drotgn, &stress[1], &c__1, ndi, nshr);
	    rotsig_(&statev[1], dr, &eplas[1], &c__1, ndi, nshr);
	    rotsig_(&statev[*ntens + 1], dr, &alpha[1], &c__1, ndi, nshr);
	    rotsig_(&statev[*nvi + 4 + *ntens * 3], dr, &stran[1], &c__2, ndi,
		     nshr);
	    i__1 = *ntens;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		statev[i__ + 3 + *nvi + *ntens * 3] = stran[i__] + dstran[i__]
			;
	    }
	}
    }
/* ---------  Deformation plastique cumulee */
    *xlambda = statev[(*ntens << 1) + 1];
/* ---------  Ecrouissage isotrope */
    *pr = statev[(*ntens << 1) + 2];
/* ---------  Damage */

	/* initialize damage even if it is not evolving */
	*phi = statev[(*ntens << 1) + 3];
	if (*phi == zero) {
	    *phi = *xip;
	} else {
	    *xip = zero;
	}
	
	/* freeze damage at maximum */
	if (*phi >= *dc) {
		*ncd = 0;
		*phi = *dc;
	}

#if 0
    if (*ncd == 0) {
	*phi = zero;
	*xip = zero;
    } else {
	*phi = statev[(*ntens << 1) + 3];
/* ------------  Initial Damage as Initial Porosity */
	if (*phi == zero) {
	    *phi = *xip;
	} else {
	    *xip = zero;
	}
    }
#endif

/* ---------  Increment de deformation plastique cumulee precedent */
    *dlan = statev[(*ntens << 1) + 4];
/* -------   FIN Eclatement du tableau STATEV(NSTATV) */
/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* -------  Comment Check Code */
    if (*m >= 4) {
	io___654.ciunit = *mpi;
	s_wsle(&io___654);
	do_lio(&c__9, &c__1, "RECOV EPLAS", (ftnlen)11);
	e_wsle();
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___655.ciunit = *mpi;
	    s_wsfe(&io___655);
	    do_fio(&c__1, (char *)&eplas[i__], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	io___656.ciunit = *mpi;
	s_wsle(&io___656);
	do_lio(&c__9, &c__1, "RECOV ALPHA", (ftnlen)11);
	e_wsle();
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___657.ciunit = *mpi;
	    s_wsfe(&io___657);
	    do_fio(&c__1, (char *)&alpha[i__], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	io___658.ciunit = *mpi;
	s_wsfe(&io___658);
	do_fio(&c__1, (char *)&(*xlambda), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*phi), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*pr), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*dlan), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___659.ciunit = *mpi;
	s_wsfe(&io___659);
	do_fio(&c__1, (char *)&(*v), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*yield), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*pf), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*rdk), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*hk), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*rsk), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*rdi), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*hi), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*rsi), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
/* -------  END Comment Check Code */
/* ----------------------------------------------------------------------* */
/* ---------   FORMATS  -------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* L401: */
/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* ======================================================================= */
    return 0;
} /* recov_ */


/* Subroutine */ int solelas_(doublereal *stress, doublereal *statev, 
	doublereal *ddsdde, doublereal *stresselas, doublereal *alphaelas, 
	doublereal *prelas, doublereal *phi, doublereal *xk, doublereal *g, 
	doublereal *xip, integer *ncd, integer *nstatv, integer *ntens, 
	integer *ndi, integer *nvi, integer *m, integer *mpi)
{
    /* Initialized data */

    static doublereal un = 1.;
    static doublereal deux = 2.;
    static doublereal trois = 3.;

    /* Format strings */
    static char fmt_101[] = "(1x,6(d20.14,2x))";

    /* System generated locals */
    integer ddsdde_dim1, ddsdde_offset, i__1, i__2;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void), s_wsfe(cilist *), do_fio(integer *, char *, ftnlen),
	     e_wsfe(void);

    /* Local variables */
    static integer i__, j;
    static doublereal d2g, xkg, unmd, unmdg;
    extern /* Subroutine */ int affect_(integer *, doublereal *, doublereal *)
	    ;
    static doublereal unmdd2g, unmdxkg;

    /* Fortran I/O blocks */
    static cilist io___672 = { 0, 0, 0, 0, 0 };
    static cilist io___673 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___674 = { 0, 0, 0, 0, 0 };
    static cilist io___675 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___676 = { 0, 0, 0, 0, 0 };
    static cilist io___677 = { 0, 0, 0, fmt_101, 0 };


/* ====================================================================== */
/* ---------------------------------------------------------------------- */
/* ------------------------  SOLUTION ELASTIQUE  ------------------------ */
/* ---------------------------------------------------------------------- */
/* ====================================================================== */

/* ENTREES */
/* ------- */
/*  ENDO		: TENSEUR D'ENDOMMAGEMENT */
/*  STRESSELAS	: TENSEUR DES CONTRAINTES AVEC PREDICTION ELASTIQUE */
/*  PHI		: DAMAGE */
/*  XK		: COEFFICIENT DE LAME */
/*  G		: MODULE DE CISAILLEMENT */
/*  XIP		: INITIAL POROSITY */
/*  NCD		: NCD =1 --> COUPLED WITH DAMAGE */
/*  NSTATV	: TAILLE DU TABLEAU STATEV */
/*  NTENS	: TAILLE DU TABLEAU DES CONTRAINTES OU DEFORMATIONS */
/*  NDI		: NOMBRE DE COMPOSANTES DIRECTES D'UN TENSEUR */
/*  NVI		: NOMBRE DE VARIABLE INTERNES STOCKEES (EE,R,X,J2(X),Y,...) */
/*  M		: COMMENTAIRES (PARAMETRE D'IMPRESSION) */
/*  MPI		: UNITE DE FICHIER DE SORTIE */
/* -------- */
/* LOCALES */
/* ------- */
/*  DMT		: TENSEUR D'EFFET D'ENDOMMAGEMENT (M) */
/*  DVDM	 	: PRODUIT Idev * Mt */
/* -------- */
/* SORTIES */
/* ------- */
/*  STRESS	: TENSEUR DES CONTRAINTES */
/*  DDSDDE 	: MATRICE JACOBIENNE ELASTIQUE */
/* ----------------------------------------------------------------------* */
/* ======================================================================= */
/* ----------------------------------------------------------------------* */
/* .1-----  Implicit, External */
/* .2-----  Parameter,Integer */
/* .3-----  Dimension */
/* ------- */
/* .4-----  Real,Complex,Double precision,Logical,Character */
/* .5-----  Common */
/* .6-----  Equivalence */
/* .7-----  Save,Data */
    /* Parameter adjustments */
    --statev;
    --alphaelas;
    --stresselas;
    ddsdde_dim1 = *ntens;
    ddsdde_offset = 1 + ddsdde_dim1;
    ddsdde -= ddsdde_offset;
    --stress;

    /* Function Body */
/* .8-----  Definition de fonctions */
/* ----------------------------------------------------------------------* */
/* ======================================================================= */
/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* -------  Store of Cauchy Stress */

    affect_(ntens, &stresselas[1], &stress[1]);

/* -------  END  Store of Cauchy Stress */
/* ----------------------------------------------------------------------* */
/* -------  Store of Kinematic Hardening */

    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
	statev[i__ + *ntens] = alphaelas[i__];
    }

/* -------  END  Store of Kinematic Hardening */
/* ----------------------------------------------------------------------* */
/* -------  Store of Isotropic Hardening */

    statev[(*ntens << 1) + 2] = *prelas;

/* -------  END  Store of Isotropic Hardening */
/* ----------------------------------------------------------------------* */
/* -------  Stockage des variables internes */

    if (*ncd == 1 || *ncd == 2) {
	statev[(*ntens << 1) + 3] = *phi;
    }

/* -------  END  Stockage des variables internes */
/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* -------  Elastic Jacobian */

    unmd = un - *phi;
    d2g = deux * *g;
    xkg = *xk - d2g / trois;
    unmdg = unmd * *g;
    unmdd2g = unmd * d2g;
    unmdxkg = unmd * xkg;
/* ------- */
    i__1 = *ndi;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *ndi;
	for (j = 1; j <= i__2; ++j) {
	    ddsdde[i__ + j * ddsdde_dim1] = unmdxkg;
	}
	ddsdde[i__ + i__ * ddsdde_dim1] += unmdd2g;
    }
    i__1 = *ntens;
    for (i__ = *ndi + 1; i__ <= i__1; ++i__) {
	ddsdde[i__ + i__ * ddsdde_dim1] = unmdg;
    }
    if (*m >= 0) {
    io___672.ciunit = *mpi;
    s_wsle(&io___672);
    do_lio(&c__9, &c__1, "SOLELAS DDSDDE", (ftnlen)14);
    e_wsle();
    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
	io___673.ciunit = *mpi;
	s_wsfe(&io___673);
	i__2 = *ntens;
	for (j = 1; j <= i__2; ++j) {
	    do_fio(&c__1, (char *)&ddsdde[i__ + j * ddsdde_dim1], (ftnlen)
		    sizeof(doublereal));
	}
	e_wsfe();
    }
	}

/* -------  END  Elastic Jacobian */
/* ----------------------------------------------------------------------* */
/* -------  Comment Check Code */
    if (*m >= 4) {
	io___674.ciunit = *mpi;
	s_wsle(&io___674);
	do_lio(&c__9, &c__1, "SOLELAS STRESS", (ftnlen)14);
	e_wsle();
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___675.ciunit = *mpi;
	    s_wsfe(&io___675);
	    do_fio(&c__1, (char *)&stress[i__], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	io___676.ciunit = *mpi;
	s_wsle(&io___676);
	do_lio(&c__9, &c__1, "SOLELAS DDSDDE", (ftnlen)14);
	e_wsle();
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___677.ciunit = *mpi;
	    s_wsfe(&io___677);
	    i__2 = *ntens;
	    for (j = 1; j <= i__2; ++j) {
		do_fio(&c__1, (char *)&ddsdde[i__ + j * ddsdde_dim1], (ftnlen)
			sizeof(doublereal));
	    }
	    e_wsfe();
	}
    }
/* -------  END Comment Check Code */
/* ----------------------------------------------------------------------* */
/* ---------   FORMATS  -------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* ======================================================================= */
    return 0;
} /* solelas_ */


/* Subroutine */ int solsyst_(doublereal *dlambda, doublereal *phi, 
	doublereal *wx, doublereal *wr, doublereal *ff, doublereal *fx, 
	doublereal *fr, doublereal *fd, doublereal *dfdl, doublereal *dfdx, 
	doublereal *dfdr, doublereal *dfdd, doublereal *dxdl, doublereal *
	dxdx, doublereal *dxdr, doublereal *dxdd, doublereal *drdl, 
	doublereal *drdr, doublereal *dhdl, doublereal *dhdx, doublereal *
	dhdr, doublereal *dhdd, integer *ncd, integer *ntens, integer *ndi, 
	integer *nshr, integer *m, integer *mpi)
{
    /* Initialized data */

    static doublereal zero = 0.;
    static doublereal one = 1.;

    /* Format strings */
    static char fmt_102[] = "(1x,a20,6(d20.14,2x))";

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static doublereal aa, ab, bb, ad, bd, dd, al, bl, dl, dr, dx, abl, dfrr, 
	    dhrr, dfxx, dhxx, dxrr, udrdr, udxdx;

	/* rescaling of update vector */
	doublereal scale = 1.0, new_phi;

    /* Fortran I/O blocks */
    static cilist io___699 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___700 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___701 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___702 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___703 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___704 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___705 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___706 = { 0, 0, 0, fmt_102, 0 };


/* ====================================================================== */
/* ---------------------------------------------------------------------- */
/* -----------  CALCUL DE LA SOLUTION CLAMBDA, CWX, CWR et CD  ---------- */
/* -----------------------------  DU SYSTEME  --------------------------- */
/* ----------  PUIS REACTUALISATION DE DLAMBDA, WX, WR et ENDO  --------- */
/* ---------------------------------------------------------------------- */
/* ====================================================================== */
/* ------------- */
/*  EN ENTREE : */
/* ------------- */
/*  FF		: VALEUR DE LA FONCTION F */
/*  FX		: VALEUR DE LA FONCTION GX */
/*  FR		: VALEUR DE LA FONCTION GR */
/*  FH		: VALEUR DE LA FONCTION H */
/*  DFDL		: DERIVEE PARTIELLE DE F PAR RAPPORT A DLAMBDA */
/*  DFDX		: DERIVEE PARTIELLE DE F PAR RAPPORT A WX */
/*  DFDR		: DERIVEE PARTIELLE DE F PAR RAPPORT A WR */
/*  DFDD		: DERIVEE PARTIELLE DE F PAR RAPPORT A ENDO */
/*  DXDL		: DERIVEE PARTIELLE DE GX PAR RAPPORT A DLAMBDA */
/*  DXDX		: DERIVEE PARTIELLE DE GX PAR RAPPORT A WX */
/*  DXDD		: DERIVEE PARTIELLE DE GX PAR RAPPORT A ENDO */
/*  DRDL		: DERIVEE PARTIELLE DE GR PAR RAPPORT A DLAMBDA */
/*  DRDR		: DERIVEE PARTIELLE DE GR PAR RAPPORT A WR */
/*  DHDL		: DERIVEE PARTIELLE DE H PAR RAPPORT A DLAMBDA */
/*  DHDX		: DERIVEE PARTIELLE DE H PAR RAPPORT A WX */
/*  DHDR		: DERIVEE PARTIELLE DE H PAR RAPPORT A WR */
/*  DHDD		: DERIVEE PARTIELLE DE H PAR RAPPORT A ENDO */
/*  NCD		: NCD=1 --> COUPLAGE ELASTOPLASTICITE-ENDOMMAGEMENT */
/*  		  NCD=0 --> ELASTOPLASTICITE SANS ENDOMMAGEMENT */
/*  NTENS	: TAILLE DU TABLEAU DES CONTRAINTES OU DEFORMATIONS */
/*  NDI		: NOMBRE DE COMPOSANTES DES CONTRAINTES DIRECTES */
/*  NSHR		: NOMBRE DE COMPOSANTES DES CONTRAINTES TANGENTIELLES */
/*  M		: COMMENTAIRES */
/*  MPI		: UNITE DE FICHIER DE SORTIE */
/* ------------- */
/*  EN SORTIE : */
/* ------------- */
/*  DLAMBDA	: INCREMENT DEFORMATION PLASTIQUE CUMULEE */
/*  WR		: ISOTROPIC DYNAMIC RECOVERY TERM: WR = DSINH[QS*R] */
/*  WX		: KINEMATIC DYNAMIC RECOVERY TERM: WX = |DMI:X| */
/*  ENDO		: TENSEUR D'ENDOMMAGEMENT */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
/* ----------------------------------------------------------------------- */
/* .1-----  Implicit, External */
/* .2-----  Parameter,Integer */
/* .3-----  Dimension */
/* .4-----  Real,Complex,Double precision,Logical,Character */
/* .5-----  Common */
/* .6-----  Equivalence */
/* .7-----  Save,Data */
/* .8-----  Definition de fonctions */
/* ---------------------------------------------------------------------* */
/* ====================================================================== */
/* ---------------------------------------------------------------------* */
/* -------  Initialisation */


/* -------  END Initialisation */
/* ---------------------------------------------------------------------* */
/* -------  Solution DD=cd (Damage) and DL=clambda (Plastic equivalent strain) */
    udxdx = one / *dxdx;
    udrdr = one / *drdr;
    dfxx = *dfdx * udxdx;
    dhxx = *dhdx * udxdx;
    dfrr = *dfdr * udrdr;
    dxrr = *dxdr * udrdr;
    dhrr = *dhdr * udrdr;
    ab = *fx - dxrr * *fr;
    abl = *dxdl - dxrr * *drdl;
/* ------- */
    aa = *ff - dfxx * ab - dfrr * *fr;
    al = *dfdl - dfxx * abl - dfrr * *drdl;
    ad = *dfdd - dfxx * *dxdd;
    bb = *fd - dhxx * ab - dhrr * *fr;
    bl = *dhdl - dhxx * abl - dhrr * *drdl;
    bd = *dhdd - dhxx * *dxdd;
/* ------------------------------------------* */
/* -------  Solution DD=cd (Damage) */
    if (*ncd == 0) {
	dd = zero;
    } else {
	dd = (aa * bl - bb * al) / (al * bd - ad * bl);
    }
    
/* ------------------------------------------* */
/* -------  Solution Dl=clambda (Plastic Equivalent Strain Increment) */
    dl = -(aa + ad * dd) / al;
/* ------------------------------------------* */
/* -------  Solutions DX=cWd (norm of Alpha) and DR=cWr (Kappa**2) */
    dx = -udxdx * (ab + abl * dl + *dxdd * dd);
    dr = -udrdr * (*fr + *drdl * dl);
/* ------------------------------------------* */
/* ---------------------------------------------------------------------* */
/* -------  Update of the solution clambda, cWx, cWr et cd */

	/* bigger than 5% change */
	if (fabs(dd) > fabs(0.05)) {
		scale = 0.05/fabs(dd);
		dl *= scale;
		dx *= scale;
		dr *= scale;
		dd *= scale;	
	}

	/* damage exceeds 1.0 */
	new_phi = *phi + dd;
	if (new_phi > 1.0) {
		scale = 0.5*(1.0 - *phi)/dd; /* half way to max */
		dl *= scale;
		dx *= scale;
		dr *= scale;
		dd *= scale;		
	}

	/* update solution */
    *dlambda += dl;
    *wx += dx;
    *wr += dr;
    *phi += dd;

/* ----------- END  Update of the solution clambda, cWx, cWr et cd */
/* ---------------------------------------------------------------------* */
/* ---------------------------------------------------------------------* */
/* -------  Comment Check Code */

    if (*m >= 4) {
	io___699.ciunit = *mpi;
	s_wsfe(&io___699);
	do_fio(&c__1, "SOLSYST DL = ", (ftnlen)13);
	do_fio(&c__1, (char *)&dl, (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___700.ciunit = *mpi;
	s_wsfe(&io___700);
	do_fio(&c__1, "SOLSYST DX = ", (ftnlen)13);
	do_fio(&c__1, (char *)&dx, (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___701.ciunit = *mpi;
	s_wsfe(&io___701);
	do_fio(&c__1, "SOLSYST DR = ", (ftnlen)13);
	do_fio(&c__1, (char *)&dr, (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___702.ciunit = *mpi;
	s_wsfe(&io___702);
	do_fio(&c__1, "SOLSYST DD = ", (ftnlen)13);
	do_fio(&c__1, (char *)&dd, (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___703.ciunit = *mpi;
	s_wsfe(&io___703);
	do_fio(&c__1, "SOLSYST DLAMBDA = ", (ftnlen)18);
	do_fio(&c__1, (char *)&(*dlambda), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___704.ciunit = *mpi;
	s_wsfe(&io___704);
	do_fio(&c__1, "SOLSYST WX      = ", (ftnlen)18);
	do_fio(&c__1, (char *)&(*wx), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___705.ciunit = *mpi;
	s_wsfe(&io___705);
	do_fio(&c__1, "SOLSYST WR      = ", (ftnlen)18);
	do_fio(&c__1, (char *)&(*wr), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___706.ciunit = *mpi;
	s_wsfe(&io___706);
	do_fio(&c__1, "SOLSYST PHI     = ", (ftnlen)18);
	do_fio(&c__1, (char *)&(*phi), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/* -------  END Comment Check Code */
/* ----------------------------------------------------------------------* */
/* ---------   FORMATS  -------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* L101: */
/* ----------------------------------------------------------------------* */
/* ======================================================================= */
    return 0;
} /* solsyst_ */


/* Subroutine */ int solution_(doublereal *stress, doublereal *eplas, 
	doublereal *alpha, doublereal *xlambda, doublereal *pr, doublereal *
	strant, doublereal *phi, doublereal *zd, doublereal *tn, doublereal *
	stres1d, doublereal *eplas1, doublereal *alpha1, doublereal *pr1, 
	doublereal *dlambda, doublereal *dlan, doublereal *tdl, doublereal *
	unmd, doublereal *untrwx, doublereal *untrwr, doublereal *rtdl, 
	doublereal *eyv, doublereal *g, doublereal *xk, doublereal *hk, 
	doublereal *hi, doublereal *xms, doublereal *teta, doublereal *dtime, 
	doublereal *stvit, integer *nitmax, integer *nstvi, integer *ntens, 
	integer *ndi, integer *nshr, integer *ncd, integer *m, integer *mpi)
{
    /* Initialized data */

    static doublereal zero = 0.;
    static doublereal one = 1.;
    static doublereal two = 2.;
    static doublereal three = 3.;

    /* Format strings */
    static char fmt_101[] = "(1x,6(d20.14,2x))";
    static char fmt_102[] = "(1x,a20,6(d20.14,2x))";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void), s_wsfe(cilist *), do_fio(integer *, char *, ftnlen),
	     e_wsfe(void);

    /* Local variables */
    static doublereal effstran, twothird;
    static integer i__;
    static doublereal x[6], eh, xh, dehk;
    extern /* Subroutine */ int aset_(doublereal *, doublereal *, integer *);
    static doublereal dltn[6];
    extern /* Subroutine */ int sinv_(doublereal *, doublereal *, doublereal *
	    , integer *, integer *);
    static doublereal xeqv;
    extern /* Subroutine */ int trace_(doublereal *, doublereal *, integer *, 
	    integer *);
    static doublereal d2gtdl;
    extern /* Subroutine */ int pdtsca_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *);

    /* Fortran I/O blocks */
    static cilist io___721 = { 0, 0, 0, 0, 0 };
    static cilist io___722 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___723 = { 0, 0, 0, 0, 0 };
    static cilist io___724 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___725 = { 0, 0, 0, 0, 0 };
    static cilist io___726 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___727 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___728 = { 0, 0, 0, fmt_102, 0 };


/* ====================================================================== */
/* ---------------------------------------------------------------------- */
/* --------------  CALCUL DE LA EPLAS(N+1), ALPHA(N+1),  ---------------- */
/* ---------------------  PR(N+1) et STRESS(N+1)  ----------------------- */
/* ---------------------------------------------------------------------- */
/* ====================================================================== */
/* ------------- */
/*  EN ENTREE : */
/* ------------- */
/*  ZD		: NORMAL DIRECTION STRESS TENSOR */
/*  TN		: NORMAL TENSOR */
/*  STRANT	: TOTAL STRAIN TENSOR */
/*  STRES1D	: KNOWN DEVIATORIC STRESS TENSOR */
/*  EPLAS1	: KNOWN PLASTIC STRAIN TENSOR */
/*  ALPHA1	: KNOWN KINEMATIC HARDENING TENSOR */
/*  PR1		: KNOWN ISOTROPIC HARDENING VARIABLE */
/*  DLAMBDA	: PLASTIC EQUIVALENT STRAIN RATE */
/*  DLAN		: PLASTIC EQUIVALENT STRAIN RATE AT PRECENDENT STEP */
/*  DM		: DAMAGE EFFECT TENSOR Mt */
/*  DMI		: DAMAGE EFFECT TENSOR M(-t) */
/*  DN		: ORIENTATION VOID OPERATOR */
/*  QD		: VOID COALESCENCE OPERATOR */
/*  TDL		: TETA*DLAMBDA */
/*  UDHM		: VALUE 1/(1-Dh)**m */
/*  UNTRWX	: 1+TETA*(PRD*DLAMBDA+PRS)*WX */
/*  UNTRWR	: 1+TETA*(RD*DLAMBDA+RS)*WR */
/*  RTDL		: R1+TETA*DLAMBDA */
/*  Y		: ENERGY RELEASE RATE Y = Ye + Yan */
/*  YT		: TRIAXIALITY YT = <Sh>/Seq */
/*  YPS		: FUNCTION YPS = (Y/S)**s */
/*  EYV		: FUNCTION Exp[3/2*<Sh>/Seq] */
/*  G		: SHEAR MODULUS */
/*  XK		: LAME COEFFICIENT */
/*  H		: ISOTROPIC HARDENING MODULUS */
/*  PH		: KINEMATIC HARDENING MODULUS */
/*  XIP		: INITIAL POROSITY */
/*  XDN		: NUCLEATION CONSTANTS XLD*XFN/(XKIC*XLF**(1/3)) */
/*  CTE		: CONSTANT OF TENSION */
/*  CCO		: CONSTANT OF COMPRESSION */
/*  CTO		: CONSTANT OF TORSION */
/*  AV		: VOID GROWTH CONSTANT */
/*  VP		: VOID GROWTH TEST (VC = 1 if Dh > Vc) */
/*  VT		: POSITIVE PART <Dh-Vc>  (Dh=1/3*Trace(Phi)) */
/*  AC		: COALESCENCE CONSTANT */
/*  TETA		: TETA-METHODE */
/*  DTIME	: TIME INCREMENT */
/*  NITMAX	: MAXIMUM ITERATION NUMBER */
/*  NSTVI	: INTERNAL VARIABLE ARRAY STVIT SIZE */
/*  NDI		: DIAGONAL COMPONENT NUMBER */
/*  NSHR		: SHEAR COMPONENT NUMBER */
/*  NTENS	: STRESS TENSOR ARRAY SIZE (NDI + NSHR). */
/*  M		: COMMENTS */
/*  MPI		: OUTPUT UNIT FILE */
/* ------------- */
/*  EN SORTIE : */
/* ------------- */
/*  STRESS	: CAUCHY STRESS TENSOR AT TIME N+1 */
/*  EPLAS	: PLASTIC STRAIN TENSOR AT TIME N+1 */
/*  EELAS	: ELASTIC STRAIN TENSOR AT TIME N+1 */
/*  ALPHA	: KINEMATIC HARDENING TENSOR AT TIME N+1 */
/*  XLAMBDA	: PLASTIC EQUIVALENT STRAIN A N+1 */
/*  DLAMBDA	: PLASTIC EQUIVALENT STRAIN INCREMENT */
/*  PR		: ISOTROPIC HARDENING TENSOR AT TIME N+1 */
/*  STVIT	: INTERNAL VARIABLE ARRAY: R, J2(X) et X */
/*  DAMAGE	: RUPTURE INDICATOR */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
/* ----------------------------------------------------------------------- */
/* .1-----  Implicit, External */
/* .2-----  Parameter,Integer */
/* .3-----  Dimension */
/* ------- */
/* .4-----  Real,Complex,Double precision,Logical,Character */
/* .5-----  Common */
/* .6-----  Equivalence */
/* .7-----  Save,Data */
    /* Parameter adjustments */
    --stvit;
    --alpha1;
    --eplas1;
    --stres1d;
    --tn;
    --zd;
    --strant;
    --alpha;
    --eplas;
    --stress;

    /* Function Body */
/* .8-----  Functions Definitions */
/* ---------------------------------------------------------------------* */
/* ====================================================================== */
/* ---------------------------------------------------------------------* */
/* -------  Initialization */

    twothird = two / three;

/* -------  END Initialization */
/* ---------------------------------------------------------------------* */
/* -------  Plastic Strain EPLAS */

    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dltn[i__ - 1] = *tdl * tn[i__];
	eplas[i__] = eplas1[i__] + dltn[i__ - 1];
    }

/* -------  END  Plastic Strain EPLAS */
/* ---------------------------------------------------------------------* */
/* -------  Kinematic Hardening Tensors ALPHA and X */

    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
	alpha[i__] = (alpha1[i__] + dltn[i__ - 1]) / *untrwx;
	x[i__ - 1] = *hk * alpha[i__];
    }
    if (*hk == zero) {
	aset_(&alpha[1], &zero, ntens);
    }

/* -------  END  Kinematic Hardening Tensors ALPHA and X */
/* ---------------------------------------------------------------------* */
/* -------  Plastic Equivalent Strain Rate DLAMBDA */

/*      DLAMBDA = (ONE-TETA)*DLAN*DTIME+TETA*DLAMBDA */

/* -------  END Plastic Equivalent Strain Rate DLAMBDA DLAMBDA */
/* ---------------------------------------------------------------------* */
/* -------  Plastic Equivalent Strain XLAMBDA */

    *xlambda = *xlambda + (one - *teta) * *dlan * *dtime + *teta * *dlambda;

/* -------  END Plastic Equivalent Strain XLAMBDA */
/* ---------------------------------------------------------------------* */
/* -------  Isotropic Hardening PR */

    *pr = *rtdl / *untrwr;
    if (*hi == zero) {
	*pr = zero;
    }

/* -------  END Isotropic Hardening PR */
/* ---------------------------------------------------------------------* */
/* -------  Cauchy Stress Tensor */

    trace_(&strant[1], &eh, ndi, ntens);
/* ------- */
    d2gtdl = two * *g * *unmd * *tdl;
    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
	stress[i__] = *unmd * stres1d[i__] - d2gtdl * tn[i__];
    }
/* ------- */
    dehk = *xk * *unmd * eh;
    i__1 = *ndi;
    for (i__ = 1; i__ <= i__1; ++i__) {
	stress[i__] += dehk;
    }

/* -------  END Cauchy Stress Tensor */
/* ---------------------------------------------------------------------* */
/* ---------------------------------------------------------------------* */
/* -------  Storage of R, J2(X), X, Y, Ye and Damage Tensors in STVIT array. */

/* ---------  Isotropic and Kinematic Hardening Stresses */
    sinv_(x, &xh, &xeqv, ndi, nshr);
    stvit[1] = *hi * *pr;
    stvit[2] = xeqv;
    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
	stvit[i__ + 2] = x[i__ - 1];
    }
/* ---------  Triaxiality <p>/Seqv */
    stvit[*ntens + 3] = *eyv;
/* ---------  Reached Iteration Number */
    stvit[*ntens + 4] = (doublereal) (*nitmax);
/* ---------  Effective Strain */
    pdtsca_(ndi, nshr, &strant[1], &strant[1], &effstran);
    stvit[*ntens + 5] = sqrt(twothird * effstran);

/* -------  END Storage of R, J2(X), X, Y, Ye and Damage Tensors in STVIT array. */
/* ---------------------------------------------------------------------* */
/* ---------------------------------------------------------------------* */
/* -------  Comment Check Code */

    if (*m >= 4) {
	io___721.ciunit = *mpi;
	s_wsle(&io___721);
	do_lio(&c__9, &c__1, "SOLUTION STRESS", (ftnlen)15);
	e_wsle();
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___722.ciunit = *mpi;
	    s_wsfe(&io___722);
	    do_fio(&c__1, (char *)&stress[i__], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	io___723.ciunit = *mpi;
	s_wsle(&io___723);
	do_lio(&c__9, &c__1, "SOLUTION EPLAS", (ftnlen)14);
	e_wsle();
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___724.ciunit = *mpi;
	    s_wsfe(&io___724);
	    do_fio(&c__1, (char *)&eplas[i__], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	io___725.ciunit = *mpi;
	s_wsle(&io___725);
	do_lio(&c__9, &c__1, "SOLUTION ALPHA", (ftnlen)14);
	e_wsle();
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___726.ciunit = *mpi;
	    s_wsfe(&io___726);
	    do_fio(&c__1, (char *)&alpha[i__], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	io___727.ciunit = *mpi;
	s_wsfe(&io___727);
	do_fio(&c__1, "SOLUTION XLAMBDA =", (ftnlen)18);
	do_fio(&c__1, (char *)&(*xlambda), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___728.ciunit = *mpi;
	s_wsfe(&io___728);
	do_fio(&c__1, "SOLUTION PR =", (ftnlen)13);
	do_fio(&c__1, (char *)&(*pr), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/* -------  END Comment Check Code */
/* ----------------------------------------------------------------------* */
/* ---------   FORMATS  -------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* ======================================================================= */
    return 0;
} /* solution_ */


/* Subroutine */ int store_(doublereal *statev, integer *nstatv, doublereal *
	eplas, doublereal *alpha, doublereal *xlambda, doublereal *pr, 
	doublereal *phi, doublereal *dlambda, doublereal *dtime, doublereal *
	stvit, integer *nstvi, integer *ntens, integer *ndi, integer *nvi, 
	integer *ncd, integer *m, integer *mpi)
{
    /* Format strings */
    static char fmt_101[] = "(1x,3(d20.14,2x))";
    static char fmt_102[] = "(1x,a20,6(d20.14,2x))";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void), s_wsfe(cilist *), do_fio(integer *, char *, ftnlen),
	     e_wsfe(void);

    /* Local variables */
    static integer i__;

    /* Fortran I/O blocks */
    static cilist io___730 = { 0, 0, 0, 0, 0 };
    static cilist io___731 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___732 = { 0, 0, 0, 0, 0 };
    static cilist io___733 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___734 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___735 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___736 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___737 = { 0, 0, 0, fmt_102, 0 };


/* ====================================================================== */
/* ---------------------------------------------------------------------- */
/* ----  ENREGISTREMENT DES VARIABLES UTILISATEURS A L'INSTANT N+1  ----- */
/* ---------------------------------------------------------------------- */
/* ====================================================================== */
/* ----------- */
/*  ENTREES : */
/* ----------- */
/*  EPLAS	: TENSEUR DES DEFORMATIONS PLASTIQUES */
/*  ALPHA	: TENSEUR DES CONTRAINTES CINEMATIQUES */
/*  XLAMBDA	: DEFORMATION PLASTIQUE CUMULEE */
/*  PR		: VARIABLE D'ECROUISSAGE ISOTROPE */
/*  PHI		: ISOTROPIC DAMAGE */
/*  DLAMBDA	: INCREMENT DE DEFORMATION PLASTIQUE CUMULEE */
/*  STVIT(NSTVI)	: TABLEAU CONTENANT LES VARIABLE (R,X,J2(X),Y,...) */
/*  NSTATV	: TAILLE DU TABLEAU STATEV */
/*  NVI		: NOMBRE DE VARIABLE INTERNES STOCKEES (R,X,J2(X),Y,...) */
/*  M		: COMMENTAIRES (PARAMETRE D'IMPRESSION) */
/*  MPI		: UNITE DE FICHIER DE SORTIE */
/* ----------- */
/*  SORTIES : */
/* ----------- */
/*  STAVEV(NSTATEV) : VARIABLES UTILISATEUR POUR LE MODELE UMAT */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
/* ----------------------------------------------------------------------- */
/* .1-----  Precision */
/* .2-----  Parameter,Integer */
/* .3-----  Dimension */
/* ------- */
/* .4-----  Real,Complex,Double precision,Logical,Character */
/* .5-----  Common */
/* .6-----  Equivalence */
/* .7-----  Save,Data */
/* .8-----  Definition de fonctions */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
/* ----------------------------------------------------------------------- */
/* -------  Rangement du tableau STATEV(NSTATV) */

/* ----------  Deformations plastiques  et elastiques */
/*            Ecrouissage cinematique (ALPHA) */
    /* Parameter adjustments */
    --statev;
    --stvit;
    --alpha;
    --eplas;

    /* Function Body */
    i__1 = *ntens;
    for (i__ = 1; i__ <= i__1; ++i__) {
	statev[i__] = eplas[i__];
	statev[i__ + *ntens] = alpha[i__];
    }
/* ---------  Deformation plastique cumulee */
    statev[(*ntens << 1) + 1] = *xlambda;
/* ---------  Ecrouissage isotrope (r) */
    statev[(*ntens << 1) + 2] = *pr;
/* ---------  Endommagement */

	/* always right-back the damage */
	statev[(*ntens << 1) + 3] = *phi;

#if 0
    if (*ncd == 1 || *ncd == 2) {
	statev[(*ntens << 1) + 3] = *phi;
    }
#endif
/* ---------  Increment de deformation plastique cumulee precedent */
	if (*dtime > 0.0)
    	statev[(*ntens << 1) + 4] = *dlambda / *dtime;
    else
    	statev[(*ntens << 1) + 4] = 0.0;
/* ------------------ */
/* -------  Test: si NVI > 0 ---> Stockage des variables internes suivantes: */
/* 		 - Isotropic Hardening (R) */
/* 		 - Kinematic Hardening Norm J2(X) if NVI > 1 */
/* 		 - Back Stress (X) if NVI > 2 */
/* 		 - Triaxiality Ye if NVI > 2 + NTENS */
/* 		 - Effective Strain Eff */
/* 		 - Number of Iterations Nitmax */
    i__1 = *nvi + 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	statev[i__ + 4 + (*ntens << 1)] = stvit[i__];
    }

/* -------   FIN Enregistrement du tableau STATEV(NSTATV) */
/* ----------------------------------------------------------------------- */
/* -------  Comment Check Code */
    if (*m >= 4) {
	io___730.ciunit = *mpi;
	s_wsle(&io___730);
	do_lio(&c__9, &c__1, "STORE EPLAS", (ftnlen)11);
	e_wsle();
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___731.ciunit = *mpi;
	    s_wsfe(&io___731);
	    do_fio(&c__1, (char *)&statev[i__], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	io___732.ciunit = *mpi;
	s_wsle(&io___732);
	do_lio(&c__9, &c__1, "STORE ALPHA", (ftnlen)11);
	e_wsle();
	i__1 = *ntens;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___733.ciunit = *mpi;
	    s_wsfe(&io___733);
	    do_fio(&c__1, (char *)&statev[i__ + *ntens], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}
	io___734.ciunit = *mpi;
	s_wsfe(&io___734);
	do_fio(&c__1, "STORE PHI     = ", (ftnlen)16);
	do_fio(&c__1, (char *)&statev[(*ntens << 1) + 3], (ftnlen)sizeof(
		doublereal));
	e_wsfe();
	io___735.ciunit = *mpi;
	s_wsfe(&io___735);
	do_fio(&c__1, "STORE XLAMBDA = ", (ftnlen)16);
	do_fio(&c__1, (char *)&statev[(*ntens << 1) + 1], (ftnlen)sizeof(
		doublereal));
	e_wsfe();
	io___736.ciunit = *mpi;
	s_wsfe(&io___736);
	do_fio(&c__1, "STORE PR      = ", (ftnlen)16);
	do_fio(&c__1, (char *)&statev[(*ntens << 1) + 2], (ftnlen)sizeof(
		doublereal));
	e_wsfe();
	io___737.ciunit = *mpi;
	s_wsfe(&io___737);
	do_fio(&c__1, "STORE DLAN    = ", (ftnlen)16);
	do_fio(&c__1, (char *)&statev[(*ntens << 1) + 4], (ftnlen)sizeof(
		doublereal));
	e_wsfe();
    }
/* -------  FIN Comment Check Code */
/* ----------------------------------------------------------------------- */
/* ---------   FORMATS  -------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
    return 0;
} /* store_ */


/* Subroutine */ int tensdev_(doublereal *tens, doublereal *tensdv, integer *
	ndi, integer *nshr)
{
    /* Initialized data */

    static doublereal zero = 0.;
    static doublereal trois = 3.;

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    static doublereal tensh, uttensh;

/* ====================================================================== */
/* ---------------------------------------------------------------------- */
/* --------------  CALCUL DU TENSEUR DEVIATORIQUE DE TENS  -------------- */
/* ---------------------------------------------------------------------- */
/* ====================================================================== */
/* ----------- */
/*  ENTREES : */
/* ----------- */
/*  TENS		: TENSEUR D'ORDRE DEUX */
/*  SINV1	: INVARIANT DE TENS SINV1=TRACE(TENS)/3 */
/*  NDI		: */
/*  NSHR		: */
/* ----------- */
/*  SORTIES : */
/* ----------- */
/*  TENSDV	: TENSEUR DEVIATORIQUE DE TENS */
/* ---------------------------------------------------------------------- */
/* ====================================================================== */
/* .1-----  Include,Precision */
/* .2-----  Parameter */
/* .3-----  Dimension */
/* .4-----  Real,Integer,Complex,Double precision,Logical,Character */
/* .5-----  Common */
/* .6-----  Equivalence */
/* .7-----  Save,Data */
    /* Parameter adjustments */
    --tensdv;
    --tens;

    /* Function Body */
/* .8-----  Definition de fonctions */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
/* ----------------------------------------------------------------------- */
    tensh = zero;
    i__1 = *ndi;
    for (i__ = 1; i__ <= i__1; ++i__) {
	tensh += tens[i__];
    }
    uttensh = tensh / trois;
    i__1 = *ndi;
    for (i__ = 1; i__ <= i__1; ++i__) {
	tensdv[i__] = tens[i__] - uttensh;
    }
    i__1 = *nshr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	tensdv[i__ + *ndi] = tens[i__ + *ndi];
    }
/* ---------------------------------------------------------------------- */
/* ====================================================================== */
    return 0;
} /* tensdev_ */


/* Subroutine */ int trace_(doublereal *tens, doublereal *tensh, integer *ndi,
	 integer *ntens)
{
    /* Initialized data */

    static doublereal zero = 0.;

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

/* ====================================================================== */
/* ---------------------------------------------------------------------- */
/* --------------------------  TRACE OF A TENSOR  ----------------------- */
/* ---------------------------------------------------------------------- */
/* ====================================================================== */
/* ------------- */
/*  EN ENTREE : */
/* ------------- */
/*  TENS		: TENSEUR */
/*  NDI		: NOMBRE DE COMPOSANTES DIAGONALES DU TENSEUR */
/*  NDI		: NOMBRE TOTAL DE COMPOSANTES DU TENSEUR */
/* ------------- */
/*  EN SORTIE : */
/* ------------- */
/*  TENSH	: TRACE DU TENSEUR */
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
    --tens;

    /* Function Body */
/* .8-----  Definition de fonctions */
/* ======================================================================= */
    *tensh = zero;
    i__1 = *ndi;
    for (i__ = 1; i__ <= i__1; ++i__) {
	*tensh += tens[i__];
    }
/* ======================================================================= */
    return 0;
} /* trace_ */


/* Subroutine */ int transpose_(integer *n, doublereal *a, doublereal *b)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j;

/* ====================================================================== */
/* ---------------------------------------------------------------------- */
/* ----------------------  TRANSPOSITION OF A MATRICE  ------------------ */
/* --------------------------   B = TRANSPOSED(A)  ---------------------- */
/* ---------------------------------------------------------------------- */
/* ====================================================================== */
/* ------------- */
/*  EN ENTREE : */
/* ------------- */
/*  A		: TENSEUR */
/*  N		: DIMENSION OF TENSOR */
/* ------------- */
/*  EN SORTIE : */
/* ------------- */
/*  B		: TRANSPOSED TENSOR */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
/* .1-----  Implicit, External */
/* .2-----  Parameter */
/* .3-----  Dimension */
/* .4-----  Real,Integer,Complex,Double precision,Logical,Character */
/* .5-----  Common */
/* .6-----  Equivalence */
/* .7-----  Data */
/* .8-----  Definition de fonctions */
/* ======================================================================= */
    /* Parameter adjustments */
    b_dim1 = *n;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    a_dim1 = *n;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    b[j + i__ * b_dim1] = a[i__ + j * a_dim1];
	}
    }
/* ======================================================================= */
    return 0;
} /* transpose_ */


/* Subroutine */ int valprop_(doublereal *a, doublereal *d__, doublereal *v, 
	integer *n, integer *np)
{
    /* Initialized data */

    static doublereal zero = 0.;
    static doublereal un = 1.;
    static doublereal deux = 2.;

    /* System generated locals */
    integer a_dim1, a_offset, v_dim1, v_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    /* Subroutine */ int s_paus(char *, ftnlen);
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);

    /* Local variables */
    static doublereal b, c__, e[3], f, g;
    static integer i__, j, k, l, m;
    static doublereal p, r__, s, dd;
    extern /* Subroutine */ int tred_(doublereal *, integer *, integer *, 
	    doublereal *, doublereal *);
    static integer iter;

/* ======================================================================= */
/* ----------------------------------------------------------------------- */
/* -------- CALCUL DES VALEURS PROPRES ET LES VECTEURS PROPRES   --------- */
/* -------  D'UNE MATRICE REELLE A CARREE SYMETRIQUE D'ORDRE M  ---------- */
/* ---------------------  EN RESOLVANT A.X=LAMDA.X  ---------------------- */
/* ----------------------------------------------------------------------- */
/* -----------------------------  METHODE  ------------------------------- */
/* --------------  REDUCTION DE A SOUS FORME TRI-DIAGONALE  -------------- */
/* ------------  RESOLUTION DU PROBLEME AUX VALEURS PROPRES  ------------- */
/* -----------------------  PAR L'ALGORITHME QR  ------------------------- */
/* ----------------------------------------------------------------------- */
/* ----------------  BIBLIOTHEQUE: NUMERICAL RECIPES  -------------------- */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
/* ----------- */
/*  ENTREES : */
/* ----------- */
/*  A 		: MATRICE REELLE SYMETRIQUE */
/*  N		: DIMENSION DE LA MATRICE A A DIAGONALISER */
/*  NP		: DIMENSION DE LA MATRICE A */
/* ----------- */
/*  SORTIES : */
/* ----------- */
/*  A 		: MATRICE UNITE */
/*  D 		: VECTEUR CONTENANT LES VALEURS PROPRES DE LA MATRICE A */
/*  V 		: TABLEAU DE MEME DIMENSION QUE LA MATRICE A ET */
/* 		  CONTENANT LES VECTEURS PROPRES (RANGES PAR COLONNE) */
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
    v_dim1 = *np;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --d__;
    a_dim1 = *np;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
/* .8-----  Definition de fonctions */
/* ----------------------------------------------------------------------* */
/* ======================================================================= */
/* ----------------------------------------------------------------------* */
/* -------  Tridiagonalisation de la matrice A */
    tred_(&a[8 * a_offset / 8], n, np, &d__[1], e);
/* -------  FIN  Tridiagonalisation de la matrice A */
/* ----------------------------------------------------------------------* */
/* -------  Initialisation de V */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    v[i__ + j * v_dim1] = zero;
	}
	v[i__ + i__ * v_dim1] = un;
    }
/* -------  FIN  Initialisation de V */
/* ----------------------------------------------------------------------* */
/* -------  Calcul des valeurs et vecteurs propres */

    if (*n > 1) {
	i__1 = *n;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    e[i__ - 2] = e[i__ - 1];
	}
	e[*n - 1] = zero;
/* ---------- */
	i__1 = *n;
	for (l = 1; l <= i__1; ++l) {
	    iter = 0;
/* ------------- */
L100:
	    i__2 = *n - 1;
	    for (m = l; m <= i__2; ++m) {
		dd = (d__1 = d__[m], abs(d__1)) + (d__2 = d__[m + 1], abs(
			d__2));
		if ((d__1 = e[m - 1], abs(d__1)) + dd == dd) {
		    goto L200;
		}
	    }
/* ------------- */
	    m = *n;
/* ------------- */
L200:
	    if (m != l) {
		if (iter == 30) {
		    s_paus("too many iterations", (ftnlen)19);
		}
		++iter;
		g = (d__[l + 1] - d__[l]) / (deux * e[l - 1]);
		r__ = sqrt(g * g + un);
		g = d__[m] - d__[l] + e[l - 1] / (g + d_sign(&r__, &g));
		s = un;
		c__ = un;
		p = zero;
/* ---------------- */
		i__2 = l;
		for (i__ = m - 1; i__ >= i__2; --i__) {
		    f = s * e[i__ - 1];
		    b = c__ * e[i__ - 1];
/* ------------------- */
		    if (abs(f) >= abs(g)) {
			c__ = g / f;
			r__ = sqrt(c__ * c__ + un);
			e[i__] = f * r__;
			s = un / r__;
			c__ *= s;
		    } else {
			s = f / g;
			r__ = sqrt(s * s + un);
			e[i__] = g * r__;
			c__ = un / r__;
			s *= c__;
		    }
/* ------------------- */
		    g = d__[i__ + 1] - p;
		    r__ = (d__[i__] - g) * s + deux * c__ * b;
		    p = s * r__;
		    d__[i__ + 1] = g + p;
		    g = c__ * r__ - b;
/* ------------------- */
		    i__3 = *n;
		    for (k = 1; k <= i__3; ++k) {
			f = v[k + (i__ + 1) * v_dim1];
			v[k + (i__ + 1) * v_dim1] = s * v[k + i__ * v_dim1] + 
				c__ * f;
			v[k + i__ * v_dim1] = c__ * v[k + i__ * v_dim1] - s * 
				f;
		    }
/* ------------------- */
		}
/* ---------------- */
		d__[l] -= p;
		e[l - 1] = g;
		e[m - 1] = zero;
		goto L100;
	    }
/* ------------- */
	}
/* ---------- */
    }

/* -------  Calcul des valeurs et vecteurs propres */
/* ----------------------------------------------------------------------* */
/* ======================================================================= */
/* ----------------------------------------------------------------------* */
    return 0;
} /* valprop_ */


/* Subroutine */ int tred_(doublereal *a, integer *n, integer *np, doublereal 
	*d__, doublereal *e)
{
    /* Initialized data */

    static doublereal zero = 0.;
    static doublereal un = 1.;

    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);

    /* Local variables */
    static doublereal f, g, h__;
    static integer i__, j, k, l;
    static doublereal hh, scale;

/* ======================================================================= */
/* ----------------------------------------------------------------------- */
/* ---------------  REDUCTION D'UNE REELLE SYMETRIQUE  ------------------- */
/* --------------------  SOUS FORME TRI-DIAGONALE  ----------------------- */
/* ----------------------------------------------------------------------- */
/* ----------------  BIBLIOTHEQUE: NUMERICAL RECIPES  -------------------- */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
/* ----------- */
/*  ENTREES : */
/* ----------- */
/*  A 		: MATRICE REELLE SYMETRIQUE */
/*  N		: DIMENSION DE LA MATRICE A A TRIDIAGONALISER */
/*  NP		: DIMENSION DE LA MATRICE A */
/* ----------- */
/*  SORTIES : */
/* ----------- */
/*  A 		: MATRICE UNITE */
/*  D 		: COMPOSANTES DIAGONALES DE LA MATRICE TRIDIAGONALES */
/*  E 		: COMPOSANTES NON-DIAGONALES AVEC E(1)=ZERO */
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
    --e;
    --d__;
    a_dim1 = *np;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
/* .8-----  Definition de fonctions */
/* ----------------------------------------------------------------------* */
/* ======================================================================= */
/* ----------------------------------------------------------------------* */
    if (*n > 1) {
	for (i__ = *n; i__ >= 2; --i__) {
	    l = i__ - 1;
	    h__ = zero;
	    scale = zero;
/* ------------- */
	    if (l > 1) {
		i__1 = l;
		for (k = 1; k <= i__1; ++k) {
		    scale += (d__1 = a[i__ + k * a_dim1], abs(d__1));
		}
/* ---------------- */
		if (scale == zero) {
		    e[i__] = a[i__ + l * a_dim1];
		} else {
		    i__1 = l;
		    for (k = 1; k <= i__1; ++k) {
			a[i__ + k * a_dim1] /= scale;
			h__ += a[i__ + k * a_dim1] * a[i__ + k * a_dim1];
		    }
/* ------------------- */
		    f = a[i__ + l * a_dim1];
		    d__1 = sqrt(h__);
		    g = -d_sign(&d__1, &f);
		    e[i__] = scale * g;
		    h__ -= f * g;
		    a[i__ + l * a_dim1] = f - g;
		    f = zero;
/* ------------------- */
		    i__1 = l;
		    for (j = 1; j <= i__1; ++j) {
			a[j + i__ * a_dim1] = a[i__ + j * a_dim1] / h__;
			g = zero;
			i__2 = j;
			for (k = 1; k <= i__2; ++k) {
			    g += a[j + k * a_dim1] * a[i__ + k * a_dim1];
			}
/* ---------------------- */
			if (l > j) {
			    i__2 = l;
			    for (k = j + 1; k <= i__2; ++k) {
				g += a[k + j * a_dim1] * a[i__ + k * a_dim1];
			    }
			}
/* ---------------------- */
			e[j] = g / h__;
			f += e[j] * a[i__ + j * a_dim1];
		    }
/* ------------------- */
		    hh = f / (h__ + h__);
/* ------------------- */
		    i__1 = l;
		    for (j = 1; j <= i__1; ++j) {
			f = a[i__ + j * a_dim1];
			g = e[j] - hh * f;
			e[j] = g;
			i__2 = j;
			for (k = 1; k <= i__2; ++k) {
			    a[j + k * a_dim1] = a[j + k * a_dim1] - f * e[k] 
				    - g * a[i__ + k * a_dim1];
			}
		    }
/* ------------------- */
		}
/* ---------------- */
	    } else {
		e[i__] = a[i__ + l * a_dim1];
	    }
/* ------------- */
	    d__[i__] = h__;
	}
    }
/* ------- */
/* ------- */
    d__[1] = zero;
    e[1] = zero;
/* ------- */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l = i__ - 1;
	if (d__[i__] != zero) {
	    i__2 = l;
	    for (j = 1; j <= i__2; ++j) {
		g = zero;
		i__3 = l;
		for (k = 1; k <= i__3; ++k) {
		    g += a[i__ + k * a_dim1] * a[k + j * a_dim1];
		}
/* ---------------- */
		i__3 = l;
		for (k = 1; k <= i__3; ++k) {
		    a[k + j * a_dim1] -= g * a[k + i__ * a_dim1];
		}
	    }
	}
/* ---------- */
	d__[i__] = a[i__ + i__ * a_dim1];
	a[i__ + i__ * a_dim1] = un;
/* ---------- */
	if (l >= 1) {
	    i__2 = l;
	    for (j = 1; j <= i__2; ++j) {
		a[i__ + j * a_dim1] = zero;
		a[j + i__ * a_dim1] = zero;
	    }
	}
    }
/* ----------------------------------------------------------------------* */
/* ======================================================================= */
/* ----------------------------------------------------------------------* */
    return 0;
} /* tred_ */


/* Subroutine */ int verifgence_(logical *convergence, doublereal *ff, 
	doublereal *fx, doublereal *fr, doublereal *fd, doublereal *epsf, 
	doublereal *epsx, doublereal *epsr, doublereal *epsd, integer *ntens, 
	integer *ndi, integer *nshr, integer *m, integer *mpi)
{
    /* Format strings */
    static char fmt_201[] = "(\002VERIFGENCE CONVERGENCE DE FF\002)";
    static char fmt_202[] = "(\002VERIFGENCE CONVERGENCE DE FGX\002)";
    static char fmt_203[] = "(\002VERIFGENCE CONVERGENCE DE FGR\002)";
    static char fmt_204[] = "(\002VERIFGENCE CONVERGENCE DE FH\002)";
    static char fmt_205[] = "(\002VERIFGENCE CONVERGENCE DE LA SOLUTION\002)";
    static char fmt_102[] = "(1x,a20,6(d20.14,2x))";

    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static logical dcrit, fcrit, rcrit, xcrit;
    static doublereal dabsfd, dabsff, dabsfr, dabsfx;

    /* Fortran I/O blocks */
    static cilist io___781 = { 0, 0, 0, fmt_201, 0 };
    static cilist io___783 = { 0, 0, 0, fmt_202, 0 };
    static cilist io___785 = { 0, 0, 0, fmt_203, 0 };
    static cilist io___787 = { 0, 0, 0, fmt_204, 0 };
    static cilist io___788 = { 0, 0, 0, fmt_205, 0 };
    static cilist io___789 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___790 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___791 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___792 = { 0, 0, 0, fmt_102, 0 };


/* ====================================================================== */
/* ---------------------------------------------------------------------- */
/* --------------  VERIFICATION DES CRITERES DE CONVERGENCE  ------------ */
/* ---------------------------------------------------------------------- */
/* ====================================================================== */
/* ----------- */
/*  ENTREES : */
/* ----------- */
/*  FF		: VALEUR DE LA FONCTION F (CRITERE) */
/*  FX		: VALEUR DE LA FONCTION GX (Wx = |Xt|) */
/*  FR		: VALEUR DE LA FONCTION (Wr = Sinh(Qs*R)) */
/*  FD		: VALEUR DE LA FONCTION TENSORIELLE D'ENDOMMAGEMENT */
/*  EPSF		: PRECISION DE LA FONCTION F (CRITERE) */
/*  EPSGX	: PRECISION DE LA FONCTION GX (Wx = |Xt|) */
/*  EPSGR	: PRECISION DE LA FONCTION (Wr = Sinh(Qs*R)) */
/*  EPSH		: PRECISION DE LA FONCTION TENSORIELLE D'ENDOMMAGEMENT */
/*  NTENS	: TAILLE DU TABLEAU DES CONTRAINTES OU DEFORMATIONS */
/*  NDI		: NOMBRE DE COMPOSANTES DES CONTRAINTES DIRECTES */
/*  NSHR		: NOMBRE DE COMPOSANTES DES CONTRAINTES TANGENTIELLES */
/*  M		: COMMENTAIRES (PARAMETRE D'IMPRESSION) */
/*  MPI		: UNITE DE FICHIER DE SORTIE */
/* ----------- */
/*  SORTIES : */
/* ----------- */
/*  CONVERGENCE	: VALEUR LOGIQUE INDIQUANT SI LES CRITERES SONT VERIFIES */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
/* ----------------------------------------------------------------------- */
/* .1-----  Precision */
/* .2-----  Parameter */
/* .3-----  Dimension */
/* ------- */
/* .4-----  Real,Integer,Complex,Double precision,Logical,Character */
/* .5-----  Common */
/* .6-----  Equivalence */
/* .7-----  Save,Data */
/* .8-----  Definition de fonctions */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */

/* -------  Initialisation */

    fcrit = FALSE_;
    xcrit = FALSE_;
    rcrit = FALSE_;
    dcrit = FALSE_;

/* -------  FIN Initialisation */
/* ----------------------------------------------------------------------- */

/* -------  Verification du critere |FF| < EPSF */

    dabsff = abs(*ff);
/* ------- */
    if (dabsff < *epsf) {
	fcrit = TRUE_;
	if (*m >= 4) {
	    io___781.ciunit = *mpi;
	    s_wsfe(&io___781);
	    e_wsfe();
	}
    }

/* ----------- FIN  Verification du critere |FF| < EPSF */
/* ----------------------------------------------------------------------- */

/* -------  Verification du critere |FX| < EPSX */

    dabsfx = abs(*fx);
/* ------- */
    if (dabsfx < *epsx) {
	xcrit = TRUE_;
	if (*m >= 4) {
	    io___783.ciunit = *mpi;
	    s_wsfe(&io___783);
	    e_wsfe();
	}
    }

/* ----------- FIN  Verification du critere |FX| < EPSX */
/* ----------------------------------------------------------------------- */
/* -------  Verification du critere |FR| < EPSR */

    dabsfr = abs(*fr);
/* ------- */
    if (dabsfr < *epsr) {
	rcrit = TRUE_;
	if (*m >= 4) {
	    io___785.ciunit = *mpi;
	    s_wsfe(&io___785);
	    e_wsfe();
	}
    }

/* ----------- FIN  Verification du critere |FR| < EPSR */
/* ----------------------------------------------------------------------- */
/* -------  Verification du critere |FD| < EPSD */

    dabsfd = abs(*fd);
/* ------- */
    if (dabsfd < *epsd) {
	dcrit = TRUE_;
	if (*m >= 4) {
	    io___787.ciunit = *mpi;
	    s_wsfe(&io___787);
	    e_wsfe();
	}
    }

/* ----------- FIN  Verification du critere |FD| < EPSD */
/* ----------------------------------------------------------------------- */

/* -------  Verification de la solution */

    if (fcrit && xcrit && rcrit && dcrit) {
	*convergence = TRUE_;
	if (*m >= 4) {
	    io___788.ciunit = *mpi;
	    s_wsfe(&io___788);
	    e_wsfe();
	}
    }

/* ----------- FIN  Verification de la solution */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* -------  Comment Check Code */
    if (*m >= 0) {
	io___789.ciunit = *mpi;
	s_wsfe(&io___789);
	do_fio(&c__1, "VERIFGENCE |FF| = ", (ftnlen)18);
	d__1 = abs(*ff);
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___790.ciunit = *mpi;
	s_wsfe(&io___790);
	do_fio(&c__1, "VERIFGENCE |FX| = ", (ftnlen)18);
	d__1 = abs(*fx);
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___791.ciunit = *mpi;
	s_wsfe(&io___791);
	do_fio(&c__1, "VERIFGENCE |FR| = ", (ftnlen)18);
	d__1 = abs(*fr);
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___792.ciunit = *mpi;
	s_wsfe(&io___792);
	do_fio(&c__1, "VERIFGENCE |FH| = ", (ftnlen)18);
	d__1 = abs(*fd);
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
/* -------  FIN Comment Check Code */
/* ----------------------------------------------------------------------- */
/* ---------   FORMATS  -------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* L101: */
/* ----------------------------------------------------------------------- */
/* ======================================================================= */
    return 0;
} /* verifgence_ */


/* Subroutine */ int disp_(doublereal *u, integer *kstep, integer *kinc, 
	doublereal *time, integer *node, integer *jdof)
{
    /* Initialized data */

    static doublereal un = 1.;

    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    static doublereal xl0, rate;

/* -----------------------------------------------------------------* */
/* -----------------------------------------------------------------* */
/* -----  CALCULATION OF IMPOSED DISPLACEMENT AND VELOCITY  --------* */
/* --------------  TO OBTAIN CONSTANT STRAIN RATE  -----------------* */
/* -----------------------------------------------------------------* */
/* -----------------------------------------------------------------* */
/*  ----------- */
/*  EN ENTREE : */
/*  ----------- */
/*     KSTEP   : NUMERO DU STEP */
/*     KINC    : NUMERO DE L'INCREMENT */
/*     TIME(1) : VALEUR ACTUELLE DU TEMPS DE STEP */
/*     TIME(2) : VALEUR ACTUELLE DU TEMPS TOTAL */
/*     NODE    : NUMERO DU NOEUD */
/*     JDOF    : DEGRE DE LIBERTE */
/*  ----------- */
/*  EN SORTIE : */
/*  ----------- */
/*     U(1)    : TOTAL VALUE OF THE PRESCRIBED VALUE AT THIS POINT. */
/*               THE VARIABLE MAY BE DISPLACEMENT, ROTATION, PORE */
/*               PRESSURE, TEMPERATURE, ETC., DEPENDING ON THE DEGREE */
/*               OF FREEDOM CONSTRAINED. */
/*     U(2)    : dU(1)/dTIME */
/*     U(3)    : d2U(1)/(dTIME)2 */
/*  =============================================================== */


    /* Parameter adjustments */
    --time;
    --u;

    /* Function Body */
/* ----------------------------------------------------------------- */
/* ----------------------------------------------------------------- */
/* ---------INTERPOLATION DU DEPLACEMENT IMPOSES */
/* ----------------------------------------------------------------- */

    xl0 = un;
    rate = 1.;
    u[1] = xl0 * exp(rate * time[2]) - xl0;
    u[2] = rate * xl0 * exp(rate * time[2]);
    u[3] = rate * u[2];

/* ====================================================================== */
    return 0;
} /* disp_ */

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* __F2C__ */

