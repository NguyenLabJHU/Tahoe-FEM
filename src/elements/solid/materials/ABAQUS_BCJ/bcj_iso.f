      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,                         
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,                
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C======================================================================
C----------------------------------------------------------------------
C-------------  MODELISATION D'UNE LOI DE COMPORTEMENT  ---------------
C-----  ELASTOPLASTIQUE AVEC ECROUISSAGE ISOTROPE ET CINEMATIQUE  -----
C-------------------  COUPLEE A L'ENDOMMAGEMENT  ----------------------
C----------------------------------------------------------------------
C======================================================================
C-------------
C  EN ENTREE :
C-------------
C  STRAN	: TABLEAU CONTENANT LES COMPOSANTES DE LA DEFORMATION
C		  TOTALE EN DEBUT D'INCREMENT.
C  DSTRAN	: TABLEAU DES INCREMENTS DE DEFORMATION.
C  TIME(1)	: VALEUR DU PAS DE TEMPS EN DEBUT D'INCREMENT.
C  TIME(2)	: VALEUR DU TEMPS TOTAL EN DEBUT D'INCREMENT.
C  DTIME	: INCREMENT DE TEMPS.
C  TEMP		: TEMPERATURE EN DEBUT D'INCREMENT.
C  DTEMP	: INCREMENT DE TEMPERATURE.
C  PREDEF	: TABLEAU CONTENANT LES VALEURS INTERPOLEES DES 
C		  VARIABLES DE CHAMPS PREDEFINIES AU POINT D'INTEGRATION
C		  EN DEBUT D'INCREMENT, BASEES SUR LES VALEURS PRISES
C		  AUX NOEUDS DANS L'OPTION *FIELD.
C  DPRED	: TABLEAU CONTENANT LES INCREMENTS DES VARIABLES DE
C		  CHAMPS PREDEFINIES.
C  CMNAME	: NOM DU MATERIAU DONNE PAR L'OPTION *MATERIAL.
C  NDI		: NOMBRE DE COMPOSANTES DES CONTRAINTES DIRECTES AU
C		  POINT D'INTEGRATION.
C  NSHR		: NOMBRE DE COMPOSANTES DES CONTRAINTES TANGENTIELLES
C		  POINT D'INTEGRATION.
C  NTENS	: TAILLE DU TABLEAU DES CONTRAINTES OU DEFORMATIONS
C		  (NDI + NSHR).
C  NSTATV	: NOMBRE DE VARIABLES D'ETAT ASSOCIEES DEPENDANT DE LA
C		  SOLUTION (TEL QU'IL EST DEFINI DANS L'OPTION *DEPVAR).
C  PROPS	: TABLEAU DES CONSTANTES CARACTERISTIQUES ENTREES DANS
C		  L'OPTION *USER MATERIAL POUR CE MATERIAU.
C  NPROPS	: NOMBRE DES CONSTANTES CARACTERISTIQUES (LA VALEUR
C		  EST DONNEE DANS LE PARAM£TRE 'CONSTANTS' DE L'OPTION
C		  *USER MATERIAL).
C  COORDS(3)	: TABLEAU CONTENANT LES COORDONNEES COURANTES DU POINT.
C  DROT(3,3)	: MATRICE INCREMENTALE DE ROTATION (RIGIDE).
C  CELENT	: LONGUEUR D'ELEMENT CARACTERISQUE.
C  DFGRD0	: TABLEAU (3,3) CONTENANT LES GRADIENTS DE DEFORMATION
C		  EN DEBUT D'INCREMENT.
C  DFGRD1	: TABLEAU (3,3) CONTENANT LES GRADIENTS DE DEFORMATION
C		  EN FIN D'INCREMENT.
C  NOEL		: NUMERO D'ELEMENT.
C  NPT		: NUMERO DU POINT D'INTEGRATION.
C  LAYER	:
C  KSPT		:
C  KSTEP	: NUMERO DU PAS.
C  KINC		: NUMERO D'INCREMENT.
C-----------
C  LOCALES :
C-----------
C  EPLAS	: TENSEUR DES DEFORMATIONS PLASTIQUES
C  ALPHA(6)	: TENSEUR DES CONTRAINTES CINEMATIQUES
C  XLAMBDA	: DEFORMATION PLASTIQUE CUMULEE
C  PR		: VARIABLE D'ECROUISSAGE ISOTROPE
C  PHI		: DAMAGE
C  E		: MODULE D'YOUNG
C  XNU		: COEFFICIENT DE POISSON
C  G		: MODULE DE CISAILLEMENT
C  XK		: COEFFICIENT DE LAME
C  YIELD	: LIMITE D'ELASTICITE
C  V		: CONSTANTE VISCOPLASTIQUE
C  PF		: CONSTANTE VISCOPLASTIQUE
C  RD		: CONSTANTE D'ECROUISSAGE ISOTROPE
C  H		: CONSTANTE D'ECROUISSAGE ISOTROPE
C  RS		: CONSTANTE D'ECROUISSAGE ISOTROPE (RECOVERY)
C  PRD		: CONSTANTE D'ECROUISSAGE CINEMATIQUE
C  PH		: CONSTANTE D'ECROUISSAGE CINEMATIQUE
C  PRS		: CONSTANTE D'ECROUISSAGE CINEMATIQUE (RECOVERY)
C  XLD		: LENGTH SCALE PARAMETER
C  XFN		: VOLUME FRACTION
C  XKIC		: ROUGHNESS PARAMETER
C  XLF		: VOLUME FRACTION
C  XDN		: CONSTANTE XLD*XFN/(XKIC*XLF**(1/3))
C  GS		: DUCTILE DAMAGE PARAMETERS
C  PS		: DUCTILE DAMAGE PARAMETERS
C  DC		: VALEUR D'ENDOMAGEMMENT CRITIQUE
C  CTE		: CONSTANTE DE TENSION
C  CCO		: CONSTANTE DE COMPRESSION
C  CTO		: CONSTANTE DE TORSION
C  VM		: VOID GROWTH INDICE M
C  DI		: VOID IMPINGEMENT CRITICAL DIAMETER RATIO
C  DS		: VOID SHEET CRITICAL DIAMETER RATIO
C  DC		: VALEUR DE L'ENDOMMAGEMENT CRITIQUE --> RUPTURE
C  TETA		: TETA-METHODE
C  EPSF		: PRECISION DE LA FONCTION DE CHARGE F
C  EPSG		: PRECISION DE LA FONCTION D'ENDOMMAGEMENT G
C  KCT		: KCT=0 --> MATRICE TANGENTE CONTINUE
C  		  KCT=1 --> MATRICE TANGENTE CONSISTENTE
C  NPLAN	: NPLAN=1 --> DEFORMATIONS PLANES
C		  NPLAN=2 --> CONTRAINTES PLANES
C  NATUR	: NATUR=0 --> PETITES DEFORMATIONS
C		  NATUR=1 --> GRANDES DEFORMATIONS
C  IROT		: IROT=0 --> PAS DE ROTATION
C		  IROT=1 --> AVEC ROTATION
C		: IROT=2 --> AVEC DERIVEE COROTATIONNELLE DE JAUMANN
C		  IROT=3 --> AVEC ROTATION D'ABAQUS (ROTSIG)
C  NVI		: NOMBRE DE VARIABLE INTERNES STOCKEES (EE,R,X,J2(X),Y,...)
C  M		: COMMENTAIRES (PARAMETRE D'IMPRESSION)
C  MPI		: UNITE DE FICHIER DE SORTIE
C---------------------------------------
C  VARIABLES DEVANT ETRE REACTUALISEES :
C---------------------------------------
C  STRESS	: TENSEUR DES CONTRAINTES.
C  STATEV	: TABLEAU CONTENANT LES VARIABLES D'ETAT DEPENDANT DE
C		  LA SOLUTION ET DEFINIES PAR L'UTILISATEUR. LA TAILLE
C		  DE CE TABLEAU EST DEFINIE PAR L'OPTION *DEPVAR.
C		  See Below for more information
C  SSE		: ENERGIE SPECIFIQUE ELASTIQUE DE DEFORMATION.
C  SPD		: DISSIPATION PLASTIQUE.
C  SCD		: "CREEP" DISSPATION
C-------------
C  EN SORTIE :
C-------------
C  DDSDDE	: MATRICE JACOBIENNE (NTENS,NTENS) DU MODELE THEORIQUE (DERIVEE
C		  DE L'INCREMENT DE DEFORMATION PAR L'INCREMENT DES CONTRAINTES).
C  PNEWDT	: COEFFICIENT DU NOUVEAU TEMPS D'INCREMENT ALLANT ETRE UTILISE.
C		  LES VARIABLES SUIVANTES SONT SEULEMENT DEFINIES SI ELLES SONT
C		  UTILISEES AVEC L'ANALYSE TEMPERATURE-DEPLACEMENT TOTALEMENT COUPLEE.
C  RPL		: GENERATION VOLUMETRIQUE DE CHALEUR PAR UNITE DE TEMPS EN FIN
C		  D'INCREMENT CAUSEE PAR LE TRAVAIL MECANIQUE DU MATERIAU.
C  DDSDDT	: VARIATION DES INCREMENTS DE CONTRAINTES PAR RAPPORT A LA
C		  TEMPERATURE.
C  DRPLDE	: VARIATION DE RPL PAR RAPPORT AUX INCREMENTS DE DEFORMATION.
C  DRPLDT	: VARIATION DE RPL PAR RAPPORT A L'INCREMENT DE TEMPERATURE.
C
C
C----------------------------------------------------------------------*
C----------------------------------------------------------------------*
C------  Internal State Variable Array STATEV :
C
C--  Internal state variables required for stress integration:
C
C	- Plastic Strain Tensor Ep
C	- Kinematic Hardening (Deformation) Alpha (b in BCJ Model)
C	- Isotropic Hardening (Deformation) r (R in BCJ Model)
C	- Anisotropic Damage Tensor Endo (Phi in BCJ Model)
C	- Plastic Equivalent Strain Rate Dlan at time t
C	- Total Strain only for Green-Naghdi derivative (Irot = 2)
C
C--  Internal state variables for information:
C
C	- Plastic Equivalent Strain Xlambda (PEEQ)
C	- Stored Variables if Nvi > 0:
C		- Isotropic Hardening R (Kappa in BCJ Model)
C		- Kinematic Hardening Norm J2(X) (J2(Alpha) in BCJ Model)
C		- Energy Release Rate Y
C		- Triaxiality Yt
C		- Reached interation number per increment 
C		- Nucleation Tensor Eta
C		- Void Growth V
C		- Coalescence Tensor c
C
C  *-------------------------------*-----------------*-----------------*
C  |     Component Definitions     |       2D        |       3D        |
C  *-------------------------------*-----------------*-----------------*
C  |      Plastic Strain Ep11      |    Statev(1)    |    Statev(1)    |
C  |      Plastic Strain Ep22      |    Statev(2)    |    Statev(2)    |
C  |      Plastic Strain Ep33      |    Statev(3)    |    Statev(3)    |
C  |      Plastic Strain Ep12      |    Statev(4)    |    Statev(4)    |
C  |      Plastic Strain Ep13      |                 |    Statev(5)    |
C  |      Plastic Strain Ep23      |                 |    Statev(6)    |
C  *-------------------------------*-----------------*-----------------*
C  |  Kinematic Hardening Alpha11  |    Statev(5)    |    Statev(7)    |
C  |  Kinematic Hardening Alpha22  |    Statev(6)    |    Statev(8)    |
C  |  Kinematic Hardening Alpha33  |    Statev(7)    |    Statev(9)    |
C  |  Kinematic Hardening Alpha12  |    Statev(8)    |    Statev(10)   |
C  |  Kinematic Hardening Alpha13  |                 |    Statev(11)   |
C  |  Kinematic Hardening Alpha23  |                 |    Statev(12)   |
C  *-------------------------------*-----------------*-----------------*
C  |   Plastic Equivalent Strain   |    Statev(9)    |    Statev(13)   |
C  *-------------------------------*-----------------*-----------------*
C  |     Isotropic Hardening r     |    Statev(10)   |    Statev(14)   |
C  *-------------------------------*-----------------*-----------------*
C  |      Istropic Damage Phi      |    Statev(11)   |    Statev(15)   |
C  *-------------------------------*-----------------*-----------------*
C  |           Precedent           |                 |                 |
C  |   Plastic Equivalent Strain   |    Statev(12)   |    Statev(16)   |
C  |          Rate  Dl/Dt          |                 |                 |
C  *-------------------------------*-----------------*-----------------*
C  | ------  Extra Internal State Variables for Post-Processing ------ |
C  *-------------------------------*-----------------*-----------------*
C  | Isotropic Hardening R (Kappa) |    Statev(13)   |    Statev(17)   |
C  *-------------------------------*-----------------*-----------------*
C  |   Kinematic Hardening Norm    |    Statev(14)   |    Statev(18)   |
C  |             J2(X)             |                 |                 |
C  *-------------------------------*-----------------*-----------------*
C  | Back Stress (Kin. Hard.) X11  |    Statev(15)   |    Statev(19)   |
C  | Back Stress (Kin. Hard.) X22  |    Statev(16)   |    Statev(20)   |
C  | Back Stress (Kin. Hard.) X33  |    Statev(17)   |    Statev(21)   |
C  | Back Stress (Kin. Hard.) X12  |    Statev(18)   |    Statev(22)   |
C  | Back Stress (Kin. Hard.) X13  |                 |    Statev(23)   |
C  | Back Stress (Kin. Hard.) X23  |                 |    Statev(24)   |
C  *-------------------------------*-----------------*-----------------*
C  |        Triaxiality Yt         |    Statev(19)   |    Statev(26)   |
C  *-------------------------------*-----------------*-----------------*
C  |   Reached Iteration Number    |    Statev(20)   |    Statev(27)   |
C  *-------------------------------*-----------------*-----------------*
C  |       Effective Strain        |    Statev(21)   |    Statev(28)   |
C  *-------------------------------*-----------------*-----------------*
C
C------  END Internal State Variable Array STATEV
C----------------------------------------------------------------------*
C----------------------------------------------------------------------*
C----------------------------------------------------------------------*
C=======================================================================
C----------------------------------------------------------------------*
C.1-----  Include,Precision
C      INCLUDE 'ABA_PARAM.INC'
      IMPLICIT NONE
C.2-----  Parameter,Integer
      INTEGER NSTATV,NPROPS,NTENS,NDI,NSHR,I
      INTEGER NOEL,NPT,LAYER,KSPT,KSTEP,KINC
      INTEGER NCD,NPLAN,NATUR,IROT,NVI,M,MPI
      INTEGER NSTVI
      PARAMETER (NSTVI=40)
C.3-----  Dimension
      DOUBLE PRECISION STRESS(NTENS),STATEV(NSTATV),                                   
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),                         
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),                   
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
      DOUBLE PRECISION SSE,SPD,SCD,RPL,DRPLDT
      DOUBLE PRECISION CELENT,DTIME,TEMP,DTEMP,PNEWDT
      DOUBLE PRECISION EPLAS(6),ALPHA(6)
      DOUBLE PRECISION STRESSELAS(6),ALPHAELAS(6)
      DOUBLE PRECISION STVIT(NSTVI)
      DOUBLE PRECISION XLAMBDA,DLAMBDA,PHI,PR,PRELAS,DLAN,CRITER
      DOUBLE PRECISION E,XNU,G,XK,YIELD,V,PF,HI,RDI,RSI,HK,RDK,RSK,
     1                 XIP,DC,XMS,TETA,EPSF,EPSX,EPSR,EPSD,VNEWDT
      DOUBLE PRECISION ZERO,UN
C.4-----  Real,Complex,Double precision,Logical,Character
      CHARACTER*8 CMNAME
      LOGICAL LONEWDT
C.5-----  Common
C.6-----  Equivalence
C.7-----  Save,Data
      DATA ZERO,UN/0.D0,1.D0/
      DATA MPI/6/
C.8-----  Definition de fonctions
C----------------------------------------------------------------------*
C=======================================================================
C----------------------------------------------------------------------*
C-------  Recuperation de la solution a l'increment precedent, des
C	  caracteristiques du materiau et des parametres de convergence.
C
      CALL RECOV(PROPS,NPROPS,STATEV,NSTATV,STRESS,STRAN,DSTRAN,
     *           EPLAS,ALPHA,DROT,DFGRD0,DFGRD1,XLAMBDA,PR,PHI,DLAN,
     *           TIME,DTIME,TEMP,DTEMP,E,XNU,G,XK,YIELD,V,PF,HI,RDI,RSI,
     *           HK,RDK,RSK,XIP,XMS,DC,TETA,EPSF,EPSX,EPSR,EPSD,
     *           VNEWDT,NCD,NPLAN,NATUR,IROT,NTENS,NDI,NSHR,NVI,M,MPI)
C
C-------  FIN Recuperation de la solution a l'increment precedent, des
C	  caracteristiques du materiau et des parametres de convergence.
C----------------------------------------------------------------------*
C------- Test : D = 1  --->  Pas d'increment et sortir de UMAT
C
      IF(PHI.EQ.UN) THEN
         CALL ASET(DDSDDE,ZERO,NTENS*NTENS)
         RETURN
      END IF
C
C---------  FIN  Test : D = 1  --->  Pas d'increment et sortir de UMAT
C----------------------------------------------------------------------*
C-------  Prediction elastique
C
C---------  Calcul elastique
C
         CALL PREDELAS(DSTRAN,STRESS,ALPHA,PR,PHI,STRESSELAS,
     *           ALPHAELAS,PRELAS,HI,RSI,HK,RSK,G,XK,
     *           NPLAN,NTENS,NDI,NSHR,M,MPI)
C-------  FIN Prediction elastique
C----------------------------------------------------------------------*
C------- Calcul du critere
         CALL CALCRIT(CRITER,STRESSELAS,ALPHAELAS,PRELAS,PHI,HK,HI,
     *                   YIELD,NTENS,NDI,NSHR,M,MPI)
C------- FIN Calcul du critere
C----------------------------------------------------------------------*
C------- Test : F < 0  --->  SOLUTION PUREMENT ELASTIQUE
C		F >=0  --->  SOLUTION AVEC CORRECTION PLASTIQUE
C
C---------  Solution elastique
         IF(CRITER.LT.ZERO) THEN
            IF(M.GE.3) WRITE(MPI,200)
            CALL SOLELAS(STRESS,STATEV,DDSDDE,STRESSELAS,ALPHAELAS,
     *           PRELAS,PHI,XK,G,XIP,NCD,NSTATV,NTENS,NDI,NVI,M,MPI)
            RETURN
         END IF
C---------  FIN Solution elastique
C----------------------------------------------------------------------*
C----------------------------------------------------------------------*
C---------  Solution plastique
C
         IF(M.GE.3) WRITE(MPI,201)
C
C-----------  Calcul de la solution
C
         CALL ELPLASTDO(DDSDDE,STRESS,STRAN,DSTRAN,EPLAS,ALPHA,
     *           XLAMBDA,DLAMBDA,PR,PHI,DLAN,STVIT,YIELD,V,PF,HK,RDK,
     *           RSK,HI,RDI,RSI,DC,XMS,G,XK,TETA,EPSF,EPSX,EPSR,EPSD,
     *           DTIME,LONEWDT,NCD,NPLAN,NATUR,
     *           NTENS,NDI,NSHR,NSTVI,M,MPI)
C
C---------  FIN Solution plastique
C----------------------------------------------------------------------*
C-----------  Cas D > 1 ---> Recommencer l'increment avec DTIME*PNEWDT
C
      IF(LONEWDT) THEN
         PNEWDT = VNEWDT
         RETURN
      END IF
C
C---------  FIN  Cas D > 1 ---> Recommencer l'increment avec DTIME*PNEWDT
C----------------------------------------------------------------------*
C----------------------------------------------------------------------*
C---------  Enregistrement de la solution plastique
C
         CALL STORE(STATEV,NSTATV,EPLAS,ALPHA,XLAMBDA,PR,PHI,
     *           DLAMBDA,DTIME,STVIT,NSTVI,NTENS,NDI,NVI,NCD,M,MPI)
C
C---------  FIN Enregistrement de la solution plastique
C----------------------------------------------------------------------*
C=======================================================================
C----------------------------------------------------------------------*
C---------   FORMATS  -------------------------------------------------*
C----------------------------------------------------------------------*
 200  FORMAT(/1X,'-------  SOLUTION ELASTIQUE  --------'/)
 201  FORMAT(/1X,'-------  SOLUTION PLASTIQUE  --------'/)
C----------------------------------------------------------------------*
C=======================================================================
      RETURN
      END
C
C $Id: bcj_iso.f,v 1.2 2004-01-05 07:39:36 paklein Exp $
      SUBROUTINE ADDTENS(NTENS,DT,T)
C======================================================================
C----------------------------------------------------------------------
C---------  EXECUTION DE L'OPERATION TENSORIELLE: T = T + DT  ---------
C----------------------------------------------------------------------
C======================================================================
C-------------
C  EN ENTREE :
C-------------
C  NTENS 	: DIMENSION DES TENSEURS
C  T 		: TENSEUR INITIAL
C  DT 		: INCREMENT DE TENSEUR A RAJOUTER A T
C-------------
C  EN SORTIE :
C-------------
C  T 		: TENSEUR REACTUALISE
C=======================================================================
C.1-----  Implicit, External
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C.2-----  Parameter
C.3-----  Dimension
      DIMENSION T(NTENS),DT(NTENS)
C.4-----  Real,Integer,Complex,Double precision,Logical,Character
C.5-----  Common
C.6-----  Equivalence
C.7-----  Save,Data
      DATA UDEMI/.5D0/
C.8-----  Definition de fonctions
C=======================================================================
      DO I=1,NTENS
         T(I)=T(I)+DT(I)
      END DO
C=======================================================================
      RETURN
      END
C
      SUBROUTINE AFFECT(N,V1,V2)
C======================================================================
C----------------------------------------------------------------------
C---------  EXECUTION DE L'OPERATION TENSORIELLE: T = T + DT  ---------
C----------------------------------------------------------------------
C======================================================================
C-------------
C  EN ENTREE :
C-------------
C  N	 	: DIMENSION DES TENSEURS
C  V1 		: TENSEUR INITIAL
C-------------
C  EN SORTIE :
C-------------
C  V2 		: TENSEUR REACTUALISE
C----------------------------------------------------------------------
C======================================================================
C----------------------------------------------------------------------
C.1-----  Implicit, External
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C.2-----  Parameter
C.3-----  Dimension
      DIMENSION V1(N),V2(N)
C.4-----  Real,Integer,Complex,Double precision,Logical,Character
C.5-----  Common
C.6-----  Equivalence
C.7-----  Save,Data
C.8-----  Definition de fonctions
C=======================================================================
      DO I=1,N
         V2(I)=V1(I)
      END DO
C=======================================================================
      RETURN
      END
C
      SUBROUTINE CALCOEFF(ZD,TN,ZDN,ZDZD,STRES1D,ESTAR,P,XKDH,DSQHMP,
     *           RYVSHDL,POLD,ALPHA1,PR1,PHI,PHIOLD,DLAMBDA,H,UNMD,
     *           UNMD1,WR,WX,DMWX,YIELD,V,PF,G,XK,HK,RDK,RSK,XMS,HI,
     *           RDI,RSI,TDL,UNTRWX,PHTRWX,PHTRWX2,UNTRWR,HTRWR,HTRWR2,
     *           TRSWR,RTDL,TETA,NCD,NITMAX,NTENS,NDI,NSHR,M,MPI)
C======================================================================
C----------------------------------------------------------------------
C-------------------  CALCUL DES TERMES REPETITIFS  -------------------
C----------------------------------------------------------------------
C======================================================================
C-------------
C  EN ENTREE :
C-------------
C  ESTAR	: TENSEUR DES DEFORMATIONS CONNUS A L'INSTANT T+DT
C  ALPHA1	: TENSEUR DES CONTRAINTES CINEMATIQUES CONNU A T
C  UNITDEV	: TENSEUR UNITE DEVIATORIQUE D'ORDRE 4
C  PR1		: VARIABLE D'ECROUISSAGE ISOTROPE CONNU A T
C  DLAMBDA	: VITESSE DE DEFORMATION PLASTIQUE CUMULEE
C  PHI		: DAMAGE VARIABLE
C  H		: DAMAGE FUNCTIONS h1, h2 and h3
C  YIELD	: LIMITE D'ELASTICITE
C  V		: CONSTANTE VISCOPLASTIQUE
C  PF		: CONSTANTE VISCOPLASTIQUE
C  RD		: CONSTANTE D'ECROUISSAGE ISOTROPE
C  H		: CONSTANTE D'ECROUISSAGE ISOTROPE
C  RS		: CONSTANTE D'ECROUISSAGE ISOTROPE (RECOVERY)
C  QS		: CONSTANTE D'ECROUISSAGE ISOTROPE (RECOVERY)
C  PRD		: CONSTANTE D'ECROUISSAGE CINEMATIQUE
C  PH		: CONSTANTE D'ECROUISSAGE CINEMATIQUE
C  PRS		: CONSTANTE D'ECROUISSAGE CINEMATIQUE (RECOVERY)
C  VM		: VOID GROWTH INDICE M
C  TETA		: TETA-METHODE
C  NDI		: NOMBRE DE COMPOSANTES DES CONTRAINTES DIRECTES
C  NSHR		: NOMBRE DE COMPOSANTES DES CONTRAINTES TANGENTIELLES
C  NTENS	: TAILLE DU TABLEAU DES CONTRAINTES OU DEFORMATIONS
C		  (NDI + NSHR)
C  M		: COMMENTAIRES
C  MPI		: UNITE DE FICHIER DE SORTIE
C-------------
C  EN SORTIE :
C-------------
C  ZD		: TENSEUR DE DIRECTION NORMALE
C  TN		: TENSEUR NORMAL
C  ZDN		: NORME DU TENSEUR ZD (1ERE METHODE)
C  ZDZD		: NORME DU TENSEUR ZD (2NDE METHODE)
C  DM		: TENSEUR EFFET DU DOMMAGE
C  DMI		: TENSEUR EFFET DU DOMMAGE INVERSE
C  DMWX		: NORME DU TENSEUR Mt:(Alpha1+Tdl*N)
C  STRESS1D	: TENSEUR DEVIATORIQUE DES CONTRAINTES CONNU A T
C  UNMD		: EXEMPLE D'UN COEFFICIENT
C  NITMAX	: NOMBRE D'ITERATIONS DU SCHEMA DE NEWTON
C-----------------------------------------------------------------------
C=======================================================================
C-----------------------------------------------------------------------
C=======================================================================
C.1-----  Implicit, External
      IMPLICIT NONE
C.2-----  Parameter,Integer
      INTEGER NITMAX,NCD,NTENS,NDI,NSHR,M,MPI,I
C.3-----  Dimension
      DOUBLE PRECISION ZD(NTENS),TN(NTENS)
      DOUBLE PRECISION ESTAR(NTENS),STRES1D(NTENS),ALPHA1(NTENS)
      DOUBLE PRECISION H(3)
C.-----
      DOUBLE PRECISION ZDN,ZDZD,DTLF,TNTN,P,POLD,XKDH,DSQHMP,RYVSHDL
      DOUBLE PRECISION DLAMBDA,PHI,PHIOLD,PR1,TDL,DLRD,UNMD,UNMD1,DPHI,
     1       WR,WX,DMWX,YIELD,V,PF,G,XK,HK,RDK,RSK,HI,RDI,RSI,XMS,
     2       UNTRWX,PHTRWX,PHTRWX2,UNTRWR,HTRWR,HTRWR2,TRSWR,RTDL,TETA
      DOUBLE PRECISION ZERO,ONE
C.4-----  Real,Complex,Double precision,Logical,Character
C.5-----  Common
C.6-----  Equivalence
C.7-----  Save,Data
      DATA ZERO,ONE/0.D0,1.D0/
C.8-----  Definition de fonctions
C---------------------------------------------------------------------*
C======================================================================
C---------------------------------------------------------------------*
C-------  Number of iterations
C
      NITMAX  = NITMAX+1
C
C-------  END Number of iterations
C---------------------------------------------------------------------*
C-------  Calculation of coefficients
C
      UNMD    = ONE-PHI
      UNMD1   = ONE/UNMD
      DPHI    = PHI-PHIOLD
C-------
      TDL     = TETA*DLAMBDA
      DLRD    = TDL*UNMD1
C-------
      DTLF    = DLAMBDA/PF
C-------
      UNTRWX  = ONE+TETA*(DLAMBDA*RDK+RSK)*WX
      PHTRWX  = HK/UNTRWX
      PHTRWX2 = PHTRWX/UNTRWX
C-------
      UNTRWR  = ONE+TETA*(DLAMBDA*RDI+RSI)*WR
      HTRWR   = HI/UNTRWR
      HTRWR2  = HTRWR/UNTRWR
      TRSWR   = TETA*RSI*WR
      RTDL    = PR1+TDL
C
C-------  END  Calculation of coefficients
C---------------------------------------------------------------------*
C-------  Tensor TN
C
      CALL CALZDTN(ZD,TN,ZDN,ZDZD,DSQHMP,RYVSHDL,STRES1D,ESTAR,
     *           ALPHA1,P,PR1,PHI,DPHI,POLD,XKDH,H,TDL,DTLF,
     *           UNMD,DMWX,PHTRWX,RTDL,HTRWR,YIELD,V,PF,G,XK,XMS,
     *           NCD,NTENS,NDI,NSHR,M,MPI)
C
C-------  END Tensor TN
C---------------------------------------------------------------------*
C-------  Comment Check Code
      IF(M.GE.4) THEN
         WRITE(*,*)'CALCOEFF TDL     = ',TDL
         WRITE(*,*)'CALCOEFF TETA    = ',TETA
         WRITE(*,*)'CALCOEFF DLAMBDA = ',DLAMBDA
         WRITE(*,*)'CALCOEFF PHTRWX  = ',PHTRWX
         WRITE(*,*)'CALCOEFF PHTRWX2 = ',PHTRWX2
         WRITE(*,*)'CALCOEFF HTRWR   = ',HTRWR
         WRITE(*,*)'CALCOEFF HTRWR2  = ',HTRWR2
         CALL PDTSCA(NDI,NSHR,TN,TN,TNTN)
         WRITE(*,*)'CALCOEFF TNTN    = ',TNTN
         WRITE(MPI,*)'CALCOEFF TN'
         DO I=1,NTENS
            WRITE(MPI,101)TN(I)
         END DO
      END IF
C-------  END Comment Check Code
C-----------------------------------------------------------------------
C---------   FORMATS  --------------------------------------------------
C-----------------------------------------------------------------------
 101  FORMAT(1X,6(D20.14,2X))
 102  FORMAT(1X,A20,6(D20.14,2X))
C-----------------------------------------------------------------------
C=======================================================================
      RETURN
      END
C
      SUBROUTINE CALCOEFFI(ZD,TN,ZDN,ZDZD,PHI,STRES1D,ESTAR,P,XKDH,
     *           DSQHMP,RYVSHDL,ALPHA1,POLD,PR1,PHIOLD,UNMD,UNMD1,H,
     *           WR,WX,V,PF,G,XK,HK,RDK,RSK,HI,RDI,RSI,XMS,DLAMBDA,
     *           TDL,DMWX,YIELD,UNTRWX,PHTRWX,PHTRWX2,UNTRWR,HTRWR,
     *           HTRWR2,TRSWR,RTDL,TETA,
     *           NITMAX,NCD,NTENS,NDI,NSHR,M,MPI)
C======================================================================
C----------------------------------------------------------------------
C-------------------  CALCUL DES TERMES REPETITIFS  -------------------
C----------------------------------------------------------------------
C======================================================================
C-------------
C  EN ENTREE :
C-------------
C  ESTAR	: TENSEUR DES DEFORMATIONS CONNUS A L'INSTANT T+DT
C  ALPHA1	: TENSEUR DES CONTRAINTES CINEMATIQUES CONNU A T
C  UNITDEV	: TENSEUR UNITE DEVIATORIQUE D'ORDRE 4
C  PR1		: VARIABLE D'ECROUISSAGE ISOTROPE CONNU A T
C  DLAMBDA	: VITESSE DE DEFORMATION PLASTIQUE CUMULEE
C  PHI		: DAMAGE VARIABLE
C  H		: DAMAGE FUNCTIONS h1, h2 and h3
C  YIELD	: LIMITE D'ELASTICITE
C  V		: CONSTANTE VISCOPLASTIQUE
C  PF		: CONSTANTE VISCOPLASTIQUE
C  RD		: CONSTANTE D'ECROUISSAGE ISOTROPE
C  H		: CONSTANTE D'ECROUISSAGE ISOTROPE
C  RS		: CONSTANTE D'ECROUISSAGE ISOTROPE (RECOVERY)
C  PRD		: CONSTANTE D'ECROUISSAGE CINEMATIQUE
C  PH		: CONSTANTE D'ECROUISSAGE CINEMATIQUE
C  PRS		: CONSTANTE D'ECROUISSAGE CINEMATIQUE (RECOVERY)
C  VM		: VOID GROWTH INDICE M
C  TETA		: TETA-METHODE
C  NDI		: NOMBRE DE COMPOSANTES DES CONTRAINTES DIRECTES
C  NSHR		: NOMBRE DE COMPOSANTES DES CONTRAINTES TANGENTIELLES
C  NTENS	: TAILLE DU TABLEAU DES CONTRAINTES OU DEFORMATIONS
C			  (NDI + NSHR)
C  M		: COMMENTAIRES
C  MPI		: UNITE DE FICHIER DE SORTIE
C-------------
C  EN SORTIE :
C-------------
C  ZD		: TENSEUR DE DIRECTION NORMALE
C  TN		: TENSEUR NORMAL
C  ZDN		: NORME DU TENSEUR ZD (1ERE METHODE)
C  ZDZD		: NORME DU TENSEUR ZD (2NDE METHODE)
C  DM		: TENSEUR EFFET DU DOMMAGE
C  DMI		: TENSEUR EFFET DU DOMMAGE INVERSE
C  DMWX		: NORME DU TENSEUR Mt:(Alpha1+Tdl*N)
C  STRESS1D	: TENSEUR DEVIATORIQUE DES CONTRAINTES CONNU A T
C  DLAMBDA	: VITESSE DE DEFORMATION PLASTIQUE CUMULEE
C  NITMAX	: INITIALISATION DU NOMBRE D'ITERATIONS (A ZERO)
C-----------------------------------------------------------------------
C=======================================================================
C-----------------------------------------------------------------------
C=======================================================================
C.1-----  Implicit, External
      IMPLICIT NONE
C.2-----  Parameter,Integer
      INTEGER NITMAX,NCD,NTENS,NDI,NSHR,M,MPI,I
C.3-----  Dimension
      DOUBLE PRECISION ZD(NTENS),TN(NTENS)
      DOUBLE PRECISION ESTAR(NTENS),STRES1D(NTENS),ALPHA1(NTENS)
      DOUBLE PRECISION H(3)
C.-----
      DOUBLE PRECISION ZDN,ZDZD,DTLF,PHI,P,POLD,XKDH,DSQHMP,RYVSHDL
      DOUBLE PRECISION PR1,PHIOLD,UNMD,UNMD1,DPHI,WR,WX,YIELD,V,PF,G,
     1       XK,HK,RDK,RSK,HI,RDI,RSI,XMS,DLAMBDA,TDL,DMWX,UNTRWX,
     2       PHTRWX,PHTRWX2,UNTRWR,HTRWR,HTRWR2,TRSWR,RTDL,TETA
      DOUBLE PRECISION ZERO,ONE
C.4-----  Real,Complex,Double precision,Logical,Character
C.5-----  Common
C.6-----  Equivalence
C.7-----  Save,Data
      DATA ZERO,ONE/0.D0,1.D0/
C.8-----  Definition de fonctions
C---------------------------------------------------------------------*
C======================================================================
C---------------------------------------------------------------------*
C-------  Number of iterations
C
      NITMAX  = 0
C
C-------  END  Number of iterations
C----------------------------------------------------------------------*
C-------  Calculation of coefficients
C
      TDL     = TETA*DLAMBDA
C-------
      DTLF    = DLAMBDA/PF
C-------
      UNTRWX  = ONE+TETA*(RDK*DLAMBDA+RSK)*WX
      PHTRWX  = HK/UNTRWX
      PHTRWX2 = PHTRWX/UNTRWX
C-------
      UNTRWR  = ONE+TETA*RSI*WR
      HTRWR   = HI/UNTRWR
      HTRWR2  = HTRWR/UNTRWR
      TRSWR   = TETA*RSI*WR
      RTDL    = PR1+TDL
C----------
      UNMD    = ONE-PHI
      UNMD1   = ONE/UNMD
      DPHI    = PHI-PHIOLD
C
C-------  END  Calculation of coefficients
C----------------------------------------------------------------------*
C-------  Tensors ZD and TN, and norms ZDN and ZDN
C
      CALL CALZDTN(ZD,TN,ZDN,ZDZD,DSQHMP,RYVSHDL,STRES1D,ESTAR,
     *           ALPHA1,P,PR1,PHI,DPHI,POLD,XKDH,H,TDL,DTLF,
     *           UNMD,DMWX,PHTRWX,RTDL,HTRWR,YIELD,V,PF,G,XK,XMS,
     *           NCD,NTENS,NDI,NSHR,M,MPI)
C
C-------  END Tensors ZD and TN, norms ZDN and ZDN
C----------------------------------------------------------------------*
C-------  Comment Check Code
      IF(M.GE.4) THEN
         WRITE(MPI,*)'CALCOEFFI TN'
         DO I=1,NTENS
            WRITE(MPI,101)TN(I)
         END DO
      END IF
C-------  END Comment Check Code
C----------------------------------------------------------------------*
C---------   FORMATS  -------------------------------------------------*
C----------------------------------------------------------------------*
 101  FORMAT(1X,6(D20.14,2X))
 102  FORMAT(1X,A20,6(D20.14,2X))
C-----------------------------------------------------------------------
C=======================================================================
      RETURN
      END
C
      SUBROUTINE CALCRIT(CRITER,STRESS,ALPHA,PR,PHI,HK,HI,
     *                   YIELD,NTENS,NDI,NSHR,M,MPI)
C======================================================================
C----------------------------------------------------------------------
C----------------  CALCUL DE LA FONCTION-SEUIL DE MISES  --------------
C----------------------------------------------------------------------
C======================================================================
C-------------
C  EN ENTREE :
C-------------
C  STRESS	: CONTRAINTES A L'INSTANT N
C  ALPHA	: TENSEUR DES CONTRAINTES CINEMATIQUES A L'INSTANT N
C  PR		: VARIABLE D'ECROUISSAGE ISOTROPE
C  PHI		: DAMAGE
C  H		: PARAMETRE D'ECROUISSAGE ISOTROPE
C  PH		: PARAMETRE D'ECROUISSAGE CINEMATIQUE
C  YIELD	: MODULE D'YOUNG
C  NTENS	: LONGUEUR DES TABLEAU ALPHA, STRESS, ...
C  NDI		: NOMBRE DE CONTRAINTES DIRECTES (S11,S22,S33)
C  NSHR		: NOMBRE DE CONTRAINTES TANGENTIELLES (S12,S13,S23)
C  M		: COMMENTAIRES
C  MPI		: UNITE DE FICHIER DE SORTIE
C-----------
C  LOCALES :
C-----------
C  SIGMXT	: CONTRAINTES EFFECTIVES A L'INSTANT N
C-------------
C  EN SORTIE :
C-------------
C  CRITER	:  FONCTION CRITERE
C----------------------------------------------------------------------*
C=======================================================================
C----------------------------------------------------------------------*
C.1-----  Implicit, External
      IMPLICIT NONE
C.2-----  Parameter,Integer
      INTEGER NTENS,NDI,NSHR,M,MPI,I
C.3-----  Dimension
      DOUBLE PRECISION STRESS(NTENS),ALPHA(NTENS)
      DOUBLE PRECISION SIGMX(6),X(6)
C-------  
      DOUBLE PRECISION CRITER,PR,HK,HI,YIELD,R,SIGMXH,SIGMXEQV
      DOUBLE PRECISION PHI,UNMD,UNMDHK
      DOUBLE PRECISION ZERO,UN,DEUX,TROIS
C.4-----  Real,Complex,Double precision,Logical,Character
C.5-----  Common
C.6-----  Equivalence
C.7-----  Save,Data
      DATA ZERO,UN,DEUX,TROIS/0.D0,1.D0,2.D0,3.D0/
C.8-----  Definition de fonctions
C=======================================================================
C----------------------------------------------------------------------*
C-------  Second Stress Invariant
C
      DO I=1,NTENS
         X(I)     = HK*ALPHA(I)
         SIGMX(I) = STRESS(I)-X(I)
      END DO
      CALL SINV(SIGMX,SIGMXH,SIGMXEQV,NDI,NSHR)
C----------------------------------------------------------------------*
C-------  Isotropic Hardening
C
      R = HI*PR
C----------------------------------------------------------------------*
C-------  Yield Function
C
      UNMD   = UN-PHI
      CRITER = SIGMXEQV/UNMD-R-YIELD
C----------------------------------------------------------------------*
C-------  Comment Check Code
C
      IF(M.GE.2)THEN
         WRITE(MPI,102)'CALCRIT CRITERE =',CRITER
         IF(M.GE.3)THEN
            WRITE(MPI,*)'CALCRIT SIGMX'
            DO I=1,NTENS
               WRITE(MPI,101)SIGMX(I)
            END DO
            WRITE(MPI,102)'CALCRIT SIGMXEQV=',SIGMXEQV
         END IF
      END IF
C-------  END Comment Check Code
C----------------------------------------------------------------------*
C---------   FORMATS  -------------------------------------------------*
C----------------------------------------------------------------------*
 101  FORMAT(1X,6(D20.14,2X))
 102  FORMAT(1X,A20,6(D20.14,2X))
C-----------------------------------------------------------------------
C=======================================================================
      RETURN
      END
C
      SUBROUTINE CALDDE(DLDE,DDDE,DNDE,UNN,TNTN,UZDNN,ALDE,DPNDE,DZDE,
     *           DXIDE,ESTAR,ALPHA1,ZD,TN,TETN,ALTN,PHI,ZDZD,DFDL,DFDX,
     *           DFDR,DFDD,DXDL,DXDX,DXDD,DRDL,DRDR,DHDL,DHDX,DHDR,
     *           DHDD,DFDS,DXDS,DRDS,DHDS,DEDT,P,EYV,DLAMBDA,TDL,UNMD,
     *           UNMD1,DMWX,WX,PHTRWX,PHTRWX2,H,G,XK,RDK,RSK,RSI,XMS,
     *           TETA,NCD,NTENS,NDI,NSHR,M,MPI)
C======================================================================
C----------------------------------------------------------------------
C-----------------------  CALCUL DES DERIVEES  ------------------------
C-------------  PAR RAPPORT AU TENSEUR DES DEFORMATIONS  --------------
C----------------------------------------------------------------------
C======================================================================
C-------------
C  EN ENTREE :
C-------------
C  TETN		: TENSEUR = ESTAR - TDLRD*TN
C  ALTN		: TENSEUR = ALPHA1 + TDLRD*TN
C  PLTN		: PRODUCT P:|dEp/dl|
C  PLTN		: PRODUCT Q:|dEp/dl|
C  ESTAR	: TENSEUR e*
C  DEDT		: TENSEUR INCREMENT DE DEFORMATION : DSTRAN
C  ALPHA1	: TENSEUR D'ECROUISSAGE CINEMATIQUE CONNU
C  TN		: TENSEUR NORMAL
C  ZDZD		: NORME DU TENSEUR ZD = RACINE (2/3*ZD:ZD)
C  DFDL		: DERIVEE PARTIELLE DE F PAR RAPPORT A DLAMBDA
C  DFDX		: DERIVEE PARTIELLE DE F PAR RAPPORT A WX
C  DFDR		: DERIVEE PARTIELLE DE F PAR RAPPORT A WR
C  DFDD		: DERIVEE PARTIELLE DE F PAR RAPPORT A ENDO
C  DXDL		: DERIVEE PARTIELLE DE GX PAR RAPPORT A DLAMBDA
C  DXDX		: DERIVEE PARTIELLE DE GX PAR RAPPORT A WX
C  DXDD		: DERIVEE PARTIELLE DE GX PAR RAPPORT A ENDO
C  DRDL		: DERIVEE PARTIELLE DE GR PAR RAPPORT A DLAMBDA
C  DRDR		: DERIVEE PARTIELLE DE GR PAR RAPPORT A WR
C  DHDL		: DERIVEE PARTIELLE DE H PAR RAPPORT A DLAMBDA
C  DHDX		: DERIVEE PARTIELLE DE H PAR RAPPORT A WX
C  DHDR		: DERIVEE PARTIELLE DE H PAR RAPPORT A WR
C  DHDD		: DERIVEE PARTIELLE DE H PAR RAPPORT A ENDO
C  DFDS		: DERIVEE PARTIELLE DE F PAR RAPPORT A DLs
C  DXDS		: DERIVEE PARTIELLE DE FGX PAR RAPPORT A DLs
C  DRDS		: DERIVEE PARTIELLE DE FGR PAR RAPPORT A DLs
C  DHDS		: DERIVEE PARTIELLE DE H PAR RAPPORT A DLs
C  DMWX		: NORME DU TENSEUR Mt:(Alpha1+Tdl*N)
C  DLAMBDA	: INCREMENT DEFORMATION PLASTIQUE CUMULEE
C  YPS1		: VALEUR DE LA FONCTION YPS = (Y/S)**(s-1)
C  TDL		: TETA*DLAMBDA
C  DLRD		: TETA*DLAMBDA/UNMD
C  WR		: ISOTROPIC DYNAMIC RECOVERY TERM: WR = DSINH[QS*R]
C  WX		: KINEMATIC DYNAMIC RECOVERY TERM: WX = |DMI:X|
C  PHTRWX	: (2/3)*PH/(1+TETA*(PRD*DLAMBDA+PRS)*WX)
C  PHTRWX2	: (2/3)*PH/(1+TETA*(PRD*DLAMBDA+PRS)*WX)**2
C  G		: MODULE DE CISAILLEMENT
C  XK		: COEFFICIENT DE LAME
C  RDK		: CONSTANTE D'ECROUISSAGE CINEMATIQUE
C  RSK		: CONSTANTE D'ECROUISSAGE CINEMATIQUE (RECOVERY)
C  RSI		: CONSTANTE D'ECROUISSAGE ISOTROPE (RECOVERY)
C  XMS		: DAMAGE RATE SENSITIVITY
C  TETA		: TETA-METHODE
C  NCD		: NCD=1 --> COUPLAGE ELASTOPLASTICITE-ENDOMMAGEMENT
C  		  NCD=0 --> ELASTOPLASTICITE SANS ENDOMMAGEMENT
C  NDI		: NOMBRE DE COMPOSANTES DES CONTRAINTES DIRECTES
C  NSHR		: NOMBRE DE COMPOSANTES DES CONTRAINTES TANGENTIELLES
C  NTENS	: NOMBRE TOTAL DE COMPOSANTES D'UN TENSEUR
C  M		: COMMENTAIRES
C  MPI		: UNITE DE FICHIER DE SORTIE
C-------------
C  EN LOCAL  :
C-------------
C  TNTN		: TENSEUR (D'ORDRE 4) = N*N
C  ALDE		: TENSEUR (D'ORDRE 4) = ALPHA1*DLDE
C  UDEVNN	: TENSEUR (D'ORDRE 4) = IDEV - 2/3 N*N
C  DPNDE	: TENSEUR (D'ORDRE 4) = 2G/ZDZD.(IDEV - 2/3 N*N)
C  DPXDE	: DERIVEE PARTIELLE DE LA FONCTION FGX PAR RAPPORT A DE
C  DPHDE	: DERIVEE PARTIELLE DE LA FONCTION H PAR RAPPORT A DE
C-------------
C  EN SORTIE :
C-------------
C  DLDE		: TENSEUR (D'ORDRE 2) DERIVEE DE DLAMBDA 
C  DDDE		: TENSEUR (D'ORDRE 2) DERIVEE DE D
C  DNDE		: TENSEUR (D'ORDRE 4) DERIVEE DU TENSEUR NORMAL N
C  UDEVNN	: TENSEUR (D'ORDRE 4) = IDEV - 2/3 N*N
C-----------------------------------------------------------------------
C=======================================================================
C-----------------------------------------------------------------------
C=======================================================================
C.1-----  Implicit, External
      IMPLICIT NONE
C.2-----  Parameter,Integer
      INTEGER NCD,NTENS,NDI,NSHR,M,MPI,I,J
C.3-----  Dimension
      DOUBLE PRECISION ESTAR(NTENS),ALPHA1(NTENS),DEDT(NTENS)
      DOUBLE PRECISION TN(NTENS),ZD(NTENS)
      DOUBLE PRECISION DLDE(NTENS),DDDE(NTENS),DNDE(NTENS,NTENS)
      DOUBLE PRECISION ALDE(NTENS,NTENS),DPNDE(NTENS,NTENS)
      DOUBLE PRECISION DZDE(NTENS,NTENS),DXIDE(NTENS,NTENS)
      DOUBLE PRECISION TETN(NTENS),ALTN(NTENS)
      DOUBLE PRECISION UZDNN(NTENS,NTENS)
      DOUBLE PRECISION UNN(NTENS,NTENS),TNTN(NTENS,NTENS)
      DOUBLE PRECISION XI(6),EEAL(6),EEALDV(6)
      DOUBLE PRECISION DXDE(6),DLDXDE(6),DRDE(6),DHDE(6)
      DOUBLE PRECISION DPFDE(6),DPXDE(6),DPHDE(6),FGXDD(6)
      DOUBLE PRECISION DSDE(6),DPXDSDE(6),DPNDET(36),DXIDET(36)
      DOUBLE PRECISION DSEFFDE(6),XIDXIDE(6)
      DOUBLE PRECISION A(6),B(6),H(3)
C-------
      DOUBLE PRECISION ZDZD,DFDL,DFDX,DFDR,DFDD,DFDS
      DOUBLE PRECISION DXDL,DXDX,DXDR,DXDD,DXDS,DRDL,DRDR,DRDS
      DOUBLE PRECISION DHDL,DHDX,DHDR,DHDD,DHDS
      DOUBLE PRECISION DLAMBDA,TDL,DLRD,DLS,ESTARH,UNMD,UNMD1
      DOUBLE PRECISION P,DMWX,WX,PHTRWX,PHTRWX2,
     1       G,XK,RDK,RSK,RSI,XMS,TETA
      DOUBLE PRECISION D2G,D2G,DQTDG,UZDZD,D2GT,XKDGT,
     1       D2GD,D2GD3,DLGDHX,DXK,H2PDXK,UNMDXK
      DOUBLE PRECISION TDLPHDMWX,TRDWX,TRSWX,TRDDLRS
      DOUBLE PRECISION SEFF,SEFF2,XIEQV,XIEQV2
      DOUBLE PRECISION EYV,XMM1,APHI,CPHI,PHI,XMS2,
     1       UNMDM,AMM,COSHA,DLCA,DLCAY,DLCAXDK,ADL,PI,ADLDKP
      DOUBLE PRECISION DH1DD,DH2DD,DH3DD
      DOUBLE PRECISION H1SEFF,H2PKSEF
C-----------------------------------------------------------*
      DOUBLE PRECISION UDXDX,UDRDR,DFXX,DFRR,DGXX,DGRR,FRXS,
     1       DXRR,DHXX,DHRR
      DOUBLE PRECISION DHS,AL,AD,BL,BD,ALBD,AB,ABL
      DOUBLE PRECISION ZERO,ONE,TWO,THREE,XNINE
      DOUBLE PRECISION HALF,ONETHIRD,TWOTHIRD,THALF
C.4-----  Real,Complex,Double precision,Logical,Character
C.5-----  Common
C.6-----  Equivalence
C.7-----  Save,Data
      DATA ZERO,ONE,TWO,THREE,XNINE/0.D0,1.D0,2.D0,3.D0,9.D0/
      DATA HALF,THALF/0.5D0,1.5D0/
C.8-----  Definition de fonctions
C---------------------------------------------------------------------*
C======================================================================
C---------------------------------------------------------------------*
C-------  Initialization
C
      ONETHIRD = ONE/THREE
      TWOTHIRD = TWO/THREE
      D2G      = TWO*G
      D2GD     = D2G*UNMD
      DQTDG    = D2GD*TWOTHIRD
      UZDZD    = ONE/ZDZD
      D2GT     = D2G/THREE
      XKDGT    = XK-D2GT
C
C-------  END Initialization
C---------------------------------------------------------------------*
C---------------------------------------------------------------------*
C-------  Calcul des tenseurs UZDNN et DPNDE
C
      CALL PDTTEN(TN,TN,TNTN,NTENS,NDI)
      DO I=1,NTENS
         DO J=1,NTENS
            UNN(I,J) = -TWOTHIRD*TNTN(I,J)
         END DO
         UNN(I,I) = UNN(I,I)+ONE
      END DO
      DO I=1,NDI
         DO J=1,NDI
	    UNN(I,J) = UNN(I,J)-ONETHIRD
	 END DO
      END DO
C-------
      DO I=1,NTENS
         DO J=1,NTENS
	    UZDNN(I,J) = UZDZD*UNN(I,J)
            DPNDE(I,J) = D2GD*UZDNN(I,J)
         END DO
      END DO
C
C-------  END Calcul des tenseurs UZDNN et DPNDE
C---------------------------------------------------------------------*
C---------------------------------------------------------------------*
C-------  Tensor DPFDE
C
      DO I=1,NTENS
         DPFDE(I) = DQTDG*TN(I)
      END DO
C
C-------  END Tensor DPFDE
C---------------------------------------------------------------------*
C---------------------------------------------------------------------*
C-------  Tensor DPXDE
C
      DO I=1,NTENS
         ALTN(I) = ALPHA1(I)+TDL*TN(I)
      END DO
      CALL TRANSPOSE(NTENS,DPNDE,DPNDET)
      CALL PDTMAT(NTENS,DPNDET,ALTN,FGXDD)
C-------
      TDLPHDMWX = TDL*PHTRWX/DMWX
      DO I=1,NTENS
         DPXDE(I) = -TDLPHDMWX*FGXDD(I)
      END DO
C
C-------  END  Tensor DPXDE
C---------------------------------------------------------------------*
C---------------------------------------------------------------------*
C-------  Definition of the Tensor DSDE (Static Recovery)
C
      IF(RSK.EQ.ZERO .AND. RSI.EQ.ZERO) THEN
         CALL ASET(DSDE,ZERO,NTENS)
      ELSE
         CALL PDTSCA(NDI,NSHR,TN,DEDT,DLS)
	 IF(DLS.EQ.ZERO) DLS = ONE
	 DO I=1,NTENS
	    DSDE(I) = TN(I)/DLS
	 END DO
      END IF
C
C-------  END  Definition of the Tensor DSDE (Static Recovery)
C---------------------------------------------------------------------*
C---------------------------------------------------------------------*
C-------  Partial Derivative DPHDE
C
      IF(NCD.EQ.0) THEN
         CALL ASET(DPHDE,ZERO,NTENS)
         GOTO 100
      END IF
C
C------------------------------------------------------*
C------------------------------------------------------*
C-------  Partial Derivative DPHDE
C
      D2GD3  = ONETHIRD*D2GD
      DLGDHX = TDL*(D2GD+PHTRWX)
      DO I=1,NTENS
         XI(I) = ZD(I)-DLGDHX*TN(I)
      END DO
C-------------------------------------------*
C---------  Norm |Xi| = Sqrt(3/2 * Xi:Xi)
      XIEQV = XI(1)*XI(1)+XI(2)*XI(2)+XI(3)*XI(3)
      DO I=1+NDI,NTENS
         XIEQV = XIEQV+TWO*XI(I)*XI(I)
      END DO
      XIEQV = DSQRT(THALF*XIEQV)
      SEFF  = DSQRT(H(1)*XIEQV*XIEQV+H(2)*P*P)
C-------------------------------------------*
C---------  Partial Derivative DXIDE = dXi/dE
      DO I=1,NTENS
         DO J=1,NTENS
            DXIDE(I,J) = -DLGDHX*DPNDE(I,J)
         END DO
         DXIDE(I,I) = DXIDE(I,I)+D2GD
      END DO
      DO I=1,NDI
         DO J=1,NDI
            DXIDE(I,J) = DXIDE(I,J)-D2GD3
         END DO
      END DO
C-------------------------------------------*
C---------  Product Xi : dXi/dE
      CALL TRANSPOSE(NTENS,DXIDE,DXIDET)
      CALL PDTMAT(NTENS,DXIDET,XI,XIDXIDE)
C-------------------------------------------*
C-------  Comment Check Code
      IF(M.GE.4) THEN
         WRITE(MPI,102)'STRESS XI'
         DO I=1,NTENS
            WRITE(MPI,101)XI(I)
         END DO
      END IF
C-----------------------------------------------------------*
C-----------------------------------------------------------*
C--------- Choice of Damage Evolution
      IF(NCD.EQ.1) THEN
C-----------------------------------------------------------*
C---------  Cocks [1989]
C-----------------------------------------------------------*
         EYV  = ONETHIRD*P/SEFF
         XMM1 = XMS/(XMS+ONE)
         AD   = TWO*(XMS+ONE)*(ONE+PHI)*UNMD**XMM1
         APHI = XNINE*XMS*PHI*UNMD*EYV/AD
C-----------------------------------------------------------*
         DH1DD  = TWOTHIRD
         DH2DD  = HALF*XMM1/(UNMD1*UNMD1)
         DH3DD  = XMM1*UNMD1**(XMM1-ONE)
         DH1DD  = ZERO
         DH2DD  = ZERO
         DH3DD  = ZERO
C-----------------------------------------------------------*
C---------  Partial Derivative of DSEFFDE = d[1/Seff]/dE
         DXK     = ONETHIRD*UNMD*XK
         XIEQV2  = XIEQV*XIEQV
         SEFF2   = SEFF*SEFF
         H1SEFF  = THALF*H(1)/SEFF2
         H2PKSEF = DXK*H(2)*P/SEFF
         DO I=1,NTENS
            DSEFFDE(I) = H1SEFF*XIDXIDE(I)
         END DO
         DO I=1,NDI
            DSEFFDE(I) = DSEFFDE(I)+H2PKSEF
         END DO
C-------------------------------------------*
C---------  Partial Derivative DPHDE = dFphi/dE
         ADL = APHI*DLAMBDA
         PI  = P
         IF(P.EQ.ZERO) PI = ONE
         ADLDKP = ADL*DXK/PI
C----------*
         DO I=1,NTENS
            DPHDE(I) = ADL*DSEFFDE(I)
         END DO
         DO I=1,NDI
            DPHDE(I) = DPHDE(I)-ADLDKP
         END DO
C-----------------------------------------------------------*
      ELSE
C-----------------------------------------------------------*
C---------  Cocks and Ashby [1980]
C-----------------------------------------------------------*
         EYV   = ONETHIRD*P/XIEQV
         XMS2  = TWO*XMS
         UNMDM = UNMD**XMS
         AMM   = TWO*(XMS2-ONE)/(XMS2+ONE)
         APHI  = AMM*EYV
         CPHI  = ONE/UNMDM-UNMD
C-------------------------------------------*
         XIEQV2  = XIEQV*XIEQV
         UNMDXK  = ONETHIRD*UNMD*XK
         COSHA   = DCOSH(APHI)
         DLCA    = DLAMBDA*CPHI*COSHA*AMM
         DLCAY   = THALF*DLCA*EYV/XIEQV2
         DLCAXDK = UNMDXK*DLCA/XIEQV
C-------------------------------------------*
C---------  Partial Derivative DHDE = dFphi/dE 
         DO I=1,NTENS
            DPHDE(I) = DLCAY*XIDXIDE(I)
         END DO
         DO I=1,NDI
            DPHDE(I) = DPHDE(I)-DLCAXDK
         END DO
C-------------------------------------------*
C-----------------------------------------------------------*
      END IF
C
C-------  END Partial Derivative DPHDE = dH/dE
C---------------------------------------------------------------------*
C---------------------------------------------------------------------*
C-------  Tensor DDDE
C
 100  UDXDX = ONE/DXDX
      UDXDX = ONE/DXDX
      UDRDR = ONE/DRDR
      DFXX  = DFDX*UDXDX
      DHXX  = DHDX*UDXDX
      DFRR  = DFDR*UDRDR
      DHRR  = DHDR*UDRDR
C-------
      DLS   = DFDS-DFXX*DXDS-DFRR*DRDS
      DHS   = DHDS-DHXX*DXDS-DHRR*DRDS
      AL    = DFDL-DFXX*DXDL-DFRR*DRDL
      AD    = DFDD-DFXX*DXDD
      BL    = DHDL-DHXX*DXDL-DHRR*DRDL
      BD    = DHDD-DHXX*DXDD
C------------------------------------------*
      DO I=1,NTENS
         A(I)  = DPFDE(I)-DFXX*DPXDE(I)+DLS*DSDE(I)
         B(I)  = DPHDE(I)-DHXX*DPXDE(I)+DHS*DSDE(I)
      END DO
C-------
      IF(NCD.EQ.0) THEN
         AD = ZERO
         CALL ASET(DDDE,ZERO,NTENS)
      ELSE
         ALBD = ONE/(AL*BD-AD*BL)
         DO I=1,NTENS
            DDDE(I) = ALBD*(BL*A(I)-AL*B(I))
         END DO
      END IF
C
C-------  END  Tensor DDDE
C---------------------------------------------------------------------*
C---------------------------------------------------------------------*
C-------  Tensor DLDE
C
      DO I=1,NTENS
         DLDE(I) = -(A(I)+AD*DDDE(I))/AL
      END DO
C
C-------  END  Tensor DLDE
C---------------------------------------------------------------------*
C---------------------------------------------------------------------*
C-------  Tensors DXDE and DRDE
C
      DO I=1,NTENS
         DXDE(I) = DXDL*DLDE(I)+DXDD*DDDE(I)+DXDS*DSDE(I)+DPXDE(I)
         DXDE(I) = -UDXDX*DXDE(I)
	 DRDE(I) = -UDRDR*(DRDL*DLDE(I)+DRDS*DSDE(I))
      END DO
C
C-------  END  Tensors DXDE and DRDE
C---------------------------------------------------------------------*
C---------------------------------------------------------------------*
C-------  Tensor DNDE
C
      TRDWX   = TETA*RDK*WX
      TRSWX   = TETA*RSK*WX
      TRDDLRS = TETA*RDK*DLAMBDA+RSK
      DO I=1,NTENS
         DLDXDE(I) = TRDWX*DLDE(I)+TRDDLRS*DXDE(I)+TRSWX*DSDE(I)
      END DO
      CALL PDTTEN(ALPHA1,DLDXDE,ALDE,NTENS,NDI)
C-------
      DO I=1,NTENS
         EEAL(I) = -D2G*ESTAR(I)
      END DO
      CALL TENSDEV(EEAL,EEALDV,NDI,NSHR)
      CALL PDTTEN(EEALDV,DDDE,DZDE,NTENS,NDI)
C-------
      DO I=1,NTENS
         DO J=1,NTENS
            DZDE(I,J) = DZDE(I,J)+PHTRWX2*ALDE(I,J)
         END DO
         DZDE(I,I) = DZDE(I,I)+D2GD
      END DO
      CALL PRDMAT(NTENS,UZDNN,DZDE,DNDE)
C
C-------  END  Tensor DNDE
C---------------------------------------------------------------------*
C-------  Comment Check Code
C
      IF(M.GE.4) THEN
         WRITE(MPI,*)'CALDDE DLDE'
         DO I=1,NTENS
            WRITE(MPI,101)DLDE(I)
         END DO
         WRITE(MPI,*)'CALDDE DXDE'
         DO I=1,NTENS
            WRITE(MPI,101)DXDE(I)
         END DO
         WRITE(MPI,*)'CALDDE DRDE'
         DO I=1,NTENS
            WRITE(MPI,101)DRDE(I)
         END DO
         WRITE(MPI,*)'CALDDE DDDE'
         DO I=1,NTENS
            WRITE(MPI,101)DDDE(I)
         END DO
         WRITE(MPI,*)'CALDDE UZDNN'
         DO I=1,NTENS
            WRITE(MPI,101)(UZDNN(I,J),J=1,NTENS)
         END DO
         WRITE(MPI,*)'CALDDE DPNDE'
         DO I=1,NTENS
            WRITE(MPI,101)(DPNDE(I,J),J=1,NTENS)
         END DO
         WRITE(MPI,*)'CALDDE ALDE'
         DO I=1,NTENS
            WRITE(MPI,101)(ALDE(I,J),J=1,NTENS)
         END DO
         WRITE(MPI,*)'CALDDE DNDE'
         DO I=1,NTENS
            WRITE(MPI,101)(DNDE(I,J),J=1,NTENS)
         END DO
      END IF
C-------  END Comment Check Code
C----------------------------------------------------------------------*
C---------   FORMATS  -------------------------------------------------*
C----------------------------------------------------------------------*
 101  FORMAT(1X,6(D20.13,2X))
 102  FORMAT(1X,A20,6(D20.14,2X))
C----------------------------------------------------------------------*
C=======================================================================
      RETURN
      END
C
      SUBROUTINE CALDPF(DFDL,DFDX,DFDR,DFDD,DFDS,DRDS,DPDD,DZDNDL,
     *           DZDNDD,ALPHZD,DFRDL,ZD,H,P,PHI,ALPHA1,STRES1D,PR1,
     *           ZDN,ZDZD,DSQHMP,PHTRWX,RYVSHDL,DLAMBDA,TDL,UNMD,
     *           UNMD1,XKDH,WX,WR,PHTRWX2,HTRWR,HTRWR2,TRSWR,RTDL,
     *           RDK,RSK,RDI,RSI,YIELD,G,XK,V,PF,XMS,XMM1,TETA,
     *           NCD,MELPST,NTENS,NDI,NSHR,M,MPI)
C======================================================================
C----------------------------------------------------------------------
C---------------  CALCUL DES DERIVEES PARTIELLES DE F  ----------------
C--------------------  PAR RAPPORT DLAMBDA ET ENDO  -------------------
C----------------------------------------------------------------------
C======================================================================
C-------------
C  EN ENTREE :
C-------------
C  ZD		: TENSEUR DE DIRECTION NORMALE
C  ENDO		: TENSEUR D'ENDOMMAGEMENT
C  ESTAR	: TENSEUR DES DEFORMATIONS CONNUS A T+DT
C  ALPHA1	: TENSEUR D'ECROUISSAGE CINEMATIQUE CONNU
C  PR1		: VARIABLE D'ECROUISSAGE ISOTROPE CONNUE
C  ZDN		: NORME DU TENSEUR DE DIRECTION NORMALE ZD
C  ZDZD		: NORME DU ZD = RACINE (2/3*ZD:ZD)
C  TDL		: TETA*DLAMBDA
C  WX		: KINEMATIC DYNAMIC RECOVERY TERM: WX = |DMI:X|
C  WR		: ISOTROPIC DYNAMIC RECOVERY TERM: WR = R
C  PHTRWX	: (2/3)*PH/(1+TETA*(RDK*DLAMBDA+RSK)*WX)
C  PHTRWX2	: (2/3)*PH/(1+TETA*(RDK*DLAMBDA+RSK)*WX)*RSK)*WX)**2
C  HTRWR	: RT(2/3)*H/(1+TETA*(RDK*DLAMBDA+RS))
C  HTRWR2	: RT(2/3)*H/(1+TETA*(RDK*DLAMBDA+RS))**2
C  TRSWR	: TETA*RS*WR
C  RTDL		: PR1+TETA*DLAMBDA
C  RD		: CONSTANTE D'ECROUISSAGE ISOTROPE
C  RS		: CONSTANTE D'ECROUISSAGE ISOTROPE (RECOVERY)
C  RDK		: CONSTANTE D'ECROUISSAGE CINEMATIQUE
C  RSK		: CONSTANTE D'ECROUISSAGE CINEMATIQUE (RECOVERY)
C  YIELD	: LIMITE D'ELASTICITE
C  G		: MODULE DE CISAILLEMENT
C  V		: CONSTANTE VISCOPLASTIQUE
C  PF		: CONSTANTE VISCOPLASTIQUE
C  TETA		: TETA-METHODE
C  NCD		: NCD=1 --> COUPLAGE ELASTOPLASTICITE-ENDOMMAGEMENT
C  		  NCD=0 --> ELASTOPLASTICITE SANS ENDOMMAGEMENT
C  MELPST	: MELPST = 0 --> MATRICE CONSISTANTE AVEC V = 1
C  NDI		: NOMBRE DE COMPOSANTES DES CONTRAINTES DIRECTES
C  NSHR		: NOMBRE DE COMPOSANTES DES CONTRAINTES TANGENTIELLES
C  NTENS	: TAILLE DU TABLEAU DES CONTRAINTES OU DEFORMATIONS
C		  (NDI + NSHR).
C  M		: COMMENTAIRES
C  MPI		: UNITE DE FICHIER DE SORTIE
C-------------
C  EN SORTIE :
C-------------
C  DFDL		: DERIVEE PARTIELLE DE F PAR RAPPORT A DLAMBDA
C  DFDX		: DERIVEE PARTIELLE DE F PAR RAPPORT A WX
C  DFDR		: DERIVEE PARTIELLE DE F PAR RAPPORT A WR
C  DFDD		: DERIVEE PARTIELLE DE F PAR RAPPORT A ENDO
C  DZDNDL	: DERIVEE PARTIELLE DE ZDN PAR RAPPORT A DLAMBDA
C  DFRDL	: DERIVEE PARTIELLE DE R PAR RAPPORT A DLAMBDA
C  DMDDESTAR	: PRODUIT CONTRACTE DMDD : ESTAR
C  ALPHZD	: PRODUIT CONTRACTE ALPHA1 : ZD
C  DMDDALPH	: PRODUIT CONTRACTE DMDD : ALPHA1
C-----------------------------------------------------------------------
C=======================================================================
C-----------------------------------------------------------------------
C.1-----  Implicit, External
      IMPLICIT NONE
C.2-----  Parameter,Integer
      INTEGER NCD,MELPST,NTENS,NDI,NSHR,M,MPI,I
C.3-----  Dimension
      DOUBLE PRECISION ZD(NTENS),ALPHA1(NTENS),STRES1D(NTENS)
      DOUBLE PRECISION H(3)
C-------
      DOUBLE PRECISION DFDL,DFDX,DFDR,DFDD
      DOUBLE PRECISION DFDS,DPDD,DZDNDL,DZDNDD,DFRDL,DRDS
      DOUBLE PRECISION ALPHZD,ZDN,ZDZD,STRESZD,DLAMBDA,PR1,PHI,TDL,
     1       UNMD,UNMD1,WX,WR,PHTRWX,PHTRWX2,HTRWR,HTRWR2,TRSWR,TRSWX,
     2       RTDL,RDK,RSK,RDI,RSI,YIELD,G,XK,V,PF,XMS,XMM1,TETA
      DOUBLE PRECISION D2G,DTALPHZDZD,UZDZD,DTD2GZD,DTPHTRWXZD,UPSD1
      DOUBLE PRECISION P,XIEQV,DSQHMP,RYVSHDL,XKDH,ONDSQHMP,DRYVSH
      DOUBLE PRECISION DH1DD,DH2DD,DH3DD,H1H1,H2H1,H3H1,DH2H1DD,DH3H1DD
      DOUBLE PRECISION ZERO,ONE,TWO,THREE,HALF,ONETHIRD,TWOTHIRD
C.4-----  Real,Integer,Complex,Double precision,Logical,Character
C.5-----  Common
C.6-----  Equivalence
C.7-----  Save,Data
      DATA ZERO,ONE,TWO,THREE/0.D0,1.D0,2.D0,3.D0/
      DATA HALF/0.5D0/
C.8-----  Definition de fonctions
C---------------------------------------------------------------------*
C======================================================================
C---------------------------------------------------------------------*
C-------  Initialisation
      ONETHIRD = ONE/THREE
      TWOTHIRD = TWO*ONETHIRD
      D2G      = TWO*G
      H2H1     = H(2)/H(1)
      H3H1     = H(3)*H(3)/H(1)
C-------  END Initialisation
C---------------------------------------------------------------------*
C-------  Partial Derivative DZNDDL = dZDN/dlambda
      IF(MELPST.EQ.0) V = ZERO
      ONDSQHMP   = ONE/DSQHMP
      DFRDL      = TETA*HTRWR2*(ONE+TRSWR-PR1*RDI*WR)
      DRYVSH     = DFRDL+V/DSQRT(PF*PF+DLAMBDA*DLAMBDA)
      DZDNDL     = TWOTHIRD*UNMD*H3H1*DRYVSH*RYVSHDL*ONDSQHMP
     *            +TETA*(D2G*UNMD+PHTRWX-TDL*RDK*WX*PHTRWX2)
C-------  END Partial Derivative DFDL = df/dlambda
C---------------------------------------------------------------------*
C-------  Partial Derivative DFDL = df/dlambda
      CALL PDTSCA(NDI,NSHR,ALPHA1,ZD,ALPHZD)
      DTALPHZDZD = TWOTHIRD*ALPHZD/ZDZD
      DFDL       = TETA*RDK*WX*PHTRWX2*DTALPHZDZD-DZDNDL
C-------  END Partial Derivative DFDL = df/dlambda
C---------------------------------------------------------------------*
C-------  Partial Derivative DFDX = df/dWx
      DFDX = TETA*(DLAMBDA*RDK+RSK)*PHTRWX2*(TDL+DTALPHZDZD)
C-------  END Partial Derivative DFDX = df/dWx
C---------------------------------------------------------------------*
C-------  Partial Derivative DFDR = df/dWr
      DFDR = TETA*(DLAMBDA*RDI+RSI)*HTRWR2*RTDL
      DFDR = TWOTHIRD*DFDR*UNMD*RYVSHDL*ONDSQHMP
C-------  END Partial Derivative DFDR = df/dWr
C---------------------------------------------------------------------*
C---------  Partial Derivative DHxDD = dhx/dphi
      IF(NCD.EQ.1) THEN
         UPSD1  = ONE/(ONE+PHI)
         DH1DD  = TWOTHIRD
         DH2DD  = HALF*XMM1*UPSD1*UPSD1
         DH3DD  = XMM1*UNMD1**(XMM1-ONE)
      ELSE IF(NCD.EQ.2) THEN
         DH1DD  = ZERO
         DH2DD  = HALF/(UNMD1*UNMD1)
         DH3DD  = -ONE
      END IF
      DH1DD  = ZERO
      DH2DD  = ZERO
      DH3DD  = ZERO
C---------
      H1H1    = ONE/(H(1)*H(1))
      H2H1    = H(2)/H(1)
      DH2H1DD = (H(1)*DH2DD-H(2)*DH1DD)/H1H1
      DH3H1DD = (TWO*H(1)*H(3)*DH3DD-H(3)*H(3)*DH1DD)/H1H1
C---------------------------------------------------------------------*
C-------  Partial Derivative DFDD = df/dD
C               Remark:  dZ/dD(t):Z = Z:dZ/dD
      IF(NCD.EQ.0) THEN
         DFDD = ZERO
      ELSE
         CALL PDTSCA(NDI,NSHR,STRES1D,ZD,STRESZD)
         DPDD   = -XKDH
         DZDNDD = (DH3H1DD*RYVSHDL*RYVSHDL-P*P*DH2H1DD-TWO*H2H1*P*DPDD)
         DZDNDD = HALF*UNMD*ONDSQHMP*DZDNDD-DSQHMP
         DZDNDD = TWOTHIRD*DZDNDD-D2G*DLAMBDA
         DFDD   = -TWOTHIRD*STRESZD/ZDZD-DZDNDD
      END IF
C-------  END Partial Derivative DFDD = df/dD
C---------------------------------------------------------------------*
C-------  Partial Derivative DRDS = dR/dLs
      TRSWR = TETA*RSI*WR
      DRDS  = UNMD*TRSWR*HTRWR2*RTDL
C-------  END Partial Derivative DRDS = dR/dLs
C---------------------------------------------------------------------*
C-------  Partial Derivative DFDS = df/dLs
      TRSWX = TETA*RSK*WX
      DFDS  = TWOTHIRD*DRDS+TRSWX*PHTRWX2*(TDL+TWOTHIRD*ALPHZD/ZDZD)
C-------  END Partial Derivative DFDS = df/dLs
C---------------------------------------------------------------------*
C---------------------------------------------------------------------*
C-------  Comment Check Code
C
      IF(M.GE.4) THEN
         WRITE(MPI,102)'CALDPF  DZDNDL = ',DZDNDL
         WRITE(MPI,102)'CALDPF  DZDNDD = ',DZDNDD
         WRITE(MPI,102)'CALDPF  DFDL   = ',DFDL
         WRITE(MPI,102)'CALDPF  DFDX   = ',DFDX
         WRITE(MPI,102)'CALDPF  DFDR   = ',DFDR
         WRITE(MPI,102)'CALDPF  DFDD   = ',DFDD
         WRITE(MPI,102)'CALDPF  DFDS   = ',DFDS
         WRITE(MPI,102)'CALDPF  DRDS   = ',DRDS
      END IF
C
C-------  END Comment Check Code
C-----------------------------------------------------------------------
C---------   FORMATS  --------------------------------------------------
C-----------------------------------------------------------------------
 101  FORMAT(1X,6(D20.14,2X))
 102  FORMAT(1X,A20,6(D20.14,2X))
C-----------------------------------------------------------------------
C=======================================================================
      RETURN
      END
C
      SUBROUTINE CALDPG(DXDL,DXDX,DXDR,DXDD,DXDS,DRDL,DRDR,
     *           DNDL,DNDX,DNDR,DNDD,DNDS,DFDR,DFRDL,DMWX,TN,ALPHA1,
     *           DLAMBDA,TDL,WX,UNTRWX,PHTRWX,PHTRWX2,RDK,RSK,
     *           TETA,NCD,NTENS,NDI,NSHR,M,MPI)
C======================================================================
C----------------------------------------------------------------------
C------------  CALCUL DES DERIVEES PARTIELLES DE GX et GR  ------------
C----------------  PAR RAPPORT DLAMBDA, WX, WR ET ENDO  ---------------
C----------------------------------------------------------------------
C======================================================================
C-------------
C  EN ENTREE :
C-------------
C  DNDL		: DERIVEE PARTIELLE DE N PAR RAPPORT A DLAMBDA
C  DNDX		: DERIVEE PARTIELLE DE N PAR RAPPORT A WX
C  DNDX		: DERIVEE PARTIELLE DE N PAR RAPPORT A WR
C  DNDD		: DERIVEE PARTIELLE DE N PAR RAPPORT A ENDO
C  DNDX		: DERIVEE PARTIELLE DE N PAR RAPPORT A LS
C  DFDR		: DERIVEE PARTIELLE DE F PAR RAPPORT A WR
C  DFRDL	: DERIVEE PARTIELLE DE R PAR RAPPORT A DLAMBDA
C  DM		: TENSEUR EFFET D'ENDOMMAGEMENT (Mt)
C  DMDDALPH	: PRODUIT CONTRACTE DMDD : ALPHA1
C  DMWX		: NORME DU TENSEUR Mt:(Alpha1+Tdl*N)
C  TN		: TENSEUR NORMAL
C  ALPHA1	: TENSEUR D'ECROUISSAGE CINEMATIQUE CONNU
C  DLAMBDA	: VITESSE DE DEFORMATION PLASTIQUE CUMULEE
C  PR1		: VARIABLE D'ECROUISSAGE ISOTROPE CONNUE
C  TDL		: TETA*DLAMBDA
C  WR		: ISOTROPIC DYNAMIC RECOVERY TERM: WR = DSINH[QS*R]
C  WX		: KINEMATIC DYNAMIC RECOVERY TERM: WX = |DMI:X|
C  UNTRWX	: 1+TETA*(PRD*DLAMBDA+PRS)*WX
C  PHTRWX	: (2/3)*PH/(1+TETA*(PRD*DLAMBDA+PRS)*WX)
C  PHTRWX2	: (2/3)*PH/(1+TETA*(PRD*DLAMBDA+PRS)*WX)**2
C  PRD		: CONSTANTE D'ECROUISSAGE CINEMATIQUE
C  PRS		: CONSTANTE D'ECROUISSAGE CINEMATIQUE (RECOVERY)
C  TETA		: TETA-METHODE
C  NCD		: NCD=1 --> COUPLAGE ELASTOPLASTICITE-ENDOMMAGEMENT
C  		  NCD=0 --> ELASTOPLASTICITE SANS ENDOMMAGEMENT
C  NDI		: NOMBRE DE COMPOSANTES DES CONTRAINTES DIRECTES
C  NSHR		: NOMBRE DE COMPOSANTES DES CONTRAINTES TANGENTIELLES
C  NTENS	: TAILLE DU TABLEAU DES CONTRAINTES OU DEFORMATIONS
C		  (NDI + NSHR).
C  M		: COMMENTAIRES
C  MPI		: UNITE DE FICHIER DE SORTIE
C-------------
C  EN SORTIE :
C-------------
C  DXDL		: DERIVEE PARTIELLE DE GX PAR RAPPORT A DLAMBDA
C  DXDX		: DERIVEE PARTIELLE DE GX PAR RAPPORT A WX
C  DXDD		: DERIVEE PARTIELLE DE GX PAR RAPPORT A ENDO
C  DRDL		: DERIVEE PARTIELLE DE GR PAR RAPPORT A DLAMBDA
C  DRDR		: DERIVEE PARTIELLE DE GR PAR RAPPORT A WR
C-----------------------------------------------------------------------
C=======================================================================
C-----------------------------------------------------------------------
C.1-----  Implicit, External
      IMPLICIT NONE
C.2-----  Parameter,Integer
      INTEGER NCD,NTENS,NDI,NSHR,M,MPI,I
C.3-----  Dimension
      DOUBLE PRECISION DNDL(NTENS),DNDX(NTENS),DNDR(NTENS),DNDD(NTENS)
      DOUBLE PRECISION DNDS(NTENS)
      DOUBLE PRECISION TN(NTENS),ALPHA1(NTENS)
      DOUBLE PRECISION ALTN(6),DNDLN(6),DNDDALTN(6)
C-------
      DOUBLE PRECISION DXDL,DXDX,DXDR,DXDD,DXDS,DRDL,DRDR
      DOUBLE PRECISION DFDR,DFRDL,DMWX,DLAMBDA,TDL,WX
      DOUBLE PRECISION UNTRWX,PHTRWX,PHTRWX2,TRSWX,RDK,RSK,TETA
      DOUBLE PRECISION UDMWX,DXDALPHA,RDDLRS,PHTRWXUW,TRSWR
      DOUBLE PRECISION ZERO,ONE,TWO,THREE,THALF,TWOTHIRD,SQTT
C.4-----  Real,Complex,Double precision,Logical,Character
C.5-----  Common
C.6-----  Equivalence
C.7-----  Save,Data
      DATA ZERO,ONE,TWO,THREE/0.D0,1.D0,2.D0,3.D0/
      DATA THALF/1.5D0/
C.8-----  Definition de fonctions
C----------------------------------------------------------------------
C======================================================================
C---------------------------------------------------------------------*
C-------  Initialisation
C
      TWOTHIRD = TWO/THREE
      SQTT  = DSQRT(TWOTHIRD)
C-------
      IF(DMWX.EQ.ZERO) THEN
         UDMWX = ONE
      ELSE
         UDMWX = ONE/DMWX
      END IF
C
C-------  FIN Initialisation
C----------------------------------------------------------------------*
C-------  Partial Derivative DXDL = dX/dlambda
C
      DO I=1,NTENS
         ALTN(I) = ALPHA1(I)+TDL*TN(I)
         DNDLN(I)  = DLAMBDA*DNDL(I)+TN(I)
      END DO
      CALL PDTSCA(NDI,NSHR,ALTN,DNDLN,DXDALPHA)
      DXDL = TETA*PHTRWX2*(RDK*WX*DMWX-UNTRWX*DXDALPHA*UDMWX)
C
C-------  FIN Partial Derivative DXDL = dGX/dlambda
C---------------------------------------------------------------------*
C-------  Partial Derivative DXDX = dX/dWx
C
      CALL PDTSCA(NDI,NSHR,ALTN,DNDX,DXDALPHA)
      RDDLRS = RDK*DLAMBDA+RSK
      DXDX   = ONE+TETA*RDDLRS*PHTRWX2*DMWX-TDL*PHTRWX*DXDALPHA*UDMWX
C
C-------  FIN Partial Derivative DXDX = dX/dWx
C---------------------------------------------------------------------*
C-------  Partial Derivative DXDR = dX/dWr
C
      CALL PDTSCA(NDI,NSHR,ALTN,DNDR,DXDALPHA)
      DXDR = -TDL*PHTRWX*DXDALPHA*UDMWX
      DXDR = ZERO
C
C-------  FIN Partial Derivative DXDR = dX/dWr
C---------------------------------------------------------------------*
C-------  Partial Derivative DXDD = dX/dD
C
      IF(NCD.EQ.0) THEN
         DXDD = ZERO
      ELSE
         CALL PDTSCA(NDI,NSHR,ALTN,DNDD,DXDALPHA)
         DXDD = -TDL*PHTRWX*DXDALPHA*UDMWX
      END IF
C
C-------  FIN Partial Derivative DXDD = dX/dD
C---------------------------------------------------------------------*
C---------------------------------------------------------------------*
C-------  Partial Derivative (Kinematic Hardening) DXDS = dX/dLs
      CALL PDTSCA(NDI,NSHR,ALTN,DNDS,DXDALPHA)
      TRSWX = TETA*RSK*WX
      DXDS  = PHTRWX2*(TRSWX*DMWX-UNTRWX*TDL*DXDALPHA*UDMWX)
C-------  END  Partial Derivative (Kinematic Hardening) DXDS = dX/dLs
C---------------------------------------------------------------------*
C-------  Partial Derivative DRDL = dR/dlambda
      DRDL = -DFRDL
C-------  FIN Partial Derivative DRDL = dR/dlambda
C---------------------------------------------------------------------*
C-------  Partial Derivative DRDR = dR/dD
      DRDR = ONE+THALF*DFDR
C-------  FIN Partial Derivative DRDR = dR/dD
C---------------------------------------------------------------------*
C---------------------------------------------------------------------*
C-------  Comment Check Code
C
      IF(M.GE.4) THEN
         WRITE(MPI,102)'CALDPG  DXDL =',DXDL
         WRITE(MPI,102)'CALDPG  DXDX =',DXDX
         WRITE(MPI,102)'CALDPG  DXDR =',DXDR
         WRITE(MPI,102)'CALDPG  DXDD =',DXDD
         WRITE(MPI,102)'CALDPG  DXDS =',DXDS
         WRITE(MPI,102)'CALDPG  DRDL =',DRDL
         WRITE(MPI,102)'CALDPG  DRDR =',DRDR
      END IF
C
C-------  FIN Comment Check Code
C----------------------------------------------------------------------*
C---------   FORMATS  -------------------------------------------------*
C----------------------------------------------------------------------*
 101  FORMAT(1X,6(D20.14,2X))
 102  FORMAT(1X,A20,6(D20.14,2X))
C----------------------------------------------------------------------*
C=======================================================================
      RETURN
      END
C
      SUBROUTINE CALDPH(DHDL,DHDX,DHDR,DHDD,DHDS,DNDL,DNDX,DNDR,DNDD,
     *           DNDS,DPDD,ZD,TN,P,PHI,STRES1D,ALPHA1,ALTN,DLAMBDA,TDL,
     *           H,UNMD,UNMD1,WX,PHTRWX,PHTRWX2,G,RDK,RSK,XMS,XMM1,
     *           TETA,NCD,NTENS,NDI,NSHR,M,MPI)
C======================================================================
C----------------------------------------------------------------------
C---------------  CALCUL DES DERIVEES PARTIELLES DE H  ----------------
C---------------  PAR RAPPORT DLAMBDA, WX, WR ET ENDO  ----------------
C----------------------------------------------------------------------
C======================================================================
C-------------
C  EN ENTREE :
C-------------
C  DNDL		: DERIVEE PARTIELLE DE N PAR RAPPORT A DLAMBDA
C  DNDX		: DERIVEE PARTIELLE DE N PAR RAPPORT A WX
C  DNDD		: DERIVEE PARTIELLE DE N PAR RAPPORT A ENDO
C  TN		: TENSEUR NORMAL
C  ENDO		: TENSEUR D'ENDOMMAGEMENT
C  ALPHA1	: TENSEUR D'ECROUISSAGE CINEMATIQUE CONNU
C  ALTN		: TENSEUR = ALPHA1 + TDL*M(-t)*TN
C  PR1		: VARIABLE D'ECROUISSAGE ISOTROPE CONNUE
C  ZDN		: NORME DU TENSEUR DE DIRECTION NORMALE ZD
C  DFDD		: DERIVEE PARTIELLE DE F PAR RAPPORT A ENDO
C  DZDNDL	: DERIVEE PARTIELLE DE ZDN PAR RAPPORT A DLAMBDA
C  ALPHZD	: PRODUIT CONTRACTE ALPHA1 : ZD
C  DM		: TENSEUR D'EFFET D'ENDOMMAGEMENT  (Mt)
C  DMI		: TENSEUR D'EFFET D'ENDOMMAGEMENT INVERSE (M(-t))
C  DN		: COMPONANTES DU VECTEUR ORIENTATION (NTENS,NTENS)
C  QD		: VOID COALESCENCE OPERATOR
C  DQD		: COALESCENCE DERIVATIVE HEAVISIDE
C  DLAMBDA	: INCREMENT DEFORMATION PLASTIQUE CUMULEE
C  YPS		: VALEUR DE LA FONCTION YPS = (Y/S)**s
C  YPS1		: VALEUR DE LA FONCTION YPS = (Y/S)**(s-1)
C  YV		: VALEUR DE LA FONCTION YV = Ye = 3/2*<Sh>/Seq
C  EYV		: VALEUR DE LA FONCTION Exp[YV] 
C  UES		: VALEUR DE LA FONCTION Exp[YV]/Seq (if <Sh>.neq.0)
C  TDL		: TETA*DLAMBDA
C  UDHM		: VALUE 1/(1-Dh)**m
C  UDHM1	: VALUE 1/(1-Dh)**(m+1)
C  WR		: ISOTROPIC DYNAMIC RECOVERY TERM: WR = DSINH[QS*R]
C  WX		: KINEMATIC DYNAMIC RECOVERY TERM: WX = |DMI:X|
C  UNTRWX	: 1+TETA*(PRD*DLAMBDA+PRS)*WX
C  PHTRWX2	: (2/3)*PH/(1+TETA*(PRD*DLAMBDA+RT(2/3)*PRS)*WX)**2
C  PHTRWX3	: (2/3)*PH/(1+TETA*(PRD*DLAMBDA+RT(2/3)*PRS)*WX)**3
C  HTRWR3	: RT(2/3)*H/(1+TETA*(RD*DLAMBDA+RS)*WR)**3
C  RTDL		: PR1+TETA*DLAMBDA
C  TRSWR	: TETA*RS*WR
C  G		: MODULE DE CISAILLEMENT
C  XK		: COEFFICIENT DE LAME
C  RD		: ISOTROPIC HARDENING DYNAMIC RECOVERY CONSTANT
C  H		: ISOTROPIC HARDENING MODULUS
C  RS		: ISOTROPIC HARDENING STATIC RECOVERY CONSTANT
C  PRD		: KINEMATIC HARDENING DYNAMIC RECOVERY CONSTANT
C  PH		: KINEMATIC HARDENING MODULUS
C  PRS		: KINEMATIC HARDENING STATIC RECOVERY CONSTANT
C  GS		: FONCTION D'ENDOMMAGEMENT DUCTILE
C  PS		: FONCTION D'ENDOMMAGEMENT DUCTILE
C  XDN		: CONSTANTE NUCLEATION XLD*XFN/(XKIC*XLF**(1/3))
C  CTE		: CONSTANTE DE TENSION
C  CCO		: CONSTANTE DE COMPRESSION
C  CTO		: CONSTANTE DE TORSION
C  AV		: VOID GROWTH CONSTANT
C  VM		: VOID GROWTH CONSTANT
C  VP		: VOID GROWTH TEST (VC = 1 if Dh > Vc)
C  AC		: COALESCENCE GROWTH CONSTANT
C  TETA		: TETA-METHODE
C  MELPST	: MELPST = 1 --> TN = ZD/ZDN (RESOLUTION)
C                 SINON TN = ZD/ZDZD (MATRICE CONSISTANTE)
C  NDI		: NOMBRE DE COMPOSANTES DES CONTRAINTES DIRECTES
C  NSHR		: NOMBRE DE COMPOSANTES DES CONTRAINTES TANGENTIELLES
C  NTENS	: TAILLE DU TABLEAU DES CONTRAINTES OU DEFORMATIONS
C		  (NDI + NSHR).
C  M		: COMMENTAIRES
C  MPI		: UNITE DE FICHIER DE SORTIE
C-------------
C  EN SORTIE :
C-------------
C  DHDL		: DERIVEE PARTIELLE DE H PAR RAPPORT A DLAMBDA
C  DHDX		: DERIVEE PARTIELLE DE H PAR RAPPORT A WX
C  DHDR		: DERIVEE PARTIELLE DE H PAR RAPPORT A WR
C  DHDD		: DERIVEE PARTIELLE DE H PAR RAPPORT A ENDO
C  UDEVNN	: TENSEUR (D'ORDRE 4) = Idev -(2/3) N*N
C  CSTRESD	: PRODUCT StresD : [ M : Idev : Mt ]
C-----------------------------------------------------------------------
C=======================================================================
C-----------------------------------------------------------------------
C.1-----  Implicit, External
      IMPLICIT NONE
C.2-----  Parameter,Integer
      INTEGER NCD,NTENS,NDI,NSHR,M,MPI,I,J
C.3-----  Dimension
      DOUBLE PRECISION DNDL(NTENS),DNDX(NTENS),DNDR(NTENS),DNDD(NTENS)
      DOUBLE PRECISION DNDS(NTENS)
      DOUBLE PRECISION ZD(NTENS),TN(NTENS)
      DOUBLE PRECISION STRES1D(NTENS),ALPHA1(NTENS),ALTN(NTENS)
      DOUBLE PRECISION STRESS(6),XI(6),DNDLN(6)
      DOUBLE PRECISION DXIDL(6),DXIDX(6),DXIDR(6),DXIDD(6),DXIDS(6)
      DOUBLE PRECISION H(3)
C-------
      DOUBLE PRECISION DHDL,DHDX,DHDR,DHDD,DHDS,DPDD,DCDD
      DOUBLE PRECISION PHI,DLAMBDA,TDL,UNMD,UNMD1,
     1       WX,PHTRWX,PHTRWX2,G,RDK,RSK,XMS,XMM1,TETA
      DOUBLE PRECISION XIEQV,P,EYV,DEYV,SEFF,DADD
      DOUBLE PRECISION UPSD,DMS2,APHI,AD,XMS2,UNMDM,AMM,SINHA,COSHA,CPHI
      DOUBLE PRECISION D2G,D2GDL,PHRDWX,PHRSWX,PHRDSDL
      DOUBLE PRECISION TGDHX,DLGDHX,TD2GPHD,T2GPHDDL
      DOUBLE PRECISION XIEQV2,SEFF2
      DOUBLE PRECISION DXIEQVDL,DXIEQVDX,DXIEQVDR,DXIEQVDD,DXIEQVDS
      DOUBLE PRECISION DSEFFDL,DSEFFDX,DSEFFDR,DSEFFDD,DSEFFDS
      DOUBLE PRECISION DH1DD,DH2DD,DH3DD,DLCAXI,DLCAY
C-----------------------------------------------------------*
      DOUBLE PRECISION ZERO,ONE,TWO,THREE,HALF,THALF,XNINE
      DOUBLE PRECISION ONETHIRD,TWOTHIRD
C.4-----  Real,Complex,Double precision,Logical,Character
C.5-----  Common
C.6-----  Equivalence
C.7-----  Save,Data
      DATA ZERO,ONE,TWO,THREE,XNINE/0.D0,1.D0,2.D0,3.D0,9.D0/
      DATA HALF,THALF/0.5D0,1.5D0/
C.8-----  Definition de fonctions
C---------------------------------------------------------------------*
C======================================================================
C---------------------------------------------------------------------*
C------- Initialization
      D2G      = TWO*G
      ONETHIRD = ONE/THREE
      TWOTHIRD = TWO/THREE
C---------------------------------------------------------------------*
      IF(NCD.EQ.1 .OR. NCD.EQ.2) THEN
C-------------------------------------------*
C---------  Tensors
         TGDHX  = TETA*(D2G*UNMD+PHTRWX)
         DLGDHX = TGDHX*DLAMBDA
         DO I=1,NTENS
            DNDLN(I)  = TDL*DNDL(I)+TN(I)
            ALTN(I)   = ALPHA1(I)+TDL*TN(I)
            XI(I)     = ZD(I)-DLGDHX*TN(I)
         END DO
C-------------------------------------------*
C---------  Norm |Xi| = Sqrt(3/2 * Xi:Xi)
         XIEQV = XI(1)*XI(1)+XI(2)*XI(2)+XI(3)*XI(3)
         DO I=1+NDI,NTENS
            XIEQV = XIEQV+TWO*XI(I)*XI(I)
         END DO
         XIEQV = DSQRT(THALF*XIEQV)
         SEFF  = DSQRT(H(1)*XIEQV*XIEQV+H(2)*P*P)
C-------------------------------------------*
C---------  Partial Derivative DXIDx = dXi/dx
         D2GDL    = D2G*TDL
         PHRDWX   = PHTRWX2*RDK*WX
         PHRSWX   = PHTRWX2*RSK*WX
         PHRDSDL  = PHTRWX2*(RDK*DLAMBDA+RSK)
         DO I=1,NTENS
            DXIDL(I) = PHRDWX*ALTN(I)-TGDHX*DNDLN(I)
            DXIDX(I) = PHRDSDL*ALTN(I)-DLGDHX*DNDX(I)
            DXIDR(I) = -DLGDHX*DNDR(I)
            DXIDD(I) = -STRES1D(I)+D2GDL*TN(I)-DLGDHX*DNDD(I)
            DXIDS(I) = PHRSWX*ALTN(I)-DLGDHX*DNDS(I)
         END DO
C-------------------------------------------*
C---------  Partial Derivative of DXIEQVDx = dXieq/dx (1st part)
         DXIEQVDL = XI(1)*DXIDL(1)+XI(2)*DXIDL(2)+XI(3)*DXIDL(3)
         DXIEQVDX = XI(1)*DXIDX(1)+XI(2)*DXIDX(2)+XI(3)*DXIDX(3)
         DXIEQVDR = XI(1)*DXIDR(1)+XI(2)*DXIDR(2)+XI(3)*DXIDR(3)
         DXIEQVDD = XI(1)*DXIDD(1)+XI(2)*DXIDD(2)+XI(3)*DXIDD(3)
         DXIEQVDS = XI(1)*DXIDS(1)+XI(2)*DXIDS(2)+XI(3)*DXIDS(3)
         DO I=1+NDI,NTENS
            DXIEQVDL = DXIEQVDL+TWO*XI(I)*DXIDL(I)
            DXIEQVDX = DXIEQVDX+TWO*XI(I)*DXIDX(I)
            DXIEQVDR = DXIEQVDR+TWO*XI(I)*DXIDR(I)
            DXIEQVDD = DXIEQVDD+TWO*XI(I)*DXIDD(I)
            DXIEQVDS = DXIEQVDS+TWO*XI(I)*DXIDS(I)
         END DO
      END IF
C-----------------------------------------------------------*
C-----------------------------------------------------------*
C--------- Choice of Damage Evolution
      IF(NCD.EQ.1) THEN
C-----------------------------------------------------------*
C---------  Cocks [1989]
C-----------------------------------------------------------*
         EYV  = ONETHIRD*P/SEFF
         XMM1 = XMS/(XMS+ONE)
         AD   = TWO*(XMS+ONE)*(ONE+PHI)*UNMD**XMM1
         APHI = XNINE*XMS*PHI*UNMD*EYV/AD
C-----------------------------------------------------------*
         DH1DD  = TWOTHIRD
         DH2DD  = HALF*XMM1/(UNMD1*UNMD1)
         DH3DD  = XMM1*UNMD1**(XMM1-ONE)
         DH1DD  = ZERO
         DH2DD  = ZERO
         DH3DD  = ZERO
C-----------------------------------------------------------*
C---------  Partial Derivative of DSEFFDx = d[1/Seff]/dx
         XIEQV2  = XIEQV*XIEQV
         SEFF2   = SEFF*SEFF
         DSEFFDL = THALF*H(1)*DXIEQVDL/SEFF2
         DSEFFDX = THALF*H(1)*DXIEQVDX/SEFF2
         DSEFFDR = THALF*H(1)*DXIEQVDR/SEFF2
         DSEFFDD = HALF*(DH1DD*XIEQV2+DH2DD*P*P)
     *                     +THALF*H(1)*DXIEQVDD+H(2)*P*DPDD
         DSEFFDD = DSEFFDD/SEFF2
         DSEFFDS = THALF*H(1)*DXIEQVDS/SEFF2
C-------------------------------------------*
C---------  Partial Derivative DHDL = dFphi/dlambda
         DHDL = -APHI*(ONE-DLAMBDA*DSEFFDL)
C-------------------------------------------*
C---------  Partial Derivative DHDK = dFphi/dWx
         DHDX = APHI*DLAMBDA*DSEFFDX
C-------------------------------------------*
C---------  Partial Derivative DHDB = dFphi/dWb
         DHDR = APHI*DLAMBDA*DSEFFDR
C-------------------------------------------*
C---------  Partial Derivative DHDD = dFphi/dphi
         DADD = (ONE+(XMM1-TWO)*PHI)/(PHI*UNMD)-ONE/(ONE+PHI)
         DEYV = ONETHIRD*DPDD/SEFF-EYV*DSEFFDD
         DHDD = APHI*DADD+XNINE*XMS*PHI*UNMD*DEYV/AD
         DHDD = ONE-DLAMBDA*DHDD
C-------------------------------------------*
C---------  Partial Derivative DHDS = dFphi/dls [Eq. 203]
         DHDS = APHI*DLAMBDA*DSEFFDS
C-----------------------------------------------------------*
      ELSE IF(NCD.EQ.2) THEN
C-----------------------------------------------------------*
C---------  Cocks and Ashby [1980]
C-----------------------------------------------------------*
         EYV   = ONETHIRD*P/XIEQV
         XMS2  = TWO*XMS
         UNMDM = UNMD**XMS
         AMM   = TWO*(XMS2-ONE)/(XMS2+ONE)
         APHI  = AMM*EYV
         SINHA = SINH(APHI)
         CPHI  = ONE/UNMDM-UNMD
C-----------------------------------------------------------*
C---------  Partial Derivative of DXIEQVDx = d[1/Xieq]/dx * Xieqv
         XIEQV2   = XIEQV*XIEQV
         DXIEQVDL = -THALF*DXIEQVDL/XIEQV2
         DXIEQVDX = -THALF*DXIEQVDX/XIEQV2
         DXIEQVDR = -THALF*DXIEQVDR/XIEQV2
         DXIEQVDD = -THALF*DXIEQVDD/XIEQV2
         DXIEQVDS = -THALF*DXIEQVDS/XIEQV2
C-------------------------------------------*
         COSHA  = DCOSH(APHI)
         DLCAY  = DLAMBDA*CPHI*COSHA*AMM*EYV
         DLCAXI = DLAMBDA*CPHI*COSHA*AMM/XIEQV
C-------------------------------------------*
C---------  Partial Derivative DHDL = dFphi/dlambda
         DHDL = -CPHI*SINHA-DLCAY*DXIEQVDL
C-------------------------------------------*
C---------  Partial Derivative DHDX = dFphi/dWx 
         DHDX = -DLCAY*DXIEQVDX
C-------------------------------------------*
C---------  Partial Derivative DHDB = dFphi/dWr
         DHDR = -DLCAY*DXIEQVDR
C-------------------------------------------*
C---------  Partial Derivative DHDD = dFphi/dphi
         DCDD = ONE+XMS*UNMD1**(XMS+ONE)
         DHDD = ONE-DLAMBDA*SINHA*DCDD-DLCAY*DXIEQVDD
     *          -ONETHIRD*DLCAXI*DPDD
C-------------------------------------------*
C---------  Partial Derivative DHDB = dFphi/dLs
         DHDS = -DLCAY*DXIEQVDS
         WRITE(MPI,102)'CALDPH DPDD = ',DPDD
C-----------------------------------------------------------*
      ELSE
C-----------------------------------------------------------*
C---------  Plasticity without Damage
C-----------------------------------------------------------*
         DHDL = ZERO
         DHDX = ZERO
         DHDR = ZERO
         DHDD = ZERO
         DHDS = ZERO
C-----------------------------------------------------------*
      END IF
C-----------------------------------------------------------*
C-------  END Partial Derivatives dFgamma/dphi  ----------------------*
C---------------------------------------------------------------------*
C-------  Comment Check Code
      IF(M.GE.4) THEN
         IF(NCD.EQ.1 .OR. NCD.EQ.2) THEN
            WRITE(MPI,102)'STRESS XI'
            DO I=1,NTENS
               WRITE(MPI,101)XI(I)
            END DO
         END IF
         WRITE(MPI,102)'CALDPH DHDL = ',DHDL
         WRITE(MPI,102)'CALDPH DHDX = ',DHDX
         WRITE(MPI,102)'CALDPH DHDR = ',DHDR
         WRITE(MPI,102)'CALDPH DHDD = ',DHDD
         WRITE(MPI,102)'CALDPH DHDS = ',DHDS
      END IF
C-------  END Comment Check Code
C-----------------------------------------------------------------------
C---------   FORMATS  --------------------------------------------------
C-----------------------------------------------------------------------
 101  FORMAT(1X,7(D20.14,2X))
 102  FORMAT(1X,A20,6(D20.14,2X))
C-----------------------------------------------------------------------
C=======================================================================
      RETURN
      END
C
      SUBROUTINE CALDPN(DNDL,DNDX,DNDR,DNDD,DNDS,ZD,ALPHA1,STRES1D,
     *           ZDZD,ZDN,DZDNDL,DZDNDD,DFDR,DFDD,ALPHZD,UNMD,
     *           DLAMBDA,TDL,WX,PHTRWX,PHTRWX2,RDK,RSK,TETA,
     *           NCD,MELPST,NTENS,NDI,NSHR,M,MPI)
C======================================================================
C----------------------------------------------------------------------
C---------------  CALCUL DES DERIVEES PARTIELLES DE N  ----------------
C---------------  PAR RAPPORT DLAMBDA, WX, XR et ENDO  ----------------
C----------------------------------------------------------------------
C======================================================================
C-------------
C  EN ENTREE :
C-------------
C  ZD		: TENSEUR DE DIRECTION NORMALE
C  TN		: TENSEUR NORMAL
C  ESTAR	: TENSEUR DES DEFORMATIONS CONNUS A T+DT
C  ALPHA1	: TENSEUR D'ECROUISSAGE CINEMATIQUE CONNU
C  TETN		: TENSEUR = ESTAR - TDLRD*TN
C  ALTN		: TENSEUR = ALPHA1 + TDLRD*TN
C  PR1		: VARIABLE D'ECROUISSAGE ISOTROPE CONNUE
C  ZDN		: NORME DU TENSEUR DE DIRECTION NORMALE ZD
C  DFDD		: DERIVEE PARTIELLE DE F PAR RAPPORT A ENDO
C  DZDNDL	: DERIVEE PARTIELLE DE ZDN PAR RAPPORT A DLAMBDA
C  ALPHZD	: PRODUIT CONTRACTE ALPHA1 : ZD
C  DM		: TENSEUR D'EFFET D'ENDOMMAGEMENT (Mt)
C  DMI		: TENSEUR D'EFFET D'ENDOMMAGEMENT INVERSE (M(-t))
C  DMDDESTAR	: PRODUIT CONTRACTE IDEV : DMDD : ESTAR
C  DMDDALPH	: PRODUIT CONTRACTE DMDD : ALPHA1
C  DLAMBDA	: INCREMENT DEFORMATION PLASTIQUE CUMULEE
C  YPS		: VALEUR DE LA FONCTION YPS = (Y/S)**s
C  YPS1		: VALEUR DE LA FONCTION YPS = (Y/S)**(s-1)
C  TDL		: TETA*DLAMBDA
C  PHTRWX	: (2/3)*PH/(1+TETA*(PRD*DLAMBDA+RT(2/3)*PRS)*WX)
C  PHTRWX2	: (2/3)*PH/(1+TETA*(PRD*DLAMBDA+RT(2/3)*PRS)*WX)**2
C  G		: MODULE DE CISAILLEMENT
C  XK		: COEFFICIENT DE LAME
C  PRD		: CONSTANTE D'ECROUISSAGE CINEMATIQUE
C  PRS		: CONSTANTE D'ECROUISSAGE CINEMATIQUE (RECOVERY)
C  NCD		: NCD=1 --> COUPLAGE ELASTOPLASTICITE-ENDOMMAGEMENT
C  		  NCD=0 --> ELASTOPLASTICITE SANS ENDOMMAGEMENT
C  TETA		: TETA-METHODE
C  MELPST	: MELPST = 1 --> TN = ZD/ZDN (RESOLUTION)
C                 SINON TN = ZD/ZDZD (MATRICE CONSISTANTE)
C  NDI		: NOMBRE DE COMPOSANTES DES CONTRAINTES DIRECTES
C  NSHR		: NOMBRE DE COMPOSANTES DES CONTRAINTES TANGENTIELLES
C  NTENS	: TAILLE DU TABLEAU DES CONTRAINTES OU DEFORMATIONS
C		  (NDI + NSHR).
C  M		: COMMENTAIRES
C  MPI		: UNITE DE FICHIER DE SORTIE
C-------------
C  EN SORTIE :
C-------------
C  DNDL		: DERIVEE PARTIELLE DE N PAR RAPPORT A DLAMBDA
C  DNDX		: DERIVEE PARTIELLE DE N PAR RAPPORT A WX
C  DNDD		: DERIVEE PARTIELLE DE N PAR RAPPORT A ENDO
C  DNDS		: DERIVEE PARTIELLE DE N PAR RAPPORT A DLs
C  DZDDD	: DERIVEE PARTIELLE DE ZD PAR RAPPORT A ENDO
C  UDEVNN	: TENSEUR (D'ORDRE 4) = Idev -(2/3) N*N
C-----------------------------------------------------------------------
C=======================================================================
C-----------------------------------------------------------------------
C=======================================================================
C.1-----  Implicit, External
      IMPLICIT NONE
C.2-----  Parameter,Integer
      INTEGER NCD,MELPST,NTENS,NDI,NSHR,M,MPI,I,J
C.3-----  Dimension
      DOUBLE PRECISION DNDL(NTENS),DNDX(NTENS),DNDR(NTENS)
      DOUBLE PRECISION DNDD(NTENS),DNDS(NTENS)
      DOUBLE PRECISION ZD(NTENS),ALPHA1(NTENS),STRES1D(NTENS)
C-------
      DOUBLE PRECISION ZDZD,ZDN,DZDNDL,DZDNDD,DFDR,DFDD,ALPHZD,DLAMBDA,
     1       TDL,UNMD,WX,PHTRWX,PHTRWX2,RDK,RSK,TETA
      DOUBLE PRECISION TRDHWX,THDLRP,STRESZD
      DOUBLE PRECISION UZDN,UZDN2,DZDNDL2,DZDNDX2,DZDNDR2,DZDNDD2
      DOUBLE PRECISION UZD,UZDZD,ALPHZD3,TRDHWXZD,TRDHWAZD3,THDLRPZD,
     1       THDLRPZD3,TRSHWX,TRSHWXZD,TRSHWAZD3,DZDZDD3
      DOUBLE PRECISION ZERO,ONE,TWO,THREE,UDEMI,TDEMI,TWOTHIRD
C.4-----  Real,Integer,Complex,Double precision,Logical,Character
C.5-----  Common
C.6-----  Equivalence
C.7-----  Save,Data
      DATA ZERO,ONE,TWO,THREE/0.D0,1.D0,2.D0,3.D0/
      DATA UDEMI,TDEMI/0.5D0,1.5D0/
C.8-----  Definition de fonctions
C---------------------------------------------------------------------*
C======================================================================
C---------------------------------------------------------------------*
C-------  Initialisation
      TWOTHIRD = TWO/THREE
C
C-------  FIN Initialisation
C---------------------------------------------------------------------*
C-------  Calcul des derivees partielles dN/dlambda, dN/dWx et DNDD=dN/dD
C
      TRDHWX = TETA*RDK*WX*PHTRWX2
      THDLRP = TETA*(DLAMBDA*RDK+RSK)*PHTRWX2
C-------
      IF(MELPST.EQ.1) THEN
C-------
         UZDN     = ONE/ZDN
	 UZDN2    = UZDN*UZDN
         DZDNDL2  = DZDNDL*UZDN2
         TRDHWXZD = UZDN*TRDHWX
	 THDLRPZD = UZDN*THDLRP
	 DZDNDX2  = UZDN*THDLRPZD*TDL
	 DZDNDR2  = DFDR*UZDN2
	 DZDNDD2  = DZDNDD*UZDN2
         DO I=1,NTENS
            DNDL(I) = TRDHWXZD*ALPHA1(I)-DZDNDL2*ZD(I)
	    DNDX(I) = THDLRPZD*ALPHA1(I)+DZDNDX2*ZD(I)
	    DNDR(I) = DZDNDR2*ZD(I)
	    DNDD(I) = -UZDN*STRES1D(I)-DZDNDD2*ZD(I)
	    DNDS(I) = ZERO
         END DO
C-------
      ELSE
C-------
         CALL PDTSCA(NDI,NSHR,STRES1D,ZD,STRESZD)
         UZD       = ONE/ZDZD
         UZDZD     = UZD*UZD
	 ALPHZD3   = TWOTHIRD*ALPHZD*UZDZD*UZDZD
         TRDHWXZD  = UZDZD*TRDHWX
         TRDHWAZD3 = TRDHWXZD*ALPHZD3
	 THDLRPZD  = UZDZD*THDLRP
	 THDLRPZD3 = THDLRP*ALPHZD3
	 TRSHWX    = TETA*RSK*WX*PHTRWX2
         TRSHWXZD  = UZDZD*TRSHWX
         TRSHWAZD3 = TRSHWXZD*ALPHZD3
         DZDZDD3   = UZDZD*UZDZD*UZDZD*STRESZD
         DO I=1,NTENS
            DNDL(I) = TRDHWXZD*ALPHA1(I)-TRDHWAZD3*ZD(I)
            DNDX(I) = THDLRPZD*ALPHA1(I)-THDLRPZD3*ZD(I)
	    DNDR(I) = ZERO
	    DNDD(I) = -UZDZD*STRES1D(I)+DZDZDD3*ZD(I)
            DNDS(I) = TRSHWXZD*ALPHA1(I)-TRSHWAZD3*ZD(I)
         END DO
C-------
      END IF
C----------  Tensor DNDD if Ncd = 0
      IF(NCD.EQ.0) CALL ASET(DNDD,ZERO,NTENS)
C
C-------  END Calcul des derivees partielles dN/dlambda, dN/dWx et dN/dD
C---------------------------------------------------------------------*
C---------------------------------------------------------------------*
C-------  Comment Check Code
C
      IF(M.GE.4) THEN
         WRITE(MPI,*)'CALDPN DNDL'
         DO I=1,NTENS
            WRITE(MPI,101)DNDL(I)
         END DO
         WRITE(MPI,*)'CALDPN DNDX'
         DO I=1,NTENS
            WRITE(MPI,101)DNDX(I)
         END DO
	 IF(MELPST.EQ.1) THEN
            WRITE(MPI,*)'CALDPN DNDR'
            DO I=1,NTENS
               WRITE(MPI,101)DNDR(I)
            END DO
	 END IF
         WRITE(MPI,*)'CALDPN DNDD'
         DO I=1,NTENS
            WRITE(MPI,101)DNDD(I)
         END DO
	 IF(MELPST.EQ.0) THEN
	    WRITE(MPI,*)'CALDPN DNDS'
	    DO I=1,NTENS
	       WRITE(MPI,101)DNDS(I)
	    END DO
	 END IF
      END IF
C
C-------  END Comment Check Code
C----------------------------------------------------------------------*
C---------   FORMATS  -------------------------------------------------*
C----------------------------------------------------------------------*
 101  FORMAT(1X,6(D20.14,2X))
 102  FORMAT(1X,A20,6(D20.14,2X))
C----------------------------------------------------------------------*
C=======================================================================
      RETURN
      END
C
      SUBROUTINE CALDSIG(DDSDDE,DNDE,TNDE,SNDDDE,DLDE,DDDE,TN,STRANT,
     *           STRES1D,DLAMBDA,UNMD,UNMD1,G,XK,
     *           NTENS,NDI,M,MPI)
C======================================================================
C----------------------------------------------------------------------
C---------------------  CALCUL DU TENSEUR DDSDDE  ---------------------
C----------------------------------------------------------------------
C======================================================================
C-------------
C  EN ENTREE :
C-------------
C  DNDE		: TENSEUR (D'ORDRE 4) = DN/DE
C  TNDE		: TENSEUR (D'ORDRE 4) = TN*DLDE
C  SNDDDE 	: TENSEUR (D'ORDRE 4) = (XK*EH*I+S1D+2G*DLRD*TN)*DDDE
C  DDDE		: TENSEUR (D'ORDRE 2) = DD/DE
C  STRANT	: TENSEUR DES DEFORMATIONS TOTALES
C  TN		: TENSEUR NORMAL
C  STRES1D	: TENSEUR DES CONTRAINTES CONNU
C  DLAMBDA	: INCREMENT DEFORMATION PLASTIQUE CUMULEE
C  UNMD		: 1-D
C  RCUNMD	: RACINE(1-D)
C  DLRD		: DLAMBDA/RACINE(1-D)
C  G		: MODULE DE CISAILLEMENT
C  XK		: COEFFICIENT DE LAME
C  NDI		: NOMBRE DE COMPOSANTES DES CONTRAINTES DIRECTES
C  NTENS	: NOMBRE TOTAL DE COMPOSANTES D'UN TENSEUR
C  M		: COMMENTAIRES
C  MPI		: UNITE DE FICHIER DE SORTIE
C-------------
C  EN LOCAL  :
C-------------
C  ALDE		: TENSEUR (D'ORDRE 4) = TN*DLDE
C-------------
C  EN SORTIE :
C-------------
C  DDSDDE	: MATRICE TANGENTE CONSISTENTE
C-----------------------------------------------------------------------
C=======================================================================
C-----------------------------------------------------------------------
C=======================================================================
C.1-----  Implicit, External
      IMPLICIT NONE
C.2-----  Parameter,Integer
      INTEGER NTENS,NDI,M,MPI,I,J
C.3-----  Dimension
      DOUBLE PRECISION DDSDDE(NTENS,NTENS),SNDDDE(NTENS,NTENS)
      DOUBLE PRECISION DNDE(NTENS,NTENS),TNDE(NTENS,NTENS)
      DOUBLE PRECISION DLDE(NTENS),DDDE(NTENS),TN(NTENS)
      DOUBLE PRECISION STRANT(NTENS),STRES1D(NTENS),S1DDN(6)
      DOUBLE PRECISION DLAMBDA,UNMD,UNMD1,G,XK
      DOUBLE PRECISION D2G,UNMD2G,UNMDXKG,D2GDL,STRANH,XKSTRANH
      DOUBLE PRECISION UDEMI,DEUX,TROIS
      DOUBLE PRECISION DE(6),SS1(6),SS2(6),SS3(6)
      DOUBLE PRECISION AA1(4,4),AA2(4,4),AA3(4,4)
C.4-----  Real,Complex,Double precision,Logical,Character
C.5-----  Common
C.6-----  Equivalence
C.7-----  Save,Data
      DATA UDEMI,DEUX,TROIS/0.5D0,2.D0,3.D0/
C.8-----  Definition de fonctions
C----------------------------------------------------------------------
C======================================================================
C----------------------------------------------------------------------
C-------  Initialization
C
      D2G     = DEUX*G
      UNMD2G  = UNMD*D2G
      UNMDXKG = UNMD*(XK-D2G/TROIS)
      D2GDL   = D2G*DLAMBDA
      CALL TRACE(STRANT,STRANH,NDI,NTENS)
      XKSTRANH = XK*STRANH
C
C-------  END Initialisation
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C-------  Consistent Tangent DDSDDE
C
      DO I=1,NTENS
         S1DDN(I) = STRES1D(I)-D2GDL*TN(I)
      END DO
      DO I=1,NDI
         S1DDN(I) = S1DDN(I)+XKSTRANH
      END DO
C-------
      CALL PDTTEN(S1DDN,DDDE,SNDDDE,NTENS,NDI)
      CALL PDTTEN(TN,DLDE,TNDE,NTENS,NDI)
C-------
      DO I=1,NTENS
         DO J=1,NTENS
            DDSDDE(I,J) = -SNDDDE(I,J)
     *                    -UNMD2G*(TNDE(I,J)+DLAMBDA*DNDE(I,J))
         END DO
         DDSDDE(I,I) = DDSDDE(I,I)+UNMD2G
      END DO
C-------
      DO I=1,NDI
         DO J=1,NDI
            DDSDDE(I,J) = DDSDDE(I,J)+UNMDXKG
         END DO
      END DO
C------
      DO I=1,NTENS
         DO J=NDI+1,NTENS
            DDSDDE(I,J) = UDEMI*DDSDDE(I,J)
         END DO
      END DO
C
C-------  END Consistent Tangent DDSDDE
C---------------------------------------------------------------------*
C---------------------------------------------------------------------*
C-------  Check Comment Code
C
      IF(M.GE.4) THEN
         WRITE(MPI,*)'CALDISG DDSDDE'
         DO I=1,NTENS
            WRITE(MPI,101)(DDSDDE(I,J),J=1,NTENS)
         END DO
      END IF
C
C-------  END Check Comment Code
C-----------------------------------------------------------------------
C---------   FORMATS  --------------------------------------------------
C-----------------------------------------------------------------------
 101  FORMAT(1X,6(D12.5,2X))
 102  FORMAT(1X,A20,6(D12.5,2X))
C-----------------------------------------------------------------------
C=======================================================================
      RETURN
      END
C
      SUBROUTINE CALFGH(FF,FX,FR,FD,ZDN,ZDZD,TN,ZD,ALPHA1,PHI,
     *           PHIOLD,UNMD,UNMD1,H,P,DLAMBDA,TDL,EYV,WX,DMWX,PHTRWX,
     *           WR,RTDL,HTRWR,G,XMS,NCD,NTENS,NDI,NSHR,M,MPI)
C======================================================================
C----------------------------------------------------------------------
C---------------  CALCUL DES FONCTIONS F, GX, GR ET H  ----------------
C----------------------------------------------------------------------
C======================================================================
C-------------
C  EN ENTREE :
C-------------
C  ENDO 	: VARIABLE D'ENDOMMAGEMENT
C  ZDN		: NORME DU TENSEUR DE DIRECTION NORMALE ZD
C  ZDZD		: NORME DU TENSEUR DE DIRECTION NORMALE ZD
C  ENDO1	: VARIABLE D'ENDOMMAGEMENT CONNUE
C  TN		: PLASTIC UNIT NORMAL TENSOR
C  DMI		: INVERSE DAMAGE EFFECT TENSOR
C  DN		: COMPOSANTES DU VECTEUR ORIENTATION
C  TDL		: TETA*DLAMBDA
C  YPS		: VALEUR DE LA FONCTION YPS = (-Y/S)**s
C  EYV		: FUNCTION Ye = exp[3/2*<Sh>/Seq]
C  UDHM		: VALUE 1/(1-Dh)**m
C  UDHM1	: VALUE 1/(1-Dh)**(m+1)
C  WX		: KINEMATIC DYNAMIC RECOVERY UNKNOWN: WX = |DMI:X|
C  DMWX		: NORME DU TENSEUR Mt:(Alpha1+Tdl*N)
C  PHTRWX	: (2/3)*PH/(1+TETA*(PRD*DLAMBDA+PRS)*WX)
C  WR		: ISOTROPIC DYNAMIC RECOVERY UNKNOWN: WR = R
C  RTDL		: PR1+TETA*DLAMDA
C  HTRWR	: RT(2/3)*H/(1+TETA*(RD*DLAMBDA+RS)*WR)
C  XDN		: CONSTANTE NUCLEATION XLD*XFN/(XKIC*XLF**(1/3))
C  CTE		: CONSTANTE DE TENSION
C  CCO		: CONSTANTE DE COMPRESSION
C  CTO		: CONSTANTE DE TORSION
C  AV		: VOID GROWTH CONSTANT
C  VM		: VOID GROWTH INDICE M
C  VC		: VOID GROWTH CRITICAL THRESHOLD
C  AC		: COALESCENCE CONSTANT
C  FKIC		: FRACTURE TOUGHNESS CONSTANT
C  DSCALE	: LENGTH SCALE PARAMETER
C  DI		: VOID IMPINGEMENT CRITICAL DIAMETER FRACTION
C  DS		: VOID SHEET CRITICAL DIAMETER FRACTION
C  NCD		: NCD=1 --> COUPLAGE ELASTOPLASTICITE-ENDOMMAGEMENT
C  		  NCD=0 --> ELASTOPLASTICITE SANS ENDOMMAGEMENT
C  NTENS	: TAILLE DU TABLEAU DES CONTRAINTES OU DEFORMATIONS
C  NDI		: NOMBRE DE CONTRAINTES DIRECTES (S11,S22,S33)
C  NSHR		: NOMBRE DE CONTRAINTES TANGENTIELLES (S12,S13,S23)
C  M		: COMMENTAIRES
C  MPI		: UNITE DE FICHIER DE SORTIE
C-------------
C  EN SORTIE :
C-------------
C  FF		: FONCTION EQUATION CRITERE
C  FGX		: FONCTION NORME DES CONTRAINTES EFFECTIVES Wx = |Xt|
C  FGR		: FONCTION SINUS HYPERBOLIQUE DE Wr = Sinh(Qs*R)
C  FH		: FONCTION EQUATION TENSORIELLE D'ENDOMMAGEMENT
C  DENS		: DENSITE DE DISTRIBUTION (TENSEUR)
C  ENDOEQ	: NORME DE LA DENSITE DE DISTRIBUTION
C  QD		: VOID COALESCENCE OPERATOR
C  DQD		: COALESCENCE DERIVATIVE HEAVISIDE
C  VP		: VOID GROWTH TEST (VC = 1 if Dh > Vc)
C  VT		: POSITIVE PART <Dh-Vc>  (Dh=1/3*Trace(Phi))
C-----------------------------------------------------------------------
C=======================================================================
C-----------------------------------------------------------------------
C.1-----  Implicit, External
      IMPLICIT NONE
C.2-----  Parameter,Integer
      INTEGER NCD,NTENS,NDI,NSHR,M,MPI,I,J
C.3-----  Dimension
      DOUBLE PRECISION TN(NTENS),ZD(NTENS),ALPHA1(NTENS),H(3)
      DOUBLE PRECISION STRESS(6),ALTN(6),XI(6),DNDLN(6)
C-------
      DOUBLE PRECISION FF,FX,FR,FD,ZDN,ZDZD,PHI,PHIOLD,UNMD,UNMD1,
     1       DLAMBDA,TDL,EYV,WX,DMWX,PHTRWX,WR,RTDL,HTRWR,G,XMS
      DOUBLE PRECISION P,XIEQV,SEFF,D2G,DLGDHX
      DOUBLE PRECISION APHI,CPHI,SINHA,UNMDM,XMS2,AMM,XMM1
      DOUBLE PRECISION ZERO,HALF,ONE,TWO,THREE,THALF,XNINE,ONETHIRD
C.4-----  Real,Complex,Double precision,Logical,Character
C.5-----  Common
C.6-----  Equivalence
C.7-----  Save,Data
      DATA ZERO,ONE,TWO,THREE,XNINE/0.D0,1.D0,2.D0,3.D0,9.D0/
      DATA HALF,THALF/0.5D0,1.5D0/
C.8-----  Definition of Functions
C---------------------------------------------------------------------*
C======================================================================
C---------------------------------------------------------------------*
C-------  Initialization
      D2G      = TWO*G
      ONETHIRD = ONE/THREE
C-------  END Initialization
C---------------------------------------------------------------------*
C---------------------------------------------------------------------*
C-------  Yield Function
      FF = ZDZD-ZDN
C-------  END Yield Function
C---------------------------------------------------------------------*
C-------  Kinematic Hardening Function
      FX = WX-PHTRWX*DMWX
C-------  END Kinematic Hardening Function
C---------------------------------------------------------------------*
C-------  Isotropic Hardening Function
      FR = WR-HTRWR*RTDL
C-------  END Isotropic Hardening Function
C---------------------------------------------------------------------*
      IF(NCD.EQ.1 .OR. NCD.EQ.2) THEN
C-------  Damage Function
C-------------------------------------------*
C---------  Stress Tensors
         DLGDHX = TDL*(D2G*UNMD+PHTRWX)
         DO I=1,NTENS
            ALTN(I)   = ALPHA1(I)+TDL*TN(I)
            XI(I)     = ZD(I)-DLGDHX*TN(I)
         END DO
C-------------------------------------------*
C---------  Norm |Xi| = Sqrt(3/2 * Xi:Xi)
         XIEQV = XI(1)*XI(1)+XI(2)*XI(2)+XI(3)*XI(3)
         DO I=1+NDI,NTENS
            XIEQV = XIEQV+TWO*XI(I)*XI(I)
         END DO
         XIEQV = SQRT(THALF*XIEQV)
C---------  Norm Seff = Sqrt(h1*Xieqv^2+h2*p^2)
         SEFF = SQRT(H(1)*XIEQV*XIEQV+H(2)*P*P)
      END IF
C-------------------------------------------*
      IF(NCD.EQ.1) THEN
         EYV  = ONETHIRD*P/SEFF
         XMM1 = XMS/(XMS+ONE)
         APHI = TWO*(XMS+ONE)*(ONE+PHI)*UNMD**XMM1
         APHI = XNINE*XMS*PHI*UNMD*EYV/APHI
         FD   = PHI-PHIOLD-TDL*APHI
      ELSE IF(NCD.EQ.2) THEN
         EYV   = ONETHIRD*P/XIEQV
         XMS2  = TWO*XMS
         UNMDM = UNMD**XMS
         AMM   = TWO*(XMS2-ONE)/(XMS2+ONE)
         APHI  = AMM*EYV
         SINHA = SINH(APHI)
         CPHI  = ONE/UNMDM-UNMD
         FD    = PHI-PHIOLD-TDL*CPHI*SINHA
      ELSE IF(NCD.EQ.0) THEN
         FD  = ZERO
      END IF
C-------  END Damage Function
C---------------------------------------------------------------------*
C---------------------------------------------------------------------*
C-------  Comment Check Code
      IF(M.GE.2) THEN
         WRITE(MPI,102)'CALFG  FF = ',FF
         WRITE(MPI,102)'CALFG  FX = ',FX
         WRITE(MPI,102)'CALFG  FR = ',FR
         WRITE(MPI,102)'CALFG  FD = ',FD
      END IF
C-------  END Comment Check Code
C----------------------------------------------------------------------*
C---------   FORMATS  -------------------------------------------------*
C----------------------------------------------------------------------*
 101  FORMAT(1X,6(D20.14,2X))
 102  FORMAT(1X,A20,6(D20.14,2X))
C----------------------------------------------------------------------*
C=======================================================================
      RETURN
      END
C
      SUBROUTINE CALFGHI(FF,FX,FR,FD,ZDN,ZDZD,PHI,PHIOLD,UNMD,WX,
     *           DMWX,PHTRWX,WR,RTDL,HTRWR,XMS,XMM1,NCD,M,MPI)
C======================================================================
C----------------------------------------------------------------------
C-------------  CALCUL DES FONCTIONS Fi, GXi, GRi ET Hi  --------------
C----------------------------------------------------------------------
C======================================================================
C-------------
C  EN ENTREE :
C-------------
C  ENDO 	: TENSEUR D'ENDOMMAGEMENT
C  ZDN		: NORME DU TENSEUR DE DIRECTION NORMALE ZD
C  ZDZD		: NORME DU TENSEUR DE DIRECTION NORMALE ZD
C  ENDO1	: TENSEUR D'ENDOMMAGEMENT CONNU
C  WX		: KINEMATIC DYNAMIC RECOVERY TERM: WX = |DMI:X|
C  DMWX		: NORME DU TENSEUR Mt:(ALPHA1+TDL*TN)
C  PHTRWX	: (2/3)*PH/(1+TETA*(PRD*DLAMBDA+PRS)*WX)
C  WR		: ISOTROPIC DYNAMIC RECOVERY TERM: WR = R
C  RTDL		: PR1+TETA*DLAMDA
C  HTRWR	: RT(2/3)*H/(1+TETA*(RD*DLAMBDA+RS)*WR)
C  VC		: VOID GROWTH CRITICAL THRESHOLD
C  FKIC		: FRACTURE TOUGHNESS CONSTANT
C  DSCALE	: LENGTH SCALE PARAMETER
C  DI		: VOID IMPINGEMENT CRITICAL DIAMETER FRACTION
C  DS		: VOID SHEET CRITICAL DIAMETER FRACTION
C  NCD		: NCD=1 --> COUPLAGE ELASTOPLASTICITE-ENDOMMAGEMENT
C  		  NCD=0 --> ELASTOPLASTICITE SANS ENDOMMAGEMENT
C  M		: COMMENTAIRES
C  MPI		: UNITE DE FICHIER DE SORTIE
C-------------
C  EN SORTIE :
C-------------
C  FF		: YIELD FUNCTION
C  FX		: KINMETIC 
C  FR		: ISOTROPIC HARDENING STRESS Wr = R
C  FD		: DAMAGE FUNCTION FH
C-----------------------------------------------------------------------
C=======================================================================
C-----------------------------------------------------------------------
C.1-----  Implicit, External
      IMPLICIT NONE
C.2-----  Parameter,Integer
      INTEGER NCD,M,MPI,I,J
C.3-----  Dimension
      DOUBLE PRECISION FF,FX,FR,FD,ZDN,ZDZD,PHI,PHIOLD,UNMD,WX,DMWX,
     1       PHTRWX,WR,RTDL,HTRWR
      DOUBLE PRECISION XMS,XMS1,XMM1,UPSD
      DOUBLE PRECISION ZERO
C.4-----  Real,Complex,Double precision,Logical,Character
C.5-----  Common
C.6-----  Equivalence
C.7-----  Save,Data
      DATA ZERO/0.D0/
C.8-----  Definition de fonctions
C---------------------------------------------------------------------*
C======================================================================
C---------------------------------------------------------------------*
C-------  Yield Function Ff
C
      FF = ZDZD-ZDN
C
C-------  END Yield Function Ff
C---------------------------------------------------------------------*
C-------  Kinematic Hardening Function Fx
C
      FX = WX-PHTRWX*DMWX
C
C-------  END Kinematic Hardening Function Fx
C---------------------------------------------------------------------*
C-------  Isotropic Hardening Function Fr
C
      FR = WR-HTRWR*RTDL
C
C-------  END Isotropic Hardening Function Fr
C---------------------------------------------------------------------*
C-------  Damage Functions Fh
C
      IF(NCD.EQ.0) THEN
         FD = ZERO
      ELSE
         FD = PHI-PHIOLD
      END IF
C
C-------  END Damage Functions Fh
C---------------------------------------------------------------------*
C-------  Comment Check Code
      IF(M.GE.4) THEN
         WRITE(MPI,102)'CALFGI  FF = ',FF
         WRITE(MPI,102)'CALFGI  FX = ',FX
         WRITE(MPI,102)'CALFGI  FR = ',FR
         WRITE(MPI,102)'CALFGI  FH = ',FD
      END IF
C-------  END Comment Check Code
C----------------------------------------------------------------------*
C---------   FORMATS  -------------------------------------------------*
C----------------------------------------------------------------------*
 101  FORMAT(1X,6(D20.14,2X))
 102  FORMAT(1X,A20,6(D20.14,2X))
C----------------------------------------------------------------------*
C=======================================================================
      RETURN
      END
C
      SUBROUTINE CALZDTN(ZD,TN,ZDN,ZDZD,DSQHMP,RYVSHDL,STRES1D,ESTAR,
     *           ALPHA1,P,PR1,PHI,DPHI,POLD,XKDH,H,TDL,DTLF,
     *           UNMD,DMWX,PHTRWX,RTDL,HTRWR,YIELD,V,PF,G,XK,XMS,
     *           NCD,NTENS,NDI,NSHR,M,MPI)
C======================================================================
C----------------------------------------------------------------------
C--------------------  CALCUL DU TENSEUR NORMAL  ----------------------
C----------------------------------------------------------------------
C======================================================================
C-------------
C  EN ENTREE :
C-------------
C  ESTAR	: TENSEUR DES DEFORMATIONS CONNUS A L'INSTANT T+DT
C  ALPHA1	: TENSEUR DES CONTRAINTES CINEMATIQUES CONNU A T
C  PR1		: VARIABLE D'ECROUISSAGE ISOTROPE CONNU A T
C  XKDEH	: Xk*Tr(d)
C  YIELD	: LIMITE D'ELASTICITE
C  TDL		: TETA*DLAMBDA
C  DTLF		: RT(2/3)*DLAMBDA/PF
C  PHTRWX	: (2/3)*PH/(1+TETA*(PRD*DLAMBDA+RT(2/3)*PRS)*WX)
C  RTDL		: PR1+TETA*DLAMDA
C  HTRWR	: RT(2/3)*H/(1+(TETA*RD*DLAMBDA+RS)*WR)
C  YIELD	: LIMITE D'ELASTICITE
C  V		: CONSTANTE VISCOPLASTIQUE
C  PF		: CONSTANTE VISCOPLASTIQUE
C  G		: MODULE DE CISAILLEMENT
C  PH		: CONSTANTE D'ECROUISSAGE CINEMATIQUE
C  NDI		: NOMBRE DE COMPOSANTES DES CONTRAINTES DIRECTES
C  NSHR		: NOMBRE DE COMPOSANTES DES CONTRAINTES TANGENTIELLES
C  NTENS	: TAILLE DU TABLEAU DES CONTRAINTES OU DEFORMATIONS
C		  (NDI + NSHR)
C  M		: COMMENTAIRES
C  MPI		: UNITE DE FICHIER DE SORTIE
C-------------
C  EN SORTIE :
C-------------
C  ZD		: TENSEUR DE DIRECTION NORMALE
C  TN		: TENSEUR NORMAL
C  ZDN		: NORME DU TENSEUR ZD
C  ZDZD		: NORME DU TENSEUR ZD EGALE A RACINE(2/3*ZD:ZD)
C  P		: HYDROSTATIC STRESS
C  POLD		: HYDROSTATIC STRESS AT TIME T
C  DMWX		: NORME DU TENSEUR Mt:(Alpha1+Tdl*N)
C  STRESS1D	: TENSEUR DEVIATORIQUE DES CONTRAINTES CONNU A T
C  H		: DAMAGE FUNCTIONS h1, h2 and h3
C-----------------------------------------------------------------------
C=======================================================================
C-----------------------------------------------------------------------
C=======================================================================
C.1-----  Implicit, External
      IMPLICIT NONE
C.2-----  Parameter,Integer
      INTEGER NTENS,NDI,NSHR,NCD,M,MPI,I,J
C.3-----  Dimension
      DOUBLE PRECISION ZD(NTENS),TN(NTENS)
      DOUBLE PRECISION ESTAR(NTENS),STRES1D(NTENS),ALPHA1(NTENS)
      DOUBLE PRECISION GESTAR(6),ALTN(6),H(3)
C-------
      DOUBLE PRECISION ZDN,ZDZD,PHI,PR1,TDL,DTLF,DMWX,PHTRWX,RTDL,HTRWR,
     1       UNMD,YIELD,V,PF,G,XK,XMS
      DOUBLE PRECISION DASINHDTLF,HMP,DSQHMP,RYVSHDL
      DOUBLE PRECISION P,POLD,XKDH,D2G,DPHI,XMS1,XMM1,UPSD
      DOUBLE PRECISION ZERO,HALF,ONE,TWO,THREE,TWOTHIRD,SQTT
C.4-----  Real,Complex,Double precision,Logical,Character
C.5-----  Common
C.6-----  Equivalence
C.7-----  Save,Data
      DATA ZERO,HALF,ONE,TWO,THREE/0.D0,0.5D0,1.D0,2.D0,3.D0/
C.8-----  Definition de fonctions
C---------------------------------------------------------------------*
C======================================================================
C---------------------------------------------------------------------*
C-------  Initialisation
      TWOTHIRD = TWO/THREE
      SQTT     = DSQRT(TWOTHIRD)
C-------  FIN  Initialisation
C---------------------------------------------------------------------*
C-------  Hydrostatic Stress
C      P = POLD+UNMD*XKDEH+XK*DPHI
      P = UNMD*XKDH
C-------  END  Hydrostatic Stress
C---------------------------------------------------------------------*
C-------  Sigma1D
      D2G = TWO*G
      CALL TENSDEV(ESTAR,GESTAR,NDI,NSHR)
C-------
      DO I=1,NTENS
         STRES1D(I) = D2G*GESTAR(I)
      END DO
C-------  END  Sigma1D
C---------------------------------------------------------------------*
C-------  Tensor ZD
      DO I=1,NTENS
         ZD(I) = UNMD*STRES1D(I)-PHTRWX*ALPHA1(I)
      END DO
C-------  END  Tensor ZD
C---------------------------------------------------------------------*
C-------  Marin and McDowell Damage Functions H
      XMS1 = XMS+ONE
      XMM1 = ONE/XMS1
      UPSD = ONE+PHI
      IF(NCD.EQ.1) THEN
         H(1) = ONE+TWOTHIRD*PHI
         H(2) = HALF*XMM1*PHI/UPSD
         H(3) = UNMD**XMM1
      ELSE IF(NCD.EQ.2) THEN
         H(1) = ONE
         H(2) = HALF*PHI/UPSD
         H(3) = UNMD
      ELSE
         H(1) = ONE
         H(2) = ZERO
         H(3) = ONE
      END IF
      H(1) = ONE
      H(2) = ZERO
      H(3) = ONE
C-------------------------------------------*
C-------  First Norm of Tensor ZD with the Sine Function
      DASINHDTLF = DLOG(DTLF+DSQRT(ONE+DTLF*DTLF))
      RYVSHDL    = HTRWR*RTDL+YIELD+V*DASINHDTLF
      HMP        = (H(3)*H(3)*RYVSHDL*RYVSHDL-H(2)*P*P)/H(1)
      DSQHMP     = DSQRT(HMP)
      ZDN        = TWOTHIRD*UNMD*DSQHMP+TDL*(UNMD*D2G+PHTRWX)
C-------  Second Norm of Tensor ZD with |Zd*Zd|
      CALL PDTSCA(NDI,NSHR,ZD,ZD,ZDZD)
      ZDZD = DSQRT(TWOTHIRD*ZDZD)
C-------  END  Calcul de la norme du tenseur ZD suivant deux methodes
C---------------------------------------------------------------------*
C-------  Tensor TN
      DO I=1,NTENS
         TN(I) = ZD(I)/ZDN
      END DO
C-------  END  Tensor TN
C---------------------------------------------------------------------*
C-------  Norm DMWX = |Mt:(Alpha1+Tdl*N)|
      DO I=1,NTENS
         ALTN(I)  = ALPHA1(I)+TDL*TN(I)
      END DO
      CALL PDTSCA(NDI,NSHR,ALTN,ALTN,DMWX)
      DMWX = DSQRT(DMWX)
C-------  END  Norm DMWX = |Mt:(Alpha1+Tdl*N)|
C---------------------------------------------------------------------*
C-------  Comment Check Code
C
      IF(M.GE.4) THEN
         WRITE(MPI,102)'CALZDTN ZDZD = ',ZDZD
         WRITE(MPI,102)'CALZDTN ZDN  = ',ZDN
         WRITE(MPI,102)'CALZDTN DMWX = ',DMWX
         WRITE(MPI,*)'CALZDTN ZD'
         DO I=1,NTENS
            WRITE(MPI,101)ZD(I)
         END DO
         WRITE(MPI,*)'CALZDTN TN'
         DO I=1,NTENS
            WRITE(MPI,101)TN(I)
         END DO
         WRITE(MPI,102)'CALZDTN H'
         DO I=1,3
            WRITE(MPI,101)H(I)
         END DO
      END IF
C
C-------  FIN Comment Check Code
C----------------------------------------------------------------------*
C---------   FORMATS  -------------------------------------------------*
C----------------------------------------------------------------------*
 101  FORMAT(1X,6(D20.14,2X))
 102  FORMAT(1X,A20,6(D20.14,2X))
C----------------------------------------------------------------------*
C=======================================================================
      RETURN
      END
C
      SUBROUTINE DPOLR(ROT,F,M,MPI)
C======================================================================
C----------------------------------------------------------------------
C---------------------  POLAR DECOMPOSITION OF F  ---------------------
C----------------------------------------------------------------------
C======================================================================
C-----------
C  ENTREES :
C-----------
C  F		: DEFORMATION GRADIENT
C  M		: PARAMETRE D'IMPRESSION
C  MPI		: UNITE DE FICHIER D'IMPRESSION
C-----------
C  SORTIES :
C-----------
C  ROT		: MATRICE ROTATION PROPRE
C----------------------------------------------------------------------*
C=======================================================================
C----------------------------------------------------------------------*
C.1-----  Precision
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C.2-----  Parameter
      PARAMETER (N=3)
C.3-----  Dimension
      DIMENSION ROT(3,3),F(3,3),TF(3,3),FTF(3,3)
      DIMENSION D(3),V(3,3),VN(3,3),DV1(3,3),TROT(3,3),RTR(3,3)
C.4-----  Real,Integer,Complex,Double precision,Logical,Character
C.5-----  Common
C.6-----  Equivalence
C.7-----  Save,Data
      DATA UDEMI,UN,DEUX/0.5D0,1.D0,2.D0/
C.8-----  Definition de fonctions
C----------------------------------------------------------------------*
C=======================================================================
C----------------------------------------------------------------------*
C-------  Transposition de la matrice F et calcul de TF*F
C
      DO I=1,N
         DO J=1,N
            TF(I,J)=F(J,I)
         END DO
      END DO
      CALL PRDMAT(N,TF,F,FTF)
C
C-------  FIN  Transposition des matrices F0 et F1
C----------------------------------------------------------------------*
C-------  Calcul des valeurs et vecteurs propres de FTF
C
      CALL VALPROP(FTF,D,V,N,N)
C-------  FIN  Calcul des valeurs et vecteurs propres de FTF
C----------------------------------------------------------------------*
C-------  Calcul de l'inverse du tenseur de deformations pure gauche DV
C
C
C-------  Calcul de l'inverse des racines des valeurs propres
C
      RD1=UN/DSQRT(D(1))
      RD2=UN/DSQRT(D(2))
      RD3=UN/DSQRT(D(3))
C
C-------  Calcul des vecteurs propres nI, nII et nIII
C
      VN(1,1)=RD1*(F(1,1)*V(1,1)+F(1,2)*V(2,1)+F(1,3)*V(3,1))
      VN(2,1)=RD1*(F(2,1)*V(1,1)+F(2,2)*V(2,1)+F(2,3)*V(3,1))
      VN(3,1)=RD1*(F(3,1)*V(1,1)+F(3,2)*V(2,1)+F(3,3)*V(3,1))
C-------
      VN(1,2)=RD2*(F(1,1)*V(1,2)+F(1,2)*V(2,2)+F(1,3)*V(3,2))
      VN(2,2)=RD2*(F(2,1)*V(1,2)+F(2,2)*V(2,2)+F(2,3)*V(3,2))
      VN(3,2)=RD2*(F(3,1)*V(1,2)+F(3,2)*V(2,2)+F(3,3)*V(3,2))
C-------
      VN(1,3)=RD3*(F(1,1)*V(1,3)+F(1,2)*V(2,3)+F(1,3)*V(3,3))
      VN(2,3)=RD3*(F(2,1)*V(1,3)+F(2,2)*V(2,3)+F(2,3)*V(3,3))
      VN(3,3)=RD3*(F(3,1)*V(1,3)+F(3,2)*V(2,3)+F(3,3)*V(3,3))
C
C-------  Calcul de DV1 inverse de DV
C
      DV1(1,1)=RD1*VN(1,1)*VN(1,1)+RD2*VN(1,2)*VN(1,2)
     *                               +RD3*VN(1,3)*VN(1,3)
      DV1(1,2)=RD1*VN(1,1)*VN(2,1)+RD2*VN(1,2)*VN(2,2)
     *                               +RD3*VN(1,3)*VN(2,3)
      DV1(1,3)=RD1*VN(1,1)*VN(3,1)+RD2*VN(1,2)*VN(3,2)
     *                               +RD3*VN(1,3)*VN(3,3)
      DV1(2,2)=RD1*VN(2,1)*VN(2,1)+RD2*VN(2,2)*VN(2,2)
     *                               +RD3*VN(2,3)*VN(2,3)
      DV1(2,3)=RD1*VN(2,1)*VN(3,1)+RD2*VN(2,2)*VN(3,2)
     *                               +RD3*VN(2,3)*VN(3,3)
      DV1(3,3)=RD1*VN(3,1)*VN(3,1)+RD2*VN(3,2)*VN(3,2)
     *                               +RD3*VN(3,3)*VN(3,3)
      DV1(2,1)=DV1(1,2)
      DV1(3,1)=DV1(1,3)
      DV1(3,2)=DV1(2,3)
C
C-------  FIN  Calcul de l'inverse du tenseur de deformations pure DV
C----------------------------------------------------------------------*
C-------  Calcul de ROT
C
      CALL PRDMAT(N,DV1,F,ROT)
      CALL TRANSPOSE(N,ROT,TROT)
      CALL PRDMAT(N,ROT,TROT,RTR)
C
C-------  FIN  Calcul de DR
C----------------------------------------------------------------------*
C----------------------------------------------------------------------*
C-------  Comment Check Code
C
      IF(M.GE.4) THEN
         WRITE(MPI,*)'DPOLR ROT'
         DO I=1,N
            WRITE(MPI,101)(ROT(I,J),J=1,N)
         END DO
      END IF
C
C-------  FIN Comment Check Code
C----------------------------------------------------------------------*
C=======================================================================
C----------------------------------------------------------------------*
C---------   FORMATS  -------------------------------------------------*
C----------------------------------------------------------------------*
 101  FORMAT(1X,3(D18.12,2X))
C-----------------------------------------------------------------------
C=======================================================================
      RETURN
      END
C
      SUBROUTINE ELPLASTDO(DDSDDE,STRESS,STRAN,DSTRAN,EPLAS,ALPHA,
     *           XLAMBDA,DLAMBDA,PR,PHI,DLAN,STVIT,YIELD,V,PF,HK,RDK,
     *           RSK,HI,RDI,RSI,DC,XMS,G,XK,TETA,EPSF,EPSX,EPSR,EPSD,
     *           DTIME,LONEWDT,NCD,NPLAN,NATUR,
     *           NTENS,NDI,NSHR,NSTVI,M,MPI)
C======================================================================
C----------------------------------------------------------------------
C---------------  SOLUTION PLASTIQUE AVEC RETOUR RADIAL  --------------
C----------------------------------------------------------------------
C======================================================================
C-----------
C  ENTREES :
C-----------
C  STRAN	: TENSEUR DES DEFORMATIONS TOTALES (ESTRAN)
C  DSTRAN	: INCREMENT DE DEFORMATION TOTALE
C  UNITDEV	: TENSEUR UNITE DEVIATORIQUE Idev=I-1/3(1*1)
C  E		: MODULE D'YOUNG
C  XNU		: COEFFICIENT DE POISSON
C  G		: MODULE DE CISAILLEMENT
C  XK		: COEFFICIENT DE LAME
C  YIELD	: LIMITE D'ELASTICITE
C  V		: CONSTANTE VISCOPLASTIQUE
C  PF		: CONSTANTE VISCOPLASTIQUE
C  RD		: CONSTANTE D'ECROUISSAGE ISOTROPE
C  H		: CONSTANTE D'ECROUISSAGE ISOTROPE
C  RS		: CONSTANTE D'ECROUISSAGE ISOTROPE (RECOVERY)
C  PRD		: CONSTANTE D'ECROUISSAGE CINEMATIQUE
C  PH		: CONSTANTE D'ECROUISSAGE CINEMATIQUE
C  PRS		: CONSTANTE D'ECROUISSAGE CINEMATIQUE (RECOVERY)
C  XIP		: INITIAL POROSITY
C  DC		: VALEUR DE L'ENDOMMAGEMENT CRITIQUE --> RUPTURE
C  XIP		: INITIAL POROSITY
C  XLF		: VOLUME FRACTION
C  XDN		: CONSTANTE NUCLEATION XLD*XFN/(XKIC*XLF**(1/3))
C  GS		: FONCTION D'ENDOMMAGEMENT DUCTILE
C  PS		: FONCTION D'ENDOMMAGEMENT DUCTILE
C  DC		: VALEUR D'ENDOMAGEMMENT CRITIQUE
C  CTE		: CONSTANTE DE TENSION
C  CCO		: CONSTANTE DE COMPRESSION
C  CTO		: CONSTANTE DE TORSION
C  AV		: VOID GROWTH CONSTANT
C  VM		: VOID GROWTH INDICE M
C  VC		: VOID GROWTH CRITICAL THRESHOLD
C  AC		: COALESCENCE CONSTANT
C  DSCALE	: LENGTH SCALE PARAMETER
C  DI		: VOID IMPINGEMENT CRITICAL DIAMETER RATIO
C  DS		: VOID SHEET CRITICAL DIAMETER RATIO
C  TETA		: TETA-METHODE
C  EPSD		: PRECISION DE LA FONCTION DE CHARGE H
C  EPSG		: PRECISION DE LA FONCTION D'ENDOMMAGEMENT G
C  KCT		: KCT=0 --> MATRICE TANGENTE CONTINUE
C  		  KCT=1 --> MATRICE TANGENTE CONSISTENTE
C  NPLAN	: NPLAN=1 --> DEFORMATIONS PLANES
C		  NPLAN=2 --> CONTRAINTES PLANES
C  NATUR	: NATUR=0 --> PETITES DEFORMATIONS
C		  NATUR=1 --> GRANDES DEFORMATIONS
C  IROT		: IROT=0 --> PAS DE ROTATION
C		  IROT=1 --> AVEC ROTATION
C		: IROT=2 --> AVEC DERIVEE COROTATIONNELLE DE JAUMANN
C		  IROT=3 --> AVEC ROTATION D'ABAQUS (ROTSIG)
C  NCD		: NCD=1 --> COUPLAGE ELASTOPLASTICITE-ENDOMMAGEMENT
C  		  NCD=0 --> ELASTOPLASTICITE SANS ENDOMMAGEMENT
C  NDI		: NOMBRE DE COMPOSANTES DES CONTRAINTES DIRECTES
C  NSHR		: NOMBRE DE COMPOSANTES DES CONTRAINTES TANGENTIELLES
C  NTENS	: TAILLE DU TABLEAU DES CONTRAINTES OU DEFORMATIONS
C		  (NDI + NSHR).
C  NSTVI	: TAILLE DU TABLEAU STVIT DE VARIABLES INTERNES
C  M		: COMMENTAIRES (PARAMETRE D'IMPRESSION)
C  MPI		: UNITE DE FICHIER DE SORTIE
C-----------
C  LOCALES :
C-----------
C  ZD		: TENSEUR DE DIRECTION NORMALE
C  TN		: TENSEUR NORMAL
C  STRANT	: TENSEUR DES DEFORMATIONS TOTALES : STRAN+DSTRAN
C  DEDT		: TENSEUR INCREMENT DE DEFORMATION : DTIME*DSTRAN
C  STRES1D	: TENSEUR DEVIATORIQUE DES CONTRAINTES CONNU A T
C  EPLAS1	: TENSEUR DES DEFORMATIONS PLASTIQUES CONNU A T
C  ALPHA1	: TENSEUR DES CONTRAINTES CINEMATIQUES CONNU A T
C  PR1		: VARIABLE D'ECROUISSAGE ISOTROPE CONNU A T
C---------------------------------------
C  VARIABLES DEVANT ETRE REACTUALISEES :
C---------------------------------------
C  STRESS	: TENSEUR DES CONTRAINTES
C  EPLAS	: TENSEUR DES DEFORMATIONS PLASTIQUES
C  ALPHA	: TENSEUR DES CONTRAINTES CINEMATIQUES
C  XLAMBDA	: DEFORMATION PLASTIQUE CUMULEE
C  PR		: VARIABLE D'ECROUISSAGE ISOTROPE
C  ENDO		: TENSEUR D'ENDOMMAGEMENT
C  DM		: TENSEUR D'EFFET D'ENDOMMAGEMENT (Mt)
C  DMI		: TENSEUR D'EFFET D'ENDOMMAGEMENT INVERSE (M(-t))
C-----------
C  SORTIES :
C-----------
C  DDSDDE	: MATRICE TANGENTE CONSISTENTE (DSTRAN/DSTRESS)
C  STVIT	: TABLEAU DES VARIABLE INTERNES (R,X,J2(X),Y,...)
C----------------------------------------------------------------------*
C=======================================================================
C----------------------------------------------------------------------*
C.1-----  Precision
      IMPLICIT NONE
C.2-----  Parameter,Integer
      INTEGER NCD,NPLAN,NATUR,NTENS,NDI,NSHR,NSTVI,M,MPI,I,NITMAX,MELPST
C.3-----  Dimension
      DOUBLE PRECISION STRESS(NTENS),STRAN(NTENS),DSTRAN(NTENS)
      DOUBLE PRECISION DDSDDE(NTENS,NTENS)
      DOUBLE PRECISION EPLAS(NTENS),ALPHA(NTENS)
      DOUBLE PRECISION EPLAS1(6),ALPHA1(6),TN(6),ZD(6),DEDT(6)
      DOUBLE PRECISION STVIT(NSTVI)
      DOUBLE PRECISION GESTAR(6),ESTAR(6),STRES1D(6),STRANT(6)
      DOUBLE PRECISION TETN(6),ALTN(6)
      DOUBLE PRECISION UDEVNN(36),TNTN(36),DZDDD(36)
      DOUBLE PRECISION DNDL(6),DNDX(6),DNDR(6),DNDD(6),DNDS(6)
      DOUBLE PRECISION H(3)
C-------
      DOUBLE PRECISION XLAMBDA,DLAMBDA,PHI,PR,DLAN,YIELD,V,PF,HK,RDK,
     1       RSK,HI,RDI,RSI,DC,XMS,G,XK,TETA,EPSF,EPSX,EPSR,EPSD,DTIME
      DOUBLE PRECISION FF,FX,FR,FD,DZDNDL,DZDNDD,ALPHZD,DFRDL
      DOUBLE PRECISION ZDN,ZDZD,P,POLD,PR1,PHI1,TDL,WR,WX,DMWX,
     1       RTDL,UNMD,UNMD1,TRSWR,UNTRWX,PHTRWX,PHTRWX2,
     2       UNTRWR,HTRWR,HTRWR2
      DOUBLE PRECISION DFDL,DFDX,DFDR,DFDD,DFDS
      DOUBLE PRECISION DXDL,DXDX,DXDR,DXDD,DXDS
      DOUBLE PRECISION DRDL,DRDR,DRDS
      DOUBLE PRECISION DHDL,DHDX,DHDR,DHDD,DHDS
      DOUBLE PRECISION EYV,XKDH,DPDD,DSQHMP,RYVSHDL,XMM1
      DOUBLE PRECISION ZERO,UN,DAMAGE
C.4-----  Real,Complex,Double precision,Logical,Character
      LOGICAL CONVERGENCE,LONEWDT
C.5-----  Common
C.6-----  Equivalence
C.7-----  Save,Data
      DATA ZERO,UN/0.D0,1.D0/
      DATA MELPST/1/
C.8-----  Definition de fonctions
C----------------------------------------------------------------------*
C=======================================================================
C----------------------------------------------------------------------*
C--------------------------   INITIALISATION  -------------------------*
C----------------------------------------------------------------------*
C-------  LoNewdt --> False si D < 1 and LoNewdt --> True si D > 1
      LONEWDT = .FALSE.
C--------- END LoNewdt --> False si D < 1 and LoNewdt --> True si D > 1
C----------------------------------------------------------------------*
C----------------------------------------------------------------------*
C-------  Calcul des constantes E*, e*, alpha1, r1, d1.
C
      CALL KTIMEN(STRESS,STRAN,DSTRAN,EPLAS,ALPHA,P,PR,PHI,XKDH,
     *           DLAMBDA,DLAN,STRANT,DEDT,ESTAR,GESTAR,STRES1D,EPLAS1,
     *           ALPHA1,PR1,PHI1,POLD,WR,WX,XK,G,HK,RDK,RSK,HI,RDI,RSI,
     *           XMS,TETA,DTIME,NCD,NTENS,NDI,NSHR,M,MPI)
C--------- END Calcul des  constantes E*, e*, alpha1, r1, d1.
C----------------------------------------------------------------------*
C----------------------------------------------------------------------*
C---------  Initial non constant repetitive coefficients
C
      CALL CALCOEFFI(ZD,TN,ZDN,ZDZD,PHI,STRES1D,ESTAR,P,XKDH,
     *           DSQHMP,RYVSHDL,ALPHA1,POLD,PR1,PHI1,UNMD,UNMD1,H,
     *           WR,WX,V,PF,G,XK,HK,RDK,RSK,HI,RDI,RSI,XMS,DLAMBDA,TDL,
     *           DMWX,YIELD,UNTRWX,PHTRWX,PHTRWX2,UNTRWR,HTRWR,
     *           HTRWR2,TRSWR,RTDL,TETA,
     *           NITMAX,NCD,NTENS,NDI,NSHR,M,MPI)
C----------- END  Initial non constant repetitive coefficients
C----------------------------------------------------------------------*
C----------------------------------------------------------------------*
C-------  Initial functions Fi, GXi, GRi and Hi
C
      CALL CALFGHI(FF,FX,FR,FD,ZDN,ZDZD,PHI,PHI1,UNMD,WX,
     *           DMWX,PHTRWX,WR,RTDL,HTRWR,XMS,XMM1,NCD,M,MPI)
C
C--------- END Initial functions Fi, GXi, GRi and Hi
C----------------------------------------------------------------------*
C------------------------   END  INITIALISATION  ----------------------*
C----------------------------------------------------------------------*
C----------------------------------------------------------------------*
C
C-------  Solution Z, DLAMBDA and ENDO
C
      CONVERGENCE=.FALSE.
      DOWHILE(.NOT.CONVERGENCE .AND. NITMAX.LT.100)
C
C----------------------------------------------------------------------*
C----------------------------------------------------------------------*
C-------  Partial derivatives of df/dlambda, df/dWx, df/dWr and df/dD
C
      CALL CALDPF(DFDL,DFDX,DFDR,DFDD,DFDS,DRDS,DPDD,DZDNDL,
     *           DZDNDD,ALPHZD,DFRDL,ZD,H,P,PHI,ALPHA1,STRES1D,PR1,
     *           ZDN,ZDZD,DSQHMP,PHTRWX,RYVSHDL,DLAMBDA,TDL,UNMD,
     *           UNMD1,XKDH,WX,WR,PHTRWX2,HTRWR,HTRWR2,TRSWR,RTDL,
     *           RDK,RSK,RDI,RSI,YIELD,G,XK,V,PF,XMS,XMM1,TETA,
     *           NCD,MELPST,NTENS,NDI,NSHR,M,MPI)
C
C----------- END  Partial derivatives of de f
C----------------------------------------------------------------------*
C----------------------------------------------------------------------*
C-------  Partial derivatives of dN/dlambda, dN/dWx and dN/dD
C
      CALL CALDPN(DNDL,DNDX,DNDR,DNDD,DNDS,ZD,ALPHA1,STRES1D,
     *           ZDZD,ZDN,DZDNDL,DZDNDD,DFDR,DFDD,ALPHZD,UNMD,
     *           DLAMBDA,TDL,WX,PHTRWX,PHTRWX2,RDK,RSK,TETA,
     *           NCD,MELPST,NTENS,NDI,NSHR,M,MPI)
C
C----------- END    Partial derivatives of de N
C----------------------------------------------------------------------*
C----------------------------------------------------------------------*
C-------  Partial derivatives of dgX/dlambda, dgX/dWx and dgX/dD
C                                        dgR/dlambda and dgR/dWr
C
      CALL CALDPG(DXDL,DXDX,DXDR,DXDD,DXDS,DRDL,DRDR,
     *           DNDL,DNDX,DNDR,DNDD,DNDS,DFDR,DFRDL,DMWX,TN,ALPHA1,
     *           DLAMBDA,TDL,WX,UNTRWX,PHTRWX,PHTRWX2,RDK,RSK,
     *           TETA,NCD,NTENS,NDI,NSHR,M,MPI)
C
C----------- END  Partial derivatives of de gX and de gR
C----------------------------------------------------------------------*
C----------------------------------------------------------------------*
C-------  Partial derivatives of dh/dlambda, dh/dWx, dh/dWr and dh/dD
C
      CALL CALDPH(DHDL,DHDX,DHDR,DHDD,DHDS,DNDL,DNDX,DNDR,DNDD,
     *           DNDS,DPDD,ZD,TN,P,PHI,STRES1D,ALPHA1,ALTN,DLAMBDA,TDL,
     *           H,UNMD,UNMD1,WX,PHTRWX,PHTRWX2,G,RDK,RSK,XMS,XMM1,
     *           TETA,NCD,NTENS,NDI,NSHR,M,MPI)
C
C----------- END    Partial derivatives of de h
C----------------------------------------------------------------------*
C----------------------------------------------------------------------*
C-------  System solution clambda, cwx, cwr and cd
C
      CALL SOLSYST(DLAMBDA,PHI,WX,WR,FF,FX,FR,FD,
     *           DFDL,DFDX,DFDR,DFDD,DXDL,DXDX,DXDR,DXDD,DRDL,DRDR,
     *           DHDL,DHDX,DHDR,DHDD,NCD,NTENS,NDI,NSHR,M,MPI)
C
C----------- END  System solution clambda, cwx, cwr and cd
C----------------------------------------------------------------------*
C----------------------------------------------------------------------*
C-------  Test: if D > 1 then decrease of loading time increment
      IF(PHI.GT.UN) THEN
         LONEWDT=.TRUE.
         WRITE(MPI,300)
	 WRITE(MPI,*)'SUB ELPLASTDO PHI = ',PHI
         WRITE(MPI,*)'ELPLASTDO NITMAX =',NITMAX
         RETURN
      END IF
C----------- END  Test: if D > 1 then decrease of loading time increment
C----------------------------------------------------------------------*
C----------------------------------------------------------------------*
C---------  Update of non constant repetitive coefficients
C
      CALL CALCOEFF(ZD,TN,ZDN,ZDZD,STRES1D,ESTAR,P,XKDH,DSQHMP,
     *           RYVSHDL,POLD,ALPHA1,PR1,PHI,PHI1,DLAMBDA,H,UNMD,
     *           UNMD1,WR,WX,DMWX,YIELD,V,PF,G,XK,HK,RDK,RSK,XMS,HI,
     *           RDI,RSI,TDL,UNTRWX,PHTRWX,PHTRWX2,UNTRWR,HTRWR,HTRWR2,
     *           TRSWR,RTDL,TETA,NCD,NITMAX,NTENS,NDI,NSHR,M,MPI)
C----------- END  Update of non constant repetitive coefficients
C----------------------------------------------------------------------*
C----------------------------------------------------------------------*
C-------  Functions F, FX, FR and FH
C
      CALL CALFGH(FF,FX,FR,FD,ZDN,ZDZD,TN,ZD,ALPHA1,PHI,
     *           PHI1,UNMD,UNMD1,H,P,DLAMBDA,TDL,EYV,WX,DMWX,PHTRWX,
     *           WR,RTDL,HTRWR,G,XMS,NCD,NTENS,NDI,NSHR,M,MPI)
C----------- END  Functions F, GX, GR and H
C----------------------------------------------------------------------*
C-------  Verifivation of convergence criteria
C
      CALL VERIFGENCE(CONVERGENCE,FF,FX,FR,FD,
     *           EPSF,EPSX,EPSR,EPSD,NTENS,NDI,NSHR,M,MPI)
C----------- END  Verifivation of convergence criteria
C----------------------------------------------------------------------*
C----------------------------------------------------------------------*
C      PAUSE
      END DO
      IF(M.GE.1) WRITE(MPI,400) NITMAX
C----------------------------------------------------------------------*
C----------------------------------------------------------------------*
C-------  Test: if NITMAX >100 then return and restart
      IF(NITMAX.GE.100) THEN
         LONEWDT=.TRUE.
         WRITE(MPI,200)
         RETURN
      END IF
C-------  END  Test: if NITMAX >100
C
C--------- END Solution Z, DLAMBDA and ENDO
C----------------------------------------------------------------------*
C----------------------------------------------------------------------*
C-------  Tensors EPLAS(n+1), ALPHA(n+1), PR(n+1) and STRESS(n+1)
C
      CALL SOLUTION(STRESS,EPLAS,ALPHA,XLAMBDA,PR,STRANT,PHI,
     *           ZD,TN,STRES1D,EPLAS1,ALPHA1,PR1,DLAMBDA,DLAN,TDL,UNMD,
     *           UNTRWX,UNTRWR,RTDL,EYV,G,XK,HK,HI,XMS,TETA,DTIME,
     *           STVIT,NITMAX,NSTVI,NTENS,NDI,NSHR,NCD,M,MPI)
C
C-------  END  Tensors EPLAS(n+1), ALPHA(n+1), PR(n+1) and STRESS(n+1)
C----------------------------------------------------------------------*
C----------------------------------------------------------------------*
C-------  Test: if D > DC then D = 1, STRESS and DDSDDE equal to zero
C
      IF(PHI.GT.DC) THEN
         PHI = UN
         CALL ASET(STRESS,ZERO,NTENS)
         CALL ASET(DDSDDE,ZERO,NTENS*NTENS)
         CALL ASET(STVIT,ZERO,NTENS+2)
         RETURN
      END IF
C
C----------- END  Test: if D > DC then D = 1, STRESS and DDSDDE equal to zero
C----------------------------------------------------------------------*
C----------------------------------------------------------------------*
C-------  Tangent consistent matrix J=DDSDDE
C
      CALL JDDSDDE(DDSDDE,TN,ZD,STRANT,STRES1D,ESTAR,ALPHA1,TETN,
     *           ALTN,DEDT,ZDN,ZDZD,P,DSQHMP,RYVSHDL,EYV,PHI,H,PR1,
     *           XKDH,DLAMBDA,TDL,UNMD,UNMD1,DMWX,WX,WR,UNTRWX,PHTRWX,
     *           PHTRWX2,HTRWR,HTRWR2,RTDL,TRSWR,YIELD,V,PF,G,XK,RDK,
     *           RSK,RDI,RSI,XMS,XMM1,TETA,NCD,NTENS,NDI,NSHR,M,MPI)
C
C-------  END  Tangent consistent matrix J=DDSDDE
C----------------------------------------------------------------------*
C----------------------------------------------------------------------*
C----------------------------------------------------------------------*
C-------  Comment Check Code
C-------  END Comment Check Code
C----------------------------------------------------------------------*
C=======================================================================
C----------------------------------------------------------------------*
C---------   FORMATS  --------------------------------------------------
C----------------------------------------------------------------------*
 101  FORMAT(1X,3(D15.8,2X))
 102  FORMAT(1X,A20,6(D20.14,2X))
 200  FORMAT(1X,'SUB ELPLASTDO: MAXIMUM ITERATION',1X,
     *          'NUMBER HAS BEEN REACHED')
 300  FORMAT(1X,'SUB ELPLASTDO: DAMAGE COMPONENT GREATER THAN 1')
 400  FORMAT(1X,'*----------------------------------------------*'/
     *       1X,'| *----*  REACHED ITERATION NUMBER: ',I3,'  *----* |',/
     *       1X,'*----------------------------------------------*')
C----------------------------------------------------------------------*
C=======================================================================
      RETURN
      END
C
      SUBROUTINE JDDSDDE(DDSDDE,TN,ZD,STRANT,STRES1D,ESTAR,ALPHA1,TETN,
     *           ALTN,DEDT,ZDN,ZDZD,P,DSQHMP,RYVSHDL,EYV,PHI,H,PR1,
     *           XKDH,DLAMBDA,TDL,UNMD,UNMD1,DMWX,WX,WR,UNTRWX,PHTRWX,
     *           PHTRWX2,HTRWR,HTRWR2,RTDL,TRSWR,YIELD,V,PF,G,XK,RDK,
     *           RSK,RDI,RSI,XMS,XMM1,TETA,NCD,NTENS,NDI,NSHR,M,MPI)
C=======================================================================
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-------------  CALCUL DE LA MATRICE TANGENTE CONSISTENTE  -------------
C--------------------  J = DDSDDE = DSTRESS/DSTRAN  --------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C=======================================================================
C-------------
C  EN ENTREE :
C-------------
C  TN		: TENSEUR NORMAL
C  ZD		: TENSEUR DE DIRECTION NORMALE
C  STRESS	: TENSEUR DES CONTRAINTES
C  STRANT	: TENSEUR DES DEFORMATIONS TOTALES
C  DEDT		: TENSEUR INCREMENT DE DEFORMATION : DTIME*DSTRAN
C  ESTAR	: TENSEUR e* DES DEFORMATIONS CONNUES
C  ENDO		: TENSEUR D'ENDOMMAGEMENT
C  ALPHA1	: TENSEUR D'ECROUISSAGE CINEMATIQUE CONNUE
C  TETN		: TENSEUR = ESTAR - TDLRD*TN
C  ALTN		: TENSEUR = ALPHA1 + TDLRD*TN
C  QD		: VOID COALESCENCE OPERATOR
C  DQD		: COALESCENCE DERIVATIVE HEAVISIDE
C  PLTN		: PRODUIT P:|dEp/dl|
C  DM		: TENSEUR D'EFFET D'ENDOMMAGEMENT (Mt)
C  DMI		: TENSEUR D'EFFET D'ENDOMMAGEMENT INVERSE (M(-t))
C  DN		: COMPONANTES DU VECTEUR ORIENTATION (NTENS,NTENS)
C  DENS		: DENSITE DE DISTRIBUTION DU DOMMAGE
C  UNITDEV	: TENSEUR UNITE DEVIATORIQUE
C  ZDN		: NORME DU TENSEUR ZD (1ERE METHODE)
C  ZDZD		: NORME DU TENSEUR ZD = RACINE (2/3*ZD:ZD)
C  PR1		: VARIABLE D'ECROUISSAGE ISOTROPE CONNUE
C  DLAMBDA	: VITESSE DE DEFORMATION PLASTIQUE CUMULEE
C  YPS		: VALEUR DE LA FONCTION YPS = (-Y/S)**s
C  YPS1		: VALEUR DE LA FONCTION YPS = (-Y/S)**(s-1)
C  YV		: VALEUR DE LA FONCTION YV = Ye = 3/2*<Sh>/Seq
C  EYV		: VALEUR DE LA FONCTION Exp[YV]
C  UES		: VALEUR DE LA FONCTION 1/Seq (if <Sh>.neq.0)
C  TDL		: TETA*DLAMBDA
C  UDHM		: VALUE 1/(1-Dh)**m
C  UDHM1	: VALUE 1/(1-Dh)**(m+1)
C  WR		: ISOTROPIC DYNAMIC RECOVERY TERM: WR = DSINH[QS*R]
C  WX		: KINEMATIC DYNAMIC RECOVERY TERM: WX = |DMI:X|
C  UNTRWX	: 1+TETA*(PRD*DLAMBDA+PRS)*WX
C  PHTRWX	: (2/3)*PH/(1+TETA*(PRD*DLAMBDA+PRS)*WX)
C  PHTRWX2	: (2/3)*PH/(1+TETA*(PRD*DLAMBDA+RT(2/3)*PRS)*WX)**2
C  HTRWR3	: RT(2/3)*H/(1+TETA*(RD*DLAMBDA+RS)*WR)**3
C  RTDL		: PR1+TETA*DLAMBDA
C  TRSWR	: TETA*RS*WR
C  G		: MODULE DE CISAILLEMENT
C  XK		: COEFFICIENT DE LAME
C  RD		: CONSTANTE D'ECROUISSAGE ISOTROPE
C  H		: CONSTANTE D'ECROUISSAGE ISOTROPE
C  RS		: CONSTANTE D'ECROUISSAGE ISOTROPE (RECOVERY)
C  PRD		: CONSTANTE D'ECROUISSAGE CINEMATIQUE
C  PH		: CONSTANTE D'ECROUISSAGE CINEMATIQUE
C  PRS		: CONSTANTE D'ECROUISSAGE CINEMATIQUE (RECOVERY)
C  GS		: FONCTION D'ENDOMMAGEMENT DUCTILE
C  PS		: FONCTION D'ENDOMMAGEMENT DUCTILE
C  XDN		: CONSTANTE NUCLEATION XLD*XFN/(XKIC*XLF**(1/3))
C  AV		: VOID GROWTH CONSTANT
C  VM		: VOID GROWTH CONSTANT
C  VC		: VOID GROWTH TEST (VC = 1 if Dh > Vc)
C  AC		: COALESCENCE CONSTANT
C  TETA		: TETA-METHODE
C  NCD		: NCD=1 --> COUPLAGE ELASTOPLASTICITE-ENDOMMAGEMENT
C  		  KCT=1 --> MATRICE TANGENTE CONSISTENTE
C  NDI		: NOMBRE DE COMPOSANTES DES CONTRAINTES DIRECTES
C  NSHR		: NOMBRE DE COMPOSANTES DES CONTRAINTES TANGENTIELLES
C  NTENS	: TAILLE DU TABLEAU DES CONTRAINTES OU DEFORMATIONS
C		  (NDI + NSHR).
C  M		: COMMENTAIRES
C  MPI		: UNITE DE FICHIER DE SORTIE
C-------------
C  EN SORTIE :
C-------------
C  DDSDDE	: MATRICE TANGENTE CONSISTANTE
C----------------------------------------------------------------------*
C=======================================================================
C----------------------------------------------------------------------*
C.1-----  Implicit, External
      IMPLICIT NONE
C.2-----  Parameter,Integer
      INTEGER NCD,NTENS,NDI,NSHR,M,MPI,MELPST,I,J
C.3-----  Dimension
      DOUBLE PRECISION DDSDDE(NTENS,NTENS)
      DOUBLE PRECISION TN(NTENS),ZD(NTENS),STRANT(NTENS)
      DOUBLE PRECISION STRES1D(6),ESTAR(NTENS),ALPHA1(NTENS)
      DOUBLE PRECISION TETN(NTENS),ALTN(NTENS),DEDT(NTENS)
      DOUBLE PRECISION DLDE(6),DDDE(6)
      DOUBLE PRECISION DNDL(6),DNDX(6),DNDR(6),DNDD(36),DNDS(6)
      DOUBLE PRECISION TNDE(36),DPNDE(36),DNDE(36),ALDE(36),SNDDDE(36)
      DOUBLE PRECISION UNN(36),UZDNN(36),UDEVNN(36),TNTN(36),DZDDD(36)
      DOUBLE PRECISION DZDE(36),DXIDE(36)
      DOUBLE PRECISION H(3)
C--------
      DOUBLE PRECISION ZDN,ZDZD,DZDNDL,DZDNDD,ALPHZD,DFRDL,P,XIEQV,XKDH
      DOUBLE PRECISION DFDL,DFDX,DFDR,DFDD,DFDS
      DOUBLE PRECISION DXDL,DXDX,DXDR,DXDD,DXDS
      DOUBLE PRECISION DRDL,DRDR,DRDS
      DOUBLE PRECISION DHDL,DHDX,DHDR,DHDD,DHDS
      DOUBLE PRECISION DLAMBDA,PHI,PR1,TDL,UNMD,UNMD1,DMWX,
     1       WX,WR,UNTRWX,PHTRWX,PHTRWX2,HTRWR,HTRWR2,RTDL,TRSWR,
     2       YIELD,V,PF,G,XK,RDK,RSK,RDI,RSI,XMS,TETA
      DOUBLE PRECISION EYV,XKDEH,DPDD,DSQHMP,RYVSHDL,XMM1
      DOUBLE PRECISION ZERO,UN,DEUX,TROIS,UDEMI
C.4-----  Real,Complex,Double precision,Logical,Character
C.5-----  Common
C.6-----  Equivalence
C.7-----  Save,Data
      DATA ZERO,UN,DEUX,TROIS,UDEMI/0.D0,1.D0,2.D0,3.D0,0.5D0/
      DATA MELPST/0/
C.8-----  Definition of fonctions
C----------------------------------------------------------------------*
C=======================================================================
C----------------------------------------------------------------------*
C----------------------------------------------------------------------*
C-------  Partial derivatives df/dlambda, df/dWx, df/dWr, df/dD and df/dLs
C
      CALL CALDPF(DFDL,DFDX,DFDR,DFDD,DFDS,DRDS,DPDD,DZDNDL,
     *           DZDNDD,ALPHZD,DFRDL,ZD,H,P,PHI,ALPHA1,STRES1D,PR1,
     *           ZDN,ZDZD,DSQHMP,PHTRWX,RYVSHDL,DLAMBDA,TDL,UNMD,
     *           UNMD1,XKDH,WX,WR,PHTRWX2,HTRWR,HTRWR2,TRSWR,RTDL,
     *           RDK,RSK,RDI,RSI,YIELD,G,XK,V,PF,XMS,XMM1,TETA,
     *           NCD,MELPST,NTENS,NDI,NSHR,M,MPI)
C
C----------- END  Partial derivatives of f
C----------------------------------------------------------------------*
C----------------------------------------------------------------------*
C-------  Partial derivatives dN/dlambda, dN/dWx, dN/dWr, dN/dD and dN/ds
C
      CALL CALDPN(DNDL,DNDX,DNDR,DNDD,DNDS,ZD,ALPHA1,STRES1D,
     *           ZDZD,ZDN,DZDNDL,DZDNDD,DFDR,DFDD,ALPHZD,UNMD,
     *           DLAMBDA,TDL,WX,PHTRWX,PHTRWX2,RDK,RSK,TETA,
     *           NCD,MELPST,NTENS,NDI,NSHR,M,MPI)
C
C----------- END Partial derivatives of N
C----------------------------------------------------------------------*
C----------------------------------------------------------------------*
C-------  Partial derivatives dgX/dlambda, dgX/dWx and dgX/dD
C                                        dgR/dlambda and dgR/dWr
C
      CALL CALDPG(DXDL,DXDX,DXDR,DXDD,DXDS,DRDL,DRDR,
     *           DNDL,DNDX,DNDR,DNDD,DNDS,DFDR,DFRDL,DMWX,TN,ALPHA1,
     *           DLAMBDA,TDL,WX,UNTRWX,PHTRWX,PHTRWX2,RDK,RSK,
     *           TETA,NCD,NTENS,NDI,NSHR,M,MPI)
C
C----------- END  Partial derivatives of gX and of gR
C----------------------------------------------------------------------*
C----------------------------------------------------------------------*
C-------  Partial derivatives dh/dlambda, dh/dWx, dh/dWr and dh/dD
C
      CALL CALDPH(DHDL,DHDX,DHDR,DHDD,DHDS,DNDL,DNDX,DNDR,DNDD,
     *           DNDS,DPDD,ZD,TN,P,PHI,STRES1D,ALPHA1,ALTN,DLAMBDA,TDL,
     *           H,UNMD,UNMD1,WX,PHTRWX,PHTRWX2,G,RDK,RSK,XMS,XMM1,
     *           TETA,NCD,NTENS,NDI,NSHR,M,MPI)
C
C----------- END  Partial derivatives of h
C----------------------------------------------------------------------*
C----------------------------------------------------------------------*
C-------  Partial derivatives in respect with  E = STRANT
C
      CALL CALDDE(DLDE,DDDE,DNDE,UNN,TNTN,UZDNN,ALDE,DPNDE,DZDE,
     *           DXIDE,ESTAR,ALPHA1,ZD,TN,TETN,ALTN,PHI,ZDZD,DFDL,DFDX,
     *           DFDR,DFDD,DXDL,DXDX,DXDD,DRDL,DRDR,DHDL,DHDX,DHDR,
     *           DHDD,DFDS,DXDS,DRDS,DHDS,DEDT,P,EYV,DLAMBDA,TDL,UNMD,
     *           UNMD1,DMWX,WX,PHTRWX,PHTRWX2,H,G,XK,RDK,RSK,RSI,XMS,
     *           TETA,NCD,NTENS,NDI,NSHR,M,MPI)
C
C-------- END  Partial derivatives in respect with  E = STRANT
C----------------------------------------------------------------------*
C----------------------------------------------------------------------*
C-------  Tangent consistent matrix DDSDDE
C
      CALL CALDSIG(DDSDDE,DNDE,TNDE,SNDDDE,DLDE,DDDE,TN,STRANT,
     *           STRES1D,DLAMBDA,UNMD,UNMD1,G,XK,
     *           NTENS,NDI,M,MPI)
C
C-------  END  Tangent consistent matrix DDSDDE
C----------------------------------------------------------------------*
C----------------------------------------------------------------------*
C-------  Comment Check Code
C
      IF(M.GE.4) THEN
         WRITE(MPI,*)'JDDSDDE DDSDDE'
         DO I=1,NTENS
            WRITE(MPI,101)(DDSDDE(I,J),J=1,NTENS)
         END DO
      END IF
C
C-------  END Comment Check Code
C-----------------------------------------------------------------------*
C---------   FORMATS  --------------------------------------------------*
C-----------------------------------------------------------------------*
 101  FORMAT(1X,6(D20.14,2X))
 102  FORMAT(1X,A20,6(D12.6,2X))
C-----------------------------------------------------------------------*
C========================================================================
      RETURN
      END
C
      SUBROUTINE KTIMEN(STRESS,STRAN,DSTRAN,EPLAS,ALPHA,P,PR,PHI,XKDH,
     *           DLAMBDA,DLAN,STRANT,DEDT,ESTAR,GESTAR,STRES1D,EPLAS1,
     *           ALPHA1,PR1,PHI1,POLD,WR,WX,XK,G,HK,RDK,RSK,HI,RDI,RSI,
     *           XMS,TETA,DTIME,NCD,NTENS,NDI,NSHR,M,MPI)
C=======================================================================
C-----------------------------------------------------------------------
C-------------  CALCUL DES TERMES E*, e*, ALPHA1, r1, D1  --------------
C-----------------------------------------------------------------------
C=======================================================================
C-------------
C  EN ENTREE :
C-------------
C  STRESS	: TENSEUR DES CONTRAINTES
C  STRAN	: DEFORMATION TOTALE
C  EPLAS	: TENSEUR DES DEFORMATIONS PLASTIQUES
C  ALPHA	: TENSEUR DES CONTRAINTES CINEMATIQUES
C  PR		: VARIABLE D'ECROUISSAGE ISOTROPE
C  PHI		: VARIABLE D'ENDOMMAGEMENT
C  DLAN		: INCREMENT DEFORMATION PLASTIQUE CUMULEE PRECEDENT
C  NTENS	: LONGUEUR DES TABLEAU DSTRAN, STRESS, ...
C  G		: MODULE DE CISAILLEMENT
C  RD		: CONSTANTE D'ECROUISSAGE ISOTROPE
C  H		: CONSTANTE D'ECROUISSAGE ISOTROPE
C  RS		: CONSTANTE D'ECROUISSAGE ISOTROPE (RECOVERY)
C  PRD		: CONSTANTE D'ECROUISSAGE CINEMATIQUE
C  PH		: CONSTANTE D'ECROUISSAGE CINEMATIQUE
C  PRS		: CONSTANTE D'ECROUISSAGE CINEMATIQUE (RECOVERY)
C  XDN		: CONSTANTE XLD*XFN/(XKIC*XLF**(1/3))
C  GS		: FONCTION D'ENDOMMAGEMENT DUCTILE
C  PS		: FONCTION D'ENDOMMAGEMENT DUCTILE
C  DC		: VALEUR D'ENDOMAGEMMENT CRITIQUE
C  CTE		: CONSTANTE DE TENSION
C  CCO		: CONSTANTE DE COMPRESSION
C  CTO		: CONSTANTE DE TORSION
C  TETA		: TETA-METHODE
C  DTIME	: INCREMENT DE TEMPS
C  NCD		: NCD=1 --> COUPLAGE ELASTOPLASTICITE-ENDOMMAGEMENT
C  		  NCD=0 --> ELASTOPLASTICITE SANS ENDOMMAGEMENT
C  NDI		: NOMBRE DE COMPOSANTES DES CONTRAINTES DIRECTES
C  NSHR		: NOMBRE DE COMPOSANTES DES CONTRAINTES TANGENTIELLES
C  NTENS	: TAILLE DU TABLEAU DES CONTRAINTES OU DEFORMATIONS
C		  (NDI + NSHR).
C  M		: COMMENTAIRES
C  MPI		: UNITE DE FICHIER DE SORTIE
C-------------
C  EN SORTIE :
C-------------
C  STRANT	: TENSEUR DES DEFORMATIONS TOTALES : STRAN+DSTRAN
C  DEDT		: TENSEUR INCREMENT DE DEFORMATION : DSTRAN
C  ESTAR	: TENSEUR DES DEFORMATIONS CONNUS A T+DT
C  GESTAR	: TENSEUR DEVIATORIQUE DE ESTAR
C  STRES1D	: TENSEUR DEVIATORIQUE DES CONTRAINTES CONNU A T
C  EPLAS1	: TENSEUR DES DEFORMATIONS PLASTIQUES CONNUS A T+DT
C  ALPHA1	: TENSEUR DES CONTRAINTES CINEMATIQUES CONNUS A T+DT
C  PR1		: VARIABLE D'ECROUISSAGE ISOTROPE CONNUS A T+DT
C  PHI1		: TENSEUR D'ENDOMMAGEMENT CONNUS A T+DT
C  P		: HYDROSTATIC STRESS
C  POLD		: HYDROSTATIC STRESS AT TIME T
C  PLTN		: PRODUIT P:|dEp/dl|
C  DN		: COMPONANTES DU VECTEUR ORIENTATION
C  DLAMBDA	: INCREMENT DE DEFORMATION PLASTIQUE EQUIVALENTE
C  WR		: ISOTROPIC DYNAMIC RECOVERY TERM: WR = R
C  WX		: KINEMATIC DYNAMIC RECOVERY TERM: WX = |DMI:X|
C-----------------------------------------------------------------------
C=======================================================================
C-----------------------------------------------------------------------
C=======================================================================
C.1-----  Implicit, External
      IMPLICIT NONE
C.2-----  Parameter,Integer
      INTEGER NCD,NTENS,NDI,NSHR,M,MPI,I
C.3-----  Dimension
      DOUBLE PRECISION STRESS(NTENS),STRES1D(NTENS)
      DOUBLE PRECISION STRAN(NTENS),DSTRAN(NTENS),DEDT(NTENS)
      DOUBLE PRECISION EPLAS(NTENS),ALPHA(NTENS)
      DOUBLE PRECISION EPLAS1(NTENS),ALPHA1(NTENS)
      DOUBLE PRECISION STRANT(NTENS),ESTAR(NTENS),GESTAR(NTENS)
      DOUBLE PRECISION Z(6),TN(6),STRESSDV(6),ALPHADV(6)
      DOUBLE PRECISION EELAS(6),DEEELAS(6)
C-------
      DOUBLE PRECISION PHI,PR,DLAMBDA,DLAN,PHI1,PR1,WR,WX,XK,G,
     1       HK,RDK,RSK,HI,RDI,RSI,XMS,TETA,DTIME
      DOUBLE PRECISION ZZ,TDZZ,P,POLD,R,RWR,RDLSWX,XKDEH,XKDH,D2G
      DOUBLE PRECISION UNMTETA,UNMTETADT,UNMTETADTDL,UNMD,UNMDHK
      DOUBLE PRECISION ZERO,ONE,TWO,THREE,HALF,THALF,ONETHIRD
C.4-----  Real,Complex,Double precision,Logical,Character
C.5-----  Common
C.6-----  Equivalence
C.7-----  Save,Data
      DATA ZERO,ONE,TWO,THREE/0.D0,1.D0,2.D0,3D0/
      DATA HALF,THALF/0.5D0,1.5D0/
C.8-----  Definition de fonctions
C---------------------------------------------------------------------*
C======================================================================
C---------------------------------------------------------------------*
C-------  Initialisation
      ONETHIRD     = ONE/THREE
      UNMTETA     = ONE-TETA
      UNMTETADT   = UNMTETA*DTIME
      UNMTETADTDL = UNMTETADT*DLAN
C-------  FIN Initialisation
C---------------------------------------------------------------------*
C-------  Prediction Increment de deformation plastique cumulee
      IF(TETA.EQ.ONE) THEN
         DLAMBDA = DLAN*DTIME
      ELSE
         DLAMBDA = ZERO
      END IF
C-------  END  Prediction Increment de deformation plastique cumulee
C---------------------------------------------------------------------*
C-------  Total Strain Tensor STRANT = STRAN + DSTRAN
C
      DO I=1,NDI
         DEDT(I)   = DSTRAN(I)
         STRANT(I) = STRAN(I)+DEDT(I)
      END DO
      DO I=1+NDI,NTENS
         DEDT(I)   = HALF*DSTRAN(I)
         STRANT(I) = HALF*STRAN(I)+DEDT(I)
      END DO
C-------  END  Total Strain Tensor STRANT = STRAN + DSTRAN
C---------------------------------------------------------------------*
C-------  Tensors Z and TN at time t
      UNMD   = ONE-PHI
      UNMDHK = UNMD*HK
      CALL TENSDEV(STRESS,STRESSDV,NDI,NSHR)
C-------
      DO I=1,NTENS
         Z(I) = STRESSDV(I)-HK*ALPHA(I)
      END DO
C-------
      CALL PDTSCA(NDI,NSHR,Z,Z,ZZ)
      IF(ZZ.EQ.ZERO) ZZ=ONE
      TDZZ = THALF/DSQRT(THALF*ZZ)
      DO I=1,NTENS
         TN(I )    = TDZZ*Z(I)
      END DO
C-------  FIN Tensors Z and TN at time t
C---------------------------------------------------------------------*
C-------  Norm WX
      CALL PDTSCA(NDI,NSHR,ALPHA,ALPHA,WX)
      WX = UNMDHK*DSQRT(WX)
C-------  END  Norm WX
C---------------------------------------------------------------------*
C-------  Tensor Ep1
      DO I=1,NTENS
         EPLAS1(I) = EPLAS(I)+UNMTETADTDL*TN(I)
      END DO
C-------  END  Tensor Ep1
C---------------------------------------------------------------------*
C-------  Tensor e*
      DO I=1,NTENS
         ESTAR(I) = STRANT(I)-EPLAS1(I)
      END DO
C-------  END  Tensor e*
C---------------------------------------------------------------------*
C-------  Tensor E* (Deviatoric tensor of e*)
      CALL TENSDEV(ESTAR,GESTAR,NDI,NSHR)
C-------  END  Tensor E*
C---------------------------------------------------------------------*
C-------  Tensor Sigma1D
      D2G = TWO*G
C-------
      DO I=1,NTENS
         STRES1D(I) = D2G*GESTAR(I)
      END DO
C-------  END  Tensor Sigma1D
C---------------------------------------------------------------------*
C-------  Tensor ALPHA1
      RDLSWX = (RDK*DLAN+RSK)*WX
C-------
      DO I=1,NTENS
         ALPHA1(I) = ALPHA(I)+UNMTETADT*(DLAN*TN(I)-RDLSWX*ALPHA(I))
      END DO
C-------  FIN Tensor ALPHA1
C---------------------------------------------------------------------*
C-------  Hydrostratic Stress
      POLD  = ONETHIRD*(STRESS(1)+STRESS(2)+STRESS(3))
      XKDEH = XK*(DSTRAN(1)+DSTRAN(2)+DSTRAN(3))
      XKDH  = XK*(STRANT(1)+STRANT(2)+STRANT(3))
      P     = POLD+UNMD*XKDEH
C-------  END  Tensor Sigma1D
C---------------------------------------------------------------------*
C-------  Calcul de r1
      R   = HI*PR
      WR  = R
      RWR = PR*WR
      PR1 = PR+UNMTETADTDL*(ONE-RDI*RWR)-UNMTETADT*RSI*RWR
C-------  END  Calcul de r1
C---------------------------------------------------------------------*
C-------  Damage Phi1
      IF(NCD.EQ.0) THEN
         PHI1 = ZERO
      ELSE
         PHI1 = PHI
      END IF
C
C-------  END  Damage Phi1
C---------------------------------------------------------------------*
C-------  Comment Check Code
      IF(M.GE.4)THEN
         WRITE(MPI,*)'KTIMEN TN'
         DO I=1,NTENS
            WRITE(MPI,101)TN(I)
         END DO
         WRITE(MPI,*)'KTIMEN STRAN'
         DO I=1,NTENS
            WRITE(MPI,101)STRAN(I)
         END DO
         WRITE(MPI,*)'KTIMEN STRANT'
         DO I=1,NTENS
            WRITE(MPI,101)STRANT(I)
         END DO
         WRITE(MPI,*)'KTIMEN ESTAR'
         DO I=1,NTENS
            WRITE(MPI,101)ESTAR(I)
         END DO
         WRITE(MPI,*)'KTIMEN GESTAR'
         DO I=1,NTENS
            WRITE(MPI,101)GESTAR(I)
         END DO
         WRITE(MPI,*)'KTIMEN EPLAS1'
         DO I=1,NTENS
            WRITE(MPI,101)EPLAS1(I)
         END DO
         WRITE(MPI,*)'KTIMEN ALPHA1'
         DO I=1,NTENS
            WRITE(MPI,101)ALPHA1(I)
         END DO
         WRITE(MPI,*)'KTIMEN STRES1D'
         DO I=1,NTENS
            WRITE(MPI,101)STRES1D(I)
         END DO
         WRITE(MPI,102)'KTIMEN PR1  =',PR1
         WRITE(MPI,102)'KTIMEN PHI1 =',PHI1
         WRITE(MPI,102)'KTIMEN WR   =',WR
         WRITE(MPI,102)'KTIMEN WX   =',WX
      END IF
C-------  FIN Comment Check Code
C-----------------------------------------------------------------------
C---------   FORMATS  --------------------------------------------------
C-----------------------------------------------------------------------
 101  FORMAT(1X,6(D20.14,2X))
 102  FORMAT(1X,A20,6(D20.14,2X))
C-----------------------------------------------------------------------
C=======================================================================
      RETURN
      END
C
      SUBROUTINE PDTMAT(N,A,B,AB)
C======================================================================
C----------------------------------------------------------------------
C-------  CALCUL DU PRODUIT CONTRACTE D'UN TENSEUR D'ORDRE 4  ---------
C--------------  AVEC UN TENSEUR D'ORDRE 2 SYMETRIQUE  ----------------
C----------------------------------------------------------------------
C======================================================================
C-------------
C  EN ENTREE :
C-------------
C  A(N,N)	: TENSEURS D'ORDRE QUATRE
C  B(N)		: TENSEURS D'ORDRE DEUX
C-------------
C  EN SORTIE :
C-------------
C  AB(N)	: PRODUIT CONTRACTE A:B
C-----------------------------------------------------------------------
C=======================================================================
C.1-----  Implicit, External
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C.2-----  Parameter
C.3-----  Dimension
      DIMENSION A(N,N)
      DIMENSION B(N),AB(N)
C.4-----  Real,Integer,Complex,Double precision,Logical,Character
C.5-----  Common
C.6-----  Equivalence
C.7-----  Data
      DATA ZERO/0.D0/
C.8-----  Definition de fonctions
C=======================================================================
      DO I=1,N
     	 AB(I)=ZERO
         DO J=1,N
            AB(I)=AB(I)+A(I,J)*B(J)
         END DO
      END DO
C=======================================================================
      RETURN
      END
C
      SUBROUTINE PDTSCA(NDI,NSHR,A,B,AB)
C======================================================================
C----------------------------------------------------------------------
C----------  CALCUL DU PRODUIT CONTRACTE DE DEUX TENSEURS   -----------
C--------------------   D'ORDRE 2 SYMETRIQUES  ------------------------
C----------------------------------------------------------------------
C======================================================================
C-------------
C  EN ENTREE :
C-------------
C  NDI		:  NOMBRE DE COMPOSANTES DES CONTRAINTES DIRECTES.
C  NSHR		:  NOMBRE DE COMPOSANTES DES CONTRAINTES TANGENTIELLES.
C  A,B		:  TENSEURS A ET B
C-------------
C  EN SORTIE :
C-------------
C  AB		:  PRODUIT CONTRACTE A:B
C
C=======================================================================
C.1-----  Implicit, External
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C.2-----  Parameter
C.3-----  Dimension
      DIMENSION A(NDI+NSHR),B(NDI+NSHR)
C.4-----  Real,Integer,Complex,Double precision,Logical,Character
C.5-----  Common
C.6-----  Equivalence
C.7-----  Data
      DATA ZERO,DEUX/0.D0,2.D0/
C.8-----  Definition de fonctions
C=======================================================================
      AB=ZERO
      DO I=1,NDI
         AB=AB+A(I)*B(I)
      END DO
      DO I=1,NSHR
         AB=AB+DEUX*A(NDI+I)*B(NDI+I)
      END DO
C=======================================================================
      RETURN
      END
C
      SUBROUTINE PDTTEN(A,B,ATB,NTENS,NDI)
C======================================================================
C----------------------------------------------------------------------
C----------  CALCUL DU PRODUIT TENSORIEL DE DEUX TENSEURS   -----------
C--------------------   D'ORDRE 2 SYMETRIQUES  ------------------------
C----------------------------------------------------------------------
C======================================================================
C-------------
C  EN ENTREE :
C-------------
C  NTENS	: DIMENSION DES TENSEURS
C  A et B	: TENSEURS D'ORDRE DEUX
C-------------
C  EN SORTIE :
C-------------
C  ATB		: ATB = A * B --->  TENSEUR D'ORDRE QUATRE
C=======================================================================
C.1-----  Implicit, External
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C.2-----  Parameter
C.3-----  Dimension
      DIMENSION A(NTENS),B(NTENS)
      DIMENSION ATB(NTENS,NTENS)
C.4-----  Real,Integer,Complex,Double precision,Logical,Character
C.5-----  Common
C.6-----  Equivalence
C.7-----  Data
      DATA DEUX/2.D0/
C.8-----  Definition de fonctions
C=======================================================================
      DO I=1,NTENS
         DO J=1,NTENS
            ATB(I,J)=A(I)*B(J)
         END DO
      END DO
C-------
      DO I=1,NTENS
         DO J=NDI+1,NTENS
            ATB(I,J)=DEUX*ATB(I,J)
         END DO
      END DO
C=======================================================================
      RETURN
      END
C
      SUBROUTINE PRDMAT(NTENS,VM1,VM2,VM3)
C======================================================================
C----------------------------------------------------------------------
C-----------------------  PRODUIT DE DEUX MATRICES  -------------------
C----------------------------------------------------------------------
C======================================================================
C-------------
C  EN ENTREE :
C-------------
C  VM1		: MATRICE OU TENSEUR D'ORDRE 4 (NDIM1,NDIM2)
C  VM2		: MATRICE OU TENSEUR D'ORDRE 4 (NDIM2,NDIM3)
C-------------
C  EN SORTIE :
C-------------
C  VM3		: MATRICE OU TENSEUR D'ORDRE 4 (NDIM1,NDIM3)
C-----------------------------------------------------------------------
C=======================================================================
C.1-----  Implicit, External
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C.2-----  Parameter
C.3-----  Dimension
      DIMENSION VM1(NTENS,NTENS),VM2(NTENS,NTENS),VM3(NTENS,NTENS)
C.4-----  Real,Integer,Complex,Double precision,Logical,Character
C.5-----  Common
C.6-----  Equivalence
C.7-----  Data
      DATA ZERO/0.D0/
C.8-----  Definition de fonctions
C=======================================================================
      DO I=1,NTENS
         DO J=1,NTENS
            R=ZERO
            DO K=1,NTENS
               R=R+VM1(I,K)*VM2(K,J)
            END DO
         VM3(I,J)=R
         END DO
      END DO
C=======================================================================
      RETURN
      END
C
      SUBROUTINE PREDELAS(DSTRAN,STRESS,ALPHA,PR,PHI,STRESSELAS,
     *           ALPHAELAS,PRELAS,HI,RSI,HK,RSK,G,XK,
     *           NPLAN,NTENS,NDI,NSHR,M,MPI)
C======================================================================
C----------------------------------------------------------------------
C-------------  CALCUL D'UNE SOLUTION PUREMENT ELASTIQUE  -------------
C----------------------------------------------------------------------
C======================================================================
C-------------
C  EN ENTREE :
C-------------
C  DSTRAN	: INCREMENT DE DEFORMATION
C  STRESS	: CONTRAINTES A L'INSTANT N
C  ALPHA	: TENSEUR DES CONTRAINTES CINEMATIQUES
C  PR		: VARIABLE D'ECROUISSAGE ISOTROPE
C  PHI		: DAMAGE
C  H		: CONSTANTE D'ECROUISSAGE ISOTROPE
C  RS		: CONSTANTE D'ECROUISSAGE ISOTROPE (RECOVERY)
C  PH		: CONSTANTE D'ECROUISSAGE CINEMATIQUE
C  PRS		: CONSTANTE D'ECROUISSAGE CINEMATIQUE (RECOVERY)
C  G		: MODULE DE CISAILLEMENT
C  XK		: COEFFICIENT DE LAME
C  XNU		: COEFFICIENT DE POISSON
C  NTENS	: LONGUEUR DES TABLEAU DSTRAN, STRESS, ...
C  NDI		: NOMBRE DE COMPOSANTES DIRECTES D'UN TENSEUR
C  NSHR		: NOMBRE DE COMPOSANTES TANGENTIELLES D'UN TENSEUR
C		  POINT D'INTEGRATION.
C  M		: COMMENTAIRES
C  MPI		: UNITE DE FICHIER DE SORTIE
C-------------
C  EN SORTIE :
C-------------
C  STRESSELAS	:  CONTRAINTES A L'INSTANT N+1 AVEC PREDICTION ELASTIQUE
C  ALPHAELAS	:  ECROUISSAGE CINEMATIQUE A L'INSTANT N+1 (TRIAL)
C  PRELAS	:  ECROUISSAGE ISOTROPE A L'INSTANT N+1 (TRIAL)
C  UNIT		:  TENSEUR UNITE 1*1
C  UNITDEV	:  TENSEUR UNITE DEVIATORIQUE Idev=I-1/3(1*1)
C-----------------------------------------------------------------------
C=======================================================================
C-----------------------------------------------------------------------
C.1-----  Implicit, External
      IMPLICIT NONE
C.2-----  Parameter,Integer
      INTEGER NPLAN,NTENS,NDI,NSHR,M,MPI,I,J
C.3-----  Dimension
      DOUBLE PRECISION DSTRAN(NTENS),STRESS(NTENS),ALPHA(NTENS)
      DOUBLE PRECISION STRESSELAS(NTENS),ALPHAELAS(NTENS)
      DOUBLE PRECISION UNITDEV(NTENS,NTENS)
      DOUBLE PRECISION DSTRAN2(6),DSTRANDV(6)
C-------
      DOUBLE PRECISION PHI,PR,PRELAS,HK,RSK,HI,RSI,G,XK
      DOUBLE PRECISION D2G,D2GD,UNMD,DSTRANH,DEHK,
     1       UNRSWX,UNMDH,XKGDH,DMWX,DMWXH
      DOUBLE PRECISION ZERO,ONE,TWO
C.4-----  Real,Complex,Double precision,Logical,Character
C.5-----  Common
C.6-----  Equivalence
C.7-----  Save,Data
      DATA ZERO,ONE,TWO/0.D0,1.D0,2.D0/
C.8-----  Definition de fonctions
C---------------------------------------------------------------------*
C======================================================================
C---------------------------------------------------------------------*
C-------  Division par deux des deformations tangentielles
C
      DO I=1,NDI
         DSTRAN2(I) = DSTRAN(I)
      END DO
      DO I=1+NDI,NTENS
         DSTRAN2(I) = DSTRAN(I)/TWO
      END DO
C
C-------  FIN Division par deux des deformations tangentielles
C---------------------------------------------------------------------*
C-------  Calcul des contraintes modifiees
C
C      IF(NPLAN.EQ.2) DSTRAN2(3)=-XNU/(ONE-XNU)*(DSTRAN2(1)+DSTRAN2(2))
C-------
      UNMD    = ONE-PHI
      D2G     = TWO*G
      D2GD    = UNMD*D2G
      DSTRANH = DSTRAN(1)+DSTRAN(2)+DSTRAN(3)
      DEHK    = UNMD*XK*DSTRANH
      CALL TENSDEV(DSTRAN2,DSTRANDV,NDI,NSHR)
      DO I=1,NTENS
         STRESSELAS(I) = STRESS(I)+D2GD*DSTRANDV(I)
      END DO
      DO I=1,NDI
         STRESSELAS(I) = STRESSELAS(I)+DEHK
      END DO
C
C-------  FIN Calcul des contraintes modifiees
C---------------------------------------------------------------------*
C-------  Calcul des contraintes cinematiques modifiees
C
      CALL SINV(ALPHA,DMWXH,DMWX,NDI,NSHR)
      DMWX = HK*DMWX
      UNRSWX = ONE-RSK*DMWX
      DO I=1,NTENS
         ALPHAELAS(I) = UNRSWX*ALPHA(I)
      END DO
C
C-------  FIN Calcul des contraintes cinematiques modifiees
C---------------------------------------------------------------------*
C-------  Calcul des contraintes isotropes modifiees
C
      PRELAS = (ONE-RSI*HI*PR)*PR
C
C-------  FIN Calcul des contraintes isotropes modifiees
C---------------------------------------------------------------------*
C-------  Comment Check Code
      IF(M.GE.4)THEN
         WRITE(MPI,*)'PREDELAS STRESSELAS'
         DO I=1,NTENS
            WRITE(MPI,101)STRESSELAS(I)
         END DO
         WRITE(MPI,*)'PREDELAS ALPHAELAS'
         DO I=1,NTENS
            WRITE(MPI,101)ALPHAELAS(I)
         END DO
         WRITE(MPI,102)'PREDELAS  PRELAS  = ',PRELAS
      END IF
C-------  FIN Comment Check Code
C----------------------------------------------------------------------*
C---------   FORMATS  -------------------------------------------------*
C----------------------------------------------------------------------*
 101  FORMAT(1X,6(D20.14,2X))
 102  FORMAT(1X,A20,6(D20.14,2X))
C----------------------------------------------------------------------*
C=======================================================================
      RETURN
      END
C
      SUBROUTINE RECOV(PROPS,NPROPS,STATEV,NSTATV,STRESS,STRAN,DSTRAN,
     *           EPLAS,ALPHA,DROT,DFGRD0,DFGRD1,XLAMBDA,PR,PHI,DLAN,
     *           TIME,DTIME,TEMP,DTEMP,E,XNU,G,XK,YIELD,V,PF,HI,RDI,RSI,
     *           HK,RDK,RSK,XIP,XMS,DC,TETA,EPSF,EPSX,EPSR,EPSD,
     *           VNEWDT,NCD,NPLAN,NATUR,IROT,NTENS,NDI,NSHR,NVI,M,MPI)
C======================================================================
C----------------------------------------------------------------------
C------  EXTRACTION DES VARIABLES UTILISATEURS ET DES CONSTANTES  -----
C----------------------------------------------------------------------
C======================================================================
C-----------
C  ENTREES :
C-----------
C  STAVEV(NSTATEV) : VARIABLES UTILISATEUR POUR LE MODELE UMAT
C  PROPS(NPROPS)   : CONSTANTES DEFINIES DANS *USER MATERIAL
C-----------
C  SORTIES :
C-----------
C  STRAN	: TENSEUR DES DEFORMATIONS TOTALES
C  DSTRAN	: TENSEUR INCREMENT DES DEFORMATIONS
C  EPLAS	: TENSEUR DES DEFORMATIONS PLASTIQUES
C  ALPHA	: TENSEUR DES CONTRAINTES CINEMATIQUES
C  DROT		: INCREMENT DE MATRICE ROTATION (WINGET)
C  DFGRD0	: GRADIENT DES DEFORMATIONS EN DEBUT D'INCREMENT
C  DFGRD1	: GRADIENT DES DEFORMATIONS EN FIN D'INCREMENT
C  XLAMBDA	: DEFORMATION PLASTIQUE CUMULEE
C  PR		: VARIABLE D'ECROUISSAGE ISOTROPE
C  ENDO		: TENSEUR D'ENDOMMAGEMENT
C  DLAN		: INCREMENT DE DEFORMATION PLASTIQUE CUMULEE PRECEDENT
C  TIME		: STEP TIME (#1) AND TOTAL TIME (#2) AT BEGINNING
C		  OF THE INCREMENT
C  DTIME	: TIME INCREMENT
C  TEMP		: TEMPERATURE AT BEGINNING OF INCREMENT
C  DTEMP	: TEMPERATURE INCREMENT
C  T0		: TEMPERATURE INITIALE
C  RHO		: DENSITE VOLUMIQUE
C  C		: COEFFICIENT ADIABATIQUE
C  E		: MODULE D'YOUNG
C  G		: MODULE DE CISAILLEMENT
C  XK		: COEFFICIENT DE LAME
C  XNU		: COEFFICIENT DE POISSON
C  YIELD	: LIMITE D'ELASTICITE
C  V		: CONSTANTE VISCOPLASTIQUE
C  PF		: CONSTANTE VISCOPLASTIQUE
C  RD		: CONSTANTE D'ECROUISSAGE ISOTROPE
C  H		: CONSTANTE D'ECROUISSAGE ISOTROPE
C  RS		: CONSTANTE D'ECROUISSAGE ISOTROPE (RECOVERY)
C  PRD		: CONSTANTE D'ECROUISSAGE CINEMATIQUE
C  PH		: CONSTANTE D'ECROUISSAGE CINEMATIQUE
C  PRS		: CONSTANTE D'ECROUISSAGE CINEMATIQUE (RECOVERY)
C  XIP		: INITIAL POROSITY
C  DC		: VALEUR DE L'ENDOMMAGEMENT CRITIQUE --> RUPTURE
C  XMS		: DAMAGE RATE SENSIVITY
C  TETA		: TETA-METHODE
C  EPSF		: PRECISION DE LA FONCTION DE CHARGE F
C  EPSX		: PRECISION DE LA FONCTION D'ENDOMMAGEMENT OMEGA-X
C  EPSR		: PRECISION DE LA FONCTION D'ENDOMMAGEMENT OMEGA-R
C  EPSD		: PRECISION DE LA FONCTION TENSORIELLE D'ENDOMMAGEMENT PHI
C  NCD		: NCD=1 --> COUPLAGE ELASTOPLASTICITE-ENDOMMAGEMENT
C  		  NCD=0 --> ELASTOPLASTICITE SANS ENDOMMAGEMENT
C  NPLAN	: NPLAN=1 --> DEFORMATIONS PLANES
C		  NPLAN=2 --> CONTRAINTES PLANES
C  NATUR	: NATUR=0 --> PETITES DEFORMATIONS
C		  NATUR=1 --> GRANDES DEFORMATIONS
C  IROT		: IROT=0 --> PAS DE ROTATION
C		  IROT=1 --> AVEC ROTATION
C		: IROT=2 --> AVEC DERIVEE COROTATIONNELLE DE JAUMANN
C		  IROT=3 --> AVEC ROTATION D'ABAQUS (ROTSIG)
C  NVI		: NOMBRE DE VARIABLE INTERNES STOCKEES (R,X,J2(X),Y,...)
C  M		: COMMENTAIRES (PARAMETRE D'IMPRESSION)
C  MPI		: UNITE DE FICHIER DE SORTIE
C----------------------------------------------------------------------*
C======================================================================*
C----------------------------------------------------------------------*
C.1-----  Precision
      IMPLICIT NONE
C.2-----  Parameter,Integer
      INTEGER NPT0,NPRHO,NPC,NPE,NPXNU,NPCTS,NPINITF,NPDC,NPMS,NPFN,
     1        NPTETA,NPEPSF,NPEPSX,NPEPSR,NPEPSD,NPNEWDT,NPNCD,
     2        NPPLAN,NPNATUR,NPIROT,NPNVI,NPMIMP,NCTS
      INTEGER NPROPS,NSTATV,NTENS,NDI,NSHR
      INTEGER NCD,NPLAN,NATUR,IROT,NVI,M,MPI,I,J
C-------
      PARAMETER(NPT0=1,NPRHO=2,NPC=3,NPE=4,NPXNU=5,NPCTS=5)
      PARAMETER(NPINITF=28,NPDC=29,NPMS=30)
      PARAMETER(NPTETA=31)
      PARAMETER(NPEPSF=32,NPEPSX=33,NPEPSR=34,NPEPSD=35,NPNEWDT=36)
      PARAMETER(NPNCD=37,NPPLAN=38,NPNATUR=39,NPIROT=40,NPNVI=41)
      PARAMETER(NPMIMP=42)
      PARAMETER(NCTS=22)
C.3-----  Dimension
      DOUBLE PRECISION STATEV(NSTATV),PROPS(NPROPS),STRESS(NTENS)
      DOUBLE PRECISION STRAN(NTENS),DSTRAN(NTENS)
      DOUBLE PRECISION EPLAS(NTENS),ALPHA(NTENS)
      DOUBLE PRECISION CTS(NCTS),TIME(2),STRESSI(6)
      DOUBLE PRECISION DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
      DOUBLE PRECISION TDROT(3,3),ROT0(3,3),ROT1(3,3)
      DOUBLE PRECISION TROT0(3,3),DR(3,3),DROTGN(3,3)
C-------
      DOUBLE PRECISION XLAMBDA,PHI,PR,DLAN,DTIME,TEMP,DTEMP
      DOUBLE PRECISION E,XNU,G,XK,YIELD,V,PF,HI,RDI,RSI,HK,RDK,RSK,
     1       XIP,DC,XMS,TETA,EPSF,EPSX,EPSR,EPSD,VNEWDT
      DOUBLE PRECISION T0,C,RHO
      DOUBLE PRECISION TUNMDXNU,UT0,YA,YB,DNORM
      DOUBLE PRECISION ZERO,UN,DEUX,TROIS,UDEMI,UTIER,DTIER,RTTD
C.4-----  Real,Complex,Double precision,Logical,Character
C.5-----  Common
C.6-----  Equivalence
C.7-----  Save,Data
      DATA ZERO,UN,DEUX,TROIS,UDEMI/0.D0,1.D0,2.D0,3.D0,0.5D0/
C.8-----  Definition de fonctions
C----------------------------------------------------------------------*
C======================================================================*
C----------------------------------------------------------------------*
C-------  Eclatement du tableau PROPS(NPROPS)
      M      = IDINT(PROPS(NPMIMP))
      T0     = PROPS(NPT0)
C      T0     = TEMP+DTEMP
      RHO    = PROPS(NPRHO)
      C      = PROPS(NPC)
      E      = PROPS(NPE)
      XNU    = PROPS(NPXNU)
      XIP    = PROPS(NPINITF)
      DC     = PROPS(NPDC)
      XMS    = PROPS(NPMS)
      TETA   = PROPS(NPTETA)
      EPSF   = PROPS(NPEPSF)
      EPSX   = PROPS(NPEPSX)
      EPSR   = PROPS(NPEPSR)
      EPSD   = PROPS(NPEPSD)
      VNEWDT = PROPS(NPNEWDT)
      NCD    = IDINT(PROPS(NPNCD))
      NPLAN  = IDINT(PROPS(NPPLAN))
      NATUR  = IDINT(PROPS(NPNATUR))
      IROT   = IDINT(PROPS(NPIROT))
      NVI    = IDINT(PROPS(NPNVI))
C----------- Constantes du modele BCJ (c1 --> c22)
      DO I=1,NCTS
         CTS(I)=PROPS(I+NPCTS)
      END DO
C-------  FIN Eclatement du tableau PROPS(NPROPS)
C----------------------------------------------------------------------*
C-------  Calcul du module de cisaillement G et de K
C
      G        = UDEMI*E/(UN+XNU)
      TUNMDXNU = TROIS*(UN-DEUX*XNU)
      XK       = E/TUNMDXNU
C
C-------  FIN  Calcul du module de cisaillement G et de K
C----------------------------------------------------------------------*
C-------  Calcul de la constante de nucleation
C
      UTIER = UN/TROIS
      DTIER = DEUX*UTIER
      RTTD  = DSQRT(TROIS*UDEMI)
C
C-------  FIN  Calcul de la constante de nucleation
C----------------------------------------------------------------------*
C-------  Calcul des constantes elastoplastiques
C
      UT0   = UN/T0
      YA    = UDEMI*CTS(3)/(CTS(21)+DEXP(-CTS(4)*UT0))
      YB    = UN+DTANH(CTS(19)*(CTS(20)-T0))
C-------
      V     = CTS(1)*DEXP(-CTS(2)*UT0)
      YIELD = YA*YB
C      YIELD = CTS(3)*DEXP(CTS(4)*UT0)
      PF    = DTIME*CTS(5)*DEXP(-CTS(6)*UT0)
      HK    = DTIER*(CTS(9)-CTS(10)*T0)
      RDK   = RTTD*CTS(7)*DEXP(-CTS(8)*UT0)
      RSK   = DTIME*RTTD*CTS(11)*DEXP(-CTS(12)*UT0)
      HI    = CTS(15)-CTS(16)*T0
      RDI   = CTS(13)*DEXP(-CTS(14)*UT0)
      RSI   = DTIME*CTS(17)*DEXP(-CTS(18)*UT0)
C-------
C      HK    = ZERO
C      RDK   = ZERO
C      RSK   = ZERO
C      HI    = ZERO
C      RDI   = ZERO
C      RSI   = ZERO
C      V     = ZERO
C
C-------  FIN  Calcul des constantes elastoplastiques
C----------------------------------------------------------------------*
C-------  Eclatement du tableau STATEV(NSTATV)
C----------  Deformations totales, plastiques et elastiques
C            Ecrouissage cinematique
      IF(NATUR.EQ.0) THEN
         DO I=1,NTENS
      	    EPLAS(I) = STATEV(I)
            ALPHA(I) = STATEV(I+NTENS)
         END DO
      ELSE 
         IF(IROT.EQ.1) THEN
            CALL ROTSIG(STATEV(1),DROT,EPLAS,1,NDI,NSHR)
            CALL ROTSIG(STATEV(1+NTENS),DROT,ALPHA,1,NDI,NSHR)
         ELSE
            DO I=1,NTENS
               STRESSI(I)=STRESS(I)
            END DO
            CALL TRANSPOSE(3,DROT,TDROT)
	    CALL DPOLR(ROT0,DFGRD0,M,MPI)
	    CALL DPOLR(ROT1,DFGRD1,M,MPI)
	    CALL TRANSPOSE(3,ROT0,TROT0)
	    CALL PRDMAT(3,ROT1,TROT0,DR)
	    DO I=1,3
	       DO J=1,3
	          DR(I,J)=UDEMI*DR(I,J)
	       END DO
	       DR(I,I)=UDEMI+DR(I,I)
	    END DO
            CALL PRDMAT(3,DR,TDROT,DROTGN)
            CALL ROTSIG(STRESSI,DROTGN,STRESS,1,NDI,NSHR)
            CALL ROTSIG(STATEV(1),DR,EPLAS,1,NDI,NSHR)
            CALL ROTSIG(STATEV(1+NTENS),DR,ALPHA,1,NDI,NSHR)
            CALL ROTSIG(STATEV(4+NVI+3*NTENS),DR,STRAN,2,NDI,NSHR)
            DO I=1,NTENS
               STATEV(I+3+NVI+3*NTENS)=STRAN(I)+DSTRAN(I)
            END DO
         END IF
      END IF
C---------  Deformation plastique cumulee
         XLAMBDA = STATEV(1+2*NTENS)
C---------  Ecrouissage isotrope
         PR = STATEV(2+2*NTENS)
C---------  Damage
      IF(NCD.EQ.0) THEN
         PHI = ZERO
         XIP = ZERO
      ELSE
         PHI = STATEV(3+2*NTENS)
C------------  Initial Damage as Initial Porosity
         IF(PHI.EQ.ZERO) THEN
            PHI = XIP
         ELSE
            XIP = ZERO
         END IF
      END IF
C---------  Increment de deformation plastique cumulee precedent
      DLAN=STATEV(4+2*NTENS)
C-------   FIN Eclatement du tableau STATEV(NSTATV)
C----------------------------------------------------------------------*
C----------------------------------------------------------------------*
C-------  Comment Check Code
      IF(M.GE.4)THEN
         WRITE(MPI,*)'RECOV EPLAS'
         DO I=1,NTENS
            WRITE(MPI,101)EPLAS(I)
         END DO
         WRITE(MPI,*)'RECOV ALPHA'
         DO I=1,NTENS
            WRITE(MPI,101)ALPHA(I)
         END DO
         WRITE(MPI,102)XLAMBDA,PHI,PR,DLAN
         WRITE(MPI,103) V,YIELD,PF,RDK,HK,RSK,RDI,HI,RSI
      END IF
C-------  END Comment Check Code
C----------------------------------------------------------------------*
C---------   FORMATS  -------------------------------------------------*
C----------------------------------------------------------------------*
 101  FORMAT(1X,3(D20.14,2X))
 401  FORMAT(1X,8(D20.14,2X))
 102  FORMAT(1X,'XLAMBDA=',D20.14,2X,'PHI=',D20.14,2X,
     *     'PR=',D20.14,2X,'DLAN=',D20.14,2X)
 103  FORMAT(
     * '*-----------------------------------------------------*'//,
     * '*-----   Material Constants   -----*'/,
     *       '  V(T)      = ',E12.5/,'  Y(T)      = ',E12.5/,
     *       '  f(T)      = ',E12.5/,'  rd(T)     = ',E12.5/,
     *       '  h(T)      = ',E12.5/,'  rs(T)     = ',E12.5/,
     *       '  Rd(T)     = ',E12.5/,'  H(T)      = ',E12.5/,
     *       '  Rs(T)     = ',E12.5/)
C----------------------------------------------------------------------*
C----------------------------------------------------------------------*
C=======================================================================
      RETURN
      END
C
      SUBROUTINE SOLELAS(STRESS,STATEV,DDSDDE,STRESSELAS,ALPHAELAS,
     *           PRELAS,PHI,XK,G,XIP,NCD,NSTATV,NTENS,NDI,NVI,M,MPI)
C======================================================================
C----------------------------------------------------------------------
C------------------------  SOLUTION ELASTIQUE  ------------------------
C----------------------------------------------------------------------
C======================================================================
C                                
C ENTREES
C -------
C  ENDO		: TENSEUR D'ENDOMMAGEMENT
C  STRESSELAS	: TENSEUR DES CONTRAINTES AVEC PREDICTION ELASTIQUE
C  PHI		: DAMAGE
C  XK		: COEFFICIENT DE LAME
C  G		: MODULE DE CISAILLEMENT
C  XIP		: INITIAL POROSITY
C  NCD		: NCD =1 --> COUPLED WITH DAMAGE
C  NSTATV	: TAILLE DU TABLEAU STATEV
C  NTENS	: TAILLE DU TABLEAU DES CONTRAINTES OU DEFORMATIONS
C  NDI		: NOMBRE DE COMPOSANTES DIRECTES D'UN TENSEUR
C  NVI		: NOMBRE DE VARIABLE INTERNES STOCKEES (EE,R,X,J2(X),Y,...)
C  M		: COMMENTAIRES (PARAMETRE D'IMPRESSION)
C  MPI		: UNITE DE FICHIER DE SORTIE
C--------	
C LOCALES
C -------
C  DMT		: TENSEUR D'EFFET D'ENDOMMAGEMENT (M)
C  DVDM	 	: PRODUIT Idev * Mt
C--------	
C SORTIES
C -------
C  STRESS	: TENSEUR DES CONTRAINTES
C  DDSDDE 	: MATRICE JACOBIENNE ELASTIQUE
C----------------------------------------------------------------------*
C=======================================================================
C----------------------------------------------------------------------*
C.1-----  Implicit, External
      IMPLICIT NONE
C.2-----  Parameter,Integer
      INTEGER NCD,NSTATV,NTENS,NDI,NVI,M,MPI,I,J
C.3-----  Dimension
      DOUBLE PRECISION DDSDDE(NTENS,NTENS)
      DOUBLE PRECISION STRESS(NTENS),STATEV(NSTATV)
      DOUBLE PRECISION STRESSELAS(NTENS),ALPHAELAS(NTENS)
C-------
      DOUBLE PRECISION PHI,PRELAS,XK,G,XIP,D2G,XKG
      DOUBLE PRECISION UNMD,UNMDG,UNMDD2G,UNMDXKG
      DOUBLE PRECISION UN,DEUX,TROIS,UDEMI
C.4-----  Real,Complex,Double precision,Logical,Character
C.5-----  Common
C.6-----  Equivalence
C.7-----  Save,Data
      DATA UN,DEUX,TROIS,UDEMI/1.D0,2.D0,3.D0,0.5D0/
C.8-----  Definition de fonctions
C----------------------------------------------------------------------*
C=======================================================================
C----------------------------------------------------------------------*
C----------------------------------------------------------------------*
C-------  Store of Cauchy Stress
C
      CALL AFFECT(NTENS,STRESSELAS,STRESS)
C
C-------  END  Store of Cauchy Stress
C----------------------------------------------------------------------*
C-------  Store of Kinematic Hardening
C
      DO I=1,NTENS
         STATEV(I+NTENS) = ALPHAELAS(I)
      END DO
C
C-------  END  Store of Kinematic Hardening
C----------------------------------------------------------------------*
C-------  Store of Isotropic Hardening
C
      STATEV(2+2*NTENS) = PRELAS
C
C-------  END  Store of Isotropic Hardening
C----------------------------------------------------------------------*
C-------  Stockage des variables internes
C
         IF(NCD.EQ.1 .OR. NCD.EQ.2) THEN
            STATEV(3+2*NTENS) = PHI
         END IF
C
C-------  END  Stockage des variables internes
C----------------------------------------------------------------------*
C----------------------------------------------------------------------*
C-------  Elastic Jacobian
C
      UNMD    = UN-PHI
      D2G     = DEUX*G
      XKG     = XK-D2G/TROIS
      UNMDG   = UNMD*G
      UNMDD2G = UNMD*D2G
      UNMDXKG = UNMD*XKG
C-------
      DO I=1,NDI
         DO J=1,NDI
            DDSDDE(I,J) = UNMDXKG
         END DO
         DDSDDE(I,I) = DDSDDE(I,I)+UNMDD2G
      END DO
      DO I=NDI+1,NTENS
         DDSDDE(I,I) = UNMDG
      END DO
         WRITE(MPI,*)'SOLELAS DDSDDE'
         DO I=1,NTENS                                                 
            WRITE(MPI,101)(DDSDDE(I,J),J=1,NTENS)
         END DO
C
C-------  END  Elastic Jacobian
C----------------------------------------------------------------------*
C-------  Comment Check Code
      IF(M.GE.4)THEN
         WRITE(MPI,*)'SOLELAS STRESS'
         DO I=1,NTENS                                                 
            WRITE(MPI,101)STRESS(I)
         END DO
         WRITE(MPI,*)'SOLELAS DDSDDE'
         DO I=1,NTENS                                                 
            WRITE(MPI,101)(DDSDDE(I,J),J=1,NTENS)
         END DO
      END IF
C-------  END Comment Check Code
C----------------------------------------------------------------------*
C---------   FORMATS  -------------------------------------------------*
C----------------------------------------------------------------------*
 101  FORMAT(1X,6(D20.14,2X))
C----------------------------------------------------------------------*
C=======================================================================
      RETURN
      END
C
      SUBROUTINE SOLSYST(DLAMBDA,PHI,WX,WR,FF,FX,FR,FD,
     *           DFDL,DFDX,DFDR,DFDD,DXDL,DXDX,DXDR,DXDD,DRDL,DRDR,
     *           DHDL,DHDX,DHDR,DHDD,NCD,NTENS,NDI,NSHR,M,MPI)
C======================================================================
C----------------------------------------------------------------------
C-----------  CALCUL DE LA SOLUTION CLAMBDA, CWX, CWR et CD  ----------
C-----------------------------  DU SYSTEME  ---------------------------
C----------  PUIS REACTUALISATION DE DLAMBDA, WX, WR et ENDO  ---------
C----------------------------------------------------------------------
C======================================================================
C-------------
C  EN ENTREE :
C-------------
C  FF		: VALEUR DE LA FONCTION F
C  FX		: VALEUR DE LA FONCTION GX
C  FR		: VALEUR DE LA FONCTION GR
C  FH		: VALEUR DE LA FONCTION H
C  DFDL		: DERIVEE PARTIELLE DE F PAR RAPPORT A DLAMBDA
C  DFDX		: DERIVEE PARTIELLE DE F PAR RAPPORT A WX
C  DFDR		: DERIVEE PARTIELLE DE F PAR RAPPORT A WR
C  DFDD		: DERIVEE PARTIELLE DE F PAR RAPPORT A ENDO
C  DXDL		: DERIVEE PARTIELLE DE GX PAR RAPPORT A DLAMBDA
C  DXDX		: DERIVEE PARTIELLE DE GX PAR RAPPORT A WX
C  DXDD		: DERIVEE PARTIELLE DE GX PAR RAPPORT A ENDO
C  DRDL		: DERIVEE PARTIELLE DE GR PAR RAPPORT A DLAMBDA
C  DRDR		: DERIVEE PARTIELLE DE GR PAR RAPPORT A WR
C  DHDL		: DERIVEE PARTIELLE DE H PAR RAPPORT A DLAMBDA
C  DHDX		: DERIVEE PARTIELLE DE H PAR RAPPORT A WX
C  DHDR		: DERIVEE PARTIELLE DE H PAR RAPPORT A WR
C  DHDD		: DERIVEE PARTIELLE DE H PAR RAPPORT A ENDO
C  NCD		: NCD=1 --> COUPLAGE ELASTOPLASTICITE-ENDOMMAGEMENT
C  		  NCD=0 --> ELASTOPLASTICITE SANS ENDOMMAGEMENT
C  NTENS	: TAILLE DU TABLEAU DES CONTRAINTES OU DEFORMATIONS
C  NDI		: NOMBRE DE COMPOSANTES DES CONTRAINTES DIRECTES
C  NSHR		: NOMBRE DE COMPOSANTES DES CONTRAINTES TANGENTIELLES
C  M		: COMMENTAIRES
C  MPI		: UNITE DE FICHIER DE SORTIE
C-------------
C  EN SORTIE :
C-------------
C  DLAMBDA	: INCREMENT DEFORMATION PLASTIQUE CUMULEE
C  WR		: ISOTROPIC DYNAMIC RECOVERY TERM: WR = DSINH[QS*R]
C  WX		: KINEMATIC DYNAMIC RECOVERY TERM: WX = |DMI:X|
C  ENDO		: TENSEUR D'ENDOMMAGEMENT
C-----------------------------------------------------------------------
C=======================================================================
C-----------------------------------------------------------------------
C.1-----  Implicit, External
      IMPLICIT NONE
C.2-----  Parameter,Integer
      INTEGER NCD,NTENS,NDI,NSHR,M,MPI,I,J
C.3-----  Dimension
      DOUBLE PRECISION DLAMBDA,PHI,WX,WR,FF,FX,FR,FD
      DOUBLE PRECISION DFDL,DFDX,DFDR,DFDD,DXDL,DXDX,DXDR,DXDD,
     1       DRDL,DRDR,DHDL,DHDX,DHDR,DHDD
      DOUBLE PRECISION UDXDX,UDRDR,DFXX,DHXX,DFRR,DXRR,DHRR,AB,ABL
      DOUBLE PRECISION AA,AL,AD,BB,BL,BD
      DOUBLE PRECISION DL,DX,DR,DD
      DOUBLE PRECISION ZERO,ONE
C.4-----  Real,Complex,Double precision,Logical,Character
C.5-----  Common
C.6-----  Equivalence
C.7-----  Save,Data
      DATA ZERO,ONE/0.D0,1.D0/
C.8-----  Definition de fonctions
C---------------------------------------------------------------------*
C======================================================================
C---------------------------------------------------------------------*
C-------  Initialisation
C
C
C-------  END Initialisation
C---------------------------------------------------------------------*
C-------  Solution DD=cd (Damage) and DL=clambda (Plastic equivalent strain)
      UDXDX = ONE/DXDX
      UDRDR = ONE/DRDR
      DFXX  = DFDX*UDXDX
      DHXX  = DHDX*UDXDX
      DFRR  = DFDR*UDRDR
      DXRR  = DXDR*UDRDR
      DHRR  = DHDR*UDRDR
      AB    = FX-DXRR*FR
      ABL   = DXDL-DXRR*DRDL
C-------
      AA    = FF-DFXX*AB-DFRR*FR
      AL    = DFDL-DFXX*ABL-DFRR*DRDL
      AD    = DFDD-DFXX*DXDD
      BB    = FD-DHXX*AB-DHRR*FR
      BL    = DHDL-DHXX*ABL-DHRR*DRDL
      BD    = DHDD-DHXX*DXDD
C------------------------------------------*
C-------  Solution DD=cd (Damage)
      IF(NCD.EQ.0) THEN
         DD = ZERO
      ELSE
         DD = (AA*BL-BB*AL)/(AL*BD-AD*BL)
      END IF
C------------------------------------------*
C-------  Solution Dl=clambda (Plastic Equivalent Strain Increment)
      DL = -(AA+AD*DD)/AL
C------------------------------------------*
C-------  Solutions DX=cWd (norm of Alpha) and DR=cWr (Kappa**2)
      DX = -UDXDX*(AB+ABL*DL+DXDD*DD)
      DR = -UDRDR*(FR+DRDL*DL)
C------------------------------------------*
C---------------------------------------------------------------------*
C-------  Update of the solution clambda, cWx, cWr et cd
C
      DLAMBDA = DLAMBDA+DL
      WX      = WX+DX
      WR      = WR+DR
      PHI     = PHI+DD
C
C----------- END  Update of the solution clambda, cWx, cWr et cd
C---------------------------------------------------------------------*
C---------------------------------------------------------------------*
C-------  Comment Check Code
C
      IF(M.GE.4) THEN
         WRITE(MPI,102)'SOLSYST DL = ',DL
         WRITE(MPI,102)'SOLSYST DX = ',DX
         WRITE(MPI,102)'SOLSYST DR = ',DR
         WRITE(MPI,102)'SOLSYST DD = ',DD
         WRITE(MPI,102)'SOLSYST DLAMBDA = ',DLAMBDA
         WRITE(MPI,102)'SOLSYST WX      = ',WX
         WRITE(MPI,102)'SOLSYST WR      = ',WR
         WRITE(MPI,102)'SOLSYST PHI     = ',PHI
      END IF
C
C-------  END Comment Check Code
C----------------------------------------------------------------------*
C---------   FORMATS  -------------------------------------------------*
C----------------------------------------------------------------------*
 101  FORMAT(1X,6(D20.14,2X))
 102  FORMAT(1X,A20,6(D20.14,2X))
C----------------------------------------------------------------------*
C=======================================================================
      RETURN
      END
C
      SUBROUTINE SOLUTION(STRESS,EPLAS,ALPHA,XLAMBDA,PR,STRANT,PHI,
     *           ZD,TN,STRES1D,EPLAS1,ALPHA1,PR1,DLAMBDA,DLAN,TDL,UNMD,
     *           UNTRWX,UNTRWR,RTDL,EYV,G,XK,HK,HI,XMS,TETA,DTIME,
     *           STVIT,NITMAX,NSTVI,NTENS,NDI,NSHR,NCD,M,MPI)
C======================================================================
C----------------------------------------------------------------------
C--------------  CALCUL DE LA EPLAS(N+1), ALPHA(N+1),  ----------------
C---------------------  PR(N+1) et STRESS(N+1)  -----------------------
C----------------------------------------------------------------------
C======================================================================
C-------------
C  EN ENTREE :
C-------------
C  ZD		: NORMAL DIRECTION STRESS TENSOR
C  TN		: NORMAL TENSOR
C  STRANT	: TOTAL STRAIN TENSOR
C  STRES1D	: KNOWN DEVIATORIC STRESS TENSOR
C  EPLAS1	: KNOWN PLASTIC STRAIN TENSOR
C  ALPHA1	: KNOWN KINEMATIC HARDENING TENSOR
C  PR1		: KNOWN ISOTROPIC HARDENING VARIABLE
C  DLAMBDA	: PLASTIC EQUIVALENT STRAIN RATE
C  DLAN		: PLASTIC EQUIVALENT STRAIN RATE AT PRECENDENT STEP
C  DM		: DAMAGE EFFECT TENSOR Mt
C  DMI		: DAMAGE EFFECT TENSOR M(-t)
C  DN		: ORIENTATION VOID OPERATOR
C  QD		: VOID COALESCENCE OPERATOR
C  TDL		: TETA*DLAMBDA
C  UDHM		: VALUE 1/(1-Dh)**m
C  UNTRWX	: 1+TETA*(PRD*DLAMBDA+PRS)*WX
C  UNTRWR	: 1+TETA*(RD*DLAMBDA+RS)*WR
C  RTDL		: R1+TETA*DLAMBDA
C  Y		: ENERGY RELEASE RATE Y = Ye + Yan
C  YT		: TRIAXIALITY YT = <Sh>/Seq
C  YPS		: FUNCTION YPS = (Y/S)**s
C  EYV		: FUNCTION Exp[3/2*<Sh>/Seq]
C  G		: SHEAR MODULUS
C  XK		: LAME COEFFICIENT
C  H		: ISOTROPIC HARDENING MODULUS
C  PH		: KINEMATIC HARDENING MODULUS
C  XIP		: INITIAL POROSITY
C  XDN		: NUCLEATION CONSTANTS XLD*XFN/(XKIC*XLF**(1/3))
C  CTE		: CONSTANT OF TENSION
C  CCO		: CONSTANT OF COMPRESSION
C  CTO		: CONSTANT OF TORSION
C  AV		: VOID GROWTH CONSTANT
C  VP		: VOID GROWTH TEST (VC = 1 if Dh > Vc)
C  VT		: POSITIVE PART <Dh-Vc>  (Dh=1/3*Trace(Phi))
C  AC		: COALESCENCE CONSTANT
C  TETA		: TETA-METHODE
C  DTIME	: TIME INCREMENT
C  NITMAX	: MAXIMUM ITERATION NUMBER
C  NSTVI	: INTERNAL VARIABLE ARRAY STVIT SIZE
C  NDI		: DIAGONAL COMPONENT NUMBER
C  NSHR		: SHEAR COMPONENT NUMBER
C  NTENS	: STRESS TENSOR ARRAY SIZE (NDI + NSHR).
C  M		: COMMENTS
C  MPI		: OUTPUT UNIT FILE
C-------------
C  EN SORTIE :
C-------------
C  STRESS	: CAUCHY STRESS TENSOR AT TIME N+1
C  EPLAS	: PLASTIC STRAIN TENSOR AT TIME N+1
C  EELAS	: ELASTIC STRAIN TENSOR AT TIME N+1
C  ALPHA	: KINEMATIC HARDENING TENSOR AT TIME N+1
C  XLAMBDA	: PLASTIC EQUIVALENT STRAIN A N+1
C  DLAMBDA	: PLASTIC EQUIVALENT STRAIN INCREMENT
C  PR		: ISOTROPIC HARDENING TENSOR AT TIME N+1
C  STVIT	: INTERNAL VARIABLE ARRAY: R, J2(X) et X
C  DAMAGE	: RUPTURE INDICATOR
C-----------------------------------------------------------------------
C=======================================================================
C-----------------------------------------------------------------------
C.1-----  Implicit, External
      IMPLICIT NONE
C.2-----  Parameter,Integer
      INTEGER NITMAX,NSTVI,NTENS,NDI,NSHR,NCD,M,MPI,I
C.3-----  Dimension
      DOUBLE PRECISION STRESS(NTENS),EPLAS(NTENS),ALPHA(NTENS)
      DOUBLE PRECISION STRANT(NTENS),EPLAS1(NTENS),ALPHA1(NTENS)
      DOUBLE PRECISION STVIT(NSTVI)
      DOUBLE PRECISION ZD(NTENS),TN(NTENS),STRES1D(NTENS)
      DOUBLE PRECISION X(6),DLTN(6)
C-------
      DOUBLE PRECISION XLAMBDA,PR,PR1,DLAMBDA,PHI,DLAN,TDL,UNMD,
     1       UNTRWX,UNTRWR,RTDL,EYV,G,XK,HK,HI,XMS,TETA,DTIME
      DOUBLE PRECISION D2GTDL,EH,DEHK,XH,XEQV,EFFSTRAN
      DOUBLE PRECISION ZERO,ONE,TWO,THREE,TWOTHIRD
C.4-----  Real,Complex,Double precision,Logical,Character
C.5-----  Common
C.6-----  Equivalence
C.7-----  Save,Data
      DATA ZERO,ONE,TWO,THREE/0.D0,1.D0,2.D0,3.D0/
C.8-----  Functions Definitions
C---------------------------------------------------------------------*
C======================================================================
C---------------------------------------------------------------------*
C-------  Initialization
C
      TWOTHIRD = TWO/THREE
C
C-------  END Initialization
C---------------------------------------------------------------------*
C-------  Plastic Strain EPLAS
C
      DO I=1,NTENS
         DLTN(I)  = TDL*TN(I)
         EPLAS(I) = EPLAS1(I)+DLTN(I)
      END DO
C
C-------  END  Plastic Strain EPLAS
C---------------------------------------------------------------------*
C-------  Kinematic Hardening Tensors ALPHA and X
C
      DO I=1,NTENS
         ALPHA(I) = (ALPHA1(I)+DLTN(I))/UNTRWX
         X(I)     = HK*ALPHA(I)
      END DO
      IF(HK.EQ.ZERO) CALL ASET(ALPHA,ZERO,NTENS)
C
C-------  END  Kinematic Hardening Tensors ALPHA and X
C---------------------------------------------------------------------*
C-------  Plastic Equivalent Strain Rate DLAMBDA
C
C      DLAMBDA = (ONE-TETA)*DLAN*DTIME+TETA*DLAMBDA
C
C-------  END Plastic Equivalent Strain Rate DLAMBDA DLAMBDA
C---------------------------------------------------------------------*
C-------  Plastic Equivalent Strain XLAMBDA
C
      XLAMBDA = XLAMBDA+(ONE-TETA)*DLAN*DTIME+TETA*DLAMBDA
C
C-------  END Plastic Equivalent Strain XLAMBDA
C---------------------------------------------------------------------*
C-------  Isotropic Hardening PR
C
      PR = RTDL/UNTRWR
      IF(HI.EQ.ZERO) PR = ZERO
C
C-------  END Isotropic Hardening PR
C---------------------------------------------------------------------*
C-------  Cauchy Stress Tensor
C
      CALL TRACE(STRANT,EH,NDI,NTENS)
C-------
      D2GTDL = TWO*G*UNMD*TDL
      DO I=1,NTENS
         STRESS(I) = UNMD*STRES1D(I)-D2GTDL*TN(I)
      END DO
C-------
      DEHK = XK*UNMD*EH
      DO I=1,NDI
         STRESS(I) = STRESS(I)+DEHK
      END DO
C
C-------  END Cauchy Stress Tensor
C---------------------------------------------------------------------*
C---------------------------------------------------------------------*
C-------  Storage of R, J2(X), X, Y, Ye and Damage Tensors in STVIT array.
C
C---------  Isotropic and Kinematic Hardening Stresses
      CALL SINV(X,XH,XEQV,NDI,NSHR)
      STVIT(1)       = HI*PR
      STVIT(2)       = XEQV
      DO I=1,NTENS
         STVIT(I+2)  = X(I)
      END DO
C---------  Triaxiality <p>/Seqv
      STVIT(3+NTENS) = EYV
C---------  Reached Iteration Number
      STVIT(4+NTENS) = DBLE(NITMAX)
C---------  Effective Strain
      CALL PDTSCA(NDI,NSHR,STRANT,STRANT,EFFSTRAN)
      STVIT(5+NTENS) = DSQRT(TWOTHIRD*EFFSTRAN)
C
C-------  END Storage of R, J2(X), X, Y, Ye and Damage Tensors in STVIT array.
C---------------------------------------------------------------------*
C---------------------------------------------------------------------*
C-------  Comment Check Code
C
      IF(M.GE.4) THEN
         WRITE(MPI,*)'SOLUTION STRESS'
         DO I=1,NTENS
            WRITE(MPI,101)STRESS(I)
         END DO
         WRITE(MPI,*)'SOLUTION EPLAS'
         DO I=1,NTENS
            WRITE(MPI,101)EPLAS(I)
         END DO
         WRITE(MPI,*)'SOLUTION ALPHA'
         DO I=1,NTENS
            WRITE(MPI,101)ALPHA(I)
         END DO
         WRITE(MPI,102)'SOLUTION XLAMBDA =',XLAMBDA
         WRITE(MPI,102)'SOLUTION PR =',PR
      END IF
C
C-------  END Comment Check Code
C----------------------------------------------------------------------*
C---------   FORMATS  -------------------------------------------------*
C----------------------------------------------------------------------*
 101  FORMAT(1X,6(D20.14,2X))
 102  FORMAT(1X,A20,6(D20.14,2X))
C----------------------------------------------------------------------*
C=======================================================================
      RETURN
      END
C
      SUBROUTINE STORE(STATEV,NSTATV,EPLAS,ALPHA,XLAMBDA,PR,PHI,
     *           DLAMBDA,DTIME,STVIT,NSTVI,NTENS,NDI,NVI,NCD,M,MPI)
C======================================================================
C----------------------------------------------------------------------
C----  ENREGISTREMENT DES VARIABLES UTILISATEURS A L'INSTANT N+1  -----
C----------------------------------------------------------------------
C======================================================================
C-----------
C  ENTREES :
C-----------
C  EPLAS	: TENSEUR DES DEFORMATIONS PLASTIQUES
C  ALPHA	: TENSEUR DES CONTRAINTES CINEMATIQUES
C  XLAMBDA	: DEFORMATION PLASTIQUE CUMULEE
C  PR		: VARIABLE D'ECROUISSAGE ISOTROPE
C  PHI		: ISOTROPIC DAMAGE
C  DLAMBDA	: INCREMENT DE DEFORMATION PLASTIQUE CUMULEE
C  STVIT(NSTVI)	: TABLEAU CONTENANT LES VARIABLE (R,X,J2(X),Y,...)
C  NSTATV	: TAILLE DU TABLEAU STATEV
C  NVI		: NOMBRE DE VARIABLE INTERNES STOCKEES (R,X,J2(X),Y,...)
C  M		: COMMENTAIRES (PARAMETRE D'IMPRESSION)
C  MPI		: UNITE DE FICHIER DE SORTIE
C-----------
C  SORTIES :
C-----------
C  STAVEV(NSTATEV) : VARIABLES UTILISATEUR POUR LE MODELE UMAT
C-----------------------------------------------------------------------
C=======================================================================
C-----------------------------------------------------------------------
C.1-----  Precision
      IMPLICIT NONE
C.2-----  Parameter,Integer
      INTEGER NSTATV,NSTVI,NTENS,NDI,NVI,NCD,M,MPI,I
C.3-----  Dimension
      DOUBLE PRECISION STATEV(NSTATV),STVIT(NSTVI)
      DOUBLE PRECISION EPLAS(NTENS),ALPHA(NTENS)
C-------
      DOUBLE PRECISION XLAMBDA,PR,PHI,DLAMBDA,DTIME
C.4-----  Real,Complex,Double precision,Logical,Character
C.5-----  Common
C.6-----  Equivalence
C.7-----  Save,Data
C.8-----  Definition de fonctions
C-----------------------------------------------------------------------
C=======================================================================
C-----------------------------------------------------------------------
C-------  Rangement du tableau STATEV(NSTATV)
C
C----------  Deformations plastiques  et elastiques 
C            Ecrouissage cinematique (ALPHA)
         DO I=1,NTENS
      	    STATEV(I)       = EPLAS(I)
            STATEV(I+NTENS) = ALPHA(I)
         END DO
C---------  Deformation plastique cumulee
         STATEV(1+2*NTENS) = XLAMBDA
C---------  Ecrouissage isotrope (r)
         STATEV(2+2*NTENS) = PR
C---------  Endommagement
         IF(NCD.EQ.1 .OR. NCD.EQ.2) THEN
            STATEV(3+2*NTENS) = PHI
         END IF
C---------  Increment de deformation plastique cumulee precedent
         STATEV(4+2*NTENS) = DLAMBDA/DTIME
C------------------
C-------  Test: si NVI > 0 ---> Stockage des variables internes suivantes:
C		 - Isotropic Hardening (R)
C		 - Kinematic Hardening Norm J2(X) if NVI > 1
C		 - Back Stress (X) if NVI > 2
C		 - Triaxiality Ye if NVI > 2 + NTENS
C		 - Effective Strain Eff
C		 - Number of Iterations Nitmax
         DO I=1,NVI+1
            STATEV(I+4+2*NTENS) = STVIT(I)
         END DO
C
C-------   FIN Enregistrement du tableau STATEV(NSTATV)
C-----------------------------------------------------------------------
C-------  Comment Check Code
      IF(M.GE.4)THEN
         WRITE(MPI,*)'STORE EPLAS'
         DO I=1,NTENS
            WRITE(MPI,101)STATEV(I)
         END DO
         WRITE(MPI,*)'STORE ALPHA'
         DO I=1,NTENS
            WRITE(MPI,101)STATEV(I+NTENS)
         END DO
         WRITE(MPI,102)'STORE PHI     = ',STATEV(3+2*NTENS)
         WRITE(MPI,102)'STORE XLAMBDA = ',STATEV(1+2*NTENS)
         WRITE(MPI,102)'STORE PR      = ',STATEV(2+2*NTENS)
         WRITE(MPI,102)'STORE DLAN    = ',STATEV(4+2*NTENS)
      END IF
C-------  FIN Comment Check Code
C-----------------------------------------------------------------------
C---------   FORMATS  --------------------------------------------------
C-----------------------------------------------------------------------
 101  FORMAT(1X,3(D20.14,2X))
 102  FORMAT(1X,A20,6(D20.14,2X))
C-----------------------------------------------------------------------
C=======================================================================
      RETURN
      END
C
      SUBROUTINE TENSDEV(TENS,TENSDV,NDI,NSHR)
C======================================================================
C----------------------------------------------------------------------
C--------------  CALCUL DU TENSEUR DEVIATORIQUE DE TENS  --------------
C----------------------------------------------------------------------
C======================================================================
C-----------
C  ENTREES :
C-----------
C  TENS		: TENSEUR D'ORDRE DEUX
C  SINV1	: INVARIANT DE TENS SINV1=TRACE(TENS)/3
C  NDI		:
C  NSHR		:
C-----------
C  SORTIES :
C-----------
C  TENSDV	: TENSEUR DEVIATORIQUE DE TENS
C----------------------------------------------------------------------
C======================================================================
C.1-----  Include,Precision
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C.2-----  Parameter
C.3-----  Dimension
      DIMENSION TENS(NDI+NSHR),TENSDV(NDI+NSHR)
C.4-----  Real,Integer,Complex,Double precision,Logical,Character
C.5-----  Common
C.6-----  Equivalence
C.7-----  Save,Data
      DATA ZERO,TROIS/0.D0,3.D0/
C.8-----  Definition de fonctions
C-----------------------------------------------------------------------
C=======================================================================
C-----------------------------------------------------------------------
      TENSH=ZERO
      DO I=1,NDI
         TENSH=TENSH+TENS(I)
      END DO
      UTTENSH=TENSH/TROIS
      DO I=1,NDI
         TENSDV(I)=TENS(I)-UTTENSH
      END DO
      DO I=1,NSHR
         TENSDV(I+NDI)=TENS(I+NDI)
      END DO
C----------------------------------------------------------------------
C======================================================================
      RETURN
      END
C
      SUBROUTINE TRACE(TENS,TENSH,NDI,NTENS)
C======================================================================
C----------------------------------------------------------------------
C--------------------------  TRACE OF A TENSOR  -----------------------
C----------------------------------------------------------------------
C======================================================================
C-------------
C  EN ENTREE :
C-------------
C  TENS		: TENSEUR
C  NDI		: NOMBRE DE COMPOSANTES DIAGONALES DU TENSEUR
C  NDI		: NOMBRE TOTAL DE COMPOSANTES DU TENSEUR
C-------------
C  EN SORTIE :
C-------------
C  TENSH	: TRACE DU TENSEUR
C-----------------------------------------------------------------------
C=======================================================================
C.1-----  Implicit, External
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C.2-----  Parameter
C.3-----  Dimension
      DIMENSION TENS(NTENS)
C.4-----  Real,Integer,Complex,Double precision,Logical,Character
C.5-----  Common
C.6-----  Equivalence
C.7-----  Data
      DATA ZERO/0.D0/
C.8-----  Definition de fonctions
C=======================================================================
      TENSH=ZERO
      DO I=1,NDI
         TENSH=TENSH+TENS(I)
      END DO
C=======================================================================
      RETURN
      END
C
      SUBROUTINE TRANSPOSE(N,A,B)
C======================================================================
C----------------------------------------------------------------------
C----------------------  TRANSPOSITION OF A MATRICE  ------------------
C--------------------------   B = TRANSPOSED(A)  ----------------------
C----------------------------------------------------------------------
C======================================================================
C-------------
C  EN ENTREE :
C-------------
C  A		: TENSEUR
C  N		: DIMENSION OF TENSOR
C-------------
C  EN SORTIE :
C-------------
C  B		: TRANSPOSED TENSOR
C-----------------------------------------------------------------------
C=======================================================================
C.1-----  Implicit, External
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C.2-----  Parameter
C.3-----  Dimension
      DIMENSION A(N,N),B(N,N)
C.4-----  Real,Integer,Complex,Double precision,Logical,Character
C.5-----  Common
C.6-----  Equivalence
C.7-----  Data
C.8-----  Definition de fonctions
C=======================================================================
      DO I=1,N
         DO J=1,N
            B(J,I)=A(I,J)
         END DO
      END DO
C=======================================================================
      RETURN
      END
C
      SUBROUTINE VALPROP(A,D,V,N,NP)
C=======================================================================
C-----------------------------------------------------------------------
C-------- CALCUL DES VALEURS PROPRES ET LES VECTEURS PROPRES   ---------
C-------  D'UNE MATRICE REELLE A CARREE SYMETRIQUE D'ORDRE M  ----------
C---------------------  EN RESOLVANT A.X=LAMDA.X  ----------------------
C-----------------------------------------------------------------------
C-----------------------------  METHODE  -------------------------------
C--------------  REDUCTION DE A SOUS FORME TRI-DIAGONALE  --------------
C------------  RESOLUTION DU PROBLEME AUX VALEURS PROPRES  -------------
C-----------------------  PAR L'ALGORITHME QR  -------------------------
C-----------------------------------------------------------------------
C----------------  BIBLIOTHEQUE: NUMERICAL RECIPES  --------------------
C-----------------------------------------------------------------------
C=======================================================================
C-----------
C  ENTREES :
C-----------
C  A 		: MATRICE REELLE SYMETRIQUE
C  N		: DIMENSION DE LA MATRICE A A DIAGONALISER
C  NP		: DIMENSION DE LA MATRICE A
C-----------
C  SORTIES :
C-----------
C  A 		: MATRICE UNITE
C  D 		: VECTEUR CONTENANT LES VALEURS PROPRES DE LA MATRICE A
C  V 		: TABLEAU DE MEME DIMENSION QUE LA MATRICE A ET
C 		  CONTENANT LES VECTEURS PROPRES (RANGES PAR COLONNE)
C----------------------------------------------------------------------*
C=======================================================================
C----------------------------------------------------------------------*
C.1-----  Implicit, External
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C.2-----  Parameter
      PARAMETER (NMAX=100)
C.3-----  Dimension
      DIMENSION A(NP,NP),V(NP,NP),D(NP),E(3)
C.4-----  Real,Integer,Complex,Double precision,Logical,Character
C.5-----  Common
C.6-----  Equivalence
C.7-----  Save,Data
      DATA ZERO,UN,DEUX/0.D0,1.D0,2.D0/
C.8-----  Definition de fonctions
C----------------------------------------------------------------------*
C=======================================================================
C----------------------------------------------------------------------*
C-------  Tridiagonalisation de la matrice A
      CALL TRED(A,N,NP,D,E)
C-------  FIN  Tridiagonalisation de la matrice A
C----------------------------------------------------------------------*
C-------  Initialisation de V
      DO I=1,N
        DO J=1,N
           V(I,J)=ZERO
        END DO
        V(I,I)=UN
      END DO
C-------  FIN  Initialisation de V
C----------------------------------------------------------------------*
C-------  Calcul des valeurs et vecteurs propres
C
      IF(N.GT.1) THEN
         DO I=2,N
            E(I-1)=E(I)
         END DO
         E(N)=ZERO
C----------
         DO L=1,N
            ITER=0
C-------------
100         DO M=L,N-1
               DD=DABS(D(M))+DABS(D(M+1))
               IF (DABS(E(M))+DD.EQ.DD) GOTO 200
            END DO
C-------------
            M=N
C-------------
200         IF(M.NE.L)THEN
               IF(ITER.EQ.30)PAUSE 'too many iterations'
               ITER=ITER+1
               G=(D(L+1)-D(L))/(DEUX*E(L))
               R=SQRT(G*G+UN)
               G=D(M)-D(L)+E(L)/(G+SIGN(R,G))
               S=UN
               C=UN
               P=ZERO
C----------------
               DO I=M-1,L,-1
                  F=S*E(I)
                  B=C*E(I)
C-------------------
                  IF(DABS(F).GE.DABS(G))THEN
                     C=G/F
                     R=DSQRT(C*C+UN)
                     E(I+1)=F*R
                     S=UN/R
                    C=C*S
                  ELSE
                     S=F/G
                     R=DSQRT(S*S+UN)
                     E(I+1)=G*R
                     C=UN/R  
                     S=S*C
                  ENDIF
C-------------------
                  G=D(I+1)-P
                  R=(D(I)-G)*S+DEUX*C*B
                  P=S*R
                  D(I+1)=G+P
                  G=C*R-B
C-------------------
                  DO K=1,N
                     F=V(K,I+1)
                     V(K,I+1)=S*V(K,I)+C*F
                     V(K,I)=C*V(K,I)-S*F
                  END DO
C-------------------
               END DO
C----------------
               D(L)=D(L)-P
               E(L)=G
               E(M)=ZERO
               GOTO 100
            ENDIF
C-------------
         END DO
C----------
      ENDIF
C
C-------  Calcul des valeurs et vecteurs propres
C----------------------------------------------------------------------*
C=======================================================================
C----------------------------------------------------------------------*
      RETURN
      END
C
      SUBROUTINE TRED(A,N,NP,D,E)
C=======================================================================
C-----------------------------------------------------------------------
C---------------  REDUCTION D'UNE REELLE SYMETRIQUE  -------------------
C--------------------  SOUS FORME TRI-DIAGONALE  -----------------------
C-----------------------------------------------------------------------
C----------------  BIBLIOTHEQUE: NUMERICAL RECIPES  --------------------
C-----------------------------------------------------------------------
C=======================================================================
C-----------
C  ENTREES :
C-----------
C  A 		: MATRICE REELLE SYMETRIQUE
C  N		: DIMENSION DE LA MATRICE A A TRIDIAGONALISER
C  NP		: DIMENSION DE LA MATRICE A
C-----------
C  SORTIES :
C-----------
C  A 		: MATRICE UNITE
C  D 		: COMPOSANTES DIAGONALES DE LA MATRICE TRIDIAGONALES
C  E 		: COMPOSANTES NON-DIAGONALES AVEC E(1)=ZERO
C----------------------------------------------------------------------*
C=======================================================================
C----------------------------------------------------------------------*
C.1-----  Implicit, External
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C.2-----  Parameter
C.3-----  Dimension
      DIMENSION A(NP,NP),D(NP),E(NP)
C.4-----  Real,Integer,Complex,Double precision,Logical,Character
C.5-----  Common
C.6-----  Equivalence
C.7-----  Save,Data
      DATA ZERO,UN/0.D0,1.D0/
C.8-----  Definition de fonctions
C----------------------------------------------------------------------*
C=======================================================================
C----------------------------------------------------------------------*
      IF(N.GT.1)THEN
         DO I=N,2,-1  
            L=I-1
            H=ZERO
            SCALE=ZERO
C-------------
            IF(L.GT.1)THEN
               DO K=1,L
                  SCALE=SCALE+DABS(A(I,K))
               END DO
C----------------
               IF(SCALE.EQ.ZERO)THEN
                  E(I)=A(I,L)
               ELSE
                  DO K=1,L
                     A(I,K)=A(I,K)/SCALE
                     H=H+A(I,K)*A(I,K)
                  END DO
C-------------------
                  F=A(I,L)
                  G=-SIGN(DSQRT(H),F)
                  E(I)=SCALE*G
                  H=H-F*G
                  A(I,L)=F-G
                  F=ZERO
C-------------------
                  DO J=1,L
                     A(J,I)=A(I,J)/H
                     G=ZERO
                     DO K=1,J
                        G=G+A(J,K)*A(I,K)
                     END DO
C----------------------
                     IF(L.GT.J)THEN
                        DO K=J+1,L
                           G=G+A(K,J)*A(I,K)
                        END DO
                     ENDIF
C----------------------
                     E(J)=G/H
                     F=F+E(J)*A(I,J)
                  END DO
C-------------------
                  HH=F/(H+H)
C-------------------
                  DO J=1,L
                     F=A(I,J)
                     G=E(J)-HH*F
                     E(J)=G
                     DO K=1,J
                        A(J,K)=A(J,K)-F*E(K)-G*A(I,K)
                     END DO
                  END DO
C-------------------
               ENDIF
C----------------
            ELSE
               E(I)=A(I,L)
            ENDIF
C-------------
            D(I)=H
         END DO
      ENDIF
C-------
C-------
      D(1)=ZERO
      E(1)=ZERO
C-------
      DO I=1,N
         L=I-1
         IF(D(I).NE.ZERO)THEN
            DO J=1,L
               G=ZERO
               DO K=1,L
                  G=G+A(I,K)*A(K,J)
               END DO
C----------------
               DO K=1,L
                  A(K,J)=A(K,J)-G*A(K,I)
               END DO
            END DO
         ENDIF
C----------
         D(I)=A(I,I)
         A(I,I)=UN
C----------
         IF(L.GE.1)THEN
            DO J=1,L
               A(I,J)=ZERO
               A(J,I)=ZERO
            END DO
         ENDIF
      END DO
C----------------------------------------------------------------------*
C=======================================================================
C----------------------------------------------------------------------*
      RETURN
      END
C
      SUBROUTINE VERIFGENCE(CONVERGENCE,FF,FX,FR,FD,
     *           EPSF,EPSX,EPSR,EPSD,NTENS,NDI,NSHR,M,MPI)
C======================================================================
C----------------------------------------------------------------------
C--------------  VERIFICATION DES CRITERES DE CONVERGENCE  ------------
C----------------------------------------------------------------------
C======================================================================
C-----------
C  ENTREES :
C-----------
C  FF		: VALEUR DE LA FONCTION F (CRITERE)
C  FX		: VALEUR DE LA FONCTION GX (Wx = |Xt|)
C  FR		: VALEUR DE LA FONCTION (Wr = Sinh(Qs*R))
C  FD		: VALEUR DE LA FONCTION TENSORIELLE D'ENDOMMAGEMENT
C  EPSF		: PRECISION DE LA FONCTION F (CRITERE)
C  EPSGX	: PRECISION DE LA FONCTION GX (Wx = |Xt|)
C  EPSGR	: PRECISION DE LA FONCTION (Wr = Sinh(Qs*R))
C  EPSH		: PRECISION DE LA FONCTION TENSORIELLE D'ENDOMMAGEMENT
C  NTENS	: TAILLE DU TABLEAU DES CONTRAINTES OU DEFORMATIONS
C  NDI		: NOMBRE DE COMPOSANTES DES CONTRAINTES DIRECTES
C  NSHR		: NOMBRE DE COMPOSANTES DES CONTRAINTES TANGENTIELLES
C  M		: COMMENTAIRES (PARAMETRE D'IMPRESSION)
C  MPI		: UNITE DE FICHIER DE SORTIE
C-----------
C  SORTIES :
C-----------
C  CONVERGENCE	: VALEUR LOGIQUE INDIQUANT SI LES CRITERES SONT VERIFIES
C-----------------------------------------------------------------------
C=======================================================================
C-----------------------------------------------------------------------
C.1-----  Precision
      IMPLICIT NONE
C.2-----  Parameter
C.3-----  Dimension
C-------
      DOUBLE PRECISION FF,FX,FR,FD,EPSF,EPSX,EPSR,EPSD
      DOUBLE PRECISION DABSFF,DABSFX,DABSFR,DABSFD
C.4-----  Real,Integer,Complex,Double precision,Logical,Character
      INTEGER NTENS,NDI,NSHR,M,MPI
      LOGICAL CONVERGENCE,FCRIT,XCRIT,RCRIT,DCRIT
C.5-----  Common
C.6-----  Equivalence
C.7-----  Save,Data
C.8-----  Definition de fonctions
C-----------------------------------------------------------------------
C=======================================================================
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
C-------  Initialisation
C
      FCRIT  = .FALSE.
      XCRIT = .FALSE.
      RCRIT = .FALSE.
      DCRIT  = .FALSE.
C
C-------  FIN Initialisation
C-----------------------------------------------------------------------
C
C-------  Verification du critere |FF| < EPSF
C
      DABSFF = DABS(FF)
C-------
      IF(DABSFF.LT.EPSF) THEN
         FCRIT = .TRUE.
         IF(M.GE.4) WRITE(MPI,201)
      END IF
C
C----------- FIN  Verification du critere |FF| < EPSF
C-----------------------------------------------------------------------
C
C-------  Verification du critere |FX| < EPSX
C
      DABSFX = DABS(FX)
C-------
      IF(DABSFX.LT.EPSX) THEN
         XCRIT = .TRUE.
         IF(M.GE.4) WRITE(MPI,202)
      END IF        
C
C----------- FIN  Verification du critere |FX| < EPSX
C-----------------------------------------------------------------------
C-------  Verification du critere |FR| < EPSR
C
      DABSFR = DABS(FR)
C-------
      IF(DABSFR.LT.EPSR) THEN
         RCRIT = .TRUE.
         IF(M.GE.4) WRITE(MPI,203)
      END IF        
C
C----------- FIN  Verification du critere |FR| < EPSR
C-----------------------------------------------------------------------
C-------  Verification du critere |FD| < EPSD
C
      DABSFD = DABS(FD)
C-------
      IF(DABSFD.LT.EPSD) THEN
         DCRIT = .TRUE.
         IF(M.GE.4) WRITE(MPI,204)
      END IF        
C
C----------- FIN  Verification du critere |FD| < EPSD
C-----------------------------------------------------------------------
C
C-------  Verification de la solution
C
      IF(FCRIT .AND. XCRIT .AND. RCRIT .AND. DCRIT) THEN
         CONVERGENCE = .TRUE.
         IF(M.GE.4) WRITE(MPI,205)
      END IF
C
C----------- FIN  Verification de la solution
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-------  Comment Check Code
      IF(M.GE.0) THEN
         WRITE(MPI,102)'VERIFGENCE |FF| = ',DABS(FF)
         WRITE(MPI,102)'VERIFGENCE |FX| = ',DABS(FX)
         WRITE(MPI,102)'VERIFGENCE |FR| = ',DABS(FR)
         WRITE(MPI,102)'VERIFGENCE |FH| = ',DABS(FD)
       END IF
C-------  FIN Comment Check Code
C-----------------------------------------------------------------------
C---------   FORMATS  --------------------------------------------------
C-----------------------------------------------------------------------
 101  FORMAT(1X,3(D20.14,2X))
 102  FORMAT(1X,A20,6(D20.14,2X))
 201  FORMAT('VERIFGENCE CONVERGENCE DE FF')
 202  FORMAT('VERIFGENCE CONVERGENCE DE FGX')
 203  FORMAT('VERIFGENCE CONVERGENCE DE FGR')
 204  FORMAT('VERIFGENCE CONVERGENCE DE FH')
 205  FORMAT('VERIFGENCE CONVERGENCE DE LA SOLUTION')
C-----------------------------------------------------------------------
C=======================================================================
      RETURN
      END
C
      SUBROUTINE DISP(U,KSTEP,KINC,TIME,NODE,JDOF)
C-----------------------------------------------------------------*
C-----------------------------------------------------------------*
C-----  CALCULATION OF IMPOSED DISPLACEMENT AND VELOCITY  --------*
C--------------  TO OBTAIN CONSTANT STRAIN RATE  -----------------*
C-----------------------------------------------------------------*
C-----------------------------------------------------------------*
C  -----------
C  EN ENTREE :
C  -----------
C     KSTEP   : NUMERO DU STEP
C     KINC    : NUMERO DE L'INCREMENT
C     TIME(1) : VALEUR ACTUELLE DU TEMPS DE STEP
C     TIME(2) : VALEUR ACTUELLE DU TEMPS TOTAL
C     NODE    : NUMERO DU NOEUD
C     JDOF    : DEGRE DE LIBERTE
C  -----------
C  EN SORTIE :
C  -----------
C     U(1)    : TOTAL VALUE OF THE PRESCRIBED VALUE AT THIS POINT.
C               THE VARIABLE MAY BE DISPLACEMENT, ROTATION, PORE
C               PRESSURE, TEMPERATURE, ETC., DEPENDING ON THE DEGREE
C               OF FREEDOM CONSTRAINED.
C     U(2)    : dU(1)/dTIME
C     U(3)    : d2U(1)/(dTIME)2
C  ===============================================================
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION U(3),TIME(2)
      DATA ZERO,UN,DEUX/0.D0,1.D0,2.D0/
C-----------------------------------------------------------------
C-----------------------------------------------------------------
C---------INTERPOLATION DU DEPLACEMENT IMPOSES
C-----------------------------------------------------------------
C
      XL0  = UN
      RATE = +1.D+00
      U(1) = XL0*DEXP(RATE*TIME(2))-XL0
      U(2) = RATE*XL0*DEXP(RATE*TIME(2))
      U(3) = RATE*U(2)
C
C======================================================================
      RETURN
      END
C
