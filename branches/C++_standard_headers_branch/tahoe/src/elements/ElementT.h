/* $Id: ElementT.h,v 1.54 2007-05-17 17:01:56 cfoster01 Exp $ */
#ifndef _ELEMENT_T_H_
#define _ELEMENT_T_H_

/* forward declarations */
#include "ios_fwd_decl.h"

namespace Tahoe {

/** class to define element type enumeration. */
class ElementT
{
public:

	/** element types */
	enum TypeT {
                    kRod = 1, /**< pair potential */
                kElastic = 2, /**< small strain solid */
           kHyperElastic = 3, /**< updated Lagragian large strain solid */
             kLocalizing = 4, /**< experimental */
                kVariTri = 5,
              kSWDiamond = 6, /**< diamond-cubic with Stillinger-Weber potentials */
         kMixedSWDiamond = 7,
         kUnConnectedRod = 8,
             kVirtualRod = 9, /**< pair potential with periodic boundary conditions */
            kVirtualSWDC = 10,
        kCohesiveSurface = 11,
         kThermalSurface = 12,
            kViscousDrag = 13,
         kPenaltyContact = 14,
             kBEMelement = 15,
          kAugLagContact = 16,
     kTotLagHyperElastic = 17, /**< total Lagragian large strain solid */
        kMeshFreeElastic = 18,
      kMeshFreeFDElastic = 19,
    kD2MeshFreeFDElastic = 20,
        kLinearDiffusion = 21,
     kMFCohesiveSurface  = 22,
           kACME_Contact = 23,
    kMultiplierContact3D = 24,
   kTotLagrExternalField = 26, /**< experimental/temporary for loosely coupled problems */
   kNonsingularContinuum = 27, /**< obsolete */
kMultiplierContactElement2D = 28,
       kSimoFiniteStrain = 29,  /**< enhanced strain element */
kPenaltyContactElement2D = 30,
    kStaggeredMultiScale = 31,
            kCoarseScale = 32,
              kFineScale = 33,
kPenaltyContactElement3D = 34,
          kBridgingScale = 35,
               kSimoQ1P0 = 36, /**< Q1P0, finite strain, mixed element */
               kAdhesion = 37, /**< adhesive tractions between surfaces */
           kParticlePair = 38,  /**< particles with pair interactions */
                    kEAM = 39,  /**< particles with EAM potental */
      kParticleThreeBody = 40,
     kNonLinearDiffusion = 41,
       kMeshfreeBridging = 45,
			 kFSMatForce = 60,    /*UpLag with material force calculation*/
			kSSMatForceD = 61,
			kSSMatForceS = 62,
		kSmallStrainQ2P1 = 64, /*small strain with mat force calculation*/		     
			   kSSQ2P1MF = 65,
		kSmallStrainQ1P0 = 66,
			   kSSQ1P0MF = 67,
				kAPSgrad = 68, /* anti-plane shear gradient plasticity */
			  kSS_SCNIMF = 70, /**< small strain stabilized, conforming nodally-integrated Galerkin MF */
			  kFS_SCNIMF = 71, /**< finite deformation ditto */
        kGradSmallStrain = 72,
			   kAPSVgrad = 80, /* anti-plane shear gradient plasticity, vector element */
		  kMeshfreeGradP = 85, /* meshfree gradient plasticity */
		kSS_EnhStrainLoc = 90, /* small strain enhanced strain embedded discontinuity element */
   kSS_EnhStrainLocCraig = 91, /* small strain enhanced strain embedded discontinuity element */
kSS_EnhStrainLocDieterich = 92, /* small strain enhanced strain embedded discontinuity element */
   kSS_EnhStrainLocOpen = 93, /* small strain enhanced strain embedded discontinuity element */
 kSS_EnhStrainLocLinear = 94, /* small strain embedded discontinuity element w/ linear jump*/
            kElasticAxi = 102, /**< small strain axisymmetric solid */
       kHyperElasticAxi = 103, /**<  updated Lagrangian, large strain axisymmetric solid */
 kTotLagHyperElasticAxi = 104, /**<  total Lagrangian, large strain axisymmetric solid */
           kSimoQ1P0Axi = 105, /**<  Q1P0 (mixed), large strain axisymmetric solid */
  kMeshFreeFDElasticAxi = 106, /**<  total Lagragian, large strain axisymmetric meshfree solid */     
   kHyperElasticInitCSE = 111, /**< large strain solid that triggers CSE */
	kPenaltyContactDrag = 114, /**< contact with constant drag traction */
kMeshfreePenaltyContact = 115, /**< contact with meshfree strikers */
kTotLagSplitIntegration = 117,
           kSimoQ1P0Inv = 136, /**< Q1P0, finite strain, mixed element with inverse dilation */
        kSimoQ1P0InvAxi = 137, /**< axisymmetric Q1P0, finite strain, mixed element with inverse dilation */
            kTotLagFlat = 217,  /**< simplified total Lagragian solid */
            	kTersoff = 138,	/**< Particles with Tersoff Potential */
        kFSSolidFluidMix = 139	/**< 3D dynamic Total Lagrangian solid fluid mixture element */
	};

	/** convert integer to ElementT::TypeT */
	static TypeT int2TypeT(int i);
};

} /* namespace Tahoe */

#endif /* _ELEMENT_T_H_ */
