/* $Id: ElementT.h,v 1.19.4.1 2003-05-22 20:32:51 paklein Exp $ */
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
           kParticlePair = 38, /**< particles with pair interactions */
                    kEAM = 39, /**< particles with EAM potental */
       kMeshfreeBridging = 45,
	     kFSMatForce = 60,     /**< UpLag with material force calculation*/
	     kSSMatForce = 61};    /**< small strain with mat force calculation*/
 /** stream extraction operator */ 
	friend istream& operator>>(istream& in, ElementT::TypeT& type);
};

} // namespace Tahoe 
#endif /* _ELEMENT_T_H_ */
