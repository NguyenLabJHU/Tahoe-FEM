/* $Id: ElementT.h,v 1.8 2002-07-02 19:55:13 cjkimme Exp $ */

#ifndef _ELEMENT_T_H_
#define _ELEMENT_T_H_

/* forward declarations */
#include "ios_fwd_decl.h"

/** class to define element type enumeration. */

namespace Tahoe {

class ElementT
{
public:

	/** element types */
	enum TypeT {
                    kRod = 1, /**< pair potential */
                kElastic = 2, /**< small strain solid */
           kHyperElastic = 3, /**< large strain solid */
             kLocalizing = 4,
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
        kAugLagContact2D = 16,
     kTotLagHyperElastic = 17,
        kMeshFreeElastic = 18,
      kMeshFreeFDElastic = 19,
    kD2MeshFreeFDElastic = 20,
        kLinearDiffusion = 21,
     kMFCohesiveSurface  = 22,
           kACME_Contact = 23,
    kMultiplierContact3D = 24,
   kTotLagrExternalField = 26, /**< experimental/temporary for loosely coupled problems */
   kNonsingularContinuum = 27, /**< nonsingular continuum element */ 
    kMultiplierContactElement2D = 28,
       kSimoFiniteStrain = 29,  /**< enhanced strain element */
kPenaltyContactElement2D = 30,
             kMultiScale = 31,
	    kCoarseScale = 32,
	      kFinePhest = 33,
kPenaltyContactElement3D = 34
	};

	/** stream extraction operator */ 
	friend istream& operator>>(istream& in, ElementT::TypeT& type);
};

} // namespace Tahoe 
#endif /* _ELEMENT_T_H_ */
