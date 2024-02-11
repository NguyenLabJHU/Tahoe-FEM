// $Id: Iso_MatlT.h,v 1.5 2003-04-23 23:34:26 creigh Exp $
#ifndef _ISO_MATLT_ 
#define _ISO_MATLT_ 

#include "VMF_MaterialT.h"

namespace Tahoe {

class Iso_MatlT : public VMF_MaterialT 
{
	public:

		 Iso_MatlT ( void ) { Allocate(); }
		~Iso_MatlT ( void ) { } 

		enum ParamT { 
						kE, 			// YoungsModulus
						kPr, 			// PoissonsRatio
						kLamda, 	// Lamda (Lame constant)
		 				kMu, 			// Shear Modulus (Lame constant)
						kBulk, 		// Bulk Modulus
						kDensity, 

						kE1, 			// 2nd YoungsModulus for Cba
						kPr1, 		// 2nd PoissonsRatio for Cba
						kLamda1,	// 2nd Lamda (Lame constant) for Cba
		 				kMu1, 		// 2nd Shear Modulus (Lame constant) for Cba
						kE2, 			// 2nd YoungsModulus for Cba
						kPr2, 		// 2nd PoissonsRatio for Cba
						kLamda2,	// 2nd Lamda (Lame constant) for Cba
		 				kMu2, 		// 2nd Shear Modulus (Lame constant) for Cba
						kPi, 			// Used for Cba 
						kRho, 		// Used for Cba
		 				kGamma_b, // Used for Cba
		 				kAlphaY,  // Used for Cba

		        kNUM_ISO_MATL_PARAMS };

		void Allocate (void) { n_mp = kNUM_ISO_MATL_PARAMS; Parameter.Dimension ( n_mp ); }

};

}
#endif
