// $Id: Iso_MatlT.h,v 1.4 2003-03-28 21:36:39 creigh Exp $
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
						kE, 		// YoungsModulus
						kPr, 		// PoissonsRatio
						kLamda, // Lamda (Lame constant)
		 				kMu, 		// Shear Modulus (Lame constant)
						kBulk, 	// Bulk Modulus
						kPi, 		// Used for Eb control 
						kRho, 	// Used for Eb control 
						kDensity, 
		        kNUM_ISO_MATL_PARAMS };

		void Allocate (void) { n_mp = kNUM_ISO_MATL_PARAMS; Parameter.Dimension ( n_mp ); }

};

}
#endif
