// $Id: Iso_MatlT.h,v 1.3 2003-03-17 22:05:33 creigh Exp $
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
		        kNUM_ISO_MATL_PARAMS };

		void Allocate (void) { n_mp = kNUM_ISO_MATL_PARAMS; Parameter.Dimension ( n_mp ); }

};

}
#endif
