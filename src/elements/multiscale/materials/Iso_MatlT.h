// $Id: Iso_MatlT.h,v 1.2 2003-02-03 04:40:28 paklein Exp $
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
						kLamda,
		 				kMu, 		// Shear Modulus
						kBulk, 	// Bulk Modulus
		        kNUM_ISO_MATL_PARAMS };

		void Allocate (void) { n_mp = kNUM_ISO_MATL_PARAMS; Parameter.Dimension ( n_mp ); }

};

}
#endif
