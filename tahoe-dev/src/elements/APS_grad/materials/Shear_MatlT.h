// $Id: Shear_MatlT.h,v 1.2 2003-09-19 00:47:07 raregue Exp $
#ifndef _SHEAR_MATLT_ 
#define _SHEAR_MATLT_ 

#include "APS_MaterialT.h"

namespace Tahoe {

class Shear_MatlT : public APS_MaterialT 
{
	public:

		 Shear_MatlT ( void ) { Allocate(); }
		~Shear_MatlT ( void ) { } 

		enum ParamT { 
		 				kMu, 			// Shear Modulus
		 				kDensity,		// Mass density
		        		kNUM_SHEAR_MATL_PARAMS };

		void Allocate (void) { n_mp = kNUM_SHEAR_MATL_PARAMS; Parameter.Dimension ( n_mp ); }

};

}
#endif
