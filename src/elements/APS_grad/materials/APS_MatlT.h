// $Id: APS_MatlT.h,v 1.1 2003-07-10 17:26:55 raregue Exp $
#ifndef _APS_MATLT_
#define _APS_MATLT_

#include "APS_MaterialT.h"

namespace Tahoe {

class APS_MatlT : public APS_MaterialT 
{
	public:

		 APS_MatlT 	( void ) { Allocate(); }
		~APS_MatlT 	( void ) { } 

		enum ParamT { 		
	 					kMu, 			// 	Shear Modulus
						kgamma0_dot,
						km_rate,
						km1,
						km2,	 					
						kl,
		 				kH,
		        		kNUM_APS_MATL_PARAMS };

		void Allocate (void) { n_mp = kNUM_APS_MATL_PARAMS; Parameter.Dimension ( n_mp ); }

};

}
#endif
