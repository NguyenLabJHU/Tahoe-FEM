// $Id: APS_MatlT.h,v 1.4 2004-02-04 00:40:52 raregue Exp $
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
						kgamma0_dot_1,
						kgamma0_dot_2,
						kgamma0_dot_3,
						km_rate,
						km1_x,
						km1_y,	
						km2_x,
						km2_y,	
						km3_x,
						km3_y,	 					
						kl,
		 				kH,
		 				kkappa0_1,
		 				kkappa0_2,
		 				kkappa0_3,
		        		kNUM_APS_MATL_PARAMS };

		void Allocate (void) { n_mp = kNUM_APS_MATL_PARAMS; Parameter.Dimension ( n_mp ); }

};

}
#endif
