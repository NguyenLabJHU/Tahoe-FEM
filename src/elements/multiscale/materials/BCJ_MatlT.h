// $Id: BCJ_MatlT.h,v 1.4 2003-03-17 22:05:33 creigh Exp $
#ifndef _BCJ_MATLT_
#define _BCJ_MATLT_

#include "VMF_MaterialT.h"

namespace Tahoe {

class BCJ_MatlT : public VMF_MaterialT 
{
	public:

		 BCJ_MatlT 	( void ) { Allocate(); }
		~BCJ_MatlT 	( void ) { } 

		enum ParamT { 
						kE, 		// 	Young's Modulus
						kPr,		//	Poisson's	Ratio
						kLamda,
		 				kMu, 		// 	Shear Modulus
						kBulk,	// 	Bulk Modulus
						kl,
					 	kc_zeta,
		 				kh, 
		 				kf,
		 				kY,
		 				kV,
		 				kCo,
		 				kc,
		 				kH,
		 				kPi,
		 				kRho,
		 				kPlastic_Modulus_K,
		        kNUM_BCJ_MATL_PARAMS };

		void Allocate (void) { n_mp = kNUM_BCJ_MATL_PARAMS; Parameter.Dimension ( n_mp ); }

};

}
#endif
