// $Id: BCJ_MatlT.h,v 1.2 2003-02-03 04:40:28 paklein Exp $
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
		        kNUM_BCJ_MATL_PARAMS };

		void Allocate (void) { n_mp = kNUM_BCJ_MATL_PARAMS; Parameter.Dimension ( n_mp ); }

};

}
#endif
