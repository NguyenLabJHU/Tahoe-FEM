// $Id: BCJ_MatlT.h,v 1.5 2003-04-23 23:34:26 creigh Exp $
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

						kE, 			// 	Young's Modulus
						kPr,			//	Poisson's	Ratio
						kLamda,
		 				kMu, 			// 	Shear Modulus

						kE1, 		  // 	for CCba Young's Modulus
						kPr1,		  //	for CCba Poisson's Ratio
						kLamda1, 	//  For CCba
		 				kMu1, 		//  For CCba	
						kE2, 		  // 	for CCba Young's Modulus
						kPr2,		  //	for CCba Poisson's Ratio
						kLamda2, 	//  For CCba
		 				kMu2, 		//  For CCba	
		 				kGamma_b, //  For CCba	
		 				kAlphaY,  //  For CCba	

						kBulk,		//  Bulk Modulus
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
