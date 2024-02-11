// $Id: E_Pr_MatlT.h,v 1.2 2003-02-03 04:40:28 paklein Exp $
#ifndef _E_PR_MATLT_ 
#define _E_PR_MATLT_ 

#include "VMF_MaterialT.h"

namespace Tahoe {

class E_Pr_MatlT : public VMF_MaterialT 
{
	public:

		 E_Pr_MatlT ( void ) { Allocate(); }
		~E_Pr_MatlT ( void ) { } 

		enum ParamT { 
						kE, 		// YoungsModulus
						kPr, 		// PoissonsRatio
		        kNUM_ISO_MATL_PARAMS };

		void Allocate (void) { n_mp = kNUM_ISO_MATL_PARAMS; Parameter.Dimension ( n_mp ); }

};

}
#endif
