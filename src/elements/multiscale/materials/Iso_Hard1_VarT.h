// $Id: Iso_Hard1_VarT.h,v 1.1 2003-03-07 22:27:54 creigh Exp $
#ifndef _ISO_HARD1_T_ 
#define _ISO_HARD1_T_ 

#include "VMF_MaterialT.h"

namespace Tahoe {

class Iso_Hard1_VarT : public VMF_VariableT 
{
	public:

		 Iso_Hard1_VarT ( void ) { Allocate(); }
		~Iso_Hard1_VarT ( void ) { } 

		enum ParamT { 
						kAlpha, 		
						kAlpha_n, 		
		        kNUM_ISO_HARD1_PARAMS };

		void Allocate (void) { n_mp = kNUM_ISO_HARD1_PARAMS; Parameter.Dimension ( n_mp ); }

};

}
#endif
