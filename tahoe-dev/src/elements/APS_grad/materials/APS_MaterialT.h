// $Id: APS_MaterialT.h,v 1.2 2003-09-03 00:02:57 raregue Exp $
#ifndef _APS_MATERIALT_
#define _APS_MATERIALT_

#include "dArrayT.h"

namespace Tahoe {

class APS_MaterialT 
{
	public:

		APS_MaterialT	( void ); 
		virtual ~APS_MaterialT() { } 
		
		virtual void Allocate ( void )	=0;

		/** ---- Normal Methods (don't use const type arg&, want to send (kE,29.0000) etc.)
		 * these methods are not called enough to affect speed. */

		void		Assign 		( int, double ); 	
		void		Assign 		( double* );
		void		Print 		( void 	);
		double& 	Retrieve 	( int 	); 				
		void 		Number 		( int	); 
		int& 		Number 		( void 	);

		enum dParamT { 
						kVolumeColor,
						kEdgeColor,
						kReflectivity,
						kGlossiness,
		        		kNUM_APS_MATERIAL_dCOLOR_PARAMS };
		
		enum iParamT { 
						kTexture,
						kDithering,
		        		kNUM_APS_MATERIAL_iCOLOR_PARAMS };
		
	//protected:

		//dArrayT iColor_Paramaters;
		//dArrayT dColor_Paramaters;
		
		dArrayT Parameter;				
		int matl_no;
		int n_mp; // Number of material parameters (set by derived class during allocation) 

};

#if 0
class APS_VariableT : public APS_MaterialT { public: };  // Exact same thing different name
#endif

}
#endif
