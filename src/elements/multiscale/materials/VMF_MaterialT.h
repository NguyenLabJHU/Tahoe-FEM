// $Id: VMF_MaterialT.h,v 1.3 2003-03-07 22:24:04 creigh Exp $
#ifndef _VMF_MATERIALT_
#define _VMF_MATERIALT_

#include "dArrayT.h"

namespace Tahoe {

class VMF_MaterialT 
{
	public:

		VMF_MaterialT	( void ); 
	  virtual ~VMF_MaterialT() { } 

		//---- Pure Virtual Functions
		
		virtual void 		Allocate ( void )	=0;

		/** ---- Normal Methods (don't use const type arg&, want to send (kE,29.0000) etc.)
		 * these methods are not called enough to affect speed. */

		void		Assign 		( int, double ); 	
		void		Assign 		( double* );
		void		Print 		( void 	);
		double& Retrieve 	( int 	); 				
		void 		Number 		( int		); 
		int& 		Number 		( void 	); 
		void    E_Nu_2_Lamda_Mu ( int iE,int iNu,int iLamda,int iMu);

		enum dParamT { 
						kVolumeColor,
						kEdgeColor,
						kReflectivity,
						kGlossiness,
		        kNUM_VMF_MATERIAL_dCOLOR_PARAMS };
		
		enum iParamT { 
						kTexture,
						kDithering,
		        kNUM_VMF_MATERIAL_iCOLOR_PARAMS };
		
	//protected:

		//dArrayT iColor_Paramaters;
		//dArrayT dColor_Paramaters;
		
		dArrayT Parameter;				
		int matl_no;
		int n_mp; // Number of material parameters (set by derived class during allocation) 

};


class VMF_VariableT : public VMF_MaterialT { public: };  // Exact same thing different name

}
#endif
