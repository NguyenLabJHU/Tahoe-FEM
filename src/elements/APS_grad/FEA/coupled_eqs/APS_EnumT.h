// $Id: APS_EnumT.h,v 1.6 2003-10-06 18:34:39 raregue Exp $
#ifndef _APS_ENUM_H_ 
#define _APS_ENUM_H_ 

namespace Tahoe {

/** APS_Variable_ListT: This class contains the names of the variables used in 
 * APS formulation. Enums 
 * can be accessed from anyone/anywhere by APS::kgrad_u etc. 
 * The APS:: is always mandatory since the enums are members of this class
 * only (there not global).
**/

class APS
{

  public:

    enum VarT_vector {  
    				//grad_u must be treated as a matrix given the GradU operation in ShapeFunctionT	
					//kgrad_u, 
					kgammap, 
					kstate,
	                kNUM_APS_VECTOR_VARS }; // <-- Keep this one last !!
	                
	enum VarT_matrix {  
					kgrad_u, 
					kgrad_u_surf, 	 
					kgrad_gammap, 
	                kNUM_APS_MATRIX_VARS }; // <-- Keep this one last !!
};

//---------------------------------------------------------------------
} // namespace Tahoe 

#endif /* _APS_ENUM_H_ */



