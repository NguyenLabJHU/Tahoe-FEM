// $Id: APS_EnumT.h,v 1.8 2004-02-04 00:40:44 raregue Exp $
#ifndef _APS_ENUM_H_ 
#define _APS_ENUM_H_ 

namespace Tahoe {

/** APS_Variable_ListT: This class contains the names of the variables used in 
 * APS formulation. Enums can be accessed from anyone/anywhere by APS::kgrad_u etc. 
 * The APS:: is always mandatory since the enums are members of this class
 * only (they're not global).
**/

class APS
{

  public:

    enum VarT_vector {  
					kgammap, 
					kgammap_surf, 
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



