// $Id: MFGP_EnumT.h
#ifndef _MFGP_ENUM_H_ 
#define _MFGP_ENUM_H_ 

namespace Tahoe 
{

/** MFGP_Variable_ListT: This class contains the names of the variables used in 
 * MFGP formulation. Enums can be accessed from anyone/anywhere by MFGP::kgrad_u etc. 
 * The MFGP:: is always mandatory since the enums are members of this class
 * only (they're not global).
**/

class MFGP
{

  public:

    enum VarT_vector {  
					//kgammap,  
					kstate,
	                kNUM_MFGP_VECTOR_VARS }; // <-- Keep this one last !!
	                
	enum VarT_matrix {  
					//kgrad_u,  	 
					//kgrad_gammap, 
	                kNUM_MFGP_MATRIX_VARS }; // <-- Keep this one last !!
};

//---------------------------------------------------------------------
} // namespace Tahoe 

#endif /* _MFGP_ENUM_H_ */



