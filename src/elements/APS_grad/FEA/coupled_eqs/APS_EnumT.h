// $Id: APS_EnumT.h,v 1.2 2003-09-19 00:47:02 raregue Exp $
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

    enum VarT {  	
					kgrad_u, 
					kgrad_gammap, 
	                kNUM_APS_VARS }; // <-- Keep this one last !!
};

//---------------------------------------------------------------------
} // namespace Tahoe 

#endif /* _APS_ENUM_H_ */



