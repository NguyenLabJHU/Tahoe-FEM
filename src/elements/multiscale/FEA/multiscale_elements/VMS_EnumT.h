//DEVELOPMENT

#ifndef _VMS_ENUM_H_ 
#define _VMS_ENUM_H_ 

namespace Tahoe {

/** VMS_Variable_ListT: This class contains the names of the variables used in 
 * VMS formulation. Enums 
 * can be accessed from anyone/anywhere by VMS::kGRAD_u etc. 
 * The VMS:: is always mandatory since the enums are members of this class
 * only (there not global).
 * Sandia National Laboratory and the University of Michigan **/

class VMS
{

  public:

    enum VarT {  	
						      kGRAD_u, 
						      kGRAD_ua, 
                 	kGRAD_ub, 
                 	kgrad_u, 
                 	kgrad_ua, 
                 	kgrad_ub, 
		 							kF,
	         				kFi, 
	         				kFa, 
	         				kFai, 
                 	kFb, 
                 	kFbi,
	                kNUM_VMS_VARS }; // <-- Keep this one last !!
};

//---------------------------------------------------------------------
} // namespace Tahoe 

#endif /* _VMS_ENUM_H_ */



