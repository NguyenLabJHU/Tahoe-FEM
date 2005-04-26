/* $Id: MFGPFractureSupportT.h,v 1.1 2005-04-26 22:32:28 kyonten Exp $ */
#ifndef _MFGP_FRACTURE_T_H_
#define _MFGP_FRACTURE_T_H_

/* base class */
#include "MFGPElementSupportT.h"
//#include "MeshFreeFractureSupportT.h"

/* direct members */
#include "dArray2DT.h"
#include "nVariArray2DT.h"
#include "dSymMatrixT.h"

namespace Tahoe {

/** support for meshfree calculations including representation of cracks using
 * cutting surfaces */
class MFGPFractureSupportT: public MFGPElementSupportT//, public MeshFreeFractureSupportT
{
public:

	/** constructor */
	MFGPFractureSupportT(void);

	/** destructor */
	//virtual ~MFGPFractureSupportT(void);
	
	/** initialization of meshless information. This method must be called once after 
	 * a call to MeshFreeElementSupportT::TakeParameterList */
	virtual void InitSupport(ostream& out, AutoArrayT<ElementCardT>& elem_cards_displ, 
		AutoArrayT<ElementCardT>& elem_cards_plast, const iArrayT& surface_nodes, 
		int numDOF_displ, int numDOF_plast, int max_node_num, ModelManagerT* model);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	//virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	//virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	//virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/
};

} // namespace Tahoe 
#endif /* _MFGP_FRACTURE_T_H_ */
