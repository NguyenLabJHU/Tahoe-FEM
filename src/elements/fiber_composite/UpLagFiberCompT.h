/* $Id: UpLagFiberCompT.h,v 1.2 2006-10-31 15:54:10 rjones Exp $ */
/* created: paklein (07/03/1996) */

#ifndef _UPLAG_FIB_COMP_T_
#define _UPLAG_FIB_COMP_T_

/* base class */
#include "UpdatedLagrangianT.h"

/* direct members */
#include "dArray2DT.h"

namespace Tahoe {

class FSFiberMatSupportT;

/** Update Lagrangian fiber composite element.  Fibers families are assumed   *
	to line in the plane of the element. Fiber orientations are specified in  *
	input using sidesets and can be specified either in global (lab) or       *
	local (parent) element coordinates.                                       */
class UpLagFiberCompT: public UpdatedLagrangianT
{
public:

	/** constructor */
	UpLagFiberCompT(const ElementSupportT& support);
	
	/**destructor**/
	~UpLagFiberCompT(void);

	/* get fiber orientation vectors*/
	const dArray2DT& FiberVec(const int elem);
	
	/*number of fiber orientations*/
	const int NumFibers(const int elem);

	/* information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/* a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

	/** extract the list of material parameters */
	virtual void CollectMaterialInfo(const ParameterListT& all_params, ParameterListT& mat_params) const;

	enum CoorSysT {kCartesian =0,
	               kPolar = 1};
protected:

	/*Rotates fiber coordinates from local coordinates to element coordinates*/ 
	void ReadFiberVec(const ParameterListT& list);

	/*reads orientations from side sets*/
	void ReadSideSetVec(const ParameterListT& fibers);

	/*reads orientations from an analytical surface*/
	void ReadAnalyticVec(const ParameterListT& fibers);


	/** construct a new material support and return a pointer. Recipient is responsible for
	 * for freeing the pointer.
	 * \param p an existing MaterialSupportT to be initialized. If NULL, allocate
	 *        a new MaterialSupportT and initialize it. */
	virtual MaterialSupportT* NewMaterialSupport(MaterialSupportT* p = NULL) const;

	/** return a pointer to a new material list. Recipient is responsible for freeing 
	 * the pointer. 
	 * \param name list identifier
	 * \param size length of the list */
	virtual MaterialListT* NewMaterialList(const StringT& name, int size);

protected:
	/*element list of fiber orientation vectors in global (lab) coordinates*/
	/*num_elem< num_fibers x nsd >*/
	ArrayT<dArray2DT> fFiber_list;
	
	/*material support*/
	FSFiberMatSupportT* fFiberSupport;

	/** element block ID */
	StringT fID;
};

inline const dArray2DT& UpLagFiberCompT::FiberVec(const int elem)
{
	if (fFiber_list[elem].MajorDim() < 1)
		ExceptionT::GeneralFail("UpLagFiberCompT::FiberVect", 
			"No fiber orientations specified for element %d", elem);

	return(fFiber_list[elem]);
}

inline const int UpLagFiberCompT::NumFibers(const int elem)
{
	return(fFiber_list[elem].MajorDim());
}

} // namespace Tahoe 
#endif /* _UPLAG_FIB_COMP_T_ */
