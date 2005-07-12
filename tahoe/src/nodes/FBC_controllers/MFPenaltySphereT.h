/* $Id: MFPenaltySphereT.h,v 1.7.20.1 2005-07-02 22:50:38 paklein Exp $ */
/* created: paklein (04/17/2000) */
#ifndef _MF_PENALTY_SPHERE_T_H_
#define _MF_PENALTY_SPHERE_T_H_

/* base class */
#include "PenaltySphereT.h"

namespace Tahoe {

/* forward declarations */
class ElementBaseT;

class MFPenaltySphereT: public PenaltySphereT
{
public:

	/* constructor */
	MFPenaltySphereT(void);

	/* system contributions */
	//virtual void ApplyLHS(void);
	//TEMP - not quite right, but leave it for now

	/** (re-)set the configuration */
	virtual GlobalT::InitStatusT UpdateConfiguration(void);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/
	
protected:

	/* accumulate the contact force vector fContactForce */
	virtual void ComputeContactForce(double kforce);

private:

	/* get element group pointer */
	void SetElementGroup(void);

protected:

	/* element group */
	const ElementBaseT* fElementGroup;
	
	/* work space */
	dArray2DT fInitCoords;
	dArray2DT fCurrCoords;

	/* need MeshFreeSupportT to do this right */
	/* quick and dirty -> request disp from element group */
	/* check if element groups has interpolant DOF's */

};

} // namespace Tahoe 
#endif /* _MF_PENALTY_SPHERE_T_H_ */
