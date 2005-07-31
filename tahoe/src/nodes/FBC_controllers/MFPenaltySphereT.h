/* $Id: MFPenaltySphereT.h,v 1.1.1.1 2001-01-29 08:20:40 paklein Exp $ */
/* created: paklein (04/17/2000)                                          */

#ifndef _MF_PENALTY_SPHERE_T_H_
#define _MF_PENALTY_SPHERE_T_H_

/* base class */
#include "PenaltySphereT.h"

/* forward declarations */
class ElementBaseT;

class MFPenaltySphereT: public PenaltySphereT
{
public:

	/* constructor */
	MFPenaltySphereT(FEManagerT& fe_manager, const iArray2DT& eqnos,
		const dArray2DT& coords, const dArray2DT* vels);

	/* input processing */
	virtual void EchoData(ifstreamT& in, ostream& out);

	/* initialize data */
	virtual void Initialize(void);

	/* system contributions */
	//virtual void ApplyLHS(void);
	//TEMP - not quite right, but leave it for now
	
protected:

	/* accumulate the contact force vector fContactForce */
	virtual void ComputeContactForce(double kforce);

private:

	/* get element group pointer */
	void SetElementGroup(void);

protected:

	/* element group */
	int fGroupNumber;
	const ElementBaseT* fElementGroup;
	
	/* work space */
	dArray2DT fCoords;
	dArray2DT fCurrCoords;

	/* need MeshFreeSupportT to do this right */
	/* quick and dirty -> request disp from element group */
	/* check if element groups has interpolant DOF's */

};

#endif /* _MF_PENALTY_SPHERE_T_H_ */
