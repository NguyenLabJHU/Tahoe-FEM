/* $Id: PenaltySphereT.h,v 1.1.1.1 2001-01-29 08:20:40 paklein Exp $ */
/* created: paklein (04/30/1998)                                          */

#ifndef _PENATLY_SPHERE_T_H_
#define _PENATLY_SPHERE_T_H_

/* base class */
#include "PenaltyRegionT.h"

/* direct members */
#include "ElementMatrixT.h"

class PenaltySphereT: public PenaltyRegionT
{
public:

	/* constructor */
	PenaltySphereT(FEManagerT& fe_manager, const iArray2DT& eqnos, const dArray2DT& coords,
		const dArray2DT* vels);

	/* input processing */
	virtual void EchoData(ifstreamT& in, ostream& out);

	/* initialize data */
	virtual void Initialize(void);

	/* form of tangent matrix */
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/* tangent term */
	virtual void ApplyLHS(void);
	
protected:

	/* accumulate the contact force vector fContactForce */
	virtual void ComputeContactForce(double kforce);

protected:

	/* sphere radius */
	double fRadius;
	
	/* center to striker distances */
	dArrayT	fDistances;
	
	/* workspace */
	dArrayT        fv_OP; //vector from center to contact node
	ElementMatrixT fLHS;  //tangent matrix
	dArrayT        fd_sh; //shallow
	iArrayT        fi_sh; //shallow
};

#endif /* _PENATLY_SPHERE_T_H_ */
