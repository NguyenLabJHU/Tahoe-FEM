/* $Id: PenaltySphereT.h,v 1.2.2.1 2002-06-27 18:03:58 cjkimme Exp $ */
/* created: paklein (04/30/1998) */

#ifndef _PENATLY_SPHERE_T_H_
#define _PENATLY_SPHERE_T_H_

/* base class */
#include "PenaltyRegionT.h"

/* direct members */
#include "ElementMatrixT.h"


namespace Tahoe {

class PenaltySphereT: public PenaltyRegionT
{
public:

	/* constructor */
	PenaltySphereT(FEManagerT& fe_manager, int group, const iArray2DT& eqnos, 
		const dArray2DT& coords, const dArray2DT* vels);

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

} // namespace Tahoe 
#endif /* _PENATLY_SPHERE_T_H_ */
