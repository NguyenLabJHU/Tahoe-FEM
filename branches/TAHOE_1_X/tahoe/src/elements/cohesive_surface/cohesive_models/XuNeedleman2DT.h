/* $Id: XuNeedleman2DT.h,v 1.10 2003-05-26 01:51:46 paklein Exp $ */
/* created: paklein (11/14/1997) */

#ifndef _XU_NEEDLE_2D_T_H_
#define _XU_NEEDLE_2D_T_H_

/* base class */
#include "SurfacePotentialT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;

/** Xu-Needleman 2D cohesive surface potential */
class XuNeedleman2DT: public SurfacePotentialT
{
public:

	/** constructor */
	XuNeedleman2DT(ifstreamT& in);

	/** return the number of state variables needed by the model */
	int NumStateVariables(void) const { return 0; };

	/** dissipated energy */
	virtual double FractureEnergy(const ArrayT<double>& state);

	/** potential energy */
	virtual double Potential(const dArrayT& jump_u, const ArrayT<double>& state);
	
	/** surface traction. Internal variables are integrated over the current
	 * time step. */	
	virtual const dArrayT& Traction(const dArrayT& jump_u, ArrayT<double>& state, const dArrayT& sigma, bool qIntegrate);

	/** tangent stiffness */
	virtual const dMatrixT& Stiffness(const dArrayT& jump_u, const ArrayT<double>& state, const dArrayT& sigma);

	/** surface status */
	virtual StatusT Status(const dArrayT& jump_u, const ArrayT<double>& state);

	/** write model name to output */
	virtual void PrintName(ostream& out) const;

	/** write model parameters */
	virtual void Print(ostream& out) const;
	
private:

	/* traction potential parameters */
	double q;     // phi_t/phi_n
	double r;     // delta_n* /d_n
	
	double d_n;   // characteristic normal opening
	double d_t;   // characteristic tangent opening  	
	double phi_n; // mode I work to fracture

	double r_fail; // d/d_(n/t) for which surface is considered failed

/* additional penetration stiffness */
double fKratio; // stiffening ratio
double fK;
};

} // namespace Tahoe 
#endif /* _XU_NEEDLE_2D_T_H_ */
