/* $Id: XuNeedleman3DT.h,v 1.3 2001-10-11 00:53:42 paklein Exp $ */
/* created: paklein (06/23/1999) */

#ifndef _XU_NEEDLE_3D_T_H_
#define _XU_NEEDLE_3D_T_H_

/* base class */
#include "SurfacePotentialT.h"

/* forward declarations */
class ifstreamT;

/** Xu-Needleman 3D cohesive surface potential */
class XuNeedleman3DT: public SurfacePotentialT
{
public:

	/** constructor */
	XuNeedleman3DT(ifstreamT& in);

	/** return the number of state variables needed by the model */
	int NumStateVariables(void) const { return 0; };

	/** dissipated energy */
	virtual double FractureEnergy(void);

	/** potential energy */
	virtual double Potential(const dArrayT& jump_u, const dArrayT& state);
	
	/** surface traction. Internal variables are integrated over the current
	 * time step. */	
	virtual const dArrayT& Traction(const dArrayT& jump_u, dArrayT& state);

	/** tangent stiffness */
	virtual const dMatrixT& Stiffness(const dArrayT& jump_u, const dArrayT& state);

	/** surface status */
	virtual StatusT Status(const dArrayT& jump_u, const dArrayT& state);

	/** write model name to output */
	virtual void PrintName(ostream& out) const;

	/** write model parameters */
	virtual void Print(ostream& out) const;
	
private:

	/* traction potential parameters */
	double q; // phi_t/phi_n
	double r; // delta_n* /d_n
	
	double d_n; // characteristic normal opening
	double d_t; // characteristic tangent opening
	
	double phi_n;  // mode I work to fracture
double r_fail; // d/d_(n/t) for which surface is considered failed

/* additional penetration stiffness */
double fKratio; // stiffening ratio
double fK;      // penetration stiffness
};

#endif /* _XU_NEEDLE_3D_T_H_ */
