/* $Id: XuNeedleman2DT.h,v 1.1.1.1 2001-01-29 08:20:38 paklein Exp $ */
/* created: paklein (11/14/1997)                                          */
/* Xu-Needleman 2D cohesive surface potential                             */

#ifndef _XU_NEEDLE_2D_T_H_
#define _XU_NEEDLE_2D_T_H_

/* base class */
#include "SurfacePotentialT.h"

/* forward declarations */
class ifstreamT;

class XuNeedleman2DT: public SurfacePotentialT
{
public:

	/* constructor */
	XuNeedleman2DT(ifstreamT& in);

	/* surface potential */
	virtual double Potential(const dArrayT& jump_u);
	
	/* traction vector given displacement jump vector */	
	virtual const dArrayT& Traction(const dArrayT& jump_u);

	/* potential stiffness */
	virtual const dMatrixT& Stiffness(const dArrayT& jump_u);

	/* surface status */
	virtual StatusT Status(const dArrayT& jump_u);

	/* print parameters to the output stream */
	virtual void PrintName(ostream& out) const;
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

#endif /* _XU_NEEDLE_2D_T_H_ */
