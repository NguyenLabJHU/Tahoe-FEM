/* $Id: TvergHutch2DT.h,v 1.1.1.1 2001-01-29 08:20:38 paklein Exp $ */
/* created: paklein (02/05/2000)                                          */
/* cohesive potential from Tvergaard and Hutchinson,                      */
/* JMPS v41, n6, 1995, 1119-1135.                                         */

#ifndef _TVERG_HUTCH_2D_T_H_
#define _TVERG_HUTCH_2D_T_H_

/* base class */
#include "SurfacePotentialT.h"

/* forward declarations */
class ifstreamT;

class TvergHutch2DT: public SurfacePotentialT
{
public:

	/* constructor */
	TvergHutch2DT(ifstreamT& in);

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

	/* returns the number of variables computed for nodal extrapolation
	 * during for element output, ie. internal variables. Returns 0
	 * by default */
	virtual int NumOutputVariables(void) const; // 0 by default
	virtual void OutputLabels(ArrayT<StringT>& labels) const; // none by default
	virtual void ComputeOutput(const dArrayT& jump_u, dArrayT& output);

protected:

	/* return true if the potential has compatible (type and sequence)
	 * nodal output - FALSE by default */
	virtual bool CompatibleOutput(const SurfacePotentialT& potential) const;
	
private:

	/* traction potential parameters */
	double fsigma_max; // cohesive stress
	double fd_c_n;     // characteristic normal opening to failure
	double fd_c_t;     // characteristic tangential opening to failure
	
	/* non-dimensional opening parameters */
	double fL_1; // opening to initial peak traction
	double fL_2; // opening to final peak traction
	double fL_fail; // opening to irreversible failure

	/* penetration stiffness */
	double fpenalty; // stiffening multiplier
	double fK;       // penetration stiffness
};

#endif /* _TVERG_HUTCH_2D_T_H_ */
