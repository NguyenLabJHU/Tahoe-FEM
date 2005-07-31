/* $Id: SurfacePotentialT.h,v 1.1.1.1 2001-01-29 08:20:38 paklein Exp $ */
/* created: paklein (06/20/1999)                                          */
/* base class for surface potential with jump vector arguments            */

#ifndef _SURFACE_POTENTIAL_T_H_
#define _SURFACE_POTENTIAL_T_H_

/* direct members */
#include "dArrayT.h"
#include "dMatrixT.h"

/* forward declarations */
#include "ios_fwd_decl.h"
class StringT;

class SurfacePotentialT
{
public:

	/* surface potential types - derived classes */
	enum CodeT {kXuNeedleman = 0,
	    kTvergaardHutchinson = 1,
	           kLinearDamage = 2};

	/* status codes */
	enum StatusT {Precritical = 0,
	                 Critical = 1,
	                   Failed = 2};

	/* constructor */
	SurfacePotentialT(int ndof);

	/* destructor */
	virtual ~SurfacePotentialT(void);
	
	/* surface potential */
	virtual double Potential(const dArrayT& jump_u) = 0;
	
	/* traction vector given displacement jump vector */	
	virtual const dArrayT& Traction(const dArrayT& jump_u) = 0;

	/* potential stiffness */
	virtual const dMatrixT& Stiffness(const dArrayT& jump_u) = 0;

	/* surface status */
	virtual StatusT Status(const dArrayT& jump_u) = 0;
	
	/* print parameters to the output stream */
	virtual void PrintName(ostream& out) const = 0;
	virtual void Print(ostream& out) const = 0;

	/* returns true if two materials have compatible nodal outputs */
	static bool CompatibleOutput(const SurfacePotentialT&, const SurfacePotentialT&);

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

protected:

	/* return values */
	dArrayT  fTraction;
	dMatrixT fStiffness;
};

#endif /* _SURFACE_POTENTIAL_T_H_ */
