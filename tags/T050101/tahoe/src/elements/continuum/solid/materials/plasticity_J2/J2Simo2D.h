/* $Id: J2Simo2D.h,v 1.1.1.1 2001-01-29 08:20:30 paklein Exp $ */
/* created: paklein (06/22/1997)                                          */

#ifndef _J2_SIMO_2D_H_
#define _J2_SIMO_2D_H_

/* base classes */
#include "SimoIso2D.h"
#include "J2SimoLinHardT.h"

/* direct members */
#include "LocalArrayT.h"

class J2Simo2D: public SimoIso2D, public J2SimoLinHardT
{
public:

	/* constructor */
	J2Simo2D(ifstreamT& in, const ElasticT& element);

	/* form of tangent matrix (symmetric by default) */
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/* update internal variables */
	virtual void UpdateHistory(void);

	/* reset internal variables to last converged solution */
	virtual void ResetHistory(void);

	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;
	
	/* modulus */
	virtual const dMatrixT& c_ijkl(void);
	
	/* stress */
	virtual const dSymMatrixT& s_ij(void);

	/* returns the strain energy density for the specified strain */
	virtual double StrainEnergyDensity(void);
	 	 	
	/* required parameter flags */
	virtual bool NeedLastDisp(void) const;

private:

	/* compute F_total and f_relative */
	void ComputeGradients(void);

private:

	/* last converged disp - needed for f_relative */
	const LocalArrayT& fLocLastDisp;
	LocalArrayT	fRelDisp;

	/* deformation gradients - 3D*/
	dMatrixT fFtot;
	dMatrixT ffrel;
	
	/* work space */
	dMatrixT fF_temp;
	dMatrixT fFtot_2D;
	dMatrixT ffrel_2D;
};

#endif /* _J2_SIMO_2D_H_ */
