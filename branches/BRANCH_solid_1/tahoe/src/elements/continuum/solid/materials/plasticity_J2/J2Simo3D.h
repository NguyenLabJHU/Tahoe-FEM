/* $Id: J2Simo3D.h,v 1.2.2.1 2001-06-13 00:08:48 paklein Exp $ */
/* created: paklein (04/30/2001)                                          */

#ifndef _J2_SIMO_3D_H_
#define _J2_SIMO_3D_H_

/* base classes */
#include "SimoIso3D.h"
//#include "J2SimoLinHardT.h"
#include "J2SimoC0HardeningT.h"

/* direct members */
#include "LocalArrayT.h"

//class J2Simo3D: public SimoIso3D, public J2SimoHardeningT
class J2Simo3D: public SimoIso3D, public J2SimoC0HardeningT
{
public:

	/** constructor */
	J2Simo3D(ifstreamT& in, const ElasticT& element);

	/** form of tangent matrix (symmetric by default) */
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/** update internal variables */
	virtual void UpdateHistory(void);

	/** reset internal variables to last converged solution */
	virtual void ResetHistory(void);

	/** print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;
	
	/** modulus */
	virtual const dMatrixT& c_ijkl(void);
	
	/** stress */
	virtual const dSymMatrixT& s_ij(void);

	/** returns the strain energy density for the specified strain */
	virtual double StrainEnergyDensity(void);
	 	 	
	/** required parameter flags */
	virtual bool NeedLastDisp(void) const;

	/** returns the number of output variables */
	virtual int NumOutputVariables(void) const;

	/** returns labels for output variables */
	virtual void OutputLabels(ArrayT<StringT>& labels) const;

	/** compute output variables */
	virtual void ComputeOutput(dArrayT& output);

private:

	/** compute F_total and f_relative */
	void ComputeGradients(void);

private:

	/* deformation gradients */
	dMatrixT fFtot;
	dMatrixT ffrel;
	
	/* work space */
	dMatrixT fF_temp;
};

#endif /* _J2_SIMO_3D_H_ */
