/* $Id: J2Simo2D.h,v 1.11 2003-10-12 01:39:03 paklein Exp $ */
/* created: paklein (06/22/1997) */
#ifndef _J2_SIMO_2D_H_
#define _J2_SIMO_2D_H_

/* base classes */
#include "SimoIso2D.h"
#include "J2SimoC0HardeningT.h"

/* direct members */
#include "LocalArrayT.h"

namespace Tahoe {

/** finite strain J2 plasticity */
class J2Simo2D: public SimoIso2D, public J2SimoC0HardeningT
{
public:

	/** constructor */
	J2Simo2D(ifstreamT& in, const FSMatSupportT& support);

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

	/** incremental heat generation */
	virtual double IncrementalHeat(void);

	/** this model does generate heat */
	virtual bool HasIncrementalHeat(void) const { return true; };
	 	 	
	/** required parameter flags */
	virtual bool Need_F_last(void) const { return true; };

	/** returns the number of output variables */
	virtual int NumOutputVariables(void) const;

	/** returns labels for output variables */
	virtual void OutputLabels(ArrayT<StringT>& labels) const;

	/** compute output variables */
	virtual void ComputeOutput(dArrayT& output);

private:

	/** flag to indicate whether material supports thermal strains.
	 * Returns true. */
	virtual bool SupportsThermalStrain(void) const { return true; };

	/** compute F_mechanical and f_relative for the current step */
	void ComputeGradients(void);

private:

	/* deformation gradients - 3D*/
	dMatrixT fFmech;
	dMatrixT ffrel;
	
	/* work space */
	dMatrixT fF_temp;
	dMatrixT fFmech_2D;
	dMatrixT ffrel_2D;
};

} // namespace Tahoe 
#endif /* _J2_SIMO_2D_H_ */
