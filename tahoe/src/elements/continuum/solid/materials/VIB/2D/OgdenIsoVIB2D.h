/* $Id: OgdenIsoVIB2D.h,v 1.9.30.1 2004-03-02 17:46:19 paklein Exp $ */
/* created: paklein (11/08/1997) */
#ifndef _OGDEN_ISO_VIB_2D_H_
#define _OGDEN_ISO_VIB_2D_H_

/* base classes */
#include "OgdenIsotropicT.h"
#include "VIB.h"

namespace Tahoe {

/* forward declarations */
class CirclePointsT;

/** 2D Isotropic VIB using Ogden's spectral formulation */
class OgdenIsoVIB2D: public OgdenIsotropicT, public VIB
{
public:

	/* constructor */
	OgdenIsoVIB2D(ifstreamT& in, const FSMatSupportT& support);

	/* destructor */
	~OgdenIsoVIB2D(void);
	
	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;

	/* strain energy density */
	virtual double StrainEnergyDensity(void);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;
	/*@}*/

protected:

	/* principal values given principal values of the stretch tensors,
	 * i.e., the principal stretches squared */
	virtual void dWdE(const dArrayT& eigenstretch2, dArrayT& eigenstress);
	virtual void ddWddE(const dArrayT& eigenstretch2, dArrayT& eigenstress,
		dSymMatrixT& eigenmod);

	/* return true of model is purely 2D, plain stress */
	virtual bool PurePlaneStress(void) const { return true; };

	/* strained lengths in terms of the Lagrangian stretch eigenvalues */
	void ComputeLengths(const dArrayT& eigs);

private:

	/* initialize angle tables */
	void Construct(void);

protected:
	
	/* integration point generator */
	CirclePointsT*	fCircle;
};

} // namespace Tahoe 
#endif /* _OGDEN_ISO_VIB_2D_H_ */
