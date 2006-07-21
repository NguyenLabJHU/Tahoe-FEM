/* $Id: OgdenIsoVIB2D.h,v 1.10 2004-07-15 08:27:45 paklein Exp $ */
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
	OgdenIsoVIB2D(void);

	/* destructor */
	~OgdenIsoVIB2D(void);
	
	/* strain energy density */
	virtual double StrainEnergyDensity(void);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
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
