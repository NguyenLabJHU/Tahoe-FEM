/* $Id: OgdenIsoVIB3D.h,v 1.3.6.1 2002-06-27 18:03:21 cjkimme Exp $ */
/* created: paklein (11/08/1997)                                          */
/* 3D Isotropic VIB using Ogden's spectral formulation                    */

#ifndef _OGDEN_ISO_VIB_3D_H_
#define _OGDEN_ISO_VIB_3D_H_

/* base classes */
#include "OgdenIsotropicT.h"
#include "VIB.h"

/* forward declarations */

namespace Tahoe {

class SpherePointsT;

class OgdenIsoVIB3D: public OgdenIsotropicT, public VIB
{
public:

	/* constructor */
	OgdenIsoVIB3D(ifstreamT& in, const FiniteStrainT& element);

	/* destructor */
	~OgdenIsoVIB3D(void);
	
	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;

	/* strain energy density */
	virtual double StrainEnergyDensity(void);

protected:

	/* principal values given principal values of the stretch tensors,
	 * i.e., the principal stretches squared */
	virtual void dWdE(const dArrayT& eigenstretch2, dArrayT& eigenstress);
	virtual void ddWddE(const dArrayT& eigenstretch2, dArrayT& eigenstress,
		dSymMatrixT& eigenmod);

	/* strained lengths in terms of the Lagrangian stretch eigenvalues */
	void ComputeLengths(const dArrayT& eigs);

private:

	/* initialize angle tables */
	void Construct(void);

protected:
	
	/* integration point generator */
	SpherePointsT* fSphere;
};

} // namespace Tahoe 
#endif /* _OGDEN_ISO_VIB_3D_H_ */
