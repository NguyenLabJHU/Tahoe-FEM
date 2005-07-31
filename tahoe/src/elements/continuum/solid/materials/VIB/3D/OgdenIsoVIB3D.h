/* $Id: OgdenIsoVIB3D.h,v 1.1.1.1 2001-01-29 08:20:25 paklein Exp $ */
/* created: paklein (11/08/1997)                                          */
/* 3D Isotropic VIB using Ogden's spectral formulation                    */

#ifndef _OGDEN_ISO_VIB_3D_H_
#define _OGDEN_ISO_VIB_3D_H_

/* base classes */
#include "OgdenIsotropicT.h"
#include "VIB.h"

/* forward declarations */
class SpherePointsT;

class OgdenIsoVIB3D: public OgdenIsotropicT, public VIB
{
public:

	/* constructor */
	OgdenIsoVIB3D(ifstreamT& in, const ElasticT& element);

	/* destructor */
	~OgdenIsoVIB3D(void);
	
	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;

	/* strain energy density */
	virtual double StrainEnergyDensity(void);

protected:

	/* principal values given principal stretches */
	virtual void dWdE(const dArrayT& eigenstretch, dArrayT& eigenstress);
	virtual void ddWddE(const dArrayT& eigenstretch, dArrayT& eigenstress,
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

#endif /* _OGDEN_ISO_VIB_3D_H_ */
