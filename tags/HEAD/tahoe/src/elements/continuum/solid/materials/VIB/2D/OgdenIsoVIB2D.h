/* $Id: OgdenIsoVIB2D.h,v 1.1.1.1 2001-01-29 08:20:24 paklein Exp $ */
/* created: paklein (11/08/1997)                                          */
/* 2D Isotropic VIB using Ogden's spectral formulation                    */

#ifndef _OGDEN_ISO_VIB_2D_H_
#define _OGDEN_ISO_VIB_2D_H_

/* base classes */
#include "OgdenIsotropicT.h"
#include "Material2DT.h"
#include "VIB.h"

/* forward declarations */
class CirclePointsT;

class OgdenIsoVIB2D: public OgdenIsotropicT, public Material2DT, public VIB
{
public:

	/* constructor */
	OgdenIsoVIB2D(ifstreamT& in, const ElasticT& element);

	/* destructor */
	~OgdenIsoVIB2D(void);
	
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

	/* return true of model is purely 2D, plain stress */
	virtual bool PurePlainStress(void) const { return true; };

	/* strained lengths in terms of the Lagrangian stretch eigenvalues */
	void ComputeLengths(const dArrayT& eigs);

private:

	/* initialize angle tables */
	void Construct(void);

protected:
	
	/* integration point generator */
	CirclePointsT*	fCircle;
};

#endif /* _OGDEN_ISO_VIB_2D_H_ */
