/* $Id: DiffusionMaterialT.h,v 1.1.1.1 2001-01-29 08:20:25 paklein Exp $ */
/* created: paklein (10/02/1999)                                          */
/* Defines the interface for materials for diffusion.                     */

#ifndef _DIFFUSION_MATERIALT_H_
#define _DIFFUSION_MATERIALT_H_

#include "GlobalT.h"

/* base class */
#include "ContinuumMaterialT.h"

/* direct members */
#include "dMatrixT.h"
#include "dArrayT.h"

/* forward declarations */
class ifstreamT;
class ShapeFunctionT;
class LocalArrayT;
class DiffusionT;

class DiffusionMaterialT: public ContinuumMaterialT
{
public:

	/* constructor */
	DiffusionMaterialT(ifstreamT& in, const DiffusionT& element);

	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;

	/* conductivity */
	const dMatrixT& k_ij(void);

	/* heat flux */
	const dArrayT& q_i(void);

	/* returns the density */
	double Density(void) const;
	double SpecificHeat(void) const;
	double Capacity(void) const;
	
protected:

	/* interpolation and gradients */
	const ShapeFunctionT& fShapes;

	/* local displacements */
	const LocalArrayT& fLocDisp;

	/* parameters */
	double   fDensity;
	double   fSpecificHeat;
	dMatrixT fConductivity;	// always symmetric?
	
	/* derived */
	double fCapacity;

	/* return value */
	dMatrixT fT_x; // temperature gradient
	dArrayT fq_i;  // heat flux
};

/* returns the density */
inline double DiffusionMaterialT::Density(void) const { return fDensity; }
inline double DiffusionMaterialT::SpecificHeat(void) const { return fSpecificHeat; }
inline double DiffusionMaterialT::Capacity(void) const { return fCapacity; }

/* conductivity */
inline const dMatrixT& DiffusionMaterialT::k_ij(void) { return fConductivity; }

#endif /* _DIFFUSION_MATERIALT_H_ */
