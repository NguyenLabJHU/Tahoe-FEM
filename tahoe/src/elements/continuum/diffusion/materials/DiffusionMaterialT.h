/* $Id: DiffusionMaterialT.h,v 1.5 2002-11-14 17:06:39 paklein Exp $ */
/* created: paklein (10/02/1999) */
#ifndef _DIFFUSION_MATERIALT_H_
#define _DIFFUSION_MATERIALT_H_

/* base class */
#include "ContinuumMaterialT.h"

/* direct members */
#include "dMatrixT.h"
#include "dArrayT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;
class DiffusionT;
class DiffusionMatSupportT;

/** interface for materials for diffusion */
class DiffusionMaterialT: public ContinuumMaterialT
{
public:

	/** constructor */
	DiffusionMaterialT(ifstreamT& in, const DiffusionMatSupportT& support);

	/** \name print parameters */
	/*@{*/
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;
	/*@}*/

	/** \name parameters at the current field point */
	/*@{*/
	/** conductivity */
	const dMatrixT& k_ij(void);

	/** heat flux */
	const dArrayT& q_i(void);

	double Density(void) const;
	double SpecificHeat(void) const;
	double Capacity(void) const;
	/*@}*/

protected:

	/** support for diffusion materials */
	const DiffusionMatSupportT& fDiffusionMatSupport;

	/* parameters */
	double   fDensity;
	double   fSpecificHeat;
	dMatrixT fConductivity;	// always symmetric?
	
	/* derived */
	double fCapacity;

	/** heat flux return value*/
	dArrayT fq_i;  
};

/* returns the density */
inline double DiffusionMaterialT::Density(void) const { return fDensity; }
inline double DiffusionMaterialT::SpecificHeat(void) const { return fSpecificHeat; }
inline double DiffusionMaterialT::Capacity(void) const { return fCapacity; }

/* conductivity */
inline const dMatrixT& DiffusionMaterialT::k_ij(void) { return fConductivity; }

} // namespace Tahoe 
#endif /* _DIFFUSION_MATERIALT_H_ */
