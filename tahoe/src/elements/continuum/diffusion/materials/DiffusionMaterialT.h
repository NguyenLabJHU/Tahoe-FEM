/* $Id: DiffusionMaterialT.h,v 1.8 2004-07-15 08:26:22 paklein Exp $ */
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
	DiffusionMaterialT(void);

	/** set support */
	void SetDiffusionMatSupport(const DiffusionMatSupportT* support);

	/** \name parameters at the current field point */
	/*@{*/
	/** conductivity */
	virtual const dMatrixT& k_ij(void);

	/** heat flux */
	virtual const dArrayT& q_i(void);

	/** change in heat flux with temperature */
	virtual const dArrayT& dq_i_dT(void) { return fdq_i; };

	double Density(void) const;
	virtual double SpecificHeat(void) const;
	double Capacity(void) const;

	/** change in capacity with temperature */
	virtual double dCapacity_dT(void) const { return 0; };
	/*@}*/

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** support for diffusion materials */
	const DiffusionMatSupportT* fDiffusionMatSupport;

	/** \name parameters */
	/*@{*/
	double   fDensity;
	double   fSpecificHeat;
	dMatrixT fConductivity;
	/*@}*/

	/** \name return values */
	/*@{*/
	/** heat flux return value */
	dArrayT fq_i;  

	/** heat flux variation return value */
	dArrayT fdq_i;  
	/*@}*/
};

/* returns the density */
inline double DiffusionMaterialT::Density(void) const { return fDensity; }
inline double DiffusionMaterialT::SpecificHeat(void) const { return fSpecificHeat; }
inline double DiffusionMaterialT::Capacity(void) const { return fDensity*SpecificHeat(); }

/* conductivity */
inline const dMatrixT& DiffusionMaterialT::k_ij(void) { return fConductivity; }

} // namespace Tahoe 
#endif /* _DIFFUSION_MATERIALT_H_ */
