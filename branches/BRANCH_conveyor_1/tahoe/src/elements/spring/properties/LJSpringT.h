/* $Id: LJSpringT.h,v 1.4 2002-11-30 16:33:42 paklein Exp $ */
/* created: paklein (11/20/1996) */
#ifndef _LJ_SPRINGT_H_
#define _LJ_SPRINGT_H_

/* base class */
#include "RodMaterialT.h"

namespace Tahoe {

/* forward declarations */
class ThermalDilatationT;
class ElementBaseT;

/** Lennard-Jones 6/12.
The unmodified Lennard Jones potential is
 \f[
	\phi_{LJ}(r) = 4 \epsilon \left[ (\sigma/r)^{12} - (\sigma/r)^{6} \right].
 \f]
 * In terms of these parameters, equilibrium length of a single, unmodified
 * Lennard-Jones bond is
 \f[
 	r_0 = 2^{1/6} \sigma,
 \f]
 and the depth of the energy well is
 \f[
	\phi_{LJ}(r_0) = -\epsilon.
 \f]
 */
class LJSpringT: public RodMaterialT
{
public:

	/* constructor */
	LJSpringT(ifstreamT& in);

	/* returns trues TRUE since the initial length is always assumed
	 * to be non-equilibrium */
	virtual int HasInternalStrain(void) const;
	
	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;
	
	/* potential function and derivatives */
	virtual double Potential(double rmag, double Rmag) const;
	virtual double DPotential(double rmag, double Rmag) const;
	virtual double DDPotential(double rmag, double Rmag) const;
	
private:

	/** \name user-defined parameters */
	/*@{*/
	/** energy scaling */
	double f_eps;

	/** length scaling */
	double f_sigma;
	/*@}*/	
};

} // namespace Tahoe 
#endif /* _LJ_SPRINGT_H_ */
