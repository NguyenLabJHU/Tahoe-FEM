/* $Id: HarmonicPairT.h,v 1.5 2004-07-15 08:29:49 paklein Exp $ */
#ifndef _HARMONIC_PAIR_T_H_
#define _HARMONIC_PAIR_T_H_

/* base class */
#include "PairPropertyT.h"

namespace Tahoe {

/** harmonic pair interaction */
class HarmonicPairT: public PairPropertyT
{
public:

	/** constructor */
	HarmonicPairT(double mass, double R0, double K);
	HarmonicPairT(void);

	/** \name return interaction functions */
	/*@{*/
	/** return a pointer to the energy function */
	virtual EnergyFunction getEnergyFunction(void);

	/** return a pointer to the force function */
	virtual ForceFunction getForceFunction(void);

	/** return a pointer to the stiffness function */
	virtual StiffnessFunction getStiffnessFunction(void);
	/*@}*/

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

private:

	/** \name interaction functions */
	/*@{*/
	static double Energy(double r_ab, double* data_a, double* data_b);
	static double Force(double r_ab, double* data_a, double* data_b);
	static double Stiffness(double r_ab, double* data_a, double* data_b);
	/*@}*/

private:

	/** equilibrium spacing */
	double fR0;

	/** stiffness */
	double fK;
	
	/** \name static parameters 
	 * There parameters are use during evaluation of the static interaction 
	 * functions. These are copied to static when a function pointer is requested */
	/*@{*/
	static double sR0;
	static double sK;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _HARMONIC_PAIR_T_H_ */
