/* $Id: LennardJonesPairT.h,v 1.7 2004-07-15 08:29:49 paklein Exp $ */
#ifndef _LENNARD_JONES_PAIR_T_H_
#define _LENNARD_JONES_PAIR_T_H_

/* base class */
#include "PairPropertyT.h"

/* define NULL */
#include <stdlib.h>

namespace Tahoe {

/** Lennard-Jones pair interaction. The Lennard-Jones 6/12 with a smooth
 * cut-off function. The potential is defined as
 \f[
	\phi(r) = \phi_{LJ}(r) - \phi_{LJ}(r_c) - (r - r_c) \phi'_{LJ}(r_c),
 \f]
 * where \f$ r_c = \alpha \sigma \f$. The unmodified Lennard Jones potential is
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
 The cut-off parameter \f$ \alpha \f$ is ignored if \f$ \alpha < 0 \f$.
 */
class LennardJonesPairT: public PairPropertyT
{
public:

	/** constructor */
	LennardJonesPairT(double mass, double eps, double sigma, double alpha);
	LennardJonesPairT(void);

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

	/** \name user-defined parameters */
	/*@{*/
	/** energy scaling */
	double f_eps;

	/** length scaling */
	double f_sigma;
	
	/** cut-off distance */
	double f_alpha;
	/*@}*/
	
	/** \name unmodified potential at the cut-off */
	/*@{*/
	/** energy */
	double f_phi_rc;

	/** force */
	double f_dphi_rc;
	/*@}*/
	
	/** \name static parameters 
	 * There parameters are use during evaluation of the static interaction 
	 * functions. These are copied to static when a function pointer is requested */
	/*@{*/
	static double s_eps;
	static double s_sigma;
	static double s_alpha;
	static double s_phi_rc;
	static double s_dphi_rc;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _HARMONIC_PAIR_T_H_ */
