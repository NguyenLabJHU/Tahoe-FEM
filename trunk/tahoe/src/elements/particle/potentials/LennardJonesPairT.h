/* $Id: LennardJonesPairT.h,v 1.1 2002-11-25 07:19:46 paklein Exp $ */
#ifndef _LENNARD_JONES_PAIR_T_H_
#define _LENNARD_JONES_PAIR_T_H_

/* base class */
#include "PairPropertyT.h"

namespace Tahoe {

/** harmonic pair interaction */
class LennardJonesPairT: public PairPropertyT
{
public:

	/** constructor */
	LennardJonesPairT(double eps, double sigma, double cut_off);

	/** \name return interaction functions */
	/*@{*/
	/** return a pointer to the energy function */
	virtual EnergyFunction getEnergyFunction(void);

	/** return a pointer to the force function */
	virtual ForceFunction getForceFunction(void);

	/** return a pointer to the stiffness function */
	virtual StiffnessFunction getStiffnessFunction(void);
	/*@}*/

private:

	/** \name interaction functions */
	/*@{*/
	static double Energy(double r_ab, double* data_a, double* data_b);
	static double Force(double r_ab, double* data_a, double* data_b);
	static double Stiffness(double r_ab, double* data_a, double* data_b);
	/*@}*/

private:

	/** energy scaling */
	double f_eps;

	/** length scaling */
	double f_sigma;
	
	/** cut-off distance */
	double f_cut;
	
	/** \name static parameters 
	 * There parameters are use during evaluation of the static interaction 
	 * functions. These are copied to static when a function pointer is requested */
	/*@{*/
	static double s_eps;
	static double s_sigma;
	static double s_cut;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _HARMONIC_PAIR_T_H_ */
