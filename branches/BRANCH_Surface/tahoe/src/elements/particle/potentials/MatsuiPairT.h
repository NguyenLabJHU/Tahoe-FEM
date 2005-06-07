/* $Id: MatsuiPairT.h,v 1.3 2004-07-15 08:29:49 paklein Exp $ */
#ifndef _MATSUI_PAIR_T_H_
#define _MATSUI_PAIR_T_H_

/* base class */
#include "PairPropertyT.h"

/* define NULL */
#include <stdlib.h>

namespace Tahoe {

/* Matsui pair interaction */
class MatsuiPairT: public PairPropertyT
{
public:

	/** constructor */
	MatsuiPairT(double mass, double sqr_q, double two_A, double two_B, double sqr_C, double f, double rc);
	MatsuiPairT(void);

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
	/** charge */
	double f_sqr_q;
	
	/** A1 + A2 */
	double f_two_A;

	/** B1 + B2 */
	double f_two_B;
	
	/** C1 * C2 */
	double f_sqr_C;

	/** standard force */
	double f_f;
	
	/** cut-off radius */
	double f_rc;
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
	static double s_sqr_q;
	static double s_two_A;
	static double s_two_B;
	static double s_sqr_C;
	static double s_f;
	static double s_rc;
	static double s_phi_rc;
	static double s_dphi_rc;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _HARMONIC_PAIR_T_H_ */
