/* $Id: StillingerWeberT.h,v 1.1 2004-11-23 01:43:20 cjkimme Exp $ */
#ifndef _STILLINGER_WEBER_T_H_
#define _STILLINGER_WEBER_T_H_

/* base classes */
#include "ThreeBodyPropertyT.h"

/* define NULL */
#include <stdlib.h>

namespace Tahoe {

class StillingerWeberT: public ThreeBodyPropertyT
{
public:

	/** constructor */
	StillingerWeberT(void);

	/** \name return interaction functions */
	/*@{*/
	/** return a pointer to the 2-body energy function */
	virtual PairPropertyT::EnergyFunction getEnergyFunction(void);

	/** return a pointer to the 2-body force function */
	virtual PairPropertyT::ForceFunction getForceFunction(void);

	/** return a pointer to the 2-body stiffness function */
	virtual PairPropertyT::StiffnessFunction getStiffnessFunction(void);

	/** return a pointer to the 3-body energy function */
	virtual ThreeBodyPropertyT::EnergyFunction getThreeBodyEnergyFunction(void);

	/** return a pointer to the 3-body force function */
	virtual ThreeBodyPropertyT::ForceFunction getThreeBodyForceFunction(void);

	/** return a pointer to the 3-body stiffness function */
	virtual ThreeBodyPropertyT::StiffnessFunction getThreeBodyStiffnessFunction(void);
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
	static double TwoBodyEnergy(double r_ab, double* data_a, double* data_b);
	static double TwoBodyForce(double r_ab, double* data_a, double* data_b);
	static double TwoBodyStiffness(double r_ab, double* data_a, double* data_b);
	static double ThreeBodyEnergy(double* ri, double* rj, double* rk);
	static double* ThreeBodyForce(const double* ri, const double* rj, const double* rk, double* fij, double *fik);
	static double* ThreeBodyStiffness(const double* ri, const double* rj, const double* rk, dMatrixT& K_ijk);
	/*@}*/

private:

	/** \name user-defined parameters */
	/*@{*/
	/** energy scaling */
	double f_eps;

	/** length scaling */
	double f_sigma;
	
	/** cut-off distance */
	double f_b;

	/** 3-body secondary energy scale */
	double f_A;

	/** repulsive coefficient */
	double f_B;

	/** repulsive exponent */
	double f_p;

	/** attractive exponent */
	double f_q;
	
	/** 3-body secondary energy scale */
	double f_lambda;
	
	/** 3-body exponential parameter */
	double f_gamma;
	
	double f_costheta_ideal;
	/*@}*/
	
	/** \name static parameters 
	 * There parameters are use during evaluation of the static interaction 
	 * functions. These are copied to static when a function pointer is requested */
	/*@{*/
	static double s_eps;
	static double s_sigma;
	static double s_b;
	static double s_A;
	static double s_B;
	static double s_p;
	static double s_q;
	static double s_lambda;
	static double s_gamma;
	static double s_costheta_ideal;
	static bool use_pow;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _STILLINGER_WEBER_T_H_ */
