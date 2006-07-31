/* $Id: TersoffPairT.h,v 1.3 2006-07-31 14:50:30 d-farrell2 Exp $ */
#ifndef _TERSOFF_PAIR_T_H_
#define _TERSOFF_PAIR_T_H_

/* base classes */
#include "TersoffPropertyT.h"
#include "iArrayT.h"
#include "AutoArrayT.h"
#include "ArrayT.h"
#include "nMatrixT.h"
#include "dArray2DT.h"

/* define NULL */
#include <stdlib.h>

namespace Tahoe {

class TersoffPairT: public TersoffPropertyT
{
public:

	/** constructor */
	/* Uses Tersoff potential as defined in Tersoff PRB 1989. 
	 * Basic form of energy for a bond i is the following:
	 * U_i = f_C(r_ij) * [f_R(r_ij) + (b_ij(r_ik,\theta_ijk) * f_A(r_ij))]
	 * 
	 * f_C	: cutoff function
	 * f_R	: repulsive part (2-body term)
	 * b_ij	: bond order parameter (angle dependance is in here)
	 * f_A	: attractive part
	 *
	 * here I use (f_C * f_R) and (f_C * b_ij * f_A) as an environment dependant 2-body term 
	 */
	
	TersoffPairT(void);

	/** \name return interaction functions */
	/*@{*/
	/** return a pointer to the energy function */
	virtual TersoffPropertyT::EnergyFunction getEnergyFunction(void);

	/** return a pointer to the force function */
	virtual TersoffPropertyT::ForceFunction getForceFunction(void);

	/** return a pointer to the stiffness function */
	virtual TersoffPropertyT::StiffnessFunction getStiffnessFunction(void);
	/*@}*/

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	
	/** Acessor to obtain cutoff function length parameter 1 */
	double GetR(void) { return f_R; };
	
	/** Acessor to obtain cutoff function length parameter 2 */
	double GetS(void) { return f_S; };
	
	/*@}*/

private:

	/** \name interaction functions */
	/*@{*/
	static double Energy(double rij, iArrayT neighbors, const int j, const AutoArrayT<int> type, ArrayT<TersoffPropertyT*> tersoff_properties, nMatrixT<int>& properties_map, const dArray2DT& coords);
	static double Force(double rij, iArrayT neighbors, const int j, const AutoArrayT<int> type, ArrayT<TersoffPropertyT*> tersoff_properties, nMatrixT<int>& properties_map, const dArray2DT& coords);
	static double Stiffness(double rij, iArrayT neighbors, const int j, const AutoArrayT<int> typ, ArrayT<TersoffPropertyT*> tersoff_propertiese, nMatrixT<int>& properties_map, const dArray2DT& coords);
	/*@}*/

private:

	/** \name user-defined parameters */
	/*@{*/
	/** Repulsive part energy scaling term. When atoms i & j are different species A = sqrt(A_i * A_j). Part of the 2-body term */
	double f_A;
	
	/** Attractive part energy scaling term. When atoms i & j are different species B = sqrt(B_i * B_j). Part of '3-body' term */
	double f_B;
	
	/** Repulsive part energy exponent term. When atoms i & j are different species lambda = 1/2 (lambda_i + lambda_j).*/
	double f_lambda;
	
	/** Attractive part energy exponent term. When atoms i & j are different species mu = 1/2 (mu_i + mu_j).*/
	double f_mu;
	
	/** Bond order parameter coefficient 1. */
	double f_beta;
	
	/** Bond order parameter exponent term. */
	double f_n;
	
	/** Bond order parameter coefficient 2. */
	double f_c;
	
	/** Bond order parameter coefficient 3. */
	double f_d;
	
	/** Bond order parameter coefficient 4. */
	double f_h;
	
	/** Bond order paramter scaling coefficient. When atoms i & j are the same species, chi = 1.0 */
	double f_chi;
	
	/** Cutoff function length parameter 1. When atoms i & j are different species R = sqrt(R_i * R_j).*/
	double f_R;
	
	/** Cutoff function length parameter 2. When atoms i & j are different species S = sqrt(S_i * S_j).*/
	double f_S;
	/*@}*/
	
	/** \name static parameters 
	 * There parameters are use during evaluation of the static interaction 
	 * functions. These are copied to static when a function pointer is requested */
	/*@{*/
	static double s_A;
	static double s_B;
	static double s_lambda;
	static double s_mu;
	static double s_beta;
	static double s_n;
	static double s_c;
	static double s_d;
	static double s_h;
	static double s_chi;
	static double s_R;
	static double s_S;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _TERSOFF_PAIR_T_H_ */
