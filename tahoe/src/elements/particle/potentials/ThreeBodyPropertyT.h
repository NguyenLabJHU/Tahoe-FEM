/* $Id: ThreeBodyPropertyT.h,v 1.2 2004-11-23 01:55:38 cjkimme Exp $ */
#ifndef _THREE_BODY_PROPERTY_T_H_
#define _THREE_BODY_PROPERTY_T_H_

/* base class */
#include "PairPropertyT.h"

namespace Tahoe {

/* forward declarations */
class BasicSupportT;
class dMatrixT;

/** defines interface for pair interactions */
class ThreeBodyPropertyT: public PairPropertyT
{
public:

	/** pair property factory. Properties that require BasicSupportT cannot be
	 * constructed unless a pointer is provided. */
	static ThreeBodyPropertyT* New(const char* name, const BasicSupportT* support = NULL);

	/** constructor */
	ThreeBodyPropertyT(void);
	
	/** \name Three-body interactions
	 * Since the interaction functions will be called many times during
	 * a calculation, we do not want these functions to be virtual. Instead,
	 * we define some virtual functions that return pointers to functions
	 * that do the evaluation of the pair interactions. In this way, a single
	 * virtual function call is needed to return a pointer to a static
	 * member function, which is then called many times. 
	 * 
	 * Three-body interactions are currently only impelemented for energies
	 * of the form phi_ijk(r_ij,r_ik,theta_ijk). Functions with arguments
	 * of ri, rj, and rk must are thusly ordered. 
	 */
	/*@{*/
	/** definition of function that returns three-body energy
	 * \param ri 
	 * \param rj
	 * \param rk
	 * \return the pair energy */
	typedef double (*EnergyFunction)(double* ri, double* rj, double* rk);

	/** definition of function that returns three-body force
	 * \param ri -- coords of outer loop particle
	 * \param rj -- coords of second loop particle
	 * \param rk -- coords of innermost loop particle
	 * \param fij -- ij force i.e. d phi/ d r_ij -- return value
	 * \param fik -- ik force i.e. d phi/ d r_ik -- return value
	 * \return NULL if beyond cutoff, f's otherwise
	 */
	typedef double* (*ForceFunction)(const double* ri, const double* rj, const double* rk, double* fij, double* fik);

	/** definition of function that returns the entire 3-body contribution to the stiffness
	 * \param ri position vector for atom i
	 * \param rj position vector for atom j
	 * \param rk position vector for atom k
	 * \param K_ijk symmetric stiffness d^2 phi/d r_l d r_m l,m in {i,j,k}
	 * \return NULL if beyond cutoff, K otherwise
	 */
	typedef double* (*StiffnessFunction)(const double* ri, const double* rj, const double* rk, 
										dMatrixT& K_ijk);

	/** return a pointer to the energy function */
	virtual EnergyFunction getThreeBodyEnergyFunction(void) = 0; 

	/** return a pointer to the force function */
	virtual ForceFunction getThreeBodyForceFunction(void) = 0;

	/** return a pointer to the stiffness function */
	virtual StiffnessFunction getThreeBodyStiffnessFunction(void) = 0;

	/*@}*/
};

} /* namespace Tahoe */

#endif /* _THREE_BODY_PROPERTY_T_H_ */
