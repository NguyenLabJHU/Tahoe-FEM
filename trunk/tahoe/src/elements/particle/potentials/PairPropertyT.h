/* $Id: PairPropertyT.h,v 1.3 2004-07-15 08:29:49 paklein Exp $ */
#ifndef _PAIR_PROPERTY_T_H_
#define _PAIR_PROPERTY_T_H_

/* base class */
#include "ParticlePropertyT.h"

namespace Tahoe {

/* forward declarations */
class BasicSupportT;

/** defines interface for pair interactions */
class PairPropertyT: public ParticlePropertyT
{
public:

	/** pair property factory. Properties that require BasicSupportT cannot be
	 * constructed unless a pointer is provided. */
	static PairPropertyT* New(const char* name, const BasicSupportT* support = NULL);

	/** constructor */
	PairPropertyT(void);
	
	/** \name pair interactions
	 * Since the interaction functions will be called many times during
	 * a calculation, we do not want these functions to be virtual. Instead,
	 * we define some virtual functions that return pointers to functions
	 * that do the evaluation of the pair interactions. In this way, a single
	 * virtual function call is needed to return a pointer to a static
	 * member function, which is then called many times. */
	/*@{*/
	/** definition of function that returns pair energy
	 * \param r_ab distance between the particles 
	 * \param data_a properties of particle a
	 * \param data_b properties of particle b
	 * \return the pair energy */
	typedef double (*EnergyFunction)(double r_ab, double* data_a, double* data_b);

	/** definition of function that returns pair force
	 * \param r_ab distance between the particles 
	 * \param data_a properties of particle a
	 * \param data_b properties of particle b
	 * \return the pair force */
	typedef double (*ForceFunction)(double r_ab, double* data_a, double* data_b);

	/** definition of function that returns pair stiffness
	 * \param r_ab distance between the particles 
	 * \param data_a properties of particle a
	 * \param data_b properties of particle b
	 * \return the stiffness of the pair interaction */
	typedef double (*StiffnessFunction)(double r_ab, double* data_a, double* data_b);

	/** return a pointer to the energy function */
	virtual EnergyFunction getEnergyFunction(void) = 0; 

	/** return a pointer to the force function */
	virtual ForceFunction getForceFunction(void) = 0;

	/** return a pointer to the stiffness function */
	virtual StiffnessFunction getStiffnessFunction(void) = 0;

	/** return Paradyn-style coefficients table.
	 * returns false if no table is available. */
	virtual bool getParadynTable(const double** coeff, double& dr, 
		int& row_size, int& num_rows) const;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _PAIR_PROPERTY_T_H_ */
