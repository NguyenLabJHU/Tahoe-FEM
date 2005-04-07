/* $Id: EAMPropertyT.h,v 1.2 2004-04-09 02:03:05 hspark Exp $ */
#ifndef _EAM_PROPERTY_T_H_
#define _EAM_PROPERTY_T_H_

/* base class */
#include "ParticlePropertyT.h"

namespace Tahoe {

/** defines interface for embedded atom method (EAM) interactions */
class EAMPropertyT: public ParticlePropertyT
{
public:

	/** constructor */
	EAMPropertyT(void);
	
	/** \name EAM 
	 * Since the interaction functions will be called many times during
	 * a calculation, we do not want these functions to be virtual. Instead,
	 * we define some virtual functions that return pointers to functions
	 * that do the evaluation of the pair interactions. In this way, a single
	 * virtual function call is needed to return a pointer to a static
	 * member function, which is then called many times. */
	/*@{*/
	typedef double (*PairEnergyFunction)(double r_ab, double* data_a, double* data_b);
	typedef double (*EmbedEnergyFunction)(double rho_ab, double* data_a, double* data_b);
	typedef double (*EDEnergyFunction)(double r_ab, double* data_a, double* data_b);

	typedef double (*PairForceFunction)(double r_ab, double* data_a, double* data_b);
	typedef double (*EmbedForceFunction)(double rho_ab, double* data_a, double* data_b);
	typedef double (*EDForceFunction)(double r_ab, double* data_a, double* data_b);

	typedef double (*PairStiffnessFunction)(double r_ab, double* data_a, double* data_b);
	typedef double (*EmbedStiffnessFunction)(double rho_ab, double* data_a, double* data_b);
	typedef double (*EDStiffnessFunction)(double r_ab, double* data_a, double* data_b);
	/*@}*/

	/** \name return interaction functions */
	/*@{*/
	/** return a pointer to the energy function */
	virtual PairEnergyFunction getPairEnergy(void) = 0;
	virtual EmbedEnergyFunction getEmbedEnergy(void) = 0;
	virtual EDEnergyFunction getElecDensEnergy(void) = 0;

	/** return a pointer to the force function */
	virtual PairForceFunction getPairForce(void) = 0;
	virtual EmbedForceFunction getEmbedForce(void) = 0;
	virtual EDForceFunction getElecDensForce(void) = 0;

	/** return a pointer to the stiffness function */
	virtual PairStiffnessFunction getPairStiffness(void) = 0;
	virtual EmbedStiffnessFunction getEmbedStiffness(void) = 0;
	virtual EDStiffnessFunction getElecDensStiffness(void) = 0;
	/*@}*/

	/** return Paradyn-style coefficients table.
	 * returns false if no table is available. */
	virtual bool getParadynTable(const double** coeff, double& dr, 
		int& row_size, int& num_rows) const = 0;
		
	virtual double GetLatticeParameter(void) const = 0;
};

} /* namespace Tahoe */

#endif /* _EAM_PROPERTY_T_H_ */
