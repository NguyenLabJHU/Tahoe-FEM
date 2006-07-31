/* $Id: TersoffPropertyT.h,v 1.3 2006-07-31 14:50:30 d-farrell2 Exp $ */
#ifndef _TERSOFF_PROPERTY_T_H_
#define _TERSOFF_PROPERTY_T_H_

/* base class */
#include "ParticlePropertyT.h"
#include "iArrayT.h"
#include "AutoArrayT.h"
#include "ArrayT.h"
#include "nMatrixT.h"
#include "dArray2DT.h"

namespace Tahoe {

/* forward declarations */
class BasicSupportT;

/** defines interface for pair interactions */
class TersoffPropertyT: public ParticlePropertyT
{
public:

	/** pair property factory. Properties that require BasicSupportT cannot be
	 * constructed unless a pointer is provided. */
	static TersoffPropertyT* New(const char* name, const BasicSupportT* support = NULL);

	/** constructor */
	TersoffPropertyT(void);
	
	/** \name Tersoff interactions
	 * Since the interaction functions will be called many times during
	 * a calculation, we do not want these functions to be virtual. Instead,
	 * we define some virtual functions that return pointers to functions
	 * that do the evaluation of the pair interactions. In this way, a single
	 * virtual function call is needed to return a pointer to a static
	 * member function, which is then called many times. 
	 * 
	 * Uses Tersoff potential as defined in Tersoff PRB 1989. 
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
	/*@{*/
	/** definition of function that returns energy
	 * \param rij distance between particles i,j
	 * \param neighbors neighbor list for particle i
	 * \param j index of jth particle
	 * \param type particle type list
	 * \param coords coordinates of particles
	 * \return the pair energy */
	typedef double (*EnergyFunction)(double rij, iArrayT neighbors, const int j, const AutoArrayT<int> type, ArrayT<TersoffPropertyT*> tersoff_properties, nMatrixT<int>& properties_map, const dArray2DT& coords);

	/** definition of function that returns force
	 * \param rij distance between particles i,j
	 * \param neighbors neighbor list for particle i
	 * \param j index of jth particle
	 * \param type particle type list
	 * \param coords coordinates of particles
	 * \return the pair force
	 */
	typedef double (*ForceFunction)(double rij, iArrayT neighbors, const int j, const AutoArrayT<int> type, ArrayT<TersoffPropertyT*> tersoff_properties, nMatrixT<int>& properties_map, const dArray2DT& coords);

	/** definition of function that returns the entire contribution to the stiffness
	 * \param rij distance between particles i,j
	 * \param neighbors neighbor list for particle i
	 * \param j index of jth particle
	 * \param type particle type list
	 * \param coords coordinates of particles
	 * \return the stiffness of the pair interaction
	 */
	typedef double (*StiffnessFunction)(double rij, iArrayT neighbors, const int j, const AutoArrayT<int> type, ArrayT<TersoffPropertyT*> tersoff_properties, nMatrixT<int>& properties_map, const dArray2DT& coords);

	/** return a pointer to the energy function */
	virtual EnergyFunction getEnergyFunction(void) = 0; 

	/** return a pointer to the force function */
	virtual ForceFunction getForceFunction(void) = 0;

	/** return a pointer to the stiffness function */
	virtual StiffnessFunction getStiffnessFunction(void) = 0;
	
	/** return an accessor to obtain cutoff function length parameter 1 */
	virtual double GetR(void) = 0;
	
	/** return an accessor to obtain cutoff function length parameter 2 */
	virtual double GetS(void) = 0;

	/*@}*/
};

} /* namespace Tahoe */

#endif /* _TERSOFF_PROPERTY_T_H_ */
