/* $Id: ViscTvergHutch2DT.h,v 1.9 2004-07-15 08:26:02 paklein Exp $ */
#ifndef _VISC_TVERG_HUTCH_2D_T_H_
#define _VISC_TVERG_HUTCH_2D_T_H_

/* base class */
#include "SurfacePotentialT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;

/** cohesive potential from Tvergaard and Hutchinson with viscous
 * dissipation. This model is described in JMPS v41, n6, 1995, 1119-1135. 
 * The model has been augmented with a simple viscous term:
\f[
	T^{visc}_i = - \eta(\lambda) \dot{\Delta_i}
\f]
 * where the viscosity goes as \f$\eta(\lambda) = \eta_0 (1 - \lambda) \f$.
 * \note Aside from the implementation of the viscous effect, this
 * implementation largely duplicates TvergHutch2DT.
 */
class ViscTvergHutch2DT: public SurfacePotentialT
{
public:

	/** constructor */
	ViscTvergHutch2DT(void);

	/** set the source of the time step */
	virtual void SetTimeStep(const double& time_step) { fTimeStep = &time_step; };

	/** return the number of state variables needed by the model.
	 * Need to store the opening displacement from the previous
	 * time increment. The incremental dissipation is also stored */
	int NumStateVariables(void) const;

	/** incremental heat. The amount of energy per unit undeformed area
	 * released as heat over the current increment */
	virtual double IncrementalHeat(const dArrayT& jump, const ArrayT<double>& state);

	/** dissipated energy */
	virtual double FractureEnergy(const ArrayT<double>& state);

	/** potential energy */
	virtual double Potential(const dArrayT& jump_u, const ArrayT<double>& state);
	
	/** surface traction. Internal variables are integrated over the current
	 * time step. */	
	virtual const dArrayT& Traction(const dArrayT& jump_u, ArrayT<double>& state, const dArrayT& sigma, bool qIntegrate);

	/** tangent stiffness */
	virtual const dMatrixT& Stiffness(const dArrayT& jump_u, const ArrayT<double>& state, const dArrayT& sigma);

	/** surface status */
	virtual StatusT Status(const dArrayT& jump_u, const ArrayT<double>& state);

	/** return the number of output variables. returns 0 by default. */
	virtual int NumOutputVariables(void) const;

	/** return labels for the output variables.
	 * \param labels returns with the labels for the output variables. Space is
	 *        allocate by the function. Returns empty by default. */
	virtual void OutputLabels(ArrayT<StringT>& labels) const;

	/** compute the output variables.
	 * \param destination of output values. Allocated by the host code */
	virtual void ComputeOutput(const dArrayT& jump, const ArrayT<double>& state, 
		dArrayT& output);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters  */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** return true if the potential has compatible (type and sequence)
	 * nodal output - FALSE by default */
	virtual bool CompatibleOutput(const SurfacePotentialT& potential) const;
	
private:

	/** index of state variables */
	enum VariableIndexT {
		    kd_n = 1,
		kIncHeat = 4
	};

	/** the time step */
	const double* fTimeStep;

	/* traction potential parameters */
	double fsigma_max; /**< cohesive stress */
	double fd_c_n;     /**< characteristic normal opening to failure */
	double fd_c_t;     /**< characteristic tangential opening to failure */
	
	/* non-dimensional opening parameters */
	double fL_1;    /**< non-dimensional opening to initial peak traction */
	double fL_2;    /**< non-dimensional opening to final peak traction */
	double fL_fail; /**< non-dimensional opening to irreversible failure */

	/** Taylor-Quinney heating factor */
	double fbeta;

	/** damping parameter */
	double feta0;
	
	/* penetration stiffness */
	double fpenalty; /**< stiffening multiplier */
	double fK;       /**< penetration stiffness calculated as a function of penalty
	                  * and the initial stiffness of the cohesive potential */
};

} // namespace Tahoe 
#endif /* _VISC_TVERG_HUTCH_2D_T_H_ */
