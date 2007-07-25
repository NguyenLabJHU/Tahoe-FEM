/* $Id: TvergHutchIrrev3DT.h,v 1.1 2007-03-12 02:20:08 cjkimme Exp $ */
/* created: paklein (02/05/2000) */
#ifndef _TVERG_HUTCH_IRREV_3D_T_H_
#define _TVERG_HUTCH_IRREV_3D_T_H_

/* base class */
#include "SurfacePotentialT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;

/** Modified cohesive potential from Tvergaard and Hutchinson. The 
 * unmodified model is
 * described in JMPS v41, n6, 1995, 1119-1135. The modification is
 * described in IJNME v44 p1267.  See SurfacePotentialT
 * for more information about the base class.*/
class TvergHutchIrrev3DT: public SurfacePotentialT
{
public:

	/** constructors */
	TvergHutchIrrev3DT(dArrayT& params);
	TvergHutchIrrev3DT(void);	

	/** return the number of state variables needed by the model */
	int NumStateVariables(void) const;

	/** Initialize state variable array */
	virtual void InitStateVariables(ArrayT<double>& state);

	/** Inform CZM of the time step */
	virtual void SetTimeStep(const double& time_step) { fTimeStep = &time_step; };

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

	/* traction potential parameters */
	double fsigma_max; /**< cohesive stress */
	double fd_c_n;     /**< characteristic normal opening to failure */
	double fd_c_t;     /**< characteristic tangential opening to failure */	
	/* non-dimensional opening parameters */
	double fL_1;    /**< non-dimensional opening to initial peak traction */
	double fL_2;    /**< non-dimensional opening to final peak traction */
	double fL_fail; /**< non-dimensional opening to irreversible failure */

	/* penetration stiffness */
	double fpenalty; /**< stiffening multiplier */
	double fK;       /**< penetration stiffness calculated as a function of penalty
	                  * and the initial stiffness of the cohesive potential */


	const double *fTimeStep;

	/** return the secant instead of the tangent stiffness */
	bool fSecantStiffness;
};

} // namespace Tahoe 
#endif /* _TVERG_HUTCH_IRREV_3D_T_H_ */
