/* $Id: InelasticDuctile2DT.h,v 1.1 2003-01-24 18:46:24 paklein Exp $ */
//DEVELOPMENT
#ifndef _INELASTIC_DUCTILE_2D_T_H_
#define _INELASTIC_DUCTILE_2D_T_H_

/* base class */
#include "SurfacePotentialT.h"

/* direct members */
#include "LAdMatrixT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;

/** Inelastic cohesive zone model for ductile fracture */
class InelasticDuctile2DT: public SurfacePotentialT
{
public:

	/** constructor.
	 * \param time_step reference to the current time step */
	InelasticDuctile2DT(ifstreamT& in, const double& time_step);

	/** return the number of state variables needed by the model */
	int NumStateVariables(void) const;

	/** dissipated energy */
	virtual double FractureEnergy(const ArrayT<double>& state);

	/** potential energy */
	virtual double Potential(const dArrayT& jump_u, const ArrayT<double>& state);
	
	/** surface traction. Internal variables are integrated over the current
	 * time step. */	
	virtual const dArrayT& Traction(const dArrayT& jump_u, ArrayT<double>& state, const dArrayT& sigma);

	/** tangent stiffness */
	virtual const dMatrixT& Stiffness(const dArrayT& jump_u, const ArrayT<double>& state, const dArrayT& sigma);

	/** surface status */
	virtual StatusT Status(const dArrayT& jump_u, const ArrayT<double>& state);

	/** write model name to output */
	virtual void PrintName(ostream& out) const;

	/** write model parameters */
	virtual void Print(ostream& out) const;

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

protected:

	/** return true if the potential has compatible (type and sequence)
	 * nodal output - FALSE by default */
	virtual bool CompatibleOutput(const SurfacePotentialT& potential) const;

	/** evaluate the rates */
	void Rates(const ArrayT<double>& q, const dArrayT& D, const dArrayT& T,
		dArrayT& dD, dArrayT& dq);

	/** evaluate the Jacobian of the local iteration */
	void Jacobian(const ArrayT<double>& q, const dArrayT& D, const dArrayT& T,
		const dArrayT& dq, dMatrixT& K);
	
private:

	/** reference to the time step */
	const double& fTimeStep;

	/** \name parameters */
	/*@{*/
	/** initial width of the localized zone */
	double fw_0;
	
	/** rate-independent strain rate limit */
	double feps_0;

	/** critical void volume fraction */
	double fphi_init;

	/** true if damage is reversible */
	bool fReversible;
	/*@}*/

	/** \name work space */
	/*@{*/
	dArrayT fdq;
	dArrayT fdD;
	
	dArrayT fR;
	LAdMatrixT fK;
	/*@}*/
	
	/** \name state variable data */
	/*@{*/
	ArrayT<double> fState;
	
	dArrayT fDelta;
	dArrayT fTraction;
	double& fkappa;
	double& fphi;
	double& fphi_s;
	double& fdissipation;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _INELASTIC_DUCTILE_2D_T_H_ */
