/* $Id: InelasticDuctile2DT.h,v 1.11 2003-11-04 17:35:11 cjkimme Exp $ */
#ifndef _INELASTIC_DUCTILE_2D_T_H_
#define _INELASTIC_DUCTILE_2D_T_H_

/* base class */
#include "SurfacePotentialT.h"
#include "TiedPotentialBaseT.h"

/* direct members */
#include "LAdMatrixT.h"

#include "GlobalT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;

/** Inelastic cohesive zone model for ductile fracture. A cohesive zone model
 * which is in complementary to the kinetic equations for the BCJ model, which
 * are implemented in BCJKineticEqn. */
class InelasticDuctile2DT: public SurfacePotentialT, public TiedPotentialBaseT
{
public:

	/** constructor.
	 * \param time_step reference to the current time step */
	InelasticDuctile2DT(ifstreamT& in, const double& time_step);

	/** return the number of state variables needed by the model */
	int NumStateVariables(void) const;

	/** initialize the state variable array. By default, initialization
	 * involves only setting the array to zero. */
	virtual void InitStateVariables(ArrayT<double>& state);

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

	/** form of stiffness matrix */
	virtual GlobalT::SystemTypeT TangentType(void) const { return GlobalT::kNonSymmetric; }

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

	/** \name implementation of the TiedPotentialBaseT interface */
	/*@{*/
	/** true if nodal release depends on bulk element groups. InelasticDuctile2DT
	 * uses nodal information if bulk groups are specified. */
	virtual bool NeedsNodalInfo(void) const { return iBulkGroups.Length() > 0; };
	
	/** release condition depends on this bulk quantity */
	virtual int NodalQuantityNeeded(void) const;

	/** nodal value is not stress, so does not need to be transformed to local frame */
	virtual bool RotateNodalQuantity(void) const { return false; };
	
	/** true if a nodal release condition is satisfied */
	virtual bool InitiationQ(const nArrayT<double>& sigma) const;

	/** true if the tied potential may ask for nodes to be retied later */
	virtual bool NodesMayRetie(void) const{ return false; };
	
	/** true if node should be retied */
	virtual bool RetieQ(const nArrayT<double>&, const ArrayT<double>&, const dArrayT&) const { return false; }; 

	/** location in state variable array of the state flag */
	virtual int TiedStatusPosition(void) const;
	/*@}*/

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

	/** evaluate the rates */
	void Rates_1(const ArrayT<double>& q, const dArrayT& D, const dArrayT& T,
		dArrayT& dD, dArrayT& dq);

	/** evaluate the Jacobian of the local iteration */
	void Jacobian_1(const ArrayT<double>& q, const dArrayT& D, const dArrayT& T,
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

	/** strength multiplication */
	double fkappa_scale;

	/** true if damage is reversible */
	bool fReversible;
	/*@}*/

	/** \name work space */
	/*@{*/
	dArrayT fdD;
	
	dArrayT fR;
	LAdMatrixT fK;
	/*@}*/
	
	/** \name state variable data */
	/*@{*/
	ArrayT<double> fState;
	
	dArrayT fDelta;
//	dArrayT fTraction;
	dArrayT fdq;
	double& fkappa;
	double& fphi;
	double& fphi_s;
	double& fdissipation;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _INELASTIC_DUCTILE_2D_T_H_ */
