/* $Id: InelasticDuctile_RP2DT.h,v 1.8 2006-05-21 17:47:59 paklein Exp $ */
#ifndef _INELASTIC_DUCTILE_RP_2D_T_H_
#define _INELASTIC_DUCTILE_RP_2D_T_H_

/* base class */
#include "SurfacePotentialT.h"
#include "TiedPotentialBaseT.h"

/* direct members */
#include "LAdMatrixT.h"
#include "nVariMatrixT.h"
#include "VariArrayT.h"

#include "GlobalT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;
class ofstreamT;

/** Inelastic cohesive zone model for ductile fracture. A cohesive zone model
 * which is in complementary to the kinetic equations for the BCJ model, which
 * are implemented in BCJKineticEqn. Although the model is derived from 
 * TiedPotentialBaseT, it enforces rigid constraints internally using a
 * penalty formulation. */
class InelasticDuctile_RP2DT: public SurfacePotentialT, public TiedPotentialBaseT
{
public:

	/** constructor.
	 * \param time_step reference to the current time step */
	InelasticDuctile_RP2DT(void);

	/** \name set data pointers */
	/*@{*/
	/** set the pointer to the number of iterations */
	void SetIterationPointer(const int* iteration) { fIteration = iteration; };

	/** set the pointer to the number of iterations */
	void SetTimeStepPointer(const double* time_step) { fTimeStep = time_step; };

	/** set the pointer to the integration point */
	void SetAreaPointer(const double* area) { fArea = area; };
	/*@}*/

	/** \name methods needed for externally enforcing constraints */
	/*@{*/
	/** create an alias to the traction vector within the state variable array */
	static void GetTraction(const ArrayT<double>& state, dArrayT& traction);

	/** create an alias to the flags indicating the evolution equations are activve 
	 * from within the state variable array. Flags should be set to 1.0 to indicate
	 * the equation is active or to 0.0 to indicate that it is inactive. */
	static void GetActiveFlags(const ArrayT<double>& state, dArrayT& active);

	/** assess whether the evolution equations will behave rigidly based on the
	 * given state variables */
	void RigidQ(const ArrayT<double>& state, ArrayT<bool>& rigid) const;

	/** update state variables with the given bulk data */
	void UpdateState(const dArrayT& bulk_nodal_data, ArrayT<double>& state) const;
	/*@}*/

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
	virtual const dArrayT& Traction_penalty(const dArrayT& jump_u, ArrayT<double>& state, const dArrayT& sigma, bool qIntegrate);

	/** tangent stiffness */
	virtual const dMatrixT& Stiffness(const dArrayT& jump_u, const ArrayT<double>& state, const dArrayT& sigma);
	virtual const dMatrixT& Stiffness_penalty(const dArrayT& jump_u, const ArrayT<double>& state, const dArrayT& sigma);

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
	/** true if nodal release depends on bulk element groups. InelasticDuctile_RP2DT
	 * uses nodal information if bulk groups are specified. */
	virtual bool NeedsNodalInfo(void) const { return iBulkGroups.Length() > 0; };
	
	/** release condition depends on this bulk quantity */
	virtual int NodalQuantityNeeded(void) const;

	/** nodal value is not stress, so does not need to be transformed to local frame */
	virtual bool RotateNodalQuantity(void) const { return false; };
	
	/** true if a nodal release condition is satisfied */
	virtual bool InitiationQ(const nArrayT<double>&) const { return true; };

	/** true if the tied potential may ask for nodes to be retied later */
	virtual bool NodesMayRetie(void) const{ return false; };
	
	/** true if node should be retied */
	virtual bool RetieQ(const nArrayT<double>&, const ArrayT<double>&, const dArrayT&) const { return false; }; 

	/** location in state variable array of the state flag */
	virtual int TiedStatusPosition(void) const;
	/*@}*/

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** return true if the potential has compatible (type and sequence)
	 * nodal output - FALSE by default */
	virtual bool CompatibleOutput(const SurfacePotentialT& potential) const;

	/** \name BCJ-like kinetic equations */
	/*@{*/
	/** evaluate the rates */
	void Rates(const ArrayT<double>& q, const dArrayT& D, const dArrayT& T,
		dArrayT& dD, dArrayT& dq, dArrayT& active);

	/** evaluate the Jacobian of the local iteration */
	void Jacobian(const ArrayT<double>& q, const dArrayT& D, const dArrayT& T,
		const dArrayT& dq, dMatrixT& K);
	/*@}*/

private:

	/** reduce local the active equations */
	void Assemble(const dArrayT& R, const dMatrixT& K, const dArrayT& active);

private:

	/** pointer to the area of the current evaluation point */
	const double* fArea;

	/** pointer to the area of the current time step */
	const double* fTimeStep;

	/** pointer to the iteration number counter */
	const int* fIteration;

	/** \name parameters */
	/*@{*/
	/** initial width of the localized zone */
	double fw_0;

	/** critical void volume fraction */
	double fphi_init;

	/** void volume fraction at failure */
	double fphi_max;

	/** damage exponent */
	double fdamage_exp;
	
	/** absolute tolerance of local iteration */
	double fabs_tol;

	/** relative tolerance of local iteration */
	double frel_tol;

	/** true if damage is reversible */
	bool fReversible;

	/** penalty stiffness for enforcing rigid behavior */
	double fConstraintStiffness;
	
	/** number of iterations between calculating an updated traction vector. For
	 * InelasticDuctile_RP2DT:: fUpdateIterations > 1, the tractions are assumed
	 * to behave like fixed external forces and therefore do not return any stiffness */
	int fUpdateIterations;
	/*@}*/

	/** \name BCJ model kinetic parameters */
	/*@{*/
	double fTemperature;
	double fC1, fC2, fC3, fC4, fC5, fC6;
	double fC19, fC20, fC21;
	double fphi_0;
	/*@}*/

	/** \name BCJ model derived parameters */
	/*@{*/
	/** rate-independent strain rate limit */
	double feps_0;
	double fV, fY;
	/*@}*/

	/** \name work space */
	/*@{*/
	int fFixedIterationCount;
	dArrayT fdD;
	
	dArrayT fR, fR_temp;
	LAdMatrixT fK_temp;
	dMatrixT fK;

	VariArrayT<double> fR_man;
	nVariMatrixT<double> fK_man;
	/*@}*/
	
	/** \name state variable data */
	/*@{*/
	ArrayT<double> fState;
	
	dArrayT fDelta;
	dArrayT fdq;

	/** 1 or 0 depending on where that equation is active */
	dArrayT feq_active;
	
	double& fkappa;
	double& fphi;
	double& fphi_s;
	double& fdissipation;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _INELASTIC_DUCTILE_RP_2D_T_H_ */
