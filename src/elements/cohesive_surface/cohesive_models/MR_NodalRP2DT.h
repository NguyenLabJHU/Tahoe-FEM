/* $Id: MR_NodalRP2DT.h,v 1.1 2006-10-08 19:22:19 regueiro Exp $ */
#ifndef _MR_NODAL_RP_2D_T_H_
#define _MR_NODAL_RP_2D_T_H_

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

/** Rigid plastic geomaterial cohesive zone model. Although the model is derived from 
 * TiedPotentialBaseT, it enforces rigid constraints internally using a
 * penalty formulation. */
class MR_NodalRP2DT: public SurfacePotentialT, public TiedPotentialBaseT
{
public:

	/** constructor.
	 * \param time_step reference to the current time step */
	MR_NodalRP2DT(void);

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
	/*@}*/

	/** return the number of state variables needed by the model */
	int NumStateVariables(void) const;
	
	/** initialize the state variable array. By default, initialization
	 * involves only setting the array to zero. */
	virtual void InitStateVariables(ArrayT<double>& state);

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
	/** true if nodal release depends on bulk element groups. MR_NodalRP2DT
	 * uses nodal information if bulk groups are specified. */
	virtual bool NeedsNodalInfo(void) const { return true; };
	
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

	/* utility function */
	double signof(double r);
	
	/** yield  */
	virtual double YFValue(const ArrayT<double>& state);
	
	/** yield function */
	double& Yield_f(const dArrayT& Sig, const dArrayT& qn, double& ff);
    
	/** hardening functions*/
	dArrayT& qbar_f(const dArrayT& Sig, const dArrayT& qn, dArrayT& qbar);
    
	/** normal to yield surface */
	dArrayT& dfdSig_f(const dArrayT& Sig, const dArrayT& qn, dArrayT& dfdSig);
    
	/** normal to plastic potential */
	dArrayT& dQdSig_f(const dArrayT& Sig, const dArrayT& qn, dArrayT& dQdSig);
    
	/** derivatives for local iteration and calculation of algorithmic tangent operator */
	dArrayT& dfdq_f(const dArrayT& Sig, const dArrayT& qn, dArrayT& dfdq);    
	dMatrixT& dQdSig2_f(const dArrayT& qn, dMatrixT& dQdSig2);
	dMatrixT& dQdSigdq_f(const dArrayT& Sig, const dArrayT& qn, dMatrixT& dQdSigdq);
	dMatrixT& dqbardSig_f(const dArrayT& Sig, const dArrayT& qn, dMatrixT& dqbardSig);
	dMatrixT& dqbardq_f(const dArrayT& Sig, const dArrayT& qn, dMatrixT& dqbardq);

private:

	/** pointer to the area of the current evaluation point */
	const double* fArea;

	/** pointer to the area of the current time step */
	const double* fTimeStep;

	/** pointer to the iteration number counter */
	const int* fIteration;

	/** \name parameters */
	/*@{*/
	double fGf_I;    /* Mode_I Fracture Energy */
	double fGf_II;   /* Mode_II Fracture Energy */

	/* Inelastic response parameters */
	double fchi_p; /* peak tensile strength*/  
	double fchi_r; /* residual tensile strength */
	double fc_p;   /* peak cohesion */
	double fc_r;   /* residual cohesion */
	double fphi_p; /* peak friction angle */
	double fphi_r; /* critical state friction angle */
	double fpsi_p; /* peak dilation angle */
	double falpha_chi; /* Coefficient of chi degredation */
	double falpha_c; /* Coefficient of c degredation */
	double falpha_phi; /*  Coefficient of phi degredation */
	double falpha_psi; /*  Coefficient of psi degredatione */
	double fTol_1;    /*  Tolerance for Yield Function */
	double fTol_2; /*  Tolerance for Residuals */
	double fchi, fc, fphi, fpsi;

	/** penalty stiffness for enforcing rigid behavior */
	double fConstraintStiffness;
	
	/** number of iterations between calculating an updated traction vector. For
	 * MR_NodalRP2DT:: fUpdateIterations > 1, the tractions are assumed
	 * to behave like fixed external forces and therefore do not return any stiffness */
	int fUpdateIterations;
	/*@}*/

	/** \name work space */
	/*@{*/
	int fFixedIterationCount;
	/*@}*/
	
	/** \name state variable data */
	/*@{*/
	ArrayT<double> fState;
	
	/** 1 or 0 depending on where that equation is active */
	dArrayT feq_active;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _MR_NODAL_RP_2D_T_H_ */
