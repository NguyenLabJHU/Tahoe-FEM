/*$Id: MR3DT.h,v 1.2 2007-03-27 21:46:31 skyu Exp $*/
/* created by manzari*/
/* Elastoplastic Cohesive Model for Geomaterials*/
#ifndef _MR_3D_T_H_
#define _MR_3D_T_H_

/* base class */
#include "SurfacePotentialT.h"

/* direct members */
#include "ifstreamT.h"
#include "ofstreamT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;
class ofstreamT;

/** An elastoplastic Traction-Displacement Model for Geometrials.*/

class MR3DT: public SurfacePotentialT
{
public:

	/** constructor */
	MR3DT(void);

	/** return the number of state variables needed by the model.
	 * Need to store the opening displacement from the previous
	 * time increment. */
	int NumStateVariables(void) const;
	
	virtual void InitStateVariables(ArrayT<double>& state);	
	
	/** dissipated energy. Total amount of energy dissipated reaching
	 * the current state. */
	virtual double FractureEnergy(const ArrayT<double>& state);
	
	/** potential energy */
	virtual double Potential(const dArrayT& jump, const ArrayT<double>& state);
	
	/** surface status */
	virtual StatusT Status(const dArrayT& jump, const ArrayT<double>& state);
	
	/** dissipated energy */
	virtual double YFValue(const ArrayT<double>& state);

	/** surface traction. Internal variables are integrated over the current
	 * time step. */	

	virtual const dArrayT& Traction(const dArrayT& jump_u, ArrayT<double>& state, const dArrayT& sigma, bool qIntegrate);

	double& Yield_f(const dArrayT& Sig, const dArrayT& qn, double& ff);
	dArrayT& qbar_f(const dArrayT& Sig, const dArrayT& qn, dArrayT& qbar);
	dArrayT& dfdSig_f(const dArrayT& Sig, const dArrayT& qn, dArrayT& dfdSig);
	dArrayT& dQdSig_f(const dArrayT& Sig, const dArrayT& qn, dArrayT& dQdSig);
	dArrayT& dfdq_f(const dArrayT& Sig, const dArrayT& qn, dArrayT& dfdq);
	dMatrixT& dQdSig2_f(const dArrayT& Sig, const dArrayT& qn, dMatrixT& dQdSig2);
	dMatrixT& dQdSigdq_f(const dArrayT& Sig, const dArrayT& qn, dMatrixT& dQdSigdq);
	dMatrixT& dqbardSig_f(const dArrayT& Sig, const dArrayT& qn, dMatrixT& dqbardSig);
	dMatrixT& dqbardq_f(const dArrayT& Sig, const dArrayT& qn, dMatrixT& dqbardq);

	/** tangent stiffness */
	virtual const dMatrixT& Stiffness(const dArrayT& jump_u, const ArrayT<double>& state, const dArrayT& sigma);

	/** form of stiffness matrix */
	virtual GlobalT::SystemTypeT TangentType(void) const { return GlobalT::kNonSymmetric; }

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

	/** For MR3DT, returns true to compute nodal tractions. */
	virtual bool NeedsNodalInfo(void);
	
	double signof(double r);
	
	virtual int NodalQuantityNeeded(void);
	//	virtual double ComputeNodalValue(const dArrayT &);
	//	virtual void UpdateStateVariables(const dArrayT &, ArrayT<double> &);
	virtual void SetElementGroupsNeeded(iArrayT& iGroups);
	
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
	
	
private:

	//bool initiationQ(const ArrayT<double>&);

	/* traction rate parameters */
	double fE_t1;		/* elastic tangential stiffness in 1 direction*/
	double fE_t2;		/* elastic tangential stiffness in 2 direction*/
	double fE_n;		/* elastic normal stiffness */
	double fGf_I;		/* Mode_I Fracture Energy */
	double fGf_II;		/* Mode_II Fracture Energy */

	/* Inelastic response parameters */
	double fchi_p;		/* peak tensile strength*/
	double fchi_r;		/* residual tensile strength */
	double fc_p;		/* peak cohesion */
	double fc_r;		/* residual cohesion */
	double fphi_p;		/* peak friction angle */
	double fphi_r;		/* critical state friction angle */
	double fpsi_p;		/* peak dilation angle */
	double falpha_chi;	/* Coefficient of chi degredation */
	double falpha_c;	/* Coefficient of c degredation */
	double falpha_phi;	/*  Coefficient of phi degredation */
	double falpha_psi;	/*  Coefficient of psi degredatione */
	double fTol_1;		/*  Tolerance for Yield Function */
	double fTol_2;		/*  Tolerance for Residuals */
	double fchi, fc, fphi, fpsi;
	
	int fGroup;		/* element group to obtain hydrostatic stress from */
	double fSteps;		/* number of steps for k_n to go to 0 after failure */


	/** write output for debugging */
	/*@{*/
	/** output file stream */
	ofstreamT mr_ep_3d_out;

	/** line output formating variables */
	int outputPrecision, outputFileWidth;
	/*@}*/
};

} // namespace Tahoe 
#endif /* _MR_3D_T_H_ */

