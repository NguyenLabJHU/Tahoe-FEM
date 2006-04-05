/*$Id: MR_RP2DT.h,v 1.14 2006-04-05 17:25:58 raregue Exp $*/
/* created by manzari*/
/* Rigid Plastic Cohesive Model for Geomaterials*/
#ifndef _MR_RP_2D_T_H_
#define _MR_RP_2D_T_H_

/* base class */
#include "SurfacePotentialT.h"
#include "TiedPotentialBaseT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;

/** An elastoplastic Traction-Displacement Model for Geometrials.*/

class MR_RP2DT: public SurfacePotentialT, public TiedPotentialBaseT
{
public:

	/** constructor */
	MR_RP2DT(ifstreamT& in);

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

	/** surface traction. Internal variables are integrated over the current
	 * time step. */	
	virtual const dArrayT& Traction(const dArrayT& jump_u, ArrayT<double>& state, const dArrayT& sigma, bool qIntegrate);
    
	/** algorithmic tangent stiffness */
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

	/** \name functions needed for using TiedNodes KBC controller */
	/*@{*/
	/** For MR_RP2DT, returns true to compute nodal stresses. */
	virtual bool NeedsNodalInfo(void) const;
	
	virtual int NodalQuantityNeeded(void) const;

	/** transform stresses to local frame */
	virtual bool RotateNodalQuantity(void) const { return true; };
	
	virtual bool InitiationQ(const nArrayT<double>& sigma) const;
	
	/** Whether or not the nodes will retie */
	virtual bool NodesMayRetie(void) const;

	/** location in state variable array of the state flag */
	virtual int TiedStatusPosition(void) const;
	/*@}*/
	
protected:
	
	/* utility function */
	double signof(double r);
	
	/** dissipated energy */
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

	/* Fracture Energy parameters */
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
	
	int fGroup; /* element group to obtain hydrostatic stress from */
	double fSteps; /* number of steps for k_n to go to 0 after failure */
};

} // namespace Tahoe 
#endif /* _MR_RP_2D_T_H_ */

