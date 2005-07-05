/* $Id: MFGPMaterialT.h,v 1.3 2005-07-05 07:14:32 paklein Exp $  */
#ifndef _MFGP_MATERIAL_T_H_
#define _MFGP_MATERIAL_T_H_

#include "Environment.h"
#include "GlobalT.h"
#include "ios_fwd_decl.h"

/* base class */
#include "ParameterInterfaceT.h"

/* direct members */
#include "MFGPMatSupportT.h"
#include "dMatrixT.h"

namespace Tahoe {

/* forward declarations */
class ElementCardT;
class dArrayT;
template <class TYPE> class ArrayT;
class StringT;

class ifstreamT;
class ElementBaseT;
class dMatrixT;
class dSymMatrixT;
class LocalArrayT;
class MFGP_AssemblyT;

/** interface for continuum materials. */
class MFGPMaterialT: virtual public ParameterInterfaceT
{
public:

	/** \name 2D constrain options */
	enum ConstraintT {
		kNoConstraint = 0, /**< no constraint, material is 3D */
		kPlaneStress = 1, /**< plane stress */
		kPlaneStrain = 2  /**< plane strain */};
	ConstraintT static int2ConstraintT(int i);
		
	/** constructor */
	MFGPMaterialT(void);

	/** destructor */
	virtual ~MFGPMaterialT(void);

	/** set the material support or pass NULL to clear */
	virtual void SetMFGPMatSupport(const MFGPMatSupportT* support);

	/** form of tangent matrix. \return symmetric by default */
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/** reference to the material support */
	const MFGPMatSupportT& MFGPMatSupport(void) const;

	/** reference to the host element */
	const MFGP_AssemblyT& MFGPAssembly(void) const;

	/** number of element nodes in the host element group */
	int NumElementNodes() const;

	/** number of degrees of freedom (per node) in the host
	 * element group. */
	int NumDOF(void) const;

	/** number of spatial dimensions in the host element group. */
	int NumSD(void) const;

	/** the total number of integration points per element in the
	 * host element group. */
	int NumIP(void) const;

	/** the current integration point within the element of evaluation. */
	int CurrIP(void) const;

	/** return the total number of elements in the host element
	 * group. */
	int NumElements(void) const;

	/** return the number of the current element of evaluation. */
	int CurrElementNumber(void) const;

	/** reference to the ElementCardT for the  specified element. */
	ElementCardT& ElementCard(int card) const;

	/** reference to the ElementCardT for the current element of
	 * evaluation */
	ElementCardT& CurrentElement(void) const;

	/** apply pre-conditions at the current time step. Called once for
	 * the model at the beginning of a time increment */
	virtual void InitStep(void);

	/** finalize the current time step. Called once for the model at 
	 * the end of a time increment */
	virtual void CloseStep(void);

	/** \name history variables */
	/*@{*/
	/** return true if the material has history variables.
	 * \return false by default. */
	virtual bool HasHistory(void) const { return false; };
	
	/** return true if model needs ContinuumMaterialT::PointInitialize
	 * to be called for every integration point of every element as
	 * part of the model initialization. \return false by default. */
	virtual bool NeedsPointInitialization(void) const;
	
	/** model initialization. Called per integration point for every
	 * element using the model. Deformation variables are available
	 * during this call. */
	virtual void PointInitialize(void);

	/** update internal variables. Called once per element for all
	 * elements using the model, hence no deformation variables are
	 * available during this call. */
	virtual void UpdateHistory(void);

	/** restore internal variables to their state at the beginning of
	 * the current time increment. Called once per element for all
	 * elements using the model, hence no deformation variables are
	 * available during this call. */
	virtual void ResetHistory(void);
	/*@}*/

	/** \name material output variables */
	/*@{*/
	/** return the number of constitutive model output parameters
	 * per evaluation point. Used by the host element group in
	 * conjunction with ContinuumMaterialT::OutputLabels and
	 * ContinuumMaterialT::ComputeOutput to collect model variables
	 * for output. \return zero by default */
	virtual int NumOutputVariables(void) const;

	/** return the labels for model output parameters
	 * per evaluation point. Used by the host element group in
	 * conjunction with ContinuumMaterialT::NumOutputVariables and
	 * ContinuumMaterialT::ComputeOutput to collect model variables
	 * for output.
	 * \param labels destination for the variable labels. Returns
	 *        with length ContinuumMaterialT::NumOutputVariables */
	virtual void OutputLabels(Tahoe::ArrayT<StringT>& labels) const;

	/** return material output variables. Used by the host element group 
	 * in conjunction with ContinuumMaterialT::NumOutputVariables and
	 * ContinuumMaterialT::OutputLabels to collect model variables
	 * for output. Called per integration point. Deformation variables
	 * are available.
	 * \param output destination for the output. Must be passed in with
	 *        length ContinuumMaterialT::NumOutputVariables */
	virtual void ComputeOutput(Tahoe::dArrayT& output);

	/** returns true if two materials have compatible output variables.
	 * Used by the host element to determine whether the two material
	 * models can be used within the same host element group when
	 * requesting model-specific, materials output. */
	static bool CompatibleOutput(const MFGPMaterialT& m1, const MFGPMaterialT& m2);
	/*@}*/
	
	/** \name spatial description */
	/*@{*/
	/** spatial tangent modulus */
	virtual const dMatrixT& c_ijkl(void) = 0;

	/** Cauchy stress */
	virtual const dSymMatrixT& s_ij(void) = 0;

	/** return the pressure associated with the last call to 
	 * SolidMaterialT::s_ij. The value is not guaranteed to
	 * persist during intervening calls to any other non-const
	 * accessor. \return 1/3 of the trace of the three-dimensional
	 * stress tensor, regardless of the dimensionality of the
	 * problem. */
	virtual double Pressure(void) const = 0;
	/*@}*/

	/** \name material description */
	/*@{*/
	/** material tangent moduli */
	virtual const dMatrixT& C_IJKL(void) = 0;

	/** 2nd Piola-Kirchhoff stress */
	virtual const dSymMatrixT& S_IJ(void) = 0;
	/*@}*/
	
	/** pass to the material model */
	virtual void Initialize(ElementCardT element, int ip, int n_ip, dSymMatrixT strain_ip, dSymMatrixT strain_lap_ip, 
					dArrayT lambdaPM_ip, dArrayT lambdaPM_lap_ip) = 0;
	
	/** yield function */
	virtual const double& YieldF(void) = 0;
	
	virtual const dMatrixT& c_UU1_ijkl(void) = 0; 
	virtual const dMatrixT& c_UU2_ijkl(void) = 0;
	virtual const dMatrixT& c_ULam1_ij(void) = 0; 
	virtual const dMatrixT& c_ULam2_ij(void) = 0;
	virtual const dMatrixT& c_LamU1_ij(void) = 0; 
	virtual const dMatrixT& c_LamU2_ij(void) = 0;
	virtual const dMatrixT& c_LamLam1(void) = 0; 
	virtual const dMatrixT& c_LamLam2(void) = 0;

	/** 2D constrain options or kNoConstraint::kNoConstraint if the material
	 * is not 2D */
	ConstraintT Constraint(void) const { return fConstraint; };
	
	virtual bool NeedDisp(void) const     { return false; };
	virtual bool NeedLastDisp(void) const { return false; };
	virtual bool NeedVel(void) const      { return false; };
	
	/** return the strain in the material. The definition of strain will be
	 * dependent on the subclass */
	virtual void Strain(dSymMatrixT& strain) = 0;
	
	/** \return mass density */
	virtual double Density(void);
	
	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list. */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/
	
protected:

	/** support from the host code */
	const MFGPMatSupportT* fMFGPMatSupport;
	
	/** number of degrees of freedom */
	int fNumDOF;

	/** number of degrees of spatial dimensions */
	int fNumSD;
	
	/** number of integration points */
	int fNumIP;
	
	/* mass density */
	double fDensity;

	/** 2D constrain option */
	ConstraintT fConstraint;
};

/* inlines */
inline int MFGPMaterialT::NumDOF(void) const { return fNumDOF; }
inline int MFGPMaterialT::NumSD(void) const { return fNumSD; }
inline int MFGPMaterialT::NumIP(void) const { return fNumIP; }

inline const MFGPMatSupportT& MFGPMaterialT::MFGPMatSupport(void) const
{ 
#if __option(extended_errorcheck)
	if (!fMFGPMatSupport)
		ExceptionT::GeneralFail("MFGPMaterialT::MFGPMatSupport", "material support not set");
#endif
	return *fMFGPMatSupport; 
}

inline int MFGPMaterialT::CurrIP(void) const { return MFGPMatSupport().CurrIP(); };

inline const MFGP_AssemblyT& MFGPMaterialT::MFGPAssembly(void) const { 
	return *(MFGPMatSupport().MFGPAssembly()); 
}

/* returns the density */
inline double MFGPMaterialT::Density(void) { return fDensity; }

} // namespace Tahoe 

#endif /* _MFGP_MATERIAL_T_H_ */
