/* $Id: GradJ2SSKStV1D.h,v 1.3 2004-07-22 21:10:23 paklein Exp $ */
#ifndef _GRAD_J2_SS_KSTV_1D_H_
#define _GRAD_J2_SS_KSTV_1D_H_

/* base classes */
#include "IsotropicT.h"
#include "HookeanMatT.h"
#include "GradJ2SSC0Hardening1DT.h"

namespace Tahoe {

/** small strain J2 plastic material */
class GradJ2SSKStV1D: public IsotropicT,
						public HookeanMatT,
						public GradJ2SSC0Hardening1DT
{
public:

	/** constructor */
	GradJ2SSKStV1D(void);

	/** \name flags */
	/*@{*/
	virtual bool HasHistory(void) const { return true; };
	virtual bool Need_Strain_last(void) const { return true; };
	/*@}*/	

	/* update internal variables */
	virtual void UpdateHistory(void);

	/* reset internal variables to last converged solution */
	virtual void ResetHistory(void);

	/** \name spatial description */
	/*@{*/
	/** spatial tangent modulus */
	virtual const dMatrixT& c_ijkl(void);

	/** \name spatial description */
	/*@{*/
	/** off diagonal modulus for Kar */
	virtual const dMatrixT& odm_bh_ij(void);

	/** \name spatial description */
	/*@{*/
	/** off diagonal modulus for Kra */
	virtual const dMatrixT& odm_hb_ij(void);

	/** \name spatial description */
	/*@{*/
	/** modulus for first term in Krr */
	virtual const dMatrixT& gm_hh(void);

	/** \name spatial description */
	/*@{*/
	/** modulus for second term in Krr */
	virtual const dMatrixT& gm_hp(void);

	/** \name spatial description */
	/*@{*/
	/** modulus for third term in Krr */
	virtual const dMatrixT& gm_hq(void);

	/** Cauchy stress */
	virtual const dSymMatrixT& s_ij(void);

	/** yield criteria moduli */
	virtual double yc(void);
	
	/** return the pressure associated with the last call to 
	 * SolidMaterialT::s_ij. See SolidMaterialT::Pressure
	 * for more information. */
	virtual double Pressure(void) const { return fStress.Trace()/3.0; };
	/*@}*/

	/* returns the strain energy density for the specified strain */
	virtual double StrainEnergyDensity(void);
	
	/* returns the number of variables computed for nodal extrapolation
	 * during for element output, ie. internal variables */
	virtual int NumOutputVariables(void) const;
	virtual void OutputLabels(ArrayT<StringT>& labels) const;
	virtual void ComputeOutput(dArrayT& output);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** return the description of the given inline subordinate parameter list */
	void DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const;
	
	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/* set modulus */
	virtual void SetModulus(dMatrixT& modulus);

private:

	/* return values */
	dSymMatrixT	fStress;
	dMatrixT	fModulus;

	/* workspaces */
	dSymMatrixT	fStress_3D;

	/** \name return values */
	/*@{*/
	dMatrixT    fOffDiagonalModulus_bh, fOffDiagonalModulus_hb;
	dMatrixT    fGradientModulus_hh, fGradientModulus_hp, fGradientModulus_hq;
	/*@}*/	

	/** \name material input parameters */
	/*@{*/
	double fk_r;           /**< nonlinear isotropic hardening coefficient */
	double fc_r;           /**< length scale for isotropic hardening */
	/*@}*/
};

} // namespace Tahoe 
#endif /* _GRAD_J2_SS_KSTV_1D_H_ */
