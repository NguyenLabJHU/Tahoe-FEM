/* $Id: GradJ2SSKStV1D.h,v 1.6 2004-08-05 23:18:59 paklein Exp $ */
#ifndef _GRAD_J2_SS_KSTV_1D_H_
#define _GRAD_J2_SS_KSTV_1D_H_

/* base classes */
#include "GradSSSolidMatT.h"
#include "IsotropicT.h"
#include "HookeanMatT.h"
#include "ParameterInterfaceT.h"

/* direct members */
#include "dSymMatrixT.h"
#include "dMatrixT.h"
#include "dArrayT.h"
#include "C1FunctionT.h"

namespace Tahoe {

/* forward declarations */
class ElementCardT;
class ifstreamT;
class dSymMatrixT;

class GradJ2SSKStV1D: public GradSSSolidMatT,
		public IsotropicT,
		public HookeanMatT
{
public:

	/** constructor */
	GradJ2SSKStV1D(void);

	/** destructor */
	virtual ~GradJ2SSKStV1D(void);

	/** \name flags */
	/*@{*/
	virtual bool HasHistory(void) const { return true; };
	virtual bool Need_Strain_last(void) const { return true; };
	/*@}*/	

	/** \name update internal variables */
	virtual void UpdateHistory(void);
	
	/** \name reset internal variables to last converged solution */
	virtual void ResetHistory(void);
	
	/** \name spatial description */
	/*@{*/
	/** spatial tangent modulus */
	virtual const dMatrixT& c_ijkl(void);
	
	/** off diagonal moduli for Kar */
	virtual const dMatrixT& odm_bh_ij(void);
	
	/** off diagonal moduli for Kra */
	virtual const dMatrixT& odm_hb_ij(void);
	
	/** modulus for first term in Krr */
	virtual const dMatrixT& gm_hh(void);
	
	/** modulus for second term in Krr */
	virtual const dMatrixT& gm_hp(void);

	/** modulus for third term in Krr */
	virtual const dMatrixT& gm_hq(void);
	
	/** Cauchy stress */
	virtual const dSymMatrixT& s_ij(void);
	
	/** yield criteria moduli */
	virtual double yc(void);
	
	/** returns 1 if the ip has weakened during the iteration, 0 otherwise */
	virtual int weakened(void);
	
	/** return the pressure associated with the last call to 
	 * SolidMaterialT::s_ij. See SolidMaterialT::Pressure
	 * for more information. */
	virtual double Pressure(void) const { return fStress.Trace()/3.0; };
	/*@}*/
	
	/** \name returns the strain energy density for the specified strain */
	virtual double StrainEnergyDensity(void);
	
	/** returns the number of variables computed for nodal extrapolation
	 * during for element output, ie. internal variables */
	virtual int  NumOutputVariables(void) const;
	virtual void OutputLabels(ArrayT<StringT>& labels) const;
	virtual void ComputeOutput(dArrayT& output);
	
	/** \name required parameter flags */
	/*@{*/
	virtual bool NeedLambda(void) const { return true; };
	virtual bool NeedLastLambda(void) const { return true; };
	/*@}*/

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

private:
	
	/* \name indexes to access internal variable (scalars) array */
	enum InternalVariablesT {kYieldCrt       = 0,     /**< yield criteria */
							 kYieldStep      = 1,     /**< boolean for plastic element */
							 kStrain         = 2,     /**< elastic strain */
							 kWeakened       = 3};    /**< elastic strain */

	/** hardening function types */
	enum HardeningFunctionT {kLinear            = 0,
							 kLinearExponential = 1,
							 kCubicSpline       = 2};

	/** incremental change in Lambda */
	virtual double del_Lambda(void);
	virtual double del_GradLambda(void);
	virtual double del_LapLambda(void);
	
	/* \name set modulus */
	virtual void SetModulus(dMatrixT& modulus);

	void AllocateAllElements(void);
	
	void LoadData(const ElementCardT& element, int ip);
	
	const dSymMatrixT& ElasticStrain(const dSymMatrixT& totalstrain,
					 const dSymMatrixT& plasticstrain );
	
	/* hardening functions and their derivatives */
	double      K(double r) const;
	double     dK(double r) const;
	double    ddK(double r) const;
	double   dddK(double r) const;
	double  ddddK(double r) const;

	/** gradients of Isotropic Hardening conjugate force */
	virtual double Grad1R(double fLambda, double fGradLambda, double fLapLambda);
	virtual double Grad2R(double fLambda, double fGradLambda, double fLapLambda);
	virtual double Grad3R(double fLambda, double fGradLambda, double fLapLambda);
	virtual double Grad4R(double fLambda, double fGradLambda, double fLapLambda);

	/** \name yield criteria */
	virtual double YieldCondition(double lambda, double laplambda);

private:
	
	/** \name number of integration points */
	int fNumIP;

	/** \name material input parameters */
	/*@{*/
	double fk_r;           /**< nonlinear isotropic hardening coefficient */
	double fc_r;           /**< length scale for isotropic hardening */
	/*@}*/	
	
	/** C1 isotropic hardening function */
	HardeningFunctionT fType;
	C1FunctionT* fK;	

	/** \name return values */
	/*@{*/
	dSymMatrixT fStress;
	dMatrixT	fModulus;
	dMatrixT    fOffDiagonalModulus_bh, fOffDiagonalModulus_hb;
	dMatrixT    fGradientModulus_hh, fGradientModulus_hp, fGradientModulus_hq;
	/*@}*/
	
	/** \name elastic strain */
	dSymMatrixT	fElasticStrain;
	
	/** \name element level internal variables at current time step*/
	dSymMatrixT fPlasticStrain_j;    /**< plastic strain */
	dSymMatrixT fNorm_j;             /**< unit normal to the stress surface */
	dArrayT     fInternal_j;         /**< internal variables */
	
	/** \name element level internal variables at end of "last" converged time step*/
	dSymMatrixT fPlasticStrain_0;    /**< plastic strain */
	dSymMatrixT fNorm_0;             /**< unit normal to the stress surface */
	dArrayT     fInternal_0;         /**< internal variables */
};
 
/* hardening functions and their 1st derivatives */
inline double GradJ2SSKStV1D::K(double a) const { return fK->Function(a); }
inline double GradJ2SSKStV1D::dK(double a) const { return fK->DFunction(a); }
inline double GradJ2SSKStV1D::ddK(double a) const { return fK->DDFunction(a); }
inline double GradJ2SSKStV1D::dddK(double a) const { return fK->DDDFunction(a); }
inline double GradJ2SSKStV1D::ddddK(double a) const { return fK->DDDDFunction(a); }

} // namespace Tahoe 
#endif /* _GRAD_J2_SS_KSTV_1D_H_ */
