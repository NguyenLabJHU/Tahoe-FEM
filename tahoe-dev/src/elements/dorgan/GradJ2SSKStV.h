/* $Id: GradJ2SSKStV.h,v 1.1 2004-09-02 18:25:04 rdorgan Exp $ */
#ifndef _GRAD_J2_SS_KSTV_H_
#define _GRAD_J2_SS_KSTV_H_

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

class GradJ2SSKStV: public GradSSSolidMatT,
		public IsotropicT,
		public HookeanMatT
{
public:

	/** constructor */
	GradJ2SSKStV(void);

	/** destructor */
	virtual ~GradJ2SSKStV(void);

	/** \name flags */
	/*@{*/
	virtual bool HasHistory(void) const { return true; };
	virtual bool Need_Strain_last(void) const { return true; };
	/*@}*/	

	/** \name update internal variables */
	virtual void UpdateHistory(void);
	
	/** \name reset internal variables to last converged solution */
	virtual void ResetHistory(void);
	
	/** \name update flag describing a weakened ip */
	virtual	void UpdateWeakened(const ElementCardT& element, int ip);
	
	/** \name reset flag describing a weakened ip */
	virtual	void ResetWeakened(const ElementCardT& element, int ip);
	
	/** \name spatial description */
	/*@{*/
	virtual const dMatrixT& c_ijkl(void);    /**< spatial tangent modulus */
	virtual const dMatrixT& odm_bh_ij(void); /**< off diagonal moduli for Kar */
	virtual const dMatrixT& odm_hb_ij(void); /**< off diagonal moduli for Kra */
	virtual const dMatrixT& gm_hh(void);     /**< modulus for first term in Krr */
	virtual const dMatrixT& gm_hp(void);     /**< modulus for second term in Krr */
	virtual const dMatrixT& gm_hq(void);     /**< modulus for third term in Krr */
	virtual const dSymMatrixT& s_ij(void);   /**< Cauchy stress */
	virtual const dSymMatrixT& n_ij(void);   /**< unit norm */
	virtual double yc(void);                 /**< yield criteria moduli */
	virtual double ys(void);                 /**< yield strenth */
	virtual int weakened(void);              /**< returns 1 if the ip has weakened during the iteration, 0 otherwise */
	/*@}*/
	
	/** return the pressure associated with the last call to 
	 * SolidMaterialT::s_ij. See SolidMaterialT::Pressure
	 * for more information. */
	/*@{*/
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

	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** return the description of the given inline subordinate parameter list */
	void DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const;
	
	/** a pointer to the ParameterInyterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/* returns elastic strain */
	virtual const dSymMatrixT& ElasticStrain(const dSymMatrixT& totalstrain,
		const ElementCardT& element, int ip);

private:
	
	/* \name indexes to access internal variable (scalars) array */
	enum InternalVariablesT {kYieldCrt       = 0,     /**< yield criteria */
							 kYieldStep      = 1,     /**< boolean for plastic element */
							 kStrain         = 2,     /**< total strain */
							 kWeakened       = 3,     /**< flag to indicate if ip is weakened*/
							 kIsoHard        = 4,     /**< isotropic hardening */
							 kLapIsoHard     = 5,     /**< laplacian of isotropic hardening */
							 kYieldStrength  = 6};    /**< yield strength = syp - R */

	/* \name hardening function types */
	enum HardeningFunctionT {kLinear            = 0,
							 kLinearExponential = 1,
							 kCubicSpline       = 2};

	/* \name* incremental change in Lambda */
	virtual double   del_Lambda(void);
	virtual dMatrixT del_GradLambda(void);
	virtual double   del_LapLambda(void);
	
	/* \name set modulus */
	virtual void SetModulus(dMatrixT& modulus);

	void AllocateAllElements(void);
	
	void LoadData(const ElementCardT& element, int ip);
	
	/* \name hardening functions and their derivatives */
	double      K(double r) const;
	double     dK(double r) const;
	double    ddK(double r) const;
	double   dddK(double r) const;
	double  ddddK(double r) const;

	/* \name gradients of Isotropic Hardening conjugate force */
	virtual dMatrixT Grad1R(double fLambda, dMatrixT fGradLambda, double fLapLambda);
	virtual double   Grad2R(double fLambda, dMatrixT fGradLambda, double fLapLambda);
	virtual dMatrixT Grad3R(double fLambda, dMatrixT fGradLambda, double fLapLambda);
	virtual double   Grad4R(double fLambda, dMatrixT fGradLambda, double fLapLambda);

	/** \name compute the consistent elastic tangent moduli */
	void TangentModulus();

	/** \name yield criteria */
	virtual double YieldCondition(double isohard, dMatrixT gradisohard, double lapisohard);

private:
	
	/** \name number of integration points */
	int fNumIP;

	/** \name material input parameters */
	/*@{*/
	double fk_r;           /**< nonlinear isotropic hardening coefficient */
	double fc_r;           /**< length scale for isotropic hardening */
	/*@}*/	
	
	/* \name C1 isotropic hardening function */
	HardeningFunctionT fType;
	C1FunctionT* fK;	

	/** \name return values */
	/*@{*/
	dSymMatrixT fStress;
	dMatrixT    fOffDiagonalModulus_bh, fOffDiagonalModulus_hb;
	dMatrixT    fGradientModulus_hh, fGradientModulus_hp, fGradientModulus_hq;
	/*@}*/
	
	/** \name elastic strain */
	dSymMatrixT	fElasticStrain;
	
	/** \name element level internal variables at current time step*/
	dSymMatrixT fUnitNorm;           /**< unit normal to the stress surface */
	dSymMatrixT fPlasticStrain_j;    /**< plastic strain */
	dArrayT     fInternal_j;         /**< internal variables */
	dMatrixT    fGradIsoHard_j;      /**< gradient of isotropic hardening */
	dMatrixT	fEModulus;            /**< pseudo-elastic stiffness operator */
	
	/** \name element level internal variables at end of "last" converged time step*/
	dSymMatrixT fPlasticStrain_0;    /**< plastic strain */
	dArrayT     fInternal_0;         /**< internal variables */
	dMatrixT    fGradIsoHard_0;      /**< gradient of isotropic hardening */

	/* work space */
	dSymMatrixT fRelStress;
	dMatrixT fTensorTemp1;	
	dMatrixT fTensorTemp2;	
	dMatrixT fTensorTemp3;	
	dSymMatrixT fTensorTemp4;	
};
 
/* hardening functions and their 1st derivatives */
inline double GradJ2SSKStV::K(double a) const { return fK->Function(a); }
inline double GradJ2SSKStV::dK(double a) const { return fK->DFunction(a); }
inline double GradJ2SSKStV::ddK(double a) const { return fK->DDFunction(a); }
inline double GradJ2SSKStV::dddK(double a) const { return fK->DDDFunction(a); }
inline double GradJ2SSKStV::ddddK(double a) const { return fK->DDDDFunction(a); }

} // namespace Tahoe 
#endif /* _GRAD_J2_SS_KSTV_H_ */
