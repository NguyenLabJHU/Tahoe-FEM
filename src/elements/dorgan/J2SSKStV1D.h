/* $Id: J2SSKStV1D.h,v 1.8 2004-07-15 08:28:12 paklein Exp $ */
#ifndef _J2_SS_KSTV_1D_H_
#define _J2_SS_KSTV_1D_H_

/* base classes */
#include "SSIsotropicMatT.h"
#include "HookeanMatT.h"

/* direct members */
#include "dSymMatrixT.h"
#include "dMatrixT.h"
#include "dArrayT.h"
#include "C1FunctionT.h"

#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
class ElementCardT;
class ifstreamT;

class J2SSKStV1D: public SSIsotropicMatT, public HookeanMatT
{
public:

	/* constructor */
	J2SSKStV1D(ifstreamT& in, const SSMatSupportT& support);
	J2SSKStV1D(void);

	/* initialization */
	virtual void Initialize(void);

	/* update internal variables */
	virtual void UpdateHistory(void);

	/* reset internal variables to last converged solution */
	virtual void ResetHistory(void);

	/** \name spatial description */
	/*@{*/
	/** spatial tangent modulus */
	virtual const dMatrixT& c_ijkl(void);

	/** Cauchy stress */
	virtual const dSymMatrixT& s_ij(void);

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

protected:

	/* set modulus */
	virtual void SetModulus(dMatrixT& modulus);

private:

	/** status flags */
	enum LoadingStatusT { kIsPlastic = 0,
						  kIsElastic = 1,
						  kReset     = 3}; // indicator not to repeat update

	/** hardening function types */
	enum HardeningFunctionT { kLinear            = 0,
							  kLinearExponential = 1,
							  kCubicSpline       = 2};
	 	 	
	/* returns elastic strain */
	virtual const dSymMatrixT& ElasticStrain(const dSymMatrixT& totalstrain,
						 const ElementCardT& element, int ip);
		
	/* return the correction to stress vector computed by the mapping the
	 * stress back to the yield surface, if needed */
	const dSymMatrixT& StressCorrection(const dSymMatrixT& trialstrain,
					    ElementCardT& element, int ip);

	/* return the correction to moduli due to plasticity (if any)
	 *
	 * Note: Return mapping occurs during the call to StressCorrection.
	 *       The element passed in is already assumed to carry current
	 *       internal variable values */
	const dMatrixT& ModuliCorrection(const ElementCardT& element, int ip);

	/* return a pointer to a new plastic element object constructed with
	 * the data from element */
	void AllocateElement(ElementCardT& element);

	enum InternalVariablesT { kalpha      = 0,  // equivalent plastic strain
							  kstressnorm = 1,  // norm of the relative stress
							  kdgamma     = 2,  // consistency parameter
							  kftrial     = 3}; // yield function value

	/* element level data */
	void Update(ElementCardT& element);
	void Reset(ElementCardT& element);

	/* hardening functions and their 1st derivatives */
	double H(double) const { return 0.0; }; // no kinematic hardening yet
	double dH(double) const { return 0.0; };
	double K(double a) const;
	double dK(double a) const;

	/** returns the value value of the yield function given the
	 * Cauchy stress vector and state variables, where beta and alpha
	 * represent the kinematic and isotropic hardening, respectively */
	double YieldCondition(double alpha) const;

	/* load element data for the specified integration point */
	void LoadData(const ElementCardT& element, int ip);

	/* returns 1 if the trial elastic strain state lies outside of the
	 * yield surface */
	int PlasticLoading(const dSymMatrixT& trialstrain, ElementCardT& element,int intpt);

	/* construct isotropic hardening function */
	void ConstructHardeningFunction(ifstreamT& in);

private:

	/* return values */
	dSymMatrixT	fStress;
	dMatrixT	fModulus;

	/* workspaces */
	dSymMatrixT	fStress_3D;

	/* element level internal variables */
	dSymMatrixT fPlasticStrain; //deviatoric part of the plastic strain
	dArrayT     fInternal;      //internal variables

	/* number of integration points */
	int fNumIP;

	/* material parameters **/
	double fYoung;
	
	/** C1 isotropic hardening function */
	HardeningFunctionT fType;
	C1FunctionT* fK;	

	/* return values */
	dSymMatrixT	fElasticStrain;
	dSymMatrixT	fStressCorr;
	dMatrixT	fModuliCorr;
		
	/* work space */
	dSymMatrixT fRelStress;
	dSymMatrixT fDevStrain; /* deviatoric part of the strain tensor */
	dMatrixT    fTensorTemp;

};

/* hardening functions and their 1st derivatives */
inline double J2SSKStV1D::K(double a) const { return fK->Function(a); }
inline double J2SSKStV1D::dK(double a) const { return fK->DFunction(a); }

} // namespace Tahoe 
#endif /* _J2_SS_KSTV_1D_H_ */
