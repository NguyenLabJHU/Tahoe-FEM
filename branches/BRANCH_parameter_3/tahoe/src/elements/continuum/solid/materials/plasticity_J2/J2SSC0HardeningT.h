/* $Id: J2SSC0HardeningT.h,v 1.6.18.1 2004-06-08 16:01:34 paklein Exp $ */
#ifndef _J2_SS_C0_HARD_T_H_
#define _J2_SS_C0_HARD_T_H_

/* base class */
#include "ParameterInterfaceT.h"

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

class J2SSC0HardeningT: virtual public ParameterInterfaceT
{
public:

	/** constructor */
	J2SSC0HardeningT(void);

	/* destructor */
	virtual ~J2SSC0HardeningT(void);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** return the description of the given inline subordinate parameter list */
	virtual void DefineInlineSub(const StringT& sub, ParameterListT::ListOrderT& order, 
		SubListT& sub_sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& list_name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** status flags */
	enum LoadingStatusT {kIsPlastic = 0,
                         kIsElastic = 1,
                             kReset = 3}; // indicator not to repeat update

	/** hardening function types */
	enum HardeningFunctionT {
               kLinear = 0,
    kLinearExponential = 1,
          kCubicSpline = 2};

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

	enum InternalVariablesT {kalpha = 0,  // equivalent plastic strain
	                    kstressnorm = 1,  // norm of the relative stress
				            kdgamma = 2,  // consistency parameter
				            kftrial = 3}; // yield function value

	/* element level data */
	void Update(ElementCardT& element);
	void Reset(ElementCardT& element);

	/* hardening functions and their 1st derivatives */
	double   H(double) const { return 0.0; }; // no kinematic hardening yet
	double  dH(double) const { return 0.0; };
	double  K(double a) const;
	double dK(double a) const;

	/** returns the value value of the yield function given the
	 * Cauchy stress vector and state variables, where beta and alpha
	 * represent the kinematic and isotropic hardening, respectively */
	double YieldCondition(const dSymMatrixT& relstress, double alpha) const;

private:

	/* load element data for the specified integration point */
	void LoadData(const ElementCardT& element, int ip);

	/* returns 1 if the trial elastic strain state lies outside of the
	 * yield surface */
	int PlasticLoading(const dSymMatrixT& trialstrain, ElementCardT& element, int intpt);

	/* computes the relative stress corresponding for the given element
	 * and elastic strain.  The functions returns a reference to the
	 * relative stress in fRelStress */
	dSymMatrixT& RelativeStress(const dSymMatrixT& totalstrain,
		const ElementCardT& element);

	/* construct isotropic hardening function */
//	void ConstructHardeningFunction(ifstreamT& in);

protected:

	/* element level internal variables */
	dSymMatrixT fPlasticStrain; //deviatoric part of the plastic strain
	dSymMatrixT fUnitNorm;      //unit normal to the stress surface
	dSymMatrixT fBeta;          //stress surface "center", kinematic hardening
	dArrayT     fInternal;      //internal variables

private:

	/* number of integration points */
	int fNumIP;

	/* material parameters **/
	double fmu;
	
	/** kinematic/isotropic hardening mixity (unused) */
	double ftheta;

	/** C1 isotropic hardening function */
	bool fIsLinear;
	C1FunctionT* fK;	

	/* return values */
	dSymMatrixT	fElasticStrain;
	dSymMatrixT	fStressCorr;
	dMatrixT	fModuliCorr;
		
	/* work space */
	dSymMatrixT fRelStress;
	dSymMatrixT fDevStrain; /* deviatoric part of the strain tensor */
	dMatrixT   fTensorTemp;
	
};

/* hardening functions and their 1st derivatives */
inline double J2SSC0HardeningT::K(double a) const { return fK->Function(a); }
inline double J2SSC0HardeningT::dK(double a) const { return fK->DFunction(a); }

} // namespace Tahoe 
#endif /* _J2_SS_C0_HARD_T_H_ */
