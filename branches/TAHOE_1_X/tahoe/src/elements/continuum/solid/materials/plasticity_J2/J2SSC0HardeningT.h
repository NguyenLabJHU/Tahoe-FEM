/* $Id: J2SSC0HardeningT.h,v 1.6.28.2 2005-04-05 23:33:30 thao Exp $ */

#ifndef _J2_SS_C0_HARD_T_H_
#define _J2_SS_C0_HARD_T_H_

/* direct members */
#include "dSymMatrixT.h"
#include "dMatrixT.h"
#include "dArrayT.h"
#include "C1FunctionT.h"
#include "iArrayT.h"

#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
class ElementCardT;
class ifstreamT;

class J2SSC0HardeningT
{
public:

	/* constructor */
	J2SSC0HardeningT(ifstreamT& in, int num_ip, double mu);

	/* destructor */
	virtual ~J2SSC0HardeningT(void);

	/* output name */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;

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
	/*accessors for internal variables*/
	const dSymMatrixT& Get_PlasticStrain(const ElementCardT& element, int ip);
	const dSymMatrixT& Get_Beta(const ElementCardT& element, int ip);
	const dArrayT& Get_Internal(const ElementCardT& element, int ip);

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
	void ConstructHardeningFunction(ifstreamT& in);

protected:

	/* element level internal variables */
	dSymMatrixT fPlasticStrain; //deviatoric part of the plastic strain
	dSymMatrixT fUnitNorm;      //unit normal to the stress surface
	dSymMatrixT fBeta;          //stress surface "center", kinematic hardening
	dArrayT     fInternal;      //internal variables
	
	/*internal variables and conjugates for material force calculations*/
	dArrayT fInternalStressVars;
	dArrayT fInternalStrainVars;
	iArrayT fInternalDOF;

private:

	/* number of integration points */
	int fNumIP;

	/* material parameters **/
	double fmu;
	
	/** kinematic/isotropic hardening mixity (unused) */
	double ftheta;

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
	dMatrixT   fTensorTemp;
	
};

/* hardening functions and their 1st derivatives */
inline double J2SSC0HardeningT::K(double a) const { return fK->Function(a); }
inline double J2SSC0HardeningT::dK(double a) const { return fK->DFunction(a); }

} // namespace Tahoe 
#endif /* _J2_SS_C0_HARD_T_H_ */
