/* $Id: J2SSC0HardeningT.h,v 1.8 2005-04-06 23:32:45 thao Exp $ */
#ifndef _J2_SS_C0_HARD_T_H_
#define _J2_SS_C0_HARD_T_H_

/* base class */
#include "J2_C0HardeningT.h"

/* direct members */
#include "dSymMatrixT.h"
#include "dMatrixT.h"
#include "dArrayT.h"
#include "C1FunctionT.h"

namespace Tahoe {

/* forward declarations */
class ElementCardT;
class ifstreamT;

/** small strain J2 plasticity */
class J2SSC0HardeningT: public J2_C0HardeningT
{
public:

	/** constructor */
	J2SSC0HardeningT(void);

protected:

	/** status flags */
	enum LoadingStatusT {kIsPlastic = 0,
                         kIsElastic = 1,
                             kReset = 3}; // indicator not to repeat update

	/* returns elastic strain */
	virtual const dSymMatrixT& ElasticStrain(const dSymMatrixT& totalstrain,
		const ElementCardT& element, int nip, int ip);
			
	/* return the correction to stress vector computed by the mapping the
	 * stress back to the yield surface, if needed */
	const dSymMatrixT& StressCorrection(const dSymMatrixT& trialstrain,
		ElementCardT& element, double mu, int nip, int ip);

	/* return the correction to moduli due to plasticity (if any)
	 *
	 * Note: Return mapping occurs during the call to StressCorrection.
	 *       The element passed in is already assumed to carry current
	 *       internal variable values */
	const dMatrixT& ModuliCorrection(const ElementCardT& element, double mu, int nip, int ip);

	/* return a pointer to a new plastic element object constructed with
	 * the data from element */
	void AllocateElement(ElementCardT& element, int nip);

	enum InternalVariablesT {kalpha = 0,  // equivalent plastic strain
	                    kstressnorm = 1,  // norm of the relative stress
				            kdgamma = 2,  // consistency parameter
				            kftrial = 3}; // yield function value

	/* element level data */
	void Update(ElementCardT& element, int nip);
	void Reset(ElementCardT& element, int nip);

	/*accessors for internal variables*/
	const dSymMatrixT& Get_PlasticStrain(const ElementCardT& element, int nip, int ip);
	const dSymMatrixT& Get_Beta(const ElementCardT& element, int nip, int ip);
	const dArrayT& Get_Internal(const ElementCardT& element, int nip, int ip);

private:

	/* load element data for the specified integration point */
	void LoadData(const ElementCardT& element, int nip, int ip);

	/* returns 1 if the trial elastic strain state lies outside of the
	 * yield surface */
	int PlasticLoading(const dSymMatrixT& trialstrain, ElementCardT& element, double mu, int nip, int ip);

	/* computes the relative stress corresponding for the given element
	 * and elastic strain.  The functions returns a reference to the
	 * relative stress in fRelStress */
	dSymMatrixT& RelativeStress(const dSymMatrixT& totalstrain,
		const ElementCardT& element, double mu);

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

	/* return values */
	dSymMatrixT	fElasticStrain;
	dSymMatrixT	fStressCorr;
	dMatrixT	fModuliCorr;
		
	/* work space */
	dSymMatrixT fRelStress;
	dSymMatrixT fDevStrain; /* deviatoric part of the strain tensor */
	dMatrixT   fTensorTemp;
	
};

} /* namespace Tahoe */

#endif /* _J2_SS_C0_HARD_T_H_ */
