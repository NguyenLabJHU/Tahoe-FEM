/* $Id: DPSSLinHardT.h,v 1.1.1.1 2001-01-29 08:20:30 paklein Exp $ */
/* created: myip (06/01/1999)                                             */
/* Interface for a elastoplastic material that is linearly                */
/* isotropically elastic subject to the Huber-von Mises yield             */
/* condition as fYield with isotropic hardening                           */
/* 	Note: all calculations are peformed in 3D.                            */

#ifndef _DP_SS_LIN_HARD_T_H_
#define _DP_SS_LIN_HARD_T_H_

/* base class */
#include "DPPrimitiveT.h"

/* direct members */
#include "dSymMatrixT.h"
#include "dMatrixT.h"
#include "dArrayT.h"

/* forward declarations */
class ElementCardT;

class DPSSLinHardT: public DPPrimitiveT
{
public:

	/* constructor */
	DPSSLinHardT(ifstreamT& in, int num_ip, double mu, double lambda);

	/* output name */
	virtual void PrintName(ostream& out) const;

protected:

	/* status flags */
	enum LoadingStatusT {kIsPlastic = 0,
                         kIsElastic = 1,
                             kReset = 3}; // indicate not to repeat update

	/* returns elastic strain (3D) */
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

	enum InternalVariablesT {kalpha_dev = 0,  // deviatoric internal variable
				             kalpha_vol = 1,  // volumetric internal variable
	                        kstressnorm = 2,  // norm of the relative stress
				                kdgamma = 3,  // consistency parameter
				                kftrial = 4}; // yield function value
	
	/* element level data */
	void Update(ElementCardT& element);
	void Reset(ElementCardT& element);
	

	//TEMP
	/* returns 1 if the trial elastic strain state lies outside of the
	 * yield surface */
	int PlasticLoading(const dSymMatrixT& trialstrain, const ElementCardT& element,
		int ip);

	/* computes the deviatoric stress corresponding for the given element
	 * and elastic strain.  The functions returns a reference to the
	 * relative stress in fDevStress */
	dSymMatrixT& DeviatoricStress(const dSymMatrixT& trialstrain,
		const ElementCardT& element);

	/* computes the hydrostatic (mean) stress. */    //**mien
	double MeanStress(const dSymMatrixT& trialstrain,
		const ElementCardT& element);
	//TEMP

private:

	/* load element data for the specified integration point */
	void LoadData(const ElementCardT& element, int ip);

	/* returns 1 if the trial elastic strain state lies outside of the
	 * yield surface */
	//	int PlasticLoading(const dSymMatrixT& trialstrain, ElementCardT& element,
	//		int ip);

	/* computes the deviatoric stress corresponding for the given element
	 * and elastic strain.  The functions returns a reference to the
	 * relative stress in fDevStress */
	//	dSymMatrixT& DeviatoricStress(const dSymMatrixT& totalstrain,
	//		const ElementCardT& element);

	/* computes the hydrostatic (mean) stress. */    //**mien
	//	double& MeanStress(const dSymMatrixT& totalstrain,
	//		const ElementCardT& element);

protected:

	/* element level internal variables */
	dSymMatrixT fPlasticStrain; //deviatoric part of the plastic strain
	dSymMatrixT fUnitNorm;      //unit normal to the stress surface
	dArrayT     fInternal;      //internal variables

private:

	/* number of integration points */
	int fNumIP;

	/* material parameters **/
	double fmu;
	double flambda;
	double fkappa;
	double fX_H;
double fMeanStress;

	/* return values */
	dSymMatrixT	fElasticStrain;
	dSymMatrixT	fStressCorr;
	dMatrixT	fModuliCorr;
		
	/* work space */
	dSymMatrixT fDevStress;
	dSymMatrixT fDevStrain; /* deviatoric part of the strain tensor */
	dSymMatrixT  IdentityTensor2;  //**mien

	dMatrixT   fTensorTemp;
	dSymMatrixT   One;  //**mien
	
};

#endif /* _DP_SS_LIN_HARD_T_H_ */
