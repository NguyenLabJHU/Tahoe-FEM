/* $Id: GRAD_MRSSNLHardT.h,v 1.7 2005-05-13 22:01:16 kyonten Exp $ */
/* created: Karma Yonten (03/04/2004)                   
   MR version modified to incorporate gradient plasticity 
   theory.
*/
/* Interface for a nonassociative, small strain,     */
/* pressure dependent gradient plasticity model      */
/* with nonlinear isotropic hardening/softening.     */

/*
 *	Note: all calculations are peformed in 3D.
 */

#ifndef _GRAD_MR_SS_NL_HARD_T_H_
#define _GRAD_MR_SS_NL_HARD_T_H_

/* base class */
#include "GRAD_MRPrimitiveT.h"

/* direct members */
#include "dSymMatrixT.h"
#include "dMatrixT.h"
#include "dArrayT.h"

namespace Tahoe 
{

/* forward declarations */
class ElementCardT;

class GRAD_MRSSNLHardT: public GRAD_MRPrimitiveT
{

public:
	/* constructor */
	GRAD_MRSSNLHardT(int num_ip, double mu, double lambda);

	/* status flags */
	enum LoadingStatusT {kIsPlastic = 0,
                         kIsElastic = 1,
                       kIsLocalized = 2, 
                             kReset = 3}; // indicate not to repeat update
                             
	/* returns elastic strain (3D) */
	virtual const dSymMatrixT& ElasticStrain(const dSymMatrixT& totalstrain, 
		const ElementCardT& element, int ip);
		
	/* returns Laplacian of elastic strain (3D) */
	virtual const dSymMatrixT& LapElasticStrain(const dSymMatrixT& lap_totalstrain, 
		const ElementCardT& element, int ip);
		
	/* return correction to stress vector computed by mapping the
	 * stress back to the yield surface, if needed */
	const dSymMatrixT& StressCorrection(const dSymMatrixT& trialstrain,
	    const dSymMatrixT& lap_trialstrain, const dArrayT& dlambda, const dArrayT& lap_dlambda,  
		ElementCardT& element, int ip); // dlam and lap_dlam at the ip     
		
	double& Yield_f(const dArrayT& Sig, const dArrayT& qn, double& ff);
    dArrayT& h_f(const dArrayT& Sig, const dArrayT& qn, dArrayT& hh);
    dArrayT& g_f(const dArrayT& Sig, const dArrayT& qn, const dArrayT& ls, dArrayT& gg);
    dArrayT& n_f(const dArrayT& Sig, const dArrayT& qn, dArrayT& nn);
    dArrayT& r_f(const dArrayT& Sig, const dArrayT& qn, dArrayT& rr);
    dArrayT& m_f(const dArrayT& Sig, const dArrayT& qn, dArrayT& mm);       
    dMatrixT& devstress_f(const dArrayT& Sig, dMatrixT& devstress);
    dMatrixT& dmdSig_f(const dArrayT& qn, dMatrixT& dmdSig);
    dMatrixT& dmdq_f(const dArrayT& Sig, const dArrayT& qn, dMatrixT& dmdq);
    dMatrixT& dhdSig_f(const dArrayT& Sig, const dArrayT& qn, dMatrixT& dhdSig);
    dMatrixT& dgdSig_f(const dArrayT& Sig, const dArrayT& qn, const dArrayT& ls, dMatrixT& dgdSig);
    dMatrixT& dhdq_f(const dArrayT& Sig, const dArrayT& qn, dMatrixT& dhdq);
    dMatrixT& dhdm_f(const dArrayT& Sig, const dArrayT& qn, dMatrixT& dhdm);
    dMatrixT& dgdq_f(const dArrayT& Sig, const dArrayT& qn, const dArrayT& ls, dMatrixT& dgdq);
    
    /* utility function */
	double signof(double& r);

	/* return the correction to moduli due to plasticity (if any)
	 *
	 * Note: Return mapping occurs during the call to StressCorrection.
	 *       The element passed in is already assumed to carry current
	 *       internal variable values */
	const dMatrixT& Moduli(const ElementCardT& element, int ip); 

    /* Modulus for checking perfectly plastic bifurcation */
	const dMatrixT& ModuliPerfPlas(const ElementCardT& element, int ip);
	
	/*@{*/
	const dMatrixT& Moduli_UU1(void) const { return fModuli_UU1; };
	const dMatrixT& Moduli_UU2(void) const { return fModuli_UU2; };
	const dMatrixT& Moduli_ULam1(void) const { return fModuli_ULam1; };
	const dMatrixT& Moduli_ULam2(void) const { return fModuli_ULam2; };
	const dMatrixT& Moduli_LamU1(void) const { return fModuli_LamU1; };
	const dMatrixT& Moduli_LamU2(void) const { return fModuli_LamU2; };
	const dMatrixT& Moduli_LamLam1(void) const { return fModuli_LamLam1; };
	const dMatrixT& Moduli_LamLam2(void) const { return fModuli_LamLam2; };
	/*@{*/
	
	/* return yield condition, f */
	const double& YieldFunction(void) const { return fYield; };

	/* return a pointer to a new plastic element object constructed with
	 * the data from element */
	void AllocateElement(ElementCardT& element);

	enum InternalVariablesT {keps11 = 6, //strains
	                         keps22 = 7,
	                         keps33 = 8,
	                         keps23 = 9,
	                         keps13 = 10,
	                         keps12 = 11,
							 kchi = 30,  // stress-like internal state variable
	                         kc   = 31,
	                      ktanphi = 32,
	                      ktanpsi = 33,
                         kdlambda = 35,  // consistency parameter
                         kplastic = 37,  // Plastic Index
                          kftrial = 34}; // yield function value

	/** internal variables */
	dArrayT& Internal(void) { return fInternal; };
	
	/* element level data */
	void Update(ElementCardT& element);
	void Reset(ElementCardT& element);

	/* returns 1 if the trial elastic strain state lies outside of the 
	 * yield surface */
	int PlasticLoading(const dSymMatrixT& trialstrain, const dSymMatrixT& del2_trialstrain, 
                        ElementCardT& element, int ip);

	/* computes the deviatoric stress corresponding to the given element
	 * and elastic strain.  The function returns a reference to the
	 * stress in fDevStress */
	dSymMatrixT& DeviatoricStress(const dSymMatrixT& trialstrain, 
	          const dSymMatrixT& lap_trialstrain, const ElementCardT& element);

	/* computes the hydrostatic (mean) stress. */
	double MeanStress(const dSymMatrixT& trialstrain, const dSymMatrixT& lap_trialstrain, 
	       const ElementCardT& element);
	
	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

  private:

	/* load element data for the specified integration point */
	void LoadData(const ElementCardT& element, int ip);

  protected:

  	/* element level internal state variables */
  	dSymMatrixT fPlasticStrain; // total plastic strain (deviatoric and volumetric)
  	dSymMatrixT fLapPlasticStrain; // Laplacian of total plastic strain (deviatoric and volumetric)
  	dArrayT     fInternal;      //internal variables

  private:

	/* number of integration points */
	int fNumIP;

  	/* material parameters **/
  	double fmu;
	double flambda;
	double fkappa;
	double fmu_ast;
	double flambda_ast;
	double fkappa_ast;
    double fMeanStress;
    double fLapMeanStress;
  
  	/* return values */
  	dSymMatrixT	fElasticStrain;
  	dSymMatrixT fLapElasticStrain;
  	dSymMatrixT	fStressCorr;
  	dMatrixT	fModuli;
    dMatrixT    fModuliPerfPlas;
    /*@}*/
    dMatrixT	fModuli_UU1;
    dMatrixT	fModuli_UU2;
    dMatrixT	fModuli_ULam1;
    dMatrixT	fModuli_ULam2;
    dMatrixT	fModuli_LamU1;
    dMatrixT	fModuli_LamU2;
    dMatrixT	fModuli_LamLam1;
    dMatrixT	fModuli_LamLam2;
    /*@}*/
    double fYield;
  		
	/* work space */
	dSymMatrixT fDevStress;
	dSymMatrixT fLapDevStress;
	dSymMatrixT fDevStrain; /* deviatoric part of the strain tensor */
	dSymMatrixT fLapDevStrain; /* deviatoric part of the Laplacian of strain tensor */
	dSymMatrixT IdentityTensor2;  

	dMatrixT      fTensorTemp;
	dSymMatrixT   One;  
  	
};


} // namespace Tahoe 
#endif /* _GRAD_MR_SS_NL_HARD_T_H_ */
