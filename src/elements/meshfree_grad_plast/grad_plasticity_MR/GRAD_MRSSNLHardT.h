/* created: Majid T. Manzari (04/16/2003)            */
/*  
 * Interface for a nonassociative, small strain,     */
/* pressure dependent gradient plasticity model      */
/* with nonlinear isotropic hardening/softening.     */

/*
 *	Note: all calculations are peformed in 3D.
 */

#ifndef _GRAD_MR_SS_LIN_HARD_T_H_
#define _GRAD_MR_SS_LIN_HARD_T_H_

/* base class */
#include "GRAD_MRPrimitiveT.h"

/* direct members */
#include "dSymMatrixT.h"
#include "dMatrixT.h"
#include "dArrayT.h"
#include "LAdMatrixT.h" //for solving Ax=b linear system of equations

namespace Tahoe 
{

/* forward declarations */
class ElementCardT;

class GRAD_MRSSNLHardT: public GRAD_MRPrimitiveT
{

public:
	/* constructor */
	GRAD_MRSSNLHardT(ifstreamT& in, int num_ip, double mu, double lambda);

  	/* output name */
	virtual void PrintName(ostream& out) const;

  protected:

	/* status flags */
	enum LoadingStatusT {kIsPlastic = 0,
                         kIsElastic = 1,
                             kReset = 3}; // indicate not to repeat update
                             
	/* return correction to stress vector computed by mapping the
	 * stress back to the yield surface, if needed */
	const dSymMatrixT& StressCorrection(const dSymMatrixT& trialstrain,
	    const dSymMatrixT& del2_trialstrain, double dlam, double del2_dlam,  
		ElementCardT& element, int ip); // dlam and del2_dlam at the ip
		
	virtual const dSymMatrixT& ElasticStrain(const dSymMatrixT& totalstrain, 
	        const dSymMatrixT& del2_totalstrain, double dlam, double del2_dlam,
	        const ElementCardT& element, int ip);
		
	double& Yield_f(const dArrayT& Sig, const dArrayT& qn, double& ff);
    dArrayT& h_f(const dArrayT& Sig, const dArrayT& qn, dArrayT& hh);
    dArrayT& g_f(const dArrayT& Sig, const dArrayT& qn, double lsp_v, double lsp_s, dArrayT& gg);
    dArrayT& nn_f(const dArrayT& Sig, const dArrayT& qn, dArrayT& nn);
    dArrayT& rr_f(const dArrayT& Sig, const dArrayT& qn, dArrayT& rr);
    dArrayT& dm_f(const dArrayT& Sig, const dArrayT& qn, dArrayT& mm);       
    dMatrixT& dmdSig_f(const dArrayT& qn, dMatrixT& dmdSig);
    dMatrixT& dmdq_f(const dArrayT& Sig, const dArrayT& qn, dMatrixT& dmdq);
    dMatrixT& dhdSig_f(const dArrayT& Sig, const dArrayT& qn, dMatrixT& dhdSig);
    dMatrixT& dgdSig_f(const dArrayT& Sig, const dArrayT& qn, double lsp_v, double lsp_s, dMatrixT& dgdSig);
    dMatrixT& dhdq_f(const dArrayT& Sig, const dArrayT& qn, dMatrixT& dhdq);
    dMatrixT& dgdq_f(const dArrayT& Sig, const dArrayT& qn, double lsp_v, double lsp_s, dMatrixT& dqbardq);
    
    LAdMatrixT dRR; //set up matrix A for solving Ax=b system
    /* utility function */
	double signof(double& r);

	/* return the correction to moduli due to plasticity (if any)
	 *
	 * Note: Return mapping occurs during the call to StressCorrection.
	 *       The element passed in is already assumed to carry current
	 *       internal variable values */
	const dMatrixT& Moduli(const ElementCardT& element, int ip); 

        /* Modulus for checking discontinuous bifurcation */

	const dMatrixT& ModuliDisc(const ElementCardT& element, int ip);

	/* return a pointer to a new plastic element object constructed with
	 * the data from element */
	void AllocateElement(ElementCardT& element);

	enum InternalVariablesT {kchi = 30,  // stress-like internal state variable
	                         kc   = 31,
	                      ktanphi = 32,
	                      ktanpsi = 33,
                         kdlambda = 35,  // consistency parameter
                         kplastic = 37,  // Plastic Index
                          kftrial = 34}; // yield function value

	/* element level data */
	void Update(ElementCardT& element);
	void Reset(ElementCardT& element);

	/* returns 1 if the trial elastic strain state lies outside of the 
	 * yield surface */
	int PlasticLoading(const dSymMatrixT& trialstrain, const dSymMatrixT& del2_trialstrain, 
                       double dlam, double del2_lam, ElementCardT& element,int ip);

	/* computes the deviatoric stress corresponding to the given element
	 * and elastic strain.  The function returns a reference to the
	 * stress in fDevStress */
	dSymMatrixT& DeviatoricStress(const dSymMatrixT& trialstrain, 
	          const dSymMatrixT& del2_trialstrain, double dlam, double del2_dlam 
		const ElementCardT& element);

	/* computes the hydrostatic (mean) stress. */
	double MeanStress(const dSymMatrixT& trialstrain, const dSymMatrixT& del2_trialstrain, 
	       double dlam, double del2_dlam, const ElementCardT& element);

  private:

	/* load element data for the specified integration point */
	void LoadData(const ElementCardT& element, int ip);

  protected:

  	/* element level internal state variables */
  	dSymMatrixT fPlasticStrain; //total plastic strain (deviatoric and volumetric)
  	dSymMatrixT fGradPlasticStrain; //gradient total plastic strain (deviatoric and volumetric)
  	dArrayT     fInternal;      //internal variables

  private:

	/* number of integration points */
	int fNumIP;

  	/* material parameters **/
  	double fmu;
	double flambda;
	double fkappa;
    double fMeanStress;
    //include length scale parameters here??
  
  	/* return values */
  	dSymMatrixT	fElasticStrain;
  	dSymMatrixT fGradElasticStrain
  	dSymMatrixT	fStressCorr;
  	dMatrixT	fModuli;
    dMatrixT    fModuliDisc;
    double fYieldFunction;
  		
	/* work space */
	dSymMatrixT fDevStress;
	dSymMatrixT fDevStrain; /* deviatoric part of the strain tensor */
	dSymMatrixT IdentityTensor2;  

	dMatrixT      fTensorTemp;
	dSymMatrixT   One;  
  	
};

} // namespace Tahoe 
#endif /* _GRAD_MR_SS_NL_HARD_T_H_ */
