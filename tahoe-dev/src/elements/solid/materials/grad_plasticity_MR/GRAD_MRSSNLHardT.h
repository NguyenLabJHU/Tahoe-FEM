/* $Id: GRAD_MRSSNLHardT.h,v 1.16 2006-08-29 21:18:41 kyonten Exp $ */
/* created: Karma Yonten (03/04/2004)                   
   Gradient Enhanced MR Model
*/
/* interface for a nonassociative, small strain,     */
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
                             kReset = 3,}; // indicate not to repeat update 
                             
	/* returns elastic strain (3D) */
	virtual const dSymMatrixT& ElasticStrain(const dSymMatrixT& totalstrain, 
		const ElementCardT& element, int ip);
		
	/* returns laplacian of elastic strain (3D) */
	virtual const dSymMatrixT& LapElasticStrain(const dSymMatrixT& lap_totalstrain, 
		const ElementCardT& element, int ip);
		
	/* return correction to stress vector computed by mapping the
	 * stress back to the yield surface, if needed */
	const dSymMatrixT& StressCorrection(const dSymMatrixT& totalstrain,
	    const dSymMatrixT& lap_totalstrain, const dArrayT& lambda, const dArrayT& lap_lambda,  
		ElementCardT& element, int ip);     
		
	void n_f(const dSymMatrixT& Sig, const dArrayT& qn, dArrayT& dfdSig);
    void r_f(const dSymMatrixT& Sig, const dArrayT& qn, dArrayT& dfdq);
    void m_f(const dSymMatrixT& Sig, const dArrayT& qn, dArrayT& dQdSig);       
    void dmdSig_f(const dSymMatrixT& Sig, const dArrayT& qn, dMatrixT& dmdSig);
    void dmdq_f(const dSymMatrixT& Sig, const dArrayT& qn, dMatrixT& dmdq);
    void h_f(const dSymMatrixT& Sig, const dArrayT& qn, dArrayT& hh);
    void dhdSig_f(const dSymMatrixT& Sig, const dArrayT& qn, dMatrixT& dhdSig);
    void dhdq_f(const dSymMatrixT& Sig, const dArrayT& qn, dMatrixT& dhdq);
    void dhdm_f(const dSymMatrixT& Sig, const dArrayT& qn, dMatrixT& dhdm);
    void g_f(const dSymMatrixT& Sig, const dArrayT& qn, dArrayT& gg);
    void dgdSig_f(const dSymMatrixT& Sig, const dArrayT& qn, dMatrixT& dgdSig);
    void dgdq_f(const dSymMatrixT& Sig, const dArrayT& qn, dMatrixT& dgdq);
    
    /* utility function */
	double signof(double r);
	
	/* double contraction of fourth order tensor with second
	 * order tensor    */
	void Contract4To2(const dMatrixT& mat,const dSymMatrixT& vec, 
	                  dSymMatrixT& res);

	/* return the correction to moduli due to plasticity (if any)
	 *
	 * Note: return mapping occurs during the call to StressCorrection.
	 *       The element passed in is already assumed to carry current
	 *       internal variable values */
	const dMatrixT& Moduli(const ElementCardT& element, int ip); 

    /* modulus for checking perfectly plastic bifurcation */
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

	enum InternalVariablesT {kchi = 0,  // stress-like internal state variable
	                         kc   = 1,
	                      ktanphi = 2,
	                      ktanpsi = 3,
                          kftrial = 4, /* yield function value */
                          klambda = 5, /* plastic multiplier */
                       klaplambda = 6,}; /* laplacian of plastic multiplier */
    
	/** internal variables */
	dArrayT& Internal(void) { return fInternal; };
	
	/* initial internal variables */
	dArrayT& IniInternal(void) { return fIniInternal; };
	
	/* element level data */
	void Update(ElementCardT& element);
	void Reset(ElementCardT& element);

	/* returns 1 if the trial elastic strain state lies outside of the 
	 * yield surface */
	int PlasticLoading(const dSymMatrixT& totalstrain, const dSymMatrixT& lap_trialstrain, 
                        const dArrayT& lambda, ElementCardT& element, int ip);

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
  	dSymMatrixT fStress;           // updated stress
  	dArrayT     fInternal;      // internal variables
  	dArrayT     fIniInternal;      // initial internal variables

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
    double fYield;
    /*@}*/
  		
	/* work space */
	dSymMatrixT fDevStress;
	dSymMatrixT fLapDevStress;
	dSymMatrixT fDevStrain; /* deviatoric part of the strain tensor */
	dSymMatrixT fLapDevStrain; /* deviatoric part of the laplacian of strain tensor */
};


} // namespace Tahoe 
#endif /* _GRAD_MR_SS_NL_HARD_T_H_ */
