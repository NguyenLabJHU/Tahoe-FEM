/* $Id: J2SimoC0HardeningT.h,v 1.1 2001-05-05 19:28:34 paklein Exp $ */
/* created: paklein (05/01/2001) */

#ifndef _J2_SIMO_C0_HARD_T_H_
#define _J2_SIMO_C0_HARD_T_H_

/* environment */
#include "Environment.h"

/* direct members */
#include "dSymMatrixT.h"
#include "dMatrixT.h"
#include "dArrayT.h"
#include "iArrayT.h"
#include "C1FunctionT.h"

/* forward declarations */
class ElementCardT;
class ifstreamT;
#include "ios_fwd_decl.h"

/** finite strain, J2 elastoplasticity following framework in
 * Simo, J.C. (1988) CMAME v66 and v68. */
class J2SimoC0HardeningT
{
public:

	/** constructor */
	J2SimoC0HardeningT(ifstreamT& in, int num_ip, double mu);

	/** destructor */
	~J2SimoC0HardeningT(void);

protected:

	/** number of stored variables */
	static const int kNumInternal;

	/** index of stored values */
	enum VariablesT {
        kalpha = 0, /* isotropic hardening         */
   kstressnorm = 1, /* norm of the relative stress */
       kdgamma = 2, /* consistency parameter       */
       kftrial = 3, /* yield function value        */
       kmu_bar = 4, 
   kmu_bar_bar = 5,
     kDetF_tot = 6};/* determinant of total F      */

	/** hardening function types */
	enum HardeningFunctionT {
               kLinear = 0,
    kLinearExponential = 1,
          kCubicSpline = 2};

	/** write parameters */
	void Print(ostream& out) const;
	void PrintName(ostream& out) const;
	
	/** compute trial elastic state - return reference to isochoric,
	 * trial elastic stretch */
	const dSymMatrixT& TrialElasticState(const dMatrixT& F_total,
		const dMatrixT& f_relative, ElementCardT& element, int ip);
			
	/** return the correction to stress vector computed by the mapping the
	 * stress back to the yield surface */
	const dSymMatrixT& StressCorrection(const dMatrixT& F_total,
		const dMatrixT& f_relative, ElementCardT& element, int ip);

	/** return the correction to moduli due to plasticity (if any)
	 *
	 * Note: Return mapping occurs during the call to StressCorrection.
	 *       The element passed in is already assumed to carry current
	 *       internal variable values. */
	const dMatrixT& ModuliCorrection(ElementCardT& element, int ip);

	/** allocate element storage */
	void AllocateElement(ElementCardT& element);

	/** element level data */
	void Update(ElementCardT& element);
	void Reset(ElementCardT& element);

private:

	/** initialize intermediate state from F_n (for ) */
	void InitIntermediate(const dMatrixT& F_total, const dMatrixT& f_relative);

	/** load element data for the specified integration point */
	void LoadData(const ElementCardT& element, int ip);

	/** returns 1 if the trial elastic strain state lies outside of the
	 * yield surface - assumes trial elastic state has been set first */
	int PlasticLoading(const dMatrixT& F_total, const dMatrixT& f_relative,
		ElementCardT& element, int ip);

	/** evaluate single parameters yield function */
	double YieldCondition(const dSymMatrixT& stress, double alpha) const;

	/* hardening functions and their 1st derivatives */
	double   H(double) const { return 0.0; }; // no kinematic hardening yet
	double  dH(double) const { return 0.0; };
	double   K(double a) const;
	double  dK(double a) const;

	/* construct isotropic hardening function */
	void ConstructHardeningFunction(ifstreamT& in);

protected:

	/* internal state variables */
	dArrayT     fInternal; //internal variables
	dSymMatrixT fb_bar;    //isochoric, elastic part of b
	dSymMatrixT fbeta_bar; //stress surface "center", kinematic hardening

private:

	/** number of integration points */
	int fNumIP;

	/** shear modulus */
	double fmu;

	/** C1 isotropic hardening function */
	HardeningFunctionT fType;
	C1FunctionT* fK;	

	/* return values */
	dSymMatrixT	fStressCorr;
	dMatrixT	fModuliCorr;

	/* work space */
	dSymMatrixT fRelStress;
	
	dMatrixT    f_f_bar;
	dSymMatrixT fb_bar_trial;
	
	/* kinematic hardening */
	double      ftrace_beta_trial;
	dSymMatrixT fbeta_bar_trial;

	dMatrixT   fMatrixTemp1;
	dMatrixT   fMatrixTemp2;

	dSymMatrixT fRed2Temp;
	dMatrixT    fRed4Temp1;
	dMatrixT    fRed4Temp2;

	/* internal state variable work space */
	dSymMatrixT fUnitNorm;        //unit normal to the stress surface
	dSymMatrixT fb_bar_trial_;    //unit normal to the stress surface
	dSymMatrixT fbeta_bar_trial_; //unit normal to the stress surface
};

/* hardening functions and their 1st derivatives */
inline double J2SimoC0HardeningT::K(double a) const { return fK->Function(a); }
inline double J2SimoC0HardeningT::dK(double a) const { return fK->DFunction(a); }

#endif /* _J2_SIMO_C1_HARD_T_H_ */
