/* $Id: J2SimoLinHardT.h,v 1.4 2002-07-02 19:56:12 cjkimme Exp $ */
/* created: paklein (06/19/1997)                                          */
/* Interface for a elastoplastic material that is linearly                */
/* isotropically elastic subject to the Huber-von Mises yield             */
/* condition as fYield with kinematic/isotropic hardening laws            */
/* given by:                                                              */
/* 		H(a) = (1 - ftheta) fH_bar a                                      */
/*      K(a) = fYield + ftheta fH_bar a                                   */
/* 		where a is the internal hardening variable                        */
/* 	Note: all calculations are peformed in 3D.                            */

#ifndef _J2_SIMO_LIN_HARD_T_H_
#define _J2_SIMO_LIN_HARD_T_H_

/* base class */
#include "J2PrimitiveT.h"

/* direct members */
#include "dSymMatrixT.h"
#include "dMatrixT.h"
#include "dArrayT.h"
#include "iArrayT.h"

/* forward declarations */

namespace Tahoe {

class ElementCardT;

/** finite strain, J2 elastoplasticity following framework in
 * Simo, J.C. (1988) CMAME v66 and v68. */
class J2SimoLinHardT: public J2PrimitiveT
{
public:

	/** constructor */
	J2SimoLinHardT(ifstreamT& in, int num_ip, double mu);

protected:

	/** elastic/plastic flag values */
	enum LoadStateT {
		kNotInit =-1,
	  kIsPlastic = 0,
	  kIsElastic = 1};

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

protected:

	/* internal state variables */
	dArrayT     fInternal; //internal variables
	dSymMatrixT fb_bar;    //isochoric, elastic part of b
	dSymMatrixT fbeta_bar; //stress surface "center", kinematic hardening

	/** shear modulus */
	double fmu;

	/* internal state variable work space */
	dSymMatrixT fUnitNorm;        //unit normal to the stress surface
	dSymMatrixT fb_bar_trial_;    //unit normal to the stress surface
	dSymMatrixT fbeta_bar_trial_; //unit normal to the stress surface

private:

	/** number of integration points */
	int fNumIP;

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
};

} // namespace Tahoe 
#endif /* _J2_SIMO_LIN_HARD_T_H_ */
