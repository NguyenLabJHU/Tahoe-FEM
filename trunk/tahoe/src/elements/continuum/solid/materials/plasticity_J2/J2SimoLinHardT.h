/* $Id: J2SimoLinHardT.h,v 1.1.1.1 2001-01-29 08:20:30 paklein Exp $ */
/* created: paklein (06/19/1997)                                          */
/* NOTE: THIS FORMULATION OF PLASTICITY RESULTS IN A NON-SYMMETRIC        */
/* CONSISTENT TANGENT MODULI MATRIX. DO NOT USE UNLESS THE                */
/* ELEMENT FORMULATION AND GLOBAL EQUATION MATRIX HAVE BEEN               */
/* CHECK FOR ASSUMPTIONS ABOUT SYMMETRY                                   */
/* Interface for a elastoplastic material that is linearly                */
/* isotropically elastic subject to the Huber-von Mises yield             */
/* condition as fYield with kinematic/isotropic hardening laws            */
/* given by:                                                              */
/* 		H(a) = (1 - ftheta) fH_bar a                                         */
/* K(a) = fYield + ftheta fH_bar a                                        */
/* 		where a is the internal hardening variable                           */
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
class ElementCardT;

class J2SimoLinHardT: public J2PrimitiveT
{
public:

	/* constructor */
	J2SimoLinHardT(ifstreamT& in, int num_ip, double mu);
	
	/* returns isochoric elastic stretch */
	const dSymMatrixT& ElasticStretch(const dMatrixT& F_total,
		const dMatrixT& f_relative, ElementCardT& element, int ip);
			
	/* return the correction to stress vector computed by the mapping the
	 * stress back to the yield surface */
	const dSymMatrixT& StressCorrection(const dMatrixT& F_total,
		const dMatrixT& f_relative, ElementCardT& element, int ip);

	/* return the correction to moduli due to plasticity (if any)
	 *
	 * Note: Return mapping occurs during the call to StressCorrection.
	 *       The element passed in is already assumed to carry current
	 *       internal variable values. */
	const dMatrixT& ModuliCorrection(ElementCardT& element, int ip);

	/* allocate element storage */
	void AllocateElement(ElementCardT& element);

protected:

	/* element level data */
	void Update(ElementCardT& element);
	void Reset(ElementCardT& element);

private:

	/* initialize intermediate state from F_n (for ) */
	void InitIntermediate(const dMatrixT& F_total, const dMatrixT& f_relative);

	/* load element data for the specified integration point */
	void LoadData(const ElementCardT& element, int ip);

	/* returns 1 if the trial elastic strain state lies outside of the
	 * yield surface */
	int PlasticLoading(const dMatrixT& F_total, const dMatrixT& f_relative,
		ElementCardT& element, int ip);

	/* computes the relative stress and returns a reference to the
	 * relative stress in fRelStress */
	dSymMatrixT& RelativeStress(const dMatrixT& F_total,
		const dMatrixT& f_relative, const ElementCardT& element);
		
private:

	/* number of integration points */
	int fNumIP;

	/* material parameters */
	double fmu;

	/* return values */
	dSymMatrixT	fbelastic;
	dSymMatrixT	fStressCorr;
	dMatrixT	fModuliCorr;

	/* work space */
	dSymMatrixT fRelStress;
	dSymMatrixT fDev_b;
	
	dMatrixT f_f_bar;
	dMatrixT f_b_trial;

	dMatrixT   fMatrixTemp1;
	dMatrixT   fMatrixTemp2;

	dSymMatrixT fRed2Temp;
	dMatrixT   fRed4Temp1;
	dMatrixT   fRed4Temp2;
	
	dSymMatrixT fb_bar;    //isochoric, elastic part of b
	dSymMatrixT fUnitNorm; //unit normal to the stress surface
	dSymMatrixT fBeta;     //stress surface "center", kinematic hardening
	dArrayT    fInternal; //internal variables
};

#endif /* _J2_SIMO_LIN_HARD_T_H_ */
