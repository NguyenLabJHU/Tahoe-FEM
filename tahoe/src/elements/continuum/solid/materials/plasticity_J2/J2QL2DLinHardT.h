/* $Id: J2QL2DLinHardT.h,v 1.1.1.1 2001-01-29 08:20:30 paklein Exp $ */
/* created: paklein (06/29/1997)                                          */
/* Interface for a elastoplastic material that is linearly                */
/* isotropically elastic subject to the Huber-von Mises yield             */
/* condition as fYield with kinematic/isotropic hardening laws            */
/* given by:                                                              */
/* 		H(a) = (1 - ftheta) fH_bar a                                         */
/* K(a) = fYield + ftheta fH_bar a                                        */
/* 		where a is the internal hardening variable                           */
/* 	Note: all calculations are peformed in 3D.                            */

#ifndef _J2_QL_LIN_HARD_2D_T_H_
#define _J2_QL_LIN_HARD_2D_T_H_

/* base classes */
#include "QuadLog2D.h"
#include "J2PrimitiveT.h"

/* direct members */
#include "dSymMatrixT.h"
#include "dMatrixT.h"
#include "dArrayT.h"
#include "iArrayT.h"

class J2QL2DLinHardT: public QuadLog2D, public J2PrimitiveT
{
public:

	/* constructor */
	J2QL2DLinHardT(ifstreamT& in, const ElasticT& element);

	/* update internal variables */
	virtual void UpdateHistory(void);

	/* reset internal variables to last converged solution */
	virtual void ResetHistory(void);

	/* print parameters */
	virtual void Print(ostream& out) const;

	/* modulus */
	virtual const dMatrixT& c_ijkl(void);
	
	/* stress */
	virtual const dSymMatrixT& s_ij(void);

	/* strain energy density */
	virtual double StrainEnergyDensity(void);

	/* required parameter flags */
	virtual bool NeedLastDisp(void) const;

	/* returns the number of variables computed for nodal extrapolation
	 * during for element output, ie. internal variables */
	virtual int NumOutputVariables(void) const;
	virtual void OutputLabels(ArrayT<StringT>& labels) const;
	virtual void ComputeOutput(dArrayT& output);
	
protected:

	/* returns the trial stretch */
	const dSymMatrixT& TrialStretch(const dMatrixT& F_total,
		const dMatrixT& f_relative, int ip);
			
	/* apply return mapping if necessary */
	void ReturnMapping(const dSymMatrixT& b_tr, dArrayT& beta, int ip);

	/*
	 * Return the correction to moduli due to plasticity (if any)
	 *
	 * Note: Return mapping occurs during the call to StressCorrection.
	 *       The element passed in is already assumed to carry current
	 *       internal variable values.
	 */
	void ElastoPlasticCorrection(dMatrixT& a_ep, dArrayT& beta, int ip);

	/* allocate element storage */
	void AllocateElement(ElementCardT& elememt);

private:

	/* compute F_total and f_relative */
	void ComputeGradients(void);

	/* initialize intermediate state from F_n (for ) */
	void InitIntermediate(const dMatrixT& F_total, const dMatrixT& f_relative);

	/* load element data for the specified integration point */
	void LoadData(const ElementCardT& element, int ip);

	/*
	 * Returns 1 if the trial elastic strain state lies outside of the
	 * yield surface.
	 */
	int PlasticLoading(const dArrayT& beta, ElementCardT& element,
		int ip);

	double YieldCondition(const dArrayT& devpstress, double alpha) const;

private:

	/* displacements from the last time step */
	const LocalArrayT& fLocLastDisp;

	/* return values */
	dSymMatrixT fb_elastic; //return value
	dMatrixT   fEPModuli;  //elastoplastic moduli in principal stress space

	/* work space */
	dMatrixT fa_inverse; //inverse of 3 x 3 moduli
	dMatrixT fMatrixTemp1;
	dMatrixT fMatrixTemp2;
	dArrayT  fdev_beta; //deviatoric part of principal stress

	dSymMatrixT fb_n;      //last converged elastic state
	dSymMatrixT fb_tr;     //trial elastic state
	dArrayT	    fbeta_tr;  //latest trial principal stresses
	dArrayT     fUnitNorm; //unit normal to principal stress surface
	dArrayT     fInternal; //internal variables

	dMatrixT fFtot;
	dMatrixT ffrel;
	dMatrixT fF_temp;
	dMatrixT fFtot_2D;
	dMatrixT ffrel_2D;	
};

#endif /* _J2_QL_LIN_HARD_2D_T_H_ */
