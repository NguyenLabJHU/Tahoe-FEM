/* $Id: J2QLLinHardT.h,v 1.6 2003-01-29 07:35:02 paklein Exp $ */
/* created: paklein (10/26/2000) */
#ifndef _J2_QL_LIN_HARD_T_H_
#define _J2_QL_LIN_HARD_T_H_

/* base classes */
#include "QuadLog3D.h"
#include "J2PrimitiveT.h"

/* direct members */
#include "dSymMatrixT.h"
#include "dMatrixT.h"
#include "dArrayT.h"
#include "iArrayT.h"

namespace Tahoe {

/** interface for a elastoplastic material that is linearly
 * isotropically elastic subject to the Huber-von Mises yield
 * condition as fYield with kinematic/isotropic hardening laws
 * given by:
 *      H(a) = (1 - ftheta) fH_bar a
 *      K(a) = fYield + ftheta fH_bar a
 * 		where a is the internal hardening variable
 */
class J2QLLinHardT: public QuadLog3D, public J2PrimitiveT
{
public:

	/* constructor */
	J2QLLinHardT(ifstreamT& in, const FSMatSupportT& support);

	/** required parameter flags */
	virtual bool Need_F_last(void) const { return true; };

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
	void ElastoPlasticCorrection(dSymMatrixT& a_ep, dArrayT& beta, int ip);

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

	/* return values */
	dSymMatrixT fb_elastic; //return value
	dSymMatrixT fEPModuli;  //elastoplastic moduli in principal stress space

	/* work space */
	dSymMatrixT fa_inverse; //inverse of 3 x 3 moduli
	dMatrixT    fMatrixTemp1;
	dMatrixT    fMatrixTemp2;
	dSymMatrixT fMatrixTemp3;
	dArrayT     fdev_beta; //deviatoric part of principal stress

	dSymMatrixT fb_n;      //last converged elastic state
	dSymMatrixT fb_tr;     //trial elastic state
	dArrayT	    fbeta_tr;  //latest trial principal stresses
	dArrayT     fUnitNorm; //unit normal to principal stress surface
	dArrayT     fInternal; //internal variables

	dMatrixT fFtot;
	dMatrixT ffrel;
	dMatrixT fF_temp;
};

} // namespace Tahoe 
#endif /* _J2_QL_LIN_HARD_T_H_ */
