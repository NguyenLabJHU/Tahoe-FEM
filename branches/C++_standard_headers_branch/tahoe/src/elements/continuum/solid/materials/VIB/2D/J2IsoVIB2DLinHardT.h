/* $Id: J2IsoVIB2DLinHardT.h,v 1.7 2004-07-15 08:27:45 paklein Exp $ */
/* created: paklein (10/18/1998) */
#ifndef _J2_ISOVIB_2D_T_H_
#define _J2_ISOVIB_2D_T_H_

/* base classes */
#include "IsoVIB3D.h"
#include "J2PrimitiveT.h"

/* direct members */
#include "dSymMatrixT.h"
#include "dMatrixT.h"
#include "dArrayT.h"
#include "iArrayT.h"

namespace Tahoe {

/** VIB plus principal stretch elasticity
 * Interface for a elastoplastic material that is linearly
 * isotropically elastic subject to the Huber-von Mises yield
 * condition as fYield with kinematic/isotropic hardening laws
 * given by:
 *      H(a) = (1 - ftheta) fH_bar a
 *      K(a) = fYield + ftheta fH_bar a
 * 		where a is the internal hardening variable
 * 	\note all calculations are peformed in 3D
 */
class J2IsoVIB2DLinHardT: public IsoVIB3D, public J2PrimitiveT
{
public:

	/* constructor */
	J2IsoVIB2DLinHardT(void);

	/* update internal variables */
	virtual void UpdateHistory(void);

	/* reset internal variables to last converged solution */
	virtual void ResetHistory(void);

	/* spatial description */
	virtual const dMatrixT& c_ijkl(void);
	virtual const dSymMatrixT& s_ij(void);

	/* material description */
	virtual const dMatrixT& C_IJKL(void); // material tangent moduli
	virtual const dSymMatrixT& S_IJ(void); // PK2 stress
//TEMP - not yet optimized for total Lagrangian formulation.
//       calls to these write error message and throw ExceptionT::xception

	/* strain energy density */
	virtual double StrainEnergyDensity(void);

	/* required parameter flags */
	virtual bool NeedLastDisp(void) const;

	/* returns the number of internal variables */
	virtual int NumOutputVariables(void) const;
	virtual void OutputLabels(ArrayT<StringT>& labels) const;
	virtual void ComputeOutput(dArrayT& output);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/
	
protected:

	/* returns the trial stretch */
	const dSymMatrixT& TrialStretch(const dMatrixT& F_total,
		const dMatrixT& f_relative, int ip);
			
	/* apply return mapping if necessary */
	void ReturnMapping(const dSymMatrixT& b_tr, const dArrayT& b_eigs,
		dArrayT& beta, int ip);

	/* return a pointer to a new plastic element object constructed with
	 * the data from element */
	void AllocateElement(ElementCardT& element);

private:

	/* compute F_total and f_relative */
	void ComputeGradients(void);

	/* initialize intermediate state from F_n (for ) */
	void InitIntermediate(const dMatrixT& F_total, const dMatrixT& f_relative);

	/* load element data for the specified integration point */
	void LoadData(const ElementCardT& element, int ip);

	/* returns 1 if the trial elastic strain state lies outside of the
	 * yield surface */
	int PlasticLoading(const dArrayT& beta, ElementCardT& element,
		int ip);

	double YieldCondition(const dArrayT& devpstress, double alpha) const;

	/* integrate Kirchhoff stress principal values */
	void ComputeBeta(const dArrayT& eigs, dArrayT& beta);
	void Computeddw(const dArrayT& eigs, dMatrixT& ddw);

private:

//TEMP - overrides IsoVIB3D::fEigs
	dArrayT    fEigs;
	
	/* isotropic stored energy function derivatives wrt to log stretch */
	dArrayT    fBeta;
	dMatrixT   fddW;

	/* return values */
	dSymMatrixT fb_elastic; //return value
	dMatrixT   fEPModuli;  //elastoplastic moduli in principal stress space

	/* work space */
	dMatrixT fMatrixTemp1;
	dMatrixT fMatrixTemp2;
	dArrayT  fdev_beta;   //deviatoric part of principal stress

	dSymMatrixT fb_n;      // last converged elastic state
	dSymMatrixT fb_tr;     // trial elastic state
	dArrayT	   fbeta_tr;  // latest trial principal stresses
	dArrayT    flog_e;    // latest log elastic stretches
	dArrayT    fUnitNorm; // unit normal to principal stress surface
	dArrayT    fInternal; // internal variables

	dMatrixT fFtot;
	dMatrixT ffrel;

	/* 2D translation */
	dMatrixT fF_temp;
	dMatrixT fFtot_2D;
	dMatrixT ffrel_2D;
	
	/* return values */
	dSymMatrixT fStress2D;
	dSymMatrixT fb_2D;	
};

} // namespace Tahoe 
#endif /* _J2_ISOVIB_2D_T_H_ */
