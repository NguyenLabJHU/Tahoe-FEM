/* $Id: VerondaWestmannT.h,v 1.1 2006-05-02 00:58:53 thao Exp $ */
#ifndef _V_W_
#define _V_W_

/* base class */
#include "OgdenIsotropicT.h"

namespace Tahoe {

/** principal stretch version of Quad Log model */
class VerondaWestmannT: public OgdenIsotropicT
{
public:

	/** constructor */
	VerondaWestmannT(void);
	
	/** strain energy density */
	virtual double StrainEnergyDensity(void);

protected:

	/* principal values given principal values of the stretch tensors,
	 * i.e., the principal stretches squared */
	virtual void dWdE(const dArrayT& eigenstretch2, dArrayT& eigenstress);
	virtual void ddWddE(const dArrayT& eigenstretch2, dArrayT& eigenstress,
		dSymMatrixT& eigenmod);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

private:
	/*derivatives of strain energy with respect to stretch invariants*/
	void Compute_Ik(const dArrayT& eigenstretch2, dArrayT& Ik);
	void Compute_Ik(dArrayT& Ik);
	void Compute_dWdIk (const dArrayT& Ik, dArrayT& dWdIk);
	void Compute_ddWddIk (const dArrayT& Ik, dArrayT& dWdIk, dSymMatrixT& ddWddIk);

private:
	/*material constants*/
	double falpha;
	double fbeta;
	double fgamma;
	double fkappa;

	/*invariants of stretch tensor*/
	dArrayT fIk;
	
	/*derivatives of strain energy density with respect to stretch tensor*/
	dArrayT fdWdIk;
	dSymMatrixT fddWddIk;

	/*work space*/
	dSymMatrixT fMat; 
};

} // namespace Tahoe 
#endif /* _V_W_ */
