/* $Id: VWPotentialT.h,v 1.1 2007-04-09 23:33:26 tdnguye Exp $ */
#ifndef _VWPotentialT_
#define _VWPotentialT_

/* base class */
#include "PotentialT.h"

namespace Tahoe {

class VWPotentialT: public PotentialT
{
  public:

	/* constructor */
	VWPotentialT(void);

	/* set parameters */
	void SetKappaMu(double kappa, double mu);
	virtual void DefineParameters(ParameterListT& list) const;
	virtual void TakeParameterList(const ParameterListT& list);

	/*free energy density*/
	virtual double Energy(const dArrayT& lambda_bar,const double& J);

	/*Kirchoff stress measures*/
	virtual void DevStress(const dArrayT& lambda_bar, dArrayT& tau);
	/*derivative of Kirchoff stress with log strain*/
	virtual void DevMod(const dArrayT& lambda_bar,dSymMatrixT& eigenmodulus);

  private:  

  	/*elastic moduli*/
	double falpha;
	double fbeta;
	double fgamma;
};
}
#endif /* _VWPotentialT_ */
