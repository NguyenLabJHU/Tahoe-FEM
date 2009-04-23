/* $Id: VWPotentialT.h,v 1.2 2009-04-23 03:22:46 tdnguye Exp $ */
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
	virtual double Energy(const dArrayT& lambda_bar,const double& J,  double temperature=1.0);

	/*Kirchoff stress measures*/
	virtual void DevStress(const dArrayT& lambda_bar, dArrayT& tau,  double temperature=1.0);
	/*derivative of Kirchoff stress with log strain*/
	virtual void DevMod(const dArrayT& lambda_bar,dSymMatrixT& eigenmodulus,  double temperature=1.0);

  private:  

  	/*elastic moduli*/
	double falpha;
	double fbeta;
	double fgamma;
};
}
#endif /* _VWPotentialT_ */
