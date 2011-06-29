/* $Id: MooneyRivlin.h,v 1.2 2009-04-23 03:22:46 tdnguye Exp $ */
/* created: TDN (01/22/2001) */
#ifndef _MooneyRivlin_
#define _MooneyRivlin_

/* base class */
#include "PotentialT.h"

namespace Tahoe {

class MooneyRivlin: public PotentialT
{
  public:

	/* constructor */
	MooneyRivlin(void);

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

  	/*inelastic moduli*/
	double fc1;
	double fc2;
};
}
#endif /* _RG_MooneyRivlin3D_ */
