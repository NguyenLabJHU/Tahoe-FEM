/* $Id: MooneyRivlin.h,v 1.1 2007-04-09 23:33:26 tdnguye Exp $ */
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
	virtual double Energy(const dArrayT& lambda_bar,const double& J);

	/*Kirchoff stress measures*/
	virtual void DevStress(const dArrayT& lambda_bar, dArrayT& tau);
	/*derivative of Kirchoff stress with log strain*/
	virtual void DevMod(const dArrayT& lambda_bar,dSymMatrixT& eigenmodulus);

  private:  

  	/*inelastic moduli*/
	double fc1;
	double fc2;
};
}
#endif /* _RG_MooneyRivlin3D_ */
