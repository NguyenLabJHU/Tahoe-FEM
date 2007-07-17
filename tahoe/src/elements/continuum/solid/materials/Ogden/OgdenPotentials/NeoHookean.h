/* $Id: NeoHookean.h,v 1.4 2007-07-17 20:12:30 tdnguye Exp $ */
/* created: TDN (01/22/2001) */
#ifndef _NeoHookean_
#define _NeoHookean_

/* base class */
#include "PotentialT.h"

namespace Tahoe {

class NeoHookean: public PotentialT
{
  public:

	/* constructor */
	NeoHookean(void);

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

	protected:
		void SetMu(double mu);
		double GetMu (void);
  private:  

};
}
#endif /* _RG_NeoHookean3D_ */
