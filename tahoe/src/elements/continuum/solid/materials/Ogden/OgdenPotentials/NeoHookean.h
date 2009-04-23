/* $Id: NeoHookean.h,v 1.5 2009-04-23 03:22:46 tdnguye Exp $ */
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
	virtual double Energy(const dArrayT& lambda_bar,const double& J,  double temperature=1.0);

	/*Kirchoff stress measures*/
	virtual void DevStress(const dArrayT& lambda_bar, dArrayT& tau,  double temperature=1.0);
	/*derivative of Kirchoff stress with log strain*/
	virtual void DevMod(const dArrayT& lambda_bar,dSymMatrixT& eigenmodulus,  double temperature=1.0);

	protected:
		void SetMu(double mu);
		double GetMu (void);
  private:  

};
}
#endif /* _RG_NeoHookean3D_ */
