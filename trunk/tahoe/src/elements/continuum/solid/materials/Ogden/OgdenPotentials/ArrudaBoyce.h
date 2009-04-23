/* $Id: ArrudaBoyce.h,v 1.2 2009-04-23 03:22:46 tdnguye Exp $ */
/* created: TDN (01/22/2001) */
#ifndef _ArrudaBoyce_
#define _ArrudaBoyce_


/* base class */
#include "PotentialT.h"
#include "InvLangevin.h"

namespace Tahoe {

class ArrudaBoyce: public PotentialT
{
  public:

	/* constructor */
	ArrudaBoyce(void);

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
		double GetMu (void);
		
  private:  
	InvLangevin fLangevin; /*Inverse langevin function for strain hardening*/
	double fmuN;		/*network stiffness*/
	double flambdaL;	/*locking stretch*/

};
}
#endif /* _RG_ArrudaBoyce3D_ */
