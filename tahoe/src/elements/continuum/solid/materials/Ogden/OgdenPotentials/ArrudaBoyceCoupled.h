/* $Id: ArrudaBoyceCoupled.h,v 1.3 2021-01-10 Zheliang Wang $ */
/* created: Zheliang Wang (2021 Jan) */
#ifndef _ArrudaBoyceCoupled_
#define _ArrudaBoyceCoupled_


/* base class */
#include "PotentialT.h"
#include "InvLangevin.h"

namespace Tahoe {

class ArrudaBoyceCoupled: public PotentialT
{
  public:

	/* constructor */
	ArrudaBoyceCoupled(void);

	/* set parameters */
	void SetKappaMu(double kappa, double mu);

	virtual void DefineParameters(ParameterListT& list) const;
	virtual void TakeParameterList(const ParameterListT& list);

	/*free energy density*/
	virtual double Energy(const dArrayT& lambda_bar,const double& J,  double temperature=1.0);

	/*Kirchoff stress measures*/
    /*This is the coupled form, although the name is deviatoric stress, actually it's not*/
	virtual void DevStress(const dArrayT& lambda_bar, dArrayT& tau,  double temperature=1.0);
    
	/*derivative of Kirchoff stress with log strain*/
	virtual void DevMod(const dArrayT& lambda_bar,dSymMatrixT& eigenmodulus,  double temperature=1.0);
    
    
	protected:
		double GetMu (void);
		
  private:  
	InvLangevin fLangevin; /*Inverse langevin function for strain hardening*/
	double fmuN;		/*network stiffness*/
	double flambdaL;	/*locking stretch*/
	double fT0;			/*reference temperature*/

};
}
#endif /* _RG_ArrudaBoyceCoupled3D_ */
