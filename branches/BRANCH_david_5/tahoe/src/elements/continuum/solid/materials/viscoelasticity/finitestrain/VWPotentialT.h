/* $Id: VWPotentialT.h,v 1.1 2005-07-14 19:44:38 regueiro Exp $ */
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

	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;
	virtual void Initialize(void);

	/*free energy density*/
	virtual double Energy(const dArrayT& lambda_bar,const double& J);

	/*Kirchoff stress measures*/
	virtual void DevStress(const dArrayT& lambda_bar, dArrayT& tau);
	virtual double MeanStress(const double& J);

	/*derivative of Kirchoff stress with log strain*/
	virtual void DevMod(const dArrayT& lambda_bar,dSymMatrixT& eigenmodulus);
	virtual double MeanMod(const double& J);

  private:  

  	/*elastic moduli*/
	double fMu;
	double fKappa;
};
}
#endif /* _VWPotentialT_ */
