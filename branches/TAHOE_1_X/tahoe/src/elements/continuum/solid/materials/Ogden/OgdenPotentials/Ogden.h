/* $Id: Ogden.h,v 1.1.2.1 2005-02-22 00:15:36 thao Exp $ */
/* created: TDN (01/22/2001) */

#ifndef _OGDEN_
#define _OGDEN_

#include "PotentialT.h"

namespace Tahoe {

class Ogden: public PotentialT
{
  public:

	/* constructor */
	Ogden(ifstreamT& in);

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

	virtual const double Kappa(void);
	virtual const double Mu(void);	

  private:  
  	/* moduli*/
	int fr;
	
	dArrayT falpha_r;
	dArrayT fmu_r;
	double fkappa;
	double fmu;
	
	const double fthird;
};
inline const double Ogden::Kappa(void) {return fkappa;}
inline const double Ogden::Mu(void) {return fmu;}
}
#endif /* _OGDEN_ */
