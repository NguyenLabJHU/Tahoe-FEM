/* $Id: NeoHookean.h,v 1.1.50.1 2005-02-22 00:15:36 thao Exp $ */
/* created: TDN (01/22/2001) */

#ifndef _NeoHookean_
#define _NeoHookean_

#include "PotentialT.h"

namespace Tahoe {

class NeoHookean: public PotentialT
{
  public:

	/* constructor */
	NeoHookean(ifstreamT& in);

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
  	/*inelastic moduli*/
	double fMu;
	double fKappa;
	
	const double fthird;
};
inline const double NeoHookean::Kappa(void) {return fKappa;}
inline const double NeoHookean::Mu(void) {return fMu;}

}
#endif /* _RG_NeoHookean3D_ */
