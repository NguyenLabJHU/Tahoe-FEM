/* $Id: RG_NeoHookean2D.h,v 1.2 2002-11-14 17:06:10 paklein Exp $ */
/* created: TDN (01/22/2001) */

#ifndef _RG_NeoHookean2D_
#define _RG_NeoHookean2D_

#include "RG_VDSplit2D.h"

namespace Tahoe {

class RG_NeoHookean2D: public RG_VDSplit2D
{
  public:

	/* constructor */
	RG_NeoHookean2D(ifstreamT& in, const FDMatSupportT& support);

	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;
	virtual void Initialize(void);
  protected:
	virtual double Phi(const dArrayT& eigenstretch_bar, const double& J,
			   const int SpringType);
	virtual void devTau(const dArrayT& eigenstretch_bar, 
			    dArrayT& eigenstress, const int SpringType);
	virtual double meanTau(const double& J, const int SpringType);
	virtual void DdevTauDepsilon(const dArrayT& eigenstretch_bar, 
				     dSymMatrixT&  eigenmodulus, 
				     const int SpringType);
	virtual double DmeanTauDepsilon(const double& J, const int SpringType);
	virtual void OutOfPlaneStretch(const dArrayT& eigenstretch_bar, 
				       const double& J, const int SpringType);

  private:  
  	/*inelastic moduli*/
	dArrayT fMu;
	dArrayT fKappa;

	/*out of plane strain*/
	double fl2_bar;
};
}
#endif /* _RG_NeoHookean2D_ */
