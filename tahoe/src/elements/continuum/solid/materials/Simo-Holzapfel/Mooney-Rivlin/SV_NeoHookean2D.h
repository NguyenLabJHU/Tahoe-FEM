/* $Id: SV_NeoHookean2D.h,v 1.5 2003-01-29 07:34:51 paklein Exp $ */
/* created:   TDN (5/31/2001) */

#ifndef _SV_NEOHOOKEAN2D_H_
#define _SV_NEOHOOKEAN2D_H_

#include "FDSimoVisco2D.h"

namespace Tahoe {

class ifstreamT;

/*Compressible Neo-Hookean2D Potential*/
class SV_NeoHookean2D: public FDSimoVisco2D
{
	public:

	/*constructor*/
        SV_NeoHookean2D(ifstreamT& in, const FSMatSupportT& support);
	
	/*print parameters*/
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;

        protected:

	virtual double Phi(const dMatrixT& Fbar, const double& J, 
			   const int SpringType);
		
	virtual void Sig_dev(const dMatrixT& Fbar, const double& J,
		      dSymMatrixT& stress, const int SpringType);
	virtual void Devmod(const dMatrixT& Fbar, const double& J, 
			    dMatrixT& mod, const int SpringType);
	virtual double dUdJ(const double& J, const int SpringType);
	virtual double ddUddJ(const double& J, const int SpringType);

	virtual void OutOfPlaneStretch(const dMatrixT& Fbar,const double& J, 
					const int SpringType);

        protected:

	/*stretch tensor*/
	dSymMatrixT fCbar;
	double fCbar33;

	/*material parameters*/
	dArrayT fKappa;
	dArrayT fMu;
};
}
#endif /*_SV_NEOHOOKEAN2D_H_*/
