/* $Id: SV_NeoHookean3D.h,v 1.3 2002-11-14 17:06:15 paklein Exp $ */
/* created: TDN (5/31/2001) */
#ifndef _SV_NEOHOOKEAN3D_H_
#define _SV_NEOHOOKEAN3D_H_

#include "FDSimoVisco3D.h"

namespace Tahoe {

class ifstreamT;

/** compressible Neo-Hookean3D Potential */
class SV_NeoHookean3D: public FDSimoVisco3D
{
	public:

	/*constructor*/
        SV_NeoHookean3D(ifstreamT& in, const FDMatSupportT& support);
	
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

        protected:

	/*stretch tensor*/
	dSymMatrixT fCbar;

	/*material parameters*/
	dArrayT fKappa;
	dArrayT fMu;
};
}
#endif /*_SV_NEOHOOKEAN3D_H_*/
