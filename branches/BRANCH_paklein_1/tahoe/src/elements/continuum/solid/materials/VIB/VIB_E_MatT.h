/* $Id: VIB_E_MatT.h,v 1.2 2002-07-02 19:55:55 cjkimme Exp $ */
/* created: paklein (11/08/1997)                                          */
/* Base class for isotropic VIB_E_MatT solvers.                           */

#ifndef _VIB_E_H_
#define _VIB_E_H_

/* base class */
#include "VIB.h"


namespace Tahoe {

class VIB_E_MatT: public VIB
{
public:

	/* constructor */
	VIB_E_MatT(ifstreamT& in, int nsd);

	/* print parameters */
	virtual void PrintName(ostream& out) const;	

protected:

	/* returns the strain energy density for the specified strain */
	double VIBEnergyDensity(const dSymMatrixT& E);

	/* compute strained lengths */
	void ComputeLengths(const dSymMatrixT& strain);

	/* convenience */
	void SetStressPointers2D(double*&,double*&,double*&);
	void SetStressPointers3D(double*&,double*&,double*&,
	                         double*&,double*&,double*&);

	void SetModuliPointers2D(double*&, double*&, double*&,
							 double*&, double*&);
	void SetModuliPointers3D(double*&, double*&, double*&, double*&, double*&,
	                         double*&, double*&, double*&, double*&, double*&,
	                         double*&, double*&, double*&, double*&, double*&);

};

} // namespace Tahoe 
#endif /* _VIB_E_H_ */
