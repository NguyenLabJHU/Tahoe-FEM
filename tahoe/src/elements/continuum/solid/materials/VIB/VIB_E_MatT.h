/* $Id: VIB_E_MatT.h,v 1.1.1.1 2001-01-29 08:20:24 paklein Exp $ */
/* created: paklein (11/08/1997)                                          */
/* Base class for isotropic VIB_E_MatT solvers.                           */

#ifndef _VIB_E_H_
#define _VIB_E_H_

/* base class */
#include "VIB.h"

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

#endif /* _VIB_E_H_ */
