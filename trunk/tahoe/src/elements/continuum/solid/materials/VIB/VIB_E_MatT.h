/* $Id: VIB_E_MatT.h,v 1.4 2005-03-16 10:20:42 paklein Exp $ */
/* created: paklein (11/08/1997) */
#ifndef _VIB_E_H_
#define _VIB_E_H_

/* base class */
#include "VIB.h"

namespace Tahoe {

/** base class for isotropic VIB Green elastic materials */
class VIB_E_MatT: public VIB
{
public:

	/* constructor */
	VIB_E_MatT(int nsd);

protected:

	/* set reference energy */
	void SetReferenceEnergy(void);

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

protected:

	/** reference energy */
	double fU_0;
};

} /* namespace Tahoe */

#endif /* _VIB_E_H_ */
