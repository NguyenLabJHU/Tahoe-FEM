/* $Id: VIB_E_MatT.h,v 1.3 2004-07-15 08:27:40 paklein Exp $ */
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

} /* namespace Tahoe */

#endif /* _VIB_E_H_ */
