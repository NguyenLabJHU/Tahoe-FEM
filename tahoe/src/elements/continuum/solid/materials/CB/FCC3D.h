/* $Id: FCC3D.h,v 1.3.12.1 2004-06-09 23:17:30 paklein Exp $ */
#ifndef _FCC_3D_H_
#define _FCC_3D_H_

/* base class */
#include "NL_E_MatT.h"

namespace Tahoe {

/* forward declarations */
class FCCLatticeT;
class PairPropertyT;

/** 3D Cauchy-Born material for FCC crystals with pair potential interactions. */
class FCC3D: public NL_E_MatT
{
public:

	/** constructor */
	FCC3D(ifstreamT& in, const FSMatSupportT& support);
	
	/** destructor */
	~FCC3D(void);
	
protected:

	/** compute the symetric Cij reduced index matrix */
	virtual void ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli);	
	
	/** symmetric 2nd Piola-Kirchhoff stress */
	virtual void ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2);
					                 					
	/** strain energy density */
	virtual double ComputeEnergyDensity(const dSymMatrixT& E);

	/** return the equi-axed stretch at which the stress is zero. This method
	 * assumes the material is isotropic when subject to equi-axed stretch. */
	double ZeroStressStretch(void);

private:

	/** nearest neighbor distance */
	double fNearestNeighbor;

	/** bond information */
	dMatrixT fQ;
	FCCLatticeT* fFCCLattice;

	/** pair interaction potential */
	PairPropertyT* fPairProperty;

	/** \name work space */
	/*@{*/
	dMatrixT fBondTensor4;
	dArrayT  fBondTensor2;
	/*@}*/

	/** atomic volume */
	double fAtomicVolume;
};

} /* namespace Tahoe */

#endif /* _FCC_3D_H_ */
