/* $Id: FCC3D.h,v 1.1.2.1 2003-02-21 01:16:32 paklein Exp $ */
#ifndef _FCC_3D_H_
#define _FCC_3D_H_

/* base class */
#include "NL_E_MatT.h"

namespace Tahoe {

/* forward declarations */
class FCCLatticeT;
class PairPropertyT;

/** plane stress hexagonal lattice */
class FCC3D: public NL_E_MatT
{
public:

	/** constructor */
	FCC3D(ifstreamT& in, const FSMatSupportT& support);
	
	/** destructor */
	~FCC3D(void);
	
	/** \name write parameters */
	/*@{*/
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;
	/*@}*/

protected:

	/** compute the symetric Cij reduced index matrix */
	virtual void ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli);	
	
	/** symmetric 2nd Piola-Kirchhoff stress */
	virtual void ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2);
					                 					
	/** strain energy density */
	virtual double ComputeEnergyDensity(const dSymMatrixT& E);

	/** return the equitriaxial stretch at which the stress is zero */
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

	/** reference volume */
	double fCellVolume;
};

} /* namespace Tahoe */

#endif /* _FCC_3D_H_ */
