/* $Id: Chain1D.h,v 1.2 2004-06-26 05:56:41 paklein Exp $ */
#ifndef _CHAIN_1D_H_
#define _CHAIN_1D_H_

/* base class */
#include "NL_E_MatT.h"

namespace Tahoe {

/* forward declarations */
class Lattice1DT;
class PairPropertyT;
class BondLatticeT;

/** 1D Cauchy-Born material with pair potential interactions. */
class Chain1D: public NL_E_MatT
{
public:

	/** constructor */
	Chain1D(ifstreamT& in, const FSMatSupportT& support);
	
	/** destructor */
	~Chain1D(void);
	
	/** \name write parameters */
	/*@{*/
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;
	/*@}*/

	/** \name Cauchy-Born parameters */
	/*@{*/
	/** return a reference to the bond lattice */
	const BondLatticeT& BondLattice(void) const;

	/** reference volume */
	double CellVolume(void) const { return fAtomicVolume; };

	/** nearest neighbor distance */
	double NearestNeighbor(void) const { return fNearestNeighbor; };
	/*@}*/

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
	Lattice1DT* fLattice1D;

	/** pair interaction potential */
	PairPropertyT* fPairProperty;

	/** \name work space */
	/*@{*/
	dMatrixT fBondTensor4;
	dArrayT  fBondTensor2;
	/*@}*/

	/** atomic volume */
	double fAtomicVolume;

	/** dummy full bond density array */
	dArrayT fFullDensity;

	/** flag to indicate whether stress calculation for output should include
	 * the full bond density */
	bool fFullDensityForStressOutput;
};

} /* namespace Tahoe */

#endif /* _FCC_3D_H_ */
