/* $Id: Hex2D.h,v 1.3 2004-06-26 06:02:24 paklein Exp $ */
#ifndef _HEX_2D_H_
#define _HEX_2D_H_

/* base class */
#include "NL_E_Mat2DT.h"

namespace Tahoe {

/* forward declarations */
class HexLattice2DT;
class PairPropertyT;
class BondLatticeT;

/** plane stress hexagonal lattice */
class Hex2D: public NL_E_Mat2DT
{
public:

	/** constructor */
	Hex2D(ifstreamT& in, const FSMatSupportT& support);
	
	/** destructor */
	~Hex2D(void);
	
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
	double CellVolume(void) const { return fCellVolume; };

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
	dMatrixT fQ;
	HexLattice2DT* fHexLattice2D;

	/** pair interaction potential */
	PairPropertyT* fPairProperty;

	/** \name work space */
	/*@{*/
	dMatrixT fBondTensor4;
	dArrayT  fBondTensor2;
	/*@}*/

	/** reference volume */
	double fCellVolume;
	
	/** dummy full bond density array */
	dArrayT fFullDensity;
	
	/** flag to indicate whether stress calculation for output should include
	 * the full bond density */
	bool fFullDensityForStressOutput;
};

} /* namespace Tahoe */

#endif /* _HEX_2D_H_ */
