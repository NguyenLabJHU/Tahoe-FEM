/* $Id: Hex2D.h,v 1.1.2.1 2003-02-19 01:12:59 paklein Exp $ */
#ifndef _HEX_2D_H_
#define _HEX_2D_H_

/* base class */
#include "NL_E_Mat2DT.h"

namespace Tahoe {

/* forward declarations */
class HexLattice2DT;
class PairPropertyT;

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

protected:

	/** compute the symetric Cij reduced index matrix */
	virtual void ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli);	
	
	/** symmetric 2nd Piola-Kirchhoff stress */
	virtual void ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2);
					                 					
	/** strain energy density */
	virtual double ComputeEnergyDensity(const dSymMatrixT& E);

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
};

} /* namespace Tahoe */

#endif /* _HEX_2D_H_ */
