/* $Id: IsoVIB3D.h,v 1.4 2002-07-05 22:28:18 paklein Exp $ */
/* created: paklein (03/15/1998)                                          */
/* 3D Isotropic VIB solver using spectral decomposition formulation       */

#ifndef _ISO_VIB_3D_H_
#define _ISO_VIB_3D_H_

/* base classes */
#include "FDStructMatT.h"
#include "VIB.h"

/* direct members */
#include "SpectralDecompT.h"

namespace Tahoe {

/* forward declarations */
class SpherePointsT;

class IsoVIB3D: public FDStructMatT, public VIB
{
public:

	/* constructor */
	IsoVIB3D(ifstreamT& in, const FiniteStrainT& element);

	/* destructor */
	virtual ~IsoVIB3D(void);
	
	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;	
	
	/* spatial description */
	virtual const dMatrixT& c_ijkl(void); // spatial tangent moduli
	virtual const dSymMatrixT& s_ij(void); // Cauchy stress

	/* material description */
	virtual const dMatrixT& C_IJKL(void); // material tangent moduli
	virtual const dSymMatrixT& S_IJ(void); // PK2 stress
//TEMP - not yet optimized for total Lagrangian formulation.
//       calls to these write error message and throw exception

	/* strain energy density */
	virtual double StrainEnergyDensity(void);

protected:

	/* strained lengths in terms of the Lagrangian stretch eigenvalues */
	void ComputeLengths(const dArrayT& eigs);

private:

	/* initialize angle tables */
	void Construct(void);

protected:

	/* spectral decomp solver */
	SpectralDecompT	fSpectral;

	/* eigenvalues */
	dArrayT	    fEigs;		//rank2 eigenvalues
	dSymMatrixT	fEigmods;	//rank4 eigenvalues  	
	
private:

	/* stretch */
	dSymMatrixT fb;

	/* integration point generator */
	SpherePointsT*	fSphere;
	
	/* return values */
	dMatrixT    fModulus;
	dSymMatrixT fStress;
};

} // namespace Tahoe 
#endif /* _ISO_VIB_3D_H_ */
