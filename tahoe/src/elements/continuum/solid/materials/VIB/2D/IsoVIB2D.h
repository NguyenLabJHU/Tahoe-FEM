/* $Id: IsoVIB2D.h,v 1.1.1.1 2001-01-29 08:20:24 paklein Exp $ */
/* created: paklein (11/08/1997)                                          */
/* 2D Isotropic VIB solver using spectral decomposition formulation       */

#ifndef _ISO_VIB_2D_H_
#define _ISO_VIB_2D_H_

/* base classes */
#include "FDStructMatT.h"
#include "Material2DT.h"
#include "VIB.h"

/* direct members */
#include "SpectralDecompT.h"

/* forward declarations */
class CirclePointsT;

class IsoVIB2D: public FDStructMatT, public Material2DT, public VIB
{
public:

	/* constructor */
	IsoVIB2D(ifstreamT& in, const ElasticT& element);

	/* destructor */
	~IsoVIB2D(void);
	
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

	//TEMP
	const dSymMatrixT& CurvatureTensor(void);

protected:

	/* strained lengths in terms of the Lagrangian stretch eigenvalues */
	void ComputeLengths(const dArrayT& eigs);

private:

	/* initialize angle tables */
	void Construct(void);

protected:
	
	/* integration point generator */
	CirclePointsT*	fCircle;
	  	
	/* eigenvalues */
	dArrayT	    fEigs;    // rank2 eigenvalues
	dSymMatrixT	fEigmods; // rank4 eigenvalues
	
	/* spectral decomp solver */
	SpectralDecompT	fSpectral;
	
private:

	/* return value */
	dMatrixT fModulus;
};

#endif /* _ISO_VIB_2D_H_ */
