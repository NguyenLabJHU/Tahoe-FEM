/* $Id: QuadLog3D.h,v 1.3.2.1 2001-06-07 03:01:20 paklein Exp $ */
/* created: paklein (06/27/1997)                                          */
/* Hyperelastic material governed by quadratic logarithmic potential.     */

#ifndef _QUAD_LOG_3D_H_
#define _QUAD_LOG_3D_H_

/* base classes */
#include "FDStructMatT.h"
#include "IsotropicT.h"
#include "SpectralDecompT.h"

class QuadLog3D: public FDStructMatT, public IsotropicT
{
public:

	/* constructor */
	QuadLog3D(ifstreamT& in, const ElasticT& element);
	
	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;

	/* modulus */
	virtual const dMatrixT& c_ijkl(void);
	
	/* stresses */
	virtual const dSymMatrixT& s_ij(void);

	/* material description */
	virtual const dMatrixT& C_IJKL(void); // material tangent moduli
	virtual const dSymMatrixT& S_IJ(void); // PK2 stress
//TEMP - no reason to use these in total Lagrangian formulation.
//       calls to these write error message and throw exception

	/* strain energy density */
	virtual double StrainEnergyDensity(void);

protected:

	/* computation routines */
	void ComputeModuli(const dSymMatrixT& b, dMatrixT& moduli);
	void ComputeCauchy(const dSymMatrixT& b, dSymMatrixT& cauchy);
	double ComputeEnergy(const dArrayT& loge);

	/* compute logarithmic stretches from the given eigenvalues */
	void LogStretches(const dArrayT& eigs);

protected:

	/* spectral decomposition solver */
	SpectralDecompT fSpectral;

	/* left stretch */
	dSymMatrixT fb;

	/* return values */
	dSymMatrixT	fStress;
	dMatrixT	fModulus;

	/* deviatoric operator */
	dSymMatrixT fDevOp3;

	/* spectral decomposition */
	dArrayT	fEigs;  //principal value of b
	dArrayT	floge;  //logarithmic stretches
	dArrayT	fBeta;  //principal stresses
	dSymMatrixT fEigMod;//modulus in principal stretches
};

#endif /* _QUAD_LOG_3D_H_ */
