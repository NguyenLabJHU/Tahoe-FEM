/* $Id: QuadLog3D.h,v 1.1.1.1 2001-01-29 08:20:25 paklein Exp $ */
/* created: paklein (06/27/1997)                                          */
/* Hyperelastic material governed by quadratic logarithmic potential.     */

#ifndef _QUAD_LOG_3D_H_
#define _QUAD_LOG_3D_H_

/* base classes */
#include "FDStructMatT.h"
#include "IsotropicT.h"

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

	/* compute spectral decomposition of b
*
* NOTE: Repeated eigenvalues are perturbed
	 */
	void SpectralDecomp(const dSymMatrixT& b, dArrayT& eigs, bool perturb_repeated);

	/* compute logarithmic stretches from the given eigenvalues */
	void LogStretches(dArrayT& eigs);

	/* components of algorithmic tangent moduli */
	void Set_b_Tensor(const dSymMatrixT& b);
	const dMatrixT& SpatialTensor(const dSymMatrixT& b, int A);
	
private:

	/* returns min */
	static double Min(double d1, double d2, double d3);

protected:

	/* return value */
	dSymMatrixT	fStress;
	dMatrixT	fModulus;

/* elastic constants */
double fmu;
double flambda;
	
	/* fixed forms */
	dSymMatrixT	fISym;
	dMatrixT	fIdentity3;
	dMatrixT	f1x1;
	dMatrixT	fDevOp3;

	/* spectral decomposition */
	dArrayT	 fEigs;  //principal value of b
	dArrayT	 floge;  //logarithmic stretches
	dArrayT	 fBeta;  //principal stresses
	dMatrixT fEigMod;//modulus in principal stretches

	ArrayT<dSymMatrixT> fm; //array of rank 1 matrices	
	dSymMatrixT& fn0xn0; //rank 1 principal direction matricies
	dSymMatrixT& fn1xn1;
	dSymMatrixT& fn2xn2;
	
	/* decomp work space */
	dSymMatrixT fm1;
	dSymMatrixT fm2;
	
	/* spatial tensor work space */
	dMatrixT   fSpatTensor;
	dMatrixT   fc_b; //part of spatial tensor dependent only on b
	dMatrixT   fRank4;
	dSymMatrixT fRank2;
	dMatrixT   fIdentity4;
};

/* inline functions */

/* returns max */
inline double QuadLog3D::Min(double d1, double d2, double d3)
{
	return ( (d1 < d2) ?
			((d1 < d3) ? d1 : d3) :
			((d2 < d3) ? d2 : d3) );
}

#endif /* _QUAD_LOG_3D_H_ */
