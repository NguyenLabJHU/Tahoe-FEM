/* $Id: SimoIso3D.h,v 1.1.1.1 2001-01-29 08:20:25 paklein Exp $ */
/* created: paklein (03/02/1997)                                          */
/* Hyperelastic material governed by Simo's split volumetric/deviatoric   */
/* stored energy function.                                                */
/* Note: This material is inherently 3D.                                  */

#ifndef _SIMO_ISO_3D_H_
#define _SIMO_ISO_3D_H_

/* base classes */
#include "FDStructMatT.h"
#include "IsotropicT.h"

class SimoIso3D: public FDStructMatT, public IsotropicT
{
public:

	/* constructor */
	SimoIso3D(ifstreamT& in, const ElasticT& element);
	
	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;

	/* spatial description */
	virtual const dMatrixT& c_ijkl(void); // spatial tangent moduli
	virtual const dSymMatrixT& s_ij(void); // Cauchy stress

	/* material description */
	virtual const dMatrixT& C_IJKL(void); // material tangent moduli
	virtual const dSymMatrixT& S_IJ(void); // PK2 stress
//TEMP - no reason to use these in total Lagrangian formulation.
//       calls to these write error message and throw exception

	/* strain energy density */
	virtual double StrainEnergyDensity(void);
	
protected:
	
	/* computation routines - split volumetric/deviatoric */
	void ComputeModuli(double J, const dSymMatrixT& b_bar, dMatrixT& moduli);
	void ComputeCauchy(double J, const dSymMatrixT& b_bar, dSymMatrixT& cauchy);
	double ComputeEnergy(double J, const dSymMatrixT& b_bar);

	/* Volumetric energy function and derivatives */
	double   U(double J) const;
	double  dU(double J) const;
	double ddU(double J) const;

protected:

	/* return value */
	dSymMatrixT	fStress;
	dMatrixT    fModulus;

private:

	/* moduli */
	double fmu; 	//shear modulus
	double fkappa;	//bulk modulus

	/* work space */
	dSymMatrixT	fb_bar;
	dSymMatrixT	fnorm;
	dMatrixT	frank4;
	
	/* fixed forms */
	dSymMatrixT	fIdentity;
	dMatrixT	fIcrossI;
	dMatrixT	fIdentity4;
	dMatrixT	fDevOp4;
};

#endif /* _SIMO_ISO_3D_H_ */
