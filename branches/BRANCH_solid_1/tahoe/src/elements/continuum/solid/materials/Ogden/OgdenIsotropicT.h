/* $Id: OgdenIsotropicT.h,v 1.2.2.1 2001-06-06 16:24:14 paklein Exp $ */
/* created: paklein (10/01/2000)                                          */
/* base class for large deformation isotropic material following          */
/* Ogden's formulation.                                                   */

#ifndef _OGDEN_ISOTROPIC_T_H_
#define _OGDEN_ISOTROPIC_T_H_

/* base classes */
#include "FDStructMatT.h"
#include "IsotropicT.h"

/* direct members */
#include "SpectralDecompT.h"

class OgdenIsotropicT: public FDStructMatT, public IsotropicT
{
public:

	/* constructor */
	OgdenIsotropicT(ifstreamT& in, const ElasticT& element);
	
	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;

	/* class specific initializations */
	virtual void Initialize(void);

	/* spatial description */
	virtual const dMatrixT& c_ijkl(void);
	virtual const dSymMatrixT& s_ij(void);

	/* material description */
	virtual const dMatrixT& C_IJKL(void); // material tangent moduli
	virtual const dSymMatrixT& S_IJ(void); // PK2 stress

protected:

	/* principal values given principal values of the stretch tensors,
	 * i.e., the principal stretches squared */
	virtual void dWdE(const dArrayT& eigenstretch2, dArrayT& eigenstress) = 0;
	virtual void ddWddE(const dArrayT& eigenstretch2, dArrayT& eigenstress,
		dSymMatrixT& eigenmod) = 0;

	/* return true of model is purely 2D, plain stress */
	virtual bool PurePlainStress(void) const { return false; };

private:

	/* construct symmetric rank-4 mixed-direction tensor (6.1.44) */
	void MixedRank4_2D(const dArrayT& a, const dArrayT& b,
		dMatrixT& rank4_ab) const;
	void MixedRank4_3D(const dArrayT& a, const dArrayT& b,
		dMatrixT& rank4_ab) const;

protected:

	/* spectral operations */
	SpectralDecompT fSpectralDecomp;

	/* work space */
	dArrayT     fEigs; //TEMP - need this??
	dArrayT     fdWdE;
	dSymMatrixT fddWddE;
	dMatrixT    fModMat;
	
	/* return values */
	dMatrixT    fModulus;
	dSymMatrixT fStress;
};

#endif /* _OGDEN_ISOTROPIC_T_H_ */
