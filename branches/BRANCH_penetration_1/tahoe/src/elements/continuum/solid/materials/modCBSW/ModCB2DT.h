/* $Id: ModCB2DT.h,v 1.7 2003-01-29 07:34:59 paklein Exp $ */
/* created: paklein (05/31/1997) */
#ifndef _MODCB_2DT_H_
#define _MODCB_2DT_H_

/* base class */
#include "NL_E_Mat2DT.h"

/* direct members */
#include "SWDataT.h"

namespace Tahoe {

/* forward declarations */
class ModCBSolverT;

class ModCB2DT: public NL_E_Mat2DT
{
public:

	/* plane codes - for crystal axes rotated wrt global axes */
	enum PlaneCodeT {kDC001 = 0,
                     kDC101 = 1,
                     kDC111 = 2};

	/* constructor */
	ModCB2DT(ifstreamT& in, const FSMatSupportT& support, bool equilibrate, 
		PlaneCodeT plane_code);

	/* destructor */
	virtual ~ModCB2DT(void);
	
	/* print Parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;

protected:

	/* symmetric Cij reduced index matrix */
	virtual void ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli);	
	
	/* symmetric 2nd Piola-Kirchhoff stress */
	virtual void ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2);

	/* strain energy density */
	virtual double ComputeEnergyDensity(const dSymMatrixT& E);
	
private:

	/*
	 * Compute the 3D stretch tensor from the 2D reduced index
	 * strain vector (assuming plane strain)
	 */
	void StrainToStretch(const dSymMatrixT& strain2D, dMatrixT& stretch3D);
	
private:
	
	/* 2D plane code */
	PlaneCodeT fPlaneCode;

	/* modified CB solver */
	ModCBSolverT*	fModCBSolver;
	
	/* work space */
	dMatrixT	fCij3D;
	dArrayT		fXsi; //internal DOF vector
	dMatrixT	fStretch3D;
	dMatrixT	fStretch2D;
	dMatrixT	fStress3D;
		
};

} // namespace Tahoe 
#endif /* _MODCB_2DT_H_ */
