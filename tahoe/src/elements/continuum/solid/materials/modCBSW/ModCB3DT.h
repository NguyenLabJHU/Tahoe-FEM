/* $Id: ModCB3DT.h,v 1.1.1.1 2001-01-29 08:20:26 paklein Exp $ */
/* created: paklein (10/14/1998)                                          */

#ifndef _MODCB_3D_T_H_
#define _MODCB_3D_T_H_

/* base class */
#include "NL_E_MatT.h"

/* direct members */
#include "SWDataT.h"

/* forward declarations */
class ModCBSolverT;

class ModCB3DT: public NL_E_MatT
{
public:

	/* constructor */
	ModCB3DT(ifstreamT& in, const ElasticT& element, bool equilibrate);

	/* destructor */
	virtual ~ModCB3DT(void);
	
	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;

protected:

	/* compute the symmetric Cij reduced index matrix */
	virtual void ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli);	
	
	/* compute the symmetric 2nd Piola-Kirchhoff reduced index vector */
	virtual void ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2);

	/* returns the strain energy density for the specified strain */
	virtual double ComputeEnergyDensity(const dSymMatrixT& E);
	
private:

	/* compute the 3D stretch tensor from the 2D reduced index
	 * strain vector (assuming plane strain) */
	void StrainToStretch(const dSymMatrixT& E, dMatrixT& C);
	
private:
	
	/* orientation code */
	int fOrientationCode;

	/* modified CB solver */
	ModCBSolverT* fModCBSolver;
	
	/* work space */
	dArrayT	 fXsi; //internal DOF vector
	dMatrixT fC;	
	dMatrixT fPK2;		
};

#endif /* _MODCB_3D_T_H_ */
