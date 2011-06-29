/* $Id: ModCB2DT.h,v 1.9 2004-09-10 22:39:22 paklein Exp $ */
/* created: paklein (05/31/1997) */
#ifndef _MODCB_2DT_H_
#define _MODCB_2DT_H_

/* base class */
#include "NL_E_MatT.h"

namespace Tahoe {

/* forward declarations */
class ModCBSolverT;

class ModCB2DT: public NL_E_MatT
{
public:

	/* constructor */
	ModCB2DT(void);

	/* destructor */
	virtual ~ModCB2DT(void);
	
	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
 	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

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
