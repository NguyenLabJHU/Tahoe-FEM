/* $Id: EAMFCC3DMatT_surf.h,v 1.2 2006-05-29 17:22:56 paklein Exp $ */
/* created: paklein (10/25/1998) */
#ifndef _EAMFCC3DMatT_SURF_H_
#define _EAMFCC3DMatT_SURF_H_

/* base class */
#include "NL_E_MatT.h"

namespace Tahoe {

/* forward declarations */
class EAMFCC3DSym_surf;

/** plane strain EAM material */
class EAMFCC3DMatT_surf: public NL_E_MatT
{
public:

	/* constructor */
	EAMFCC3DMatT_surf(void);

	/* destructor */
	virtual ~EAMFCC3DMatT_surf(void);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;
	
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/* compute the symetric Cij reduced index matrix */
	virtual void ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli);	
	
	/* compute the symetric 2nd Piola-Kirchhoff reduced index vector */
	virtual void ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2);

	/* returns the strain energy density for the specified strain */
	virtual double ComputeEnergyDensity(const dSymMatrixT& E);
	
protected:
	
	/** Cauchy-Born EAM solver */
	EAMFCC3DSym_surf* fEAM;
};

} /* namespace Tahoe */

#endif /* _EAMFCC3DMatT_SURF_H_ */
