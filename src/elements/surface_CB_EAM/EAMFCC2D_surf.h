/* $Id: EAMFCC2D_surf.h,v 1.1 2006-05-21 15:55:19 hspark Exp $ */
/* created: paklein (12/09/1996) */
#ifndef _EAMFCC2D_SURF_H_
#define _EAMFCC2D_SURF_H_

/* base class */
#include "NL_E_MatT.h"

namespace Tahoe {

/* forward declarations */
class EAMFCC3DSym_surf;

/** plane strain EAM material */
class EAMFCC2D_surf: public NL_E_MatT
{
public:

	/* constructor */
	EAMFCC2D_surf(void);

	/* destructor */
	virtual ~EAMFCC2D_surf(void);
	
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

} // namespace Tahoe 
#endif /* _EAMFCC2D_SURF_H_ */
