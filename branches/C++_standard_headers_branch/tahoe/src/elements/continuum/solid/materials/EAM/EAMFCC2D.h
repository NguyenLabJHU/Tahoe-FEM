/* $Id: EAMFCC2D.h,v 1.9 2004-09-10 22:38:52 paklein Exp $ */
/* created: paklein (12/09/1996) */
#ifndef _EAMFCC2D_H_
#define _EAMFCC2D_H_

/* base class */
#include "NL_E_MatT.h"

namespace Tahoe {

/* forward declarations */
class EAMFCC3DSym;

/** plane strain EAM material */
class EAMFCC2D: public NL_E_MatT
{
public:

	/* constructor */
	EAMFCC2D(void);

	/* destructor */
	virtual ~EAMFCC2D(void);
	
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
	EAMFCC3DSym* fEAM;
};

} // namespace Tahoe 
#endif /* _EAMFCC2D_H_ */
