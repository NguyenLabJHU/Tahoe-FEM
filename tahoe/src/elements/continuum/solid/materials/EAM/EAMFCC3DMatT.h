/* $Id: EAMFCC3DMatT.h,v 1.7 2004-07-15 08:26:47 paklein Exp $ */
/* created: paklein (10/25/1998) */
#ifndef _EAMFCC3DMatT_H_
#define _EAMFCC3DMatT_H_

/* base class */
#include "NL_E_MatT.h"

namespace Tahoe {

/* forward declarations */
class EAMFCC3DSym;

/** plane strain EAM material */
class EAMFCC3DMatT: public NL_E_MatT
{
public:

	/* constructor */
	EAMFCC3DMatT(void);

	/* destructor */
	virtual ~EAMFCC3DMatT(void);

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

} /* namespace Tahoe */

#endif /* _EAMFCC3DMatT_H_ */
