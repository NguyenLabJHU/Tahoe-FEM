/* $Id: HookeanMatT.h,v 1.7 2004-07-22 21:09:32 paklein Exp $ */
/* created: paklein (06/09/1997) */
#ifndef _HOOKEAN_MAT_H_
#define _HOOKEAN_MAT_H_

/* base class */
#include "ParameterInterfaceT.h"

/* direct members */
#include "dMatrixT.h"

namespace Tahoe {

/* forward declarations */
class dSymMatrixT;

/** base class for all Hookean materials, defined as:
 * 	stress_ij = moduli_ijkl strain_kl */
class HookeanMatT: virtual public ParameterInterfaceT 
{
public:

	/** constructors */
	HookeanMatT(int nsd);
	HookeanMatT(void);

	/** destructor */
	virtual ~HookeanMatT(void);

	/** dimension work space before calling HookeanMatT::TakeParameterList */
	void Dimension(int nsd);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** return the description of the given inline subordinate parameter list */
	virtual void DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/* set (material) tangent modulus */
	virtual void SetModulus(dMatrixT& modulus);
	const dMatrixT& Modulus(void) const { return fModulus; };

	/* symmetric stress */
	void HookeanStress(const dSymMatrixT& strain, dSymMatrixT& stress) const;

	/* strain energy density */
	double HookeanEnergy(const dSymMatrixT& strain) const;
		
private:

	/* (constant) modulus */
	dMatrixT fModulus;
};

} // namespace Tahoe 
#endif /* _HOOKEAN_MAT_H_ */
