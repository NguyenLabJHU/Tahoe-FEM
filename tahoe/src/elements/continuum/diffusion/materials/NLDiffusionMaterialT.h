/* $Id: NLDiffusionMaterialT.h,v 1.2 2003-12-10 07:14:28 paklein Exp $ */
#ifndef _NL_DIFFUSION_MATERIALT_H_
#define _NL_DIFFUSION_MATERIALT_H_

/* base class */
#include "DiffusionMaterialT.h"

namespace Tahoe {

/* forward declarations */
class C1FunctionT;

/** interface for nonlinear materials for diffusion. Materials may
 * have temperature dependent conductivity and thermal capacity and
 * must return the associated derivatives */
class NLDiffusionMaterialT: public DiffusionMaterialT
{
public:

	/** constructor */
	NLDiffusionMaterialT(ifstreamT& in, const DiffusionMatSupportT& support);
	NLDiffusionMaterialT(void);

	/** destructor */
	~NLDiffusionMaterialT(void);

	/** \name print parameters */
	/*@{*/
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;
	/*@}*/

	/** \name parameters at the current field point */
	/*@{*/
	/** conductivity */
	virtual const dMatrixT& k_ij(void);

	/** heat flux */
	virtual const dArrayT& q_i(void);
	
	/** change in heat flux with temperature */
	virtual const dArrayT& dq_i_dT(void);

	/** specific heat */
	virtual double SpecificHeat(void) const;

	/** change in capacity with temperature */
	virtual double dCapacity_dT(void) const;
	/*@}*/

private:

	/** \name temperature dependence */
	/*@{*/
	/** variation of conductivity with temperature */
	C1FunctionT* fConductivityScaleFunction;

	/** variation of specific heat with temperature */
	C1FunctionT* fCpScaleFunction;
	/*@}*/

	/** temperature varying conductivity return value */
	dMatrixT fScaledConductivity;
};

} /* namespace Tahoe */

#endif /* _NL_DIFFUSION_MATERIALT_H_ */
