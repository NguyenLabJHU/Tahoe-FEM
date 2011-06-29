/* $Id: VIB3D.h,v 1.7 2004-07-15 08:27:51 paklein Exp $ */
/* created: paklein (04/20/1997) */
#ifndef _VIB_3D_H_
#define _VIB_3D_H_

/* base class */
#include "NL_E_MatT.h"
#include "VIB_E_MatT.h"
#include "SpherePointsT.h"

namespace Tahoe {

/* forward declarations */
class dMatrixT;

/** 3D isotropic VIB solver */
class VIB3D: public NL_E_MatT, public VIB_E_MatT
{
public:

	/* constructor */
	VIB3D(void);

	/* destructor */
	~VIB3D(void);
	
	/* set angle offset - for testing onset of amorphous behavior
	 * Angles given in degrees */
	void SetAngles(double phi, double theta);

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

	/* compute the symetric Cij reduced index matrix */
	virtual void ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli);	
	
	/* compute the symetric 2nd Piola-Kirchhoff reduced index vector */
	virtual void ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2);

	/* returns the strain energy density for the specified strain */
	virtual double ComputeEnergyDensity(const dSymMatrixT& E);
		
private:

	/* integration point generator */
	SpherePointsT*	fSphere;	
		
};

} // namespace Tahoe 
#endif /* _VIB_3D_H_ */
