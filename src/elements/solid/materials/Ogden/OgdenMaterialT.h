/* $Id: OgdenMaterialT.h,v 1.1 2003-03-19 19:00:54 thao Exp $ */
/* created: tdn (3/17/2003) */
#ifndef _OGDEN_Material_T_H_
#define _OGDEN_Material_T_H_

/* base classes */
#include "OgdenIsotropicT.h"

namespace Tahoe {

/*forward declaration*/
class PotentialT;

/**Ogden material model*/
class OgdenMaterialT: public OgdenIsotropicT
{
public:

	/* constructor/destructor */
	OgdenMaterialT(ifstreamT& in, const FSMatSupportT& support);
	~OgdenMaterialT(void);
	
	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;

protected:

	/* principal values of the PK2 stress given principal values of the stretch 
	 * tensors, i.e., the principal stretches squared */
	virtual double StrainEnergyDensity(void);
	virtual void dWdE(const dArrayT& eigenstretch, dArrayT& eigenstress);
	virtual void ddWddE(const dArrayT& eigenstretch, dArrayT& eigenstress,
		dSymMatrixT& eigenmod);

 private:
	PotentialT* fPot;
	const double fthird;

};

} // namespace Tahoe 
#endif /* _OGDEN_Material_T_H_ */
