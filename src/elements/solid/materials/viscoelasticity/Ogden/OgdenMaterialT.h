/* $Id: OgdenMaterialT.h,v 1.1 2006-10-30 23:32:05 thao Exp $ */
/* created: tdn (3/17/2003) */
#ifndef _OGDEN_Material_T_H_
#define _OGDEN_Material_T_H_

/* base classes */
#include "OgdenBaseT.h"

namespace Tahoe {

/*forward declaration*/
class PotentialT;

/**Ogden material model*/
class OgdenMaterialT: public OgdenBaseT
{
public:

	/* constructor/destructor */
	OgdenMaterialT(ifstreamT& in, const FSMatSupportT& support);
	~OgdenMaterialT(void);

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
