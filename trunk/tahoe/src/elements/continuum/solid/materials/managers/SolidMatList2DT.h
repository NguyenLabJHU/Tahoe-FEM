/* $Id: SolidMatList2DT.h,v 1.8 2002-10-05 20:04:16 paklein Exp $ */
/* created: paklein (02/14/1997) */

#ifndef _MATLIST_2D_T_H_
#define _MATLIST_2D_T_H_

/* base classes */
#include "StructuralMatListT.h"
#include "MaterialT.h"

namespace Tahoe {

/* forward declaration */
class ElasticT;
class SmallStrainT;
class FiniteStrainT;
class MultiScaleT;

class SolidMatList2DT: public StructuralMatListT, public MaterialT
{
public:

	/* constructor */
	SolidMatList2DT(int length, const ElasticT& element_group);

	/* read material data from the input stream */
	virtual void ReadMaterialData(ifstreamT& in);

	/** return true if the list contains plane stress models */
	virtual bool HasPlaneStress(void) const;

private:
	
	/* errror messages */
	void Error_no_small_strain(ostream& out, int matcode) const;
	void Error_no_finite_strain(ostream& out, int matcode) const;
	void Error_no_multi_scale(ostream& out, int matcode) const;

private:

	const ElasticT&      fElementGroup;
	const SmallStrainT*  fSmallStrain;
	const FiniteStrainT* fFiniteStrain;
	const MultiScaleT*   fMultiScale;
};

} // namespace Tahoe 
#endif /* _MATLIST_2D_T_H_ */
