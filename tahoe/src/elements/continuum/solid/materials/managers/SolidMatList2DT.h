/* $Id: SolidMatList2DT.h,v 1.2.2.1 2001-06-22 14:18:15 paklein Exp $ */
/* created: paklein (02/14/1997)                                          */

#ifndef _MATLIST_2D_T_H_
#define _MATLIST_2D_T_H_

/* base classes */
#include "SolidMatListT.h"
#include "MaterialT.h"

/* forward declaration */
class ElasticT;
class SmallStrainT;
class FiniteStrainT;

class SolidMatList2DT: public SolidMatListT, public MaterialT
{
public:

	/* constructor */
	SolidMatList2DT(int length, const ElasticT& element_group);

	/* read material data from the input stream */
	virtual void ReadMaterialData(ifstreamT& in);

private:
	
	/* errror messages */
	void Error_no_small_strain(ostream& out, int matcode) const;
	void Error_no_finite_strain(ostream& out, int matcode) const;

private:

	const ElasticT&      fElementGroup;
	const SmallStrainT*  fSmallStrain;
	const FiniteStrainT* fFiniteStrain;
};

#endif /* _MATLIST_2D_T_H_ */
