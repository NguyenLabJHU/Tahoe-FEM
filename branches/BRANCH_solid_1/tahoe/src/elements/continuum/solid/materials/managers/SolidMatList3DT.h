/* $Id: SolidMatList3DT.h,v 1.2.2.1 2001-06-22 14:18:16 paklein Exp $ */
/* created: paklein (02/14/1997)                                          */

#ifndef _MATLIST_3D_T_H_
#define _MATLIST_3D_T_H_

/* base class */
#include "SolidMatListT.h"
#include "MaterialT.h"

/* forward declaration */
class ElasticT;
class SmallStrainT;
class FiniteStrainT;

class SolidMatList3DT: public SolidMatListT, public MaterialT
{
public:

	/* constructors */
	SolidMatList3DT(int length, const ElasticT& element_group);

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

#endif /* _MATLIST_3D_T_H_ */
