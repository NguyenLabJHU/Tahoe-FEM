/* $Id: SolidMatList2DT.h,v 1.2 2001-04-27 10:53:29 paklein Exp $ */
/* created: paklein (02/14/1997)                                          */

#ifndef _MATLIST_2D_T_H_
#define _MATLIST_2D_T_H_

/* base classes */
#include "SolidMatListT.h"
#include "MaterialT.h"

/* forward declaration */
class ElasticT;

class SolidMatList2DT: public SolidMatListT, public MaterialT
{
public:

	/* constructor */
	SolidMatList2DT(int length, const ElasticT& element_group);

	/* read material data from the input stream */
	virtual void ReadMaterialData(ifstreamT& in);
	
private:

	const ElasticT& fElementGroup;
};

#endif /* _MATLIST_2D_T_H_ */
