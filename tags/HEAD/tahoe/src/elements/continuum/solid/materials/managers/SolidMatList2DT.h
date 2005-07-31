/* $Id: SolidMatList2DT.h,v 1.1.1.1 2001-01-29 08:20:25 paklein Exp $ */
/* created: paklein (02/14/1997)                                          */

#ifndef _MATLIST_2D_T_H_
#define _MATLIST_2D_T_H_

/* base class */
#include "SolidMatListT.h"

/* forward declaration */
class ElasticT;

class SolidMatList2DT: public SolidMatListT
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
