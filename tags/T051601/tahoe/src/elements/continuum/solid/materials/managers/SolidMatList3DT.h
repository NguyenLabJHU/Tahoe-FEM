/* $Id: SolidMatList3DT.h,v 1.2 2001-04-27 10:53:30 paklein Exp $ */
/* created: paklein (02/14/1997)                                          */

#ifndef _MATLIST_3D_T_H_
#define _MATLIST_3D_T_H_

/* base class */
#include "SolidMatListT.h"
#include "MaterialT.h"

/* forward declaration */
class ElasticT;

class SolidMatList3DT: public SolidMatListT, public MaterialT
{
public:

	/* constructors */
	SolidMatList3DT(int length, const ElasticT& element_group);

	/* read material data from the input stream */
	virtual void ReadMaterialData(ifstreamT& in);

private:

	const ElasticT& fElementGroup;
};

#endif /* _MATLIST_3D_T_H_ */
