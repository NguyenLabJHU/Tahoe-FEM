/* $Id: DiffusionMatListT.h,v 1.3.8.1 2002-10-28 06:49:16 paklein Exp $ */
/* created: paklein (10/02/1999) */
#ifndef _DIFFUSE_MAT_LIST_T_H_
#define _DIFFUSE_MAT_LIST_T_H_

/* base class */
#include "MaterialListT.h"

namespace Tahoe {

/* forward declarations */
class DiffusionT;
class DiffusionMatSupportT;

class DiffusionMatListT: public MaterialListT
{
public:

	/* constructors */
	DiffusionMatListT(int length, const DiffusionT& element_group);

	/* read material data from the input stream */
	virtual void ReadMaterialData(ifstreamT& in);
	
private:

	/* element group */
	const DiffusionT& fElementGroup;
	const DiffusionMatSupportT* fDiffMatSupport;
};

} // namespace Tahoe 
#endif /* _DIFFUSE_MAT_LIST_T_H_ */
