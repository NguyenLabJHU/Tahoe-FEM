/* $Id: DiffusionMatListT.h,v 1.3.8.2 2002-11-13 08:44:20 paklein Exp $ */
/* created: paklein (10/02/1999) */
#ifndef _DIFFUSE_MAT_LIST_T_H_
#define _DIFFUSE_MAT_LIST_T_H_

/* base class */
#include "MaterialListT.h"

namespace Tahoe {

/* forward declarations */
class DiffusionMatSupportT;

/** list of materials for diffusion analysis */
class DiffusionMatListT: public MaterialListT
{
public:

	/** constructors */
	DiffusionMatListT(int length, const DiffusionMatSupportT& support);

	/** read material data from the input stream */
	virtual void ReadMaterialData(ifstreamT& in);
	
private:

	/** support for diffusion materials */
	const DiffusionMatSupportT& fDiffusionMatSupport;
};

} // namespace Tahoe 
#endif /* _DIFFUSE_MAT_LIST_T_H_ */
