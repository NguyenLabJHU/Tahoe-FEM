/* $Id: DiffusionMatListT.h,v 1.5 2003-06-09 06:53:11 paklein Exp $ */
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

	/** enum defining material types */
	enum TypeT {
        kLinear = 1,
     kNonLinear = 2};

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
