/* $Id: SIERRA_Material_Data.h,v 1.1 2003-03-05 02:27:52 paklein Exp $ */
#ifndef _SIERRA_MAT_DATA_H_
#define _SIERRA_MAT_DATA_H_

/* direct members */
#include "StringT.h"

namespace Tahoe {

/** data for a single Sierra material */
class SIERRA_Material_Data
{
public:

	/** constructor */
	SIERRA_Material_Data(const StringT& name);

private:

	/** material name */
	StringT fName;
};

} /* namespace Tahoe */

#endif /* _SIERRA_MAT_DATA_H_ */
