/* $Id: SIERRA_Material_DB.h,v 1.1 2003-03-05 02:27:52 paklein Exp $ */
#ifndef _SIERRA_MAT_DB_H_
#define _SIERRA_MAT_DB_H_

/* direct members */
#include "MapT.h"
#include "StringT.h"

namespace Tahoe {

/* forward declarations */
class SIERRA_Material_Data;

/** database for Sierra materials parameters */
class SIERRA_Material_DB
{
public:

	/** constructor */
	SIERRA_Material_DB(void);

	/** destructor */
	~SIERRA_Material_DB(void);

private:

	/** array of material data cards. Maps material names to
	 * pointers to material data cards. */
	MapT<StringT, SIERRA_Material_Data*> fMaterialData;
};

} /* namespace Tahoe */

#endif /* _SIERRA_MAT_DB_H_ */
