/* $Id: FSSolidMatList3DT.h,v 1.1.2.1 2004-01-21 19:10:18 paklein Exp $ */
#ifndef _FS_MATLIST_3D_T_H_
#define _FS_MATLIST_3D_T_H_

/* base class */
#include "SolidMatListT.h"
#include "SolidT.h"

namespace Tahoe {

/** materials list for 3D structural analysis */
class FSSolidMatList3DT: public SolidMatListT, public SolidT
{
public:

	/** constructors */
	FSSolidMatList3DT(int length, const FSMatSupportT& support);
	FSSolidMatList3DT(void);

	/** read material data from the input stream */
	virtual void ReadMaterialData(ifstreamT& in);

private:

	/** support for finite strain materials */
	const FSMatSupportT* fFSMatSupport;
};

} /* namespace Tahoe */

#endif /* _FS_MATLIST_3D_T_H_ */
