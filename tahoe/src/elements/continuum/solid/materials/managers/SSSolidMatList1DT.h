/* $Id: SSSolidMatList1DT.h,v 1.1.2.1 2004-01-21 19:10:18 paklein Exp $ */
#ifndef _SS_MATLIST_1D_T_H_
#define _SS_MATLIST_1D_T_H_

/* base classes */
#include "SolidMatListT.h"
#include "SolidT.h"

namespace Tahoe {

/** materials list for 1D structural analysis */
class SSSolidMatList1DT: public SolidMatListT, public SolidT
{
public:

	/** constructor */
	SSSolidMatList1DT(int length, const SSMatSupportT& support);
	SSSolidMatList1DT(void);

	/** read material data from the input stream */
	virtual void ReadMaterialData(ifstreamT& in);

private:

	/** support for small strain materials */
	const SSMatSupportT* fSSMatSupport;

	/** support for gradient enhanced small strain materials */
	const GradSSMatSupportT* fGradSSMatSupport;
};

} /* namespace Tahoe */

#endif /* _SS_MATLIST_1D_T_H_ */
