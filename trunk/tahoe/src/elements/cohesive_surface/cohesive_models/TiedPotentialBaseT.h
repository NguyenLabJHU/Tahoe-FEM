/* $Id: TiedPotentialBaseT.h,v 1.2 2003-04-17 20:11:33 cjkimme Exp $ */
/* created: cjkimme (04/15/2002) */

#ifndef _TIED_POTENTIAL_BASE_T_H_
#define _TIED_POTENTIAL_BASE_T_H_

/* direct member */
#include "iArrayT.h"

#include "ios_fwd_decl.h"

namespace Tahoe {

/* constants for state variable flags */
const double kTiedNode = -100.;
const double kReleaseNextStep = -10;
const double kFirstFreeStep = -1.;
const double kFreeNode = 1.;

/** A base class for potentials using the TiedNodes KBC controller. */
class TiedPotentialBaseT
{
public:
	
	enum sbntmaT {kAverageCode = 2};

	/** constructor */
	TiedPotentialBaseT(void);
	
	~TiedPotentialBaseT(void);
	
	virtual bool NeedsNodalInfo(void) = 0;
	virtual int NodalQuantityNeeded(void) = 0;
	
	virtual bool InitiationQ(const double* sigma) = 0;

	virtual iArrayT& BulkGroups(void);
	
protected:

    iArrayT iBulkGroups;
};

} // namespace Tahoe 
#endif /* _TIED_POTENTIAL_BASE_T_H_ */
