/* $Id: TiedPotentialBaseT.h,v 1.4 2003-04-22 19:02:06 cjkimme Exp $ */
/* created: cjkimme (04/15/2002) */

#ifndef _TIED_POTENTIAL_BASE_T_H_
#define _TIED_POTENTIAL_BASE_T_H_

/* direct member */
#include "iArrayT.h"

#include "ios_fwd_decl.h"

namespace Tahoe {

class dArrayT;

/* constants for state variable flags */
const double kTiedNode = -100.;
const double kReleaseNextStep = -10;
const double kFirstFreeStep = -1.;
const double kFreeNode = 1.;
const double kTieNextStep = 10.;

/** A base class for potentials using the TiedNodes KBC controller. */
class TiedPotentialBaseT
{
public:
	
	enum sbntmaT {kAverageCode = 2};

	/** constructor */
	TiedPotentialBaseT(void);
	
	~TiedPotentialBaseT(void);
	
	/* true if nodal release depends on bulk element groups */
	virtual bool NeedsNodalInfo(void) = 0;
	
	/* release condition depends on this bulk quantity */
	virtual int NodalQuantityNeeded(void) = 0;
	
	/* True if a nodal release condition is satisfied */
	virtual bool InitiationQ(const double* sigma) = 0;

	/* Bulk element groups needed for calculation of nodal release conditions */
	virtual iArrayT& BulkGroups(void);
	
	/* True if the tied potential may ask for nodes to be retied later */
	virtual bool NodesMayRetie(void) = 0;
	
	/* True if node should be retied */
	virtual bool RetieQ(const double* sigma, const ArrayT<double>& state, const dArrayT& jump_u); 
	
protected:

    iArrayT iBulkGroups;

};

inline bool TiedPotentialBaseT::RetieQ(const double* sigma, const ArrayT<double>& state, const dArrayT& jump_u)
{
#pragma unused(sigma)
#pragma unused(state)
#pragma unused(jump_u)
	return false;
}

} // namespace Tahoe 
#endif /* _TIED_POTENTIAL_BASE_T_H_ */
