/* $Id: TiedPotentialBaseT.h,v 1.7 2003-12-28 23:37:01 paklein Exp $ */
/* created: cjkimme (04/15/2002) */
#ifndef _TIED_POTENTIAL_BASE_T_H_
#define _TIED_POTENTIAL_BASE_T_H_

/* direct member */
#include "iArrayT.h"

#include "ios_fwd_decl.h"

namespace Tahoe {

class dArrayT;

/** a base class for potentials using the TiedNodes KBC controller. */
class TiedPotentialBaseT
{
public:
	
	enum sbntmaT {kAverageCode = 2};

	/** constructor */
	TiedPotentialBaseT(void);
	
	/** destructor */
	virtual ~TiedPotentialBaseT(void);
	
	/** true if nodal release depends on bulk element groups */
	virtual bool NeedsNodalInfo(void) const = 0;
	
	/** release condition depends on this bulk quantity */
	virtual int NodalQuantityNeeded(void) const = 0;

	/** true if nodal values need to be transformed to local frame. This should only 
	 * return true if the nodal value given by TiedPotentialBaseT::NodalQuantityNeeded
	 * is a stress tensor */
	virtual bool RotateNodalQuantity(void) const = 0;

	/** true if a nodal release condition is satisfied */
	virtual bool InitiationQ(const nArrayT<double>& sigma) const = 0;

	/** bulk element groups needed for calculation of nodal release conditions */
	virtual const iArrayT& BulkGroups(void) const;
	
	/** true if the tied potential may ask for nodes to be retied later */
	virtual bool NodesMayRetie(void) const = 0;
	
	/** true if node should be retied */
	virtual bool RetieQ(const nArrayT<double>& sigma, const ArrayT<double>& state, const dArrayT& jump_u) const;

	/** \name constants for state variable flags */
	/*@{*/
	
	/** location in state variable array of the state flag */
	virtual int TiedStatusPosition(void) const = 0;
	
	static const double kTiedNode;
	static const double kReleaseNextStep;
	static const double kFirstFreeStep;
	static const double kFreeNode;
	static const double kTieNextStep;
	/*@}*/
	
protected:

    iArrayT iBulkGroups;
};

inline bool TiedPotentialBaseT::RetieQ(const nArrayT<double>& sigma, const ArrayT<double>& state, const dArrayT& jump_u) const
{
#pragma unused(sigma)
#pragma unused(state)
#pragma unused(jump_u)
	return false;
}

} // namespace Tahoe 
#endif /* _TIED_POTENTIAL_BASE_T_H_ */
