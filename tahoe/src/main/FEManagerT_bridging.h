/* $Id: FEManagerT_bridging.h,v 1.1.2.1 2003-02-06 02:39:08 paklein Exp $ */
#ifndef _FE_MANAGER_BRIDGING_H_
#define _FE_MANAGER_BRIDGING_H_

/* base class */
#include "FEManagerT.h"

namespace Tahoe {

class FEManagerT_bridging: public FEManagerT
{
public:

	/** constructor */
	FEManagerT_bridging(ifstreamT& input, ofstreamT& output, CommunicatorT& comm);

	/** return list of ghost nodes */
	const iArrayT& GhostNodes(void) const { return fGhostNodes; };

	/** return list of ghost nodes */
	const iArrayT& NonGhostNodes(void) const { return fNonGhostNodes; };

	/** initialize nodes that follow the field computed by this instance */
	void SetFollowers(const iArrayT& nodes, const StringT& field,
		NodeManagerT& node_manager);

	/** initialize data for the driving field */
	void SetExactSolution(const iArrayT& nodes, const StringT& field,
		NodeManagerT& node_manager);

protected:

	/** read parameters. Collect ghost nodes */
	virtual void ReadParameters(InitCodeT init);

	/** write parameters to main out */
	virtual void WriteParameters(void) const;

	/** construct node manager */
	virtual void SetNodeManager(void);

private:

	/** list of my ghost nodes */
	iArrayT fGhostNodes;

	/** list of my non-ghost nodes */
	iArrayT fNonGhostNodes;
};

} /* namespace Tahoe */

#endif /* _FE_MANAGER_BRIDGING_H_ */
