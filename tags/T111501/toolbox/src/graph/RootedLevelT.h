/* $Id: RootedLevelT.h,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
/* created: paklein (08/05/1996)                                          */

#ifndef _ROOTEDLEVEL_T_H_
#define _ROOTEDLEVEL_T_H_

#include "Environment.h"
#include "Constants.h"

/* forward declarations */
class GraphT;
class iArrayT;

class RootedLevelT
{
public:

	/* constructor */
	RootedLevelT(void);

	/* destructor */
	~RootedLevelT(void);
	
	/* generate the level rooted structure from the given graph
	 * rooted at the given node */
	void MakeRootedLevel(const GraphT& graph, int rootnode);
	void MakePartialRootedLevel(const GraphT& graph, int rootnode, bool clear);
		// only those connected to rootnode

	/* returns the root node number of last structure */
	int RootNode(void) const;
		
	/* returns the number of levels */
	int Depth(void) const;
	
	/* returns the width of the specified level */
	int Width(int level) const;
	int MaxWidth(void) const;
	
	/* return the node numbers on the specified level */
	void NodesOnLevel(iArrayT& nodes, int level) const;
	
	/* return list of nodes used */
	void NodesUsed(iArrayT& nodes_used) const;
	
	/* initialize the priorities as specified in Sloan.  Assumes DIRECT
	 * INDEXING in priorities */
	void InitializePriorities(iArrayT& priorities, int W2) const;

	/* copy in the labelled branches of the rooted level structure */
	void CopyIn(const RootedLevelT& source);

private:

	/* compute the number of nodes in each level */
	void ComputeWidths(void);
	
private:

	int	 fNumNodes;
	int* fLevels;

	int	 fNumLevels;
	int	 fMaxWidth;
	int* fWidths;
	int  fRoot;

	/* return value */
	int* fNodesOnLevel;
};

/* inlines */
inline int RootedLevelT::RootNode(void) const { return fRoot; }
inline int RootedLevelT::Depth(void) const { return fNumLevels; }
inline int RootedLevelT::Width(int level) const { return fWidths[level]; }
inline int RootedLevelT::MaxWidth(void) const { return fMaxWidth; }

#endif /* _ROOTEDLEVEL_T_H_ */
