/* $Id: LocalArrayGroupT.h,v 1.3 2003-01-27 06:42:44 paklein Exp $ */
/* created: paklein (09/11/1998) */
#ifndef _LOCALARRAY_GROUP_T_H_
#define _LOCALARRAY_GROUP_T_H_

/* base class */
#include "MemoryGroupT.h"

/* direct members */
#include "LocalArrayT.h"

namespace Tahoe {

/** Class to manage a list of group of dynamically changing size
 * LocalArrayT's. All arrays must initially be of the same dimension
 * \note all registered arrays will be shallow.
 */
class LocalArrayGroupT: private MemoryGroupT<double>
{
public:

	/* constructor */
	LocalArrayGroupT(int headroom);

	/* add array to list of managed */
	void Register(LocalArrayT& localarray);
	bool IsRegistered(const LocalArrayT& localarray) const;

	/* set number of nodes */
	void SetNumberOfNodes(int numnodes);

	/* accessors */	
	int NumberOfNodes(void) const;
	int MinorDim(void) const;
	
private:

	/* dimensions */
	int fNumNodes;
	int fMinorDim;	
};

/* inlines */

/* check */
inline bool LocalArrayGroupT::IsRegistered(const LocalArrayT& localarray) const
{
	/* inherited */
	return MemoryGroupT<double>::IsRegistered(localarray);
}

/* accessors */	
inline int LocalArrayGroupT::NumberOfNodes(void) const { return fNumNodes; }
inline int LocalArrayGroupT::MinorDim(void) const  { return fMinorDim; }

} // namespace Tahoe 
#endif /* _LOCALARRAY_GROUP_T_H_ */
