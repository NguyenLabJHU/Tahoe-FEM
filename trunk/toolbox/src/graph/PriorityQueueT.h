/* $Id: PriorityQueueT.h,v 1.4 2002-07-05 22:26:30 paklein Exp $ */
/* created: paklein (8/06/1996) */

#ifndef _PRIORITYQUEUET_H_
#define _PRIORITYQUEUET_H_

/* direct members */
#include "AutoArrayT.h"

namespace Tahoe {

/* forward declarations */
class iArrayT;

/* size parameters */
const int  kPriorityQueueSize = 25;

class PriorityQueueT
{
public:

	/* constructor */
	PriorityQueueT(iArrayT& priorities, int size = kPriorityQueueSize);
	PriorityQueueT(const iArrayT& values, iArrayT& priorities);

	/* destructor */
	~PriorityQueueT(void);
	
	/* add a value */
	void Add(int value);
	
	/* remove values - both return 0 if the queue is empty */
	int	PullHighest(int& value);
	int	PullLowest(int& value);
	
	/* take "half" */
	void ShrinkToHighest(void);
	void ShrinkToLowest(void);
	
private:

	/* shift all values down, overwriting the value at position dex */
	void ShiftDown(int dex);
	
	/* direct-mapped priorities:
	 *
	 *	priority = fPriorities[fQueue[i]] */
	void DirectHighest(int& value, int& dex) const;
	void DirectLowest(int& value, int& dex) const;
	
	/* index-mapped priorities
	 *
	 *	priority = fPriorities[i] */
	void IndexHighest(int& value, int& dex) const;
	void IndexLowest(int& value, int& dex) const;
	  	 	
private:

	int fMode;
	AutoArrayT<int> fQueue;

	iArrayT& fPriorities;
	
};

} // namespace Tahoe 
#endif /* _PRIORITYQUEUET_H_ */
