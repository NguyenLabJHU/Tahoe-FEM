/* $Id: PriorityQueueT.cpp,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
/* created: paklein (8/06/1996)                                           */

#include "PriorityQueueT.h"
#include <math.h>
#include <string.h>
#include <limits.h>
#include "iArrayT.h"

/* priority modes */
const int kDirectMapped = 0;
const int kIndexMapped  = 1;

/* constructor */
PriorityQueueT::PriorityQueueT(iArrayT& priorities, int size):
	fMode(kDirectMapped), fLogicalSize(size), fCurrSize(0),
	fPriorities(priorities)
{
	fQueue = new int[fLogicalSize];
	if (!fQueue) throw eOutOfMemory;
}

PriorityQueueT::PriorityQueueT(const iArrayT& values,
	iArrayT& priorities): fMode(kIndexMapped), fPriorities(priorities)
{
	fLogicalSize = values.Length();
	fQueue = new int[fLogicalSize];
	if (!fQueue) throw eOutOfMemory;
	
	/* copy over */
	memcpy(fQueue, values.Pointer(), sizeof(int)*fLogicalSize);
	
	fCurrSize = fLogicalSize;
}

/* destructor */
PriorityQueueT::~PriorityQueueT(void)
{
	delete[] fQueue;
}
	
/* add a value */
void PriorityQueueT::Add(int value)
{	
	/* check */
	if (value >= fPriorities.Length()) throw eGeneralFail;

	/* need to allocate more space */
	if (fCurrSize == fLogicalSize)
	{
		int* temp     = fQueue;
		fLogicalSize *= 2;
		
		fQueue = new int[fLogicalSize];
		if (!fQueue) throw eOutOfMemory;
		
		/* byte copy */
		memcpy(fQueue, temp, sizeof(int)*fCurrSize);
		
		delete[] temp;
	}
		
	/* add value to the end */
	fQueue[fCurrSize] = value;
	
	fCurrSize++;	
}
	
/* remove values - both return 0 if the queue is empty */
int	PriorityQueueT::PullHighest(int& value)
{
	/* queue is empty */
	if (fCurrSize == 0) return 0;

	/* find maximum */
	int dex;
	if (fMode == kDirectMapped)
		DirectHighest(value, dex);
	else
		IndexHighest(value, dex);

	/* update queue */
	ShiftDown(dex);
		
	return 1;
}

int	PriorityQueueT::PullLowest(int& value)
{
	/* queue is empty */
	if (fCurrSize == 0)
		return(0);
	
	/* find minumum */	
	int dex;
	if (fMode == kDirectMapped)
		DirectLowest(value, dex);
	else
		IndexLowest(value, dex);

	/* update queue */
	ShiftDown(dex);
		
	return 1;
}

/* take "half" */
void PriorityQueueT::ShrinkToHighest(void)
{
	int newsize = int(floor((fCurrSize+2)/2.0));
	
	if (newsize < fCurrSize)
	{
		int* newqueue = new int[fLogicalSize];
		if (!newqueue) throw eOutOfMemory;
		
		for (int i = 0; i < newsize; i++)
			PullHighest(newqueue[i]);
			
		/* replace */
		int* temp = fQueue;
		fQueue    = newqueue;
		fCurrSize = newsize;
		
		/* free space */
		delete[] temp;
	}
}

void PriorityQueueT::ShrinkToLowest(void)
{
	int newsize = int(floor((fCurrSize+2)/2.0));
	
	if (newsize < fCurrSize)
	{
		int* newqueue = new int[fLogicalSize];
		if (!newqueue) throw eOutOfMemory;
	
		for (int i = 0; i < newsize; i++)
			PullLowest(newqueue[i]);
			
		/* replace */
		int* temp = fQueue;
		fQueue    = newqueue;
		fCurrSize = newsize;
		
		/* free space */
		delete[] temp;
	}
}


/************************************************************************
* Private
************************************************************************/

/* shift all values down, overwriting the value at position dex */
void PriorityQueueT::ShiftDown(int dex)
{
	/* check  */
	if (dex < 0) throw eGeneralFail;

	int* p = fQueue + dex;
	
	/* copy bytes */
	memmove(p, p+1, sizeof(int)*(fCurrSize - dex));
	
	/* update size */
	fCurrSize--;
}

/* direct-mapped priorities:
*
*	priority = fPriorities[ fQueue[i] - 1] */
void PriorityQueueT::DirectHighest(int& value, int& dex) const
{
	int highest = INT_MIN;
	dex = -1;
	
	/* find maximum */	
	for (int i = 0; i < fCurrSize; i++)
	{
		int currpriority = fPriorities[fQueue[i]];
	
		/* returns the first occurence of the highest priority */
		if (currpriority > highest)
		{
			highest = currpriority;
			dex     = i;
		}	
	}

#if __option (extended_errorcheck)
	if (dex < 0 || dex >= fCurrSize) throw eGeneralFail;
#endif

	value = fQueue[dex];
}

void PriorityQueueT::DirectLowest(int& value, int& dex) const
{
	int lowest = INT_MAX;
	dex = -1;
	
	/* find minumum */	
	for (int i = 0; i < fCurrSize; i++)
	{
		int currpriority = fPriorities[fQueue[i]];
	
		/* returns the first occurence of the lowest priority */
		if (currpriority < lowest)
		{
			lowest = currpriority;
			dex    = i;
		}	
	}

#if __option (extended_errorcheck)
	if (dex < 0 || dex >= fCurrSize) throw eGeneralFail;
#endif

	value = fQueue[dex];
}

/* index-mapped priorities
*
*	priority = fPriorities[i] */
void PriorityQueueT::IndexHighest(int& value, int& dex) const
{
	int highest = INT_MIN;
	dex = -1;
	
	/* find maximum */	
	for (int i = 0; i < fCurrSize; i++)
	{
		int currpriority = fPriorities[i];
	
		/* returns the first occurence of the highest priority */
		if (currpriority > highest)
		{
			highest = currpriority;
			dex     = i;
		}	
	}

#if __option (extended_errorcheck)
	if (dex < 0 || dex >= fCurrSize) throw eGeneralFail;
#endif

	value = fQueue[dex];
}

void PriorityQueueT::IndexLowest(int& value, int& dex) const
{
	int lowest = INT_MAX;
	dex = -1;
	
	/* find minumum */	
	for (int i = 0; i < fCurrSize; i++)
	{
		int currpriority = fPriorities[i];
	
		/* returns the first occurence of the lowest priority */
		if (currpriority < lowest)
		{
			lowest = currpriority;
			dex    = i;
		}	
	}

#if __option (extended_errorcheck)
	if (dex < 0 || dex >= fCurrSize) throw eGeneralFail;
#endif

	value = fQueue[dex];
}
