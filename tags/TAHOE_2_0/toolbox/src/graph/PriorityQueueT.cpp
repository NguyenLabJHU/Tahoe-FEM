/* $Id: PriorityQueueT.cpp,v 1.4 2002-10-20 22:39:01 paklein Exp $ */
/* created: paklein (8/06/1996) */

#include "PriorityQueueT.h"
#include <math.h>
#include <string.h>
#include <limits.h>
#include "iArrayT.h"

/* priority modes */

using namespace Tahoe;

const int kDirectMapped = 0;
const int kIndexMapped  = 1;

/* constructor */
PriorityQueueT::PriorityQueueT(iArrayT& priorities, int size):
	fMode(kDirectMapped),
	fQueue(50),
	fPriorities(priorities)
{
#pragma unused(size)
}

PriorityQueueT::PriorityQueueT(const iArrayT& values, iArrayT& priorities): 
	fMode(kIndexMapped), 
	fQueue(values),
	fPriorities(priorities)
{

}

/* destructor */
PriorityQueueT::~PriorityQueueT(void)
{

}
	
/* add a value */
void PriorityQueueT::Add(int value)
{	
	/* check */
	if (value >= fPriorities.Length()) throw ExceptionT::kGeneralFail;
	
	/* append to the queue */
	fQueue.Append(value);
}
	
/* remove values - both return 0 if the queue is empty */
int	PriorityQueueT::PullHighest(int& value)
{
	/* queue is empty */
	if (fQueue.Length() == 0) return 0;

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
	if (fQueue.Length() == 0) return(0);
	
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
	int newsize = int(floor((fQueue.Length() + 2)/2.0));	
	if (newsize < fQueue.Length())
	{
		/* fill into temp */
		ArrayT<int> new_queue(newsize);
		for (int i = 0; i < newsize; i++)
			PullHighest(new_queue[i]);
			
		/* copy into the queue */
		fQueue = new_queue;
	}
}

void PriorityQueueT::ShrinkToLowest(void)
{
	int newsize = int(floor((fQueue.Length() + 2)/2.0));
	if (newsize < fQueue.Length())
	{
		/* fill into temp */
		ArrayT<int> new_queue(newsize);
		for (int i = 0; i < newsize; i++)
			PullLowest(new_queue[i]);
			
		/* copy into the queue */
		fQueue = new_queue;
	}
}


/************************************************************************
* Private
************************************************************************/

/* shift all values down, overwriting the value at position dex */
void PriorityQueueT::ShiftDown(int dex)
{
	fQueue.DeleteAt(dex);
}

/* direct-mapped priorities:
*
*	priority = fPriorities[fQueue[i]] */
void PriorityQueueT::DirectHighest(int& value, int& dex) const
{
	int highest = INT_MIN;
	dex = -1;
	
	/* find maximum */
	int curr_size = fQueue.Length();	
	for (int i = 0; i < curr_size; i++)
	{
		int currpriority = fPriorities[fQueue[i]];
	
		/* returns the first occurence of the highest priority */
		if (currpriority > highest)
		{
			highest = currpriority;
			dex     = i;
		}	
	}
	value = fQueue[dex];
}

void PriorityQueueT::DirectLowest(int& value, int& dex) const
{
	int lowest = INT_MAX;
	dex = -1;
	
	/* find minumum */	
	int curr_size = fQueue.Length();	
	for (int i = 0; i < curr_size; i++)
	{
		int currpriority = fPriorities[fQueue[i]];
	
		/* returns the first occurence of the lowest priority */
		if (currpriority < lowest)
		{
			lowest = currpriority;
			dex    = i;
		}	
	}
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
	int curr_size = fQueue.Length();
	for (int i = 0; i < curr_size; i++)
	{
		int currpriority = fPriorities[i];
	
		/* returns the first occurence of the highest priority */
		if (currpriority > highest)
		{
			highest = currpriority;
			dex     = i;
		}	
	}
	value = fQueue[dex];
}

void PriorityQueueT::IndexLowest(int& value, int& dex) const
{
	int lowest = INT_MAX;
	dex = -1;
	
	/* find minumum */
	int curr_size = fQueue.Length();
	for (int i = 0; i < curr_size; i++)
	{
		int currpriority = fPriorities[i];
	
		/* returns the first occurence of the lowest priority */
		if (currpriority < lowest)
		{
			lowest = currpriority;
			dex    = i;
		}	
	}
	value = fQueue[dex];
}
