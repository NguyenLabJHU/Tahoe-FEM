/* $Id: ReLabellerT.h,v 1.3 2002-07-05 22:26:30 paklein Exp $ */
/* created: paklein (08/05/1996)                                          */
/* relabels sequences based on connectivities registered with AddGroup(). */
/* connectivies can have an arbitrary MinorDim(), but the labels in       */
/* the range of 1...[max label] must all be used in the sequence(s)       */
/* sent to Renumber()                                                     */

#ifndef _RELABELLERT_H_
#define _RELABELLERT_H_

/* direct members */
#include "iArrayT.h"
#include "GraphT.h"
#include "RootedLevelT.h"
#include "PriorityQueueT.h"

namespace Tahoe {

/* forward declarations */
class iArray2DT;

class ReLabellerT
{
public:

	/* constructor */
	ReLabellerT(void);

	/* add a group to the graph */
	void AddGroup(const iArray2DT& moredata);
	void AddGroup(const RaggedArray2DT<int>& moredata);
	
	/* renumber the positive values in oldsequence - returns the number
	 * of re-sequenced labels */
	int Renumber(iArray2DT& oldsequence);

	/* renumber the positive values in oldsequences - returns the number
	 * of re-sequenced labels. labels are assumed to increase contiguously
	 * through the oldsequences */
	int Renumber(ArrayT<iArray2DT*>& oldsequences);

private:

	/* initialize data before renumbering */
	void Initialize(void);

	/* fill internal rooted level even for disconnected graphs */
	void BuildRootedLevel(void);

	/* select pseudo-diameter nodes */
	void SelectNodes(void);
	
	/* relabel */
	void NewSequence(void);
	void ReverseSequence(void);

	/* estimate size of the system in terms of the
	 * bandwidth and profile using the sequence */
	void ComputeSize(const iArrayT& sequence, int& bandwidth, int& profile);
	
	/* initialize priorities as specified by Sloan */
	void InitializePriorities(void);
	
	/* add a node to the labelling queue */
	void Queue(int nodenum);
	
	/* assign new number and update adjacent nodes */
	void NewNumber(int nodenum);
	
	/* make the specified node active and update adjacent nodes */
	void MakeActive(int nodenum);
	
private:

	iArrayT fSequence;
	GraphT  fGraph;
	
	/* workspace */
	iArrayT fStatus;   /* direct indexing */
	iArrayT fPriority; /* direct indexing */
	
	RootedLevelT   fRootedLevel;
	PriorityQueueT fLabelQueue;
	
	int fStartNode;
	int fEndNode;
	int	fCurrLabel;
	
	/* Sloan priority queue weights */
	int fWeight1; // default = 2
	int fWeight2; // default = 1
};

} // namespace Tahoe 
#endif /* _RELABELLERT_H_ */
