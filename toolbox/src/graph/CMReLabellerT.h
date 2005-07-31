/* $Id: CMReLabellerT.h,v 1.1.1.1 2001-01-25 20:56:27 paklein Exp $ */
/* created: paklein (08/05/1996)                                          */
/* relabels sequences based on connectivities registered with AddGroup(). */
/* connectivies can have an arbitrary MinorDim(), but the labels in       */
/* the range of 1...[max label] must all be used in the sequence(s)       */
/* sent to Renumber()                                                     */

#ifndef _CM_RELABELLERT_H_
#define _CM_RELABELLERT_H_

/* direct members */
#include "iArrayT.h"
#include "GraphT.h"
#include "RootedLevelT.h"
#include "PriorityQueueT.h"

/* forward declarations */
class iArray2DT;

class CMReLabellerT
{
public:

	/* constructor */
	CMReLabellerT(void);

	/* add a group to the graph */
	void AddGroup(const iArray2DT& moredata);
	void AddGroup(const RaggedArray2DT<int>& moredata);
	
	/* renumber the positive values in oldsequence - returns the number
	 * of re-sequenced labels */
	int Renumber(iArray2DT& oldsequence, bool is_disconnected);

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
	
	/* relabel with Sloan */
	void NewSequence(void);
	
	/* relabel with Reverse CutHill-McKee */
	void CMSequence(void);
	
	/* computes bandwidth and profile */
	void ComputeSize(const iArrayT& sequence, int& bandwidth, int& profile);
	
	/* reverse sequence */
	void ReverseSequence(void);
	
	/* add node to map */
	void AddNode(int node, int count);
	
	/* sort edges by minimum degree */
	void SortByMinDegree(ArrayT<int>& edges);
	
	/* AZ_sort: sorter taken from AZTEC az_sort.c file */
	void AZ_sort(int list[], int N, int list2[], double list3[]);
	
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
	
	GraphT fGraph;
	
	/* workspace */
	iArrayT fStatus;   /* direct indexing */
	iArrayT fPriority; /* direct indexing */
	
	RootedLevelT   fRootedLevel;
	PriorityQueueT fLabelQueue;
	
	int fStartNode;
	int fEndNode;
	int	fCurrLabel;
	
	int kWeight1;
	int kWeight2;
};

#endif /* _CM_RELABELLERT_H_ */
