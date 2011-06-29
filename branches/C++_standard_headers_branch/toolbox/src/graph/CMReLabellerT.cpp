/* $Id: CMReLabellerT.cpp,v 1.4 2003-11-21 22:41:54 paklein Exp $ */
/* created: paklein (08/05/1996)                                          */

#include "CMReLabellerT.h"
#include <limits.h>
#include <string.h>
#include "iArray2DT.h"
#include "AutoArrayT.h"

/* status codes */

using namespace Tahoe;

const int kInActive   = 0;
const int kPreActive  = 1;
const int kActive     = 2;
const int kPostActive = 3;

/* constructor */
CMReLabellerT::CMReLabellerT(void): fLabelQueue(fPriority),
	fStartNode(0),
	fEndNode(0),
	fCurrLabel(0)
{

}

/* add a group to the graph */
void CMReLabellerT::AddGroup(const iArray2DT& moredata)
{
	fGraph.AddGroup(moredata);
}

void CMReLabellerT::AddGroup(const RaggedArray2DT<int>& moredata)
{
	fGraph.AddGroup(moredata);
}

/* renumber the oldsequence */
int CMReLabellerT::Renumber(iArray2DT& oldsequence, bool is_disconnected)
{
	/* make graph */
	fGraph.MakeGraph();

	/* initialize */
	Initialize();

	/* pick start and end nodes */
	if (is_disconnected)
		BuildRootedLevel();
	else
		SelectNodes();
		
	/* relabel */
	CMSequence();
	ReverseSequence();
	
	/* re-sequence positive values */	
	int subdim = oldsequence.MinorDim();
	int	label = 0;
	for (int i = 0; i < fGraph.NumNodes(); i++)
	{
		int* pseq = oldsequence(fSequence[i]);
		for (int j = 0; j < subdim; j++)
		{	
			int& oldlabel = *pseq++;
			if (oldlabel > 0) oldlabel = ++label;
		}
	}
	
	return label;
}

/* renumber the positive values in oldsequences - returns the number
* of re-sequenced labels. labels are assumed to increase contiguously
* through the rows of the oldsequences */
int CMReLabellerT::Renumber(ArrayT<iArray2DT*>& oldsequences)
{
	/* first row handled by each sequence */
	iArrayT maxrow(oldsequences.Length());
	for (int k = 0; k < oldsequences.Length(); k++)
	{
		maxrow[k] = oldsequences[k]->MajorDim();

		/* offset from previous sequence */
		if (k > 0) maxrow[k] += maxrow[k-1];
	}

	/* make graph */
	fGraph.MakeGraph();

	/* initialize */
	Initialize();

	/* pick start and end nodes */
	SelectNodes();

	/* relabel */
	CMSequence();
	ReverseSequence();
	
	/* re-sequence positive values */	
	int	label = 0;
	for (int i = 0; i < fGraph.NumNodes(); i++)
	{
		/* next node in line for new labels */
		int next = fSequence[i];
		
		/* get pointer to the correct sequence */
		int seqnum = 0;
		while (next >= maxrow[seqnum]) seqnum++;
		if (seqnum > 0) next -= maxrow[seqnum - 1];
		
		iArray2DT& currsequence = *oldsequences[seqnum];
	
		int* pseq  = currsequence(next);
		int subdim = currsequence.MinorDim();
		for (int j = 0; j < subdim; j++)
		{	
			int& oldlabel = *pseq++;
		
			if (oldlabel > 0) oldlabel = ++label;
		}
	}

	return label;
}

/************************************************************************
* Private
************************************************************************/

/* initialize data before renumbering */
void CMReLabellerT::Initialize(void)
{
	int numnodes = fGraph.NumNodes();

	/* allocate space */
	fSequence.Dimension(numnodes);
	fStatus.Dimension(numnodes);
	fPriority.Dimension(numnodes);

	/* initialize */
	fSequence =-1;
	fStatus   =-1;
	fPriority =-1;
}

/* fill internal rooted level even for disconnected graphs with
* the root for each branch selected by minimum width */
void CMReLabellerT::BuildRootedLevel(void)
{
	/* dimensions */
	int nnd = fGraph.NumNodes(); // graph needs to be set first

	/* work space */
	RootedLevelT rooted_level;
	iArrayT nodes_used;
	ArrayT<bool> used(nnd);

	fStartNode = -1;
	fEndNode   = -1;
	used = false;
	int count = 0;
	while (count < nnd)
	{
		/* select root node */
		int root_node = 0;
		while (root_node < nnd && used[root_node] == true)
			root_node++;
	
		/* first cut */
		rooted_level.MakePartialRootedLevel(fGraph, root_node, true);

		/* find node with mininum width and mark */
		rooted_level.NodesUsed(nodes_used);
		int min_degree = fGraph.Degree(root_node);
		for (int i = 0; i < nodes_used.Length(); i++)
		{
			int   node_i = nodes_used[i];
			int degree_i = fGraph.Degree(node_i);
			
			/* check degree */
			if (degree_i < min_degree)
			{
				min_degree = degree_i;
				root_node  = node_i;
			}
			
			/* mark as used */
			used[node_i] = true;
		}
		
		/* second cut */
		rooted_level.MakePartialRootedLevel(fGraph, root_node, true);
		
		/* set start node */
		if (fStartNode == -1) fStartNode = root_node;

		/* copy new structure in */	
		fRootedLevel.CopyIn(rooted_level);
		count += nodes_used.Length();
	}
	
	/* check all nodes used */
	if (count != nnd) throw ExceptionT::kGeneralFail;
	
	/* set end node from last level */
	fRootedLevel.NodesOnLevel(nodes_used, fRootedLevel.Depth() - 1);
	fEndNode = nodes_used[0];
}


/* select pseudo-diameter nodes */
void CMReLabellerT::SelectNodes(void)
{
	fStartNode = fGraph.MinDegreeNode();
	fRootedLevel.MakeRootedLevel(fGraph, fStartNode);
	
	cout << "degree of start node = "<< fGraph.Degree(fStartNode) << '\n';
	cout << "degree of node 608   = "<<  fGraph.Degree(608) << '\n';
	
	iArrayT	topnodes, degrees;
	fEndNode = -1;
	
	while (fEndNode == -1)
	{
		int h_max = fRootedLevel.Depth();
		int	w_min = INT_MAX;

		/* collect top level info */
		fRootedLevel.NodesOnLevel(topnodes, h_max - 1);
		degrees.Dimension(topnodes.Length()); //really want to Dimension() every time?
		fGraph.ReturnDegrees(topnodes, degrees);
		
		/* order and halve */
		PriorityQueueT topqueue(topnodes, degrees);
		topqueue.ShrinkToLowest();
		
		int currnode;
		while ( topqueue.PullLowest(currnode) )
		{
			fRootedLevel.MakeRootedLevel(fGraph, currnode);
			
			int currwidth = fRootedLevel.MaxWidth();
			if (currwidth < w_min)
			{
				if (fRootedLevel.Depth() > h_max) /* quick exit */
				{
					fStartNode = currnode;
					fEndNode   = -1;
					break;
				}
				else
				{
					w_min    = currwidth;
					fEndNode = currnode;
				}	
			}
		}
	}
	fRootedLevel.MakeRootedLevel(fGraph, fStartNode);	
}

/* relabel */
void CMReLabellerT::NewSequence(void)
{
	/* initialize priorities and queue */
	InitializePriorities();
	
	/* start the queue */
	fCurrLabel = 0;
	Queue(fStartNode);

	/* relabel until queue is empty */	
	int currnode;
	while ( fLabelQueue.PullHighest(currnode) )
		NewNumber(currnode);
	
	/* check that all nodes got renumbered */
	if (fCurrLabel != fPriority.Length()) throw ExceptionT::kGeneralFail;
}

/* initialize priorities as specified by Sloan */
void CMReLabellerT::InitializePriorities(void)
{
	/* need distances to end */
	fRootedLevel.MakeRootedLevel(fGraph, fEndNode);
	
	/* initialize priority and queue */
	fPriority = 0;
	fGraph.InitializePriorities(fPriority, kWeight1);
	fRootedLevel.InitializePriorities(fPriority, kWeight2);
}

/* add a node to the labelling queue */
void CMReLabellerT::Queue(int nodenum)
{
	fPriority[nodenum] += kWeight1;	
	
	/* add to queue */
	if (fStatus[nodenum] == kInActive)
	{
		fLabelQueue.Add(nodenum);
		fStatus[nodenum] = kPreActive;
	}
}

/* assign new number and check adjacent nodes */
void CMReLabellerT::NewNumber(int nodenum)
{	
	const int* adj_i = fGraph.Edges(nodenum);
	int length_j = fGraph.Degree(nodenum);

	/* queue any adjacent nodes */
	if (fStatus[nodenum] == kPreActive)
		for (int j = 0; j < length_j; j++)
			Queue( adj_i[j] );
	
	/* assign new number and status */
	fSequence[fCurrLabel++] = nodenum;  	
	fStatus[nodenum]        = kPostActive;
	
	/* update adjacent nodes */
	for (int j = 0; j < length_j; j++)
		if (fStatus[ adj_i[j] ] == kPreActive)
			MakeActive( adj_i[j] );
}

/* make the specified node active and update adjacent nodes */
void CMReLabellerT::MakeActive(int nodenum)
{
	fPriority[nodenum] += kWeight1;	
	fStatus[nodenum]    = kActive;

	/* queue any adjacent nodes */
	const int* adj_i = fGraph.Edges(nodenum);
	int length_j = fGraph.Degree(nodenum);

	for (int j = 0; j < length_j; j++)
		Queue( adj_i[j] );
}


/* relabel */
void CMReLabellerT::CMSequence(void)
{
	// start nodes at first level
	iArrayT startnodes(fRootedLevel.Width(0));
	fRootedLevel.NodesOnLevel(startnodes, 0);
	for( int ii = 0; ii < startnodes.Length(); ii++)
	{
		//fSequence[startnodes[ii]] = ii + 1;
		fSequence[ii] = startnodes[ii];
	}
	
	// count how many nodes have been renumbered
	int count = startnodes.Length();
	
	AutoArrayT<int> Liordered, Liminus1;
	Liordered = startnodes;
	int widthiminus1 = Liordered.Length();
	
	// loop through levels starting at level 1 (the second)
	AutoArrayT<int> Liactive;
	iArrayT active_tmp;
	for(int i = 1; i < fRootedLevel.Depth(); i++)
	{
		// array of nodes on level i-1
		Liminus1 = Liordered;
		widthiminus1 = Liordered.Length();
		
		// array of active nodes (nodes on level i)
		fRootedLevel.NodesOnLevel(active_tmp, i);
		Liactive = active_tmp;
		
		// array carrying ordering of Liactive
		Liordered.Dimension(0);
		
		// loop through nodes of level i-1
		for(int j = 0; j < widthiminus1; j++)
		{
			// look at edges of nodes in L_(i-1)
			int postactivenode = Liminus1[j];
			iArrayT edges;
			edges.Alias(fGraph.Degree(postactivenode), fGraph.Edges(postactivenode));
			SortByMinDegree(edges);
			
			// while active nodes still exist in L_i, loop through edges of each node
			if(Liactive.Length() != 0)
			{
				// loop through edges of node j and see if they are in L_i
				for(int k = 0; k < edges.Length(); k++)
				{
					int edgenode = edges[k];
					int positioninLi = Liactive.PositionOf(edgenode);
					if( positioninLi != -1 )
					{
						Liactive.DeleteAt(positioninLi);
						AddNode(edgenode, count);
						count++;
						
						Liordered.Append(edgenode);
					}
					
				}
			}
		}
		
		// join disconnected rooted level
		if (Liactive.Length() > 0)
		{
			SortByMinDegree(Liactive);

			for(int k = 0; k < Liactive.Length(); k++)
			{
				AddNode(Liactive[k], count);
				count++;		
				Liordered.Append(Liactive[k]);
			}
		}
	}
}


/* compute bandwidth and profile size */
void CMReLabellerT::ComputeSize(const iArrayT& sequence, int& bandwidth, int& profile)
{
#if __option(extended_errorcheck)
	if (fGraph.NumNodes() != sequence.Length() - 1)
		throw ExceptionT::kSizeMismatch;
#endif

	// create map
	iArrayT nodemap(sequence.Length());
	for (int ii = 0; ii < sequence.Length(); ii++)
		nodemap[sequence[ii]] = ii;

	iArrayT heights(sequence.Length());
	heights = 0;

	bandwidth = 0;
	for (int i = 0; i < fGraph.NumNodes(); i++)
	{
		const int* edges = fGraph.Edges(i);
		int degree = fGraph.Degree(i);
	
		int itag = nodemap[i];
		int maxheight = 0;
		for (int j = 0; j < degree; j++)
		{
			int jtag = nodemap[edges[j]];
			int dist = itag - jtag;
			
			maxheight = (dist > maxheight) ? dist : maxheight;
		}
		
		bandwidth = (maxheight > bandwidth) ? maxheight : bandwidth;
		heights[i] = maxheight;
	}

	/* compute profile */
	profile = heights.Sum() + fGraph.NumNodes();
	
	/* bandwidth definition including diagonal */
	bandwidth++;
}


/* reverse sequence */
void CMReLabellerT::ReverseSequence(void)
{
	//fSequence *= -1;
	//fSequence += fSequence.Length();
	int* ptop = &fSequence[0];
	int* pbottom = &fSequence[fSequence.Length()-1];
	
	while( pbottom > ptop )
	{
		int temp = *pbottom;
		*pbottom = *ptop;
		*ptop = temp;
		ptop++;
		pbottom--;
	}
}


/* add node to map */
void CMReLabellerT::AddNode(int node, int count)
{
	//map ordering: fSequence[node] = count;
	fSequence[count] = node;
}


/* sort edges by minimum degree */
void CMReLabellerT::SortByMinDegree(ArrayT<int>& edges)
{
	int numedges = edges.Length();
	
	// array with degrees of edges
	iArrayT edgedegrees(numedges);
	
	for( int i = 0; i < numedges; i++)
	{
		edgedegrees[i] = fGraph.Degree(edges[i]);
	}
	
	// sort edges by degree in ascending order
	AZ_sort(edgedegrees.Pointer(), numedges, edges.Pointer(), NULL);
}


/* AZ_sort: sorter taken from AZTEC az_sort.c file */
void CMReLabellerT::AZ_sort(int list[], int N, int list2[], double list3[])

/*******************************************************************************

This routine was taken from Knuth: Sorting and Searching. It puts the input
data list into a heap and then sorts it.

Author:          Ray Tuminaro, SNL, 1422
=======

Return code:     double, maximum value in vector 'vec'
============

Parameter list:
===============

list:            On input, values to be sorted. On output, sorted values
(i.e., list[i] <= list[i+1]).

N:               length of vector 'vec'.

list2:           If on input,
a) list2 = NULL: it is unchanged on output,
b) list2 is a list associated with 'list':
on output, if list[k] on input is now element 'j' on output,
list2[j] on output is list2[k].

list3:           If on input,
a) list3 = NULL: it is unchanged on output,
b) list3 is a list associated with 'list':
on output, list3[j] is assigned the input value of list3[k],
if list[j] has been assigned the input value of list[k].

*******************************************************************************/

{

/* local variables */

int    l, r, RR, K, j, i, flag;
int    RR2;
double RR3;

/**************************** execution begins ******************************/

if (N <= 1) return;

l   = N / 2 + 1;
r   = N - 1;
l   = l - 1;
RR  = list[l - 1];
K   = list[l - 1];

if ((list2 != NULL) && (list3 != NULL)) {
RR2 = list2[l - 1];
RR3 = list3[l - 1];
while (r != 0) {
j = l;
flag = 1;

while (flag == 1) {
i = j;
j = j + j;

if (j > r + 1)
flag = 0;
else {
if (j < r + 1)
if (list[j] > list[j - 1]) j = j + 1;

if (list[j - 1] > K) {
list[ i - 1] = list[ j - 1];
list2[i - 1] = list2[j - 1];
list3[i - 1] = list3[j - 1];
}
else {
flag = 0;
}
}
}

list[ i - 1] = RR;
list2[i - 1] = RR2;
list3[i - 1] = RR3;

if (l == 1) {
RR  = list [r];
RR2 = list2[r];
RR3 = list3[r];

K = list[r];
list[r ] = list[0];
list2[r] = list2[0];
list3[r] = list3[0];
r = r - 1;
}
else {
l   = l - 1;
RR  = list[ l - 1];
RR2 = list2[l - 1];
RR3 = list3[l - 1];
K   = list[l - 1];
}
}

list[ 0] = RR;
list2[0] = RR2;
list3[0] = RR3;
}
else if (list2 != NULL) {
RR2 = list2[l - 1];
while (r != 0) {
j = l;
flag = 1;

while (flag == 1) {
i = j;
j = j + j;

if (j > r + 1)
flag = 0;
else {
if (j < r + 1)
if (list[j] > list[j - 1]) j = j + 1;

if (list[j - 1] > K) {
list[ i - 1] = list[ j - 1];
list2[i - 1] = list2[j - 1];
}
else {
flag = 0;
}
}
}

list[ i - 1] = RR;
list2[i - 1] = RR2;

if (l == 1) {
RR  = list [r];
RR2 = list2[r];

K = list[r];
list[r ] = list[0];
list2[r] = list2[0];
r = r - 1;
}
else {
l   = l - 1;
RR  = list[ l - 1];
RR2 = list2[l - 1];
K   = list[l - 1];
}
}

list[ 0] = RR;
list2[0] = RR2;
}
else if (list3 != NULL) {
RR3 = list3[l - 1];
while (r != 0) {
j = l;
flag = 1;

while (flag == 1) {
i = j;
j = j + j;

if (j > r + 1)
flag = 0;
else {
if (j < r + 1)
if (list[j] > list[j - 1]) j = j + 1;

if (list[j - 1] > K) {
list[ i - 1] = list[ j - 1];
list3[i - 1] = list3[j - 1];
}
else {
flag = 0;
}
}
}

list[ i - 1] = RR;
list3[i - 1] = RR3;

if (l == 1) {
RR  = list [r];
RR3 = list3[r];

K = list[r];
list[r ] = list[0];
list3[r] = list3[0];
r = r - 1;
}
else {
l   = l - 1;
RR  = list[ l - 1];
RR3 = list3[l - 1];
K   = list[l - 1];
}
}

list[ 0] = RR;
list3[0] = RR3;

}
else {
while (r != 0) {
j = l;
flag = 1;

while (flag == 1) {
i = j;
j = j + j;

if (j > r + 1)
flag = 0;
else {
if (j < r + 1)
if (list[j] > list[j - 1]) j = j + 1;

if (list[j - 1] > K) {
list[ i - 1] = list[ j - 1];
}
else {
flag = 0;
}
}
}

list[ i - 1] = RR;

if (l == 1) {
RR  = list [r];

K = list[r];
list[r ] = list[0];
r = r - 1;
}
else {
l   = l - 1;
RR  = list[ l - 1];
K   = list[l - 1];
}
}

list[ 0] = RR;
}

} /* AZ_sort */
