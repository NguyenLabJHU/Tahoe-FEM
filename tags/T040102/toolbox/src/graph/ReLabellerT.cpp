/* $Id: ReLabellerT.cpp,v 1.1.1.1 2001-01-25 20:56:27 paklein Exp $ */
/* created: paklein (08/05/1996)                                          */

#include "ReLabellerT.h"
#include <limits.h>
#include <string.h>
#include <fstream.h>
#include "iArray2DT.h"
#include "LocalParabolaT.h"
#include "AutoArrayT.h"

/* status codes */
const int kInActive   = 0;
const int kPreActive  = 1;
const int kActive     = 2;
const int kPostActive = 3;

/* constructor */
ReLabellerT::ReLabellerT(void): fLabelQueue(fPriority),
	fStartNode(0),
	fEndNode(0),
	fCurrLabel(0)
{

}

/* add a group to the graph */
void ReLabellerT::AddGroup(const iArray2DT& moredata)
{
	fGraph.AddGroup(moredata);
}

void ReLabellerT::AddGroup(const RaggedArray2DT<int>& moredata)
{
	fGraph.AddGroup(moredata);
}

/* renumber the oldsequence */
int ReLabellerT::Renumber(iArray2DT& oldsequence)
{
	/* make graph */
	fGraph.MakeGraph();

	/* initialize */
	Initialize();

	/* pick start and end nodes */
	SelectNodes();

	/* relabel */
	NewSequence();
	
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
int ReLabellerT::Renumber(ArrayT<iArray2DT*>& oldsequences)
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
	fSequence.SetValueToPosition();

	int bandwidth, profile;
	ComputeSize(fSequence, bandwidth, profile);
	cout << " before:\n";
	cout << "   bandwidth = " << bandwidth << '\n';
	cout << "     profile = " << profile   << '\n';
	cout << endl;

	/* find optimal parameters with bisection */
//	ofstream itout("relab.out");
	AutoArrayT<int> W1_list;
	AutoArrayT<int> profile_list;
	fWeight2 = 5;

	int W1L, W1U, profileL, profileU;

	fWeight1 = W1L = 0;
	NewSequence();
	ComputeSize(fSequence, bandwidth, profileL);
	//itout << fWeight1 << '\t' << profileL << '\n';
	W1_list.Append(fWeight1);
	profile_list.Append(profileL);

	fWeight1 = W1U = 100;
//	fWeight1 = W1U = 10;
	NewSequence();
	ComputeSize(fSequence, bandwidth, profileU);
	//itout << fWeight1 << '\t' << profileU << '\n';
	W1_list.Append(fWeight1);
	profile_list.Append(profileU);
	
	int newW1 = (W1L + W1U)/2;
	int count = 0;
	int newprofile;
	while ((W1L != newW1 || W1U != newW1) && count++ < 20)
	{
		fWeight1 = newW1;
		int dex = W1_list.PositionOf(fWeight1);
		if (dex > -1)
			newprofile = profile_list[dex];
		else
		{
			NewSequence();
			ComputeSize(fSequence, bandwidth, newprofile);
			//itout << fWeight1 << '\t' << newprofile << '\n';
			W1_list.Append(fWeight1);
			profile_list.Append(newprofile);
		}
	
		/* min lies outside of range */
		if (profileL > newprofile && newprofile > profileU)
		{
			profileL = newprofile;
			W1L      = fWeight1;
			fWeight1 = W1U = W1U + (W1U - fWeight1);

			int dex = W1_list.PositionOf(fWeight1);
			if (dex > -1)
				profileU = profile_list[dex];
			else
			{
				NewSequence();
				ComputeSize(fSequence, bandwidth, profileU);
				//itout << fWeight1 << '\t' << profileU << '\n';
				W1_list.Append(fWeight1);
				profile_list.Append(profileU);
			}
		}
		else if (profileL < newprofile && newprofile < profileU)
		{
			profileU = newprofile;
			W1U      = fWeight1;			
			fWeight1 = W1L = W1L - (fWeight1 - W1L);
			fWeight1 = (fWeight1 < 0) ? 0 : fWeight1;
					
			int dex = W1_list.PositionOf(fWeight1);
			if (dex > -1)
				profileL = profile_list[dex];
			else
			{
				NewSequence();
				ComputeSize(fSequence, bandwidth, profileL);		
				//itout << fWeight1 << '\t' << profileL << '\n';
				W1_list.Append(fWeight1);
				profile_list.Append(profileL);
			}
		}
		else if (profileL < newprofile && newprofile > profileU)
		{
			if (profileL < profileU)
			{
				profileU = newprofile;
				W1U      = fWeight1;
				fWeight1 = W1L = W1L/2;		

				int dex = W1_list.PositionOf(fWeight1);
				if (dex > -1)
					profileL = profile_list[dex];
				else
				{
					NewSequence();
					ComputeSize(fSequence, bandwidth, profileL);		
					//itout << fWeight1 << '\t' << profileL << '\n';
					W1_list.Append(fWeight1);
					profile_list.Append(profileL);
				}
			}
			else
			{
				profileL = newprofile;
				W1L      = fWeight1;	
				fWeight1 = W1U = 2*W1U;
				
				int dex = W1_list.PositionOf(fWeight1);
				if (dex > -1)
					profileU = profile_list[dex];
				else
				{
					NewSequence();
					ComputeSize(fSequence, bandwidth, profileU);
					//itout << fWeight1 << '\t' << profileU << '\n';
					W1_list.Append(fWeight1);
					profile_list.Append(profileU);				
				}
			}
		}
		else
		{
			if (newprofile < profileL)
			{
				profileL = newprofile;
				W1L      = fWeight1;
			}
			else
			{
				profileU = newprofile;
				W1U      = fWeight1;
			}
		}
	
		newW1    = (W1L + W1U)/2;
				
		cout << setw(kIntWidth) << W1L;
		cout << setw(kIntWidth) << W1U;
		cout << setw(kIntWidth) << newW1 << endl;
	}

	/* print results */
	cout << "\n after: " << count << " iterations\n";
	cout << "   bandwidth = " << bandwidth  << '\n';
	cout << "     profile = " << newprofile << '\n';
	cout << endl;	

//TEMP - old default parameters
//	fWeight1 = 2;
//	fWeight2 = 1;
//	NewSequence();
//	ComputeSize(fSequence, bandwidth, newprofile);
//	cout << "\n after: " << count << " iterations\n";
//	cout << "   bandwidth = " << bandwidth  << '\n';
//	cout << "     profile = " << newprofile << '\n';
//	cout << endl;	

//TEMP: output profile(w1)
//	int dW = 1;
//	for (fWeight1 = 0; fWeight1 < 500; fWeight1 += dW)
//	{
//		int bandwidth, profile;
//	
//		NewSequence();
//		ComputeSize(fSequence, bandwidth, profile);		
//		cout << setw(kIntWidth) << fWeight1;
//		cout << setw(kIntWidth) << bandwidth;
//		cout << setw(kIntWidth) << profile << endl;
//	
//		if (fWeight1 > 20)  dW = 5;
//		if (fWeight1 > 120) dW = 10;
//	}
//	throw;

	/* re-sequence positive values */	
	int	label = 0;
	for (int i = 0; i < fSequence.Length(); i++)
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
void ReLabellerT::Initialize(void)
{
	int numnodes = fGraph.NumNodes();

	/* allocate space */
	fSequence.Allocate(numnodes);
	fStatus.Allocate(numnodes);
	fPriority.Allocate(numnodes),

	/* initialize */
	fSequence =-1;
	fStatus   =-1;
	fPriority =-1;
}

/* fill internal rooted level even for disconnected graphs with
* the root for each branch selected by minimum width */
void ReLabellerT::BuildRootedLevel(void)
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
	if (count != nnd) throw eGeneralFail;
	
	/* set end node from last level */
	fRootedLevel.NodesOnLevel(nodes_used, fRootedLevel.Depth() - 1);
	fEndNode = nodes_used[0];
}

/* select pseudo-diameter nodes */
void ReLabellerT::SelectNodes(void)
{
	fStartNode = fGraph.MinDegreeNode();
	fRootedLevel.MakeRootedLevel(fGraph, fStartNode);
	
	iArrayT	topnodes, degrees;
	fEndNode = -1;
	
	while (fEndNode == -1)
	{
		int h_max = fRootedLevel.Depth();
		int	w_min = INT_MAX;

		/* collect top level info */
		fRootedLevel.NodesOnLevel(topnodes, h_max - 1);
		degrees.Allocate(topnodes.Length()); //really want to Allocate() every time?
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

	/* set final rooted level structure */
	if (fRootedLevel.RootNode() != fStartNode)
		fRootedLevel.MakeRootedLevel(fGraph, fStartNode);
}

/* relabel */
void ReLabellerT::NewSequence(void)
{
	/* reset status */
	fStatus = 0;

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
	if (fCurrLabel != fPriority.Length()) throw eGeneralFail;
}

/* estimate size of the system in terms of the
* bandwidth and profile using the sequence */
void ReLabellerT::ComputeSize(const iArrayT& sequence, int& bandwidth, int& profile)
{
#if __option(extended_errorcheck)
	if (fGraph.NumNodes() != sequence.Length())
		throw eSizeMismatch;
#endif

	/* generate node map */
	iArrayT nodemap(sequence.Length());
	for (int ii = 0; ii < fSequence.Length(); ii++)
		nodemap[sequence[ii]] = ii;

	iArrayT heights(sequence.Length());
	heights = 0;

	bandwidth = 0;
	for (int i = 0; i < fGraph.NumNodes(); i++)
	{
		int* edges  = fGraph.Edges(i);
		int  degree = fGraph.Degree(i);
	
		int itag = nodemap[i];
		int maxheight = 0;
		for (int j = 0; j < degree; j++)
		{
			int jtag = nodemap[edges[j]];
			int dist = itag - jtag;
			
			maxheight = (dist > maxheight) ? dist : maxheight;
		}
		
		bandwidth = (maxheight > bandwidth) ? maxheight : bandwidth;
		heights[itag] = maxheight;
	}

	/* compute profile final values */
	bandwidth += 1;
	profile = heights.Sum() + fGraph.NumNodes();
}

void ReLabellerT::ReverseSequence(void)
{
	int* ptop = fSequence.Pointer();
	int* pbot = fSequence.Pointer(fSequence.Length() - 1);
	
	while (pbot > ptop)
	{
		int temp = *pbot;
		
		*pbot = *ptop;
		*ptop = temp;

		ptop++; pbot--;
	}
}

/* initialize priorities as specified by Sloan */
void ReLabellerT::InitializePriorities(void)
{
	/* need distances to end */
	fRootedLevel.MakeRootedLevel(fGraph, fEndNode);
	
	/* initialize priority and queue */
	fPriority = 0;
	fGraph.InitializePriorities(fPriority, fWeight1);
	fRootedLevel.InitializePriorities(fPriority, fWeight2);
}

/* add a node to the labelling queue */
void ReLabellerT::Queue(int nodenum)
{
	fPriority[nodenum] += fWeight1;	
	
	/* add to queue */
	if (fStatus[nodenum] == kInActive)
	{
		fLabelQueue.Add(nodenum);
		fStatus[nodenum] = kPreActive;
	}
}

/* assign new number and check adjacent nodes */
void ReLabellerT::NewNumber(int nodenum)
{	
	int* adj_i    = fGraph.Edges(nodenum);
	int  length_j = fGraph.Degree(nodenum);

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
void ReLabellerT::MakeActive(int nodenum)
{
	fPriority[nodenum] += fWeight1;	
	fStatus[nodenum]    = kActive;

	/* queue any adjacent nodes */
	int* adj_i    = fGraph.Edges(nodenum);
	int  length_j = fGraph.Degree(nodenum);

	for (int j = 0; j < length_j; j++)
		Queue( adj_i[j] );
}
