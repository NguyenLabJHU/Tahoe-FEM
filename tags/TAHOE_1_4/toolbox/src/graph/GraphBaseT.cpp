/* $Id: GraphBaseT.cpp,v 1.16 2003-11-21 23:32:11 paklein Exp $ */
/* created: paklein (04/13/1999) */
#include "GraphBaseT.h"

#include <iostream.h>
#include <fstream.h>
#include <time.h>

#include "iArrayT.h"
#include "RaggedArray2DT.h"
#include "AutoFill2DT.h"
#include "AutoArrayT.h"
#include "iArray2DT.h"

#ifdef __METIS__
/* partitioning package */
#include "metis.h"
#endif

using namespace Tahoe;

/* rounding floating point numbers */
static inline int rnd(double number) { return int((2.0*number + 1.0)/2); }

/* array behavior */
namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<GraphBaseT*>::fByteCopy = true;
} /* namespace Tahoe */

/* constructor */
GraphBaseT::GraphBaseT(bool verbose):
	fVerbose(verbose),
	fMaxNodeNum(-1),	// uninitialized
	fMinNodeNum(-1),    // uninitialized
	fMinDegree(-1),     // uninitialized
	fMinDegreeNode(-1), // uninitialized
	fShift(0)
{

}

/* destructor */
GraphBaseT::~GraphBaseT(void) { }

/* set adjacency list */
void GraphBaseT::SetAdjacency(const RaggedArray2DT<int>& edge_list)
{
	/* alias */
	fEdgeList.Alias(edge_list);
	
	/* range */
	iArrayT tmp(edge_list.Length(), edge_list.Pointer());
	SetRange(tmp, true);

	/* dimensions */
	fMinDegree = fEdgeList.MinMinorDim(fMinDegreeNode, 0);
}

/* partition: simple partitioning scheme: contract grid until number
	of vertices approaches the partition size and then just assign
	partitions and propagate it back out */
void GraphBaseT::Partition(const iArrayT& config, const iArrayT& weight,
	iArrayT& partition, bool verbose)
{
	/* dimension check */
	if (weight.Length() != fEdgeList.MajorDim()) throw ExceptionT::kSizeMismatch;

	/* initialize partition map */
	partition.Dimension(fEdgeList.MajorDim());
	partition = 0;
	
	/* total number of partitions */
	int num_parts = config.Product();

	/* time */
	clock_t t0 = clock();
	
	/* contract down */
	AutoArrayT<GraphBaseT*> graphs;
	AutoArrayT<iArrayT*> maps;
	AutoArrayT<iArrayT*> weights;
	const GraphBaseT* graph = this;
	int depth = 0;	
	while (graph->NumNodes() > num_parts)
	{
		depth++;
	
		/* next contracted graph */
		GraphBaseT* new_graph = new GraphBaseT(fVerbose);
		if (!new_graph) throw ExceptionT::kOutOfMemory;
		graphs.Append(new_graph);
		
		iArrayT* new_map = new iArrayT();
		if (!new_map) throw ExceptionT::kOutOfMemory;
		maps.Append(new_map);
	
		/* generate contracted graph */	
		new_graph->Contract(*graph, *new_map);
		if (!new_graph->Verify(cout)) throw ExceptionT::kGeneralFail; //TEMP

		/* generated contracted nodal weights */
		iArrayT* new_weight = new iArrayT(new_graph->NumNodes());
		weights.Append(new_weight);
		const iArrayT* last_weight = (depth == 1) ? &weight : weights[depth - 2];
		(*new_weight) = 0;
		for (int i = 0; i < new_map->Length(); i++)
			(*new_weight)[(*new_map)[i]] += (*last_weight)[i];

		/* next */
		graph = new_graph;
	}

	/* time */
	clock_t t1 = clock();
	if (verbose)
		cout << setw(kDoubleWidth) << double(t1 - t0)/CLOCKS_PER_SEC
		     << " sec: GraphBaseT::Partition: contract graph" << endl;
	
	/* initialize partitions */
	int curr_num_parts = graphs[depth-1]->NumNodes();
	partition.Dimension(curr_num_parts);
	partition.SetValueToPosition();

	/* header */
	double* junk = NULL;
	int d_width = OutputWidth(cout, junk);
	if (verbose)
	{
		cout << "    Optimization:\n";
		cout << setw(kIntWidth) << "depth"
		     << setw(kIntWidth) << "reps."
		     << setw(int(1.5*kIntWidth)) << "cuts"
		     << setw(d_width)   << "time" << '\n';
	}

	iArrayT last_part_map;
	for (int i = depth-1; i > -1; i--)
	{	
		const GraphBaseT& next_graph = (i == 0) ? *this : *graphs[i-1];
		const iArrayT& map = *maps[i];
		const iArrayT& next_weight = (i == 0) ? weight: *weights[i-1];
		if (map.Length() != next_graph.NumNodes()) throw ExceptionT::kSizeMismatch; //TEMP
		
		int dim = map.Length();
		last_part_map.Swap(partition);
		partition.Dimension(map.Length());
	
		/* propagate partitions out */
		for (int j = 0; j < dim; j++)
			partition[j] = last_part_map[map[j]];
			
		/* optimize partitions */
		if (1 || i == 0)
		{			
			int reps, cuts;
			clock_t ta = clock();
			OptimizeParts(next_graph, next_weight, partition, num_parts, reps, cuts);
			clock_t tb = clock();
			if (verbose)
				cout << setw(kIntWidth) << i
			    	 << setw(kIntWidth) << reps
			    	 << setw(int(1.5*kIntWidth)) << cuts
			    	 << setw(d_width)   << double(tb - ta)/CLOCKS_PER_SEC << endl;

			if (fVerbose) cout << "cuts = " << cuts << endl; //TEMP			
		}
	}

	/* time */
	clock_t t2 = clock();
	if (verbose)
		cout << setw(kDoubleWidth) << double(t2 - t1)/CLOCKS_PER_SEC
		     << " sec: GraphBaseT::Partition: balance partitions" << endl;
	
	/* free memory */
	for (int j = 0; j < graphs.Length(); j++)
	{
		delete graphs[j];
		delete maps[j];
		delete weights[j];
	}
}	

/* generate partition using METIS */
void GraphBaseT::Partition_METIS(int num_partitions, const iArrayT& weight,
	iArrayT& partition, int volume_or_edgecut)
{
#ifndef __METIS__
#pragma unused(num_partitions)
#pragma unused(weight)
#pragma unused(partition)
#pragma unused(volume_or_edgecut)

	/* error message */
	cout << "\n GraphBaseT::Partition_METIS: requires metis module" << endl;
	throw ExceptionT::kGeneralFail;
#else
	
	/* NOTE: based on "kmetis.c" */

	/* dimension check */
	if (weight.Length() != fEdgeList.MajorDim()) throw ExceptionT::kSizeMismatch;

	/* options check */
	if (volume_or_edgecut != 0 && volume_or_edgecut != 1)
	{
		cout << "\n GraphBaseT::Partition_METIS: volume_or_edgecut must be 0 or 1: " 
		     << volume_or_edgecut << endl;
		throw ExceptionT::kGeneralFail;
	}

	/* initialize partition map */
	partition.Dimension(fEdgeList.MajorDim());
	partition = 0;
	if (num_partitions < 2) return;

	/* timing info */
	timer TOTALTmr, METISTmr;    
  	cleartimer(TOTALTmr);
  	cleartimer(METISTmr);

	starttimer(TOTALTmr);
	cout << "**********************************************************************\n";
	cout << METISTITLE;
	cout << "Graph Information ---------------------------------------------------\n";
	cout << "#Vertices: " << fEdgeList.MajorDim() << '\n';
	cout << "   #Edges: " << fEdgeList.Length()/2 << '\n';
	cout << "   #Parts: " << num_partitions << '\n';
	cout << "  Balancing Constraints: " << 0 << '\n';
	cout << "\nK-way Partitioning... -----------------------------------------------\n" << endl;

	/* offset vector for adjacency list */
	iArrayT offsets;
	fEdgeList.GenerateOffsetVector(offsets);
	
	/* partition weights - even */
	ArrayT<float> tpwgts(num_partitions);
	tpwgts = 1.0/num_partitions;

  	/* METIS options */
	iArrayT options(5);
	options[0] = 0; /* use all default values */

	/* other arguments */
	int num_vertices = fEdgeList.MajorDim();
	int edgecut; /* returns with number of edges cut by the partitioning */
	int num_flag = 0;
	int weight_flag = 2;
	int* adjwgt = NULL;

	/* partitioning method */
	starttimer(METISTmr);
	if (volume_or_edgecut == 0)
		METIS_WPartGraphVKway(&num_vertices, offsets.Pointer(), fEdgeList.Pointer(), (int*) weight.Pointer(), 
			adjwgt, &weight_flag, &num_flag, &num_partitions, 
			tpwgts.Pointer(), options.Pointer(), &edgecut, partition.Pointer());
	else if (volume_or_edgecut == 1)
		METIS_WPartGraphKway(&num_vertices, offsets.Pointer(), fEdgeList.Pointer(), (int*) weight.Pointer(), 
			adjwgt, &weight_flag, &num_flag, &num_partitions, 
			tpwgts.Pointer(), options.Pointer(), &edgecut, partition.Pointer());
	else throw;
	stoptimer(METISTmr);

	/* assess partition quality */
	GraphType* graph = CreateGraph();
	int ncon = 1; /* just 1 constraint per vertex - the weight */
	SetUpGraph(graph, OP_KMETIS, num_vertices, ncon, offsets.Pointer(), fEdgeList.Pointer(), 
		(int*) weight.Pointer(), NULL, 2);	
	ComputePartitionInfo(graph, num_partitions, partition.Pointer());
  	FreeGraph(graph);

	/* write timing info */
	stoptimer(TOTALTmr);
	cout << "\nTiming Information --------------------------------------------------\n";
	cout << "  Partitioning: \t\t " << gettimer(METISTmr) << "   (KMETIS time)\n";
	cout << "  Total:        \t\t " << gettimer(TOTALTmr) << '\n';
	cout << "**********************************************************************\n";
	cout.flush();
#endif
}	

/* fill in the degrees for the specified nodes */
void GraphBaseT::ReturnDegrees(const ArrayT<int>& nodes, ArrayT<int>& degrees) const
{
#if __option(extended_errorcheck)
	if (nodes.Length() != degrees.Length()) throw ExceptionT::kSizeMismatch;
#endif

	const int* pnodes = nodes.Pointer();
	int* pdegrees = degrees.Pointer();
	int length = nodes.Length();
	for (int i = 0; i < length; i++)
		*pdegrees++ = fEdgeList.MinorDim(*pnodes++ - fShift); 	
}

/* verify graph data is consistent - returns 1 if consistent, 0 otherwise */
int GraphBaseT::Verify(ostream& err) const
{
	int OK = 1;
	iArrayT row, row_j;
	for (int i = 0; i < fEdgeList.MajorDim(); i++)
	{
		fEdgeList.RowAlias(i, row);
		for (int j = 0; j < row.Length(); j++)
		{
			fEdgeList.RowAlias(row[j] - fShift, row_j);
			if (!row_j.HasValue(i + fShift))
			{
				err << " {node,miss} = {" << setw(kIntWidth) << row[j] + 1;
				err << "," << setw(kIntWidth) << i + fShift + 1 << "}" << endl;
				OK = 0;
			}
		}
	}
	return OK;
}

/* output graph to stream */
void GraphBaseT::Write(ostream& out) const
{
	iArrayT temp;
	for (int i = 0; i < fEdgeList.MajorDim(); i++)
	{
		temp.Alias(fEdgeList.MinorDim(i), fEdgeList(i));
		
		out << setw(kIntWidth) << fShift + i;
		temp++;
		out << temp.wrap(10, 8) << '\n';
		temp--;
	}
}

/* initialize the priorities as specified in Sloan.  Assumes DIRECT
* INDEXING in priorities */
void GraphBaseT::InitializePriorities(iArrayT& priorities, int W1) const
{
	//TEMP
	if (fShift != 0)
	{
		cout << "\n GraphBaseT::InitializePriorities: not expecting non-zero\n"
		     <<   "     node number shift: " << fShift << endl;
		throw ExceptionT::kGeneralFail;
	}

	for (int i = 0; i < fMaxNodeNum; i++)
		priorities[i] += (fMaxNodeNum - Degree(i))*W1;
}

/* set the node range and shift parameters */
void GraphBaseT::SetRange(const nArrayT<int>& node_list, bool reset)
{
	if (node_list.Length() > 0)
	{
		int min, max;
		node_list.MinMax(min, max);
		
		if (reset)
		{
			fMaxNodeNum = max;
			fMinNodeNum = min;
		}
		else
		{	
			fMaxNodeNum = (max > fMaxNodeNum) ? max : fMaxNodeNum;
			fMinNodeNum = (fMinNodeNum == -1 || min < fMinNodeNum) ? min : fMinNodeNum;
		}
		fShift = fMinNodeNum;
	}
}

/************************************************************************
* Private
************************************************************************/

/* contract the graph */
void GraphBaseT::Contract(const GraphBaseT& parent, iArrayT& map)
{
	//TEMP
	if (fShift != 0 && parent.fShift)
	{
		cout << "\n GraphBaseT::Contract: not expecting a node number shift\n"
		     <<   "     in this graph " << fShift << " or the parent "
		     << parent.fShift << endl;
		throw ExceptionT::kGeneralFail;
	}

	/* dimensions */
	int dim = parent.NumNodes();

	/* temp work space */
	AutoFill2DT<int> edgedata(dim, 1, 10, parent.MinDegree()); // under allocate
	AutoArrayT<int> degrees;
	iArrayT edges, edges2;

	/* set sequence by min to max degree */
	iArrayT all_degrees(dim), sequence(dim);
	int row_shift;
	parent.Degrees(all_degrees, row_shift);
	sequence.SetValueToPosition();
	AZ_sort(all_degrees.Pointer(), dim, sequence.Pointer(), NULL);
	
	/* mark vertices to keep/merge */
	iArrayT flags(dim);
	flags = -2; // mark all as active
	for (int i = 0; i < dim; i++)
	{
		int next = sequence[i]; // node number
		if (flags[next] == -2)
		{	
			/* get edges */
			parent.GetEdges(next, edges);
				
			/* get degrees for edges adjacent vertex */
			degrees.Dimension(edges.Length());
			parent.ReturnDegrees(edges, degrees);		
		
			/* mark edge to collapse */
			int node = SelectCollapse(edges, degrees, flags);
			if (node > -1)
			{
				flags[next] = node; // vertex to collapse
				flags[node] = -1;   // won't be retained
			}
			else
				flags[next] = next; // retain, but no collapse
		}
	}

	/* generate map map[old_num] = new_num */
	map.Dimension(dim);
	map = -1;
	int dex = 0;
	for (int j = 0; j < dim; j++) /* retained nodes */
		if (flags[j] != -1)
		{
			map[j] = dex; // retained nodes
			if (flags[j] != j)
				map[flags[j]] = dex; // collapsed nodes
			dex++;
		}

	/* perform collapse */
	for (int k = 0; k < dim; k++)
		if (flags[k] != -1)
		{
			/* get edges */
			parent.GetEdges(k, edges);

			/* retained vertices from self */
			for (int i = 0; i < edges.Length(); i++)
			{
				int dex = edges[i];
				if (flags[k] != dex) // skip duplicate
					edgedata.AppendUnique(k, map[dex]);
			}
			
			/* append retained vertices from node (not self, no repeats) */
			parent.GetEdges(flags[k], edges);
			for (int j = 0; j < edges.Length(); j++)
			{
				int dex = edges[j];
				if (dex != k) // skip duplicate
					edgedata.AppendUnique(k, map[dex]);
			}
		}
	
	/* compress/copy in */
	fEdgeList.CopyCompressed(edgedata);

	/* set graph parameters */
	fMinNodeNum = 0;
	fMaxNodeNum = fEdgeList.MajorDim();
	fMinDegree  = fEdgeList.MinMinorDim(fMinDegreeNode, 0);
}

/* return the "best" edge to collapse, -1 if none of adjacent nodes are still
* active */
int GraphBaseT::SelectCollapse(const ArrayT<int>& edges, const ArrayT<int>& degrees,
	const ArrayT<int>& flags) const
{
	//TEMP
	if (fShift != 0)
	{
		cout << "\n GraphBaseT::SelectCollapse: not expecting non-zero\n"
		     <<   "     node number shift: " << fShift << endl;
		throw ExceptionT::kGeneralFail;
	}

	int node = -1, degree = -1;
	for (int i = 0; i < edges.Length(); i++)
	{
		int node_i = edges[i];
		if (flags[node_i] == -2) // still in active set
		{
			if (node == -1 || (node > -1 && degrees[i] > degree))
			{
				node   = node_i;
				degree = degrees[i];
			}
		}
	}
	
	return node;
}

/* set gains */
int GraphBaseT::SetGains(const GraphBaseT& graph, const iArrayT& partition,
	iArray2DT& gain) const
{
	//TEMP
	if (fShift != 0)
	{
		cout << "\n GraphBaseT::SetGains: not expecting non-zero\n"
		     <<   "     node number shift: " << fShift << endl;
		throw ExceptionT::kGeneralFail;
	}

	/* dimensions */
	int nnd = gain.MajorDim();
	int dim = gain.MinorDim();

	/* connections by partition */
	gain = 0;
	for (int i = 0; i < nnd; i++)
	{
		int degree = graph.Degree(i);
		const int* edges = graph.Edges(i);
		int* gain_i = gain(i);
		for (int j = 0; j < degree; j++)
			gain_i[partition[*edges++]]++;
	}

	/* set gains and compute cut edges/select candidates */
	int cuts = 0;
	for (int k = 0; k < nnd; k++)
	{
		int* gain_k = gain(k);
		int  part_k = partition[k];
		int   dim_k = gain_k[part_k];
		for (int j = 0; j < dim; j++)
		{
			/* compute cuts */
			if (j != part_k) cuts += gain_k[j];		
		
			/* connects -> gain */
			gain_k[j] -= dim_k;
				// this gain does not include a
				// size term. don't know how to
				// update the size term yet. size
				// driving force for switching partitions
				// is applied when selecting moves.
		}
	}
	
	return cuts;
}	

/* map of destination partitions for each node */
void GraphBaseT::SetMoves(const iArray2DT& gain, const iArrayT& partition,
	const iArrayT& size, iArrayT& move) const
{
	//TEMP
	if (fShift != 0)
	{
		cout << "\n GraphBaseT::SetMoves: not expecting non-zero\n"
		     <<   "     node number shift: " << fShift << endl;
		throw ExceptionT::kGeneralFail;
	}

	/* dimensions */
	int nnd = move.Length();
	int dim = size.Length();

	/* select candidates */
	iArrayT max(nnd);
	move =-1;
	max  =-1;
	for (int k = 0; k < nnd; k++)
	{
		const int* gain_k = gain(k);
		int part_k = partition[k];
//		int   dim_k = gain_k[part_k];
		for (int j = 0; j < dim; j++)
		{
			int gain_kj = gain_k[j];
						
			/* store max gain for allowable move */
			if (part_k != j)
			{
				/* add size balancing driving force */
				// int size_diff = size[part_k] - size[j];
				int size_diff;
				if (size[j] == 0)
					size_diff = size[part_k];
				else
					size_diff = rnd(double(size[part_k])/double(size[j])) - 1;

				int size_weight = 1;
				size_diff *= size_weight;
				gain_kj += size_diff;
				if (gain_kj > max[k])
				if (size[part_k] >= size[j] && gain_kj > max[k])
				{
					 max[k] = gain_kj;
					move[k] = j;
				}
				/* randomness */
				//else if (double(rand())/double(RAND_MAX) > 0.5 && gain_kj > 0)
				//{
				//	 max[k] = gain_kj;
				//	move[k] = j;
				//}				
				//else if (size[part_k] < size[j])
				//{	
				//	max[k] = gain_kj;
				//	move[k] = j;
				//}
			}
		}
	}
}

/* apply moves/update data */
int GraphBaseT::ApplyMoves(const GraphBaseT& graph, const iArrayT& move,
	const iArrayT& weight, iArrayT& partition, iArray2DT& gain, iArrayT& size)
{
	//TEMP
	if (fShift != 0)
	{
		cout << "\n GraphBaseT::ApplyMoves: not expecting non-zero\n"
		     <<   "     node number shift: " << fShift << endl;
		throw ExceptionT::kGeneralFail;
	}

	/* dimensions */
	int nnd = move.Length();
	int dim = size.Length();

	/* change in cuts */
	int cuts = 0;

	/* collect moves */
	int count = nnd - move.Count(-1);
	if (count > 0)
	{
		iArrayT topgains(count);
		iArrayT topnodes(count);
		int dex = 0;
		for (int m = 0; m < nnd; m++)
			if (move[m] > -1)
			{
				topgains[dex] = gain(m, move[m]);
				topnodes[dex] = m;
				dex++;
			}

		/* sort ascending */
		AZ_sort(topgains.Pointer(), count, topnodes.Pointer(), NULL);
	
		/* move highest gains first */
//		int num_moves = count/2;
		int num_moves = count;
		dex = count - 1;
		for (int l = 0; l < num_moves; l++)
		{
			/* candidate node */
			int node = topnodes[dex--];

			/* affected paritions */
			int& part_i = partition[node];
			int  part_f = move[node];

			/* allowable move */
			int* gain_l = gain(node);
//			if (size[part_i] > size[part_f] && gain_l[part_f] > 0)
			if (size[part_i] > 1 && size[part_i] >= size[part_f])
//			double r = (size[part_f] > 0) ? double(size[part_i])/double(size[part_f]) : 5.0;
//			if (size[part_i] > 1 && r > 0.99)
			{
				/* global data */
				size[part_i] -= weight[node];
				size[part_f] += weight[node];

				cuts -= gain_l[part_f];
		
				/* local data */
				int degree = graph.Degree(node);
				const int* edges = graph.Edges(node);
				for (int i = 0; i < degree; i++)
				{
					int nd = *edges++;
					int* gain_i = gain(nd);

					/* transfer */
					gain_i[part_i]--;
					gain_i[part_f]++;
					
					/* reset to gain */
					int self = partition[nd];
					if (part_i == self || part_f == self)
					{
						int shift = gain_i[self];
						for (int j = 0; j < dim; j++)
							gain_i[j] -= shift;
					}
				}

				/* self */
				int shift = gain_l[part_f];
				for (int j = 0; j < dim; j++)
					gain_l[j] -= shift;
				
				/* update partition */
				part_i = part_f;		
			}
		}
	}
	
	return cuts;
}

/* optimize partitions */
void GraphBaseT::OptimizeParts(const GraphBaseT& graph, const iArrayT& weight,
	iArrayT& partition, int dim, int& repetitions, int& cuts)	
{
#if __option(extended_errorcheck)
	if (graph.NumNodes() != partition.Length()) throw ExceptionT::kSizeMismatch;
	if (graph.NumNodes() != weight.Length()) throw ExceptionT::kSizeMismatch;
#endif	

	//TEMP
	if (fShift != 0)
	{
		cout << "\n GraphBaseT::OptimizeParts: not expecting non-zero\n"
		     <<   "     node number shift: " << fShift << endl;
		throw ExceptionT::kGeneralFail;
	}

	/* dimensions */
	int nnd = graph.NumNodes();
	if (fVerbose) cout << "num_nodes = " << nnd << endl;

	/* compute partition weighted sizes */
	iArrayT size(dim);
	size = 0;
	for (int i = 0; i < nnd; i++)
		size[partition[i]] += weight[i];
		//NOTE: this measure of the size is based on the collapsed partition
		//      size not the actual partition size.

	/* assure specified number of partitions */
//NOTE: this did not assure that all partitions have non-zero size.
//      SetMoves now applies a size-balancing criteria that tends to
//      fill empty partitions. PAK (06/06/2000)
#if 0
	int curr_dim = dim - size.Count(0);
	while (curr_dim < dim && nnd > curr_dim)
	{
		int max_size = size.Max();
		int max_part;
		size.HasValue(max_size, max_part);
		
		/* use 1st occurrence */
		int seed_node;
		partition.HasValue(max_part, seed_node);
		
		/* empty part */
		int non_part;
		size.HasValue(0, non_part);
		
		/* seed part */
		size[non_part] += weight[seed_node];
		size[max_part] -= weight[seed_node];
		partition[seed_node] = non_part;
		curr_dim++;
	}	
#endif

	
	/* initialize gains */
	iArray2DT gain(nnd, dim);
	cuts = SetGains(graph, partition, gain);
	if (fVerbose) cout << " initial cuts = " << cuts << '\n';

	/* select candidates */
	iArrayT move(nnd);
	SetMoves(gain, partition, size, move);

	/* collect moves */
	int dcuts = ApplyMoves(graph, move, weight, partition, gain, size);
	if (fVerbose) cout << " dcuts = " << dcuts << '\n';

	repetitions = 0;
	int zero_count = 0; // probably twice zero means nothing will happen
	//	int dcuts_last = 1;
	double rcuts, rcuts_last = 1.0;
	int max_reps = 50;
	int max_zero_count = 5;
	while (zero_count < max_zero_count && repetitions < max_reps)
	{
		repetitions++;
		cuts += dcuts;
		
		SetMoves(gain, partition, size, move);
//############################################
#if 0
		dcuts_last = dcuts;
#endif
//############################################
		dcuts = ApplyMoves(graph, move, weight, partition, gain, size);		
		if (fVerbose) cout << " dcuts = " << dcuts << '\n';
		
//############################################
#if 0
		/* absolute change criterion */
		if (dcuts_last == 0 && dcuts == 0)
			zero_count++;
		else
			zero_count = 0;
#endif
//############################################

//############################################
#if 1
		/* relative change criterion */
		rcuts = double(dcuts)/(cuts + 1.0);
		if (fabs(rcuts_last) < 0.001 && fabs(rcuts) < 0.001)
			zero_count++;
		else
			zero_count = 0;
		rcuts_last = rcuts;
#endif
//############################################
	}

	if (fVerbose)
	{
		cout << " reps = " << repetitions << '\n';
		cout << " size = ";
		cout << size.wrap(10) << '\n';
		cout << endl;
	}

// QUESTIONS:
// (1) support arbitrary partition numbers?
// (2) what happens if dim is wrong?
//		-> empty partitions will stay empty b/c there won't
//      -> be any driving force to fill them
// (3) correct number of partitions here?
// while end condition not satisfied
	// compute partition sizes
	// compute nodal deltas	
	// select nodes to move
	// move nodes/update sizes (and deltas?)
}

/* AZ_sort: sorter taken from AZTEC az_sort.c file */
void GraphBaseT::AZ_sort(int list[], int N, int list2[], double list3[])

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
