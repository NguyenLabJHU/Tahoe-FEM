/* $Id: GraphT.cpp,v 1.1.1.1 2001-01-25 20:56:27 paklein Exp $ */
/* created: paklein (08/05/1996)                                          */

#include "GraphT.h"

#include <time.h>

#include "iArrayT.h"
#include "iArray2DT.h"
#include "RaggedArray2DT.h"
#include "AutoFill2DT.h"
#include "AutoArrayT.h"
#include "RootedLevelT.h"
#include "PartitionT.h"

/* inlines */
/* return true if node is in range and active */
inline bool GraphT::Active(int node, const iArrayT& active) const
{
	return node > -1 &&
	       node < active.Length() &&
	       active[node] == 1;
}

/* constructor */
GraphT::GraphT(bool verbose): GraphBaseT(verbose) { }

/* add a group to the graph */
void GraphT::AddGroup(const iArray2DT& groupdata)
{
	/* reset maximum node number */
	SetRange(groupdata);

	/* do not allow repeated registering of groups */
	if (!fGroupData_1.AppendUnique(&groupdata)) throw eGeneralFail;
}

void GraphT::AddGroup(const RaggedArray2DT<int>& groupdata)
{
	/* reset maximum node number*/
	iArrayT temp(groupdata.Length(), groupdata.Pointer());
	SetRange(temp);

	/* do not allow repeated registering of groups */
	if (!fGroupData_2.AppendUnique(&groupdata)) throw eGeneralFail;
}

void GraphT::ClearGroups(void)
{
	fGroupData_1.Clear();
	fGroupData_2.Clear();
}
	
/* make the graph using the current data */
/* NOTE: assumes all ien >= 0, i.e. these are connectivities */
void GraphT::MakeGraph(void)
{
	if (fMaxNodeNum == -1) return;

	/* this version does not support any shifted nodes, but
	 * does skip ien < 0 */
	if (fShift > 0)
	{
		cout << " GraphT::MakeGraph: active connectivities must be >= 0, mininum\n"
		     << "     found " << fShift << endl;
		throw eGeneralFail;
	}
	/* ignore ien < 0 */
	else if (fShift < 0)
		fShift = fMinNodeNum = 0;

	/* temp work space */
	int range = fMaxNodeNum - fMinNodeNum + 1;
	AutoFill2DT<int> edgedata(range, 25, 5);

	/* create adjacency lists */	
	const iArray2DT* currgroup;
	fGroupData_1.Top();
	while ( fGroupData_1.Next(currgroup) )
	{
		int  nel = currgroup->MajorDim();
		int  nen = currgroup->MinorDim();
		int* ien = currgroup->Pointer();
		for (int k = 0; k < nel; k++)
		{
			for (int i = 0; i < nen; i++)
			{	
				int r_i = ien[i];
				if (r_i > -1)
				{
					for (int j = i+1; j < nen; j++)
					{
						int r_j = ien[j];
						if (r_j > -1 && r_i != r_j && // no refs to self
						    edgedata.AppendUnique(r_j - fShift, r_i))
							edgedata.Append(r_i - fShift, r_j);
							// NOTE: could also AppendUnique to the shorter of
							//       r_i/r_j, but didn't seem to be faster
					}
				}
			}
			ien += nen;
		}
	}

	const RaggedArray2DT<int>* raggroup;
	fGroupData_2.Top();
	while ( fGroupData_2.Next(raggroup) )
	{
		int  nel = raggroup->MajorDim();
		for (int k = 0; k < nel; k++)
		{
			int  nen = raggroup->MinorDim(k);
			int* ien = (*raggroup)(k);
			for (int i = 0; i < nen; i++)
			{	
				int r_i = ien[i];
				if (r_i > -1)
				{
					for (int j = i+1; j < nen; j++)
					{
						int r_j = ien[j];
						if (r_j > -1 && r_i != r_j && // no refs to self
						    edgedata.AppendUnique(r_j - fShift, r_i))
							edgedata.Append(r_i - fShift, r_j);
					}
				}
			}
			ien += nen;
		}
	}
		
	/* copy/compress */
	fEdgeList.Copy(edgedata);

	/* find node with smallest degree */
	fMinDegree = fEdgeList.MinMinorDim(fMinDegreeNode, 0);
	fMinDegreeNode += fShift;
}

/* make the graph using the current data
* NOTE: assumes these are equation numbers, i.e.  1,... and skip
*       all < 1, but edge data needs to be 0,... */
void GraphT::MakeGraph(const iArrayT& active_rows, bool add_self)
{
	if (fMaxNodeNum == -1) return;
	
	/* take the range from active rows */
	SetRange(active_rows, true);

	/* generate active flags */
	int range = fMaxNodeNum - fMinNodeNum + 1;
	iArrayT active(range);
	active = 0;
	int* pactiverow = active_rows.Pointer();
	int num_active = active_rows.Length();
	for (int i = 0; i < num_active; i++)
		active[*pactiverow++ - fShift] = 1;		

	/* temp work space */
	AutoFill2DT<int> edgedata(range, 25, 5);
	
	/* initialize all lists with self */
	if (add_self)
	{
		int* pactiverow = active_rows.Pointer();
		for (int i = 0; i < num_active; i++)
		{
			int r_i = *pactiverow++;
			edgedata.Append(r_i - fShift, r_i - 1);
		}
	}

	/* create adjacency lists */	
	const iArray2DT* currgroup;
	fGroupData_1.Top();
	while (fGroupData_1.Next(currgroup))
	{
		int  nel = currgroup->MajorDim();
		int  nen = currgroup->MinorDim();
		int* ien = currgroup->Pointer(); //OFFSET 1,...
		for (int k = 0; k < nel; k++)
		{
			for (int i = 0; i < nen; i++)
			{	
				int r_i = ien[i];
				if (r_i > 0)
				{
					if (Active(r_i - fShift, active))
					{
						for (int j = i+1; j < nen; j++)
						{
							int r_j = ien[j];
							if (r_j > 0 &&
							    edgedata.AppendUnique(r_i - fShift, r_j - 1) &&
							    Active(r_j - fShift, active))
								edgedata.Append(r_j - fShift, r_i - 1);
						}
					}
					else
					{
						for (int j = i+1; j < nen; j++)
						{
							int r_j = ien[j];
							if (r_j > 0 && Active(r_j - fShift, active))
								edgedata.AppendUnique(r_j - fShift, r_i - 1);
						}
					}
				}	
			}
			ien += nen;
		}
	}

	const RaggedArray2DT<int>* raggroup;
	fGroupData_2.Top();
	while ( fGroupData_2.Next(raggroup) )
	{
		int  nel = raggroup->MajorDim();
		for (int k = 0; k < nel; k++)
		{
			int  nen = raggroup->MinorDim(k);
			int* ien = (*raggroup)(k); //OFFSET 1,...

			for (int i = 0; i < nen; i++)
			{	
				int r_i = ien[i];
				if (r_i > 0)
				{
					if (Active(r_i - fShift, active))
					{
						for (int j = i+1; j < nen; j++)
						{
							int r_j = ien[j];
							if (r_j > 0 &&
							    edgedata.AppendUnique(r_i - fShift, r_j - 1) &&
							    Active(r_j - fShift, active))
								edgedata.Append(r_j - fShift, r_i - 1);
						}
					}
					else
					{
						for (int j = i+1; j < nen; j++)
						{
							int r_j = ien[j];
							if (r_j > 0 && Active(r_j - fShift, active))
								edgedata.AppendUnique(r_j - fShift, r_i - 1);
						}
					}
				}
			}
			ien += nen;
		}
	}

	/* copy/compress */
	fEdgeList.Copy(edgedata);
	
	/* find node with smallest degree */
	fMinDegree = fEdgeList.MinMinorDim(fMinDegreeNode, 0);
	fMinDegreeNode += fShift;
}

void GraphT::MakeGraph(const iArrayT& active_rows, bool add_self, bool upper_only)
{
	if (fMaxNodeNum == -1) return;

	/* full graph */
	if (!upper_only)
		MakeGraph(active_rows, add_self);
	else
	{
		/* take the range from active rows */
		SetRange(active_rows, true);

		/* generate active flags */
		int range = fMaxNodeNum - fMinNodeNum + 1;
		iArrayT active(range);
		active = 0;
		int* pactiverow = active_rows.Pointer();
		int num_active = active_rows.Length();
		for (int i = 0; i < num_active; i++)
			active[*pactiverow++ - fShift] = 1;		

		/* temp work space */
		AutoFill2DT<int> edgedata(range, 25, 5);
	
		/* initialize all lists with self */
		if (add_self)
		{
			int* pactiverow = active_rows.Pointer();
			for (int i = 0; i < num_active; i++)
			{
				int r_i = *pactiverow++;
				edgedata.Append(r_i - fShift, r_i - 1);
			}
		}

		/* create adjacency lists */	
		const iArray2DT* currgroup;
		fGroupData_1.Top();
		while ( fGroupData_1.Next(currgroup) )
		{
			int  nel = currgroup->MajorDim();
			int  nen = currgroup->MinorDim();
			int* ien = currgroup->Pointer(); //OFFSET 1,...
			for (int k = 0; k < nel; k++)
			{
				for (int i = 0; i < nen; i++)
				{	
					int r_i = ien[i];
					if (r_i > 0)
						for (int j = i+1; j < nen; j++)
						{
							int r_j = ien[j];
							if (r_j > 0)
							{
//								if (r_j > r_i && Active(r_i - fShift, active))
//									edgedata.AppendUnique(r_i - fShift, r_j - 1);
//								else if (Active(r_j - fShift, active))
//									edgedata.AppendUnique(r_j - fShift, r_i - 1);
								if (r_j > r_i)
								{
									if (Active(r_i - fShift, active))
										edgedata.AppendUnique(r_i - fShift, r_j - 1);
								}
								else if (Active(r_j - fShift, active))
									edgedata.AppendUnique(r_j - fShift, r_i - 1);
							}
						}
				}
				ien += nen;
			}
		}

		const RaggedArray2DT<int>* raggroup;
		fGroupData_2.Top();
		while ( fGroupData_2.Next(raggroup) )
		{
			int  nel = raggroup->MajorDim();		
			for (int k = 0; k < nel; k++)
			{
				int  nen = raggroup->MinorDim(k);
				int* ien = (*raggroup)(k); //OFFSET 1,...
				for (int i = 0; i < nen; i++)
				{	
					int r_i = ien[i];
					if (r_i > 0)
						for (int j = i+1; j < nen; j++)
						{
							int r_j = ien[j];
							if (r_j > 0)
							{
								if (r_j > r_i)
								{
									if (Active(r_i - fShift, active))
										edgedata.AppendUnique(r_i - fShift, r_j - 1);
								}
								else if (Active(r_j - fShift, active))
									edgedata.AppendUnique(r_j - fShift, r_i - 1);
							}
						}
				}				
				ien += nen;
			}
		}

		/* copy/compress */
		fEdgeList.Copy(edgedata);
		
		/* this is not calculated */
		fMinDegree = -1;
		fMinDegreeNode = 1;
	}
}

/* return list of unconnected nodes */
void GraphT::UnconnectedNodes(iArrayT& nodes) const
{
	if (fMinDegree < 1)
	{
		/* collect */
		AutoArrayT<int> stray;
		stray.Allocate(0);
		for (int i = 0; i < fEdgeList.MajorDim(); i++)
			if (fEdgeList.MinorDim(i) == 0)
				stray.Append(i + fShift);
	
		/* set return value */
		nodes.Allocate(stray.Length());
		stray.CopyInto(nodes);
	}
	else
		nodes.Allocate(0);
}

/* label nodes by branch of graph */
void GraphT::LabelBranches(const iArrayT& nodes, iArrayT& branch_map)
{
	//TEMP
	if (fShift != 0)
	{
		cout << "\n GraphT::LabelBranches: not expecting non-zero\n"
		     <<   "     node number shift: " << fShift << endl;
		throw eGeneralFail;
	}

	/* rooted level structure */
	RootedLevelT rootedlevel;
	
	/* surface set data */
	iArrayT level_nodes;
	branch_map.Allocate(nodes.Max() + 1);
	branch_map = -1;
	int branch = 0;
	for (int i = 0; i < nodes.Length(); i++)
	{
		int nd = nodes[i];
		if (branch_map[nd] == -1)
		{
			/* mark set nodes */
			rootedlevel.MakePartialRootedLevel(*this, nd, true);
	
			/* mark root */
			branch_map[nd] = branch;
			
			/* mark level structure */
			for (int j = 1; j < rootedlevel.Depth(); j++)
			{
				rootedlevel.NodesOnLevel(level_nodes, j);
			
				int* pnodes = level_nodes.Pointer();
				int  dim = level_nodes.Length();
				for (int k = 0; k < dim; k++)
				{
					if (*pnodes >= branch_map.Length())
						cout << "hello" << endl;
				
					int& map = branch_map[*pnodes++];
#if __option(extended_errorcheck)
					if (map != -1) throw eGeneralFail;
#endif
					map = branch;
				}
			}
			
			/* next */
			branch++;	
		}
	}
}

void GraphT::Partition(const iArrayT& config, const iArrayT& weight,
	ArrayT<PartitionT>& partition, bool verbose)
{
	//TEMP
	if (fShift != 0)
	{
		cout << "\n GraphT::Partition: not expecting non-zero\n"
		     <<   "     node number shift: " << fShift << endl;
		throw eGeneralFail;
	}

	/* dimensions */
	int nnd = weight.Length();

	/* generate partition */
	iArrayT part_map(nnd);
	GraphBaseT::Partition(config, weight, part_map, verbose);

	/* time */
	clock_t t0 = clock();
	
	/* resolve internal/boundary nodes */
	partition.Allocate(config.Sum());
	for (int i = 0; i < partition.Length(); i++)
		partition[i].Set(partition.Length(), i, part_map, *this);
		
	/* set outgoing communication maps */
	for (int j = 0; j < partition.Length(); j++)
	{
		int ID = partition[j].ID();
		const iArrayT& commID = partition[j].CommID();
		int comm_size = commID.Length();
		
		/* collect nodes */
		ArrayT<iArrayT> nodes_out(comm_size);
		for (int i = 0; i < comm_size; i++)
		{
			const iArrayT* nodes_i = partition[commID[i]].NodesIn(ID);
			if (!nodes_i)
				throw eGeneralFail;
			else
				nodes_out[i].Alias(*nodes_i);
		}
		
		/* set */
		partition[j].SetOutgoing(nodes_out);
	}		

	/* time */
	clock_t t1 = clock();
	if (verbose)
		cout << setw(kDoubleWidth) << double(t1 - t0)/CLOCKS_PER_SEC
		     << " sec: GraphT::Partition: generate node maps" << endl;
}

/* using external graph to classify nodes */
void GraphT::Partition(const iArrayT& config, const iArrayT& weight,
	const GraphT& node_graph, ArrayT<PartitionT>& partition,
	bool verbose)
{
	//TEMP
	if (fShift != 0)
	{
		cout << "\n GraphT::Partition: not expecting non-zero\n"
		     <<   "     node number shift: " << fShift << endl;
		throw eGeneralFail;
	}

	/* dimensions */
	int nnd = weight.Length();

	/* generate partition */
	iArrayT part_map(nnd);
	GraphBaseT::Partition(config, weight, part_map, verbose);

	/* time */
	clock_t t0 = clock();
	
	/* resolve internal/boundary nodes */
	partition.Allocate(config.Sum());
	for (int i = 0; i < partition.Length(); i++)
		partition[i].Set(partition.Length(), i, part_map, node_graph);
		
	/* set outgoing communication maps */
	for (int j = 0; j < partition.Length(); j++)
	{
		int ID = partition[j].ID();
		const iArrayT& commID = partition[j].CommID();
		int comm_size = commID.Length();
		
		/* collect nodes */
		ArrayT<iArrayT> nodes_out(comm_size);
		for (int i = 0; i < comm_size; i++)
		{
			const iArrayT* nodes_i = partition[commID[i]].NodesIn(ID);
			if (!nodes_i)
				throw eGeneralFail;
			else
				nodes_out[i].Alias(*nodes_i);
		}
		
		/* set */
		partition[j].SetOutgoing(nodes_out);
	}		

	/* time */
	clock_t t1 = clock();
	if (verbose)
		cout << setw(kDoubleWidth) << double(t1 - t0)/CLOCKS_PER_SEC
		     << " sec: GraphT::Partition: generate node maps" << endl;
}
