/* $Id: GraphT.cpp,v 1.15 2003-11-21 22:41:54 paklein Exp $ */
#include "GraphT.h"

#include <time.h>

#include "iArrayT.h"
#include "iArray2DT.h"
#include "RaggedArray2DT.h"
#include "RowAutoFill2DT.h"
#include "AutoArrayT.h"
#include "RootedLevelT.h"
#include "PartitionT.h"

using namespace Tahoe;

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
	if (!fGroupData_1.AppendUnique(&groupdata)) throw ExceptionT::kGeneralFail;
}

void GraphT::AddGroup(const RaggedArray2DT<int>& groupdata)
{
	/* reset maximum node number*/
	iArrayT temp(groupdata.Length(), groupdata.Pointer());
	SetRange(temp);

	/* do not allow repeated registering of groups */
	if (!fGroupData_2.AppendUnique(&groupdata)) throw ExceptionT::kGeneralFail;
}

void GraphT::ClearGroups(void)
{
	fGroupData_1.Clear();
	fGroupData_2.Clear();
}

/* add a group to the graph */
void GraphT::AddEquivalentNodes(const iArray2DT& equivalentNodes)
{
	/* reset maximum node number */
	//SetRange(equivalentNodes);

	/* do not allow repeated registering of groups */
	if (!fEquivalentData.AppendUnique(&equivalentNodes)) throw ExceptionT::kGeneralFail;
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
		throw ExceptionT::kGeneralFail;
	}
	/* ignore ien < 0 */
	else if (fShift < 0)
		fShift = fMinNodeNum = 0;

	/* temp work space */
	int range = fMaxNodeNum - fMinNodeNum + 1;
	//AutoFill2DT<int> edgedata(range, 25, 5);
	RowAutoFill2DT<int> edgedata(range, 25, 5);
	//edgedata.SetFlushSize(8*100000000); //flush every 100MB

	/* create adjacency lists */	
	const iArray2DT* currgroup;
	fGroupData_1.Top();
	while ( fGroupData_1.Next(currgroup) )
	{
		int  nel = currgroup->MajorDim();
		int  nen = currgroup->MinorDim();
		const int* ien = currgroup->Pointer();
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
			int nen = raggroup->MinorDim(k);
			const int* ien = (*raggroup)(k);
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
	
#if 1
	const iArray2DT* currEquiv;
	fEquivalentData.Top();
	while ( fEquivalentData.Next(currEquiv) )
	{	
		int  nel = currEquiv->MajorDim();
		int  nen = currEquiv->MinorDim();
		const int* ien = currEquiv->Pointer();
		iArrayT row_i, row_j;
		for (int k = 0; k < nel; k++)
		{
			for (int i = 0; i < nen; i++)
			{	
				int r_i = ien[i];
				if (r_i > -1)
				{
					for (int j = i+1; j < nen; j++)
					{
						/* Copy all equivalent rows of edgedata to
						 * all the other equivalent rows 
						 */
						int r_j = ien[j];
						if (r_j > -1 && r_i != r_j) 	
					    {
					    	row_i.Dimension(edgedata.MinorDim(r_i - fShift));
					    	row_i.Copy(edgedata(r_i - fShift));
					    	row_j.Dimension(edgedata.MinorDim(r_j - fShift));
					    	row_j.Copy(edgedata(r_j - fShift));
							
							/* Loop over all entries to eliminate references to self */
							int *rz = row_i.Pointer();
							for (int z = 0; z < row_i.Length(); z++, rz++)
								if (*rz != r_j)
									edgedata.AppendUnique(r_j - fShift, *rz);
									
							rz = row_j.Pointer();
							for (int z = 0; z < row_j.Length(); z++, rz++)
								if (*rz != r_i)
									edgedata.AppendUnique(r_i - fShift, *rz);
						}
					}
				}
			}
			ien += nen;
		}
	
		/* Graph is no longer consistent, row_i and row_j
		 * are identical, but nodes k that were adjacent to j do not
		 * necessarily know about i.
		 */
		
		ien = currEquiv->Pointer();
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
						if (r_j > -1 && r_i != r_j) 	
					    {
					    	
					    	row_i.Dimension(edgedata.MinorDim(r_i - fShift));
					    	row_i.Copy(edgedata(r_i - fShift));
					    	row_j.Dimension(edgedata.MinorDim(r_j - fShift));
					    	row_j.Copy(edgedata(r_j - fShift));
									
					    	for (int z = 0; z < row_i.Length(); z++)
					    		if (r_j != row_i[z])
									edgedata.AppendUnique(row_i[z] - fShift,r_j);
							for (int z = 0; z < row_j.Length(); z++)
								if (r_i != row_j[z])
									edgedata.AppendUnique(row_j[z] - fShift,r_i);
						}
					}
				}
			}
			ien += nen;
		}
	}
	
#endif
		
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
	const int* pactiverow = active_rows.Pointer();
	int num_active = active_rows.Length();
	for (int i = 0; i < num_active; i++)
		active[*pactiverow++ - fShift] = 1;		

	/* temp work space */
	//AutoFill2DT<int> edgedata(range, 25, 5);
	RowAutoFill2DT<int> edgedata(range, 25, 5);
	
	/* initialize all lists with self */
	if (add_self)
	{
		const int* pactiverow = active_rows.Pointer();
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
		const int* ien = currgroup->Pointer(); //OFFSET 1,...
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
			int nen = raggroup->MinorDim(k);
			const int* ien = (*raggroup)(k); //OFFSET 1,...

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
		const int* pactiverow = active_rows.Pointer();
		int num_active = active_rows.Length();
		for (int i = 0; i < num_active; i++)
			active[*pactiverow++ - fShift] = 1;		

		/* temp work space */
		//AutoFill2DT<int> edgedata(range, 25, 5);
		RowAutoFill2DT<int> edgedata(range, 25, 5);
	
		/* initialize all lists with self */
		if (add_self)
		{
			const int* pactiverow = active_rows.Pointer();
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
			const int* ien = currgroup->Pointer(); //OFFSET 1,...
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
				int nen = raggroup->MinorDim(k);
				const int* ien = (*raggroup)(k); //OFFSET 1,...
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
		stray.Dimension(0);
		for (int i = 0; i < fEdgeList.MajorDim(); i++)
			if (fEdgeList.MinorDim(i) == 0)
				stray.Append(i + fShift);
	
		/* set return value */
		nodes.Dimension(stray.Length());
		stray.CopyInto(nodes);
	}
	else
		nodes.Dimension(0);
}

/* label nodes by branch of graph */
void GraphT::LabelBranches(const iArrayT& nodes, iArrayT& branch_map)
{
	//TEMP
	if (fShift != 0)
	{
		cout << "\n GraphT::LabelBranches: not expecting non-zero\n"
		     <<   "     node number shift: " << fShift << endl;
		throw ExceptionT::kGeneralFail;
	}

	/* rooted level structure */
	RootedLevelT rootedlevel;
	
	/* surface set data */
	iArrayT level_nodes;
	branch_map.Dimension(nodes.Max() + 1);
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
					if (map != -1) throw ExceptionT::kGeneralFail;
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
	ArrayT<PartitionT>& partition, bool verbose, int method)
{
	//TEMP
	if (fShift != 0)
	{
		cout << "\n GraphT::Partition: not expecting non-zero\n"
		     <<   "     node number shift: " << fShift << endl;
		throw ExceptionT::kGeneralFail;
	}

	/* dimensions */
	int nnd = weight.Length();

	/* generate partition */
	iArrayT part_map(nnd);
	if (method == 0)
		GraphBaseT::Partition(config, weight, part_map, verbose);
	else if (method == 1)
	{
		int num_partitions = config.Sum();
		int volume_or_edgecut = 1;
		Partition_METIS(num_partitions, weight, part_map, volume_or_edgecut);	
	}
	else throw;

	/* time */
	clock_t t0 = clock();
	
	/* resolve internal/boundary nodes */
	if (verbose) cout << " GraphT::Partition: classifying nodes" << endl;
	partition.Dimension(config.Sum());
	for (int i = 0; i < partition.Length(); i++)
		partition[i].Set(partition.Length(), i, part_map, *this);
		
	/* set outgoing communication maps */
	if (verbose) cout << " GraphT::Partition: setting communication maps" << endl;
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
				throw ExceptionT::kGeneralFail;
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
	bool verbose, int method)
{
	//TEMP
	if (fShift != 0)
	{
		cout << "\n GraphT::Partition: not expecting non-zero\n"
		     <<   "     node number shift: " << fShift << endl;
		throw ExceptionT::kGeneralFail;
	}

	/* dimensions */
	int nnd = weight.Length();

	/* generate partition */
	iArrayT part_map(nnd);
	if (method == 0)
		GraphBaseT::Partition(config, weight, part_map, verbose);
	else if (method == 1)
	{
		int num_partitions = config.Sum();
		int volume_or_edgecut = 1;
		Partition_METIS(num_partitions, weight, part_map, volume_or_edgecut);	
	}
	else throw;
	
	/* time */
	clock_t t0 = clock();
	
	/* resolve internal/boundary nodes */
	if (verbose) cout << " GraphT::Partition: classifying nodes" << endl;
	partition.Dimension(config.Sum());
	for (int i = 0; i < partition.Length(); i++)
		partition[i].Set(partition.Length(), i, part_map, node_graph);
		
	/* set outgoing communication maps */
	if (verbose) cout << " GraphT::Partition: setting communication maps" << endl;
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
				throw ExceptionT::kGeneralFail;
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

void GraphT::Partition(const iArrayT& config, const iArrayT& weight,
	const ArrayT<const iArray2DT*>& connects_1, const ArrayT<const RaggedArray2DT<int>*>& connects_2, 
	ArrayT<PartitionT>& partition, bool verbose, int method)
{
	//TEMP
	if (fShift != 0)
	{
		cout << "\n GraphT::Partition: not expecting non-zero\n"
		     <<   "     node number shift: " << fShift << endl;
		throw ExceptionT::kGeneralFail;
	}

	/* dimensions */
	int nnd = weight.Length();

	/* generate partition */
	iArrayT part_map(nnd);
	if (method == 0)
		GraphBaseT::Partition(config, weight, part_map, verbose);
	else if (method == 1)
	{
		int num_partitions = config.Sum();
		int volume_or_edgecut = 1;
		Partition_METIS(num_partitions, weight, part_map, volume_or_edgecut);	
	}
	else throw;

	/* time */
	clock_t t0 = clock();
	
	/* resolve internal/boundary nodes */
	if (verbose) cout << " GraphT::Partition: classifying nodes" << endl;
	partition.Dimension(config.Sum());
	for (int i = 0; i < partition.Length(); i++)
		partition[i].Set(partition.Length(), i, part_map, connects_1, connects_2);
		
	/* set outgoing communication maps */
	if (verbose) cout << " GraphT::Partition: setting communication maps" << endl;
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
				throw ExceptionT::kGeneralFail;
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
