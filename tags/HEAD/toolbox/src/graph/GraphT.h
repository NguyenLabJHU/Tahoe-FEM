/* $Id: GraphT.h,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
/* created: paklein (08/05/1996)                                          */
/* generates graphs for the connectivities registered with AddGroup().    */
/* connectivies can have an arbitrary MinorDim(), but the labels in       */
/* the range of 0...[max label - 1] must all be used at least once.       */

#ifndef _GRAPH_T_H_
#define _GRAPH_T_H_

/* base class */
#include "GraphBaseT.h"

/* direct members */
#include "LinkedListT.h"

/* forward declarations */
class iArray2DT;
class PartitionT;

class GraphT: public GraphBaseT
{
public:

	/* constructor */
	GraphT(bool verbose = false);

	/* add a group to the graph */
	void AddGroup(const iArray2DT& groupdata);
	void AddGroup(const RaggedArray2DT<int>& groupdata);
	void ClearGroups(void);
	
	/* make the graph using the current data */
	void MakeGraph(void);
	
	/* these not compatible with partitioning functions */
	void MakeGraph(const iArrayT& active_rows, bool add_self);
	void MakeGraph(const iArrayT& active_rows, bool add_self, bool upper_only);
	
	/* return list of unconnected nodes */
	void UnconnectedNodes(iArrayT& nodes) const;
	
	/* label nodes by branch of graph */
	void LabelBranches(const iArrayT& nodes, iArrayT& branch_map);

	/* partition -
	 *      config: i x j x k x ...  rectangular partition dimensions
	 *      weight: nodal weights used for load balancing
	 *   partition: decomposition data per partition */
	void Partition(const iArrayT& config, const iArrayT& weight,
		ArrayT<PartitionT>& partition, bool verbose);

	void Partition(const iArrayT& config, const iArrayT& weight,
		iArrayT& partition, bool verbose);

	/* using external graph to classify nodes */
	void Partition(const iArrayT& config, const iArrayT& weight,
		const GraphT& node_graph, ArrayT<PartitionT>& partition,
		bool verbose);

private:

	/* return true if node is in range and active */
	bool Active(int node, const iArrayT& active) const;	
	
private:

	/* connectivity groups */
	LinkedListT<const iArray2DT*>           fGroupData_1;
	LinkedListT<const RaggedArray2DT<int>*> fGroupData_2;
};

/* inlines */
inline void GraphT::Partition(const iArrayT& config, const iArrayT& weight,
	iArrayT& partition, bool verbose)
{
	/* inherited */
	GraphBaseT::Partition(config, weight, partition, verbose);
}

#endif /* _GRAPH_T_H_ */
