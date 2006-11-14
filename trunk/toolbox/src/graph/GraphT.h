/* $Id: GraphT.h,v 1.7 2006-11-14 03:27:38 paklein Exp $ */
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

namespace Tahoe {

/* forward declarations */
class iArray2DT;
class PartitionT;

class GraphT: public GraphBaseT
{
public:

	/* constructor */
	GraphT(bool verbose = false);

	/** \name structure data */
	/*@{*/
	/** add a regular group to the graph */
	void AddGroup(const iArray2DT& groupdata);

	/** add a irregular group to the graph */
	void AddGroup(const RaggedArray2DT<int>& groupdata);

	/** reset structure information */
	void ClearGroups(void);

	void GetGroups(ArrayT<const iArray2DT*>& group_data);
	void GetGroups(ArrayT<const RaggedArray2DT<int>*>& group_data);
	/*@}*/
	
	/* add nodes whose connectivities should be copied to each other */
	void AddEquivalentNodes(const iArray2DT& equivalentNodes);
	
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
		ArrayT<PartitionT>& partition, bool verbose, int method);

	/* forward base class partitioning function */
	void Partition(const iArrayT& config, const iArrayT& weight,
		iArrayT& partition, bool verbose);

	/* using external graph to classify nodes */
	void Partition(const iArrayT& config, const iArrayT& weight,
		const GraphT& node_graph, ArrayT<PartitionT>& partition,
		bool verbose, int method);

	/* using raw connectivities to classify nodes */
	void Partition(const iArrayT& config, const iArrayT& weight,
		const ArrayT<const iArray2DT*>& connects_1,
		const ArrayT<const RaggedArray2DT<int>*>& connects_2,
		ArrayT<PartitionT>& partition,
		bool verbose, int method);

private:

	/* return true if node is in range and active */
	bool Active(int node, const iArrayT& active) const;	
	
private:

	/* connectivity groups */
	LinkedListT<const iArray2DT*>           fGroupData_1;
	LinkedListT<const RaggedArray2DT<int>*> fGroupData_2;
	LinkedListT<const iArray2DT*>			fEquivalentData;
};

/* inlines */
inline void GraphT::Partition(const iArrayT& config, const iArrayT& weight,
	iArrayT& partition, bool verbose)
{
	/* inherited */
	GraphBaseT::Partition(config, weight, partition, verbose);
}

} // namespace Tahoe 
#endif /* _GRAPH_T_H_ */
