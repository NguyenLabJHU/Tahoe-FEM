/* $Id: GraphBaseT.h,v 1.1.1.1 2001-01-25 20:56:27 paklein Exp $ */
/* created: paklein (04/13/1999)                                          */
/* base class for graph class. does not handle making the graph.          */

#ifndef _GRAPHBASE_T_H_
#define _GRAPHBASE_T_H_

/* direct members */
#include "RaggedArray2DT.h"

/* forward declarations */
class iArrayT;
class iArray2DT;

class GraphBaseT
{
public:

	/* constructor */
	GraphBaseT(bool verbose);

	/* destructor */
	virtual ~GraphBaseT(void);
		
	/* accessors */
	int  Degree(int nodenum) const; // number of edges incident on the node
	void Degrees(ArrayT<int>& degrees, int& row_shift) const; // all of them
	int  NumNodes(void) const;
	int* Edges(int nodenum) const;  // adjacency list for the node
	int  MinDegree(void) const;
	int  MinDegreeNode(void) const;
	void GetEdges(int nodenum, ArrayT<int>& edges) const; // shallow copy

	/* set adjacency list */
	void SetAdjacency(const RaggedArray2DT<int>& edge_list);

	/* partition -
	 *      config: i x j x k x ...  rectangular partition dimensions
	 *      weight: nodal weights used for load balancing
	 *   partition: partition[node] = partition of the node */
	void Partition(const iArrayT& config, const iArrayT& weight,
		iArrayT& partition, bool verbose);

	/* fill in the degrees for the specified nodes */
	void ReturnDegrees(const ArrayT<int>& nodes, ArrayT<int>& degrees) const;

	/* verify graph data is consistent - returns 1 if consistent, 0 otherwise */
	int Verify(ostream& err) const;
	
	/* output graph to stream */
	void Write(ostream& out) const;

	/* initialize the priorities as specified in Sloan.  Assumes DIRECT
	 * INDEXING in priorities */
	void InitializePriorities(iArrayT& priorities, int W1) const;

protected:

	/* set the node range and shift parameters */
	void SetRange(const nArrayT<int>& node_list, bool reset = false);

	/* graph data - row numbers are shifted down */
	const RaggedArray2DT<int>& EdgeList(int& row_shift) const;

private:

	/* contract the graph */
	void Contract(const GraphBaseT& parent, iArrayT& map);

	/* return the "best" edge to collapse */
	int SelectCollapse(const ArrayT<int>& edges, const ArrayT<int>& degrees,
		const ArrayT<int>& flags) const;

	/* set gains - returns number of cut edges */
	int SetGains(const GraphBaseT& graph, const iArrayT& partition,
		iArray2DT& gain) const;
	
	/* map of destination partitions for each node */
	void SetMoves(const iArray2DT& gain, const iArrayT& partition, const iArrayT& size,
		iArrayT& move) const;
	
	/* apply moves/update data - returns change in cuts */
	int ApplyMoves(const GraphBaseT& graph, const iArrayT& move, const iArrayT& weight,
		iArrayT& partition, iArray2DT& gain, iArrayT& size);
	
	/* optimize partitions - returns number of exchanges */
	void OptimizeParts(const GraphBaseT& graph,  const iArrayT& weights,
		iArrayT& partition, int dim, int& repetitions, int& cuts);

	/* sort into ascending order (from the Aztec library) */
	void AZ_sort(int list[], int N, int list2[], double list3[]);
	
protected:

	//TEMP
	bool fVerbose;

	/* runtime parameters */
	int fMaxNodeNum;
	int fMinNodeNum;
	int fMinDegree;
	int fMinDegreeNode;
	
	/* adjacency lists */
	int fShift;
	RaggedArray2DT<int> fEdgeList;
};

/* in-lines */

/* accessors */
inline int GraphBaseT::Degree(int nodenum) const
{
	return fEdgeList.MinorDim(nodenum - fShift);
}

inline void GraphBaseT::Degrees(ArrayT<int>& degrees, int& row_shift) const
{
	row_shift = fShift;
	fEdgeList.MinorDim(degrees);
}

inline int GraphBaseT::NumNodes(void) const { return fEdgeList.MajorDim(); }
inline int* GraphBaseT::Edges(int nodenum) const
{
	return fEdgeList(nodenum - fShift);
}

inline int GraphBaseT::MinDegree(void) const { return fMinDegree; }
inline int GraphBaseT::MinDegreeNode(void) const { return fMinDegreeNode; }

inline void GraphBaseT::GetEdges(int nodenum, ArrayT<int>& edges) const // shallow copy
{
	fEdgeList.RowAlias(nodenum - fShift, edges);
}

/* graph data */
inline const RaggedArray2DT<int>& GraphBaseT::EdgeList(int& row_shift) const
{
	row_shift = fShift;
	return fEdgeList;
}

#endif /* _GRAPHBASE_T_H_ */
