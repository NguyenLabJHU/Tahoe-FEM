/* $Id: PartitionT.h,v 1.1.1.1 2001-01-25 20:56:27 paklein Exp $ */
/* created: paklein (11/16/1999)                                          */
/* graph partition information (following NEMESIS data model)             */
/* class generates complete decomposition information using               */
/* a node-to-partition map and a graph. initialization is in              */
/* 3 stages:                                                              */
/* (1) Set(,,) to set internal nodal data                                 */
/* (2) [cross-link] to complete communication maps                        */
/* (3) [tolocal] to convert to partition-local numbering                  */

#ifndef _PARTITION_T_H_
#define _PARTITION_T_H_

/* direct members */
#include "iArrayT.h"

/* forward declarations */
class GraphT;
class iArray2DT;
class ifstreamT;
class StringT;

class PartitionT
{
public:

	enum NumberScopeT {kUnSet = 0,
	                   kLocal = 1,
	                  kGlobal = 2};

	/* constructor */
	PartitionT(void);

	/* returns true if version if current */
	static bool CheckVersion(const StringT& version);
	
	/* accessors */
	int ID(void) const;
	NumberScopeT NumberScope(void) const;

	const iArrayT& Nodes_Internal(void) const;
	const iArrayT& Nodes_Border(void) const;
	const iArrayT& Nodes_External(void) const;
	void PartitionNodes(iArrayT& nodes, NumberScopeT scope) const;

	const iArrayT& CommID(void) const;
	const iArrayT* NodesIn(int commID) const;  // NULL if ID not found
	const iArrayT* NodesOut(int commID) const; // NULL if ID not found

	/* set partition data */
	void Set(int num_parts, int id, const iArrayT& part_map, const GraphT& graph);
	void SetOutgoing(const ArrayT<iArrayT>& nodes_in); // in sequence of CommID
	void SetScope(NumberScopeT scope);

	int NumElementBlocks(void) const;
	void InitElementBlocks(const iArrayT& blockID);	
	void SetElements(int blockID, const iArray2DT& connects);

	/* check cross-references - returns 1 if OK */
	int CrossCheck(const PartitionT& that) const;

	/* I/O */
	friend ifstreamT& operator>>(ifstreamT& in, PartitionT& partition) { return partition.Read(in); };
	friend ostream& operator<<(ostream& out, const PartitionT& partition);

	/* operator support */
	ifstreamT& Read(ifstreamT& in);

	/* maps */
	const iArrayT& NodeMap(void) const;
	const iArrayT& InverseNodeMap(int& index_shift) const;
	const iArrayT& ElementMap(int blockID) const;

	/* returns indeces of global nodes that lie within the partition */
	void ReturnPartitionNodes(const iArrayT& global_nodes,
		iArrayT& partition_indices) const;

	/* returns indeces of (block) global elements that lie within
	 * the partition */
	void ReturnPartitionElements(int blockID, const iArrayT& global_elements,
		iArrayT& partition_indices) const;

	/* mapping functions (assumes scope is currently the opposite) */
	void SetNodeScope(NumberScopeT scope, ArrayT<int>& nodes) const;
	void SetElementScope(NumberScopeT scope, int blockID, ArrayT<int>& elements) const;

	/* input operator for scope */
	friend istream& operator>>(istream& in, PartitionT::NumberScopeT& scope);

private:

	/* node and element classifications */
	enum StatusT {kInternal = 0,
	                kBorder = 1,
	              kExternal = 2};
	//there are actually 4 types -> internal,
	//                              border-internal
	//                              border-external
	//                              external

	/* resolve element block ID to index */
	int ElementBlockIndex(int blockID, const char* caller = NULL) const;

	/* number transformations */
	void MapValues(const iArrayT& map, int shift, ArrayT<int>& values) const;

	/* make inverse map (filled with -1) */
	void MakeInverseMap(const iArrayT& map, iArrayT& inv_map, int& shift) const;

	/* classify set nodes as label nodes as internal, external, or border */
	void ClassifyNodes(const iArrayT& part_map, const GraphT& graph);

	/* set receiving nodes/partition information */
	void SetReceive(const iArrayT& part_map);

	/* map status of (in range) parts into status_map */
	void MapStatus(StatusT status, const iArrayT& part, ArrayT<StatusT>& status_map,
		int offset);

	/* set numbering maps */
	void SetNodeMap(NumberScopeT scope, iArrayT& map, int& shift) const;
	void SetElementMap(NumberScopeT scope, int blockID, iArrayT& map, int& shift) const;

private:
	
	int fNumPartitions; // total number of partitions
	int fID;            // partition number
	NumberScopeT fScope; // local or global numbering
	
	// nodal information
	iArrayT fNodes_i; // internal nodes	
	iArrayT fNodes_b; // border nodes	
	iArrayT fNodes_e; // external nodes
	
	// receive/send information
	iArrayT fCommID; // ID's of communicating partitions (only)
	ArrayT<iArrayT> fNodes_in;  // nodes per comm part
	ArrayT<iArrayT> fNodes_out; // nodes per comm part
	
	// element information
	iArrayT fElementBlockID;
	ArrayT<iArrayT> fElements_i; // internal elements per block
	ArrayT<iArrayT> fElements_b; // border elements per block

	// global node map
	iArrayT fNodeMap; // global[local] for all _i, _b, _e nodes
	int fNodeMapShift;
	iArrayT fInvNodeMap;
	
	// block global element numbering (number within global blocks)
	ArrayT<iArrayT> fElementMap; // block_global[block_local]
	iArrayT fElementMapShift;
	ArrayT<iArrayT> fInvElementMap;	
};

/* inlines */
inline int PartitionT::ID(void) const { return fID; }
inline PartitionT::NumberScopeT PartitionT::NumberScope(void) const { return fScope; }

inline const iArrayT& PartitionT::Nodes_Internal(void) const { return fNodes_i; }
inline const iArrayT& PartitionT::Nodes_Border(void) const { return fNodes_b; }
inline const iArrayT& PartitionT::Nodes_External(void) const { return fNodes_e; }

inline 	const iArrayT& PartitionT::CommID(void) const { return fCommID; }

inline int PartitionT::NumElementBlocks(void) const { return fElementMap.Length(); }

/* maps */
inline const iArrayT& PartitionT::NodeMap(void) const { return fNodeMap; }
inline const iArrayT& PartitionT::ElementMap(int blockID) const
{
	return fElementMap[ElementBlockIndex(blockID, "ElementMap")];
}

#endif /* _PARTITION_T_H_ */
