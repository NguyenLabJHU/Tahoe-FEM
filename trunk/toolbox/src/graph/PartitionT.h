/* $Id: PartitionT.h,v 1.12 2004-11-17 23:21:29 paklein Exp $ */
/* created: paklein (11/16/1999) */
#ifndef _PARTITION_T_H_
#define _PARTITION_T_H_

/* direct members */
#include "iArrayT.h"
#include "StringT.h"

namespace Tahoe {

/* forward declarations */
class GraphT;
class iArray2DT;
class dArray2DT;
class ifstreamT;
template <class TYPE> class RaggedArray2DT;

/** graph partition information (following NEMESIS data model)
 * class generates complete decomposition information using
 * a node-to-partition map and a graph. Initialization is performed
 * in the following stages:\n
 * -# labels nodes using PartitionT::Set, which also sets the incoming node lists
 * -# set the outgoing node lists using PartitionT::SetOutgoing
 * -# set the list of element blocks using PartitionT::InitElementBlocks
 * -# set element blocks PartitionT::SetElements
 * -# convert the node and element numbering using PartitionT::SetScope 
 * -# optionally, set the decomposition type using PartitionT::SetDecomp 
 * Partition information can also be set by reading the information from
 * a file. */
class PartitionT
{
public:

	/** scope of node and element numbers */
	enum NumberScopeT {kUnSet = 0, /**< not defined */
	                   kLocal = 1, /**< numbering within the partition */
	                  kGlobal = 2  /**< global numbering */ };

	/** partition type */
	enum DecompTypeT {
	kUndefined =-1, /**< undefined partition type */
        kGraph = 0, /**< partition based on connectivities */
        kIndex = 1, /**< partition based on node index */
      kSpatial = 2  /**< partition based on position */
	};

	/** constructor */
	PartitionT(void);

	/** returns true if version if current */
	static bool CheckVersion(const StringT& version);

	/** \name accessors */
	/*@{*/
	/** comment character used for partition files */
	static char CommentMarker(void) { return '#'; };

	/** ID of this partition */
	int ID(void) const;

	/** current scope for the node and element numbers. The numbering
	 * scope is changed using PartitionT::SetScope. */
	NumberScopeT NumberScope(void) const;
	
	/** return the decomposition type */
	DecompTypeT DecompType(void) const { return fDecompType; };

	/** number of element blocks */
	int NumElementBlocks(void) const;

	/** id's of the element blocks */
	const ArrayT<StringT>& BlockID(void) const { return fElementBlockID; };
	/*@}*/

	/** \name nodes by classification 
	 * Node numbers are given by the current PartitionT::NumberScope */
	/*@{*/
	/** nodes owned by this partition which interact only with other nodes 
	 * owned by this partition */
	const iArrayT& Nodes_Internal(void) const;

	/** nodes owned by this partition that interact with at least one node
	 * which is not owned by this partition */
	const iArrayT& Nodes_Border(void) const;

	/** nodes not owned by this partition that interact with at least one
	 * node that is owned by this partition */
	const iArrayT& Nodes_External(void) const;
	
	/** return the number of nodes owned by this partition */
	int NumPartitionNodes(void) const;
	/*@}*/

	/** collect the node ids owned by this partition. There are the internal
	 * and border nodes.
	 * \param nodes returns with the nodes owned by this partition
	 * \param scope numbering scope for the node numbers
	 */
	void PartitionNodes(iArrayT& nodes, NumberScopeT scope) const;

	/** \name communicated nodes 
	 * The ids in the current PartitionT::NumberScope of the nodes with
	 * incoming and outgoing information. Requesting nodes to and from
	 * nodes not in the communication list returns NULL. */
	/*@{*/
	/** list of partitions communicating with this one */
	const iArrayT& CommID(void) const;

	/** list of nodes with information coming from the given partition */
	const iArrayT* NodesIn(int commID) const;

	/** list of nodes with information going to the given partition */
	const iArrayT* NodesOut(int commID) const;
	/*@}*/

	/** \name set partition node data
	 * Labels nodes as internal, border, and external and set communication
	 * lists for incoming and outgoing nodes.
     * \param num_parts total number of partitions
	 * \param id id of this partition
	 * \param part_map partition of all nodes in global numbering
	 */
	/*@{*/
	/** labels nodes as internal, border, and external using a graph of
	 * of the global model. Method also determines which partitions
	 * provide the information for all external nodes and sets the list
	 * of partitions returned by PartitionT::CommID. */
	void Set(int num_parts, int id, const iArrayT& part_map, const GraphT& graph);

	/** labels nodes as internal, border, and external using the connectivities
	 * of the global model. Method also determines which partitions
	 * provide the information for all external nodes and sets the list
	 * of partitions returned by PartitionT::CommID. */
	void Set(int num_parts, int id, const iArrayT& part_map, const ArrayT<const iArray2DT*>& connects_1,
		const ArrayT<const RaggedArray2DT<int>*>& connects_2);

	/** set the lists of nodes communicated to other partitions. Each partition could
	 * determine this based on the information provided by PartitionT::Set, but it is
	 * not. This information can be collected from other partitions using PartitionT::NodesIn.
	 * The order of the outgoing node lists must be the same as the order of
	 * PartitionT::CommID */
	void SetOutgoing(const ArrayT<iArrayT>& nodes_in);
	/*@}*/

	/** set partition element data
	 * These methods should be called after the methods used to set
	 * the partition node data */
	/*@{*/
	/** set the list of element block id's for this partition */
	void InitElementBlocks(const ArrayT<StringT>& blockID);	

	/** set the connectivities for the given block id. The is must
	 * match one of the id's set by PartitionT::InitElementBlocks. 
	 * Elements are */
	void SetElements(const StringT& blockID, const iArray2DT& connects);
	/*@}*/

	/** set the decomposition type */
	void SetDecompType(DecompTypeT t) { fDecompType = t; };

	/** change the node and element numbering scope. Nothing happens if the
	 * current scope is already the requested scope. */
	void SetScope(NumberScopeT scope);

	/** check cross-references - returns 1 if OK */
	int CrossCheck(const PartitionT& that) const;

	/** \name parameters for spatial decomposition only
	 * Methods throw exception if decomposition type is not PartitionT::kSpatial */
	/*@{*/
	void SetGridDimensions(const iArrayT& grid_dims) { fGridDims = grid_dims; };
	void SetGridPosition(const iArrayT& grid_position) { fGridPosition = grid_position; };
	
	const iArrayT& GridDimensions(void) const { return fGridDims; };
	const iArrayT& GridPosition(void) const { return fGridPosition; };
	/*@}*/

	/** \name I/O */
	/*@{*/
	/** stream extraction */
	friend ifstreamT& operator>>(ifstreamT& in, PartitionT& partition) { return partition.Read(in); };

	/** stream insertion */
	friend ostream& operator<<(ostream& out, const PartitionT& partition);

	/** stream extraction for the enum PartitionT::NumberScopeT */
	friend istream& operator>>(istream& in, PartitionT::NumberScopeT& scope);

	/** stream extraction for the enum PartitionT::DecompTypeT */
	friend istream& operator>>(istream& in, PartitionT::DecompTypeT& t);

	/** non-operator version of the read function */
	ifstreamT& Read(ifstreamT& in);
	/*@}*/

	/** \name maps */
	/*@{*/
	const iArrayT& NodeMap(void) const;
	const iArrayT& InverseNodeMap(int& index_shift) const;
	const iArrayT& ElementMap(const StringT& blockID) const;
	/*@}*/

	/** returns indeces of global nodes that lie within the partition. Includes
	 * the global nodes that appear in the internal, border, and external nodes,
	 * that is, \e all nodes that appear in this partition not just those nodes
	 * \e owned by this partition. */
	void ReturnPartitionNodes(const iArrayT& global_nodes,
		iArrayT& partition_indices) const;

	/** returns indeces of (block) global elements that lie within
	 * the partition */
	void ReturnPartitionElements(const StringT& blockID, const iArrayT& global_elements,
		iArrayT& partition_indices) const;
		
	/** return the node to processor map for nodes listed in PartitionT::NodeMap */
	void ReturnProcessorMap(ArrayT<int>& n2p) const;

	/** \name mapping functions
	 * The mapping functions assumes the scope of the node and elements
	 * being passed on is currently the opposite of the requested scope */
	/*@{*/
	void SetNodeScope(NumberScopeT scope, ArrayT<int>& nodes) const;
	void SetElementScope(NumberScopeT scope, const StringT& blockID, ArrayT<int>& elements) const;
	/*@}*/

private:

	/** node and element classifications */
	enum StatusT {   kUnset =-1,
	              kInternal = 0,
	                kBorder = 1,
	              kExternal = 2};
	//there are actually 4 types -> internal,
	//                              border-internal
	//                              border-external
	//                              external

	/** resolve element block ID to index */
	int ElementBlockIndex(const StringT& blockID, const char* caller = NULL) const;

	/** number transformations */
	void MapValues(const iArrayT& map, int shift, ArrayT<int>& values) const;

	/** make inverse map (filled with -1) */
	void MakeInverseMap(const iArrayT& map, iArrayT& inv_map, int& shift) const;

	/** classify set nodes as label nodes as internal, external, or border */
	/*@{*/
	void ClassifyNodes(const iArrayT& part_map, const GraphT& graph);
	void ClassifyNodes(const iArrayT& part_map, const ArrayT<const iArray2DT*>& connects_1,
		const ArrayT<const RaggedArray2DT<int>*>& connects_2);
	/*@}*/

	/* set receiving nodes/partition information */
	void SetReceive(const iArrayT& part_map);

	/** map status of (in range) parts into status_map */
	void MapStatus(StatusT status, const iArrayT& part, ArrayT<StatusT>& status_map,
		int offset);

	/** \name set numbering maps */
	/*@{*/
	void SetNodeMap(NumberScopeT scope, iArrayT& map, int& shift) const;
	void SetElementMap(NumberScopeT scope, const StringT& blockID, iArrayT& map, int& shift) const;
	/*@}*/
	
private:

	/** \name basic parameters */
	/*@{*/	
	int fNumPartitions;  /**< total number of partitions */
	int fID;             /**< partition number           */
	NumberScopeT fScope; /**< local or global numbering  */
	DecompTypeT  fDecompType; /**< decomposition type    */
	/*@}*/	
	
	/** \name nodal information */
	/*@{*/
	iArrayT fNodes_i; /**< internal nodes */
	iArrayT fNodes_b; /**< border nodes	  */
	iArrayT fNodes_e; /**< external nodes */
	/*@}*/
	
	/** \name receive/send information */
	/*@{*/
	iArrayT fCommID; /**< ID's of communicating partitions (only) */
	ArrayT<iArrayT> fNodes_in;  /**< nodes per comm part */
	ArrayT<iArrayT> fNodes_out; /**< nodes per comm part */
	/*@}*/
	
	/** \name element information */
	/*@{*/
	ArrayT<StringT> fElementBlockID; /**< element block ID's */
	ArrayT<iArrayT> fElements_i; /**< internal elements per block */
	ArrayT<iArrayT> fElements_b; /**< border elements per block   */
	/*@}*/

	/** \name node numbering maps */
	/*@{*/
	/** map from local to global node numbers
	 * global[local] for all _i, _b, _e nodes */
	iArrayT fNodeMap; 
	int fNodeMapShift;
	iArrayT fInvNodeMap;
	/*@}*/
	
	/** \name element numbering maps by block */
	/*@{*/
	/** block global element numbering (number within global blocks) */
	ArrayT<iArrayT> fElementMap; // block_global[block_local]
	iArrayT fElementMapShift;
	ArrayT<iArrayT> fInvElementMap;
	/*@}*/

	/** \name additional information for spatial decomposition */
	/*@{*/
	/** number of grid cells in each coordinate dimension */
	iArrayT fGridDims;
	
	/** location of this partition in the grid */
	iArrayT fGridPosition;
	/*@}*/
};

/* inlines */
inline int PartitionT::ID(void) const { return fID; }
inline PartitionT::NumberScopeT PartitionT::NumberScope(void) const { return fScope; }

inline const iArrayT& PartitionT::Nodes_Internal(void) const { return fNodes_i; }
inline const iArrayT& PartitionT::Nodes_Border(void) const { return fNodes_b; }
inline const iArrayT& PartitionT::Nodes_External(void) const { return fNodes_e; }

inline int PartitionT::NumPartitionNodes(void) const
{
	return fNodes_i.Length() + fNodes_b.Length();
}

inline 	const iArrayT& PartitionT::CommID(void) const { return fCommID; }

inline int PartitionT::NumElementBlocks(void) const { return fElementMap.Length(); }

/* maps */
inline const iArrayT& PartitionT::NodeMap(void) const { return fNodeMap; }
inline const iArrayT& PartitionT::ElementMap(const StringT& blockID) const
{
	return fElementMap[ElementBlockIndex(blockID, "ElementMap")];
}

} // namespace Tahoe 
#endif /* _PARTITION_T_H_ */
