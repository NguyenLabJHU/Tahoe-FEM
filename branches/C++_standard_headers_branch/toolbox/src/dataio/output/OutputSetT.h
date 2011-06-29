/* $Id: OutputSetT.h,v 1.21 2005-06-06 06:38:24 paklein Exp $ */
/* created: paklein (03/07/2000) */

#ifndef _OUTPUTSET_T_H_
#define _OUTPUTSET_T_H_

/* direct members */
#include "GeometryT.h"
#include "StringT.h"
#include "iArrayT.h"
#include "iArray2DT.h"
#include "VariArrayT.h"

namespace Tahoe {

/* forward declarations */
class iArray2DT;

/** class to act as specifier for output data. A class that would
 * like to write output, needs to use the following procedure:
 * -# construct an OutputSetT with its output data specifications.
 * -# register this output data specification with an output formatter.
 *     Output formatters are classes derived from OutputBaseT. Output
 *     is registered using OutputBaseT::AddElementSet. Note that this
 *     function is accessed in Tahoe through FEManagerT::RegisterOutput
 *     which passes the call through IOManager::AddElementSet to
 *     OutputBaseT::AddElementSet. Registering output returns an integer
 *     output ID (OID) which must be used in step (3). Note that the OID 
 *     obtained during registration is generally not the same as the ID 
 *     passed into the OutputSetT constructor. A class with a single ID 
 *     may register an arbitary number of output sets, receiving a unique
 *     OID for each set.
 * -# write output by sending output data to the output formatter using
 *     OutputBaseT::WriteOutput with the ID obtained during registration
 *     in step (2). Note that this function is accessed in Tahoe through 
 *     FEManagerT::WriteOutput which passes the call through 
 *     IOManager::AddElementSet to the output formatter.
 *
 * Outside of the constructor, the remainder of the interface is intended
 * for use by the output formatter, not the class sending the data for
 * output. 
 *
 * There are three modes for output sets, denoted by the OutputSetT::ModeT
 * enum:
 * - OutputSetT::kElementBlock - data written over element blocks defined
 *     in the global geometry database. Both nodal and element output is
 *     supported in this mode since local element numbers can be mapped to
 *     global element numbers.
 * - OutputSetT::kFreeSet - data written over an arbitrary set of nodes.
 *     Only nodal output is supported in this mode since no map from local
 *     to global elements is available in the PartitionT files.
 * - OutputSetT::kElementFromSideSet - data written over element blocks not
 *     defined in the global geometry database but derived from a sideset
 *     of a specific element block that IS in the global geometry database
 */
class OutputSetT
{
public:

	/** set mode */
	enum ModeT {kElementBlock = 0, /**< writing data over element blocks defined in the
	                                * global geometry file */
	            kFreeSet = 1, /**< writing data over arbitrary sets of nodes */
	            kElementFromSideSet = 2 /**< writing data over element blocks created from side sets */
	            };

	/** generate output data record.
	 * \param geometry_code GeometryT::CodeT defining the geometry associated
	 *        with the connectivities.
	 * \param block_ID list of element block ID's comprising the connectivities.
	 *        These ID's must correspond with the block ID's of the element blocks
	 *        in the global geometry database.
	 * \param connectivities elements over which data will be written. The
	 *        output formatters retain a reference to these connectivities
	 *        for use during output. These connectivities are not copied.
	 * \param n_labels list of labels for the nodal variables. The length of
	 *        this list defines the number of nodal output variables.
	 * \param e_labels list of labels for the element variables. The length of
	 *        this list defines the number of element output variables.
	 * \param changing flag to indicate whether the connectivities may change
	 *        from output step to output step. */
	OutputSetT(GeometryT::CodeT geometry_code,
		const ArrayT<StringT>& block_ID, const ArrayT<const iArray2DT*>& connectivities, 
		const ArrayT<StringT>& n_labels, const ArrayT<StringT>& e_labels, 
		bool changing);

	/** generate output data record.
	 * \param geometry_code GeometryT::CodeT defining the geometry associated
	 *        with the connectivities.
	 * \param connectivities elements over which data will be written. The
	 *        output formatters retain a reference to these connectivities
	 *        for use during output. These connectivities are not copied.
	 * \param n_labels list of labels for the nodal variables. The length of
	 *        this list defines the number of nodal output variables
	 * \param changing flag to indicate whether the connectivities may change
	 *        from output step to output step. */
	OutputSetT(GeometryT::CodeT geometry_code,
		const iArray2DT& connectivities, const ArrayT<StringT>& n_labels, bool changing = false);

	/** output data record for a set of points
	 * \param points points over which data will be written. The
	 *        output formatters retain a reference to these points
	 *        for use during output. This array is not copied.
	 * \param n_labels list of labels for the nodal variables. The length of
	 *        this list defines the number of nodal output variables
	 * \param changing flag to indicate whether the connectivities may change
	 *        from output step to output step. */
	OutputSetT(const iArrayT& points, const ArrayT<StringT>& n_labels, bool changing = false);

	/** generate output data record.
	 * \param geometry_code GeometryT::CodeT defining the geometry associated
	 *        with the connectivities.
	 * \param block_ID list of element block ID's. Presumably, these are NOT 
	 *        in the global geometry database.
	 * \param sideset_ID list of side set IDs that, presumably, ARE in the global
	 *        geometry database
	 * \param connectivities elements over which data will be written. The
	 *        output formatters retain a reference to these connectivities
	 *        for use during output. These connectivities are not copied.
	 * \param n_labels list of labels for the nodal variables. The length of
	 *        this list defines the number of nodal output variables.
	 * \param e_labels list of labels for the element variables. The length of
	 *        this list defines the number of element output variables.
	 * \param changing flag to indicate whether the connectivities may change
	 *        from output step to output step. */
	OutputSetT(GeometryT::CodeT geometry_code,
		const ArrayT<StringT>& block_ID, const ArrayT<StringT>& sideset_ID,
		const ArrayT<const iArray2DT*>& connectivities, 
		const ArrayT<StringT>& n_labels, const ArrayT<StringT>& e_labels, 
		bool changing);

	/** copy constructor */
	OutputSetT(const OutputSetT& source);
	
	/** output set mode */
	ModeT Mode(void) const { return fMode; }; 
	
	/* dimensions */
	int NumNodes(void) const; /**< return the number of nodes used by the set */
	int NumBlocks (void) const;	/**< return the number of connectivity blocks in the set */
	int NumBlockElements(const StringT& ID) const; /**< return the number of elements in the specified block */
	int NumElements(void) const; /**< return the total number of elements */
	int NumNodeValues(void) const; /** return the number of nodal output variables */
	int NumElementValues(void) const; /** return the number of element output variables */

	/* print step counter */
	int PrintStep(void) const;
	void ResetPrintStep(void);
	void IncrementPrintStep(void);

	/** return the ID for the output set */
	void SetID(const StringT& id);
	const StringT& ID(void) const { return fID; };

	/** \name changing geometry flag.
	 * set/get the changing geometry flag */
	/*@{*/
	/** return true if the set has changing geometry, false otherwise */
	bool Changing(void) const;

	/** set the flag */
	void SetChanging(bool changing) { fChanging = changing; };
	/*@}*/

	/** return the GeometryT::CodeT for the output set */
	GeometryT::CodeT Geometry(void) const;

	/** return the list of element block ID's used by the set */
	const ArrayT<StringT>& BlockID(void) const { return fBlockID; };
	
	/** return the list of side set ID's used by the set */
	const ArrayT<StringT>& SideSetID(void) const { return fSSID; };
	
	/** return the ID of the specified block */
	const StringT& BlockID(int index) const;
	
	/** return the ID of the specified side set */
	const StringT& SideSetID(int index) const;

	/** return a pointer to the connectivities for the specified block */
	const iArray2DT* Connectivities(const StringT& ID) const;

	/** return the labels for the nodal output variables */
	const ArrayT<StringT>& NodeOutputLabels(void) const;

	/** return the labels for the element output variables */
	const ArrayT<StringT>& ElementOutputLabels(void) const;

	/** return the nodes used by the output set. If the geometry is
	 * changing, the nodes used are recalculated with every call. 
	 * \note this function isn't really const since the nodes used
	 * may be recalculated, but this is the price of lazy evaluation. */
	const iArrayT& NodesUsed(void) const;

	/** return the nodes used by the given block. If the geometry is
	 * changing, the nodes used are recalculated with every call.
	 * \note this function isn't really const since the nodes used
	 * may be recalculated, but this is the price of lazy evaluation. */
	const iArrayT& BlockNodesUsed(const StringT& ID) const;

	/** block index to set index node map. For cases with sets with
	 * changing geometry, this map corresponds to the configuration
	 * during the last call to OutputSetT::BlockNodesUsed. */
	const iArrayT& BlockIndexToSetIndexMap(const StringT& ID) const;

	/** returns the index for the element block for the given */
	int BlockIndex(const StringT& ID) const;

private:

	/** make private to avoid accidental use. */
	OutputSetT& operator=(OutputSetT&) { return *this; }

	/** determine the nodes used.
	 * \param connects list of nodes
	 * \param nodes_used destination for nodes used 
	 * \param nodes_used_man memory manager for the nodes used array */
	void SetNodesUsed(const iArray2DT& connects, iArrayT& nodes_used, 
		VariArrayT<int>& nodes_used_man) const;

	/** determine the nodes used.
	 * \param connects_list list of lists of nodes
	 * \param nodes_used destination for nodes used 
	 * \param nodes_used_man memory manager for the nodes used array */
	void SetNodesUsed(const ArrayT<const iArray2DT*>& connects_list, 
		iArrayT& nodes_used, VariArrayT<int>& nodes_used_man);
	
private:

	/** output set mode */
	ModeT fMode;

	/** count of number of output steps */
	int fPrintStep;

	/** set ID */
	StringT fID;
	
	/** true if set has changing geometry, false otherwise */
	bool fChanging;
	
	/** geometry for the connectivities in the set */
	GeometryT::CodeT fGeometry;

	/** list of ID's for the connectivities */
	ArrayT<StringT> fBlockID;

	/** pointers to the connectivity data */
	ArrayT<const iArray2DT*> fConnectivities;
	
	/** list of ID's for the side sets of each element block */
	ArrayT<StringT> fSSID;
	
	/** labels for nodal output variables */
	ArrayT<StringT> fNodeOutputLabels;

	/** labels for element output variables */
	ArrayT<StringT> fElementOutputLabels;
	
	/** \name cached */
	/*@{*/
	iArrayT fNodesUsed; /**< nodes used by the whole output set */
	ArrayT<iArrayT> fBlockNodesUsed; /**< nodes used by element block */
	ArrayT<iArrayT> fBlockIndexToSetIndexMap; 
	/**< map telling for the ith block\n
	 * fBlockNodesUsed_i[j] = fNodesUsed[fBlockIndexToSetIndexMap_i[j]] */
	/*@}*/

	/** \name writing over set of points */
	/*@{*/
	const iArrayT* fPoints;
	
	/** dummy 2D connectivity used when writing data over an iArrayT */
	iArray2DT fConnects2D;
	/*@}*/

	/** \name memory managers */
	/*@{*/
	VariArrayT<int> fNodesUsed_man;
	ArrayT<VariArrayT<int> > fBlockNodesUsed_man;
	ArrayT<VariArrayT<int> > fBlockIndexToSetIndexMap_man;
	/*@}*/
};

/* inlines */
inline int OutputSetT::PrintStep(void) const { return fPrintStep; }
inline void OutputSetT::ResetPrintStep(void) { fPrintStep = -1; }
inline void OutputSetT::IncrementPrintStep(void) { fPrintStep++; }

inline bool OutputSetT::Changing(void) const { return fChanging; }
inline GeometryT::CodeT OutputSetT::Geometry(void) const { return fGeometry; }
inline int OutputSetT::NumBlocks (void) const { return fBlockID.Length(); }
inline const iArrayT& OutputSetT::NodesUsed(void) const
{
	if (fChanging) {
		
		/* not so const */
		OutputSetT* non_const_this = (OutputSetT*) this;

		/* reset alias */
		if (fPoints) non_const_this->fConnects2D.Alias(fPoints->Length(), 1, fPoints->Pointer());

		/* reset nodes used */
		non_const_this->SetNodesUsed(fConnectivities, non_const_this->fNodesUsed, non_const_this->fNodesUsed_man);
	}
	return fNodesUsed;
}

inline const ArrayT<StringT>& OutputSetT::NodeOutputLabels(void) const { return fNodeOutputLabels; }

inline const ArrayT<StringT>& OutputSetT::ElementOutputLabels(void) const { return fElementOutputLabels; }

/* dimensions */
inline int OutputSetT::NumNodes(void) const { return fNodesUsed.Length(); }
inline int OutputSetT::NumNodeValues(void) const { return fNodeOutputLabels.Length(); }
inline int OutputSetT::NumElementValues(void) const { return fElementOutputLabels.Length(); }

/* block index to set index node map */
inline const iArrayT& OutputSetT::BlockIndexToSetIndexMap(const StringT& ID) const
{ 
	return fBlockIndexToSetIndexMap[BlockIndex(ID)]; 
};

/* return the ID of the specified block */
inline const StringT& OutputSetT::BlockID(int index) const
{
	if (index < 0 || index >= fBlockID.Length())
		ExceptionT::OutOfRange("OutputSetT::BlockID", "%d out of range %d", index);
	return fBlockID[index];

}

/* return the ID of the specified side set */
inline const StringT& OutputSetT::SideSetID(int index) const
{
	if (index < 0 || index >= fSSID.Length())
		ExceptionT::OutOfRange("OutputSetT::SideSetID", "%d out of range", index);
	return fSSID[index];

}

} // namespace Tahoe 
#endif /* _OUTPUTSET_T_H_ */
