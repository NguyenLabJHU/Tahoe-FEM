/* $Id: OutputSetT.h,v 1.8 2002-02-07 23:28:33 paklein Exp $ */
/* created: paklein (03/07/2000) */

#ifndef _OUTPUTSET_T_H_
#define _OUTPUTSET_T_H_

/* direct members */
#include "GeometryT.h"
#include "StringT.h"
#include "iArrayT.h"

/* forward declarations */
class iArray2DT;

/** class to act as specifier for output data. A class that would
 * like to write output, needs to use the following procedure:\n
 * (1) construct an OutputSetT with its output data specifications\n
 * (2) register this output data specification with an output formatter.
 *     Output formatters are classes derived from OutputBaseT. Output
 *     is registered using OutputBaseT::AddElementSet. Note that this
 *     function is accessed in Tahoe through FEManagerT::RegisterOutput
 *     which passes the call through IOManager::AddElementSet to
 *     OutputBaseT::AddElementSet. Registering output returns an integer
 *     output ID (OID) which must be used in step (3). Note that the OID 
 *     obtained during registration is generally not the same as the ID 
 *     passed into the OutputSetT constructor. A class with a single ID 
 *     may register an arbitary number of output sets, receiving a unique
 *     OID for each set\n
 * (3) write output by sending output data to the output formatter using
 *     OutputBaseT::WriteOutput with the ID obtained during registration
 *     in step (2). Note that this function is accessed in Tahoe through 
 *     FEManagerT::WriteOutput which passes the call through 
 *     IOManager::AddElementSet to the output formatter\n
 *
 * Outside of the constructor, the remainder of the interface is intended
 * for use by the output formatter, not the class sending the data for
 * output\n */
class OutputSetT
{
public:

	/** generate output data record.
	 * \param ID identifier to denote the class generating the OutputSetT
	 * \param geometry_code GeometryT::CodeT defining the geometry associated
	 *        with the connectivities.
	 * \param block_ID list of ID's comprising the connectivities. This list is 
	 *        used to declare that the connectivities passed in are actually a union 
	 *        of connectivities defined by ID's declared in an associated geometry 
	 *        database file. These values are used internally to refer to the data 
	 *        associated with each member of the connectivities list. If data is only
	 *        written per node, these ID's have no other use. If there is data per
	 *        element, these values are also used by JoinOutputT when joining "element"
	 *        data, in which case the ID must correspond to the ID for the element
	 *        block in the geometry database.
	 * \param connectivities elements over which data will be written. The
	 *        output formatters retain a reference to these connectivities
	 *        for use during output.
	 * \param n_labels list of labels for the nodal variables. The length of
	 *        this list defines the number of nodal output variables.
	 * \param e_labels list of labels for the element variables. The length of
	 *        this list defines the number of element output variables.
	 * \param changing flag to indicate whether the connectivities may change
	 *        from output step to output step. */
	OutputSetT(const StringT& ID, GeometryT::CodeT geometry_code,
		const ArrayT<StringT>& block_ID, const ArrayT<const iArray2DT*>& connectivities, 
		const ArrayT<StringT>& n_labels, const ArrayT<StringT>& e_labels, 
		bool changing);

	/** copy constructor */
	OutputSetT(const OutputSetT& source);
	
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
	const StringT& ID(void) const;

	/** return true if the set has changing geometry, false otherwise */
	bool Changing(void) const;

	/** return the GeometryT::CodeT for the output set */
	GeometryT::CodeT Geometry(void) const;

	/** return the list of element block ID's used by the set */
	const ArrayT<StringT>& BlockID(void) const { return fBlockID; };
	
	/** return the ID of the specified block */
	const StringT& BlockID(int index) const;

	/** return a pointer to the connectivities for the specified block */
	const iArray2DT* Connectivities(const StringT& ID) const;

//TEMP - used to write all set connectivities at once
#if 0
	void AllConnectivities(iArray2DT& connects) const;
#endif

	/** return the labels for the nodal output variables */
	const ArrayT<StringT>& NodeOutputLabels(void) const;

	/** return the labels for the element output variables */
	const ArrayT<StringT>& ElementOutputLabels(void) const;

	/** return the nodes used by the output set. If the geometry is
	 * changing, the nodes used are recalculated with every call. */
	const iArrayT& NodesUsed(void);

	/** return the nodes used by the given block. If the geometry is
	 * changing, the nodes used are recalculated with every call. */
	const iArrayT& BlockNodesUsed(const StringT& ID);

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
	 * \param destination for nodes used */
	void SetNodesUsed(const iArray2DT& connects, iArrayT& nodes_used);

	/** determine the nodes used.
	 * \param connects_list list of lists of nodes
	 * \param destination for nodes used */
	void SetNodesUsed(const ArrayT<const iArray2DT*>& connects_list, 
		iArrayT& nodes_used);
	
private:

	/** count of number of output steps */
	int fPrintStep;

	/** set ID */
	const StringT fID;
	
	/** true if set has changing geometry, false otherwise */
	bool fChanging;
	
	/** geometry for the connectivities in the set */
	GeometryT::CodeT fGeometry;

	/** list of ID's for the connectivities */
	const ArrayT<StringT> fBlockID;

	/** pointers to the connectivity data */
	const ArrayT<const iArray2DT*> fConnectivities;
	
	/** labels for nodal output variables */
	ArrayT<StringT> fNodeOutputLabels;

	/** labels for element output variables */
	ArrayT<StringT> fElementOutputLabels;
	
	/* cached */
	iArrayT fNodesUsed; /**< nodes used by the whole output set */
	ArrayT<iArrayT> fBlockNodesUsed; /**< nodes used by element block */
	ArrayT<iArrayT> fBlockIndexToSetIndexMap; 
	/**< map telling for the ith block\n
	 * fBlockNodesUsed_i[j] = fNodesUsed[fBlockIndexToSetIndexMap_i[j]] */
};

/* inlines */
inline int OutputSetT::PrintStep(void) const { return fPrintStep; }
inline void OutputSetT::ResetPrintStep(void) { fPrintStep = -1; }
inline void OutputSetT::IncrementPrintStep(void) { fPrintStep++; }

inline const StringT& OutputSetT::ID(void) const { return fID; }
inline bool OutputSetT::Changing(void) const { return fChanging; }
inline GeometryT::CodeT OutputSetT::Geometry(void) const { return fGeometry; }
inline int OutputSetT::NumBlocks (void) const { return fBlockID.Length(); }
inline const iArrayT& OutputSetT::NodesUsed(void)
{
	if (fChanging) SetNodesUsed(fConnectivities, fNodesUsed);
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
	if (index < 0 || index >= fBlockID.Length()) {
		cout << "\n OutputSetT::BlockID: index out of range: " << index << endl;
		throw eOutOfRange;
	}
	return fBlockID[index];

}

#endif /* _OUTPUTSET_T_H_ */
