/* $Id: OutputSetT.h,v 1.5.2.3 2001-10-29 21:10:27 sawimme Exp $ */
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
	 * \param block_ID list of element block ID's comprising the connectivities.
	 *        This list is used to declare that the connectivities passed in are
	 *        actually a union of connectivities defined by ID's declared in an 
	 *        associated geometry database file. This information is needed to
	 *        join element-based data from parallel calculations. This list may
	 *        be empty.
	 * \param connectivities elements over which data will be written. The
	 *        output formatters retain a reference to these connectivities
	 *        for use during output.
	 * \param n_labels list of labels for the nodal variables. The length of
	 *        this list defines the number of nodal output variables.
	 * \param e_labels list of labels for the element variables. The length of
	 *        this list defines the number of element output variables.
	 * \param changing flag to indicate whether the connectivities may change
	 *        from output step to output step. */
	OutputSetT(int ID, GeometryT::CodeT geometry_code,
		const iArrayT& block_ID, const ArrayT<const iArray2DT*> connectivities, 
		const ArrayT<StringT>& n_labels, const ArrayT<StringT>& e_labels, 
		bool changing);

	/** copy constructor */
	OutputSetT(const OutputSetT& source);

	/* dimensions */
	int NumNodes(void) const; /**< return the number of nodes used by the set */
	int NumBlocks (void) const;	/**< return the number of connectivity blocks in the set */
	int NumBlockElements(int index) const; /**< return the number of elements in the specified block */
	int NumElements(void) const; /**< return the total number of elements */
	int NumNodeValues(void) const; /** return the number of nodal output variables */
	int NumElementValues(void) const; /** return the number of element output variables */

	/* print step counter */
	int PrintStep(void) const;
	void ResetPrintStep(void);
	void IncrementPrintStep(void);

	/** return the ID for the output set */
	int ID(void) const;

	/** return true if the set has changing geometry, false otherwise */
	bool Changing(void) const;

	/** return the GeometryT::CodeT for the output set */
	GeometryT::CodeT Geometry(void) const;

	/** return the list of element block ID's used by the set */
	const iArrayT& BlockID(void) const { return fBlockID; };

	/** return a pointer to the connectivities for the specified block */
	const iArray2DT* Connectivities(int index) const;

//TEMP - used to write all set connectivities at once
	//#if 0
	void AllConnectivities(iArray2DT& connects) const;
	//#endif

	/** return the labels for the nodal output variables */
	const ArrayT<StringT>& NodeOutputLabels(void) const;

	/** return the labels for the element output variables */
	const ArrayT<StringT>& ElementOutputLabels(void) const;

	/** return the nodes used by the output set. If the geometry if
	 * changing, the nodes used are recalculated with every call. */
	const iArrayT& NodesUsed(void);

private:

	/** make private to avoid accidental use. */
	OutputSetT& operator=(OutputSetT&) { return *this; }

	/** determine the nodes used in the output set */
	void SetNodesUsed(void);

private:

	/** count of number of output steps */
	int fPrintStep;

	int  fID;
	bool fChanging;
	GeometryT::CodeT fGeometry;

	/* connectivities */
	const iArrayT fBlockID;
	const ArrayT<const iArray2DT*> fConnectivities;
	
	/* output labels */
	ArrayT<StringT> fNodeOutputLabels;
	ArrayT<StringT> fElementOutputLabels;
	
	/* derived */
	iArrayT fNodesUsed;
};

/* inlines */
inline int OutputSetT::PrintStep(void) const { return fPrintStep; }
inline void OutputSetT::ResetPrintStep(void) { fPrintStep = -1; }
inline void OutputSetT::IncrementPrintStep(void) { fPrintStep++; }

inline int OutputSetT::ID(void) const { return fID; }
inline bool OutputSetT::Changing(void) const { return fChanging; }
inline GeometryT::CodeT OutputSetT::Geometry(void) const { return fGeometry; }
inline int OutputSetT::NumBlocks (void) const { return fBlockID.Length(); }
inline const iArrayT& OutputSetT::NodesUsed(void)
{
	if (fChanging) SetNodesUsed();
	return fNodesUsed;
}

inline const ArrayT<StringT>& OutputSetT::NodeOutputLabels(void) const { return fNodeOutputLabels; }

inline const ArrayT<StringT>& OutputSetT::ElementOutputLabels(void) const { return fElementOutputLabels; }

/* dimensions */
inline int OutputSetT::NumNodes(void) const { return fNodesUsed.Length(); }
inline int OutputSetT::NumNodeValues(void) const { return fNodeOutputLabels.Length(); }
inline int OutputSetT::NumElementValues(void) const { return fElementOutputLabels.Length(); }

#endif /* _OUTPUTSET_T_H_ */
