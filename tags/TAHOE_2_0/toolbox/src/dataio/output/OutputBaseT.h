/* $Id: OutputBaseT.h,v 1.19 2003-10-28 07:50:41 paklein Exp $ */
/* created: sawimme (05/18/1999) */
#ifndef _OUTPUTBASE_T_H_
#define _OUTPUTBASE_T_H_

/* base class */
#include "IOBaseT.h"

/* direct members */
#include "StringT.h"
#include "iArrayT.h"
#include "iAutoArrayT.h"
#include "GeometryT.h"

namespace Tahoe {

/* forward declarations */
class dArray2DT;
class iArray2DT;
class OutputSetT;

/** base class for output formatter.
 * initialization:\n
 * -# construct
 * -# SetCoordinates
 * -# AddElementSet 
 */
class OutputBaseT: public IOBaseT
{
public:
	
	/** constructor */
	OutputBaseT(ostream& out, const ArrayT<StringT>& outstrings);

	/** destructor */
	~OutputBaseT(void);

	/** \name accessors */
	/*@{*/
	const StringT& Title(void) const;
	const StringT& CodeName(void) const;
	const StringT& Version(void) const;
	const StringT& OutputRoot(void) const;
	const OutputSetT& OutputSet(int ID) const;
	/*@}*/
	
	/** return the array of nodes used by the specified output set
	 * \param ID set ID returned from the call to OutputBaseT::AddElementSet */
	 const iArrayT& NodesUsed(int ID) const;

	/** increment sequence, create new output file series */
	virtual void NextTimeSequence(int sequence_number);

	/** set nodal coordinates for the output set
	 * \param coordinate array of nodal coordinates
	 * \param node_id list of ids for the rows in the coordinate array, passing NULL
	 *        implies the row number is the id. The number of ids must match the
	 *        number of rows in the coordinate array */
	virtual void SetCoordinates(const dArray2DT& coordinates, const iArrayT* node_id);

	virtual void SetBounds(const dArray2DT& bounds);
	virtual void SetTypes(const iArrayT& types);
	virtual void SetParts(const iArrayT& parts);

	/** return the node id list
	 * \return point to the id list. If NULL implies the row number in the coordinate
	 *         array is the node id */
	const iArrayT* NodeID(void) const;

	/** return a reference to the coordinates */
	const dArray2DT& Coordinates(void) const;

	/** register the output for an element set. returns the output ID */
	virtual int AddElementSet(const OutputSetT& output_set);
	const ArrayT<OutputSetT*>& ElementSets(void) const;
	int NumElements(void) const;

	void AddNodeSet(const iArrayT& nodeset, const StringT& setID);
	void AddSideSet(const iArray2DT& sideset, const StringT& setID, const StringT& group_ID);

	/** \name output methods */
	/*@{*/
	virtual void WriteGeometry(void) = 0;
	void WriteGeometryFile(const StringT& file_name, IOBaseT::FileTypeT format) const;

	/** send data for output. Order of the nodal values is set by OutputSetT::NodesUsed. */
	virtual void WriteOutput(double time, int ID, const dArray2DT& n_values, const dArray2DT& e_values);

	/** send data for output */
	virtual void WriteOutput(double time, int ID, const ArrayT<int>& nodes, const dArray2DT& n_values, 
		const ArrayT<int>& elements, const dArray2DT& e_values);
	/*@}*/

protected:

	enum DataTypeT {kNode = 0,
	             kElement = 1};

	/** convert group node numbering to block node numbering within a connectivity */
	void LocalConnectivity(const iArrayT& node_map, const iArray2DT& connects, iArray2DT& local_connects) const;

	/** separate group element data into block element data */
	void ElementBlockValues (int ID, int block, const dArray2DT& allvalues, dArray2DT& blockvalues) const;

	/** separate group nodal data into block nodal data */
	void NodalBlockValues (int ID, int block, const dArray2DT& allvalues, dArray2DT& blockvalues, iArrayT& block_nodes) const;

	/** search for the block string ID
	 * \param  s string block ID
	 * \return g index to the group containing the block 
	 * \return b index to the block within the group */
	void ElementGroupBlockIndex (const StringT& s, int& g, int& b) const;

	/** create fElementBlockIDs by attempting to convert string block IDs directly
	    to integer IDs. If IDs are not unique, then blocks are globally numbered
	    starting at 1. */
	void CreateElementBlockIDs (void);

	/** converts side set or node set string ID values to integer values if possible,
	    if not unique, blocks are globally numbered starting at 1. */
	void String2IntIDs (const ArrayT<StringT>& s, iArrayT& i) const;

protected:

	StringT fTitle;    /**< title: description of problem */
	StringT fCodeName; /**<	qa_record codename and version */
	StringT fVersion;  /**<	qa_record inputfile version */
	StringT fOutroot;  /**<	root of all output files including the path */

	/* output data */
	const dArray2DT* fCoordinates; /**< pointer to coordinates */
	const iArrayT*   fNodeID; /**< list of node IDs */

	const dArray2DT* fBounds; /**< pointer to bounds, paradyn format */
	const iArrayT* fTypes; /**< pointer to types, paradyn format */
	const iArrayT* fParts; /**< pointer to parts, paradyn format */

	AutoArrayT<OutputSetT*> fElementSets; /**< list of output groups */
	ArrayT<iArrayT> fElementBlockIDs; /**< integer blocks IDs created from string IDs */
	
	AutoArrayT<const iArrayT*>   fNodeSets; /**< links to node sets */
	AutoArrayT<const iArray2DT*> fSideSets; /**< links to side sets */
	AutoArrayT<StringT>          fNodeSetNames; /**< node set string IDs */
	AutoArrayT<StringT>          fSideSetNames; /**< side set string IDs */
	AutoArrayT<StringT>          fSSGroupNames; /**< block name string ID that side set is in */

	int fSequence; /**< solution sequence number */
};

inline const ArrayT<OutputSetT*>& OutputBaseT::ElementSets(void) const
{
	return fElementSets;
}

/* accessors */
inline const StringT& OutputBaseT::Version(void) const { return fVersion; }
inline const StringT& OutputBaseT::CodeName(void) const { return fCodeName; }
inline const StringT& OutputBaseT::Title(void) const { return fTitle; }
inline const StringT& OutputBaseT::OutputRoot(void) const { return fOutroot; }
inline const dArray2DT& OutputBaseT::Coordinates(void) const
{
	if (!fCoordinates)
		ExceptionT::GeneralFail("OutputBaseT::Coordinates", "pointer to coordinates not set");
	return *fCoordinates;
}

inline const iArrayT* OutputBaseT::NodeID(void) const { return fNodeID; }

} // namespace Tahoe 
#endif /* _OUTPUTMANAGER_H_ */
