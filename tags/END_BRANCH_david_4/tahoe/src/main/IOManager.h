/* $Id: IOManager.h,v 1.20 2005-03-12 08:41:35 paklein Exp $ */
/* created: sawimme (10/12/1999) */
#ifndef _IOMANAGER_H_
#define _IOMANAGER_H_

/* language support */
#include <iostream.h>
#include <fstream.h>

/* direct members */
#include "StringT.h"
#include "iArrayT.h"
#include "iAutoArrayT.h"
#include "IOBaseT.h"
#include "GeometryT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;
class OutputBaseT;
class iArray2DT;
class dArray2DT;
class dArrayT;
class OutputSetT;

/** class to handle format-independent output. Data for output
 * is registered through IOManager::AddElementSet and written with
 * successive calls to IOManager::WriteOutput */
class IOManager
{
public:

	/** changing geometry */
	enum ChangingFlagT {
		kNoChangingFlag,
		kForceChanging,
		kForceNotChanging
	};

	/** constructor */
	IOManager(ostream& log, const StringT& program_name, const StringT& version,
		const StringT& title, const StringT& input_file, IOBaseT::FileTypeT output_format);

	/** destructor */
	virtual ~IOManager(void);

	/** the output formatter */
	const OutputBaseT& Output(void) const { return *fOutput; };

	/** the output format */
	IOBaseT::FileTypeT OutputFormat(void) const { return fOutputFormat; };
	
	/** log stream */
	ostream& Log(void) const { return fLog; };

	/** how to override changing geometry flags in registered output sets */
	void SetChangingFlag(ChangingFlagT changing_flag);

	/** echo format file type */
	void EchoData (ostream& o) const;

	/** increment the time sequence */
	void NextTimeSequence(int sequence_number);

	/** set nodal coordinates
	 * \param coordinate array of nodal coordinates
	 * \param node_id list of ids for the rows in the coordinate array, passing NULL
	 *        implies the row number is the id. The number of ids must match the
	 *        number of rows in the coordinate array */
	void SetCoordinates(const dArray2DT& coordinates, const iArrayT* node_id);
	
	/** register the output for an element set. returns the output ID */
	int AddElementSet(const OutputSetT& output_set);

	/** return the array of output set specifiers */
	const ArrayT<OutputSetT*>& ElementSets(void) const;

	/** add a node set to the results file. Not all output formats support
	 * defining node sets in the results. No data is written to the node
	 * set. Data can be written to a node set by registering it with
	 * IOManager::AddElementSet */
	void AddNodeSet(const iArrayT& nodeset, const StringT& setID);

	/** add a side set to the results file. Not all output formats support
	 * defining side sets in the results. No data is written to the side
	 * set. Data can be written to a side set by registering it with
	 * IOManager::AddElementSet */
	void AddSideSet(const iArray2DT& sideset, const StringT& setID, const StringT& group_ID);

	/* output functions */
	void WriteGeometry(void);
	void WriteGeometryFile(const StringT& file_name, IOBaseT::FileTypeT format) const;
	
	/** define the time for the current output step */
	void SetOutputTime(double time);

	/** write data to output
	 * \param ID set ID returned from the call to IOManager::AddElementSet
	 * \param n_values nodal output values ordered as given by IOManager::NodesUsed
	 * \param e_values element output values in the the ordered defined in the connectivities */
	virtual void WriteOutput(int ID, const dArray2DT& n_values, const dArray2DT& e_values);

	/** write data to output
	 * \param ID set ID returned from the call to IOManager::AddElementSet
	 * \param n_values nodal output values ordered as given by IOManager::NodesUsed */
	virtual void WriteOutput(int ID, const dArray2DT& n_values);

	/** write a snapshot */
	void WriteOutput(const StringT& file, const dArray2DT& coords, const iArrayT& node_map,
		const dArray2DT& values, const ArrayT<StringT>& labels) const;

	/** return the list of nodes used by the output set
	 * \param ID set ID returned from the call to IOManager::AddElementSet */
	const iArrayT& NodesUsed(int ID) const;

	/** temporarily re-route output to a database with the given filename */
	virtual void DivertOutput(const StringT& outfile);

	/** restore output to the default stream */
	void RestoreOutput(void);

	/** return the output set with the given ID
	 * \param ID set ID returned from the call to IOManager::AddElementSet */
	const OutputSetT& OutputSet(int ID) const;

protected:

	ostream& fLog;
	StringT fTitle;

	/* output formatter */
	IOBaseT::FileTypeT fOutputFormat;
	OutputBaseT* fOutput;

	/* echo interactive data to input file */
	ofstream fEchoInput;
	bool fEcho;

	/** flag for how to override changing geometry flags */
	ChangingFlagT fChangingFlag;
	
private:

	/** time label associated with the data written to IOManager::fOutput.
	 * This value is set using IOManager::SetOutputTime. This value is
	 * not directly tied to the simulation run time because this causes
	 * problems with multiple output for the same time, i.e., iteration
	 * output. In this case, some post-processors choke on the fact that
	 * the time stamp is duplicated. The remedy is set some changing
	 * stamp for the time. */
	double fOutputTime;

	/** store main out during diversions set with IOManager::DivertOutput */
	OutputBaseT* fOutput_tmp;	
};

/* inlines */
inline void IOManager::SetOutputTime(double time) { fOutputTime = time; }

} /* namespace Tahoe */

#endif
