/* $Id: TextOutputT.h,v 1.4 2005-06-05 06:22:00 paklein Exp $ */
/* created: sawimme (05/20/1999) */
#ifndef _TEXT_OUTPUT_T_H_
#define _TEXT_OUTPUT_T_H_

/* base class */
#include "OutputBaseT.h"

/* direct members */
#include <fstream.h>

namespace Tahoe {

/** text file output */
class TextOutputT: public OutputBaseT
{
public:

  /** constructor
   * \param out error stream
   * \param external flag to write array data in external files
   * \param out_strings see OutputBaseT::OutputBaseT */
	TextOutputT(ostream& out, bool external, const ArrayT<StringT>& out_strings);

	/** register the output for an element set. returns the output ID */
	virtual int AddElementSet(const OutputSetT& output_set);

	/** increment sequence, create new output file series */
	virtual void NextTimeSequence(int sequence_number);

	/** write geometry for all output sets */
	virtual void WriteGeometry(void);

	/** write geometry and node variables for output set ID */
	virtual void WriteOutput(double time, int ID, const dArray2DT& n_values,
		const dArray2DT& e_values);

	/** \name writing values */
	/*@{*/
	/** writing node array data */
	static void WriteNodeValues(ostream& out, const ArrayT<int>& node_numbers, const dArray2DT& values);

	/** writing element array data */
	static void WriteElementValues(ostream& out, const dArray2DT& values);
	/*@}*/

private:

	/** initialize the results file */
	void InitResultsFile(ostream& out, int ID);

	/** set-by-set geometry output */
	void WriteGeometryData(ostream& out, int ID);

	/** set-by-set variable output */
	void WriteOutputData(ostream& out, int ID, const dArray2DT& n_values,
		const dArray2DT& e_values);

		/** write node variable headers */
	void WriteNodeHeader(ostream& out, int num_output_nodes,
		const ArrayT<StringT>& labels) const;
		/** write element variable headers */
	void WriteElementHeader(ostream& out, int num_output_elems,
		const ArrayT<StringT>& labels) const;

private:

	bool fExternTahoeII; /**< flag to write array data to external files */
	AutoArrayT<bool> fInitGeom; /**< flags to determine if appending or creating a file */
	AutoArrayT<bool> fInitRun; /**< flags to determine if appending or creating a file */
	iArrayT fNodeSetIntIDs;
	iArrayT fSideSetIntIDs;
};

} // namespace Tahoe 
#endif // _TEXT_OUTPUT_T_H_
