/* $Id: FE_ASCIIT.h,v 1.1.1.1.2.2 2001-11-06 20:38:33 sawimme Exp $ */
/* created: sawimme (05/20/1999)                                          */

#ifndef _FE_ASCII_T_H_
#define _FE_ASCII_T_H_

/* base class */
#include "OutputBaseT.h"

/* direct members */
#include <fstream.h>

class FE_ASCIIT: public OutputBaseT
{
public:

  /** constructor
   * \param out error stream
   * \param external flag to write array data in external files
   * \param out_strings see OutputBaseT::OutputBaseT */
	FE_ASCIIT(ostream& out, bool external, const ArrayT<StringT>& out_strings);

	/** increment sequence, create new output file series */
	virtual void NextTimeSequence(int sequence_number);

	/** write geometry for all output sets */
	virtual void WriteGeometry(void);

	/** write geometry and node variables for output set ID */
	virtual void WriteOutput(double time, int ID, const dArray2DT& n_values,
		const dArray2DT& e_values);

private:

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

	/** writing node array data */
	void WriteNodeValues(ostream& out, const iArrayT& node_numbers,
		const dArray2DT& values) const;
	/** writing element array data */
	void WriteElementValues(ostream& out, const dArray2DT& values) const;

private:

	bool fExternTahoeII; /**< flag to write array data to external files */
	bool fInitGeom; /**< flag to determine if appending or creating a file */
	bool fInitRun; /**< flag to determine if appending or creating a file */
};

#endif // _FE_ASCII_T_H_
