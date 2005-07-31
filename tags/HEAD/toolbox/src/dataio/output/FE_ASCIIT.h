/* $Id: FE_ASCIIT.h,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
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

	/* constructor */
	FE_ASCIIT(ostream& out, bool external, const ArrayT<StringT>& out_strings);

	/* increment sequence, create new output file series */
	virtual void NextTimeSequence(int sequence_number);

	/* output functions */
	virtual void WriteGeometry(void);
	virtual void WriteOutput(double time, int ID, const dArray2DT& n_values,
		const dArray2DT& e_values);

private:

	/* set-by-set output */
	void WriteGeometryData(ostream& out, int ID);
	void WriteOutputData(ostream& out, int ID, const dArray2DT& n_values,
		const dArray2DT& e_values);

	/* variable headers */
	void WriteNodeHeader(ostream& out, int num_output_nodes,
		const ArrayT<StringT>& labels) const;
	void WriteElementHeader(ostream& out, int num_output_elems,
		const ArrayT<StringT>& labels) const;

	/* writing data */
	void WriteNodeValues(ostream& out, const iArrayT& node_numbers,
		const dArray2DT& values) const;
	void WriteElementValues(ostream& out, const dArray2DT& values) const;

private:

	bool fExternTahoeII;
	bool fInitGeom;
	bool fInitRun;
};

#endif // _FE_ASCII_T_H_
