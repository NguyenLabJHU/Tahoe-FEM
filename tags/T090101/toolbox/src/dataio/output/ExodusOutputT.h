/* $Id: ExodusOutputT.h,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
/* created: sawimme (05/18/1999)                                          */

#ifndef _EXODUSOUTPUT_T_H_
#define _EXODUSOUTPUT_T_H_

/* base class */
#include "OutputBaseT.h"

/* forward declarations */
class ExodusT;

class ExodusOutputT: public OutputBaseT
{
public:
	ExodusOutputT(ostream& out, const ArrayT<StringT>& out_strings);

	/* output functions */
	virtual void WriteGeometry(void);
	virtual void WriteOutput(double time, int ID, const dArray2DT& n_values,
		const dArray2DT& e_values);

private:

	/* generate database file name for the given ID */
	void FileName(int ID, StringT& filename) const;

	/* create files */
	void CreateResultsFile(int ID, ExodusT& exo);
	void CreateGeometryFile(ExodusT& exo);

	// items common to CreateResultsFile and CreateGeometryFile
	void AssembleQA (ArrayT<StringT>& qa) const;
	void WriteCoordinates (ExodusT& exo, iArrayT& nodes_used);
	void WriteConnectivity (int ID, ExodusT& exo, const iArrayT& nodes_used);
};

#endif
