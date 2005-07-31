/* $Id: AbaqusOutputT.h,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
/* created: sawimme (05/31/2000)                                          */

#ifndef _ABAQUSOUTPUT_T_H_
#define _ABAQUSOUTPUT_T_H_

/* base class */
#include "OutputBaseT.h"
#include "AbaqusT.h"

/* forward declarations */

class AbaqusOutputT: public OutputBaseT
{
public:
	AbaqusOutputT(ostream& out, const ArrayT<StringT>& out_strings, bool binary);

	/* output functions */
	virtual void WriteGeometry(void);
	virtual void WriteOutput(double time, int ID, const dArray2DT& n_values, const dArray2DT& e_values);

private:

	/* generate database file name for the given ID */
	void FileName(int ID, StringT& filename) const;
	bool OpenFile (ofstream& out, const StringT& filename, AbaqusT& aba);
	
	void CreateResultsFile (int ID, AbaqusT& aba, ofstream& out);
	void SetRecordKey (const ArrayT<StringT>& labels, ArrayT<AbaqusT::VariableKeyT>& keys) const;
	
private:
	bool fBinary;
	int fBufferWritten;
	double fOldTime;
};

#endif
