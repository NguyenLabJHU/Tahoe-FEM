/* $Id: TecPlotOutputT.h,v 1.1.1.1.2.1 2001-11-01 19:44:36 sawimme Exp $ */
/* created: sawimme (06/06/2000)                                          */

#ifndef _TECPLOTOUTPUT_T_H_
#define _TECPLOTOUTPUT_T_H_

/* base class */
#include "OutputBaseT.h"

/* forward declarations */
class TecPlotT;

class TecPlotOutputT: public OutputBaseT
{
public:
	TecPlotOutputT(ostream& out, const ArrayT<StringT>& out_strings, int digits);

	/* output functions */
	virtual void WriteGeometry(void);
	virtual void WriteOutput(double time, int ID, const dArray2DT& n_values, const dArray2DT& e_values);

private:

	/* generate database file name for the given ID */
	void FileName(int ID, StringT& filename, int printstep) const;

private:
	bool fBinary;
	int fNumDigits;
};

#endif
