/* $Id: TecPlotOutputT.h,v 1.2.4.1 2002-06-27 18:01:04 cjkimme Exp $ */
/* created: sawimme (06/06/2000)                                          */

#ifndef _TECPLOTOUTPUT_T_H_
#define _TECPLOTOUTPUT_T_H_

/* base class */
#include "OutputBaseT.h"

/* forward declarations */

namespace Tahoe {

class TecPlotT;

class TecPlotOutputT: public OutputBaseT
{
public:
  /** constructor
   * \param out error stream 
   * \param out_strings see OutputBaseT::OutputBaseT
   * \param digits number of digits to use in file name print increment */
	TecPlotOutputT(ostream& out, const ArrayT<StringT>& out_strings, int digits);

	/** write geometry for all output sets */
	virtual void WriteGeometry(void);

	/** write geometry and node variables for output set ID */
	virtual void WriteOutput(double time, int ID, const dArray2DT& n_values, const dArray2DT& e_values);

private:

	/** generate database file name for the given ID */
	void FileName(int ID, StringT& filename, int printstep) const;

private:
	/** not currently used */
	bool fBinary;

	/** number of digits in file name for print increment */
	int fNumDigits;
};

} // namespace Tahoe 
#endif