/* $Id: ExodusOutputT.h,v 1.6 2002-07-05 22:26:27 paklein Exp $ */
/* created: sawimme (05/18/1999)                                          */

#ifndef _EXODUSOUTPUT_T_H_
#define _EXODUSOUTPUT_T_H_

/* base class */
#include "OutputBaseT.h"

namespace Tahoe {

/* forward declarations */
class ExodusT;

class ExodusOutputT: public OutputBaseT
{
public:
  /** constructor
   * \param out error stream
   * \param out_strings see OutputBaseT::OutputBaseT */
	ExodusOutputT(ostream& out, const ArrayT<StringT>& out_strings);

	/** write geometry for all output sets */
	virtual void WriteGeometry(void);
	/** write geometry and node variables for output set ID */
	virtual void WriteOutput(double time, int ID, const dArray2DT& n_values,
		const dArray2DT& e_values);

private:

	/** generate database file name for the given ID */
	void FileName(int ID, StringT& filename) const;

	/** create results files */
	void CreateResultsFile(int ID, ExodusT& exo);

	/** creat geometry file */
	void CreateGeometryFile(ExodusT& exo);

	// items common to CreateResultsFile and CreateGeometryFile
	void AssembleQA (ArrayT<StringT>& qa) const; /**< create QA records */
	void WriteCoordinates (ExodusT& exo, const iArrayT& nodes_used); /**< write coordinates */
	void WriteConnectivity (int ID, ExodusT& exo, const iArrayT& nodes_used); /**< write connectivity blocks for this output set */

 private:
	iArrayT fNodeSetIntIDs; /**< integer ID values for node sets */
	iArrayT fSideSetIntIDs; /**< integer ID values for side sets */
};

} // namespace Tahoe 
#endif
