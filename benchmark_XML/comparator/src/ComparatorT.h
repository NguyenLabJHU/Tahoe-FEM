/* $Id: ComparatorT.h,v 1.4 2001-06-14 20:52:09 paklein Exp $ */

#ifndef _COMPARATOR_T_H_
#define _COMPARATOR_T_H_

/* base class */
#include "FileCrawlerT.h"

/* direct members */
#include "AutoArrayT.h"

/* forward declarations */
class dArray2DT;

/* TODO:
 * (1) flexible tolerance specifications
 *     - by point absolute/relative error tolerance
 *     - by (field) norm absolute/relative error tolerance
 * (2) ability to (re-)set tolerances with local tolerance
 *     specifications files.
 * (3) ability to tolerance by variable
 */

class ComparatorT: public FileCrawlerT
{
public:

	/* constructors */
	ComparatorT(int argc, char* argv[], char job_char, char batch_char);

protected:

	/* MUST be overloaded */
	virtual void RunJob(ifstreamT& in, ostream& status);
	
	/* batch file processing */
	virtual void RunBatch(ifstreamT& in, ostream& status);
	
private:

	/* compare results against benchmarks */
	bool PassOrFail(ifstreamT& in); //const;
	// cannot be const until since tolerances are class data that can change

	/* read data block header */
	bool ReadDataBlockInfo(ifstreamT& in, int& group, double& time) const;

	/* read block of nodal data */
	bool ReadNodalData(ifstreamT& in, ArrayT<StringT>& labels, dArray2DT& data) const;

	/* read block of element data */
	bool ReadElementData(ifstreamT& in, ArrayT<StringT>& labels, dArray2DT& data) const;

	/* compare blocks - normalized by set 1 */
	bool CompareDataBlocks(const ArrayT<StringT>& labels_1, const dArray2DT& data_1,
		const ArrayT<StringT>& labels_2, const dArray2DT& data_2) const;

private:

	/* tolerances */
	double fAbsTol;	
	double fRelTol;

	/* history */
	bool    fIsRoot;
	AutoArrayT<StringT> fFiles;
	AutoArrayT<bool>    fPassFail;
};

#endif /* _COMPARATOR_T_H_ */
