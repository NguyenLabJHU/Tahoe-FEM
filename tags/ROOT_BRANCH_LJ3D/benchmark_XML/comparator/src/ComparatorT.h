/* $Id: ComparatorT.h,v 1.13 2004-02-26 08:54:34 paklein Exp $ */

#ifndef _COMPARATOR_T_H_
#define _COMPARATOR_T_H_

/* base class */
#include "FileCrawlerT.h"

/* direct members */
#include "AutoArrayT.h"

namespace Tahoe {

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

	/* prompt input files until "quit" */
	virtual void Run(void);

protected:

	/** handle batch file command */
	virtual void BatchFileCommand(const StringT& command, ifstreamT& batch);

	/* MUST be overloaded */
	virtual void RunJob(ifstreamT& in, ostream& status);
	
	/* batch file processing */
	virtual void RunBatch(ifstreamT& in, ostream& status);

	/** recursive dispatch */
	virtual void JobOrBatch(ifstreamT& in, ostream& status);
	
private:

	/** compare results against benchmarks. The input stream is expected
	 * to be a Tahoe parameters file. The results are expected to be in
	 * the same directory, while the benchmark results are expected to be
	 * in a "benchmark" subdirectory. */
	bool PassOrFail(ifstreamT& in); //const;
	// cannot be const until since tolerances are class data that can change

	/** compare results against benchmarks. The input stream is expected
	 * to be a Tahoe parameters file. The results are expected to be in
	 * the same directory, while the benchmark results are expected to be
	 * in a "benchmark" subdirectory. 
	 * However, the results and benchmark results are expected to have
	 * an extension to the root (supplied in the command line). 
	 * Example: root.geom then becomes root{extenstion}.geom */
	bool PassOrFail_Extension (ifstreamT& in, const StringT& extension);

	/** set the tolerances, read tolerances from file if provided. */
	void Tolerance (ifstreamT& in, bool& do_relative, bool& do_absolute);

	/** compare results */
	bool PassOrFail(const StringT& file_1, const StringT& file_2, 
		bool do_rel, bool do_abs);

	/** deprecated version not using the ModelManagerT class to read data */
	bool PassOrFail_old(const StringT& file_1, const StringT& file_2, 
		bool do_rel, bool do_abs);

	/* read data block header */
	bool ReadDataBlockInfo(ifstreamT& in, double& time, int& num_blocks) const;

	/* read block of nodal data */
	bool ReadNodalData(ifstreamT& in, ArrayT<StringT>& labels, dArray2DT& data) const;

	/* read block of element data */
	bool ReadElementData(ifstreamT& in, ArrayT<StringT>& labels, dArray2DT& data, StringT& block_ID) const;

	/* compare blocks - normalized by set 1 */
	bool CompareDataBlocks(const ArrayT<StringT>& labels_1, const dArray2DT& data_1,
		const ArrayT<StringT>& labels_2, const dArray2DT& data_2,
		bool do_rel, bool do_abs) const;

private:

	/* tolerances */
	double fAbsTol;	
	double fRelTol;

	/** \name tolerances passed in the batch file
	 * These value will be different from -1.0 if tolerances were
	 * passed as command in the batch file. */
	/*@{*/
	double fAbsTol_batch;
	double fRelTol_batch;
	/*@}*/

	/* history */
	bool    fIsRoot;
	AutoArrayT<StringT> fFiles;
	AutoArrayT<bool>    fPassFail;
	
	/* labels to skip */
	AutoArrayT<StringT> fSkipLabels;
};

} //namespace Tahoe

#endif /* _COMPARATOR_T_H_ */
