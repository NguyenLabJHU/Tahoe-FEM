/* $Id: FEExecutionManagerT.h,v 1.27.2.4 2004-09-15 02:14:15 d-farrell2 Exp $ */
/* created: paklein (09/21/1997) */
#ifndef _FE_EXECMAN_T_H_
#define _FE_EXECMAN_T_H_

/* base class */
#include "ExecutionManagerT.h"
#include "FEDecomposeT.h" // DEF 4 Aug 04

/* direct members */
#include "IOBaseT.h"

namespace Tahoe {

/* forward declarations */
template <class TYPE> class ArrayT;
class iArrayT;
class OutputSetT;
class IOManager;
class FEManagerT;
class PartitionT;
class ModelManagerT;
class dArray2DT;
class dArray2DT;
class StringT;
class ParameterListT;

/** class to handle file driven finite element simulations */
class FEExecutionManagerT: public ExecutionManagerT
{
public:

	/** constructor */
	FEExecutionManagerT(int argc, char* argv[], char job_char, char batch_char,
		CommunicatorT& comm);

protected:

	/** add the command line option to the list. \returns true if the option was
	 * added, false otherwise. */
	virtual bool AddCommandLineOption(const char* str);

	/** remove the command line option to the list. \returns true if the option was
	 * removed, false otherwise. */
	virtual bool RemoveCommandLineOption(const char* str);

	/** overloaded */
	virtual void RunJob(ifstreamT& in, ostream& status);

	/** Recursive dispatch */
	virtual void JobOrBatch(ifstreamT& in, ostream& status);

private:

	/** \name execution modes */
	/*@{*/
	/** enum for execution modes */
	enum ModeT {
        kJob = 0,
  kDecompose = 1,
       kJoin = 2,
   kBridging = 3,
        kTHK = 4,
        kDTD = 5
	};
	
	/** serial and parallel driver */
	void RunJob_analysis(const StringT& input_file, ostream& status) const;

	/** generate decomposition files */
	void RunDecomp_serial(const StringT& input_file, ostream& status, CommunicatorT& comm, int size = -1) const;

	/** join parallel results files */
	void RunJoin_serial(const StringT& input_file, ostream& status, int size = -1) const;

	/** dump current parameter description file */
	void RunWriteDescription(int doc_type) const;
	/*@}*/

	/** construct and return the local IOManager */
	IOManager* NewLocalIOManager(const FEManagerT* global_FEman,
		const iArrayT& output_map) const;
};

} // namespace Tahoe 
#endif /* _FE_EXECMAN_T_H_ */
