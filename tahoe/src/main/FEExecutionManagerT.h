/* $Id: FEExecutionManagerT.h,v 1.4 2002-01-07 00:56:22 paklein Exp $ */
/* created: paklein (09/21/1997) */

#ifndef _FE_EXECMAN_T_H_
#define _FE_EXECMAN_T_H_

/* base class */
#include "ExecutionManagerT.h"

/* direct members */
#include "IOBaseT.h"

/* forward declarations */
template <class TYPE> class ArrayT;
class iArrayT;
class OutputSetT;
class IOManager;
class FEManagerT;
class PartitionT;

/** class to handle file driven finite element simulations */
class FEExecutionManagerT: public ExecutionManagerT
{
public:

	/* Constructor */
	FEExecutionManagerT(int argc, char* argv[], char job_char, char batch_char);

	/* Prompt input files until "quit" */
	virtual void Run(void);

protected:

	/** add the command line option to the list. \returns true if the option was
	 * added, false otherwise. */
	virtual bool AddCommandLineOption(const char* str);

	/** overloaded */
	virtual void RunJob(ifstreamT& in, ostream& status);

	/* basic MP support */
	int Rank(void) const;
	int Size(void) const;

private:

	/** standard serial driver */
	void RunJob_serial(ifstreamT& in, ostream& status) const;
	
	/** parallel driver */
	void RunJob_parallel(ifstreamT& in, ostream& status) const;

	/** generate decomposition files */
	void RunDecomp_serial(ifstreamT& in, ostream& status) const;

	/** join parallel results files */
	void RunJoin_serial(ifstreamT& in, ostream& status) const;

	/** print message on exception */
	void Rewind(ifstreamT& in, ostream& status) const;

	/** extract the model file name from the stream */
	void GetModelFile(ifstreamT& in, StringT& model_file,
		IOBaseT::FileTypeT& format) const;

	/** generate decomposition data */
	void Decompose(ifstreamT& in, int size, const StringT& model_file,
		const StringT& global_model_file, IOBaseT::FileTypeT format,
		const StringT& output_map_file) const;

	/** returns true if a new decomposition is needed */
	bool NeedDecomposition(ifstreamT& in, const StringT& model_file,
		int size) const;

	/** returns true if the global output model file is not found */
	bool NeedModelFile(const StringT& model_file, IOBaseT::FileTypeT format) const;

	/** returns true if a new decomposition is needed */
	bool NeedOutputMap(ifstreamT& in, const StringT& map_file,
		int size) const;
	void ReadOutputMap(ifstreamT& in, const StringT& map_file,
		iArrayT& map) const;

	/** set output map based on length of map. The map defines the output prcoessor
	 * for each OutputSetT.
	 * \param output_sets list of OutputSetT's
	 * \param output_map returns with the output processor for each OutputSetT
	 * \param size number of output processors. */
	void SetOutputMap(const ArrayT<OutputSetT*>& output_sets,
		iArrayT& output_map, int size) const;

	/** construct and return the local IOManager */
	IOManager* NewLocalIOManager(const FEManagerT* global_FEman,
		const iArrayT& output_map) const;
		
	/** write partial geometry files */
	void EchoPartialGeometry(const PartitionT& partition, const StringT& model_file,
		const StringT& partial_file, IOBaseT::FileTypeT format) const;
	void EchoPartialGeometry_ExodusII(const PartitionT& partition,
		const StringT& model_file, const StringT& partial_file) const;
	void EchoPartialGeometry_TahoeII(const PartitionT& partition,
		const StringT& model_file, const StringT& partial_file) const;
};

#endif /* _FE_EXECMAN_T_H_ */
