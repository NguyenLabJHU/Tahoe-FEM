/* $Id: FEManagerT_mpi.h,v 1.21 2004-07-15 08:31:03 paklein Exp $ */
/* created: paklein (01/12/2000) */
#ifndef _FE_MANAGER_MPI_H_
#define _FE_MANAGER_MPI_H_

/* base class */
#include "FEManagerT.h"

/* direct members */
#include "PartitionT.h"
#include "dArray2DT.h"
#include "ofstreamT.h"

namespace Tahoe {

/* forward declarations */
class IOManager_mpi;
class ParitionT;

class FEManagerT_mpi: public FEManagerT
{
public:

	/** task codes */
	enum TaskT {kDecompose = 0,
	                  kRun = 1};

	/* constructor - partition can be NULL for decomposition */
	FEManagerT_mpi(const StringT& input, ofstreamT& output, CommunicatorT& comm,
		const ArrayT<StringT>& argv, PartitionT* partition, TaskT task);

	/* destructor */
	virtual ~FEManagerT_mpi(void);

	/* system relaxation */
	virtual GlobalT::RelaxCodeT RelaxSystem(int group) const;

	/** update solution */
	virtual void Update(int group, const dArrayT& update);

	/** initiate the process of writing output from all output sets 
	 * \param time time label associated with the output data */
	virtual void WriteOutput(double time);
	
	/** write results for a single output set
	 * \param ID output set ID for the given data
	 * \param n_values nodal output values
	 * \param e_values element output values */
	virtual void WriteOutput(int ID, const dArray2DT& n_values, const dArray2DT& e_values) const;

	/* (temporarily) direct output away from main out */
	virtual void DivertOutput(const StringT& outfile);
	virtual void RestoreOutput(void);

	/* domain decomposition (graph is returned) */
	void Decompose(ArrayT<PartitionT>& partition, GraphT& graph, bool verbose, int method);

	/* return reference to partition data */
	const PartitionT& Partition(void) const;

	/* set the external IOManager */
	void SetExternalIO(IOManager_mpi* externIO);

	/* debugging */
	virtual const iArrayT* ElementMap(const StringT& block_ID) const;

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** \name solution steps 
	 * All steps return ExceptionT::kNoError = 0 unless a problem occurs. */
	/*@{*/
	/** initialize the new time interval */
	virtual ExceptionT::CodeT InitStep(void);

	/** execute the solution procedure */
	virtual ExceptionT::CodeT SolveStep(void);
	/*@}*/

	/* (re-)set system to initial conditions */
	virtual ExceptionT::CodeT InitialCondition(void);

	/** \name equation system information */
	/*@{*/
	virtual int GetGlobalEquationStart(int group, int start_eq_shift) const;
	virtual int GetGlobalNumEquations(int group) const;
	/*@}*/

	/** construct a new CommManagerT. Should be called some time after the
	 * ModelManagerT has been constructed */
	virtual CommManagerT* New_CommManager(void) const;

private:

	/** write time stamp to log file */
	void TimeStamp(const char* message) const;

	/** collect computation effort for each node */
	void WeightNodalCost(iArrayT& weight) const;

	/** \name decomposition methods */
	/*@{*/
	void DoDecompose_1(ArrayT<PartitionT>& partition, GraphT& graph, bool verbose, int method);
	void DoDecompose_2(ArrayT<PartitionT>& partition, GraphT& graph, bool verbose, int method);
	/*@}*/

private:

	/* task */
	TaskT fTask;

	/* external IO */
	IOManager_mpi* fExternIOManager;
	IOBaseT::FileTypeT fInputFormat;
	StringT fModelFile;

	/* partition information */
	PartitionT* fPartition;
	
	/** log file */
	ofstreamT flog;
};

/* return reference to partition data */
inline const PartitionT& FEManagerT_mpi::Partition(void) const
{
	if (!fPartition) throw ExceptionT::kGeneralFail;
	return *fPartition;
}

/* set the external IOManager */
inline void FEManagerT_mpi::SetExternalIO(IOManager_mpi* externIO)
{
	fExternIOManager = externIO;
}

/* debugging */
inline const iArrayT* FEManagerT_mpi::ElementMap(const StringT& block_ID) const
{
	if (fTask == kDecompose)
		return NULL; // assume no map
	else
	{
		if (!fPartition) throw ExceptionT::kGeneralFail;		
		return &(fPartition->ElementMap(block_ID));
	}
}

} // namespace Tahoe 
#endif /* _FE_MANAGER_H_ */
