/* $Id: FEManagerT_mpi.h,v 1.9 2002-07-02 19:55:30 cjkimme Exp $ */
/* created: paklein (01/12/2000) */

#ifndef _FE_MANAGER_MPI_H_
#define _FE_MANAGER_MPI_H_

/* base class */
#include "FEManagerT.h"

//TEMP
#include <time.h>
#include "fstreamT.h"

/* direct members */
#include "PartitionT.h"
#include "dArray2DT.h"
#ifdef __MPI__
#include "mpi.h"
#else

namespace Tahoe {

typedef int MPI_Request;
typedef int MPI_Op;

} // namespace Tahoe 
#endif

namespace Tahoe {

/* forward declarations */
class IOManager_mpi;
class ParitionT;

class FEManagerT_mpi: public FEManagerT
{
public:

	/* task codes */
	enum TaskT {kDecompose = 0,
	                  kRun = 1};

	/* constructor - partition can be NULL for decomposition */
	FEManagerT_mpi(ifstreamT& input, ofstreamT& output, const PartitionT* partition,
		TaskT task);

	/* destructor */
	virtual ~FEManagerT_mpi(void);
	
	/* exception handling */
//	virtual void HandleException(int exception);

	/* time sequence messaging */
//	virtual bool Step(void);

	/* solution update */
	virtual void Update(int group, const dArrayT& update);

	/* system relaxation */
	virtual GlobalT::RelaxCodeT RelaxSystem(int group) const;

	/** initiate the process of writing output from all output sets 
	 * \param time time label associated with the output data */
	virtual void WriteOutput(double time, IOBaseT::OutputModeT mode);
	
	/** write results for a single output set
	 * \param ID output set ID for the given data
	 * \param n_values nodal output values
	 * \param e_values element output values */
	virtual void WriteOutput(int ID, const dArray2DT& n_values, const dArray2DT& e_values);

	/** reference to the coordinate array */
	const dArray2DT& Coordinates(void) const;

	/* (temporarily) direct output away from main out */
	virtual void DivertOutput(const StringT& outfile);
	virtual void RestoreOutput(void);

	/* return list of ID's of external nodes */
	virtual void IncomingNodes(iArrayT& nodes_in) const;
	virtual void OutgoingNodes(iArrayT& nodes_out) const;

	/* get external nodal values */
	virtual void SendExternalData(const dArray2DT& all_out_data);
	virtual void RecvExternalData(dArray2DT& external_data);
	virtual void SendRecvExternalData(const iArray2DT& all_out_data, iArray2DT& external_data);

	/* synchronize */
	virtual void Wait(void);

	/* domain decomposition (graph is returned) */
	void Decompose(ArrayT<PartitionT>& partition, GraphT& graph, bool verbose, int method);

	/* return reference to partition data */
	const PartitionT& Partition(void) const;

	/* set the external IOManager */
	void SetExternalIO(IOManager_mpi* externIO);

	/* debugging */
	virtual const iArrayT* NodeMap(void) const;
	virtual const iArrayT* ElementMap(const StringT& block_ID) const;

	/* basic MP support */
	virtual int Rank(void) const;
	virtual int Size(void) const;

protected:

	/** initialize the new time interval */
	virtual int InitStep(void);

	/* initialization functions */
	virtual void ReadParameters(InitCodeT init);
	virtual void SetNodeManager(void);
	virtual void SetElementGroups(void);  	

	/* (re-)set system to initial conditions */
	virtual void InitialCondition(void);

	/* reduce single value */
	int AllReduce(MPI_Op operation, int value);

	/* equation system information */
	virtual int GetGlobalEquationStart(int group) const;
	virtual int GetGlobalNumEquations(int group) const;

private:

	/* allocate buffers for all-to-all communications */
	void AllocateBuffers(int minor_dim, ArrayT<dArray2DT>& recv, ArrayT<dArray2DT>& send);
	void AllocateBuffers(int minor_dim, ArrayT<iArray2DT>& recv, ArrayT<iArray2DT>& send);

	/* write time stamp to log file */
	void TimeStamp(const char* message, bool flush_stream = false) const;

	/* returns the time string */
	const char* WallTime(void) const;

	/* collect computation effort for each node */
	void WeightNodalCost(iArrayT& weight) const;

	/* decomposition methods */
	void DoDecompose_1(ArrayT<PartitionT>& partition, GraphT& graph, bool verbose, int method);
	void DoDecompose_2(ArrayT<PartitionT>& partition, GraphT& graph, bool verbose, int method);

private:

	/* task */
	TaskT fTask;

	/* external IO */
	IOManager_mpi* fExternIOManager;
	IOBaseT::FileTypeT fInputFormat;
	StringT fModelFile;
	
	/* partition information */
	int fRank;
	const PartitionT* fPartition;
	ArrayT<dArray2DT> fRecvBuffer;
	ArrayT<dArray2DT> fSendBuffer;
	ArrayT<MPI_Request> fRecvRequest;
	ArrayT<MPI_Request> fSendRequest;

	//TEMP?
	ofstreamT flog;
	clock_t  flast_time;
};

/* return reference to partition data */
inline const PartitionT& FEManagerT_mpi::Partition(void) const
{
	if (!fPartition) throw eGeneralFail;
	return *fPartition;
}

/* set the external IOManager */
inline void FEManagerT_mpi::SetExternalIO(IOManager_mpi* externIO)
{
	fExternIOManager = externIO;
}

/* debugging */
inline const iArrayT* FEManagerT_mpi::NodeMap(void) const
{
	if (fTask == kDecompose)
		return NULL; // assume no map
	else
	{
		if (!fPartition) throw eGeneralFail;
		return &(fPartition->NodeMap());
	}
}

inline const iArrayT* FEManagerT_mpi::ElementMap(const StringT& block_ID) const
{
	if (fTask == kDecompose)
		return NULL; // assume no map
	else
	{
		if (!fPartition) throw eGeneralFail;		
		return &(fPartition->ElementMap(block_ID));
	}
}

} // namespace Tahoe 
#endif /* _FE_MANAGER_H_ */
